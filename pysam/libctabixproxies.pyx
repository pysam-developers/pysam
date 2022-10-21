from cpython cimport PyBytes_FromStringAndSize

from libc.stdio cimport printf, feof, fgets
from libc.string cimport strcpy, strlen, memcmp, memcpy, memchr, strstr, strchr
from libc.stdlib cimport free, malloc, calloc, realloc
from libc.stdlib cimport atoi, atol, atof

from pysam.libcutils cimport force_bytes, force_str, charptr_to_str
from pysam.libcutils cimport encode_filename, from_string_and_size

import collections
import copy


cdef char *StrOrEmpty(char * buffer):
     if buffer == NULL:
         return ""
     else: return buffer


cdef int isNew(char * p, char * buffer, size_t nbytes):
    """return True if `p` is located within `buffer` of size
    `nbytes`
    """
    if p == NULL:
        return 0
    
    return not (buffer <= p <= buffer + nbytes)


cdef class TupleProxy:
    '''Proxy class for access to parsed row as a tuple.

    This class represents a table row for fast read-access.

    Access to individual fields is via the [] operator.
    
    Only read-only access is implemented.

    '''

    def __cinit__(self, encoding="ascii"): 
        self.data = NULL
        self.fields = NULL
        self.nbytes = 0
        self.is_modified = 0
        self.nfields = 0
        # start counting at field offset
        self.offset = 0
        self.encoding = encoding

    def __dealloc__(self):
        cdef int x
        if self.is_modified:
            for x from 0 <= x < self.nfields:
                if isNew(self.fields[x], self.data, self.nbytes):
                    free(self.fields[x])
                    self.fields[x] = NULL

        if self.data != NULL:
            free(self.data)
        if self.fields != NULL:
            free(self.fields)

    def __copy__(self):
        if self.is_modified:
            raise NotImplementedError(
                "copying modified tuples is not implemented")
        cdef TupleProxy n = type(self)()
        n.copy(self.data, self.nbytes, reset=True)
        return n

    def compare(self, TupleProxy other):
        '''return -1,0,1, if contents in this are binary
        <,=,> to *other*

        '''
        if self.is_modified or other.is_modified:
            raise NotImplementedError(
                'comparison of modified TupleProxies is not implemented')
        if self.data == other.data:
            return 0

        if self.nbytes < other.nbytes:
            return -1
        elif self.nbytes > other.nbytes:
            return 1
        return memcmp(self.data, other.data, self.nbytes)

    def __richcmp__(self, TupleProxy other, int op):
        if op == 2:  # == operator
            return self.compare(other) == 0
        elif op == 3:  # != operator
            return self.compare(other) != 0
        else:
            err_msg = "op {0} isn't implemented yet".format(op)
            raise NotImplementedError(err_msg)

    cdef take(self, char * buffer, size_t nbytes):
        '''start presenting buffer.

        Take ownership of the pointer.
        '''
        self.data = buffer
        self.nbytes = nbytes
        self.update(buffer, nbytes)

    cdef present(self, char * buffer, size_t nbytes):
        '''start presenting buffer.

        Do not take ownership of the pointer.
        '''
        self.update(buffer, nbytes)

    cdef copy(self, char * buffer, size_t nbytes, bint reset=False):
        '''start presenting buffer of size *nbytes*.

        Buffer is a '\0'-terminated string without the '\n'.

        Take a copy of buffer.
        '''
        # +1 for '\0'
        cdef int s = sizeof(char) *  (nbytes + 1)
        self.data = <char*>malloc(s)
        if self.data == NULL:
            raise ValueError("out of memory in TupleProxy.copy()")
        memcpy(<char*>self.data, buffer, s)

        if reset:
            for x from 0 <= x < nbytes:
                if self.data[x] == b'\0':
                    self.data[x] = b'\t'

        self.update(self.data, nbytes)

    cpdef int getMinFields(self):
        '''return minimum number of fields.'''
        # 1 is not a valid tabix entry, but TupleProxy
        # could be more generic.
        return 1

    cpdef int getMaxFields(self):
        '''return maximum number of fields. Return 
        0 for unknown length.'''
        return 0

    cdef update(self, char * buffer, size_t nbytes):
        '''update internal data.

        *buffer* is a \0 terminated string.

        *nbytes* is the number of bytes in buffer (excluding
        the \0)

        Update starts work in buffer, thus can be used
        to collect any number of fields until nbytes
        is exhausted.

        If max_fields is set, the number of fields is initialized to
        max_fields.

        '''
        cdef char * pos
        cdef char * old_pos
        cdef int field
        cdef int max_fields, min_fields, x

        assert strlen(buffer) == nbytes, \
            "length of buffer (%i) != number of bytes (%i)" % (
            strlen(buffer), nbytes)

        if buffer[nbytes] != 0:
            raise ValueError("incomplete line at %s" % buffer)

        #################################
        # remove line breaks and feeds and update number of bytes
        x = nbytes - 1
        while x > 0 and (buffer[x] == b'\n' or buffer[x] == b'\r'): 
            buffer[x] = b'\0'
            x -= 1
        self.nbytes = x + 1

        #################################
        # clear data
        if self.fields != NULL:
            free(self.fields)
 
        for field from 0 <= field < self.nfields:
            if isNew(self.fields[field], self.data, self.nbytes):
                free(self.fields[field])
                
        self.is_modified = self.nfields = 0

        #################################
        # allocate new
        max_fields = self.getMaxFields()
        # pre-count fields - better would be
        # to guess or dynamically grow
        if max_fields == 0:
            for x from 0 <= x < nbytes:
                if buffer[x] == b'\t':
                    max_fields += 1
            max_fields += 1

        self.fields = <char **>calloc(max_fields, sizeof(char *)) 
        if self.fields == NULL:
            raise ValueError("out of memory in TupleProxy.update()")

        #################################
        # start filling
        field = 0
        self.fields[field] = pos = buffer
        field += 1
        old_pos = pos
        while 1:
            
            pos = <char*>memchr(pos, b'\t', nbytes)
            if pos == NULL:
                break
            if field >= max_fields:
                raise ValueError(
                    "parsing error: more than %i fields in line: %s" %
                    (max_fields, buffer))

            pos[0] = b'\0'
            pos += 1
            self.fields[field] = pos
            field += 1
            nbytes -= pos - old_pos
            if nbytes < 0:
                break
            old_pos = pos
        self.nfields = field
        if self.nfields < self.getMinFields():
            raise ValueError(
                "parsing error: fewer than %i fields in line: %s" %
                (self.getMinFields(), buffer))

    def _getindex(self, int index):
        '''return item at idx index'''
        cdef int i = index
        if i < 0:
            i += self.nfields
        if i < 0:
            raise IndexError("list index out of range")
        # apply offset - separating a fixed number 
        # of fields from a variable number such as in VCF
        i += self.offset
        if i >= self.nfields:
            raise IndexError(
                "list index out of range %i >= %i" %
                (i, self.nfields))
        return force_str(self.fields[i], self.encoding)

    def __getitem__(self, key):
        if type(key) == int:
            return self._getindex(key)
        # slice object
        start, end, step = key.indices(self.nfields)
        result = []
        for index in range(start, end, step):
            result.append(self._getindex(index))
        return result

    def _setindex(self, index, value):
        '''set item at idx index.'''
        cdef int idx = index
        if idx < 0:
            raise IndexError("list index out of range")
        if idx >= self.nfields:
            raise IndexError("list index out of range")

        if isNew(self.fields[idx], self.data, self.nbytes):
            free(self.fields[idx])

        self.is_modified = 1

        if value is None:
            self.fields[idx] = NULL
            return

        # conversion with error checking
        value = force_bytes(value)
        cdef char * tmp = <char*>value
        self.fields[idx] = <char*>malloc((strlen( tmp ) + 1) * sizeof(char))
        if self.fields[idx] == NULL:
            raise ValueError("out of memory" )
        strcpy(self.fields[idx], tmp)

    def __setitem__(self, index, value):
        '''set item at *index* to *value*'''
        cdef int i = index
        if i < 0:
            i += self.nfields
        i += self.offset
        
        self._setindex(i, value)

    def __len__(self):
        return self.nfields

    def __iter__(self):
        return TupleProxyIterator(self)

    def __str__(self):
        '''return original data'''
        # copy and replace \0 bytes with \t characters
        cdef char * cpy
        if self.is_modified:
            # todo: treat NULL values
            result = []
            for x in xrange(0, self.nfields):
                result.append(StrOrEmpty(self.fields[x]).decode(self.encoding))
            return "\t".join(result)
        else:
            cpy = <char*>calloc(sizeof(char), self.nbytes+1)
            if cpy == NULL:
                raise ValueError("out of memory")
            memcpy(cpy, self.data, self.nbytes+1)
            for x from 0 <= x < self.nbytes:
                if cpy[x] == b'\0':
                    cpy[x] = b'\t'
            result = cpy[:self.nbytes]
            free(cpy)
            r = result.decode(self.encoding)
            return r


cdef class TupleProxyIterator:
    def __init__(self, proxy):
        self.proxy = proxy
        self.index = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self.index >= self.proxy.nfields:
            raise StopIteration
        cdef char *retval = self.proxy.fields[self.index]
        self.index += 1
        return force_str(retval, self.proxy.encoding) if retval != NULL else None


def toDot(v):
    '''convert value to '.' if None'''
    if v is None:
        return "." 
    else:
        return str(v)

def quote(v):
    '''return a quoted attribute.'''
    if isinstance(v, str):
        return '"%s"' % v
    else: 
        return str(v)


cdef class NamedTupleProxy(TupleProxy):

    map_key2field = {}

    def __setattr__(self, key, value):
        '''set attribute.'''
        cdef int idx
        idx, f = self.map_key2field[key]
        if self.nfields < idx:
            raise KeyError("field %s not set" % key)
        TupleProxy.__setitem__(self, idx, str(value))

    def __getattr__(self, key):
        cdef int idx
        idx, f = self.map_key2field[key]
        if self.nfields < idx:
            raise KeyError("field %s not set" % key)
        if f == str:
            return force_str(self.fields[idx],
                              self.encoding)
        return f(self.fields[idx])


cdef dot_or_float(v):
    if v == "" or v == b".":
        return None
    else:
        try:
            return int(v)
        except ValueError:
            return float(v)


cdef dot_or_int(v):
    if v == "" or v == b".":
        return None
    else:
        return int(v)


cdef dot_or_str(v):
    if v == "" or v == b".":
        return None
    else:
        return force_str(v)


cdef int from1based(v):
    return atoi(v) - 1


cdef str to1based(int v):
    return str(v + 1)


cdef class GTFProxy(NamedTupleProxy):
    '''Proxy class for access to GTF fields.

    This class represents a GTF entry for fast read-access.
    Write-access has been added as well, though some care must
    be taken. If any of the string fields (contig, source, ...)
    are set, the new value is tied to the lifetime of the
    argument that was supplied.

    The only exception is the attributes field when set from
    a dictionary - this field will manage its own memory.

    '''
    separator = "; "

    # first value is field index, the tuple contains conversion
    # functions for getting (converting internal string representation
    # to pythonic value) and setting (converting pythonic value to
    # interval string representation)
    map_key2field = {
        'contig' : (0, (str, str)),
        'source' : (1, (dot_or_str, str)),
        'feature': (2, (dot_or_str, str)),
        'start' : (3, (from1based, to1based)),
        'end' : (4, (int, int)),
        'score' : (5, (dot_or_float, toDot)),
        'strand' : (6, (dot_or_str, str)),
        'frame' : (7, (dot_or_int, toDot)),
        'attributes': (8, (str, str))}
    
    def __cinit__(self): 
        # automatically calls TupleProxy.__cinit__
        self.attribute_dict = None
        
    cpdef int getMinFields(self):
        '''return minimum number of fields.'''
        return 9

    cpdef int getMaxFields(self):
        '''return max number of fields.'''
        return 9

    def to_dict(self):
        """parse attributes - return as dict

        The dictionary can be modified to update attributes.
        """
        if not self.attribute_dict:
            self.attribute_dict = self.attribute_string2dict(
                self.attributes)
            self.is_modified = True
        return self.attribute_dict

    def as_dict(self):
        """deprecated: use :meth:`to_dict`
        """
        return self.to_dict()

    def from_dict(self, d):
        '''set attributes from a dictionary.'''
        self.attribute_dict = None
        attribute_string = force_bytes(
            self.attribute_dict2string(d),
            self.encoding)
        self._setindex(8, attribute_string)

    def __str__(self):
        cdef char * cpy
        cdef int x

        if self.is_modified:
            return "\t".join( 
                (self.contig, 
                 toDot(self.source), 
                 toDot(self.feature), 
                 str(self.start + 1),
                 str(self.end),
                 toDot(self.score),
                 toDot(self.strand),
                 toDot(self.frame),
                 self.attributes))
        else: 
            return TupleProxy.__str__(self)

    def invert(self, int lcontig):
        '''invert coordinates to negative strand coordinates
        
        This method will only act if the feature is on the
        negative strand.'''

        if self.strand[0] == '-':
            start = min(self.start, self.end)
            end = max(self.start, self.end)
            self.start, self.end = lcontig - end, lcontig - start

    def keys(self):
        '''return a list of attributes defined in this entry.'''
        if not self.attribute_dict:
            self.attribute_dict = self.attribute_string2dict(
                self.attributes)
        return self.attribute_dict.keys()

    def __getitem__(self, key):
        return self.__getattr__(key)

    def setAttribute(self, name, value):
        '''convenience method to set an attribute.
        '''
        if not self.attribute_dict:
            self.attribute_dict = self.attribute_string2dict(
                self.attributes)
        self.attribute_dict[name] = value
        self.is_modified = True

    def attribute_string2dict(self, s):
        return collections.OrderedDict(
            self.attribute_string2iterator(s))
    
    def __cmp__(self, other):
        return (self.contig, self.strand, self.start) < \
            (other.contig, other.strand, other.start)

    # python 3 compatibility
    def __richcmp__(GTFProxy self, GTFProxy other, int op):
        if op == 0:
            return (self.contig, self.strand, self.start) < \
                (other.contig, other.strand, other.start)
        elif op == 1:
            return (self.contig, self.strand, self.start) <= \
                (other.contig, other.strand, other.start)
        elif op == 2:
            return self.compare(other) == 0
        elif op == 3:
            return self.compare(other) != 0
        else:
            err_msg = "op {0} isn't implemented yet".format(op)
            raise NotImplementedError(err_msg)

    def dict2attribute_string(self, d):
        """convert dictionary to attribute string in GTF format.

        """
        aa = []
        for k, v in d.items():
            if isinstance(v, str):
                aa.append('{} "{}"'.format(k, v))
            else:
                aa.append("{} {}".format(k, str(v)))

        return self.separator.join(aa) + ";"

    def attribute_string2iterator(self, s):
        """convert attribute string in GTF format to records
        and iterate over key, value pairs.
        """
        
        # remove comments
        attributes = force_str(s, encoding=self.encoding)

        # separate into fields
        # Fields might contain a ";", for example in ENSEMBL GTF file
        # for mouse, v78:
        # ...; transcript_name "TXNRD2;-001"; ....
        # The current heuristic is to split on a semicolon followed by a
        # space, see also http://mblab.wustl.edu/GTF22.html

        # Remove white space to prevent a last empty field.
        fields = [x.strip() for x in attributes.strip().split("; ")]
        for f in fields:

            # strip semicolon (GTF files without a space after the last semicolon)
            if f.endswith(";"):
                f = f[:-1]

            # split at most once in order to avoid separating
            # multi-word values
            d = [x.strip() for x in f.split(" ", 1)]

            n, v = d[0], d[1]
            if len(d) > 2:
                v = d[1:]

            if v[0] == '"' and v[-1] == '"':
                v = v[1:-1]
            else:
                ## try to convert to a value
                try:
                    v = float(v)
                    v = int(v)
                except ValueError:
                    pass
                except TypeError:
                    pass
                
            yield n, v
       
    def __getattr__(self, key):
        """Generic lookup of attribute from GFF/GTF attributes 
        """

        # Only called if there *isn't* an attribute with this name
        cdef int idx
        idx, f = self.map_key2field.get(key, (-1, None))
        if idx >= 0:
            # deal with known attributes (fields 0-8)
            if idx == 8:
                # flush attributes if requested
                if self.is_modified and self.attribute_dict is not None:
                    s = self.dict2attribute_string(self.attribute_dict)
                    TupleProxy._setindex(self, idx, s)
                    self.attribute_dict = None
                    return s
                                         
            if f[0] == str:
                return force_str(self.fields[idx],
                                 self.encoding)
            else:
                return f[0](self.fields[idx])
        else:
            # deal with generic attributes (gene_id, ...)
            if self.attribute_dict is None:
                self.attribute_dict = self.attribute_string2dict(
                    self.attributes)
            return self.attribute_dict[key]

    def __setattr__(self, key, value):
        '''set attribute.'''

        # Note that __setattr__ is called before properties, so __setattr__ and
        # properties don't mix well. This is different from __getattr__ which is
        # called after any properties have been resolved.
        cdef int idx
        idx, f = self.map_key2field.get(key, (-1, None))

        if idx >= 0:
            if value is None:
                s = "."
            elif f[1] == str:
                s = force_bytes(value,
                                self.encoding)
            else:
                s = str(f[1](value))
            TupleProxy._setindex(self, idx, s)
        else:
            if self.attribute_dict is None:
                self.attribute_dict = self.attribute_string2dict(
                    self.attributes)
            self.attribute_dict[key] = value
            self.is_modified = True

    # for backwards compatibility
    def asDict(self, *args, **kwargs):
        return self.to_dict(*args, **kwargs)

    def fromDict(self, *args, **kwargs):
        return self.from_dict(*args, **kwargs)
    

cdef class GFF3Proxy(GTFProxy):

    def dict2attribute_string(self, d):
        """convert dictionary to attribute string."""
        return ";".join(["{}={}".format(k, v) for k, v in d.items()])
        
    def attribute_string2iterator(self, s):
        """convert attribute string in GFF3 format to records
        and iterate over key, value pairs.
        """
        
        for f in (x.strip() for x in s.split(";")):
            if not f:
                continue

            key, value = f.split("=", 1)
            value = value.strip()
            
            ## try to convert to a value
            try:
                value = float(value)
                value = int(value)
            except ValueError:
                pass
            except TypeError:
                pass
                
            yield key.strip(), value
   

cdef class BedProxy(NamedTupleProxy):
    '''Proxy class for access to Bed fields.

    This class represents a BED entry for fast read-access.
    '''
    map_key2field = {
        'contig' : (0, str),
        'start' : (1, int),
        'end' : (2, int),
        'name' : (3, str),
        'score' : (4, float),
        'strand' : (5, str),
        'thickStart' : (6, int),
        'thickEnd' : (7, int),
        'itemRGB' : (8, str),
        'blockCount': (9, int),
        'blockSizes': (10, str),
        'blockStarts': (11, str), } 

    cpdef int getMinFields(self):
        '''return minimum number of fields.'''
        return 3

    cpdef int getMaxFields(self):
        '''return max number of fields.'''
        return 12

    cdef update(self, char * buffer, size_t nbytes):
        '''update internal data.

        nbytes does not include the terminal '\0'.
        '''
        TupleProxy.update(self, buffer, nbytes)

        if self.nfields < 3:
            raise ValueError(
                "bed format requires at least three columns")

        # determines bed format
        self.bedfields = self.nfields

        # do automatic conversion
        self.contig = self.fields[0]
        self.start = atoi(self.fields[1]) 
        self.end = atoi(self.fields[2])

    # __setattr__ in base class seems to take precedence
    # hence implement setters in __setattr__
    #property start:
    #    def __get__( self ): return self.start
    #property end:
    #    def __get__( self ): return self.end

    def __str__(self):

        cdef int save_fields = self.nfields
        # ensure fields to use correct format
        self.nfields = self.bedfields
        retval = TupleProxy.__str__(self)
        self.nfields = save_fields
        return retval

    def __setattr__(self, key, value):
        '''set attribute.'''
        if key == "start":
            self.start = value
        elif key == "end":
            self.end = value

        cdef int idx
        idx, f = self.map_key2field[key]
        TupleProxy._setindex(self, idx, str(value))
        

cdef class VCFProxy(NamedTupleProxy):
    '''Proxy class for access to VCF fields.

    The genotypes are accessed via a numeric index.
    Sample headers are not available.
    '''
    map_key2field = { 
        'contig' : (0, str),
        'pos' : (1, int),
        'id' : (2, str),
        'ref' : (3, str),
        'alt' : (4, str),
        'qual' : (5, str),
        'filter' : (6, str),
        'info' : (7, str),
        'format' : (8, str) }

    def __cinit__(self): 
        # automatically calls TupleProxy.__cinit__
        # start indexed access at genotypes
        self.offset = 9

    cdef update(self, char * buffer, size_t nbytes):
        '''update internal data.
        
        nbytes does not include the terminal '\0'.
        '''
        TupleProxy.update(self, buffer, nbytes)

        self.contig = self.fields[0]
        # vcf counts from 1 - correct here
        self.pos = atoi(self.fields[1]) - 1
                             
    def __len__(self):
        '''return number of genotype fields.'''
        return max(0, self.nfields - 9)

    property pos:
       '''feature end (in 0-based open/closed coordinates).'''
       def __get__(self): 
           return self.pos

    def __setattr__(self, key, value):
        '''set attribute.'''
        if key == "pos": 
            self.pos = value
            value += 1

        cdef int idx
        idx, f = self.map_key2field[key]
        TupleProxy._setindex(self, idx, str(value))


__all__ = [
    "TupleProxy",
    "NamedTupleProxy",
    "GTFProxy",
    "GFF3Proxy",
    "BedProxy",
    "VCFProxy"]
