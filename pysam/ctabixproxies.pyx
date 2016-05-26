from cpython cimport PyBytes_FromStringAndSize

from libc.stdio cimport printf, feof, fgets
from libc.string cimport strcpy, strlen, memcmp, memcpy, memchr, strstr, strchr
from libc.stdlib cimport free, malloc, calloc, realloc
from libc.stdlib cimport atoi, atol, atof

from pysam.cutils cimport force_bytes, force_str, charptr_to_str
from pysam.cutils cimport encode_filename, from_string_and_size

import collections

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
    return not (buffer <= p < buffer + nbytes)


cdef class TupleProxy:
    '''Proxy class for access to parsed row as a tuple.

    This class represents a table row for fast read-access.

    Access to individual fields is via the [] operator.
    
    Only read-only access is implemented.

    '''

    def __cinit__(self, encoding="ascii"): 
        self.data = NULL
        self.fields = NULL
        self.index = 0
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
                if self.data[x] == '\0':
                    self.data[x] = '\t'

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
        while x > 0 and (buffer[x] == '\n' or buffer[x] == '\r'): 
            buffer[x] = '\0'
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
                if buffer[x] == '\t':
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
            
            pos = <char*>memchr(pos, '\t', nbytes)
            if pos == NULL:
                break
            if field >= max_fields:
                raise ValueError(
                    "parsing error: more than %i fields in line: %s" %
                    (max_fields, buffer))

            pos[0] = '\0'
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
                "parsing error: fewer that %i fields in line: %s" %
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
            free(self.fields[idx] )

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
        self.index = 0
        return self

    def __next__(self): 
        """python version of next().
        """
        if self.index >= self.nfields:
            raise StopIteration
        cdef char * retval = self.fields[self.index]
        self.index += 1
        if retval == NULL:
            return None
        else:
            return force_str(retval, self.encoding)

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
                if cpy[x] == '\0':
                    cpy[x] = '\t'
            result = cpy[:self.nbytes]
            free(cpy)
            r = result.decode(self.encoding)
            return r

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


cdef class GTFProxy(TupleProxy):
    '''Proxy class for access to GTF fields.

    This class represents a GTF entry for fast read-access.
    Write-access has been added as well, though some care must
    be taken. If any of the string fields (contig, source, ...)
    are set, the new value is tied to the lifetime of the
    argument that was supplied.

    The only exception is the attributes field when set from
    a dictionary - this field will manage its own memory.
    '''

    def __cinit__(self): 
        # automatically calls TupleProxy.__cinit__
        self.hasOwnAttributes = False
        self._attributes = NULL

    def __dealloc__(self):
        # automatically calls TupleProxy.__dealloc__
        if self.hasOwnAttributes:
            free(self._attributes)

    cpdef int getMinFields(self):
        '''return minimum number of fields.'''
        return 9

    cpdef int getMaxFields(self):
        '''return max number of fields.'''
        return 9

    property contig:
        '''contig of feature.'''
        def __get__(self):
            return self._getindex(0)
        def __set__(self, value):
            self._setindex(0, value)

    property source:
        '''feature source.'''
        def __get__(self):
            return self._getindex(1)
        def __set__(self, value):
            if value is None:
                value = "."
            self._setindex(1, value)

    property feature:
        '''feature name.'''
        def __get__(self):
            return self._getindex(2)
        def __set__(self, value):
            if value is None:
                value = "."
            self._setindex(2, value)

    property start:
        '''feature start (in 0-based open/closed coordinates).'''
        def __get__(self ):
            return int( self._getindex(3)) - 1
        def __set__(self, value ):
            self._setindex(3, str(value+1))

    property end:
        '''feature end (in 0-based open/closed coordinates).'''
        def __get__(self):
            return int(self._getindex(4))
        def __set__(self, value):
            self._setindex(4, str(value))

    property score:
        '''feature score.'''
        def __get__(self): 
            v = self._getindex(5)
            if v == "" or v[0] == '.':
                return None
            else:
                return float(v)

        def __set__(self, value):
            if value is None:
                value = "."
            self._setindex(5, str(value))

    property strand:
        '''feature strand.'''
        def __get__(self):
           return self._getindex(6)
        def __set__(self, value ):
            if value is None:
                value = "."
            self._setindex(6, value)

    property frame:
       '''feature frame.'''
       def __get__(self):
            v = self._getindex(7)
            if v == "" or v[0] == '.':
                return v
            else:
                return int(v)

       def __set__(self, value):
            if value is None:
                value = "."
            self._setindex(7, str(value))

    property attributes:
        '''feature attributes (as a string).'''
        def __get__(self): 
            if self.hasOwnAttributes:
                return force_str(self._attributes)
            else:
                return force_str(self._getindex(8))
        def __set__( self, value): 
            if self.hasOwnAttributes:
                free(self._attributes)
                self._attributes = NULL
                self.hasOwnAttributes = False
            self._setindex(8, value)

    cdef char * getAttributes(self):
        '''return pointer to attributes.'''
        cdef char * attributes
        if self.hasOwnAttributes:
            attributes = self._attributes
        else:
            attributes = self.fields[8]
        if attributes == NULL:
            raise KeyError("no attributes defined GTF entry")
        return attributes

    def asDict(self):
        """parse attributes - return as dict
        """

        # remove comments
        attributes = self.attributes

        # separate into fields
        # Fields might contain a ";", for example in ENSEMBL GTF file
        # for mouse, v78:
        # ...; transcript_name "TXNRD2;-001"; ....
        # The current heuristic is to split on a semicolon followed by a
        # space, see also http://mblab.wustl.edu/GTF22.html

        # Remove white space to prevent a last empty field.
        fields = [x.strip() for x in attributes.strip().split("; ")]
        
        result = collections.OrderedDict()

        for f in fields:

            # strip semicolon (GTF files without a space after the last semicolon)
            if f.endswith(";"):
                f = f[:-1]

            # split at most once in order to avoid separating
            # multi-word values
            d = [x.strip() for x in f.split(" ", 1)]

            n,v = d[0], d[1]
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

            result[n] = v
        
        return result
    
    def fromDict(self, d):
        '''set attributes from a dictionary.'''
        cdef char * p
        cdef int l

        # clean up if this field is set twice
        if self.hasOwnAttributes: 
            free(self._attributes)

        aa = []
        for k,v in d.items():
            if isinstance(v, str):
                aa.append( '%s "%s"' % (k,v) )
            else:
                aa.append( '%s %s' % (k,str(v)) )

        a = force_bytes("; ".join(aa) + ";")
        p = a
        l = len(a)
        self._attributes = <char *>calloc(l + 1, sizeof(char))
        if self._attributes == NULL:
            raise ValueError("out of memory")
        memcpy(self._attributes, p, l)

        self.hasOwnAttributes = True
        self.is_modified = True

    def __str__(self):
        cdef char * cpy
        cdef int x

        if self.is_modified:
            return "\t".join( 
                (self.contig, 
                 self.source, 
                 self.feature, 
                 str(self.start+1),
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
        r = self.attributes
        return [x.strip().split(" ")[0]
                # separator is ';' followed by space
                for x in r.split("; ") if x.strip() != '']

    def __getitem__(self, key):
        return self.__getattr__(key)

    def __getattr__(self, item):
        """Generic lookup of attribute from GFF/GTF attributes 
        Only called if there *isn't* an attribute with this name
        """
        cdef char * start
        cdef char * query
        cdef char * cpy
        cdef char * end
        cdef int l

        #
        # important to use the getAttributes function.
        # Using the self.attributes property to access
        # the attributes caused a hard-to-trace bug
        # in which fields in the attribute string were
        # set to 0.
        # Running through valgrind complained that
        # memory was accessed in the memory field
        # that has been released. It is not clear
        # why this happened and might be a cython bug
        # (Version 0.16). The valgrind warnings
        # disappeard after accessing the C data structures
        # directly and so did the bug.
        cdef char * attributes = self.getAttributes()
        if attributes == NULL:
            raise KeyError("key %s not found, no attributes" % item)

        # add space in order to make sure
        # to not pick up a field that is a prefix of another field
        r = force_bytes(item + " ")
        query = r
        start = strstr(attributes, query)

        if start == NULL:
            raise AttributeError("'GTFProxy' has no attribute '%s'" % item)

        start += strlen(query)
        # skip gaps before
        while start[0] == ' ':
            start += 1

        if start[0] == '"':
            start += 1
            end = start
            while end[0] != '\0' and end[0] != '"':
                end += 1
            l = end - start
            result = force_str(PyBytes_FromStringAndSize(start, l),
                                self.encoding)
            return result
        else:
            return force_str(start, self.encoding)

    def setAttribute(self, name, value):
        '''convenience method to set an attribute.'''
        r = self.asDict()
        r[name] = value
        self.fromDict(r)

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

    def __setattr__(self, key, value ):
        '''set attribute.'''
        if key == "start":
            self.start = value
        elif key == "end":
            self.end = value

        cdef int idx
        idx, f = self.map_key2field[key]
        TupleProxy._setindex(self, idx, str(value) )

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

