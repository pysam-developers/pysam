import types
from cpython cimport PyString_FromStringAndSize, PyString_AsString, PyString_AS_STRING

cdef char * nextItem( char * buffer ):
    cdef char * pos
    pos = strchr( buffer, '\t' )
    if pos == NULL: raise ValueError( "malformatted entry at %s" % buffer )
    pos[0] = '\0'
    pos += 1
    return pos

cdef char *StrOrEmpty( char * buffer ):
     if buffer == NULL: return ""
     else: return buffer

cdef int isNew( char * p, char * buffer, size_t nbytes ):
     if p == NULL: return 0
     return not (buffer <= p < buffer + nbytes )

cdef class TupleProxy:
    '''Proxy class for access to parsed row as a tuple.

    This class represents a table row for fast read-access.

    Access to individual fields is via the [] operator.
    
    Only read-only access is implemented.
    '''

    def __cinit__(self ): 
        self.data = NULL
        self.fields = NULL
        self.index = 0
        self.nbytes = 0
        self.is_modified = 0
        self.nfields = 0
        # start counting at field offset
        self.offset = 0

    def __dealloc__(self):
        cdef int x
        if self.is_modified:
            for x from 0 <= x < self.nfields:
                if isNew( self.fields[x], self.data, self.nbytes ):
                    free( self.fields[x] )
                    self.fields[x] = NULL

        if self.data != NULL: free(self.data)
        if self.fields != NULL: free( self.fields )

    cdef take( self, char * buffer, size_t nbytes ):
        '''start presenting buffer.

        Take ownership of the pointer.
        '''
        self.data = buffer
        self.nbytes = nbytes
        self.update( buffer, nbytes )

    cdef present( self, char * buffer, size_t nbytes ):
        '''start presenting buffer.

        Do not take ownership of the pointer.
        '''
        self.update( buffer, nbytes )

    cdef copy( self, char * buffer, size_t nbytes ):
        '''start presenting buffer.

        Take a copy of buffer.
        '''
        cdef int s
        # +1 for '\0'
        s = sizeof(char) *  (nbytes + 1)
        self.data = <char*>malloc( s ) 
        if self.data == NULL:
            raise ValueError("out of memory" )
        self.nbytes = nbytes
        memcpy( <char*>self.data, buffer, s )
        self.update( self.data, nbytes )

    cdef int getMaxFields( self, size_t nbytes ):
        '''initialize fields.'''
        return nbytes / 2

    cdef update( self, char * buffer, size_t nbytes ):
        '''update internal data.

        Update starts work in buffer, thus can be used
        to collect any number of fields until nbytes
        is exhausted.

        If max_fields is set, the number of fields is initialized to max_fields.

        '''
        cdef char * pos
        cdef char * old_pos
        cdef int field
        cdef int max_fields, x

        if buffer[nbytes] != 0:
            raise ValueError( "incomplete line at %s" % buffer )

        #################################
        # clear data
        if self.fields != NULL: free(self.fields)
        
        for field from 0 <= field < self.nfields:
            if isNew( self.fields[field], self.data, self.nbytes ):
                free( self.fields[field] )
                
        self.is_modified = self.nfields = 0

        #################################
        # allocate new
        max_fields = self.getMaxFields( nbytes )
        self.fields = <char **>calloc( max_fields, sizeof(char *) ) 
        if self.fields == NULL:
            raise ValueError("out of memory" )

        #################################
        # start filling
        field = 0
        self.fields[field] = pos = buffer
        field += 1
        old_pos = pos
        
        while 1:

            pos = <char*>memchr( pos, '\t', nbytes )
            if pos == NULL: break
            pos[0] = '\0'
            pos += 1
            self.fields[field] = pos
            field += 1
            if field >= max_fields:
                raise ValueError("row too large - more than %i fields" % max_fields )
            nbytes -= pos - old_pos
            if nbytes < 0: break
            old_pos = pos

        self.nfields = field

    def _getindex( self, int index ):
        '''return item at idx index'''
        cdef int i = index
        if i < 0: i += self.nfields
        if i < 0: raise IndexError( "list index out of range" )
        i += self.offset
        if i >= self.nfields:
            raise IndexError( "list index out of range %i >= %i" % (i, self.nfields ))
        return self.fields[i]

    def __getitem__( self, key ):
        if type(key) == int: return self._getindex( key )
        # slice object
        start, end, step = key.indices( self.nfields )
        result = []
        for index in range( start, end, step ):
            result.append( self._getindex( index ) )
        return result

    def _setindex( self, index, value ):
        '''set item at idx index.'''
        cdef int idx = index
        if idx < 0: raise IndexError( "list index out of range" )        
        if idx >= self.nfields:
            raise IndexError( "list index out of range" )

        if isNew( self.fields[idx], self.data, self.nbytes ):
            free( self.fields[idx] )

        self.is_modified = 1

        if value == None:
            self.fields[idx] = NULL
            return

        # conversion with error checking
        cdef char * tmp = PyString_AsString( value )
        self.fields[idx] = <char*>malloc( (strlen( tmp ) + 1) * sizeof(char) )
        if self.fields[idx] == NULL:
            raise ValueError("out of memory" )
        strcpy( self.fields[idx], tmp )

    def __setitem__(self, index, value ):
        '''set item at *index* to *value*'''
        cdef int i = index
        if i < 0: i += self.nfields
        i += self.offset
        
        self._setindex( i, value )

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
        if retval == NULL: return None
        else: return retval

    def __str__(self):
        '''return original data'''
        # copy and replace \0 bytes with \t characters
        if self.is_modified:
            # todo: treat NULL values
            return "\t".join( [StrOrEmpty( self.fields[x]) for x in xrange(0, self.nfields ) ] )
        else:
            cpy = <char*>calloc( sizeof(char), self.nbytes+1 )
            if cpy == NULL:
                raise ValueError("out of memory" )
            memcpy( cpy, self.data, self.nbytes+1)
            for x from 0 <= x < self.nbytes:
                if cpy[x] == '\0': cpy[x] = '\t'
            result = PyString_FromStringAndSize(cpy, self.nbytes)
            free(cpy)
            return result

def toDot( v ):
    '''convert value to '.' if None'''
    if v == None: return "." 
    else: return str(v)

def quote( v ):
    '''return a quoted attribute.'''
    if type(v) in types.StringTypes:
        return '"%s"' % v
    else: 
        return str(v)

cdef class GTFProxy( TupleProxy ):
    '''Proxy class for access to GTF fields.

    This class represents a GTF entry for fast read-access.
    Write-access has been added as well, though some care must
    be taken. If any of the string fields (contig, source, ...)
    are set, the new value is tied to the lifetime of the
    argument that was supplied.

    The only exception is the attributes field when set from
    a dictionary - this field will manage its own memory.

    '''

    def __cinit__(self ): 
        # automatically calls TupleProxy.__cinit__
        self.hasOwnAttributes = False

    def __dealloc__(self):
        # automatically calls TupleProxy.__dealloc__
        if self.hasOwnAttributes:
            free(self.attributes)

    cdef int getMaxFields( self, size_t nbytes ):
        '''return max number of fields.'''
        return 9

    cdef update( self, char * buffer, size_t nbytes ):
        '''update internal data.

        nbytes does not include the terminal '\0'.
        '''
        cdef int end
        cdef char * cstart, * cend, * cscore
        self.contig = buffer
        cdef char * pos

        if buffer[nbytes] != 0:
            raise ValueError( "incomplete line at %s" % buffer )
        
        self.source = pos = nextItem( buffer )
        self.feature = pos = nextItem( pos )
        cstart = pos = nextItem( pos )
        cend = pos = nextItem( pos )
        self.score = pos = nextItem( pos )
        self.strand = pos = nextItem( pos )
        self.frame = pos = nextItem( pos )
        self.attributes = pos = nextItem( pos )

        self.start = atoi( cstart ) - 1
        self.end = atoi( cend )
        self.nfields = 9
       
    property contig:
       '''contig of feature.'''
       def __get__( self ): return self.contig
       def __set__( self, value ): 
           self.is_modified = True
           self.contig = value

    property feature:
       '''feature name.'''
       def __get__( self ): return self.feature
       def __set__( self, value ): 
           self.is_modified = True
           self.feature = value

    property source:
       '''feature source.'''
       def __get__( self ): return self.source
       def __set__( self, value ): 
           self.is_modified = True
           self.source = value

    property start:
       '''feature start (in 0-based open/closed coordinates).'''
       def __get__( self ): return self.start
       def __set__( self, value ): 
           self.is_modified = True
           self.start = value

    property end:
       '''feature end (in 0-based open/closed coordinates).'''
       def __get__( self ): return self.end
       def __set__( self, value ): 
           self.is_modified = True
           self.end = value

    property score:
       '''feature score.'''
       def __get__( self ): 
           if self.score[0] == '.' and self.score[1] == '\0' :
               return None
           else:
               return atof(self.score)
       def __set__( self, value ): 
           self.is_modified = True
           self.score = value

    property strand:
       '''feature strand.'''
       def __get__( self ): return self.strand
       def __set__( self, value ): 
           self.is_modified = True
           self.strand = value

    property frame:
       '''feature frame.'''
       def __get__( self ): return self.frame
       def __set__( self, value ): 
           self.is_modified = True
           self.frame = value

    property attributes:
       '''feature attributes (as a string).'''
       def __get__( self ): return self.attributes
       def __set__( self, value ): 
           self.is_modified = True
           self.attributes = value

    def asDict( self ):
        """parse attributes - return as dict
        """

        # remove comments
        attributes = self.attributes

        # separate into fields
        fields = [ x.strip() for x in attributes.split(";")[:-1]]
        
        result = {}

        for f in fields:
            
            d = [ x.strip() for x in f.split(" ")]
            
            n,v = d[0], d[1]
            if len(d) > 2: v = d[1:]

            if v[0] == '"' and v[-1] == '"':
                v = v[1:-1]
            else:
                ## try to convert to a value
                try:
                    v = float( v )
                    v = int( v )
                except ValueError:
                    pass
                except TypeError:
                    pass

            result[n] = v
        
        return result
    
    def fromDict( self, d ):
        '''set attributes from a dictionary.'''
        cdef char * p
        cdef int l

        # clean up if this field is set twice
        if self.hasOwnAttributes:
            free(self.attributes)

        aa = []
        for k,v in d.items():
            if type(v) == types.StringType:
                aa.append( '%s "%s"' % (k,v) )
            else:
                aa.append( '%s %s' % (k,str(v)) )

        a = "; ".join( aa ) + ";"
        p = a
        l = len(a)
        self.attributes = <char *>calloc( l + 1, sizeof(char) )
        if self.attributes == NULL:
            raise ValueError("out of memory" )
        memcpy( self.attributes, p, l )

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
                 self.strand,
                 self.frame,
                 self.attributes ) )
        else: 
            return TupleProxy.__str__(self)

    def invert( self, int lcontig ):
        '''invert coordinates to negative strand coordinates
        
        This method will only act if the feature is on the
        negative strand.'''

        if self.strand[0] == '-':
            start = min(self.start, self.end)
            end = max(self.start, self.end)
            self.start, self.end = lcontig - end, lcontig - start

    def keys( self ):
        '''return a list of attributes defined in this entry.'''
        r = self.attributes
        return [ x.strip().split(" ")[0] for x in r.split(";") if x.strip() != '' ]

    def __getitem__(self, key):
        return self.__getattr__( key )

    def __getattr__(self, item ):
        """Generic lookup of attribute from GFF/GTF attributes 
        Only called if there *isn't* an attribute with this name
        """
        cdef char * start
        cdef char * query 
        cdef char * cpy
        cdef char * end
        cdef int l
        query = item
        
        start = strstr( self.attributes, query)
        if start == NULL:
            raise AttributeError("'GTFProxy' has no attribute '%s'" % item )

        start += strlen(query) + 1
        # skip gaps before
        while start[0] == ' ': start += 1
        if start[0] == '"':
            start += 1
            end = start
            while end[0] != '\0' and end[0] != '"': end += 1
            l = end - start
            result = PyString_FromStringAndSize( start, l )
            return result
        else:
            return start

    def setAttribute( self, name, value ):
        '''convenience method to set an attribute.'''
        r = self.asDict()
        r[name] = value
        self.fromDict( r )

cdef class NamedTupleProxy( TupleProxy ):

    map_key2field = {}

    def __setattr__(self, key, value ):
        '''set attribute.'''
        cdef int idx
        idx, f = self.map_key2field[key]
        if self.nfields < idx:
            raise KeyError( "field %s not set" % key )
        TupleProxy.__setitem__(self, idx, str(value) )

    def __getattr__(self, key ):
        cdef int idx
        idx, f = self.map_key2field[key]
        if self.nfields < idx:
            raise KeyError( "field %s not set" % key )
        return f( self.fields[idx] )

cdef class BedProxy( NamedTupleProxy ):
    '''Proxy class for access to Bed fields.

    This class represents a GTF entry for fast read-access.
    '''
    map_key2field = { 
        'contig' : (0, str),
        'start' : (1, int),
        'end' : (2, int),
        'name' : (3, str),
        'score' : (4, float),
        'strand' : (5, str),
        'thickStart' : (6,int ),
        'thickEnd' : (7,int),
        'itemRGB' : (8,str),
        'blockCount': (9,int),
        'blockSizes': (10,str),
        'blockStarts': (11,str), } 

    cdef int getMaxFields( self, size_t nbytes ):
        '''return max number of fields.'''
        return 12

    cdef update( self, char * buffer, size_t nbytes ):
        '''update internal data.

        nbytes does not include the terminal '\0'.
        '''
        TupleProxy.update( self, buffer, nbytes )

        if self.nfields < 3:
            raise ValueError( "bed format requires at least three columns" )

        # determines bed format
        self.bedfields = self.nfields

        # do automatic conversion
        self.contig = self.fields[0]
        self.start = atoi( self.fields[1] ) 
        self.end = atoi( self.fields[2] )

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
        retval = TupleProxy.__str__( self )
        self.nfields = save_fields
        return retval

    def __setattr__(self, key, value ):
        '''set attribute.'''
        if key == "start": self.start = value
        elif key == "end": self.end = value

        cdef int idx
        idx, f = self.map_key2field[key]
        TupleProxy._setindex(self, idx, str(value) )

cdef class VCFProxy( NamedTupleProxy ):
    '''Proxy class for access to VCF fields.

    The genotypes are accessed via index.
    '''
    map_key2field = { 
        'contig' : (0, str),
        'pos' : (1, int),
        'id' : (2, str),
        'ref' : (3, str),
        'alt' : (4, str),
        'qual' : (5, str),
        'filter' : (6,str),
        'info' : (7,str),
        'format' : (8,str) }

    def __cinit__(self ): 
        # automatically calls TupleProxy.__cinit__
        # start indexed access at genotypes
        self.offset = 9

    cdef update( self, char * buffer, size_t nbytes ):
        '''update internal data.
        
        nbytes does not include the terminal '\0'.
        '''
        TupleProxy.update( self, buffer, nbytes )

        self.contig = self.fields[0]
        # vcf counts from 1 - correct here
        self.pos = atoi( self.fields[1] ) - 1
                             
    def __len__(self):
        return max(0, self.nfields - 9)

    property pos:
       '''feature end (in 0-based open/closed coordinates).'''
       def __get__( self ): 
           return self.pos

    def __setattr__(self, key, value ):
        '''set attribute.'''
        if key == "pos": 
            self.pos = value
            value += 1

        cdef int idx
        idx, f = self.map_key2field[key]
        TupleProxy._setindex(self, idx, str(value) )

