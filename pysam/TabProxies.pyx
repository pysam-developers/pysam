import types
from cpython cimport PyString_FromStringAndSize, PyString_AS_STRING

cdef char * nextItem( char * buffer ):
    cdef char * pos
    pos = strchr( buffer, '\t' )
    if pos == NULL: raise ValueError( "malformatted entry at %s" % buffer )
    pos[0] = '\0'
    pos += 1
    return pos

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

    def __dealloc__(self):
        if self.data != NULL:
            free(self.data)

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
        self.nbytes = nbytes
        memcpy( <char*>self.data, buffer, s )
        self.update( self.data, nbytes )

    cdef update( self, char * buffer, size_t nbytes ):
        '''update internal data.'''
        cdef char * pos
        cdef char * old_pos
        cdef int field
        cdef int max_fields
        field = 0

        if buffer[nbytes] != 0:
            raise ValueError( "incomplete line at %s" % buffer )
        
        if self.fields != NULL:
            free(self.fields)
        
        max_fields = nbytes / 4
        self.fields = <char **>calloc( max_fields, sizeof(char *) ) 
        
        pos = buffer
        self.fields[0] = pos
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

    def __getitem__( self, key ):

        cdef int i
        i = key
        if i < 0: i += self.nfields
        if i >= self.nfields or i < 0:
            raise IndexError( "list index out of range" )
        return self.fields[i]

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
        self.index += 1
        return self.fields[self.index-1]

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
        self.isModified = False
        self.hasOwnAttributes = False

    def __dealloc__(self):
        # automatically calls TupleProxy.__dealloc__
        if self.hasOwnAttributes:
            free(self.attributes)

    cdef take( self, char * buffer, size_t nbytes ):
        '''start presenting buffer.

        Take ownership of the pointer.
        '''
        self.isModified = False
        TupleProxy.take( self, buffer, nbytes )

    cdef present( self, char * buffer, size_t nbytes ):
        '''start presenting buffer.

        Do not take ownership of the pointer.
        '''
        self.isModified = False
        TupleProxy.present( self, buffer, nbytes )

    cdef copy( self, char * buffer, size_t nbytes ):
        '''start presenting buffer.

        Take a copy of buffer.
        '''
        self.isModified = False
        TupleProxy.copy( self, buffer, nbytes )

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
        
        pos = strchr( buffer, '\t' )
        if pos == NULL: raise ValueError( "malformatted entry at %s" % buffer )
        pos[0] = '\0'
        pos += 1
        self.source = pos

        pos = strchr( pos, '\t' )
        if pos == NULL: raise ValueError( "malformatted entry at %s" % buffer )
        pos[0] = '\0'
        pos += 1
        self.feature = pos

        pos = strchr( pos, '\t' )
        if pos == NULL: raise ValueError( "malformatted entry at %s" % buffer )
        pos[0] = '\0'
        pos += 1
        cstart = pos

        pos = strchr( pos, '\t' )
        if pos == NULL: raise ValueError( "malformatted entry at %s" % buffer )
        pos[0] = '\0'
        pos += 1
        cend = pos

        pos = strchr( pos, '\t' )
        if pos == NULL: raise ValueError( "malformatted entry at %s" % buffer )
        pos[0] = '\0'
        pos += 1
        self.score = pos

        pos = strchr( pos, '\t' )
        if pos == NULL: raise ValueError( "malformatted entry at %s" % buffer )
        pos[0] = '\0'
        pos += 1
        self.strand = pos

        pos = strchr( pos, '\t' )
        if pos == NULL: raise ValueError( "malformatted entry at %s" % buffer )
        pos[0] = '\0'
        pos += 1
        self.frame = pos

        pos = strchr( pos, '\t' )
        if pos == NULL: raise ValueError( "malformatted entry at %s" % buffer )
        pos[0] = '\0'
        pos += 1
        self.attributes = pos
        self.start = atoi( cstart ) - 1
        self.end = atoi( cend )
                      
    property contig:
       '''contig of feature.'''
       def __get__( self ): return self.contig
       def __set__( self, value ): 
           self.isModified = True
           self.contig = value

    property feature:
       '''feature name.'''
       def __get__( self ): return self.feature
       def __set__( self, value ): 
           self.isModified = True
           self.feature = value

    property source:
       '''feature source.'''
       def __get__( self ): return self.source
       def __set__( self, value ): 
           self.isModified = True
           self.source = value

    property start:
       '''feature start (in 0-based open/closed coordinates).'''
       def __get__( self ): return self.start
       def __set__( self, value ): 
           self.isModified = True
           self.start = value

    property end:
       '''feature end (in 0-based open/closed coordinates).'''
       def __get__( self ): return self.end
       def __set__( self, value ): 
           self.isModified = True
           self.end = value

    property score:
       '''feature score.'''
       def __get__( self ): 
           if self.score[0] == '.' and self.score[1] == '\0' :
               return None
           else:
               return atof(self.score)
       def __set__( self, value ): 
           self.isModified = True
           self.score = value

    property strand:
       '''feature strand.'''
       def __get__( self ): return self.strand
       def __set__( self, value ): 
           self.isModified = True
           self.strand = value

    property frame:
       '''feature frame.'''
       def __get__( self ): return self.frame
       def __set__( self, value ): 
           self.isModified = True
           self.frame = value

    property attributes:
       '''feature attributes (as a string).'''
       def __get__( self ): return self.attributes
       def __set__( self, value ): 
           self.isModified = True
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
        memcpy( self.attributes, p, l )

        self.hasOwnAttributes = True
        self.isModified = True

    def __str__(self):
        cdef char * cpy
        cdef int x

        if self.isModified:
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
            cpy = <char*>calloc( sizeof(char), self.nbytes+1 )
            assert cpy != NULL
            memcpy( cpy, self.data, self.nbytes+1)
            for x from 0 <= x < self.nbytes:
                if cpy[x] == '\0': cpy[x] = '\t'
            result = PyString_FromStringAndSize(cpy, self.nbytes)
            free(cpy)
            return result

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

    def __getitem__(self, item):
        return self.__getattr__( item )


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
        while start[0] == " ": start += 1
        if start[0] == '"':
            start += 1
            end = start
            while end[0] != '\0' and end[0] != '"': end += 1
            l = end - start + 1
            cpy = <char*>calloc( l, sizeof(char ) )
            memcpy( cpy, start, l )
            cpy[l-1] = '\0'
            result = cpy
            free(cpy)
            return result
        else:
            return start

    def setAttribute( self, name, value ):
        '''convenience method to set an attribute.'''
        r = self.asDict()
        r[name] = value
        self.fromDict( r )

cdef class VCFProxy( TupleProxy ):
    '''Proxy class for access to VCF fields.
    '''

    def __cinit__(self ): 
        # automatically calls TupleProxy.__cinit__
        pass

    def __dealloc__(self):
        # automatically calls TupleProxy.__dealloc__
        pass

    cdef update( self, char * buffer, size_t nbytes ):
        '''update internal data.

        nbytes does not include the terminal '\0'.
        '''
        cdef int end
        cdef char * cpos
        self.contig = buffer
        cdef char * pos

        if buffer[nbytes] != 0:
            raise ValueError( "incomplete line at %s" % buffer )
        
        cpos = pos = nextItem( buffer )
        self.id = pos = nextItem( pos )
        self.ref = pos = nextItem( pos )
        self.alt = pos = nextItem( pos )
        self.qual = pos = nextItem( pos )
        self.filter = pos = nextItem( pos )
        self.info = pos = nextItem( pos )
        self.format = pos = nextItem( pos )
        self.genotypes = pos = nextItem( pos )
        # vcf counts from 1 - correct here
        self.pos = atoi( cpos ) - 1

    def __str__(self):
        cdef char * cpy
        cdef int x

        if self.isModified:
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
            cpy = <char*>calloc( sizeof(char), self.nbytes+1 )
            assert cpy != NULL
            memcpy( cpy, self.data, self.nbytes+1)
            for x from 0 <= x < self.nbytes:
                if cpy[x] == '\0': cpy[x] = '\t'
            result = PyString_FromStringAndSize(cpy, self.nbytes)
            free(cpy)
            return result


