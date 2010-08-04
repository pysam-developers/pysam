# cython: embedsignature=True
# adds doc-strings for sphinx

import tempfile, os, sys, types, itertools, struct, ctypes

cdef class Tabixfile:
    '''*(filename, mode='r')*

    opens a :term:`tabix file` for reading. A missing
    index (*filename* + ".tbi") will raise an exception.
    '''

    cdef char * filename

    # pointer to tabixfile
    cdef tabix_t * tabixfile

    def __cinit__(self, *args, **kwargs ):
        self.tabixfile = NULL
        self._open( *args, **kwargs )

    def _isOpen( self ):
        '''return true if samfile has been opened.'''
        return self.tabixfile != NULL

    def _open( self, 
               char * filename, 
               mode ='r',
              ):
        '''open a :term:`tabix file` for reading.
        '''

        assert mode in ( "r",), "invalid file opening mode `%s`" % mode

        # close a previously opened file
        if self.tabixfile != NULL: self.close()
        self.tabixfile = NULL

        self.filename = filename
        filename_index = filename + ".tbi"

        if mode[0] == 'w':
            # open file for writing
            pass

        elif mode[0] == "r":
            # open file for reading
            if not os.path.exists( self.filename ):
                raise IOError( "file `%s` not found" % self.filename)

            if not os.path.exists( filename_index ):
                raise IOError( "index `%s` not found" % filename_index)

            # open file and load index
            self.tabixfile = ti_open( self.filename, filename_index )

        if self.tabixfile == NULL:
            raise IOError("could not open file `%s`" % filename )

    def _parseRegion( self, 
                      reference = None, 
                      start = None, 
                      end = None, 
                      region = None ):
        '''parse region information.

        raise ValueError for for invalid regions.

        returns a tuple of region, tid, start and end. Region
        is a valid samtools :term:`region` or None if the region
        extends over the whole file.

        Note that regions are 1-based, while start,end are python coordinates.
        '''
        ti_lazy_index_load( self.tabixfile )

        cdef int rtid
        cdef int rstart
        cdef int rend
        cdef int max_pos
        max_pos = 2 << 29

        rtid = rstart = rend = 0

        # translate to a region
        if reference:
            if start != None and end != None:
                region = "%s:%i-%i" % (reference, start+1, end)
            elif start == None and end != None:
                region = "%s:%i-%i" % (reference, 1, end)
            elif end == None and start != None:
                region = "%s:%i-%i" % (reference, start+1, max_pos-1)
            else:
                region = reference

        if region:
            ti_parse_region( self.tabixfile.idx, region, &rtid, &rstart, &rend)        
            if rtid < 0: raise ValueError( "invalid region `%s`" % region )
            if rstart > rend: raise ValueError( 'invalid region: start (%i) > end (%i)' % (rstart, rend) )
            if not 0 <= rstart < max_pos: raise ValueError( 'start out of range (%i)' % rstart )
            if not 0 <= rend < max_pos: raise ValueError( 'end out of range (%i)' % rend )

        return region, rtid, rstart, rend

    def fetch( self, 
               reference = None,
               start = None, 
               end = None, 
               region = None,
               parser = None ):
        '''
               
        fetch one or more rows in a :term:`region` using 0-based indexing. The region is specified by
        :term:`reference`, *start* and *end*. Alternatively, a samtools :term:`region` string can be supplied.

        Without *reference* or *region* all entries will be fetched. 
        
        If only *reference* is set, all reads matching on *reference* will be fetched.

        If *parser* is None, the results are returned as an unparsed string.
        Otherwise, *parser* is assumed to be a functor that will return parsed 
        data (see for example :meth:`asTuple` and :meth:`asGTF`).
        '''
        ti_lazy_index_load( self.tabixfile )

        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )

        region, rtid, rstart, rend = self._parseRegion( reference, start, end, region )

        if parser == None:
            if region:
                return TabixIterator( self, rtid, rstart, rend )
            else:
                return TabixIterator( self, -1, 0, 0 )
        else:
            if region:
                return TabixIteratorParsed( self, rtid, rstart, rend, parser )
            else:
                return TabixIteratorParsed( self, -1, 0, 0, parser )

    property contigs:
       '''chromosome names'''
       def __get__(self):
           cdef char ** sequences
           cdef int nsequences
           
           ti_lazy_index_load( self.tabixfile )
           sequences = ti_seqname( self.tabixfile.idx, &nsequences ) 
           cdef int x
           result = []
           for x from 0 <= x < nsequences:
               result.append( sequences[x] )
           return result
            
cdef class TabixIterator:
    """iterates over rows in *tabixfile* in region
    given by *tid*, *start* and *end*.
    """
    
    cdef ti_iter_t iterator
    cdef tabix_t * tabixfile

    def __cinit__(self, Tabixfile tabixfile, 
                  int tid, int start, int end ):
        
        assert tabixfile._isOpen()
        
        # makes sure that samfile stays alive as long as the
        # iterator is alive.
        self.tabixfile = tabixfile.tabixfile

        if tid < 0:
            # seek to start of file to ensure iteration is over
            # all entries.
            bgzf_seek( self.tabixfile.fp, 0, 0)
            self.iterator = ti_iter_first()
        else:
            self.iterator = ti_queryi(self.tabixfile, tid, start, end) 

        if <void*>self.iterator == NULL:
            raise ValueError("malformatted query or wrong sequence name.\n")

    def __iter__(self):
        return self 

    def __next__(self): 
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """
    
        cdef char * s
        cdef int len
        s = ti_read(self.tabixfile, self.iterator, &len)
        if s == NULL: raise StopIteration
        return s

    def __dealloc__(self):
        if <void*>self.iterator != NULL:
            ti_iter_destroy(self.iterator)

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

cdef class TupleProxy:
    '''Proxy class for access to parsed row as a tuple.

    This class represents a table row for fast read-access.
    '''

    cdef:
        char * data
        char ** fields
        int nfields
        int index

    def __cinit__(self ): 

        self.data = NULL
        self.fields = NULL
        self.index = 0

    cdef take( self, char * buffer, size_t nbytes ):
        '''start presenting buffer.

        Take ownership of the pointer.
        '''
        self.data = buffer
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

    def __dealloc__(self):
        if self.data != NULL:
            free(self.data)

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

cdef class GTFProxy:
    '''Proxy class for access to GTF fields.

    This class represents a GTF entry for fast read-access.
    Write-access has been added as well, though some care must
    be taken. If any of the string fields (contig, source, ...)
    are set, the new value is tied to the lifetime of the
    argument that was supplied.

    The only exception is the attributes field when set from
    a dictionary - this field will manage its own memory.

    '''

    cdef:
        char * contig
        char * source
        char * feature
        uint32_t start
        uint32_t end
        char * score
        char * strand
        char * frame
        char * attributes
        int nbytes
        char * data
        cdef bint isModified
        cdef bint hasOwnAttributes

    def __cinit__(self ): 
        self.data = NULL
        self.isModified = False
        self.hasOwnAttributes = False

    cdef take( self, char * buffer, size_t nbytes ):
        '''start presenting buffer.

        Take ownership of the pointer.
        '''
        self.data = buffer
        self.update( buffer, nbytes )
        self.isModified = False

    cdef present( self, char * buffer, size_t nbytes ):
        '''start presenting buffer.

        Do not take ownership of the pointer.
        '''
        self.update( buffer, nbytes )
        self.isModified = False

    cdef copy( self, char * buffer, size_t nbytes ):
        '''start presenting buffer.

        Take a copy of buffer.
        '''
        cdef int s
        # +1 for '\0'
        s = sizeof(char) *  (nbytes + 1)
        self.data = <char*>malloc( s ) 
        memcpy( <char*>self.data, buffer, s )
        self.update( self.data, nbytes )
        self.isModified = False

    cdef update( self, char * buffer, size_t nbytes ):
        '''update internal data.

        nbytes does not include the terminal '\0'.
        '''
        cdef int end
        cdef char * cstart, * cend, * cscore
        self.contig = buffer
        self.nbytes = nbytes
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
            memcpy( cpy, self.data, self.nbytes+1)
            for x from 0 <= x < self.nbytes:
                if cpy[x] == '\0': cpy[x] = '\t'
            result = cpy
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

    def __dealloc__(self):
        if self.data != NULL:
            free(self.data)
        if self.hasOwnAttributes:
            free(self.attributes)

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

cdef class Parser:
    pass

cdef class asTuple(Parser):
    '''converts a :term:`tabix row` into a python tuple.''' 
    def __call__(self, char * buffer, int len):
        cdef TupleProxy r
        r = TupleProxy()
        # need to copy - there were some
        # persistence issues with "present"
        r.copy( buffer, len )
        return r

cdef class asGTF(Parser):
    '''converts a :term:`tabix row` into a GTF record.''' 
    def __call__(self, char * buffer, int len):
        cdef GTFProxy r
        r = GTFProxy()
        r.copy( buffer, len )
        return r

cdef class TabixIteratorParsed:
    """iterates over mapped reads in a region.
    """
    
    cdef ti_iter_t iterator
    cdef tabix_t * tabixfile
    cdef Parser parser

    def __cinit__(self, 
                  Tabixfile tabixfile, 
                  int tid, 
                  int start, 
                  int end,
                  Parser parser ):

        assert tabixfile._isOpen()
        self.parser = parser

        # makes sure that samfile stays alive as long as the
        # iterator is alive.
        self.tabixfile = tabixfile.tabixfile

        if tid < 0:
            # seek to start of file to ensure iteration is over
            # all entries.
            bgzf_seek( self.tabixfile.fp, 0, 0)
            self.iterator = ti_iter_first()
        else:
            self.iterator = ti_queryi(self.tabixfile, tid, start, end) 

        if <void*>self.iterator == NULL:
            raise ValueError("malformatted query or wrong sequence name.\n")

    def __iter__(self):
        return self 

    def __next__(self): 
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """
    
        cdef char * s
        cdef int len
        s = ti_read(self.tabixfile, self.iterator, &len)
        if s == NULL: raise StopIteration
        return self.parser(s, len)

    def __dealloc__(self):
        if <void*>self.iterator != NULL:
            ti_iter_destroy(self.iterator)
        
def tabix_compress( filename_in, 
              filename_out,
              force = False ):

    '''
    compress *filename_in* writing the output to *filename_out*.
    
    Raise an IOError if *filename_out* already exists, unless *force* is set.
    '''

    if not force and os.path.exists(filename_out ):
        raise IOError( "Filename '%s' already exists, use *force* to overwrite" % filename_out)

    cdef int WINDOW_SIZE
    cdef int c, r
    cdef void * buffer
    cdef BGZF * fp
    cdef int fd_src

    cdef int O_RDONLY
    O_RDONLY = os.O_RDONLY

    WINDOW_SIZE = 64 * 1024

    fp = bgzf_open( filename_out, "w")
    if fp == NULL:
        raise IOError( "could not open '%s' for writing" )

    fd_src = open(filename_in, O_RDONLY)
    if fd_src == 0:
        raise IOError( "could not open '%s' for reading" )

    buffer = malloc(WINDOW_SIZE)

    while c > 0:
        c = read(fd_src, buffer, WINDOW_SIZE)
        r = bgzf_write(fp, buffer, c)
        if r < 0:
            free( buffer )
            raise OSError("writing failed")
        
    free( buffer )
    r = bgzf_close(fp)
    if r < 0: raise OSError("writing failed")

def tabix_index( filename, 
                 force = False,
                 seq_col = None, 
                 start_col = None, 
                 end_col = None,
                 preset = None,
                 meta_char = "#",
                 zerobased = False,
                ):
    '''
    index tab-separated *filename* using tabix.

    An existing index will not be overwritten unless
    *force* is set.

    The index will be built from coordinates
    in columns *seq_col*, *start_col* and *end_col*.

    The contents of *filename* have to be sorted by 
    contig and position - the method does not check
    if the file is sorted.

    Column indices are 0-based. Coordinates in the file
    are assumed to be 1-based.

    If *preset* is provided, the column coordinates
    are taken from a preset. Valid values for preset
    are "gff", "bed", "sam", "vcf", psltbl", "pileup".
    
    Lines beginning with *meta_char* and the first
    *line_skip* lines will be skipped.
    
    If *filename* does not end in ".gz", it will be automatically
    compressed. The original file will be removed and only the 
    compressed file will be retained. 

    If *filename* ends in *gz*, the file is assumed to be already
    compressed with bgzf.

    returns the filename of the compressed data
    '''
    
    if not os.path.exists(filename): raise IOError("No such file '%s'" % filename)

    if not filename.endswith(".gz"): 
        
        tabix_compress( filename, filename + ".gz", force = force )
        os.unlink( filename )
        filename += ".gz"

    if not force and os.path.exists(filename + ".tbi" ):
        raise IOError( "Filename '%s.tbi' already exists, use *force* to overwrite" )

    # columns (1-based)
    # preset-code, contig, start, end, metachar for commends, lines to ignore at beginning
    # 0 is a missing column
    preset2conf = {
        'gff' : ( 0, 1, 4, 5, ord('#'), 0 ),
        'bed' : ( 0x10000, 1, 2, 3, ord('#'), 0 ),
        'psltbl' : ( 0x10000, 15, 17, 18, ord('#'), 0 ),
        'sam' : ( 1, 3, 4, 0, ord('#'), 0 ),
        'vcf' : ( 2, 1, 2, 0, ord('#'), 0 ),
        'pileup': (3, 1, 2, 0, ord('#'), 0 ),
        }

    if preset:
        try:
            conf_data = preset2conf[preset]
        except KeyError:
            raise KeyError( "unknown preset '%s', valid presets are '%s'" % (preset, ",".join(preset2conf.keys() )))
    else:
        if end_col == None: end_col = -1
        preset = 0

        # note that tabix internally works with 0-based coordinates and open/closed intervals.
        # When using a preset, conversion is automatically taken care of.
        # Otherwise, the coordinates are assumed to be 1-based closed intervals and 
        # -1 is subtracted from the start coordinate. To avoid doing this, set
        # the TI_FLAG_UCSC=0x10000 flag:
        if zerobased: preset = preset | 0x10000

        conf_data = (preset, seq_col+1, start_col+1, end_col+1, ord(meta_char), 0)
                
    cdef ti_conf_t conf
    conf.preset, conf.sc, conf.bc, conf.ec, conf.meta_char, conf.line_skip = conf_data

    ti_index_build( filename, &conf)
    
    return filename
    
__all__ = ["tabix_index", 
           "tabix_compress",
           "Tabixfile", 
           "asTuple",
           "asGTF",
           ]
