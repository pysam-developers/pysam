# cython: embedsignature=True
# adds doc-strings for sphinx

# Helper functions for python 3 compatibility - taken
# from csamtools.pyx
import tempfile, os, sys, types, itertools, struct, ctypes, gzip
import io
cimport TabProxies

from cpython cimport PyErr_SetString, PyBytes_Check, \
    PyUnicode_Check, PyBytes_FromStringAndSize, \
    PyObject_AsFileDescriptor

PYTHON3 = PY_MAJOR_VERSION >= 3

# from cpython cimport PyString_FromStringAndSize, PyString_AS_STRING
from cpython.version cimport PY_MAJOR_VERSION

cdef from_string_and_size(char* s, size_t length):
    if PY_MAJOR_VERSION < 3:
        return s[:length]
    else:
        return s[:length].decode("ascii")

# filename encoding (copied from lxml.etree.pyx)
cdef str _FILENAME_ENCODING
_FILENAME_ENCODING = sys.getfilesystemencoding()
if _FILENAME_ENCODING is None:
    _FILENAME_ENCODING = sys.getdefaultencoding()
if _FILENAME_ENCODING is None:
    _FILENAME_ENCODING = 'ascii'

#cdef char* _C_FILENAME_ENCODING
#_C_FILENAME_ENCODING = <char*>_FILENAME_ENCODING

cdef bytes _my_encodeFilename(object filename):
    u"""Make sure a filename is 8-bit encoded (or None).
    """
    if filename is None:
        return None
    elif PyBytes_Check(filename):
        return filename
    elif PyUnicode_Check(filename):
        return filename.encode(_FILENAME_ENCODING)
    else:
        raise TypeError, u"Argument must be string or unicode."

cdef bytes _force_bytes(object s):
    u"""convert string or unicode object to bytes, assuming ascii encoding.
    """
    if PY_MAJOR_VERSION < 3:
        return s
    elif s is None:
        return None
    elif PyBytes_Check(s):
        return s
    elif PyUnicode_Check(s):
        return s.encode('ascii')
    else:
        raise TypeError, u"Argument must be string, bytes or unicode."

cdef inline bytes _force_cmdline_bytes(object s):
    return _force_bytes(s)

cdef _charptr_to_str(char* s):
    if PY_MAJOR_VERSION < 3:
        return s
    else:
        return s.decode("ascii")

cdef _force_str(object s):
    """Return s converted to str type of current Python (bytes in Py2, unicode in Py3)"""
    if s is None:
        return None
    if PY_MAJOR_VERSION < 3:
        return s
    elif PyBytes_Check(s):
        return s.decode('ascii')
    else:
        # assume unicode
        return s


cdef class Tabixfile:
    '''*(filename, mode='r')*

    opens a :term:`tabix file` for reading. A missing
    index (*filename* + ".tbi") will raise an exception.
    '''
    def __cinit__(self, filename, mode = 'r', *args, **kwargs ):
        self.tabixfile = NULL
        self._open( filename, mode, *args, **kwargs )

    def _isOpen( self ):
        '''return true if samfile has been opened.'''
        return self.tabixfile != NULL

    def _open( self, 
               filename,
               mode ='r',
              ):
        '''open a :term:`tabix file` for reading.
        '''

        assert mode in ( "r",), "invalid file opening mode `%s`" % mode

        # close a previously opened file
        if self.tabixfile != NULL: self.close()
        self.tabixfile = NULL

        filename_index = filename + ".tbi"
        self.isremote = filename.startswith( "http:") or filename.startswith( "ftp:" )

        # encode all the strings
        filename = _my_encodeFilename(filename)
        filename_index = _my_encodeFilename(filename_index)
        cdef bytes bmode = mode.encode('ascii')

        if self._filename != NULL: free(self._filename )

        self._filename = strdup(filename)

        if mode[0] == 'w':
            # open file for writing
            raise NotImplementedError("writing to tabix files not implemented" )

        elif mode[0] == "r":
            # open file for reading
            
            if not self.isremote:
                if not os.path.exists( filename ):
                    raise IOError( "file `%s` not found" % filename)

                if not os.path.exists( filename_index ):
                    raise IOError( "index `%s` not found" % filename_index)

            # open file and load index
            self.tabixfile = ti_open( filename, filename_index )

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
            region = _force_bytes(region)
            ti_parse_region( self.tabixfile.idx, region, 
                             &rtid, &rstart, &rend)        
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

    ###############################################################
    ###############################################################
    ###############################################################
    ## properties
    ###############################################################
    property filename:
        '''filename associated with this object.'''
        def __get__(self):
            if not self._isOpen(): raise ValueError( "I/O operation on closed file" )
            return self._filename

    property header:
        '''the file header.
          
        .. note::
            The header is returned as an iterator over lines without the
            newline character.
        '''
        
        def __get__( self ):
            return TabixHeaderIterator( self )

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
            
    def close( self ):
        '''
        closes the :class:`pysam.Tabixfile`.'''
        if self.tabixfile != NULL:
            ti_close( self.tabixfile )
            self.tabixfile = NULL

    def __dealloc__( self ):
        # remember: dealloc cannot call other python methods
        # note: no doc string
        # note: __del__ is not called.
        if self.tabixfile != NULL:
            ti_close( self.tabixfile )
            self.tabixfile = NULL
        if self._filename != NULL: free( self._filename )

cdef class TabixIterator:
    """iterates over rows in *tabixfile* in region
    given by *tid*, *start* and *end*.
    """
    
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
        # metachar filtering does not work within tabix 
        # though it should. Getting the metachar is a pain
        # as ti_index_t is incomplete type.

        # simply use '#' for now.
        while 1:
            s = ti_read(self.tabixfile, self.iterator, &len)
            if s == NULL: raise StopIteration
            if s[0] != '#': break

        retval = _charptr_to_str( s )
        return retval

    def __dealloc__(self):
        if <void*>self.iterator != NULL:
            ti_iter_destroy(self.iterator)

cdef class TabixHeaderIterator:
    """return header lines.
    """
    
    def __cinit__(self, Tabixfile tabixfile ):
        
        assert tabixfile._isOpen()
        
        # makes sure that samfile stays alive as long as the
        # iterator is alive.
        self.tabixfile = tabixfile.tabixfile

        self.iterator = ti_query(self.tabixfile, NULL, 0, 0) 

        if <void*>self.iterator == NULL:
            raise ValueError("can't open header.\n")

    def __iter__(self):
        return self 

    def __next__(self): 
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """
    
        cdef char * s
        cdef int len

        # Getting the metachar is a pain as ti_index_t is incomplete type.
        # simply use '#' for now.
        s = ti_read(self.tabixfile, self.iterator, &len)
        if s == NULL: raise StopIteration
        # stop at first non-header line
        if s[0] != '#': raise StopIteration

        return s

    def __dealloc__(self):
        if <void*>self.iterator != NULL:
            ti_iter_destroy(self.iterator)


#########################################################
#########################################################
#########################################################
cdef class Parser:
    pass

cdef class asTuple(Parser):
    '''converts a :term:`tabix row` into a python tuple.

    Access is by numeric index.
    ''' 
    def __call__(self, char * buffer, int len):
        cdef TabProxies.TupleProxy r
        r = TabProxies.TupleProxy()
        # need to copy - there were some
        # persistence issues with "present"
        r.copy( buffer, len )
        return r

cdef class asGTF(Parser):
    '''converts a :term:`tabix row` into a GTF record with the following 
    fields:

    contig
       contig
    feature
       feature
    source
       source
    start
       genomic start coordinate (0-based)
    end
       genomic end coordinate plus one (0-based)
    score
       feature score
    strand
       strand
    frame
       frame
    attributes
       attribute string.

    GTF formatted entries also defined the attributes:

    gene_id
       the gene identifier
    transcript_ind
       the transcript identifier
    
    ''' 
    def __call__(self, char * buffer, int len):
        cdef TabProxies.GTFProxy r
        r = TabProxies.GTFProxy()
        r.copy( buffer, len )
        return r

cdef class asBed( Parser ):
    '''converts a :term:`tabix row` into a bed record
    with the following fields:

    contig
       contig
    start
       genomic start coordinate (zero-based)
    end
       genomic end coordinate plus one (zero-based)
    name
       name of feature.
    score
       score of feature
    strand
       strand of feature
    thickStart
       thickStart
    thickEnd
       thickEnd
    itemRGB
       itemRGB
    blockCount
       number of bocks
    blockSizes
       ',' separated string of block sizes
    blockStarts
       ',' separated string of block genomic start positions

    Only the first three fields are required. Additional
    fields are optional, but if one is defined, all the preceeding
    need to be defined as well.

    ''' 
    def __call__(self, char * buffer, int len):
        cdef TabProxies.BedProxy r
        r = TabProxies.BedProxy()
        r.copy( buffer, len )
        return r

cdef class asVCF( Parser ): 
    '''converts a :term:`tabix row` into a VCF record with
    the following fields:
    
    contig
       contig
    pos
       chromosomal position, zero-based
    id 
       id
    ref
       reference
    alt
       alt
    qual
       qual
    filter
       filter
    info
       info
    format
       format specifier.

    Access to genotypes is via index::

        contig = vcf.contig
        first_sample_genotype = vcf[0]
        second_sample_genotype = vcf[1]

    '''
    def __call__(self, char * buffer, int len ):
        cdef TabProxies.VCFProxy r
        r = TabProxies.VCFProxy()
        r.copy( buffer, len )
        return r
    
#########################################################
#########################################################
#########################################################
cdef class TabixIteratorParsed:
    """iterates over mapped reads in a region.

    Returns parsed data.
    """

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
        while 1:
            s = ti_read(self.tabixfile, self.iterator, &len)
            if s == NULL: raise StopIteration
            # todo: read metachar from configuration
            if s[0] != '#': break
            
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

    fn = _force_bytes(filename_out)
    fp = bgzf_open( fn, "w")
    if fp == NULL:
        raise IOError( "could not open '%s' for writing" )

    fn = _force_bytes(filename_in)
    fd_src = open(fn, O_RDONLY)
    if fd_src == 0:
        raise IOError( "could not open '%s' for reading" )

    buffer = malloc(WINDOW_SIZE)
    c = 1

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

    if preset == None and (seq_col == None or start_col == None or end_col == None):
        raise ValueError("neither preset nor seq_col,start_col and end_col given" )

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

    fn = _my_encodeFilename( filename )
    ti_index_build( fn, &conf)
    
    return filename

#########################################################
#########################################################
#########################################################
## Iterators for parsing through unindexed files.
#########################################################
ctypedef class tabix_inplace_iterator:
    '''iterate over ``infile``.

    This iterator is not safe. If the :meth:`__next__()` method is called 
    after ``infile`` is closed, the result is undefined (see ``fclose()``).

    The iterator might either raise a StopIteration or segfault.
    '''


    def __cinit__(self, infile, int buffer_size = 65536 ):

        cdef int fd = PyObject_AsFileDescriptor( infile )
        if fd == -1: raise ValueError( "I/O operation on closed file." )
        self.infile = fdopen( fd, 'r')

        if self.infile == NULL: raise ValueError( "I/O operation on closed file." )

        self.buffer = <char*>malloc( buffer_size )        
        self.size = buffer_size

    def __iter__(self):
        return self

    cdef __cnext__(self):

        cdef char * b
        cdef size_t nbytes
        b = self.buffer
        r = self.Parser()

        while not feof( self.infile ):
            nbytes = getline( &b, &self.size, self.infile)

            # stop at first error or eof
            if (nbytes == -1): break
            # skip comments
            if (b[0] == '#'): continue

            # skip empty lines
            if b[0] == '\0' or b[0] == '\n' or b[0] == '\r': continue

            # make sure that entry is complete
            if b[nbytes-1] != '\n' and b[nbytes-1] != '\r':
                result = b
                raise ValueError( "incomplete line at %s" % result )

            # make sure that this goes fully through C
            # otherwise buffer is copied to/from a
            # Python object causing segfaults as
            # the wrong memory is freed
            r.present( b, nbytes )
            return r 

        raise StopIteration

    def __dealloc__(self):
        free(self.buffer)

    def __next__(self):
        return self.__cnext__()

ctypedef class tabix_copy_iterator:
    '''iterate over ``infile``.

    This iterator is not save. If the :meth:`__next__()` method is called 
    after ``infile`` is closed, the result is undefined (see ``fclose()``).

    The iterator might either raise a StopIteration or segfault.
    '''

    def __cinit__(self, infile, Parser parser ):

        cdef int fd = PyObject_AsFileDescriptor( infile )
        if fd == -1: raise ValueError( "I/O operation on closed file." )
        self.infile = fdopen( fd, 'r')
        if self.infile == NULL: raise ValueError( "I/O operation on closed file." )
        self.parser = parser

    def __iter__(self):
        return self

    cdef __cnext__(self):

        cdef char * b
        cdef size_t nbytes
        cdef int x

        b = NULL        

        while not feof( self.infile ):

            # getline allocates on demand
            # return number of characters read excluding null byte
            nbytes = getline( &b, &nbytes, self.infile)
            # stop at first error
            if (nbytes == -1): break
            # skip comments
            if (b[0] == '#'): continue

            # skip empty lines
            if b[0] == '\0' or b[0] == '\n' or b[0] == '\r': continue

            # make sure that entry is complete
            if b[nbytes-1] != '\n' and b[nbytes-1] != '\r':
                result = b
                free(b)
                raise ValueError( "incomplete line at %s" % result )

            # make sure that this goes fully through C
            # otherwise buffer is copied to/from a
            # Python object causing segfaults as
            # the wrong memory is freed
            # -1 to remove the new-line character
            return self.parser(b, nbytes)

        free(b)
        raise StopIteration

    def __next__(self):
        return self.__cnext__()
    
class tabix_generic_iterator:
    '''iterate over ``infile``.
    
    Permits the use of file-like objects for example from the gzip module.
    '''
    def __init__(self, infile, parser ):

        self.infile = infile
        if self.infile.closed: raise ValueError( "I/O operation on closed file." )
        self.parser = parser

    def __iter__(self):
        return self

    # cython version - required for python 3
    def __next__(self):
        
        cdef char * b, * cpy
        cdef size_t nbytes
        while 1:

            line = self.infile.readline()
            if not line: break
            
            s = _force_bytes( line )
            b = s
            nbytes = len( line )
            assert b[nbytes] == '\0'

            # skip comments
            if (b[0] == '#'): continue

            # skip empty lines
            if b[0] == '\0' or b[0] == '\n' or b[0] == '\r': continue
            
            # make sure that entry is complete
            if b[nbytes-1] != '\n' and b[nbytes-1] != '\r':
                raise ValueError( "incomplete line at %s" % line )
            
            # create a copy
            cpy = <char*>malloc(nbytes+1)        
            if cpy == NULL: raise MemoryError()
            memcpy( cpy, b, nbytes+1)

            return self.parser(cpy, nbytes)            

        raise StopIteration

    # python version - required for python 2.7
    def next(self):
        return self.__next__()
    
def tabix_iterator( infile, parser ):
    """return an iterator over all entries in a file."""
    return tabix_generic_iterator( infile, parser )
    # file objects can use C stdio
    # used to be: isinstance( infile, file):
    # if PYTHON3:
    #     if isinstance( infile, io.IOBase ):
    #         return tabix_copy_iterator( infile, parser )
    #     else:
    #         return tabix_generic_iterator( infile, parser )
    # else:
#        if isinstance( infile, file ):
#            return tabix_copy_iterator( infile, parser )
#        else:
#            return tabix_generic_iterator( infile, parser )
    
__all__ = ["tabix_index", 
           "tabix_compress",
           "Tabixfile", 
           "asTuple",
           "asGTF",
           "asVCF",
           "asBed",
           "tabix_iterator", 
           "tabix_inplace_iterator"
           ]
