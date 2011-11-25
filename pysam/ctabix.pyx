# cython: embedsignature=True
# adds doc-strings for sphinx

import tempfile, os, sys, types, itertools, struct, ctypes, gzip
from cpython cimport PyString_FromStringAndSize, PyString_AS_STRING
cimport TabProxies

cdef class Tabixfile:
    '''*(filename, mode='r')*

    opens a :term:`tabix file` for reading. A missing
    index (*filename* + ".tbi") will raise an exception.
    '''

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

        if self._filename != NULL: free(self._filename )
        self._filename = strdup( filename )

        filename_index = filename + ".tbi"

        if mode[0] == 'w':
            # open file for writing
            pass

        elif mode[0] == "r":
            # open file for reading
            if not os.path.exists( self._filename ):
                raise IOError( "file `%s` not found" % self._filename)

            if not os.path.exists( filename_index ):
                raise IOError( "index `%s` not found" % filename_index)

            # open file and load index
            self.tabixfile = ti_open( self._filename, filename_index )

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
        # metachar filtering does not work within tabix 
        # though it should. Getting the metachar is a pain
        # as ti_index_t is incomplete type.

        # simply use '#' for now.
        while 1:
            s = ti_read(self.tabixfile, self.iterator, &len)
            if s == NULL: raise StopIteration
            if s[0] != '#': break

        return s

    def __dealloc__(self):
        if <void*>self.iterator != NULL:
            ti_iter_destroy(self.iterator)

cdef class TabixHeaderIterator:
    """return header lines.
    """
    
    cdef ti_iter_t iterator
    cdef tabix_t * tabixfile

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
           "asVCF",
           "asBed",
           ]
