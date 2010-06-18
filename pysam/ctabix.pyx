# cython: embedsignature=True
# adds doc-strings for sphinx

import tempfile, os, sys, types, itertools, struct, ctypes
 
cdef class Tabixfile:
    '''*(filename, mode='r', template = None, referencenames = None, referencelengths = None, text = NULL, header = None)*
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
        '''open a sam/bam file.

        If _open is called on an existing bamfile, the current file will be
        closed and a new file will be opened.
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
               region = None ):
        '''*(reference = None, start = None, end = None, region = None, callback = None, until_eof = False)*
               
        fetch :meth:`AlignedRead` objects in a :term:`region` using 0-based indexing. The region is specified by
        :term:`reference`, *start* and *end*. Alternatively, a samtools :term:`region` string can be supplied.

        Without *reference* or *region* all entries will be fetched. 
        
        If only *reference* is set, all reads matching on *reference* will be fetched.
        '''
        ti_lazy_index_load( self.tabixfile )

        if not self._isOpen():
            raise ValueError( "I/O operation on closed file" )

        region, rtid, rstart, rend = self._parseRegion( reference, start, end, region )

        if region:
            return TabixIterator( self, rtid, rstart, rend )
        else:
            return TabixIterator( self, -1, 0, 0 )

cdef class TabixIterator:
    """iterates over mapped reads in a region.
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
        

def tabix_compress( filename_in, 
              filename_out,
              force = False ):

    '''compress a file with bgzf.'''

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
              
                ):
    '''index tab-separated *filename* using tabix.

    An existing index will not be overwritten unless
    *force* is set.

    The index will be built from coordinates
    in columns *seq_col*, *start_col* and *end_col*.

    The contents of *filename* have to be sorted by 
    contig and position - the method does not check
    if the file is sorted.

    Column indices are 0-based.

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
        conf_data = (0, seq_col+1, start_col+1, end_col+1, ord(meta_char), 0)
                
    cdef ti_conf_t conf
    conf.preset, conf.sc, conf.bc, conf.ec, conf.meta_char, conf.line_skip = conf_data

    ti_index_build( filename, &conf)
    
    return filename
    
