# cython: embedsignature=True
# adds doc-strings for sphinx
import os
import sys

from libc.stdio cimport printf, fprintf, stderr
from libc.string cimport strerror
from libc.errno cimport errno
from posix.unistd cimport dup

from cpython cimport PyErr_SetString, PyBytes_Check, \
    PyUnicode_Check, PyBytes_FromStringAndSize, \
    PyObject_AsFileDescriptor

from cpython.version cimport PY_MAJOR_VERSION

cimport TabProxies

from chtslib cimport htsFile, hts_open, hts_close, HTS_IDX_START,\
    BGZF, bgzf_open, bgzf_close, bgzf_write, \
    ks_init, ks_destroy, gzFile, ks_getuntil, kstring_t, \
    tbx_index_build, tbx_index_load, tbx_itr_queryi, tbx_itr_querys, \
    tbx_conf_t, tbx_seqnames, tbx_itr_next, tbx_itr_destroy, \
    tbx_destroy, gzopen, gzclose, gzerror, gzdopen

PYTHON3 = PY_MAJOR_VERSION >= 3

# filename encoding (copied from lxml.etree.pyx)
cdef str _FILENAME_ENCODING
_FILENAME_ENCODING = sys.getfilesystemencoding()
if _FILENAME_ENCODING is None:
    _FILENAME_ENCODING = sys.getdefaultencoding()
if _FILENAME_ENCODING is None:
    _FILENAME_ENCODING = 'ascii'

#cdef char* _C_FILENAME_ENCODING
#_C_FILENAME_ENCODING = <char*>_FILENAME_ENCODING

cdef inline bytes _encodeFilename(object filename):
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

cdef inline bytes _force_bytes(object s, encoding="ascii"):
    u"""convert string or unicode object to bytes, assuming ascii encoding.
    """
    if PY_MAJOR_VERSION < 3:
        return s
    elif s is None:
        return None
    elif PyBytes_Check(s):
        return s
    elif PyUnicode_Check(s):
        return s.encode(encoding)
    else:
        raise TypeError, u"Argument must be string, bytes or unicode."

cdef inline _charptr_to_str(char* s, encoding="ascii"):
    if PY_MAJOR_VERSION < 3:
        return s
    else:
        return s.decode(encoding)

cdef _force_str(object s, encoding="ascii"):
    """Return s converted to str type of current Python
    (bytes in Py2, unicode in Py3)"""
    if s is None:
        return None
    if PY_MAJOR_VERSION < 3:
        return s
    elif PyBytes_Check(s):
        return s.decode(encoding)
    else:
        # assume unicode
        return s


cdef class Parser:

    def __init__(self, encoding="ascii"):
        self.encoding = encoding

    def set_encoding(self, encoding):
        self.encoding = encoding

    def get_encoding(self):
        return self.encoding

    cdef parse(self, char * buffer, int length):
        raise NotImplementedError(
            'parse method of %s not implemented' % str(self))

    def __call__(self, char * buffer, int length):
        return self.parse(buffer, length)


cdef class asTuple(Parser):
    '''converts a :term:`tabix row` into a python tuple.

    A field in a row is accessed by numeric index.
    ''' 
    cdef parse(self, char * buffer, int len):
        cdef TabProxies.TupleProxy r
        r = TabProxies.TupleProxy(self.encoding)
        # need to copy - there were some
        # persistence issues with "present"
        r.copy(buffer, len)
        return r


cdef class asGTF(Parser):
    '''converts a :term:`tabix row` into a GTF record with the following
    fields:
   
    +----------+----------+-------------------------------+
    |*Column*  |*Name*    |*Content*                      |
    +----------+----------+-------------------------------+
    |1         |contig    |the chromosome name            |
    +----------+----------+-------------------------------+
    |2         |feature   |The feature type               |
    +----------+----------+-------------------------------+
    |3         |source    |The feature source             |
    +----------+----------+-------------------------------+
    |4         |start     |genomic start coordinate       |
    |          |          |(0-based)                      |
    +----------+----------+-------------------------------+
    |5         |end       |genomic end coordinate         |
    |          |          |(0-based)                      |
    +----------+----------+-------------------------------+
    |6         |score     |feature score                  |
    +----------+----------+-------------------------------+
    |7         |strand    |strand                         |
    +----------+----------+-------------------------------+
    |8         |frame     |frame                          |
    +----------+----------+-------------------------------+
    |9         |attributes|the attribute field            |
    +----------+----------+-------------------------------+

    GTF formatted entries also define the following fields that
    are derived from the attributes field:

    +--------------------+------------------------------+
    |*Name*              |*Content*                     |
    +--------------------+------------------------------+
    |gene_id             |the gene identifier           |
    +--------------------+------------------------------+
    |transcript_id       |the transcript identifier     |
    +--------------------+------------------------------+

    ''' 
    cdef parse(self, char * buffer, int len):
        cdef TabProxies.GTFProxy r
        r = TabProxies.GTFProxy(self.encoding)
        r.copy(buffer, len)
        return r


cdef class asBed(Parser):
    '''converts a :term:`tabix row` into a bed record
    with the following fields:

    +-----------+-----------+------------------------------------------+
    |*Column*   |*Field*    |*Contents*                                |
    |           |           |                                          |
    +-----------+-----------+------------------------------------------+
    |1          |contig     |contig                                    |
    |           |           |                                          |
    +-----------+-----------+------------------------------------------+
    |2          |start      |genomic start coordinate (zero-based)     |
    +-----------+-----------+------------------------------------------+
    |3          |end        |genomic end coordinate plus one           |
    |           |           |(zero-based)                              |
    +-----------+-----------+------------------------------------------+
    |4          |name       |name of feature.                          |
    +-----------+-----------+------------------------------------------+
    |5          |score      |score of feature                          |
    +-----------+-----------+------------------------------------------+
    |6          |strand     |strand of feature                         |
    +-----------+-----------+------------------------------------------+
    |7          |thickStart |thickStart                                |
    +-----------+-----------+------------------------------------------+
    |8          |thickEnd   |thickEnd                                  |
    +-----------+-----------+------------------------------------------+
    |9          |itemRGB    |itemRGB                                   |
    +-----------+-----------+------------------------------------------+
    |10         |blockCount |number of bocks                           |
    +-----------+-----------+------------------------------------------+
    |11         |blockSizes |',' separated string of block sizes       |
    +-----------+-----------+------------------------------------------+
    |12         |blockStarts|',' separated string of block genomic     |
    |           |           |start positions                           |
    +-----------+-----------+------------------------------------------+

    Only the first three fields are required. Additional
    fields are optional, but if one is defined, all the preceeding
    need to be defined as well.

    ''' 
    cdef parse(self, char * buffer, int len):
        cdef TabProxies.BedProxy r
        r = TabProxies.BedProxy(self.encoding)
        r.copy(buffer, len)
        return r


cdef class asVCF(Parser): 
    '''converts a :term:`tabix row` into a VCF record with
    the following fields:
    
    +----------+---------+------------------------------------+
    |*Column*  |*Field*  |*Contents*                          |
    |          |         |                                    |
    +----------+---------+------------------------------------+
    |1         |contig   |chromosome                          |
    +----------+---------+------------------------------------+
    |2         |pos      |chromosomal position, zero-based    |
    +----------+---------+------------------------------------+
    |3         |id       |id                                  |
    +----------+---------+------------------------------------+
    |4         |ref      |reference allele                    |
    +----------+---------+------------------------------------+
    |5         |alt      |alternate alleles                   |
    +----------+---------+------------------------------------+
    |6         |qual     |quality                             |
    +----------+---------+------------------------------------+
    |7         |filter   |filter                              |
    +----------+---------+------------------------------------+
    |8         |info     |info                                |
    +----------+---------+------------------------------------+
    |9         |format   |format specifier.                   |
    +----------+---------+------------------------------------+

    Access to genotypes is via index::

        contig = vcf.contig
        first_sample_genotype = vcf[0]
        second_sample_genotype = vcf[1]

    '''
    cdef parse(self, char * buffer, int len):
        cdef TabProxies.VCFProxy r
        r = TabProxies.VCFProxy(self.encoding)
        r.copy(buffer, len)
        return r


cdef class TabixFile:
    '''*(filename, mode='r', parser = None)*

    opens a :term:`tabix file` for reading. A missing
    index (*filename* + ".tbi") will raise an exception. *index*
    specifies an alternative name of the index.

    *parser* sets the default parser for this tabix file. If *parser*
    is None, the results are returned as an unparsed string.
    Otherwise, *parser* is assumed to be a functor that will return
    parsed data (see for example :class:`~pysam.asTuple` and
    :class:`~pysam.asGTF`).

    '''
    def __cinit__(self,
                  filename,
                  mode = 'r',
                  parser=None,
                  index=None,
                  encoding="ascii",
                  *args,
                  **kwargs ):

        self.tabixfile = NULL
        self.parser = parser
        self._open(filename, mode, index, *args, **kwargs)
        self.encoding = encoding

    def _open( self, 
               filename,
               mode='r',
               index=None,
              ):
        '''open a :term:`tabix file` for reading.
        '''

        assert mode in ("r",), "invalid file opening mode `%s`" % mode

        if self.tabixfile != NULL:
            self.close()
        self.tabixfile = NULL

        filename_index = index or (filename + ".tbi")
        self.isremote = filename.startswith("http:") or filename.startswith("ftp:")

        if not self.isremote:
            if not os.path.exists(filename):
                raise IOError("file `%s` not found" % filename)

            if not os.path.exists(filename_index):
                raise IOError("index `%s` not found" % filename_index)

        self._filename = filename
        self._filename_index = filename_index

        # encode all the strings to pass to tabix
        _encoded_filename = _encodeFilename(filename)
        _encoded_index = _encodeFilename(filename_index)

        # open file
        self.tabixfile = hts_open(_encoded_filename, 'r')
        if self.tabixfile == NULL:
            raise IOError("could not open file `%s`" % filename)
        
        self.index = tbx_index_load(_encoded_index)
        if self.index == NULL:
            raise IOError("could not open index for `%s`" % filename)

    def _dup(self):
        '''return a copy of this tabix file.
        
        The file is being re-opened.
        '''
        return TabixFile(self._filename,
                         mode="r", 
                         parser=self.parser,
                         index=self._filename_index,
                         encoding=self.encoding)

    def _isOpen(self):
        '''return true if samfile has been opened.'''
        return self.tabixfile != NULL


    def fetch(self, 
              reference=None,
              start=None, 
              end=None, 
              region=None,
              parser=None,
              multiple_iterators=False):
        '''fetch one or more rows in a :term:`region` using 0-based
        indexing. The region is specified by :term:`reference`,
        *start* and *end*. Alternatively, a samtools :term:`region`
        string can be supplied.

        Without *reference* or *region* all entries will be fetched. 
        
        If only *reference* is set, all reads matching on *reference*
        will be fetched.

        If *parser* is None, the default parser will be used for
        parsing.
        
        Set *multiple_iterators* to true if you will be using multiple
        iterators on the same file at the same time. The iterator
        returned will receive its own copy of a filehandle to the file
        effectively re-opening the file. Re-opening a file creates
        some overhead, so beware.

        '''
        if not self._isOpen():
            raise ValueError("I/O operation on closed file")

        # convert coordinates to region string
        if reference:
            if end is not None:
                if start is None:
                    start = 0
                region = '%s:%i-%i' % (reference, start + 1, end)
                if start > end:
                    raise ValueError(
                        'start (%i) > end (%i)' % (start, end))
            elif start is not None:
                region = '%s:%i' % (reference, start + 1)
            else:
                region = reference

        # get iterator
        cdef hts_itr_t * iter
        cdef TabixFile fileobj

        # reopen the same file if necessary
        if multiple_iterators:
            fileobj = self._dup()
        else:
            fileobj = self

        if region is None:
            # without region or reference - iterate from start
            iter = tbx_itr_queryi(fileobj.index,
                                  HTS_IDX_START,
                                  0,
                                  0)
        else:
            s = _force_bytes(region, encoding=fileobj.encoding)
            iter = tbx_itr_querys(fileobj.index, s)

        if iter == NULL:
            if region is None:
                # possible reason is that the file is empty -
                # return an empty iterator
                return EmptyIterator()
            else:
                raise ValueError(
                    "could not create iterator for region '%s'" %
                    region)
            
        # use default parser if no parser is specified
        if parser is None:
            parser = fileobj.parser

        cdef TabixIterator a
        if parser is None: 
            a = TabixIterator(encoding=fileobj.encoding)
        else:
            parser.set_encoding(fileobj.encoding)
            a = TabixIteratorParsed(parser)

        a.tabixfile = fileobj
        a.iterator = iter

        return a

    ###############################################################
    ###############################################################
    ###############################################################
    ## properties
    ###############################################################
    property filename:
        '''filename associated with this object.'''
        def __get__(self):
            if not self._isOpen():
                raise ValueError("I/O operation on closed file")
            return self._filename

    property header:
        '''the file header.
          
        .. note::
            The header is returned as an iterator presenting lines without the
            newline character.
        '''
        
        def __get__(self):
            return GZIteratorHead(self.filename)

    property contigs:
        '''list of chromosome names'''
        def __get__(self):
            cdef char ** sequences
            cdef int nsequences
            
            sequences = tbx_seqnames(self.index, &nsequences) 
            cdef int x
            result = []
            for x from 0 <= x < nsequences:
                result.append(sequences[x])
            
            # htslib instructions:
            # only free container, not the sequences themselves
            free(sequences)

            return result
            
    def close(self):
        '''
        closes the :class:`pysam.TabixFile`.'''
        if self.tabixfile != NULL:
            hts_close(self.tabixfile)
            self.tabixfile = NULL
        if self.index != NULL:
            tbx_destroy(self.index)
            self.index = NULL

    def __dealloc__( self ):
        # remember: dealloc cannot call other python methods
        # note: no doc string
        # note: __del__ is not called.
        if self.tabixfile != NULL:
            hts_close(self.tabixfile)
            self.tabixfile = NULL
        if self.index != NULL:
            tbx_destroy(self.index)


cdef class TabixIterator:
    """iterates over rows in *tabixfile* in region
    given by *tid*, *start* and *end*.
    """

    def __init__(self, encoding="ascii"):
        self.encoding = encoding
    
    def __iter__(self):
        self.buffer.s = NULL
        self.buffer.l = 0
        self.buffer.m = 0

        return self 

    cdef int __cnext__(self):
        '''iterate to next element.
        
        Return -5 if file has been closed when this function
        was called.
        '''
        if self.tabixfile.tabixfile == NULL:
            return -5

        cdef int retval

        while 1:
                
            retval = tbx_itr_next(
                self.tabixfile.tabixfile,
                self.tabixfile.index,
                self.iterator,
                &self.buffer)
            if retval < 0:
                break

            if self.buffer.s[0] != '#':
                break

        return retval

    def __next__(self): 
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """
        
        cdef int retval = self.__cnext__()
        if retval == -5:
            raise IOError("iteration on closed file")
        elif retval < 0:
            raise StopIteration

        return _charptr_to_str(self.buffer.s, self.encoding)

    def next(self):
        return self.__next__()

    def __dealloc__(self):
        if <void*>self.iterator != NULL:
            tbx_itr_destroy(self.iterator)
        if self.buffer.s != NULL:
            free(self.buffer.s)


class EmptyIterator:
    '''empty iterator'''

    def __iter__(self):
        return self

    def next(self):
        raise StopIteration()

    def __next__(self):
        raise StopIteration()


cdef class TabixIteratorParsed(TabixIterator):
    """iterates over mapped reads in a region.

    The *parser* determines the encoding.

    Returns parsed data.
    """

    def __init__(self, 
                 Parser parser):
        
        TabixIterator.__init__(self)
        self.parser = parser

    def __next__(self): 
        """python version of next().

        pyrex uses this non-standard name instead of next()
        """
        
        cdef int retval = self.__cnext__()
        if retval == -5:
            raise IOError("iteration on closed file")
        elif retval < 0:
            raise StopIteration

        return self.parser.parse(self.buffer.s,
                                 self.buffer.l)


cdef class GZIterator:
    def __init__(self, filename, int buffer_size=65536, encoding="ascii"):
        '''iterate line-by-line through gzip (or bgzip)
        compressed file.
        '''
        if not os.path.exists(filename):
            raise IOError("No such file or directory: %s" % filename)

        filename = _encodeFilename(filename)
        self.gzipfile = gzopen(filename, "r")
        self._filename = filename
        self.kstream = ks_init(self.gzipfile)
        self.encoding = encoding

        self.buffer.l = 0
        self.buffer.m = 0
        self.buffer.s = <char*>malloc(buffer_size)

    def __dealloc__(self):
        '''close file.'''
        if self.gzipfile != NULL:
            gzclose(self.gzipfile)
            self.gzipfile = NULL
        if self.buffer.s != NULL:
            free(self.buffer.s)
        ks_destroy(self.kstream)

    def __iter__(self):
        return self

    cdef int __cnext__(self):
        cdef int dret = 0
        cdef int retval = 0
        while 1:
            retval = ks_getuntil(self.kstream, '\n', &self.buffer, &dret)
            
            if retval < 0: 
                break

            return dret
        return -1

    def __next__(self):
        """python version of next().
        """
        cdef int retval = self.__cnext__()
        if retval < 0:
            raise StopIteration
        return _force_str(self.buffer.s, self.encoding)


cdef class GZIteratorHead(GZIterator):
    '''iterate line-by-line through gzip (or bgzip)
    compressed file returning comments at top of file.
    '''

    def __next__(self):
        """python version of next().
        """
        cdef int retval = self.__cnext__()
        if retval < 0:
            raise StopIteration
        if self.buffer.s[0] == '#':
            return self.buffer.s
        else:
            raise StopIteration


cdef class GZIteratorParsed(GZIterator):
    '''iterate line-by-line through gzip (or bgzip)
    compressed file returning comments at top of file.
    '''

    def __init__(self, parser):
        self.parser = parser

    def __next__(self):
        """python version of next().
        """
        cdef int retval = self.__cnext__()
        if retval < 0:
            raise StopIteration

        return self.parser.parse(self.buffer.s,
                                 self.buffer.l)


def tabix_compress(filename_in, 
                   filename_out,
                   force=False):
    '''compress *filename_in* writing the output to *filename_out*.
    
    Raise an IOError if *filename_out* already exists, unless *force*
    is set.
    '''

    if not force and os.path.exists(filename_out ):
        raise IOError(
            "Filename '%s' already exists, use *force* to overwrite" % filename_out)

    cdef int WINDOW_SIZE
    cdef int c, r
    cdef void * buffer
    cdef BGZF * fp
    cdef int fd_src

    cdef int O_RDONLY
    O_RDONLY = os.O_RDONLY

    WINDOW_SIZE = 64 * 1024

    fn = _encodeFilename(filename_out)
    fp = bgzf_open( fn, "w")
    if fp == NULL:
        raise IOError("could not open '%s' for writing" % (filename_out, ))

    fn = _encodeFilename(filename_in)
    fd_src = open(fn, O_RDONLY)
    if fd_src == 0:
        raise IOError("could not open '%s' for reading" % (filename_in, ))

    buffer = malloc(WINDOW_SIZE)
    c = 1

    while c > 0:
        c = read(fd_src, buffer, WINDOW_SIZE)
        r = bgzf_write(fp, buffer, c)
        if r < 0:
            free(buffer)
            raise OSError("writing failed")
        
    free(buffer)
    r = bgzf_close(fp)
    if r < 0:
        raise OSError("writing to file %s failed" % filename_out)

    r = close(fd_src)
    if r < 0:
        raise OSError("error when closing file %s" % filename_in)

def tabix_index( filename, 
                 force = False,
                 seq_col = None, 
                 start_col = None, 
                 end_col = None,
                 preset = None,
                 meta_char = "#",
                 zerobased = False,
                 min_shift = -1,
                ):
    '''index tab-separated *filename* using tabix.

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

    *min-shift* sets the minimal interval size to 1<<INT; 0 for the
    old tabix index. The default of -1 is changed inside htslib to 
    the old tabix default of 0.

    returns the filename of the compressed data

    '''
    
    if not os.path.exists(filename):
        raise IOError("No such file '%s'" % filename)

    if preset is None and \
       (seq_col is None or start_col is None or end_col is None):
        raise ValueError(
            "neither preset nor seq_col,start_col and end_col given")

    if not filename.endswith(".gz"): 
        tabix_compress(filename, filename + ".gz", force=force)
        os.unlink( filename )
        filename += ".gz"

    if not force and os.path.exists(filename + ".tbi"):
        raise IOError(
            "Filename '%s.tbi' already exists, use *force* to overwrite")

    # columns (1-based):
    #   preset-code, contig, start, end, metachar for
    #     comments, lines to ignore at beginning
    # 0 is a missing column
    preset2conf = {
        'gff' : (0, 1, 4, 5, ord('#'), 0),
        'bed' : (0x10000, 1, 2, 3, ord('#'), 0),
        'psltbl' : (0x10000, 15, 17, 18, ord('#'), 0),
        'sam' : (1, 3, 4, 0, ord('@'), 0),
        'vcf' : (2, 1, 2, 0, ord('#'), 0),
        'pileup': (3, 1, 2, 0, ord('#'), 0),
        }

    if preset:
        try:
            conf_data = preset2conf[preset]
        except KeyError:
            raise KeyError(
                "unknown preset '%s', valid presets are '%s'" %
                (preset, ",".join(preset2conf.keys())))
    else:
        if end_col == None:
            end_col = -1
        preset = 0

        # note that tabix internally works with 0-based coordinates
        # and open/closed intervals.  When using a preset, conversion
        # is automatically taken care of.  Otherwise, the coordinates
        # are assumed to be 1-based closed intervals and -1 is
        # subtracted from the start coordinate. To avoid doing this,
        # set the TI_FLAG_UCSC=0x10000 flag:
        if zerobased:
            preset = preset | 0x10000

        conf_data = (preset, seq_col+1, start_col+1, end_col+1, ord(meta_char), 0)
                
    cdef tbx_conf_t conf
    conf.preset, conf.sc, conf.bc, conf.ec, conf.meta_char, conf.line_skip = conf_data


    fn = _encodeFilename(filename)
    tbx_index_build(fn, min_shift, &conf)
    
    return filename

# #########################################################
# cdef class tabix_file_iterator_old:
#     '''iterate over ``infile``.

#     This iterator is not safe. If the :meth:`__next__()` method is called 
#     after ``infile`` is closed, the result is undefined (see ``fclose()``).

#     The iterator might either raise a StopIteration or segfault.
#     '''


#     def __cinit__(self, 
#                   infile, 
#                   Parser parser,
#                   int buffer_size = 65536 ):

#         cdef int fd = PyObject_AsFileDescriptor( infile )
#         if fd == -1: raise ValueError( "I/O operation on closed file." )
#         self.infile = fdopen( fd, 'r')

#         if self.infile == NULL: raise ValueError( "I/O operation on closed file." )

#         self.buffer = <char*>malloc( buffer_size )        
#         self.size = buffer_size
#         self.parser = parser

#     def __iter__(self):
#         return self

#     cdef __cnext__(self):

#         cdef char * b
#         cdef size_t nbytes
#         b = self.buffer

#         while not feof( self.infile ):
#             nbytes = getline( &b, &self.size, self.infile)

#             # stop at first error or eof
#             if (nbytes == -1): break
#             # skip comments
#             if (b[0] == '#'): continue

#             # skip empty lines
#             if b[0] == '\0' or b[0] == '\n' or b[0] == '\r': continue

#             # make sure that entry is complete
#             if b[nbytes-1] != '\n' and b[nbytes-1] != '\r':
#                 result = b
#                 raise ValueError( "incomplete line at %s" % result )

#             # make sure that this goes fully through C
#             # otherwise buffer is copied to/from a
#             # Python object causing segfaults as
#             # the wrong memory is freed
#             return self.parser.parse( b, nbytes )

#         raise StopIteration

#     def __dealloc__(self):
#         free(self.buffer)

#     def __next__(self):
#         return self.__cnext__()

#########################################################
#########################################################
#########################################################
## Iterators for parsing through unindexed files.
#########################################################
cdef buildGzipError(void *gzfp):
    cdef int errnum = 0
    cdef char *s = gzerror(gzfp, &errnum)
    return "error (%d): %s (%d: %s)" % (errno, strerror(errno), errnum, s)


cdef class tabix_file_iterator:
    '''iterate over a compressed or uncompressed ``infile``.
    '''

    def __cinit__(self, 
                  infile, 
                  Parser parser,
                  int buffer_size=65536):

        if infile.closed:
            raise ValueError("I/O operation on closed file.")

        self.infile = infile

        cdef int fd = PyObject_AsFileDescriptor(infile)
        if fd == -1:
            raise ValueError("I/O operation on closed file.")

        self.duplicated_fd = dup(fd)

        # From the manual:
        # gzopen can be used to read a file which is not in gzip format; 
        # in this case gzread will directly read from the file without decompression. 
        # When reading, this will be detected automatically by looking 
        # for the magic two-byte gzip header. 
        self.fh = gzdopen(self.duplicated_fd, 'r')

        if self.fh == NULL: 
            raise IOError('%s' % strerror(errno))

        self.kstream = ks_init(self.fh) 
        
        self.buffer.s = <char*>malloc(buffer_size)
        #if self.buffer == NULL:
        #    raise MemoryError( "tabix_file_iterator: could not allocate %i bytes" % buffer_size)
        #self.size = buffer_size
        self.parser = parser

    def __iter__(self):
        return self

    cdef __cnext__(self):

        cdef char * b
        cdef int dret = 0
        cdef int retval = 0
        while 1:
            
            retval = ks_getuntil(self.kstream, '\n', &self.buffer, &dret)
            
            if retval < 0: 
                break
                #raise IOError('gzip error: %s' % buildGzipError( self.fh ))

            b = self.buffer.s
            
            # skip comments
            if (b[0] == '#'):
                continue

            # skip empty lines
            if b[0] == '\0' or b[0] == '\n' or b[0] == '\r':
                continue

            # gzgets terminates at \n, no need to test

            # parser creates a copy
            return self.parser.parse( b, self.buffer.l)

        raise StopIteration

    def __dealloc__(self):
        free(self.buffer.s)
        ks_destroy(self.kstream)
        gzclose(self.fh)
        
    def __next__(self):
        return self.__cnext__()

    def next(self):
        return self.__cnext__()
    

class tabix_generic_iterator:
    '''iterate over ``infile``.
    
    Permits the use of file-like objects for example from the gzip module.
    '''
    def __init__(self, infile, parser):

        self.infile = infile
        if self.infile.closed:
            raise ValueError("I/O operation on closed file.")
        self.parser = parser

    def __iter__(self):
        return self

    # cython version - required for python 3
    def __next__(self):
        
        cdef char * b
        cdef char * cpy
        cdef size_t nbytes

        encoding = self.parser.get_encoding()

        # note that GzipFile.close() does not close the file
        # reading is still possible.
        if self.infile.closed:
            raise ValueError("I/O operation on closed file.")

        while 1:

            line = self.infile.readline()
            if not line:
                break
            
            s = _force_bytes(line, encoding)
            b = s
            nbytes = len(line)
            assert b[nbytes] == '\0'

            # skip comments
            if b[0] == '#':
                continue

            # skip empty lines
            if b[0] == '\0' or b[0] == '\n' or b[0] == '\r':
                continue
            
            # make sure that entry is complete
            if b[nbytes-1] != '\n' and b[nbytes-1] != '\r':
                raise ValueError("incomplete line at %s" % line)
            
            bytes_cpy = <bytes> b
            cpy = <char *> bytes_cpy

            return self.parser(cpy, nbytes)            

        raise StopIteration

    # python version - required for python 2.7
    def next(self):
        return self.__next__()

def tabix_iterator(infile, parser):
    """return an iterator over all entries in a file.
    
    Results are returned parsed as specified by the *parser*. If
    *parser* is None, the results are returned as an unparsed string.
    Otherwise, *parser* is assumed to be a functor that will return
    parsed data (see for example :class:`~pysam.asTuple` and
    :class:`~pysam.asGTF`).

    """
    if PYTHON3:
        return tabix_generic_iterator(infile, parser)
    else:
        return tabix_file_iterator(infile, parser)
        
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
    
cdef class Tabixfile(TabixFile):
    pass


__all__ = [
    "tabix_index", 
    "tabix_compress",
    "TabixFile",
    "Tabixfile",
    "asTuple",
    "asGTF",
    "asVCF",
    "asBed",
    "GZIterator",
    "GZIteratorHead",
    "tabix_iterator", 
    "tabix_generic_iterator", 
    "tabix_file_iterator", 
]
