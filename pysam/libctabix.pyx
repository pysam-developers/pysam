# cython: language_level=3
# cython: embedsignature=True
# cython: profile=True
###############################################################################
###############################################################################
# Cython wrapper for access to tabix indexed files in bgzf format
###############################################################################
# The principal classes and functions defined in this module are:
#
# class TabixFile  class wrapping tabix indexed files in bgzf format
#
# class asTuple  Parser class for tuples
# class asGTF    Parser class for GTF formatted rows
# class asGFF3   Parser class for GFF3 formatted rows
# class asBed    Parser class for Bed formatted rows
# class asVCF    Parser class for VCF formatted rows
#
# class tabix_generic_iterator  Streamed iterator of bgzf formatted files
#
# Additionally this module defines several additional classes that are part
# of the internal API. These are:
#
# class Parser                base class for parsers of tab-separated rows
# class tabix_file_iterator
# class TabixIterator         iterator class over rows in bgzf file
# class EmptyIterator
#
# For backwards compatibility, the following classes are also defined:
#
# class Tabixfile   equivalent to TabixFile
#
###############################################################################
#
# The MIT License
#
# Copyright (c) 2015 Andreas Heger
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.
#
###############################################################################
import os
import sys

from libc.stdio cimport printf, fprintf, stderr
from libc.string cimport strerror
from libc.errno cimport errno
from posix.unistd cimport dup

from cpython cimport PyErr_SetString, PyBytes_Check, \
    PyUnicode_Check, PyBytes_FromStringAndSize, \
    PyObject_AsFileDescriptor

cimport pysam.libctabixproxies as ctabixproxies

from pysam.libchtslib cimport htsFile, hts_open, hts_close, HTS_IDX_START,\
    BGZF, bgzf_open, bgzf_dopen, bgzf_close, bgzf_write, \
    tbx_index_build2, tbx_index_load2, tbx_itr_queryi, tbx_itr_querys, \
    tbx_conf_t, tbx_seqnames, tbx_itr_next, tbx_itr_destroy, \
    tbx_destroy, hisremote, region_list, hts_getline, \
    TBX_GENERIC, TBX_SAM, TBX_VCF, TBX_UCSC, hts_get_format, htsFormat, \
    no_compression, bcf, bcf_index_build2

from pysam.libcutils cimport force_bytes, force_str, charptr_to_str
from pysam.libcutils cimport encode_filename, from_string_and_size

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
        cdef ctabixproxies.TupleProxy r
        r = ctabixproxies.TupleProxy(self.encoding)
        # need to copy - there were some
        # persistence issues with "present"
        r.copy(buffer, len)
        return r


cdef class asGFF3(Parser):
    '''converts a :term:`tabix row` into a GFF record with the following
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

    ''' 
    cdef parse(self, char * buffer, int len):
        cdef ctabixproxies.GFF3Proxy r
        r = ctabixproxies.GFF3Proxy(self.encoding)
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
        cdef ctabixproxies.GTFProxy r
        r = ctabixproxies.GTFProxy(self.encoding)
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
    fields are optional, but if one is defined, all the preceding
    need to be defined as well.

    ''' 
    cdef parse(self, char * buffer, int len):
        cdef ctabixproxies.BedProxy r
        r = ctabixproxies.BedProxy(self.encoding)
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
        cdef ctabixproxies.VCFProxy r
        r = ctabixproxies.VCFProxy(self.encoding)
        r.copy(buffer, len)
        return r


cdef class TabixFile:
    """Random access to bgzf formatted files that
    have been indexed by :term:`tabix`.

    The file is automatically opened. The index file of file
    ``<filename>`` is expected to be called ``<filename>.tbi``
    by default (see parameter `index`).
    
    Parameters
    ----------
    
    filename : string
        Filename of bgzf file to be opened.

    index : string
        The filename of the index. If not set, the default is to
        assume that the index is called ``filename.tbi``

    mode : char
        The file opening mode. Currently, only ``r`` is permitted.
        
    parser : :class:`pysam.Parser`
    
        sets the default parser for this tabix file. If `parser`
        is None, the results are returned as an unparsed string.
        Otherwise, `parser` is assumed to be a functor that will return
        parsed data (see for example :class:`~pysam.asTuple` and
        :class:`~pysam.asGTF`).

    encoding : string

        The encoding passed to the parser

    threads: integer
        Number of threads to use for decompressing Tabix files.
        (Default=1)


    Raises
    ------
    
    ValueError
        if index file is missing.

    IOError
        if file could not be opened
    """
    def __cinit__(self,
                  filename,
                  mode='r',
                  parser=None,
                  index=None,
                  encoding="ascii",
                  threads=1,
                  *args,
                  **kwargs ):

        self.htsfile = NULL
        self.is_remote = False
        self.is_stream = False
        self.parser = parser
        self.threads = threads
        self._open(filename, mode, index, *args, **kwargs)
        self.encoding = encoding

    def _open( self,
               filename,
               mode='r',
               index=None,
               threads=1,
              ):
        '''open a :term:`tabix file` for reading.'''

        if mode != 'r':
            raise ValueError("invalid file opening mode `%s`" % mode)

        if self.htsfile != NULL:
            self.close()
        self.htsfile = NULL
        self.threads=threads

        filename_index = index or (filename + ".tbi")
        # encode all the strings to pass to tabix
        self.filename = encode_filename(filename)
        self.filename_index = encode_filename(filename_index)

        self.is_stream = self.filename == b'-'
        self.is_remote = hisremote(self.filename)

        if not self.is_remote:
            if not os.path.exists(filename):
                raise IOError("file `%s` not found" % filename)

            if not os.path.exists(filename_index):
                raise IOError("index `%s` not found" % filename_index)

        # open file
        cdef char *cfilename = self.filename
        cdef char *cfilename_index = self.filename_index
        with nogil:
            self.htsfile = hts_open(cfilename, 'r')

        if self.htsfile == NULL:
            raise IOError("could not open file `%s`" % filename)
        
        #if self.htsfile.format.category != region_list:
        #    raise ValueError("file does not contain region data")

        with nogil:
            self.index = tbx_index_load2(cfilename, cfilename_index)

        if self.index == NULL:
            raise IOError("could not open index for `%s`" % filename)

        if not self.is_stream:
            self.start_offset = self.tell()

    def _dup(self):
        '''return a copy of this tabix file.
        
        The file is being re-opened.
        '''
        return TabixFile(self.filename,
                         mode="r",
                         threads=self.threads,
                         parser=self.parser,
                         index=self.filename_index,
                         encoding=self.encoding)

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
        if not self.is_open():
            raise ValueError("I/O operation on closed file")

        # convert coordinates to region string, which is one-based
        if reference:
            if end is not None:
                if end < 0:
                    raise ValueError("end out of range (%i)" % end)
                if start is None:
                    start = 0
                    
                if start < 0:
                    raise ValueError("start out of range (%i)" % end)
                elif start > end:
                    raise ValueError(
                        'start (%i) >= end (%i)' % (start, end))
                elif start == end:
                    return EmptyIterator()
                else:
                    region = '%s:%i-%i' % (reference, start + 1, end)
            elif start is not None:
                if start < 0:
                    raise ValueError("start out of range (%i)" % end)
                region = '%s:%i' % (reference, start + 1)
            else:
                region = reference

        # get iterator
        cdef hts_itr_t * itr
        cdef char *cstr
        cdef TabixFile fileobj

        # reopen the same file if necessary
        if multiple_iterators:
            fileobj = self._dup()
        else:
            fileobj = self

        if region is None:
            # without region or reference - iterate from start
            with nogil:
                itr = tbx_itr_queryi(fileobj.index,
                                     HTS_IDX_START,
                                     0,
                                     0)
        else:
            s = force_bytes(region, encoding=fileobj.encoding)
            cstr = s
            with nogil:
                itr = tbx_itr_querys(fileobj.index, cstr)

        if itr == NULL:
            if region is None:
                if len(self.contigs) > 0:
                    # when accessing a tabix file created prior tabix 1.0
                    # the full-file iterator is empty.
                    raise ValueError(
                        "could not create iterator, possible "
                        "tabix version mismatch")
                else:
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
        a.iterator = itr

        return a

    ###############################################################
    ###############################################################
    ###############################################################
    ## properties
    ###############################################################
    property header:
        '''the file header.

        The file header consists of the lines at the beginning of a
        file that are prefixed by the comment character ``#``.
       
        .. note::
            The header is returned as an iterator presenting lines
            without the newline character.
        '''
        
        def __get__(self):

            cdef char *cfilename = self.filename
            cdef char *cfilename_index = self.filename_index
            
            cdef kstring_t buffer
            buffer.l = buffer.m = 0
            buffer.s = NULL
            
            cdef htsFile * fp = NULL
            cdef int KS_SEP_LINE = 2
            cdef tbx_t * tbx = NULL
            lines = []
            with nogil:
                fp = hts_open(cfilename, 'r')
                
            if fp == NULL:
                raise OSError("could not open {} for reading header".format(self.filename))

            with nogil:
                tbx = tbx_index_load2(cfilename, cfilename_index)
                
            if tbx == NULL:
                raise OSError("could not load .tbi/.csi index of {}".format(self.filename))

            while hts_getline(fp, KS_SEP_LINE, &buffer) >= 0:
                if not buffer.l or buffer.s[0] != tbx.conf.meta_char:
                    break
                lines.append(force_str(buffer.s, self.encoding))

            with nogil:
                hts_close(fp)
                free(buffer.s)

            return lines

    property contigs:
        '''list of chromosome names'''
        def __get__(self):
            cdef const char ** sequences
            cdef int nsequences
            
            with nogil:
                sequences = tbx_seqnames(self.index, &nsequences)
            cdef int x
            result = []
            for x from 0 <= x < nsequences:
                result.append(force_str(sequences[x]))
            
            # htslib instructions:
            # only free container, not the sequences themselves
            free(sequences)

            return result
            
    def close(self):
        '''
        closes the :class:`pysam.TabixFile`.'''
        if self.htsfile != NULL:
            hts_close(self.htsfile)
            self.htsfile = NULL
        if self.index != NULL:
            tbx_destroy(self.index)
            self.index = NULL

    def __dealloc__( self ):
        # remember: dealloc cannot call other python methods
        # note: no doc string
        # note: __del__ is not called.
        if self.htsfile != NULL:
            hts_close(self.htsfile)
            self.htsfile = NULL
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
        if self.tabixfile.htsfile == NULL:
            return -5

        cdef int retval

        while 1:
            with nogil:
                retval = tbx_itr_next(
                    self.tabixfile.htsfile,
                    self.tabixfile.index,
                    self.iterator,
                    &self.buffer)

            if retval < 0:
                break

            if self.buffer.s[0] != b'#':
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

        return charptr_to_str(self.buffer.s, self.encoding)

    def __dealloc__(self):
        if <void*>self.iterator != NULL:
            tbx_itr_destroy(self.iterator)
        if self.buffer.s != NULL:
            free(self.buffer.s)


class EmptyIterator:
    '''empty iterator'''

    def __iter__(self):
        return self

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

        filename = encode_filename(filename)
        cdef char *cfilename = filename
        with nogil:
            self.gzipfile = bgzf_open(cfilename, "r")
        self._filename = filename
        self.kstream = ks_init(self.gzipfile)
        self.encoding = encoding

        self.buffer.l = 0
        self.buffer.m = 0
        self.buffer.s = <char*>malloc(buffer_size)

    def __dealloc__(self):
        '''close file.'''
        if self.gzipfile != NULL:
            bgzf_close(self.gzipfile)
            self.gzipfile = NULL
        if self.buffer.s != NULL:
            free(self.buffer.s)
        if self.kstream != NULL:
            ks_destroy(self.kstream)

    def __iter__(self):
        return self

    cdef int __cnext__(self):
        cdef int dret = 0
        cdef int retval = 0
        while 1:
            with nogil:
                retval = ks_getuntil(self.kstream, b'\n', &self.buffer, &dret)
            
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
        return force_str(self.buffer.s, self.encoding)


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
        if self.buffer.s[0] == b'#':
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

    if not force and os.path.exists(filename_out):
        raise IOError(
            "Filename '%s' already exists, use *force* to "
            "overwrite" % filename_out)

    cdef int WINDOW_SIZE
    cdef int c, r
    cdef void * buffer
    cdef BGZF * fp
    cdef int fd_src
    cdef bint is_empty = True
    cdef int O_RDONLY
    O_RDONLY = os.O_RDONLY

    WINDOW_SIZE = 64 * 1024

    fn = encode_filename(filename_out)
    cdef char *cfn = fn
    with nogil:
        fp = bgzf_open(cfn, "w")
    if fp == NULL:
        raise IOError("could not open '%s' for writing" % filename_out)

    fn = encode_filename(filename_in)
    fd_src = open(fn, O_RDONLY)
    if fd_src == 0:
        raise IOError("could not open '%s' for reading" % filename_in)

    buffer = malloc(WINDOW_SIZE)
    c = 1
    
    while c > 0:
        with nogil:
            c = read(fd_src, buffer, WINDOW_SIZE)
            if c > 0:
                is_empty = False
            r = bgzf_write(fp, buffer, c)
        if r < 0:
            free(buffer)
            raise IOError("writing failed")
        
    free(buffer)
    r = bgzf_close(fp)
    if r < 0:
        raise IOError("error %i when writing to file %s" % (r, filename_out))

    r = close(fd_src)
    # an empty file will return with -1, thus ignore this.
    if r < 0:
        if not (r == -1 and is_empty):
            raise IOError("error %i when closing file %s" % (r, filename_in))


def tabix_index(filename,
                force=False,
                seq_col=None,
                start_col=None,
                end_col=None,
                preset=None,
                meta_char="#",
                int line_skip=0,
                zerobased=False,
                int min_shift=-1,
                index=None,
                keep_original=False,
                csi=False,
                ):
    '''index tab-separated *filename* using tabix.

    An existing index will not be overwritten unless *force* is set.

    The index will be built from coordinates in columns *seq_col*,
    *start_col* and *end_col*.

    The contents of *filename* have to be sorted by contig and
    position - the method does not check if the file is sorted.

    Column indices are 0-based. Note that this is different from the
    tabix command line utility where column indices start at 1.
    
    Coordinates in the file are assumed to be 1-based unless
    *zerobased* is set.

    If *preset* is provided, the column coordinates are taken from a
    preset. Valid values for preset are "gff", "bed", "sam", "vcf",
    psltbl", "pileup".
    
    Lines beginning with *meta_char* and the first *line_skip* lines
    will be skipped.

    If *filename* is not detected as a gzip file it will be automatically
    compressed. The original file will be removed and only the compressed
    file will be retained.

    By default or when *min_shift* is 0, creates a TBI index. If *min_shift*
    is greater than zero and/or *csi* is True, creates a CSI index with a
    minimal interval size of 1<<*min_shift* (1<<14 if only *csi* is set).

    *index* controls the filename which should be used for creating the index.
    If not set, the default is to append ``.tbi`` to *filename*.

    When automatically compressing files, if *keep_original* is set the
    uncompressed file will not be deleted.

    returns the filename of the compressed data

    '''
    
    if preset is None and \
       (seq_col is None or start_col is None or end_col is None):
        raise ValueError(
            "neither preset nor seq_col,start_col and end_col given")

    fn = encode_filename(filename)
    cdef char *cfn = fn

    cdef htsFile *fp = hts_open(cfn, "r")
    if fp == NULL:
        raise IOError("Could not open file '%s': %s" % (filename, force_str(strerror(errno))))

    cdef htsFormat fmt = hts_get_format(fp)[0]
    hts_close(fp)

    if fmt.compression == no_compression:
        tabix_compress(filename, filename + ".gz", force=force)
        if not keep_original:
            os.unlink(filename)
        filename += ".gz"
        fn = encode_filename(filename)
        cfn = fn

    # columns (1-based):
    #   preset-code, contig, start, end, metachar for
    #     comments, lines to ignore at beginning
    # 0 is a missing column
    preset2conf = {
        'gff' : (TBX_GENERIC, 1, 4, 5, ord('#'), 0),
        'bed' : (TBX_UCSC, 1, 2, 3, ord('#'), 0),
        'psltbl' : (TBX_UCSC, 15, 17, 18, ord('#'), 0),
        'sam' : (TBX_SAM, 3, 4, 0, ord('@'), 0),
        'vcf' : (TBX_VCF, 1, 2, 0, ord('#'), 0),
        }
    
    conf_data = None
    if preset == "bcf" or fmt.format == bcf:
        csi = True
    elif preset:
        try:
            conf_data = preset2conf[preset]
        except KeyError:
            raise KeyError(
                "unknown preset '%s', valid presets are '%s'" %
                (preset, ",".join(preset2conf.keys())))
    else:
        if end_col is None:
            end_col = -1
            
        preset = 0
        # tabix internally works with 0-based coordinates and
        # open/closed intervals.  When using a preset, conversion is
        # automatically taken care of.  Otherwise, the coordinates are
        # assumed to be 1-based closed intervals and -1 is subtracted
        # from the start coordinate. To avoid doing this, set the
        # TI_FLAG_UCSC=0x10000 flag:
        if zerobased:
            preset = preset | TBX_UCSC

        conf_data = (preset, seq_col + 1, start_col + 1, end_col + 1, ord(meta_char), line_skip)

    cdef tbx_conf_t conf
    if conf_data:
        conf.preset, conf.sc, conf.bc, conf.ec, conf.meta_char, conf.line_skip = conf_data

    if csi or min_shift > 0:
        suffix = ".csi"
        if min_shift <= 0: min_shift = 14
    else:
        suffix = ".tbi"
        min_shift = 0

    index = index or filename + suffix    
    fn_index = encode_filename(index)

    if not force and os.path.exists(index):
        raise IOError(
            "filename '%s' already exists, use *force* to overwrite" % index)
    
    cdef char *fnidx = fn_index
    cdef int retval = 0

    if csi and fmt.format == bcf:
        with nogil:
            retval = bcf_index_build2(cfn, fnidx, min_shift)
    else:
        with nogil:
            retval = tbx_index_build2(cfn, fnidx, min_shift, &conf)
            
    if retval != 0:
        raise OSError("building of index for {} failed".format(filename))
    
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
# cdef buildGzipError(void *gzfp):
#     cdef int errnum = 0
#     cdef char *s = gzerror(gzfp, &errnum)
#     return "error (%d): %s (%d: %s)" % (errno, strerror(errno), errnum, s)


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
        self.fh = bgzf_dopen(self.duplicated_fd, 'r')

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
            with nogil:
                retval = ks_getuntil(self.kstream, b'\n', &self.buffer, &dret)
            
            if retval < 0: 
                break
                #raise IOError('gzip error: %s' % buildGzipError( self.fh ))

            b = self.buffer.s
            
            # skip comments
            if (b[0] == b'#'):
                continue

            # skip empty lines
            if b[0] == b'\0' or b[0] == b'\n' or b[0] == b'\r':
                continue

            # gzgets terminates at \n, no need to test

            # parser creates a copy
            return self.parser.parse(b, self.buffer.l)

        raise StopIteration

    def __dealloc__(self):
        free(self.buffer.s)
        ks_destroy(self.kstream)
        bgzf_close(self.fh)
        
    def __next__(self):
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
            
            s = force_bytes(line, encoding)
            b = s
            nbytes = len(line)
            assert b[nbytes] == b'\0'

            # skip comments
            if b[0] == b'#':
                continue

            # skip empty lines
            if b[0] == b'\0' or b[0] == b'\n' or b[0] == b'\r':
                continue
            
            # make sure that entry is complete
            if b[nbytes-1] != b'\n' and b[nbytes-1] != b'\r':
                raise ValueError("incomplete line at %s" % line)
            
            bytes_cpy = <bytes> b
            cpy = <char *> bytes_cpy

            return self.parser(cpy, nbytes)            

        raise StopIteration


def tabix_iterator(infile, parser):
    """return an iterator over all entries in a file.
    
    Results are returned parsed as specified by the *parser*. If
    *parser* is None, the results are returned as an unparsed string.
    Otherwise, *parser* is assumed to be a functor that will return
    parsed data (see for example :class:`~pysam.asTuple` and
    :class:`~pysam.asGTF`).

    """
    return tabix_generic_iterator(infile, parser)
        
    
cdef class Tabixfile(TabixFile):
    """Tabixfile is deprecated: use TabixFile instead"""
    pass


__all__ = [
    "tabix_index", 
    "tabix_compress",
    "TabixFile",
    "Tabixfile",
    "asTuple",
    "asGTF",
    "asGFF3",
    "asVCF",
    "asBed",
    "GZIterator",
    "GZIteratorHead",
    "tabix_iterator", 
    "tabix_generic_iterator", 
    "tabix_file_iterator", 
]
