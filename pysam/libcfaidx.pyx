# cython: embedsignature=True
# cython: profile=True
###############################################################################
###############################################################################
# Cython wrapper for SAM/BAM/CRAM files based on htslib
###############################################################################
# The principal classes defined in this module are:
#
# class FastaFile   random read read/write access to faidx indexd files
# class FastxFile   streamed read/write access to fasta/fastq files
#
# Additionally this module defines several additional classes that are part
# of the internal API. These are:
#
# class FastqProxy
# class FastxRecord
#
# For backwards compatibility, the following classes are also defined:
#
# class Fastafile   equivalent to FastaFile
# class FastqFile   equivalent to FastxFile
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
import sys
import os
import re


from libc.errno  cimport errno
from libc.string cimport strerror

from cpython cimport array

from cpython cimport PyErr_SetString, \
    PyBytes_Check, \
    PyUnicode_Check, \
    PyBytes_FromStringAndSize

from pysam.libchtslib cimport \
    faidx_nseq, fai_load, fai_load3, fai_destroy, fai_fetch, \
    faidx_seq_len, faidx_iseq, faidx_seq_len, \
    faidx_fetch_seq, hisremote, \
    bgzf_open, bgzf_close

from pysam.libcutils cimport force_bytes, force_str, charptr_to_str
from pysam.libcutils cimport encode_filename, from_string_and_size
from pysam.libcutils cimport qualitystring_to_array, parse_region

cdef class FastqProxy
cdef makeFastqProxy(kseq_t * src):
    '''enter src into AlignedRead.'''
    cdef FastqProxy dest = FastqProxy.__new__(FastqProxy)
    dest._delegate = src
    return dest

## TODO:
##        add automatic indexing.
##        add function to get sequence names.
cdef class FastaFile:
    """Random access to fasta formatted files that
    have been indexed by :term:`faidx`.

    The file is automatically opened. The index file of file
    ``<filename>`` is expected to be called ``<filename>.fai``.

    Parameters
    ----------

    filename : string
        Filename of fasta file to be opened.

    filepath_index : string
        Optional, filename of the index. By default this is
        the filename + ".fai".

    filepath_index_compressed : string
        Optional, filename of the index if fasta file is. By default this is
        the filename + ".gzi".

    Raises
    ------

    ValueError
        if index file is missing

    IOError
        if file could not be opened

    """

    def __cinit__(self, *args, **kwargs):
        self.fastafile = NULL
        self._filename = None
        self._references = None
        self._lengths = None
        self.reference2length = None
        self._open(*args, **kwargs)

    def is_open(self):
        '''return true if samfile has been opened.'''
        return self.fastafile != NULL

    def __len__(self):
        if self.fastafile == NULL:
            raise ValueError("calling len() on closed file")

        return faidx_nseq(self.fastafile)

    def _open(self, filename, filepath_index=None, filepath_index_compressed=None):
        '''open an indexed fasta file.

        This method expects an indexed fasta file.
        '''

        # close a previously opened file
        if self.fastafile != NULL:
            self.close()

        self._filename = encode_filename(filename)
        cdef char *cfilename = self._filename
        cdef char *cindexname = NULL
        cdef char *cindexname_compressed = NULL
        self.is_remote = hisremote(cfilename)
        
        # open file for reading
        if (self._filename != b"-"
            and not self.is_remote
            and not os.path.exists(filename)):
            raise IOError("file `%s` not found" % filename)

        # 3 modes to open:
        # compressed fa: fai_load3 with filename, index_fai and index_gzi
        # uncompressed fa: fai_load3 with filename and index_fai
        # uncompressed fa: fai_load with default index name
        if filepath_index:
            # when opening, set flags to 0 - do not automatically
            # build index if it does not exist.

            if not os.path.exists(filepath_index):
                raise IOError("filename {} does not exist".format(filepath_index))
            cindexname = bindex_filename = encode_filename(filepath_index)
            
            if filepath_index_compressed:
                if not os.path.exists(filepath_index_compressed):
                    raise IOError("filename {} does not exist".format(filepath_index_compressed))
                cindexname_compressed = bindex_filename_compressed = encode_filename(filepath_index_compressed)
                with nogil:
                    self.fastafile = fai_load3(cfilename, cindexname, cindexname_compressed, 0)
            else:
                with nogil:
                    self.fastafile = fai_load3(cfilename, cindexname, NULL, 0)
        else:
            with nogil:
                self.fastafile = fai_load(cfilename)

        if self.fastafile == NULL:
            raise IOError("error when opening file `%s`" % filename)

        cdef int nreferences = faidx_nseq(self.fastafile)
        cdef int x
        cdef const char * s
        self._references = []
        self._lengths = []
        for x from 0 <= x < nreferences:
            s = faidx_iseq(self.fastafile, x)
            ss = force_str(s)
            self._references.append(ss)
            self._lengths.append(faidx_seq_len(self.fastafile, s))
        self.reference2length = dict(zip(self._references, self._lengths))

    def close(self):
        """close the file."""
        if self.fastafile != NULL:
            fai_destroy(self.fastafile)
            self.fastafile = NULL

    def __dealloc__(self):
        if self.fastafile != NULL:
            fai_destroy(self.fastafile)
            self.fastafile = NULL

    # context manager interface
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    property closed:
        """bool indicating the current state of the file object.
        This is a read-only attribute; the close() method changes the value.
        """
        def __get__(self):
            return not self.is_open()

    property filename:
        """filename associated with this object. This is a read-only attribute."""
        def __get__(self):
            return self._filename

    property references:
        '''tuple with the names of :term:`reference` sequences.'''
        def __get__(self):
            return self._references

    property nreferences:
        """int with the number of :term:`reference` sequences in the file.
        This is a read-only attribute."""
        def __get__(self):
            return len(self._references) if self.references else None

    property lengths:
        """tuple with the lengths of :term:`reference` sequences."""
        def __get__(self):
            return self._lengths

    def fetch(self,
              reference=None,
              start=None,
              end=None,
              region=None):
        """fetch sequences in a :term:`region`.

        A region can
        either be specified by :term:`reference`, `start` and
        `end`. `start` and `end` denote 0-based, half-open
        intervals.

        Alternatively, a samtools :term:`region` string can be
        supplied.

        If any of the coordinates are missing they will be replaced by the
        minimum (`start`) or maximum (`end`) coordinate.

        Note that region strings are 1-based, while `start` and `end` denote
        an interval in python coordinates.
        The region is specified by :term:`reference`, `start` and `end`.

        Returns
        -------

        string : a string with the sequence specified by the region.

        Raises
        ------

        IndexError
            if the coordinates are out of range

        ValueError
            if the region is invalid

        """

        if not self.is_open():
            raise ValueError("I/O operation on closed file" )

        cdef int length
        cdef char *seq
        cdef char *ref
        cdef int rstart, rend

        contig, rstart, rend = parse_region(reference, start, end, region)

        if contig is None:
            raise ValueError("no sequence/region supplied.")

        if rstart == rend:
            return ""

        contig_b = force_bytes(contig)
        ref = contig_b
        with nogil:
            length = faidx_seq_len(self.fastafile, ref)
        if length == -1:
            raise KeyError("sequence '%s' not present" % contig)
        if rstart >= length:
            return ""

        # fai_fetch adds a '\0' at the end
        with nogil:
            seq = faidx_fetch_seq(self.fastafile,
                                  ref,
                                  rstart,
                                  rend-1,
                                  &length)

        if not seq:
            if errno:
                raise IOError(errno, strerror(errno))
            else:
                raise ValueError("failure when retrieving sequence on '%s'" % contig)

        try:
            return charptr_to_str(seq)
        finally:
            free(seq)

    cdef char *_fetch(self, char *reference, int start, int end, int *length) except? NULL:
        '''fetch sequence for reference, start and end'''

        cdef char *seq
        with nogil:
            seq = faidx_fetch_seq(self.fastafile,
                                  reference,
                                  start,
                                  end-1,
                                  length)

        if not seq:
            if errno:
                raise IOError(errno, strerror(errno))
            else:
                raise ValueError("failure when retrieving sequence on '%s'" % reference)

        return seq

    def get_reference_length(self, reference):
        '''return the length of reference.'''
        return self.reference2length[reference]

    def __getitem__(self, reference):
        return self.fetch(reference)

    def __contains__(self, reference):
        '''return true if reference in fasta file.'''
        return reference in self.reference2length


cdef class FastqProxy:
    """A single entry in a fastq file."""
    def __init__(self):
        raise ValueError("do not instantiate FastqProxy directly")

    property name:
        """The name of each entry in the fastq file."""
        def __get__(self):
            return charptr_to_str(self._delegate.name.s)

    property sequence:
        """The sequence of each entry in the fastq file."""
        def __get__(self):
            return charptr_to_str(self._delegate.seq.s)

    property comment:
        def __get__(self):
            if self._delegate.comment.l:
                return charptr_to_str(self._delegate.comment.s)
            else:
                return None

    property quality:
        """The quality score of each entry in the fastq file, represented as a string."""
        def __get__(self):
            if self._delegate.qual.l:
                return charptr_to_str(self._delegate.qual.s)
            else:
                return None

    cdef cython.str to_string(self):
        if self.comment is None:
            comment = ""
        else:
            comment = " %s" % self.comment

        if self.quality is None:
            return ">%s%s\n%s" % (self.name, comment, self.sequence)
        else:
            return "@%s%s\n%s\n+\n%s" % (self.name, comment,
                                         self.sequence, self.quality)
    
    cdef cython.str tostring(self):
        """deprecated : use :meth:`to_string`"""
        return self.to_string()
    
    def __str__(self):
        return self.to_string()

    cpdef array.array get_quality_array(self, int offset=33):
        '''return quality values as integer array after subtracting offset.'''
        if self.quality is None:
            return None
        return qualitystring_to_array(force_bytes(self.quality),
                                      offset=offset)

cdef class FastxRecord:
    """A fasta/fastq record.

    A record must contain a name and a sequence. If either of them are
    None, a ValueError is raised on writing.

    """
    def __init__(self,
                 name=None,
                 comment=None,
                 sequence=None,
                 quality=None,
                 FastqProxy proxy=None):
        if proxy is not None:
            self.comment = proxy.comment
            self.quality = proxy.quality
            self.sequence = proxy.sequence
            self.name = proxy.name
        else:
            self.comment = comment
            self.quality = quality
            self.sequence = sequence
            self.name = name

    def __copy__(self):
        return FastxRecord(self.name, self.comment, self.sequence, self.quality)

    def __deepcopy__(self, memo):
        return FastxRecord(self.name, self.comment, self.sequence, self.quality)

    cdef cython.str to_string(self):
        if self.name is None:
            raise ValueError("can not write record without name")

        if self.sequence is None:
            raise ValueError("can not write record without a sequence")
        
        if self.comment is None:
            comment = ""
        else:
            comment = " %s" % self.comment

        if self.quality is None:
            return ">%s%s\n%s" % (self.name, comment, self.sequence)
        else:
            return "@%s%s\n%s\n+\n%s" % (self.name, comment,
                                         self.sequence, self.quality)
        
    cdef cython.str tostring(self):
        """deprecated : use :meth:`to_string`"""
        return self.to_string()

    def set_name(self, name):
        if name is None:
            raise ValueError("FastxRecord must have a name and not None")
        self.name = name

    def set_comment(self, comment):
        self.comment = comment    
        
    def set_sequence(self, sequence, quality=None):
        """set sequence of this record.

        """
        self.sequence = sequence
        if quality is not None:
            if len(sequence) != len(quality):
                raise ValueError("sequence and quality length do not match: {} vs {}".format(
                    len(sequence), len(quality)))

            self.quality = quality
        else:
            self.quality = None

    def __str__(self):
        return self.to_string()

    cpdef array.array get_quality_array(self, int offset=33):
        '''return quality values as array after subtracting offset.'''
        if self.quality is None:
            return None
        return qualitystring_to_array(force_bytes(self.quality),
                                      offset=offset)


cdef class FastxFile:
    r"""Stream access to :term:`fasta` or :term:`fastq` formatted files.

    The file is automatically opened.

    Entries in the file can be both fastq or fasta formatted or even a
    mixture of the two.

    This file object permits iterating over all entries in the
    file. Random access is not implemented. The iteration returns
    objects of type :class:`FastqProxy`

    Parameters
    ----------

    filename : string
        Filename of fasta/fastq file to be opened.

    persist : bool

        If True (default) make a copy of the entry in the file during
        iteration. If set to False, no copy will be made. This will
        permit much faster iteration, but an entry will not persist
        when the iteration continues and an entry is read-only.

    Notes
    -----
    Prior to version 0.8.2, this class was called FastqFile.

    Raises
    ------

    IOError
        if file could not be opened


    Examples
    --------
    >>> with pysam.FastxFile(filename) as fh:
    ...    for entry in fh:
    ...        print(entry.name)
    ...        print(entry.sequence)
    ...        print(entry.comment)
    ...        print(entry.quality)
    >>> with pysam.FastxFile(filename) as fin, open(out_filename, mode='w') as fout:
    ...    for entry in fin:
    ...        fout.write(str(entry) + '\n')

    """
    def __cinit__(self, *args, **kwargs):
        # self.fastqfile = <gzFile*>NULL
        self._filename = None
        self.entry = NULL
        self._open(*args, **kwargs)

    def is_open(self):
        '''return true if samfile has been opened.'''
        return self.entry != NULL

    def _open(self, filename, persist=True):
        '''open a fastq/fasta file in *filename*

        Paramentes
        ----------

        persist : bool

            if True return a copy of the underlying data (default
            True).  The copy will persist even if the iteration
            on the file continues.

        '''
        if self.fastqfile != NULL:
            self.close()

        self._filename = encode_filename(filename)
        cdef char *cfilename = self._filename
        self.is_remote = hisremote(cfilename)

        # open file for reading
        if (self._filename != b"-"
            and not self.is_remote
            and not os.path.exists(filename)):
            raise IOError("file `%s` not found" % filename)

        self.persist = persist

        with nogil:
            self.fastqfile = bgzf_open(cfilename, "r")
            self.entry = kseq_init(self.fastqfile)
        self._filename = filename

    def close(self):
        '''close the file.'''
        if self.fastqfile != NULL:
            bgzf_close(self.fastqfile)
            self.fastqfile = NULL
        if self.entry != NULL:
            kseq_destroy(self.entry)
            self.entry = NULL

    def __dealloc__(self):
        if self.fastqfile != NULL:
            bgzf_close(self.fastqfile)
        if self.entry:
            kseq_destroy(self.entry)

    # context manager interface
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
        return False

    property closed:
        """bool indicating the current state of the file object.
        This is a read-only attribute; the close() method changes the value.
        """
        def __get__(self):
            return not self.is_open()

    property filename:
        """string with the filename associated with this object."""
        def __get__(self):
            return self._filename

    def __iter__(self):
        if not self.is_open():
            raise ValueError("I/O operation on closed file")
        return self

    cdef kseq_t * getCurrent(self):
        return self.entry

    cdef int cnext(self):
        '''C version of iterator
        '''
        with nogil:
            return kseq_read(self.entry)

    def __next__(self):
        """
        python version of next().
        """
        cdef int l
        with nogil:
            l = kseq_read(self.entry)
        if (l >= 0):
            if self.persist:
                return FastxRecord(proxy=makeFastqProxy(self.entry))
            return makeFastqProxy(self.entry)
        elif (l == -1):
            raise StopIteration
        elif (l == -2):
            raise ValueError('truncated quality string in {0}'
                             .format(self._filename))
        else:
            raise ValueError('unknown problem parsing {0}'
                             .format(self._filename))

# Compatibility Layer for pysam 0.8.1
cdef class FastqFile(FastxFile):
    """FastqFile is deprecated: use FastxFile instead"""
    pass

# Compatibility Layer for pysam < 0.8
cdef class Fastafile(FastaFile):
    """Fastafile is deprecated: use FastaFile instead"""
    pass

__all__ = ["FastaFile",
           "FastqFile",
           "FastxFile",
           "Fastafile",
           "FastxRecord",
           "FastqProxy"]
