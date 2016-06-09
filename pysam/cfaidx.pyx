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
# class PersistentFastqProxy
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
from cpython cimport array

from cpython cimport PyErr_SetString, \
    PyBytes_Check, \
    PyUnicode_Check, \
    PyBytes_FromStringAndSize

from cpython.version cimport PY_MAJOR_VERSION

from pysam.chtslib cimport \
    faidx_nseq, fai_load, fai_destroy, fai_fetch, \
    faidx_seq_len, \
    faidx_fetch_seq, hisremote, \
    bgzf_open, bgzf_close

from pysam.cutils cimport force_bytes, force_str, charptr_to_str
from pysam.cutils cimport encode_filename, from_string_and_size
from pysam.cutils cimport qualitystring_to_array, parse_region

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

    def _open(self, filename, filepath_index=None):
        '''open an indexed fasta file.

        This method expects an indexed fasta file.
        '''

        # close a previously opened file
        if self.fastafile != NULL:
            self.close()

        self._filename = encode_filename(filename)
        cdef char *cfilename = self._filename
        self.is_remote = hisremote(cfilename)

        if filepath_index is not None:
            raise NotImplementedError(
                "setting an explicit path for the index "
                "is not implemented")

        # open file for reading
        if (self._filename != b"-"
            and not self.is_remote
            and not os.path.exists(filename)):
            raise IOError("file `%s` not found" % filename)

        with nogil:
            self.fastafile = fai_load(cfilename)

        if self.fastafile == NULL:
            raise IOError("could not open file `%s`" % filename)

        if self.is_remote:
            filepath_index = os.path.basename(
                re.sub("[^:]+:[/]*", "", filename)) + ".fai"
        elif filepath_index is None:
            filepath_index = filename + ".fai"

        if not os.path.exists(filepath_index):
            raise ValueError("could not locate index file {}".format(
                filepath_index))

        with open(filepath_index) as inf:
            data = [x.split("\t") for x in inf]
            self._references = tuple(x[0] for x in data)
            self._lengths = tuple(int(x[1]) for x in data)
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
        """"bool indicating the current state of the file object.
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
        """"int with the number of :term:`reference` sequences in the file.
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

        reference, rstart, rend = parse_region(reference, start, end, region)

        if reference is None:
            raise ValueError("no sequence/region supplied.")

        if rstart == rend:
            return ""

        ref = reference
        with nogil:
            length = faidx_seq_len(self.fastafile, ref)
        if length == -1:
            raise KeyError("sequence '%s' not present" % reference)
        if rstart >= length:
            return ""

        # fai_fetch adds a '\0' at the end
        with nogil:
            seq = faidx_fetch_seq(self.fastafile,
                                  ref,
                                  rstart,
                                  rend-1,
                                  &length)

        if seq == NULL:
            raise ValueError(
                "failure when retrieving sequence on '%s'" % reference)

        try:
            return charptr_to_str(seq)
        finally:
            free(seq)

    cdef char * _fetch(self, char * reference, int start, int end, int * length):
        '''fetch sequence for reference, start and end'''

        with nogil:
            return faidx_fetch_seq(self.fastafile,
                                   reference,
                                   start,
                                   end-1,
                                   length)

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
    def __init__(self): pass

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

    cdef cython.str tostring(self):
        if self.comment is None:
            comment = ""
        else:
            comment = " %s" % self.comment

        if self.quality is None:
            return ">%s%s\n%s" % (self.name, comment, self.sequence)
        else:
            return "@%s%s\n%s\n+\n%s" % (self.name, comment,
                                         self.sequence, self.quality)

    def __str__(self):
        return self.tostring()

    cpdef array.array get_quality_array(self, int offset=33):
        '''return quality values as integer array after subtracting offset.'''
        if self.quality is None:
            return None
        return qualitystring_to_array(force_bytes(self.quality),
                                      offset=offset)

cdef class PersistentFastqProxy:
    """
    Python container for pysam.cfaidx.FastqProxy with persistence.
    Needed to compare multiple fastq records from the same file.
    """
    def __init__(self, FastqProxy FastqRead):
        self.comment = FastqRead.comment
        self.quality = FastqRead.quality
        self.sequence = FastqRead.sequence
        self.name = FastqRead.name

    cdef cython.str tostring(self):
        if self.comment is None:
            comment = ""
        else:
            comment = " %s" % self.comment

        if self.quality is None:
            return ">%s%s\n%s" % (self.name, comment, self.sequence)
        else:
            return "@%s%s\n%s\n+\n%s" % (self.name, comment,
                                         self.sequence, self.quality)

    def __str__(self):
        return self.tostring()

    cpdef array.array get_quality_array(self, int offset=33):
        '''return quality values as array after subtracting offset.'''
        if self.quality is None:
            return None
        return qualitystring_to_array(force_bytes(self.quality),
                                      offset=offset)


cdef class FastxFile:
    """Stream access to :term:`fasta` or :term:`fastq` formatted files.

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
        permit faster iteration, but an entry will not persist when
        the iteration continues.

    Notes
    -----
    Prior to version 0.8.2, this was called FastqFile.

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
        """"bool indicating the current state of the file object.
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
                return PersistentFastqProxy(makeFastqProxy(self.entry))
            return makeFastqProxy(self.entry)
        else:
            raise StopIteration

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
           "Fastafile"]
