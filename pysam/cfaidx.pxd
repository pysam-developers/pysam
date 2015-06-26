from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.string cimport memcpy, memcmp, strncpy, strlen, strdup
from libc.stdio cimport FILE, printf
from cpython cimport array
cimport cython

from chtslib cimport faidx_t, kseq_t, gzFile

ctypedef pFastqProxy pFastqProxy_t
ctypedef pFastqFile pFastqFile_t

cdef extern from "htslib/kstring.h" nogil:
    ctypedef struct kstring_t:
        size_t l, m
        char *s

cdef class FastaFile:
    cdef object _filename, _references, _lengths, reference2length
    cdef faidx_t* fastafile
    cdef char* _fetch(self, char* reference,
                      int start, int end, int* length)


cdef class FastqProxy:
    cdef kseq_t * _delegate


cdef class FastxFile:
    cdef object _filename
    cdef gzFile fastqfile
    cdef kseq_t * entry

    cdef kseq_t * getCurrent(self)
    cdef int cnext(self)

# Compatibility Layer for pysam 0.8.1
cdef class FastqFile(FastxFile):
    pass

# Compatibility Layer for pysam < 0.8
cdef class Fastafile(FastaFile):
    pass

cdef class Fastqfile(FastxFile):
    pass

cdef class pFastqFile:
    cdef public FastxFile handle
    cpdef close(self)

cdef class pFastqProxy:
    """
    Python container for pysam.cfaidx.FastqProxy with persistence.
    """
    cdef public bytes comment, quality, sequence, name
    cdef cython.str tostring(self)
    cpdef array.array getQualArray(self)

cdef array.array cs_to_ph(bytes input_str)
