# cython: language_level=3
from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.string cimport memcpy, memcmp, strncpy, strlen, strdup
from libc.stdio cimport FILE, printf
cimport cython

from cpython cimport array
from pysam.libchtslib cimport faidx_t, kstring_t, BGZF

cdef extern from "htslib/kseq.h" nogil:
    """
    struct __kstream_t;
    #define kstream_t  struct __kstream_t
    __KSEQ_TYPE(type_t_unused_here)
    #undef kstream_t
    """
    ctypedef struct kseq_t:
        kstring_t name
        kstring_t comment
        kstring_t seq
        kstring_t qual


cdef class FastaFile:
    cdef bint is_remote
    cdef object _filename, _references, _lengths, reference2length
    cdef faidx_t* fastafile
    cdef char* _fetch(self, char* reference,
                      int start, int end, int* length) except? NULL


cdef class FastqProxy:
    cdef kseq_t * _delegate
    cdef cython.str to_string(self)
    cdef cython.str tostring(self)
    cpdef array.array get_quality_array(self, int offset=*)


cdef class FastxRecord:
    """
    Python container for pysam.libcfaidx.FastqProxy with persistence.
    """
    cdef public str comment, quality, sequence, name
    cdef cython.str to_string(self)
    cdef cython.str tostring(self)
    cpdef array.array get_quality_array(self, int offset=*)

cdef class FastxFile:
    cdef object _filename
    cdef BGZF * fastqfile
    cdef kseq_t * entry
    cdef bint persist
    cdef bint is_remote

    cdef kseq_t * getCurrent(self)
    cdef int cnext(self)


# Compatibility Layer for pysam 0.8.1
cdef class FastqFile(FastxFile):
    pass


# Compatibility Layer for pysam < 0.8
cdef class Fastafile(FastaFile):
    pass

