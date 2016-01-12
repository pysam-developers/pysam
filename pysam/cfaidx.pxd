from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.string cimport memcpy, memcmp, strncpy, strlen, strdup
from libc.stdio cimport FILE, printf
cimport cython

from cpython cimport array
from pysam.chtslib cimport faidx_t, gzFile, kstring_t

# These functions are put here and not in chtslib.pxd in order
# to avoid warnings for unused functions.
cdef extern from "pysam_stream.h" nogil:

    ctypedef struct kstream_t:
        pass

    ctypedef struct kseq_t:
        kstring_t name
        kstring_t comment
        kstring_t seq
        kstring_t qual

    gzFile gzopen(char *, char *)
    kseq_t *kseq_init(gzFile)
    int kseq_read(kseq_t *)
    void kseq_destroy(kseq_t *)
    int gzclose(gzFile)

    kstream_t *ks_init(gzFile)
    void ks_destroy(kstream_t *)

    # Retrieve characters from stream until delimiter
    # is reached placing results in str.
    int ks_getuntil(kstream_t *,
                    int delimiter,
                    kstring_t * str,
                    int * dret)

cdef class FastaFile:
    cdef object _filename, _references, _lengths, reference2length
    cdef faidx_t* fastafile
    cdef char* _fetch(self, char* reference,
                      int start, int end, int* length)


cdef class FastqProxy:
    cdef kseq_t * _delegate
    cdef cython.str tostring(self)
    cpdef array.array get_quality_array(self, int offset=*)


cdef class PersistentFastqProxy:
    """
    Python container for pysam.cfaidx.FastqProxy with persistence.
    """
    cdef public str comment, quality, sequence, name
    cdef cython.str tostring(self)
    cpdef array.array get_quality_array(self, int offset=*)


cdef class FastxFile:
    cdef object _filename
    cdef gzFile fastqfile
    cdef kseq_t * entry
    cdef bint persist

    cdef kseq_t * getCurrent(self)
    cdef int cnext(self)


# Compatibility Layer for pysam 0.8.1
cdef class FastqFile(FastxFile):
    pass


# Compatibility Layer for pysam < 0.8
cdef class Fastafile(FastaFile):
    pass

