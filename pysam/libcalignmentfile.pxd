from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.string cimport memcpy, memcmp, strncpy, strlen, strdup
from libc.stdio cimport FILE, printf

from pysam.libcfaidx cimport faidx_t, FastaFile
from pysam.libcalignedsegment cimport AlignedSegment
from pysam.libchtslib cimport *

from cpython cimport array
cimport cython

cdef extern from *:
    ctypedef char* const_char_ptr "const char*"

cdef extern from "htslib_util.h":

    char * pysam_bam_get_qname(bam1_t * b)

####################################################################
# Utility types

ctypedef struct __iterdata:
    htsFile * htsfile
    bam_hdr_t * header
    hts_itr_t * iter
    faidx_t * fastafile
    int tid
    char * seq
    int seq_len
    int min_mapping_quality
    int flag_require
    int flag_filter
    bint compute_baq
    bint redo_baq
    bint ignore_orphans
    int adjust_capq_threshold


cdef class AlignmentHeader(object):
    cdef bam_hdr_t *ptr

cdef class AlignmentFile(HTSFile):
    cdef readonly object reference_filename
    cdef readonly AlignmentHeader header

    # pointer to index
    cdef hts_idx_t *index

    # current read within iteration
    cdef bam1_t * b

    cdef bam1_t * getCurrent(self)
    cdef int cnext(self)

    # write an aligned read
    cpdef int write(self, AlignedSegment read) except -1


cdef class IteratorRow:
    cdef int retval
    cdef bam1_t * b
    cdef AlignmentFile samfile
    cdef htsFile * htsfile
    cdef hts_idx_t * index
    cdef AlignmentHeader header
    cdef int owns_samfile


cdef class IteratorRowRegion(IteratorRow):
    cdef hts_itr_t * iter
    cdef bam1_t * getCurrent(self)
    cdef int cnext(self)


cdef class IteratorRowHead(IteratorRow):
    cdef int max_rows
    cdef int current_row
    cdef bam1_t * getCurrent(self)
    cdef int cnext(self)


cdef class IteratorRowAll(IteratorRow):
    cdef bam1_t * getCurrent(self)
    cdef int cnext(self)


cdef class IteratorRowAllRefs(IteratorRow):
    cdef int         tid
    cdef IteratorRowRegion rowiter


cdef class IteratorRowSelection(IteratorRow):
    cdef int current_pos
    cdef positions
    cdef bam1_t * getCurrent(self)
    cdef int cnext(self)


cdef class IteratorColumn:

    # result of the last plbuf_push
    cdef IteratorRowRegion iter
    cdef int tid
    cdef int pos
    cdef int n_plp
    cdef uint32_t min_base_quality
    cdef const bam_pileup1_t * plp
    cdef bam_mplp_t pileup_iter
    cdef __iterdata iterdata
    cdef AlignmentFile samfile
    cdef FastaFile fastafile
    cdef stepper
    cdef int max_depth
    cdef bint ignore_overlaps

    cdef int cnext(self)
    cdef char * get_sequence(self)
    cdef _setup_iterator(self,
                         int tid,
                         int start,
                         int stop,
                         int multiple_iterators=?)
    cdef _setup_raw_rest_iterator(self)

    cdef reset(self, tid, start, stop)
    cdef _free_pileup_iter(self)
    # backwards compatibility
    cdef char * getSequence(self)


cdef class IteratorColumnRegion(IteratorColumn):
    cdef int start
    cdef int stop
    cdef int truncate


cdef class IteratorColumnAllRefs(IteratorColumn):
    pass


cdef class IteratorColumnAll(IteratorColumn):
    pass


cdef class IndexedReads:
    cdef AlignmentFile samfile
    cdef htsFile * htsfile
    cdef object index
    cdef int owns_samfile
    cdef AlignmentHeader header
