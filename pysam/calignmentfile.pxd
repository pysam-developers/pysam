from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libc.stdlib cimport malloc, calloc, realloc, free
from libc.string cimport memcpy, memcmp, strncpy, strlen, strdup
from libc.stdio cimport FILE, printf

from cfaidx cimport faidx_t, Fastafile
from chtslib cimport *

cdef extern from *:
    ctypedef char* const_char_ptr "const char*"

cdef extern from "htslib_util.h":

    int hts_set_verbosity(int verbosity)
    int hts_get_verbosity()

    # add *nbytes* into the variable length data of *src* at *pos*
    bam1_t * pysam_bam_update(bam1_t * b,
                              size_t nbytes_old,
                              size_t nbytes_new,
                              uint8_t * pos)

    # now: static
    int aux_type2size(int)

    char * pysam_bam_get_qname(bam1_t * b)
    uint32_t * pysam_bam_get_cigar(bam1_t * b)
    uint8_t * pysam_bam_get_seq(bam1_t * b)
    uint8_t * pysam_bam_get_qual(bam1_t * b)
    uint8_t * pysam_bam_get_aux(bam1_t * b)
    int pysam_bam_get_l_aux(bam1_t * b)
    char pysam_bam_seqi(uint8_t * s, int i)

    uint16_t pysam_get_bin(bam1_t * b)
    uint8_t pysam_get_qual(bam1_t * b)
    uint8_t pysam_get_l_qname(bam1_t * b)
    uint16_t pysam_get_flag(bam1_t * b)
    uint16_t pysam_get_n_cigar(bam1_t * b)
    void pysam_set_bin(bam1_t * b, uint16_t v)
    void pysam_set_qual(bam1_t * b, uint8_t v)
    void pysam_set_l_qname(bam1_t * b, uint8_t v)
    void pysam_set_flag(bam1_t * b, uint16_t v)
    void pysam_set_n_cigar(bam1_t * b, uint16_t v)
    void pysam_update_flag(bam1_t * b, uint16_t v, uint16_t flag)


cdef extern from "samfile_util.h":

    int bam_cap_mapQ(bam1_t *b, char *ref, int thres)
    int bam_prob_realn(bam1_t *b, const char *ref)

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

# Exposing pysam extension classes
#
# Note: need to declare all C fields and methods here
cdef class AlignedSegment:

    # object that this AlignedSegment represents
    cdef bam1_t * _delegate

    # add an alignment tag with value to the AlignedSegment
    # an existing tag of the same name will be replaced.
    cpdef set_tag(self, tag, value, value_type=?, replace=?)

    # add an alignment tag with value to the AlignedSegment
    # an existing tag of the same name will be replaced.
    cpdef get_tag(self, tag)

    # return true if tag exists
    cpdef has_tag(self, tag)

cdef class AlignmentFile:

    cdef object _filename

    # pointer to htsFile structure
    cdef htsFile * htsfile

    # pointer to index
    cdef hts_idx_t *index
    # header structure
    cdef bam_hdr_t * header
    # true if file is bam format
    cdef readonly bint is_bam
    # true if file is bam format
    cdef readonly bint is_cram
    # true if not a file but a stream
    cdef readonly bint is_stream
    # true if file is not on the local filesystem
    cdef readonly bint is_remote
    # current read within iteration
    cdef bam1_t * b
    # file opening mode
    cdef char * mode

    # beginning of read section
    cdef int64_t start_offset

    cdef bam_hdr_t * _buildHeader(self, new_header)
    cdef bam1_t * getCurrent(self)
    cdef int cnext(self)

    # write an aligned read
    cpdef int write(self, AlignedSegment read) except -1

    cdef char * _getrname(self, int tid)

cdef class PileupColumn:
    cdef bam_pileup1_t ** plp
    cdef int tid
    cdef int pos
    cdef int n_pu

cdef class PileupRead:
    cdef AlignedSegment _alignment
    cdef int32_t  _qpos
    cdef int _indel
    cdef int _level
    cdef uint32_t _is_del
    cdef uint32_t _is_head
    cdef uint32_t _is_tail
    cdef uint32_t _is_refskip

cdef class IteratorRow:
    cdef int retval
    cdef bam1_t * b
    cdef AlignmentFile samfile
    cdef htsFile * htsfile
    cdef bam_hdr_t * header
    cdef int owns_samfile

cdef class IteratorRowRegion(IteratorRow):
    cdef hts_itr_t * iter
    cdef bam1_t * getCurrent( self )
    cdef int cnext(self)

cdef class IteratorRowHead(IteratorRow):
    cdef int max_rows
    cdef int current_row
    cdef bam1_t * getCurrent(self)
    cdef int cnext(self)

cdef class IteratorRowAll(IteratorRow):
    cdef bam1_t * getCurrent( self )
    cdef int cnext(self)

cdef class IteratorRowAllRefs(IteratorRow):
    cdef int         tid
    cdef IteratorRowRegion rowiter

cdef class IteratorRowSelection(IteratorRow):
    cdef int current_pos
    cdef positions
    cdef bam1_t * getCurrent( self )
    cdef int cnext(self)

cdef class IteratorColumn:

    # result of the last plbuf_push
    cdef IteratorRowRegion iter
    cdef int tid
    cdef int pos
    cdef int n_plp
    cdef int mask
    cdef bam_pileup1_t * plp
    cdef bam_plp_t pileup_iter
    cdef __iterdata iterdata
    cdef AlignmentFile samfile
    cdef Fastafile fastafile
    cdef stepper
    cdef int max_depth

    cdef int cnext(self)
    cdef char * getSequence( self )
    cdef setMask(self, mask)
    cdef setupIteratorData(self,
                           int tid,
                           int start,
                           int end,
                           int multiple_iterators = ?)

    cdef reset(self, tid, start, end)
    cdef _free_pileup_iter(self)

cdef class IteratorColumnRegion(IteratorColumn):
    cdef int start
    cdef int end
    cdef int truncate

cdef class IteratorColumnAllRefs(IteratorColumn):
    pass

cdef class IndexedReads:
    cdef AlignmentFile samfile
    cdef htsFile * htsfile
    cdef index
    cdef int owns_samfile
    cdef bam_hdr_t * header
