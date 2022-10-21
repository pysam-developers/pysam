# cython: language_level=3
from pysam.libchtslib cimport *

cdef extern from "htslib_util.h":

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

    uint8_t pysam_get_qual(bam1_t * b)
    uint32_t pysam_get_n_cigar(bam1_t * b)
    void pysam_set_qual(bam1_t * b, uint8_t v)
    void pysam_set_n_cigar(bam1_t * b, uint32_t v)
    void pysam_update_flag(bam1_t * b, uint16_t v, uint16_t flag)


from pysam.libcalignmentfile cimport AlignmentFile, AlignmentHeader
ctypedef AlignmentFile AlignmentFile_t


# Note: need to declare all C fields and methods here
cdef class AlignedSegment:

    # object that this AlignedSegment represents
    cdef bam1_t * _delegate

    # the header that a read is associated with
    cdef readonly AlignmentHeader header

    # caching of array properties for quick access
    cdef object cache_query_qualities
    cdef object cache_query_alignment_qualities
    cdef object cache_query_sequence
    cdef object cache_query_alignment_sequence

    # add an alignment tag with value to the AlignedSegment
    # an existing tag of the same name will be replaced.
    cpdef set_tag(self, tag, value, value_type=?, replace=?)

    # get an alignment tag from the AlignedSegment
    cpdef get_tag(self, tag, with_value_type=?)

    # return true if tag exists
    cpdef has_tag(self, tag)

    # returns a valid sam alignment string
    cpdef to_string(self)

    # returns a valid sam alignment string (deprecated)
    cpdef tostring(self, htsfile=*)


cdef class PileupColumn:
    cdef const bam_pileup1_t ** plp
    cdef int tid
    cdef int pos
    cdef int n_pu
    cdef AlignmentHeader header
    cdef uint32_t min_base_quality
    cdef kstring_t buf
    cdef char * reference_sequence

cdef class PileupRead:
    cdef int32_t  _qpos
    cdef AlignedSegment _alignment
    cdef int _indel
    cdef int _level
    cdef uint32_t _is_del
    cdef uint32_t _is_head
    cdef uint32_t _is_tail
    cdef uint32_t _is_refskip

# factory methods
cdef AlignedSegment makeAlignedSegment(
    bam1_t * src,
    AlignmentHeader header)

cdef PileupColumn makePileupColumn(
    const bam_pileup1_t ** plp,
    int tid,
    int pos,
    int n_pu,
    uint32_t min_base_quality,
    char * reference_sequence,
    AlignmentHeader header)

cdef PileupRead makePileupRead(const bam_pileup1_t * src,
		               AlignmentHeader header)

cdef uint32_t get_alignment_length(bam1_t * src)
