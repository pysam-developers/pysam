from pysam.chtslib cimport *

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


from pysam.calignmentfile cimport AlignmentFile
ctypedef AlignmentFile AlignmentFile_t


# Note: need to declare all C fields and methods here
cdef class AlignedSegment:

    # object that this AlignedSegment represents
    cdef bam1_t * _delegate

    # the file from which this AlignedSegment originates (can be None)
    cdef AlignmentFile _alignment_file

    # caching of array properties for quick access
    cdef object cache_query_qualities
    cdef object cache_query_alignment_qualities
    cdef object cache_query_sequence
    cdef object cache_query_alignment_sequence

    # add an alignment tag with value to the AlignedSegment
    # an existing tag of the same name will be replaced.
    cpdef set_tag(self, tag, value, value_type=?, replace=?)

    # add an alignment tag with value to the AlignedSegment
    # an existing tag of the same name will be replaced.
    cpdef get_tag(self, tag, with_value_type=?)

    # return true if tag exists
    cpdef has_tag(self, tag)

    # returns a valid sam alignment string
    cpdef tostring(self, AlignmentFile_t handle)


cdef class PileupColumn:
    cdef bam_pileup1_t ** plp
    cdef int tid
    cdef int pos
    cdef int n_pu
    cdef AlignmentFile _alignment_file


cdef class PileupRead:
    cdef AlignedSegment _alignment
    cdef int32_t  _qpos
    cdef int _indel
    cdef int _level
    cdef uint32_t _is_del
    cdef uint32_t _is_head
    cdef uint32_t _is_tail
    cdef uint32_t _is_refskip

# factor methods
cdef makeAlignedSegment(bam1_t * src, AlignmentFile alignment_file)
cdef makePileupColumn(bam_pileup1_t ** plp, int tid, int pos, int n_pu, AlignmentFile alignment_file)
cdef inline makePileupRead(bam_pileup1_t * src, AlignmentFile alignment_file)
cdef inline uint32_t get_alignment_length(bam1_t * src)
