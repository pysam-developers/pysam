# cython: language_level=3
from pysam.libcalignmentfile cimport AlignedSegment, AlignmentFile

#################################################
# Compatibility Layer for pysam < 0.8

# import all declarations from htslib
from pysam.libchtslib cimport *

cdef class AlignedRead(AlignedSegment):
    pass

cdef class Samfile(AlignmentFile):
    pass

# import the conversion functions
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
    uint32_t pysam_get_n_cigar(bam1_t * b)
    void pysam_set_bin(bam1_t * b, uint16_t v)
    void pysam_set_qual(bam1_t * b, uint8_t v)
    void pysam_set_l_qname(bam1_t * b, uint8_t v)
    void pysam_set_flag(bam1_t * b, uint16_t v)
    void pysam_set_n_cigar(bam1_t * b, uint32_t v)
    void pysam_update_flag(bam1_t * b, uint16_t v, uint16_t flag)
