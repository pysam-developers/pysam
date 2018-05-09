from pysam.libchtslib cimport bam1_t

cdef extern from "bam_md.h" nogil:
    void bam_fillmd1_core(bam1_t *b, char *ref, int ref_len, int flag, int max_nm)
