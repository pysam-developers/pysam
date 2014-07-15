from libc.stdlib cimport calloc, free

cdef extern from "pysam_util.h":
    int pysam_dispatch(int argc, char *argv[])
    void pysam_set_stderr(int fd)
    void pysam_unset_stderr()


cdef extern from "sam.h":

    ctypedef struct bam1_t

    # functions not actually declared in sam.h, but available
    # as extern
    # 
    # implemented in samtools/bam_md.c
    int bam_prob_realn(bam1_t *b, char *ref)
    int bam_cap_mapQ(bam1_t *b, char *ref, int thres)

