from libc.stdlib cimport calloc, free


cdef extern from "pysam_util.h":
    int pysam_dispatch(int argc, char *argv[])
    void pysam_set_stderr(int fd)
    void pysam_unset_stderr()
