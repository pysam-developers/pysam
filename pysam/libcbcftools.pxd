# cython: language_level=3
cdef extern from "bcftools.pysam.h":

    int bcftools_dispatch(int argc, char *argv[])
    void bcftools_set_stderr(int fd)
    void bcftools_close_stderr()
    void bcftools_set_stdout(int fd)
    void bcftools_set_stdout_fn(const char *)
    void bcftools_close_stdout()
