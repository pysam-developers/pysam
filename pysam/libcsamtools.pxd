# cython: language_level=3
cdef extern from "samtools.pysam.h":

    int samtools_dispatch(int argc, char *argv[])
    void samtools_set_stderr(int fd)
    void samtools_close_stderr()
    void samtools_set_stdout(int fd)
    void samtools_set_stdout_fn(const char *)
    void samtools_close_stdout()
