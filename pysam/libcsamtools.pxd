cdef extern from "samtools.pysam.h":

    int samtools_main(int argc, char *argv[])
    void samtools_set_stderr(int fd)
    void samtools_unset_stderr()
    void samtools_set_stdout(int fd)
    void samtools_set_stdout_fn(const char *)
    void samtools_unset_stdout()
    void samtools_set_optind(int)
