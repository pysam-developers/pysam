cdef extern from "bcftools.pysam.h":

    int bcftools_main(int argc, char *argv[])
    void bcftools_set_stderr(int fd)
    void bcftools_unset_stderr()
    void bcftools_set_stdout(int fd)
    void bcftools_set_stdout_fn(const char *)
    void bcftools_unset_stdout()
    void bcftools_set_optind(int)
