#########################################################################
# Utility functions used across pysam
#########################################################################
cimport cython
from cpython cimport array as c_array

cpdef parse_region(reference=*, start=*, end=*, region=*)

#########################################################################
# Utility functions for quality string conversions

cpdef c_array.array qualitystring_to_array(input_str, int offset=*)
cpdef array_to_qualitystring(c_array.array arr, int offset=*)
cpdef qualities_to_qualitystring(qualities, int offset=*)

########################################################################
########################################################################
########################################################################
## Python 3 compatibility functions
########################################################################
cdef charptr_to_str(const char *s, encoding=*)
cdef bytes charptr_to_bytes(const char *s, encoding=*)
cdef charptr_to_str_w_len(const char* s, size_t n, encoding=*)
cdef force_str(object s, encoding=*)
cdef bytes force_bytes(object s, encoding=*)
cdef bytes encode_filename(object filename)
cdef from_string_and_size(const char *s, size_t length)

cdef extern from "pysam_util.h":

    int samtools_main(int argc, char *argv[])
    int bcftools_main(int argc, char *argv[])
    void pysam_set_stderr(int fd)
    void pysam_unset_stderr()
    void pysam_set_stdout(int fd)
    void pysam_set_stdout_fn(const char *)
    void pysam_unset_stdout()
    void set_optind(int)
