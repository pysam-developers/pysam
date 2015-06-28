cdef class TestClass:
     cdef int x

from cpython cimport array

cdef array.array _chars_to_array(bytes input_str, int offset=*)

########################################################################
########################################################################
########################################################################
## Python 3 compatibility functions
########################################################################
cdef charptr_to_str(char *s, encoding=*)
cdef force_str(object s, encoding=*)
cdef bytes force_bytes(object s, encoding=*)
cdef bytes force_cmdline_bytes(object s, encoding=*)
cdef bytes encode_filename(object filename)
cdef from_string_and_size(char *s, size_t length)