cdef class TestClass:
     cdef int x

from cpython cimport array

cdef array.array _chars_to_array(bytes input_str, int offset=*)

########################################################################
########################################################################
########################################################################
## Python 3 compatibility functions
########################################################################
cdef _charptr_to_str(char *s, encoding=*)
cdef _force_str(object s, encoding=*)
cdef bytes _force_bytes(object s, encoding=*)
cdef bytes _force_cmdline_bytes(object s, encoding=*)
cdef bytes _encode_filename(object filename)
cdef from_string_and_size(char *s, size_t length)