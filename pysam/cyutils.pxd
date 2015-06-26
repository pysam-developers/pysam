cdef class TestClass:
     cdef int x

from cpython cimport array

cdef array.array _chars_to_array(bytes input_str)