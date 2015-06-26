import types
import sys
import string

from cpython.version cimport PY_MAJOR_VERSION
from cpython cimport PyErr_SetString, PyBytes_Check
from cpython cimport PyUnicode_Check, PyBytes_FromStringAndSize

from libc.stdio cimport printf

from cpython cimport array

cdef class TestClass:
    def __init__(self):
        pass

cdef array.array _chars_to_array(bytes input_str, int offset=0):
    '''convert a buffer of characters to a byte array'''
    cdef char i
    return array.array('B', [i - offset for i in input_str])
