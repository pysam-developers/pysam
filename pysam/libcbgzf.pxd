# cython: language_level=3

from pysam.libchtslib cimport BGZF

cdef class BGZFile(object):
    cdef BGZF *bgzf
    cdef readonly object name, index
