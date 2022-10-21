# cython: embedsignature=True
# cython: profile=True
# adds doc-strings for sphinx
import tempfile
import os
import sys
import types
import itertools
import struct
import ctypes
import collections
import re
import platform
import warnings
from cpython cimport PyErr_SetString, \
    PyBytes_Check, \
    PyUnicode_Check, \
    PyBytes_FromStringAndSize

from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment


cdef class Samfile(AlignmentFile):
    '''Deprecated alternative for :class:`~pysam.AlignmentFile`

    Added for backwards compatibility with pysam <= 0.8.0
    '''
    pass


cdef class AlignedRead(AlignedSegment):
    '''Deprecated alternative for :class:`~pysam.AlignedSegment`

    Added for backwards compatibility with pysam <= 0.8.0
    '''
    pass


__all__ = ['Samfile', 'AlignedRead']


