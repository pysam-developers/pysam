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

from cpython.version cimport PY_MAJOR_VERSION

########################################################################
########################################################################
########################################################################
## Python 3 compatibility functions
########################################################################
IS_PYTHON3 = PY_MAJOR_VERSION >= 3

cdef from_string_and_size(char* s, size_t length):
    if PY_MAJOR_VERSION < 3:
        return s[:length]
    else:
        return s[:length].decode("ascii")

# filename encoding (copied from lxml.etree.pyx)
cdef str _FILENAME_ENCODING
_FILENAME_ENCODING = sys.getfilesystemencoding()
if _FILENAME_ENCODING is None:
    _FILENAME_ENCODING = sys.getdefaultencoding()
if _FILENAME_ENCODING is None:
    _FILENAME_ENCODING = 'ascii'

#cdef char* _C_FILENAME_ENCODING
#_C_FILENAME_ENCODING = <char*>_FILENAME_ENCODING

cdef bytes _encodeFilename(object filename):
    """Make sure a filename is 8-bit encoded (or None)."""
    if filename is None:
        return None
    elif PyBytes_Check(filename):
        return filename
    elif PyUnicode_Check(filename):
        return filename.encode(_FILENAME_ENCODING)
    else:
        raise TypeError, u"Argument must be string or unicode."

cdef _forceStr(object s):
    """Return s converted to str type of current Python
    (bytes in Py2, unicode in Py3)"""
    if s is None:
        return None
    if PY_MAJOR_VERSION < 3:
        return s
    elif PyBytes_Check(s):
        return s.decode('ascii')
    else:
        # assume unicode
        return s

cdef bytes _forceBytes(object s):
    u"""convert string or unicode object to bytes, assuming ascii encoding.
    """
    if PY_MAJOR_VERSION < 3:
        return s
    elif s is None:
        return None
    elif PyBytes_Check(s):
        return s
    elif PyUnicode_Check(s):
        return s.encode('ascii')
    else:
        raise TypeError, u"Argument must be string, bytes or unicode."

cdef inline bytes _forceCmdlineBytes(object s):
    return _forceBytes(s)

cdef _charptr_to_str(char* s):
    if PY_MAJOR_VERSION < 3:
        return s
    else:
        return s.decode("ascii")

__all__ = []

