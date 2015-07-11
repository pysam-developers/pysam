import types
import sys
import string

from cpython.version cimport PY_MAJOR_VERSION
from cpython cimport PyBytes_Check, PyUnicode_Check

from cpython cimport array as c_array
cimport cython


#################################################################
# Utility functions for quality string conversions

PHRED_OFFSET_STRING64 = string.maketrans(
        "".join(chr(x) for x in xrange(
            0, 63)),
        "".join(chr(x + 64) for x in xrange(
            0, 63)))

PHRED_OFFSET_STRING33 = string.maketrans(
        "".join(chr(x) for x in xrange(
            0, 94)),
        "".join(chr(x + 33) for x in xrange(
            0, 94)))

cpdef c_array.array qualitystring_to_array(bytes input_str, int offset=33):
    """convert a qualitystring to an array of quality values."""
    if input_str == None:
        return None
    cdef char i
    return c_array.array('B', [i - offset for i in input_str])


cpdef array_to_qualitystring(c_array.array qualities, int offset=33):
    """convert an array of quality values to a string."""
    if qualities is None:
        return None
    cdef char i
    if offset == 33:
        return qualities.tostring().translate(PHRED_OFFSET_STRING33)
    elif offset == 64:
        return qualities.tostring().translate(PHRED_OFFSET_STRING64)
    else:
        return "".join([i + offset for i in qualities])


cpdef qualities_to_qualitystring(qualities, int offset=33):
    """convert a list or array of quality scores to the string
    representation used in the SAM format.

    Parameters
    ----------
    offset : int
        offset to be added to the quality scores to arrive at
        the characters of the quality string (default=33).

    Returns
    -------
    string
         a quality string

    """
    cdef char x
    if qualities is None:
        return None
    elif isinstance(qualities, c_array.array):
        return array_to_qualitystring(qualities, offset=offset)
    else:
        # tuples and lists
        return "".join([chr(x + offset) for x in qualities])


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

cdef bytes encode_filename(object filename):
    """Make sure a filename is 8-bit encoded (or None)."""
    if filename is None:
        return None
    elif PyBytes_Check(filename):
        return filename
    elif PyUnicode_Check(filename):
        return filename.encode(_FILENAME_ENCODING)
    else:
        raise TypeError, u"Argument must be string or unicode."

cdef bytes force_bytes(object s, encoding="ascii"):
    u"""convert string or unicode object to bytes, assuming
    ascii encoding.
    """
    if PY_MAJOR_VERSION < 3:
        return s
    elif s is None:
        return None
    elif PyBytes_Check(s):
        return s
    elif PyUnicode_Check(s):
        return s.encode(encoding)
    else:
        raise TypeError, u"Argument must be string, bytes or unicode."

cdef bytes force_cmdline_bytes(object s, encoding="ascii"):
    return force_bytes(s)

cdef charptr_to_str(char* s, encoding="ascii"):
    if PY_MAJOR_VERSION < 3:
        return s
    else:
        return s.decode(encoding)

cdef force_str(object s, encoding="ascii"):
    """Return s converted to str type of current Python
    (bytes in Py2, unicode in Py3)"""
    if s is None:
        return None
    if PY_MAJOR_VERSION < 3:
        return s
    elif PyBytes_Check(s):
        return s.decode(encoding)
    else:
        # assume unicode
        return s

__all__ = ["qualitystring_to_array",
           "array_to_qualitystring",
           "qualities_to_qualitystring"]
