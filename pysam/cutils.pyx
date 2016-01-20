import types
import sys
import string
import re

from cpython.version cimport PY_MAJOR_VERSION
from cpython cimport PyBytes_Check, PyUnicode_Check

from cpython cimport array as c_array
cimport cython

#####################################################################
# hard-coded constants
cdef int MAX_POS = 2 << 29

#################################################################
# Utility functions for quality string conversions
cpdef c_array.array qualitystring_to_array(input_str, int offset=33):
    """convert a qualitystring to an array of quality values."""
    if input_str is None:
        return None
    qs = force_bytes(input_str)
    cdef char i
    return c_array.array('B', [i - offset for i in qs])


cpdef array_to_qualitystring(c_array.array qualities, int offset=33):
    """convert an array of quality values to a string."""
    if qualities is None:
        return None
    cdef int x
    
    cdef c_array.array result
    result = c_array.clone(qualities, len(qualities), zero=False)
    
    for x from 0 <= x < len(qualities):
        result[x] = qualities[x] + offset
    return force_str(result.tostring())


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
        return force_str("".join([chr(x + offset) for x in qualities]))


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
        raise TypeError(u"Argument must be string or unicode.")

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
        raise TypeError(u"Argument must be string, bytes or unicode.")

cdef bytes force_cmdline_bytes(object s, encoding="ascii"):
    return force_bytes(s)

cdef charptr_to_str(char* s, encoding="ascii"):
    if s == NULL:
        return None
    if PY_MAJOR_VERSION < 3:
        return s
    else:
        return s.decode(encoding)

cdef bytes charptr_to_bytes(char* s, encoding="ascii"):
    if s == NULL:
        return None
    else:
        return s

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

cpdef parse_region(reference=None,
                   start=None,
                   end=None,
                   region=None):
    """parse alternative ways to specify a genomic region. A region can
    either be specified by :term:`reference`, `start` and
    `end`. `start` and `end` denote 0-based, half-open
    intervals.

    Alternatively, a samtools :term:`region` string can be
    supplied.

    If any of the coordinates are missing they will be replaced by the
    minimum (`start`) or maximum (`end`) coordinate.

    Note that region strings are 1-based, while `start` and `end` denote
    an interval in python coordinates.

    Returns
    -------

    tuple :  a tuple of `reference`, `start` and `end`.

    Raises
    ------

    ValueError
       for invalid or out of bounds regions.

    """
    cdef int rtid
    cdef long long rstart
    cdef long long rend

    rtid = -1
    rstart = 0
    rend = MAX_POS
    if start != None:
        try:
            rstart = start
        except OverflowError:
            raise ValueError('start out of range (%i)' % start)

    if end != None:
        try:
            rend = end
        except OverflowError:
            raise ValueError('end out of range (%i)' % end)

    if region:
        region = force_str(region)
        parts = re.split("[:-]", region)
        reference = parts[0]
        if len(parts) >= 2:
            rstart = int(parts[1]) - 1
        if len(parts) >= 3:
            rend = int(parts[2])

    if not reference:
        return None, 0, 0

    if not 0 <= rstart < MAX_POS:
        raise ValueError('start out of range (%i)' % rstart)
    if not 0 <= rend <= MAX_POS:
        raise ValueError('end out of range (%i)' % rend)
    if rstart > rend:
        raise ValueError(
            'invalid region: start (%i) > end (%i)' % (rstart, rend))

    return force_bytes(reference), rstart, rend


__all__ = ["qualitystring_to_array",
           "array_to_qualitystring",
           "qualities_to_qualitystring"]
