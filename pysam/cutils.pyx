import types
import sys
import string
import re
import tempfile
import os
import io
from contextlib import contextmanager

from cpython.version cimport PY_MAJOR_VERSION
from cpython cimport PyBytes_Check, PyUnicode_Check
from cpython cimport array as c_array
from libc.stdlib cimport calloc, free
from libc.string cimport strncpy
from libc.stdio cimport fprintf, stderr, fflush
from libc.stdio cimport stdout as c_stdout
from posix.fcntl cimport open as c_open, O_WRONLY

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
cdef bint IS_PYTHON3 = PY_MAJOR_VERSION >= 3

cdef from_string_and_size(const char* s, size_t length):
    if IS_PYTHON3:
        return s[:length].decode("ascii")
    else:
        return s[:length]

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
    if not IS_PYTHON3:
        return s
    elif s is None:
        return None
    elif PyBytes_Check(s):
        return s
    elif PyUnicode_Check(s):
        return s.encode(encoding)
    else:
        raise TypeError(u"Argument must be string, bytes or unicode.")

cdef charptr_to_str(const char* s, encoding="ascii"):
    if s == NULL:
        return None
    if PY_MAJOR_VERSION < 3:
        return s
    else:
        return s.decode(encoding)

cdef charptr_to_str_w_len(const char* s, size_t n, encoding="ascii"):
    if s == NULL:
        return None
    if PY_MAJOR_VERSION < 3:
        return s[:n]
    else:
        return s[:n].decode(encoding)

cdef bytes charptr_to_bytes(const char* s, encoding="ascii"):
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


def _pysam_dispatch(collection,
                    method,
                    args=None,
                    catch_stdout=True,
                    save_stdout=None):
    '''call ``method`` in samtools/bcftools providing arguments in args.
    
    Catching of stdout can be turned off by setting *catch_stdout* to
    False.

    '''

    if method == "index":
        if not os.path.exists(args[0]):
            raise IOError("No such file or directory: '%s'" % args[0])
            
    if args is None:
        args = []
    else:
        args = list(args)

    # redirect stderr to file
    stderr_h, stderr_f = tempfile.mkstemp()
    pysam_set_stderr(stderr_h)

    # redirect stdout to file
    if save_stdout:
        stdout_f = save_stdout
        stdout_h = c_open(force_bytes(stdout_f),
                          O_WRONLY)
        if stdout_h == -1:
            raise OSError("error while opening {} for writing".format(stdout_f))

        pysam_set_stdout_fn(force_bytes(stdout_f))
        pysam_set_stdout(stdout_h)
    elif catch_stdout:
        stdout_h, stdout_f = tempfile.mkstemp()

        MAP_STDOUT_OPTIONS = {
            "samtools": {
                "view": "-o {}",
                "mpileup": "-o {}",
                "depad": "-o {}",
                "calmd": "",  # uses pysam_stdout_fn
            },
            "bcftools": {}
        }

        stdout_option = None
        if collection == "bcftools":
            # in bcftools, most methods accept -o, the exceptions
            # are below:
            if method not in ("index", "roh", "stats"):
                stdout_option = "-o {}"
        elif method in MAP_STDOUT_OPTIONS[collection]:
            stdout_option = MAP_STDOUT_OPTIONS[collection][method]

        if stdout_option is not None:
            os.close(stdout_h)
            pysam_set_stdout_fn(force_bytes(stdout_f))
            args.extend(stdout_option.format(stdout_f).split(" "))
        else:
            pysam_set_stdout(stdout_h)
    else:
        pysam_set_stdout_fn("-")

    # setup the function call to samtools/bcftools main
    cdef char ** cargs
    cdef int i, n, retval, l
    n = len(args)
    method = force_bytes(method)
    collection = force_bytes(collection)
    args = [force_bytes(a) for a in args]

    # allocate two more for first (dummy) argument (contains command)
    cdef int extra_args = 0
    if method == b"index":
        extra_args = 1
    # add extra arguments for commands accepting optional arguments
    # such as 'samtools index x.bam [out.index]'
    cargs = <char**>calloc(n + 2 + extra_args, sizeof(char *))
    cargs[0] = collection
    cargs[1] = method

    # create copies of strings - getopt for long options permutes
    # arguments
    for i from 0 <= i < n:
        l = len(args[i])
        cargs[i + 2] = <char *>calloc(l + 1, sizeof(char))
        strncpy(cargs[i + 2], args[i], l)
    
    # reset getopt. On OsX there getopt reset is different
    # between getopt and getopt_long
    if method in [b'index', b'cat', b'quickcheck',
                  b'faidx', b'kprobaln']:
        set_optind(1)
    else:
        set_optind(0)

    # call samtools/bcftools
    if collection == b"samtools":
        retval = samtools_main(n + 2, cargs)
    elif collection == b"bcftools":
        retval = bcftools_main(n + 2, cargs)

    for i from 0 <= i < n:
        free(cargs[i + 2])
    free(cargs)

    # get error messages
    def _collect(fn):
        out = []
        try:
            with open(fn, "r") as inf:
                out = inf.read()
        except UnicodeDecodeError:
            with open(fn, "rb") as inf:
                # read binary output
                out = inf.read()
        finally:
            os.remove(fn)
        return out

    pysam_unset_stderr()
    out_stderr = _collect(stderr_f)

    if save_stdout:
        pysam_unset_stdout()
        out_stdout = None
    elif catch_stdout:
        pysam_unset_stdout()
        out_stdout = _collect(stdout_f)
    else:
        out_stdout = None

    return retval, out_stderr, out_stdout


__all__ = ["qualitystring_to_array",
           "array_to_qualitystring",
           "qualities_to_qualitystring"]
