# cython: language_level=3
import types
import sys
import string
import re
import tempfile
import os
import io
from contextlib import contextmanager
from codecs import register_error

from cpython.version cimport PY_MAJOR_VERSION, PY_MINOR_VERSION
from cpython cimport PyBytes_Check, PyUnicode_Check
from cpython cimport array as c_array
from libc.stdlib cimport calloc, free
from libc.string cimport strncpy
from libc.stdint cimport INT32_MAX, int32_t
from libc.stdio cimport fprintf, stderr, fflush
from libc.stdio cimport stdout as c_stdout
from posix.fcntl cimport open as c_open, O_WRONLY

from libcsamtools cimport samtools_dispatch, samtools_set_stdout, samtools_set_stderr, \
    samtools_close_stdout, samtools_close_stderr, samtools_set_stdout_fn

from libcbcftools cimport bcftools_dispatch, bcftools_set_stdout, bcftools_set_stderr, \
    bcftools_close_stdout, bcftools_close_stderr, bcftools_set_stdout_fn

#####################################################################
# hard-coded constants
cdef int MAX_POS = (1 << 31) - 1

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
    return force_str(result.tobytes())


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
## String encoding configuration facilities
########################################################################

# Codec error handler that just interprets each bad byte as ISO-8859-1.
def latin1_replace(exception):
    return (chr(exception.object[exception.start]), exception.end)

register_error('pysam.latin1replace', latin1_replace)


cdef str ERROR_HANDLER = 'strict'

cpdef get_encoding_error_handler():
    return ERROR_HANDLER

cpdef set_encoding_error_handler(name):
    global ERROR_HANDLER
    previous = ERROR_HANDLER
    ERROR_HANDLER = name
    return previous

########################################################################
## Python 3 compatibility functions
########################################################################

cdef from_string_and_size(const char* s, size_t length):
    return s[:length].decode('utf-8', ERROR_HANDLER)

# filename encoding (adapted from lxml.etree.pyx)
cdef str FILENAME_ENCODING = sys.getfilesystemencoding() or sys.getdefaultencoding() or 'ascii'
cdef str TEXT_ENCODING = 'utf-8'

cdef bytes encode_filename(object filename):
    """Make sure a filename is 8-bit encoded (or None)."""
    if filename is None:
        return None
    return os.fsencode(filename)


cdef bytes force_bytes(object s, encoding=None, errors=None):
    """convert string or unicode object to bytes, assuming
    utf8 encoding.
    """
    if s is None:
        return None
    elif PyBytes_Check(s):
        return s
    elif PyUnicode_Check(s):
        return s.encode(encoding or TEXT_ENCODING, errors or ERROR_HANDLER)
    else:
        raise TypeError("Argument must be string, bytes or unicode.")


cdef charptr_to_str(const char* s, encoding=None, errors=None):
    if s == NULL:
        return None
    return s.decode(encoding or TEXT_ENCODING, errors or ERROR_HANDLER)


cdef charptr_to_str_w_len(const char* s, size_t n, encoding=None, errors=None):
    if s == NULL:
        return None
    return s[:n].decode(encoding or TEXT_ENCODING, errors or ERROR_HANDLER)


cdef bytes charptr_to_bytes(const char* s, encoding=None, errors=None):
    if s == NULL:
        return None
    else:
        return s


cdef force_str(object s, encoding=None, errors=None):
    """Return s converted to str type of current Python
    (bytes in Py2, unicode in Py3)"""
    if s is None:
        return None
    if PyBytes_Check(s):
        return s.decode(encoding or TEXT_ENCODING, errors or ERROR_HANDLER)
    # assume unicode
    return s


cdef decode_bytes(bytes s, encoding=None, errors=None):
    """Return s converted to current Python's str type,
    always decoding even in Python 2"""
    if s is None:
        return None
    else:
        return s.decode(encoding or TEXT_ENCODING, errors or ERROR_HANDLER)


cpdef parse_region(contig=None,
                   start=None,
                   stop=None,
                   region=None,
                   reference=None,
                   end=None):
    """parse alternative ways to specify a genomic region. A region can
    either be specified by :term:`reference`, `start` and
    `end`. `start` and `end` denote 0-based, half-open intervals.
    
    :term:`reference` and `end` are also accepted for backward
    compatibility as synonyms for :term:`contig` and `stop`,
    respectively.

    Alternatively, a samtools :term:`region` string can be supplied.

    If any of the coordinates are missing they will be replaced by the
    minimum (`start`) or maximum (`end`) coordinate.

    Note that region strings are 1-based, while `start` and `end`
    denote an interval in python coordinates.

    Returns
    -------

    tuple :  a tuple of `reference`, `start` and `end`.

    Raises
    ------

    ValueError
       for invalid or out of bounds regions.

    """
    cdef int32_t rstart
    cdef int32_t rstop

    
    if reference is not None:
        if contig is not None:
           raise ValueError('contig and reference should not both be specified')
        contig = reference

    if contig is not None and region is not None:
        raise ValueError('contig/reference and region should not both be specified')
        
    if end is not None:
        if stop is not None:
            raise ValueError('stop and end should not both be specified')
        stop = end

    if contig is None and region is None:
        raise ValueError("neither contig nor region are given")

    rstart = 0
    rstop = MAX_POS
    if start is not None:
        try:
            rstart = start
        except OverflowError:
            raise ValueError('start out of range (%i)' % start)

    if stop is not None:
        try:
            rstop = stop
        except OverflowError:
            raise ValueError('stop out of range (%i)' % stop)

    if region:
        if ":" in region:
            contig, coord = region.split(":")
            parts = coord.split("-")
            rstart = int(parts[0]) - 1
            if len(parts) >= 1:
                rstop = int(parts[1])
        else:
            contig = region

    if rstart > rstop:
        raise ValueError('invalid coordinates: start (%i) > stop (%i)' % (rstart, rstop))
    if not 0 <= rstart < MAX_POS:
        raise ValueError('start out of range (%i)' % rstart)
    if not 0 <= rstop <= MAX_POS:
        raise ValueError('stop out of range (%i)' % rstop)

    return contig, rstart, rstop


def _pysam_dispatch(collection,
                    method,
                    args=None,
                    catch_stdout=True,
                    is_usage=False,
                    save_stdout=None):
    '''call ``method`` in samtools/bcftools providing arguments in args.
    
    By default, stdout is redirected to a temporary file using the patched
    C sources except for a few commands that have an explicit output option
    (typically: -o). In these commands (such as samtools view), this explicit
    option is used. If *is_usage* is True, then these explicit output options
    will not be used.

    Catching of stdout can be turned off by setting *catch_stdout* to
    False.
    '''

    if method == "index" and args:
        # We make sure that at least the first specified input file exists,
        # and if it doesn't we raise an IOError.
        ARGUMENTS = ['-m', '--min-shift', '-o', '--output', '--output-file', '-@', '--threads']
        skip_next = False
        for arg in args:
            if skip_next:
                skip_next = False
                continue
            if arg.startswith('-'):
                # Skip next argument for e.g. '--min-shift' '12' or '-m' '12' but not '-m12'
                if arg in ARGUMENTS:
                    skip_next = True
                continue
            if not os.path.exists(arg):
                raise IOError("No such file or directory: '%s'" % arg)
            else:
                break
            
    if args is None:
        args = []
    else:
        args = list(args)

    # redirect stderr to file
    stderr_h, stderr_f = tempfile.mkstemp()
        
    # redirect stdout to file
    if save_stdout:
        stdout_f = save_stdout
        stdout_h = c_open(force_bytes(stdout_f),
                          O_WRONLY)
        if stdout_h == -1:
            raise IOError("error while opening {} for writing".format(stdout_f))

        samtools_set_stdout_fn(force_bytes(stdout_f))
        bcftools_set_stdout_fn(force_bytes(stdout_f))
            
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
            # special case - samtools view -c outputs on stdout
            if not(method == "view" and "-c" in args):
                stdout_option = MAP_STDOUT_OPTIONS[collection][method]

        if stdout_option is not None and not is_usage:
            os.close(stdout_h)
            samtools_set_stdout_fn(force_bytes(stdout_f))
            bcftools_set_stdout_fn(force_bytes(stdout_f))
            args.extend(stdout_option.format(stdout_f).split(" "))
            stdout_h = c_open(b"/dev/null", O_WRONLY)
    else:
        samtools_set_stdout_fn("-")
        bcftools_set_stdout_fn("-")
        stdout_h = c_open(b"/dev/null", O_WRONLY)

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

    # call samtools/bcftools
    if collection == b"samtools":
        samtools_set_stdout(stdout_h)
        samtools_set_stderr(stderr_h)
        retval = samtools_dispatch(n + 2, cargs)
        samtools_close_stdout()
        samtools_close_stderr()
    elif collection == b"bcftools":
        bcftools_set_stdout(stdout_h)
        bcftools_set_stderr(stderr_h)
        retval = bcftools_dispatch(n + 2, cargs)
        bcftools_close_stdout()
        bcftools_close_stderr()
    else:
        # unknown -- just return a Unix shell's "command not found" exit status
        retval = 127

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

    out_stderr = _collect(stderr_f)
    if save_stdout:
        out_stdout = None
    elif catch_stdout:
        out_stdout = _collect(stdout_f)
    else:
        out_stdout = None
        
    return retval, out_stderr, out_stdout


__all__ = [
    "qualitystring_to_array",
    "array_to_qualitystring",
    "qualities_to_qualitystring",
    "get_encoding_error_handler",
    "set_encoding_error_handler",
]
