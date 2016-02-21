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


@contextmanager
def stdout_redirector(to=os.devnull):
    '''
    import os

    with stdout_redirected(to=filename):
        print("from Python")
        os.system("echo non-Python applications are also supported")

    see http://stackoverflow.com/questions/5081657/how-do-i-prevent-a-c-shared-library-to-print-on-stdout-in-python/17954769#17954769
    '''
    fd = sys.stdout.fileno()

    def _redirect_stdout(to):
        # flush C-level stdout
        try:
            fflush(c_stdout)
            sys.stdout.close()
        except (OSError, IOError):
            # some tools close stdout
            # Py3: OSError
            # Py2: IOError
            pass

        # fd writes to 'to' file
        os.dup2(to.fileno(), fd)
        # Python writes to fd
        if IS_PYTHON3:
            sys.stdout = io.TextIOWrapper(
                os.fdopen(fd, 'wb'))
        else:
            sys.stdout = os.fdopen(fd, 'w')
        
    with os.fdopen(os.dup(fd), 'w') as old_stdout:
        _redirect_stdout(to)
        try:
            yield # allow code to be run with the redirected stdout
        finally:
            _redirect_stdout(old_stdout)
            # restore stdout.
            # buffering and flags may be different

# def stdout_redirector(stream):
#     """
#     See discussion in:

#     http://eli.thegreenplace.net/2015/redirecting-all-kinds-of-stdout-in-python/
#     """

#     # The original fd stdout points to. Usually 1 on POSIX systems.
#     original_stdout_fd = sys.stdout.fileno()
#     print ("original_fd=", original_stdout_fd)
#     def _redirect_stdout(to_fd):
#         """Redirect stdout to the given file descriptor."""
#         # Flush the C-level buffer stdout
#         fflush(c_stdout)
#         # Flush and close sys.stdout - also closes the file descriptor
#         # (fd)
#         sys.stdout.close()
#         # Make original_stdout_fd point to the same file as to_fd
#         os.dup2(to_fd, original_stdout_fd)
#         # Create a new sys.stdout that points to the redirected fd
#         if IS_PYTHON3:
#             sys.stdout = io.TextIOWrapper(
#                 os.fdopen(original_stdout_fd, 'wb'))

#     # Save a copy of the original stdout fd in saved_stdout_fd
#     saved_stdout_fd = os.dup(original_stdout_fd)
#     try:
#         # Create a temporary file and redirect stdout to it
#         tfile = tempfile.TemporaryFile(mode='w+b')
#         _redirect_stdout(tfile.fileno())
#         # Yield to caller, then redirect stdout back to the saved fd
#         yield
#         _redirect_stdout(saved_stdout_fd)
#         # Copy contents of temporary file to the given stream
#         tfile.flush()
#         tfile.seek(0, io.SEEK_SET)
#         stream.write(tfile.read())
#     finally:
#         tfile.close()
#         os.close(saved_stdout_fd)


def _pysam_dispatch(collection,
                    method,
                    args=(),
                    catch_stdout=True):
    '''call ``method`` in samtools/bcftools providing arguments in args.
    
    .. note:: 
       This method redirects stdout to capture it 
       from samtools. If for some reason stdout disappears
       the reason might be in this method.

    .. note::
       This method captures stdout and stderr using temporary files,
       which are then read into memory in their entirety. This method
       is slow and might cause large memory overhead.

    Catching of stdout can be turned of by setting *catch_stdout* to
    False.

    See http://bytes.com/topic/c/answers/487231-how-capture-stdout-temporarily
    on the topic of redirecting stderr/stdout.

    '''

    # note that debugging this module can be a problem
    # as stdout/stderr will not appear on the terminal
    # some special cases
    if method == "index":
        if not os.path.exists(args[0]):
            raise IOError("No such file or directory: '%s'" % args[0])

    # redirect stderr and stdout to file
    stderr_h, stderr_f = tempfile.mkstemp()
    pysam_set_stderr(stderr_h)

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
    if catch_stdout:
        with tempfile.TemporaryFile(mode='w+b') as tfile:
            with stdout_redirector(tfile):
                if collection == b"samtools":
                    retval = samtools_main(n + 2, cargs)
                elif collection == b"bcftools":
                    retval = bcftools_main(n + 2, cargs)
            tfile.flush()
            tfile.seek(0)
            # do not force str, as output might be binary,
            # for example BAM, VCF.gz, etc.
            out_stdout = tfile.read()
    else:
        if collection == b"samtools":
            retval = samtools_main(n + 2, cargs)
        elif collection == b"bcftools":
            retval = bcftools_main(n + 2, cargs)
        out_stdout = None

    for i from 0 <= i < n:
        free(cargs[i + 2])
    free(cargs)

    # get error messages
    pysam_unset_stderr()
    out_stderr = []
    try:
        with open(stderr_f, "r") as inf:
            out_stderr = inf.readlines()
    except UnicodeDecodeError:
        with open( stderr_f, "rb") as inf:
            # read binary output
            out_stderr = inf.read()
    finally:
        os.remove(stderr_f)

    return retval, out_stderr, out_stdout


__all__ = ["qualitystring_to_array",
           "array_to_qualitystring",
           "qualities_to_qualitystring"]
