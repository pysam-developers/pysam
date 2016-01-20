import types
import sys
import string
import re
import tempfile
import os

from cpython.version cimport PY_MAJOR_VERSION
from cpython cimport PyBytes_Check, PyUnicode_Check
from cpython cimport array as c_array
from libc.stdlib cimport calloc, free
from libc.string cimport strncpy

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


class Outs:
    '''http://mail.python.org/pipermail/python-list/2000-June/038406.html'''
    def __init__(self, id = 1):
        self.streams = []
        self.id = id

    def setdevice(self, filename):
        '''open an existing file, like "/dev/null"'''
        fd = os.open(filename, os.O_WRONLY)
        self.setfd(fd)

    def setfile(self, filename):
        '''open a new file.'''
        fd = os.open(filename, os.O_WRONLY|os.O_CREAT, 0660)
        self.setfd(fd)

    def setfd(self, fd):
        ofd = os.dup(self.id)      #  Save old stream on new unit.
        self.streams.append(ofd)
        sys.stdout.flush()          #  Buffered data goes to old stream.
        sys.stderr.flush()          #  Buffered data goes to old stream.
        os.dup2(fd, self.id)        #  Open unit 1 on new stream.
        os.close(fd)                #  Close other unit (look out, caller.)

    def restore(self):
        '''restore previous output stream'''
        if self.streams:
            # the following was not sufficient, hence flush both stderr and stdout
            # os.fsync( self.id )
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(self.streams[-1], self.id)
            os.close(self.streams[-1])
            del self.streams[-1]


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

    Catching of stdout can be turned of by setting *catch_stdout* to False.

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
                                        
    if catch_stdout:
        stdout_h, stdout_f = tempfile.mkstemp()
        try:
            stdout_save = Outs(sys.stdout.fileno())
            stdout_save.setfd(stdout_h)
        except AttributeError:
            # stdout has already been redirected
            catch_stdout = False

        # patch for `samtools view`
        # samtools `view` closes stdout, from which I can not
        # recover. Thus redirect output to file with -o option.
        if collection == "samtools" and method == "view":
            if "-o" in args:
                raise ValueError("option -o is forbidden in samtools view")
            args =  ("-o", stdout_f) + args

    # do the function call to samtools
    cdef char ** cargs
    cdef int i, n, retval, l

    n = len(args)
    method = force_cmdline_bytes(method)
    collection = force_cmdline_bytes(collection)
    args = [force_cmdline_bytes(a) for a in args]

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
    
    # reset getopt
    reset_getopt()

    # call samtools
    if collection == b"samtools":
        retval = samtools_main(n + 2, cargs)
    elif collection == b"bcftools":
        retval = bcftools_main(n + 2, cargs)

    for i from 0 <= i < n:
        free(cargs[i + 2])
    free(cargs)
    # restore stdout/stderr. This will also flush, so
    # needs to be before reading back the file contents
    if catch_stdout:
        stdout_save.restore()
        try:
            with open(stdout_f, "r") as inf:
                out_stdout = inf.readlines()
        except UnicodeDecodeError:
            with open(stdout_f, "rb") as inf:
                # read binary output
                out_stdout = inf.read()
        os.remove(stdout_f)
    else:
        out_stdout = []

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
