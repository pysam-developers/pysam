# cython: embedsignature=True
# cython: profile=True
# adds doc-strings for sphinx
import tempfile
import os
import sys
import platform
from cpython cimport PyBytes_Check, PyUnicode_Check
from cpython.version cimport PY_MAJOR_VERSION

########################################################################
########################################################################
########################################################################
## Python 3 compatibility functions
########################################################################
IS_PYTHON3 = PY_MAJOR_VERSION >= 3

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


def _samtools_dispatch(method,
                       args = (),
                       catch_stdout = True):
    '''call ``method`` in samtools providing arguments in args.
    
    .. note:: 
       This method redirects stdout to capture it 
       from samtools. If for some reason stdout disappears
       the reason might be in this method.

    .. note::
       The current implementation might only work on linux.

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
            stdout_save = Outs( sys.stdout.fileno() )
            stdout_save.setfd( stdout_h )
        except AttributeError:
            # stdout has already been redirected
            catch_stdout = False

        # patch for `samtools view`
        # samtools `view` closes stdout, from which I can not
        # recover. Thus redirect output to file with -o option.
        if method == "view":
            if "-o" in args:
                raise ValueError("option -o is forbidden in samtools view")
            args = ( "-o", stdout_f ) + args

    # do the function call to samtools
    cdef char ** cargs
    cdef int i, n, retval

    n = len(args)
    method = _forceCmdlineBytes(method)
    args = [ _forceCmdlineBytes(a) for a in args ]

    # allocate two more for first (dummy) argument (contains command)
    cargs = <char**>calloc( n+2, sizeof( char *) )
    cargs[0] = "samtools"
    cargs[1] = method
    for i from 0 <= i < n: cargs[i+2] = args[i]
    
    retval = pysam_dispatch(n+2, cargs)
    free( cargs )
    
    # restore stdout/stderr. This will also flush, so
    # needs to be before reading back the file contents
    if catch_stdout:
        stdout_save.restore()
        try:
            with open( stdout_f, "r") as inf:
                out_stdout = inf.readlines()
        except UnicodeDecodeError:
            with open( stdout_f, "rb") as inf:
                # read binary output
                out_stdout = inf.read()
        os.remove( stdout_f )
    else:
        out_stdout = []

    # get error messages
    pysam_unset_stderr()
    try:
        with open( stderr_f, "r") as inf:
            out_stderr = inf.readlines()
    except UnicodeDecodeError:
        with open( stderr_f, "rb") as inf:
            # read binary output
            out_stderr = inf.read()
    else:
        out_stderr = []
    finally:
        os.remove( stderr_f )

    return retval, out_stderr, out_stdout

