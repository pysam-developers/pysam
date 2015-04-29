from pysam.libchtslib import *

import pysam.ctabix as ctabix
from pysam.ctabix import *
import pysam.csamfile as csamfile
from pysam.csamfile import *
import pysam.calignmentfile as calignmentfile
from pysam.calignmentfile import *
import pysam.cfaidx as cfaidx
from pysam.cfaidx import *
import pysam.cvcf as cvcf
from pysam.cvcf import *
import pysam.cbcf as cbcf
from pysam.cbcf import *
import pysam.csamtools as csamtools

import pysam.Pileup as Pileup
import os


class SamtoolsError(Exception):
    '''exception raised in case of an error incurred in the samtools
    library.'''

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class SamtoolsDispatcher(object):
    '''samtools dispatcher.

    Emulates the samtools command line as module calls.

    Captures stdout and stderr.

    Raises a :class:`pysam.SamtoolsError` exception in case
    samtools exits with an error code other than 0.

    Some command line options are associated with parsers.  For
    example, the samtools command "pileup -c" creates a tab-separated
    table on standard output. In order to associate parsers with
    options, an optional list of parsers can be supplied. The list
    will be processed in order checking for the presence of each
    option.

    If no parser is given or no appropriate parser is found, the
    stdout output of samtools commands will be returned.

    '''

    dispatch = None
    parsers = None

    def __init__(self, dispatch, parsers):
        self.dispatch = dispatch
        self.parsers = parsers
        self.stderr = []

    def __call__(self, *args, **kwargs):
        '''execute a samtools command
        '''
        retval, stderr, stdout = csamtools._samtools_dispatch(
            self.dispatch, args, catch_stdout=kwargs.get("catch_stdout", True))

        if retval:
            raise SamtoolsError(
                'csamtools returned with error %i: %s' %
                (retval, "\n".join(stderr)))
        self.stderr = stderr
        # samtools commands do not propagate the return code correctly.
        # I have thus added this patch to throw if there is output on stderr.
        # Note that there is sometimes output on stderr that is not an error,
        # for example: [sam_header_read2] 2 sequences loaded.
        # Ignore messages like these
        stderr = [x for x in stderr
                  if not (x.startswith("[sam_header_read2]") or
                          x.startswith("[bam_index_load]") or
                          x.startswith("[bam_sort_core]") or
                          x.startswith("[samopen] SAM header is present"))]
        if stderr:
            raise SamtoolsError("\n".join(stderr))

        # call parser for stdout:
        if not kwargs.get("raw") and stdout and self.parsers:
            for options, parser in self.parsers:
                for option in options:
                    if option not in args:
                        break
                else:
                    return parser(stdout)

        return stdout

    def get_messages(self):
        return self.stderr

    def usage(self):
        '''return the samtools usage information for this command'''
        retval, stderr, stdout = csamtools._samtools_dispatch(
            self.dispatch)
        return "".join(stderr)

#
# samtools command line options to export in python
#
# import is a python reserved word.
SAMTOOLS_DISPATCH = {
    # samtools 'documented' commands
    "view": ("view", None),
    "sort": ("sort", None),
    "mpileup": ("mpileup", None),
    "depth": ("depth", None),
    "faidx": ("faidx", None),
    "tview": ("tview", None),
    "index": ("index", None),
    "idxstats": ("idxstats", None),
    "fixmate": ("fixmate", None),
    "flagstat": ("flagstat", None),
    "calmd": ("calmd", None),
    "merge": ("merge", None),
    "rmdup": ("rmdup", None),
    "reheader": ("reheader", None),
    "cat": ("cat", None),
    "targetcut": ("targetcut", None),
    "phase": ("phase", None),
    # others
    "samimport": ("import", None),
    "bam2fq": ("bam2fq", None),
    "pad2unpad": ("pad2unpad", None),
    "depad": ("pad2unpad", None),
    "bedcov": ("bedcov", None),
    "bamshuf": ("bamshuf", None),
    # obsolete
    # "pileup": "pileup", ( (("-c",), Pileup.iterate),),),
}

# instantiate samtools commands as python functions
for key, options in SAMTOOLS_DISPATCH.items():
    cmd, parser = options
    globals()[key] = SamtoolsDispatcher(cmd, parser)

# hack to export all the symbols from separate modules
__all__ = \
    libchtslib.__all__ + \
    ctabix.__all__ + \
    cvcf.__all__ +\
    cbcf.__all__ +\
    cfaidx.__all__ +\
    calignmentfile.__all__ +\
    csamfile.__all__ +\
    ["SamtoolsError", "SamtoolsDispatcher"] +\
    list(SAMTOOLS_DISPATCH) +\
    ["Pileup"]

from pysam.version import __version__, __samtools_version__


def get_include():
    '''return a list of include directories.'''
    dirname = os.path.abspath(os.path.join(os.path.dirname(__file__)))
    return [dirname,
            os.path.join(dirname, 'include', 'htslib'),
            os.path.join(dirname, 'include', 'samtools')]


def get_defines():
    '''return a list of defined compilation parameters.'''
    return [('_FILE_OFFSET_BITS', '64'),
            ('_USE_KNETFILE', '')]


def get_libraries():
    '''return a list of libraries to link against.'''
    # Note that this list does not include csamtools.so as there are
    # numerous name conflicts with libchtslib.so.
    dirname = os.path.abspath(os.path.join(os.path.dirname(__file__)))
    return [os.path.join(dirname, x) for x in (
        'libchtslib.so',
        'TabProxies.so',
        'cfaidx.so',
        'csamfile.so',
        'cvcf.so',
        'cbcf.so',
        'ctabix.so')]
