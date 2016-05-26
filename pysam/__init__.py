import os
import sys
import sysconfig

from pysam.libchtslib import *
from pysam.cutils import *
import pysam.cutils as cutils
import pysam.cfaidx as cfaidx
from pysam.cfaidx import *
import pysam.ctabix as ctabix
from pysam.ctabix import *
import pysam.csamfile as csamfile
from pysam.csamfile import *
import pysam.calignmentfile as calignmentfile
from pysam.calignmentfile import *
import pysam.calignedsegment as calignedsegment
from pysam.calignedsegment import *
import pysam.cvcf as cvcf
from pysam.cvcf import *
import pysam.cbcf as cbcf
from pysam.cbcf import *
from pysam.utils import SamtoolsError
import pysam.Pileup as Pileup
from pysam.samtools import *
import pysam.config


# export all the symbols from separate modules
__all__ = \
    libchtslib.__all__ +\
    cutils.__all__ +\
    ctabix.__all__ +\
    cvcf.__all__ +\
    cbcf.__all__ +\
    cfaidx.__all__ +\
    calignmentfile.__all__ +\
    calignedsegment.__all__ +\
    csamfile.__all__ +\
    ["SamtoolsError"] +\
    ["Pileup"]

from pysam.version import __version__, __samtools_version__


def get_include():
    '''return a list of include directories.'''
    dirname = os.path.abspath(os.path.join(os.path.dirname(__file__)))

    #
    # Header files may be stored in different relative locations
    # depending on installation mode (e.g., `python setup.py install`,
    # `python setup.py develop`. The first entry in each list is
    # where develop-mode headers can be found.
    #
    htslib_possibilities = [os.path.join(dirname, '..', 'htslib'),
                            os.path.join(dirname, 'include', 'htslib')]
    samtool_possibilities = [os.path.join(dirname, '..', 'samtools'),
                             os.path.join(dirname, 'include', 'samtools')]

    includes = [dirname]
    for header_locations in [htslib_possibilities, samtool_possibilities]:
        for header_location in header_locations:
            if os.path.exists(header_location):
                includes.append(os.path.abspath(header_location))
                break

    return includes


def get_defines():
    '''return a list of defined compilation parameters.'''
    return [] #('_FILE_OFFSET_BITS', '64'),
    # ('_USE_KNETFILE', '')]


def get_libraries():
    '''return a list of libraries to link against.'''
    # Note that this list does not include csamtools.so as there are
    # numerous name conflicts with libchtslib.so.
    dirname = os.path.abspath(os.path.join(os.path.dirname(__file__)))
    pysam_libs = ['ctabixproxies',
                  'cfaidx',
                  'csamfile',
                  'cvcf',
                  'cbcf',
                  'ctabix']
    if pysam.config.HTSLIB == "builtin":
        pysam_libs.append('libchtslib')

    if sys.version_info.major >= 3:
        if sys.version_info.minor >= 5:
            return [os.path.join(dirname, x + ".{}.so".format(
                sysconfig.get_config_var('SOABI'))) for x in pysam_libs]
        else:
            return [os.path.join(dirname, x + ".{}{}.so".format(
                sys.implementation.cache_tag,
                sys.abiflags)) for x in pysam_libs]
    else:
        return [os.path.join(dirname, x + ".so") for x in pysam_libs]
