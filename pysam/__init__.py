import os
import sysconfig

from pysam.libchtslib import *
import pysam.libchtslib as libchtslib
from pysam.libcsamtools import *
from pysam.libcbcftools import *
from pysam.libcutils import *
import pysam.libcutils as libcutils
import pysam.libcfaidx as libcfaidx
from pysam.libcfaidx import *
import pysam.libctabix as libctabix
from pysam.libctabix import *
import pysam.libctabixproxies as libctabixproxies
from pysam.libctabixproxies import *
import pysam.libcsamfile as libcsamfile
from pysam.libcsamfile import *
import pysam.libcalignmentfile as libcalignmentfile
from pysam.libcalignmentfile import *
import pysam.libcalignedsegment as libcalignedsegment
from pysam.libcalignedsegment import *
import pysam.libcvcf as libcvcf
from pysam.libcvcf import *
import pysam.libcbcf as libcbcf
from pysam.libcbcf import *
import pysam.libcbgzf as libcbgzf
from pysam.libcbgzf import *
from pysam.utils import SamtoolsError
import pysam.Pileup as Pileup
from pysam.samtools import *
import pysam.config


# export all the symbols from separate modules
__all__ = (
    libchtslib.__all__ +  # type: ignore
    libcutils.__all__ +  # type: ignore
    libctabix.__all__ +  # type: ignore
    libcvcf.__all__ +  # type: ignore
    libcbcf.__all__ +  # type: ignore
    libcbgzf.__all__ +  # type: ignore
    libcfaidx.__all__ +  # type: ignore
    libctabixproxies.__all__ +  # type: ignore
    libcalignmentfile.__all__ +  # type: ignore
    libcalignedsegment.__all__ +  # type: ignore
    libcsamfile.__all__ +  # type: ignore
    ["SamtoolsError"] +
    ["Pileup"]
)
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
    # ('_FILE_OFFSET_BITS', '64'),
    # ('_USE_KNETFILE', '')]
    return []


def get_libraries():
    '''return a list of libraries to link against.'''
    # Note that this list does not include libcsamtools.so as there are
    # numerous name conflicts with libchtslib.so.
    dirname = os.path.abspath(os.path.join(os.path.dirname(__file__)))
    pysam_libs = ['libctabixproxies',
                  'libcfaidx',
                  'libcsamfile',
                  'libcvcf',
                  'libcbcf',
                  'libctabix']
    if pysam.config.HTSLIB == "builtin":
        pysam_libs.append('libchtslib')

    so = sysconfig.get_config_var('EXT_SUFFIX') or sysconfig.get_config_var('SO')
    return [os.path.join(dirname, x + so) for x in pysam_libs]
