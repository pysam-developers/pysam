#! /usr/bin/python

'''pysam - a python module for reading, manipulating and writing
genomic data sets.

pysam is a lightweight wrapper of the htslib C-API and provides
facilities to read and write SAM/BAM/VCF/BCF/BED/GFF/GTF/FASTA/FASTQ
files as well as access to the command line functionality of the
samtools and bcftools packages. The module supports compression and
random access through indexing.

This module provides a low-level wrapper around the htslib C-API as
using cython and a high-level API for convenient access to the data
within standard genomic file formats.

The current version wraps htslib-1.3.1, samtools-1.3.1 and bcftools-1.3.1.

See:
http://www.htslib.org
https://github.com/pysam-developers/pysam
http://pysam.readthedocs.org/en/stable

'''

import collections
import glob
import os
import platform
import re
import subprocess
import sys
import sysconfig
from contextlib import contextmanager
from setuptools import Extension, setup

IS_PYTHON3 = sys.version_info.major >= 3


@contextmanager
def changedir(path):
    save_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(save_dir)


def run_configure(option):
    try:
        retcode = subprocess.call(
            " ".join(("./configure", option)),
            shell=True)
        if retcode != 0:
            return False
        else:
            return True
    except OSError as e:
        return False


def run_make_print_config():
    stdout = subprocess.check_output(["make", "print-config"])
    if IS_PYTHON3:
        stdout = stdout.decode("ascii")

    result = dict([[x.strip() for x in line.split("=")]
                   for line in stdout.splitlines()])
    return result


def configure_library(library_dir, env_options=None, options=[]):

    configure_script = os.path.join(library_dir, "configure")

    if not os.path.exists(configure_script):
        raise ValueError(
            "configure script {} does not exist".format(configure_script))

    with changedir(library_dir):
        if env_options is not None:
            if run_configure(env_options):
                return env_options

        for option in options:
            if run_configure(option):
                return option

    return None


def distutils_dir_name(dname):
    """Returns the name of a distutils build directory
    see: http://stackoverflow.com/questions/14320220/
               testing-python-c-libraries-get-build-path
    """
    f = "{dirname}.{platform}-{version[0]}.{version[1]}"
    return f.format(dirname=dname,
                    platform=sysconfig.get_platform(),
                    version=sys.version_info)

# How to link against HTSLIB
# separate: use included htslib and include in each extension
#           module. No dependencies between modules and works
#           with setup.py install, but wasteful in terms of
#           memory and compilation time.
# shared: share chtslib across extension modules. This would be
#         the ideal method, but currently requires
#         LD_LIBRARY_PATH to be set correctly when using
#         pysam.
# external: use shared libhts.so compiled outside of
#           pysam
HTSLIB_MODE = os.environ.get("HTSLIB_MODE", "shared")
HTSLIB_LIBRARY_DIR = os.environ.get("HTSLIB_LIBRARY_DIR", None)
HTSLIB_INCLUDE_DIR = os.environ.get("HTSLIB_INCLUDE_DIR", None)
HTSLIB_CONFIGURE_OPTIONS = os.environ.get("HTSLIB_CONFIGURE_OPTIONS", None)
HTSLIB_SOURCE = None

package_list = ['pysam',
                'pysam.include',
                'pysam.include.samtools',
                'pysam.include.bcftools',
                'pysam.include.samtools.win32']
package_dirs = {'pysam': 'pysam',
                'pysam.include.samtools': 'samtools',
                'pysam.include.bcftools': 'bcftools'}
config_headers = ["samtools/config.h"]

from cy_build import CyExtension as Extension, cy_build_ext as build_ext

cmdclass = {'build_ext': build_ext}

# Check if cython is available
#
# If cython is available, the pysam will be built using cython from
# the .pyx files. If no cython is available, the C-files included in the
# distribution will be used.
try:
    import cython
    HAVE_CYTHON = True
    print ("# pysam: cython is available - using cythonize if necessary")
    source_pattern = "pysam/c%s.pyx"
    if HTSLIB_MODE != "external":
        HTSLIB_MODE = "shared"
except ImportError:
    HAVE_CYTHON = False
    print ("# pysam: no cython available - using pre-compiled C")
    # no Cython available - use existing C code
    source_pattern = "pysam/c%s.c"
    if HTSLIB_MODE != "external":
        HTSLIB_MODE = "shared"

# collect pysam version
sys.path.insert(0, "pysam")
import version
version = version.__version__

# exclude sources that contain a main function
EXCLUDE = {
    "samtools": (
        "razip.c", "bgzip.c", "main.c",
        "calDepth.c", "bam2bed.c", "wgsim.c",
        "md5fa.c", "md5sum-lite.c", "maq2sam.c",
        "bamcheck.c", "chk_indel.c", "vcf-miniview.c",
        "htslib-1.3",   # do not import twice
        "hfile_irods.c",  # requires irods library
    ),
    "bcftools": (
        "test", "plugins", "peakfit.c",
        "peakfit.h",
        # needs to renamed, name conflict with samtools reheader
        "reheader.c",
        "polysomy.c"),
    "htslib": (
        'htslib/tabix.c', 'htslib/bgzip.c',
        'htslib/htsfile.c', 'htslib/hfile_irods.c'),
}

print ("# pysam: htslib mode is {}".format(HTSLIB_MODE))
print ("# pysam: HTSLIB_CONFIGURE_OPTIONS={}".format(
    HTSLIB_CONFIGURE_OPTIONS))
htslib_configure_options = None

if HTSLIB_MODE in ['shared', 'separate']:
    package_list += ['pysam.include.htslib',
                     'pysam.include.htslib.htslib']
    package_dirs.update({'pysam.include.htslib':'htslib'})

    htslib_configure_options = configure_library(
        "htslib",
        HTSLIB_CONFIGURE_OPTIONS,
        ["--enable-libcurl",
         "--disable-libcurl"])

    HTSLIB_SOURCE = "builtin"
    print ("# pysam: htslib configure options: {}".format(
        str(htslib_configure_options)))

    config_headers += ["htslib/config.h"]
    if htslib_configure_options is None:
        # create empty config.h file
        with open("htslib/config.h", "w") as outf:
            outf.write(
                "/* empty config.h created by pysam */\n")
            outf.write(
                "/* conservative compilation options */\n")

    with changedir("htslib"):
        htslib_make_options = run_make_print_config()

    for key, value in htslib_make_options.items():
        print ("# pysam: htslib_config {}={}".format(key, value))

    external_htslib_libraries = ['z']
    if "LIBS" in htslib_make_options:
        external_htslib_libraries.extend(
            [re.sub("^-l", "", x) for x in htslib_make_options["LIBS"].split(" ") if x.strip()])

    shared_htslib_sources = [re.sub("\.o", ".c", os.path.join("htslib", x))
                             for x in
                             htslib_make_options["LIBHTS_OBJS"].split(" ")]

    htslib_sources = []

if HTSLIB_LIBRARY_DIR:
    # linking against a shared, externally installed htslib version, no
    # sources required for htslib
    htslib_sources = []
    shared_htslib_sources = []
    chtslib_sources = []
    htslib_library_dirs = [HTSLIB_LIBRARY_DIR]
    htslib_include_dirs = [HTSLIB_INCLUDE_DIR]
    internal_htslib_libraries = []
    external_htslib_libraries = ['z', 'hts']

elif HTSLIB_MODE == 'separate':
    # add to each pysam component a separately compiled
    # htslib
    htslib_sources = shared_htslib_sources
    shared_htslib_sources = htslib_sources
    htslib_library_dirs = []
    htslib_include_dirs = ['htslib']
    internal_htslib_libraries = []

elif HTSLIB_MODE == 'shared':
    # link each pysam component against the same
    # htslib built from sources included in the pysam
    # package.
    htslib_library_dirs = [
        'pysam',
        ".",
        os.path.join("build",
                     distutils_dir_name("lib"),
                     "pysam")]

    htslib_include_dirs = ['htslib']

    if IS_PYTHON3:
        if sys.version_info.minor >= 5:
            internal_htslib_libraries = ["chtslib.{}".format(
                sysconfig.get_config_var('SOABI'))]
        else:
            if sys.platform == "darwin":
                # On OSX, python 3.3 and 3.4 Libs have no platform tags.
                internal_htslib_libraries = ["chtslib"]
            else:
                internal_htslib_libraries = ["chtslib.{}{}".format(
                    sys.implementation.cache_tag,
                    sys.abiflags)]
    else:
        internal_htslib_libraries = ["chtslib"]

else:
    raise ValueError("unknown HTSLIB value '%s'" % HTSLIB_MODE)

# build config.py
with open(os.path.join("pysam", "config.py"), "w") as outf:
    outf.write('HTSLIB = "{}"\n'.format(HTSLIB_SOURCE))
    config_values = collections.defaultdict(int)

    if HTSLIB_SOURCE == "builtin":
        with open(os.path.join("htslib", "config.h")) as inf:
            for line in inf:
                if line.startswith("#define"):
                    key, value = re.match(
                        "#define (\S+)\s+(\S+)", line).groups()
                    config_values[key] = int(value)
            for key in ["ENABLE_PLUGINS",
                        "HAVE_COMMONCRYPTO",
                        "HAVE_GMTIME_R",
                        "HAVE_HMAC",
                        "HAVE_IRODS",
                        "HAVE_LIBCURL",
                        "HAVE_MMAP"]:
                outf.write("{} = {}\n".format(key, config_values[key]))
                print ("# pysam: config_option: {}={}".format(key, config_values[key]))

# create empty config.h files if they have not been created automatically
# or created by the user:
for fn in config_headers:
    if not os.path.exists(fn):
        with open(fn, "w") as outf:
            outf.write(
                "/* empty config.h created by pysam */\n")
            outf.write(
                "/* conservative compilation options */\n")

parts = ["samtools",
         "bcftools",
         "htslib",
         "tabix",
         "faidx",
         "samfile",
         "utils",
         "alignmentfile",
         "tabixproxies",
         "vcf",
         "bcf"]

# Exit if there are no pre-compiled files and no cython available
fn = source_pattern % "htslib"
if not os.path.exists(fn):
    raise ValueError(
        "no cython installed, but can not find {}."
        "Make sure that cython is installed when building "
        "from the repository"
        .format(fn))


#######################################################
classifiers = """
Development Status :: 3 - Beta
Operating System :: MacOS :: MacOS X
Operating System :: POSIX
Operating System :: POSIX :: Linux
Operating System :: Unix
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bioinformatics
"""

#######################################################

#######################################################
# Windows compatibility - untested
if platform.system() == 'Windows':
    include_os = ['win32']
    os_c_files = ['win32/getopt.c']
    extra_compile_args = []
else:
    include_os = []
    os_c_files = []
    # for python 3.4, see for example
    # http://stackoverflow.com/questions/25587039/
    # error-compiling-rpy2-on-python3-4-due-to-werror-
    # declaration-after-statement
    extra_compile_args = [
        "-Wno-unused",
        "-Wno-strict-prototypes",
        "-Wno-sign-compare",
        "-Wno-error=declaration-after-statement"]

define_macros = []

chtslib = Extension(
    "pysam.libchtslib",
    [source_pattern % "htslib",
     "pysam/htslib_util.c"] +
    shared_htslib_sources +
    os_c_files,
    library_dirs=htslib_library_dirs,
    runtime_library_dirs=htslib_library_dirs,
    include_dirs=["pysam", "."] + include_os + htslib_include_dirs,
    libraries=external_htslib_libraries,
    language="c",
    extra_compile_args=extra_compile_args,
    define_macros=define_macros
)

# samfile requires functions defined in bam_md.c
# for __advance_samtools method.
# Selected ones have been copied into samfile_utils.c
# Needs to be devolved somehow.
csamfile = Extension(
    "pysam.csamfile",
    [source_pattern % "samfile",
     "pysam/htslib_util.c",
     "pysam/samfile_util.c",
     "samtools/kprobaln.c"] +
    htslib_sources +
    os_c_files,
    library_dirs=htslib_library_dirs,
    include_dirs=["pysam", "samtools", "."] + include_os + htslib_include_dirs,
    libraries=external_htslib_libraries + internal_htslib_libraries,
    language="c",
    extra_compile_args=extra_compile_args,
    define_macros=define_macros
)

# alignmentfile requires functions defined in bam_md.c
# for __advance_samtools method.
# Selected ones have been copied into samfile_utils.c
# Needs to be devolved somehow.
calignmentfile = Extension(
    "pysam.calignmentfile",
    [source_pattern % "alignmentfile",
     "pysam/htslib_util.c",
     "pysam/samfile_util.c",
     "samtools/kprobaln.c"] +
    htslib_sources +
    os_c_files,
    library_dirs=htslib_library_dirs,
    include_dirs=["pysam", "samtools"] + include_os + htslib_include_dirs,
    libraries=external_htslib_libraries + internal_htslib_libraries,
    language="c",
    extra_compile_args=extra_compile_args,
    define_macros=define_macros
)

# alignmentfile requires functions defined in bam_md.c
# for __advance_samtools method.
# Selected ones have been copied into samfile_utils.c
# Needs to be devolved somehow.
calignedsegment = Extension(
    "pysam.calignedsegment",
    [source_pattern % "alignedsegment",
     "pysam/htslib_util.c",
     "pysam/samfile_util.c",
     "samtools/kprobaln.c"] +
    htslib_sources +
    os_c_files,
    library_dirs=htslib_library_dirs,
    include_dirs=["pysam", "samtools", "."] + include_os + htslib_include_dirs,
    libraries=external_htslib_libraries + internal_htslib_libraries,
    language="c",
    extra_compile_args=extra_compile_args,
    define_macros=define_macros
)

ctabix = Extension(
    "pysam.ctabix",
    [source_pattern % "tabix",
     "pysam/tabix_util.c"] +
    htslib_sources +
    os_c_files,
    library_dirs=["pysam"] + htslib_library_dirs,
    include_dirs=["pysam", "."] + include_os + htslib_include_dirs,
    libraries=external_htslib_libraries + internal_htslib_libraries,
    language="c",
    extra_compile_args=extra_compile_args,
    define_macros=define_macros
)

cutils = Extension(
    "pysam.cutils",
    [source_pattern % "utils", "pysam/pysam_util.c"] +
    glob.glob(os.path.join("samtools", "*.pysam.c")) +
    # glob.glob(os.path.join("samtools", "*", "*.pysam.c")) +
    glob.glob(os.path.join("bcftools", "*.pysam.c")) +
    # glob.glob(os.path.join("bcftools", "*", "*.pysam.c")) +
    htslib_sources +
    os_c_files,
    library_dirs=["pysam"] + htslib_library_dirs,
    include_dirs=["samtools", "bcftools", "pysam", "."] +
    include_os + htslib_include_dirs,
    libraries=external_htslib_libraries + internal_htslib_libraries,
    language="c",
    extra_compile_args=extra_compile_args,
    define_macros=define_macros
)

cfaidx = Extension(
    "pysam.cfaidx",
    [source_pattern % "faidx"] +
    htslib_sources +
    os_c_files,
    library_dirs=["pysam"] + htslib_library_dirs,
    include_dirs=["pysam", "."] + include_os + htslib_include_dirs,
    libraries=external_htslib_libraries + internal_htslib_libraries,
    language="c",
    extra_compile_args=extra_compile_args,
    define_macros=define_macros
)

ctabixproxies = Extension(
    "pysam.ctabixproxies",
    [source_pattern % "tabixproxies"] +
    os_c_files,
    library_dirs=htslib_library_dirs,
    include_dirs=include_os,
    libraries=external_htslib_libraries + internal_htslib_libraries,
    language="c",
    extra_compile_args=extra_compile_args,
    define_macros=define_macros
)

cvcf = Extension(
    "pysam.cvcf",
    [source_pattern % "vcf"] +
    os_c_files,
    library_dirs=htslib_library_dirs,
    include_dirs=["htslib", "."] + include_os + htslib_include_dirs,
    libraries=external_htslib_libraries + internal_htslib_libraries,
    language="c",
    extra_compile_args=extra_compile_args,
    define_macros=define_macros
)

cbcf = Extension(
    "pysam.cbcf",
    [source_pattern % "bcf"] +
    htslib_sources +
    os_c_files,
    library_dirs=htslib_library_dirs,
    include_dirs=["htslib", "."] + include_os + htslib_include_dirs,
    libraries=external_htslib_libraries + internal_htslib_libraries,
    language="c",
    extra_compile_args=extra_compile_args,
    define_macros=define_macros
)

metadata = {
    'name': "pysam",
    'version': version,
    'description': "pysam",
    'long_description': __doc__,
    'author': "Andreas Heger",
    'author_email': "andreas.heger@gmail.com",
    'license': "MIT",
    'platforms': "ALL",
    'url': "https://github.com/pysam-developers/pysam",
    'packages': package_list,
    'requires': ['cython (>=0.21)'],
    'ext_modules': [chtslib,
                    csamfile,
                    calignmentfile,
                    calignedsegment,
                    ctabix,
                    ctabixproxies,
                    cvcf,
                    cbcf,
                    cfaidx,
                    cutils],
    'cmdclass': cmdclass,
    'package_dir': package_dirs,
    'package_data': {'': ['*.pxd', '*.h'], },
    # do not pack in order to permit linking to csamtools.so
    'zip_safe': False,
    'use_2to3': True,
}

if __name__ == '__main__':
    dist = setup(**metadata)
