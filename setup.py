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

The current version wraps htslib-1.7, samtools-1.7 and bcftools-1.6.

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
from cy_build import CyExtension as Extension, cy_build_ext as build_ext
try:
    import cython
    HAVE_CYTHON = True
except ImportError:
    HAVE_CYTHON = False

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
    stdout = subprocess.check_output(["make", "-s", "print-config"])
    if IS_PYTHON3:
        stdout = stdout.decode("ascii")

    make_print_config = {}
    for line in stdout.splitlines():
        if "=" in line:
            row = line.split("=")
            if len(row) == 2:
                make_print_config.update(
                    {row[0].strip(): row[1].strip()})
    return make_print_config


def configure_library(library_dir, env_options=None, options=[]):

    configure_script = os.path.join(library_dir, "configure")

    on_rtd = os.environ.get("READTHEDOCS") == "True"
    # RTD has no bzip2 development libraries installed:
    if on_rtd:
        env_options = "--disable-bz2"

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


def get_pysam_version():
    sys.path.insert(0, "pysam")
    import version
    return version.__version__
    

# How to link against HTSLIB
# shared:   build shared chtslib from builtin htslib code.
# external: use shared libhts.so compiled outside of
#           pysam
# separate: use included htslib and include in each extension
#           module. No dependencies between modules and works with
#           setup.py install, but wasteful in terms of memory and
#           compilation time. Fallback if shared module compilation
#           fails.

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

# list of config files that will be automatically generated should
# they not already exist or be created by configure scripts in the
# subpackages.
config_headers = ["samtools/config.h",
                  "bcftools/config.h"]

cmdclass = {'build_ext': build_ext}

# If cython is available, the pysam will be built using cython from
# the .pyx files. If no cython is available, the C-files included in the
# distribution will be used.
if HAVE_CYTHON:
    print ("# pysam: cython is available - using cythonize if necessary")
    source_pattern = "pysam/libc%s.pyx"
else:
    print ("# pysam: no cython available - using pre-compiled C")
    source_pattern = "pysam/libc%s.c"

# Exit if there are no pre-compiled files and no cython available
fn = source_pattern % "htslib"
if not os.path.exists(fn):
    raise ValueError(
        "no cython installed, but can not find {}."
        "Make sure that cython is installed when building "
        "from the repository"
        .format(fn))

# exclude sources that contain a main function
EXCLUDE = {
    "samtools": (
    ),
    "bcftools": (
        "test", "plugins", "peakfit.c",
        "peakfit.h",
        # needs to renamed, name conflict with samtools reheader
        "reheader.c",
        "polysomy.c"),
    "htslib": (
        'htslib/tabix.c',
        'htslib/bgzip.c',
        'htslib/htsfile.c'),
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
    external_htslib_libraries = ['z', 'hts']
elif HTSLIB_MODE == 'separate':
    # add to each pysam component a separately compiled
    # htslib
    htslib_sources = shared_htslib_sources
    shared_htslib_sources = htslib_sources
    htslib_library_dirs = []
    htslib_include_dirs = ['htslib']
elif HTSLIB_MODE == 'shared':
    # link each pysam component against the same
    # htslib built from sources included in the pysam
    # package.
    htslib_library_dirs = [
        "pysam",  # when using setup.py develop?
        ".",  # when using setup.py develop?
        os.path.join("build", distutils_dir_name("lib"), "pysam")]

    htslib_include_dirs = ['htslib']
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
                    config_values[key] = value
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

suffix = sysconfig.get_config_var('EXT_SUFFIX')
if not suffix:
    suffix = sysconfig.get_config_var('SO')

internal_htslib_libraries = [
    os.path.splitext("chtslib{}".format(suffix))[0]]
internal_samtools_libraries = [
    os.path.splitext("csamtools{}".format(suffix))[0],
    os.path.splitext("cbcftools{}".format(suffix))[0],
    ]
internal_pysamutil_libraries = [
    os.path.splitext("cutils{}".format(suffix))[0]]

libraries_for_pysam_module = external_htslib_libraries + internal_htslib_libraries + internal_pysamutil_libraries

# Order of modules matters in order to make sure that dependencies are resolved.
# The structures of dependencies is as follows:
# libchtslib: htslib utility functions and htslib itself if builtin is set.
# libcsamtools: samtools code (builtin)
# libcbcftools: bcftools code (builtin)
# libcutils: General utility functions, depends on all of the above
# libcXXX (pysam module): depends on libchtslib and libcutils

# The list below uses the union of include_dirs and library_dirs for
# reasons of simplicity.

modules = [
    dict(name="pysam.libchtslib",
         sources=[source_pattern % "htslib", "pysam/htslib_util.c"] + shared_htslib_sources + os_c_files,
         libraries=external_htslib_libraries),
    dict(name="pysam.libcsamtools",
         sources=[source_pattern % "samtools"] + glob.glob(os.path.join("samtools", "*.pysam.c")) +
         [os.path.join("samtools", "lz4", "lz4.c")] + htslib_sources + os_c_files,
         libraries=external_htslib_libraries + internal_htslib_libraries),
    dict(name="pysam.libcbcftools",
         sources=[source_pattern % "bcftools"] + glob.glob(os.path.join("bcftools", "*.pysam.c")) + htslib_sources + os_c_files,
         libraries=external_htslib_libraries + internal_htslib_libraries),
    dict(name="pysam.libcutils",
         sources=[source_pattern % "utils", "pysam/pysam_util.c"] + htslib_sources + os_c_files,
         libraries=external_htslib_libraries + internal_htslib_libraries + internal_samtools_libraries),
    dict(name="pysam.libcalignmentfile",
         sources=[source_pattern % "alignmentfile"] + htslib_sources + os_c_files,
         libraries=libraries_for_pysam_module),
    dict(name="pysam.libcsamfile",
         sources=[source_pattern % "samfile"] + htslib_sources + os_c_files,
         libraries=libraries_for_pysam_module),
    dict(name="pysam.libcalignedsegment",
         sources=[source_pattern % "alignedsegment"] + htslib_sources + os_c_files,
         libraries=libraries_for_pysam_module),
    dict(name="pysam.libctabix",
         sources=[source_pattern % "tabix"] + htslib_sources + os_c_files,
         libraries=libraries_for_pysam_module),
    dict(name="pysam.libcfaidx",
         sources=[source_pattern % "faidx"] + htslib_sources + os_c_files,
         libraries=libraries_for_pysam_module),
    dict(name="pysam.libcbcf",
         sources=[source_pattern % "bcf"] + htslib_sources + os_c_files,
         libraries=libraries_for_pysam_module),
    dict(name="pysam.libcbgzf",
         sources=[source_pattern % "bgzf"] + htslib_sources + os_c_files,
         libraries=libraries_for_pysam_module),
    dict(name="pysam.libctabixproxies",
         sources=[source_pattern % "tabixproxies"] + htslib_sources + os_c_files,
         libraries=libraries_for_pysam_module),
    dict(name="pysam.libcvcf",
         sources=[source_pattern % "vcf"] + htslib_sources + os_c_files,
         libraries=libraries_for_pysam_module),
]

common_options = dict(
    language="c",
    extra_compile_args=extra_compile_args,
    define_macros=define_macros,
    # for out-of-tree compilation, use absolute paths
    library_dirs=[os.path.abspath(x) for x in ["pysam"] + htslib_library_dirs],
    include_dirs=[os.path.abspath(x) for x in htslib_include_dirs + \
                  ["samtools", "samtools/lz4", "bcftools", "pysam", "."] + include_os])

# add common options (in python >3.5, could use n = {**a, **b}
for module in modules:
    module.update(**common_options)

classifiers = """
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""
    
metadata = {
    'name': "pysam",
    'version': get_pysam_version(),
    'description': "pysam",
    'long_description': __doc__,
    'author': "Andreas Heger",
    'author_email': "andreas.heger@gmail.com",
    'license': "MIT",
    'platforms': ["POSIX", "UNIX", "MacOS"],
    'classifiers': [_f for _f in classifiers.split("\n") if _f],
    'url': "https://github.com/pysam-developers/pysam",
    'packages': package_list,
    'requires': ['cython (>=0.21)'],
    'ext_modules': [Extension(**opts) for opts in modules],
    'cmdclass': cmdclass,
    'package_dir': package_dirs,
    'package_data': {'': ['*.pxd', '*.h'], },
    # do not pack in order to permit linking to csamtools.so
    'zip_safe': False,
    'use_2to3': True,
}

if __name__ == '__main__':
    dist = setup(**metadata)
