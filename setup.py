#! /usr/bin/python

'''The SAM/BAM/CRAM format is a way to store efficiently large numbers
of alignments, such as those routinely are created by next-generation
sequencing methods.

This module provides a low-level wrapper around the htslib C-API as
using cython and a high-level API for convenient access to the data in
SAM/BAM formatted files. Also included is an interface to the samtools
command line utilities and the tabix C-API for reading compressed and
indexed tabular data.

The current version wraps htslib-1.3 and samtools-1.3.

See:
http://www.htslib.org
https://github.com/pysam-developers/pysam
http://pysam.readthedocs.org/en/stable

'''

import collections
import fnmatch
import glob
import hashlib
import os
import platform
import re
import shutil
import subprocess
import sys
from contextlib import contextmanager
from setuptools import Extension, setup


@contextmanager
def changedir(path):
    save_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(save_dir)


def configure_library(library_dir, env_options=None, options=[]):

    configure_script = os.path.join(library_dir, "configure")

    if not os.path.exists(configure_script):
        raise ValueError(
            "configure script {} does not exist".format(configure_script))

    def run_configure(option):
        try:
            retcode = subprocess.call(
                " ".join(("./configure", option)),
                shell=True)
            if retcode != 0:
                return False
            else:
                print "# successful configure run with options {}".format(
                    option)
                return True
        except OSError as e:
            return False

    with changedir(library_dir):
        if env_options is not None:
            if run_configure(env_options):
                return True

        for option in options:
            if run_configure(option):
                break


def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
    supplied root directory.
    '''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)


def _update_pysam_files(cf, destdir):
    '''update pysam files applying redirection of ouput'''
    for filename in cf:
        if not filename:
            continue
        dest = filename + ".pysam.c"
        with open(filename) as infile:
            with open(dest, "w") as outfile:
                outfile.write('#include "pysam.h"\n\n')
                outfile.write(
                    re.sub("stderr", "pysamerr", "".join(infile.readlines())))
            with open(os.path.join(destdir, "pysam.h"), "w")as outfile:
                outfile.write("""#ifndef PYSAM_H
#define PYSAM_H
#include "stdio.h"
extern FILE * pysamerr;
#endif
""")

IS_PYTHON3 = sys.version_info[0] >= 3

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
HTSLIB_MODE = "separate"
HTSLIB_LIBRARY_DIR = os.environ.get('HTSLIB_LIBRARY_DIR', None)
HTSLIB_INCLUDE_DIR = os.environ.get('HTSLIB_INCLUDE_DIR', None)

# collect pysam version
sys.path.insert(0, "pysam")
import version
version = version.__version__

# exclude sources that contains a main function
EXCLUDE = {
    "samtools": ("razip.c",
                 "bgzip.c",
                 "main.c",
                 "calDepth.c",
                 "bam2bed.c",
                 "wgsim.c",
                 "md5fa.c",
                 "md5sum-lite.c",
                 "maq2sam.c",
                 "bamcheck.c",
                 "chk_indel.c",
                 "vcf-miniview.c",
                 "htslib-1.3",   # do not import twice
                 "hfile_irods.c",  # requires irods library
             ),
    "bcftools": ("test",
                 "plugins",
                 "peakfit.c",
                 "peakfit.h",
                 # needs to renamed, name conflict with samtools reheader
                 "reheader.c",
                 "polysomy.c"),
    "htslib": ('htslib/tabix.c',
               'htslib/bgzip.c',
               'htslib/htsfile.c',
               'htslib/hfile_irods.c'),
    }

# destination directories for import of samtools and tabix
samtools_dest = os.path.abspath("samtools")

if HTSLIB_LIBRARY_DIR:
    # linking against a shared, externally installed htslib version, no
    # sources required for htslib
    htslib_sources = []
    shared_htslib_sources = []
    chtslib_sources = []
    htslib_library_dirs = [HTSLIB_LIBRARY_DIR]
    htslib_include_dirs = [HTSLIB_INCLUDE_DIR]
    htslib_libraries = ['hts']

elif HTSLIB_MODE == 'separate':
    # add to each pysam component a separately compiled
    # htslib
    htslib_sources = [
        x for x in
        glob.glob(os.path.join("htslib", "*.c")) +
        glob.glob(os.path.join("htslib", "cram", "*.c"))
        if x not in EXCLUDE["htslib"]]
    shared_htslib_sources = htslib_sources
    htslib_library_dirs = []
    htslib_include_dirs = ['htslib']
    htslib_libraries = []

elif HTSLIB_MODE == 'shared':
    # link each pysam component against the same
    # htslib built from sources included in the pysam
    # package.
    htslib_sources = []
    shared_htslib_sources = [
        x for x in
        glob.glob(os.path.join("htslib", "*.c")) +
        glob.glob(os.path.join("htslib", "cram", "*.c"))
        if x not in EXCLUDE["htslib"]]
    htslib_library_dirs = ['pysam']
    htslib_include_dirs = ['htslib']
    htslib_libraries = ['chtslib']
else:
    raise ValueError("unknown HTSLIB value '%s'" % HTSLIB_MODE)


print "# htslib mode is {}".format(HTSLIB_MODE)


if HTSLIB_MODE in ['shared', 'separate']:

    configure_library(
        "htslib",
        os.environ.get('HTSLIB_COMPILE_OPTIONS', None),
        ["--enable-libcurl --enable-plugins",
         "--enable-plugins",
         ""])

    HTSLIB_MODE = "builtin"

# build config.py
with open(os.path.join("pysam", "config.py"), "w") as outf:
    outf.write('HTSLIB_MODE = "{}"\n'.format(HTSLIB_MODE))
    config_values = collections.defaultdict(int)

    if HTSLIB_MODE == "builtin":
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



#################################################################
# Importing samtools and htslib
#
# For htslib, simply copy the whole release tar-ball
# into the directory "htslib" and recreate the file version.h
#
# rm -rf htslib
# mv download/htslib htslib
# git checkout -- htslib/version.h
# Edit the file htslib/version.h to set the right version number.
#
# For samtools, type:
# rm -rf samtools
# python setup.py import samtools download/samtools
# Manually, then:
# modify config.h to set compatibility flags
# change bamtk.c.pysam.c/main to bamtk.c.pysam.c/samtools_main
#
# For bcftools, type:
# rm -rf bedtools
# python setup.py import bedtools download/bedtools

if len(sys.argv) >= 2 and sys.argv[1] == "import":
    if len(sys.argv) != 4:
        raise ValueError("import requires dest src")

    dest, srcdir = sys.argv[2:4]
    if dest not in EXCLUDE:
        raise ValueError("import expected one of %s" %
                         ",".join(EXCLUDE.keys()))
    exclude = EXCLUDE[dest]
    destdir = os.path.abspath(dest)
    srcdir = os.path.abspath(srcdir)
    if not os.path.exists(srcdir):
        raise IOError(
            "source directory `%s` does not exist." % srcdir)

    cfiles = locate("*.c", srcdir)
    hfiles = locate("*.h", srcdir)

    # remove unwanted files and htslib subdirectory.
    cfiles = [x for x in cfiles if os.path.basename(x) not in exclude
              and not re.search("htslib-", x)]

    hfiles = [x for x in hfiles if os.path.basename(x) not in exclude
              and not re.search("htslib-", x)]

    ncopied = 0

    def _compareAndCopy(src, srcdir, destdir, exclude):

        d, f = os.path.split(src)
        common_prefix = os.path.commonprefix((d, srcdir))
        subdir = re.sub(common_prefix, "", d)[1:]
        targetdir = os.path.join(destdir, subdir)
        if not os.path.exists(targetdir):
            os.makedirs(targetdir)
        old_file = os.path.join(targetdir, f)
        if os.path.exists(old_file):
            md5_old = hashlib.md5(
                "".join(open(old_file, "r").readlines())).digest()
            md5_new = hashlib.md5(
                "".join(open(src, "r").readlines())).digest()
            if md5_old != md5_new:
                raise ValueError(
                    "incompatible files for %s and %s" %
                    (old_file, src))

        shutil.copy(src, targetdir)
        return old_file

    for src_file in hfiles:
        _compareAndCopy(src_file, srcdir, destdir, exclude)
        ncopied += 1

    cf = []
    for src_file in cfiles:
        cf.append(_compareAndCopy(src_file,
                                  srcdir,
                                  destdir,
                                  exclude))
        ncopied += 1

    sys.stdout.write(
        "installed latest source code from %s: "
        "%i files copied\n" % (srcdir, ncopied))
    # redirect stderr to pysamerr and replace bam.h with a stub.
    sys.stdout.write("applying stderr redirection\n")

    _update_pysam_files(cf, destdir)

    sys.exit(0)


if len(sys.argv) >= 2 and sys.argv[1] == "refresh":
    sys.stdout.write("refreshing latest source code from .c to .pysam.c")
    # redirect stderr to pysamerr and replace bam.h with a stub.
    sys.stdout.write("applying stderr redirection")
    for destdir in ('samtools', ):
        pysamcfiles = locate("*.pysam.c", destdir)
        for f in pysamcfiles:
            os.remove(f)
        cfiles = locate("*.c", destdir)
        _update_pysam_files(cfiles, destdir)

    sys.exit(0)


###################
# populate headers
# mkdir pysam/include pysam/include/win32
# touch pysam/include/__init__.py pysam/include/win32/__init__.py
# cp samtools/*.h pysam/*.h pysam/include
# cp samtools/win32/*.h pysam/include/win32



#######################################################
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

try:
    from cy_build import CyExtension as Extension, cy_build_ext as build_ext
except ImportError:
    # no Cython available - use existing C code
    cmdclass = {}
    source_pattern = "pysam/c%s.c"
else:
    # remove existing files to recompute
    # necessary to be both compatible for python 2.7 and 3.3
    if IS_PYTHON3:
        for part in parts:
            try:
                os.unlink("pysam/c%s.c" % part)
            except:
                pass
    source_pattern = "pysam/c%s.pyx"
    cmdclass = {'build_ext': build_ext}

#######################################################
classifiers = """
Development Status :: 3 - Beta
Operating System :: MacOS :: MacOS X
Operating System :: OS Independent
Operating System :: POSIX
Operating System :: POSIX :: Linux
Operating System :: Unix
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bioinformatics
"""

#######################################################
# Windows compatibility
if platform.system() == 'Windows':
    include_os = ['win32']
    os_c_files = ['win32/getopt.c']
else:
    include_os = []
    os_c_files = []

#######################################################
extra_compile_args = ["-Wno-error=declaration-after-statement",
                      "-DSAMTOOLS=1"]
define_macros = [('_FILE_OFFSET_BITS', '64'),
                 ('_USE_KNETFILE', '')]

chtslib = Extension(
    "pysam.libchtslib",
    [source_pattern % "htslib",
     "pysam/htslib_util.c"] +
    shared_htslib_sources +
    os_c_files,
    library_dirs=htslib_library_dirs,
    include_dirs=["pysam"] + include_os + htslib_include_dirs,
    libraries=["z"] + htslib_libraries,
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
    include_dirs=["pysam", "samtools"] + include_os + htslib_include_dirs,
    libraries=["z"] + htslib_libraries,
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
    libraries=["z"] + htslib_libraries,
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
    include_dirs=["pysam", "samtools"] + include_os + htslib_include_dirs,
    libraries=["z"] + htslib_libraries,
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
    include_dirs=["pysam"] + include_os + htslib_include_dirs,
    libraries=["z"] + htslib_libraries,
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
    include_dirs=["samtools", "bcftools", "pysam"] +
    include_os + htslib_include_dirs,
    libraries=["z"] + htslib_libraries,
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
    include_dirs=["pysam"] + include_os + htslib_include_dirs,
    libraries=["z"] + htslib_libraries,
    language="c",
    extra_compile_args=extra_compile_args,
    define_macros=define_macros
)

ctabixproxies = Extension(
    "pysam.ctabixproxies",
    [source_pattern % "tabixproxies"] + 
    os_c_files,
    library_dirs=[],
    include_dirs=include_os,
    libraries=["z"],
    language="c",
    extra_compile_args=extra_compile_args,
    define_macros=define_macros
)

cvcf = Extension(
    "pysam.cvcf",
    [source_pattern % "vcf"] + 
    os_c_files,
    library_dirs=[],
    include_dirs=["htslib"] + include_os + htslib_include_dirs,
    libraries=["z"],
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
    include_dirs=["htslib"] + include_os + htslib_include_dirs,
    libraries=["z"] + htslib_libraries,
    language="c",
    extra_compile_args=extra_compile_args,
    define_macros=define_macros
)

metadata = {
    'name': "psyam",
    'version': version,
    'description': "pysam",
    'long_description': __doc__,
    'author': "Andreas Heger",
    'author_email': "andreas.heger@gmail.com",
    'license': "MIT",
    'platforms': "ALL",
    'url': "https://github.com/pysam-developers/pysam",
    'packages': ['pysam',
                 'pysam.include',
                 'pysam.include.htslib',
                 'pysam.include.htslib.htslib',
                 'pysam.include.samtools',
                 'pysam.include.bcftools',
                 'pysam.include.samtools.win32'],
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
    'package_dir': {'pysam': 'pysam',
                    'pysam.include.htslib': 'htslib',
                    'pysam.include.samtools': 'samtools',
                    'pysam.include.bcftools': 'bcftools'},
    'package_data': {'': ['*.pxd', '*.h'], },
    # do not pack in order to permit linking to csamtools.so
    'zip_safe': False,
    'use_2to3': True,
}

if __name__ == '__main__':
    dist = setup(**metadata)
