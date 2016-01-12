#! /usr/bin/python

'''The SAM/BAM/CRAM format is a way to store efficiently large numbers
of alignments, such as those routinely are created by next-generation
sequencing methods.

This module provides a low-level wrapper around the htslib C-API as
using cython and a high-level API for convenient access to the data in
SAM/BAM formatted files. Also included is an interface to the samtools
command line utilities and the tabix C-API for reading compressed and
indexed tabular data.

The current version wraps htslib-1.2.1 and samtools-1.2.

See:
http://www.htslib.org
https://github.com/pysam-developers/pysam
http://pysam.readthedocs.org/en/stable

'''

import os
import sys
import glob
import shutil
import hashlib
import re
import fnmatch
import platform

name = "pysam"

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
samtools_exclude = ("bamtk.c",
                    "razip.c",
                    "bgzip.c",
                    "main.c",
                    "calDepth.c",
                    "bam2bed.c",
                    "wgsim.c",
                    "md5fa.c",
                    "maq2sam.c",
                    "bamcheck.c",
                    "chk_indel.c",
                    "vcf-miniview.c",
                    "htslib-1.2.1",   # do not import twice
                    "hfile_irods.c",  # requires irods library
                    )

htslib_exclude = ('htslib/tabix.c',
                  'htslib/bgzip.c',
                  'htslib/htsfile.c',
                  'htslib/hfile_irods.c')

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
        if x not in htslib_exclude]
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
        if x not in htslib_exclude]
    htslib_library_dirs = ['pysam']
    htslib_include_dirs = ['htslib']
    htslib_libraries = ['chtslib']
else:
    raise ValueError("unknown HTSLIB value '%s'" % HTSLIB_MODE)


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
# python setup.py import download/samtools
#
if len(sys.argv) >= 2 and sys.argv[1] == "import":
    if len(sys.argv) < 3:
        raise ValueError("missing PATH to samtools source directory")

    destdir = samtools_dest
    srcdir = sys.argv[2]
    exclude = samtools_exclude

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

from setuptools import Extension, setup

#######################################################
parts = ["samtools", "htslib", "tabix",
         "faidx", "samfile", "utils",
         "alignmentfile", "tabixproxies",
         "vcf", "bcf"]

try:
    from Cython.Distutils import build_ext
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
Development Status :: 2 - Alpha
Operating System :: MacOS :: MacOS X
Operating System :: Microsoft :: Windows :: Windows NT/2000
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

csamtools = Extension(
    "pysam.csamtools",
    [source_pattern % "samtools",
     "pysam/pysam_util.c"] +
    glob.glob(os.path.join("samtools", "*.pysam.c")) +
    glob.glob(os.path.join("samtools", "*", "*.pysam.c")) +
    os_c_files +
    htslib_sources,
    library_dirs=htslib_library_dirs,
    include_dirs=["samtools", "pysam"] + include_os + htslib_include_dirs,
    libraries=["z"] + htslib_libraries,
    language="c",
    extra_compile_args=extra_compile_args,
    define_macros=define_macros
)

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
    [source_pattern % "utils"] +
    htslib_sources +
    os_c_files,
    library_dirs=["pysam"] + htslib_library_dirs,
    include_dirs=["pysam"] + include_os + htslib_include_dirs,
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
    'name': name,
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
                 # 'pysam.include.samtools.bcftools',
                 'pysam.include.samtools.win32'],
    'requires': ['cython (>=0.21)'],
    'ext_modules': [csamtools,
                    chtslib,
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
                    'pysam.include.samtools': 'samtools'},
    'package_data': {'': ['*.pxd', '*.h'], },
    # do not pack in order to permit linking to csamtools.so
    'zip_safe': False,
    'use_2to3': True,
}

if __name__ == '__main__':
    dist = setup(**metadata)
