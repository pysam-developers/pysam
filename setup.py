#!/usr/bin/python
'''

pysam
*****

'''

import os
import sys
import glob
import shutil
import hashlib
import re
import fnmatch
import platform
import subprocess

name = "pysam"

IS_PYTHON3 = sys.version_info[0] >= 3

# collect pysam version
sys.path.insert(0, "pysam")
import version

version = version.__version__

samtools_exclude = ("bamtk.c", "razip.c", "bgzip.c",
                    "main.c", "calDepth.c", "bam2bed.c",
                    "wgsim.c", "md5fa.c", "maq2sam.c",
                    "bamcheck.c",
                    "chk_indel.c")
htslib_exclude = ('htslib/tabix.c',)
samtools_dest = os.path.abspath("samtools")
tabix_exclude = ("main.c",)
tabix_dest = os.path.abspath("tabix")


def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
    supplied root directory.'''
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

# copy samtools source
if len(sys.argv) >= 2 and sys.argv[1] == "import":
    if len(sys.argv) < 3:
        raise ValueError("missing PATH to samtools source directory")
    if len(sys.argv) < 4:
        raise ValueError("missing PATH to tabix source directory")

    for destdir, srcdir, exclude in zip(
            (samtools_dest, tabix_dest),
            sys.argv[2:4],
            (samtools_exclude, tabix_exclude)):

        srcdir = os.path.abspath(srcdir)
        if not os.path.exists(srcdir):
            raise IOError("samtools src dir `%s` does not exist." % srcdir)

        cfiles = locate("*.c", srcdir)
        hfiles = locate("*.h", srcdir)
        ncopied = 0

        def _compareAndCopy(src, srcdir, destdir, exclude):

            d, f = os.path.split(src)
            if f in exclude:
                return None
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
                        "incompatible files for %s and %s" % (old_file, src))

            shutil.copy(src, targetdir)
            return old_file

        for src_file in hfiles:
            _compareAndCopy(src_file, srcdir, destdir, exclude)
            ncopied += 1

        cf = []
        for src_file in cfiles:
            cf.append(_compareAndCopy(src_file, srcdir, destdir, exclude))
            ncopied += 1

        sys.stdout.write(
            "installed latest source code from %s: "
            "%i files copied" % (srcdir, ncopied))
        # redirect stderr to pysamerr and replace bam.h with a stub.
        sys.stdout.write("applying stderr redirection")

        _update_pysam_files(cf, destdir)

    sys.exit(0)


if len(sys.argv) >= 2 and sys.argv[1] == "refresh":
    sys.stdout.write("refreshing latest source code from .c to .pysam.c")
    # redirect stderr to pysamerr and replace bam.h with a stub.
    sys.stdout.write("applying stderr redirection")
    for destdir in ('samtools', 'tabix'):
        pysamcfiles = locate("*.pysam.c", destdir)
        for f in pysamcfiles:
            os.remove(f)
        cfiles = locate("*.c", destdir)
        _update_pysam_files(cfiles, destdir)

    sys.exit(0)


# checkout latest version of htslib
if len(sys.argv) == 2 and sys.argv[1] == "htslib":
    if not os.path.exists("htslib"):
        subprocess.call(["git", "clone", "git@github.com:samtools/htslib.git"])
        with open(os.path.join("htslib", "version.h"), "w") as outfile:
            outfile.write('#define HTS_VERSION "0.0.1"')
    else:
        os.chdir('htslib')
        subprocess.call(["git", "pull"])

    sys.exit(0)

###################
# populate headers
# mkdir pysam/include pysam/include/win32
# touch pysam/include/__init__.py pysam/include/win32/__init__.py
# cp samtools/*.h pysam/*.h pysam/include
# cp samtools/win32/*.h pysam/include/win32

try:
    from setuptools import Extension, setup, find_packages
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()
    from setuptools import Extension, setup, find_packages

#######################################################
#######################################################
try:
    from Cython.Distutils import build_ext
except ImportError:
    # no Cython available - use existing C code
    cmdclass = {}
    csamtools_sources = ["pysam/csamtools.c"]
    chtslib_sources = ["pysam/chtslib.c"]
    tabix_sources = ["pysam/ctabix.c"]
    tabproxies_sources = ["pysam/TabProxies.c"]
    cvcf_sources = ["pysam/cvcf.c"]
else:
    # remove existing files to recompute
    # necessary to be both compatible for python 2.7 and 3.3
    if IS_PYTHON3:
        for f in ("pysam/csamtools.c",
                  "pysam/chtslib.c",
                  "pysam/ctabix.c",
                  "pysam/TabProxies.c",
                  "pysam/cvcf.c"):
            try:
                os.unlink(f)
            except:
                pass

    cmdclass = {'build_ext': build_ext}
    csamtools_sources = ["pysam/csamtools.pyx"]
    chtslib_sources = ["pysam/chtslib.pyx"]
    tabix_sources = ["pysam/ctabix.pyx"]
    tabproxies_sources = ["pysam/TabProxies.pyx"]
    cvcf_sources = ["pysam/cvcf.pyx"]

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
samtools = Extension(
    "pysam.csamtools",
    csamtools_sources +
    ["pysam/%s" % x for x in (
        "pysam_util.c", )] +
    glob.glob(os.path.join("samtools", "*.pysam.c")) +
    os_c_files +
    glob.glob(os.path.join("samtools", "*", "*.pysam.c")),
    library_dirs=[],
    include_dirs=["samtools", "pysam"] + include_os,
    libraries=["z", ],
    language="c",
    extra_compile_args=["-Wno-error=declaration-after-statement"],
    define_macros=[('_FILE_OFFSET_BITS', '64'),
                   ('_USE_KNETFILE', '')]
)

#######################################################
htslib = Extension(
    "pysam.chtslib",
    chtslib_sources +
    ["pysam/%s" % x for x in (
        "htslib_util.c", )] +
    [x for x in glob.glob(
        os.path.join("htslib", "*.c")) +
     glob.glob(
         os.path.join("htslib", "cram", "*.c"))
     if x not in htslib_exclude] +
    os_c_files,
    library_dirs=[], # "/home/andreas/devel/htslib"],
    include_dirs=["htslib",
                  "pysam"] + include_os,
    # at later stage, to include libhts.so, add "hts",
    libraries=["z"],
    language="c",
    extra_compile_args=["-Wno-error=declaration-after-statement",
                        "-DSAMTOOLS=1"],
    define_macros=[('_FILE_OFFSET_BITS', '64'),
                   ('_USE_KNETFILE', '')]
)

tabix = Extension(
    "pysam.ctabix",
    tabix_sources +
    ["pysam/%s" % x for x in ("tabix_util.c", )] +
    [x for x in glob.glob(
        os.path.join("htslib", "*.c")) +
     glob.glob(
         os.path.join("htslib", "cram", "*.c"))
     if x not in htslib_exclude] +
    os_c_files,
    # glob.glob(os.path.join("tabix", "*.pysam.c")),
    library_dirs=[],
    include_dirs=["htslib", "pysam"] + include_os,
    libraries=["z", ],
    language="c",
    extra_compile_args=["-Wno-error=declaration-after-statement",
                        "-DSAMTOOLS=1"],
    define_macros=[('_FILE_OFFSET_BITS', '64'),
                   ('_USE_KNETFILE', '')],
)

tabproxies = Extension(
    "pysam.TabProxies",
    tabproxies_sources + os_c_files,
    library_dirs=[],
    include_dirs=include_os,
    libraries=["z", ],
    language="c",
    extra_compile_args=["-Wno-error=declaration-after-statement"],
)

cvcf = Extension(
    "pysam.cvcf",
    cvcf_sources + os_c_files,
    library_dirs=[],
    include_dirs=["tabix", "htslib"] + include_os,
    libraries=["z", ],
    language="c",
    extra_compile_args=["-Wno-error=declaration-after-statement"],
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
    'url': "http://code.google.com/p/pysam/",
    'packages': ['pysam',
                 'pysam.include',
                 'pysam.include.htslib',
                 'pysam.include.samtools',
                 'pysam.include.samtools.bcftools',
                 'pysam.include.samtools.win32'],
                 #'pysam.include.tabix'],
    'requires': ['cython (>=0.17)'],
    'ext_modules': [samtools, htslib, tabix, tabproxies, cvcf],
    'cmdclass': cmdclass,
    'install_requires': ['cython>=0.17', ],
    'package_dir': {'pysam': 'pysam',
                    'pysam.include.htslib': 'htslib',
                    'pysam.include.samtools': 'samtools'},
                    #'pysam.include.tabix': 'tabix'},
    'package_data': {'': ['*.pxd', '*.h'], },
    # do not pack in order to permit linking to csamtools.so
    'zip_safe': False,
    'use_2to3': True,
}

if __name__ == '__main__':
    dist = setup(**metadata)
