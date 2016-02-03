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
# python import.py samtools download/samtools
# Manually, then:
# modify config.h to set compatibility flags
# change bamtk.c.pysam.c/main to bamtk.c.pysam.c/samtools_main
#
# For bcftools, type:
# rm -rf bedtools
# python import.py bedtools download/bedtools
import os
import sys
import fnmatch


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


if len(sys.argv) >= 1:
    if len(sys.argv) != 3:
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


# if len(sys.argv) >= 2 and sys.argv[1] == "refresh":
#     sys.stdout.write("refreshing latest source code from .c to .pysam.c")
#     # redirect stderr to pysamerr and replace bam.h with a stub.
#     sys.stdout.write("applying stderr redirection")
#     for destdir in ('samtools', ):
#         pysamcfiles = locate("*.pysam.c", destdir)
#         for f in pysamcfiles:
#             os.remove(f)
#         cfiles = locate("*.c", destdir)
#         _update_pysam_files(cfiles, destdir)

#     sys.exit(0)

