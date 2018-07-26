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
# git checkout -- samtools/version.h
#
# Manually, then:
# modify config.h to set compatibility flags
#
# For bcftools, type:
# rm -rf bcftools
# python import.py bcftools download/bedtools
# git checkout -- bcftools/version.h
# rm -rf bedtools/test bedtools/plugins

import fnmatch
import os
import re
import itertools
import shutil
import sys
import hashlib


EXCLUDE = {
    "samtools": (
        "razip.c",
        "bgzip.c",
        "main.c",
        "calDepth.c",
        "bam2bed.c",
        "wgsim.c",
        "bam_tview.c",
        "bam_tview.h",
        "bam_tview_html.c",
        "bam_tview_curses.c",
        "md5fa.c",
        "md5sum-lite.c",
        "maq2sam.c",
        "bamcheck.c",
        "chk_indel.c",
        "vcf-miniview.c",
        "hfile_irods.c",  # requires irods library
    ),
    "bcftools": (
        "test", "plugins", "peakfit.c",
        "peakfit.h",
        # needs to renamed, name conflict with samtools reheader
        # "reheader.c",
        "polysomy.c"),
    "htslib": (
        'htslib/tabix.c', 'htslib/bgzip.c',
        'htslib/htsfile.c', 'htslib/hfile_irods.c'),
}


MAIN = {
    "samtools": "bamtk",
    "bcftools": "main"
}



def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
    supplied root directory.
    '''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)


def _update_pysam_files(cf, destdir):
    '''update pysam files applying redirection of ouput'''
    basename = os.path.basename(destdir)
    for filename in cf:
        if not filename:
            continue
        dest = filename + ".pysam.c"
        with open(filename, encoding="utf-8") as infile:
            lines = "".join(infile.readlines())

            with open(dest, "w", encoding="utf-8") as outfile:
                outfile.write('#include "{}.pysam.h"\n\n'.format(basename))
                subname, _ = os.path.splitext(os.path.basename(filename))
                if subname in MAIN.get(basename, []):
                    lines = re.sub("int main\(", "int {}_main(".format(
                        basename), lines)
                else:
                    lines = re.sub("int main\(", "int {}_{}_main(".format(
                        basename, subname), lines)
                lines = re.sub("stderr", "{}_stderr".format(basename), lines)
                lines = re.sub("stdout", "{}_stdout".format(basename), lines)
                lines = re.sub(" printf\(", " fprintf({}_stdout, ".format(basename), lines)
                lines = re.sub("([^kf])puts\(", r"\1{}_puts(".format(basename), lines)
                lines = re.sub("putchar\(([^)]+)\)",
                               r"fputc(\1, {}_stdout)".format(basename), lines)

                fn = os.path.basename(filename)
                # some specific fixes:
                SPECIFIC_SUBSTITUTIONS = {
                    "bam_md.c": (
                        'sam_open_format("-", mode_w',
                        'sam_open_format({}_stdout_fn, mode_w'.format(basename)),
                    "phase.c": (
                        'putc("ACGT"[f->seq[j] == 1? (c&3, {}_stdout) : (c>>16&3)]);'.format(basename),
                        'putc("ACGT"[f->seq[j] == 1? (c&3) : (c>>16&3)], {}_stdout);'.format(basename)),
                    "cut_target.c": (
                        'putc(33 + (cns[j]>>8>>2, {}_stdout));'.format(basename),
                        'putc(33 + (cns[j]>>8>>2), {}_stdout);'.format(basename))
                    }
                if fn in SPECIFIC_SUBSTITUTIONS:
                    lines = lines.replace(
                        SPECIFIC_SUBSTITUTIONS[fn][0],
                        SPECIFIC_SUBSTITUTIONS[fn][1])
                outfile.write(lines)

    with open(os.path.join("import", "pysam.h")) as inf, \
         open(os.path.join(destdir, "{}.pysam.h".format(basename)), "w") as outf:
        outf.write(re.sub("@pysam@", basename, inf.read()))

    with open(os.path.join("import", "pysam.c")) as inf, \
         open(os.path.join(destdir, "{}.pysam.c".format(basename)), "w") as outf:
        outf.write(re.sub("@pysam@", basename, inf.read()))


if len(sys.argv) >= 1:
    if len(sys.argv) != 3:
        raise ValueError("import requires dest src")

    dest, srcdir = sys.argv[1:3]
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
    mfiles = itertools.chain(locate("README", srcdir), locate("LICENSE", srcdir))
    
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
                "".join(open(old_file, "r", encoding="utf-8").readlines()).encode()).digest()
            md5_new = hashlib.md5(
                "".join(open(src, "r", encoding="utf-8").readlines()).encode()).digest()
            if md5_old != md5_new:
                raise ValueError(
                    "incompatible files for %s and %s" %
                    (old_file, src))

        shutil.copy(src, targetdir)
        return old_file

    for src_file in hfiles:
        _compareAndCopy(src_file, srcdir, destdir, exclude)
        ncopied += 1

    for src_file in mfiles:
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

