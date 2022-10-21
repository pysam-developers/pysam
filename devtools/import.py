#################################################################
# Importing samtools, bcftools, and htslib
#
# For each package PKG:
#
# rm -rf PKG
# python import.py PKG path/to/download/PKG-X.Y

import fnmatch
import os
import re
import itertools
import shutil
import sys
import hashlib


EXCLUDE = {
    "samtools": (
        "test", "misc",
        "bgzip.c",
        "main.c",
        "calDepth.c",
        "bam2bed.c",
        "bam_tview.c",
        "bam_tview.h",
        "bam_tview_html.c",
        "bam_tview_curses.c",
        "bam2bcf.c",
        "bam2bcf.h",
        "vcf-miniview.c",
    ),
    "bcftools": (
        "test", "plugins", "peakfit.c",
        "peakfit.h",
        'regidx.c',  # duplicated from htslib
        "polysomy.c"),
    "htslib": (
        'htslib/tabix.c',
        'htslib/bgzip.c',
        'htslib/htsfile.c',
        "test", "tests"),
}


MAIN = {
    "samtools": "bamtk",
    "bcftools": "main"
}


C_VERSION = {
    "htslib":   "HTS_VERSION_TEXT",
    "samtools": "SAMTOOLS_VERSION",
    "bcftools": "BCFTOOLS_VERSION"
}


def locate(pattern, root=os.curdir, exclude=[], exclude_htslib=False):
    '''Locate all files matching supplied filename pattern (but not listed
    in exclude) in and below the supplied root directory. Omit any files under
    directories listed in exclude or (if exclude_htslib=True) matching /htslib-*/.
    '''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            if filename not in exclude:
                yield os.path.join(path, filename)

        for dirname in exclude:
            if dirname in dirs: dirs.remove(dirname)
        if exclude_htslib:
            for dirname in [d for d in dirs if re.match(r"htslib-", d)]:
                dirs.remove(dirname)


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
                    lines = re.sub(r"int main\(", "int {}_main(".format(
                        basename), lines)
                else:
                    lines = re.sub(r"int main\(", "int {}_{}_main(".format(
                        basename, subname), lines)
                lines = re.sub(r"\b({}_stdout)\b".format(basename), r"\1_internal", lines)
                lines = re.sub(r"\bexit\(", "{}_exit(".format(basename), lines)
                lines = re.sub(r"\bstderr\b", "{}_stderr".format(basename), lines)
                lines = re.sub(r"\bstdout\b", "{}_stdout".format(basename), lines)
                lines = re.sub(r" printf\(", " fprintf({}_stdout, ".format(basename), lines)
                lines = re.sub(r"([^kf])puts\(", r"\1{}_puts(".format(basename), lines)
                lines = re.sub(r"putchar\(([^)]+)\)",
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
                if fn == "bamtk.c":
                    lines = re.sub(r'(#include "version.h")', r'\1\n#include "samtools_config_vars.h"', lines)
                    lines = re.sub(r'(else if.*"tview")', r'//\1', lines)

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

    cfiles = locate("*.c", srcdir, exclude=exclude, exclude_htslib=True)
    hfiles = locate("*.h", srcdir, exclude=exclude, exclude_htslib=True)
    mfiles = itertools.chain(locate("README", srcdir), locate("LICENSE", srcdir),
                             locate("version.sh", srcdir, exclude_htslib=True))

    if dest == "htslib":
        # Add build files, including *.ac *.in *.mk *.m4
        mfiles = itertools.chain(mfiles, locate("Makefile", srcdir),
                                 locate("configure", srcdir),
                                 locate("*.[aim][cnk4]", srcdir))

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

    if dest in MAIN:
        # redirect stderr to pysamerr and replace bam.h with a stub.
        sys.stdout.write("applying stderr redirection\n")

        _update_pysam_files(cf, destdir)

    def _getVersion(srcdir):
        with open(os.path.join(srcdir, "version.sh"), encoding="utf-8") as inf:
            for line in inf:
                m = re.match(r"VERSION=(\S+)", line)
                if m: return m.group(1)
            raise ValueError("no VERSION line in version.sh")

    def _update_version_file(key, value, filename):
        tmpfilename = filename + ".tmp"
        with open(filename, encoding="utf-8") as inf:
            with open(tmpfilename, "w", encoding="utf-8") as outf:
                for line in inf:
                    if key in line:
                        line = re.sub(r'"[^"]*"', '"{}"'.format(value), line)
                    outf.write(line)
        os.rename(tmpfilename, filename)

    def _update_version_doc_file(dest, value, filename):
        tmpfilename = filename + ".tmp"
        with open(filename, encoding="utf-8") as inf:
            with open(tmpfilename, "w", encoding="utf-8") as outf:
                for line in inf:
                    if " wraps " in line:
                        # hide the sentence's fullstop from the main regexp
                        line = re.sub(r'\.$', ',DOT', line)
                        line = re.sub(r'{}-[^*,]*'.format(dest),
                                      '{}-{}'.format(dest, value), line)
                        line = re.sub(',DOT', '.', line)
                    outf.write(line)
        os.rename(tmpfilename, filename)

    version = _getVersion(srcdir)
    _update_version_file("__{}_version__".format(dest), version, "pysam/version.py")
    _update_version_file(C_VERSION[dest], version + " (pysam)", "pysam/version.h")
    _update_version_doc_file(dest, version, "README.rst")
    _update_version_doc_file(dest, version, "doc/index.rst")

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

