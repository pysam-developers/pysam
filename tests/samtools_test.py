#!/usr/bin/env python
'''unit testing code for pysam.

Execute in the :file:`tests` directory as it requires the Makefile
and data files located there.
'''

import pysam
import pysam.samtools
import pysam.bcftools
import unittest
import os
import re
import glob
import sys
import subprocess
import shutil
from TestUtils import checkBinaryEqual, check_lines_equal, \
    check_samtools_view_equal, get_temp_filename, force_bytes

IS_PYTHON3 = sys.version_info[0] >= 3

WORKDIR = "pysam_test_work"
DATADIR = "pysam_data"


def run_command(cmd):
    '''run a samtools command'''
    try:
        retcode = subprocess.call(cmd, shell=True,
                                  stderr=subprocess.PIPE)
        if retcode < 0:
            print("Child was terminated by signal", -retcode)
    except OSError as e:
        print("Execution failed:", e)


def get_version(executable):
    '''return samtools/bcftools version'''

    with subprocess.Popen(executable, shell=True,
                          stderr=subprocess.PIPE).stderr as pipe:
        lines = b"".join(pipe.readlines())

    if IS_PYTHON3:
        lines = lines.decode('ascii')
    try:
        x = re.search("Version:\s+(\S+)", lines).groups()[0]
    except AttributeError:
        raise ValueError("could not get version from %s" % lines)
    return x


class SamtoolsTest(unittest.TestCase):

    '''test samtools command line commands and compare
    against pysam commands.

    Tests fail, if the output is not binary identical.
    '''

    requisites = [
        "ex1.fa", "ex1.fa.fai",
        "ex1.sam.gz",
        "ex1.bam", "ex1.bam.bai",
        "ex1.sam", "ex2.bam",
        "ex1.bed"]

    # a list of statements to test
    # should contain at least one %(out)s component indicating
    # an output file.
    statements = [
        "view ex1.bam > %(out)s_ex1.view",
        # ("view -bT ex1.fa -o %(out)s_ex1.view2 ex1.sam",
        "sort ex1.bam -o %(out)s_ex1.sort.bam",
        "mpileup ex1.bam > %(out)s_ex1.pileup",
        "depth ex1.bam > %(out)s_ex1.depth",
        # TODO: issues with file naming
        # "faidx ex1.fa; %(out)s_ex1.fa.fai",
        "index ex1.bam %(out)s_ex1.bam.fai",
        "idxstats ex1.bam > %(out)s_ex1.idxstats",
        "fixmate ex1.bam %(out)s_ex1.fixmate.bam",
        "flagstat ex1.bam > %(out)s_ex1.flagstat",
        # Fails python 3.3 on linux, passes on OsX and when
        # run locally
        "calmd ex1.bam ex1.fa > %(out)s_ex1.calmd.bam",
        # use -s option, otherwise the following error in samtools 1.2:
        # Samtools-htslib-API: bam_get_library() not yet implemented
        # causes downstream problems
        # TODO: The following cause subsequent commands to fail
        # unknow option
        # "rmdup -s ex1.bam %(out)s_ex1.rmdup.bam",
        # "merge -f %(out)s_ex1.merge.bam ex1.bam ex1.bam",
        "reheader ex1.sam ex1.bam > %(out)s_ex1.reheader",
        "cat -o %(out)s_ex1.cat.bam ex1.bam ex1.bam",
        "targetcut ex1.bam > %(out)s_ex1.targetcut",
        "phase ex1.bam > %(out)s_ex1.phase",
        "import ex1.fa.fai ex1.sam.gz %(out)s_ex1.bam",
        "bam2fq ex1.bam > %(out)s_ex1.bam2fq",
        # TODO: not the same
        # "pad2unpad -T ex1.fa ex2.bam > %(out)s_ex2.unpad",
        # TODO: command line option problem
        # "bamshuf ex1.bam -O --output-fmt SAM > %(out)s_ex1.bamshuf.sam",
        # "collate ex1.bam %(out)s_ex1.collate",
        "bedcov ex1.bed ex1.bam > %(out)s_ex1.bedcov",
        "stats ex1.bam > %(out)s_ex1.stats",
        "dict ex1.bam > %(out)s_ex1.dict",
        # TODO: not the same
        # ("addreplacerg -r 'RG\tID:ga\tSM:hs' ex1.bam > %(out)s_ex1.addreplacerg",
    ]

    map_command = {
        "import": "samimport"}

    executable = "samtools"

    module = pysam.samtools

    def check_version(self):

        samtools_version = get_version(self.executable)
        def _r(s):
            # patch - remove any of the alpha/beta suffixes, i.e., 0.1.12a ->
            # 0.1.12
            if s.count('-') > 0:
                s = s[0:s.find('-')]
            return re.sub("[^0-9.]", "", s)

        if _r(samtools_version) != _r(pysam.__samtools_version__):
            raise ValueError(
                "versions of pysam.%s and %s differ: %s != %s" %
                (self.executable,
                 self.executable,
                 pysam.__samtools_version__,
                 samtools_version))

    def setUp(self):
        '''setup tests.

        For setup, all commands will be run before the first test is
        executed. Individual tests will then just compare the output
        files.

        '''

        self.check_version()

        if not os.path.exists(WORKDIR):
            os.makedirs(WORKDIR)

        for f in self.requisites:
            shutil.copy(os.path.join(DATADIR, f),
                        os.path.join(WORKDIR, f))

        self.savedir = os.getcwd()
        os.chdir(WORKDIR)

        return

    def check_statement(self, statement):

        parts = statement.split(" ")
        r_samtools = {"out": self.executable}
        r_pysam = {"out": "pysam"}

        command = parts[0]
        command = self.map_command.get(command, command)
        # self.assertTrue(command in pysam.SAMTOOLS_DISPATCH)

        targets = [x for x in parts if "%(out)s" in x]
        samtools_targets = [x % r_samtools for x in targets]
        pysam_targets = [x % r_pysam for x in targets]

        pysam_method = getattr(self.module, command)
        # run samtools
        full_statement = re.sub("%\(out\)s", self.executable, statement)
        run_command(" ".join((self.executable, full_statement)))
        # sys.stdout.write("%s %s ok" % (command, self.executable))

        # run pysam
        if ">" in statement:
            assert parts[-2] == ">"
            parts = parts[:-2]

        # avoid interpolation to preserve string quoting, tab chars, etc.
        pysam_parts = [re.sub("%\(out\)s", "pysam", x) for x in parts[1:]]
        output = pysam_method(*pysam_parts,
                              raw=True,
                              catch_stdout=True)
        # sys.stdout.write(" pysam ok\n")
        if ">" in statement:
            with open(pysam_targets[-1], "wb") as outfile:
                if output is not None:
                    outfile.write(force_bytes(output))

        for samtools_target, pysam_target in zip(samtools_targets,
                                                 pysam_targets):
            if os.path.isdir(samtools_target):
                samtools_files = glob.glob(os.path.join(
                    samtools_target, "*"))
                pysam_files = glob.glob(os.path.join(pysam_target, "*"))
                self.assertEqual(len(samtools_files), len(pysam_files))
                # need to be able to exclude files like README, etc.
                continue
            else:
                samtools_files = [samtools_target]
                pysam_files = [pysam_target]

            for s, p in zip(samtools_files, pysam_files):
                binary_equal = checkBinaryEqual(s, p)
                error_msg = "%s failed: files %s and %s are not the same" % (command, s, p)
                if binary_equal:
                    continue
                if s.endswith(".bam"):
                    self.assertTrue(
                        check_samtools_view_equal(
                            s, p, without_header=True),
                        error_msg)
                check_lines_equal(
                    self, s, p,
                    filter_f=lambda x: x.startswith("#"),
                    msg=error_msg)

    def testStatements(self):
        for statement in self.statements:
            if (statement.startswith("calmd") and 
                list(sys.version_info[:2]) == [3, 3]):
                # skip calmd test, fails only on python 3.3.5
                # in linux (empty output). Works in OsX and passes
                # for 3.4 and 3.5, see issue #293
                continue
            self.check_statement(statement)

    def tearDown(self):
        if os.path.exists(WORKDIR):
            shutil.rmtree(WORKDIR)
        os.chdir(self.savedir)


class EmptyIndexTest(unittest.TestCase):

    def testEmptyIndex(self):
        self.assertRaises(IOError, pysam.samtools.index,
                          "exdoesntexist.bam")

class TestReturnType(unittest.TestCase):
    
    def testReturnValueString(self):
        retval = pysam.idxstats(os.path.join(DATADIR, "ex1.bam"))
        if IS_PYTHON3:
            self.assertFalse(isinstance(retval, bytes))
            self.assertTrue(isinstance(retval, str))
        else:
            self.assertTrue(isinstance(retval, bytes))
            self.assertTrue(isinstance(retval, basestring))

    def testReturnValueData(self):
        args = "-O BAM {}".format(os.path.join(DATADIR, "ex1.bam")).split(" ")
        retval = pysam.view(*args)

        if IS_PYTHON3:
            self.assertTrue(isinstance(retval, bytes))
            self.assertFalse(isinstance(retval, str))
        else:
            self.assertTrue(isinstance(retval, bytes))
            self.assertTrue(isinstance(retval, basestring))


class StdoutTest(unittest.TestCase):
    '''test if stdout can be redirected.'''

    def testWithRedirectedStdout(self):
        r = pysam.samtools.flagstat(
            os.path.join(DATADIR, "ex1.bam"))
        self.assertTrue(len(r) > 0)

    def testWithoutRedirectedStdout(self):
        r = pysam.samtools.flagstat(
            os.path.join(DATADIR, "ex1.bam"),
            catch_stdout=False)
        self.assertEqual(r, None)

    def testDoubleCalling(self):
        # The following would fail if there is an
        # issue with stdout being improperly caught.
        retvals = pysam.idxstats(
            os.path.join(DATADIR, "ex1.bam"))
        retvals = pysam.idxstats(
            os.path.join(DATADIR, "ex1.bam"))

    def testSaveStdout(self):
        outfile = get_temp_filename(suffix=".tsv")
        r = pysam.samtools.flagstat(
            os.path.join(DATADIR, "ex1.bam"),
            save_stdout=outfile)
        self.assertEqual(r, None)
        with open(outfile) as inf:
            r = inf.read()
        self.assertTrue(len(r) > 0)


class PysamTest(SamtoolsTest):
    """check access to samtools command in the pysam 
    main package.

    This is for backwards capability.
    """

    module = pysam


class BcftoolsTest(SamtoolsTest):

    requisites = [
        "ex1.fa",
        "ex1.vcf.gz",
        "ex1.vcf.gz.tbi",
    ]
    # a list of statements to test
    # should contain at least one %(out)s component indicating
    # an output file.
    statements = [
        # "index -n ex1.vcf.gz > %(out)s_ex1.index",

        "annotate -x ID ex1.vcf.gz > %(out)s_ex1.annotate",
        "concat -a ex1.vcf.gz ex1.vcf.gz > %(out)s_ex1.concat",
        "isec -p %(out)s_ex1.isec ex1.vcf.gz ex1.vcf.gz",
        "merge --force-samples ex1.vcf.gz ex1.vcf.gz > %(out)s_ex1.norm",
        "norm -m +both ex1.vcf.gz > %(out)s_ex1.norm",

        # "plugin",
        # "query -f '%CHROM\n' ex1.vcf.gz > %(out)s_ex1.query",
        # "reheader -s A > %(out)s_ex1.reheader",
        # "view ex1.vcf.gz > %(out)s_ex1.view",
        # "call -m ex1.vcf.gz > %(out)s_ex1.call",
        # bad file descriptor
        # "consensus -f ex1.fa ex1.vcf.gz  > %(out)s_ex1.consensus"
        # need appropriate VCF file
        # "cnv",
        # segfault
        # "filter -s A ex1.vcf.gz  > %(out)s_ex1.filter",
        # exit
        # "gtcheck -s A ex1.vcf.gz  > %(out)s_ex1.gtcheck",
        "roh -s A ex1.vcf.gz > %(out)s_ex1.roh",
        "stats ex1.vcf.gz > %(out)s_ex1.stats",
    ]

    map_command = {
        "import": "samimport"}

    executable = "bcftools"

    module = pysam.bcftools


if __name__ == "__main__":
    # build data files
    print ("building data files")
    subprocess.call("make -C %s" % DATADIR, shell=True)
    print ("starting tests")
    unittest.main()
    print ("completed tests")
