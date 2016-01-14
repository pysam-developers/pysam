#!/usr/bin/env python
'''unit testing code for pysam.

Execute in the :file:`tests` directory as it requires the Makefile
and data files located there.
'''

import pysam
import unittest
import os
import re
import sys
import subprocess
import shutil
from TestUtils import checkBinaryEqual

IS_PYTHON3 = sys.version_info[0] >= 3

SAMTOOLS = "samtools"
WORKDIR = "pysam_test_work"
DATADIR = "pysam_data"


def run_samtools(cmd):
    '''run a samtools command'''
    try:
        retcode = subprocess.call(cmd, shell=True,
                                  stderr=subprocess.PIPE)
        if retcode < 0:
            print("Child was terminated by signal", -retcode)
    except OSError as e:
        print("Execution failed:", e)


def getSamtoolsVersion():
    '''return samtools version'''

    with subprocess.Popen(SAMTOOLS, shell=True,
                          stderr=subprocess.PIPE).stderr as pipe:
        lines = b"".join(pipe.readlines())

    if IS_PYTHON3:
        lines = lines.decode('ascii')
    return re.search("Version:\s+(\S+)", lines).groups()[0]


class BinaryTest(unittest.TestCase):

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
        ("view ex1.bam > %(out)s_ex1.view"),
        # ("view -bT ex1.fa -o %(out)s_ex1.view2 ex1.sam"),
        ("sort ex1.bam -o %(out)s_ex1.sort.bam"),
        ("mpileup ex1.bam > %(out)s_ex1.pileup"),
        ("depth ex1.bam > %(out)s_ex1.depth"),
        # TODO: issues with file naming
        # ("faidx ex1.fa; %(out)s_ex1.fa.fai"),
        ("index ex1.bam %(out)s_ex1.bam.fai"),
        ("idxstats ex1.bam > %(out)s_ex1.idxstats"),
        ("fixmate ex1.bam %(out)s_ex1.fixmate.bam"),
        ("flagstat ex1.bam > %(out)s_ex1.flagstat"),
        ("calmd ex1.bam ex1.fa > %(out)s_ex1.calmd.bam"),
        # use -s option, otherwise the following error in samtools 1.2:
        # Samtools-htslib-API: bam_get_library() not yet implemented
        # causes downstream problems
        # TODO: The following cause subsequent commands to fail
        # unknow option
        # ("rmdup -s ex1.bam %(out)s_ex1.rmdup.bam"),
        # ("merge -f %(out)s_ex1.merge.bam ex1.bam ex1.bam"),
        ("reheader ex1.sam ex1.bam > %(out)s_ex1.reheader"),
        ("cat -o %(out)s_ex1.cat.bam ex1.bam ex1.bam"),
        ("targetcut ex1.bam > %(out)s_ex1.targetcut"),
        ("phase ex1.bam > %(out)s_ex1.phase"),
        ("import ex1.fa.fai ex1.sam.gz %(out)s_ex1.bam"),
        ("bam2fq ex1.bam > %(out)s_ex1.bam2fq"),
        # TODO: not the same
        # ("pad2unpad -T ex1.fa ex2.bam > %(out)s_ex2.unpad"),
        # TODO: command line option problem
        # ("bamshuf ex1.bam -O --output-fmt SAM > %(out)s_ex1.bamshuf.sam"),
        # ("collate ex1.bam %(out)s_ex1.collate"),
        ("bedcov ex1.bed ex1.bam > %(out)s_ex1.bedcov"),
        ("stats ex1.bam > %(out)s_ex1.stats"),
        ("dict ex1.bam > %(out)s_ex1.dict"),
        # TODO: not the same
        # ("addreplacerg -r 'RG\tID:ga\tSM:hs' ex1.bam > %(out)s_ex1.addreplacerg"),
    ]

    map_command = {
        "import": "samimport"}

    def setUp(self):
        '''setup tests. 

        For setup, all commands will be run before the first test is
        executed. Individual tests will then just compare the output
        files.

        '''
        if not os.path.exists(WORKDIR):
            os.makedirs(WORKDIR)

        for f in self.requisites:
            shutil.copy(os.path.join(DATADIR, f),
                        os.path.join(WORKDIR, f))

        samtools_version = getSamtoolsVersion()

        def _r(s):
            # patch - remove any of the alpha/beta suffixes, i.e., 0.1.12a ->
            # 0.1.12
            if s.count('-') > 0:
                s = s[0:s.find('-')]
            return re.sub("[^0-9.]", "", s)

        if _r(samtools_version) != _r(pysam.__samtools_version__):
            raise ValueError(
                "versions of pysam/samtools and samtools differ: %s != %s" %
                (pysam.__samtools_version__,
                 samtools_version))

        self.savedir = os.getcwd()
        os.chdir(WORKDIR)

        return

    def check_statement(self, statement):

        parts = statement.split(" ")
        r_samtools = {"out": "samtools"}
        r_pysam = {"out": "pysam"}

        command = parts[0]
        command = self.map_command.get(command, command)
        self.assertTrue(command in pysam.SAMTOOLS_DISPATCH)

        targets = [x for x in parts if "%(out)s" in x]
        samtools_targets = [x % r_samtools for x in targets]
        pysam_targets = [x % r_pysam for x in targets]

        pysam_method = getattr(pysam.samtools, command)

        # run samtools
        samtools_statement = statement % {"out": "samtools"}
        run_samtools(" ".join((SAMTOOLS, samtools_statement)))
        sys.stdout.write(" samtools ok")

        # run pysam
        if ">" in statement:
            assert parts[-2] == ">"
            parts = parts[:-2]

        # avoid interpolation to preserve string quoting, tab chars, etc.
        pysam_parts = [re.sub("%\(out\)s", "pysam", x) for x in parts[1:]]
        output = pysam_method(*pysam_parts,
                              raw=True,
                              catch_stdout=True)
        
        sys.stdout.write(" pysam ok\n")

        if ">" in statement:
            with open(pysam_targets[-1], "wb") as outfile:
                if type(output) == list:
                    if IS_PYTHON3:
                        for line in output:
                            outfile.write(line.encode('ascii'))
                    else:
                        for line in output:
                            outfile.write(line)
                else:
                    outfile.write(output)

        for samtools_target, pysam_target in zip(samtools_targets,
                                                 pysam_targets):
            self.assertTrue(
                checkBinaryEqual(samtools_target, pysam_target),
                "%s failed: files %s and %s are not the same" %
                (command, samtools_target, pysam_target))

    def testStatements(self):
        for statement in self.statements:
            self.check_statement(statement)
        
    def tearDown(self):
        if os.path.exists(WORKDIR):
            shutil.rmtree(WORKDIR)
        os.chdir(self.savedir)

class EmptyIndexTest(unittest.TestCase):

    def testEmptyIndex(self):
        self.assertRaises(IOError, pysam.samtools.index,
                          "exdoesntexist.bam")


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
        self.assertTrue(len(r) == 0)

if __name__ == "__main__":
    # build data files
    print ("building data files")
    subprocess.call("make -C %s" % DATADIR, shell=True)
    print ("starting tests")
    unittest.main()
    print ("completed tests")
