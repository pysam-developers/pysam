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


def runSamtools(cmd):
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

    first_time = True

    # a dictionary of commands to test
    # first entry: (samtools output file, samtools command)
    # second entry: (pysam output file, (pysam function, pysam options) )
    commands = \
        {
            "view":
            (
                ("ex1.view", "view ex1.bam > ex1.view"),
                ("pysam_ex1.view", (pysam.view, "ex1.bam")),
            ),
            "view2":
            (
                ("ex1.view", "view -bT ex1.fa -o ex1.view2 ex1.sam"),
                # note that -o ex1.view2 throws exception.
                ("pysam_ex1.view",
                 (pysam.view, "-bT ex1.fa -oex1.view2 ex1.sam")),
            ),
            "sort":
            (
                ("ex1.sort.bam", "sort ex1.bam ex1.sort"),
                ("pysam_ex1.sort.bam", (pysam.sort, "ex1.bam pysam_ex1.sort")),
            ),
            "mpileup":
            (
                ("ex1.pileup", "mpileup ex1.bam > ex1.pileup"),
                ("pysam_ex1.mpileup", (pysam.mpileup, "ex1.bam")),
            ),
            "depth":
            (
                ("ex1.depth", "depth ex1.bam > ex1.depth"),
                ("pysam_ex1.depth", (pysam.depth, "ex1.bam")),
            ),
            "faidx":
            (
                ("ex1.fa.fai", "faidx ex1.fa"),
                ("pysam_ex1.fa.fai", (pysam.faidx, "ex1.fa")),
            ),
            "index":
            (
                ("ex1.bam.bai", "index ex1.bam"),
                ("pysam_ex1.bam.bai", (pysam.index, "pysam_ex1.bam")),
            ),
            "idxstats":
            (
                ("ex1.idxstats", "idxstats ex1.bam > ex1.idxstats"),
                ("pysam_ex1.idxstats", (pysam.idxstats, "pysam_ex1.bam")),
            ),
            "fixmate":
            (
                ("ex1.fixmate", "fixmate ex1.bam ex1.fixmate"),
                ("pysam_ex1.fixmate",
                 (pysam.fixmate, "pysam_ex1.bam pysam_ex1.fixmate")),
            ),
            "flagstat":
            (
                ("ex1.flagstat", "flagstat ex1.bam > ex1.flagstat"),
                ("pysam_ex1.flagstat", (pysam.flagstat, "pysam_ex1.bam")),
            ),
            "calmd":
            (
                ("ex1.calmd", "calmd ex1.bam ex1.fa > ex1.calmd"),
                ("pysam_ex1.calmd", (pysam.calmd, "pysam_ex1.bam ex1.fa")),
            ),
            "merge":
            (
                ("ex1.merge", "merge -f ex1.merge ex1.bam ex1.bam"),
                # -f option does not work - following command will cause the subsequent
                # command to fail
                ("pysam_ex1.merge",
                 (pysam.merge, "pysam_ex1.merge pysam_ex1.bam pysam_ex1.bam")),
            ),
            "rmdup":
            (
                ("ex1.rmdup", "rmdup ex1.bam ex1.rmdup"),
                ("pysam_ex1.rmdup",
                 (pysam.rmdup, "pysam_ex1.bam pysam_ex1.rmdup")),
            ),
            "reheader":
            (
                ("ex1.reheader", "reheader ex1.bam ex1.bam > ex1.reheader"),
                ("pysam_ex1.reheader", (pysam.reheader, "ex1.bam ex1.bam")),
            ),
            "cat":
            (
                ("ex1.cat", "cat ex1.bam ex1.bam > ex1.cat"),
                ("pysam_ex1.cat", (pysam.cat, "ex1.bam ex1.bam")),
            ),
            "targetcut":
            (
                ("ex1.targetcut", "targetcut ex1.bam > ex1.targetcut"),
                ("pysam_ex1.targetcut", (pysam.targetcut, "pysam_ex1.bam")),
            ),
            "phase":
            (
                ("ex1.phase", "phase ex1.bam > ex1.phase"),
                ("pysam_ex1.phase", (pysam.phase, "pysam_ex1.bam")),
            ),
            "import":
            (
                ("ex1.bam", "import ex1.fa.fai ex1.sam.gz ex1.bam"),
                ("pysam_ex1.bam",
                 (pysam.samimport, "ex1.fa.fai ex1.sam.gz pysam_ex1.bam")),
            ),
            "bam2fq":
            (
                ("ex1.bam2fq", "bam2fq ex1.bam > ex1.bam2fq"),
                ("pysam_ex1.bam2fq", (pysam.bam2fq, "pysam_ex1.bam")),
            ),
            "pad2unpad":
            (
                ("ex2.unpad", "pad2unpad -T ex1.fa ex2.bam > ex2.unpad"),
                ("pysam_ex2.unpad", (pysam.pad2unpad, "-T ex1.fa ex2.bam")),
            ),
            "bamshuf":
            (
                ("ex1.bamshuf.bam", "bamshuf ex1.bam ex1.bamshuf"),
                ("pysam_ex1.bamshuf.bam",
                 (pysam.bamshuf, "ex1.bam pysam_ex1.bamshuf")),
            ),
            "bedcov":
            (
                ("ex1.bedcov", "bedcov ex1.bed ex1.bam > ex1.bedcov"),
                ("pysam_ex1.bedcov", (pysam.bedcov, "ex1.bed ex1.bam")),
            ),
        }

    # some tests depend on others. The order specifies in which order
    # the samtools commands are executed.
    # The first three (faidx, import, index) need to be in that order,
    # the rest is arbitrary.
    order = ('faidx', 'import', 'index',
             # 'pileup1', 'pileup2', deprecated
             # 'glfview', deprecated
             'view', 'view2',
             'sort',
             'mpileup',
             'depth',
             'idxstats',
             # 'fixmate',
             'flagstat',
             # 'calmd',
             'merge',
             # 'rmdup',
             'reheader',
             'cat',
             'bedcov',
             'targetcut',
             'phase',
             # 'bamshuf',
             'bam2fq',
             # 'pad2unpad',
             )

    def setUp(self):
        '''setup tests. 

        For setup, all commands will be run before the first test is
        executed. Individual tests will then just compare the output
        files.
        '''
        if BinaryTest.first_time:

            # remove previous files
            if os.path.exists(WORKDIR):
                shutil.rmtree(WORKDIR)
                pass

            # copy the source files to WORKDIR
            os.makedirs(WORKDIR)

            for f in ("ex1.fa", "ex1.sam.gz",
                      "ex1.sam", "ex2.bam",
                      "ex1.bed"):
                shutil.copy(os.path.join(DATADIR, f),
                            os.path.join(WORKDIR, f))

            # cd to workdir
            savedir = os.getcwd()
            os.chdir(WORKDIR)
            for label in self.order:
                # print ("command=", label)
                command = self.commands[label]
                # build samtools command and target and run
                samtools_target, samtools_command = command[0]
                runSamtools(" ".join((SAMTOOLS, samtools_command)))

                # get pysam command and run
                try:
                    pysam_target, pysam_command = command[1]
                except ValueError as msg:
                    raise ValueError("error while setting up %s=%s: %s" %
                                     (label, command, msg))

                pysam_method, pysam_options = pysam_command

                try:
                    output = pysam_method(*pysam_options.split(" "), raw=True)
                except pysam.SamtoolsError as msg:
                    raise pysam.SamtoolsError(
                        "error while executing %s: options=%s: msg=%s" %
                        (label, pysam_options, msg))

                if ">" in samtools_command:
                    with open(pysam_target, "wb") as outfile:
                        if type(output) == list:
                            if IS_PYTHON3:
                                for line in output:
                                    outfile.write(line.encode('ascii'))
                            else:
                                for line in output:
                                    outfile.write(line)
                        else:
                            outfile.write(output)

            os.chdir(savedir)
            BinaryTest.first_time = False

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

    def checkCommand(self, command):
        if command:
            samtools_target, pysam_target = self.commands[
                command][0][0], self.commands[command][1][0]
            samtools_target = os.path.join(WORKDIR, samtools_target)
            pysam_target = os.path.join(WORKDIR, pysam_target)
            self.assertTrue(
                checkBinaryEqual(samtools_target, pysam_target),
                "%s failed: files %s and %s are not the same" %
                (command, samtools_target, pysam_target))

    def testImport(self):
        self.checkCommand("import")

    def testIndex(self):
        self.checkCommand("index")

    def testSort(self):
        self.checkCommand("sort")

    def testMpileup(self):
        self.checkCommand("mpileup")

    def testDepth(self):
        self.checkCommand("depth")

    def testIdxstats(self):
        self.checkCommand("idxstats")

    # def testFixmate(self):
    #     self.checkCommand("fixmate")

    def testFlagstat(self):
        self.checkCommand("flagstat")

    def testMerge(self):
        self.checkCommand("merge")

    # def testRmdup(self):
    #     self.checkCommand("rmdup")

    def testReheader(self):
        self.checkCommand("reheader")

    def testCat(self):
        self.checkCommand("cat")

    def testTargetcut(self):
        self.checkCommand("targetcut")

    def testPhase(self):
        self.checkCommand("phase")

    def testBam2fq(self):
        self.checkCommand("bam2fq")

    def testBedcov(self):
        self.checkCommand("bedcov")

    # def testBamshuf(self):
    #     self.checkCommand("bamshuf")

    # def testPad2Unpad(self):
    #     self.checkCommand("pad2unpad")

    # def testPileup1( self ):
    #     self.checkCommand( "pileup1" )

    # def testPileup2( self ):
    #     self.checkCommand( "pileup2" )

    # deprecated
    # def testGLFView( self ):
    #     self.checkCommand( "glfview" )

    def testView(self):
        self.checkCommand("view")

    def testEmptyIndex(self):
        self.assertRaises(IOError, pysam.index, "exdoesntexist.bam")

    def __del__(self):
        if os.path.exists(WORKDIR):
            pass
        # shutil.rmtree( WORKDIR )


class StdoutTest(unittest.TestCase):
    '''test if stdout can be redirected.'''

    def testWithRedirectedStdout(self):
        r = pysam.flagstat(os.path.join(DATADIR, "ex1.bam"))
        self.assertTrue(len(r) > 0)

    def testWithoutRedirectedStdout(self):
        r = pysam.flagstat(os.path.join(DATADIR, "ex1.bam"),
                           catch_stdout=False)
        self.assertTrue(len(r) == 0)

if __name__ == "__main__":
    # build data files
    print ("building data files")
    subprocess.call("make -C %s" % DATADIR, shell=True)
    print ("starting tests")
    unittest.main()
    print ("completed tests")
