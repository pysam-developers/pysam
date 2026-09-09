import warnings
import os
import re
import glob
import pytest
import sys
import subprocess
import shutil

import pysam
import pysam.samtools
import pysam.bcftools
from TestUtils import checkBinaryEqual, slurp_file, \
    check_samtools_view_equal, force_bytes, \
    make_data_files, BAM_DATADIR


def setUpModule():
    make_data_files(BAM_DATADIR)


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
        lines = b"".join(pipe.readlines()).decode("ascii")

    try:
        x = re.search(r"Version:\s+(\S+)", lines).groups()[0]
    except AttributeError:
        raise ValueError("could not get version from %s" % lines)
    return x


class TestSamtools:

    '''test samtools command line commands and compare
    against pysam commands.

    Tests fail, if the output is not binary identical.
    '''

    requisites = [
        "ex1.fa", "ex1.fa.fai",
        "ex1.sam.gz",
        "ex1.bam", "ex1.bam.bai",
        "ex1.sam",
        "ex1.sam",
        "ex2.bam",
        "ex2.sam",
        "ex1.bed"]

    # a list of statements to test
    # should contain at least one %(out)s component indicating
    # an output file.
    statements = [
        "view ex1.bam > %(out)s_ex1.view",
        "view -c ex1.bam > %(out)s_ex1.count",
        # ("view -bT ex1.fa -o %(out)s_ex1.view2 ex1.sam",
        "sort ex1.bam -o %(out)s_ex1.sort.bam",
        "mpileup ex1.bam > %(out)s_ex1.pileup",
        "depth ex1.bam > %(out)s_ex1.depth",
        # TODO: issues with file naming
        # "faidx ex1.fa; %(out)s_ex1.fa.fai",
        "index ex1.bam %(out)s_ex1.bam.fai",
        "index -@2 ex1.bam %(out)s_ex1.bam.fai",
        "idxstats ex1.bam > %(out)s_ex1.idxstats",
        # TODO: fixmate behaviour changed in 1.21
        # "fixmate ex1.bam %(out)s_ex1.fixmate.bam",
        "flagstat ex1.bam > %(out)s_ex1.flagstat",
        "calmd ex1.bam ex1.fa > %(out)s_ex1.calmd.bam",
        # use -s option, otherwise the following error in samtools 1.2:
        # Samtools-htslib-API: bam_get_library() not yet implemented
        # causes downstream problems
        # TODO: The following cause subsequent commands to fail
        # unknown option
        # "rmdup -s ex1.bam %(out)s_ex1.rmdup.bam",
        # "merge -f %(out)s_ex1.merge.bam ex1.bam ex1.bam",
        "reheader ex2.sam ex1.bam > %(out)s_ex1.reheader.bam",
        "cat -o %(out)s_ex1.cat.bam ex1.bam ex1.bam",
        "targetcut ex1.bam > %(out)s_ex1.targetcut",
        "phase ex1.bam > %(out)s_ex1.phase",
        "view -bt ex1.fa.fai ex1.sam.gz > %(out)s_ex1.bam",
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
    }

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
            warnings.warn(
                "versions of pysam.%s and %s differ: %s != %s" %
                (self.executable,
                 self.executable,
                 pysam.__samtools_version__,
                 samtools_version))

    @pytest.fixture(autouse=True)
    def setup_method(self, tmp_path):
        '''setup tests.

        For setup, all commands will be run before the first test is
        executed. Individual tests will then just compare the output
        files.

        '''
        self.check_version()

        self.workdir = str(tmp_path)

        for f in self.requisites:
            shutil.copy(os.path.join(BAM_DATADIR, f),
                        os.path.join(self.workdir, f))

        self.savedir = os.getcwd()
        os.chdir(self.workdir)

        yield

        # Unlike most other tests, do remove these voluminous output files
        for f in glob.glob(os.path.join(self.workdir, "*")):
            os.unlink(f)

        os.chdir(self.savedir)

    def get_command(self, statement, map_to_internal=True):
        """return samtools command from statement"""
        parts = statement.split(" ")
        command = parts[0]
        if map_to_internal:
            return self.map_command.get(command, command)
        else:
            return command

    def check_statement(self, statement):

        parts = statement.split(" ")
        r_samtools = {"out": self.executable}
        r_pysam = {"out": "pysam"}

        command = self.get_command(statement)

        targets = [x for x in parts if "%(out)s" in x]
        samtools_targets = [x % r_samtools for x in targets]
        pysam_targets = [x % r_pysam for x in targets]

        pysam_method = getattr(self.module, command)

        # run samtools
        full_statement = re.sub(r"%\(out\)s", self.executable, statement)
        run_command(" ".join((self.executable, full_statement)))
        # sys.stdout.write("%s %s ok" % (command, self.executable))

        # run pysam
        if ">" in statement:
            assert parts[-2] == ">"
            parts = parts[:-2]

        # avoid interpolation to preserve string quoting, tab chars, etc.
        pysam_parts = [re.sub(r"%\(out\)s", "pysam", x) for x in parts[1:]]
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
                assert len(samtools_files) == len(pysam_files)
                # need to be able to exclude files like README, etc.
                continue
            else:
                samtools_files = [samtools_target]
                pysam_files = [pysam_target]

            for s, p in zip(samtools_files, pysam_files):
                binary_equal = checkBinaryEqual(s, p)
                error_msg = f"{command} failed: files {s} and {p} are not the same"
                if binary_equal:
                    continue
                elif s.endswith(".bam"):
                    assert check_samtools_view_equal(s, p, without_header=True), error_msg
                else:
                    lines_s = slurp_file(s, omit=lambda x: x.startswith("#"))
                    lines_p = slurp_file(p, omit=lambda x: x.startswith("#"))
                    assert lines_s == lines_p, error_msg

    def testStatements(self):
        for statement in self.statements:
            command = self.get_command(statement, map_to_internal=False)
            # bam2fq differs between version 1.5 and 1.6 - re-enable if
            # bioconda samtools will be available.
            # flagstat differs between version <=1.12 and >=1.13
            if command in ("bedcov", "stats", "dict", "bam2fq", "flagstat"):
                continue

            self.check_statement(statement)

    @pytest.mark.skipif(not sys.stdin.isatty(), reason="skipping usage tests, stdin is not a tty")
    def testUsage(self):
        if self.executable == "bcftools":
            # bcftools usage messages end with exit(1)
            return

        for statement in self.statements:
            command = self.get_command(statement, map_to_internal=False)
            # ignore commands that exit or cause other failures
            # TODO: check - if reheader or phase is run in testStatements, sort fails
            # here
            if command in ("view", "sort", "bam2fq", "flagstat", "reheader",
                           "stats", "idxstats"):
                continue
            mapped_command = self.get_command(statement, map_to_internal=True)
            pysam_method = getattr(self.module, mapped_command)
            usage_msg = pysam_method.usage()
            expected = r"Usage:\s+{} {}".format(self.executable, command)
            assert re.search(expected, usage_msg) is not None


class TestEmptyIndex:
    def testEmptyIndex(self):
        with pytest.raises(IOError):
            pysam.samtools.index("exdoesntexist.bam")

    def testEmptyIndexWithExtraFlag(self):
        with pytest.raises(IOError):
            pysam.samtools.index("-c", "exdoesntexist.bam")

    def testEmptyIndexWithExtraArg(self):
        with pytest.raises(IOError):
            pysam.samtools.index("-c", "-m", "14", "exdoesntexist.bam")


class TestSubcommands:
    def testFailingSamtools(self):
        with pytest.raises(pysam.SamtoolsError):
            pysam.samtools.view("nonexistent.bam")

    def testFailingBCFtools(self):
        with pytest.raises(pysam.SamtoolsError):
            pysam.bcftools.view("nonexistent.vcf")


class TestReturnType:
    def testReturnValueString(self):
        retval = pysam.idxstats(os.path.join(BAM_DATADIR, "ex1.bam"))
        assert not isinstance(retval, bytes)
        assert isinstance(retval, str)

    def testReturnValueData(self):
        retval = pysam.view("-O", "BAM", os.path.join(BAM_DATADIR, "ex1.bam"))
        assert isinstance(retval, bytes)
        assert not isinstance(retval, str)


class TestStdout:
    '''test if stdout can be redirected.'''

    def testWithRedirectedStdout(self):
        r = pysam.samtools.flagstat(
            os.path.join(BAM_DATADIR, "ex1.bam"))
        assert len(r) > 0

    def testWithoutRedirectedStdout(self):
        r = pysam.samtools.flagstat(
            os.path.join(BAM_DATADIR, "ex1.bam"),
            catch_stdout=False)
        assert r is None

    def testDoubleCalling(self):
        # The following would fail if there is an
        # issue with stdout being improperly caught.
        retvals = pysam.idxstats(
            os.path.join(BAM_DATADIR, "ex1.bam"))
        retvals = pysam.idxstats(
            os.path.join(BAM_DATADIR, "ex1.bam"))

    def testSaveStdout(self, tmp_path):
        outfile = str(tmp_path / "flagstat.tsv")
        r = pysam.samtools.flagstat(
            os.path.join(BAM_DATADIR, "ex1.bam"),
            save_stdout=outfile)
        assert r is None
        with open(outfile) as inf:
            r = inf.read()
        assert len(r) > 0


class TestBcftoolsOutput:
    @pytest.fixture
    def norm_inputs(self, tmp_path):
        reference = tmp_path / "reference.fa"
        reference.write_text(">chr1\nAAAAAAAAAA\n")
        pysam.faidx(str(reference))

        vcf = tmp_path / "input.vcf"
        vcf.write_text(
            "##fileformat=VCFv4.2\n"
            "##contig=<ID=chr1,length=10>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t2\t.\tAA\tA\t.\tPASS\t.\n"
        )
        compressed_vcf = tmp_path / "input.vcf.gz"
        pysam.tabix_compress(str(vcf), str(compressed_vcf))
        pysam.tabix_index(str(compressed_vcf), preset="vcf")
        return reference, compressed_vcf

    @pytest.mark.parametrize("output_option", ["-o", "--output"])
    def test_norm_honours_explicit_output_path(self, norm_inputs, tmp_path, output_option):
        reference, compressed_vcf = norm_inputs
        output = tmp_path / "normalized.vcf"

        expected = pysam.bcftools.norm(
            "-Ov", "--no-version", "-f", str(reference), str(compressed_vcf), raw=True
        )
        returned = pysam.bcftools.norm(
            "-Ov",
            "--no-version",
            output_option,
            str(output),
            "-f",
            str(reference),
            str(compressed_vcf),
            raw=True,
        )

        assert expected
        assert returned in (None, "")
        assert output.exists()
        assert output.read_text() == expected

    def test_norm_honours_equals_output_path(self, norm_inputs, tmp_path):
        reference, compressed_vcf = norm_inputs
        output = tmp_path / "normalized.vcf"

        pysam.bcftools.norm(
            "-Ov",
            "--no-version",
            f"--output={output}",
            "-f",
            str(reference),
            str(compressed_vcf),
            raw=True,
        )

        assert output.exists()


class TestPysam(TestSamtools):
    """check access to samtools command in the pysam
    main package.

    This is for backwards capability.
    """

    module = pysam


# class TestBcftools(TestSamtools):

#     requisites = [
#         "ex1.fa",
#         "ex1.vcf.gz",
#         "ex1.vcf.gz.tbi",
#     ]
#     # a list of statements to test
#     # should contain at least one %(out)s component indicating
#     # an output file.
#     statements = [
#         # "index -n ex1.vcf.gz > %(out)s_ex1.index",

#         "annotate -x ID ex1.vcf.gz > %(out)s_ex1.annotate",
#         "concat -a ex1.vcf.gz ex1.vcf.gz > %(out)s_ex1.concat",
#         "isec -p %(out)s_ex1.isec ex1.vcf.gz ex1.vcf.gz",
#         "merge --force-samples ex1.vcf.gz ex1.vcf.gz > %(out)s_ex1.norm",
#         "norm -m +both ex1.vcf.gz > %(out)s_ex1.norm",

#         # "plugin",
#         # "query -f '%CHROM\n' ex1.vcf.gz > %(out)s_ex1.query",
#         # "reheader -s A > %(out)s_ex1.reheader",
#         # "view ex1.vcf.gz > %(out)s_ex1.view",
#         # "call -m ex1.vcf.gz > %(out)s_ex1.call",
#         # bad file descriptor
#         # "consensus -f ex1.fa ex1.vcf.gz  > %(out)s_ex1.consensus"
#         # need appropriate VCF file
#         # "cnv",
#         # segfault
#         # "filter -s A ex1.vcf.gz  > %(out)s_ex1.filter",
#         # exit
#         # "gtcheck -s A ex1.vcf.gz  > %(out)s_ex1.gtcheck",
#         # segfault, used to work with bcftools 1.3
#         # "roh -s A ex1.vcf.gz > %(out)s_ex1.roh",
#         "stats ex1.vcf.gz > %(out)s_ex1.stats",
#     ]

#     map_command = {
#     }

#     executable = "bcftools"

#     module = pysam.bcftools
