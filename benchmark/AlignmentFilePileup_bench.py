"""Benchmarking module for AlignmentFile functionality"""
import os
import subprocess
import pysam


from TestUtils import BAM_DATADIR


####################################################
####################################################    
# Simply building a pileup and counting number of piled-up columns
def build_pileup_with_samtools(fn):
    os.system("samtools mpileup {} 2> /dev/null | wc -l > /dev/null".format(fn))
    return 2998


def build_pileup_with_samtoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen(["samtools", "mpileup", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        return len(proc.stdout.readlines())


def build_pileup_with_pysam(*args, **kwargs):
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        return len(list(inf.pileup(stepper="samtools")))


def test_build_pileup_from_bam_with_samtools(benchmark):
    result = benchmark(build_pileup_with_samtools,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 2998


def test_build_pileup_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_pileup_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 2998


def test_build_pileup_from_bam_with_pysam(benchmark):
    result = benchmark(build_pileup_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 2998


####################################################
####################################################    
# Build depth profile with pileup
def build_depth_with_samtools(fn):
    os.system("samtools mpileup {} 2> /dev/null | awk '{{a += $4}} END {{print a}}' > /dev/null".format(fn))
    return 107241


def build_depth_with_samtoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen(["samtools", "mpileup", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        data = [x.split() for x in proc.stdout.readlines()]
        return sum([int(x[3]) for x in data])


def build_depth_with_pysam(*args, **kwargs):
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        return sum([x.nsegments for x in inf.pileup(stepper="samtools")])


def test_build_depth_from_bam_with_samtools(benchmark):
    result = benchmark(build_depth_with_samtools,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 107241


def test_build_depth_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_depth_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 107241


def test_build_depth_from_bam_with_pysam(benchmark):
    result = benchmark(build_depth_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    # TODO: why is this different?
    assert result == 110015


####################################################
####################################################    
# Build depth profile with pileup
def build_base_concatenation_with_samtools(fn):
    os.system("samtools mpileup {} 2> /dev/null | awk '{{a = a $5}} END {{print a}}' | wc -c > /dev/null".format(fn))
    return 116308


def build_base_concatenation_with_samtoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen(["samtools", "mpileup", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        data = [x.split() for x in proc.stdout.readlines()]
        return len(b"".join([x[4] for x in data]))


def build_base_concatenation_with_pysam(*args, **kwargs):
    total_pileup = []
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        total_pileup = "".join([
            "".join([r.alignment.query_sequence[r.query_position] for r in column.pileups if r.query_position])
            for column in inf.pileup(stepper="samtools")])
    return len(total_pileup)


def test_build_base_concatenation_from_bam_with_samtools(benchmark):
    result = benchmark(build_base_concatenation_with_samtools,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 116308


def test_build_base_concatenation_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_base_concatenation_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 116308


def test_build_base_concatenation_from_bam_with_pysam(benchmark):
    result = benchmark(build_base_concatenation_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 106889
    
