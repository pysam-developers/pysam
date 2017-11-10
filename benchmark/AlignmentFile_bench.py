"""Benchmarking module for AlignmentFile functionality"""
import os
import subprocess
import pysam


from TestUtils import BAM_DATADIR


def count_number_lines_with_samtools(fn):
    os.system("samtools view {} | wc -l > /dev/null".format(fn))
    return 3270


def count_number_lines_with_samtoolspipe(fn):
    with subprocess.Popen(["samtools", "view", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE) as proc:
        return len(proc.stdout.readlines())


def count_number_lines_with_pysam(*args, **kwargs):
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        return len(list(inf.fetch()))


def test_count_number_lines_from_sam_with_samtools(benchmark):
    result = benchmark(count_number_lines_with_samtools,
                       os.path.join(BAM_DATADIR, "ex2.sam"))
    assert result == 3270


def test_count_number_lines_from_sam_with_samtoolspipe(benchmark):
    result = benchmark(count_number_lines_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.sam"))
    assert result == 3270


def test_count_number_lines_from_sam_with_pysam(benchmark):
    result = benchmark(count_number_lines_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.sam"), "r")
    assert result == 3270


def test_count_number_lines_from_bam_with_samtools(benchmark):
    result = benchmark(count_number_lines_with_samtools,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 3270


def test_count_number_lines_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(count_number_lines_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 3270


def test_count_number_lines_from_bam_with_pysam(benchmark):
    result = benchmark(count_number_lines_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 3270
