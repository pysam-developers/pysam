"""Benchmarking module for AlignmentFile functionality"""
import os

from TestUtils import BAM_DATADIR, force_str, flatten_nested_list
from PileupTestUtils import *


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


def test_build_depth_from_bam_with_samtools(benchmark):
    result = benchmark(build_depth_with_samtools,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 107248


def test_build_depth_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_depth_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert sum(result) == 107248


def test_build_depth_from_bam_with_pysam(benchmark):
    result = benchmark(build_depth_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    # different value, as samtools filters with a minimum
    # base quality of 13
    assert sum(result) == 110015


def test_build_depth_with_filter_from_bam_with_pysam(benchmark):
    result = benchmark(build_depth_with_filter_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert sum(result) == 107248


def test_build_query_bases_from_bam_with_samtools(benchmark):
    result = benchmark(build_query_bases_with_samtools,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 116314


def test_build_query_bases_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_query_bases_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len("".join(flatten_nested_list(result))) == 116314


def test_build_query_bases_from_bam_with_pysam_pileups(benchmark):
    result = benchmark(build_query_bases_with_pysam_pileups,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len("".join(flatten_nested_list(result))) == 116314


def test_build_query_bases_from_bam_with_pysam(benchmark):
    result = benchmark(build_query_bases_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len("".join(flatten_nested_list(result))) == 116314


def test_build_query_qualities_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_query_qualities_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len("".join(result)) == 107248


def test_build_query_qualities_from_bam_with_pysam(benchmark):
    result = benchmark(build_query_qualities_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert sum([len(x) for x in result]) == 107248


def test_build_query_names_from_bam_with_Pysam(benchmark):
    result = benchmark(build_query_names_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len("".join([x for column in result for x in column])) == 2307497


def test_build_mapping_qualities_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_mapping_qualities_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len("".join(result)) == 107248


def test_build_mapping_qualities_from_bam_with_pysam(benchmark):
    result = benchmark(build_mapping_qualities_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert sum([len(x) for x in result]) == 107248


def test_build_query_positions_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_query_positions_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    # positions output by samtools are 1-based
    assert sum([sum(x) - len(x) for x in result]) == 1841736


def test_build_query_positions_from_bam_with_pysam(benchmark):
    result = benchmark(build_query_positions_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert sum([sum(x) for x in result]) == 1841736