"""Benchmarking module for AlignmentFile functionality"""
import os
import pytest


from TestUtils import BAM_DATADIR, force_str, flatten_nested_list
from AlignmentFileFetchTestUtils import *


def test_build_fetch_from_bam_with_samtoolsshell(benchmark):
    result = benchmark(build_fetch_with_samtoolsshell,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 3270


def test_build_fetch_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_fetch_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 3270


def test_build_fetch_from_bam_with_pysam(benchmark):
    result = benchmark(build_fetch_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 3270


def test_build_query_sequences_from_bam_with_samtoolsshell(benchmark):
    result = benchmark(build_query_sequences_with_samtoolsshell,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len(result) == 3270


def test_build_query_sequences_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_query_sequences_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len(result) == 3270


def test_build_query_sequences_from_bam_with_pysam(benchmark):
    result = benchmark(build_query_sequences_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len(result) == 3270


def test_build_query_qualities_from_bam_with_pysam(benchmark):
    result = benchmark(build_query_qualities_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len(result) == 3270


def test_build_query_sequences_from_bam_flagfilter_with_samtoolsshell(benchmark):
    result = benchmark(build_query_sequences_flagfilter_with_samtoolsshell,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len(result) == 3124


def test_build_query_sequences_from_bam_flagfilter_with_samtoolspipe(benchmark):
    result = benchmark(build_query_sequences_flagfilter_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len(result) == 3124


def test_build_query_sequences_from_bam_flagfilter_with_pysam(benchmark):
    result = benchmark(build_query_sequences_flagfilter_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len(result) == 3124


def test_build_query_sequences_from_bam_directflagfilter_with_pysam(benchmark):
    result = benchmark(build_query_sequences_flagfilter_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len(result) == 3124


@pytest.mark.aligned_pairs
def test_build_aligned_pairs_default_with_pysam(benchmark):
    result = benchmark(build_aligned_pairs_with_pysam,
                       os.path.join(BAM_DATADIR, "with_md.bam"))
    assert len(result) == 3235


@pytest.mark.aligned_pairs
def test_build_aligned_pairs_matchesonly_with_pysam(benchmark):
    result = benchmark(build_aligned_pairs_with_pysam,
                       os.path.join(BAM_DATADIR, "with_md.bam"),
                       matches_only=True)
    assert len(result) == 3235


@pytest.mark.aligned_pairs    
def test_build_aligned_pairs_withseq_with_pysam(benchmark):
    result = benchmark(build_aligned_pairs_with_pysam,
                       os.path.join(BAM_DATADIR, "with_md.bam"),
                       with_seq=True)
    assert len(result) == 3235
    

