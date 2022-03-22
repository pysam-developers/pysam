"""Benchmarking module for AlignmentFile functionality"""
from TestUtils import BAM_DATADIR, flatten_nested_list
from PileupTestUtils import *


def test_build_pileup_from_bam_with_samtoolsshell(benchmark):
    result = benchmark(build_pileup_with_samtoolsshell,
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


def test_build_depth_from_bam_with_samtoolsshell(benchmark):
    result = benchmark(build_depth_with_samtoolsshell,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 107241


def test_build_depth_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_depth_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert sum(result) == 107241


def test_build_depth_from_bam_with_pysam(benchmark):
    result = benchmark(build_depth_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    # different value, as samtools filters with a minimum
    # base quality of 13
    assert sum(result) == 110015


def test_build_depth_with_filter_from_bam_with_pysam(benchmark):
    result = benchmark(build_depth_with_filter_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert sum(result) == 107241


def test_build_query_bases_from_bam_with_samtoolsshell(benchmark):
    result = benchmark(build_query_bases_with_samtoolsshell,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 116308


def test_build_query_bases_from_bam_with_samtoolspysam(benchmark):
    result = benchmark(build_query_bases_with_samtoolspysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len("".join(flatten_nested_list(result))) == 116308
    

def test_build_query_bases_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_query_bases_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len("".join(flatten_nested_list(result))) == 116308


def test_build_query_bases_from_bam_with_pysam_pileups(benchmark):
    # note that there is no overlap detection here
    result = benchmark(build_query_bases_with_pysam_pileups,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len("".join(flatten_nested_list(result))) == 107241


def test_build_query_bases_from_bam_with_pysam(benchmark):
    result = benchmark(build_query_bases_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len("".join(flatten_nested_list(result))) == 116308


# note that pileups with/without reference sequence will differ due to
# realignment.
def test_build_query_bases_with_reference_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_query_bases_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"),
                       "-f", os.path.join(BAM_DATADIR, "ex1.fa"))
    assert len("".join(flatten_nested_list(result))) == 115924
    

def test_build_query_bases_with_reference_from_bam_with_pysam(benchmark):
    with pysam.FastaFile(os.path.join(BAM_DATADIR, "ex1.fa")) as fasta:
        result = benchmark(build_query_bases_with_pysam,
                           os.path.join(BAM_DATADIR, "ex2.bam"),
                           fastafile=fasta)
    assert len("".join(flatten_nested_list(result))) == 115924
    

def test_build_query_bases_with_reference_from_bam_with_samtoolspysam(benchmark):
    result = benchmark(build_query_bases_with_samtoolspysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"),
                       "-f", os.path.join(BAM_DATADIR, "ex1.fa"))
    assert len("".join(flatten_nested_list(result))) == 115924


def test_build_query_qualities_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_query_qualities_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len("".join(result)) == 107241


def test_build_query_qualities_from_bam_with_pysam(benchmark):
    result = benchmark(build_query_qualities_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert sum([len(x) for x in result]) == 107241


def test_build_query_names_from_bam_with_pysam(benchmark):
    result = benchmark(build_query_names_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len("".join([x for column in result for x in column])) == 2307343


def test_build_mapping_qualities_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_mapping_qualities_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len("".join(result)) == 107241


def test_build_mapping_qualities_from_bam_with_pysam(benchmark):
    result = benchmark(build_mapping_qualities_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert sum([len(x) for x in result]) == 107241


def test_build_query_positions_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_query_positions_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    # positions output by samtools are 1-based
    assert sum([sum(x) - len(x) for x in result]) == 1841699


def test_build_query_positions_from_bam_with_pysam(benchmark):
    result = benchmark(build_query_positions_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert sum([sum(x) for x in result]) == 1841699
