"""Benchmarking the cfaidx module. Usage::

pytest benchmark/faidx_bench.py
"""

import pysam


from TestUtils import BAM_DATADIR


def iterate_over_fastx(fn, persist=True):
    return len(list(pysam.FastxFile(fn, persist=persist)))


def iterate_over_fastx_as_file(fn):
    with open(fn) as inf:
        return len(inf.read())
    

def test_fasta_iteration_short_sequences(benchmark):
    result = benchmark(iterate_over_fastx, os.path.join(BAM_DATADIR, "faidx_ex1.fa"))
    assert result == 3270


def test_fasta_iteration_long_sequences(benchmark):
    result = benchmark(iterate_over_fastx, os.path.join(BAM_DATADIR, "ex1.fa"))
    assert result == 2


def test_fasta_iteration_short_sequences_without_persistence(benchmark):
    result = benchmark(iterate_over_fastx, os.path.join(BAM_DATADIR, "faidx_ex1.fa"), persist=False)
    assert result == 3270


def test_fasta_iteration_long_sequences_without_persistence(benchmark):
    result = benchmark(iterate_over_fastx, os.path.join(BAM_DATADIR, "ex1.fa"), persist=False)
    assert result == 2


def test_fasta_iteration_short_sequences_as_file(benchmark):
    result = benchmark(iterate_over_fastx_as_file, os.path.join(BAM_DATADIR, "faidx_ex1.fa"))
    assert result == 195399


def test_fasta_iteration_long_sequences_as_file(benchmark):
    result = benchmark(iterate_over_fastx_as_file, os.path.join(BAM_DATADIR, "ex1.fa"))
    assert result == 3225


def test_fastq_iteration_short_sequences(benchmark):
    result = benchmark(iterate_over_fastx, os.path.join(BAM_DATADIR, "faidx_ex1.fq"))
    assert result == 3270


def test_fastq_iteration_short_sequences_without_persistence(benchmark):
    result = benchmark(iterate_over_fastx, os.path.join(BAM_DATADIR, "faidx_ex1.fq"), persist=False)
    assert result == 3270
    

def test_fastq_iteration_short_sequences_as_file(benchmark):
    result = benchmark(iterate_over_fastx_as_file, os.path.join(BAM_DATADIR, "faidx_ex1.fq"))
    assert result == 320458
