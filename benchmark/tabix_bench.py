import gzip
import os
import pysam

from TestUtils import TABIX_DATADIR

FN_COMPRESSED = "example.bed.gz"
FN_UNCOMPRESSED = "example.bed"
FN_LARGE_COMPRESSED = "example_large.bed.gz"
FN_LARGE_UNCOMPRESSED = "example_large.bed"


def read_python_compressed(fn):
    '''iterate through with python.'''
    with gzip.open(fn, mode="r") as f:
        return len([x.split(b"\t") for x in f])


def read_python_uncompressed(fn):
    with open(fn) as f:
        return len([x.split("\t") for x in f])


def fetch_plain(fn):
    with pysam.Tabixfile(fn) as f:
        return len(list(f.fetch()))


def fetch_parsed(fn):
    with pysam.Tabixfile(fn) as f:
        return len(list(f.fetch(parser=pysam.asBed())))


def iterate_generic_compressed(fn):
    with gzip.open(fn) as f:
        return len(list(pysam.tabix_generic_iterator(f, parser=pysam.asBed())))


def iterate_generic_uncompressed(fn):
    with open(fn) as f:
        return len(list(pysam.tabix_generic_iterator(f, parser=pysam.asBed())))


def iterate_parsed_compressed(fn):
    with gzip.open(fn) as f:
        return len(list(pysam.tabix_iterator(f, parser=pysam.asBed())))


def iterate_parsed_uncompressed(fn):
    with open(fn) as f:
        return len(list(pysam.tabix_iterator(f, parser=pysam.asBed())))


def iterate_file_compressed(fn):
    with gzip.open(fn) as f:
        return len(list(pysam.tabix_file_iterator(f, parser=pysam.asBed())))


def iterate_file_uncompressed(fn):
    with open(fn) as f:
        return len(list(pysam.tabix_file_iterator(f, parser=pysam.asBed())))


def test_read_python_compressed(benchmark):
    result = benchmark(read_python_compressed, os.path.join(TABIX_DATADIR, FN_COMPRESSED))
    assert result == 164


def test_read_python_uncompressed(benchmark):
    result = benchmark(read_python_uncompressed, os.path.join(TABIX_DATADIR, FN_UNCOMPRESSED))
    assert result == 164


def test_fetch_plain(benchmark):
    result = benchmark(fetch_plain, os.path.join(TABIX_DATADIR, FN_COMPRESSED))
    assert result == 164


def test_fetch_parsed(benchmark):
    result = benchmark(fetch_parsed, os.path.join(TABIX_DATADIR, FN_COMPRESSED))
    assert result == 164


def test_iterate_generic_compressed(benchmark):
    result = benchmark(iterate_generic_compressed, os.path.join(TABIX_DATADIR, FN_COMPRESSED))
    assert result == 164


def test_iterate_generic_uncompressed(benchmark):
    result = benchmark(iterate_generic_uncompressed, os.path.join(TABIX_DATADIR, FN_UNCOMPRESSED))
    assert result == 164


def test_iterate_parsed_compressed(benchmark):
    result = benchmark(iterate_parsed_compressed, os.path.join(TABIX_DATADIR, FN_COMPRESSED))
    assert result == 164


def test_iterate_parsed_uncompressed(benchmark):
    result = benchmark(iterate_parsed_uncompressed, os.path.join(TABIX_DATADIR, FN_UNCOMPRESSED))
    assert result == 164


def test_iterate_file_compressed(benchmark):
    result = benchmark(iterate_file_compressed, os.path.join(TABIX_DATADIR, FN_COMPRESSED))
    assert result == 164


def test_iterate_file_uncompressed(benchmark):
    result = benchmark(iterate_file_uncompressed, os.path.join(TABIX_DATADIR, FN_UNCOMPRESSED))
    assert result == 164


def test_read_python_large_compressed(benchmark):
    result = benchmark(read_python_compressed, os.path.join(TABIX_DATADIR, FN_LARGE_COMPRESSED))
    assert result == 100000


def test_read_python_large_uncompressed(benchmark):
    result = benchmark(read_python_uncompressed, os.path.join(TABIX_DATADIR, FN_LARGE_UNCOMPRESSED))
    assert result == 100000


def test_fetch_plain(benchmark):
    result = benchmark(fetch_plain, os.path.join(TABIX_DATADIR, FN_LARGE_COMPRESSED))
    assert result == 100000


def test_fetch_parsed(benchmark):
    result = benchmark(fetch_parsed, os.path.join(TABIX_DATADIR, FN_LARGE_COMPRESSED))
    assert result == 100000


def test_iterate_generic_large_compressed(benchmark):
    result = benchmark(iterate_generic_compressed, os.path.join(TABIX_DATADIR, FN_LARGE_COMPRESSED))
    assert result == 100000


def test_iterate_generic_large_uncompressed(benchmark):
    result = benchmark(iterate_generic_uncompressed, os.path.join(TABIX_DATADIR, FN_LARGE_UNCOMPRESSED))
    assert result == 100000


def test_iterate_parsed_large_compressed(benchmark):
    result = benchmark(iterate_parsed_compressed, os.path.join(TABIX_DATADIR, FN_LARGE_COMPRESSED))
    assert result == 100000


def test_iterate_parsed_large_uncompressed(benchmark):
    result = benchmark(iterate_parsed_uncompressed, os.path.join(TABIX_DATADIR, FN_LARGE_UNCOMPRESSED))
    assert result == 100000


def test_iterate_file_large_compressed(benchmark):
    result = benchmark(iterate_file_compressed, os.path.join(TABIX_DATADIR, FN_LARGE_COMPRESSED))
    assert result == 100000


def test_iterate_file_large_uncompressed(benchmark):
    result = benchmark(iterate_file_uncompressed, os.path.join(TABIX_DATADIR, FN_LARGE_UNCOMPRESSED))
    assert result == 100000
