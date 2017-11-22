"""Benchmarking module for AlignedSegment functionality"""
import os
import array
import pysam


from TestUtils import BAM_DATADIR


def set_binary_tag():
    read = pysam.AlignedSegment()
    read.set_tag('FZ', array.array('H', range(1000)))
    return len(read.get_tag('FZ'))


def read_binary_tag(fn):
    with pysam.AlignmentFile(fn) as inf:
        read = next(inf.fetch())
    return len(read.get_tag('FZ'))


def test_set_binary_tag(benchmark):
    result = benchmark(set_binary_tag)
    assert result == 1000


def test_read_binary_tag(benchmark):
    result = benchmark(read_binary_tag, os.path.join(
        BAM_DATADIR, "example_btag.bam"))
    assert result == 260
