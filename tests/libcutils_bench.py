"""Benchmarking the libcutils module. Usage::

pytest tests/libcutils_bench.py
"""
import pysam


def test_qualitystring_to_array_long_sequences(benchmark):
    result = benchmark(pysam.array_to_qualitystring, pysam.qualitystring_to_array("123") * 500)
    assert result == "123" * 500
