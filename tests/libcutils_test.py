"""Benchmarking the libcutils module. Usage::

pytest tests/libcutils_bench.py
"""
import pysam


def test_qualitystring_to_array_empty():
    result = pysam.array_to_qualitystring(pysam.qualitystring_to_array(""))
    assert result == ""
