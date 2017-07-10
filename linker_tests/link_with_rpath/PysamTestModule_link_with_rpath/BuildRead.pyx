from __future__ import absolute_import

from pysam.libcalignedsegment cimport AlignedSegment
import pysam


cpdef build_read():
    cdef AlignedSegment read = pysam.AlignedSegment()
    read.query_name = "hello"
    read.query_sequence = "ACGT"
    return read
