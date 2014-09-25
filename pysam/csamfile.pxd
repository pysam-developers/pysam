from calignmentfile cimport AlignedSegment, AlignmentFile

# Compatibility Layer for pysam < 0.8
cdef class AlignedRead(AlignedSegment):
    pass

cdef class Samfile(AlignmentFile):
    pass
