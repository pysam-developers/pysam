# cython: language_level=3

from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from pysam.libctabix cimport Tabixfile

cdef AlignmentFile samfile
cdef Tabixfile tabixfile


def testCountBAM(AlignmentFile samfile):
    '''test reading from a BAM file accessing
    the flag field directly.'''

    cdef AlignedSegment read
    cdef int n = 0
    
    for read in samfile.fetch():
        flag = read._delegate.core.flag
        n += 1
            
    return n

def testCountGTF(Tabixfile tabixfile):
    '''test reading from a tabixfile.'''
    
    cdef int n = 0

    for entry in tabixfile.fetch():
        n += 1

    return n
