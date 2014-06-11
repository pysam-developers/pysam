from pysam.chtslib cimport Samfile, AlignedRead
from pysam.ctabix cimport Tabixfile

cdef Samfile samfile
cdef Tabixfile tabixfile


def testCountBAM( Samfile samfile ):
    '''test reading from a BAM file accessing
    the flag field directly.'''

    cdef AlignedRead read
    cdef int n = 0
    
    for read in samfile.fetch():
        flag = read._delegate.core.flag
        n += 1
            
    return n

def testCountGTF( Tabixfile tabixfile ):
    '''test reading from a tabixfile.'''
    
    cdef int n = 0

    for entry in tabixfile.fetch():
        n += 1

    return n
