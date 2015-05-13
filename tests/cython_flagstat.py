import pysam

import pyximport
pyximport.install()
import _cython_flagstat

is_paired, is_proper = _cython_flagstat.count(
    pysam.AlignmentFile("ex1.bam", "rb"))

print ("there are alignments of %i paired reads" % is_paired)
print ("there are %i proper paired alignments" % is_proper)
