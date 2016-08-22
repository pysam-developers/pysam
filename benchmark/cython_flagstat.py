"""compute number of reads/alignments from BAM file
===================================================

This is a benchmarking utility script with limited functionality.

Compute simple flag stats on a BAM-file using
the pysam cython interface.

"""

import sys
import pysam
import pyximport
pyximport.install()
import _cython_flagstat

assert len(sys.argv) == 2, "USAGE: {} filename.bam".format(sys.argv[0])

is_paired, is_proper = _cython_flagstat.count(
    pysam.AlignmentFile(sys.argv[1], "rb"))

print ("there are alignments of %i paired reads" % is_paired)
print ("there are %i proper paired alignments" % is_proper)
