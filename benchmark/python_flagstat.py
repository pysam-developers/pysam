"""compute number of reads/alignments from BAM file
===================================================

This is a benchmarking utility script with limited functionality.

Compute simple flag stats on a BAM-file using
the pysam python interface.
"""

import sys
import pysam

assert len(sys.argv) == 2, "USAGE: {} filename.bam".format(sys.argv[0])

is_paired = 0
is_proper = 0

for read in pysam.AlignmentFile(sys.argv[1], "rb"):
    is_paired += read.is_paired
    is_proper += read.is_proper_pair

print ("there are alignments of %i paired reads" % is_paired)
print ("there are %i proper paired alignments" % is_proper)
