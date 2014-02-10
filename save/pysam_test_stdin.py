# usage: cat pysam_ex1.bam | python pysam_test_stdin.pyx

import pysam

samfile = pysam.Samfile( "-", "rb" )

# set up the modifying iterators
l = list(samfile.fetch( until_eof = True ))
assert len(l) == 3270

