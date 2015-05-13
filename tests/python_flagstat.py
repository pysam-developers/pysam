import pysam

is_paired = 0
is_proper = 0

for read in pysam.AlignmentFile("ex1.bam", "rb"):
    is_paired += read.is_paired
    is_proper += read.is_proper_pair

print ("there are alignments of %i paired reads" % is_paired)
print ("there are %i proper paired alignments" % is_proper)
