## This script contains some example code 
## illustrating ways to to use the pysam 
## interface to samtools.
##
## The unit tests in the script pysam_test.py
## contain more examples.
##

import pysam

samfile = pysam.Samfile( "ex1.bam", "rb" )

print "###################"
# check different ways to iterate
print len(list(samfile.fetch()))
print len(list(samfile.fetch( "chr1", 10, 200 )))
print len(list(samfile.fetch( region="chr1:10-200" )))
print len(list(samfile.fetch( "chr1" )))
print len(list(samfile.fetch( region="chr1")))
print len(list(samfile.fetch( "chr2" )))
print len(list(samfile.fetch( region="chr2")))
print len(list(samfile.fetch()))
print len(list(samfile.fetch( "chr1" )))
print len(list(samfile.fetch( region="chr1")))
print len(list(samfile.fetch()))

print len(list(samfile.pileup( "chr1", 10, 200 )))
print len(list(samfile.pileup( region="chr1:10-200" )))
print len(list(samfile.pileup( "chr1" )))
print len(list(samfile.pileup( region="chr1")))
print len(list(samfile.pileup( "chr2" )))
print len(list(samfile.pileup( region="chr2")))
print len(list(samfile.pileup()))
print len(list(samfile.pileup()))

print "########### fetch with callback ################"
def my_fetch_callback( alignment ): print str(alignment)
samfile.fetch( region="chr1:10-200", callback=my_fetch_callback )

print "########## pileup with callback ################"
def my_pileup_callback( column ): print str(column)
samfile.pileup( region="chr1:10-200", callback=my_pileup_callback )


print "########### Using a callback object ###### "

class Counter:
    mCounts = 0
    def __call__(self, alignment):
        self.mCounts += 1

c = Counter()
samfile.fetch( region = "chr1:10-200", callback = c )
print "counts=", c.mCounts

print "########### Calling a samtools command line function ############"

for p in pysam.mpileup( "-c", "ex1.bam" ):
    print str(p)

print pysam.mpileup.getMessages()

print "########### Investigating headers #######################"

# playing arount with headers
samfile = pysam.Samfile( "ex3.sam", "r" )
print samfile.references
print samfile.lengths
print samfile.text
print samfile.header
header = samfile.header
samfile.close()

header["HD"]["SO"] = "unsorted"
outfile = pysam.Samfile( "out.sam", "wh", 
                         header = header )

outfile.close()

