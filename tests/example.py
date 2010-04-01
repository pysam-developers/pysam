import sys
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

print "##########iterator row #################"
iter = pysam.IteratorRow( samfile, 0, 10, 200)
for x in iter: print str(x)

print "##########iterator col #################"
iter = pysam.IteratorColumn( samfile, 0, 10, 200 )
for x in iter: print str(x)

print "#########row all##################"
iter = pysam.IteratorRowAll( samfile )
for x in iter: print str(x)


print "###################"

class Counter:
    mCounts = 0
    def __call__(self, alignment):
        self.mCounts += 1

c = Counter()
samfile.fetch( "chr1:10-200", c )
print "counts=", c.mCounts

sys.exit(0)
print samfile.getTarget( 0 )
print samfile.getTarget( 1 )

for p in pysam.pileup( "-c", "ex1.bam" ):
    print str(p)

print pysam.pileup.getMessages()

for p in pysam.pileup( "-c", "ex1.bam", raw=True ):
    print str(p),



print "###########################"

samfile = pysam.Samfile( "ex2.sam.gz", "r" )

print "num targets=", samfile.getNumTargets()

iter = pysam.IteratorRowAll( samfile )
for x in iter: print str(x)

samfile.close()

print "###########################"
samfile = pysam.Samfile( "ex2.sam.gz", "r" )
def my_fetch_callback( alignment ):
    print str(alignment)

try:
    samfile.fetch( "chr1:10-20", my_fetch_callback )
except AssertionError:
    print "caught fetch exception"

samfile.close()

print "###########################"
samfile = pysam.Samfile( "ex2.sam.gz", "r" )
def my_pileup_callback( pileups ):
    print str(pileups)
try:
    samfile.pileup( "chr1:10-20", my_pileup_callback )
except NotImplementedError:
    print "caught pileup exception"

# playing arount with headers
samfile = pysam.Samfile( "ex3.sam", "r" )
print samfile.targets
print samfile.lengths
print samfile.text
print samdile.header
header = samfile.header
samfile.close()

header["HD"]["SO"] = "unsorted"
outfile = pysam.Samfile( "out.sam", "wh", 
                         header = header )

outfile.close()

