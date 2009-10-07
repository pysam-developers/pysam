import pysam

samfile = pysam.Samfile( "ex1.bam", "rb" )

print "###########################"
def my_fetch_callback( alignment ):
    print str(alignment)
samfile.fetch( "seq1:10-20", my_fetch_callback )

print "###########################"
def my_pileup_callback( pileups ):    print str(pileups)
samfile.pileup( "seq1:10-20", my_pileup_callback )

print "###########################"
iter = pysam.IteratorRow( samfile, "seq1:10-20")
for x in iter: print str(x)

print "###########################"
iter = pysam.IteratorColumn( samfile, "seq1:10-20")
for x in iter: print str(x)

print "###########################"
iter = pysam.IteratorRowAll( samfile )
for x in iter: print str(x)

class Counter:
    mCounts = 0
    def __call__(self, alignment):
        self.mCounts += 1

c = Counter()
samfile.fetch( "seq1:10-20", c )
print "counts=", c.mCounts

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
    samfile.fetch( "seq1:10-20", my_fetch_callback )
except AssertionError:
    print "caught fetch exception"

samfile.close()

print "###########################"
samfile = pysam.Samfile( "ex2.sam.gz", "r" )
def my_pileup_callback( pileups ):
    print str(pileups)
try:
    samfile.pileup( "seq1:10-20", my_pileup_callback )
except NotImplementedError:
    print "caught pileup exception"

