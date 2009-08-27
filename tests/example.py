import pysam

samfile = pysam.Samfile()
samfile.open( "ex1.bam", "rb" )

def my_fetch_callback( alignment ):
    print str(alignment)
samfile.fetch( "seq1:10-20", my_fetch_callback )

def my_pileup_callback( pileups ):
    print str(pileups)
samfile.pileup( "seq1:10-20", my_pileup_callback )

print "###########################"
iter = pysam.IteratorRow( samfile, "seq1:10-20")
for x in iter: print str(x)

print "###########################"
iter = pysam.IteratorColumn( samfile, "seq1:10-20")
for x in iter: print str(x)

