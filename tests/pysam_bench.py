import pysam
import timeit 

iterations = 10

print "Samfile.fetch"
s = '''
f = pysam.Samfile( "ex1.bam", "rb" )
results = list(f.fetch())
'''
print timeit.repeat( s, number = iterations, setup="from __main__ import pysam" )

print "pysam.view"
s = '''
f = pysam.view( "ex1.bam" )
'''
print timeit.repeat( s, number = iterations, setup="from __main__ import pysam" )

print "Samfile.pileup"
s = '''
f = pysam.Samfile( "ex1.bam", "rb" )
results = list(f.pileup())
'''
print timeit.repeat( s, number = iterations, setup="from __main__ import pysam" )

print "Samfile.pileup with coverage retrieval"
s = '''
f = pysam.Samfile( "ex1.bam", "rb" )
results = [ x.n for x in f.pileup() ]
'''
print timeit.repeat( s, number = iterations, setup="from __main__ import pysam" )


print "Samfile.pileup with full retrieval"
s = '''
f = pysam.Samfile( "ex1.bam", "rb" )
results = [ x.pileups for x in f.pileup() ]
'''
print timeit.repeat( s, number = iterations, setup="from __main__ import pysam" )

print "pysam.pileup"
s = '''
f = pysam.view( "ex1.bam" )
'''
print timeit.repeat( s, number = iterations, setup="from __main__ import pysam" )


