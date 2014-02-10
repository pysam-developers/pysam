'''benchmark pysam BAM/SAM access with the samtools commandline tools.

samtools functions are called via the pysam interface to avoid the over-head
of starting additional processes.
'''

import pysam
import timeit 

iterations = 10

def runBenchmark( test, 
                  pysam_way,
                  samtools_way = None):
    print test
    print timeit.repeat( pysam_way, number = iterations, setup="from __main__ import pysam" )
    if samtools_way:
        print timeit.repeat( samtools_way, number = iterations, setup="from __main__ import pysam" )

runBenchmark( "Samfile.fetch",
'''
f = pysam.Samfile( "ex1.bam", "rb" )
results = list(f.fetch())
''',
'''
f = pysam.view( "ex1.bam" )
'''
)

runBenchmark( "Samfile.pileup",
'''
f = pysam.Samfile( "ex1.bam", "rb" )
results = list(f.pileup())
''',
'''
f = pysam.pileup( "ex1.bam" )
''')

runBenchmark( "Samfile.pileup with coverage retrieval",
'''
f = pysam.Samfile( "ex1.bam", "rb" )
results = [ x.n for x in f.pileup() ]
''' )

runBenchmark( "Samfile.pileup with full retrieval",
'''
f = pysam.Samfile( "ex1.bam", "rb" )
results = [ x.pileups for x in f.pileup() ]
''' )

runBenchmark( "Samfile.pileup - many references",
'''
f = pysam.Samfile( "manyrefs.bam", "rb" )
results = [ x.pileups for x in f.pileup() ]
''',
'''
f = pysam.pileup( "manyrefs.bam" )
'''
 )




