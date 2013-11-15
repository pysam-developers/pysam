import gzip
import pysam
import timeit

def test_python():
    '''iterate through with python.'''
    f = gzip.open( "windows_small.bed.gz")
    l = len( [x.split("\t") for x in f])

def test_fetch_plain():
    """Stupid test function"""
    f = pysam.Tabixfile("windows_small.bed.gz")
    l = len( list(f.fetch()) )

def test_fetch_parsed():
    """Stupid test function"""
    f = pysam.Tabixfile("windows_small.bed.gz")
    l = len( list(f.fetch( parser = pysam.asBed())) )

def test_iterator_generic_compressed():
    f = gzip.open("windows_small.bed.gz")
    l = len( list( pysam.tabix_generic_iterator( f, parser = pysam.asBed() )))

def test_iterator_generic_uncompressed():
    f = open("windows_small.bed")
    l = len( list( pysam.tabix_generic_iterator( f, parser = pysam.asBed() )))

def test_iterator_parsed_compressed():
    f = gzip.open("windows_small.bed.gz")
    l = len( list( pysam.tabix_iterator( f, parser = pysam.asBed() )))

def test_iterator_parsed_uncompressed():
    f = open("windows_small.bed")
    l = len( list( pysam.tabix_iterator( f, parser = pysam.asBed() )))

def test_iterator_file_compressed():
    f = gzip.open("windows_small.bed")
    l = len( list( pysam.tabix_file_iterator( f, parser = pysam.asBed() )))

def test_iterator_file_uncompressed():
    f = open("windows_small.bed")
    l = len( list( pysam.tabix_file_iterator( f, parser = pysam.asBed() )))

iterations = 100

print "iterations=", iterations

tests = ( test_python, 
          test_fetch_plain, 
          test_fetch_parsed,
          test_iterator_generic_compressed, 
          test_iterator_generic_uncompressed, 
          test_iterator_parsed_compressed,
          test_iterator_parsed_uncompressed,
          test_iterator_file_compressed,
          test_iterator_file_uncompressed )

for test in tests:
    t = timeit.timeit( test, number = iterations )
    print "%5.2f\t%s" % (t,str(test))


