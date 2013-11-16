import gzip
import pysam
import timeit

iterations = 5
repeats = 100
print "repeats=", repeats, "iterations=", iterations

fn_compressed = '/tmp/windows_small.bed.gz'
fn_uncompressed = '/tmp/windows_small.bed'

def test_python_compressed():
    '''iterate through with python.'''
    f = gzip.open( fn_compressed)
    l = len( [x.split("\t") for x in f])
 
def test_python_uncompressed():
    '''iterate through with python.'''
    f = open( "windows_small.bed")
    l = len( [x.split("\t") for x in f])
 
def test_fetch_plain():
    """Stupid test function"""
    f = pysam.Tabixfile(fn_compressed)
    l = len( list(f.fetch()) )

def test_fetch_parsed():
    """Stupid test function"""
    f = pysam.Tabixfile(fn_compressed)
    l = len( list(f.fetch( parser = pysam.asBed())) )

def test_iterator_generic_compressed():
    f = gzip.open(fn_compressed)
    l = len( list( pysam.tabix_generic_iterator( f, parser = pysam.asBed() )))

def test_iterator_generic_uncompressed():
    f = open("windows_small.bed")
    l = len( list( pysam.tabix_generic_iterator( f, parser = pysam.asBed() )))

def test_iterator_parsed_compressed():
    f = gzip.open(fn_compressed)
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

tests = ( test_python_compressed, 
          test_python_uncompressed, 
          test_fetch_plain, 
          test_fetch_parsed,
          test_iterator_generic_compressed, 
          test_iterator_generic_uncompressed, 
          test_iterator_parsed_compressed,
          test_iterator_parsed_uncompressed,
          test_iterator_file_compressed,
          test_iterator_file_uncompressed )

for repeat in range( repeats ):
    print "# repeat=", repeat
    for test in tests:
        try:
            t = timeit.timeit( test, number = iterations )
        except AttributeError:
            continue
        print "%5.2f\t%s" % (t,str(test))


