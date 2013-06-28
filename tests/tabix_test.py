#!/usr/bin/env python
'''unit testing code for pysam.

Execute in the :file:`tests` directory as it requires the Makefile
and data files located there.
'''

import sys, os, shutil, gzip
import pysam
import unittest
import itertools
import subprocess
import glob
import re

IS_PYTHON3 = sys.version_info[0] >= 3

def myzip_open( infile, mode = "r" ):
    '''open compressed file and decode.'''

    def _convert(f):
        for l in f:
            yield l.decode("ascii")

    if IS_PYTHON3:
        if mode == "r":
            return _convert(gzip.open(infile,"r"))
    else:
        return gzip.open( mode )

def loadAndConvert( infile ):
    '''load and convert all fields to bytes'''
    data = []
    if infile.endswith(".gz"):
        for line in gzip.open( infile ):
            line = line.decode("ascii")
            if line.startswith("#"): continue
            d = line.strip().split("\t")
            data.append( [x.encode("ascii") for x in d ] )
    else:
        with open(infile) as f:
            for line in f:
                if line.startswith("#"): continue
                d = line.strip().split("\t")
                data.append( [x.encode("ascii") for x in d ] )

    return data

def splitToBytes( s ):
    '''split string and return list of bytes.'''
    return [x.encode("ascii") for x in s.split("\t")]

def checkBinaryEqual( filename1, filename2 ):
    '''return true if the two files are binary equal.'''
    if os.path.getsize( filename1 ) !=  os.path.getsize( filename2 ):
        return False

    infile1 = open(filename1, "rb")
    infile2 = open(filename2, "rb")

    d1, d2 = infile1.read(), infile2.read()
    found = False
    for c1,c2 in zip( d1, d2 ):
        if c1 != c2: break
    else:
        found = True

    infile1.close()
    infile2.close()
    return found

class TestIndexing(unittest.TestCase):
    filename = "example.gtf.gz" 
    filename_idx = "example.gtf.gz.tbi" 

    def setUp( self ):
        
        self.tmpfilename = "tmp_%i.gtf.gz" % id(self)
        shutil.copyfile( self.filename, self.tmpfilename )

    def testIndexPreset( self ):
        '''test indexing via preset.'''

        pysam.tabix_index( self.tmpfilename, preset = "gff" )
        checkBinaryEqual( self.tmpfilename + ".tbi", self.filename_idx )

    def tearDown( self ):
        os.unlink( self.tmpfilename )
        os.unlink( self.tmpfilename + ".tbi" )

class TestCompression(unittest.TestCase):
    filename = "example.gtf.gz" 
    filename_idx = "example.gtf.gz.tbi" 

    def setUp( self ):
        
        self.tmpfilename = "tmp_%i.gtf" % id(self)
        infile = gzip.open( self.filename, "rb")
        outfile = open( self.tmpfilename, "wb" )
        outfile.write( infile.read() )
        outfile.close()
        infile.close()

    def testIndexPreset( self ):
        '''test indexing via preset.'''
        
        pysam.tabix_index( self.tmpfilename, preset = "gff" )
        checkBinaryEqual( self.tmpfilename + ".gz", self.filename )
        checkBinaryEqual( self.tmpfilename + ".gz.tbi", self.filename_idx )

    def testCompression( self ):
        '''see also issue 106'''
        pysam.tabix_compress( self.tmpfilename, self.tmpfilename + ".gz" )
        checkBinaryEqual( self.tmpfilename, self.tmpfilename + ".gz" )
        
    def tearDown( self ):
        os.unlink( self.tmpfilename + ".gz" )
        if os.path.exists( self.tmpfilename + ".gz.tbi" ):
            os.unlink( self.tmpfilename + ".gz.tbi" )

class TestIteration( unittest.TestCase ):

    filename = "example.gtf.gz" 

    def setUp( self ):

        self.tabix = pysam.Tabixfile( self.filename )
        lines = []
        inf = gzip.open( self.filename, "rb")
        for line in inf:
            line = line.decode('ascii')
            if line.startswith("#"): continue
            lines.append( line )
        inf.close()
        # creates index of contig, start, end, adds content without newline.
        self.compare = [ 
            (x[0][0], int(x[0][3]), int(x[0][4]), x[1])
            for x in [ (y.split("\t"), y[:-1]) for y in lines ] ]
                         
    def getSubset( self, contig = None, start = None, end = None):
        
        if contig == None:
            # all lines
            subset = [ x[3] for x in self.compare ]
        else:
            if start != None and end == None:
                # until end of contig
                subset = [ x[3] for x in self.compare if x[0] == contig and x[2] > start ]
            elif start == None and end != None:
                # from start of contig
                subset = [ x[3] for x in self.compare if x[0] == contig and x[1] <= end ]
            elif start == None and end == None:
                subset = [ x[3] for x in self.compare if x[0] == contig ]
            else:
                # all within interval
                subset = [ x[3] for x in self.compare if x[0] == contig and \
                               min( x[2], end) - max(x[1], start) > 0 ]
            
        return subset

    def checkPairwise( self, result, ref ):
        '''check pairwise results.
        '''
        result.sort()
        ref.sort()

        a = set(result)
        b = set(ref)

        self.assertEqual( len(result), len(ref),
                          "unexpected number of results: result=%i, expected ref=%i, differences are %s: %s" \
                              % (len(result), len(ref),
                                 a.difference(b), 
                                 b.difference(a) ))

        for x, d in enumerate( list(zip( result, ref ))):
            self.assertEqual( d[0], d[1],
                              "unexpected results in pair %i:\n'%s', expected\n'%s'" % \
                                  (x, 
                                   d[0], 
                                   d[1]) )


    def testAll( self ):
        result = list(self.tabix.fetch())
        ref = self.getSubset( )
        self.checkPairwise( result, ref )

    def testPerContig( self ):
        for contig in ("chr1", "chr2", "chr1", "chr2" ):
            result = list(self.tabix.fetch( contig ))
            ref = self.getSubset( contig )
            self.checkPairwise( result, ref )
            
    def testPerContigToEnd( self ):
        
        end = None
        for contig in ("chr1", "chr2", "chr1", "chr2" ):
            for start in range( 0, 200000, 1000):
                result = list(self.tabix.fetch( contig, start, end ))
                ref = self.getSubset( contig, start, end )
                self.checkPairwise( result, ref )

    def testPerContigFromStart( self ):
        
        start = None
        for contig in ("chr1", "chr2", "chr1", "chr2" ):
            for end in range( 0, 200000, 1000):
                result = list(self.tabix.fetch( contig, start, end ))
                ref = self.getSubset( contig, start, end )
                self.checkPairwise( result, ref )

    def testPerContig( self ):
        
        start, end  = None, None
        for contig in ("chr1", "chr2", "chr1", "chr2" ):
            result = list(self.tabix.fetch( contig, start, end ))
            ref = self.getSubset( contig, start, end )
            self.checkPairwise( result, ref )
                
    def testPerInterval( self ):
        
        start, end  = None, None
        for contig in ("chr1", "chr2", "chr1", "chr2" ):
            for start in range( 0, 200000, 2000):
                for end in range( start, start + 2000, 500):
                    result = list(self.tabix.fetch( contig, start, end ))
                    ref = self.getSubset( contig, start, end )
                    self.checkPairwise( result, ref )
                

    def testInvalidIntervals( self ):
        
        self.assertRaises( ValueError, self.tabix.fetch, "chr1", 0, -10)
        self.assertRaises( ValueError, self.tabix.fetch, "chr1", -10, 200)
        self.assertRaises( ValueError, self.tabix.fetch, "chr1", 200, 0)
        self.assertRaises( ValueError, self.tabix.fetch, "chr1", -10, -20)
        self.assertRaises( ValueError, self.tabix.fetch, "chrUn" )

    def testGetContigs( self ):
        self.assertEqual( sorted(self.tabix.contigs), [b"chr1", b"chr2"] )
        # check that contigs is read-only
        self.assertRaises( AttributeError, setattr, self.tabix, "contigs", ["chr1", "chr2"] )

    def testHeader( self ):
        ref = []
        inf = gzip.open( self.filename )
        for x in inf:
            x = x.decode("ascii")
            if not x.startswith("#"): break
            ref.append( x[:-1].encode('ascii') )
        inf.close()

        header = list( self.tabix.header )
        self.assertEqual( ref, header )

    def testReopening( self ):
        '''test repeated opening of the same file.'''
        def func1():
            # opens any tabix file
            inf = pysam.Tabixfile(self.filename)
            return

        for i in range(10000):
            func1()

class TestParser( unittest.TestCase ):

    filename = "example.gtf.gz" 

    def setUp( self ):

        self.tabix = pysam.Tabixfile( self.filename )
        self.compare = loadAndConvert( self.filename )

    def testRead( self ):

        for x, r in enumerate(self.tabix.fetch( parser = pysam.asTuple() )):
            self.assertEqual( self.compare[x], list(r) )
            self.assertEqual( len(self.compare[x]), len(r) )

            # test indexing
            for c in range(0,len(r)):
                self.assertEqual( self.compare[x][c], r[c] )

            # test slicing access
            for c in range(0, len(r)-1):
                for cc in range(c+1, len(r)):
                    self.assertEqual( self.compare[x][c:cc],
                                      r[c:cc] )

    def testWrite( self ):
        
        for x, r in enumerate(self.tabix.fetch( parser = pysam.asTuple() )):
            self.assertEqual( self.compare[x], list(r) )
            c = list(r)
            for y in range(len(r)):
                r[y] = "test_%05i" % y
                c[y] = "test_%05i" % y
            self.assertEqual( [x.encode("ascii") for x in c], list(r) )
            self.assertEqual( "\t".join( c ), str(r) )
            # check second assignment
            for y in range(len(r)):
                r[y] = "test_%05i" % y
            self.assertEqual( [x.encode("ascii") for x in c], list(r) )
            self.assertEqual( "\t".join( c ), str(r) )

    def testUnset( self ):
        for x, r in enumerate(self.tabix.fetch( parser = pysam.asTuple() )):
            self.assertEqual( self.compare[x], list(r) )
            c = list(r)
            e = [ x.decode('ascii') for x in r ]
            for y in range(len(r)):
                r[y] = None
                c[y] = None
                e[y] = ""
                self.assertEqual( c, list(r) )
                self.assertEqual( "\t".join(e), str(r) )

    def testIteratorCompressed( self ):
        '''test iteration from compressed file.'''
        with gzip.open( self.filename ) as infile:
            for x, r in enumerate(pysam.tabix_iterator( infile, pysam.asTuple() )):
                self.assertEqual( self.compare[x], list(r) )
                self.assertEqual( len(self.compare[x]), len(r) )

                # test indexing
                for c in range(0,len(r)):
                    self.assertEqual( self.compare[x][c], r[c] )

                # test slicing access
                for c in range(0, len(r)-1):
                    for cc in range(c+1, len(r)):
                        self.assertEqual( self.compare[x][c:cc],
                                          r[c:cc] )

    def testIteratorUncompressed( self ):
        '''test iteration from uncompressed file.'''
        tmpfilename = 'tmp_testIteratorUncompressed'
        infile = gzip.open( self.filename, "rb")
        outfile = open( tmpfilename, "wb" )
        outfile.write( infile.read() )
        outfile.close()
        infile.close()

        with open( tmpfilename ) as infile:
            for x, r in enumerate(pysam.tabix_iterator( infile, pysam.asTuple() )):
                self.assertEqual( self.compare[x], list(r) )
                self.assertEqual( len(self.compare[x]), len(r) )

                # test indexing
                for c in range(0,len(r)):
                    self.assertEqual( self.compare[x][c], r[c] )

                # test slicing access
                for c in range(0, len(r)-1):
                    for cc in range(c+1, len(r)):
                        self.assertEqual( self.compare[x][c:cc],
                                          r[c:cc] )

        os.unlink( tmpfilename )

class TestGTF( TestParser ):

    def testRead( self ):

        for x, r in enumerate(self.tabix.fetch( parser = pysam.asGTF() )):
            c = self.compare[x]
            self.assertEqual( len(c), len(r) )
            self.assertEqual( list(c), list(r) )
            self.assertEqual( c, splitToBytes( str(r) ) )
            self.assertTrue( r.gene_id.startswith("ENSG") )
            if r.feature != b'gene':
                self.assertTrue( r.transcript_id.startswith("ENST") )
            self.assertEqual( c[0], r.contig )

class TestBed( unittest.TestCase ):
    filename = "example.bed.gz"

    def setUp( self ):

        self.tabix = pysam.Tabixfile( self.filename)
        self.compare = loadAndConvert( self.filename )

    def testRead( self ):

        for x, r in enumerate(self.tabix.fetch( parser = pysam.asBed() )):
            c = self.compare[x]
            self.assertEqual( len(c), len(r) )
            self.assertEqual( c, splitToBytes( str(r) ) )
            self.assertEqual( list(c), list(r) )
            self.assertEqual( c[0], r.contig)
            self.assertEqual( int(c[1]), r.start)
            self.assertEqual( int(c[2]), r.end)

    def testWrite( self ):

        for x, r in enumerate(self.tabix.fetch( parser = pysam.asBed() )):
            c = self.compare[x]
            self.assertEqual( c, splitToBytes(str(r) ))
            self.assertEqual( list(c), list(r) )

            r.contig = "test"
            self.assertEqual( b"test", r.contig)
            self.assertEqual( b"test", r[0])

            r.start += 1
            self.assertEqual( int(c[1]) + 1, r.start )
            self.assertEqual( str(int(c[1]) + 1), r[1].decode("ascii" ))

            r.end += 1
            self.assertEqual( int(c[2]) + 1, r.end )
            self.assertEqual( str(int(c[2]) + 1), r[2].decode("ascii") )

class TestVCF( unittest.TestCase ):

    filename = "example.vcf40"

    def setUp( self ):
        
        self.tmpfilename = "tmp_%s.vcf" % id(self)
        shutil.copyfile( self.filename, self.tmpfilename )
        pysam.tabix_index( self.tmpfilename, preset = "vcf" )

    def tearDown( self ):
        os.unlink( self.tmpfilename + ".gz" )
        if os.path.exists( self.tmpfilename + ".gz.tbi" ):
            os.unlink( self.tmpfilename + ".gz.tbi" )

class TestVCFFromTabix( TestVCF ):

    columns = ("contig", "pos", "id", 
               "ref", "alt", "qual", 
               "filter", "info", "format" )


    def setUp( self ):

        TestVCF.setUp( self )

        self.tabix = pysam.Tabixfile( self.tmpfilename + ".gz" )
        self.compare = loadAndConvert( self.filename )

    def testRead( self ):
        
        ncolumns = len(self.columns) 

        for x, r in enumerate(self.tabix.fetch( parser = pysam.asVCF() )):
            c = self.compare[x]
            for y, field in enumerate( self.columns ):
                # it is ok to have a missing format column
                if y == 8 and y == len(c): continue
                if field == "pos":
                    self.assertEqual( int(c[y]) - 1, getattr( r, field ) )
                    self.assertEqual( int(c[y]) - 1, r.pos )
                else:
                    self.assertEqual( c[y], getattr( r, field ), 
                                      "mismatch in field %s: %s != %s" %\
                                          ( field,c[y], getattr( r, field ) ) )
            if len(c) == 8:
                self.assertEqual( 0, len(r) )
            else:
                self.assertEqual( len(c), len( r ) + ncolumns )

            for y in range(len(c) - ncolumns):
                self.assertEqual( c[ncolumns+y], r[y] )
   
    def testWrite( self ):

        ncolumns = len(self.columns) 

        for x, r in enumerate(self.tabix.fetch( parser = pysam.asVCF() )):
            c = self.compare[x]
            # check unmodified string
            cmp_string = str(r)
            ref_string = "\t".join( [x.decode() for x in c] )

            self.assertEqual( ref_string, cmp_string )
            
            # set fields and compare field-wise
            for y, field in enumerate( self.columns ):
                # it is ok to have a missing format column
                if y == 8 and y == len(c): continue
                if field == "pos":
                    rpos = getattr( r, field )
                    self.assertEqual( int(c[y]) - 1, rpos )
                    self.assertEqual( int(c[y]) - 1, r.pos )
                    # increment pos by 1
                    setattr( r, field, rpos + 1 )
                    self.assertEqual( getattr( r, field ), rpos + 1 )
                    c[y] = str(int(c[y]) + 1 ) 
                else:
                    setattr( r, field, "test_%i" % y)
                    c[y] = ("test_%i" % y).encode('ascii')
                    self.assertEqual( c[y], getattr( r, field ), 
                                      "mismatch in field %s: %s != %s" %\
                                          ( field,c[y], getattr( r, field ) ) )

            if len(c) == 8:
                self.assertEqual( 0, len(r) )
            else:
                self.assertEqual( len(c), len( r ) + ncolumns )
            
            for y in range(len(c) - ncolumns):
                c[ncolumns+y] = ("test_%i" % y).encode('ascii')
                r[y] = ("test_%i" % y).encode('ascii')
                self.assertEqual( c[ncolumns+y], r[y] )

class TestVCFFromVCF( TestVCF ):

    columns = ("chrom", "pos", "id", 
               "ref", "alt", "qual", 
               "filter", "info", "format" )

    # tests failing while parsing
    fail_on_parsing = ( (5, "Flag fields should not have a value"),
                       (9, "aouao" ),
                       (13, "aoeu" ),
                       (18, "Error BAD_NUMBER_OF_PARAMETERS" ),
                       (24, "Error HEADING_NOT_SEPARATED_BY_TABS" ) )

    # tests failing on opening
    fail_on_opening = ( (24, "Error HEADING_NOT_SEPARATED_BY_TABS" ),
                     )
    def setUp( self ):
        
        TestVCF.setUp( self )

        self.vcf = pysam.VCF()
        self.compare = loadAndConvert( self.filename )

    def testParsing( self ):

        # self.vcf.connect( self.tmpfilename + ".gz" )
        ncolumns = len(self.columns) 

        fn = os.path.basename( self.filename )
        with open(self.filename) as f:

            for x, msg in self.fail_on_opening:
                if "%i.vcf" % x == fn:
                    self.assertRaises( ValueError, self.vcf.parse, f )
                    return
            else:
                iter = self.vcf.parse(f)

            for x, msg in self.fail_on_parsing:
                if "%i.vcf" % x == fn:
                    self.assertRaises( ValueError, list, iter )
                    break
                    # python 2.7
                    # self.assertRaisesRegexp( ValueError, re.compile(msg), self.vcf.parse, f )
            else:
                # do the actual parsing
                for x, r in enumerate(iter):
                    c = self.compare[x]
                    for y, field in enumerate( self.columns ):
                        # it is ok to have a missing format column
                        if y == 8 and y == len(c): continue

                        val = r[field] 
                        if field == "pos":
                            self.assertEqual( int(c[y]) - 1, val )
                        elif field == "alt":
                            if c[y] == ".":
                                # convert . to empty list
                                self.assertEqual( [], val, 
                                                  "mismatch in field %s: expected %s, got %s" %\
                                                      ( field,c[y], val ) )
                            else:
                                # convert to list
                                self.assertEqual( c[y].split(","), val, 
                                                  "mismatch in field %s: expected %s, got %s" %\
                                                      ( field,c[y], val ) )

                        elif field == "filter":
                            if c[y] == "PASS" or c[y] == ".":
                                # convert PASS to empty list
                                self.assertEqual( [], val, 
                                                  "mismatch in field %s: expected %s, got %s" %\
                                                      ( field,c[y], val ) )
                            else:
                                # convert to list
                                self.assertEqual( c[y].split(";"), val, 
                                                  "mismatch in field %s: expected %s, got %s" %\
                                                      ( field,c[y], val ) )

                        elif field == "info":
                            # tests for info field not implemented
                            pass
                        elif field == "qual" and c[y] == ".":
                            self.assertEqual( -1, val,
                                               "mismatch in field %s: expected %s, got %s" %\
                                                   ( field,c[y], val ) )
                        elif field == "format":
                            # format field converted to list
                            self.assertEqual( c[y].split(":"), val,
                                              "mismatch in field %s: expected %s, got %s" %\
                                                  ( field,c[y], val ) )

                        elif type(val) in (int, float):
                            if c[y] == ".":
                                self.assertEqual( None, val, 
                                                  "mismatch in field %s: expected %s, got %s" %\
                                                      ( field,c[y], val ) )

                            else:
                                self.assertEqual( float( c[y]), float(val), 
                                                  "mismatch in field %s: expected %s, got %s" %\
                                                      ( field,c[y], val ) )

                        else:
                            self.assertEqual( c[y], val, 
                                              "mismatch in field %s: expected %s, got %s" %\
                                                  ( field,c[y], val ) )

############################################################################ 
# create a test class for each example vcf file.
# Two samples are created - 
# 1. Testing pysam/tabix access
# 2. Testing the VCF class
vcf_files = glob.glob( "vcf-examples/*.vcf" )

for vcf_file in vcf_files:
    n = "VCFFromTabixTest_%s" % os.path.basename( vcf_file[:-4] )
    globals()[n] = type( n, (TestVCFFromTabix,), dict( filename=vcf_file,) )
    n = "VCFFromVCFTest_%s" % os.path.basename( vcf_file[:-4] )
    globals()[n] = type( n, (TestVCFFromVCF,), dict( filename=vcf_file,) )

############################################################################                   
class TestRemoteFileHTTP( unittest.TestCase):

    url = "http://genserv.anat.ox.ac.uk/downloads/pysam/test/example.gtf.gz"
    region = "chr1:1-1000"
    local = "example.gtf.gz"

    def testFetchAll( self ):
        remote_file = pysam.Tabixfile(self.url, "r")  
        remote_result = list(remote_file.fetch())
        local_file = pysam.Tabixfile(self.local, "r")  
        local_result = list(local_file.fetch())

        self.assertEqual( len(remote_result), len(local_result) )
        for x, y in zip(remote_result, local_result):
            self.assertEqual( x, y )


if __name__ == "__main__":

    unittest.main()


