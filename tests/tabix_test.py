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

def checkBinaryEqual( filename1, filename2 ):
    '''return true if the two files are binary equal.'''
    if os.path.getsize( filename1 ) !=  os.path.getsize( filename2 ):
        return False

    infile1 = open(filename1, "rb")
    infile2 = open(filename2, "rb")

    def chariter( infile ):
        while 1:
            c = infile.read(1)
            if c == "": break
            yield c

    found = False
    for c1,c2 in itertools.izip( chariter( infile1), chariter( infile2) ):
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
        infile = gzip.open( self.filename, "r")
        outfile = open( self.tmpfilename, "w" )
        outfile.write( "".join(infile.readlines()) )
        outfile.close()
        infile.close()

    def testIndexPreset( self ):
        '''test indexing via preset.'''

        pysam.tabix_index( self.tmpfilename, preset = "gff" )
        checkBinaryEqual( self.tmpfilename + ".gz", self.filename )
        checkBinaryEqual( self.tmpfilename + ".gz.tbi", self.filename_idx )

    def tearDown( self ):
        os.unlink( self.tmpfilename + ".gz" )
        os.unlink( self.tmpfilename + ".gz.tbi" )

class TestIteration( unittest.TestCase ):

    filename = "example.gtf.gz" 

    def setUp( self ):

        self.tabix = pysam.Tabixfile( self.filename )
        lines = [ x for x in gzip.open(self.filename).readlines() if not x.startswith("#") ]
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

        result.sort()
        ref.sort()

        a = set(result)
        b = set(ref)

        self.assertEqual( len(result), len(ref),
                          "unexpected number of results: %i, expected %i, differences are %s: %s" \
                              % (len(result), len(ref),
                                 a.difference(b), 
                                 b.difference(a) ))

        for x, d in enumerate( zip( result, ref )):
            self.assertEqual( d[0], d[1],
                              "unexpected results in pair %i: '%s', expected '%s'" % \
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
        self.assertEqual( sorted(self.tabix.contigs), ["chr1", "chr2"] )
        # check that contigs is read-only
        self.assertRaises( AttributeError, setattr, self.tabix, "contigs", ["chr1", "chr2"] )

    def testHeader( self ):
        ref = []
        for x in gzip.open( self.filename ):
            if not x.startswith("#"): break
            ref.append( x[:-1] )
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
        self.compare = [ x[:-1].split("\t") for x in gzip.open( self.filename, "r") if not x.startswith("#") ]

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
            self.assertEqual( c, list(r) )
            self.assertEqual( "\t".join( c ), str(r) )
            # check second assignment
            for y in range(len(r)):
                r[y] = "test_%05i" % y
            self.assertEqual( c, list(r) )
            self.assertEqual( "\t".join( c ), str(r) )

    def testUnset( self ):
        for x, r in enumerate(self.tabix.fetch( parser = pysam.asTuple() )):
            self.assertEqual( self.compare[x], list(r) )
            c = list(r)
            e = list(r)
            for y in range(len(r)):
                r[y] = c[y] = None
                e[y] = ""
                self.assertEqual( c, list(r) )
                self.assertEqual( "\t".join(e), str(r) )

class TestGTF( TestParser ):

    def testRead( self ):

        for x, r in enumerate(self.tabix.fetch( parser = pysam.asGTF() )):
            c =  self.compare[x]
            
            self.assertEqual( len(c), len(r) )
            self.assertEqual( "\t".join(c), str(r) )
            self.assertTrue( r.gene_id.startswith("ENSG") )
            if r.feature != "gene":
                self.assertTrue( r.transcript_id.startswith("ENST") )
            self.assertEqual( c[0], r.contig )

class TestBed( unittest.TestCase ):
    filename = "example.bed.gz"

    def setUp( self ):

        self.tabix = pysam.Tabixfile( self.filename )
        self.compare = [ x[:-1].split("\t") for x in gzip.open( self.filename, "r") if not x.startswith("#") ]

    def testRead( self ):

        for x, r in enumerate(self.tabix.fetch( parser = pysam.asBed() )):
            c = self.compare[x]
            self.assertEqual( "\t".join( c ), str(r) )
            self.assertEqual( list(c), list(r) )
            self.assertEqual( c[0], r.contig)
            self.assertEqual( int(c[1]), r.start)
            self.assertEqual( int(c[2]), r.end)

    def testWrite( self ):

        for x, r in enumerate(self.tabix.fetch( parser = pysam.asBed() )):
            c = self.compare[x]
            self.assertEqual( "\t".join( c ), str(r) )
            self.assertEqual( list(c), list(r) )

            r.contig = "test"
            self.assertEqual( "test", r.contig)
            self.assertEqual( "test", r[0])

            r.start += 1
            self.assertEqual( int(c[1]) + 1, r.start )
            self.assertEqual( str(int(c[1]) + 1), r[1] )

            r.end += 1
            self.assertEqual( int(c[2]) + 1, r.end )
            self.assertEqual( str(int(c[2]) + 1), r[2] )

class TestVCF( TestParser ):

    filename = "example.vcf40.gz"
    columns = ("contig", "pos", "id", 
               "ref", "alt", "qual", 
               "filter", "info", "format" )

    def testRead( self ):
        
        ncolumns = len(self.columns) 

        for x, r in enumerate(self.tabix.fetch( parser = pysam.asVCF() )):
            c = self.compare[x]
            for y, field in enumerate( self.columns ):
                if field == "pos":
                    self.assertEqual( int(c[y]) - 1, getattr( r, field ) )
                    self.assertEqual( int(c[y]) - 1, r.pos )
                else:
                    self.assertEqual( c[y], getattr( r, field ), 
                                      "mismatch in field %s: %s != %s" %\
                                          ( field,c[y], getattr( r, field ) ) )
            self.assertEqual( len(c), len( r ) + ncolumns )
            
            for y in range(len(c) - ncolumns):
                self.assertEqual( c[ncolumns+y], r[y] )
                
    def testWrite( self ):

        ncolumns = len(self.columns) 

        for x, r in enumerate(self.tabix.fetch( parser = pysam.asVCF() )):

            c = self.compare[x]

            # check unmodified string
            ref_string = "\t".join( c )
            cmp_string = str(r)
            self.assertEqual( ref_string, cmp_string )

            # set fields and compare field-wise
            for y, field in enumerate( self.columns ):
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
                    c[y] = "test_%i" % y
                    self.assertEqual( c[y], getattr( r, field ), 
                                      "mismatch in field %s: %s != %s" %\
                                          ( field,c[y], getattr( r, field ) ) )

            self.assertEqual( len(c), len( r ) + ncolumns )
            
            for y in range(len(c) - ncolumns):
                c[ncolumns+y] = "test_%i" % y
                r[y] = "test_%i" % y
                self.assertEqual( c[ncolumns+y], r[y] )

class TestVCF( TestParser ):

    filename = "example.vcf40.gz"

    def testOpening( self ):
        while 1:
            infile = pysam.Tabixfile( self.filename )
            infile.close()

                
            # check strings
            ref_string = "\t".join( c )
            cmp_string = str(r)
            
            self.assertEqual( ref_string, cmp_string )

if __name__ == "__main__":

    unittest.main()


