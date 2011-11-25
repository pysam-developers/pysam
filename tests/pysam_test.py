#!/usr/bin/env python
'''unit testing code for pysam.

Execute in the :file:`tests` directory as it requires the Makefile
and data files located there.
'''

import pysam
import unittest
import os, re, sys
import itertools, collections
import subprocess
import shutil
import logging

SAMTOOLS="samtools"
WORKDIR="pysam_test_work"

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

def runSamtools( cmd ):
    '''run a samtools command'''

    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode < 0:
            print >>sys.stderr, "Child was terminated by signal", -retcode
    except OSError, e:
        print >>sys.stderr, "Execution failed:", e

def getSamtoolsVersion():
    '''return samtools version'''

    pipe = subprocess.Popen(SAMTOOLS, shell=True, stderr=subprocess.PIPE).stderr
    lines = "".join(pipe.readlines())
    return re.search( "Version:\s+(\S+)", lines).groups()[0]

class BinaryTest(unittest.TestCase):
    '''test samtools command line commands and compare
    against pysam commands.

    Tests fail, if the output is not binary identical.
    '''

    first_time = True

    # a list of commands to test
    commands = \
        { 
          "view" :
              (
                ("ex1.view", "view ex1.bam > ex1.view"),
                ("pysam_ex1.view", (pysam.view, "ex1.bam" ) ),
                ),
          "view2" :
              (
                ("ex1.view", "view -bT ex1.fa -o ex1.view2 ex1.sam"),
                # note that -o ex1.view2 throws exception.
                ("pysam_ex1.view", (pysam.view, "-bT ex1.fa -oex1.view2 ex1.sam" ) ),
                ),
          "sort" :
              (
                ( "ex1.sort.bam", "sort ex1.bam ex1.sort" ),
                ( "pysam_ex1.sort.bam", (pysam.sort, "ex1.bam pysam_ex1.sort" ) ),
                ),
          "mpileup" :
              (
                ("ex1.pileup", "mpileup ex1.bam > ex1.pileup" ),
                ("pysam_ex1.mpileup", (pysam.mpileup, "ex1.bam" ) ),
                ),
          "depth" :
              (
                ("ex1.depth", "depth ex1.bam > ex1.depth" ),
                ("pysam_ex1.depth", (pysam.depth, "ex1.bam" ) ),
                ),
          "faidx" : 
              ( 
                ("ex1.fa.fai", "faidx ex1.fa"), 
                ("pysam_ex1.fa.fai", (pysam.faidx, "ex1.fa") ),
                ),
          "index":
              (
                ("ex1.bam.bai", "index ex1.bam" ),
                ("pysam_ex1.bam.bai", (pysam.index, "pysam_ex1.bam" ) ),
                ),
          "idxstats" :
              ( 
                ("ex1.idxstats", "idxstats ex1.bam > ex1.idxstats" ),
                ("pysam_ex1.idxstats", (pysam.idxstats, "pysam_ex1.bam" ) ),
                ),
          "fixmate" :
              (
                ("ex1.fixmate", "fixmate ex1.bam ex1.fixmate" ),
                ("pysam_ex1.fixmate", (pysam.fixmate, "pysam_ex1.bam pysam_ex1.fixmate") ),
                ),
          "flagstat" :
              (
                ("ex1.flagstat", "flagstat ex1.bam > ex1.flagstat" ),
                ("pysam_ex1.flagstat", (pysam.flagstat, "pysam_ex1.bam") ),
                ),
          "calmd" :
              (
                ("ex1.calmd", "calmd ex1.bam ex1.fa > ex1.calmd" ),
                ("pysam_ex1.calmd", (pysam.calmd, "pysam_ex1.bam ex1.fa") ),
                ),
          "merge" :
              (
                ("ex1.merge", "merge -f ex1.merge ex1.bam ex1.bam" ),
                # -f option does not work - following command will cause the subsequent
                # command to fail
                ("pysam_ex1.merge", (pysam.merge, "pysam_ex1.merge pysam_ex1.bam pysam_ex1.bam") ),
                ),
          "rmdup" :
              (
                ("ex1.rmdup", "rmdup ex1.bam ex1.rmdup" ),
                ("pysam_ex1.rmdup", (pysam.rmdup, "pysam_ex1.bam pysam_ex1.rmdup" )),
                ),
          "reheader" :
              (
                ( "ex1.reheader", "reheader ex1.bam ex1.bam > ex1.reheader"),
                ( "pysam_ex1.reheader", (pysam.reheader, "ex1.bam ex1.bam" ) ),
                ),
          "cat":
              (
                ( "ex1.cat", "cat ex1.bam ex1.bam > ex1.cat"),
                ( "pysam_ex1.cat", (pysam.cat, "ex1.bam ex1.bam" ) ),
                ),
          "targetcut":
              (
                ("ex1.targetcut", "targetcut ex1.bam > ex1.targetcut" ),
                ("pysam_ex1.targetcut", (pysam.targetcut, "pysam_ex1.bam") ),
                ),
          "phase":
              (
                ("ex1.phase", "phase ex1.bam > ex1.phase" ),
                ("pysam_ex1.phase", (pysam.phase, "pysam_ex1.bam") ),
                ),
          "import" :
              (
                ("ex1.bam", "import ex1.fa.fai ex1.sam.gz ex1.bam" ),
                ("pysam_ex1.bam", (pysam.samimport, "ex1.fa.fai ex1.sam.gz pysam_ex1.bam") ),
                ),
          "bam2fq":
              (
                ("ex1.bam2fq", "bam2fq ex1.bam > ex1.bam2fq" ),
                ("pysam_ex1.bam2fq", (pysam.bam2fq, "pysam_ex1.bam") ),
                ),
        }

    # some tests depend on others. The order specifies in which order
    # the samtools commands are executed.
    # The first three (faidx, import, index) need to be in that order, 
    # the rest is arbitrary.
    order = ('faidx', 'import', 'index', 
              # 'pileup1', 'pileup2', deprecated
              # 'glfview', deprecated
              'view', 'view2',
              'sort',
              'mpileup',
              'depth',
              'idxstats',
              'fixmate',
              'flagstat',
              # 'calmd',
              'merge',
              'rmdup',
              'reheader',
              'cat',
              'targetcut',
              'phase',
              'bam2fq',
              )

    def setUp( self ):
        '''setup tests. 

        For setup, all commands will be run before the first test is
        executed. Individual tests will then just compare the output
        files.
        '''
        if BinaryTest.first_time:

            # remove previous files
            if os.path.exists( WORKDIR ):
                shutil.rmtree( WORKDIR )
                
            # copy the source files to WORKDIR
            os.makedirs( WORKDIR )

            shutil.copy( "ex1.fa", os.path.join( WORKDIR, "pysam_ex1.fa" ) )
            shutil.copy( "ex1.fa", os.path.join( WORKDIR, "ex1.fa" ) )
            shutil.copy( "ex1.sam.gz", os.path.join( WORKDIR, "ex1.sam.gz" ) )
            shutil.copy( "ex1.sam", os.path.join( WORKDIR, "ex1.sam" ) )

            # cd to workdir
            savedir = os.getcwd()
            os.chdir( WORKDIR )
            
            for label in self.order:
                command = self.commands[label]
                samtools_target, samtools_command = command[0]
                try:
                    pysam_target, pysam_command = command[1]
                except ValueError, msg:
                    raise ValueError( "error while setting up %s=%s: %s" %\
                                          (label, command, msg) )
                runSamtools( " ".join( (SAMTOOLS, samtools_command )))
                pysam_method, pysam_options = pysam_command
                try:
                    output = pysam_method( *pysam_options.split(" "), raw=True)
                except pysam.SamtoolsError, msg:
                    raise pysam.SamtoolsError( "error while executing %s: options=%s: msg=%s" %\
                                                   (label, pysam_options, msg) )
                if ">" in samtools_command:
                    outfile = open( pysam_target, "wb" )
                    for line in output: outfile.write( line )
                    outfile.close()
                    
            os.chdir( savedir )
            BinaryTest.first_time = False

            

        samtools_version = getSamtoolsVersion()

        
        def _r( s ):
            # patch - remove any of the alpha/beta suffixes, i.e., 0.1.12a -> 0.1.12
            if s.count('-') > 0: s = s[0:s.find('-')]
            return re.sub( "[^0-9.]", "", s )

        if _r(samtools_version) != _r( pysam.__samtools_version__):
            raise ValueError("versions of pysam/samtools and samtools differ: %s != %s" % \
                                 (pysam.__samtools_version__,
                                  samtools_version ))

    def checkCommand( self, command ):
        if command:
            samtools_target, pysam_target = self.commands[command][0][0], self.commands[command][1][0]
            samtools_target = os.path.join( WORKDIR, samtools_target )
            pysam_target = os.path.join( WORKDIR, pysam_target )
            self.assertTrue( checkBinaryEqual( samtools_target, pysam_target ), 
                             "%s failed: files %s and %s are not the same" % (command, samtools_target, pysam_target) )
            
    def testImport( self ):
        self.checkCommand( "import" )

    def testIndex( self ):
        self.checkCommand( "index" )

    def testSort( self ):
        self.checkCommand( "sort" )

    def testMpileup( self ):
        self.checkCommand( "mpileup" )

    def testDepth( self ):
        self.checkCommand( "depth" )

    def testIdxstats( self ):
        self.checkCommand( "idxstats" )

    def testFixmate( self ):
        self.checkCommand( "fixmate" )

    def testFlagstat( self ):
        self.checkCommand( "flagstat" )
        
    def testMerge( self ):
        self.checkCommand( "merge" )

    def testRmdup( self ):
        self.checkCommand( "rmdup" )

    def testReheader( self ):
        self.checkCommand( "reheader" )

    def testCat( self ):
        self.checkCommand( "cat" )

    def testTargetcut( self ):
        self.checkCommand( "targetcut" )

    def testPhase( self ):
        self.checkCommand( "phase" )

    def testBam2fq( self ):
        self.checkCommand( "bam2fq" )

    # def testPileup1( self ):
    #     self.checkCommand( "pileup1" )
    
    # def testPileup2( self ):
    #     self.checkCommand( "pileup2" )

    # deprecated
    # def testGLFView( self ):
    #     self.checkCommand( "glfview" )

    def testView( self ):
        self.checkCommand( "view" )

    def testEmptyIndex( self ):
        self.assertRaises( IOError, pysam.index, "exdoesntexist.bam" )

    def __del__(self):
        if os.path.exists( WORKDIR ):
            shutil.rmtree( WORKDIR )

class IOTest(unittest.TestCase):
    '''check if reading samfile and writing a samfile are consistent.'''

    def checkEcho( self, input_filename, reference_filename, 
                   output_filename, 
                   input_mode, output_mode, use_template = True ):
        '''iterate through *input_filename* writing to *output_filename* and
        comparing the output to *reference_filename*. 
        
        The files are opened according to the *input_mode* and *output_mode*.

        If *use_template* is set, the header is copied from infile using the
        template mechanism, otherwise target names and lengths are passed 
        explicitely. 

        '''

        infile = pysam.Samfile( input_filename, input_mode )
        if use_template:
            outfile = pysam.Samfile( output_filename, output_mode, template = infile )
        else:
            outfile = pysam.Samfile( output_filename, output_mode, 
                                     referencenames = infile.references,
                                     referencelengths = infile.lengths,
                                     add_sq_text = False )

        iter = infile.fetch()
        for x in iter: outfile.write( x )
        infile.close()
        outfile.close()

        self.assertTrue( checkBinaryEqual( reference_filename, output_filename), 
                         "files %s and %s are not the same" % (reference_filename, output_filename) )

    def testReadWriteBam( self ):
        
        input_filename = "ex1.bam"
        output_filename = "pysam_ex1.bam"
        reference_filename = "ex1.bam"

        self.checkEcho( input_filename, reference_filename, output_filename,
                        "rb", "wb" )

    def testReadWriteBamWithTargetNames( self ):
        
        input_filename = "ex1.bam"
        output_filename = "pysam_ex1.bam"
        reference_filename = "ex1.bam"

        self.checkEcho( input_filename, reference_filename, output_filename,
                        "rb", "wb", use_template = False )

    def testReadWriteSamWithHeader( self ):
        
        input_filename = "ex2.sam"
        output_filename = "pysam_ex2.sam"
        reference_filename = "ex2.sam"

        self.checkEcho( input_filename, reference_filename, output_filename,
                        "r", "wh" )

    def testReadWriteSamWithoutHeader( self ):
        
        input_filename = "ex2.sam"
        output_filename = "pysam_ex2.sam"
        reference_filename = "ex1.sam"

        self.checkEcho( input_filename, reference_filename, output_filename,
                        "r", "w" )

    def testReadSamWithoutHeaderWriteSamWithoutHeader( self ):
        
        input_filename = "ex1.sam"
        output_filename = "pysam_ex1.sam"
        reference_filename = "ex1.sam"

        # disabled - reading from a samfile without header
        # is not implemented.
        
        # self.checkEcho( input_filename, reference_filename, output_filename,
        #                 "r", "w" )

    def testFetchFromClosedFile( self ):

        samfile = pysam.Samfile( "ex1.bam", "rb" )
        samfile.close()
        self.assertRaises( ValueError, samfile.fetch, 'chr1', 100, 120)

    def testClosedFile( self ):
        '''test that access to a closed samfile raises ValueError.'''

        samfile = pysam.Samfile( "ex1.bam", "rb" )
        samfile.close()
        self.assertRaises( ValueError, samfile.fetch, 'chr1', 100, 120)
        self.assertRaises( ValueError, samfile.pileup, 'chr1', 100, 120)
        self.assertRaises( ValueError, samfile.getrname, 0 )
        self.assertRaises( ValueError, samfile.tell )
        self.assertRaises( ValueError, samfile.seek, 0 )
        self.assertRaises( ValueError, getattr, samfile, "nreferences" )
        self.assertRaises( ValueError, getattr, samfile, "references" )
        self.assertRaises( ValueError, getattr, samfile, "lengths" )
        self.assertRaises( ValueError, getattr, samfile, "text" )
        self.assertRaises( ValueError, getattr, samfile, "header" )

        # write on closed file 
        self.assertEqual( 0, samfile.write(None) )

    def testAutoDetection( self ):
        '''test if autodetection works.'''

        samfile = pysam.Samfile( "ex3.sam" )
        self.assertRaises( ValueError, samfile.fetch, 'chr1' )
        samfile.close()

        samfile = pysam.Samfile( "ex3.bam" )
        samfile.fetch('chr1')
        samfile.close()

    def testReadingFromSamFileWithoutHeader( self ):
        '''read from samfile without header.
        '''
        samfile = pysam.Samfile( "ex7.sam" )
        self.assertRaises( NotImplementedError, samfile.__iter__ )

    def testReadingFromFileWithoutIndex( self ):
        '''read from bam file without index.'''

        assert not os.path.exists( "ex2.bam.bai" )
        samfile = pysam.Samfile( "ex2.bam", "rb" )
        self.assertRaises( ValueError, samfile.fetch )
        self.assertEqual( len(list( samfile.fetch(until_eof = True) )), 3270 )

    def testReadingUniversalFileMode( self ):
        '''read from samfile without header.
        '''

        input_filename = "ex2.sam"
        output_filename = "pysam_ex2.sam"
        reference_filename = "ex1.sam"

        self.checkEcho( input_filename, reference_filename, output_filename,
                        "rU", "w" )

class TestFloatTagBug( unittest.TestCase ):
    '''see issue 71'''

    def testFloatTagBug( self ): 
        '''a float tag before another exposed a parsing bug in bam_aux_get - expected to fail

        This test is expected to fail until samtools is fixed.
        '''
        samfile = pysam.Samfile("tag_bug.bam")
        read = samfile.fetch(until_eof=True).next()
        self.assertTrue( ('XC',1) in read.tags )
        self.assertEqual(read.opt('XC'), 1)

class TestTagParsing( unittest.TestCase ):
    '''tests checking the accuracy of tag setting and retrieval.'''

    def makeRead( self ):
        a = pysam.AlignedRead()
        a.qname = "read_12345"
        a.tid = 0
        a.seq="ACGT" * 3
        a.flag = 0
        a.rname = 0
        a.pos = 1
        a.mapq = 20
        a.cigar = ( (0,10), (2,1), (0,25) )
        a.mrnm = 0
        a.mpos=200
        a.isize = 0
	a.qual ="1234" * 3
        # todo: create tags
        return a

    def testNegativeIntegers( self ):
        x = -2
        aligned_read = self.makeRead()
        aligned_read.tags = [("XD", int(x) ) ]
        print aligned_read.tags

    def testNegativeIntegers2( self ):
        x = -2
        r = self.makeRead()
        r.tags = [("XD", int(x) ) ]
        outfile = pysam.Samfile( "test.bam",
                                 "wb",
                                 referencenames = ("chr1",),
                                 referencelengths = (1000,) )
        outfile.write (r )
        outfile.close()

class TestIteratorRow(unittest.TestCase):

    def setUp(self):
        self.samfile=pysam.Samfile( "ex1.bam","rb" )

    def checkRange( self, rnge ):
        '''compare results from iterator with those from samtools.'''
        ps = list(self.samfile.fetch(region=rnge))
        sa = list(pysam.view( "ex1.bam", rnge, raw = True) )
        self.assertEqual( len(ps), len(sa), "unequal number of results for range %s: %i != %i" % (rnge, len(ps), len(sa) ))
        # check if the same reads are returned and in the same order
        for line, pair in enumerate( zip( ps, sa ) ):
            a,b = pair
            d = b.split("\t")
            self.assertEqual( a.qname, d[0], "line %i: read id mismatch: %s != %s" % (line, a.rname, d[0]) )
            self.assertEqual( a.pos, int(d[3])-1, "line %i: read position mismatch: %s != %s, \n%s\n%s\n" % \
                                  (line, a.pos, int(d[3])-1,
                                   str(a), str(d) ) )
            self.assertEqual( a.qual, d[10], "line %i: quality mismatch: %s != %s, \n%s\n%s\n" % \
                                  (line, a.qual, d[10],
                                   str(a), str(d) ) )

    def testIteratePerContig(self):
        '''check random access per contig'''
        for contig in self.samfile.references:
            self.checkRange( contig )

    def testIterateRanges(self):
        '''check random access per range'''
        for contig, length in zip(self.samfile.references, self.samfile.lengths):
            for start in range( 1, length, 90):
                self.checkRange( "%s:%i-%i" % (contig, start, start + 90) ) # this includes empty ranges

    def tearDown(self):
        self.samfile.close()

class TestIteratorRowAll(unittest.TestCase):

    def setUp(self):
        self.samfile=pysam.Samfile( "ex1.bam","rb" )

    def testIterate(self):
        '''compare results from iterator with those from samtools.'''
        ps = list(self.samfile.fetch())
        sa = list(pysam.view( "ex1.bam", raw = True) )
        self.assertEqual( len(ps), len(sa), "unequal number of results: %i != %i" % (len(ps), len(sa) ))
        # check if the same reads are returned
        for line, pair in enumerate( zip( ps, sa ) ):
            data = pair[1].split("\t")
            self.assertEqual( pair[0].qname, data[0], "read id mismatch in line %i: %s != %s" % (line, pair[0].rname, data[0]) )

    def tearDown(self):
        self.samfile.close()

class TestIteratorColumn(unittest.TestCase):
    '''test iterator column against contents of ex3.bam.'''
    
    # note that samfile contains 1-based coordinates
    # 1D means deletion with respect to reference sequence
    # 
    mCoverages = { 'chr1' : [ 0 ] * 20 + [1] * 36 + [0] * (100 - 20 -35 ),
                   'chr2' : [ 0 ] * 20 + [1] * 35 + [0] * (100 - 20 -35 ),
                   }

    def setUp(self):
        self.samfile=pysam.Samfile( "ex4.bam","rb" )

    def checkRange( self, rnge ):
        '''compare results from iterator with those from samtools.'''
        # check if the same reads are returned and in the same order
        for column in self.samfile.pileup(region=rnge):
            thiscov = len(column.pileups)
            refcov = self.mCoverages[self.samfile.getrname(column.tid)][column.pos]
            self.assertEqual( thiscov, refcov, "wrong coverage at pos %s:%i %i should be %i" % (self.samfile.getrname(column.tid), column.pos, thiscov, refcov))

    def testIterateAll(self):
        '''check random access per contig'''
        self.checkRange( None )

    def testIteratePerContig(self):
        '''check random access per contig'''
        for contig in self.samfile.references:
            self.checkRange( contig )

    def testIterateRanges(self):
        '''check random access per range'''
        for contig, length in zip(self.samfile.references, self.samfile.lengths):
            for start in range( 1, length, 90):
                self.checkRange( "%s:%i-%i" % (contig, start, start + 90) ) # this includes empty ranges

    def testInverse( self ):
        '''test the inverse, is point-wise pileup accurate.'''
        for contig, refseq in self.mCoverages.items():
            refcolumns = sum(refseq)
            for pos, refcov in enumerate( refseq ):
                columns = list(self.samfile.pileup( contig, pos, pos+1) )
                if refcov == 0:
                    # if no read, no coverage
                    self.assertEqual( len(columns), refcov, "wrong number of pileup columns returned for position %s:%i, %i should be %i" %(contig,pos,len(columns), refcov) )
                elif refcov == 1:
                    # one read, all columns of the read are returned
                    self.assertEqual( len(columns), refcolumns, "pileup incomplete at position %i: got %i, expected %i " %\
                                          (pos, len(columns), refcolumns))

                    
    
    def tearDown(self):
        self.samfile.close()
    
class TestAlignedReadFromBam(unittest.TestCase):

    def setUp(self):
        self.samfile=pysam.Samfile( "ex3.bam","rb" )
        self.reads=list(self.samfile.fetch())

    def testARqname(self):
        self.assertEqual( self.reads[0].qname, "read_28833_29006_6945", "read name mismatch in read 1: %s != %s" % (self.reads[0].qname, "read_28833_29006_6945") )
        self.assertEqual( self.reads[1].qname, "read_28701_28881_323b", "read name mismatch in read 2: %s != %s" % (self.reads[1].qname, "read_28701_28881_323b") )

    def testARflag(self):
        self.assertEqual( self.reads[0].flag, 99, "flag mismatch in read 1: %s != %s" % (self.reads[0].flag, 99) )
        self.assertEqual( self.reads[1].flag, 147, "flag mismatch in read 2: %s != %s" % (self.reads[1].flag, 147) )

    def testARrname(self):
        self.assertEqual( self.reads[0].rname, 0, "chromosome/target id mismatch in read 1: %s != %s" % (self.reads[0].rname, 0) )
        self.assertEqual( self.reads[1].rname, 1, "chromosome/target id mismatch in read 2: %s != %s" % (self.reads[1].rname, 1) )

    def testARpos(self):
        self.assertEqual( self.reads[0].pos, 33-1, "mapping position mismatch in read 1: %s != %s" % (self.reads[0].pos, 33-1) )
        self.assertEqual( self.reads[1].pos, 88-1, "mapping position mismatch in read 2: %s != %s" % (self.reads[1].pos, 88-1) )

    def testARmapq(self):
        self.assertEqual( self.reads[0].mapq, 20, "mapping quality mismatch in read 1: %s != %s" % (self.reads[0].mapq, 20) )
        self.assertEqual( self.reads[1].mapq, 30, "mapping quality mismatch in read 2: %s != %s" % (self.reads[1].mapq, 30) )

    def testARcigar(self):
        self.assertEqual( self.reads[0].cigar, [(0, 10), (2, 1), (0, 25)], "read name length mismatch in read 1: %s != %s" % (self.reads[0].cigar, [(0, 10), (2, 1), (0, 25)]) )
        self.assertEqual( self.reads[1].cigar, [(0, 35)], "read name length mismatch in read 2: %s != %s" % (self.reads[1].cigar, [(0, 35)]) )

    def testARmrnm(self):
        self.assertEqual( self.reads[0].mrnm, 0, "mate reference sequence name mismatch in read 1: %s != %s" % (self.reads[0].mrnm, 0) )
        self.assertEqual( self.reads[1].mrnm, 1, "mate reference sequence name mismatch in read 2: %s != %s" % (self.reads[1].mrnm, 1) )
        self.assertEqual( self.reads[0].rnext, 0, "mate reference sequence name mismatch in read 1: %s != %s" % (self.reads[0].rnext, 0) )
        self.assertEqual( self.reads[1].rnext, 1, "mate reference sequence name mismatch in read 2: %s != %s" % (self.reads[1].rnext, 1) )

    def testARmpos(self):
        self.assertEqual( self.reads[0].mpos, 200-1, "mate mapping position mismatch in read 1: %s != %s" % (self.reads[0].mpos, 200-1) )
        self.assertEqual( self.reads[1].mpos, 500-1, "mate mapping position mismatch in read 2: %s != %s" % (self.reads[1].mpos, 500-1) )
        self.assertEqual( self.reads[0].pnext, 200-1, "mate mapping position mismatch in read 1: %s != %s" % (self.reads[0].pnext, 200-1) )
        self.assertEqual( self.reads[1].pnext, 500-1, "mate mapping position mismatch in read 2: %s != %s" % (self.reads[1].pnext, 500-1) )

    def testARisize(self):
        self.assertEqual( self.reads[0].isize, 167, "insert size mismatch in read 1: %s != %s" % (self.reads[0].isize, 167) )
        self.assertEqual( self.reads[1].isize, 412, "insert size mismatch in read 2: %s != %s" % (self.reads[1].isize, 412) )
        self.assertEqual( self.reads[0].tlen, 167, "insert size mismatch in read 1: %s != %s" % (self.reads[0].tlen, 167) )
        self.assertEqual( self.reads[1].tlen, 412, "insert size mismatch in read 2: %s != %s" % (self.reads[1].tlen, 412) )

    def testARseq(self):
        self.assertEqual( self.reads[0].seq, "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG", "sequence mismatch in read 1: %s != %s" % (self.reads[0].seq, "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG") )
        self.assertEqual( self.reads[1].seq, "ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA", "sequence size mismatch in read 2: %s != %s" % (self.reads[1].seq, "ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA") )
        self.assertEqual( self.reads[3].seq, "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG", "sequence mismatch in read 4: %s != %s" % (self.reads[3].seq, "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG") )

    def testARqual(self):
        self.assertEqual( self.reads[0].qual, "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<", "quality string mismatch in read 1: %s != %s" % (self.reads[0].qual, "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<") )
        self.assertEqual( self.reads[1].qual, "<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<", "quality string mismatch in read 2: %s != %s" % (self.reads[1].qual, "<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<") )
        self.assertEqual( self.reads[3].qual, "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<", "quality string mismatch in read 3: %s != %s" % (self.reads[3].qual, "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<") )

    def testARquery(self):
        self.assertEqual( self.reads[0].query, "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG", "query mismatch in read 1: %s != %s" % (self.reads[0].query, "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG") )
        self.assertEqual( self.reads[1].query, "ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA", "query size mismatch in read 2: %s != %s" % (self.reads[1].query, "ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA") )
        self.assertEqual( self.reads[3].query, "TAGCTAGCTACCTATATCTTGGTCTT", "query mismatch in read 4: %s != %s" % (self.reads[3].query, "TAGCTAGCTACCTATATCTTGGTCTT") )

    def testARqqual(self):
        self.assertEqual( self.reads[0].qqual, "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<", "qquality string mismatch in read 1: %s != %s" % (self.reads[0].qqual, "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<") )
        self.assertEqual( self.reads[1].qqual, "<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<", "qquality string mismatch in read 2: %s != %s" % (self.reads[1].qqual, "<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<") )
        self.assertEqual( self.reads[3].qqual, "<<<<<<<<<<<<<<<<<:<9/,&,22", "qquality string mismatch in read 3: %s != %s" % (self.reads[3].qqual, "<<<<<<<<<<<<<<<<<:<9/,&,22") )

    def testPresentOptionalFields(self):
        self.assertEqual( self.reads[0].opt('NM'), 1, "optional field mismatch in read 1, NM: %s != %s" % (self.reads[0].opt('NM'), 1) )
        self.assertEqual( self.reads[0].opt('RG'), 'L1', "optional field mismatch in read 1, RG: %s != %s" % (self.reads[0].opt('RG'), 'L1') )
        self.assertEqual( self.reads[1].opt('RG'), 'L2', "optional field mismatch in read 2, RG: %s != %s" % (self.reads[1].opt('RG'), 'L2') )
        self.assertEqual( self.reads[1].opt('MF'), 18, "optional field mismatch in read 2, MF: %s != %s" % (self.reads[1].opt('MF'), 18) )

    def testPairedBools(self):
        self.assertEqual( self.reads[0].is_paired, True, "is paired mismatch in read 1: %s != %s" % (self.reads[0].is_paired, True) )
        self.assertEqual( self.reads[1].is_paired, True, "is paired mismatch in read 2: %s != %s" % (self.reads[1].is_paired, True) )
        self.assertEqual( self.reads[0].is_proper_pair, True, "is proper pair mismatch in read 1: %s != %s" % (self.reads[0].is_proper_pair, True) )
        self.assertEqual( self.reads[1].is_proper_pair, True, "is proper pair mismatch in read 2: %s != %s" % (self.reads[1].is_proper_pair, True) )

    def testTags( self ):
        self.assertEqual( self.reads[0].tags, 
                          [('NM', 1), ('RG', 'L1'), 
                           ('PG', 'P1'), ('XT', 'U')] )
        self.assertEqual( self.reads[1].tags, 
                          [('MF', 18), ('RG', 'L2'), 
                           ('PG', 'P2'),('XT', 'R') ] )

    def testOpt( self ):
        self.assertEqual( self.reads[0].opt("XT"), "U" )
        self.assertEqual( self.reads[1].opt("XT"), "R" )

    def testMissingOpt( self ):
        self.assertRaises( KeyError, self.reads[0].opt, "XP" )

    def testEmptyOpt( self ):
        self.assertRaises( KeyError, self.reads[2].opt, "XT" )

    def tearDown(self):
        self.samfile.close()

class TestAlignedReadFromSam(TestAlignedReadFromBam):

    def setUp(self):
        self.samfile=pysam.Samfile( "ex3.sam","r" )
        self.reads=list(self.samfile.fetch())

# needs to be implemented 
# class TestAlignedReadFromSamWithoutHeader(TestAlignedReadFromBam):
#
#     def setUp(self):
#         self.samfile=pysam.Samfile( "ex7.sam","r" )
#         self.reads=list(self.samfile.fetch())

class TestHeaderSam(unittest.TestCase):

    header = {'SQ': [{'LN': 1575, 'SN': 'chr1'}, 
                     {'LN': 1584, 'SN': 'chr2'}], 
              'RG': [{'LB': 'SC_1', 'ID': 'L1', 'SM': 'NA12891', 'PU': 'SC_1_10', "CN":"name:with:colon"}, 
                     {'LB': 'SC_2', 'ID': 'L2', 'SM': 'NA12891', 'PU': 'SC_2_12', "CN":"name:with:colon"}],
              'PG': [{'ID': 'P1', 'VN': '1.0'}, {'ID': 'P2', 'VN': '1.1'}], 
              'HD': {'VN': '1.0'},
              'CO' : [ 'this is a comment', 'this is another comment'],
              }

    def compareHeaders( self, a, b ):
        '''compare two headers a and b.'''
        for ak,av in a.iteritems():
            self.assertTrue( ak in b, "key '%s' not in '%s' " % (ak,b) )
            self.assertEqual( av, b[ak] )

    def setUp(self):
        self.samfile=pysam.Samfile( "ex3.sam","r" )

    def testHeaders(self):
        self.compareHeaders( self.header, self.samfile.header )
        self.compareHeaders( self.samfile.header, self.header )

    def testNameMapping( self ):
        for x, y in enumerate( ("chr1", "chr2")):
            tid = self.samfile.gettid( y )
            ref = self.samfile.getrname( x )
            self.assertEqual( tid, x )
            self.assertEqual( ref, y )

        self.assertEqual( self.samfile.gettid("chr?"), -1 )
        self.assertRaises( ValueError, self.samfile.getrname, 2 )

    def tearDown(self):
        self.samfile.close()

class TestHeaderBam(TestHeaderSam):

    def setUp(self):
        self.samfile=pysam.Samfile( "ex3.bam","rb" )

class TestUnmappedReads(unittest.TestCase):

    def testSAM(self):
        samfile=pysam.Samfile( "ex5.sam","r" )
        self.assertEqual( len(list(samfile.fetch( until_eof = True))), 2 ) 
        samfile.close()

    def testBAM(self):
        samfile=pysam.Samfile( "ex5.bam","rb" )
        self.assertEqual( len(list(samfile.fetch( until_eof = True))), 2 ) 
        samfile.close()

class TestPileupObjects(unittest.TestCase):

    def setUp(self):
        self.samfile=pysam.Samfile( "ex1.bam","rb" )

    def testPileupColumn(self):
        for pcolumn1 in self.samfile.pileup( region="chr1:105" ):
            if pcolumn1.pos == 104:
                self.assertEqual( pcolumn1.tid, 0, "chromosome/target id mismatch in position 1: %s != %s" % (pcolumn1.tid, 0) )
                self.assertEqual( pcolumn1.pos, 105-1, "position mismatch in position 1: %s != %s" % (pcolumn1.pos, 105-1) )
                self.assertEqual( pcolumn1.n, 2, "# reads mismatch in position 1: %s != %s" % (pcolumn1.n, 2) )
        for pcolumn2 in self.samfile.pileup( region="chr2:1480" ):
            if pcolumn2.pos == 1479:
                self.assertEqual( pcolumn2.tid, 1, "chromosome/target id mismatch in position 1: %s != %s" % (pcolumn2.tid, 1) )
                self.assertEqual( pcolumn2.pos, 1480-1, "position mismatch in position 1: %s != %s" % (pcolumn2.pos, 1480-1) )
                self.assertEqual( pcolumn2.n, 12, "# reads mismatch in position 1: %s != %s" % (pcolumn2.n, 12) )

    def testPileupRead(self):
        for pcolumn1 in self.samfile.pileup( region="chr1:105" ):
            if pcolumn1.pos == 104:
                self.assertEqual( len(pcolumn1.pileups), 2, "# reads aligned to column mismatch in position 1: %s != %s" % (len(pcolumn1.pileups), 2) )
#                self.assertEqual( pcolumn1.pileups[0]  # need to test additional properties here

    def tearDown(self):
        self.samfile.close()

class TestContextManager(unittest.TestCase):

    def testManager( self ):
        with pysam.Samfile('ex1.bam', 'rb') as samfile:
            samfile.fetch()
        self.assertEqual( samfile._isOpen(), False )

class TestExceptions(unittest.TestCase):

    def setUp(self):
        self.samfile=pysam.Samfile( "ex1.bam","rb" )

    def testMissingFile(self):

        self.assertRaises( IOError, pysam.Samfile, "exdoesntexist.bam", "rb" )
        self.assertRaises( IOError, pysam.Samfile, "exdoesntexist.sam", "r" )
        self.assertRaises( IOError, pysam.Samfile, "exdoesntexist.bam", "r" )
        self.assertRaises( IOError, pysam.Samfile, "exdoesntexist.sam", "rb" )

    def testBadContig(self):
        self.assertRaises( ValueError, self.samfile.fetch, "chr88" )

    def testMeaninglessCrap(self):
        self.assertRaises( ValueError, self.samfile.fetch, "skljf" )

    def testBackwardsOrderNewFormat(self):
        self.assertRaises( ValueError, self.samfile.fetch, 'chr1', 100, 10 )

    def testBackwardsOrderOldFormat(self):
        self.assertRaises( ValueError, self.samfile.fetch, region="chr1:100-10")
        
    def testOutOfRangeNegativeNewFormat(self):
        self.assertRaises( ValueError, self.samfile.fetch, "chr1", 5, -10 )
        self.assertRaises( ValueError, self.samfile.fetch, "chr1", 5, 0 )
        self.assertRaises( ValueError, self.samfile.fetch, "chr1", -5, -10 )

        self.assertRaises( ValueError, self.samfile.count, "chr1", 5, -10 )
        self.assertRaises( ValueError, self.samfile.count, "chr1", 5, 0 )        
        self.assertRaises( ValueError, self.samfile.count, "chr1", -5, -10 )

    def testOutOfRangeNegativeOldFormat(self):
        self.assertRaises( ValueError, self.samfile.fetch, region="chr1:-5-10" )
        self.assertRaises( ValueError, self.samfile.fetch, region="chr1:-5-0" )
        self.assertRaises( ValueError, self.samfile.fetch, region="chr1:-5--10" )

        self.assertRaises( ValueError, self.samfile.count, region="chr1:-5-10" )
        self.assertRaises( ValueError, self.samfile.count, region="chr1:-5-0" )
        self.assertRaises( ValueError, self.samfile.count, region="chr1:-5--10" )

    def testOutOfRangNewFormat(self):
        self.assertRaises( ValueError, self.samfile.fetch, "chr1", 9999999999, 99999999999 )
        self.assertRaises( ValueError, self.samfile.count, "chr1", 9999999999, 99999999999 )

    def testOutOfRangeLargeNewFormat(self):
        self.assertRaises( ValueError, self.samfile.fetch, "chr1", 9999999999999999999999999999999, 9999999999999999999999999999999999999999 )
        self.assertRaises( ValueError, self.samfile.count, "chr1", 9999999999999999999999999999999, 9999999999999999999999999999999999999999 )

    def testOutOfRangeLargeOldFormat(self):
        self.assertRaises( ValueError, self.samfile.fetch, "chr1:99999999999999999-999999999999999999" )
        self.assertRaises( ValueError, self.samfile.count, "chr1:99999999999999999-999999999999999999" )

    def testZeroToZero(self):        
        '''see issue 44'''
        self.assertEqual( len(list(self.samfile.fetch('chr1', 0, 0))), 0)

    def tearDown(self):
        self.samfile.close()

class TestWrongFormat(unittest.TestCase):
    '''test cases for opening files not in bam/sam format.'''

    def testOpenSamAsBam( self ):
        self.assertRaises( ValueError, pysam.Samfile, 'ex1.sam', 'rb' )

    def testOpenBamAsSam( self ):
        # test fails, needs to be implemented.
        # sam.fetch() fails on reading, not on opening
        # self.assertRaises( ValueError, pysam.Samfile, 'ex1.bam', 'r' )
        pass

    def testOpenFastaAsSam( self ):
        # test fails, needs to be implemented.
        # sam.fetch() fails on reading, not on opening
        # self.assertRaises( ValueError, pysam.Samfile, 'ex1.fa', 'r' )
        pass

    def testOpenFastaAsBam( self ):
        self.assertRaises( ValueError, pysam.Samfile, 'ex1.fa', 'rb' )

class TestFastaFile(unittest.TestCase):

    mSequences = { 'chr1' :
                       "CACTAGTGGCTCATTGTAAATGTGTGGTTTAACTCGTCCATGGCCCAGCATTAGGGAGCTGTGGACCCTGCAGCCTGGCTGTGGGGGCCGCAGTGGCTGAGGGGTGCAGAGCCGAGTCACGGGGTTGCCAGCACAGGGGCTTAACCTCTGGTGACTGCCAGAGCTGCTGGCAAGCTAGAGTCCCATTTGGAGCCCCTCTAAGCCGTTCTATTTGTAATGAAAACTATATTTATGCTATTCAGTTCTAAATATAGAAATTGAAACAGCTGTGTTTAGTGCCTTTGTTCAACCCCCTTGCAACAACCTTGAGAACCCCAGGGAATTTGTCAATGTCAGGGAAGGAGCATTTTGTCAGTTACCAAATGTGTTTATTACCAGAGGGATGGAGGGAAGAGGGACGCTGAAGAACTTTGATGCCCTCTTCTTCCAAAGATGAAACGCGTAACTGCGCTCTCATTCACTCCAGCTCCCTGTCACCCAATGGACCTGTGATATCTGGATTCTGGGAAATTCTTCATCCTGGACCCTGAGAGATTCTGCAGCCCAGCTCCAGATTGCTTGTGGTCTGACAGGCTGCAACTGTGAGCCATCACAATGAACAACAGGAAGAAAAGGTCTTTCAAAAGGTGATGTGTGTTCTCATCAACCTCATACACACACATGGTTTAGGGGTATAATACCTCTACATGGCTGATTATGAAAACAATGTTCCCCAGATACCATCCCTGTCTTACTTCCAGCTCCCCAGAGGGAAAGCTTTCAACGCTTCTAGCCATTTCTTTTGGCATTTGCCTTCAGACCCTACACGAATGCGTCTCTACCACAGGGGGCTGCGCGGTTTCCCATCATGAAGCACTGAACTTCCACGTCTCATCTAGGGGAACAGGGAGGTGCACTAATGCGCTCCACGCCCAAGCCCTTCTCACAGTTTCTGCCCCCAGCATGGTTGTACTGGGCAATACATGAGATTATTAGGAAATGCTTTACTGTCATAACTATGAAGAGACTATTGCCAGATGAACCACACATTAATACTATGTTTCTTATCTGCACATTACTACCCTGCAATTAATATAATTGTGTCCATGTACACACGCTGTCCTATGTACTTATCATGACTCTATCCCAAATTCCCAATTACGTCCTATCTTCTTCTTAGGGAAGAACAGCTTAGGTATCAATTTGGTGTTCTGTGTAAAGTCTCAGGGAGCCGTCCGTGTCCTCCCATCTGGCCTCGTCCACACTGGTTCTCTTGAAAGCTTGGGCTGTAATGATGCCCCTTGGCCATCACCCAGTCCCTGCCCCATCTCTTGTAATCTCTCTCCTTTTTGCTGCATCCCTGTCTTCCTCTGTCTTGATTTACTTGTTGTTGGTTTTCTGTTTCTTTGTTTGATTTGGTGGAAGACATAATCCCACGCTTCCTATGGAAAGGTTGTTGGGAGATTTTTAATGATTCCTCAATGTTAAAATGTCTATTTTTGTCTTGACACCCAACTAATATTTGTCTGAGCAAAACAGTCTAGATGAGAGAGAACTTCCCTGGAGGTCTGATGGCGTTTCTCCCTCGTCTTCTTA",
                   'chr2' :
                       "TTCAAATGAACTTCTGTAATTGAAAAATTCATTTAAGAAATTACAAAATATAGTTGAAAGCTCTAACAATAGACTAAACCAAGCAGAAGAAAGAGGTTCAGAACTTGAAGACAAGTCTCTTATGAATTAACCCAGTCAGACAAAAATAAAGAAAAAAATTTTAAAAATGAACAGAGCTTTCAAGAAGTATGAGATTATGTAAAGTAACTGAACCTATGAGTCACAGGTATTCCTGAGGAAAAAGAAAAAGTGAGAAGTTTGGAAAAACTATTTGAGGAAGTAATTGGGGAAAACCTCTTTAGTCTTGCTAGAGATTTAGACATCTAAATGAAAGAGGCTCAAAGAATGCCAGGAAGATACATTGCAAGACAGACTTCATCAAGATATGTAGTCATCAGACTATCTAAAGTCAACATGAAGGAAAAAAATTCTAAAATCAGCAAGAGAAAAGCATACAGTCATCTATAAAGGAAATCCCATCAGAATAACAATGGGCTTCTCAGCAGAAACCTTACAAGCCAGAAGAGATTGGATCTAATTTTTGGACTTCTTAAAGAAAAAAAAACCTGTCAAACACGAATGTTATGCCCTGCTAAACTAAGCATCATAAATGAAGGGGAAATAAAGTCAAGTCTTTCCTGACAAGCAAATGCTAAGATAATTCATCATCACTAAACCAGTCCTATAAGAAATGCTCAAAAGAATTGTAAAAGTCAAAATTAAAGTTCAATACTCACCATCATAAATACACACAAAAGTACAAAACTCACAGGTTTTATAAAACAATTGAGACTACAGAGCAACTAGGTAAAAAATTAACATTACAACAGGAACAAAACCTCATATATCAATATTAACTTTGAATAAAAAGGGATTAAATTCCCCCACTTAAGAGATATAGATTGGCAGAACAGATTTAAAAACATGAACTAACTATATGCTGTTTACAAGAAACTCATTAATAAAGACATGAGTTCAGGTAAAGGGGTGGAAAAAGATGTTCTACGCAAACAGAAACCAAATGAGAGAAGGAGTAGCTATACTTATATCAGATAAAGCACACTTTAAATCAACAACAGTAAAATAAAACAAAGGAGGTCATCATACAATGATAAAAAGATCAATTCAGCAAGAAGATATAACCATCCTACTAAATACATATGCACCTAACACAAGACTACCCAGATTCATAAAACAAATACTACTAGACCTAAGAGGGATGAGAAATTACCTAATTGGTACAATGTACAATATTCTGATGATGGTTACACTAAAAGCCCATACTTTACTGCTACTCAATATATCCATGTAACAAATCTGCGCTTGTACTTCTAAATCTATAAAAAAATTAAAATTTAACAAAAGTAAATAAAACACATAGCTAAAACTAAAAAAGCAAAAACAAAAACTATGCTAAGTATTGGTAAAGATGTGGGGAAAAAAGTAAACTCTCAAATATTGCTAGTGGGAGTATAAATTGTTTTCCACTTTGGAAAACAATTTGGTAATTTCGTTTTTTTTTTTTTCTTTTCTCTTTTTTTTTTTTTTTTTTTTGCATGCCAGAAAAAAATATTTACAGTAACT",
                   }

    def setUp(self):
        self.file=pysam.Fastafile( "ex1.fa" )

    def testFetch(self):
        for id, seq in self.mSequences.items():
            self.assertEqual( seq, self.file.fetch( id ) )
            for x in range( 0, len(seq), 10):
                self.assertEqual( seq[x:x+10], self.file.fetch( id, x, x+10) )
                # test x:end
                self.assertEqual( seq[x:], self.file.fetch( id, x) )
                # test 0:x
                self.assertEqual( seq[:x], self.file.fetch( id, None, x) )

        
        # unknown sequence returns ""
        self.assertEqual( "", self.file.fetch("chr12") )

    def testOutOfRangeAccess( self ):
        '''test out of range access.'''
        # out of range access returns an empty string
        for contig, s in self.mSequences.iteritems():
            self.assertEqual( self.file.fetch( contig, len(s), len(s)+1), "" )

        self.assertEqual( self.file.fetch( "chr3", 0 , 100), "" ) 

    def testFetchErrors( self ):
        self.assertRaises( ValueError, self.file.fetch )
        self.assertRaises( ValueError, self.file.fetch, "chr1", -1, 10 )
        self.assertRaises( ValueError, self.file.fetch, "chr1", 20, 10 )

    def testLength( self ):
        self.assertEqual( len(self.file), 2 )
        
    def tearDown(self):
        self.file.close()

class TestAlignedRead(unittest.TestCase):
    '''tests to check if aligned read can be constructed
    and manipulated.
    '''

    def checkFieldEqual( self, read1, read2, exclude = []):
        '''check if two reads are equal by comparing each field.'''

        for x in ("qname", "seq", "flag",
                  "rname", "pos", "mapq", "cigar",
                  "mrnm", "mpos", "isize", "qual",
                  "is_paired", "is_proper_pair",
                  "is_unmapped", "mate_is_unmapped",
                  "is_reverse", "mate_is_reverse",
                  "is_read1", "is_read2",
                  "is_secondary", "is_qcfail",
                  "is_duplicate", "bin"):
            if x in exclude: continue
            self.assertEqual( getattr(read1, x), getattr(read2,x), "attribute mismatch for %s: %s != %s" % 
                              (x, getattr(read1, x), getattr(read2,x)))
    
    def testEmpty( self ):
        a = pysam.AlignedRead()
        self.assertEqual( a.qname, None )
        self.assertEqual( a.seq, None )
        self.assertEqual( a.qual, None )
        self.assertEqual( a.flag, 0 )
        self.assertEqual( a.rname, 0 )
        self.assertEqual( a.mapq, 0 )
        self.assertEqual( a.cigar, None )
        self.assertEqual( a.tags, [] )
        self.assertEqual( a.mrnm, 0 )
        self.assertEqual( a.mpos, 0 )
        self.assertEqual( a.isize, 0 )

    def buildRead( self ):
        '''build an example read.'''
        
        a = pysam.AlignedRead()
        a.qname = "read_12345"
        a.seq="ACGT" * 3
        a.flag = 0
        a.rname = 0
        a.pos = 33
        a.mapq = 20
        a.cigar = ( (0,10), (2,1), (0,25) )
        a.mrnm = 0
        a.mpos=200
        a.isize=167
	a.qual="1234" * 3
        # todo: create tags
        return a

    def testUpdate( self ):
        '''check if updating fields affects other variable length data
        '''
        a = self.buildRead()
        b = self.buildRead()

        # check qname
        b.qname = "read_123"
        self.checkFieldEqual( a, b, "qname" )
        b.qname = "read_12345678"
        self.checkFieldEqual( a, b, "qname" )
        b.qname = "read_12345"
        self.checkFieldEqual( a, b)

        # check cigar
        b.cigar = ( (0,10), )
        self.checkFieldEqual( a, b, "cigar" )
        b.cigar = ( (0,10), (2,1), (0,25), (2,1), (0,25) )
        self.checkFieldEqual( a, b, "cigar" )
        b.cigar = ( (0,10), (2,1), (0,25) )
        self.checkFieldEqual( a, b)

        # check seq 
        b.seq = "ACGT"
        self.checkFieldEqual( a, b, ("seq", "qual") )
        b.seq = "ACGT" * 10
        self.checkFieldEqual( a, b, ("seq", "qual") )
        b.seq = "ACGT" * 3
        self.checkFieldEqual( a, b, ("qual",))

        # reset qual
        b = self.buildRead()

        # check flags:
        for x in (
            "is_paired", "is_proper_pair",
            "is_unmapped", "mate_is_unmapped",
            "is_reverse", "mate_is_reverse",
            "is_read1", "is_read2",
            "is_secondary", "is_qcfail",
            "is_duplicate"):
            setattr( b, x, True )
            self.assertEqual( getattr(b, x), True )
            self.checkFieldEqual( a, b, ("flag", x,) )
            setattr( b, x, False )
            self.assertEqual( getattr(b, x), False )
            self.checkFieldEqual( a, b )

    def testLargeRead( self ):
        '''build an example read.'''
        
        a = pysam.AlignedRead()
        a.qname = "read_12345"
        a.seq="ACGT" * 200
        a.flag = 0
        a.rname = 0
        a.pos = 33
        a.mapq = 20
        a.cigar = ( (0,10), (2,1), (0,25) )
        a.mrnm = 0
        a.mpos=200
        a.isize=167
	a.qual="1234" * 200

        return a

    def testTagParsing( self ):
        '''test for tag parsing

        see http://groups.google.com/group/pysam-user-group/browse_thread/thread/67ca204059ea465a
        '''
        samfile=pysam.Samfile( "ex8.bam","rb" )

        for entry in samfile:
            before = entry.tags
            entry.tags = entry.tags
            after = entry.tags
            self.assertEqual( after, before )

class TestDeNovoConstruction(unittest.TestCase):
    '''check BAM/SAM file construction using ex3.sam
    
    (note these are +1 coordinates):
    
    read_28833_29006_6945	99	chr1	33	20	10M1D25M	=	200	167	AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG	<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<	NM:i:1	RG:Z:L1
    read_28701_28881_323b	147	chr2	88	30	35M	=	500	412	ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA	<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<	MF:i:18	RG:Z:L2
    '''

    header = { 'HD': {'VN': '1.0'},
               'SQ': [{'LN': 1575, 'SN': 'chr1'}, 
                      {'LN': 1584, 'SN': 'chr2'}], }

    bamfile = "ex6.bam"
    samfile = "ex6.sam"

    def checkFieldEqual( self, read1, read2, exclude = []):
        '''check if two reads are equal by comparing each field.'''

        for x in ("qname", "seq", "flag",
                  "rname", "pos", "mapq", "cigar",
                  "mrnm", "mpos", "isize", "qual",
                  "bin",
                  "is_paired", "is_proper_pair",
                  "is_unmapped", "mate_is_unmapped",
                  "is_reverse", "mate_is_reverse",
                  "is_read1", "is_read2",
                  "is_secondary", "is_qcfail",
                  "is_duplicate"):
            if x in exclude: continue
            self.assertEqual( getattr(read1, x), getattr(read2,x), "attribute mismatch for %s: %s != %s" % 
                              (x, getattr(read1, x), getattr(read2,x)))

    def setUp( self ):

        
        a = pysam.AlignedRead()
        a.qname = "read_28833_29006_6945"
        a.seq="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
        a.flag = 99
        a.rname = 0
        a.pos = 32
        a.mapq = 20
        a.cigar = ( (0,10), (2,1), (0,25) )
        a.mrnm = 0
        a.mpos=199
        a.isize=167
	a.qual="<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<"
	a.tags = ( ("NM", 1),
                   ("RG", "L1") )

        b = pysam.AlignedRead()
        b.qname = "read_28701_28881_323b"
        b.seq="ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA"
        b.flag = 147
        b.rname = 1
        b.pos = 87
        b.mapq = 30
        b.cigar = ( (0,35), )
        b.mrnm = 1
        b.mpos=499
        b.isize=412
	b.qual="<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<"
	b.tags = ( ("MF", 18),
                   ("RG", "L2") )

        self.reads = (a,b)

    def testSAMWholeFile( self ):
        
        tmpfilename = "tmp_%i.sam" % id(self)

        outfile = pysam.Samfile( tmpfilename, "wh", header = self.header )

        for x in self.reads: outfile.write( x )
        outfile.close()
        
        self.assertTrue( checkBinaryEqual( tmpfilename, self.samfile ),
                         "mismatch when construction SAM file, see %s %s" % (tmpfilename, self.samfile))
        
        os.unlink( tmpfilename )

    def testBAMPerRead( self ):
        '''check if individual reads are binary equal.'''
        infile = pysam.Samfile( self.bamfile, "rb")

        others = list(infile)
        for denovo, other in zip( others, self.reads):
            self.checkFieldEqual( other, denovo )
            self.assertEqual( other.compare( denovo ), 0 )

    def testSAMPerRead( self ):
        '''check if individual reads are binary equal.'''
        infile = pysam.Samfile( self.samfile, "r")

        others = list(infile)
        for denovo, other in zip( others, self.reads):
            self.checkFieldEqual( other, denovo )
            self.assertEqual( other.compare( denovo), 0 )
            
    def testBAMWholeFile( self ):
        
        tmpfilename = "tmp_%i.bam" % id(self)

        outfile = pysam.Samfile( tmpfilename, "wb", header = self.header )

        for x in self.reads: outfile.write( x )
        outfile.close()
        
        self.assertTrue( checkBinaryEqual( tmpfilename, self.bamfile ),
                         "mismatch when construction BAM file, see %s %s" % (tmpfilename, self.bamfile))
        
        os.unlink( tmpfilename )


class TestDoubleFetch(unittest.TestCase):
    '''check if two iterators on the same bamfile are independent.'''
    
    def testDoubleFetch( self ):

        samfile1 = pysam.Samfile('ex1.bam', 'rb')

        for a,b in zip(samfile1.fetch(), samfile1.fetch()):
            self.assertEqual( a.compare( b ), 0 )

    def testDoubleFetchWithRegion( self ):

        samfile1 = pysam.Samfile('ex1.bam', 'rb')
        chr, start, stop = 'chr1', 200, 3000000
        self.assertTrue(len(list(samfile1.fetch ( chr, start, stop))) > 0) #just making sure the test has something to catch

        for a,b in zip(samfile1.fetch( chr, start, stop), samfile1.fetch( chr, start, stop)):
            self.assertEqual( a.compare( b ), 0 ) 

    def testDoubleFetchUntilEOF( self ):

        samfile1 = pysam.Samfile('ex1.bam', 'rb')

        for a,b in zip(samfile1.fetch( until_eof = True), 
                       samfile1.fetch( until_eof = True )):
            self.assertEqual( a.compare( b), 0 )

class TestRemoteFileFTP(unittest.TestCase):
    '''test remote access.

    '''

    # Need to find an ftp server without password on standard
    # port.

    url = "ftp://ftp.sanger.ac.uk/pub/rd/humanSequences/CV.bam"
    region = "1:1-1000"

    def testFTPView( self ):
        return
        result = pysam.view( self.url, self.region )
        self.assertEqual( len(result), 36 )
        
    def testFTPFetch( self ):
        return
        samfile = pysam.Samfile(self.url, "rb")  
        result = list(samfile.fetch( region = self.region ))
        self.assertEqual( len(result), 36 )

class TestRemoteFileHTTP( unittest.TestCase):

    url = "http://genserv.anat.ox.ac.uk/downloads/pysam/test/ex1.bam"
    region = "chr1:1-1000"
    local = "ex1.bam"

    def testView( self ):
        samfile_local = pysam.Samfile(self.local, "rb")  
        ref = list(samfile_local.fetch( region = self.region ))
        
        result = pysam.view( self.url, self.region )
        self.assertEqual( len(result), len(ref) )
        
    def testFetch( self ):
        samfile = pysam.Samfile(self.url, "rb")  
        result = list(samfile.fetch( region = self.region ))
        samfile_local = pysam.Samfile(self.local, "rb")  
        ref = list(samfile_local.fetch( region = self.region ))

        self.assertEqual( len(ref), len(result) )
        for x, y in zip(result, ref):
            self.assertEqual( x.compare( y ), 0 )

    def testFetchAll( self ):
        samfile = pysam.Samfile(self.url, "rb")  
        result = list(samfile.fetch())
        samfile_local = pysam.Samfile(self.local, "rb")  
        ref = list(samfile_local.fetch() )

        self.assertEqual( len(ref), len(result) )
        for x, y in zip(result, ref):
            self.assertEqual( x.compare( y ), 0 )

class TestLargeOptValues( unittest.TestCase ):

    ints = ( 65536, 214748, 2147484, 2147483647 )
    floats = ( 65536.0, 214748.0, 2147484.0 )

    def check( self, samfile ):
        
        i = samfile.fetch()
        for exp in self.ints:
            rr = i.next()
            obs = rr.opt("ZP")
            self.assertEqual( exp, obs, "expected %s, got %s\n%s" % (str(exp), str(obs), str(rr)))

        for exp in [ -x for x in self.ints ]:
            rr = i.next()
            obs = rr.opt("ZP")
            self.assertEqual( exp, obs, "expected %s, got %s\n%s" % (str(exp), str(obs), str(rr)))

        for exp in self.floats:
            rr = i.next()
            obs = rr.opt("ZP")
            self.assertEqual( exp, obs, "expected %s, got %s\n%s" % (str(exp), str(obs), str(rr)))

        for exp in [ -x for x in self.floats ]:
            rr = i.next()
            obs = rr.opt("ZP")
            self.assertEqual( exp, obs, "expected %s, got %s\n%s" % (str(exp), str(obs), str(rr)))

    def testSAM( self ):
        samfile = pysam.Samfile("ex10.sam", "r")
        self.check( samfile )

    def testBAM( self ):
        samfile = pysam.Samfile("ex10.bam", "rb")
        self.check( samfile )

# class TestSNPCalls( unittest.TestCase ):
#     '''test pysam SNP calling ability.'''

#     def checkEqual( self, a, b ):
#         for x in ("reference_base", 
#                   "pos",
#                   "genotype",
#                   "consensus_quality",
#                   "snp_quality",
#                   "mapping_quality",
#                   "coverage" ):
#             self.assertEqual( getattr(a, x), getattr(b,x), "%s mismatch: %s != %s\n%s\n%s" % 
#                               (x, getattr(a, x), getattr(b,x), str(a), str(b)))

#     def testAllPositionsViaIterator( self ):
#         samfile = pysam.Samfile( "ex1.bam", "rb")  
#         fastafile = pysam.Fastafile( "ex1.fa" )
#         try: 
#             refs = [ x for x in pysam.pileup( "-c", "-f", "ex1.fa", "ex1.bam" ) if x.reference_base != "*"]
#         except pysam.SamtoolsError:
#             pass

#         i = samfile.pileup( stepper = "samtools", fastafile = fastafile )
#         calls = list(pysam.IteratorSNPCalls(i))
#         for x,y in zip( refs, calls ):
#             self.checkEqual( x, y )

#     def testPerPositionViaIterator( self ):
#         # test pileup for each position. This is a slow operation
#         # so this test is disabled 
#         return
#         samfile = pysam.Samfile( "ex1.bam", "rb")  
#         fastafile = pysam.Fastafile( "ex1.fa" )
#         for x in pysam.pileup( "-c", "-f", "ex1.fa", "ex1.bam" ):
#             if x.reference_base == "*": continue
#             i = samfile.pileup( x.chromosome, x.pos, x.pos+1,
#                                 fastafile = fastafile,
#                                 stepper = "samtools" )
#             z = [ zz for zz in pysam.IteratorSamtools(i) if zz.pos == x.pos ]
#             self.assertEqual( len(z), 1 )
#             self.checkEqual( x, z[0] )

#     def testPerPositionViaCaller( self ):
#         # test pileup for each position. This is a fast operation
#         samfile = pysam.Samfile( "ex1.bam", "rb")  
#         fastafile = pysam.Fastafile( "ex1.fa" )
#         i = samfile.pileup( stepper = "samtools", fastafile = fastafile )
#         caller = pysam.SNPCaller( i )

#         for x in pysam.pileup( "-c", "-f", "ex1.fa", "ex1.bam" ):
#             if x.reference_base == "*": continue
#             call = caller.call( x.chromosome, x.pos )
#             self.checkEqual( x, call )

# class TestIndelCalls( unittest.TestCase ):
#     '''test pysam indel calling.'''

#     def checkEqual( self, a, b ):

#         for x in ("pos",
#                   "genotype",
#                   "consensus_quality",
#                   "snp_quality",
#                   "mapping_quality",
#                   "coverage",
#                   "first_allele",
#                   "second_allele",
#                   "reads_first",
#                   "reads_second",
#                   "reads_diff"):
#             if b.genotype == "*/*" and x == "second_allele":
#                 # ignore test for second allele (positions chr2:439 and chr2:1512)
#                 continue
#             self.assertEqual( getattr(a, x), getattr(b,x), "%s mismatch: %s != %s\n%s\n%s" % 
#                               (x, getattr(a, x), getattr(b,x), str(a), str(b)))

#     def testAllPositionsViaIterator( self ):

#         samfile = pysam.Samfile( "ex1.bam", "rb")  
#         fastafile = pysam.Fastafile( "ex1.fa" )
#         try: 
#             refs = [ x for x in pysam.pileup( "-c", "-f", "ex1.fa", "ex1.bam" ) if x.reference_base == "*"]
#         except pysam.SamtoolsError:
#             pass

#         i = samfile.pileup( stepper = "samtools", fastafile = fastafile )
#         calls = [ x for x in list(pysam.IteratorIndelCalls(i)) if x != None ]
#         for x,y in zip( refs, calls ):
#             self.checkEqual( x, y )

#     def testPerPositionViaCaller( self ):
#         # test pileup for each position. This is a fast operation
#         samfile = pysam.Samfile( "ex1.bam", "rb")  
#         fastafile = pysam.Fastafile( "ex1.fa" )
#         i = samfile.pileup( stepper = "samtools", fastafile = fastafile )
#         caller = pysam.IndelCaller( i )

#         for x in pysam.pileup( "-c", "-f", "ex1.fa", "ex1.bam" ):
#             if x.reference_base != "*": continue
#             call = caller.call( x.chromosome, x.pos )
#             self.checkEqual( x, call )

class TestLogging( unittest.TestCase ):
    '''test around bug issue 42,

    failed in versions < 0.4
    '''

    def check( self, bamfile, log ):

        if log:
            logger = logging.getLogger('franklin')
            logger.setLevel(logging.INFO)
            formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
            log_hand  = logging.FileHandler('log.txt')
            log_hand.setFormatter(formatter)
            logger.addHandler(log_hand)

        bam  = pysam.Samfile(bamfile, 'rb')
        cols = bam.pileup()
        self.assert_( True )

    def testFail1( self ):
        self.check( "ex9_fail.bam", False )
        self.check( "ex9_fail.bam", True )

    def testNoFail1( self ):
        self.check( "ex9_nofail.bam", False )
        self.check( "ex9_nofail.bam", True )

    def testNoFail2( self ):
        self.check( "ex9_nofail.bam", True )
        self.check( "ex9_nofail.bam", True )
        
# TODOS
# 1. finish testing all properties within pileup objects
# 2. check exceptions and bad input problems (missing files, optional fields that aren't present, etc...)
# 3. check: presence of sequence

class TestSamfileUtilityFunctions( unittest.TestCase ):

    def testCount( self ):

        samfile = pysam.Samfile( "ex1.bam", "rb" )

        for contig in ("chr1", "chr2" ):
            for start in xrange( 0, 2000, 100 ):
                end = start + 1
                self.assertEqual( len( list( samfile.fetch( contig, start, end ) ) ),
                                  samfile.count( contig, start, end ) )

                # test empty intervals
                self.assertEqual( len( list( samfile.fetch( contig, start, start ) ) ),
                                  samfile.count( contig, start, start ) )

                # test half empty intervals
                self.assertEqual( len( list( samfile.fetch( contig, start ) ) ),
                                  samfile.count( contig, start ) )

    def testMate( self ):
        '''test mate access.'''

        readnames = [ x.split("\t")[0] for x in open( "ex1.sam", "rb" ).readlines() ]
        counts = collections.defaultdict( int )
        for x in readnames: counts[x] += 1

        samfile = pysam.Samfile( "ex1.bam", "rb" )
        for read in samfile.fetch():
            if not read.is_paired:
                self.assertRaises( ValueError, samfile.mate, read )
            elif read.mate_is_unmapped:
                self.assertRaises( ValueError, samfile.mate, read )
            else:
                if counts[read.qname] == 1:
                    self.assertRaises( ValueError, samfile.mate, read )
                else:
                    mate = samfile.mate( read )
                    self.assertEqual( read.qname, mate.qname )
                    self.assertEqual( read.is_read1, mate.is_read2 )
                    self.assertEqual( read.is_read2, mate.is_read1 )
                    self.assertEqual( read.pos, mate.mpos )
                    self.assertEqual( read.mpos, mate.pos )

    def testIndexStats( self ):
        '''test if total number of mapped/unmapped reads is correct.'''

        samfile = pysam.Samfile( "ex1.bam", "rb" )
        self.assertEqual( samfile.mapped, 3235 )
        self.assertEqual( samfile.unmapped, 35 )

class TestSamtoolsProxy( unittest.TestCase ):
    '''tests for sanity checking access to samtools functions.'''

    def testIndex( self ):
        self.assertRaises( IOError, pysam.index, "missing_file" )

    def testView( self ):
        # note that view still echos "open: No such file or directory"
        self.assertRaises( pysam.SamtoolsError, pysam.view, "missing_file" )

    def testSort( self ):
        self.assertRaises( pysam.SamtoolsError, pysam.sort, "missing_file" )

class TestSamfileIndex( unittest.TestCase):
    
    def testIndex( self ):
        samfile = pysam.Samfile( "ex1.bam", "rb" )
        index = pysam.IndexedReads( samfile )
        index.build()

        reads = collections.defaultdict( int )

        for read in samfile: reads[read.qname] += 1
            
        for qname, counts in reads.iteritems():
            found = list(index.find( qname ))
            self.assertEqual( len(found), counts )
            for x in found: self.assertEqual( x.qname, qname )
            

if __name__ == "__main__":
    # build data files
    print "building data files"
    subprocess.call( "make", shell=True)
    print "starting tests"
    unittest.main()
