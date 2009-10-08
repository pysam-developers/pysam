#!/usr/bin/env python
'''unit testing code for pysam.'''

import pysam
import unittest
import os
import itertools
import subprocess
import shutil

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

class BinaryTest(unittest.TestCase):
    '''test samtools command line commands and compare
    against pysam commands.

    Tests fail, if the output is not binary identical.
    '''

    first_time = True

    # a list of commands to test
    mCommands = \
        { "faidx" : \
        ( 
            ("ex1.fa.fai", "samtools faidx ex1.fa"), 
            ("pysam_ex1.fa.fai", (pysam.faidx, "ex1.fa") ),
            ),
          "import" :
              (
                ("ex1.bam", "samtools import ex1.fa.fai ex1.sam.gz ex1.bam" ),
                ("pysam_ex1.bam", (pysam.samimport, "ex1.fa.fai ex1.sam.gz pysam_ex1.bam") ),
                ),
          "index":
              (
                ("ex1.bam.bai", "samtools index ex1.bam" ),
                ("pysam_ex1.bam.bai", (pysam.index, "pysam_ex1.bam" ) ),
                ),
          "pileup1" :
          (
                ("ex1.pileup", "samtools pileup -cf ex1.fa ex1.bam > ex1.pileup" ),
                ("pysam_ex1.pileup", (pysam.pileup, "-c -f ex1.fa ex1.bam" ) )
                ),
          "pileup2" :
          (
                ("ex1.glf", "samtools pileup -gf ex1.fa ex1.bam > ex1.glf" ),
                ("pysam_ex1.glf", (pysam.pileup, "-g -f ex1.fa ex1.bam" ) )
                ),
          "glfview" :
        (
                ("ex1.glfview", "samtools glfview ex1.glf > ex1.glfview"),
                ("pysam_ex1.glfview", (pysam.glfview, "ex1.glf" ) ),
                ),
          "view" :
        (
                ("ex1.view", "samtools view ex1.bam > ex1.view"),
                ("pysam_ex1.view", (pysam.view, "ex1.bam" ) ),
                ),
        }

    # some tests depend on others. The order specifies in which order
    # the samtools commands are executed.
    mOrder = ('faidx', 'import', 'index', 'pileup1', 'pileup2', 'glfview', 'view' )

    def setUp( self ):
        '''setup tests. 

        For setup, all commands will be run before the first test is
        executed. Individual tests will then just compare the output
        files.
        '''

        if BinaryTest.first_time:
            # copy the source 
            shutil.copy( "ex1.fa", "pysam_ex1.fa" )

            for label in self.mOrder:
                command = self.mCommands[label]
                samtools_target, samtools_command = command[0]
                pysam_target, pysam_command = command[1]
                runSamtools( samtools_command )
                pysam_method, pysam_options = pysam_command
                output = pysam_method( *pysam_options.split(" "), raw=True)
                if ">" in samtools_command:
                    outfile = open( pysam_target, "w" )
                    for line in output: outfile.write( line )
                    outfile.close()

            BinaryTest.first_time = False

    def checkCommand( self, command ):
        if command:
            samtools_target, pysam_target = self.mCommands[command][0][0], self.mCommands[command][1][0]
            self.assertTrue( checkBinaryEqual( samtools_target, pysam_target ), 
                             "%s failed: files %s and %s are not the same" % (command, samtools_target, pysam_target) )

    def testImport( self ):
        self.checkCommand( "import" )

    def testIndex( self ):
        self.checkCommand( "index" )
        
    def testPileup1( self ):
        self.checkCommand( "pileup1" )
    
    def testPileup2( self ):
        self.checkCommand( "pileup2" )

    def testGLFView( self ):
        self.checkCommand( "glfview" )

    def testView( self ):
        self.checkCommand( "view" )

    def __del__(self):

        for label, command in self.mCommands.iteritems():
            samtools_target, samtools_command = command[0]
            pysam_target, pysam_command = command[1]
            if os.path.exists( samtools_target): os.remove( samtools_target )
            if os.path.exists( pysam_target): os.remove( pysam_target )
        if os.path.exists( "pysam_ex1.fa" ): os.remove( "pysam_ex1.fa" )

class IOTest(unittest.TestCase):
    '''check if reading samfile and writing a samfile are consistent.'''

    def checkEcho( self, input_filename, reference_filename, 
                   output_filename, 
                   input_mode, output_mode, use_template = True):
        '''iterate through *input_filename* writing to *output_filename* and
        comparing the output to *reference_filename*. 
        
        The files are opened according to the *input_mode* and *output_mode*.

        If *use_template* is set, the header is copied from infile using the
        template mechanism, otherwise target names and lengths are passed explicitely. 
        '''

        infile = pysam.Samfile( input_filename, input_mode )
        if use_template:
            outfile = pysam.Samfile( output_filename, output_mode, template = infile )
        else:
            outfile = pysam.Samfile( output_filename, output_mode, 
                                     targetnames = infile.targets,
                                     targetlengths = infile.lengths )

        iter = pysam.IteratorRowAll( infile )
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

class TestIteratorRow(unittest.TestCase):

    def setUp(self):
        self.samfile=pysam.Samfile( "ex1.bam","rb" )

    def checkRange( self, rnge ):
        '''compare results from iterator with those from samtools.'''
        ps = list(pysam.IteratorRow(self.samfile, rnge))
        sa = list(pysam.view( "ex1.bam", rnge , raw = True) )
        self.assertEqual( len(ps), len(sa), "unequal number of results for range %s: %i != %i" % (rnge, len(ps), len(sa) ))
        # check if the same reads are returned and in the same order
        for line, pair in enumerate( zip( ps, sa ) ):
            data = pair[1].split("\t")
            self.assertEqual( pair[0].qname, data[0], "read id mismatch in line %i: %s != %s" % (line, pair[0].rname, data[0]) )

    def testIteratePerContig(self):
        '''check random access per contig'''
        for contig in self.samfile.targets:
            self.checkRange( contig )

    def testIterateRanges(self):
        '''check random access per range'''
        for contig, length in zip(self.samfile.targets, self.samfile.lengths):
            for start in range( 0, length, 100):
                self.checkRange( "%s:%i-%i" % (contig, start, start + 100) )

    def tearDown(self):
        self.samfile.close()

class TestIteratorRowAll(unittest.TestCase):

    def setUp(self):
        self.samfile=pysam.Samfile( "ex1.bam","rb" )

    def testIterate(self):
        '''compare results from iterator with those from samtools.'''
        ps = list(pysam.IteratorRowAll(self.samfile))
        sa = list(pysam.view( "ex1.bam", raw = True) )
        self.assertEqual( len(ps), len(sa), "unequal number of results: %i != %i" % (len(ps), len(sa) ))
        # check if the same reads are returned
        for line, pair in enumerate( zip( ps, sa ) ):
            data = pair[1].split("\t")
            self.assertEqual( pair[0].qname, data[0], "read id mismatch in line %i: %s != %s" % (line, pair[0].rname, data[0]) )

    def tearDown(self):
        self.samfile.close()

if __name__ == "__main__":
    unittest.main()
