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
        }

    # some tests depend on others. The order specifies in which order
    # the samtools commands are executed.
    mOrder = ('faidx', 'import', 'index', 'pileup1', 'pileup2', 'glfview' )

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
                self.runSamtools( samtools_command )
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
            self.assertTrue( self.checkBinaryEqual( samtools_target, pysam_target ), 
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

    def __del__(self):

        for label, command in self.mCommands.iteritems():
            samtools_target, samtools_command = command[0]
            pysam_target, pysam_command = command[1]
            if os.path.exists( samtools_target): os.remove( samtools_target )
            if os.path.exists( pysam_target): os.remove( pysam_target )
        if os.path.exists( "pysam_ex1.fa" ): os.remove( "pysam_ex1.fa" )

if __name__ == "__main__":
    unittest.main()
