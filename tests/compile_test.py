'''test script for checking if compilation against
pysam and tabix works.'''

import pyximport
pyximport.install( build_in_temp = False )
import _compile_test

import unittest
import pysam

class BAMTest(unittest.TestCase):

    input_filename = "ex1.bam"
    
    def testCount( self ):
        
        nread = _compile_test.testCountBAM( pysam.Samfile( self.input_filename ) )
        self.assertEqual( nread, 3270 )

class GTFTest(unittest.TestCase):

    input_filename = "example.gtf.gz"
    
    def testCount( self ):
        
        nread = _compile_test.testCountGTF( pysam.Tabixfile( self.input_filename ) )


if __name__ == "__main__":
    unittest.main()
