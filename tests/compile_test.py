'''
compile_test.py - check pyximport functionality with pysam
==========================================================

test script for checking if compilation against
pysam and tabix works.
'''

# clean up previous compilation
import os
import unittest
import pysam
from TestUtils import BAM_DATADIR, TABIX_DATADIR

try:
    os.unlink('tests/_compile_test.c')
    os.unlink('tests/_compile_test.pyxbldc')
except OSError:
    pass

import pyximport
pyximport.install(build_in_temp=False)
import _compile_test


class BAMTest(unittest.TestCase):

    input_filename = os.path.join(BAM_DATADIR, "ex1.bam")

    def testCount(self):

        nread = _compile_test.testCountBAM(
            pysam.Samfile(self.input_filename))
        self.assertEqual(nread, 3270)


class GTFTest(unittest.TestCase):

    input_filename = os.path.join(TABIX_DATADIR, "example.gtf.gz")

    def testCount(self):
        nread = _compile_test.testCountGTF(
            pysam.Tabixfile(self.input_filename))
        self.assertEqual(nread, 237)


if __name__ == "__main__":
    unittest.main()
