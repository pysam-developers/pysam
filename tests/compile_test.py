'''
compile_test.py - check pyximport
=================================

test script for checking if compilation against
pysam and tabix works.
'''
# clean up previous compilation
import os
try:
    os.unlink('_compile_test.c')
    os.unlink('_compile_test.pyxbldc')
except OSError:
    pass


import pyximport
pyximport.install(build_in_temp=False)
import _compile_test

import unittest
import pysam


class BAMTest(unittest.TestCase):

    input_filename = "pysam_data/ex1.bam"

    def testCount(self):

        nread = _compile_test.testCountBAM(
            pysam.Samfile(self.input_filename))
        self.assertEqual(nread, 3270)


class GTFTest(unittest.TestCase):

    input_filename = "tabix_data/example.gtf.gz"

    def testCount(self):
        nread = _compile_test.testCountGTF(
            pysam.Tabixfile(self.input_filename))
        self.assertEqual(nread, 237)

if __name__ == "__main__":
    unittest.main()
