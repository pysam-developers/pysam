#!/usr/bin/env python
'''unit testing code for pysam.'''

import pysam
import unittest
import os
import itertools
import subprocess
import shutil

class TestExceptions(unittest.TestCase):

    def setUp(self):
        self.samfile=pysam.Samfile( "ex1.bam","rb" )

    def testOutOfRangeNegativeNewFormat(self):
        self.assertRaises( ValueError, self.samfile.fetch, "chr1", 5, -10 )
        self.assertRaises( ValueError, self.samfile.fetch, "chr1", 5, 0 )
        self.assertRaises( ValueError, self.samfile.fetch, "chr1", -5, -10 )

    def testOutOfRangeNegativeOldFormat(self):
        self.assertRaises( ValueError, self.samfile.fetch, "chr1:-5-10" )
        self.assertRaises( ValueError, self.samfile.fetch, "chr1:-5-0" )
        self.assertRaises( ValueError, self.samfile.fetch, "chr1:-5--10" )

    def testOutOfRangeLargeNewFormat(self):
        self.assertRaises( ValueError, self.samfile.fetch, "chr1", 99999999999999999, 999999999999999999 )

    def testOutOfRangeLargeOldFormat(self):
        self.assertRaises( ValueError, self.samfile.fetch, "chr1:99999999999999999-999999999999999999" )

    def tearDown(self):
        self.samfile.close()

if __name__ == "__main__":
    unittest.main()

