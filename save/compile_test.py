'''test script for checking if compilation against
pysam and tabix works.'''

import pyximport
pyximport.install( build_in_temp=False)
import _compile_test

import unittest

if __name__ == "__main__":
    unittest.main()
