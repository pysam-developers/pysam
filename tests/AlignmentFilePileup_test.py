"""Benchmarking module for AlignmentFile functionality"""
import os
import unittest
from TestUtils import BAM_DATADIR, force_str, flatten_nested_list
from PileupTestUtils import *


class PileUpColumnTests(unittest.TestCase):

    fn = os.path.join(BAM_DATADIR, "ex2.bam")

    def test_pileup_depths_are_equal(self):
        samtools_result = build_depth_with_samtoolspipe(self.fn)
        pysam_result = build_depth_with_filter_with_pysam(self.fn)
        self.assertEqual(pysam_result, samtools_result)

    def test_pileup_query_bases_are_equal(self):
        samtools_result = build_query_bases_with_samtoolspipe(self.fn)
        pysam_result = build_query_bases_with_pysam(self.fn)
        self.assertEqual(["".join(x) for x in pysam_result], samtools_result)

    def test_pileup_query_qualities_are_equal(self):
        samtools_result = build_query_qualities_with_samtoolspipe(self.fn)
        pysam_result = build_query_qualities_with_pysam(self.fn)
        pysam_result = [
            [chr(min(126, x + 33)) for x in l] for l in pysam_result]
        self.assertEqual("".join(flatten_nested_list(pysam_result)),
                         "".join(flatten_nested_list(samtools_result)))

    def test_pileup_mapping_qualities_are_equal(self):
        samtools_result = build_mapping_qualities_with_samtoolspipe(self.fn)
        pysam_result = build_mapping_qualities_with_pysam(self.fn)
        # convert to chars
        pysam_result = [
            [chr(min(126, x + 33)) for x in l] for l in pysam_result]
                
        self.assertEqual("".join(flatten_nested_list(pysam_result)),
                         "".join(flatten_nested_list(samtools_result)))
        
        
        

if __name__ == "__main__":
    unittest.main()
