import glob
import os
import unittest
import pysam

DATADIR="cbcf_data"
from tabix_test import loadAndConvert

class TestMissingGenotypes(unittest.TestCase):

    filename = "missing_genotypes.vcf"

    def setUp(self):
        self.compare = loadAndConvert(os.path.join(DATADIR, self.filename),
                                      encode=False)

    def check(self, filename):
        """see issue 203 - check for segmentation fault"""
        fn = os.path.join(DATADIR, filename)
        self.assertEqual(True, os.path.exists(fn))
        v = pysam.VariantFile(fn)
        for site in v:
            for ss,rec in site.samples.items():
                a, b = ss, rec

        v = pysam.VariantFile(fn)
        for x, site in enumerate(v):
            for ss,rec in site.samples.items():
                a, b = ss, rec.alleles
                a, b = ss, rec.allele_indices

    def testVCF(self):
        self.check(self.filename)

    def testVCFGZ(self):
        self.check(self.filename + ".gz")


class TestOpening(unittest.TestCase):

    def testMissingFile(self):
        self.assertRaises(IOError, pysam.VariantFile,
                          "missing_file.vcf")

    def testMissingFileVCFGZ(self):
        self.assertRaises(IOError, pysam.VariantFile,
                          "missing_file.vcf.gz")
        
    def testEmptyFileVCF(self):
        with open("tmp_testEmptyFile.vcf", "w"):
            pass

        self.assertRaises(ValueError, pysam.VariantFile,
                          "tmp_testEmptyFile.vcf")

        os.unlink("tmp_testEmptyFile.vcf")

    def testEmptyFileVCFGZ(self):
        with open("tmp_testEmptyFile.vcf", "w"):
            pass

        pysam.tabix_compress("tmp_testEmptyFile.vcf",
                             "tmp_testEmptyFile.vcf.gz")

        self.assertRaises(ValueError, pysam.VariantFile,
                          "tmp_testEmptyFile.vcf.gz")

        os.unlink("tmp_testEmptyFile.vcf")
        os.unlink("tmp_testEmptyFile.vcf.gz")

if __name__ == "__main__":
    unittest.main()
