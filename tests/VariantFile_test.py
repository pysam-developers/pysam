import glob
import os
import unittest
import pysam

DATADIR="cbcf_data"
from tabix_test import loadAndConvert


# class TestMissingGenotypes(unittest.TestCase):

#     filename = "missing_genotypes.vcf"

#     def setUp(self):
#         self.compare = loadAndConvert(os.path.join(DATADIR, self.filename),
#                                       encode=False)

#     def check(self, filename):
#         return 
#         fn = os.path.join(DATADIR, filename)
#         self.assertEqual(True, os.path.exists(fn))
#         # v = pysam.VariantFile(fn)
#         # for site in v:
#         #     for ss,rec in site.samples.items():
#         #         print ss, rec

#         v = pysam.VariantFile(fn)
#         for x, site in enumerate(v):
#             import pdb; pdb.set_trace()
#             print "======"
#             print site
#             for ss,rec in site.samples.items():
#                 print ss
#                 print rec.alleles

#     def testVCF(self):
#         self.check(self.filename)

#     def testVCFGZ(self):
#         self.check(self.filename + ".gz")


if __name__ == "__main__":
    unittest.main()
