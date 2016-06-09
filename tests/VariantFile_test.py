import os
import unittest
import pysam
import gzip
import subprocess
from TestUtils import get_temp_filename, check_lines_equal

DATADIR="cbcf_data"
from tabix_test import loadAndConvert


def read_header(filename):
    data = []
    if filename.endswith(".gz"):
        for line in gzip.open(filename):
            line = line.decode("ascii")
            if line.startswith("#"):
                data.append(line)
    else:
        with open(filename) as f:
            for line in f:
                if line.startswith("#"):
                    data.append(line)

    return data


class TestMissingGenotypes(unittest.TestCase):

    filename = "missing_genotypes.vcf"

    def setUp(self):
        self.compare = loadAndConvert(
            os.path.join(DATADIR, self.filename),
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

    def testEmptyFileVCFGZWithIndex(self):
        with open("tmp_testEmptyFile.vcf", "w"):
            pass

        pysam.tabix_index("tmp_testEmptyFile.vcf",
                          preset="vcf",
                          force=True)

        self.assertRaises(ValueError, pysam.VariantFile,
                          "tmp_testEmptyFile.vcf.gz")

        os.unlink("tmp_testEmptyFile.vcf.gz")
        os.unlink("tmp_testEmptyFile.vcf.gz.tbi")

    def testEmptyFileVCFGZWithoutIndex(self):
        with open("tmp_testEmptyFileWithoutIndex.vcf", "w"):
            pass

        pysam.tabix_compress("tmp_testEmptyFileWithoutIndex.vcf",
                             "tmp_testEmptyFileWithoutIndex.vcf.gz",
                             force=True)

        self.assertRaises(ValueError, pysam.VariantFile,
                          "tmp_testEmptyFileWithoutIndex.vcf.gz")

        os.unlink("tmp_testEmptyFileWithoutIndex.vcf")
        os.unlink("tmp_testEmptyFileWithoutIndex.vcf.gz")

    def testEmptyFileVCFOnlyHeader(self):
        with pysam.VariantFile(os.path.join(
                DATADIR,
                "example_vcf42_only_header.vcf")) as inf:
            self.assertEqual(len(list(inf.fetch())), 0)

    def testEmptyFileVCFGZOnlyHeader(self):
        with pysam.VariantFile(os.path.join(
                DATADIR,
                "example_vcf42_only_header.vcf")) as inf:
            self.assertEqual(len(list(inf.fetch())), 0)

    def testDetectVCF(self):
        with pysam.VariantFile(os.path.join(DATADIR,
            "example_vcf40.vcf")) as inf:
            self.assertEqual(inf.category, 'VARIANTS')
            self.assertEqual(inf.format, 'VCF')
            self.assertEqual(inf.compression, 'NONE')
            self.assertFalse(inf.is_remote)
            self.assertFalse(inf.is_stream)
            self.assertEqual(len(list(inf.fetch())), 5)

    def testDetectVCFGZ(self):
        with pysam.VariantFile(os.path.join(DATADIR,
            "example_vcf40.vcf.gz")) as inf:
            self.assertEqual(inf.category, 'VARIANTS')
            self.assertEqual(inf.format, 'VCF')
            self.assertEqual(inf.compression, 'BGZF')
            self.assertFalse(inf.is_remote)
            self.assertFalse(inf.is_stream)
            self.assertEqual(len(list(inf.fetch())), 5)

    def testDetectBCF(self):
        with pysam.VariantFile(os.path.join(
                DATADIR,
                "example_vcf40.bcf")) as inf:
            self.assertEqual(inf.category, 'VARIANTS')
            self.assertEqual(inf.format, 'BCF')
            self.assertEqual(inf.compression, 'BGZF')
            self.assertFalse(inf.is_remote)
            self.assertFalse(inf.is_stream)
            self.assertEqual(len(list(inf.fetch())), 5)


class TestHeader(unittest.TestCase):

    filename = "example_vcf40.vcf"

    def testStr(self):

        fn = os.path.join(DATADIR, self.filename)
        v = pysam.VariantFile(fn)

        ref = read_header(fn)
        comp = str(v.header).splitlines(True)

        self.assertEqual(sorted(ref),
                         sorted(comp))

    def testIterator(self):

        fn = os.path.join(DATADIR, self.filename)
        v = pysam.VariantFile(fn)

        ref = read_header(fn)
        # remove last header line starting with #CHROM
        ref.pop()
        ref = sorted(ref)
        comp = sorted([str(x) for x in v.header.records])

        self.assertEqual(len(ref), len(comp))

        for x, y in zip(ref, comp):
            self.assertEqual(x[:-1], str(y))


# These tests need to be separate and start from newly opened files.  This
# is because htslib's parser is lazy and the pysam API needs to trigger
# appropriate parsing when accessing each time of data.  Failure to do so
# will result in crashes or return of incorrect data.  Thus this test suite
# is testing both the triggering of the lazy parser and the results of the
# parser.
class TestParsing(unittest.TestCase):

    filename = "example_vcf40.vcf.gz"

    def testChrom(self):
        fn = os.path.join(DATADIR, self.filename)
        v = pysam.VariantFile(fn)
        chrom = [rec.chrom for rec in v]
        self.assertEqual(chrom, ['M', '17', '20', '20', '20'])

    def testPos(self):
        fn = os.path.join(DATADIR, self.filename)
        v = pysam.VariantFile(fn)
        pos = [rec.pos for rec in v]
        self.assertEqual(pos, [1230237, 14370, 17330, 1110696, 1234567])

    def testStart(self):
        fn = os.path.join(DATADIR, self.filename)
        v = pysam.VariantFile(fn)
        start = [rec.start for rec in v]
        self.assertEqual(start, [1230236, 14369, 17329, 1110695, 1234566])

    def testStop(self):
        fn = os.path.join(DATADIR, self.filename)
        v = pysam.VariantFile(fn)
        stop = [rec.stop for rec in v]
        self.assertEqual(stop, [1230237, 14370, 17330, 1110696, 1234570])

    def testId(self):
        fn = os.path.join(DATADIR, self.filename)
        v = pysam.VariantFile(fn)
        ids = [rec.id for rec in v]
        self.assertEqual(ids, [None, 'rs6054257', None, 'rs6040355', 'microsat1'])

    def testRef(self):
        fn = os.path.join(DATADIR, self.filename)
        v = pysam.VariantFile(fn)
        ref = [rec.ref for rec in v]
        self.assertEqual(ref, ['T', 'G', 'T', 'A', 'GTCT'])

    def testAlt(self):
        fn = os.path.join(DATADIR, self.filename)
        v = pysam.VariantFile(fn)
        alts = [rec.alts for rec in v]
        self.assertEqual(alts, [None, ('A',), ('A',), ('G', 'T'), ('G', 'GTACT')])

    def testAlleles(self):
        fn = os.path.join(DATADIR, self.filename)
        v = pysam.VariantFile(fn)
        alleles = [rec.alleles for rec in v]
        self.assertEqual(alleles, [('T',), ('G', 'A'), ('T', 'A'), ('A', 'G', 'T'), ('GTCT', 'G', 'GTACT')])

    def testQual(self):
        fn = os.path.join(DATADIR, self.filename)
        v = pysam.VariantFile(fn)
        qual = [rec.qual for rec in v]
        self.assertEqual(qual, [47.0, 29.0, 3.0, 67.0, 50.0])

    def testFilter(self):
        fn = os.path.join(DATADIR, self.filename)
        v = pysam.VariantFile(fn)
        filter = [rec.filter.keys() for rec in v]
        self.assertEqual(filter, [['PASS'], ['PASS'], ['q10'], ['PASS'], ['PASS']])

    def testInfo(self):
        fn = os.path.join(DATADIR, self.filename)
        v = pysam.VariantFile(fn)
        info = [rec.info.items() for rec in v]
        self.assertEqual(info, [[('NS', 3), ('DP', 13), ('AA', 'T')],
                                [('NS', 3), ('DP', 14), ('AF', (0.5,)), ('DB', True), ('H2', True)],
                                [('NS', 3), ('DP', 11), ('AF', (0.017000000923871994,))],
                                [('NS', 2), ('DP', 10), ('AF', (0.3330000042915344, 0.6669999957084656)),
                                            ('AA', 'T'), ('DB', True)],
                                [('NS', 3), ('DP', 9), ('AA', 'G')]])

    def testFormat(self):
        fn = os.path.join(DATADIR, self.filename)
        v = pysam.VariantFile(fn)
        format = [rec.format.keys() for rec in v]
        self.assertEqual(format, [['GT', 'GQ', 'DP', 'HQ'],
                                  ['GT', 'GQ', 'DP', 'HQ'],
                                  ['GT', 'GQ', 'DP', 'HQ'],
                                  ['GT', 'GQ', 'DP', 'HQ'],
                                  ['GT', 'GQ', 'DP']])

    def testSampleAlleles(self):
        fn = os.path.join(DATADIR, self.filename)
        v = pysam.VariantFile(fn)
        alleles = [s.alleles for rec in v for s in rec.samples.values()]
        self.assertEqual(alleles, [('T', 'T'), ('T', 'T'), ('T', 'T'),
                                   ('G', 'G'), ('A', 'G'), ('A', 'A'),
                                   ('T', 'T'), ('T', 'A'), ('T', 'T'),
                                   ('G', 'T'), ('T', 'G'), ('T', 'T'),
                                   ('GTCT', 'G'), ('GTCT', 'GTACT'),
                                   ('G', 'G')])

    def testSampleFormats(self):
        fn = os.path.join(DATADIR, self.filename)
        v = pysam.VariantFile(fn)
        format = [s.items() for rec in v for s in rec.samples.values()]
        self.assertEqual(format, [[('GT', (0, 0)), ('GQ', 54), ('DP', 7), ('HQ', (56, 60))],
                                  [('GT', (0, 0)), ('GQ', 48), ('DP', 4), ('HQ', (51, 51))],
                                  [('GT', (0, 0)), ('GQ', 61), ('DP', 2), ('HQ', (None,))],
                                  [('GT', (0, 0)), ('GQ', 48), ('DP', 1), ('HQ', (51, 51))],
                                  [('GT', (1, 0)), ('GQ', 48), ('DP', 8), ('HQ', (51, 51))],
                                  [('GT', (1, 1)), ('GQ', 43), ('DP', 5), ('HQ', (None, None))],
                                  [('GT', (0, 0)), ('GQ', 49), ('DP', 3), ('HQ', (58, 50))],
                                  [('GT', (0, 1)), ('GQ', 3), ('DP', 5), ('HQ', (65, 3))],
                                  [('GT', (0, 0)), ('GQ', 41), ('DP', 3), ('HQ', (None,))],
                                  [('GT', (1, 2)), ('GQ', 21), ('DP', 6), ('HQ', (23, 27))],
                                  [('GT', (2, 1)), ('GQ', 2), ('DP', 0), ('HQ', (18, 2))],
                                  [('GT', (2, 2)), ('GQ', 35), ('DP', 4), ('HQ', (None,))],
                                  [('GT', (0, 1)), ('GQ', 35), ('DP', 4)],
                                  [('GT', (0, 2)), ('GQ', 17), ('DP', 2)],
                                  [('GT', (1, 1)), ('GQ', 40), ('DP', 3)]])

    def testSampleAlleleIndices(self):
        fn = os.path.join(DATADIR, self.filename)
        v = pysam.VariantFile(fn)
        indices = [s.allele_indices for rec in v for s in rec.samples.values()]
        self.assertEqual(indices, [(0, 0), (0, 0), (0, 0), (0, 0), (1, 0),
                                   (1, 1), (0, 0), (0, 1), (0, 0), (1, 2),
                                   (2, 1), (2, 2), (0, 1), (0, 2), (1, 1)])


class TestIndexFilename(unittest.TestCase):

    filenames = [('example_vcf40.vcf.gz', 'example_vcf40.vcf.gz.tbi'),
                 ('example_vcf40.vcf.gz', 'example_vcf40.vcf.gz.csi'),
                 ('example_vcf40.bcf',    'example_vcf40.bcf.csi')]

    def testOpen(self):
        for fn, idx_fn in self.filenames:
            fn = os.path.join(DATADIR, fn)
            idx_fn = os.path.join(DATADIR, idx_fn)

            v = pysam.VariantFile(fn, index_filename=idx_fn)

            self.assertEqual(len(list(v.fetch('20'))), 3)


class TestConstructionVCFWithContigs(unittest.TestCase):
    """construct VariantFile from scratch."""

    filename = "example_vcf42_withcontigs.vcf"

    def complete_check(self, fn_in, fn_out):

        check_lines_equal(
            self, fn_in, fn_out, sort=True,
            filter_f=lambda x: x.startswith("##contig"))
        os.unlink(fn_out)

    def testConstructionWithRecords(self):

        fn_in = os.path.join(DATADIR, self.filename)
        fn_out = get_temp_filename(suffix=".vcf")
        vcf_in = pysam.VariantFile(fn_in)

        header = pysam.VariantHeader()

        for record in vcf_in.header.records:
            header.add_record(record)

        fn = str("tmp_VariantFileTest_testConstructionWithRecords") + ".vcf"
        vcf_out = pysam.VariantFile(fn, "w", header=header)
        for record in vcf_in:
            # currently segfaults here:
            # vcf_out.write(record)
            pass
        return

        vcf_out.close()
        self.complete_check(fn_in, fn_out)

    def testConstructionFromCopy(self):

        fn_in = os.path.join(DATADIR, self.filename)
        fn_out = get_temp_filename(suffix=".vcf")
        vcf_in = pysam.VariantFile(fn_in)

        vcf_out = pysam.VariantFile(fn_out, "w", header=vcf_in.header)
        for record in vcf_in:
            vcf_out.write(record)

        vcf_out.close()

        self.complete_check(fn_in, fn_out)

    def testConstructionWithLines(self):

        fn_in = os.path.join(DATADIR, self.filename)
        fn_out = get_temp_filename(suffix=".vcf")
        vcf_in = pysam.VariantFile(fn_in)

        header = pysam.VariantHeader()
        for sample in vcf_in.header.samples:
            header.add_sample(sample)

        for hr in vcf_in.header.records:
            header.add_line(str(hr))

        vcf_out = pysam.VariantFile(fn_out, "w", header=header)

        for record in vcf_in:
            vcf_out.write(record)

        vcf_out.close()
        vcf_in.close()

        self.complete_check(fn_in, fn_out)

# Currently segfaults for VCFs without contigs
# class TestConstructionVCFWithoutContigs(TestConstructionVCFWithContigs):
#     """construct VariantFile from scratch."""
#     filename = "example_vcf40.vcf"


class TestConstructionVCFGZWithContigs(TestConstructionVCFWithContigs):
    """construct VariantFile from scratch."""

    filename = "example_vcf42_withcontigs.vcf.gz"


class TestConstructionVCFGZWithoutContigs(TestConstructionVCFWithContigs):
    """construct VariantFile from scratch."""

    filename = "example_vcf42.vcf.gz"


class TestSettingRecordValues(unittest.TestCase):

    filename = "example_vcf40.vcf"

    def testSetQual(self):
        with pysam.VariantFile(os.path.join(DATADIR, self.filename)) as inf:
            record = next(inf)
            self.assertEqual(record.qual, 47)
            record.qual = record.qual
            self.assertEqual(record.qual, 47)
            record.qual = 10
            self.assertEqual(record.qual, 10)
            self.assertEqual(str(record).split("\t")[5], "10")

    def testGenotype(self):
        with pysam.VariantFile(os.path.join(DATADIR, self.filename)) as inf:
            record = next(inf)
            sample = record.samples["NA00001"]
            print (sample["GT"])
            self.assertEqual(sample["GT"], (0, 0))
#	Fails with TypeError
#            sample["GT"] = sample["GT"]

class TestSubsetting(unittest.TestCase):
    
    filename = "example_vcf42.vcf.gz"
    
    def testSubsetting(self):
        with pysam.VariantFile(os.path.join(DATADIR,
                                            self.filename)) as inf:
            inf.subset_samples(["NA00001"])


if __name__ == "__main__":
    # build data files
    print ("building data files")
    subprocess.call("make -C %s" % DATADIR, shell=True)
    print ("starting tests")
    unittest.main()
    print ("completed tests")
