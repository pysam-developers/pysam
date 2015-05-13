#!/usr/bin/env python
'''unit testing code for pysam.

Execute in the :file:`tests` directory as it requires the Makefile
and data files located there.
'''

import pysam
import unittest
import os
import shutil
import sys
import collections
import subprocess
import logging
from TestUtils import checkBinaryEqual, checkURL

IS_PYTHON3 = sys.version_info[0] >= 3

SAMTOOLS = "samtools"
WORKDIR = "pysam_test_work"
DATADIR = "pysam_data"


class BasicTestBAMFetch(unittest.TestCase):

    '''basic first test - detailed testing
    if information in file is consistent
    with information in AlignedRead object.'''

    def setUp(self):
        self.samfile = pysam.Samfile(
            os.path.join(DATADIR, "ex3.bam"),
            "rb")
        self.reads = list(self.samfile.fetch())

    def testARqname(self):
        self.assertEqual(
            self.reads[0].qname,
            "read_28833_29006_6945",
            "read name mismatch in read 1: %s != %s" % (
                self.reads[0].qname, "read_28833_29006_6945"))
        self.assertEqual(
            self.reads[1].qname,
            "read_28701_28881_323b",
            "read name mismatch in read 2: %s != %s" % (
                self.reads[1].qname, "read_28701_28881_323b"))

    def testARflag(self):
        self.assertEqual(
            self.reads[0].flag, 99,
            "flag mismatch in read 1: %s != %s" % (
                self.reads[0].flag, 99))
        self.assertEqual(
            self.reads[1].flag, 147,
            "flag mismatch in read 2: %s != %s" % (
                self.reads[1].flag, 147))

    def testARrname(self):
        self.assertEqual(
            self.reads[0].rname, 0,
            "chromosome/target id mismatch in read 1: %s != %s" %
            (self.reads[0].rname, 0))
        self.assertEqual(
            self.reads[1].rname, 1,
            "chromosome/target id mismatch in read 2: %s != %s" %
            (self.reads[1].rname, 1))

    def testARpos(self):
        self.assertEqual(
            self.reads[0].pos, 33 - 1,
            "mapping position mismatch in read 1: %s != %s" %
            (self.reads[0].pos, 33 - 1))
        self.assertEqual(
            self.reads[1].pos, 88 - 1,
            "mapping position mismatch in read 2: %s != %s" %
            (self.reads[1].pos, 88 - 1))

    def testARmapq(self):
        self.assertEqual(
            self.reads[0].mapq, 20,
            "mapping quality mismatch in read 1: %s != %s" %
            (self.reads[0].mapq, 20))
        self.assertEqual(
            self.reads[1].mapq, 30,
            "mapping quality mismatch in read 2: %s != %s" % (
                self.reads[1].mapq, 30))

    def testARcigar(self):
        self.assertEqual(
            self.reads[0].cigar,
            [(0, 10), (2, 1), (0, 25)],
            "read name length mismatch in read 1: %s != %s" %
            (self.reads[0].cigar, [(0, 10), (2, 1), (0, 25)]))
        self.assertEqual(
            self.reads[1].cigar, [(0, 35)],
            "read name length mismatch in read 2: %s != %s" %
            (self.reads[1].cigar, [(0, 35)]))

    def testARcigarstring(self):
        self.assertEqual(self.reads[0].cigarstring, '10M1D25M')
        self.assertEqual(self.reads[1].cigarstring, '35M')

    def testARmrnm(self):
        self.assertEqual(
            self.reads[0].mrnm, 0,
            "mate reference sequence name mismatch in read 1: %s != %s" %
            (self.reads[0].mrnm, 0))
        self.assertEqual(
            self.reads[1].mrnm, 1,
            "mate reference sequence name mismatch in read 2: %s != %s" %
            (self.reads[1].mrnm, 1))
        self.assertEqual(
            self.reads[0].rnext, 0,
            "mate reference sequence name mismatch in read 1: %s != %s" %
            (self.reads[0].rnext, 0))
        self.assertEqual(
            self.reads[1].rnext, 1,
            "mate reference sequence name mismatch in read 2: %s != %s" %
            (self.reads[1].rnext, 1))

    def testARmpos(self):
        self.assertEqual(self.reads[
                         0].mpos, 200 - 1, "mate mapping position mismatch in read 1: %s != %s" % (self.reads[0].mpos, 200 - 1))
        self.assertEqual(self.reads[
                         1].mpos, 500 - 1, "mate mapping position mismatch in read 2: %s != %s" % (self.reads[1].mpos, 500 - 1))
        self.assertEqual(self.reads[
                         0].pnext, 200 - 1, "mate mapping position mismatch in read 1: %s != %s" % (self.reads[0].pnext, 200 - 1))
        self.assertEqual(self.reads[
                         1].pnext, 500 - 1, "mate mapping position mismatch in read 2: %s != %s" % (self.reads[1].pnext, 500 - 1))

    def testARisize(self):
        self.assertEqual(self.reads[0].isize, 167, "insert size mismatch in read 1: %s != %s" % (
            self.reads[0].isize, 167))
        self.assertEqual(self.reads[1].isize, 412, "insert size mismatch in read 2: %s != %s" % (
            self.reads[1].isize, 412))
        self.assertEqual(self.reads[0].tlen, 167, "insert size mismatch in read 1: %s != %s" % (
            self.reads[0].tlen, 167))
        self.assertEqual(self.reads[1].tlen, 412, "insert size mismatch in read 2: %s != %s" % (
            self.reads[1].tlen, 412))

    def testARseq(self):
        self.assertEqual(self.reads[0].seq, "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG", "sequence mismatch in read 1: %s != %s" % (
            self.reads[0].seq, "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"))
        self.assertEqual(self.reads[1].seq, "ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA", "sequence size mismatch in read 2: %s != %s" % (
            self.reads[1].seq, "ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA"))
        self.assertEqual(self.reads[3].seq, "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG", "sequence mismatch in read 4: %s != %s" % (
            self.reads[3].seq, "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"))

    def testARqual(self):
        self.assertEqual(self.reads[0].qual, "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<",
                         "quality string mismatch in read 1: %s != %s" % (self.reads[0].qual, "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<"))
        self.assertEqual(self.reads[1].qual, "<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<", "quality string mismatch in read 2: %s != %s" % (
            self.reads[1].qual, "<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<"))
        self.assertEqual(self.reads[3].qual, "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<",
                         "quality string mismatch in read 3: %s != %s" % (self.reads[3].qual, "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<"))

    def testARquery(self):
        self.assertEqual(self.reads[0].query, "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG", "query mismatch in read 1: %s != %s" % (
            self.reads[0].query, "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"))
        self.assertEqual(self.reads[1].query, "ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA", "query size mismatch in read 2: %s != %s" % (
            self.reads[1].query, "ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA"))
        self.assertEqual(self.reads[3].query, "TAGCTAGCTACCTATATCTTGGTCTT", "query mismatch in read 4: %s != %s" % (
            self.reads[3].query, "TAGCTAGCTACCTATATCTTGGTCTT"))

    def testARqqual(self):
        self.assertEqual(
            self.reads[0].qqual, "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<",
            "qquality string mismatch in read 1: %s != %s" %
            (self.reads[0].qqual, "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<"))
        self.assertEqual(
            self.reads[1].qqual, "<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<",
            "qquality string mismatch in read 2: %s != %s" %
            (self.reads[1].qqual, "<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<"))
        self.assertEqual(
            self.reads[3].qqual, "<<<<<<<<<<<<<<<<<:<9/,&,22",
            "qquality string mismatch in read 3: %s != %s" %
            (self.reads[3].qqual, "<<<<<<<<<<<<<<<<<:<9/,&,22"))

    def testPresentOptionalFields(self):
        self.assertEqual(
            self.reads[0].opt('NM'), 1,
            "optional field mismatch in read 1, NM: %s != %s" %
            (self.reads[0].opt('NM'), 1))
        self.assertEqual(
            self.reads[0].opt('RG'), 'L1',
            "optional field mismatch in read 1, RG: %s != %s" %
            (self.reads[0].opt('RG'), 'L1'))
        self.assertEqual(
            self.reads[1].opt('RG'), 'L2',
            "optional field mismatch in read 2, RG: %s != %s" %
            (self.reads[1].opt('RG'), 'L2'))
        self.assertEqual(
            self.reads[1].opt('MF'), 18,
            "optional field mismatch in read 2, MF: %s != %s" %
            (self.reads[1].opt('MF'), 18))

    def testPairedBools(self):
        self.assertEqual(self.reads[0].is_paired, True,
                         "is paired mismatch in read 1: %s != %s" % (
            self.reads[0].is_paired, True))
        self.assertEqual(self.reads[1].is_paired, True,
                         "is paired mismatch in read 2: %s != %s" % (
            self.reads[1].is_paired, True))
        self.assertEqual(self.reads[0].is_proper_pair, True,
                         "is proper pair mismatch in read 1: %s != %s" % (
            self.reads[0].is_proper_pair, True))
        self.assertEqual(self.reads[1].is_proper_pair, True,
                         "is proper pair mismatch in read 2: %s != %s" % (
            self.reads[1].is_proper_pair, True))

    def testTags(self):
        self.assertEqual(self.reads[0].tags,
                         [('NM', 1), ('RG', 'L1'),
                          ('PG', 'P1'), ('XT', 'U')])
        self.assertEqual(self.reads[1].tags,
                         [('MF', 18), ('RG', 'L2'),
                          ('PG', 'P2'), ('XT', 'R')])

    def testAddTags(self):
        self.assertEqual(sorted(self.reads[0].tags),
                         sorted([('NM', 1), ('RG', 'L1'),
                                 ('PG', 'P1'), ('XT', 'U')]))

        self.reads[0].setTag('X1', 'C')
        self.assertEqual(sorted(self.reads[0].tags),
                         sorted([('X1', 'C'), ('NM', 1), ('RG', 'L1'),
                                 ('PG', 'P1'), ('XT', 'U'), ]))
        self.reads[0].setTag('X2', 5)
        self.assertEqual(sorted(self.reads[0].tags),
                         sorted([('X2', 5), ('X1', 'C'),
                                 ('NM', 1), ('RG', 'L1'),
                                 ('PG', 'P1'), ('XT', 'U'), ]))
        # add with replacement
        self.reads[0].setTag('X2', 10)
        self.assertEqual(sorted(self.reads[0].tags),
                         sorted([('X2', 10), ('X1', 'C'),
                                 ('NM', 1), ('RG', 'L1'),
                                 ('PG', 'P1'), ('XT', 'U'), ]))

        # add without replacement
        self.reads[0].setTag('X2', 5, replace=False)
        self.assertEqual(sorted(self.reads[0].tags),
                         sorted([('X2', 10), ('X1', 'C'),
                                 ('X2', 5),
                                 ('NM', 1), ('RG', 'L1'),
                                 ('PG', 'P1'), ('XT', 'U'), ]))

    def testAddTagsType(self):
        self.reads[0].tags = None
        self.assertEqual(self.reads[0].tags, [])

        self.reads[0].setTag('X1', 5.0)
        self.reads[0].setTag('X2', "5.0")
        self.reads[0].setTag('X3', 5)

        self.assertEqual(sorted(self.reads[0].tags),
                         sorted([('X1', 5.0),
                                 ('X2', "5.0"),
                                 ('X3', 5)]))

        # test setting float for int value
        self.reads[0].setTag('X4', 5, value_type='d')
        self.assertEqual(sorted(self.reads[0].tags),
                         sorted([('X1', 5.0),
                                 ('X2', "5.0"),
                                 ('X3', 5),
                                 ('X4', 5.0)]))

        # test setting int for float value - the
        # value will be rounded.
        self.reads[0].setTag('X5', 5.2, value_type='i')
        self.assertEqual(sorted(self.reads[0].tags),
                         sorted([('X1', 5.0),
                                 ('X2', "5.0"),
                                 ('X3', 5),
                                 ('X4', 5.0),
                                 ('X5', 5)]))

        # test setting invalid type code
        self.assertRaises(ValueError, self.reads[0].setTag, 'X6', 5.2, 'g')

    def testTagsUpdatingFloat(self):
        self.assertEqual(self.reads[0].tags,
                         [('NM', 1), ('RG', 'L1'),
                          ('PG', 'P1'), ('XT', 'U')])
        self.reads[0].tags += [('XC', 5.0)]
        self.assertEqual(self.reads[0].tags,
                         [('NM', 1), ('RG', 'L1'),
                          ('PG', 'P1'), ('XT', 'U'), ('XC', 5.0)])

    def testOpt(self):
        self.assertEqual(self.reads[0].opt("XT"), "U")
        self.assertEqual(self.reads[1].opt("XT"), "R")

    def testMissingOpt(self):
        self.assertRaises(KeyError, self.reads[0].opt, "XP")

    def testEmptyOpt(self):
        self.assertRaises(KeyError, self.reads[2].opt, "XT")

    def tearDown(self):
        self.samfile.close()


class BasicTestBAMFile(BasicTestBAMFetch):

    def setUp(self):
        self.samfile = pysam.Samfile(
            os.path.join(DATADIR, "ex3.sam"),
            "r")
        self.reads = [r for r in self.samfile]


class BasicTestSAMFile(BasicTestBAMFetch):

    def setUp(self):
        self.samfile = pysam.Samfile(
            os.path.join(DATADIR, "ex3.sam"),
            "r")
        self.reads = [r for r in self.samfile]


class BasicTestSAMFetch(BasicTestBAMFetch):

    def setUp(self):
        self.samfile = pysam.Samfile(
            os.path.join(DATADIR, "ex3.sam"),
            "r")
        self.reads = list(self.samfile.fetch())


# needs to be implemented
# class TestAlignedReadFromSamWithoutHeader(TestAlignedReadFromBam):
#
#     def setUp(self):
#         self.samfile=pysam.Samfile( "ex7.sam","r" )
#         self.reads=list(self.samfile.fetch())


class TestIO(unittest.TestCase):

    '''check if reading samfile and writing a samfile are consistent.'''

    def checkEcho(self,
                  input_filename,
                  reference_filename,
                  output_filename,
                  input_mode, output_mode,
                  use_template=True):
        '''iterate through *input_filename* writing to *output_filename* and
        comparing the output to *reference_filename*.

        The files are opened according to the *input_mode* and *output_mode*.

        If *use_template* is set, the header is copied from infile
        using the template mechanism, otherwise target names and
        lengths are passed explicitely.

        '''

        infile = pysam.Samfile(os.path.join(DATADIR, input_filename),
                               input_mode)
        if use_template:
            outfile = pysam.Samfile(output_filename,
                                    output_mode,
                                    template=infile)
        else:
            outfile = pysam.Samfile(output_filename,
                                    output_mode,
                                    referencenames=infile.references,
                                    referencelengths=infile.lengths,
                                    add_sq_text=False)

        iter = infile.fetch()

        for x in iter:
            outfile.write(x)
        infile.close()
        outfile.close()

        self.assertTrue(
            checkBinaryEqual(os.path.join(DATADIR, reference_filename),
                             output_filename),
            "files %s and %s are not the same" % (reference_filename,
                                                  output_filename))

    def testReadWriteBam(self):

        input_filename = "ex1.bam"
        output_filename = "pysam_ex1.bam"
        reference_filename = "ex1.bam"

        self.checkEcho(input_filename, reference_filename, output_filename,
                       "rb", "wb", use_template=True)

    # Disabled - should work, files are not binary equal, but are
    # non-binary equal:
    # diff <(samtools view pysam_ex1.bam) <(samtools view pysam_data/ex1.bam)
    # def testReadWriteBamWithTargetNames(self):
    #     input_filename = "ex1.bam"
    #     output_filename = "pysam_ex1.bam"
    #     reference_filename = "ex1.bam"

    #     self.checkEcho(input_filename, reference_filename, output_filename,
    #                    "rb", "wb", use_template=False)

    def testReadWriteSamWithHeader(self):

        input_filename = "ex2.sam"
        output_filename = "pysam_ex2.sam"
        reference_filename = "ex2.sam"

        self.checkEcho(input_filename,
                       reference_filename,
                       output_filename,
                       "r", "wh")

    # Release 0.8.0
    # no samfiles without header
    def testReadWriteSamWithoutHeader(self):

        input_filename = "ex2.sam"
        output_filename = "pysam_ex2.sam"
        reference_filename = "ex1.sam"

        self.checkEcho(input_filename,
                       reference_filename,
                       output_filename,
                       "r", "w")

    def testReadSamWithoutTargetNames(self):
        '''see issue 104.'''
        input_filename = os.path.join(DATADIR,
                                      "example_unmapped_reads_no_sq.sam")

        # raise exception in default mode
        self.assertRaises(ValueError, pysam.Samfile, input_filename, "r")

        # raise exception if no SQ files
        self.assertRaises(ValueError, pysam.Samfile,
                          input_filename, "r",
                          check_header=True)

        infile = pysam.Samfile(
            input_filename,
            check_header=False,
            check_sq=False)

        # TODO
        # result = list(infile.fetch(until_eof=True))
        # self.assertEqual(2, len(result))

    def testReadBamWithoutTargetNames(self):
        '''see issue 104.'''
        input_filename = os.path.join(
            DATADIR, "example_unmapped_reads_no_sq.bam")

        # raise exception in default mode
        self.assertRaises(ValueError, pysam.Samfile, input_filename, "r")

        # raise exception if no SQ files
        self.assertRaises(ValueError, pysam.Samfile, input_filename, "r",
                          check_header=True)

        infile = pysam.Samfile(
            input_filename, check_header=False, check_sq=False)
        result = list(infile.fetch(until_eof=True))

    # TODO
    def testReadSamWithoutHeader(self):
        input_filename = os.path.join(DATADIR, "ex1.sam")

        # reading from a samfile without header is not
        # implemented
        self.assertRaises(ValueError,
                          pysam.Samfile,
                          input_filename,
                          "r")

        # TODO
        # without check_header header is no read
        # leading to segfault
        # self.assertRaises(ValueError,
        #                   pysam.Samfile,
        #                   input_filename,
        #                   "r",
        #                   check_header=False)

    # TODO
    # def testReadUnformattedFile(self):
    #     '''test reading from a file that is not bam/sam formatted'''
    #     input_filename = os.path.join(DATADIR, 'Makefile')

    #     # bam - file raise error
    #     self.assertRaises(ValueError,
    #                       pysam.Samfile,
    #                       input_filename,
    #                       "rb")

    #     # sam - file error, but can't fetch
    #     self.assertRaises(ValueError,
    #                       pysam.Samfile,
    #                       input_filename,
    #                       "r")

    #     self.assertRaises(ValueError,
    #                       pysam.Samfile,
    #                       input_filename,
    #                       "r",
    #                       check_header=False)

    def testBAMWithoutAlignedReads(self):
        '''see issue 117'''
        input_filename = os.path.join(DATADIR, "test_unaligned.bam")
        samfile = pysam.Samfile(input_filename, "rb", check_sq=False)
        samfile.fetch(until_eof=True)

    def testBAMWithShortBAI(self):
        '''see issue 116'''
        input_filename = os.path.join(DATADIR, "example_bai.bam")
        samfile = pysam.Samfile(input_filename, "rb", check_sq=False)
        samfile.fetch('chr2')

    def testFetchFromClosedFile(self):

        samfile = pysam.Samfile(os.path.join(DATADIR, "ex1.bam"),
                                "rb")
        samfile.close()
        self.assertRaises(ValueError, samfile.fetch, 'chr1', 100, 120)

    def testClosedFile(self):
        '''test that access to a closed samfile raises ValueError.'''

        samfile = pysam.Samfile(os.path.join(DATADIR, "ex1.bam"),
                                "rb")
        samfile.close()
        self.assertRaises(ValueError, samfile.fetch, 'chr1', 100, 120)
        self.assertRaises(ValueError, samfile.pileup, 'chr1', 100, 120)
        self.assertRaises(ValueError, samfile.getrname, 0)
        # TODO
        self.assertRaises(ValueError, samfile.tell)
        self.assertRaises(ValueError, samfile.seek, 0)
        self.assertRaises(ValueError, getattr, samfile, "nreferences")
        self.assertRaises(ValueError, getattr, samfile, "references")
        self.assertRaises(ValueError, getattr, samfile, "lengths")
        self.assertRaises(ValueError, getattr, samfile, "text")
        self.assertRaises(ValueError, getattr, samfile, "header")

        # write on closed file
        self.assertEqual(0, samfile.write(None))

    def testAutoDetection(self):
        '''test if autodetection works.'''

        # TODO
        # samfile = pysam.Samfile(os.path.join(DATADIR, "ex3.sam"))
        # self.assertRaises(ValueError, samfile.fetch, 'chr1')
        # samfile.close()

        samfile = pysam.Samfile(os.path.join(DATADIR, "ex3.bam"))
        samfile.fetch('chr1')
        samfile.close()

    # TOOD
    # def testReadingFromSamFileWithoutHeader(self):
    #     '''read from samfile without header.
    #     '''
    #     samfile = pysam.Samfile(os.path.join(DATADIR, "ex7.sam"),
    #                             check_header=False,
    #                             check_sq=False)
    #     self.assertRaises(NotImplementedError, samfile.__iter__)

    def testReadingFromFileWithoutIndex(self):
        '''read from bam file without index.'''

        shutil.copyfile(os.path.join(DATADIR, "ex2.bam"), 'tmp_ex2.bam')
        samfile = pysam.Samfile('tmp_ex2.bam',
                                "rb")
        self.assertRaises(ValueError, samfile.fetch)
        self.assertEqual(len(list(samfile.fetch(until_eof=True))),
                         3270)
        os.unlink('tmp_ex2.bam')

    # def testReadingUniversalFileMode(self):
    #     '''read from samfile without header.
    #     '''

    #     input_filename = "ex2.sam"
    #     output_filename = "pysam_ex2.sam"
    #     reference_filename = "ex1.sam"

    #     self.checkEcho(input_filename,
    #                    reference_filename,
    #                    output_filename,
    #                    "rU", "w")

    def testHead(self):
        '''test IteratorRowHead'''
        samfile = pysam.Samfile(os.path.join(DATADIR, "ex1.bam"),
                                "rb")
        l10 = list(samfile.head(10))
        l100 = list(samfile.head(100))
        self.assertEqual(len(l10), 10)
        self.assertEqual(len(l100), 100)
        self.assertEqual(list(map(str, l10)),
                         list(map(str, l100[:10])))


class TestFloatTagBug(unittest.TestCase):

    '''see issue 71'''

    def testFloatTagBug(self):
        '''a float tag before another exposed a parsing bug in bam_aux_get.

        Fixed in 0.1.19
        '''
        samfile = pysam.Samfile(os.path.join(DATADIR, "tag_bug.bam"))
        read = next(samfile.fetch(until_eof=True))
        self.assertTrue(('XC', 1) in read.tags)
        self.assertEqual(read.opt('XC'), 1)


class TestLargeFieldBug(unittest.TestCase):

    '''see issue 100'''

    def testLargeFileBug(self):
        '''when creating a read with a large entry in the tag field
        causes an errror:
            NotImplementedError: tags field too large
        '''
        samfile = pysam.Samfile(os.path.join(DATADIR, "issue100.bam"))
        read = next(samfile.fetch(until_eof=True))
        new_read = pysam.AlignedRead()
        new_read.tags = read.tags
        self.assertEqual(new_read.tags, read.tags)


class TestTagParsing(unittest.TestCase):

    '''tests checking the accuracy of tag setting and retrieval.'''

    def makeRead(self):
        a = pysam.AlignedRead()
        a.qname = "read_12345"
        a.tid = 0
        a.seq = "ACGT" * 3
        a.flag = 0
        a.rname = 0
        a.pos = 1
        a.mapq = 20
        a.cigar = ((0, 10), (2, 1), (0, 25))
        a.mrnm = 0
        a.mpos = 200
        a.isize = 0
        a.qual = "1234" * 3
        # todo: create tags
        return a

    def testNegativeIntegers(self):
        x = -2
        aligned_read = self.makeRead()
        aligned_read.tags = [("XD", int(x))]
        # print (aligned_read.tags)

    def testNegativeIntegers2(self):
        x = -2
        r = self.makeRead()
        r.tags = [("XD", int(x))]
        outfile = pysam.Samfile("test.bam",
                                "wb",
                                referencenames=("chr1",),
                                referencelengths = (1000,))
        outfile.write(r)
        outfile.close()

    def testCigarString(self):
        r = self.makeRead()
        self.assertEqual(r.cigarstring, "10M1D25M")
        r.cigarstring = "20M10D20M"
        self.assertEqual(r.cigar, [(0, 20), (2, 10), (0, 20)])
        # unsetting cigar string
        r.cigarstring = None
        self.assertEqual(r.cigarstring, None)

    def testCigar(self):
        r = self.makeRead()
        self.assertEqual(r.cigar, [(0, 10), (2, 1), (0, 25)])
        # unsetting cigar string
        r.cigar = None
        self.assertEqual(r.cigar, [])

    def testLongTags(self):
        '''see issue 115'''

        r = self.makeRead()
        rg = 'HS2000-899_199.L3'
        tags = [('XC', 85), ('XT', 'M'), ('NM', 5),
                ('SM', 29), ('AM', 29), ('XM', 1),
                ('XO', 1), ('XG', 4), ('MD', '37^ACCC29T18'),
                ('XA', '5,+11707,36M1I48M,2;21,-48119779,46M1I38M,2;hs37d5,-10060835,40M1D45M,3;5,+11508,36M1I48M,3;hs37d5,+6743812,36M1I48M,3;19,-59118894,46M1I38M,3;4,-191044002,6M1I78M,3;')]

        r.tags = tags
        r.tags += [("RG", rg)] * 100
        tags += [("RG", rg)] * 100

        self.assertEqual(tags, r.tags)


class TestClipping(unittest.TestCase):

    def testClipping(self):

        self.samfile = pysam.Samfile(os.path.join(DATADIR, "softclip.bam"),
                                     "rb")
        for read in self.samfile:

            if read.qname == "r001":
                self.assertEqual(read.seq, 'AAAAGATAAGGATA')
                self.assertEqual(read.query, 'AGATAAGGATA')
                self.assertEqual(read.qual, None)
                self.assertEqual(read.qqual, None)

            elif read.qname == "r002":

                self.assertEqual(read.seq, 'GCCTAAGCTAA')
                self.assertEqual(read.query, 'AGCTAA')
                self.assertEqual(read.qual, '01234567890')
                self.assertEqual(read.qqual, '567890')

            elif read.qname == "r003":

                self.assertEqual(read.seq, 'GCCTAAGCTAA')
                self.assertEqual(read.query, 'GCCTAA')
                self.assertEqual(read.qual, '01234567890')
                self.assertEqual(read.qqual, '012345')

            elif read.qname == "r004":

                self.assertEqual(read.seq, 'TAGGC')
                self.assertEqual(read.query, 'TAGGC')
                self.assertEqual(read.qual, '01234')
                self.assertEqual(read.qqual, '01234')


class TestIteratorRow(unittest.TestCase):

    def setUp(self):
        self.samfile = pysam.Samfile(os.path.join(DATADIR, "ex1.bam"),
                                     "rb")

    def checkRange(self, rnge):
        '''compare results from iterator with those from samtools.'''
        ps = list(self.samfile.fetch(region=rnge))
        sa = list(pysam.view(os.path.join(DATADIR, "ex1.bam"),
                             rnge,
                             raw=True))
        self.assertEqual(len(ps), len(
            sa), "unequal number of results for range %s: %i != %i" % (rnge, len(ps), len(sa)))
        # check if the same reads are returned and in the same order
        for line, (a, b) in enumerate(list(zip(ps, sa))):
            d = b.split("\t")
            self.assertEqual(
                a.qname, d[0], "line %i: read id mismatch: %s != %s" % (line, a.rname, d[0]))
            self.assertEqual(a.pos, int(d[3]) - 1, "line %i: read position mismatch: %s != %s, \n%s\n%s\n" %
                             (line, a.pos, int(d[3]) - 1,
                              str(a), str(d)))
            qual = d[10]
            self.assertEqual(a.qual, qual, "line %i: quality mismatch: %s != %s, \n%s\n%s\n" %
                             (line, a.qual, qual,
                              str(a), str(d)))

    def testIteratePerContig(self):
        '''check random access per contig'''
        for contig in self.samfile.references:
            self.checkRange(contig)

    def testIterateRanges(self):
        '''check random access per range'''
        for contig, length in zip(self.samfile.references, self.samfile.lengths):
            for start in range(1, length, 90):
                # this includes empty ranges
                self.checkRange("%s:%i-%i" % (contig, start, start + 90))

    def tearDown(self):
        self.samfile.close()


class TestIteratorRowAll(unittest.TestCase):

    def setUp(self):
        self.samfile = pysam.Samfile(os.path.join(DATADIR, "ex1.bam"),
                                     "rb")

    def testIterate(self):
        '''compare results from iterator with those from samtools.'''
        ps = list(self.samfile.fetch())
        sa = list(pysam.view(os.path.join(DATADIR, "ex1.bam"),
                             raw=True))
        self.assertEqual(
            len(ps), len(sa), "unequal number of results: %i != %i" % (len(ps), len(sa)))
        # check if the same reads are returned
        for line, pair in enumerate(list(zip(ps, sa))):
            data = pair[1].split("\t")
            self.assertEqual(pair[0].qname, data[
                             0], "read id mismatch in line %i: %s != %s" % (line, pair[0].rname, data[0]))

    def tearDown(self):
        self.samfile.close()


class TestIteratorColumn(unittest.TestCase):

    '''test iterator column against contents of ex4.bam.'''

    # note that samfile contains 1-based coordinates
    # 1D means deletion with respect to reference sequence
    #
    mCoverages = {'chr1': [0] * 20 + [1] * 36 + [0] * (100 - 20 - 35),
                  'chr2': [0] * 20 + [1] * 35 + [0] * (100 - 20 - 35),
                  }

    def setUp(self):
        self.samfile = pysam.Samfile(os.path.join(DATADIR, "ex4.bam"),
                                     "rb")

    def checkRange(self, contig, start=None, end=None, truncate=False):
        '''compare results from iterator with those from samtools.'''
        # check if the same reads are returned and in the same order
        for column in self.samfile.pileup(contig, start, end,
                                          truncate=truncate):
            if truncate:
                self.assertGreaterEqual(column.pos, start)
                self.assertLess(column.pos, end)
            thiscov = len(column.pileups)
            refcov = self.mCoverages[
                self.samfile.getrname(column.tid)][column.pos]
            self.assertEqual(
                thiscov, refcov, "wrong coverage at pos %s:%i %i should be %i" % (
                    self.samfile.getrname(column.tid), column.pos, thiscov, refcov))

    def testIterateAll(self):
        '''check random access per contig'''
        self.checkRange(None)

    def testIteratePerContig(self):
        '''check random access per contig'''
        for contig in self.samfile.references:
            self.checkRange(contig)

    def testIterateRanges(self):
        '''check random access per range'''
        for contig, length in zip(self.samfile.references, self.samfile.lengths):
            for start in range(1, length, 90):
                # this includes empty ranges
                self.checkRange(contig, start, start + 90)

    def testInverse(self):
        '''test the inverse, is point-wise pileup accurate.'''
        for contig, refseq in list(self.mCoverages.items()):
            refcolumns = sum(refseq)
            for pos, refcov in enumerate(refseq):
                columns = list(self.samfile.pileup(contig, pos, pos + 1))
                if refcov == 0:
                    # if no read, no coverage
                    self.assertEqual(
                        len(columns),
                        refcov,
                        "wrong number of pileup columns returned for position %s:%i, %i should be %i" % (
                            contig, pos,
                            len(columns), refcov))
                elif refcov == 1:
                    # one read, all columns of the read are returned
                    self.assertEqual(
                        len(columns),
                        refcolumns,
                        "pileup incomplete at position %i: got %i, expected %i " %
                        (pos, len(columns), refcolumns))

    def testIterateTruncate(self):
        '''check random access per range'''
        for contig, length in zip(self.samfile.references, self.samfile.lengths):
            for start in range(1, length, 90):
                # this includes empty ranges
                self.checkRange(contig, start, start + 90, truncate=True)

    def tearDown(self):
        self.samfile.close()


class TestIteratorColumn2(unittest.TestCase):

    '''test iterator column against contents of ex1.bam.'''

    def setUp(self):
        self.samfile = pysam.Samfile(os.path.join(DATADIR, "ex1.bam"),
                                     "rb")

    def testStart(self):
        # print self.samfile.fetch().next().pos
        # print self.samfile.pileup().next().pos
        pass

    def testTruncate(self):
        '''see issue 107.'''
        # note that ranges in regions start from 1
        p = self.samfile.pileup(region='chr1:170:172', truncate=True)
        columns = [x.pos for x in p]
        self.assertEqual(len(columns), 3)
        self.assertEqual(columns, [169, 170, 171])

        p = self.samfile.pileup('chr1', 169, 172, truncate=True)
        columns = [x.pos for x in p]

        self.assertEqual(len(columns), 3)
        self.assertEqual(columns, [169, 170, 171])

    def testAccessOnClosedIterator(self):
        '''see issue 131

        Accessing pileup data after iterator has closed.
        '''
        pcolumn = self.samfile.pileup('chr1', 170, 180).__next__()
        self.assertRaises(ValueError, getattr, pcolumn, "pileups")


class TestHeaderSam(unittest.TestCase):

    header = {'SQ': [{'LN': 1575, 'SN': 'chr1'},
                     {'LN': 1584, 'SN': 'chr2'}],
              'RG': [{'LB': 'SC_1', 'ID': 'L1', 'SM': 'NA12891', 'PU': 'SC_1_10', "CN": "name:with:colon"},
                     {'LB': 'SC_2', 'ID': 'L2', 'SM': 'NA12891', 'PU': 'SC_2_12', "CN": "name:with:colon"}],
              'PG': [{'ID': 'P1', 'VN': '1.0'}, {'ID': 'P2', 'VN': '1.1'}],
              'HD': {'VN': '1.0'},
              'CO': ['this is a comment', 'this is another comment'],
              }

    def compareHeaders(self, a, b):
        '''compare two headers a and b.'''
        for ak, av in a.items():
            self.assertTrue(ak in b, "key '%s' not in '%s' " % (ak, b))
            self.assertEqual(av, b[ak])

    def setUp(self):
        self.samfile = pysam.Samfile(os.path.join(DATADIR, "ex3.sam"),
                                     "r")

    def testHeaders(self):
        self.compareHeaders(self.header, self.samfile.header)
        self.compareHeaders(self.samfile.header, self.header)

    def testNameMapping(self):
        for x, y in enumerate(("chr1", "chr2")):
            tid = self.samfile.gettid(y)
            ref = self.samfile.getrname(x)
            self.assertEqual(tid, x)
            self.assertEqual(ref, y)

        self.assertEqual(self.samfile.gettid("chr?"), -1)
        self.assertRaises(ValueError, self.samfile.getrname, 2)

    def tearDown(self):
        self.samfile.close()


class TestHeaderBam(TestHeaderSam):

    def setUp(self):
        self.samfile = pysam.Samfile(os.path.join(DATADIR, "ex3.bam"),
                                     "rb")


class TestHeaderFromRefs(unittest.TestCase):

    '''see issue 144

    reference names need to be converted to string for python 3
    '''

    # def testHeader( self ):
    #     refs = ['chr1', 'chr2']
    #     tmpfile = "tmp_%i" % id(self)
    #     s = pysam.Samfile(tmpfile, 'wb',
    #                       referencenames=refs,
    #                       referencelengths=[100]*len(refs))
    #     s.close()

    #     self.assertTrue( checkBinaryEqual( 'issue144.bam', tmpfile ),
    #                      'bam files differ')
    #     os.unlink( tmpfile )


class TestHeader1000Genomes(unittest.TestCase):

    '''see issue 110'''
    # bamfile = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/phase2b_alignment/data/NA07048/exome_alignment/NA07048.unmapped.ILLUMINA.bwa.CEU.exome.20120522_p2b.bam"
    bamfile = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/phase3_EX_or_LC_only_alignment/data/HG00104/alignment/HG00104.chrom11.ILLUMINA.bwa.GBR.low_coverage.20130415.bam"

    def testRead(self):

        if not checkURL(self.bamfile):
            return

        f = pysam.Samfile(self.bamfile, "rb")
        data = f.header.copy()
        self.assertTrue(data)


class TestUnmappedReads(unittest.TestCase):

    # TODO
    # def testSAM(self):
    #     samfile = pysam.Samfile(os.path.join(DATADIR, "ex5.sam"),
    #                             "r")
    #     self.assertEqual(len(list(samfile.fetch(until_eof=True))), 2)
    #     samfile.close()

    def testBAM(self):
        samfile = pysam.Samfile(os.path.join(DATADIR, "ex5.bam"),
                                "rb")
        self.assertEqual(len(list(samfile.fetch(until_eof=True))), 2)
        samfile.close()


class TestPileupObjects(unittest.TestCase):

    def setUp(self):
        self.samfile = pysam.Samfile(os.path.join(DATADIR, "ex1.bam"),
                                     "rb")

    def testPileupColumn(self):
        for pcolumn1 in self.samfile.pileup(region="chr1:105"):
            if pcolumn1.pos == 104:
                self.assertEqual(
                    pcolumn1.tid, 0, "chromosome/target id mismatch in position 1: %s != %s" % (pcolumn1.tid, 0))
                self.assertEqual(
                    pcolumn1.pos, 105 - 1, "position mismatch in position 1: %s != %s" % (pcolumn1.pos, 105 - 1))
                self.assertEqual(
                    pcolumn1.n, 2, "# reads mismatch in position 1: %s != %s" % (pcolumn1.n, 2))
        for pcolumn2 in self.samfile.pileup(region="chr2:1480"):
            if pcolumn2.pos == 1479:
                self.assertEqual(
                    pcolumn2.tid, 1, "chromosome/target id mismatch in position 1: %s != %s" % (pcolumn2.tid, 1))
                self.assertEqual(
                    pcolumn2.pos, 1480 - 1, "position mismatch in position 1: %s != %s" % (pcolumn2.pos, 1480 - 1))
                self.assertEqual(
                    pcolumn2.n, 12, "# reads mismatch in position 1: %s != %s" % (pcolumn2.n, 12))

    def testPileupRead(self):
        for pcolumn1 in self.samfile.pileup(region="chr1:105"):
            if pcolumn1.pos == 104:
                self.assertEqual(
                    len(pcolumn1.pileups), 2,
                    "# reads aligned to column mismatch in position 1"
                    ": %s != %s" %
                    (len(pcolumn1.pileups), 2))


# self.assertEqual( pcolumn1.pileups[0]  # need to test additional
# properties here

    def tearDown(self):
        self.samfile.close()

    def testIteratorOutOfScope(self):
        '''test if exception is raised if pileup col is accessed after
        iterator is exhausted.'''

        for pileupcol in self.samfile.pileup():
            pass

        self.assertRaises(ValueError, getattr, pileupcol, "pileups")


class TestContextManager(unittest.TestCase):

    def testManager(self):
        with pysam.Samfile(os.path.join(DATADIR, 'ex1.bam'),
                           'rb') as samfile:
            samfile.fetch()
        self.assertEqual(samfile._isOpen(), False)


class TestExceptions(unittest.TestCase):

    def setUp(self):
        self.samfile = pysam.Samfile(os.path.join(DATADIR, "ex1.bam"),
                                     "rb")

    def testMissingFile(self):

        self.assertRaises(IOError, pysam.Samfile, "exdoesntexist.bam", "rb")
        self.assertRaises(IOError, pysam.Samfile, "exdoesntexist.sam", "r")
        self.assertRaises(IOError, pysam.Samfile, "exdoesntexist.bam", "r")
        self.assertRaises(IOError, pysam.Samfile, "exdoesntexist.sam", "rb")

    def testBadContig(self):
        self.assertRaises(ValueError, self.samfile.fetch, "chr88")

    def testMeaninglessCrap(self):
        self.assertRaises(ValueError, self.samfile.fetch, "skljf")

    def testBackwardsOrderNewFormat(self):
        self.assertRaises(ValueError, self.samfile.fetch, 'chr1', 100, 10)

    def testBackwardsOrderOldFormat(self):
        self.assertRaises(ValueError, self.samfile.fetch, region="chr1:100-10")

    def testOutOfRangeNegativeNewFormat(self):
        self.assertRaises(ValueError, self.samfile.fetch, "chr1", 5, -10)
        self.assertRaises(ValueError, self.samfile.fetch, "chr1", 5, 0)
        self.assertRaises(ValueError, self.samfile.fetch, "chr1", -5, -10)

        self.assertRaises(ValueError, self.samfile.count, "chr1", 5, -10)
        self.assertRaises(ValueError, self.samfile.count, "chr1", 5, 0)
        self.assertRaises(ValueError, self.samfile.count, "chr1", -5, -10)

    def testOutOfRangeNegativeOldFormat(self):
        self.assertRaises(ValueError, self.samfile.fetch, region="chr1:-5-10")
        self.assertRaises(ValueError, self.samfile.fetch, region="chr1:-5-0")
        self.assertRaises(ValueError, self.samfile.fetch, region="chr1:-5--10")

        self.assertRaises(ValueError, self.samfile.count, region="chr1:-5-10")
        self.assertRaises(ValueError, self.samfile.count, region="chr1:-5-0")
        self.assertRaises(ValueError, self.samfile.count, region="chr1:-5--10")

    def testOutOfRangNewFormat(self):
        self.assertRaises(
            ValueError, self.samfile.fetch, "chr1", 9999999999, 99999999999)
        self.assertRaises(
            ValueError, self.samfile.count, "chr1", 9999999999, 99999999999)

    def testOutOfRangeLargeNewFormat(self):
        self.assertRaises(ValueError, self.samfile.fetch, "chr1",
                          9999999999999999999999999999999, 9999999999999999999999999999999999999999)
        self.assertRaises(ValueError, self.samfile.count, "chr1",
                          9999999999999999999999999999999, 9999999999999999999999999999999999999999)

    def testOutOfRangeLargeOldFormat(self):
        self.assertRaises(
            ValueError, self.samfile.fetch, "chr1:99999999999999999-999999999999999999")
        self.assertRaises(
            ValueError, self.samfile.count, "chr1:99999999999999999-999999999999999999")

    def testZeroToZero(self):
        '''see issue 44'''
        self.assertEqual(len(list(self.samfile.fetch('chr1', 0, 0))), 0)

    def tearDown(self):
        self.samfile.close()


class TestWrongFormat(unittest.TestCase):

    '''test cases for opening files not in bam/sam format.'''

    def testOpenSamAsBam(self):
        self.assertRaises(ValueError,
                          pysam.Samfile,
                          os.path.join(DATADIR, 'ex1.sam'),
                          'rb')

    def testOpenBamAsSam(self):
        # test fails, needs to be implemented.
        # sam.fetch() fails on reading, not on opening
        # self.assertRaises( ValueError, pysam.Samfile, 'ex1.bam', 'r' )
        pass

    def testOpenFastaAsSam(self):
        # test fails, needs to be implemented.
        # sam.fetch() fails on reading, not on opening
        # self.assertRaises( ValueError, pysam.Samfile, 'ex1.fa', 'r' )
        pass

    def testOpenFastaAsBam(self):
        self.assertRaises(ValueError,
                          pysam.Samfile,
                          os.path.join(DATADIR, 'ex1.fa'),
                          'rb')


class ReadTest(unittest.TestCase):

    def checkFieldEqual(self, read1, read2, exclude=[]):
        '''check if two reads are equal by comparing each field.'''

        # add the . for refactoring purposes.
        for x in (".qname", ".seq", ".flag",
                  ".rname", ".pos", ".mapq", ".cigar",
                  ".mrnm", ".mpos", ".isize",
                  ".qual",
                  ".bin",
                  ".is_paired", ".is_proper_pair",
                  ".is_unmapped", ".mate_is_unmapped",
                  ".is_reverse", ".mate_is_reverse",
                  ".is_read1", ".is_read2",
                  ".is_secondary", ".is_qcfail",
                  ".is_duplicate"):
            n = x[1:]
            if n in exclude:
                continue
            self.assertEqual(getattr(read1, n), getattr(read2, n),
                             "attribute mismatch for %s: %s != %s" %
                             (n, getattr(read1, n), getattr(read2, n)))


class TestAlignedRead(ReadTest):

    '''tests to check if aligned read can be constructed
    and manipulated.
    '''

    def testEmpty(self):
        a = pysam.AlignedRead()
        self.assertEqual(a.qname, None)
        self.assertEqual(a.seq, None)
        self.assertEqual(a.qual, None)
        self.assertEqual(a.flag, 0)
        self.assertEqual(a.rname, 0)
        self.assertEqual(a.mapq, 0)
        self.assertEqual(a.cigar, [])
        self.assertEqual(a.tags, [])
        self.assertEqual(a.mrnm, 0)
        self.assertEqual(a.mpos, 0)
        self.assertEqual(a.isize, 0)

    def testStrOfEmptyRead(self):
        a = pysam.AlignedRead()
        s = str(a)
        self.assertEqual(
            "None\t0\t0\t0\t0\tNone\t0\t0\t0\tNone\tNone\t[]",
            s)

    def buildRead(self):
        '''build an example read.'''

        a = pysam.AlignedRead()
        a.qname = "read_12345"
        a.seq = "ACGT" * 10
        a.flag = 0
        a.rname = 0
        a.pos = 20
        a.mapq = 20
        a.cigar = ((0, 10), (2, 1), (0, 9), (1, 1), (0, 20))
        a.mrnm = 0
        a.mpos = 200
        a.isize = 167
        a.qual = "1234" * 10
        # todo: create tags
        return a

    def testUpdate(self):
        '''check if updating fields affects other variable length data
        '''
        a = self.buildRead()
        b = self.buildRead()

        # check qname
        b.qname = "read_123"
        self.checkFieldEqual(a, b, "qname")
        b.qname = "read_12345678"
        self.checkFieldEqual(a, b, "qname")
        b.qname = "read_12345"
        self.checkFieldEqual(a, b)

        # check cigar
        b.cigar = ((0, 10), )
        self.checkFieldEqual(a, b, "cigar")
        b.cigar = ((0, 10), (2, 1), (0, 10))
        self.checkFieldEqual(a, b, "cigar")
        b.cigar = ((0, 10), (2, 1), (0, 9), (1, 1), (0, 20))
        self.checkFieldEqual(a, b)

        # check seq
        b.seq = "ACGT"
        self.checkFieldEqual(a, b, ("seq", "qual"))
        b.seq = "ACGT" * 3
        self.checkFieldEqual(a, b, ("seq", "qual"))
        b.seq = "ACGT" * 10
        self.checkFieldEqual(a, b, ("qual",))

        # reset qual
        b = self.buildRead()

        # check flags:
        for x in (
                "is_paired", "is_proper_pair",
                "is_unmapped", "mate_is_unmapped",
                "is_reverse", "mate_is_reverse",
                "is_read1", "is_read2",
                "is_secondary", "is_qcfail",
                "is_duplicate"):
            setattr(b, x, True)
            self.assertEqual(getattr(b, x), True)
            self.checkFieldEqual(a, b, ("flag", x,))
            setattr(b, x, False)
            self.assertEqual(getattr(b, x), False)
            self.checkFieldEqual(a, b)

    def testUpdate2(self):
        '''issue 135: inplace update of sequence and quality score.

        This does not work as setting the sequence will erase
        the quality scores.
        '''
        a = self.buildRead()
        a.seq = a.seq[5:10]
        self.assertEqual(a.qual, None)

        a = self.buildRead()
        s = a.qual
        a.seq = a.seq[5:10]
        a.qual = s[5:10]

        self.assertEqual(a.qual, s[5:10])

    def testLargeRead(self):
        '''build an example read.'''

        a = pysam.AlignedRead()
        a.qname = "read_12345"
        a.seq = "ACGT" * 200
        a.flag = 0
        a.rname = 0
        a.pos = 20
        a.mapq = 20
        a.cigar = ((0, 4 * 200), )
        a.mrnm = 0
        a.mpos = 200
        a.isize = 167
        a.qual = "1234" * 200

        return a

    def testTagParsing(self):
        '''test for tag parsing

        see http://groups.google.com/group/pysam-user-group/browse_thread/thread/67ca204059ea465a
        '''
        samfile = pysam.Samfile(os.path.join(DATADIR, "ex8.bam"),
                                "rb")

        for entry in samfile:
            before = entry.tags
            entry.tags = entry.tags
            after = entry.tags
            self.assertEqual(after, before)

    def testUpdateTlen(self):
        '''check if updating tlen works'''
        a = self.buildRead()
        oldlen = a.tlen
        oldlen *= 2
        a.tlen = oldlen
        self.assertEqual(a.tlen, oldlen)

    def testPositions(self):
        a = self.buildRead()
        self.assertEqual(a.positions,
                         [20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                          31, 32, 33, 34, 35, 36, 37, 38, 39,
                          40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
                          50, 51, 52, 53, 54, 55, 56, 57, 58, 59])

        self.assertEqual(a.aligned_pairs,
                         [(0, 20), (1, 21), (2, 22), (3, 23), (4, 24),
                          (5, 25), (6, 26), (7, 27), (8, 28), (9, 29),
                          (None, 30),
                          (10, 31), (11, 32), (12, 33), (13, 34), (14, 35),
                          (15, 36), (16, 37), (17, 38), (18, 39), (19, None),
                          (20, 40), (21, 41), (22, 42), (23, 43), (24, 44),
                          (25, 45), (26, 46), (27, 47), (28, 48), (29, 49),
                          (30, 50), (31, 51), (32, 52), (33, 53), (34, 54),
                          (35, 55), (36, 56), (37, 57), (38, 58), (39, 59)])

        self.assertEqual(
            a.positions,
            [x[1] for x in a.aligned_pairs
             if x[0] is not None and x[1] is not None])
        # alen is the length of the aligned read in genome
        self.assertEqual(a.alen, a.aligned_pairs[-1][0] + 1)
        # aend points to one beyond last aligned base in ref
        self.assertEqual(a.positions[-1], a.aend - 1)

    def testBlocks(self):
        a = self.buildRead()
        self.assertEqual(a.blocks,
                         [(20, 30), (31, 40), (40, 60)])

    # Disabled as not backwards compatible
    # def testFancyStr(self):
    #     a = self.buildRead()
    #     output = a.fancy_str()
    #     self.assertEqual(len(output), 9)


class TestDeNovoConstruction(ReadTest):

    '''check BAM/SAM file construction using ex6.sam

    (note these are +1 coordinates):

    read_28833_29006_6945	99	chr1	33	20	10M1D25M	=	200	167	AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG	<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<	NM:i:1	RG:Z:L1
    read_28701_28881_323b	147	chr2	88	30	35M	=	500	412	ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA	<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<	MF:i:18	RG:Z:L2
    '''

    header = {'HD': {'VN': '1.0'},
              'SQ': [{'LN': 1575, 'SN': 'chr1'},
                     {'LN': 1584, 'SN': 'chr2'}], }

    bamfile = os.path.join(DATADIR, "ex6.bam")
    samfile = os.path.join(DATADIR, "ex6.sam")

    def setUp(self):

        a = pysam.AlignedRead()
        a.qname = "read_28833_29006_6945"
        a.seq = "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
        a.flag = 99
        a.rname = 0
        a.pos = 32
        a.mapq = 20
        a.cigar = ((0, 10), (2, 1), (0, 25))
        a.mrnm = 0
        a.mpos = 199
        a.isize = 167
        a.qual = "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<"
        a.tags = (("NM", 1),
                  ("RG", "L1"))

        b = pysam.AlignedRead()
        b.qname = "read_28701_28881_323b"
        b.seq = "ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA"
        b.flag = 147
        b.rname = 1
        b.pos = 87
        b.mapq = 30
        b.cigar = ((0, 35), )
        b.mrnm = 1
        b.mpos = 499
        b.isize = 412
        b.qual = "<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<"
        b.tags = (("MF", 18),
                  ("RG", "L2"))

        self.reads = (a, b)

    # TODO
    # def testSAMWholeFile(self):

    #     tmpfilename = "tmp_%i.sam" % id(self)

    #     outfile = pysam.Samfile(tmpfilename,
    #                             "wh",
    #                             header=self.header)

    #     for x in self.reads:
    #         outfile.write(x)
    #     outfile.close()
    #     self.assertTrue(checkBinaryEqual(tmpfilename, self.samfile),
    #                     "mismatch when construction SAM file, see %s %s" % (tmpfilename, self.samfile))

    #     os.unlink(tmpfilename)

    def testBAMPerRead(self):
        '''check if individual reads are binary equal.'''
        infile = pysam.Samfile(self.bamfile, "rb")

        others = list(infile)
        for denovo, other in zip(others, self.reads):
            self.checkFieldEqual(other, denovo)
            self.assertEqual(other.compare(denovo), 0)

    # TODO
    # def testSAMPerRead(self):
    #     '''check if individual reads are binary equal.'''
    #     infile = pysam.Samfile(self.samfile, "r")

    #     others = list(infile)
    #     for denovo, other in zip(others, self.reads):
    #         self.checkFieldEqual(other, denovo)
    #         self.assertEqual(other.compare(denovo), 0)

    def testBAMWholeFile(self):

        tmpfilename = "tmp_%i.bam" % id(self)

        outfile = pysam.Samfile(tmpfilename, "wb", header=self.header)

        for x in self.reads:
            outfile.write(x)
        outfile.close()

        self.assertTrue(checkBinaryEqual(tmpfilename, self.bamfile),
                        "mismatch when construction BAM file, see %s %s" % (tmpfilename, self.bamfile))

        os.unlink(tmpfilename)


class TestDeNovoConstructionUserTags(TestDeNovoConstruction):

    '''test de novo construction with a header that contains lower-case tags.'''

    header = {'HD': {'VN': '1.0'},
              'SQ': [{'LN': 1575, 'SN': 'chr1'},
                     {'LN': 1584, 'SN': 'chr2'}],
              'x1': {'A': 2, 'B': 5},
              'x3': {'A': 6, 'B': 5},
              'x2': {'A': 4, 'B': 5}}

    bamfile = os.path.join(DATADIR, "example_user_header.bam")
    samfile = os.path.join(DATADIR, "example_user_header.sam")


class TestEmptyHeader(unittest.TestCase):

    '''see issue 84.'''

    def testEmptyHeader(self):

        s = pysam.Samfile(os.path.join(DATADIR, 'example_empty_header.bam'))
        self.assertEqual(s.header, {'SQ': [{'LN': 1000, 'SN': 'chr1'}]})


class TestBTagSam(unittest.TestCase):

    '''see issue 81.'''

    compare = [[100, 1, 91, 0, 7, 101, 0, 201, 96, 204, 0, 0, 87, 109, 0, 7, 97, 112, 1, 12, 78, 197, 0, 7, 100, 95, 101, 202, 0, 6, 0, 1, 186, 0, 84, 0, 244, 0, 0, 324, 0, 107, 195, 101, 113, 0, 102, 0, 104, 3, 0, 101, 1, 0, 212, 6, 0, 0, 1, 0, 74, 1, 11, 0, 196, 2, 197, 103, 0, 108, 98, 2, 7, 0, 1, 2, 194, 0, 180, 0, 108, 0, 203, 104, 16, 5, 205, 0, 0, 0, 1, 1, 100, 98, 0, 0, 204, 6, 0, 79, 0, 0, 101, 7, 109, 90, 265, 1, 27, 10, 109, 102, 9, 0, 292, 0, 110, 0, 0, 102, 112, 0, 0, 84, 100, 103, 2, 81, 126, 0, 2, 90, 0, 15, 96, 15, 1, 0, 2, 0, 107, 92, 0, 0, 101, 3, 98, 15, 102, 13, 116, 116, 90, 93, 198, 0, 0, 0, 199, 92, 26, 495, 100, 5, 0, 100, 5, 209, 0, 92, 107, 90, 0, 0, 0, 0, 109, 194, 7, 94, 200, 0, 40, 197, 0, 11, 0, 0, 112, 110, 6, 4, 200, 28, 0, 196, 0, 203, 1, 129, 0, 0, 1, 0, 94, 0, 1, 0, 107, 5, 201, 3, 3, 100, 0, 121, 0, 7, 0, 1, 105, 306, 3, 86, 8, 183, 0, 12, 163, 17, 83, 22, 0, 0, 1, 8, 109, 103, 0, 0, 295, 0, 200, 16, 172, 3, 16, 182, 3, 11, 0, 0, 223, 111, 103, 0, 5, 225, 0, 95],
               [-100, 200, -300, -400],
               [-100, 12],
               [12, 15],
               [-1.0, 5.0, 2.5]]

    filename = os.path.join(DATADIR, 'example_btag.sam')

    def testRead(self):

        s = pysam.Samfile(self.filename)
        for x, read in enumerate(s):
            if x == 0:
                self.assertEqual(read.tags, [('RG', 'QW85I'), ('PG', 'tmap'), ('MD', '140'), ('NM', 0), ('AS', 140), ('FZ', [100, 1, 91, 0, 7, 101, 0, 201, 96, 204, 0, 0, 87, 109, 0, 7, 97, 112, 1, 12, 78, 197, 0, 7, 100, 95, 101, 202, 0, 6, 0, 1, 186, 0, 84, 0, 244, 0, 0, 324, 0, 107, 195, 101, 113, 0, 102, 0, 104, 3, 0, 101, 1, 0, 212, 6, 0, 0, 1, 0, 74, 1, 11, 0, 196, 2, 197, 103, 0, 108, 98, 2, 7, 0, 1, 2, 194, 0, 180, 0, 108, 0, 203, 104, 16, 5, 205, 0, 0, 0, 1, 1, 100, 98, 0, 0, 204, 6, 0, 79, 0, 0, 101, 7, 109, 90, 265, 1, 27, 10, 109, 102, 9, 0, 292, 0, 110, 0, 0, 102, 112, 0, 0, 84, 100, 103, 2, 81, 126, 0, 2, 90, 0, 15, 96, 15, 1, 0, 2, 0, 107, 92, 0, 0, 101, 3, 98, 15, 102, 13, 116, 116, 90, 93, 198, 0, 0, 0, 199, 92, 26, 495, 100, 5, 0, 100, 5, 209, 0, 92, 107, 90, 0, 0, 0, 0, 109, 194, 7, 94, 200, 0, 40, 197, 0, 11, 0, 0, 112, 110, 6, 4, 200, 28, 0, 196, 0, 203, 1, 129, 0, 0, 1, 0, 94, 0, 1, 0, 107, 5, 201, 3, 3, 100, 0, 121, 0, 7, 0, 1, 105, 306, 3, 86, 8, 183, 0, 12, 163, 17, 83, 22, 0, 0, 1, 8, 109, 103, 0, 0, 295, 0, 200, 16, 172, 3, 16, 182, 3, 11, 0, 0, 223, 111, 103, 0, 5, 225, 0, 95]), ('XA', 'map2-1'), ('XS', 53), ('XT', 38), ('XF', 1), ('XE', 0)]
                                 )

            fz = dict(read.tags)["FZ"]
            self.assertEqual(fz, self.compare[x])
            self.assertEqual(read.opt("FZ"), self.compare[x])

    def testWrite(self):

        s = pysam.Samfile(self.filename)
        for read in s:
            before = read.tags
            read.tags = read.tags
            after = read.tags
            self.assertEqual(after, before)


class TestBTagBam(TestBTagSam):
    filename = os.path.join(DATADIR, 'example_btag.bam')


class TestDoubleFetch(unittest.TestCase):

    '''check if two iterators on the same bamfile are independent.'''

    filename = os.path.join(DATADIR, 'ex1.bam')

    def testDoubleFetch(self):

        samfile1 = pysam.Samfile(self.filename, 'rb')

        for a, b in zip(samfile1.fetch(multiple_iterators=True),
                        samfile1.fetch(multiple_iterators=True)):
            self.assertEqual(a.compare(b), 0)

    def testDoubleFetchWithRegion(self):

        samfile1 = pysam.Samfile(self.filename, 'rb')
        chr, start, stop = 'chr1', 200, 3000000
        # just making sure the test has something to catch
        self.assertTrue(len(list(samfile1.fetch(chr, start, stop))) > 0)

        for a, b in zip(samfile1.fetch(chr, start, stop),
                        samfile1.fetch(chr, start, stop,
                                       multiple_iterators=True)):
            self.assertEqual(a.compare(b), 0)

    def testDoubleFetchUntilEOF(self):

        samfile1 = pysam.Samfile(self.filename, 'rb')

        for a, b in zip(samfile1.fetch(until_eof=True),
                        samfile1.fetch(until_eof=True,
                                       multiple_iterators=True)):
            self.assertEqual(a.compare(b), 0)


class TestRemoteFileFTP(unittest.TestCase):

    '''test remote access.

    '''

    # Need to find an ftp server without password on standard
    # port.

    url = "ftp://ftp.sanger.ac.uk/pub/rd/humanSequences/CV.bam"
    region = "1:1-1000"

    def testFTPView(self):
        return
        if not checkURL(self.url):
            return

        result = pysam.view(self.url, self.region)
        self.assertEqual(len(result), 36)

    def testFTPFetch(self):
        return
        if not checkURL(self.url):
            return

        samfile = pysam.Samfile(self.url, "rb")
        result = list(samfile.fetch(region=self.region))
        self.assertEqual(len(result), 36)


class TestRemoteFileHTTP(unittest.TestCase):

    url = "http://genserv.anat.ox.ac.uk/downloads/pysam/test/ex1.bam"
    region = "chr1:1-1000"
    local = os.path.join(DATADIR, "ex1.bam")

    def testView(self):
        if not checkURL(self.url):
            return

        samfile_local = pysam.Samfile(self.local, "rb")
        ref = list(samfile_local.fetch(region=self.region))

        result = pysam.view(self.url, self.region)
        self.assertEqual(len(result), len(ref))

    def testFetch(self):
        if not checkURL(self.url):
            return

        samfile = pysam.Samfile(self.url, "rb")
        result = list(samfile.fetch(region=self.region))
        samfile_local = pysam.Samfile(self.local, "rb")
        ref = list(samfile_local.fetch(region=self.region))

        self.assertEqual(len(ref), len(result))
        for x, y in zip(result, ref):
            self.assertEqual(x.compare(y), 0)

    def testFetchAll(self):
        if not checkURL(self.url):
            return

        samfile = pysam.Samfile(self.url, "rb")
        result = list(samfile.fetch())
        samfile_local = pysam.Samfile(self.local, "rb")
        ref = list(samfile_local.fetch())

        self.assertEqual(len(ref), len(result))
        for x, y in zip(result, ref):
            self.assertEqual(x.compare(y), 0)


class TestLargeOptValues(unittest.TestCase):

    ints = (65536, 214748, 2147484, 2147483647)
    floats = (65536.0, 214748.0, 2147484.0)

    def check(self, samfile):

        i = samfile.fetch()
        for exp in self.ints:
            rr = next(i)
            obs = rr.opt("ZP")
            self.assertEqual(exp, obs,
                             "expected %s, got %s\n%s" %
                             (str(exp), str(obs), str(rr)))

        for exp in [-x for x in self.ints]:
            rr = next(i)
            obs = rr.opt("ZP")
            self.assertEqual(exp, obs,
                             "expected %s, got %s\n%s" %
                             (str(exp), str(obs), str(rr)))

        for exp in self.floats:
            rr = next(i)
            obs = rr.opt("ZP")
            self.assertEqual(exp, obs,
                             "expected %s, got %s\n%s" %
                             (str(exp), str(obs), str(rr)))

        for exp in [-x for x in self.floats]:
            rr = next(i)
            obs = rr.opt("ZP")
            self.assertEqual(exp, obs, "expected %s, got %s\n%s" %
                             (str(exp), str(obs), str(rr)))

    def testSAM(self):
        samfile = pysam.Samfile(
            os.path.join(DATADIR, "ex10.sam"),
            "r")
        self.check(samfile)

    def testBAM(self):
        samfile = pysam.Samfile(
            os.path.join(DATADIR, "ex10.bam"),
            "rb")
        self.check(samfile)


class TestPileup(unittest.TestCase):

    '''test pileup functionality.'''

    samfilename = "pysam_data/ex1.bam"
    fastafilename = "pysam_data/ex1.fa"

    def setUp(self):

        self.samfile = pysam.Samfile(self.samfilename)
        self.fastafile = pysam.Fastafile(self.fastafilename)

    def checkEqual(self, references, iterator):

        for x, column in enumerate(iterator):
            (contig, pos, reference_base,
             read_bases, read_qualities, alignment_mapping_qualities) \
                = references[x][:-1].split("\t")
            self.assertEqual(int(pos) - 1, column.pos)

    def testSamtoolsStepper(self):
        refs = pysam.mpileup(
            "-f", self.fastafilename,
            self.samfilename)
        iterator = self.samfile.pileup(
            stepper="samtools",
            fastafile=self.fastafile)
        self.checkEqual(refs, iterator)

    def testAllStepper(self):
        refs = pysam.mpileup(
            "-f", self.fastafilename,
            "-A", "-B",
            self.samfilename)

        iterator = self.samfile.pileup(
            stepper="all",
            fastafile=self.fastafile)
        self.checkEqual(refs, iterator)


class TestLogging(unittest.TestCase):

    '''test around bug issue 42,

    failed in versions < 0.4
    '''

    def check(self, bamfile, log):

        if log:
            logger = logging.getLogger('franklin')
            logger.setLevel(logging.INFO)
            formatter = logging.Formatter(
                '%(asctime)s %(levelname)s %(message)s')
            log_hand = logging.FileHandler('log.txt')
            log_hand.setFormatter(formatter)
            logger.addHandler(log_hand)

        bam = pysam.Samfile(bamfile, 'rb')
        cols = bam.pileup()
        self.assertTrue(True)

    def testFail1(self):
        self.check(os.path.join(DATADIR, "ex9_fail.bam"),
                   False)
        self.check(os.path.join(DATADIR, "ex9_fail.bam"),
                   True)

    def testNoFail1(self):
        self.check(os.path.join(DATADIR, "ex9_nofail.bam"),
                   False)
        self.check(os.path.join(DATADIR, "ex9_nofail.bam"),
                   True)

    def testNoFail2(self):
        self.check(os.path.join(DATADIR, "ex9_nofail.bam"),
                   True)
        self.check(os.path.join(DATADIR, "ex9_nofail.bam"),
                   True)

# TODOS
# 1. finish testing all properties within pileup objects
# 2. check exceptions and bad input problems (missing files, optional fields that aren't present, etc...)
# 3. check: presence of sequence


class TestSamfileUtilityFunctions(unittest.TestCase):

    def testCount(self):

        samfile = pysam.Samfile(os.path.join(DATADIR, "ex1.bam"),
                                "rb")

        for contig in ("chr1", "chr2"):
            for start in range(0, 2000, 100):
                end = start + 1
                self.assertEqual(
                    len(list(samfile.fetch(contig, start, end))),
                    samfile.count(contig, start, end),
                    'number mismatch for %s:%i-%i %i != %i' % (
                        contig, start, end,
                        len(list(samfile.fetch(contig, start, end))),
                        samfile.count(contig, start, end)))

                # test empty intervals
                self.assertEqual(
                    len(list(samfile.fetch(contig, start, start))),
                    samfile.count(contig, start, start),
                    'number mismatch for %s:%i-%i %i != %i' % (
                        contig, start, start,
                        len(list(samfile.fetch(contig, start, start))),
                        samfile.count(contig, start, start)))

                # test half empty intervals
                self.assertEqual(len(list(samfile.fetch(contig, start))),
                                 samfile.count(contig, start))

                self.assertEqual(
                    len(list(samfile.fetch(contig, start))),
                    samfile.count(contig, start),
                    'number mismatch for %s:%i %i != %i' % (
                        contig, start,
                        len(list(samfile.fetch(contig, start))),
                        samfile.count(contig, start)))

    def testMate(self):
        '''test mate access.'''

        with open(os.path.join(DATADIR, "ex1.sam"), "rb") as inf:
            readnames = [x.split(b"\t")[0] for x in inf.readlines()]
        if sys.version_info[0] >= 3:
            readnames = [name.decode('ascii') for name in readnames]

        counts = collections.defaultdict(int)
        for x in readnames:
            counts[x] += 1

        samfile = pysam.Samfile(os.path.join(DATADIR, "ex1.bam"),
                                "rb")

        for read in samfile.fetch():
            if not read.is_paired:
                self.assertRaises(ValueError, samfile.mate, read)
            elif read.mate_is_unmapped:
                self.assertRaises(ValueError, samfile.mate, read)
            else:
                if counts[read.qname] == 1:
                    self.assertRaises(ValueError, samfile.mate, read)
                else:
                    mate = samfile.mate(read)
                    self.assertEqual(read.qname, mate.qname)
                    self.assertEqual(read.is_read1, mate.is_read2)
                    self.assertEqual(read.is_read2, mate.is_read1)
                    self.assertEqual(read.pos, mate.mpos)
                    self.assertEqual(read.mpos, mate.pos)

    def testIndexStats(self):
        '''test if total number of mapped/unmapped reads is correct.'''

        samfile = pysam.Samfile(os.path.join(DATADIR, "ex1.bam"),
                                "rb")
        self.assertEqual(samfile.mapped, 3235)
        self.assertEqual(samfile.unmapped, 35)
        self.assertEqual(samfile.nocoordinate, 0)


class TestSamtoolsProxy(unittest.TestCase):

    '''tests for sanity checking access to samtools functions.'''

    def testIndex(self):
        self.assertRaises(IOError, pysam.index, "missing_file")

    def testView(self):
        # note that view still echos "open: No such file or directory"
        self.assertRaises(pysam.SamtoolsError, pysam.view, "missing_file")

    def testSort(self):
        self.assertRaises(pysam.SamtoolsError, pysam.sort, "missing_file")


class TestSamfileIndex(unittest.TestCase):

    def testIndex(self):
        samfile = pysam.Samfile(os.path.join(DATADIR, "ex1.bam"),
                                "rb")
        index = pysam.IndexedReads(samfile)
        index.build()
        reads = collections.defaultdict(int)

        for read in samfile:
            reads[read.qname] += 1

        for qname, counts in reads.items():
            found = list(index.find(qname))
            self.assertEqual(len(found), counts)
            for x in found:
                self.assertEqual(x.qname, qname)


if __name__ == "__main__":
    # build data files
    print ("building data files")
    subprocess.call("make -C %s" % DATADIR, shell=True)
    print ("starting tests")
    unittest.main()
    print ("completed tests")
