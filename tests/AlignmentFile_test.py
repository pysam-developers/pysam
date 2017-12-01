#!/usr/bin/env python
'''unit testing code for pysam.

Execute in the :file:`tests` directory as it requires the Makefile
and data files located there.
'''

import unittest
import os
import shutil
import sys
import re
import copy
import collections
import subprocess
import logging
import array
if sys.version_info.major >= 3:
    from io import StringIO
else:
    from StringIO import StringIO

from functools import partial

import pysam
import pysam.samtools
from TestUtils import checkBinaryEqual, checkURL, \
    check_samtools_view_equal, checkFieldEqual, force_str, \
    get_temp_filename, BAM_DATADIR


##################################################
#
# Detailed test of file contents
#
# Data are read either through file based iterator
# access (BasicTestBAMFromFile) or by calling fetch
# without coordinates (BasicTestBAMFromFetch)
##################################################
class BasicTestBAMFromFetch(unittest.TestCase):

    '''basic first test - detailed testing
    if information in file is consistent
    with information in AlignedSegment object.'''

    def setUp(self):
        self.samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex3.bam"),
            "rb")
        self.reads = list(self.samfile.fetch())

    def tearDown(self):
        self.samfile.close()

    def testARqname(self):
        self.assertEqual(
            self.reads[0].query_name,
            "read_28833_29006_6945",
            "read name mismatch in read 1: %s != %s" % (
                self.reads[0].query_name, "read_28833_29006_6945"))
        self.assertEqual(
            self.reads[1].query_name,
            "read_28701_28881_323b",
            "read name mismatch in read 2: %s != %s" % (
                self.reads[1].query_name, "read_28701_28881_323b"))

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
            self.reads[0].reference_id, 0,
            "chromosome/target id mismatch in read 1: %s != %s" %
            (self.reads[0].reference_id, 0))
        self.assertEqual(
            self.reads[1].reference_id, 1,
            "chromosome/target id mismatch in read 2: %s != %s" %
            (self.reads[1].reference_id, 1))

    def testARpos(self):
        self.assertEqual(
            self.reads[0].reference_start, 33 - 1,
            "mapping position mismatch in read 1: %s != %s" %
            (self.reads[0].reference_start, 33 - 1))
        self.assertEqual(
            self.reads[1].reference_start, 88 - 1,
            "mapping position mismatch in read 2: %s != %s" %
            (self.reads[1].reference_start, 88 - 1))

    def testARmapq(self):
        self.assertEqual(
            self.reads[0].mapping_quality, 20,
            "mapping quality mismatch in read 1: %s != %s" %
            (self.reads[0].mapping_quality, 20))
        self.assertEqual(
            self.reads[1].mapping_quality, 30,
            "mapping quality mismatch in read 2: %s != %s" % (
                self.reads[1].mapping_quality, 30))

    def testARcigar(self):
        self.assertEqual(
            self.reads[0].cigartuples,
            [(0, 10), (2, 1), (0, 25)],
            "read name length mismatch in read 1: %s != %s" %
            (self.reads[0].cigartuples, [(0, 10), (2, 1), (0, 25)]))
        self.assertEqual(
            self.reads[1].cigartuples, [(0, 35)],
            "read name length mismatch in read 2: %s != %s" %
            (self.reads[1].cigartuples, [(0, 35)]))

    def testARcigarstring(self):
        self.assertEqual(self.reads[0].cigarstring, '10M1D25M')
        self.assertEqual(self.reads[1].cigarstring, '35M')

    def testARmrnm(self):
        self.assertEqual(
            self.reads[0].next_reference_id, 0,
            "mate reference sequence name mismatch in read 1: %s != %s" %
            (self.reads[0].next_reference_id, 0))
        self.assertEqual(
            self.reads[1].next_reference_id, 1,
            "mate reference sequence name mismatch in read 2: %s != %s" %
            (self.reads[1].next_reference_id, 1))
        self.assertEqual(
            self.reads[0].next_reference_id, 0,
            "mate reference sequence name mismatch in read 1: %s != %s" %
            (self.reads[0].next_reference_id, 0))
        self.assertEqual(
            self.reads[1].next_reference_id, 1,
            "mate reference sequence name mismatch in read 2: %s != %s" %
            (self.reads[1].next_reference_id, 1))

    def testARmpos(self):
        self.assertEqual(self.reads[
                         0].next_reference_start, 200 - 1, "mate mapping position mismatch in read 1: %s != %s" % (self.reads[0].next_reference_start, 200 - 1))
        self.assertEqual(self.reads[
                         1].next_reference_start, 500 - 1, "mate mapping position mismatch in read 2: %s != %s" % (self.reads[1].next_reference_start, 500 - 1))
        self.assertEqual(self.reads[
                         0].next_reference_start, 200 - 1, "mate mapping position mismatch in read 1: %s != %s" % (self.reads[0].next_reference_start, 200 - 1))
        self.assertEqual(self.reads[
                         1].next_reference_start, 500 - 1, "mate mapping position mismatch in read 2: %s != %s" % (self.reads[1].next_reference_start, 500 - 1))

    def testARQueryLength(self):
        self.assertEqual(
            self.reads[0].query_length, 35,
            "insert size mismatch in read 1: %s != %s" %
            (self.reads[0].query_length, 35))
        self.assertEqual(
            self.reads[1].query_length, 35,
            "insert size mismatch in read 2: %s != %s" %
            (self.reads[1].query_length, 35))
        self.assertEqual(
            self.reads[0].query_length, 35,
            "insert size mismatch in read 1: %s != %s" %
            (self.reads[0].query_length, 35))
        self.assertEqual(
            self.reads[1].query_length, 35,
            "insert size mismatch in read 2: %s != %s" %
            (self.reads[1].query_length, 35))

    def testARseq(self):
        self.assertEqual(self.reads[0].query_sequence, "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG", "sequence mismatch in read 1: %s != %s" % (
            self.reads[0].query_sequence, "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"))
        self.assertEqual(self.reads[1].query_sequence, "ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA", "sequence size mismatch in read 2: %s != %s" % (
            self.reads[1].query_sequence, "ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA"))
        self.assertEqual(self.reads[3].query_sequence, "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG", "sequence mismatch in read 4: %s != %s" % (
            self.reads[3].query_sequence, "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"))

    def testARqual(self):
        self.assertEqual(
            pysam.qualities_to_qualitystring(self.reads[0].query_qualities),
            "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
        self.assertEqual(
            pysam.qualities_to_qualitystring(self.reads[1].query_qualities),
            "<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<")
        self.assertEqual(
            pysam.qualities_to_qualitystring(self.reads[3].query_qualities),
            "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")

    def testARquery(self):
        self.assertEqual(
            self.reads[0].query_alignment_sequence,
            "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG",
            "query mismatch in read 1: %s != %s" %
            (self.reads[0].query_alignment_sequence,
             "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"))
        self.assertEqual(
            self.reads[1].query_alignment_sequence,
            "ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA",
            "query size mismatch in read 2: %s != %s" %
            (self.reads[1].query_alignment_sequence,
             "ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA"))
        self.assertEqual(
            self.reads[3].query_alignment_sequence,
            "TAGCTAGCTACCTATATCTTGGTCTT",
            "query mismatch in read 4: %s != %s" %
            (self.reads[3].query_alignment_sequence,
             "TAGCTAGCTACCTATATCTTGGTCTT"))

    def testARqqual(self):
        self.assertEqual(
            pysam.qualities_to_qualitystring(self.reads[0].query_alignment_qualities),
            "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<",
            "qquality string mismatch in read 1: %s != %s" %
            (pysam.qualities_to_qualitystring(self.reads[0].query_alignment_qualities),
             "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<"))
        self.assertEqual(
            pysam.qualities_to_qualitystring(self.reads[1].query_alignment_qualities),
            "<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<",
            "qquality string mismatch in read 2: %s != %s" %
            (pysam.qualities_to_qualitystring(self.reads[1].query_alignment_qualities),
             "<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<"))
        self.assertEqual(
            pysam.qualities_to_qualitystring(self.reads[3].query_alignment_qualities),
            "<<<<<<<<<<<<<<<<<:<9/,&,22",
            "qquality string mismatch in read 3: %s != %s" %
            (pysam.qualities_to_qualitystring(self.reads[3].query_alignment_qualities),
             "<<<<<<<<<<<<<<<<<:<9/,&,22"))

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
        self.assertEqual(self.reads[0].is_paired, True, "is paired mismatch in read 1: %s != %s" % (
            self.reads[0].is_paired, True))
        self.assertEqual(self.reads[1].is_paired, True, "is paired mismatch in read 2: %s != %s" % (
            self.reads[1].is_paired, True))
        self.assertEqual(self.reads[0].is_proper_pair, True, "is proper pair mismatch in read 1: %s != %s" % (
            self.reads[0].is_proper_pair, True))
        self.assertEqual(self.reads[1].is_proper_pair, True, "is proper pair mismatch in read 2: %s != %s" % (
            self.reads[1].is_proper_pair, True))

    def testTags(self):
        self.assertEqual(self.reads[0].tags,
                         [('NM', 1), ('RG', 'L1'),
                          ('PG', 'P1'), ('XT', 'U')])
        self.assertEqual(self.reads[1].tags,
                         [('MF', 18), ('RG', 'L2'),
                          ('PG', 'P2'), ('XT', 'R')])

    def testOpt(self):
        self.assertEqual(self.reads[0].opt("XT"), "U")
        self.assertEqual(self.reads[1].opt("XT"), "R")


class BasicTestSAMFromFetch(BasicTestBAMFromFetch):

    def setUp(self):
        self.samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex3.sam"),
            "r")
        self.reads = list(self.samfile.fetch())


class BasicTestCRAMFromFetch(BasicTestBAMFromFetch):

    def setUp(self):
        self.samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex3.cram"),
            "rc")
        self.reads = list(self.samfile.fetch())

    def testTags(self):
        self.assertEqual(
            sorted(self.reads[0].tags),
            sorted([('RG', 'L1'),
                    ('NM', 22),
                    ('MD', '0C0T1G1C0C0A1G0^G0C1C1G1A0T2G0G0G0A1C1G1G1A2C0'),
                    ('PG', 'P1'),
                    ('XT', 'U'),
                    ]))
        self.assertEqual(
            sorted(self.reads[1].tags),
            sorted([('RG', 'L2'),
                    ('NM', 26),
                    ('MD',
                     '1G0A0A1G1G0G2C0A0G0A0A0C0T0T0G0A0A0G0A0C0A0A1T2C0T0T1'),
                    ('MF', 18),
                    ('PG', 'P2'),
                    ('XT', 'R')]))

    def testPresentOptionalFields(self):
        self.assertEqual(
            self.reads[0].opt('NM'), 22,
            "optional field mismatch in read 1, NM: %s != %s" %
            (self.reads[0].opt('NM'), 22))
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


class BasicTestSAMFromFilename(BasicTestBAMFromFetch):

    def setUp(self):
        self.samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex3.sam"),
            "r")
        self.reads = [r for r in self.samfile]


class BasicTestCRAMFromFilename(BasicTestCRAMFromFetch):

    def setUp(self):
        self.samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex3.cram"),
            "rc")
        self.reads = [r for r in self.samfile]


class BasicTestBAMFromFilename(BasicTestBAMFromFetch):

    def setUp(self):
        self.samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex3.bam"),
            "rb")
        self.reads = [r for r in self.samfile]


class BasicTestBAMFromFile(BasicTestBAMFromFetch):

    def setUp(self):
        with open(os.path.join(BAM_DATADIR, "ex3.bam")) as f:
            self.samfile = pysam.AlignmentFile(
                f, "rb")
        self.reads = [r for r in self.samfile]


class BasicTestBAMFromFileNo(BasicTestBAMFromFetch):

    def setUp(self):
        with open(os.path.join(BAM_DATADIR, "ex3.bam")) as f:
            self.samfile = pysam.AlignmentFile(
                f.fileno(), "rb")
        self.reads = [r for r in self.samfile]


class BasicTestSAMFromFile(BasicTestBAMFromFetch):

    def setUp(self):
        with open(os.path.join(BAM_DATADIR, "ex3.sam")) as f:
            self.samfile = pysam.AlignmentFile(
                f, "r")
        self.reads = [r for r in self.samfile]


class BasicTestSAMFromFileNo(BasicTestBAMFromFetch):

    def setUp(self):
        with open(os.path.join(BAM_DATADIR, "ex3.sam")) as f:
            self.samfile = pysam.AlignmentFile(
                f.fileno(), "r")
        self.reads = [r for r in self.samfile]


class BasicTestCRAMFromFile(BasicTestCRAMFromFetch):

    def setUp(self):
        with open(os.path.join(BAM_DATADIR, "ex3.cram")) as f:
            self.samfile = pysam.AlignmentFile(f, "rc")
        self.reads = [r for r in self.samfile]


class BasicTestCRAMFromFileNo(BasicTestCRAMFromFetch):

    def setUp(self):
        with open(os.path.join(BAM_DATADIR, "ex3.cram")) as f:
            self.samfile = pysam.AlignmentFile(
                f.fileno(), "rc")
        self.reads = [r for r in self.samfile]


class BasicTestSAMFromStringIO(BasicTestBAMFromFetch):

    def testRaises(self):
        statement = "samtools view -h {}".format(
                os.path.join(BAM_DATADIR, "ex3.bam"))
        stdout = subprocess.check_output(statement.split(" "))
        bam = StringIO()
        if sys.version_info.major >= 3:
            bam.write(stdout.decode('ascii'))
        else:
            bam.write(stdout)
        bam.seek(0)
        self.assertRaises(NotImplementedError,
                          pysam.AlignmentFile, bam)
        # self.reads = [r for r in samfile]


##################################################
#
# Test of basic File I/O
#
# * format conversions
# * reading with/without index
# * reading from closed files
#
##################################################
class TestIO(unittest.TestCase):

    '''check if reading samfile and writing a samfile
    are consistent.'''

    def checkEcho(self,
                  input_filename,
                  reference_filename,
                  output_filename,
                  input_mode,
                  output_mode,
                  sequence_filename=None,
                  use_template=True,
                  checkf=checkBinaryEqual,
                  **kwargs):
        '''iterate through *input_filename* writing to
        *output_filename* and comparing the output to
        *reference_filename*.

        The files are opened according to the *input_mode* and
        *output_mode*.

        If *use_template* is set, the header is copied from infile
        using the template mechanism, otherwise target names and
        lengths are passed explicitly.

        The *checkf* is used to determine if the files are
        equal.
        '''

        with pysam.AlignmentFile(
                os.path.join(BAM_DATADIR, input_filename),
                input_mode) as infile:

            if "b" in input_mode:
                self.assertTrue(infile.is_bam)
                self.assertFalse(infile.is_cram)
            elif "c" in input_mode:
                self.assertFalse(infile.is_bam)
                self.assertTrue(infile.is_cram)
            else:
                self.assertFalse(infile.is_cram)
                self.assertFalse(infile.is_bam)

            if use_template:
                outfile = pysam.AlignmentFile(
                    output_filename,
                    output_mode,
                    reference_filename=sequence_filename,
                    template=infile, **kwargs)
            else:
                outfile = pysam.AlignmentFile(
                    output_filename,
                    output_mode,
                    reference_names=infile.references,
                    reference_lengths=infile.lengths,
                    reference_filename=sequence_filename,
                    add_sq_text=False,
                    **kwargs)

            iter = infile.fetch()

            for x in iter:
                outfile.write(x)

            outfile.close()

        self.assertTrue(checkf(
            os.path.join(BAM_DATADIR, reference_filename),
            output_filename),
            "files %s and %s are not the same" %
            (reference_filename,
             output_filename))

        os.unlink(output_filename)

    def testSAM2SAM(self):
        self.checkEcho("ex2.sam",
                       "ex2.sam",
                       "tmp_ex2.sam",
                       "r", "wh")

    def testSAM2SAMWithoutHeader(self):
        self.checkEcho("ex2.sam",
                       "ex1.sam",
                       "tmp_ex2.sam",
                       "r", "w",
                       add_sam_header=False)

    def testBAM2BAM(self):
        self.checkEcho("ex2.bam",
                       "ex2.bam",
                       "tmp_ex2.bam",
                       "rb", "wb")

    def testCRAM2CRAM(self):
        self.checkEcho("ex2.cram",
                       "ex2.cram",
                       "tmp_ex2.cram",
                       "rc", "wc",
                       sequence_filename=os.path.join(BAM_DATADIR, "ex1.fa"),
                       checkf=check_samtools_view_equal)

    def testSAM2BAM(self):
        self.checkEcho("ex2.sam",
                       "ex2.bam",
                       "tmp_ex2.bam",
                       "r", "wb")

    def testBAM2SAM(self):
        self.checkEcho("ex2.bam",
                       "ex2.sam",
                       "tmp_ex2.sam",
                       "rb", "wh")

    def testBAM2CRAM(self):
        # ignore header (md5 sum)
        self.checkEcho("ex2.bam",
                       "ex2.cram",
                       "tmp_ex2.cram",
                       "rb", "wc",
                       sequence_filename=os.path.join(BAM_DATADIR, "ex1.fa"),
                       checkf=partial(
                           check_samtools_view_equal,
                           without_header=True))

    def testCRAM2BAM(self):
        # ignore header (md5 sum)
        self.checkEcho("ex2.cram",
                       "ex2.bam",
                       "tmp_ex2.bam",
                       "rc", "wb",
                       sequence_filename=os.path.join(BAM_DATADIR, "ex1.fa"),
                       checkf=partial(
                           check_samtools_view_equal,
                           without_header=True))

    def testSAM2CRAM(self):
        self.checkEcho("ex2.sam",
                       "ex2.cram",
                       "tmp_ex2.cram",
                       "r", "wc",
                       sequence_filename=os.path.join(BAM_DATADIR, "ex1.fa"),
                       checkf=partial(
                           check_samtools_view_equal,
                           without_header=True))

    def testCRAM2SAM(self):
        self.checkEcho("ex2.cram",
                       "ex2.sam",
                       "tmp_ex2.sam",
                       "rc", "wh",
                       sequence_filename=os.path.join(BAM_DATADIR, "ex1.fa"),
                       checkf=partial(
                           check_samtools_view_equal,
                           without_header=True))

    # Disabled - should work, files are not binary equal, but are
    # non-binary equal:
    # diff <(samtools view pysam_ex1.bam) <(samtools view pysam_data/ex1.bam)
    # def testReadWriteBamWithTargetNames(self):
    #     input_filename = "ex1.bam"
    #     output_filename = "pysam_ex1.bam"
    #     reference_filename = "ex1.bam"

    #     self.checkEcho(input_filename, reference_filename, output_filename,
    #                    "rb", "wb", use_template=False)

    def testReadSamWithoutTargetNames(self):
        '''see issue 104.'''
        input_filename = os.path.join(
            BAM_DATADIR,
            "example_unmapped_reads_no_sq.sam")

        # raise exception in default mode
        self.assertRaises(ValueError,
                          pysam.AlignmentFile,
                          input_filename,
                          "r")

        # raise exception if no SQ files
        self.assertRaises(ValueError,
                          pysam.AlignmentFile,
                          input_filename, "r",
                          check_header=True)

        with pysam.AlignmentFile(
            input_filename,
            check_header=False,
            check_sq=False) as infile:
            result = list(infile.fetch(until_eof=True))
            self.assertEqual(2, len(result))

    def testReadBamWithoutTargetNames(self):
        '''see issue 104.'''
        input_filename = os.path.join(
            BAM_DATADIR, "example_unmapped_reads_no_sq.bam")

        # raise exception in default mode
        self.assertRaises(ValueError,
                          pysam.AlignmentFile,
                          input_filename,
                          "r")

        # raise exception if no SQ files
        self.assertRaises(ValueError,
                          pysam.AlignmentFile,
                          input_filename,
                          "r",
                          check_header=True)

        with pysam.AlignmentFile(
            input_filename, check_sq=False) as infile:
            result = list(infile.fetch(until_eof=True))

    def test_fail_read_sam_without_header(self):
        input_filename = os.path.join(BAM_DATADIR, "ex1.sam")

        self.assertRaises(ValueError,
                          pysam.AlignmentFile,
                          input_filename,
                          "r")

    def test_pass_read_sam_without_header_with_refs(self):
        with pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex1.sam"),
                                 "r",
                                 reference_names=["chr1", "chr2"],
                                 reference_lengths=[1575, 1584]) as samfile:
            self.assertEqual(len(list(samfile.fetch(until_eof=True))), 3270)

    def test_pass_read_sam_with_header_without_header_check(self):
        with pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex2.sam"),
                                 "r", check_header=False) as samfile:
            self.assertEqual(len(list(samfile.fetch(until_eof=True))), 3270)

    def test_fail_when_reading_unformatted_files(self):
        '''test reading from a file that is not bam/sam formatted'''
        input_filename = os.path.join(BAM_DATADIR, 'Makefile')

        self.assertRaises(ValueError,
                          pysam.AlignmentFile,
                          input_filename,
                          "rb")

        self.assertRaises(ValueError,
                          pysam.AlignmentFile,
                          input_filename,
                          "r")

    def testBAMWithoutAlignedSegments(self):
        '''see issue 117'''
        input_filename = os.path.join(BAM_DATADIR, "test_unaligned.bam")
        samfile = pysam.AlignmentFile(input_filename,
                                      "rb",
                                      check_sq=False)
        samfile.fetch(until_eof=True)

    def testBAMWithShortBAI(self):
        '''see issue 116'''
        input_filename = os.path.join(BAM_DATADIR, "example_bai.bam")
        samfile = pysam.AlignmentFile(input_filename,
                                      "rb",
                                      check_sq=False)
        samfile.fetch('chr2')

    def testFetchFromClosedFile(self):

        samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex1.bam"),
            "rb")
        samfile.close()
        self.assertRaises(ValueError, samfile.fetch, 'chr1', 100, 120)

    def testFetchFromClosedFileObject(self):

        f = open(os.path.join(BAM_DATADIR, "ex1.bam"))
        samfile = pysam.AlignmentFile(f, "rb")
        f.close()
        self.assertTrue(f.closed)
        # access to Samfile still works
        self.checkEcho("ex1.bam",
                       "ex1.bam",
                       "tmp_ex1.bam",
                       "rb", "wb")

        f = open(os.path.join(BAM_DATADIR, "ex1.bam"))
        samfile = pysam.AlignmentFile(f, "rb")
        self.assertFalse(f.closed)
        samfile.close()
        # python file needs to be closed separately
        self.assertFalse(f.closed)

    def testClosedFile(self):
        '''test that access to a closed samfile raises ValueError.'''

        samfile = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex1.bam"),
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


    # TOOD
    # def testReadingFromSamFileWithoutHeader(self):
    #     '''read from samfile without header.
    #     '''
    #     samfile = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex7.sam"),
    #                             check_header=False,
    #                             check_sq=False)
    #     self.assertRaises(NotImplementedError, samfile.__iter__)

    def testReadingFromFileWithoutIndex(self):
        '''read from bam file without index.'''

        shutil.copyfile(os.path.join(BAM_DATADIR, "ex2.bam"),
                        'tests/tmp_ex2.bam')
        samfile = pysam.AlignmentFile('tests/tmp_ex2.bam',
                                      "rb")
        self.assertRaises(ValueError, samfile.fetch)
        self.assertEqual(
            len(list(samfile.fetch(until_eof=True))),
            3270)
        os.unlink('tests/tmp_ex2.bam')

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
        samfile = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex1.bam"),
                                      "rb")
        l10 = list(samfile.head(10))
        l100 = list(samfile.head(100))
        self.assertEqual(len(l10), 10)
        self.assertEqual(len(l100), 100)
        self.assertEqual(list(map(str, l10)),
                         list(map(str, l100[:10])))

    def testWriteUncompressedBAMFile(self):
        '''read from uncompressed BAM file, see issue #43'''

        input_filename = "ex2.sam"
        output_filename = "pysam_uncompressed.bam"
        reference_filename = "uncompressed.bam"

        self.checkEcho(input_filename,
                       reference_filename,
                       output_filename,
                       "r", "wb0")

        self.checkEcho(input_filename,
                       reference_filename,
                       output_filename,
                       "r", "wbu")

    def testEmptyBAM(self):
        samfile = pysam.Samfile(os.path.join(BAM_DATADIR, "empty.bam"),
                                "rb")
        self.assertEqual(samfile.mapped, 0)
        self.assertEqual(samfile.unmapped, 0)
        self.assertEqual(samfile.nocoordinate, 0)

    def testEmptyWithHeaderBAM(self):
        self.assertRaises(
            ValueError,
            pysam.Samfile,
            os.path.join(BAM_DATADIR, "example_empty_with_header.bam"),
            "rb")

        samfile = pysam.Samfile(
            os.path.join(BAM_DATADIR, "example_empty_with_header.bam"),
            "rb",
            check_sq=False)
        self.assertEqual(samfile.mapped, 0)
        self.assertEqual(samfile.unmapped, 0)
        self.assertEqual(samfile.nocoordinate, 0)
        self.assertEqual([], list(samfile.fetch()))

    def testOpenFromFilename(self):

        samfile = pysam.AlignmentFile(
            filename=os.path.join(BAM_DATADIR, "ex1.bam"),
            mode="rb")
        self.assertEqual(len(list(samfile.fetch())), 3270)

    def testBAMWithCSIIndex(self):
        '''see issue 116'''
        input_filename = os.path.join(BAM_DATADIR, "ex1_csi.bam")
        samfile = pysam.AlignmentFile(input_filename,
                                      "rb",
                                      check_sq=False)
        samfile.fetch('chr2')

    def test_fetch_by_tid(self):
        with pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex1.bam"), "rb") as samfile:
            self.assertEqual(len(list(samfile.fetch('chr1'))),
                             len(list(samfile.fetch(tid=0))))
            self.assertEqual(len(list(samfile.fetch('chr2'))),
                             len(list(samfile.fetch(tid=1))))
            self.assertRaises(
                IndexError,
                samfile.fetch,
                tid=2)
            self.assertRaises(
                IndexError,
                samfile.fetch,
                tid=-1)
            self.assertEqual(len(list(samfile.fetch('chr1',start=1000, end=2000))),
                             len(list(samfile.fetch(tid=0, start=1000, end=2000))))


class TestAutoDetect(unittest.TestCase):

    def testSAM(self):
        """test SAM autodetection."""

        with pysam.AlignmentFile(
                os.path.join(BAM_DATADIR, "ex3.sam")) as inf:
            self.assertFalse(inf.is_bam)
            self.assertFalse(inf.is_cram)

            self.assertRaises(ValueError, inf.fetch, 'chr1')

    def testBAM(self):
        """test BAM autodetection."""

        with pysam.AlignmentFile(
                os.path.join(BAM_DATADIR, "ex3.bam")) as inf:
            self.assertTrue(inf.is_bam)
            self.assertFalse(inf.is_cram)
            self.assertEqual(len(list(inf.fetch('chr1'))), 1)
            self.assertEqual(len(list(inf.fetch('chr2'))), 3)

    def testCRAM(self):
        """test CRAM autodetection."""

        with pysam.AlignmentFile(
                os.path.join(BAM_DATADIR, "ex3.cram")) as inf:
            self.assertFalse(inf.is_bam)
            self.assertTrue(inf.is_cram)
            self.assertEqual(len(list(inf.fetch('chr1'))), 1)
            self.assertEqual(len(list(inf.fetch('chr2'))), 3)


##################################################
#
# Random access iterator tests
#
##################################################
class TestIteratorRowBAM(unittest.TestCase):

    filename = os.path.join(BAM_DATADIR, "ex2.bam")
    mode = "rb"
    reference_filename = None

    def setUp(self):
        self.samfile = pysam.AlignmentFile(
            self.filename,
            self.mode,
            reference_filename=self.reference_filename,
        )

    def tearDown(self):
        self.samfile.close()

    def checkRange(self, rnge):
        '''compare results from iterator with those from samtools.'''
        ps = list(self.samfile.fetch(region=rnge))
        sa = force_str(
            pysam.samtools.view(
                self.filename,
                rnge,
                raw=True)).splitlines(True)
        self.assertEqual(
            len(ps), len(sa),
            "unequal number of results for range %s: %i != %i" %
            (rnge, len(ps), len(sa)))
        # check if the same reads are returned and in the same order
        for line, (a, b) in enumerate(list(zip(ps, sa))):
            d = b.split("\t")
            self.assertEqual(
                a.query_name, d[0],
                "line %i: read id mismatch: %s != %s" %
                (line, a.reference_id, d[0]))
            self.assertEqual(
                a.reference_start,
                int(d[3]) - 1,
                "line %i: read position mismatch: %s != %s, \n%s\n%s\n" %
                (line, a.reference_start, int(d[3]) - 1,
                 str(a), str(d)))
            qual = d[10]
            self.assertEqual(
                pysam.qualities_to_qualitystring(a.query_qualities),
                qual,
                "line %i: quality mismatch: %s != %s, \n%s\n%s\n" %
                (line, pysam.qualities_to_qualitystring(a.query_qualities), qual,
                 str(a), str(d)))

    def testIteratePerContig(self):
        '''check random access per contig'''
        for contig in self.samfile.references:
            self.checkRange(contig)

    def testIterateRanges(self):
        '''check random access per range'''
        for contig, length in zip(self.samfile.references,
                                  self.samfile.lengths):
            for start in range(1, length, 90):
                # this includes empty ranges
                self.checkRange("%s:%i-%i" %
                                (contig, start, start + 90))


class TestIteratorRowAllBAM(unittest.TestCase):

    filename = os.path.join(BAM_DATADIR, "ex2.bam")
    mode = "rb"

    def setUp(self):
        self.samfile = pysam.AlignmentFile(
            self.filename,
            self.mode)

    def testIterate(self):
        '''compare results from iterator with those from samtools.'''
        ps = list(self.samfile.fetch())
        sa = pysam.samtools.view(self.filename,
                                 raw=True).splitlines()
        self.assertEqual(
            len(ps), len(sa),
            "unequal number of results: %i != %i" %
            (len(ps), len(sa)))
        # check if the same reads are returned
        for line, pair in enumerate(list(zip(ps, sa))):
            data = force_str(pair[1]).split("\t")
            self.assertEqual(
                pair[0].query_name,
                data[0],
                "read id mismatch in line %i: %s != %s" %
                (line, pair[0].reference_id, data[0]))

    def tearDown(self):
        self.samfile.close()


class TestIteratorColumnBAM(unittest.TestCase):

    '''test iterator column against contents of ex4.bam.'''

    # note that samfile contains 1-based coordinates
    # 1D means deletion with respect to reference sequence
    #
    mCoverages = {'chr1': [0] * 20 + [1] * 36 + [0] * (100 - 20 - 35),
                  'chr2': [0] * 20 + [1] * 35 + [0] * (100 - 20 - 35),
                  }

    def setUp(self):
        self.samfile = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex4.bam"),
                                           "rb")

    def checkRange(self, contig, start=None, end=None, truncate=False):
        '''compare results from iterator with those from samtools.'''
        # check if the same reads are returned and in the same order
        for column in self.samfile.pileup(
                contig, start, end, truncate=truncate):
            if truncate:
                self.assertGreaterEqual(column.reference_pos, start)
                self.assertLess(column.reference_pos, end)
            thiscov = len(column.pileups)
            refcov = self.mCoverages[
                self.samfile.getrname(column.reference_id)][column.reference_pos]
            self.assertEqual(thiscov, refcov,
                             "wrong coverage at pos %s:%i %i should be %i" % (
                                 self.samfile.getrname(column.reference_id),
                                 column.reference_pos, thiscov, refcov))

    def testIterateAll(self):
        '''check random access per contig'''
        self.checkRange(None)

    def testIteratePerContig(self):
        '''check random access per contig'''
        for contig in self.samfile.references:
            self.checkRange(contig)

    def testIterateRanges(self):
        '''check random access per range'''
        for contig, length in zip(
                self.samfile.references, self.samfile.lengths):
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
        for contig, length in zip(self.samfile.references,
                                  self.samfile.lengths):
            for start in range(1, length, 90):
                # this includes empty ranges
                self.checkRange(contig, start, start + 90, truncate=True)

    def tearDown(self):
        self.samfile.close()


class TestIteratorRowCRAM(TestIteratorRowBAM):
    filename = os.path.join(BAM_DATADIR, "ex2.cram")
    mode = "rc"


class TestIteratorRowCRAMWithReferenceFilename(TestIteratorRowCRAM):
    reference_filename = os.path.join(BAM_DATADIR, "ex1.fa")


##########################################################
##########################################################
##########################################################
# needs to be implemented
# class TestAlignedSegmentFromSamWithoutHeader(TestAlignedSegmentFromBam):
#
#     def setUp(self):
#         self.samfile=pysam.AlignmentFile( "ex7.sam","r" )
#         self.reads=list(self.samfile.fetch())


class TestIteratorColumn2(unittest.TestCase):

    '''test iterator column against contents of ex1.bam.'''

    def setUp(self):
        self.samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex1.bam"),
            "rb")

    def testStart(self):
        # print self.samfile.fetch().next().reference_start
        # print self.samfile.pileup().next().reference_start
        pass

    def testTruncate(self):
        '''see issue 107.'''
        # note that ranges in regions start from 1
        p = self.samfile.pileup(region='chr1:170:172', truncate=True)
        columns = [x.reference_pos for x in p]
        self.assertEqual(len(columns), 3)
        self.assertEqual(columns, [169, 170, 171])

        p = self.samfile.pileup('chr1', 169, 172, truncate=True)
        columns = [x.reference_pos for x in p]

        self.assertEqual(len(columns), 3)
        self.assertEqual(columns, [169, 170, 171])

    def testAccessOnClosedIterator(self):
        '''see issue 131

        Accessing pileup data after iterator has closed.
        '''
        pcolumn = self.samfile.pileup('chr1', 170, 180).__next__()
        self.assertRaises(ValueError, getattr, pcolumn, "pileups")

    def testStr(self):
        '''test if PileupRead can be printed.'''
        iter = self.samfile.pileup('chr1', 170, 180)
        pcolumn = iter.__next__()
        s = str(pcolumn)
        self.assertEqual(len(s.split("\n")), 2)


class TestFloatTagBug(unittest.TestCase):

    '''see issue 71'''

    def testFloatTagBug(self):
        '''a float tag before another exposed a parsing bug in bam_aux_get.

        Fixed in 0.1.19
        '''
        samfile = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "tag_bug.bam"))
        read = next(samfile.fetch(until_eof=True))
        self.assertTrue(('XC', 1) in read.tags)
        self.assertEqual(read.opt('XC'), 1)


class TestLargeFieldBug(unittest.TestCase):

    '''see issue 100'''

    def testLargeFileBug(self):
        '''when creating a read with a large entry in the tag field
        causes an error:
            NotImplementedError: tags field too large
        '''
        samfile = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "issue100.bam"))
        read = next(samfile.fetch(until_eof=True))
        new_read = pysam.AlignedSegment()
        new_read.tags = read.tags
        self.assertEqual(new_read.tags, read.tags)


class TestTagParsing(unittest.TestCase):

    '''tests checking the accuracy of tag setting and retrieval.'''

    def makeRead(self):
        a = pysam.AlignedSegment()
        a.query_name = "read_12345"
        a.reference_id = 0
        a.query_sequence = "ACGT" * 3
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 1
        a.mapping_quality = 20
        a.cigartuples = ((0, 10), (2, 1), (0, 25))
        a.next_reference_id = 0
        a.next_reference_start = 200
        a.template_length = 0
        a.query_qualities = pysam.qualitystring_to_array("1234") * 3
        # todo: create tags
        return a

    def testNegativeIntegers(self):
        x = -2
        aligned_read = self.makeRead()
        aligned_read.tags = [("XD", int(x))]
        self.assertEqual(aligned_read.opt('XD'), x)
        # print (aligned_read.tags)

    def testNegativeIntegers2(self):
        x = -2
        r = self.makeRead()
        r.tags = [("XD", x)]
        outfile = pysam.AlignmentFile(
            "tests/test.bam",
            "wb",
            referencenames=("chr1",),
            referencelengths = (1000,))
        outfile.write(r)
        outfile.close()
        infile = pysam.AlignmentFile("tests/test.bam")
        r = next(infile)
        self.assertEqual(r.tags, [("XD", x)])
        infile.close()
        os.unlink("tests/test.bam")

    def testCigarString(self):
        r = self.makeRead()
        self.assertEqual(r.cigarstring, "10M1D25M")
        r.cigarstring = "20M10D20M"
        self.assertEqual(r.cigartuples, [(0, 20), (2, 10), (0, 20)])
        # unsetting cigar string
        r.cigarstring = None
        self.assertEqual(r.cigarstring, None)

    def testCigar(self):
        r = self.makeRead()
        self.assertEqual(r.cigartuples, [(0, 10), (2, 1), (0, 25)])
        # unsetting cigar string
        r.cigartuples = None
        self.assertEqual(r.cigartuples, None)

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

    def testArrayTags(self):

        r = self.makeRead()

        def c(r, l):
            r.tags = [('ZM', l)]
            self.assertEqual(list(r.opt("ZM")), list(l))

        # signed integers
        c(r, (-1, 1))
        c(r, (-1, 100))
        c(r, (-1, 200))
        c(r, (-1, 1000))
        c(r, (-1, 30000))
        c(r, (-1, 50000))
        c(r, (1, -1))
        c(r, (1, -100))
        c(r, (1, -200))
        c(r, (1, -1000))
        c(r, (1, -30000))
        c(r, (1, -50000))

        # unsigned integers
        c(r, (1, 100))
        c(r, (1, 1000))
        c(r, (1, 10000))
        c(r, (1, 100000))

        # floats
        c(r, (1.0, 100.0))


class TestClipping(unittest.TestCase):

    def testClipping(self):

        self.samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "softclip.bam"),
            "rb")

        for read in self.samfile:

            if read.query_name == "r001":
                self.assertEqual(read.query_sequence, 'AAAAGATAAGGATA')
                self.assertEqual(read.query_alignment_sequence, 'AGATAAGGATA')
                self.assertEqual(pysam.qualities_to_qualitystring(read.query_qualities),
                                 None)
                self.assertEqual(
                    pysam.qualities_to_qualitystring(read.query_alignment_qualities),
                    None)

            elif read.query_name == "r002":

                self.assertEqual(read.query_sequence, 'GCCTAAGCTAA')
                self.assertEqual(read.query_alignment_sequence, 'AGCTAA')
                self.assertEqual(
                    pysam.qualities_to_qualitystring(read.query_qualities),
                    '01234567890')
                self.assertEqual(
                    pysam.qualities_to_qualitystring(read.query_alignment_qualities),
                    '567890')

            elif read.query_name == "r003":

                self.assertEqual(read.query_sequence, 'GCCTAAGCTAA')
                self.assertEqual(read.query_alignment_sequence, 'GCCTAA')
                self.assertEqual(
                    pysam.qualities_to_qualitystring(read.query_qualities),
                    '01234567890')
                self.assertEqual(
                    pysam.qualities_to_qualitystring(read.query_alignment_qualities),
                    '012345')

            elif read.query_name == "r004":

                self.assertEqual(read.query_sequence, 'TAGGC')
                self.assertEqual(read.query_alignment_sequence, 'TAGGC')
                self.assertEqual(
                    pysam.qualities_to_qualitystring(read.query_qualities),
                    '01234')
                self.assertEqual(
                    pysam.qualities_to_qualitystring(read.query_alignment_qualities),
                    '01234')


class TestHeaderSAM(unittest.TestCase):

    """testing header manipulation"""

    header = {'SQ': [{'LN': 1575, 'SN': 'chr1', 'AH': 'chr1:5000000-5010000'},
                     {'LN': 1584, 'SN': 'chr2', 'AH': '*'}],
              'RG': [{'LB': 'SC_1', 'ID': 'L1', 'SM': 'NA12891',
                      'PU': 'SC_1_10', "CN": "name:with:colon"},
                     {'LB': 'SC_2', 'ID': 'L2', 'SM': 'NA12891',
                      'PU': 'SC_2_12', "CN": "name:with:colon"}],
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
        self.samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex3.sam"),
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


class TestHeaderBAM(TestHeaderSAM):

    def setUp(self):
        self.samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex3.bam"),
            "rb")


class TestHeaderCRAM(TestHeaderSAM):

    def setUp(self):
        self.samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex3.cram"),
            "rc")

    def compareHeaders(self, a, b):
        '''compare two headers a and b.'''
        def _strip(dd):
            for x in dd:
                for y in ("M5", "UR"):
                    if y in x:
                        del x[y]

        for ak, av in a.items():
            _strip(av)
            self.assertTrue(ak in b, "key '%s' not in '%s' " % (ak, b))
            _strip(b[ak])

            self.assertEqual(av, b[ak])


class TestHeaderFromRefs(unittest.TestCase):

    '''see issue 144

    reference names need to be converted to string for python 3
    '''

    # def testHeader( self ):
    #     refs = ['chr1', 'chr2']
    #     tmpfile = "tmp_%i" % id(self)
    #     s = pysam.AlignmentFile(tmpfile, 'wb',
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

        f = pysam.AlignmentFile(self.bamfile, "rb")
        data = f.header.copy()
        self.assertTrue(data)


class TestHeaderWriteRead(unittest.TestCase):
    header = {'SQ': [{'LN': 1575, 'SN': 'chr1'},
                     {'LN': 1584, 'SN': 'chr2'}],
              'RG': [{'LB': 'SC_1', 'ID': 'L1', 'SM': 'NA12891',
                      'PU': 'SC_1_10', "CN": "name:with:colon"},
                     {'LB': 'SC_2', 'ID': 'L2', 'SM': 'NA12891',
                      'PU': 'SC_2_12', "CN": "name:with:colon"}],
              'PG': [{'ID': 'P1', 'VN': '1.0', 'CL': 'tool'},
                     {'ID': 'P2', 'VN': '1.1', 'CL': 'tool with in option -R a\tb',
                      'PP': 'P1'}],
              'HD': {'VN': '1.0'},
              'CO': ['this is a comment', 'this is another comment'],
    }

    def compare_headers(self, a, b):
        '''compare two headers a and b.

        Ignore M5 and UR field as they are set application specific.
        '''
        for ak, av in a.items():
            self.assertTrue(ak in b, "key '%s' not in '%s' " % (ak, b))
            self.assertEqual(
                len(av), len(b[ak]),
                "unequal number of entries for key {}: {} vs {}"
                .format(ak, av, b[ak]))

            for row_a, row_b in zip(av, b[ak]):
                if isinstance(row_b, dict):
                    for x in ["M5", "UR"]:
                        try:
                            del row_b[x]
                        except KeyError:
                            pass
                self.assertEqual(row_a, row_b)

    def check_read_write(self, flag_write, header):

        fn = get_temp_filename()
        with pysam.AlignmentFile(
                fn,
                flag_write,
                header=header,
                reference_filename=os.path.join(BAM_DATADIR, "ex1.fa")) as outf:
            a = pysam.AlignedSegment()
            a.query_name = "abc"
            outf.write(a)

        with pysam.AlignmentFile(fn) as inf:
            read_header = inf.header

        os.unlink(fn)
        self.compare_headers(header, read_header)

    def test_SAM(self):
        self.check_read_write("wh", self.header)

    def test_BAM(self):
        self.check_read_write("wb", self.header)

    def test_CRAM(self):
        header = copy.copy(self.header)
        # for CRAM, \t needs to be quoted:
        header['PG'][1]['CL'] = re.sub(r"\t", r"\\\\t", header['PG'][1]['CL'])
        self.check_read_write("wc", header)


class TestUnmappedReads(unittest.TestCase):

    # TODO
    # def testSAM(self):
    #     samfile = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex5.sam"),
    #                             "r")
    #     self.assertEqual(len(list(samfile.fetch(until_eof=True))), 2)
    #     samfile.close()

    def testBAM(self):
        samfile = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex5.bam"),
                                      "rb")
        self.assertEqual(len(list(samfile.fetch(until_eof=True))), 2)
        samfile.close()


class TestPileupObjects(unittest.TestCase):

    def setUp(self):
        self.samfile = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex1.bam"),
                                           "rb")

    def testPileupColumn(self):
        for pcolumn1 in self.samfile.pileup(region="chr1:105"):
            if pcolumn1.reference_pos == 104:
                self.assertEqual(
                    pcolumn1.reference_id, 0,
                    "chromosome/target id mismatch in position 1: %s != %s" %
                    (pcolumn1.reference_id, 0))
                self.assertEqual(
                    pcolumn1.reference_pos, 105 - 1,
                    "position mismatch in position 1: %s != %s" %
                    (pcolumn1.reference_pos, 105 - 1))
                self.assertEqual(
                    pcolumn1.nsegments, 2,
                    "# reads mismatch in position 1: %s != %s" %
                    (pcolumn1.nsegments, 2))
        for pcolumn2 in self.samfile.pileup(region="chr2:1480"):
            if pcolumn2.reference_pos == 1479:
                self.assertEqual(
                    pcolumn2.reference_id, 1,
                    "chromosome/target id mismatch in position 1: %s != %s" %
                    (pcolumn2.reference_id, 1))
                self.assertEqual(
                    pcolumn2.reference_pos, 1480 - 1,
                    "position mismatch in position 1: %s != %s" %
                    (pcolumn2.reference_pos, 1480 - 1))
                self.assertEqual(
                    pcolumn2.nsegments, 12,
                    "# reads mismatch in position 1: %s != %s" %
                    (pcolumn2.nsegments, 12))

    def testPileupRead(self):
        for pcolumn1 in self.samfile.pileup(region="chr1:105"):
            if pcolumn1.reference_pos == 104:
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
        with pysam.AlignmentFile(os.path.join(BAM_DATADIR, 'ex1.bam'),
                                 'rb') as samfile:
            samfile.fetch()
        self.assertEqual(samfile.closed, True)


class TestExceptions(unittest.TestCase):

    def setUp(self):
        self.samfile = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex1.bam"),
                                           "rb")

    def testMissingFile(self):

        self.assertRaises(
            IOError, pysam.AlignmentFile, "exdoesntexist.bam", "rb")
        self.assertRaises(
            IOError, pysam.AlignmentFile, "exdoesntexist.sam", "r")
        self.assertRaises(
            IOError, pysam.AlignmentFile, "exdoesntexist.bam", "r")
        self.assertRaises(
            IOError, pysam.AlignmentFile, "exdoesntexist.sam", "rb")

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
                          pysam.AlignmentFile,
                          os.path.join(BAM_DATADIR, 'ex1.sam'),
                          'rb')

    def testOpenBamAsSam(self):
        # test fails, needs to be implemented.
        # sam.fetch() fails on reading, not on opening
        #self.assertRaises(ValueError, pysam.AlignmentFile,
        #                  os.path.join(BAM_DATADIR, 'ex1.bam'),
        #                  'r')
        pass

    def testOpenFastaAsSam(self):
        # test fails, needs to be implemented.
        # sam.fetch() fails on reading, not on opening
        # self.assertRaises( ValueError, pysam.AlignmentFile, 'ex1.fa', 'r' )
        pass

    def testOpenFastaAsBam(self):
        self.assertRaises(ValueError,
                          pysam.AlignmentFile,
                          os.path.join(BAM_DATADIR, 'ex1.fa'),
                          'rb')


class TestDeNovoConstruction(unittest.TestCase):

    '''check BAM/SAM file construction using ex6.sam

    (note these are +1 coordinates):

    read_28833_29006_6945	99	chr1	33	20	10M1D25M	=	200	167	AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG	<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<	NM:i:1	RG:Z:L1
    read_28701_28881_323b	147	chr2	88	30	35M	=	500	412	ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA	<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<	MF:i:18	RG:Z:L2
    '''

    header = {'HD': {'VN': '1.0'},
              'SQ': [{'LN': 1575, 'SN': 'chr1'},
                     {'LN': 1584, 'SN': 'chr2'}], }

    bamfile = os.path.join(BAM_DATADIR, "ex6.bam")
    samfile = os.path.join(BAM_DATADIR, "ex6.sam")

    def setUp(self):

        a = pysam.AlignedSegment()
        a.query_name = "read_28833_29006_6945"
        a.query_sequence = "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
        a.flag = 99
        a.reference_id = 0
        a.reference_start = 32
        a.mapping_quality = 20
        a.cigartuples = ((0, 10), (2, 1), (0, 25))
        a.next_reference_id = 0
        a.next_reference_start = 199
        a.template_length = 167
        a.query_qualities = pysam.qualitystring_to_array(
            "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<")
        a.tags = (("NM", 1),
                  ("RG", "L1"))

        b = pysam.AlignedSegment()
        b.query_name = "read_28701_28881_323b"
        b.query_sequence = "ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA"
        b.flag = 147
        b.reference_id = 1
        b.reference_start = 87
        b.mapping_quality = 30
        b.cigartuples = ((0, 35), )
        b.next_reference_id = 1
        b.next_reference_start = 499
        b.template_length = 412
        b.query_qualities = pysam.qualitystring_to_array(
            "<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<")
        b.tags = (("MF", 18),
                  ("RG", "L2"))

        self.reads = (a, b)

    # TODO
    # def testSAMWholeFile(self):

    #     tmpfilename = "tmp_%i.sam" % id(self)

    #     outfile = pysam.AlignmentFile(tmpfilename,
    #                             "wh",
    #                             header=self.header)

    #     for x in self.reads:
    #         outfile.write(x)
    #     outfile.close()
    #     self.assertTrue(checkBinaryEqual(tmpfilename, self.samfile),
    #                     "mismatch when construction SAM file, see %s %s" % (tmpfilename, self.samfile))

    #     os.unlink(tmpfilename)

    def test_pass_if_reads_binary_equal(self):
        '''check if individual reads are binary equal.'''
        infile = pysam.AlignmentFile(self.bamfile, "rb")

        others = list(infile)
        for denovo, other in zip(others, self.reads):
            checkFieldEqual(self, other, denovo)
            self.assertEqual(other.compare(denovo), 0)

    # TODO
    # def testSAMPerRead(self):
    #     '''check if individual reads are binary equal.'''
    #     infile = pysam.AlignmentFile(self.samfile, "r")

    #     others = list(infile)
    #     for denovo, other in zip(others, self.reads):
    #         checkFieldEqual(self, other, denovo)
    #         self.assertEqual(other.compare(denovo), 0)

    def testBAMWholeFile(self):

        tmpfilename = "tmp_%i.bam" % id(self)

        outfile = pysam.AlignmentFile(tmpfilename, "wb", header=self.header)

        for x in self.reads:
            outfile.write(x)
        outfile.close()

        self.assertTrue(
            checkBinaryEqual(tmpfilename, self.bamfile),
            "mismatch when construction BAM file, see %s %s" %
            (tmpfilename, self.bamfile))

        os.unlink(tmpfilename)


class TestDeNovoConstructionUserTags(TestDeNovoConstruction):

    '''test de novo construction with a header that contains lower-case tags.'''

    header = {'HD': {'VN': '1.0'},
              'SQ': [{'LN': 1575, 'SN': 'chr1'},
                     {'LN': 1584, 'SN': 'chr2'}],
              'x1': {'A': 2, 'B': 5},
              'x3': {'A': 6, 'B': 5},
              'x2': {'A': 4, 'B': 5}}

    bamfile = os.path.join(BAM_DATADIR, "example_user_header.bam")
    samfile = os.path.join(BAM_DATADIR, "example_user_header.sam")


class TestEmptyHeader(unittest.TestCase):

    '''see issue 84.'''

    def testEmptyHeader(self):
        s = pysam.AlignmentFile(os.path.join(BAM_DATADIR,
                                             'example_empty_header.bam'))
        self.assertEqual(s.header, {'SQ': [{'LN': 1000, 'SN': 'chr1'}]})


class TestHeaderWithProgramOptions(unittest.TestCase):

    '''see issue 39.'''

    def testHeader(self):
        s = pysam.AlignmentFile(os.path.join(BAM_DATADIR,
                                             'rg_with_tab.bam'))
        self.assertEqual(
            s.header,
            {'SQ': [{'LN': 1575, 'SN': 'chr1'},
                    {'LN': 1584, 'SN': 'chr2'}],
             'PG': [{'PN': 'bwa',
                     'ID': 'bwa',
                     'VN': '0.7.9a-r786',
                     'CL': 'bwa mem -p -t 8 -M -R '
                     '@RG\tID:None\tSM:None\t/mnt/data/hg19.fa\t'
                     '/mnt/analysis/default-0.fastq'}]})


class TestTruncatedBAM(unittest.TestCase):

    '''see pull request 50.'''

    def testTruncatedBam2(self):
        self.assertRaises(IOError,
                          pysam.AlignmentFile,
                          os.path.join(BAM_DATADIR, 'ex2_truncated.bam'))

    def testTruncatedBam2(self):
        s = pysam.AlignmentFile(os.path.join(BAM_DATADIR, 'ex2_truncated.bam'),
                                ignore_truncation=True)
        iterall = lambda x: len([a for a in x])
        self.assertRaises(IOError, iterall, s)


COMPARE_BTAG = [100, 1, 91, 0, 7, 101, 0, 201, 96, 204,
                0, 0, 87, 109, 0, 7, 97, 112, 1, 12, 78,
                197, 0, 7, 100, 95, 101, 202, 0, 6, 0, 1,
                186, 0, 84, 0, 244, 0, 0, 324, 0, 107, 195,
                101, 113, 0, 102, 0, 104, 3, 0, 101, 1, 0,
                212, 6, 0, 0, 1, 0, 74, 1, 11, 0, 196, 2,
                197, 103, 0, 108, 98, 2, 7, 0, 1, 2, 194,
                0, 180, 0, 108, 0, 203, 104, 16, 5, 205,
                0, 0, 0, 1, 1, 100, 98, 0, 0, 204, 6, 0,
                79, 0, 0, 101, 7, 109, 90, 265, 1, 27, 10,
                109, 102, 9, 0, 292, 0, 110, 0, 0, 102,
                112, 0, 0, 84, 100, 103, 2, 81, 126, 0, 2,
                90, 0, 15, 96, 15, 1, 0, 2, 0, 107, 92, 0,
                0, 101, 3, 98, 15, 102, 13, 116, 116, 90, 93,
                198, 0, 0, 0, 199, 92, 26, 495, 100, 5, 0,
                100, 5, 209, 0, 92, 107, 90, 0, 0, 0, 0, 109,
                194, 7, 94, 200, 0, 40, 197, 0, 11, 0, 0, 112,
                110, 6, 4, 200, 28, 0, 196, 0, 203, 1, 129,
                0, 0, 1, 0, 94, 0, 1, 0, 107, 5, 201, 3, 3, 100,
                0, 121, 0, 7, 0, 1, 105, 306, 3, 86, 8, 183, 0,
                12, 163, 17, 83, 22, 0, 0, 1, 8, 109, 103, 0, 0,
                295, 0, 200, 16, 172, 3, 16, 182, 3, 11, 0, 0,
                223, 111, 103, 0, 5, 225, 0, 95]

class TestBTagSam(unittest.TestCase):

    '''see issue 81.'''

    compare = [COMPARE_BTAG,
               [-100, 200, -300, -400],
               [-100, 12],
               [12, 15],
               [-1.0, 5.0, 2.5]]

    filename = os.path.join(BAM_DATADIR, 'example_btag.sam')

    read0 = [('RG', 'QW85I'),
             ('PG', 'tmap'),
             ('MD', '140'),
             ('NM', 0),
             ('AS', 140),
             ('FZ', array.array('H', COMPARE_BTAG)),
             ('XA', 'map2-1'),
             ('XS', 53),
             ('XT', 38),
             ('XF', 1),
             ('XE', 0)]

    def testReadTags(self):

        s = pysam.AlignmentFile(self.filename)
        for x, read in enumerate(s):
            tags = read.tags
            if x == 0:
                self.assertEqual(tags, self.read0)

            fz = list(dict(tags)["FZ"])
            self.assertEqual(fz, self.compare[x])
            self.assertEqual(list(read.opt("FZ")), self.compare[x])
            self.assertEqual(tags, read.get_tags())
            for tag, value in tags:
                self.assertEqual(value, read.get_tag(tag))

    def testReadWriteTags(self):

        s = pysam.AlignmentFile(self.filename)
        for read in s:
            before = read.tags
            read.tags = before
            self.assertEqual(read.tags, before)

            read.set_tags(before)
            self.assertEqual(read.tags, before)

            for tag, value in before:
                read.set_tag(tag, value)
                self.assertEqual(value, read.get_tag(tag))


class TestBTagBam(TestBTagSam):
    filename = os.path.join(BAM_DATADIR, 'example_btag.bam')


class TestDoubleFetchBAM(unittest.TestCase):
    '''check if two iterators on the same bamfile are independent.'''

    filename = os.path.join(BAM_DATADIR, 'ex1.bam')
    mode = "rb"

    def testDoubleFetch(self):

        with pysam.AlignmentFile(self.filename, self.mode) as samfile1:
            for a, b in zip(samfile1.fetch(multiple_iterators=True),
                            samfile1.fetch(multiple_iterators=True)):
                self.assertEqual(a.compare(b), 0)

    def testDoubleFetchWithRegion(self):

        with pysam.AlignmentFile(self.filename, self.mode) as samfile1:
            contig, start, stop = 'chr1', 200, 3000000
            # just making sure the test has something to catch
            self.assertTrue(len(list(samfile1.fetch(contig, start, stop))) > 0)

            # see Issue #293
            # The following fails for CRAM files, but works for BAM
            # files when the first is multiple_iterators=False:
            for a, b in zip(samfile1.fetch(contig, start, stop,
                                           multiple_iterators=True),
                            samfile1.fetch(contig, start, stop,
                                           multiple_iterators=True)):
                self.assertEqual(a.compare(b), 0)

    def testDoubleFetchUntilEOF(self):

        with pysam.AlignmentFile(self.filename, self.mode) as samfile1:

            for a, b in zip(samfile1.fetch(until_eof=True),
                            samfile1.fetch(until_eof=True,
                                           multiple_iterators=True)):
                self.assertEqual(a.compare(b), 0)


class TestDoubleFetchCRAM(TestDoubleFetchBAM):
    filename = os.path.join(BAM_DATADIR, 'ex2.cram')
    mode = "rc"


class TestDoubleFetchCRAMWithReference(TestDoubleFetchBAM):
    filename = os.path.join(BAM_DATADIR, 'ex2.cram')
    mode = "rc"
    reference_filename = os.path.join(BAM_DATADIR, 'ex1.fa')


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

        result = pysam.samtools.view(self.url, self.region)
        self.assertEqual(len(result), 36)

    def testFTPFetch(self):
        return
        if not checkURL(self.url):
            return

        samfile = pysam.AlignmentFile(self.url, "rb")
        result = list(samfile.fetch(region=self.region))
        self.assertEqual(len(result), 36)


class TestRemoteFileHTTP(unittest.TestCase):

    url = "http://genserv.anat.ox.ac.uk/downloads/pysam/test/ex1.bam"
    region = "chr1:1-1000"
    local = os.path.join(BAM_DATADIR, "ex1.bam")

    def testView(self):
        if not checkURL(self.url):
            return

        samfile_local = pysam.AlignmentFile(self.local, "rb")
        ref = list(samfile_local.fetch(region=self.region))

        result = pysam.samtools.view(self.url, self.region)
        self.assertEqual(len(result.splitlines()), len(ref))

    def testFetch(self):
        if not checkURL(self.url):
            return

        with pysam.AlignmentFile(self.url, "rb") as samfile:
            result = list(samfile.fetch(region=self.region))

        with pysam.AlignmentFile(self.local, "rb") as samfile_local:
            ref = list(samfile_local.fetch(region=self.region))

        self.assertEqual(len(ref), len(result))
        for x, y in zip(result, ref):
            self.assertEqual(x.compare(y), 0)

    def testFetchAll(self):
        if not checkURL(self.url):
            return

        with pysam.AlignmentFile(self.url, "rb") as samfile:
            result = list(samfile.fetch())

        with pysam.AlignmentFile(self.local, "rb") as samfile_local:
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
        samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex10.sam"),
            "r")
        self.check(samfile)

    def testBAM(self):
        samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex10.bam"),
            "rb")
        self.check(samfile)


class TestPileup(unittest.TestCase):

    '''test pileup functionality.'''

    samfilename = os.path.join(BAM_DATADIR, "ex1.bam")
    fastafilename = os.path.join(BAM_DATADIR, "ex1.fa")

    def setUp(self):

        self.samfile = pysam.AlignmentFile(self.samfilename)
        self.fastafile = pysam.Fastafile(self.fastafilename)

    def tearDown(self):
        self.samfile.close()
        self.fastafile.close()

    def checkEqual(self, references, iterator):

        for x, column in enumerate(iterator):
            v = references[x][:-1].split("\t")
            self.assertEqual(
                len(v), 6,
                "expected 6 values, got {}".format(v))
            (contig, pos, reference_base,
             read_bases, read_qualities, alignment_mapping_qualities) \
                = v
            self.assertEqual(int(pos) - 1, column.reference_pos)

    def testSamtoolsStepper(self):
        refs = force_str(
            pysam.samtools.mpileup(
                "-f", self.fastafilename,
                self.samfilename)).splitlines(True)
        iterator = self.samfile.pileup(
            stepper="samtools",
            fastafile=self.fastafile)
        self.checkEqual(refs, iterator)

    def testAllStepper(self):
        refs = force_str(
            pysam.samtools.mpileup(
                "-f", self.fastafilename,
                "-A", "-B",
                self.samfilename)).splitlines(True)

        iterator = self.samfile.pileup(
            stepper="all",
            fastafile=self.fastafile)
        self.checkEqual(refs, iterator)


class TestCountCoverage(unittest.TestCase):

    samfilename = os.path.join(BAM_DATADIR, "ex1.bam")
    fastafilename = os.path.join(BAM_DATADIR, "ex1.fa")

    def setUp(self):

        self.samfile = pysam.AlignmentFile(self.samfilename)
        self.fastafile = pysam.Fastafile(self.fastafilename)

        samfile = pysam.AlignmentFile(
            "tests/test_count_coverage_read_all.bam", 'wb',
            template=self.samfile)
        for ii, read in enumerate(self.samfile.fetch()):
            # if ii % 2 == 0: # setting BFUNMAP makes no sense...
                #read.flag = read.flag | 0x4
            if ii % 3 == 0:
                read.flag = read.flag | 0x100
            if ii % 5 == 0:
                read.flag = read.flag | 0x200
            if ii % 7 == 0:
                read.flag = read.flag | 0x400
            samfile.write(read)
        samfile.close()
        pysam.samtools.index("tests/test_count_coverage_read_all.bam")

    def tearDown(self):
        self.samfile.close()
        self.fastafile.close()
        os.unlink("tests/test_count_coverage_read_all.bam")
        os.unlink("tests/test_count_coverage_read_all.bam.bai")

    def count_coverage_python(self, bam, chrom, start, stop,
                              read_callback,
                              quality_threshold=15):
        l = stop - start
        count_a = array.array('L', [0] * l)
        count_c = array.array('L', [0] * l)
        count_g = array.array('L', [0] * l)
        count_t = array.array('L', [0] * l)
        for p in bam.pileup(chrom, start, stop, truncate=True,
                            stepper='nofilter'):
            rpos = p.reference_pos - start
            for read in p.pileups:
                if not read.is_del and not read.is_refskip and \
                   read_callback(read.alignment):
                    try:
                        if read.alignment.query_qualities[read.query_position] \
                           >= quality_threshold:
                            letter = read.alignment.query[read.query_position]
                            if letter == 'A':
                                count_a[rpos] += 1
                            elif letter == 'C':
                                count_c[rpos] += 1
                            elif letter == 'G':
                                count_g[rpos] += 1
                            elif letter == 'T':
                                count_t[rpos] += 1
                    except IndexError:
                        pass
        return count_a, count_c, count_g, count_t

    def test_count_coverage(self):
        chrom = 'chr1'
        start = 0
        stop = 2000
        manual_counts = self.count_coverage_python(
            self.samfile, chrom, start, stop,
            lambda read: True,
            quality_threshold=0)
        fast_counts = self.samfile.count_coverage(
            chrom, start, stop,
            read_callback=lambda read: True,
            quality_threshold=0)

        self.assertEqual(list(fast_counts[0]), list(manual_counts[0]))
        self.assertEqual(list(fast_counts[1]), list(manual_counts[1]))
        self.assertEqual(list(fast_counts[2]), list(manual_counts[2]))
        self.assertEqual(list(fast_counts[3]), list(manual_counts[3]))

    def test_count_coverage_quality_filter(self):
        chrom = 'chr1'
        start = 0
        stop = 2000
        manual_counts = self.count_coverage_python(
            self.samfile, chrom, start, stop,
            lambda read: True,
            quality_threshold=0)
        fast_counts = self.samfile.count_coverage(
            chrom, start, stop,
            read_callback=lambda read: True,
            quality_threshold=15)
        # we filtered harder, should be less
        for i in range(4):
            for r in range(start, stop):
                self.assertTrue(fast_counts[i][r] <= manual_counts[i][r])

    def test_count_coverage_read_callback(self):
        chrom = 'chr1'
        start = 0
        stop = 2000
        manual_counts = self.count_coverage_python(
            self.samfile, chrom, start, stop,
            lambda read: read.flag & 0x10,
            quality_threshold=0)
        fast_counts = self.samfile.count_coverage(
            chrom, start, stop,
            read_callback=lambda read: True,
            quality_threshold=0)
        for i in range(4):
            for r in range(start, stop):
                self.assertTrue(fast_counts[i][r] >= manual_counts[i][r])
        fast_counts = self.samfile.count_coverage(
            chrom, start, stop,
            read_callback=lambda read: read.flag & 0x10,
            quality_threshold=0)

        self.assertEqual(fast_counts[0], manual_counts[0])
        self.assertEqual(fast_counts[1], manual_counts[1])
        self.assertEqual(fast_counts[2], manual_counts[2])
        self.assertEqual(fast_counts[3], manual_counts[3])

    def test_count_coverage_read_all(self):

        chrom = 'chr1'
        start = 0
        stop = 2000

        def filter(read):
            return not (read.flag & (0x4 | 0x100 | 0x200 | 0x400))

        with pysam.AlignmentFile("tests/test_count_coverage_read_all.bam") as samfile:

            fast_counts = samfile.count_coverage(
                chrom, start, stop,
                read_callback='all',
                #read_callback = lambda read: ~(read.flag & (0x4 | 0x100 | 0x200 | 0x400)),
                quality_threshold=0)
            manual_counts = samfile.count_coverage(
                chrom, start, stop,
                read_callback=lambda read: not(
                    read.flag & (0x4 | 0x100 | 0x200 | 0x400)),
                quality_threshold=0)

        self.assertEqual(fast_counts[0], manual_counts[0])
        self.assertEqual(fast_counts[1], manual_counts[1])
        self.assertEqual(fast_counts[2], manual_counts[2])
        self.assertEqual(fast_counts[3], manual_counts[3])

    def test_count_coverage_nofilter(self):
        samfile = pysam.AlignmentFile(
            "tests/test_count_coverage_nofilter.bam", 'wb', template=self.samfile)
        for ii, read in enumerate(self.samfile.fetch()):
            # if ii % 2 == 0: # setting BFUNMAP makes no sense...
                #read.flag = read.flag | 0x4
            if ii % 3 == 0:
                read.flag = read.flag | 0x100
            if ii % 5 == 0:
                read.flag = read.flag | 0x200
            if ii % 7 == 0:
                read.flag = read.flag | 0x400
            samfile.write(read)
        samfile.close()
        pysam.samtools.index("tests/test_count_coverage_nofilter.bam")
        chr = 'chr1'
        start = 0
        stop = 2000

        with pysam.AlignmentFile("tests/test_count_coverage_nofilter.bam") as samfile:

            fast_counts = samfile.count_coverage(chr, start, stop,
                                                 read_callback='nofilter',
                                                 quality_threshold=0)

            manual_counts = self.count_coverage_python(samfile, chr, start, stop,
                                                       read_callback=lambda x: True,
                                                       quality_threshold=0)

        os.unlink("tests/test_count_coverage_nofilter.bam")
        os.unlink("tests/test_count_coverage_nofilter.bam.bai")
        self.assertEqual(fast_counts[0], manual_counts[0])
        self.assertEqual(fast_counts[1], manual_counts[1])
        self.assertEqual(fast_counts[2], manual_counts[2])
        self.assertEqual(fast_counts[3], manual_counts[3])


class TestPileupQueryPosition(unittest.TestCase):

    filename = "test_query_position.bam"

    def testPileup(self):
        last = {}
        with pysam.AlignmentFile(os.path.join(BAM_DATADIR, self.filename)) as inf:
            for col in inf.pileup():
                for r in col.pileups:
                    # print r.alignment.query_name
                    # print r.query_position, r.query_position_or_next, r.is_del
                    if r.is_del:
                        self.assertEqual(r.query_position, None)
                        self.assertEqual(r.query_position_or_next,
                                         last[r.alignment.query_name] + 1)
                    else:
                        self.assertNotEqual(r.query_position, None)
                        last[r.alignment.query_name] = r.query_position


class TestFindIntrons(unittest.TestCase):
    samfilename = os.path.join(BAM_DATADIR, "ex_spliced.bam")

    def setUp(self):
        self.samfile = pysam.AlignmentFile(self.samfilename)

    def tearDown(self):
        self.samfile.close()

    def test_total(self):
        all_read_counts = self.samfile.count()
        splice_sites = self.samfile.find_introns(self.samfile.fetch())
        self.assertEqual(sum(splice_sites.values()), all_read_counts -1)  # there is a single unspliced read in there

    def test_first(self):
        reads = list(self.samfile.fetch())[:10]
        splice_sites = self.samfile.find_introns(reads)
        starts = [14792+38 - 1]
        stops = [14792+38 + 140 - 1]
        self.assertEqual(len(splice_sites), 1)
        self.assertTrue((starts[0], stops[0]) in splice_sites)
        self.assertEqual(splice_sites[(starts[0], stops[0])], 9) # first one is the unspliced read

    def test_all(self):
        reads = list(self.samfile.fetch())
        splice_sites = self.samfile.find_introns(reads)
        should = collections.Counter({
            (14829, 14969): 33,
            (15038, 15795): 24,
            (15947, 16606): 3,
            (16765, 16857): 9,
            (16765, 16875): 1,
            (17055, 17232): 19,
            (17055, 17605): 3,
            (17055, 17914): 1,
            (17368, 17605): 7,
            })
        self.assertEqual(should,  splice_sites)


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

        with pysam.AlignmentFile(bamfile, 'rb') as bam:
            cols = bam.pileup()
        self.assertTrue(True)

    def testFail1(self):
        self.check(os.path.join(BAM_DATADIR, "ex9_fail.bam"),
                   False)
        self.check(os.path.join(BAM_DATADIR, "ex9_fail.bam"),
                   True)

    def testNoFail1(self):
        self.check(os.path.join(BAM_DATADIR, "ex9_nofail.bam"),
                   False)
        self.check(os.path.join(BAM_DATADIR, "ex9_nofail.bam"),
                   True)

    def testNoFail2(self):
        self.check(os.path.join(BAM_DATADIR, "ex9_nofail.bam"),
                   True)
        self.check(os.path.join(BAM_DATADIR, "ex9_nofail.bam"),
                   True)

# TODOS
# 1. finish testing all properties within pileup objects
# 2. check exceptions and bad input problems (missing files, optional fields that aren't present, etc...)
# 3. check: presence of sequence


class TestAlignmentFileUtilityFunctions(unittest.TestCase):

    def testCount(self):

        with pysam.AlignmentFile(
                os.path.join(BAM_DATADIR, "ex1.bam"),
                "rb") as samfile:

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

        with open(os.path.join(BAM_DATADIR, "ex1.sam"), "rb") as inf:
            readnames = [x.split(b"\t")[0] for x in inf.readlines()]
        if sys.version_info[0] >= 3:
            readnames = [name.decode('ascii') for name in readnames]

        counts = collections.defaultdict(int)
        for x in readnames:
            counts[x] += 1

        with pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex1.bam"),
                                 "rb") as samfile:

            for read in samfile.fetch():
                if not read.is_paired:
                    self.assertRaises(ValueError, samfile.mate, read)
                elif read.mate_is_unmapped:
                    self.assertRaises(ValueError, samfile.mate, read)
                else:
                    if counts[read.query_name] == 1:
                        self.assertRaises(ValueError, samfile.mate, read)
                    else:
                        mate = samfile.mate(read)
                        self.assertEqual(read.query_name, mate.query_name)
                        self.assertEqual(read.is_read1, mate.is_read2)
                        self.assertEqual(read.is_read2, mate.is_read1)
                        self.assertEqual(
                            read.reference_start, mate.next_reference_start)
                        self.assertEqual(
                            read.next_reference_start, mate.reference_start)

    def testIndexStats(self):
        '''test if total number of mapped/unmapped reads is correct.'''

        with pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex1.bam"),
                                 "rb") as samfile:
            self.assertEqual(samfile.mapped, 3235)
            self.assertEqual(samfile.unmapped, 35)
            self.assertEqual(samfile.nocoordinate, 0)


class TestMappedUnmapped(unittest.TestCase):
    filename = "test_mapped_unmapped.bam"

    def test_counts_of_mapped_and_unmapped_are_correct(self):

        with pysam.AlignmentFile(os.path.join(BAM_DATADIR,
                                              self.filename)) as inf:
            unmapped_flag = 0
            unmapped_nopos = 0
            mapped_flag = 0
            for x in inf.fetch(until_eof=True):
                if x.is_unmapped:
                    if x.reference_id < 0:
                        unmapped_nopos += 1
                    else:
                        unmapped_flag += 1
                else:
                    mapped_flag += 1

            self.assertEqual(inf.mapped, mapped_flag)
            self.assertEqual(inf.unmapped, unmapped_flag + unmapped_nopos)

            inf.reset()
            self.assertEqual(inf.count(),
                             inf.mapped + unmapped_flag)

            inf.reset()
            self.assertEqual(inf.count(until_eof=True),
                             inf.mapped + unmapped_flag + unmapped_nopos)

            inf.reset()
            self.assertEqual(inf.count(read_callback="all"),
                             inf.mapped)

            inf.reset()
            self.assertEqual(inf.count(until_eof=True, read_callback="all"),
                             inf.mapped)

    def test_counts_of_mapped_and_unmapped_are_correct_per_chromosome(self):

        with pysam.AlignmentFile(os.path.join(BAM_DATADIR,
                                              self.filename)) as inf:

            counts = inf.get_index_statistics()

            counts_contigs = [x.contig for x in counts]
            self.assertEqual(sorted(counts_contigs),
                             sorted(inf.references))
            
            for contig in inf.references:
                unmapped_flag = 0
                unmapped_nopos = 0
                mapped_flag = 0
                for x in inf.fetch(contig=contig):
                    if x.is_unmapped:
                        unmapped_flag += 1
                    else:
                        mapped_flag += 1

                cc = [c for c in counts if c.contig == contig][0]
                self.assertEqual(cc.mapped, mapped_flag)
                self.assertEqual(cc.unmapped, unmapped_flag)
                self.assertEqual(cc.total, mapped_flag + unmapped_flag)

class TestSamtoolsProxy(unittest.TestCase):

    '''tests for sanity checking access to samtools functions.'''

    def testIndex(self):
        self.assertRaises(IOError, pysam.samtools.index, "missing_file")

    def testView(self):
        # note that view still echos "open: No such file or directory"
        self.assertRaises(pysam.SamtoolsError,
                          pysam.samtools.view,
                          "missing_file")

    def testSort(self):
        self.assertRaises(pysam.SamtoolsError,
                          pysam.samtools.sort,
                          "missing_file")


class TestAlignmentFileIndex(unittest.TestCase):

    def testIndex(self):
        samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "ex1.bam"),
            "rb")
        index = pysam.IndexedReads(samfile)
        index.build()
        reads = collections.defaultdict(int)

        for read in samfile:
            reads[read.query_name] += 1

        for qname, counts in reads.items():
            found = list(index.find(qname))
            self.assertEqual(len(found), counts)
            for x in found:
                self.assertEqual(x.query_name, qname)


class TestExplicitIndex(unittest.TestCase):

    def testExplicitIndexBAM(self):
        with pysam.AlignmentFile(
                os.path.join(BAM_DATADIR, "explicit_index.bam"),
                "rb",
                filepath_index=os.path.join(BAM_DATADIR, 'ex1.bam.bai')) as samfile:
            samfile.fetch("chr1")

    def testExplicitIndexCRAM(self):
        with pysam.AlignmentFile(
                os.path.join(BAM_DATADIR, "explicit_index.cram"),
                "rc",
                filepath_index=os.path.join(BAM_DATADIR, 'ex1.cram.crai')) as samfile:
            samfile.fetch("chr1")

    def testRemoteExplicitIndexBAM(self):
        if not checkURL(
                "http://genserv.anat.ox.ac.uk/downloads/pysam/test/noindex.bam"):
            return

        with pysam.AlignmentFile(
                "http://genserv.anat.ox.ac.uk/downloads/pysam/test/noindex.bam",
                "rb",
                filepath_index=os.path.join(BAM_DATADIR, 'ex1.bam.bai')) as samfile:
            samfile.fetch("chr1")


class TestVerbosity(unittest.TestCase):

    '''test if setting/getting of verbosity works.'''

    def testVerbosity(self):
        self.assertEqual(pysam.get_verbosity(), 3)
        old = pysam.set_verbosity(0)
        self.assertEqual(pysam.get_verbosity(), 0)
        new = pysam.set_verbosity(old)
        self.assertEqual(new, 0)
        self.assertEqual(pysam.get_verbosity(), 3)


class TestSanityCheckingBAM(unittest.TestCase):

    mode = "wb"

    def check_write(self, read):

        fn = "tmp_test_sanity_check.bam"
        names = ["chr1"]
        lengths = [10000]
        with pysam.AlignmentFile(
                fn,
                self.mode,
                reference_names=names,
                reference_lengths=lengths) as outf:
            outf.write(read)

        if os.path.exists(fn):
            os.unlink(fn)

    def test_empty_read_gives_value_error(self):
        read = pysam.AlignedSegment()
        self.check_write(read)

# SAM writing fails, as query length is 0
# class TestSanityCheckingSAM(TestSanityCheckingSAM):
#     mode = "w"


if __name__ == "__main__":
    # build data files
    print ("building data files")
    subprocess.call("make -C %s" % BAM_DATADIR, shell=True)
    print ("starting tests")
    unittest.main()
    print ("completed tests")
