#!/usr/bin/env python
'''unit testing code for pysam.

Execute in the :file:`tests` directory as it requires the Makefile
and data files located there.
'''

import unittest
import pytest
import os
import shutil
import sys
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
from TestUtils import checkBinaryEqual, checkGZBinaryEqual, check_url, \
    check_samtools_view_equal, checkFieldEqual, force_str, \
    get_temp_filename, make_data_files, BAM_DATADIR


def setUpModule():
    make_data_files(BAM_DATADIR)


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
        self.assertEqual(self.reads[0].query_name,
                         "read_28833_29006_6945")
        self.assertEqual(self.reads[1].query_name,
                         "read_28701_28881_323b")

    def testARflag(self):
        self.assertEqual(self.reads[0].flag, 99)
        self.assertEqual(self.reads[1].flag, 147)

    def testARrname(self):
        self.assertEqual(self.reads[0].reference_id, 0)
        self.assertEqual(self.reads[1].reference_id, 1)

    def testARpos(self):
        self.assertEqual(self.reads[0].reference_start, 33 - 1)
        self.assertEqual(self.reads[1].reference_start, 88 - 1)

    def testARmapq(self):
        self.assertEqual(self.reads[0].mapping_quality, 20)
        self.assertEqual(self.reads[1].mapping_quality, 30)

    def testARcigar(self):
        self.assertEqual(self.reads[0].cigartuples, [(0, 10), (2, 1), (0, 25)])
        self.assertEqual(self.reads[1].cigartuples, [(0, 35)])

    def testARcigarstring(self):
        self.assertEqual(self.reads[0].cigarstring, '10M1D25M')
        self.assertEqual(self.reads[1].cigarstring, '35M')

    def testARmrnm(self):
        self.assertEqual(self.reads[0].next_reference_id, 0)
        self.assertEqual(self.reads[1].next_reference_id, 1)
        self.assertEqual(self.reads[0].next_reference_id, 0)
        self.assertEqual(self.reads[1].next_reference_id, 1)

    def testARmpos(self):
        self.assertEqual(self.reads[0].next_reference_start, 200 - 1)
        self.assertEqual(self.reads[1].next_reference_start, 500 - 1)
        self.assertEqual(self.reads[0].next_reference_start, 200 - 1)
        self.assertEqual(self.reads[1].next_reference_start, 500 - 1)

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
        self.assertEqual(
            self.reads[0].query_sequence,
            "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG")
        self.assertEqual(
            self.reads[1].query_sequence,
            "ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA")
        self.assertEqual(
            self.reads[3].query_sequence,
            "AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG")

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
            pysam.qualities_to_qualitystring(
                self.reads[0].query_alignment_qualities),
            "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<",
            "qquality string mismatch in read 1: %s != %s" %
            (pysam.qualities_to_qualitystring(self.reads[0].query_alignment_qualities),
             "<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<"))
        self.assertEqual(
            pysam.qualities_to_qualitystring(
                self.reads[1].query_alignment_qualities),
            "<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<",
            "qquality string mismatch in read 2: %s != %s" %
            (pysam.qualities_to_qualitystring(self.reads[1].query_alignment_qualities),
             "<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<"))
        self.assertEqual(
            pysam.qualities_to_qualitystring(
                self.reads[3].query_alignment_qualities),
            "<<<<<<<<<<<<<<<<<:<9/,&,22",
            "qquality string mismatch in read 3: %s != %s" %
            (pysam.qualities_to_qualitystring(self.reads[3].query_alignment_qualities),
             "<<<<<<<<<<<<<<<<<:<9/,&,22"))

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
                         [('NM', 22), ('RG', 'L1'),
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
                       "rb", "wb",
                       checkf=checkGZBinaryEqual)

    def testCRAM2CRAM(self):
        # in some systems different reference sequence paths might be
        # embedded in the CRAM files which will result in different headers
        # see #542
        self.checkEcho("ex2.cram",
                       "ex2.cram",
                       "tmp_ex2.cram",
                       "rc", "wc",
                       sequence_filename=os.path.join(BAM_DATADIR, "ex1.fa"),
                       checkf=partial(
                           check_samtools_view_equal,
                           without_header=True))

    def testSAM2BAM(self):
        self.checkEcho("ex2.sam",
                       "ex2.bam",
                       "tmp_ex2.bam",
                       "r", "wb",
                       checkf=checkGZBinaryEqual)

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
                       "rb", "wb",
                       checkf=checkGZBinaryEqual)

        f = open(os.path.join(BAM_DATADIR, "ex1.bam"))
        samfile = pysam.AlignmentFile(f, "rb")
        self.assertFalse(f.closed)
        samfile.close()
        # python file needs to be closed separately
        self.assertFalse(f.closed)

    def test_accessing_attributes_in_closed_file_raises_errors(self):
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

        self.assertEqual(samfile.header, None)
        # write on closed file
        self.assertEqual(0, samfile.write(None))

    def test_header_available_after_closing_file(self):

        def load_bam():
            with pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex1.bam"), "rb") as inf:
                header = inf.header
            return header

        header = load_bam()
        self.assertTrue(header)
        self.assertEqual(header.nreferences, 2)
        self.assertEqual(header.references, ("chr1", "chr2"))

    def test_reference_name_available_after_closing_file(self):
        """read tids can be mapped to references after AlignmentFile has been closed.
        
        see issue #517"""
        
        def load_bam():
            with pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex1.bam"), "rb") as inf:
                read = next(iter(inf))
            return read

        read = load_bam()
        self.assertEqual(read.reference_name, "chr1")
        
    # TODO
    # def testReadingFromSamFileWithoutHeader(self):
    #     '''read from samfile without header.
    #     '''
    #     samfile = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex7.sam"),
    #                             check_header=False,
    #                             check_sq=False)
    #     self.assertRaises(NotImplementedError, samfile.__iter__)

    def testReadingFromFileWithoutIndex(self):
        '''read from bam file without index.'''

        dest = get_temp_filename("tmp_ex2.bam")
        shutil.copyfile(os.path.join(BAM_DATADIR, "ex2.bam"),
                        dest)
        samfile = pysam.AlignmentFile(dest,
                                      "rb")
        self.assertRaises(ValueError, samfile.fetch)
        self.assertEqual(
            len(list(samfile.fetch(until_eof=True))),
            3270)
        os.unlink(dest)

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
            self.assertEqual(len(list(samfile.fetch('chr1', start=1000, end=2000))),
                             len(list(samfile.fetch(tid=0, start=1000, end=2000))))

    def test_write_bam_to_unknown_path_fails(self):
        '''see issue 116'''
        input_filename = os.path.join(BAM_DATADIR, "ex1.bam")
        with pysam.AlignmentFile(input_filename) as inf:
            self.assertRaises(IOError,
                              pysam.AlignmentFile,
                              "missing_directory/new_file.bam",
                              "wb",
                              template=inf)
        

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


class TestIteratorRowCRAM(TestIteratorRowBAM):
    filename = os.path.join(BAM_DATADIR, "ex2.cram")
    mode = "rc"


class TestIteratorRowCRAMWithReferenceFilename(TestIteratorRowCRAM):
    reference_filename = os.path.join(BAM_DATADIR, "ex1.fa")


# needs to be implemented
# class TestAlignedSegmentFromSamWithoutHeader(TestAlignedSegmentFromBam):
#
#     def setUp(self):
#         self.samfile=pysam.AlignmentFile( "ex7.sam","r" )
#         self.reads=list(self.samfile.fetch())


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
        samfile = pysam.AlignmentFile(
            os.path.join(BAM_DATADIR, "issue100.bam"))
        read = next(samfile.fetch(until_eof=True))
        new_read = pysam.AlignedSegment()
        new_read.tags = read.tags
        self.assertEqual(new_read.tags, read.tags)


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
                    pysam.qualities_to_qualitystring(
                        read.query_alignment_qualities),
                    None)

            elif read.query_name == "r002":

                self.assertEqual(read.query_sequence, 'GCCTAAGCTAA')
                self.assertEqual(read.query_alignment_sequence, 'AGCTAA')
                self.assertEqual(
                    pysam.qualities_to_qualitystring(read.query_qualities),
                    '01234567890')
                self.assertEqual(
                    pysam.qualities_to_qualitystring(
                        read.query_alignment_qualities),
                    '567890')

            elif read.query_name == "r003":

                self.assertEqual(read.query_sequence, 'GCCTAAGCTAA')
                self.assertEqual(read.query_alignment_sequence, 'GCCTAA')
                self.assertEqual(
                    pysam.qualities_to_qualitystring(read.query_qualities),
                    '01234567890')
                self.assertEqual(
                    pysam.qualities_to_qualitystring(
                        read.query_alignment_qualities),
                    '012345')

            elif read.query_name == "r004":

                self.assertEqual(read.query_sequence, 'TAGGC')
                self.assertEqual(read.query_alignment_sequence, 'TAGGC')
                self.assertEqual(
                    pysam.qualities_to_qualitystring(read.query_qualities),
                    '01234')
                self.assertEqual(
                    pysam.qualities_to_qualitystring(
                        read.query_alignment_qualities),
                    '01234')


class TestUnmappedReadsRetrieval(unittest.TestCase):

    def test_fetch_from_sam_with_until_eof_reads_unmapped_reads(self):
        with pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex5.sam"),
                                 "rb") as samfile:
            self.assertEqual(len(list(samfile.fetch(until_eof=True))), 2)

    def test_fetch_from_bam_with_until_eof_reads_unmapped_reads(self):
        with pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex5.bam"),
                                 "rb") as samfile:
            self.assertEqual(len(list(samfile.fetch(until_eof=True))), 2)

    def test_fetch_with_asterisk_only_returns_unmapped_reads(self):
        with pysam.AlignmentFile(os.path.join(BAM_DATADIR, "test_mapped_unmapped.bam"),
                                 "rb") as samfile:
            self.assertEqual(len(list(samfile.fetch(region="*"))), 4)

    def test_fetch_with_asterisk_only_returns_unmapped_reads_by_contig(self):
        with pysam.AlignmentFile(os.path.join(BAM_DATADIR, "test_mapped_unmapped.bam"),
                                 "rb") as samfile:
            self.assertEqual(len(list(samfile.fetch(contig="*"))), 4)


class TestContextManager(unittest.TestCase):

    def testManager(self):
        with pysam.AlignmentFile(os.path.join(BAM_DATADIR, 'ex1.bam'),
                                 'rb') as samfile:
            samfile.fetch()
        self.assertEqual(samfile.closed, True)


class TestMultiThread(unittest.TestCase):

    def testSingleThreadEqualsMultithread(self):
        input_bam = os.path.join(BAM_DATADIR, 'ex1.bam')
        single_thread_out = get_temp_filename("tmp_single.bam")
        multi_thread_out = get_temp_filename("tmp_multi.bam")
        with pysam.AlignmentFile(input_bam,
                                 'rb') as samfile:
            reads = [r for r in samfile]
            with pysam.AlignmentFile(single_thread_out,
                                     mode='wb',
                                     template=samfile,
                                     threads=1) as single_out:
                [single_out.write(r) for r in reads]
            with pysam.AlignmentFile(multi_thread_out,
                                     mode='wb',
                                     template=samfile,
                                     threads=2) as multi_out:
                [single_out.write(r) for r in reads]
        with pysam.AlignmentFile(input_bam) as inp, \
            pysam.AlignmentFile(single_thread_out) as single, \
            pysam.AlignmentFile(multi_thread_out) as multi:
            for r1, r2, r3 in zip(inp, single, multi):
                assert r1.to_string == r2.to_string == r3.to_string

    def TestNoMultiThreadingWithIgnoreTruncation(self):
        self.assertRaises(
            ValueError, pysam.AlignmentFile(os.path.join(BAM_DATADIR, 'ex1.bam'),
                                            threads=2,
                                            ignore_truncation=True)
        )


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
        # self.assertRaises(ValueError, pysam.AlignmentFile,
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


class TestRegionParsing(unittest.TestCase):

    def test_dash_in_chr(self):
        with pysam.AlignmentFile(
                os.path.join(BAM_DATADIR, "example_dash_in_chr.bam")) as inf:
            self.assertEqual(len(list(inf.fetch(contig="chr-1"))), 1)
            self.assertEqual(len(list(inf.fetch(contig="chr2"))), 1)
            self.assertEqual(len(list(inf.fetch(region="chr-1"))), 1)
            self.assertEqual(len(list(inf.fetch(region="chr2"))), 1)


class TestDeNovoConstruction(unittest.TestCase):

    '''check BAM/SAM file construction using ex6.sam

    (note these are +1 coordinates):

    read_28833_29006_6945	99	chr1	33	20	10M1D25M	=	200	167	AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG	<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<	NM:i:1	RG:Z:L1
    read_28701_28881_323b	147	chr2	88	30	35M	=	500	412	ACCTATATCTTGGCCTTGGCCGATGCGGCCTTGCA	<<<<<;<<<<7;:<<<6;<<<<<<<<<<<<7<<<<	MF:i:18	RG:Z:L2
    '''  # noqa

    header = {'HD': {'VN': '1.0'},
              'SQ': [{'LN': 1575, 'SN': 'chr1'},
                     {'LN': 1584, 'SN': 'chr2'}], }

    bamfile = os.path.join(BAM_DATADIR, "ex6.bam")
    samfile = os.path.join(BAM_DATADIR, "ex6.sam")

    def setUp(self):

        header = pysam.AlignmentHeader.from_dict(self.header)

        a = pysam.AlignedSegment(header)
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

        b = pysam.AlignedSegment(header)
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

        references = list(infile)
        for denovo, reference in zip(references, self.reads):
            checkFieldEqual(self, reference, denovo)
            self.assertEqual(reference.compare(denovo), 0)

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
            checkGZBinaryEqual(tmpfilename, self.bamfile),
            "mismatch when construction BAM file, see %s %s" %
            (tmpfilename, self.bamfile))

        os.unlink(tmpfilename)


class TestEmptyHeader(unittest.TestCase):
    '''see issue 84.'''

    def testEmptyHeader(self):
        s = pysam.AlignmentFile(os.path.join(BAM_DATADIR,
                                             'example_empty_header.bam'))
        self.assertEqual(s.header.to_dict(), {'SQ': [{'LN': 1000, 'SN': 'chr1'}]})

    def test_bam_without_seq_in_header(self):
        s = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "0example_no_seq_in_header.bam"))
        self.assertTrue("SQ" in s.header.to_dict())
        self.assertTrue("@SQ" in str(s.header))

    def test_bam_without_seq_with_null_bytes_in_header(self):
        s = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "0example_no_seq_in_header_null_bytes.bam"))
        self.assertTrue("SQ" in s.header.to_dict())
        self.assertTrue("@SQ" in str(s.header))


class TestMismatchingHeader(unittest.TestCase):
    '''see issue 716.'''

    def testMismatchingHeader(self):
        # Note: no chr2
        header = {
            'SQ': [{'SN': 'chr1', 'LN': 1575}],
            'PG': [{'ID': 'bwa', 'PN': 'bwa', 'VN': '0.7.15', 'CL': 'bwa mem xx -'}],
        }

        dest = get_temp_filename("tmp_ex3.bam")
        with pysam.AlignmentFile(os.path.join(BAM_DATADIR, 'ex3.bam')) as inf:
            with pysam.AlignmentFile(dest, mode="wb", header=header) as outf:
                for read in inf:
                    if read.reference_name == "chr1":
                        outf.write(read)
                    else:
                        self.assertRaises(ValueError,
                                          outf.write,
                                          read)
        os.unlink(dest)


class TestHeaderWithProgramOptions(unittest.TestCase):

    '''see issue 39.'''

    def testHeader(self):
        s = pysam.AlignmentFile(os.path.join(BAM_DATADIR,
                                             'rg_with_tab.bam'))
        self.assertEqual(
            s.header.to_dict(),
            {'SQ': [{'LN': 1575, 'SN': 'chr1'},
                    {'LN': 1584, 'SN': 'chr2'}],
             'PG': [{'PN': 'bwa',
                     'ID': 'bwa',
                     'VN': '0.7.9a-r786',
                     'CL': 'bwa mem -p -t 8 -M -R '
                     '@RG ID:None SM:None /mnt/data/hg19.fa '
                     '/mnt/analysis/default-0.fastq'}]})


class TestTruncatedBAM(unittest.TestCase):

    '''see pull request 50.'''

    def testTruncatedBam2(self):
        self.assertRaises(IOError,
                          pysam.AlignmentFile,
                          os.path.join(BAM_DATADIR, 'ex2_truncated.bam'))

    @pytest.mark.filterwarnings('ignore:no BGZF EOF marker')
    def testTruncatedBamIterator(self):
        s = pysam.AlignmentFile(os.path.join(BAM_DATADIR, 'ex2_truncated.bam'),
                                ignore_truncation=True)

        def iterall(x):
            return len([a for a in x])
        self.assertRaises(IOError, iterall, s)

        # Ignore closing errors, as s is now in an error state
        try:
            s.close()
        except IOError:
            pass


class TestCorruptBAM(unittest.TestCase):
    """See pull request 1035."""

    def testCorruptBamIterator(self):
        s = pysam.AlignmentFile(os.path.join(BAM_DATADIR, "ex2_corrupt.bam"))

        def iterall(x):
            return len([a for a in x])

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

    def testDoubleFetchWithRegionTrueTrue(self):

        with pysam.AlignmentFile(self.filename, self.mode) as samfile1:
            contig, start, stop = 'chr2', 200, 3000000
            # just making sure the test has something to catch
            self.assertTrue(len(list(samfile1.fetch(contig, start, stop))) > 0)

            for a, b in zip(samfile1.fetch(contig, start, stop,
                                           multiple_iterators=True),
                            samfile1.fetch(contig, start, stop,
                                           multiple_iterators=True)):
                self.assertEqual(a.compare(b), 0)

    def testDoubleFetchWithRegionFalseTrue(self):
        with pysam.AlignmentFile(self.filename, self.mode) as samfile1:
            contig, start, stop = 'chr2', 200, 3000000
            # just making sure the test has something to catch
            self.assertTrue(len(list(samfile1.fetch(contig, start, stop))) > 0)
        
            for a, b in zip(samfile1.fetch(contig, start, stop,
                                           multiple_iterators=False),
                            samfile1.fetch(contig, start, stop,
                                           multiple_iterators=True)):
                self.assertEqual(a.compare(b), 0)

    def testDoubleFetchWithRegionTrueFalse(self):
        with pysam.AlignmentFile(self.filename, self.mode) as samfile1:
            contig, start, stop = 'chr2', 200, 3000000
            # just making sure the test has something to catch
            self.assertTrue(len(list(samfile1.fetch(contig, start, stop))) > 0)
        
            for a, b in zip(samfile1.fetch(contig, start, stop,
                                           multiple_iterators=True),
                            samfile1.fetch(contig, start, stop,
                                           multiple_iterators=False)):
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
        if not check_url(self.url):
            return

        result = pysam.samtools.view(self.url, self.region)
        self.assertEqual(len(result), 36)

    def testFTPFetch(self):
        return
        if not check_url(self.url):
            return

        samfile = pysam.AlignmentFile(self.url, "rb")
        result = list(samfile.fetch(region=self.region))
        self.assertEqual(len(result), 36)


class TestRemoteFileHTTP(unittest.TestCase):

    url = "http://genserv.anat.ox.ac.uk/downloads/pysam/test/ex1.bam"
    region = "chr1:1-1000"
    local = os.path.join(BAM_DATADIR, "ex1.bam")

    def testView(self):
        if not check_url(self.url):
            return

        samfile_local = pysam.AlignmentFile(self.local, "rb")
        ref = list(samfile_local.fetch(region=self.region))

        result = pysam.samtools.view(self.url, self.region)
        self.assertEqual(len(result.splitlines()), len(ref))

    def testFetch(self):
        if not check_url(self.url):
            return

        with pysam.AlignmentFile(self.url, "rb") as samfile:
            result = list(samfile.fetch(region=self.region))

        with pysam.AlignmentFile(self.local, "rb") as samfile_local:
            ref = list(samfile_local.fetch(region=self.region))

        self.assertEqual(len(ref), len(result))
        for x, y in zip(result, ref):
            self.assertEqual(x.compare(y), 0)

    def testFetchAll(self):
        if not check_url(self.url):
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


class TestCountCoverage(unittest.TestCase):

    samfilename = os.path.join(BAM_DATADIR, "ex1.bam")
    fastafilename = os.path.join(BAM_DATADIR, "ex1.fa")

    def setUp(self):

        self.fastafile = pysam.Fastafile(self.fastafilename)
        self.tmpfilename = get_temp_filename(".bam")

        with pysam.AlignmentFile(self.samfilename) as inf:
            with pysam.AlignmentFile(
                    self.tmpfilename,
                    'wb',
                    template=inf) as outf:
                for ii, read in enumerate(inf.fetch()):
                    # if ii % 2 == 0: # setting BFUNMAP makes no sense...
                        #read.flag = read.flag | 0x4
                    if ii % 3 == 0:
                        read.flag = read.flag | 0x100
                    if ii % 5 == 0:
                        read.flag = read.flag | 0x200
                    if ii % 7 == 0:
                        read.flag = read.flag | 0x400
                    outf.write(read)
        pysam.samtools.index(self.tmpfilename)

    def tearDown(self):
        self.fastafile.close()
        os.unlink(self.tmpfilename)
        os.unlink(self.tmpfilename + ".bai")

    def count_coverage_python(self,
                              bam, chrom, start, stop,
                              read_callback,
                              quality_threshold=15):
        stop = min(stop, bam.get_reference_length(chrom))
        l = stop - start
        count_a = array.array('L', [0] * l)
        count_c = array.array('L', [0] * l)
        count_g = array.array('L', [0] * l)
        count_t = array.array('L', [0] * l)
        for p in bam.pileup(chrom, start, stop,
                            truncate=True,
                            stepper='nofilter',
                            min_base_quality=quality_threshold,
                            ignore_overlaps=False):
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

    def test_count_coverage_with_coordinates_works(self):

        with pysam.AlignmentFile(self.samfilename) as inf:
            c = inf.count_coverage("chr1")
            self.assertEqual(len(c[0]), inf.get_reference_length("chr1"))
            self.assertEqual(len(c[0]), 1575)

            c = inf.count_coverage("chr1", 100)
            self.assertEqual(len(c[0]), inf.get_reference_length("chr1") - 100)

            c = inf.count_coverage("chr1", 100, 200)
            self.assertEqual(len(c[0]), 200 - 100)
            
            c = inf.count_coverage("chr1", None, 200)
            self.assertEqual(len(c[0]), 200)

            c = inf.count_coverage("chr1", None, inf.get_reference_length("chr1") + 10000)
            self.assertEqual(len(c[0]), inf.get_reference_length("chr1"))

    def test_count_coverage_without_coordinates_fails(self):
        with pysam.AlignmentFile(self.samfilename) as inf:
            self.assertRaises(TypeError, inf.count_coverage)

    def test_count_coverage_wrong_coordinates_fails(self):
        with pysam.AlignmentFile(self.samfilename) as inf:
            self.assertRaises(ValueError, inf.count_coverage, "chr1", 200, 100)
            self.assertRaises(KeyError, inf.count_coverage, "chrUnknown", 100, 200)
            
    def test_counting_the_same_region_works(self):

        with pysam.AlignmentFile(self.samfilename) as inf:
            c1 = inf.count_coverage("chr1")
            c2 = inf.count_coverage("chr1")
            self.assertEqual(c1, c2)
            
    def test_count_coverage_counts_as_expected(self):
        chrom = 'chr1'
        start = 0
        stop = 1000

        with pysam.AlignmentFile(self.samfilename) as inf:
            manual_counts = self.count_coverage_python(
                inf, chrom, start, stop,
                lambda read: True,
                quality_threshold=0)

            fast_counts = inf.count_coverage(
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
        stop = 1000
        with pysam.AlignmentFile(self.samfilename) as inf:
            manual_counts = self.count_coverage_python(
                inf, chrom, start, stop,
                lambda read: True,
                quality_threshold=0)
            fast_counts = inf.count_coverage(
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
        stop = 1000
        with pysam.AlignmentFile(self.samfilename) as inf:
            manual_counts = self.count_coverage_python(
                inf, chrom, start, stop,
                lambda read: read.flag & 0x10,
                quality_threshold=0)
            fast_counts = inf.count_coverage(
                chrom, start, stop,
                read_callback=lambda read: True,
                quality_threshold=0)
            
            for i in range(4):
                for r in range(start, stop):
                    self.assertTrue(fast_counts[i][r] >= manual_counts[i][r])
                    
            fast_counts = inf.count_coverage(
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
        stop = 1000

        def filter(read):
            return not (read.flag & (0x4 | 0x100 | 0x200 | 0x400))

        with pysam.AlignmentFile(self.tmpfilename) as samfile:
            fast_counts = samfile.count_coverage(
                chrom, start, stop,
                read_callback='all',
                # read_callback = lambda read: ~(read.flag & (0x4 | 0x100 | 0x200 | 0x400)),
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

        with pysam.AlignmentFile(self.samfilename) as inf:
            with pysam.AlignmentFile(
                    self.tmpfilename, 'wb', template=inf) as outf:
        
                for ii, read in enumerate(inf.fetch()):
                    # if ii % 2 == 0: # setting BFUNMAP makes no sense...
                    # read.flag = read.flag | 0x4
                    if ii % 3 == 0:
                        read.flag = read.flag | 0x100
                    if ii % 5 == 0:
                        read.flag = read.flag | 0x200
                    if ii % 7 == 0:
                        read.flag = read.flag | 0x400
                    outf.write(read)
        
        pysam.samtools.index(self.tmpfilename)
        chr = 'chr1'
        start = 0
        stop = 1000

        with pysam.AlignmentFile(self.tmpfilename) as inf:
            fast_counts = inf.count_coverage(chr, start, stop,
                                             read_callback='nofilter',
                                             quality_threshold=0)

            manual_counts = self.count_coverage_python(inf, chr, start, stop,
                                                       read_callback=lambda x: True,
                                                       quality_threshold=0)

        self.assertEqual(fast_counts[0], manual_counts[0])
        self.assertEqual(fast_counts[1], manual_counts[1])
        self.assertEqual(fast_counts[2], manual_counts[2])
        self.assertEqual(fast_counts[3], manual_counts[3])
        

class TestFindIntrons(unittest.TestCase):
    samfilename = os.path.join(BAM_DATADIR, "ex_spliced.bam")

    def setUp(self):
        self.samfile = pysam.AlignmentFile(self.samfilename)

    def tearDown(self):
        self.samfile.close()

    def test_total(self):
        all_read_counts = self.samfile.count()
        splice_sites = self.samfile.find_introns(self.samfile.fetch())
        # there is a single unspliced read and a single unmapped read in there
        self.assertEqual(sum(splice_sites.values()), all_read_counts - 2)

    def test_first(self):
        reads = list(self.samfile.fetch())[:10]
        splice_sites = self.samfile.find_introns(reads)
        starts = [14792 + 38 - 1]
        stops = [14792 + 38 + 140 - 1]
        self.assertEqual(len(splice_sites), 1)
        self.assertTrue((starts[0], stops[0]) in splice_sites)
        # first one is the unspliced read
        self.assertEqual(splice_sites[(starts[0], stops[0])], 9)

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
        self.assertEqual(should, splice_sites)


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

    def tearDown(self):
        if os.path.exists('log.txt'):
            os.unlink('log.txt')

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
        if not check_url(
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


class TestLargeCigar(unittest.TestCase):

    def setUp(self):
        self.read_length = 70000
        self.header = pysam.AlignmentHeader.from_references(
            ["chr1", "chr2"],
            [self.read_length * 2, self.read_length * 2])

    def build_read(self):
        '''build an example read.'''

        a = pysam.AlignedSegment(self.header)
        l = self.read_length
        a.query_name = "read_12345"
        a.query_sequence = "A" * (l + 1)
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 20
        a.mapping_quality = 20
        a.cigarstring = "1M1D" * l + "1M"
        self.assertEqual(len(a.cigartuples), 2 * l + 1)
        a.next_reference_id = 0
        a.next_reference_start = 0
        a.template_length = l
        a.query_qualities = pysam.qualitystring_to_array("1") * (l + 1)
        return a

    def check_read(self, read, mode="bam"):
        fn = get_temp_filename("tmp_largecigar.{}".format(mode))
        fn_reference = get_temp_filename("tmp_largecigar.fa")

        nrows = int(self.read_length * 2 / 80)

        s = "\n".join(["A" * 80 for x in range(nrows)])
        with open(fn_reference, "w") as outf:
            outf.write(">chr1\n{seq}\n>chr2\n{seq}\n".format(
                seq=s))

        if mode == "bam":
            write_mode = "wb"
        elif mode == "sam":
            write_mode = "w"
        elif mode == "cram":
            write_mode = "wc"

        with pysam.AlignmentFile(fn, write_mode,
                                 header=self.header,
                                 reference_filename=fn_reference) as outf:
            outf.write(read)
        with pysam.AlignmentFile(fn) as inf:
            ref_read = next(iter(inf))

        if mode == "cram":
            # in CRAM, the tag field is kept, while it is emptied by the BAM/SAM reader
            self.assertEqual(read.cigarstring, ref_read.cigarstring)
        else:
            self.assertEqual(read, ref_read)

        os.unlink(fn)
        os.unlink(fn_reference)

    def test_reading_writing_sam(self):
        read = self.build_read()
        self.check_read(read, mode="sam")

    def test_reading_writing_bam(self):
        read = self.build_read()
        self.check_read(read, mode="bam")

    @unittest.skip("fails on linux - https issue?")
    def test_reading_writing_cram(self):
        # The following test fails with htslib 1.9, but worked with previous versions.
        # Error is:
        # [W::find_file_url] Failed to open reference "https://www.ebi.ac.uk/ena/cram/md5/ac9fac7c3e9c476f74f1d0e47abc8be2": Input/output error
        # Error can be reproduced using samtools 1.9 command line.
        # Could be a conda configuration issue, see
        # https://github.com/bioconda/bioconda-recipes/issues/9056
        return
        read = self.build_read()
        self.check_read(read, mode="cram")
        
# SAM writing fails, as query length is 0
# class TestSanityCheckingSAM(TestSanityCheckingSAM):
#     mode = "w"

if __name__ == "__main__":
    print("starting tests")
    unittest.main()
    print("completed tests")
