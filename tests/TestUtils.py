import sys
import os
import pysam
import difflib


IS_PYTHON3 = sys.version_info[0] >= 3

if IS_PYTHON3:
    from itertools import zip_longest
    from urllib.request import urlopen
else:
    from itertools import izip as zip_longest
    from urllib2 import urlopen


def checkBinaryEqual(filename1, filename2):
    '''return true if the two files are binary equal.
    '''
    if os.path.getsize(filename1) != os.path.getsize(filename2):
        return False

    infile1 = open(filename1, "rb")
    infile2 = open(filename2, "rb")

    def chariter(infile):
        while 1:
            c = infile.read(1)
            if c == b"":
                break
            yield c

    found = False
    for c1, c2 in zip_longest(chariter(infile1), chariter(infile2)):
        if c1 != c2:
            break
    else:
        found = True

    infile1.close()
    infile2.close()
    return found


def checkSamtoolsViewEqual(filename1, filename2,
                           without_header=False):
    '''return true if the two files are equal in their
    content through samtools view.
    '''

    # strip MD and NM tags, as not preserved in CRAM files
    args = ["-x", "MD", "-x", "NM"]
    if not without_header:
        args.append("-h")

    lines1 = pysam.view(*(args + [filename1]))
    lines2 = pysam.view(*(args + [filename2]))

    if len(lines1) != len(lines2):
        return False

    if lines1 != lines2:
        # line by line comparison
        # sort each line, as tags get rearranged between
        # BAM/CRAM
        for n, pair in enumerate(zip(lines1, lines2)):
            l1, l2 = pair
            l1 = sorted(l1[:-1].split("\t"))
            l2 = sorted(l2[:-1].split("\t"))
            if l1 != l2:
                print "mismatch in line %i" % n
                print l1
                print l2
                return False
        else:
            return False

    return True


def checkURL(url):
    '''return True if URL is available.

    A URL might not be available if it is the wrong URL
    or there is no connection to the URL.
    '''
    try:
        urlopen(url, timeout=1)
        return True
    except:
        return False


def checkFieldEqual(cls, read1, read2, exclude=[]):
    '''check if two reads are equal by comparing each field.'''

    # add the . for refactoring purposes.
    for x in (".query_name",
              ".query_sequence",
              ".flag",
              ".reference_id",
              ".reference_start",
              ".mapping_quality",
              ".cigartuples",
              ".next_reference_id",
              ".next_reference_start",
              ".template_length",
              ".query_length",
              ".query_qualities",
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
        cls.assertEqual(getattr(read1, n), getattr(read2, n),
                        "attribute mismatch for %s: %s != %s" %
                        (n, getattr(read1, n), getattr(read2, n)))
