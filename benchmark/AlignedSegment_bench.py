"""Benchmarking module for AlignedSegment functionality"""

import timeit

iterations = 10000
repeats = 5

setup_binary_tag = """
import pysam
import array
read = pysam.AlignedSegment()
read.set_tag('FZ', array.array('H', range(1000)))
"""

setup_binary_tag_from_file = """
import pysam
with pysam.AlignmentFile("../tests/pysam_data/example_btag.bam", "rb") as inf:
    read = inf.fetch().next()
"""

def test_read_binary_get_tag(read):
    tags = read.get_tag('FZ')

def test_read_and_process_binary_get_tag(read):
    tags = sum(read.get_tag('FZ'))

tests = (
    ("test_read_binary_get_tag", "setup_binary_tag"),
    ("test_read_binary_get_tag", "setup_binary_tag_from_file"),
    ("test_read_and_process_binary_get_tag", "setup_binary_tag"),
    )

for repeat in range(repeats):
    print ("# repeat=", repeat)
    for testf, setup_name in tests:
        setup = locals()[setup_name]
        setup += """\nfrom __main__ import %s""" % testf
        #try:
        t = timeit.timeit("%s(read)" % testf, number=iterations, setup=setup)
        #except AttributeError, msg:
        #    print msg
        #    continue
        print ("%5.2f\t%s\t%s" % (t,testf, setup_name))
