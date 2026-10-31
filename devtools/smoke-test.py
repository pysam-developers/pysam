#!/usr/bin/env python3
"""
Quick check that an installed pysam works, for platforms where the full
test suite does not yet run. Writes test.bam and test.vcf.gz in the
current directory. If a BAM file URL is given, also reads it via libcurl.
"""
import sys

import pysam, pysam.bcftools
with pysam.AlignmentFile("test.bam", "wb", reference_names=["chr1"], reference_lengths=[100]) as outf:
    read = pysam.AlignedSegment(outf.header)
    read.query_name, read.reference_id, read.reference_start = "read1", 0, 10
    read.cigarstring, read.query_sequence = "4M", "ACGT"
    outf.write(read)
pysam.samtools.index("test.bam")
with pysam.AlignmentFile("test.bam") as inf:
    assert [r.query_name for r in inf.fetch("chr1", 0, 50)] == ["read1"]
assert pysam.samtools.view("-c", "test.bam", "chr1").strip() == "1"
header = pysam.VariantHeader()
header.contigs.add("chr1", length=100)
with pysam.VariantFile("test.vcf.gz", "wz", header=header) as outf:
    outf.write(outf.new_record(contig="chr1", start=10, alleles=("A", "T")))
assert pysam.bcftools.view("-H", "test.vcf.gz").count("\n") == 1
if len(sys.argv) > 1:
    with pysam.AlignmentFile(sys.argv[1]) as inf:
        print(sys.argv[1], "has", sum(1 for _ in inf), "reads")
print(pysam.__version__, pysam.config.HTSLIB, "OK")
