import pysam

def test_idxstats_parse():
    bam_filename = "./pysam_data/ex2.bam"
    lines = pysam.idxstats(bam_filename)
    for line in lines:
        _seqname, _seqlen, nmapped, _nunmapped = line.split()

def test_bedcov():
    bam_filename = "./pysam_data/ex1.bam"
    bed_filename = "./pysam_data/ex1.bed"
    lines = pysam.bedcov(bed_filename, bam_filename)
    for line in lines:
        fields = line.split('\t')
        assert len(fields) in [4, 5], "bedcov should give tab delimited output with 4 or 5 fields.  Split line (%s) gives %d fields." % (fields, len(fields))
