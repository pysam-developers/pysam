import pysam

def test_idxstats_parse_split_lines():
    bam_filename = "./pysam_data/ex2.bam"
    lines = pysam.idxstats(bam_filename, split_lines=True)  # Test pysam 0.8.X style output, which returns a list of lines
    for line in lines:
        _seqname, _seqlen, nmapped, _nunmapped = line.split()


def test_bedcov_split_lines():
    bam_filename = "./pysam_data/ex1.bam"
    bed_filename = "./pysam_data/ex1.bed"
    lines = pysam.bedcov(bed_filename, bam_filename, split_lines=True)  # Test pysam 0.8.X style output, which returns a list of lines
    for line in lines:
        fields = line.split('\t')
        assert len(fields) in [4, 5], "bedcov should give tab delimited output with 4 or 5 fields.  Split line (%s) gives %d fields." % (fields, len(fields))


def test_idxstats_parse():
    bam_filename = "./pysam_data/ex2.bam"
    idxstats_string = pysam.idxstats(bam_filename, split_lines=False)  # Test pysam 0.9.X style output, which returns a string that needs to be split by \n
    lines = idxstats_string.splitlines()
    for line in lines:
        splt = line.split("\t")
        _seqname, _seqlen, nmapped, _nunmapped = splt


def test_bedcov():
    bam_filename = "./pysam_data/ex1.bam"
    bed_filename = "./pysam_data/ex1.bed"
    bedcov_string = pysam.bedcov(bed_filename, bam_filename, split_lines=False)  # Test pysam 0.9.X style output, which returns a string that needs to be split by \n
    lines = bedcov_string.splitlines()
    for line in lines:
        fields = line.split('\t')
        assert len(fields) in [4, 5], "bedcov should give tab delimited output with 4 or 5 fields.  Split line (%s) gives %d fields." % (fields, len(fields))
