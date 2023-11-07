import pysam
import pytest

from TestUtils import make_data_files, CBCF_DATADIR


def setUpModule():
    make_data_files(CBCF_DATADIR)


@pytest.fixture
def vcf_header():
    vcf_header = pysam.VariantHeader()
    vcf_header.add_samples("sample1", "sample2")
    vcf_header.contigs.add("1")
    return vcf_header

## segfault without coordinates

def test_ascii_annotation_can_be_added(vcf_header):
    vcf_header.formats.add("AN", 1, "String", "An annotation")
    record = vcf_header.new_record(
        contig="1",
        start=12,
        stop=13,
        samples=[
            {"AN": "anno1"},
            {"AN": "anno2"}])
    assert str(record)[:-1].split("\t")[-2:] == ["anno1", "anno2"]


def test_ascii_annotation_with_variable_length_can_be_added(vcf_header):
    vcf_header.formats.add("AN", 1, "String", "An annotation")
    record = vcf_header.new_record(
        contig="1",
        start=12,
        stop=13,
        samples=[
            {"AN": "anno1b"},
            {"AN": "anno1"}])
    assert str(record)[:-1].split("\t")[-2:] == ["anno1b", "anno1"]
    record = vcf_header.new_record(
        contig="1",
        start=12,
        stop=13,
        samples=[
            {"AN": "anno2"},
            {"AN": "anno2b"}])
    assert str(record)[:-1].split("\t")[-2:] == ["anno2", "anno2b"]
    

def test_unicode_annotation_can_be_added(vcf_header):
    vcf_header.formats.add("AN", 1, "String", "An annotation")
    record = vcf_header.new_record(
        contig="1",
        start=12,
        stop=13,
        samples=[
            {"AN": "anno1"},
            {"AN": "Friedrich-Alexander-Universit\u00E4t_Erlangen-N\u00FCrnberg"}])
    assert str(record)[:-1].split("\t")[-2:] == [
        "anno1",
        "Friedrich-Alexander-Universit\u00E4t_Erlangen-N\u00FCrnberg"]

def test_set_sample_alleles(vcf_header):
    vcf_header.formats.add('GT',1,'String',"Genotype") # id, number, type, description
    record = vcf_header.new_record(
        contig="1",
        start=20,
        stop=21,
        alleles=('A','T')
        )

    record.samples['sample1'].alleles = ('T', 'A')
    assert record.samples['sample1'].alleles  == ('T','A')

    # Empty record:
    record.samples['sample1'].alleles = (None, )
    assert record.samples['sample1'].alleles   == tuple()
    record.samples['sample1'].alleles = None
    assert record.samples['sample1'].alleles   == tuple()
    record.samples['sample1'].alleles = tuple()
    assert record.samples['sample1'].alleles   == tuple()

    # check error conditions:
    with pytest.raises(ValueError, match='One or more of the supplied sample alleles are not defined'):
        record.samples['sample1'].alleles = ('C', 'A')

    with pytest.raises(ValueError, match='Use .allele_indices to set integer allele indices'):
        record.samples['sample1'].alleles = (1, 0)
