import pysam
import pytest

try:
    from pathlib import Path
except ImportError:
    Path = None

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
