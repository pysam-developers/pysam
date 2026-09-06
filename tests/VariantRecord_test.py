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

# segfault without coordinates


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
    vcf_header.formats.add('GT', 1, 'String', "Genotype")  # id, number, type, description
    record = vcf_header.new_record(
        contig="1",
        start=20,
        stop=21,
        alleles=('A', 'T')
        )

    record.samples['sample1'].alleles = ('T', 'A')
    assert record.samples['sample1'].alleles  == ('T', 'A')

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


def test_update_samples(vcf_header):
    vcf_header.formats.add("DP", 1, "Integer", "Read depth")
    record = vcf_header.new_record(contig="1", start=20, stop=30, alleles=("A", "T"), samples=[{"DP": 10}, {"DP": 20}])
    with pytest.raises(TypeError): record.samples.update({"sample1": "foo"})
    with pytest.raises(TypeError): record.samples.pop("sample2")


def test_sample_update_dict(vcf_header):
    vcf_header.formats.add("DA", 1, "Integer", "Misc")
    vcf_header.formats.add("DB", 1, "String", "Misc")
    rec = vcf_header.new_record(contig="1", start=30, stop=40)

    rec.samples["sample1"].update({"DA": 28, "DB": "test"})
    assert dict(rec.samples["sample1"]) == {"DA": 28, "DB": "test"}


def test_sample_update_iterable(vcf_header):
    vcf_header.formats.add("YA", 1, "Integer", "Misc")
    vcf_header.formats.add("YB", 1, "String", "Misc")
    rec = vcf_header.new_record(contig="1", start=30, stop=40)

    def yield_formats():
        yield ("YA", 28)
        yield ("YB", "test")

    rec.samples["sample1"].update(yield_formats())
    assert dict(rec.samples["sample1"]) == {"YA": 28, "YB": "test"}


def test_sample_update_keywords(vcf_header):
    vcf_header.formats.add("KA", 1, "Integer", "Misc")
    vcf_header.formats.add("KB", 1, "String", "Misc")
    rec = vcf_header.new_record(contig="1", start=30, stop=40)

    rec.samples["sample1"].update(KA=28, KB="test")
    assert dict(rec.samples["sample1"]) == {"KA": 28, "KB": "test"}


def test_sample_update_dict_and_keywords(vcf_header):
    vcf_header.formats.add("DA", 1, "Integer", "Misc")
    vcf_header.formats.add("DB", 1, "String", "Misc")
    vcf_header.formats.add("KA", 1, "Integer", "Misc")
    vcf_header.formats.add("KB", 1, "String", "Misc")
    rec = vcf_header.new_record(contig="1", start=30, stop=40)

    rec.samples["sample1"].update({"DA": 28, "DB": "test"}, KA=28, KB="test")
    assert dict(rec.samples["sample1"]) == {"DA": 28, "DB": "test", "KA": 28, "KB": "test"}


def test_repeated_new_record(vcf_header):
    vcf_header.formats.add('GT', 1, 'String', "Genotype")
    vcf_header.formats.add("AA", 1, "String", "An annotation")
    vcf_header.formats.add("BB", 1, "String", "Another annotation")

    data = {'id': 'INS_1', 'contig': '1', 'start': 10, 'stop': 15, 'alleles': ['A', 'TCGA'],
            'samples': [{'AA': ('one'), 'GT': (0, 1), 'BB': ('two')},
                        {'GT': (1, 0), 'BB': ('three')}]}

    record1 = vcf_header.new_record(**data)
    assert '\tGT:' in str(record1)  # Verify that GT is output first
    assert record1.samples['sample1'].alleles == ('A', 'TCGA')
    assert record1.samples['sample2'].alleles == ('TCGA', 'A')

    record2 = vcf_header.new_record(**data)
    assert '\tGT:' in str(record2)  # Verify that GT is actually emitted and is output first
    assert record2.samples['sample1'].alleles == ('A', 'TCGA')
    assert record2.samples['sample2'].alleles == ('TCGA', 'A')


def test_format_number_P(vcf_header):
    vcf_header.formats.add("GT", 1, "String", "Genotype")
    vcf_header.formats.add("X", "P", "Integer", "Per-genotype-allele value")

    samples=[{"GT": (0, 1), "X": (20, 43)}, {"GT": (1,), "X": 37}]
    rec = vcf_header.new_record(contig="1", start=10, alleles=["A", "T"], samples=samples)

    assert rec.samples["sample1"]["X"] == (20, 43)
    assert rec.samples["sample2"]["X"] == (37,)


def test_format_number_LA_LR_LG(vcf_header):
    vcf_header.formats.add("GT", 1, "String", "Genotype")
    vcf_header.formats.add("LAA", ".", "Integer", "Local alleles")
    vcf_header.formats.add("LAD", "LR", "Integer", "Local-allele representation of AD")
    vcf_header.formats.add("LEC", "LA", "Integer", "Local EC")
    vcf_header.formats.add("LPL", "LG", "Integer", "Local PL")

    samples=[{"GT": (0, 1), "LAA": (1, 2),    "LAD": (50, 10, 20),     "LEC": (11, 21),     "LPL": (1,2,3,4,5,6)},
             {"GT": (1,),   "LAA": (1, 3, 4), "LAD": (50, 10, 30, 40), "LEC": (11, 31, 41), "LPL": (1,2,3,4)}]
    rec = vcf_header.new_record(contig="1", start=10, alleles=["A", "C", "G", "T"], samples=samples)

    assert len(rec.samples["sample1"]["LAA"]) == 2
    #assert rec.samples["sample1"]["LAD"] == (50, 10, 20)
    #assert rec.samples["sample1"]["LEC"] == (11, 21)
    assert rec.samples["sample1"]["LPL"] == (1,2,3,4,5,6)

    assert len(rec.samples["sample2"]["LAA"]) == 3
    assert rec.samples["sample2"]["LAD"] == (50, 10, 30, 40)
    assert rec.samples["sample2"]["LEC"] == (11, 31, 41)
    assert rec.samples["sample2"]["LPL"] == (1,2,3,4)
