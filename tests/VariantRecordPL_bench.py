import pytest

from pysam import VariantFile, VariantHeader


@pytest.mark.benchmark(min_rounds=1)
def test_access_pl_values_10_samples(benchmark, vcf_with_n_samples):
    vcf = vcf_with_n_samples(10)
    result = benchmark(access_pl_values, vcf)
    assert result == (0, 1, 2)


@pytest.mark.benchmark(min_rounds=1)
def test_access_pl_values_100_samples(benchmark, vcf_with_n_samples):
    vcf = vcf_with_n_samples(100)
    result = benchmark(access_pl_values, vcf)
    assert result == (0, 1, 2)


@pytest.mark.benchmark(min_rounds=1)
def test_access_pl_values_1000_samples(benchmark, vcf_with_n_samples):
    vcf = vcf_with_n_samples(1000)
    result = benchmark(access_pl_values, vcf)
    assert result == (0, 1, 2)


@pytest.mark.benchmark(min_rounds=1)
def test_access_pl_values_10000_samples(benchmark, vcf_with_n_samples):
    vcf = vcf_with_n_samples(10000)
    result = benchmark(access_pl_values, vcf)
    assert result == (0, 1, 2)


@pytest.mark.benchmark(min_rounds=1)
def test_access_pl_values_100000_samples(benchmark, vcf_with_n_samples):
    vcf = vcf_with_n_samples(100000)
    result = benchmark(access_pl_values, vcf)
    assert result == (0, 1, 2)


def access_pl_values(data):
    pl = None
    with VariantFile(data) as vcf:
        for record in vcf:
            for sample in record.samples.values():
                pl = sample["PL"]
    return pl


@pytest.fixture()
def vcf_with_n_samples(tmp_path):
    def vcf_with_n_samples(n_samples):
        vcfh = VariantHeader()
        for s in (f"s{i}" for i in range(n_samples)):
            vcfh.add_sample(s)
        vcfh.add_meta("contig", items=[("ID", "chr20")])
        vcfh.add_meta(
            "FORMAT",
            items=dict(ID="PL", Number="G", Type="Integer", Description="Phred Scaled Likelihood").items(),
        )
        vcfh.add_meta(
            "FORMAT",
            items=dict(ID="GT", Number="1", Type="String", Description="True Genotype").items(),
        )

        vcf = tmp_path / "large.vcf.gz"
        with VariantFile(vcf, "w", header=vcfh) as vcf_out:
            r = vcf_out.new_record(contig="chr20", start=1, stop=1, alleles=["C", "T"])
            for n, samp in enumerate(r.samples.values()):
                samp["GT"] = (0, 0)
                samp["PL"] = [0, 1, 2]
            vcf_out.write(r)
        return vcf
    return vcf_with_n_samples
