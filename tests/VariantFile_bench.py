"""Benchmarking module for AlignmentFile functionality"""
import pytest


from VariantFileFetchTestUtils import *


GENOMES_URL = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

CHROM = 22


@pytest.fixture
def genomes_data():
    url = GENOMES_URL.format(chrom=CHROM)
    fn = os.path.basename(url)
    print(fn)
    if not os.path.exists(fn):
        os.system("wget {}".format(url))
    if not os.path.exists(fn + ".tbi"):
        os.system("wget {}".format(url + ".tbi"))

    fn_small = "small.vcf.gz"
    if not os.path.exists(fn_small):
        os.system("bcftools view {} | head -n 10000 | bgzip > {}".format(fn, fn_small))
        os.system("tabix -f -p vcf {}".format(fn_small))
        
    return fn_small


@pytest.mark.benchmark(min_rounds=1)
def test_build_filter_from_vcf_with_bcftoolsshell(benchmark, genomes_data):
    result = benchmark(build_filter_from_vcf_with_samtoolsshell, genomes_data)
    assert result == 9120


@pytest.mark.benchmark(min_rounds=1)
def test_build_filter_from_vcf_with_bcftoolpipe(benchmark, genomes_data):
    result = benchmark(build_filter_from_vcf_with_bcftoolspipe, genomes_data)
    assert result == 9120


@pytest.mark.benchmark(min_rounds=1)
def test_build_filter_from_vcf_with_cyvcf2(benchmark, genomes_data):
    result = benchmark(build_filter_from_vcf_with_cyvcf2, genomes_data)
    # note: inconsistent with bcftools
    assert result == 9114


@pytest.mark.benchmark(min_rounds=1)
def test_build_filter_from_vcf_with_pysam(benchmark, genomes_data):
    result = benchmark(build_filter_from_vcf_with_pysam, genomes_data)
    # note: inconsistent with bcftools
    assert result == 9114
    

    
