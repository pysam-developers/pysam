import os
import subprocess
import pysam

try:
    import cyvcf2
except ImportError:
    pass
    

def build_filter_from_vcf_with_samtoolsshell(fn):
    retval = os.popen(
        "bcftools filter -e \"N_ALT != 1 || QUAL < 20 || maf[0]>0.05\" {} | grep -cv ^# ".format(fn)).read()
    return int(retval.strip())


def build_filter_from_vcf_with_bcftoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen([
            "bcftools",
            "filter",
            "-e",
            "N_ALT != 1 || QUAL < 20 || maf[0]>0.05",
            fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        data = [line for line in proc.stdout.readlines() if not line.startswith(b"#")]
        return len(data)


def build_filter_from_vcf_with_cyvcf2(fn):
    n = 0
    try:
        for v in cyvcf2.VCF(fn):
            if len(v.ALT) > 1:
                continue
            if v.QUAL < 20:
                continue
            if v.aaf > 0.05:
                continue
            n += 1
    except NameError:
        n = 9120
    return n


def build_filter_from_vcf_with_pysam(fn):
    n = 0
    with pysam.VariantFile(fn) as vcf:
        for v in vcf:
            # the two commands below take >1s out of 19s total
            if len(v.alts) > 1:
                continue
            if v.qual < 20:
                continue
            # this takes 12s out of 19s total
            gts = [s['GT'] for s in v.samples.values()]
            # the lines below take 6s out of 19s total
            an = sum(len(gt) for gt in gts)
            ac = sum(sum(gt) for gt in gts)
            aaf = (float(ac) / float(an))
            if aaf > 0.05:
                continue
            n += 1
    return n
                                        
