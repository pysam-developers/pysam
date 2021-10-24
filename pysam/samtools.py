from pysam.utils import PysamDispatcher

# samtools command line options to export in python
SAMTOOLS_DISPATCH = {
    # samtools 'documented' commands
    "view": ("view", None),
    "sort": ("sort", None),
    "mpileup": ("mpileup", None),
    "depth": ("depth", None),
    "faidx": ("faidx", None),
    "fqidx": ("fqidx", None),
    "tview": ("tview", None),
    "index": ("index", None),
    "idxstats": ("idxstats", None),
    "fixmate": ("fixmate", None),
    "flagstat": ("flagstat", None),
    "calmd": ("calmd", None),
    "merge": ("merge", None),
    "markdup": ("markdup", None),
    "rmdup": ("rmdup", None),
    "reheader": ("reheader", None),
    "cat": ("cat", None),
    "targetcut": ("targetcut", None),
    "phase": ("phase", None),
    "bam2fq": ("bam2fq", None),
    "dict": ("dict", None),
    "addreplacerg": ("addreplacerg", None),
    "pad2unpad": ("pad2unpad", None),
    "depad": ("pad2unpad", None),
    "bedcov": ("bedcov", None),
    "coverage": ("coverage", None),
    "bamshuf": ("bamshuf", None),
    "collate": ("collate", None),
    "stats": ("stats", None),
    "fasta": ("fasta", None),
    "fastq": ("fastq", None),
    "quickcheck": ("quickcheck", None),
    "split": ("split", None),
    "flags": ("flags", None),
    "ampliconclip": ("ampliconclip", None),
    "ampliconstats": ("ampliconstats", None),
    "version": ("version", None),
    "fqimport": ("import", None),
    "samples": ("samples", None),
}

# instantiate samtools commands as python functions
for key, options in SAMTOOLS_DISPATCH.items():
    cmd, parser = options
    globals()[key] = PysamDispatcher("samtools", cmd, parser)

__all__ = list(SAMTOOLS_DISPATCH)
