from pysam.utils import PysamDispatcher

BCFTOOLS_DISPATCH = [
    "index",
    "annotate",
    "concat",
    "convert",
    "isec",
    "merge",
    "norm",
    "plugin",
    "query",
    "reheader",
    "sort",
    "view",
    "call",
    "consensus",
    "cnv",
    "csq",
    "filter",
    "gtcheck",
    "mpileup",
    "roh",
    "stats"]

# instantiate bcftools commands as python functions
for cmd in BCFTOOLS_DISPATCH:
    globals()[cmd] = PysamDispatcher("bcftools", cmd, None)
