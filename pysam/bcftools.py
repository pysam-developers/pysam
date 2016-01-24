from utils import PysamDispatcher

BCFTOOLS_DISPATCH = [
    "index",
    "annotate",
    "concat",
    "isec",
    "merge",
    "norm",
    "plugin",
    "query",
    "reheader",
    "view",
    "call",
    "consensus",
    "cnv",
    "filter",
    "gtcheck",
    "roh",
    "stats"]

# instantiate bcftools commands as python functions
for cmd in BCFTOOLS_DISPATCH:
    globals()[cmd] = PysamDispatcher("bcftools", cmd, None)
