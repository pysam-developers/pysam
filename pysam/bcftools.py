from utils import PysamDispatcher

BCFTOOLS_DISPATCH = {
    "stats": ("stats", None),
}


# instantiate samtools commands as python functions
for key, options in BCFTOOLS_DISPATCH.items():
    cmd, parser = options
    globals()[key] = PysamDispatcher("bcftools", cmd, parser)
