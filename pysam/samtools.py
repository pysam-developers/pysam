import platform
from typing import (
    Callable,
    List,
    Tuple,
    Iterable,
    Union,
)
try:
    from typing import Final
    HAVE_FINAL = True
except ImportError:
    HAVE_FINAL = False

from pysam.utils import PysamDispatcher


# samtools command line options to export in python
_SAMTOOLS_DISPATCH = {
    # samtools 'documented' commands
    "view": ("view", ()),
    "head": ("head", ()),
    "sort": ("sort", ()),
    "mpileup": ("mpileup", ()),
    "consensus": ("consensus", ()),
    "depth": ("depth", ()),
    "faidx": ("faidx", ()),
    "fqidx": ("fqidx", ()),
    "tview": ("tview", ()),
    "index": ("index", ()),
    "idxstats": ("idxstats", ()),
    "fixmate": ("fixmate", ()),
    "flagstat": ("flagstat", ()),
    "calmd": ("calmd", ()),
    "merge": ("merge", ()),
    "markdup": ("markdup", ()),
    "rmdup": ("rmdup", ()),
    "reference": ("reference", ()),
    "reheader": ("reheader", ()),
    "reset": ("reset", ()),
    "cat": ("cat", ()),
    "targetcut": ("targetcut", ()),
    "phase": ("phase", ()),
    "bam2fq": ("bam2fq", ()),
    "dict": ("dict", ()),
    "addreplacerg": ("addreplacerg", ()),
    "pad2unpad": ("pad2unpad", ()),
    "depad": ("pad2unpad", ()),
    "bedcov": ("bedcov", ()),
    "coverage": ("coverage", ()),
    "bamshuf": ("bamshuf", ()),
    "collate": ("collate", ()),
    "stats": ("stats", ()),
    "fasta": ("fasta", ()),
    "fastq": ("fastq", ()),
    "cram_size": ("cram-size", ()),
    "quickcheck": ("quickcheck", ()),
    "split": ("split", ()),
    "flags": ("flags", ()),
    "ampliconclip": ("ampliconclip", ()),
    "ampliconstats": ("ampliconstats", ()),
    "version": ("version", ()),
    "fqimport": ("import", ()),
    "import_": ("import", ()),
    "samples": ("samples", ()),
}


def _wrap_command(
    dispatch: str,
    parsers: Iterable[Tuple[str, Callable[[Union[str, List[str]]], Union[str, List[str]]]]],
) -> PysamDispatcher:
    return PysamDispatcher("samtools", dispatch, parsers)


if not HAVE_FINAL:
    # python 3.7
    for key, options in _SAMTOOLS_DISPATCH.items():
        cmd, parser = options
        globals()[key] = PysamDispatcher("samtools", cmd, parser)

    __all__ = list(_SAMTOOLS_DISPATCH)
else:
    # python >=3.8
    view: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["view"][0], _SAMTOOLS_DISPATCH["view"][1])

    head: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["head"][0], _SAMTOOLS_DISPATCH["head"][1])

    sort: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["sort"][0], _SAMTOOLS_DISPATCH["sort"][1])

    mpileup: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["mpileup"][0], _SAMTOOLS_DISPATCH["mpileup"][1])

    consensus: Final[PysamDispatcher] = _wrap_command(
        _SAMTOOLS_DISPATCH["consensus"][0],
        _SAMTOOLS_DISPATCH["consensus"][1],
    )

    depth: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["depth"][0], _SAMTOOLS_DISPATCH["depth"][1])

    faidx: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["faidx"][0], _SAMTOOLS_DISPATCH["faidx"][1])

    fqidx: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["fqidx"][0], _SAMTOOLS_DISPATCH["fqidx"][1])

    tview: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["tview"][0], _SAMTOOLS_DISPATCH["tview"][1])

    index: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["index"][0], _SAMTOOLS_DISPATCH["index"][1])

    idxstats: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["idxstats"][0], _SAMTOOLS_DISPATCH["idxstats"][1])

    fixmate: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["fixmate"][0], _SAMTOOLS_DISPATCH["fixmate"][1])

    flagstat: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["flagstat"][0], _SAMTOOLS_DISPATCH["flagstat"][1])

    calmd: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["calmd"][0], _SAMTOOLS_DISPATCH["calmd"][1])

    merge: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["merge"][0], _SAMTOOLS_DISPATCH["merge"][1])

    markdup: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["markdup"][0], _SAMTOOLS_DISPATCH["markdup"][1])

    rmdup: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["rmdup"][0], _SAMTOOLS_DISPATCH["rmdup"][1])

    reference: Final[PysamDispatcher] = _wrap_command(
        _SAMTOOLS_DISPATCH["reference"][0],
        _SAMTOOLS_DISPATCH["reference"][1],
    )

    reheader: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["reheader"][0], _SAMTOOLS_DISPATCH["reheader"][1])

    reset: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["reset"][0], _SAMTOOLS_DISPATCH["reset"][1])

    cat: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["cat"][0], _SAMTOOLS_DISPATCH["cat"][1])

    targetcut: Final[PysamDispatcher] = _wrap_command(
        _SAMTOOLS_DISPATCH["targetcut"][0],
        _SAMTOOLS_DISPATCH["targetcut"][1],
    )

    phase: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["phase"][0], _SAMTOOLS_DISPATCH["phase"][1])

    bam2fq: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["bam2fq"][0], _SAMTOOLS_DISPATCH["bam2fq"][1])

    dict: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["dict"][0], _SAMTOOLS_DISPATCH["dict"][1])

    addreplacerg: Final[PysamDispatcher] = _wrap_command(
        _SAMTOOLS_DISPATCH["addreplacerg"][0],
        _SAMTOOLS_DISPATCH["addreplacerg"][1],
    )

    pad2unpad: Final[PysamDispatcher] = _wrap_command(
        _SAMTOOLS_DISPATCH["pad2unpad"][0],
        _SAMTOOLS_DISPATCH["pad2unpad"][1],
    )

    depad: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["depad"][0], _SAMTOOLS_DISPATCH["depad"][1])

    bedcov: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["bedcov"][0], _SAMTOOLS_DISPATCH["bedcov"][1])

    coverage: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["coverage"][0], _SAMTOOLS_DISPATCH["coverage"][1])

    bamshuf: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["bamshuf"][0], _SAMTOOLS_DISPATCH["bamshuf"][1])

    collate: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["collate"][0], _SAMTOOLS_DISPATCH["collate"][1])

    stats: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["stats"][0], _SAMTOOLS_DISPATCH["stats"][1])

    fasta: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["fasta"][0], _SAMTOOLS_DISPATCH["fasta"][1])

    fastq: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["fastq"][0], _SAMTOOLS_DISPATCH["fastq"][1])

    cram_size: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["cram_size"][0], _SAMTOOLS_DISPATCH["cram_size"][1])

    quickcheck: Final[PysamDispatcher] = _wrap_command(
        _SAMTOOLS_DISPATCH["quickcheck"][0],
        _SAMTOOLS_DISPATCH["quickcheck"][1],
    )

    split: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["split"][0], _SAMTOOLS_DISPATCH["split"][1])

    flags: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["flags"][0], _SAMTOOLS_DISPATCH["flags"][1])

    ampliconclip: Final[PysamDispatcher] = _wrap_command(
        _SAMTOOLS_DISPATCH["ampliconclip"][0],
        _SAMTOOLS_DISPATCH["ampliconclip"][1],
    )

    ampliconstats: Final[PysamDispatcher] = _wrap_command(
        _SAMTOOLS_DISPATCH["ampliconstats"][0],
        _SAMTOOLS_DISPATCH["ampliconstats"][1],
    )

    version: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["version"][0], _SAMTOOLS_DISPATCH["version"][1])

    fqimport: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["fqimport"][0], _SAMTOOLS_DISPATCH["fqimport"][1])

    import_: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["import_"][0], _SAMTOOLS_DISPATCH["import_"][1])

    samples: Final[PysamDispatcher] = _wrap_command(_SAMTOOLS_DISPATCH["samples"][0], _SAMTOOLS_DISPATCH["samples"][1])
