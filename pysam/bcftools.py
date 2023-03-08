try:
    from typing import Final
    HAVE_FINAL = True
except ImportError:
    HAVE_FINAL = False

from pysam.utils import PysamDispatcher

_BCFTOOLS_DISPATCH = [
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
    "head",
    "call",
    "consensus",
    "cnv",
    "csq",
    "filter",
    "gtcheck",
    "mpileup",
    "roh",
    "stats"]


def _wrap_command(dispatch: str) -> PysamDispatcher:
    return PysamDispatcher("bcftools", dispatch, ())


if not HAVE_FINAL:
    # instantiate bcftools commands as python functions
    for cmd in _BCFTOOLS_DISPATCH:
        globals()[cmd] = PysamDispatcher("bcftools", cmd, None)
else:
    # python >=3.8
    index: Final[PysamDispatcher] = _wrap_command("index")
    annotate: Final[PysamDispatcher] = _wrap_command("annotate")
    concat: Final[PysamDispatcher] = _wrap_command("concat")
    convert: Final[PysamDispatcher] = _wrap_command("convert")
    isec: Final[PysamDispatcher] = _wrap_command("isec")
    merge: Final[PysamDispatcher] = _wrap_command("merge")
    norm: Final[PysamDispatcher] = _wrap_command("norm")
    plugin: Final[PysamDispatcher] = _wrap_command("plugin")
    query: Final[PysamDispatcher] = _wrap_command("query")
    reheader: Final[PysamDispatcher] = _wrap_command("reheader")
    sort: Final[PysamDispatcher] = _wrap_command("sort")
    view: Final[PysamDispatcher] = _wrap_command("view")
    head: Final[PysamDispatcher] = _wrap_command("head")
    call: Final[PysamDispatcher] = _wrap_command("call")
    consensus: Final[PysamDispatcher] = _wrap_command("consensus")
    cnv: Final[PysamDispatcher] = _wrap_command("cnv")
    csq: Final[PysamDispatcher] = _wrap_command("csq")
    filter: Final[PysamDispatcher] = _wrap_command("filter")
    gtcheck: Final[PysamDispatcher] = _wrap_command("gtcheck")
    mpileup: Final[PysamDispatcher] = _wrap_command("mpileup")
    roh: Final[PysamDispatcher] = _wrap_command("roh")
    stats: Final[PysamDispatcher] = _wrap_command("stats")
