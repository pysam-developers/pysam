# cython: embedsignature=True
# cython: profile=True
# adds doc-strings for sphinx
from pysam.chtslib cimport *

cpdef set_verbosity(int verbosity):
    u"""Set htslib's hts_verbose global variable to the specified value.
    """
    return hts_set_verbosity(verbosity)

cpdef get_verbosity():
    u"""Return the value of htslib's hts_verbose global variable.
    """
    return hts_get_verbosity()

__all__ = [
    "get_verbosity",
    "set_verbosity"]

