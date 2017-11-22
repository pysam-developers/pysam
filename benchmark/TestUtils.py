import sys
import os

BAM_DATADIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                           "..",
                                           "tests",
                                           "pysam_data"))

TABIX_DATADIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                             "..",
                                             "tests",
                                             "tabix_data"))

IS_PYTHON3 = sys.version_info[0] >= 3


if IS_PYTHON3:
    def force_str(s):
        try:
            return s.decode('ascii')
        except AttributeError:
            return s
    def force_bytes(s):
        try:
            return s.encode('ascii')
        except AttributeError:
            return s
else:
    def force_str(s):
        return s
    def force_bytes(s):
        return s

def flatten_nested_list(l):
    return [i for ll in l for i in ll]
