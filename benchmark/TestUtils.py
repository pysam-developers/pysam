import os

BAM_DATADIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                           "..",
                                           "tests",
                                           "pysam_data"))

TABIX_DATADIR = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                             "..",
                                             "tests",
                                             "tabix_data"))
