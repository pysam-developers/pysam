#!/usr/bin/python
"""\

NCL
***

Nested contained lists are a way to index segment data. See
Alekseyenko & Lee (2007) 
(http://bioinformatics.oxfordjournals.org/cgi/content/full/23/11/1386)

The following code was taken from the author's implemetation
in pygr (http://sourceforge.net/projects/pygr) and modified. The
modifications include:

1. remove target coordinates, only target_id is kept.

"""\

import os
import sys
from distutils.core import setup, Extension
from Pyrex.Distutils import build_ext

name = "pysam"
version = "0.1"

classifiers = """
Development Status :: 1 - Pre Alpha
Operating System :: MacOS :: MacOS X
Operating System :: Microsoft :: Windows :: Windows NT/2000
Operating System :: OS Independent
Operating System :: POSIX
Operating System :: POSIX :: Linux
Operating System :: Unix
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Scientific/Engineering :: Bioinformatics
"""

pysam = Extension(
    "pysam/csamtools",                   # name of extension
    [ "pysam/csamtools.pyx",] + \
    [ "pysam/%s" % x for x in [ 
      "bgzf.c",
      "kstring.c",
      "bam_aux.c",
      "bam.c",
      "bam_import.c",
      "sam.c",
      "bam_index.c",
      "bam_pileup.c",
      "bam_lpileup.c",
      "bam_md.c",
      "glf.c",
      "razf.c",
      "faidx.c",
      "knetfile.c",
      "bam_sort.c" ]],
      library_dirs=[],
      libraries=[ "z", ],
      language="c",
    )

metadata = {
    'name': name,
    'version': version,
    'description': "NCL", 
    'long_description': __doc__,
    'author': "Andreas Heger",
    'author_email': "andreas.heger@gmail.com",
    'license': "GPL",
    'platforms': "ALL",
    'url': "TBD",
    'py_modules': [
      "pysam/__init__", ],
    'ext_modules': [pysam,],
    'cmdclass' : {'build_ext': build_ext} }

if __name__=='__main__':
   dist = setup(**metadata)
