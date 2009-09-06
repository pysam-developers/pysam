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

import os, sys, glob, shutil

samtools_exclude = ( "bamtk.c", "razip.c", "bgzip.c" )
samtools_dest = os.path.abspath( "samtools" )

# copy samtools source
if len(sys.argv) >= 2 and sys.argv[1] == "import":
   if len(sys.argv) < 3: raise ValueError("missing PATH to samtools source directory")
   samtools_src = os.path.abspath( sys.argv[2] )
   if not os.path.exists( samtools_src ): raise IOError( "samtools src dir `%s` does not exist." % samtools_src )

   cfiles = glob.glob( os.path.join( samtools_src, "*.c" ) )
   hfiles = glob.glob( os.path.join( samtools_src, "*.h" ) )
   ncopied = 0
   for p in cfiles + hfiles:
      f = os.path.basename(p)
      if f in samtools_exclude: continue
      if os.path.exists( os.path.join( samtools_dest, f )): continue
      shutil.copy( p, samtools_dest )
      ncopied += 1
   print "installed latest source code from %s: %i files copied" % (samtools_src, ncopied)
   sys.exit(0)

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
    [ "pysam/csamtools.pyx",]  +\
       [ "pysam/%s" % x for x in (
             "pysam_util.c", )] +\
       glob.glob( os.path.join( "samtools", "*.c" ) ),
    library_dirs=[],
    include_dirs=[ "samtools", ],
    libraries=[ "z", ],
    language="c",
    )

metadata = {
    'name': name,
    'version': version,
    'description': "pysam", 
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
