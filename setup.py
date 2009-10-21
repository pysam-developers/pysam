#!/usr/bin/python
'''

pysam
*****

'''

import os, sys, glob, shutil

name = "pysam"
version = "0.1.1"

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

classifiers = """
Development Status :: 2 - Alpha
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
    [ "pysam/csamtools.pyx" ]  +\
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
    'license': "MIT",
    'platforms': "ALL",
    'url': "http://code.google.com/p/pysam/",
    'py_modules': [
      "pysam/__init__", "pysam/Pileup"],
    'ext_modules': [pysam,],
    'cmdclass' : {'build_ext': build_ext} }

if __name__=='__main__':
   dist = setup(**metadata)
