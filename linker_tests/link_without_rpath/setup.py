import glob
import sys
import os

from setuptools import setup, find_packages, Extension
from Cython.Distutils import build_ext

import pysam

test_module_suffix = os.path.dirname(os.path.abspath(__file__)).split(os.sep)[-1]
test_module_name = "PysamTestModule_{}".format(test_module_suffix)

pysam_libraries = pysam.get_libraries()
pysam_libdirs, pysam_libs = zip(*[os.path.split(x) for x in pysam_libraries])
pysam_libdir = pysam_libdirs[0]
# remove lib and .so
pysam_libs = [x[3:-3] for x in pysam_libs]

TestModule = Extension(
    "{}.BuildRead".format(test_module_name),
    ["{}/BuildRead.pyx".format(test_module_name)],
    include_dirs=pysam.get_include(),
    library_dirs=[pysam_libdir],
    libraries=pysam_libs,
    language="C",
)

setup(
    name=test_module_name,
    version='0.1',
    packages=find_packages(),
    package_dir={test_module_name: test_module_name},
    ext_modules=[TestModule],
    cmdclass={'build_ext': build_ext},
)
