import glob
import sys
import os

from setuptools import setup, find_packages, Extension
from cy_build import CyExtension as Extension, cy_build_ext as build_ext

import pysam

test_module_suffix = os.path.dirname(os.path.abspath(__file__)).split(os.sep)[-1]
test_module_name = "PysamTestModule_{}".format(test_module_suffix)

TestModule = Extension(
    "{}.BuildRead".format(test_module_name),
    ["{}/BuildRead.pyx".format(test_module_name)],
    include_dirs=pysam.get_include(),
    extra_link_args=pysam.get_libraries(),
    define_macros=pysam.get_defines(),
)

setup(
    name=test_module_name,
    version='0.1',
    packages=find_packages(),
    package_dir={test_module_name: test_module_name},
    ext_modules=[TestModule],
    cmdclass={'build_ext': build_ext},
)
