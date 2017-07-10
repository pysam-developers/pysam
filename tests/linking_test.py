"""test linking against pysam.
"""

import unittest
import os
import subprocess

from TestUtils import LINKDIR

class TestPackage(unittest.TestCase):

    package_name = "link_with_rpath"

    def setUp(self):
        self.workdir = os.path.join(LINKDIR, self.package_name)
        
    def test_package_can_be_installed(self):
        subprocess.check_output(
            "cd {} && rm -rf build && python setup.py install".format(self.workdir),
                shell=True)

    def test_package_tests_work(self):
        subprocess.check_output(
            "cd {} && python test_module.py".format(os.path.join(self.workdir, "tests")),
                shell=True)
            

if __name__ == "__main__":
    unittest.main()

                                
            
