"""test linking against pysam.
"""

import unittest
import os
import subprocess
import pysam

from TestUtils import LINKDIR


def check_import(statement):
    try:
        output = subprocess.check_output(
            statement, stderr=subprocess.STDOUT, shell=True)
    except subprocess.CalledProcessError as exc:
        if b"ImportError" in exc.output:
            raise ImportError(
                "module could not be imported: {}".format(str(exc.output)))
        else:
            raise


def check_pass(statement):
    try:
        output = subprocess.check_output(
            statement, stderr=subprocess.STDOUT, shell=True)
    except subprocess.CalledProcessError as exc:
        raise ValueError("{}: {}".format(exc, exc.output))
    if b"FAILED" in output:
        raise ValueError("module tests failed")
    return True


@unittest.skipUnless(
    os.environ.get("PYSAM_LINKING_TESTS", None),
    "enable linking tests by setting PYSAM_LINKING_TESTS environment variable")
class TestLinking(unittest.TestCase):

    package_name = "link_with_rpath"

    def setUp(self):
        self.workdir = os.path.join(LINKDIR, self.package_name)

    def test_package_can_be_installed(self):
        subprocess.check_output(
            "cd {} && rm -rf build && python setup.py install".format(
                self.workdir),
            shell=True)


@unittest.skipUnless(
    os.environ.get("PYSAM_LINKING_TESTS", None),
    "enable linking tests by setting PYSAM_LINKING_TESTS environment variable")
class TestLinkWithRpath(TestLinking):

    package_name = "link_with_rpath"

    def test_package_tests_pass(self):
        self.assertTrue(check_pass(
            "cd {} && python test_module.py".format(os.path.join(self.workdir, "tests"))))


@unittest.skipUnless(
    os.environ.get("PYSAM_LINKING_TESTS", None),
    "enable linking tests by setting PYSAM_LINKING_TESTS environment variable")
class TestLinkWithoutRpath(TestLinking):

    package_name = "link_without_rpath"

    def test_package_tests_fail_on_import(self):

        self.assertRaises(
            ImportError,
            check_import,
            "cd {} && python test_module.py".format(os.path.join(self.workdir, "tests")))

    def test_package_tests_pass_if_ld_library_path_set(self):

        pysam_libraries = pysam.get_libraries()
        pysam_libdirs, pysam_libs = zip(
            *[os.path.split(x) for x in pysam_libraries])
        pysam_libdir = pysam_libdirs[0]

        self.assertTrue(check_pass(
            "export LD_LIBRARY_PATH={}:$PATH && cd {} && python test_module.py".format(
                pysam_libdir,
                os.path.join(self.workdir, "tests"))))


if __name__ == "__main__":
    unittest.main()
