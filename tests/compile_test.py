'''
compile_test.py - check pyximport functionality with pysam
==========================================================

test script for checking if compilation against
pysam and tabix works.
'''

# clean up previous compilation
import os
import pytest
import pysam
from TestUtils import make_data_files, BAM_DATADIR, TABIX_DATADIR


def setUpModule():
    make_data_files(BAM_DATADIR)
    make_data_files(TABIX_DATADIR)


try:
    os.unlink('tests/_compile_test.c')
    os.unlink('tests/_compile_test.pyxbldc')
except OSError:
    pass

NO_PYXIMPORT = False
try:
    import pyximport
    pyximport.install(build_in_temp=False)
    import _compile_test
except:
    NO_PYXIMPORT = True


@pytest.mark.skipif(NO_PYXIMPORT, reason="no pyximport")
def test_bam():

    input_filename = os.path.join(BAM_DATADIR, "ex1.bam")
    nread = _compile_test.testCountBAM(
        pysam.Samfile(input_filename))
    assert nread == 3270


@pytest.mark.skipif(NO_PYXIMPORT, "no pyximport")
def test_gtf():

    input_filename = os.path.join(TABIX_DATADIR, reason="example.gtf.gz")

    nread = _compile_test.testCountGTF(
        pysam.Tabixfile(input_filename))
    assert nread == 237
