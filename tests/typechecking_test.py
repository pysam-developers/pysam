import inspect
import os
import pathlib
import pytest
import re
from typing import TYPE_CHECKING

import pysam
import pysam.samtools
import pysam.bcftools

from TestUtils import BAM_DATADIR, CBCF_DATADIR, make_data_files

try:
    import mypy.api
except ImportError:
    pytest.skip('mypy API not available', allow_module_level=True)

PREAMBLE = """
import os
import pathlib
import re
from typing import TYPE_CHECKING

import pysam
import pysam.samtools
import pysam.bcftools

def typecheck(check_locals = True): pass
"""

MYPY_OPTIONS = ['--no-incremental', '--no-error-summary', '--follow-imports', 'silent']


def typecheck(check_locals: bool = True):
    myframe = inspect.currentframe()
    if not myframe: pytest.skip('current stack frame not available')
    caller = myframe.f_back
    if not caller: pytest.skip('caller stack frame not available')

    code = inspect.getsource(caller)
    stdout, stderr, status = mypy.api.run(MYPY_OPTIONS + ['--command', PREAMBLE + code])
    assert status == 0, f'mypy failed:\n{stdout}{stderr}'

    def plain_mypy(s):
        s = re.sub(r'builtins\.', '', s)
        return s

    types = {}
    for line in stdout.splitlines():
        m = re.search(r'note:   *(\w+): ([\w.]*)', line)
        if m: types[m.group(1)] = plain_mypy(m.group(2))

    def plain_repr(s):
        s = re.sub(r"<class '([^']*)'>", r'\1', s)
        return s

    if check_locals:
        for var, vartype in types.items():
            assert vartype == plain_repr(repr(type(caller.f_locals[var]))), f'Incorrect type for {var!r}'

    return types


def setUpModule():
    make_data_files(BAM_DATADIR)


@pytest.fixture
def aln():
    header = pysam.AlignmentHeader.from_references(['chr1', 'chr2'], [50000, 20000])
    a = pysam.AlignedSegment(header)
    a.query_name = 'read_one'
    a.flag = pysam.FPAIRED | pysam.FREAD1
    a.reference_id = 0
    a.reference_start = 1000
    a.mapping_quality = 20
    a.next_reference_id = 1
    a.next_reference_start = 5000
    a.template_length = 0
    a.cigartuples = ((pysam.CMATCH, 6), (pysam.CINS, 4), (pysam.CMATCH, 2))
    a.query_sequence = 'ATGCATGCATGC'
    a.query_qualities = pysam.qualitystring_to_array('abcdefghijkl')
    return a


@pytest.fixture
def sam_fname():
    return os.path.join(BAM_DATADIR, 'ex3.sam')


@pytest.fixture
def bam_fname():
    return os.path.join(BAM_DATADIR, 'ex1.bam')


@pytest.fixture
def vcf_fname():
    # This file is a committed source file, so we don't need to make_data_files() here
    return os.path.join(CBCF_DATADIR, 'example_vcf43.vcf')


def test_AlignedSegment_is_mapped(aln: pysam.AlignedSegment) -> None:
    rmap = aln.is_mapped
    mmap = aln.mate_is_mapped
    rfwd = aln.is_forward
    mfwd = aln.mate_is_forward

    if TYPE_CHECKING: reveal_locals()
    types = typecheck()


def test_AlignmentHeader_get_reference_length(sam_fname: str) -> None:
    inf = pysam.AlignmentFile(sam_fname)
    n1 = inf.get_reference_length('chr1')
    hdr = inf.header
    n2 = hdr.get_reference_length('chr1')

    if TYPE_CHECKING: reveal_locals()
    types = typecheck()


def test_pileup_iterator_column(bam_fname: str) -> None:
    inf = pysam.AlignmentFile(bam_fname)
    pu = inf.pileup('chr1')
    for p in pu:
        pid = p.reference_id
        ppos = p.reference_pos

    if TYPE_CHECKING: reveal_locals()
    types = typecheck(check_locals=False)
    assert re.search(r'\.IteratorColumn(|All|AllRefs|Region)$', types['pu'])
    assert types['p'].endswith('PileupColumn')
    assert types['pid'] == 'int'
    assert types['ppos'] == 'int'


def test_samtools_subcommands() -> None:
    p1_samtools_faidx = pysam.samtools.faidx
    p2_faidx = pysam.faidx
    p3_view = pysam.view
    p4_bcftools_view = pysam.bcftools.view

    if TYPE_CHECKING: reveal_locals()
    types = typecheck()
    for var, vartype in types.items():
        assert vartype.endswith('PysamDispatcher'), f'{var!r} is not a dispatcher'

def test_AlignmentFile_filenames(sam_fname: str) -> None:
    fp1 = pysam.AlignmentFile(sam_fname)
    fp2 = pysam.AlignmentFile(sam_fname.encode())
    fp3 = pysam.AlignmentFile(pathlib.Path(sam_fname))

    fd = os.open(sam_fname, os.O_RDONLY)
    fp4 = pysam.AlignmentFile(fd)
    os.close(fd)

    class WrappedFd:
        def __init__(self, fname): self.fd = os.open(fname, os.O_RDONLY)
        def fileno(self): return self.fd
        def __del__(self): os.close(self.fd)

    fp5 = pysam.AlignmentFile(WrappedFd(sam_fname))

    fio = open(sam_fname)
    fp6 = pysam.AlignmentFile(fio)

    if TYPE_CHECKING: reveal_locals()
    types = typecheck()

def test_VariantFile_filenames(vcf_fname: str) -> None:
    fp1 = pysam.VariantFile(vcf_fname)
    fp2 = pysam.VariantFile(vcf_fname.encode())
    fp3 = pysam.VariantFile(pathlib.Path(vcf_fname))

    fd = os.open(vcf_fname, os.O_RDONLY)
    fp4 = pysam.VariantFile(fd)
    os.close(fd)

    class WrappedFd:
        def __init__(self, fname): self.fd = os.open(fname, os.O_RDONLY)
        def fileno(self): return self.fd
        def __del__(self): os.close(self.fd)

    fp5 = pysam.VariantFile(WrappedFd(vcf_fname))

    fio = open(vcf_fname)
    fp6 = pysam.VariantFile(fio)

    if TYPE_CHECKING: reveal_locals()
    types = typecheck()
