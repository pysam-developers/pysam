import os
import subprocess
import sys

import pytest

import pysam.samtools  # noqa: F401  (dispatcher must be importable)

from TestUtils import BAM_DATADIR

# The samtools/bcftools wrappers run in-process and read file descriptor 0, so
# a terminal on stdin leaves them waiting for input that never arrives.  Each
# test runs the dispatch in a child process with a controlled stdin and fails
# if it does not return.
TIMEOUT = 60

IGNORE_ERROR = ("import pysam\n"
                "try:\n"
                "    pysam.samtools.view('-')\n"
                "except pysam.SamtoolsError:\n"
                "    pass\n")


def run(code, stdin, data=None):
    '''run *code* with *stdin*, feeding it *data*, and return its stdout.'''
    proc = subprocess.Popen([sys.executable, "-c", code], stdin=stdin,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        out, err = proc.communicate(input=data, timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.communicate()
        raise AssertionError("dispatch did not return within %i seconds" %
                             TIMEOUT)
    assert proc.returncode == 0, err.decode()
    return out


@pytest.mark.skipif(not hasattr(os, "openpty"),
                    reason="requires a pseudo terminal")
class TestStdinDoesNotBlock:
    '''samtools and bcftools must not wait forever for stdin.'''

    def test_stdin_from_terminal_returns(self):
        '''a terminal as stdin must not make the dispatch hang'''
        master, slave = os.openpty()
        try:
            run(IGNORE_ERROR, slave)
        finally:
            os.close(master)
            os.close(slave)

    def test_stdin_from_pipe_still_reads(self):
        '''a pipe carrying real input must keep being read'''
        with open(os.path.join(BAM_DATADIR, "ex_spliced.sam"), "rb") as inf:
            data = inf.read()
        n_records = sum(1 for line in data.split(b"\n")
                        if line and not line.startswith(b"@"))
        out = run("import sys, pysam.samtools\n"
                  "sys.stdout.write("
                  "str(len(pysam.samtools.view('-').splitlines())))\n",
                  subprocess.PIPE, data)
        assert int(out) == n_records, "streamed input was not passed through"
