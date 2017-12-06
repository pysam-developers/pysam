import os
import subprocess
import pysam

from TestUtils import BAM_DATADIR, force_str

def build_fetch_with_samtoolsshell(fn):
    retval = os.popen("samtools view {} 2> /dev/null | wc -l".format(fn)).read()
    return int(retval.strip())


def build_fetch_with_samtoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen(["samtools", "view", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        return len(proc.stdout.readlines())

    
def build_fetch_with_pysam(*args, **kwargs):
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        return len(list(inf.fetch()))


def build_query_sequences_with_samtoolsshell(fn):
    retval = os.popen("samtools view {} 2> /dev/null | cut -f 11".format(fn)).read()
    return force_str(retval).splitlines()


def build_query_sequences_with_samtoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen(["samtools", "view", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        data = [force_str(x).split()[10] for x in proc.stdout.readlines()]
    return data

    
def build_query_sequences_with_pysam(*args, **kwargs):
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        data = [x.query_sequence for x in inf]
    return data


def build_query_qualities_with_pysam(*args, **kwargs):
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        data = [x.query_qualities for x in inf]
    return data


def build_query_sequences_flagfilter_with_samtoolsshell(fn):
    retval = os.popen("samtools view -f 2 {} 2> /dev/null | cut -f 11".format(fn)).read()
    return force_str(retval).splitlines()


def build_query_sequences_flagfilter_with_samtoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen(["samtools", "view", "-f", "2", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        data = [force_str(x).split()[10] for x in proc.stdout.readlines()]
    return data

def build_query_sequences_flagfilter_with_pysam(*args, **kwargs):
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        data = [x.query_sequence for x in inf if x.is_proper_pair]
    return data


def build_query_sequences_directflagfilter_with_pysam(*args, **kwargs):
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        data = [x.query_sequence for x in inf if x.flag & 2]
    return data


def build_aligned_pairs_with_pysam(*args, **kwargs):
    matches_only = kwargs.pop("matches_only", False)
    with_seq = kwargs.pop("with_seq", False)
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        data = [x.get_aligned_pairs(matches_only=matches_only, with_seq=with_seq)
                for x in inf if not x.is_unmapped]
    return data
    
