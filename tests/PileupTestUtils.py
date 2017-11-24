import os
import subprocess
import pysam

from TestUtils import BAM_DATADIR, force_str

def build_pileup_with_samtoolsshell(fn):
    os.system("samtools mpileup {} 2> /dev/null | wc -l > /dev/null".format(fn))
    return 2998


def build_pileup_with_samtoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen(["samtools", "mpileup", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        return len(proc.stdout.readlines())


def build_pileup_with_pysam(*args, **kwargs):
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        return len(list(inf.pileup(stepper="samtools")))


def build_depth_with_samtoolsshell(fn):
    os.system(
        "samtools mpileup {} 2> /dev/null | awk '{{a += $4}} END {{print a}}' > /dev/null".format(fn))
    return 107241


def build_depth_with_samtoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen(["samtools", "mpileup", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        data = [x.split() for x in proc.stdout.readlines()]
        return [int(x[3]) for x in data]


def build_depth_with_filter_with_pysam(*args, **kwargs):
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        return [x.get_num_aligned() for x in inf.pileup(stepper="samtools")]


def build_depth_with_pysam(*args, **kwargs):
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        return [x.nsegments for x in inf.pileup(stepper="samtools")]


def build_query_bases_with_samtoolsshell(fn):
    os.system("samtools mpileup {} 2> /dev/null | awk '{{a = a $5}} END {{print a}}' | wc -c > /dev/null".format(fn))
    return 116308


def build_query_bases_with_samtoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen(["samtools", "mpileup", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        stdout = proc.stdout.read().decode()
        return [x.split()[4] for x in stdout.splitlines()]


def build_query_bases_with_samtoolspysam(fn):
    return [x.split()[4] for x in pysam.samtools.mpileup(fn).splitlines()]
    

def build_query_bases_with_pysam_pileups(*args, **kwargs):
    total_pileup = []
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        total_pileup = [
            [r.alignment.query_sequence[r.query_position_or_next]
                for r in column.pileups if r.query_position_or_next is not None]
            for column in inf.pileup(stepper="samtools")]
    return total_pileup


def build_query_qualities_with_pysam_pileups(*args, **kwargs):
    total_pileup = []
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        total_pileup = [
            [r.alignment.query_qualities[r.query_position_or_next]
                for r in column.pileups if r.query_position_or_next is not None]
            for column in inf.pileup(stepper="samtools")]
    return total_pileup


def build_query_bases_with_pysam(*args, **kwargs):
    total_pileup = []
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        total_pileup = [column.get_query_sequences(add_markers=True, add_indels=True) for column in
                        inf.pileup(stepper="samtools")]
    return total_pileup


def build_query_names_with_pysam(*args, **kwargs):
    total_pileup = []
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        total_pileup = [column.get_query_names() for column in
                        inf.pileup(stepper="samtools")]
    return total_pileup


def build_query_qualities_with_samtoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen(["samtools", "mpileup", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        data = [force_str(x).split()[5] for x in proc.stdout.readlines()]
    return data


def build_query_qualities_with_pysam(*args, **kwargs):
    total_pileup = []
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        total_pileup = [column.get_query_qualities() for column in
                        inf.pileup(stepper="samtools")]
    return total_pileup


def build_mapping_qualities_with_samtoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen(["samtools", "mpileup", "-s", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        data = [force_str(x).split()[6] for x in proc.stdout.readlines()]
    return data


def build_mapping_qualities_with_pysam(*args, **kwargs):
    total_pileup = []
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        total_pileup = [column.get_mapping_qualities() for column in
                        inf.pileup(stepper="samtools")]
    return total_pileup


def build_query_positions_with_samtoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen(["samtools", "mpileup", "-O", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        data = [list(map(int, force_str(x).split()[6].split(",")))
                for x in proc.stdout.readlines()]
    return data


def build_query_positions_with_pysam(*args, **kwargs):
    total_pileup = []
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        total_pileup = [column.get_query_positions() for column in
                        inf.pileup(stepper="samtools")]
    return total_pileup
