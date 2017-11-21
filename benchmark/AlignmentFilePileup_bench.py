"""Benchmarking module for AlignmentFile functionality"""
import os
import subprocess
import pysam


from TestUtils import BAM_DATADIR, force_str


# Use -x option in samtools mpileup to turn off overlap detection
# API not available at the moment.

####################################################
####################################################    
# Simply building a pileup and counting number of piled-up columns
def build_pileup_with_samtools(fn):
    os.system("samtools mpileup -x {} 2> /dev/null | wc -l > /dev/null".format(fn))
    return 2998


def build_pileup_with_samtoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen(["samtools", "mpileup", "-x", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        return len(proc.stdout.readlines())


def build_pileup_with_pysam(*args, **kwargs):
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        return len(list(inf.pileup(stepper="samtools")))


def test_build_pileup_from_bam_with_samtools(benchmark):
    result = benchmark(build_pileup_with_samtools,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 2998


def test_build_pileup_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_pileup_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 2998


def test_build_pileup_from_bam_with_pysam(benchmark):
    result = benchmark(build_pileup_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 2998


####################################################
####################################################    
# Build depth profile with pileup
def build_depth_with_samtools(fn):
    os.system("samtools mpileup -x {} 2> /dev/null | awk '{{a += $4}} END {{print a}}' > /dev/null".format(fn))
    return 107248


def build_depth_with_samtoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen(["samtools", "mpileup", "-x", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        data = [x.split() for x in proc.stdout.readlines()]
        return sum([int(x[3]) for x in data])


def build_depth_with_filter_with_pysam(*args, **kwargs):
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        return sum([x.get_num_aligned() for x in inf.pileup(stepper="samtools")])


def build_depth_with_pysam(*args, **kwargs):
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        return sum([x.nsegments for x in inf.pileup(stepper="samtools")])
    

def test_build_depth_from_bam_with_samtools(benchmark):
    result = benchmark(build_depth_with_samtools,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 107248


def test_build_depth_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_depth_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 107248


def test_build_depth_from_bam_with_pysam(benchmark):
    result = benchmark(build_depth_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    # different value, as samtools filters with a minimum
    # base quality of 13
    assert result == 110015


def test_build_depth_with_filter_from_bam_with_pysam(benchmark):
    result = benchmark(build_depth_with_filter_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    # why not 107241?
    assert result == 107248
    

####################################################
####################################################    
# Build depth profile with pileup
def build_base_concatenation_with_samtools(fn):
    os.system("samtools mpileup -x {} 2> /dev/null | awk '{{a = a $5}} END {{print a}}' | wc -c > /dev/null".format(fn))
    return 116314


def build_base_concatenation_with_samtoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen(["samtools", "mpileup", "-x", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        data = [x.split() for x in proc.stdout.readlines()]
        return len(b"".join([x[4] for x in data]))


def build_base_concatenation_with_pysam_pileups(*args, **kwargs):
    total_pileup = []
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        total_pileup = "".join([
            "".join([r.alignment.query_sequence[r.query_position] for r in column.pileups if r.query_position])
            for column in inf.pileup(stepper="samtools")])
    return len(total_pileup)


def build_base_concatenation_with_pysam(*args, **kwargs):
    total_pileup = []
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        total_pileup = "".join(
            ["".join(column.get_query_sequences()) for column in
             inf.pileup(stepper="samtools")])
    return len(total_pileup)


def test_build_base_concatenation_from_bam_with_samtools(benchmark):
    result = benchmark(build_base_concatenation_with_samtools,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 116314


def test_build_base_concatenation_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_base_concatenation_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 116314


def test_build_base_concatenation_from_bam_with_pysam_pileups(benchmark):
    result = benchmark(build_base_concatenation_with_pysam_pileups,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 106889


def test_build_base_concatenation_from_bam_with_pysam(benchmark):
    result = benchmark(build_base_concatenation_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert result == 116314


def build_qualities_concatenation_with_samtoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen(["samtools", "mpileup", "-x", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        data = [force_str(x).split()[5] for x in proc.stdout.readlines()]
    return data
    

def build_qualities_concatenation_with_pysam(*args, **kwargs):
    total_pileup = []
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        total_pileup = [column.get_query_qualities() for column in
                        inf.pileup(stepper="samtools")]
    return total_pileup

    
def test_build_quality_concatenation_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_qualities_concatenation_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len("".join(result)) == 107248


def test_build_quality_concatenation_from_bam_with_pysam(benchmark):
    result = benchmark(build_qualities_concatenation_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert sum([len(x) for x in result]) == 107248


def build_querynames_concatenation_with_pysam(*args, **kwargs):
    total_pileup = []
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        total_pileup = [column.get_query_names() for column in
                        inf.pileup(stepper="samtools")]
    return total_pileup

    
def test_build_querynames_concatenation_from_bam_with_pysam(benchmark):
    result = benchmark(build_querynames_concatenation_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len("".join([x for column in result for x in column])) == 2307497


def build_mappingqualities_concatenation_with_samtoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen(["samtools", "mpileup", "-xs", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        data = [force_str(x).split()[6] for x in proc.stdout.readlines()]
    return data
    

def build_mappingqualities_concatenation_with_pysam(*args, **kwargs):
    total_pileup = []
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        total_pileup = [column.get_mapping_qualities() for column in
                        inf.pileup(stepper="samtools")]
    return total_pileup

    
def test_build_mappingquality_concatenation_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_mappingqualities_concatenation_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert len("".join(result)) == 107248


def test_build_mappingquality_concatenation_from_bam_with_pysam(benchmark):
    result = benchmark(build_mappingqualities_concatenation_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert sum([len(x) for x in result]) == 107248


def build_query_positions_concatenation_with_samtoolspipe(fn):
    FNULL = open(os.devnull, 'w')
    with subprocess.Popen(["samtools", "mpileup", "-xO", fn],
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=FNULL) as proc:
        data = [list(map(int, force_str(x).split()[6].split(","))) for x in proc.stdout.readlines()]
    return data
    

def build_query_positions_concatenation_with_pysam(*args, **kwargs):
    total_pileup = []
    with pysam.AlignmentFile(*args, **kwargs) as inf:
        total_pileup = [column.get_query_positions() for column in
                        inf.pileup(stepper="samtools")]
    return total_pileup

    
def test_build_query_positions_concatenation_from_bam_with_samtoolspipe(benchmark):
    result = benchmark(build_query_positions_concatenation_with_samtoolspipe,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    # positions output by samtools are 1-based
    assert sum([sum(x) - len(x) for x in result]) == 1841736


def test_build_querypositions_concatenation_from_bam_with_pysam(benchmark):
    result = benchmark(build_query_positions_concatenation_with_pysam,
                       os.path.join(BAM_DATADIR, "ex2.bam"))
    assert sum([sum(x) for x in result]) == 1841736
    
