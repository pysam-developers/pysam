import pp
import pysam

def sortBam(bam_filepath, out_prefix):
    pysam.sort(bam_filepath, out_prefix)

job_server = pp.Server(ncpus=2, ppservers = (), secret='secret')

job1 = job_server.submit(sortBam, ('ex1.bam', 'tmp1_'), modules=('pysam',))
job2 = job_server.submit(sortBam, ('ex2.bam', 'tmp2_'), modules=('pysam',))

result1 = job1()
result2 = job2()
