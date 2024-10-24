import pysam.utils

addreplacerg = pysam.utils.PysamDispatcher('samtools', 'addreplacerg')
ampliconclip = pysam.utils.PysamDispatcher('samtools', 'ampliconclip')
ampliconstats = pysam.utils.PysamDispatcher('samtools', 'ampliconstats')
bam2fq = pysam.utils.PysamDispatcher('samtools', 'bam2fq')
bamshuf = pysam.utils.PysamDispatcher('samtools', 'bamshuf')
bedcov = pysam.utils.PysamDispatcher('samtools', 'bedcov')
calmd = pysam.utils.PysamDispatcher('samtools', 'calmd')
cat = pysam.utils.PysamDispatcher('samtools', 'cat')
collate = pysam.utils.PysamDispatcher('samtools', 'collate')
consensus = pysam.utils.PysamDispatcher('samtools', 'consensus')
coverage = pysam.utils.PysamDispatcher('samtools', 'coverage')
cram_size = pysam.utils.PysamDispatcher('samtools', 'cram-size')
depad = pysam.utils.PysamDispatcher('samtools', 'depad')
depth = pysam.utils.PysamDispatcher('samtools', 'depth')
dict = pysam.utils.PysamDispatcher('samtools', 'dict')
faidx = pysam.utils.PysamDispatcher('samtools', 'faidx')
fasta = pysam.utils.PysamDispatcher('samtools', 'fasta')
fastq = pysam.utils.PysamDispatcher('samtools', 'fastq')
fixmate = pysam.utils.PysamDispatcher('samtools', 'fixmate')
flags = pysam.utils.PysamDispatcher('samtools', 'flags')
flagstat = pysam.utils.PysamDispatcher('samtools', 'flagstat')
fqidx = pysam.utils.PysamDispatcher('samtools', 'fqidx')
fqimport = pysam.utils.PysamDispatcher('samtools', 'import')
head = pysam.utils.PysamDispatcher('samtools', 'head')
idxstats = pysam.utils.PysamDispatcher('samtools', 'idxstats')
index = pysam.utils.PysamDispatcher('samtools', 'index')
markdup = pysam.utils.PysamDispatcher('samtools', 'markdup')
merge = pysam.utils.PysamDispatcher('samtools', 'merge')
mpileup = pysam.utils.PysamDispatcher('samtools', 'mpileup')
pad2unpad = pysam.utils.PysamDispatcher('samtools', 'pad2unpad')
phase = pysam.utils.PysamDispatcher('samtools', 'phase')
quickcheck = pysam.utils.PysamDispatcher('samtools', 'quickcheck')
reference = pysam.utils.PysamDispatcher('samtools', 'reference')
reheader = pysam.utils.PysamDispatcher('samtools', 'reheader')
reset = pysam.utils.PysamDispatcher('samtools', 'reset')
rmdup = pysam.utils.PysamDispatcher('samtools', 'rmdup')
samples = pysam.utils.PysamDispatcher('samtools', 'samples')
sort = pysam.utils.PysamDispatcher('samtools', 'sort')
split = pysam.utils.PysamDispatcher('samtools', 'split')
stats = pysam.utils.PysamDispatcher('samtools', 'stats')
targetcut = pysam.utils.PysamDispatcher('samtools', 'targetcut')
tview = pysam.utils.PysamDispatcher('samtools', 'tview')
version = pysam.utils.PysamDispatcher('samtools', 'version')
view = pysam.utils.PysamDispatcher('samtools', 'view')

__all__ = [
    'addreplacerg', 'ampliconclip', 'ampliconstats',
    'bam2fq', 'bamshuf', 'bedcov', 'calmd', 'cat',
    'collate', 'consensus', 'coverage', 'cram_size',
    'depad', 'depth', 'dict', 'faidx', 'fasta',
    'fastq', 'fixmate', 'flags', 'flagstat', 'fqidx',
    'fqimport', 'head', 'idxstats', 'index',
    'markdup', 'merge', 'mpileup', 'pad2unpad',
    'phase', 'quickcheck', 'reference', 'reheader',
    'reset', 'rmdup', 'samples', 'sort', 'split',
    'stats', 'targetcut', 'tview', 'version', 'view',
]
