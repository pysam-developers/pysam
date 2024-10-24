import pysam.utils

annotate = pysam.utils.PysamDispatcher('bcftools', 'annotate')
call = pysam.utils.PysamDispatcher('bcftools', 'call')
cnv = pysam.utils.PysamDispatcher('bcftools', 'cnv')
concat = pysam.utils.PysamDispatcher('bcftools', 'concat')
consensus = pysam.utils.PysamDispatcher('bcftools', 'consensus')
convert = pysam.utils.PysamDispatcher('bcftools', 'convert')
csq = pysam.utils.PysamDispatcher('bcftools', 'csq')
filter = pysam.utils.PysamDispatcher('bcftools', 'filter')
gtcheck = pysam.utils.PysamDispatcher('bcftools', 'gtcheck')
head = pysam.utils.PysamDispatcher('bcftools', 'head')
index = pysam.utils.PysamDispatcher('bcftools', 'index')
isec = pysam.utils.PysamDispatcher('bcftools', 'isec')
merge = pysam.utils.PysamDispatcher('bcftools', 'merge')
mpileup = pysam.utils.PysamDispatcher('bcftools', 'mpileup')
norm = pysam.utils.PysamDispatcher('bcftools', 'norm')
plugin = pysam.utils.PysamDispatcher('bcftools', 'plugin')
query = pysam.utils.PysamDispatcher('bcftools', 'query')
reheader = pysam.utils.PysamDispatcher('bcftools', 'reheader')
roh = pysam.utils.PysamDispatcher('bcftools', 'roh')
sort = pysam.utils.PysamDispatcher('bcftools', 'sort')
stats = pysam.utils.PysamDispatcher('bcftools', 'stats')
view = pysam.utils.PysamDispatcher('bcftools', 'view')

__all__ = [
    'annotate', 'call', 'cnv', 'concat', 'consensus',
    'convert', 'csq', 'filter', 'gtcheck', 'head',
    'index', 'isec', 'merge', 'mpileup', 'norm',
    'plugin', 'query', 'reheader', 'roh', 'sort',
    'stats', 'view',
]
