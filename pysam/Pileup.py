'''Tools for working with files in the samtools pileup -c format.'''
import collections
import pysam
import pysam.VCF

PileupSubstitution = collections.namedtuple( "PileupSubstitution",
                                    " ".join( (\
            "chromosome", 
            "pos", 
            "reference_base", 
            "genotype",
            "consensus_quality",
            "snp_quality",
            "mapping_quality",
            "coverage",
            "read_bases",
            "base_qualities" ) ) )

PileupIndel = collections.namedtuple( "PileupIndel",
                                      " ".join( (\
            "chromosome", 
            "pos", 
            "reference_base", 
            "genotype",
            "consensus_quality",
            "snp_quality",
            "mapping_quality",
            "coverage",
            "first_allele",
            "second_allele",
            "reads_first",
            "reads_second",
            "reads_diff" ) ) )

def iterate( infile ):
    '''iterate over ``samtools pileup -c`` formatted file.

    *infile* can be any iterator over a lines.

    The function yields named tuples of the type :class:`pysam.Pileup.PileupSubstitution`
    or :class:`pysam.Pileup.PileupIndel`.

    .. note:: 
       The parser converts to 0-based coordinates
    '''
    
    conv_subst = (str,lambda x: int(x)-1,str,str,int,int,int,int,str,str)
    conv_indel = (str,lambda x: int(x)-1,str,str,int,int,int,int,str,str,int,int,int)

    for line in infile:
        d = line[:-1].split()
        if d[2] == "*":
            try:
                yield PileupIndel( *[x(y) for x,y in zip(conv_indel,d) ] )
            except TypeError:
                raise pysam.SamtoolsError( "parsing error in line: `%s`" % line)
        else:
            try:
                yield PileupSubstitution( *[x(y) for x,y in zip(conv_subst,d) ] )
            except TypeError:
                raise pysam.SamtoolsError( "parsing error in line: `%s`" % line)

ENCODE_GENOTYPE = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'AA': 'A', 'CC': 'C', 'GG': 'G', 'TT': 'T', 'UU': 'U',
    'AG': 'r', 'GA': 'R',
    'CT': 'y', 'TC': 'Y',
    'AC': 'm', 'CA': 'M',
    'GT': 'k', 'TG': 'K',
    'CG': 's', 'GC': 'S',
    'AT': 'w', 'TA': 'W',
    }        

DECODE_GENOTYPE = {
    'A': 'AA',
    'C': 'CC',
    'G': 'GG',
    'T': 'TT',
    'r': 'AG', 'R': 'AG',
    'y': 'CT', 'Y': 'CT',
    'm': 'AC', 'M': 'AC',
    'k': 'GT', 'K': 'GT',
    's': 'CG', 'S': 'CG',
    'w': 'AT', 'W': 'AT',
    }

##------------------------------------------------------------
def encodeGenotype( code ):
    '''encode genotypes like GG, GA into a one-letter code.
    The returned code is lower case if code[0] < code[1], otherwise
    it is uppercase.
    '''
    return ENCODE_GENOTYPE[ code.upper() ]

def decodeGenotype( code ):
    '''decode single letter genotypes like m, M into two letters.
    This is the reverse operation to :meth:`encodeGenotype`.
    '''
    return DECODE_GENOTYPE[ code ] 


def iterate_from_vcf( infile, sample ):
    '''iterate over a vcf-formatted file.

    *infile* can be any iterator over a lines.

    The function yields named tuples of the type :class:`pysam.Pileup.PileupSubstitution`
    or :class:`pysam.Pileup.PileupIndel`.

    Positions without a snp will be skipped.

    .. note::
        indels are not implemented yet
    '''

    
    vcf = pysam.VCF.VCFFile()
    
    iterator = vcf.parse( infile )

    for data in iterator:
        
        if sample not in data: continue
        v = data[sample]
        chromosome = data["chrom"]
        pos = data["pos"]
        reference_base = data["ref"]
        bases = [reference_base] + data["alt"]

        if len(reference_base) > 1:
            # indels
            # raise NotImplementedError( "indels not implemented" )
            pass
        elif max([len(x) for x in data["alt"]] ) > 1:
            # indels
            # raise NotImplementedError( "indels not implemented" )
            pass
        else:
            # get genotype
            genotypes = v["GT"]
            if len(genotypes) > 1:
                raise ValueError( "only single genotype per position, %s" % (str(data)))
            genotypes = [ x for x in genotypes[0] if x != "/"]

            # not a snp
            if genotypes[0] == ".": continue

            genotype = encodeGenotype( "".join([ bases[int(x)] for x in genotypes if x != "/" ]) )
            
            snp_quality = 0
            consensus_quality = v.get( "GQ", [0])[0]
            mapping_quality = data["info"].get( "MQ", [0])[0]
            coverage = v.get( "DP", 0)
            read_bases = ""
            base_qualities = ""

            yield PileupSubstitution( chromosome, pos, reference_base, genotype,
                                      consensus_quality, 
                                      snp_quality, 
                                      mapping_quality,
                                      coverage, read_bases, base_qualities ) 

    
