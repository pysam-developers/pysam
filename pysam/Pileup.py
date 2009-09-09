'''Tools for working with files in the samtools pileup -c format.'''
import collections
import pysam

PileupEntry = collections.namedtuple( "PileupEntry",
                                      " ".join( (\
            "chromosome", 
            "position", 
            "reference_base", 
            "consensus_base",
            "consensus_quality",
            "snp_quality",
            "rms_mapping_quality",
            "coverage",
            "read_bases",
            "base_qualities" ) ) )

def iterate( infile ):
    '''iterate over ``samtools pileup -c`` formatted file.

    *infile* can be any iterator over a lines.

    The function yields objects of the type :class:`pysam.Pileup.PileupEntry`.
    '''
    
    conv = (str,int,str,str,int,int,int,int,str,str)

    for line in infile:
        try:
            yield PileupEntry( *[x(y) for x,y in zip(conv,line[:-1].split()) ] )
        except TypeError:
            raise pysam.SamtoolsError( "parsing error in line: `%s`" % line)
