'''Tools for working with files in the samtools pileup -c format.'''
import collections
import pysam

PileupSubstitution = collections.namedtuple("PileupSubstitution",
                                            " ".join((
                                                "chromosome",
                                                "pos",
                                                "reference_base",
                                                "genotype",
                                                "consensus_quality",
                                                "snp_quality",
                                                "mapping_quality",
                                                "coverage",
                                                "read_bases",
                                                "base_qualities")))

PileupIndel = collections.namedtuple("PileupIndel",
                                     " ".join((
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
                                         "reads_diff")))


def iterate(infile):
    '''iterate over ``samtools pileup -c`` formatted file.

    *infile* can be any iterator over a lines.

    The function yields named tuples of the type :class:`pysam.Pileup.PileupSubstitution`
    or :class:`pysam.Pileup.PileupIndel`.

    .. note::

       The parser converts to 0-based coordinates
    '''

    conv_subst = (str, lambda x: int(x) - 1, str,
                  str, int, int, int, int, str, str)
    conv_indel = (str, lambda x: int(x) - 1, str, str, int,
                  int, int, int, str, str, int, int, int)

    for line in infile:
        d = line[:-1].split()
        if d[2] == "*":
            try:
                yield PileupIndel(*[x(y) for x, y in zip(conv_indel, d)])
            except TypeError:
                raise pysam.SamtoolsError("parsing error in line: `%s`" % line)
        else:
            try:
                yield PileupSubstitution(*[x(y) for x, y in zip(conv_subst, d)])
            except TypeError:
                raise pysam.SamtoolsError("parsing error in line: `%s`" % line)


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

# ------------------------------------------------------------


def encodeGenotype(code):
    '''encode genotypes like GG, GA into a one-letter code.
    The returned code is lower case if code[0] < code[1], otherwise
    it is uppercase.
    '''
    return ENCODE_GENOTYPE[code.upper()]


def decodeGenotype(code):
    '''decode single letter genotypes like m, M into two letters.
    This is the reverse operation to :meth:`encodeGenotype`.
    '''
    return DECODE_GENOTYPE[code]


def translateIndelGenotypeFromVCF(vcf_genotypes, ref):
    '''translate indel from vcf to pileup format.'''

    # indels
    def getPrefix(s1, s2):
        '''get common prefix of strings s1 and s2.'''
        n = min(len(s1), len(s2))
        for x in range(n):
            if s1[x] != s2[x]:
                return s1[:x]
        return s1[:n]

    def getSuffix(s1, s2):
        '''get common sufix of strings s1 and s2.'''
        n = min(len(s1), len(s2))
        if s1[-1] != s2[-1]:
            return ""
        for x in range(-2, -n - 1, -1):
            if s1[x] != s2[x]:
                return s1[x + 1:]
        return s1[-n:]

    def getGenotype(variant, ref):

        if variant == ref:
            return "*", 0

        if len(ref) > len(variant):
            # is a deletion
            if ref.startswith(variant):
                return "-%s" % ref[len(variant):], len(variant) - 1
            elif ref.endswith(variant):
                return "-%s" % ref[:-len(variant)], -1
            else:
                prefix = getPrefix(ref, variant)
                suffix = getSuffix(ref, variant)
                shared = len(prefix) + len(suffix) - len(variant)
                # print "-", prefix, suffix, ref, variant, shared, len(prefix), len(suffix), len(ref)
                if shared < 0:
                    raise ValueError()
                return "-%s" % ref[len(prefix):-(len(suffix) - shared)], len(prefix) - 1

        elif len(ref) < len(variant):
            # is an insertion
            if variant.startswith(ref):
                return "+%s" % variant[len(ref):], len(ref) - 1
            elif variant.endswith(ref):
                return "+%s" % variant[:len(ref)], 0
            else:
                prefix = getPrefix(ref, variant)
                suffix = getSuffix(ref, variant)
                shared = len(prefix) + len(suffix) - len(ref)
                if shared < 0:
                    raise ValueError()

                return "+%s" % variant[len(prefix):-(len(suffix) - shared)], len(prefix)
        else:
            assert 0, "snp?"

    # in pileup, the position refers to the base
    # after the coordinate, hence subtract 1
    # pos -= 1

    genotypes, offsets = [], []
    is_error = True

    for variant in vcf_genotypes:
        try:
            g, offset = getGenotype(variant, ref)
        except ValueError:
            break

        genotypes.append(g)
        if g != "*":
            offsets.append(offset)

    else:
        is_error = False

    if is_error:
        raise ValueError()

    assert len(set(offsets)) == 1, "multiple offsets for indel"
    offset = offsets[0]

    genotypes = "/".join(genotypes)
    return genotypes, offset


def vcf2pileup(vcf, sample):
    '''convert vcf record to pileup record.'''

    chromosome = vcf.contig
    pos = vcf.pos
    reference = vcf.ref
    allelles = [reference] + vcf.alt

    data = vcf[sample]

    # get genotype
    genotypes = data["GT"]
    if len(genotypes) > 1:
        raise ValueError("only single genotype per position, %s" % (str(vcf)))

    genotypes = genotypes[0]

    # not a variant
    if genotypes[0] == ".":
        return None

    genotypes = [allelles[int(x)] for x in genotypes if x != "/"]

    # snp_quality is "genotype quality"
    snp_quality = consensus_quality = data.get("GQ", [0])[0]
    mapping_quality = vcf.info.get("MQ", [0])[0]
    coverage = data.get("DP", 0)

    if len(reference) > 1 or max([len(x) for x in vcf.alt]) > 1:
        # indel
        genotype, offset = translateIndelGenotypeFromVCF(genotypes, reference)

        return PileupIndel(chromosome,
                           pos + offset,
                           "*",
                           genotype,
                           consensus_quality,
                           snp_quality,
                           mapping_quality,
                           coverage,
                           genotype,
                           "<" * len(genotype),
                           0,
                           0,
                           0)

    else:
        genotype = encodeGenotype("".join(genotypes))
        read_bases = ""
        base_qualities = ""

        return PileupSubstitution(chromosome, pos, reference,
                                  genotype, consensus_quality,
                                  snp_quality, mapping_quality,
                                  coverage, read_bases,
                                  base_qualities)


def iterate_from_vcf(infile, sample):
    '''iterate over a vcf-formatted file.

    *infile* can be any iterator over a lines.

    The function yields named tuples of the type
    :class:`pysam.Pileup.PileupSubstitution` or
    :class:`pysam.Pileup.PileupIndel`.

    Positions without a snp will be skipped.

    This method is wasteful and written to support same legacy code
    that expects samtools pileup output.

    Better use the vcf parser directly.

    '''
    vcf = pysam.VCF()
    vcf.connect(infile)

    if sample not in vcf.getsamples():
        raise KeyError("sample %s not vcf file")

    for row in vcf.fetch():
        result = vcf2pileup(row, sample)
        if result:
            yield result
