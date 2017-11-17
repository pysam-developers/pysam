/* The MIT License

   Copyright (c) 2016 Genome Research Ltd.

   Author: Petr Danecek <pd3@sanger.ac.uk>
   
   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:
   
   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.
   
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

 */
/*
    Things that would be nice to have
        - for stop-lost events (also in frameshifts) report the number of truncated aa's
        - memory could be greatly reduced by indexing gff (but it is quite compact already)
        - deletions that go beyond transcript boundaries are not checked at sequence level
            - alloc tscript->ref in hap_finalize, introduce fa_off_beg:16,fa_off_end:16
            - see test/csq/ENST00000573314/insertion-overlap.vcf #1476288882

    Read about transcript types here
        http://vega.sanger.ac.uk/info/about/gene_and_transcript_types.html
        http://www.ensembl.org/info/genome/variation/predicted_data.html
        http://www.gencodegenes.org/gencode_biotypes.html

    List of supported biotypes
        antisense
        IG_C_gene
        IG_D_gene
        IG_J_gene
        IG_LV_gene
        IG_V_gene
        lincRNA
        macro_lncRNA
        miRNA
        misc_RNA
        Mt_rRNA
        Mt_tRNA
        polymorphic_pseudogene
        processed_transcript
        protein_coding
        ribozyme
        rRNA
        sRNA
        scRNA
        scaRNA
        sense_intronic
        sense_overlapping
        snRNA
        snoRNA
        TR_C_gene
        TR_D_gene
        TR_J_gene
        TR_V_gene

    The gff parsing logic
        We collect features such by combining gff lines A,B,C as follows:
            A .. gene line with a supported biotype
                    A.ID=~/^gene:/

            B .. transcript line referencing A
                    B.ID=~/^transcript:/ && B.Parent=~/^gene:A.ID/

            C .. corresponding CDS, exon, and UTR lines:
                    C[3] in {"CDS","exon","three_prime_UTR","five_prime_UTR"} && C.Parent=~/^transcript:B.ID/ 

        For coding biotypes ("protein_coding" or "polymorphic_pseudogene") the
        complete chain link C -> B -> A is required. For the rest, link B -> A suffices.
        
                
    The supported consequence types, sorted by impact:
        splice_acceptor_variant .. end region of an intron changed (2bp at the 3' end of an intron)
        splice_donor_variant    .. start region of an intron changed (2bp at the 5' end of an intron)
        stop_gained             .. DNA sequence variant resulting in a stop codon
        frameshift_variant      .. number of inserted/deleted bases not a multiple of three, disrupted translational frame
        stop_lost               .. elongated transcript, stop codon changed
        start_lost              .. the first codon changed
        inframe_altering        .. combination of indels leading to unchanged reading frame and length
        inframe_insertion       .. inserted coding sequence, unchanged reading frame
        inframe_deletion        .. deleted coding sequence, unchanged reading frame
        missense_variant        .. amino acid (aa) change, unchanged length
        splice_region_variant   .. change within 1-3 bases of the exon or 3-8 bases of the intron
        synonymous_variant      .. DNA sequence variant resulting in no amino acid change
        stop_retained_variant   .. different stop codon
        non_coding_variant      .. variant in non-coding sequence, such as RNA gene
        5_prime_UTR_variant
        3_prime_UTR_variant
        intron_variant          .. reported only if none of the above
        intergenic_variant      .. reported only if none of the above


    The annotation algorithm.
        The algorithm checks if the variant falls in a region of a supported type. The
        search is performed in the following order, until a match is found:
            1. idx_cds(gf_cds_t) - lookup CDS by position, create haplotypes, call consequences
            2. idx_utr(gf_utr_t) - check UTR hits
            3. idx_exon(gf_exon_t) - check for splice variants
            4. idx_tscript(tscript_t) - check for intronic variants, RNAs, etc.

        These regidx indexes are created by parsing a gff3 file as follows:
            1.  create the array "ftr" of all UTR, CDS, exons. This will be
            processed later and pruned based on transcript types we want to keep.
            In the same go, create the hash "id2tr" of transcripts to keep
            (based on biotype) which maps from transcript_id to a transcript. At
            the same time also build the hash "gid2gene" which maps from gene_id to
            gf_gene_t pointer.
            
            2.  build "idx_cds", "idx_tscript", "idx_utr" and "idx_exon" indexes.
            Use only features from "ftr" which are present in "id2tr".

            3.  clean data that won't be needed anymore: ftr, id2tr, gid2gene.
        
    Data structures.
        idx_cds, idx_utr, idx_exon, idx_tscript:
            as described above, regidx structures for fast lookup of exons/transcripts
            overlapping a region, the payload is a pointer to tscript.cds
*/
 
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/khash.h>
#include <htslib/khash_str2int.h>
#include <htslib/kseq.h>
#include <htslib/faidx.h>
#include <errno.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include "bcftools.h"
#include "filter.h"
#include "regidx.h"
#include "kheap.h"
#include "smpl_ilist.h"
#include "rbuf.h"

#ifndef __FUNCTION__
#  define __FUNCTION__ __func__
#endif

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

// Definition of splice_region, splice_acceptor and splice_donor
#define N_SPLICE_DONOR         2      
#define N_SPLICE_REGION_EXON   3 
#define N_SPLICE_REGION_INTRON 8 

#define N_REF_PAD 10    // number of bases to avoid boundary effects

#define STRAND_REV 0
#define STRAND_FWD 1

#define TRIM_NONE   0
#define TRIM_5PRIME 1
#define TRIM_3PRIME 2

// How to treat phased/unphased genotypes
#define PHASE_REQUIRE 0     // --phase r
#define PHASE_MERGE   1     // --phase m
#define PHASE_AS_IS   2     // --phase a
#define PHASE_SKIP    3     // --phase s
#define PHASE_NON_REF 4     // --phase R
#define PHASE_DROP_GT 5     // --samples -

// Node types in the haplotype tree
#define HAP_CDS   0
#define HAP_ROOT  1 
#define HAP_SSS   2     // start/stop/splice

#define CSQ_PRINTED_UPSTREAM    (1<<0)
#define CSQ_SYNONYMOUS_VARIANT  (1<<1)
#define CSQ_MISSENSE_VARIANT    (1<<2)
#define CSQ_STOP_LOST           (1<<3)
#define CSQ_STOP_GAINED         (1<<4)
#define CSQ_INFRAME_DELETION    (1<<5)
#define CSQ_INFRAME_INSERTION   (1<<6)
#define CSQ_FRAMESHIFT_VARIANT  (1<<7)
#define CSQ_SPLICE_ACCEPTOR     (1<<8)
#define CSQ_SPLICE_DONOR        (1<<9)
#define CSQ_START_LOST          (1<<10)
#define CSQ_SPLICE_REGION       (1<<11)
#define CSQ_STOP_RETAINED       (1<<12)
#define CSQ_UTR5                (1<<13)
#define CSQ_UTR3                (1<<14)
#define CSQ_NON_CODING          (1<<15)
#define CSQ_INTRON              (1<<16)
//#define CSQ_INTERGENIC          (1<<17)
#define CSQ_INFRAME_ALTERING    (1<<18)
#define CSQ_UPSTREAM_STOP       (1<<19)     // adds * in front of the csq string
#define CSQ_INCOMPLETE_CDS      (1<<20)     // to remove START/STOP in incomplete CDS, see ENSG00000173376/synon.vcf
#define CSQ_CODING_SEQUENCE     (1<<21)     // cannot tell exactly what it is, but it does affect the coding sequence

// Haplotype-aware consequences, printed in one vcf record only, the rest has a reference @12345
#define CSQ_COMPOUND (CSQ_SYNONYMOUS_VARIANT|CSQ_MISSENSE_VARIANT|CSQ_STOP_LOST|CSQ_STOP_GAINED| \
                      CSQ_INFRAME_DELETION|CSQ_INFRAME_INSERTION|CSQ_FRAMESHIFT_VARIANT| \
                      CSQ_START_LOST|CSQ_STOP_RETAINED|CSQ_INFRAME_ALTERING|CSQ_INCOMPLETE_CDS| \
                      CSQ_UPSTREAM_STOP)
#define CSQ_START_STOP          (CSQ_STOP_LOST|CSQ_STOP_GAINED|CSQ_STOP_RETAINED|CSQ_START_LOST)

#define CSQ_PRN_STRAND(csq)     ((csq)&CSQ_COMPOUND && !((csq)&(CSQ_SPLICE_ACCEPTOR|CSQ_SPLICE_DONOR|CSQ_SPLICE_REGION)))
#define CSQ_PRN_TSCRIPT         (~(CSQ_INTRON|CSQ_NON_CODING))
#define CSQ_PRN_BIOTYPE         CSQ_NON_CODING

// see kput_vcsq()
const char *csq_strings[] = 
{
    NULL, 
    "synonymous", 
    "missense", 
    "stop_lost", 
    "stop_gained", 
    "inframe_deletion", 
    "inframe_insertion", 
    "frameshift", 
    "splice_acceptor", 
    "splice_donor", 
    "start_lost", 
    "splice_region", 
    "stop_retained", 
    "5_prime_utr", 
    "3_prime_utr", 
    "non_coding", 
    "intron", 
    "intergenic",
    "inframe_altering",
    NULL,
    NULL,
    "coding_sequence"
};


// GFF line types
#define GFF_TSCRIPT_LINE 1
#define GFF_GENE_LINE    2


/* 
    Genomic features, for fast lookup by position to overlapping features
*/
#define GF_coding_bit 6
#define GF_is_coding(x) ((x) & (1<<GF_coding_bit))
#define GF_MT_rRNA                       1                      // non-coding: 1, 2, ...
#define GF_MT_tRNA                       2
#define GF_lincRNA                       3
#define GF_miRNA                         4
#define GF_MISC_RNA                      5
#define GF_rRNA                          6
#define GF_snRNA                         7
#define GF_snoRNA                        8
#define GF_PROCESSED_TRANSCRIPT          9
#define GF_ANTISENSE                    10
#define GF_macro_lncRNA                 11
#define GF_ribozyme                     12
#define GF_sRNA                         13
#define GF_scRNA                        14
#define GF_scaRNA                       15
#define GF_SENSE_INTRONIC               16
#define GF_SENSE_OVERLAPPING            17
#define GF_PSEUDOGENE                   18
#define GF_PROCESSED_PSEUDOGENE         19
#define GF_ARTIFACT                     20
#define GF_IG_PSEUDOGENE                21
#define GF_IG_C_PSEUDOGENE              22
#define GF_IG_J_PSEUDOGENE              23
#define GF_IG_V_PSEUDOGENE              24
#define GF_TR_V_PSEUDOGENE              25
#define GF_TR_J_PSEUDOGENE              26
#define GF_MT_tRNA_PSEUDOGENE           27
#define GF_misc_RNA_PSEUDOGENE          28
#define GF_miRNA_PSEUDOGENE             29
#define GF_RIBOZYME                     30
#define GF_RETAINED_INTRON              31
#define GF_RETROTRANSPOSED              32
#define GF_tRNA_PSEUDOGENE              33
#define GF_TRANSCRIBED_PROCESSED_PSEUDOGENE     34
#define GF_TRANSCRIBED_UNPROCESSED_PSEUDOGENE   35
#define GF_TRANSCRIBED_UNITARY_PSEUDOGENE       36
#define GF_TRANSLATED_UNPROCESSED_PSEUDOGENE    37
#define GF_TRANSLATED_PROCESSED_PSEUDOGENE      38
#define GF_KNOWN_NCRNA                          39
#define GF_UNITARY_PSEUDOGENE                   40
#define GF_UNPROCESSED_PSEUDOGENE               41
#define GF_LRG_GENE                             42
#define GF_3PRIME_OVERLAPPING_ncRNA             43
#define GF_DISRUPTED_DOMAIN                     44
#define GF_vaultRNA                             45
#define GF_BIDIRECTIONAL_PROMOTER_lncRNA        46
#define GF_AMBIGUOUS_ORF                        47
#define GF_PROTEIN_CODING               (1|(1<<GF_coding_bit))  // coding: 65, 66, ...
#define GF_POLYMORPHIC_PSEUDOGENE       (2|(1<<GF_coding_bit))
#define GF_IG_C                         (3|(1<<GF_coding_bit))
#define GF_IG_D                         (4|(1<<GF_coding_bit))
#define GF_IG_J                         (5|(1<<GF_coding_bit))
#define GF_IG_LV                        (6|(1<<GF_coding_bit))
#define GF_IG_V                         (7|(1<<GF_coding_bit))
#define GF_TR_C                         (8|(1<<GF_coding_bit))
#define GF_TR_D                         (9|(1<<GF_coding_bit))
#define GF_TR_J                        (10|(1<<GF_coding_bit))
#define GF_TR_V                        (11|(1<<GF_coding_bit))
#define GF_NMD                         (12|(1<<GF_coding_bit))
#define GF_NON_STOP_DECAY              (13|(1<<GF_coding_bit))
#define GF_CDS      ((1<<(GF_coding_bit+1))+1)                  // special types: 129, 130, ...
#define GF_EXON     ((1<<(GF_coding_bit+1))+2)
#define GF_UTR3     ((1<<(GF_coding_bit+1))+3)
#define GF_UTR5     ((1<<(GF_coding_bit+1))+4)
// GF_MAX = (1<<30)-1, see hap_node_t

typedef struct _tscript_t tscript_t;
typedef struct
{
    tscript_t *tr;      // transcript
    uint32_t beg;       // the start coordinate of the CDS (on the reference strand, 0-based)
    uint32_t pos;       // 0-based index of the first exon base within the transcript (only to
                        //  update hap_node_t.sbeg in hap_init, could be calculated on the fly)
    uint32_t len;       // exon length
    uint32_t icds:30,   // exon index within the transcript
             phase:2;   // offset of the CDS
}
gf_cds_t;
typedef struct
{
    char *name;           // human readable name, e.g. ORF45
    uint8_t iseq;
}
gf_gene_t;
typedef struct
{
    uint32_t beg,end;
    tscript_t *tr;
}
gf_exon_t;
typedef enum { prime3, prime5 } utr_t;
typedef struct
{
    utr_t which;
    uint32_t beg,end;
    tscript_t *tr;
}
gf_utr_t;


/*
    Structures related to VCF output:

    vcsq_t
        information required to assemble consequence lines such as "inframe_deletion|XYZ|ENST01|+|5TY>5I|121ACG>A+124TA>T"

    vcrec_t 
        single VCF record and csq tied to this record. (Haplotype can have multiple
        consequences in several VCF records. Each record can have multiple consequences
        from multiple haplotypes.)

    csq_t
        a top-level consequence tied to a haplotype

    vbuf_t
    pos2vbuf
        VCF records with the same position clustered together for a fast lookup via pos2vbuf
*/
typedef struct _vbuf_t vbuf_t;
typedef struct _vcsq_t vcsq_t;
struct _vcsq_t
{
    uint32_t strand:1,
             type:31;   // one of CSQ_* types
    uint32_t trid;
    uint32_t biotype;   // one of GF_* types
    char *gene;         // gene name
    bcf1_t *ref;        // if type&CSQ_PRINTED_UPSTREAM, ref consequence "@1234"
    kstring_t vstr;     // variant string, eg 5TY>5I|121ACG>A+124TA>T
};
typedef struct
{
    bcf1_t *line;
    uint32_t *smpl;     // bitmask of sample consequences with first/second haplotype interleaved
    uint32_t nfmt:4, nvcsq:28, mvcsq;
    vcsq_t *vcsq;       // there can be multiple consequences for a single VCF record
}
vrec_t;
typedef struct
{
    uint32_t pos;
    vrec_t *vrec;   // vcf line that this csq is tied to; needed when printing haplotypes (hap_stage_vcf)
    int idx;        // 0-based index of the csq at the VCF line, for FMT/BCSQ
    vcsq_t type;
}
csq_t;
struct _vbuf_t
{
    vrec_t **vrec;   // buffer of VCF lines with the same position
    int n, m;
};
KHASH_MAP_INIT_INT(pos2vbuf, vbuf_t*)


/*
    Structures related to haplotype-aware consequences in coding regions

    hap_node_t
        node of a haplotype tree. Each transcript has one tree

    tscript_t
        despite its general name, it is intended for coding transcripts only

    hap_t
    hstack_t
        for traversal of the haplotype tree and braking combined
        consequences into independent parts
*/
typedef struct _hap_node_t hap_node_t;
struct _hap_node_t
{
    char *seq;          // cds segment [parent_node,this_node)
    char *var;          // variant "ref>alt"
    uint32_t type:2,    // HAP_ROOT or HAP_CDS
             csq:30;    // this node's consequence
    int dlen;           // alt minus ref length: <0 del, >0 ins, 0 substitution
    uint32_t rbeg;      // variant's VCF position (0-based, inclusive)
    int32_t rlen;       // variant's rlen; alen=rlen+dlen; fake for non CDS types
    uint32_t sbeg;      // variant's position on the spliced reference transcript (0-based, inclusive, N_REF_PAD not included)
    uint32_t icds;      // which exon does this node's variant overlaps
    hap_node_t **child, *prev;  // children haplotypes and previous coding node
    int nchild, mchild;
    bcf1_t *cur_rec, *rec;      // current VCF record and node's VCF record
    uint32_t nend;              // number of haplotypes ending in this node
    int *cur_child, mcur_child; // mapping from the allele to the currently active child
    csq_t *csq_list;            // list of haplotype's consequences, broken by position
    int ncsq_list, mcsq_list;
};
struct _tscript_t
{
    uint32_t id;        // transcript id
    uint32_t beg,end;   // transcript's beg and end coordinate (ref strand, 0-based, inclusive)
    uint32_t strand:1,  // STRAND_REV or STRAND_FWD
             ncds:31,   // number of exons
             mcds;
    gf_cds_t **cds;     // ordered list of exons
    char *ref;          // reference sequence, padded with N_REF_PAD bases on both ends
    char *sref;         // spliced reference sequence, padded with N_REF_PAD bases on both ends
    hap_node_t *root;   // root of the haplotype tree
    hap_node_t **hap;   // pointer to haplotype leaves, two for each sample
    int nhap, nsref;    // number of haplotypes and length of sref, including 2*N_REF_PAD
    uint32_t trim:2,    // complete, 5' or 3' trimmed, see TRIM_* types
             type:30;   // one of GF_* types
    gf_gene_t *gene;
};
static inline int cmp_tscript(tscript_t **a, tscript_t **b)
{
    return ( (*a)->end  < (*b)->end ) ? 1 : 0;
}
KHEAP_INIT(trhp, tscript_t*, cmp_tscript)
typedef khp_trhp_t tr_heap_t;
typedef struct
{
    hap_node_t *node;   // current node
    int ichild;         // current child in the active node
    int dlen;           // total dlen, from the root to the active node
    size_t slen;        // total sequence length, from the root to the active node
}
hstack_t;
typedef struct
{
    int mstack;
    hstack_t *stack;
    tscript_t *tr;      // tr->ref: spliced transcript on ref strand
    kstring_t sseq;     // spliced haplotype sequence on ref strand
    kstring_t tseq;     // the variable part of translated haplotype transcript, coding strand
    kstring_t tref;     // the variable part of translated reference transcript, coding strand
    uint32_t sbeg;      // stack's sbeg, for cases first node's type is HAP_SSS
    int upstream_stop;
}
hap_t;


/*
    Helper structures, only for initialization
    
    ftr_t
        temporary list of all exons, CDS, UTRs 
*/
KHASH_MAP_INIT_INT(int2tscript, tscript_t*)
KHASH_MAP_INIT_INT(int2gene, gf_gene_t*)
typedef struct
{
    int type;       // GF_CDS, GF_EXON, GF_5UTR, GF_3UTR
    uint32_t beg;
    uint32_t end;
    uint32_t trid;
    uint32_t strand:1;   // STRAND_REV,STRAND_FWD
    uint32_t phase:2;    // 0, 1 or 2
    uint32_t iseq:29;
}
ftr_t;
/*
    Mapping from GFF ID string (such as ENST00000450305 or Zm00001d027230_P001)
    to integer id.  To keep the memory requirements low, the original version
    relied on IDs in the form of a string prefix and a numerical id.  However,
    it turns out that this assumption is not valid for some ensembl GFFs, see
    for example Zea_mays.AGPv4.36.gff3.gz
 */
typedef struct
{
    void *str2id;       // khash_str2int
    int nstr, mstr;
    char **str;         // numeric id to string
}
id_tbl_t;
typedef struct
{
    // all exons, CDS, UTRs
    ftr_t *ftr;
    int nftr, mftr;

    // mapping from gene id to gf_gene_t
    kh_int2gene_t *gid2gene;

    // mapping from transcript id to tscript, for quick CDS anchoring
    kh_int2tscript_t *id2tr;

    // sequences
    void *seq2int;  // str2int hash
    char **seq;
    int nseq, mseq;

    // ignored biotypes
    void *ignored_biotypes;

    id_tbl_t gene_ids;   // temporary table for mapping between gene id (eg. Zm00001d027245) and a numeric idx
}
aux_t;

typedef struct _args_t
{
    // the main regidx lookups, from chr:beg-end to overlapping features and
    // index iterator
    regidx_t *idx_cds, *idx_utr, *idx_exon, *idx_tscript;
    regitr_t *itr;

    // temporary structures, deleted after initializtion
    aux_t init;

    // text tab-delimited output (out) or vcf/bcf output (out_fh)
    FILE *out;
    htsFile *out_fh;

    // vcf
    bcf_srs_t *sr;
    bcf_hdr_t *hdr;
    int hdr_nsmpl;          // actual number of samples in the vcf, for bcf_update_format_values()

    // include or exclude sites which match the filters
    filter_t *filter;
    char *filter_str;
    int filter_logic;       // FLT_INCLUDE or FLT_EXCLUDE

    // samples to process
    int sample_is_file;
    char *sample_list;
    smpl_ilist_t *smpl;

    char *outdir, **argv, *fa_fname, *gff_fname, *output_fname;
    char *bcsq_tag;
    int argc, output_type;
    int phase, quiet, local_csq;
    int ncsq_max, nfmt_bcsq;    // maximum number of csq per site that can be accessed from FORMAT/BCSQ
    int ncsq_small_warned;
    
    int rid;                    // current chromosome
    tr_heap_t *active_tr;       // heap of active transcripts for quick flushing
    hap_t *hap;                 // transcript haplotype recursion
    vbuf_t **vcf_buf;           // buffered VCF lines to annotate with CSQ and flush
    rbuf_t vcf_rbuf;            // round buffer indexes to vcf_buf
    kh_pos2vbuf_t *pos2vbuf;    // fast lookup of buffered lines by position
    tscript_t **rm_tr;          // buffer of transcripts to clean
    int nrm_tr, mrm_tr;
    csq_t *csq_buf;             // pool of csq not managed by hap_node_t, i.e. non-CDS csqs
    int ncsq_buf, mcsq_buf;
    id_tbl_t tscript_ids;       // mapping between transcript id (eg. Zm00001d027245_T001) and a numeric idx

    faidx_t *fai;
    kstring_t str, str2;
    int32_t *gt_arr, mgt_arr;
}
args_t;

// AAA, AAC, ...
const char *gencode = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF";
const uint8_t nt4[] =
{
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 3,4,4,4, 4,4,4,4, 4,4,4,4,
    4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 3
};
const uint8_t cnt4[] =
{
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,
    4,3,4,2, 4,4,4,1, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 0,4,4,4, 4,4,4,4, 4,4,4,4,
    4,3,4,2, 4,4,4,1, 4,4,4,4, 4,4,4,4,
    4,4,4,4, 0
};
#define dna2aa(x)  gencode[  nt4[(uint8_t)(x)[0]]<<4 |  nt4[(uint8_t)(x)[1]]<<2 |  nt4[(uint8_t)(x)[2]] ]
#define cdna2aa(x) gencode[ cnt4[(uint8_t)(x)[2]]<<4 | cnt4[(uint8_t)(x)[1]]<<2 | cnt4[(uint8_t)(x)[0]] ]

static const char *gf_strings_noncoding[] = 
{ 
    "MT_rRNA", "MT_tRNA", "lincRNA", "miRNA", "misc_RNA", "rRNA", "snRNA", "snoRNA", "processed_transcript",
    "antisense", "macro_lncRNA", "ribozyme", "sRNA", "scRNA", "scaRNA", "sense_intronic", "sense_overlapping",
    "pseudogene", "processed_pseudogene", "artifact", "IG_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene", 
    "IG_V_pseudogene", "TR_V_pseudogene", "TR_J_pseudogene", "MT_tRNA_pseudogene", "misc_RNA_pseudogene", 
    "miRNA_pseudogene", "ribozyme", "retained_intron", "retrotransposed", "Trna_pseudogene", "transcribed_processed_pseudogene", 
    "transcribed_unprocessed_pseudogene", "transcribed_unitary_pseudogene",    "translated_unprocessed_pseudogene",
    "translated_processed_pseudogene", "known_ncRNA", "unitary_pseudogene", "unprocessed_pseudogene",
    "LRG_gene", "3_prime_overlapping_ncRNA", "disrupted_domain", "vaultRNA", "bidirectional_promoter_lncRNA", "ambiguous_orf"
};
static const char *gf_strings_coding[] = { "protein_coding", "polymorphic_pseudogene", "IG_C", "IG_D", "IG_J", "IG_LV", "IG_V", "TR_C", "TR_D", "TR_J", "TR_V", "NMD", "non_stop_decay"};
static const char *gf_strings_special[] = { "CDS", "exon", "3_prime_UTR", "5_prime_UTR" };

const char *gf_type2gff_string(int type)
{
    if ( !GF_is_coding(type) )
    {
        if ( type < (1<<GF_coding_bit) ) return gf_strings_noncoding[type-1];
        type &= (1<<(GF_coding_bit+1)) - 1;
        return gf_strings_special[type - 1];
    }
    type &= (1<<GF_coding_bit) - 1;
    return gf_strings_coding[type - 1];
}

/*
    gff parsing functions
*/
static inline int feature_set_seq(args_t *args, char *chr_beg, char *chr_end)
{
    aux_t *aux = &args->init;
    char c = chr_end[1];
    chr_end[1] = 0;
    int iseq;
    if ( khash_str2int_get(aux->seq2int, chr_beg, &iseq)!=0 )
    {
        hts_expand(char*, aux->nseq+1, aux->mseq, aux->seq);
        aux->seq[aux->nseq] = strdup(chr_beg);
        iseq = khash_str2int_inc(aux->seq2int, aux->seq[aux->nseq]);
        aux->nseq++;
        assert( aux->nseq < 256 );  // see gf_gene_t.iseq
    }
    chr_end[1] = c;
    return iseq;
}
static inline char *gff_skip(const char *line, char *ss)
{
    while ( *ss && *ss!='\t' ) ss++;
    if ( !*ss ) error("[%s:%d %s] Could not parse the line: %s\n",__FILE__,__LINE__,__FUNCTION__,line);
    return ss+1;
}
static inline void gff_parse_chr(const char *line, char **chr_beg, char **chr_end)
{
    char *se = (char*) line;
    while ( *se && *se!='\t' ) se++;
    if ( !*se ) error("[%s:%d %s] Could not parse the line: %s\n",__FILE__,__LINE__,__FUNCTION__,line);
    *chr_beg = (char*) line;
    *chr_end = se-1;
}
static inline char *gff_parse_beg_end(const char *line, char *ss, uint32_t *beg, uint32_t *end)
{
    char *se = ss;
    *beg = strtol(ss, &se, 10) - 1;
    if ( ss==se ) error("[%s:%d %s] Could not parse the line:\n\t%s\n\t%s\n",__FILE__,__LINE__,__FUNCTION__,line,ss);
    ss = se+1;
    *end = strtol(ss, &se, 10) - 1;
    if ( ss==se ) error("[%s:%d %s] Could not parse the line: %s\n",__FILE__,__LINE__,__FUNCTION__,line);
    return se+1;
}
static void gff_id_init(id_tbl_t *tbl)
{
    memset(tbl, 0, sizeof(*tbl));
    tbl->str2id = khash_str2int_init();
}
static void gff_id_destroy(id_tbl_t *tbl)
{
    khash_str2int_destroy_free(tbl->str2id);
    free(tbl->str);
}
static inline uint32_t gff_id_parse(id_tbl_t *tbl, const char *line, const char *needle, char *ss)
{
    ss = strstr(ss,needle);     // e.g. "ID=transcript:"
    if ( !ss ) error("[%s:%d %s] Could not parse the line, \"%s\" not present: %s\n",__FILE__,__LINE__,__FUNCTION__,needle,line);
    ss += strlen(needle);

    char *se = ss;
    while ( *se && *se!=';' && !isspace(*se) ) se++;
    char tmp = *se;
    *se = 0;

    int id;
    if ( khash_str2int_get(tbl->str2id, ss, &id) < 0 )
    {
        id = tbl->nstr++;
        hts_expand(char*, tbl->nstr, tbl->mstr, tbl->str);
        tbl->str[id] = strdup(ss);
        int ret = khash_str2int_set(tbl->str2id, tbl->str[id], id);
    }
    *se = tmp;

    return id;
}
static inline int gff_parse_type(char *line)
{
    line = strstr(line,"ID=");
    if ( !line ) return -1;
    line += 3;
    if ( !strncmp(line,"transcript:",11) ) return GFF_TSCRIPT_LINE;
    else if ( !strncmp(line,"gene:",5) ) return GFF_GENE_LINE;
    return -1;
}
static inline int gff_parse_biotype(char *_line)
{
    char *line = strstr(_line,"biotype=");
    if ( !line ) return -1;

    line += 8;
    switch (*line)
    {
        case 'p': 
            if ( !strncmp(line,"protein_coding",14) ) return GF_PROTEIN_CODING;
            else if ( !strncmp(line,"pseudogene",10) ) return GF_PSEUDOGENE;
            else if ( !strncmp(line,"processed_transcript",20) ) return GF_PROCESSED_TRANSCRIPT;
            else if ( !strncmp(line,"processed_pseudogene",20) ) return GF_PROCESSED_PSEUDOGENE;
            else if ( !strncmp(line,"polymorphic_pseudogene",22) ) return GF_POLYMORPHIC_PSEUDOGENE;
            break;
        case 'a':
            if ( !strncmp(line,"artifact",8) ) return GF_ARTIFACT;
            else if ( !strncmp(line,"antisense",9) ) return GF_ANTISENSE;
            else if ( !strncmp(line,"ambiguous_orf",13) ) return GF_AMBIGUOUS_ORF;
            break;
        case 'I':
            if ( !strncmp(line,"IG_C_gene",9) ) return GF_IG_C;
            else if ( !strncmp(line,"IG_D_gene",9) ) return GF_IG_D;
            else if ( !strncmp(line,"IG_J_gene",9) ) return GF_IG_J;
            else if ( !strncmp(line,"IG_LV_gene",10) ) return GF_IG_LV;
            else if ( !strncmp(line,"IG_V_gene",9) ) return GF_IG_V;
            else if ( !strncmp(line,"IG_pseudogene",13) ) return GF_IG_PSEUDOGENE;
            else if ( !strncmp(line,"IG_C_pseudogene",15) ) return GF_IG_C_PSEUDOGENE;
            else if ( !strncmp(line,"IG_J_pseudogene",15) ) return GF_IG_J_PSEUDOGENE;
            else if ( !strncmp(line,"IG_V_pseudogene",15) ) return GF_IG_V_PSEUDOGENE;
            break;
        case 'T':
            if ( !strncmp(line,"TR_C_gene",9) ) return GF_TR_C;
            else if ( !strncmp(line,"TR_D_gene",9) ) return GF_TR_D;
            else if ( !strncmp(line,"TR_J_gene",9) ) return GF_TR_J;
            else if ( !strncmp(line,"TR_V_gene",9) ) return GF_TR_V;
            else if ( !strncmp(line,"TR_V_pseudogene",15) ) return GF_TR_V_PSEUDOGENE;
            else if ( !strncmp(line,"TR_J_pseudogene",15) ) return GF_TR_J_PSEUDOGENE;
            break;
        case 'M':
            if ( !strncmp(line,"Mt_tRNA_pseudogene",18) ) return GF_MT_tRNA_PSEUDOGENE;
            else if ( !strncmp(line,"Mt_tRNA",7) ) return GF_MT_tRNA;
            else if ( !strncmp(line,"Mt_rRNA",7) ) return GF_MT_tRNA;
            break;
        case 'l':
            if ( !strncmp(line,"lincRNA",7) ) return GF_lincRNA;
            break;
        case 'm':
            if ( !strncmp(line,"macro_lncRNA",12) ) return GF_macro_lncRNA;
            else if ( !strncmp(line,"misc_RNA_pseudogene",19) ) return GF_misc_RNA_PSEUDOGENE;
            else if ( !strncmp(line,"miRNA_pseudogene",16) ) return GF_miRNA_PSEUDOGENE;
            else if ( !strncmp(line,"miRNA",5) ) return GF_miRNA;
            else if ( !strncmp(line,"misc_RNA",8) ) return GF_MISC_RNA;
            break;
        case 'r':
            if ( !strncmp(line,"rRNA",4) ) return GF_rRNA;
            else if ( !strncmp(line,"ribozyme",8) ) return GF_RIBOZYME;
            else if ( !strncmp(line,"retained_intron",15) ) return GF_RETAINED_INTRON;
            else if ( !strncmp(line,"retrotransposed",15) ) return GF_RETROTRANSPOSED;
            break;
        case 's':
            if ( !strncmp(line,"snRNA",5) ) return GF_snRNA;
            else if ( !strncmp(line,"sRNA",4) ) return GF_sRNA;
            else if ( !strncmp(line,"scRNA",5) ) return GF_scRNA;
            else if ( !strncmp(line,"scaRNA",6) ) return GF_scaRNA;
            else if ( !strncmp(line,"snoRNA",6) ) return GF_snoRNA;
            else if ( !strncmp(line,"sense_intronic",14) ) return GF_SENSE_INTRONIC;
            else if ( !strncmp(line,"sense_overlapping",17) ) return GF_SENSE_OVERLAPPING;
            break;
        case 't':
            if ( !strncmp(line,"tRNA_pseudogene",15) ) return GF_tRNA_PSEUDOGENE;
            else if ( !strncmp(line,"transcribed_processed_pseudogene",32) ) return GF_TRANSCRIBED_PROCESSED_PSEUDOGENE;
            else if ( !strncmp(line,"transcribed_unprocessed_pseudogene",34) ) return GF_TRANSCRIBED_UNPROCESSED_PSEUDOGENE; 
            else if ( !strncmp(line,"transcribed_unitary_pseudogene",30) ) return GF_TRANSCRIBED_UNITARY_PSEUDOGENE;
            else if ( !strncmp(line,"translated_unprocessed_pseudogene",33) ) return GF_TRANSLATED_UNPROCESSED_PSEUDOGENE;
            else if ( !strncmp(line,"translated_processed_pseudogene",31) ) return GF_TRANSLATED_PROCESSED_PSEUDOGENE;
            break;
        case 'n':
            if ( !strncmp(line,"nonsense_mediated_decay",23) ) return GF_NMD;
            else if ( !strncmp(line,"non_stop_decay",14) ) return GF_NON_STOP_DECAY;
            break;
        case 'k':
            if ( !strncmp(line,"known_ncrna",11) ) return GF_KNOWN_NCRNA;
            break;
        case 'u':
            if ( !strncmp(line,"unitary_pseudogene",18) ) return GF_UNITARY_PSEUDOGENE;
            else if ( !strncmp(line,"unprocessed_pseudogene",22) ) return GF_UNPROCESSED_PSEUDOGENE;
            break;
        case 'L':
            if ( !strncmp(line,"LRG_gene",8) ) return GF_LRG_GENE;
            break;
        case '3':
            if ( !strncmp(line,"3prime_overlapping_ncRNA",24) ) return GF_3PRIME_OVERLAPPING_ncRNA;
            break;
        case 'd':
            if ( !strncmp(line,"disrupted_domain",16) ) return GF_DISRUPTED_DOMAIN;
            break;
        case 'v':
            if ( !strncmp(line,"vaultRNA",8) ) return GF_vaultRNA;
            break;
        case 'b':
            if ( !strncmp(line,"bidirectional_promoter_lncRNA",29) ) return GF_BIDIRECTIONAL_PROMOTER_lncRNA;
            break;
    }
    return 0;
}
static inline int gff_ignored_biotype(args_t *args, char *ss)
{
    ss = strstr(ss,"biotype=");
    if ( !ss ) return 0;

    ss += 8;
    char *se = ss, tmp;
    while ( *se && *se!=';' ) se++;
    tmp = *se;
    *se = 0;

    char *key = ss;
    int n = 0;
    if ( khash_str2int_get(args->init.ignored_biotypes, ss, &n)!=0 ) key = strdup(ss);
    khash_str2int_set(args->init.ignored_biotypes, key, n+1);

    *se = tmp;
    return 1;
}
gf_gene_t *gene_init(aux_t *aux, uint32_t gene_id)
{
    khint_t k = kh_get(int2gene, aux->gid2gene, (int)gene_id);
    gf_gene_t *gene = (k == kh_end(aux->gid2gene)) ? NULL : kh_val(aux->gid2gene, k);
    if ( !gene )
    {
        gene = (gf_gene_t*) calloc(1,sizeof(gf_gene_t));
        int ret;
        k = kh_put(int2gene, aux->gid2gene, (int)gene_id, &ret);
        kh_val(aux->gid2gene,k) = gene;
    }
    return gene;
}
void gff_parse_transcript(args_t *args, const char *line, char *ss, ftr_t *ftr)
{
    aux_t *aux = &args->init;
    int biotype = gff_parse_biotype(ss);
    if ( biotype <= 0 )
    {
        if ( !gff_ignored_biotype(args, ss) && args->quiet<2 ) fprintf(stderr,"ignored transcript: %s\n",line);
        return;
    }

    // create a mapping from transcript_id to gene_id
    uint32_t trid = gff_id_parse(&args->tscript_ids, line, "ID=transcript:", ss);
    uint32_t gene_id = gff_id_parse(&args->init.gene_ids, line, "Parent=gene:", ss);

    tscript_t *tr = (tscript_t*) calloc(1,sizeof(tscript_t));
    tr->id     = trid;
    tr->strand = ftr->strand;
    tr->gene   = gene_init(aux, gene_id);
    tr->type   = biotype;
    tr->beg    = ftr->beg;
    tr->end    = ftr->end;

    khint_t k;
    int ret;
    k = kh_put(int2tscript, aux->id2tr, (int)trid, &ret);
    kh_val(aux->id2tr,k) = tr;
}
void gff_parse_gene(args_t *args, const char *line, char *ss, char *chr_beg, char *chr_end, ftr_t *ftr)
{
    int biotype = gff_parse_biotype(ss);
    if ( biotype <= 0 )
    {
        if ( !gff_ignored_biotype(args, ss) && args->quiet<2 ) fprintf(stderr,"ignored gene: %s\n",line);
        return;
    }

    aux_t *aux = &args->init;

    // substring search for "ID=gene:ENSG00000437963"
    uint32_t gene_id = gff_id_parse(&aux->gene_ids, line, "ID=gene:", ss);
    gf_gene_t *gene = gene_init(aux, gene_id);
    assert( !gene->name );      // the gene_id should be unique

    gene->iseq = feature_set_seq(args, chr_beg,chr_end);

    // substring search for "Name=OR4F5"
    ss = strstr(chr_end+2,"Name=");
    if ( ss )
    {
        ss += 5;
        char *se = ss;
        while ( *se && *se!=';' && !isspace(*se) ) se++;
        gene->name = (char*) malloc(se-ss+1);
        memcpy(gene->name,ss,se-ss);
        gene->name[se-ss] = 0;
    }
    else
        gene->name = strdup(aux->gene_ids.str[gene_id]); // Name=<GeneName> field is not present, use the gene ID instead
}
int gff_parse(args_t *args, char *line, ftr_t *ftr)
{
    // - skip empty lines and commented lines
    // - columns 
    //      1.      chr
    //      2.      <skip>
    //      3.      CDS, transcript, gene, ...
    //      4-5.    beg,end
    //      6.      <skip>
    //      7.      strand
    //      8.      phase
    //      9.      Parent=transcript:ENST(\d+);ID=... etc

    char *ss = line;
    if ( !*ss ) return -1;      // skip blank lines
    if ( *ss=='#' ) return -1;  // skip comments

    char *chr_beg, *chr_end;
    gff_parse_chr(line, &chr_beg, &chr_end);
    ss = gff_skip(line, chr_end + 2);

    // 3. column: is this a CDS, transcript, gene, etc.
    if ( !strncmp("exon\t",ss,5) ) { ftr->type = GF_EXON; ss += 5; }
    else if ( !strncmp("CDS\t",ss,4) ) { ftr->type = GF_CDS; ss += 4; }
    else if ( !strncmp("three_prime_UTR\t",ss,16) ) { ftr->type = GF_UTR3; ss += 16; }
    else if ( !strncmp("five_prime_UTR\t",ss,15) ) { ftr->type = GF_UTR5; ss += 15; }
    else
    {
        ss = gff_skip(line, ss);
        ss = gff_parse_beg_end(line, ss, &ftr->beg,&ftr->end);
        ss = gff_skip(line, ss);
        int type = gff_parse_type(ss);
        if ( type!=GFF_TSCRIPT_LINE && type!=GFF_GENE_LINE ) 
        {
            // we ignore these, debug print to see new types:
            ss = strstr(ss,"ID=");
            if ( !ss ) return -1;   // no ID, ignore the line
            if ( !strncmp("chromosome",ss+3,10) ) return -1;
            if ( !strncmp("supercontig",ss+3,11) ) return -1;
            if ( args->quiet<2 ) fprintf(stderr,"ignored: %s\n", line);
            return -1;
        }

        // 7. column: strand
        if ( *ss == '+' ) ftr->strand = STRAND_FWD;
        else if ( *ss == '-' ) ftr->strand = STRAND_REV;
        else error("Unknown strand: %c .. %s\n", *ss,ss);

        if ( type==GFF_TSCRIPT_LINE )
            gff_parse_transcript(args, line, ss, ftr);
        else
            gff_parse_gene(args, line, ss, chr_beg, chr_end, ftr);

        return -1;
    }
    ss = gff_parse_beg_end(line, ss, &ftr->beg,&ftr->end);
    ss = gff_skip(line, ss);

    // 7. column: strand
    if ( *ss == '+' ) ftr->strand = STRAND_FWD;
    else if ( *ss == '-' ) ftr->strand = STRAND_REV;
    else { if ( args->quiet<2 ) fprintf(stderr,"Skipping unknown strand: %c\n", *ss); return -1; }
    ss += 2;

    // 8. column: phase (codon offset)
    if ( *ss == '0' ) ftr->phase = 0;
    else if ( *ss == '1' ) ftr->phase = 1;
    else if ( *ss == '2' ) ftr->phase = 2;
    else if ( *ss == '.' ) ftr->phase = 0;      // exons do not have phase
    else { if ( args->quiet<2 ) fprintf(stderr,"Skipping unknown phase: %c, %s\n", *ss, line); return -1; }
    ss += 2;

    // substring search for "Parent=transcript:ENST00000437963"
    ftr->trid = gff_id_parse(&args->tscript_ids, line, "Parent=transcript:", ss);
    ftr->iseq = feature_set_seq(args, chr_beg,chr_end);
    return 0;
}

static int cmp_cds_ptr(const void *a, const void *b)
{
    // comparison function for qsort of transcripts's CDS
    if ( (*((gf_cds_t**)a))->beg < (*((gf_cds_t**)b))->beg ) return -1;
    if ( (*((gf_cds_t**)a))->beg > (*((gf_cds_t**)b))->beg ) return 1;
    return 0;
}

static inline void chr_beg_end(aux_t *aux, int iseq, char **chr_beg, char **chr_end)
{
    *chr_beg = *chr_end = aux->seq[iseq];
    while ( (*chr_end)[1] ) (*chr_end)++;
}
tscript_t *tscript_init(aux_t *aux, uint32_t trid)
{
    khint_t k = kh_get(int2tscript, aux->id2tr, (int)trid);
    tscript_t *tr = (k == kh_end(aux->id2tr)) ? NULL : kh_val(aux->id2tr, k);
    assert( tr );
    return tr;
}
void register_cds(args_t *args, ftr_t *ftr)
{
    // Make the CDS searchable via idx_cds. Note we do not malloc tr->cds just yet.
    //  ftr is the result of parsing a gff CDS line
    aux_t *aux = &args->init;

    tscript_t *tr = tscript_init(aux, ftr->trid);
    if ( tr->strand != ftr->strand ) error("Conflicting strand in transcript %"PRIu32" .. %d vs %d\n",ftr->trid,tr->strand,ftr->strand);
    
    gf_cds_t *cds = (gf_cds_t*) malloc(sizeof(gf_cds_t));
    cds->tr    = tr;
    cds->beg   = ftr->beg;
    cds->len   = ftr->end - ftr->beg + 1;
    cds->icds  = 0;     // to keep valgrind on mac happy
    cds->phase = ftr->phase;
    
    hts_expand(gf_cds_t*,tr->ncds+1,tr->mcds,tr->cds);
    tr->cds[tr->ncds++] = cds;
}
void register_utr(args_t *args, ftr_t *ftr)
{
    aux_t *aux = &args->init;
    gf_utr_t *utr = (gf_utr_t*) malloc(sizeof(gf_utr_t));
    utr->which = ftr->type==GF_UTR3 ? prime3 : prime5;
    utr->beg   = ftr->beg;
    utr->end   = ftr->end;
    utr->tr    = tscript_init(aux, ftr->trid);

    char *chr_beg, *chr_end;
    chr_beg_end(&args->init, utr->tr->gene->iseq, &chr_beg, &chr_end);
    regidx_push(args->idx_utr, chr_beg,chr_end, utr->beg,utr->end, &utr);
}
void register_exon(args_t *args, ftr_t *ftr)
{
    aux_t *aux = &args->init;
    gf_exon_t *exon = (gf_exon_t*) malloc(sizeof(gf_exon_t));
    exon->beg = ftr->beg;
    exon->end = ftr->end;
    exon->tr  = tscript_init(aux, ftr->trid);

    char *chr_beg, *chr_end;
    chr_beg_end(&args->init, exon->tr->gene->iseq, &chr_beg, &chr_end);
    regidx_push(args->idx_exon, chr_beg,chr_end, exon->beg - N_SPLICE_REGION_INTRON, exon->end + N_SPLICE_REGION_INTRON, &exon);
}

void tscript_init_cds(args_t *args)
{
    aux_t *aux = &args->init;

    // Sort CDS in all transcripts, set offsets, check their phase, length, create index (idx_cds)
    khint_t k;
    for (k=0; k<kh_end(aux->id2tr); k++)
    {
        if ( !kh_exist(aux->id2tr, k) ) continue;
        tscript_t *tr = (tscript_t*) kh_val(aux->id2tr, k);

        // position-to-tscript lookup
        char *chr_beg, *chr_end;
        chr_beg_end(aux, tr->gene->iseq, &chr_beg, &chr_end);
        regidx_push(args->idx_tscript, chr_beg, chr_end, tr->beg, tr->end, &tr);

        if ( !tr->ncds ) continue;      // transcript with no CDS

        // sort CDs
        qsort(tr->cds, tr->ncds, sizeof(gf_cds_t*), cmp_cds_ptr);

        // trim non-coding start
        int i, len = 0;
        if ( tr->strand==STRAND_FWD )
        {
            if ( tr->cds[0]->phase ) tr->trim |= TRIM_5PRIME;
            tr->cds[0]->beg += tr->cds[0]->phase;
            tr->cds[0]->len -= tr->cds[0]->phase;
            tr->cds[0]->phase = 0;

            // sanity check phase
            for (i=0; i<tr->ncds; i++)
            {
                int phase = tr->cds[i]->phase ? 3 - tr->cds[i]->phase : 0;
                if ( phase!=len%3)
                    error("GFF3 assumption failed for transcript %s, CDS=%d: phase!=len%%3 (phase=%d, len=%d)\n",args->tscript_ids.str[tr->id],tr->cds[i]->beg+1,phase,len);
                assert( phase == len%3 );
                len += tr->cds[i]->len; 
            }
        }
        else
        {
            // Check that the phase is not bigger than CDS length. Curiously, this can really happen,
            // see Mus_musculus.GRCm38.85.gff3.gz, transcript:ENSMUST00000163141
            // todo: the same for the fwd strand
            i = tr->ncds - 1;
            int phase = tr->cds[i]->phase;
            if ( phase ) tr->trim |= TRIM_5PRIME;
            while ( i>=0 && phase > tr->cds[i]->len )
            {
                phase -= tr->cds[i]->len;
                tr->cds[i]->phase = 0;
                tr->cds[i]->len   = 0;
                i--;
            }
            tr->cds[i]->len  -= tr->cds[i]->phase;
            tr->cds[i]->phase = 0;

            // sanity check phase
            for (i=tr->ncds-1; i>=0; i--)
            {
                int phase = tr->cds[i]->phase ? 3 - tr->cds[i]->phase : 0;
                if ( phase!=len%3)
                    error("GFF3 assumption failed for transcript %s, CDS=%d: phase!=len%%3 (phase=%d, len=%d)\n",args->tscript_ids.str[tr->id],tr->cds[i]->beg+1,phase,len);
                len += tr->cds[i]->len;
            }
        }

        // set len. At the same check that CDS within a transcript do not overlap
        len = 0;
        for (i=0; i<tr->ncds; i++)
        {
            tr->cds[i]->icds = i;
            len += tr->cds[i]->len; 
            if ( !i ) continue;

            gf_cds_t *a = tr->cds[i-1];
            gf_cds_t *b = tr->cds[i];
            if ( a->beg + a->len - 1 >= b->beg ) 
                error("Error: CDS overlap in the transcript %"PRIu32": %"PRIu32"-%"PRIu32" and %"PRIu32"-%"PRIu32"\n", 
                    kh_key(aux->id2tr, k), a->beg+1,a->beg+a->len, b->beg+1,b->beg+b->len);
        }
        if ( len%3 != 0 )
        {
            // There are 13k transcripts with incomplete 3' CDS. See for example ENST00000524289
            //  http://sep2015.archive.ensembl.org/Homo_sapiens/Transcript/Sequence_cDNA?db=core;g=ENSG00000155868;r=5:157138846-157159019;t=ENST00000524289
            // Also, the incomplete CDS can be too short (1 or 2bp), so it is not enough to trim the last one.

            tr->trim |= TRIM_3PRIME;
            if ( tr->strand==STRAND_FWD )
            {
                i = tr->ncds - 1;
                while ( i>=0 && len%3 )
                {
                    int dlen = tr->cds[i]->len >= len%3 ? len%3 : tr->cds[i]->len;
                    tr->cds[i]->len -= dlen;
                    len -= dlen;
                    i--;
                }
            }
            else
            {
                i = 0;
                while ( i<tr->ncds && len%3 )
                {
                    int dlen = tr->cds[i]->len >= len%3 ? len%3 : tr->cds[i]->len;
                    tr->cds[i]->len -= dlen;
                    tr->cds[i]->beg += dlen;
                    len -= dlen;
                    i++;
                }
            }
        }

        // set CDS offsets and insert into regidx
        len=0;
        for (i=0; i<tr->ncds; i++)
        {
            tr->cds[i]->pos = len;
            len += tr->cds[i]->len;
            regidx_push(args->idx_cds, chr_beg,chr_end, tr->cds[i]->beg,tr->cds[i]->beg+tr->cds[i]->len-1, &tr->cds[i]);
        }
    }
}

void regidx_free_gf(void *payload) { free(*((gf_cds_t**)payload)); }
void regidx_free_tscript(void *payload) { tscript_t *tr = *((tscript_t**)payload); free(tr->cds); free(tr); }

void init_gff(args_t *args)
{
    aux_t *aux = &args->init;
    aux->seq2int   = khash_str2int_init();   // chrom's numeric id
    aux->gid2gene  = kh_init(int2gene);      // gene id to gf_gene_t, for idx_gene
    aux->id2tr     = kh_init(int2tscript);   // transcript id to tscript_t
    args->idx_tscript = regidx_init(NULL, NULL, regidx_free_tscript, sizeof(tscript_t*), NULL);
    aux->ignored_biotypes = khash_str2int_init();
    gff_id_init(&aux->gene_ids);
    gff_id_init(&args->tscript_ids);

    // parse gff
    kstring_t str = {0,0,0};
    htsFile *fp = hts_open(args->gff_fname,"r");
    if ( !fp ) error("Failed to read %s\n", args->gff_fname);
    while ( hts_getline(fp, KS_SEP_LINE, &str) > 0 )
    {
        hts_expand(ftr_t, aux->nftr+1, aux->mftr, aux->ftr);
        int ret = gff_parse(args, str.s, aux->ftr + aux->nftr);
        if ( !ret ) aux->nftr++;
    }
    free(str.s);
    if ( hts_close(fp)!=0 ) error("Close failed: %s\n", args->gff_fname);


    // process gff information: connect CDS and exons to transcripts
    args->idx_cds  = regidx_init(NULL, NULL, regidx_free_gf, sizeof(gf_cds_t*), NULL);
    args->idx_utr  = regidx_init(NULL, NULL, regidx_free_gf, sizeof(gf_utr_t*), NULL);
    args->idx_exon = regidx_init(NULL, NULL, regidx_free_gf, sizeof(gf_exon_t*), NULL);
    args->itr      = regitr_init(NULL);

    int i;
    for (i=0; i<aux->nftr; i++)
    {
        ftr_t *ftr = &aux->ftr[i];

        // check whether to keep this feature: is there a mapping trid -> gene_id -> gene?
        khint_t k = kh_get(int2tscript, aux->id2tr, (int)ftr->trid);
        if ( k==kh_end(aux->id2tr) ) continue;       // no such transcript

        tscript_t *tr = kh_val(aux->id2tr,k);
        if ( !tr->gene->name )
        {
            // not a supported biotype (e.g. gene:pseudogene, transcript:processed_transcript)
            regidx_free_tscript(&tr);
            kh_del(int2tscript, aux->id2tr,k);
            continue;
        }

        // populate regidx by category: 
        //      ftr->type   .. GF_CDS, GF_EXON, GF_UTR3, GF_UTR5
        //      gene->type  .. GF_PROTEIN_CODING, GF_MT_rRNA, GF_IG_C, ...
        if ( ftr->type==GF_CDS ) register_cds(args, ftr);
        else if ( ftr->type==GF_EXON ) register_exon(args, ftr);
        else if ( ftr->type==GF_UTR5 ) register_utr(args, ftr);
        else if ( ftr->type==GF_UTR3 ) register_utr(args, ftr);
        else
            error("something: %s\t%d\t%d\t%s\t%s\n", aux->seq[ftr->iseq],ftr->beg+1,ftr->end+1,args->tscript_ids.str[ftr->trid],gf_type2gff_string(ftr->type));
    }
    tscript_init_cds(args);

    if ( !args->quiet )
    {
        fprintf(stderr,"Indexed %d transcripts, %d exons, %d CDSs, %d UTRs\n", 
                regidx_nregs(args->idx_tscript),
                regidx_nregs(args->idx_exon),
                regidx_nregs(args->idx_cds),
                regidx_nregs(args->idx_utr));
    }

    free(aux->ftr);
    khash_str2int_destroy_free(aux->seq2int);
    // keeping only to destroy the genes at the end: kh_destroy(int2gene,aux->gid2gene);
    kh_destroy(int2tscript,aux->id2tr);
    free(aux->seq);
    gff_id_destroy(&aux->gene_ids);

    if ( args->quiet<2 && khash_str2int_size(aux->ignored_biotypes) )
    {
        khash_t(str2int) *ign = (khash_t(str2int)*)aux->ignored_biotypes;
        fprintf(stderr,"Ignored the following biotypes:\n");
        for (i = kh_begin(ign); i < kh_end(ign); i++)
        {
            if ( !kh_exist(ign,i)) continue;
            fprintf(stderr,"\t%dx\t.. %s\n", kh_value(ign,i), kh_key(ign,i));
        }
    }
    khash_str2int_destroy_free(aux->ignored_biotypes);
}

void init_data(args_t *args)
{
    args->nfmt_bcsq = 1 + (args->ncsq_max - 1) / 32; 

    if ( !args->quiet ) fprintf(stderr,"Parsing %s ...\n", args->gff_fname);
    init_gff(args);

    args->rid = -1;

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);

    args->fai = fai_load(args->fa_fname);
    if ( !args->fai ) error("Failed to load the fai index: %s\n", args->fa_fname);

    args->pos2vbuf  = kh_init(pos2vbuf);
    args->active_tr = khp_init(trhp);
    args->hap = (hap_t*) calloc(1,sizeof(hap_t));

    // init samples
    if ( !bcf_hdr_nsamples(args->hdr) ) args->phase = PHASE_DROP_GT;
    if ( args->sample_list && !strcmp("-",args->sample_list) )
    {
        // ignore all samples
        if ( args->output_type==FT_TAB_TEXT ) 
        {
            // significant speedup for plain VCFs
            bcf_hdr_set_samples(args->hdr,NULL,0);
        }
        args->phase = PHASE_DROP_GT;
    }
    else
        args->smpl = smpl_ilist_init(args->hdr, args->sample_list, args->sample_is_file, SMPL_STRICT);
    args->hdr_nsmpl = args->phase==PHASE_DROP_GT ? 0 : bcf_hdr_nsamples(args->hdr);

    if ( args->output_type==FT_TAB_TEXT )
    {
        args->out = args->output_fname ? fopen(args->output_fname,"w") : stdout;
        if ( !args->out ) error("Failed to open %s: %s\n", args->output_fname,strerror(errno));

        fprintf(args->out,"# This file was produced by: bcftools +csq(%s+htslib-%s)\n", bcftools_version(),hts_version());
        fprintf(args->out,"# The command line was:\tbcftools +%s", args->argv[0]);
        int i;
        for (i=1; i<args->argc; i++)
            fprintf(args->out," %s",args->argv[i]);
        fprintf(args->out,"\n");
        fprintf(args->out,"# LOG\t[2]Message\n");
        fprintf(args->out,"# CSQ"); i = 1;
        fprintf(args->out,"\t[%d]Sample", ++i);
        fprintf(args->out,"\t[%d]Haplotype", ++i);
        fprintf(args->out,"\t[%d]Chromosome", ++i);
        fprintf(args->out,"\t[%d]Position", ++i);
        fprintf(args->out,"\t[%d]Consequence", ++i);
        fprintf(args->out,"\n");
    }
    else
    {
        args->out_fh = hts_open(args->output_fname? args->output_fname : "-",hts_bcf_wmode(args->output_type));
        if ( args->out_fh == NULL ) error("Can't write to %s: %s\n", args->output_fname? args->output_fname : "standard output", strerror(errno));
        bcf_hdr_append_version(args->hdr,args->argc,args->argv,"bcftools/csq");
        bcf_hdr_printf(args->hdr,"##INFO=<ID=%s,Number=.,Type=String,Description=\"%s consequence annotation from BCFtools/csq. Format: '[*]consequence|gene|transcript|biotype[|strand|amino_acid_change|dna_change]' or, for consequences of variants split across multiple sites, a pointer to the record storing the consequences '@position'. '*' prefix indicates a consequence downstream from a stop \">",args->bcsq_tag, args->local_csq ? "Local" : "Haplotype-aware");
        if ( args->hdr_nsmpl ) 
            bcf_hdr_printf(args->hdr,"##FORMAT=<ID=%s,Number=.,Type=Integer,Description=\"Bitmask of indexes to INFO/BCSQ, with interleaved first/second haplotype. Use \\\"bcftools query -f'[%%CHROM\\t%%POS\\t%%SAMPLE\\t%%TBCSQ\\n]'\\\" to translate.\">",args->bcsq_tag);
        bcf_hdr_write(args->out_fh, args->hdr);
    }
    if ( !args->quiet ) fprintf(stderr,"Calling...\n");
}

void destroy_data(args_t *args)
{
    regidx_destroy(args->idx_cds);
    regidx_destroy(args->idx_utr);
    regidx_destroy(args->idx_exon);
    regidx_destroy(args->idx_tscript);
    regitr_destroy(args->itr);

    khint_t k,i,j;
    for (k=0; k<kh_end(args->init.gid2gene); k++)
    {
        if ( !kh_exist(args->init.gid2gene, k) ) continue;
        gf_gene_t *gene = (gf_gene_t*) kh_val(args->init.gid2gene, k);
        free(gene->name);
        free(gene);
    }
    kh_destroy(int2gene,args->init.gid2gene);

    if ( args->filter )
        filter_destroy(args->filter);

    khp_destroy(trhp,args->active_tr);
    kh_destroy(pos2vbuf,args->pos2vbuf);
    if ( args->smpl ) smpl_ilist_destroy(args->smpl);
    int ret;
    if ( args->out_fh )
        ret = hts_close(args->out_fh);
    else
        ret = fclose(args->out);
    if ( ret ) error("Error: close failed .. %s\n", args->output_fname?args->output_fname:"stdout");
    for (i=0; i<args->vcf_rbuf.m; i++)
    {
        vbuf_t *vbuf = args->vcf_buf[i];
        if ( !vbuf ) continue;
        for (j=0; j<vbuf->m; j++)
        {
            if ( !vbuf->vrec[j] ) continue;
            if ( vbuf->vrec[j]->line ) bcf_destroy(vbuf->vrec[j]->line);
            free(vbuf->vrec[j]->smpl);
            free(vbuf->vrec[j]->vcsq);
            free(vbuf->vrec[j]);
        }
        free(vbuf->vrec);
        free(vbuf);
    }
    free(args->vcf_buf);
    free(args->rm_tr);
    free(args->csq_buf);
    free(args->hap->stack);
    free(args->hap->sseq.s);
    free(args->hap->tseq.s);
    free(args->hap->tref.s);
    free(args->hap);
    fai_destroy(args->fai);
    free(args->gt_arr);
    free(args->str.s);
    free(args->str2.s);
    gff_id_destroy(&args->tscript_ids);
}

/*
    The splice_* functions are for consquences around splice sites: start,stop,splice_*
 */
#define SPLICE_VAR_REF 0   // ref: ACGT>ACGT, csq not applicable, skip completely
#define SPLICE_OUTSIDE 1   // splice acceptor or similar; csq set and is done, does not overlap the region
#define SPLICE_INSIDE  2   // overlaps coding region; csq can be set but coding prediction is needed 
#define SPLICE_OVERLAP 3   // indel overlaps region boundary, csq set but could not determine csq
typedef struct
{
    tscript_t *tr;
    struct {
        int32_t pos, rlen, alen;
        char *ref, *alt;
        bcf1_t *rec;
    } vcf;
    uint16_t check_acceptor:1,  // check distance from exon start (fwd) or end (rev)
             check_start:1,     // this is the first coding exon (relative to transcript orientation), check first (fwd) or last (rev) codon 
             check_stop:1,      // this is the last coding exon (relative to transcript orientation), check last (fwd) or first (rev) codon
             check_donor:1,     // as with check_acceptor
             check_region_beg:1,    // do/don't check for splices at this end, eg. in the first or last exon
             check_region_end:1,    // 
             check_utr:1,           // check splice sites (acceptor/donor/region_*) only if not in utr
             set_refalt:1;          // set kref,kalt, if set, check also for synonymous events
    uint32_t csq;
    int tbeg, tend;             // number of trimmed bases from beg and end of ref,alt allele
    uint32_t ref_beg,           // ref coordinates with spurious bases removed, ACC>AC can become AC>A or CC>C, whichever gives 
             ref_end;           // a more conservative csq (the first and last base in kref.s)
    kstring_t kref, kalt;       // trimmed alleles, set only with SPLICE_OLAP
}
splice_t;
void splice_init(splice_t *splice, bcf1_t *rec)
{
    memset(splice,0,sizeof(*splice));
    splice->vcf.rec  = rec;
    splice->vcf.pos  = rec->pos;
    splice->vcf.rlen = rec->rlen;
    splice->vcf.ref  = rec->d.allele[0];
}
static inline void splice_build_hap(splice_t *splice, uint32_t beg, int len)
{
    // len>0 .. beg is the first base, del filled from right
    // len<0 .. beg is the last base, del filled from left

    int rlen, alen, rbeg, abeg;     // first base to include (ref coordinates)
    if ( len<0 )
    {
        rlen = alen = -len;
        rbeg = beg - rlen + 1;
        int dlen = splice->vcf.alen - splice->vcf.rlen;
        if ( dlen<0 && beg < splice->ref_end ) // incomplete del, beg is in the middle
            dlen += splice->ref_end - beg;
        abeg = rbeg + dlen;
    }
    else
    {
        rbeg = abeg = beg;
        rlen = alen = len;
        // check for incomplete del as above??
    }

#define XDBG 0
#if XDBG
fprintf(stderr,"build_hap:  rbeg=%d + %d    abeg=%d \n",rbeg,rlen,abeg);
#endif 
    splice->kref.l = 0;
    splice->kalt.l = 0;

    // add the part before vcf.ref, in the vcf.ref and after vcf.ref
    int roff;   // how many vcf.ref bases already used
    if ( rbeg < splice->vcf.pos )
    {
        assert( splice->tr->beg <= rbeg );  // this can be extended thanks to N_REF_PAD
        kputsn(splice->tr->ref + N_REF_PAD + rbeg - splice->tr->beg, splice->vcf.pos - rbeg, &splice->kref);
        roff = 0;
    }
    else
        roff = rbeg - splice->vcf.pos;
#if XDBG
fprintf(stderr,"r1: %s  roff=%d\n",splice->kref.s,roff);
#endif

    if ( roff < splice->vcf.rlen && splice->kref.l < rlen )
    {
        int len = splice->vcf.rlen - roff;  // len still available in vcf.ref
        if ( len > rlen - splice->kref.l ) len = rlen - splice->kref.l; // how much of ref allele is still needed
        kputsn(splice->vcf.ref + roff, len, &splice->kref);
    }
#if XDBG
fprintf(stderr,"r2: %s\n",splice->kref.s);
#endif

    uint32_t end = splice->vcf.pos + splice->vcf.rlen;    // position just after the ref allele
    if ( splice->kref.l < rlen )
    {
        if ( end + rlen - splice->kref.l - 1 > splice->tr->end ) // trim, the requested sequence is too long (could be extended, see N_REF_PAD)
            rlen -= end + rlen - splice->kref.l - 1 - splice->tr->end;
        if ( splice->kref.l < rlen )
            kputsn(splice->tr->ref + N_REF_PAD + end - splice->tr->beg, rlen - splice->kref.l, &splice->kref);
    }
#if XDBG
fprintf(stderr,"r3: %s\n",splice->kref.s);
#endif


    int aoff;
    if ( abeg < splice->vcf.pos )
    {
        assert( splice->tr->beg <= abeg );
        kputsn(splice->tr->ref + N_REF_PAD + abeg - splice->tr->beg, splice->vcf.pos - abeg, &splice->kalt);
        aoff = 0;
    }
    else
        aoff = abeg - splice->vcf.pos;
#if XDBG
fprintf(stderr,"a1: %s  aoff=%d\n",splice->kalt.s,aoff);
#endif

    if ( aoff < splice->vcf.alen && splice->kalt.l < alen )
    {
        int len = splice->vcf.alen - aoff;  // len still available in vcf.alt
        if ( len > alen - splice->kalt.l ) len = alen - splice->kalt.l; // how much of alt allele is still needed
        kputsn(splice->vcf.alt + aoff, len, &splice->kalt);
        aoff -= len;
    }
    if ( aoff < 0 ) aoff = 0;
    else aoff--;
#if XDBG
fprintf(stderr,"a2: %s  aoff=%d\n",splice->kalt.s,aoff);
#endif

    end = splice->vcf.pos + splice->vcf.rlen;    // position just after the ref allele
    if ( splice->kalt.l < alen )
    {
        if ( end + alen + aoff - splice->kalt.l - 1 > splice->tr->end ) // trim, the requested sequence is too long
            alen -= end + alen + aoff - splice->kalt.l - 1 - splice->tr->end;
        if ( alen > 0 && alen > splice->kalt.l )
            kputsn(splice->tr->ref + aoff + N_REF_PAD + end - splice->tr->beg, alen - splice->kalt.l, &splice->kalt);
    }
#if XDBG
fprintf(stderr,"a3: %s\n",splice->kalt.s);
fprintf(stderr," [%s]\n [%s]\n\n",splice->kref.s,splice->kalt.s);
#endif
}
void csq_stage(args_t *args, csq_t *csq, bcf1_t *rec);
static inline int csq_stage_utr(args_t *args, regitr_t *itr, bcf1_t *rec, uint32_t trid)
{
    while ( regitr_overlap(itr) )
    {
        gf_utr_t *utr = regitr_payload(itr, gf_utr_t*);
        tscript_t *tr = utr->tr;
        if ( tr->id != trid ) continue;
        csq_t csq; 
        memset(&csq, 0, sizeof(csq_t));
        csq.pos          = rec->pos;
        csq.type.type    = utr->which==prime5 ? CSQ_UTR5 : CSQ_UTR3;
        csq.type.biotype = tr->type;
        csq.type.strand  = tr->strand;
        csq.type.trid    = tr->id;
        csq.type.gene    = tr->gene->name;
        csq_stage(args, &csq, rec);
        return csq.type.type;
    }
    return 0;
}
static inline void csq_stage_splice(args_t *args, bcf1_t *rec, tscript_t *tr, uint32_t type)
{
#if XDBG
fprintf(stderr,"csq_stage_splice %d: type=%d\n",rec->pos+1,type);
#endif
    if ( !type ) return;
    csq_t csq; 
    memset(&csq, 0, sizeof(csq_t));
    csq.pos          = rec->pos;
    csq.type.type    = type;
    csq.type.biotype = tr->type;
    csq.type.strand  = tr->strand;
    csq.type.trid    = tr->id;
    csq.type.gene    = tr->gene->name;
    csq_stage(args, &csq, rec);
}
static inline int splice_csq_ins(args_t *args, splice_t *splice, uint32_t ex_beg, uint32_t ex_end)
{
    // coordinates that matter for consequences, eg AC>ACG trimmed to C>CG, 1bp
    // before and after the inserted bases
    if ( splice->tbeg || splice->vcf.ref[0]!=splice->vcf.alt[0] )
    {
        splice->ref_beg = splice->vcf.pos + splice->tbeg - 1;
        splice->ref_end = splice->vcf.pos + splice->vcf.rlen - splice->tend;
    }
    else
    {
        if ( splice->tend ) splice->tend--;
        splice->ref_beg = splice->vcf.pos;
        splice->ref_end = splice->vcf.pos + splice->vcf.rlen - splice->tend;
    }
#if XDBG
fprintf(stderr,"ins: %s>%s .. ex=%d,%d  beg,end=%d,%d  tbeg,tend=%d,%d  check_utr=%d start,stop,beg,end=%d,%d,%d,%d\n", splice->vcf.ref,splice->vcf.alt,ex_beg,ex_end,splice->ref_beg,splice->ref_end,splice->tbeg,splice->tend,splice->check_utr,splice->check_start,splice->check_stop,splice->check_region_beg,splice->check_region_end);
#endif

    int ret;
    if ( splice->ref_beg >= ex_end )   // fully outside, beyond the exon
    {
        if ( splice->check_utr )
        {
            regitr_t *itr = regitr_init(NULL);
            const char *chr = bcf_seqname(args->hdr,splice->vcf.rec);
            if ( regidx_overlap(args->idx_utr,chr,splice->ref_beg+1,splice->ref_beg+1, itr) )     // adjacent utr
            {
                ret = csq_stage_utr(args, itr, splice->vcf.rec, splice->tr->id);
                if ( ret!=0 ) 
                {
                    regitr_destroy(itr);
                    return SPLICE_OUTSIDE; // overlaps utr
                }
            }
            regitr_destroy(itr);
        }
        if ( !splice->check_region_end ) return SPLICE_OUTSIDE;
        char *ref = NULL, *alt = NULL;
        if ( splice->set_refalt )   // seq identity is checked only when tr->ref is available
        {
            splice_build_hap(splice, ex_end+1, N_SPLICE_REGION_INTRON);
            ref = splice->kref.s, alt = splice->kalt.s;
        }
        if ( splice->ref_beg < ex_end + N_SPLICE_REGION_INTRON && splice->ref_end > ex_end + N_SPLICE_DONOR )
        {
            splice->csq |= CSQ_SPLICE_REGION;
            if ( ref && !strncmp(ref,alt,N_SPLICE_REGION_INTRON) ) splice->csq |= CSQ_SYNONYMOUS_VARIANT;
        }
        if ( splice->ref_beg < ex_end + N_SPLICE_DONOR )
        {
            if ( splice->check_donor && splice->tr->strand==STRAND_FWD ) splice->csq |= CSQ_SPLICE_DONOR;
            if ( splice->check_acceptor && splice->tr->strand==STRAND_REV ) splice->csq |= CSQ_SPLICE_ACCEPTOR;
            if ( ref && !strncmp(ref,alt,N_SPLICE_DONOR) ) splice->csq |= CSQ_SYNONYMOUS_VARIANT;
        }
        csq_stage_splice(args, splice->vcf.rec, splice->tr, splice->csq);
        return SPLICE_OUTSIDE;
    }
    if ( splice->ref_end < ex_beg || (splice->ref_end == ex_beg && !splice->check_region_beg) )    // fully outside, before the exon
    {
        if ( splice->check_utr )
        {
            regitr_t *itr = regitr_init(NULL);
            const char *chr = bcf_seqname(args->hdr,splice->vcf.rec);
            if ( regidx_overlap(args->idx_utr,chr,splice->ref_end-1,splice->ref_end-1, itr) )     // adjacent utr
            {
                ret = csq_stage_utr(args, itr, splice->vcf.rec, splice->tr->id);
                if ( ret!=0 )
                {
                    regitr_destroy(itr);
                    return SPLICE_OUTSIDE; // overlaps utr
                }
            }
            regitr_destroy(itr);
        }
        if ( !splice->check_region_beg ) return SPLICE_OUTSIDE;
        char *ref = NULL, *alt = NULL;
        if ( splice->set_refalt )   // seq identity is checked only when tr->ref is available
        {
            splice_build_hap(splice, ex_beg - N_SPLICE_REGION_INTRON, N_SPLICE_REGION_INTRON);
            ref = splice->kref.s, alt = splice->kalt.s;
        }
        if ( splice->ref_end > ex_beg - N_SPLICE_REGION_INTRON && splice->ref_beg < ex_beg - N_SPLICE_DONOR )
        {
            splice->csq |= CSQ_SPLICE_REGION;
            if ( ref && !strncmp(ref,alt,N_SPLICE_REGION_INTRON) ) splice->csq |= CSQ_SYNONYMOUS_VARIANT;
        }
        if ( splice->ref_end > ex_beg - N_SPLICE_DONOR )
        {
            if ( splice->check_donor && splice->tr->strand==STRAND_REV ) splice->csq |= CSQ_SPLICE_DONOR;
            if ( splice->check_acceptor && splice->tr->strand==STRAND_FWD ) splice->csq |= CSQ_SPLICE_ACCEPTOR;
            if ( ref && !strncmp(ref+N_SPLICE_REGION_INTRON-N_SPLICE_DONOR,alt+N_SPLICE_REGION_INTRON-N_SPLICE_DONOR,N_SPLICE_DONOR) ) splice->csq |= CSQ_SYNONYMOUS_VARIANT;
        }
        csq_stage_splice(args, splice->vcf.rec, splice->tr, splice->csq);
        return SPLICE_OUTSIDE;
    }
    // overlaps the exon or inside the exon
    // possible todo: find better alignment for frameshifting variants?
    if ( splice->ref_beg <= ex_beg + 2 )    // in the first 3bp
    {
        if ( splice->check_region_beg ) splice->csq |= CSQ_SPLICE_REGION;
        if ( splice->tr->strand==STRAND_FWD ) { if ( splice->check_start ) splice->csq |= CSQ_START_LOST; }
        else { if ( splice->check_stop ) splice->csq |= CSQ_STOP_LOST; }
    }
    if ( splice->ref_end > ex_end - 2 )
    {
        if ( splice->check_region_end ) splice->csq |= CSQ_SPLICE_REGION;
        if ( splice->tr->strand==STRAND_REV ) { if ( splice->check_start ) splice->csq |= CSQ_START_LOST; }
        else { if ( splice->check_stop ) splice->csq |= CSQ_STOP_LOST; }
    }
    if ( splice->set_refalt )
    {
        // Make sure the variant will not end up left aligned to avoid overlapping vcf records
        //      splice_build_hap(splice, splice->ref_beg, splice->vcf.alen - splice->tend - splice->tbeg + 1);
        //      splice->vcf.rlen -= splice->tbeg + splice->tend - 1;
        //      if ( splice->kref.l > splice->vcf.rlen ) { splice->kref.l = splice->vcf.rlen;  splice->kref.s[splice->kref.l] = 0; }
        if ( splice->ref_beg < splice->vcf.pos )    // this must have been caused by too much trimming from right
        {
            int dlen = splice->vcf.pos - splice->ref_beg;
            assert( dlen==1 );
            splice->tbeg += dlen;
            if ( splice->tbeg + splice->tend == splice->vcf.rlen ) splice->tend -= dlen;
            splice->ref_beg = splice->vcf.pos;
        }
        if ( splice->ref_end==ex_beg ) splice->tend--;  // prevent zero-length ref allele
        splice_build_hap(splice, splice->ref_beg, splice->vcf.alen - splice->tend - splice->tbeg + 1);
        splice->vcf.rlen -= splice->tbeg + splice->tend - 1;
        if ( splice->kref.l > splice->vcf.rlen ) { splice->kref.l = splice->vcf.rlen;  splice->kref.s[splice->kref.l] = 0; }
    }
    csq_stage_splice(args, splice->vcf.rec, splice->tr, splice->csq);
    return SPLICE_INSIDE;
}

static inline int splice_csq_del(args_t *args, splice_t *splice, uint32_t ex_beg, uint32_t ex_end)
{
    // coordinates that matter for consequences, eg AC>ACG trimmed to C>CG
    splice->ref_beg = splice->vcf.pos + splice->tbeg - 1;                       // 1b before the deleted base
    splice->ref_end = splice->vcf.pos + splice->vcf.rlen - splice->tend - 1;    // the last deleted base

#if XDBG
fprintf(stderr,"del: %s>%s .. ex=%d,%d  beg,end=%d,%d  tbeg,tend=%d,%d  check_utr=%d start,stop,beg,end=%d,%d,%d,%d\n", splice->vcf.ref,splice->vcf.alt,ex_beg,ex_end,splice->ref_beg,splice->ref_end,splice->tbeg,splice->tend,splice->check_utr,splice->check_start,splice->check_stop,splice->check_region_beg,splice->check_region_end);
#endif

    if ( splice->ref_beg + 1 < ex_beg )     // the part before the exon; ref_beg is off by -1
    {
        if ( splice->check_region_beg )
        {
            int csq = 0;
            if ( splice->check_utr )
            {
                regitr_t *itr = regitr_init(NULL);
                const char *chr = bcf_seqname(args->hdr,splice->vcf.rec);
                if ( regidx_overlap(args->idx_utr,chr,splice->ref_beg,ex_beg-1, itr) )     // adjacent utr
                    csq = csq_stage_utr(args, itr, splice->vcf.rec, splice->tr->id);
                regitr_destroy(itr);
            }
            if ( !csq )
            {
                char *ref = NULL, *alt = NULL;
                if ( splice->set_refalt )   // seq identity is checked only when tr->ref is available
                {
                    // filling from the left does not work for ENST00000341065/frame3.vcf
                    //    CAG.GTGGCCAG      CAG.GTGGCCAG
                    //    CA-.--GGCCAG  vs  CAG.---GCCAG
                    //  splice_build_hap(splice, ex_beg-1, -N_SPLICE_REGION_INTRON);
                    //
                    // filling from the right:
                    splice_build_hap(splice, ex_beg - N_SPLICE_REGION_INTRON, N_SPLICE_REGION_INTRON);
                    ref = splice->kref.s, alt = splice->kalt.s;
                }
                if ( splice->ref_end >= ex_beg - N_SPLICE_REGION_INTRON && splice->ref_beg < ex_beg - N_SPLICE_DONOR )
                {
                    splice->csq |= CSQ_SPLICE_REGION;
                    if ( ref && alt && !strncmp(ref,alt,N_SPLICE_REGION_INTRON) ) splice->csq |= CSQ_SYNONYMOUS_VARIANT;
                }
                if ( splice->ref_end >= ex_beg - N_SPLICE_DONOR )
                {
                    if ( splice->check_donor && splice->tr->strand==STRAND_REV ) splice->csq |= CSQ_SPLICE_DONOR;
                    if ( splice->check_acceptor && splice->tr->strand==STRAND_FWD ) splice->csq |= CSQ_SPLICE_ACCEPTOR;
                    if ( ref && alt && !strncmp(ref+N_SPLICE_REGION_INTRON-N_SPLICE_DONOR,alt+N_SPLICE_REGION_INTRON-N_SPLICE_DONOR,N_SPLICE_DONOR) ) splice->csq |= CSQ_SYNONYMOUS_VARIANT;
                }
            }
        }
        if ( splice->ref_end >= ex_beg ) 
        {
            splice->tbeg = splice->ref_beg - splice->vcf.pos + 1;
            splice->ref_beg = ex_beg - 1;
            if ( splice->tbeg + splice->tend == splice->vcf.alen )
            {
                // the deletion overlaps ex_beg and cannot be easily realigned to the right
                if ( !splice->tend )
                {
                    splice->csq |= CSQ_CODING_SEQUENCE;
                    return SPLICE_OVERLAP;
                }
                splice->tend--;
            }
        }
    }
    if ( ex_end < splice->ref_end )     // the part after the exon
    {
        if ( splice->check_region_end )
        {
            int csq = 0;
            if ( splice->check_utr )
            {
                regitr_t *itr = regitr_init(NULL);
                const char *chr = bcf_seqname(args->hdr,splice->vcf.rec);
                if ( regidx_overlap(args->idx_utr,chr,ex_end+1,splice->ref_end, itr) )     // adjacent utr
                    csq = csq_stage_utr(args, itr, splice->vcf.rec, splice->tr->id);
                regitr_destroy(itr);
            }
            if ( !csq )
            {
                char *ref = NULL, *alt = NULL;
                if ( splice->set_refalt )   // seq identity is checked only when tr->ref is available
                {
                    splice_build_hap(splice, ex_end+1, N_SPLICE_REGION_INTRON);  // ref,alt positioned at the first intron base
                    ref = splice->kref.s, alt = splice->kalt.s;
                }
                if ( splice->ref_beg < ex_end + N_SPLICE_REGION_INTRON && splice->ref_end > ex_end + N_SPLICE_DONOR )
                {
                    splice->csq |= CSQ_SPLICE_REGION;
                    if ( ref && alt && !strncmp(ref,alt,N_SPLICE_REGION_INTRON) ) splice->csq |= CSQ_SYNONYMOUS_VARIANT;
                }
                if ( splice->ref_beg < ex_end + N_SPLICE_DONOR )
                {
                    if ( splice->check_donor && splice->tr->strand==STRAND_FWD ) splice->csq |= CSQ_SPLICE_DONOR;
                    if ( splice->check_acceptor && splice->tr->strand==STRAND_REV ) splice->csq |= CSQ_SPLICE_ACCEPTOR;
                    if ( ref && alt && !strncmp(ref+N_SPLICE_REGION_INTRON-N_SPLICE_DONOR,alt+N_SPLICE_REGION_INTRON-N_SPLICE_DONOR,N_SPLICE_DONOR) ) splice->csq |= CSQ_SYNONYMOUS_VARIANT;
                }
            }
        }
        if ( splice->ref_beg < ex_end ) 
        {
            splice->tend = splice->vcf.rlen - (splice->ref_end - splice->vcf.pos + 1);
            splice->ref_end = ex_end;
        }
    }
    if ( splice->ref_end < ex_beg || splice->ref_beg >= ex_end )
    {
        csq_stage_splice(args, splice->vcf.rec, splice->tr, splice->csq);
        return SPLICE_OUTSIDE;
    }

    if ( splice->ref_beg < ex_beg + 2 ) // ref_beg is off by -1
    {
        if ( splice->check_region_beg ) splice->csq |= CSQ_SPLICE_REGION;
        if ( splice->tr->strand==STRAND_FWD ) { if ( splice->check_start ) splice->csq |= CSQ_START_LOST; }
        else { if ( splice->check_stop ) splice->csq |= CSQ_STOP_LOST; }
    }
    if ( splice->ref_end > ex_end - 3 )
    {
        if ( splice->check_region_end ) splice->csq |= CSQ_SPLICE_REGION;
        if ( splice->tr->strand==STRAND_REV ) { if ( splice->check_start ) splice->csq |= CSQ_START_LOST; }
        else { if ( splice->check_stop ) splice->csq |= CSQ_STOP_LOST; }
    }
    if ( splice->set_refalt )
    {
        if ( splice->tbeg>0 ) splice->tbeg--;  //why is this?
        if ( splice->vcf.rlen > splice->tbeg + splice->tend && splice->vcf.alen > splice->tbeg + splice->tend )
        {
            splice->vcf.rlen -= splice->tbeg + splice->tend;
            splice->vcf.alen -= splice->tbeg + splice->tend;
        }
        splice->kref.l = 0; kputsn(splice->vcf.ref + splice->tbeg, splice->vcf.rlen, &splice->kref); 
        splice->kalt.l = 0; kputsn(splice->vcf.alt + splice->tbeg, splice->vcf.alen, &splice->kalt); 
        if ( (splice->ref_beg+1 < ex_beg && splice->ref_end >= ex_beg) || (splice->ref_beg+1 < ex_end && splice->ref_end >= ex_end) ) // ouch, ugly ENST00000409523/long-overlapping-del.vcf
        {
            splice->csq |= (splice->ref_end - splice->ref_beg + 1)%3 ? CSQ_FRAMESHIFT_VARIANT : CSQ_INFRAME_DELETION;
            return SPLICE_OVERLAP;
        }
    }
    csq_stage_splice(args, splice->vcf.rec, splice->tr, splice->csq);
    return SPLICE_INSIDE;
}

static inline int splice_csq_mnp(args_t *args, splice_t *splice, uint32_t ex_beg, uint32_t ex_end)
{
    // not a real variant, can be ignored: eg ACGT>ACGT
    if ( splice->tbeg + splice->tend == splice->vcf.rlen ) return SPLICE_VAR_REF;

    splice->ref_beg = splice->vcf.pos + splice->tbeg;
    splice->ref_end = splice->vcf.pos + splice->vcf.rlen - splice->tend - 1;

#if XDBG
fprintf(stderr,"mnp: %s>%s .. ex=%d,%d  beg,end=%d,%d  tbeg,tend=%d,%d  check_utr=%d start,stop,beg,end=%d,%d,%d,%d\n", splice->vcf.ref,splice->vcf.alt,ex_beg,ex_end,splice->ref_beg,splice->ref_end,splice->tbeg,splice->tend,splice->check_utr,splice->check_start,splice->check_stop,splice->check_region_beg,splice->check_region_end);
#endif

    if ( splice->ref_beg < ex_beg )     // the part before the exon
    {
        if ( splice->check_region_beg )
        {
            int csq = 0;
            if ( splice->check_utr )
            {
                regitr_t *itr = regitr_init(NULL);
                const char *chr = bcf_seqname(args->hdr,splice->vcf.rec);
                if ( regidx_overlap(args->idx_utr,chr,splice->ref_beg,ex_beg-1, itr) )     // adjacent utr
                    csq = csq_stage_utr(args, itr, splice->vcf.rec, splice->tr->id);
                regitr_destroy(itr);
            }
            if ( !csq )
            {
                if ( splice->ref_end >= ex_beg - N_SPLICE_REGION_INTRON && splice->ref_beg < ex_beg - N_SPLICE_DONOR )
                    splice->csq |= CSQ_SPLICE_REGION;
                if ( splice->ref_end >= ex_beg - N_SPLICE_DONOR )
                {
                    if ( splice->check_donor && splice->tr->strand==STRAND_REV ) splice->csq |= CSQ_SPLICE_DONOR;
                    if ( splice->check_acceptor && splice->tr->strand==STRAND_FWD ) splice->csq |= CSQ_SPLICE_ACCEPTOR;
                }
            }
        }
        if ( splice->ref_end >= ex_beg ) 
        {
            splice->tbeg = splice->ref_beg - splice->vcf.pos;
            splice->ref_beg = ex_beg;
        }
    }
    if ( ex_end < splice->ref_end )     // the part after the exon
    {
        if ( splice->check_region_end )
        {
            int csq = 0;
            if ( splice->check_utr )
            {
                regitr_t *itr = regitr_init(NULL);
                const char *chr = bcf_seqname(args->hdr,splice->vcf.rec);
                if ( regidx_overlap(args->idx_utr,chr,ex_end+1,splice->ref_end, itr) )     // adjacent utr
                    csq = csq_stage_utr(args, itr, splice->vcf.rec, splice->tr->id);
                regitr_destroy(itr);
            }
            if ( !csq )
            {
                if ( splice->ref_beg <= ex_end + N_SPLICE_REGION_INTRON && splice->ref_end > ex_end + N_SPLICE_DONOR )
                    splice->csq |= CSQ_SPLICE_REGION;
                if ( splice->ref_beg <= ex_end + N_SPLICE_DONOR )
                {
                    if ( splice->check_donor && splice->tr->strand==STRAND_FWD ) splice->csq |= CSQ_SPLICE_DONOR;
                    if ( splice->check_acceptor && splice->tr->strand==STRAND_REV ) splice->csq |= CSQ_SPLICE_ACCEPTOR;
                }
            }
        }
        if ( splice->ref_beg <= ex_end ) 
        {
            splice->tend = splice->vcf.rlen - (splice->ref_end - splice->vcf.pos + 1);
            splice->ref_end = ex_end;
        }
    }
    if ( splice->ref_end < ex_beg || splice->ref_beg > ex_end )
    {
        csq_stage_splice(args, splice->vcf.rec, splice->tr, splice->csq);
        return SPLICE_OUTSIDE;
    }

    if ( splice->ref_beg < ex_beg + 3 )
    {
        if ( splice->check_region_beg ) splice->csq |= CSQ_SPLICE_REGION;
        if ( splice->tr->strand==STRAND_FWD ) { if ( splice->check_start ) splice->csq |= CSQ_START_LOST; }
        else { if ( splice->check_stop ) splice->csq |= CSQ_STOP_LOST; }
    }
    if ( splice->ref_end > ex_end - 3 )
    {
        if ( splice->check_region_end ) splice->csq |= CSQ_SPLICE_REGION;
        if ( splice->tr->strand==STRAND_REV ) { if ( splice->check_start ) splice->csq |= CSQ_START_LOST; }
        else { if ( splice->check_stop ) splice->csq |= CSQ_STOP_LOST; }
    }
    if ( splice->set_refalt )
    {
        splice->vcf.rlen -= splice->tbeg + splice->tend;
        splice->kref.l = 0; kputsn(splice->vcf.ref + splice->tbeg, splice->vcf.rlen, &splice->kref); 
        splice->kalt.l = 0; kputsn(splice->vcf.alt + splice->tbeg, splice->vcf.rlen, &splice->kalt); 
    }
    csq_stage_splice(args, splice->vcf.rec, splice->tr, splice->csq);
    return SPLICE_INSIDE;
}
static inline int splice_csq(args_t *args, splice_t *splice, uint32_t ex_beg, uint32_t ex_end)
{
    splice->csq = 0;
    splice->vcf.alen = strlen(splice->vcf.alt);

    int rlen1 = splice->vcf.rlen - 1, alen1 = splice->vcf.alen - 1, i = 0;
    splice->tbeg = 0, splice->tend = 0;

    // trim from the right, then from the left
    while ( i<=rlen1 && i<=alen1 )
    {
        if ( splice->vcf.ref[rlen1-i] != splice->vcf.alt[alen1-i] ) break;
        i++;
    }
    splice->tend = i;
    rlen1 -= i, alen1 -= i, i = 0;
    while ( i<=rlen1 && i<=alen1 )
    {
        if ( splice->vcf.ref[i] != splice->vcf.alt[i] ) break;
        i++;
    }
    splice->tbeg = i;

    // The mnp, ins and del code was split into near-identical functions for clarity and debugging;
    // possible todo: generalize once stable
    if ( splice->vcf.rlen==splice->vcf.alen ) return splice_csq_mnp(args, splice, ex_beg, ex_end);
    if ( splice->vcf.rlen < splice->vcf.alen ) return splice_csq_ins(args, splice, ex_beg, ex_end);
    if ( splice->vcf.rlen > splice->vcf.alen ) return splice_csq_del(args, splice, ex_beg, ex_end);

    return 0;
}

// return value: 0 added, 1 overlapping variant, 2 silent discard (intronic,alt=ref)
int hap_init(args_t *args, hap_node_t *parent, hap_node_t *child, gf_cds_t *cds, bcf1_t *rec, int ial)
{
    int i;
    kstring_t str = {0,0,0};
    tscript_t *tr = cds->tr;
    child->icds = cds->icds;     // index of cds in the tscript's list of exons

    splice_t splice;
    splice_init(&splice, rec);
    splice.tr = tr;
    splice.vcf.alt  = rec->d.allele[ial];
    splice.check_acceptor = splice.check_donor = splice.set_refalt = splice.check_utr = 1;
    if ( !(tr->trim & TRIM_5PRIME) )
    {
        if ( tr->strand==STRAND_FWD ) { if ( child->icds==0 ) splice.check_start = 1; }
        else { if ( child->icds==tr->ncds-1 ) splice.check_start = 1; }
    }
    if ( !(tr->trim & TRIM_3PRIME) )
    {
        if ( tr->strand==STRAND_FWD ) { if ( child->icds==tr->ncds-1 ) splice.check_stop = 1; }
        else { if ( child->icds==0 ) splice.check_stop = 1; }
    }
    if ( splice.check_start )   // do not check starts in incomplete CDS, defined as not starting with M
    {
        if ( tr->strand==STRAND_FWD ) { if ( dna2aa(tr->ref+N_REF_PAD+cds->beg-tr->beg) != 'M' ) splice.check_start = 0; }
        else { if ( cdna2aa(tr->ref+N_REF_PAD+cds->beg-tr->beg+cds->len-3) != 'M' ) splice.check_start = 0; }
    }
    if ( child->icds!=0 ) splice.check_region_beg = 1;
    if ( child->icds!=tr->ncds-1 ) splice.check_region_end = 1;

#if XDBG
fprintf(stderr,"\n%d [%s][%s]   check start:%d,stop:%d\n",splice.vcf.pos+1,splice.vcf.ref,splice.vcf.alt,splice.check_start,splice.check_stop);
#endif
    int ret = splice_csq(args, &splice, cds->beg, cds->beg + cds->len - 1);
#if XDBG
fprintf(stderr,"cds splice_csq: %d [%s][%s] .. beg,end=%d %d, ret=%d, csq=%d\n\n",splice.vcf.pos+1,splice.kref.s,splice.kalt.s,splice.ref_beg+1,splice.ref_end+1,ret,splice.csq);
#endif

    if ( ret==SPLICE_VAR_REF ) return 2;  // not a variant, eg REF=CA ALT=CA
    if ( ret==SPLICE_OUTSIDE || ret==SPLICE_OVERLAP )  // not a coding csq
    {
        free(splice.kref.s);
        free(splice.kalt.s);

        if ( !splice.csq ) return 2;        // fully intronic, no csq

        // splice_region/acceptor/donor
        child->seq  = NULL;
        child->sbeg = 0;
        child->rbeg = rec->pos;
        child->rlen = 0;
        child->dlen = 0;
        kputs(rec->d.allele[0],&str);
        kputc('>',&str);
        kputs(rec->d.allele[ial],&str);
        child->var  = str.s;
        child->type = HAP_SSS;
        child->csq  = splice.csq;
        child->prev = parent->type==HAP_SSS ? parent->prev : parent;
        child->rec  = rec;
        return 0;
    }
    if ( splice.csq & CSQ_SYNONYMOUS_VARIANT ) splice.csq &= ~CSQ_SYNONYMOUS_VARIANT;   // synonymous&splice,frame could become synonymous&frame,splice

    int dbeg = 0;
    if ( splice.ref_beg < cds->beg )
    {
        // The vcf record overlaps the exon boundary, but the variant itself
        // should fit inside since we are here. This will need more work.
        // #1475227917
        dbeg = cds->beg - splice.ref_beg;
        splice.kref.l -= dbeg;
        splice.ref_beg = cds->beg;
        assert( dbeg <= splice.kalt.l );
    }

    if ( parent->type==HAP_SSS ) parent = parent->prev;
    if ( parent->type==HAP_CDS )    
    {
        i = parent->icds;
        if ( i!=cds->icds )
        {
            // the variant is on a new exon, finish up the previous
            int len = tr->cds[i]->len - parent->rbeg - parent->rlen + tr->cds[i]->beg;
            if ( len > 0 )
                kputsn_(tr->ref + N_REF_PAD + parent->rbeg + parent->rlen - tr->beg, len, &str);
        }

        // append any skipped non-variant exons
        while ( ++i < cds->icds )
            kputsn_(tr->ref + N_REF_PAD + tr->cds[i]->beg - tr->beg, tr->cds[i]->len, &str);

        if ( parent->icds==child->icds )
        {
            int len = splice.ref_beg - parent->rbeg - parent->rlen;
            if ( len < 0 )   // overlapping variants
            {
                free(str.s);
                return 1;
            }
            kputsn_(tr->ref + N_REF_PAD + parent->rbeg + parent->rlen - tr->beg, len, &str);
        }
        else
            kputsn_(tr->ref + N_REF_PAD + cds->beg - tr->beg, splice.ref_beg - cds->beg, &str);
    }
    kputs(splice.kalt.s + dbeg, &str);

    child->seq  = str.s;
    child->sbeg = cds->pos + (splice.ref_beg - cds->beg);
    child->rbeg = splice.ref_beg;
    child->rlen = splice.kref.l;
    child->type = HAP_CDS;
    child->prev = parent;
    child->rec  = rec;
    child->csq  = splice.csq;

    // set vlen and the "ref>alt" string
    {
        int rlen = strlen(rec->d.allele[0]);
        int alen = strlen(rec->d.allele[ial]);
        child->dlen = alen - rlen;
        child->var  = (char*) malloc(rlen+alen+2);
        memcpy(child->var,rec->d.allele[0],rlen);
        child->var[rlen] = '>';
        memcpy(child->var+rlen+1,rec->d.allele[ial],alen);
        child->var[rlen+alen+1] = 0;
    }

    // yuck, the whole CDS is modified/deleted, not ready for this, todo.
    if ( child->rbeg + child->rlen > cds->beg + cds->len )
    {
        child->type = HAP_SSS;
        if ( !child->csq ) child->csq |= CSQ_CODING_SEQUENCE;  // hack, specifically for ENST00000390520/deletion-overlap.vcf
    }

    free(splice.kref.s);
    free(splice.kalt.s);
    return 0;
}
void hap_destroy(hap_node_t *hap)
{
    int i;
    for (i=0; i<hap->nchild; i++)
        if ( hap->child[i] ) hap_destroy(hap->child[i]);
    for (i=0; i<hap->mcsq_list; i++) free(hap->csq_list[i].type.vstr.s);
    free(hap->csq_list);
    free(hap->child);
    free(hap->cur_child);
    free(hap->seq);
    free(hap->var);
    free(hap);
}


/*
    ref:    spliced reference and its length (ref.l)
    seq:    part of the spliced query transcript on the reference strand to translate, its 
                length (seq.l) and the total length of the complete transcript (seq.m)
    sbeg:   seq offset within the spliced query transcript
    rbeg:   seq offset within ref, 0-based
    rend:   last base of seq within ref, plus one. If seq does not contain indels, it is rend=rbeg+seq->l
    strand: coding strand - 0:rev, 1:fwd
    tseq:   translated sequence (aa)
    fill:   frameshift, fill until the end (strand=fwd) or from the start (strand=rev)
 */
void cds_translate(kstring_t *_ref, kstring_t *_seq, uint32_t sbeg, uint32_t rbeg, uint32_t rend, int strand, kstring_t *tseq, int fill)
{
#if XDBG
fprintf(stderr,"translate: %d %d %d  fill=%d  seq.l=%d\n",sbeg,rbeg,rend,fill,(int)_seq->l);
#endif
    char tmp[3], *codon, *end;
    int i, len, npad;

    kstring_t ref = *_ref;
    kstring_t seq = *_seq;

    tseq->l = 0;
    if ( !seq.l )
    {
        kputc('?', tseq);
        return;
    }

#define DBG 0
#if DBG
 fprintf(stderr,"translate: sbeg,rbeg,rend=%d %d %d  fill=%d  seq.l=%d\n",sbeg,rbeg,rend,fill,(int)_seq->l);
 fprintf(stderr,"    ref: l=%d %s\n", (int)ref.l,ref.s);
 fprintf(stderr,"    seq: l=%d m=%d ", (int)seq.l,(int)seq.m);
 for (i=0; i<seq.l; i++) fprintf(stderr,"%c",seq.s[i]); fprintf(stderr,"\n");
 fprintf(stderr,"    sbeg,rbeg,rend: %d,%d,%d\n", sbeg,rbeg,rend);
 fprintf(stderr,"    strand,fill: %d,%d\n", strand,fill);
#endif

    if ( strand==STRAND_FWD )
    {
        // left padding
        npad = sbeg % 3;
#if DBG>1
        fprintf(stderr,"    npad: %d\n",npad);
#endif
        assert( npad<=rbeg );

        for (i=0; i<npad; i++)
            tmp[i] = ref.s[rbeg+i-npad+N_REF_PAD];
        for (; i<3 && i-npad<seq.l; i++)
            tmp[i] = seq.s[i-npad];
        len = seq.l - i + npad;    // the remaining length of padded sseq
#if DBG>1
        fprintf(stderr,"\t i=%d\n", i);
#endif
        if ( i==3 )
        {
            kputc_(dna2aa(tmp), tseq);
#if DBG>1
            fprintf(stderr,"[1]%c%c%c\n",tmp[0],tmp[1],tmp[2]);
#endif
            codon = seq.s + 3 - npad;        // next codon
            end   = codon + len - 1 - (len % 3);    // last position of a valid codon
            while ( codon < end )
            {
                kputc_(dna2aa(codon), tseq);
#if DBG>1
                fprintf(stderr,"[2]%c%c%c\n",codon[0],codon[1],codon[2]);
#endif
                codon += 3;
            }
            end = seq.s + seq.l - 1;
            for (i=0; codon+i<=end; i++) tmp[i] = codon[i];
        }

        // right padding
        codon = ref.s + rend + N_REF_PAD;
        if ( i>0 )
        {
#if DBG>1
            if(i==1)fprintf(stderr,"[3]%c\n",tmp[0]);
            if(i==2)fprintf(stderr,"[3]%c%c\n",tmp[0],tmp[1]);
#endif
            for (; i<3; i++)
            {
                tmp[i] = *codon;
                codon++;
            }
            kputc_(dna2aa(tmp), tseq);
#if DBG>1
            fprintf(stderr,"[4]%c%c%c\n",tmp[0],tmp[1],tmp[2]);
#endif
        }
        if ( fill!=0 )
        {
            end = ref.s + ref.l - N_REF_PAD;
            while ( codon+3 <= end )
            {
                kputc_(dna2aa(codon), tseq);
#if DBG>1
                fprintf(stderr,"[5]%c%c%c\t%c\n",codon[0],codon[1],codon[2],dna2aa(codon));
#endif
                codon += 3;
            }
        }
    }
    else    // STRAND_REV
    {
        // right padding - number of bases to take from ref
        npad = (seq.m - (sbeg + seq.l)) % 3; 
#if DBG>1
        fprintf(stderr,"    npad: %d\n",npad);
#endif
if ( !(npad>=0 && sbeg+seq.l+npad<=seq.m) ) fprintf(stderr,"sbeg=%d  seq.l=%d seq.m=%d\n",sbeg,(int)seq.l,(int)seq.m);
        assert( npad>=0 && sbeg+seq.l+npad<=seq.m );  // todo: first codon on the rev strand

        if ( npad==2 )
        {
            tmp[1] = ref.s[rend+N_REF_PAD];
            tmp[2] = ref.s[rend+N_REF_PAD+1];
            i = 0;
        }
        else if ( npad==1 )
        {
            tmp[2] = ref.s[rend+N_REF_PAD];
            i = 1;
        }
        else
            i = 2;

        end = seq.s + seq.l;
        for (; i>=0 && end>seq.s; i--) tmp[i] = *(--end);
#if DBG>1
        fprintf(stderr,"\t i=%d\n", i);
        if(i==1)fprintf(stderr,"[0]    %c\n",tmp[2]);
        if(i==0)fprintf(stderr,"[0]  %c%c\n",tmp[1],tmp[2]);
#endif
        if ( i==-1 )
        {
#if DBG>1
            fprintf(stderr,"[1]%c%c%c\t%c\n",tmp[0],tmp[1],tmp[2], cdna2aa(tmp));
#endif
            kputc_(cdna2aa(tmp), tseq);
            codon = end - 3;
            while ( codon >= seq.s )
            {
                kputc_(cdna2aa(codon), tseq);
#if DBG>1
                fprintf(stderr,"[2]%c%c%c\t%c\n",codon[0],codon[1],codon[2], cdna2aa(codon));
#endif
                codon -= 3;
            }
            if ( seq.s-codon==2 )
            {
                tmp[2] = seq.s[0]; 
                i = 1;
            }
            else if ( seq.s-codon==1 )
            {
                tmp[1] = seq.s[0]; 
                tmp[2] = seq.s[1];
                i = 0;
            }
            else
                i = -1;
#if DBG>1
            if(i==1)fprintf(stderr,"[3]   %c\n",tmp[2]);
            if(i==0)fprintf(stderr,"[3] %c%c\n",tmp[1],tmp[2]);
#endif
        }
        // left padding
        end = ref.s + N_REF_PAD + rbeg;
        if ( i>=0 )
        {
            for (; i>=0 && end>=ref.s; i--) tmp[i] = *(--end);
            kputc_(cdna2aa(tmp), tseq);
#if DBG>1
            fprintf(stderr,"[4]%c%c%c\t%c\n",tmp[0],tmp[1],tmp[2],cdna2aa(tmp));
#endif
        }
        if ( fill!=0 )
        {
            codon = end - 3;
            while ( codon >= ref.s + N_REF_PAD )
            {
                kputc_(cdna2aa(codon), tseq);
#if DBG>1
                fprintf(stderr,"[5]%c%c%c\t%c\n",codon[0],codon[1],codon[2],cdna2aa(codon));
#endif
                codon -= 3;
            }
        }
    }
    kputc_(0,tseq); tseq->l--;
#if DBG
 fprintf(stderr,"    tseq: %s\n", tseq->s);
#endif
}

void tscript_splice_ref(tscript_t *tr)
{
    int i, len = 0;
    for (i=0; i<tr->ncds; i++) 
        len += tr->cds[i]->len;

    tr->nsref = len + 2*N_REF_PAD;
    tr->sref  = (char*) malloc(len + 1 + 2*N_REF_PAD);
    len = 0;

    memcpy(tr->sref, tr->ref + tr->cds[0]->beg - tr->beg, N_REF_PAD);
    len += N_REF_PAD;

    for (i=0; i<tr->ncds; i++)
    {
        memcpy(tr->sref + len, tr->ref + N_REF_PAD + tr->cds[i]->beg - tr->beg, tr->cds[i]->len);
        len += tr->cds[i]->len;
    }
    memcpy(tr->sref + len, tr->ref + N_REF_PAD + tr->cds[tr->ncds-1]->beg - tr->beg, N_REF_PAD);
    len += N_REF_PAD;

    tr->sref[len] = 0;
}

// returns: 0 if consequence was added, 1 if it already exists or could not be added
int csq_push(args_t *args, csq_t *csq, bcf1_t *rec)
{
#if XDBG
fprintf(stderr,"csq_push: %d .. %d\n",rec->pos+1,csq->type.type);
#endif
    khint_t k = kh_get(pos2vbuf, args->pos2vbuf, (int)csq->pos);
    vbuf_t *vbuf = (k == kh_end(args->pos2vbuf)) ? NULL : kh_val(args->pos2vbuf, k);
    if ( !vbuf ) error("This should not happen. %s:%d  %s\n",bcf_seqname(args->hdr,rec),csq->pos+1,csq->type.vstr);

    int i;
    for (i=0; i<vbuf->n; i++)
        if ( vbuf->vrec[i]->line==rec ) break;
    if ( i==vbuf->n ) error("This should not happen.. %s:%d  %s\n", bcf_seqname(args->hdr,rec),csq->pos+1,csq->type.vstr);
    vrec_t *vrec = vbuf->vrec[i];

    // if the variant overlaps donor/acceptor and also splice region, report only donor/acceptor
    if ( csq->type.type & CSQ_SPLICE_REGION && csq->type.type & (CSQ_SPLICE_DONOR|CSQ_SPLICE_ACCEPTOR) ) 
        csq->type.type &= ~CSQ_SPLICE_REGION;

    if ( csq->type.type & CSQ_PRINTED_UPSTREAM )
    {
        for (i=0; i<vrec->nvcsq; i++)
        {
            // Same as below, to avoid records like
            //      3630 .. @3632,stop_lost|AL627309.1|ENST00000423372|protein_coding|-
            //      3632 .. stop_lost|AL627309.1|ENST00000423372|protein_coding|-|260*>260G|3630T>A+3632A>C
            if ( csq->type.type&CSQ_START_STOP && vrec->vcsq[i].type&CSQ_START_STOP )
            {
                vrec->vcsq[i] = csq->type;
                goto exit_duplicate;
            }
            if ( !(vrec->vcsq[i].type & CSQ_PRINTED_UPSTREAM) ) continue;
            if ( csq->type.ref != vrec->vcsq[i].ref ) continue;
            goto exit_duplicate;
        }
    }
    else if ( csq->type.type & CSQ_COMPOUND )
    {
        for (i=0; i<vrec->nvcsq; i++)
        {
            if ( csq->type.trid != vrec->vcsq[i].trid && (csq->type.type|vrec->vcsq[i].type)&CSQ_PRN_TSCRIPT ) continue;
            if ( csq->type.biotype != vrec->vcsq[i].biotype ) continue;
            if ( csq->type.gene != vrec->vcsq[i].gene ) continue;
            if ( csq->type.vstr.s || vrec->vcsq[i].vstr.s ) 
            {
                // This is a bit hacky, but we want a simpler and more predictable output. The splice_csq() function
                // can trigger stop/start events based on indel overlap, then another stop/start event can be triggered
                // from add_csq() or test_cds_local() based on sequence comparison, and on output we could find two
                // consequences:
                //      stop_lost|AL627309.1|ENST00000423372|protein_coding|-
                //      stop_lost&inframe_insertion|AL627309.1|ENST00000423372|protein_coding|-|260*>260CL|3630T>TAAA
                if ( !csq->type.vstr.s || !vrec->vcsq[i].vstr.s ) 
                {
                    if ( csq->type.type&CSQ_START_STOP && vrec->vcsq[i].type&CSQ_START_STOP )
                    {
                        vrec->vcsq[i].type |= csq->type.type;

                        // remove stop_lost&synonymous if stop_retained set
                        if ( vrec->vcsq[i].type&CSQ_STOP_RETAINED ) 
                            vrec->vcsq[i].type &= ~(CSQ_STOP_LOST|CSQ_SYNONYMOUS_VARIANT);

                        if ( !vrec->vcsq[i].vstr.s ) vrec->vcsq[i].vstr = csq->type.vstr;
                        goto exit_duplicate;
                    }
                    continue;
                }
                if ( strcmp(csq->type.vstr.s,vrec->vcsq[i].vstr.s) ) continue;
            }
            vrec->vcsq[i].type |= csq->type.type; 
            goto exit_duplicate;
        }
    }
    else
    {
        for (i=0; i<vrec->nvcsq; i++)
        {
            if ( csq->type.trid != vrec->vcsq[i].trid && (csq->type.type|vrec->vcsq[i].type)&CSQ_PRN_TSCRIPT) continue;
            if ( csq->type.biotype != vrec->vcsq[i].biotype ) continue;
            if ( !(vrec->vcsq[i].type & CSQ_COMPOUND) ) 
            {
                vrec->vcsq[i].type |= csq->type.type;
                goto exit_duplicate;
            }
            if ( vrec->vcsq[i].type==(vrec->vcsq[i].type|csq->type.type) ) goto exit_duplicate;
        }
    }
    // no such csq yet in this vcf record
    csq->vrec = vrec;
    csq->idx  = i;
    vrec->nvcsq++;
    hts_expand0(vcsq_t, vrec->nvcsq, vrec->mvcsq, vrec->vcsq);
    vrec->vcsq[i] = csq->type;
    return 0;

exit_duplicate:
    csq->vrec = vrec;
    csq->idx  = i;
    return 1;
}

//  soff .. position of the variant within the trimmed query transcript
//  sbeg .. position of the variant within the query transcript
//  rbeg .. position on the reference transcript (if there are no indels, then rbeg=send)
//  rpos .. VCF position
#define node2soff(i) (hap->stack[i].slen - (hap->stack[i].node->rlen + hap->stack[i].node->dlen))
#define node2sbeg(i) (hap->sbeg + node2soff(i))
#define node2send(i) (hap->sbeg + hap->stack[i].slen)
#define node2rbeg(i) (hap->stack[i].node->sbeg)
#define node2rend(i) (hap->stack[i].node->sbeg + hap->stack[i].node->rlen)
#define node2rpos(i) (hap->stack[i].node->rec->pos)

void kput_vcsq(args_t *args, vcsq_t *csq, kstring_t *str)
{
    // Remove start/stop from incomplete CDS, but only if there is another
    // consequence as something must be reported
    if ( csq->type & CSQ_INCOMPLETE_CDS && (csq->type & ~(CSQ_START_STOP|CSQ_INCOMPLETE_CDS|CSQ_UPSTREAM_STOP)) ) csq->type &= ~(CSQ_START_STOP|CSQ_INCOMPLETE_CDS);

    // Remove missense from start/stops
    if ( csq->type & CSQ_START_STOP && csq->type & CSQ_MISSENSE_VARIANT ) csq->type &= ~CSQ_MISSENSE_VARIANT;

    if ( csq->type & CSQ_PRINTED_UPSTREAM && csq->ref )
    {
        kputc_('@',str);
        kputw(csq->ref->pos+1, str);
        return;
    }
    if ( csq->type & CSQ_UPSTREAM_STOP )
        kputc_('*',str);

    int i, n = sizeof(csq_strings)/sizeof(char*);
    for (i=1; i<n; i++)
        if ( csq_strings[i] && csq->type&(1<<i) ) { kputs(csq_strings[i],str); break; }
    i++;
    for (; i<n; i++)
        if ( csq_strings[i] && csq->type&(1<<i) ) { kputc_('&',str); kputs(csq_strings[i],str); }

    kputc_('|', str);
    if ( csq->gene ) kputs(csq->gene , str);

    kputc_('|', str);
    if ( csq->type & CSQ_PRN_TSCRIPT ) kputs(args->tscript_ids.str[csq->trid], str);

    kputc_('|', str);
    kputs(gf_type2gff_string(csq->biotype), str);

    if ( CSQ_PRN_STRAND(csq->type) || csq->vstr.l )
        kputs(csq->strand==STRAND_FWD ? "|+" : "|-", str);

    if ( csq->vstr.l )
        kputs(csq->vstr.s, str);
}

void hap_add_csq(args_t *args, hap_t *hap, hap_node_t *node, int tlen, int ibeg, int iend, int dlen, int indel)
{
    int i;
    tscript_t *tr = hap->tr;
    int ref_node = tr->strand==STRAND_FWD ? ibeg : iend;

    int icsq = node->ncsq_list++;
    hts_expand0(csq_t,node->ncsq_list,node->mcsq_list,node->csq_list);
    csq_t *csq = &node->csq_list[icsq];
    csq->pos  = hap->stack[ref_node].node->rec->pos;
    csq->type.trid    = tr->id;
    csq->type.gene    = tr->gene->name;
    csq->type.strand  = tr->strand;
    csq->type.biotype = tr->type;

    // only now we see the translated sequence and can determine if the stop/start changes are real
    int rm_csq = 0; 
    csq->type.type = 0;
    for (i=ibeg; i<=iend; i++)
        csq->type.type |= hap->stack[i].node->csq & CSQ_COMPOUND;
    if ( dlen==0 && indel ) csq->type.type |= CSQ_INFRAME_ALTERING;

    int has_upstream_stop = hap->upstream_stop;
    if ( hap->stack[ibeg].node->type != HAP_SSS )
    {
        // check for truncating stops
        for (i=0; i<hap->tref.l; i++)
            if ( hap->tref.s[i]=='*' ) break;
        if ( i!=hap->tref.l )
        {
            hap->tref.l = i+1;
            hap->tref.s[i+1] = 0;
        }
        for (i=0; i<hap->tseq.l; i++)
            if ( hap->tseq.s[i]=='*' ) break;
        if ( i!=hap->tseq.l )
        {
            hap->tseq.l = i+1;
            hap->tseq.s[i+1] = 0;
            hap->upstream_stop = 1;
        }
        if ( csq->type.type & CSQ_STOP_LOST )
        {
            if ( hap->tref.s[hap->tref.l-1]=='*' && hap->tref.s[hap->tref.l-1] == hap->tseq.s[hap->tseq.l-1] ) 
            {
                rm_csq |= CSQ_STOP_LOST;
                csq->type.type |= CSQ_STOP_RETAINED;
            }
            else if ( hap->tref.s[hap->tref.l-1]!='*' )
            {
                // This is CDS 3' incomplete ENSG00000173376/synon.vcf, can also be missense
                // We observe in real data a change to a stop, ENST00000528237/retained-stop-incomplete-cds.vcf
                if ( hap->tseq.s[hap->tseq.l-1] == '*' )
                {
                    rm_csq |= CSQ_STOP_GAINED;
                    csq->type.type |= CSQ_STOP_RETAINED;
                }
                else
                    csq->type.type |= CSQ_INCOMPLETE_CDS;
            }
        }
        if ( csq->type.type & CSQ_START_LOST && hap->tref.s[0]!='M' )
        {
            rm_csq |= CSQ_START_LOST;
            csq->type.type &= ~CSQ_START_LOST;
        }
        if ( dlen!=0 )
        {
            if ( dlen%3 )
                csq->type.type |= CSQ_FRAMESHIFT_VARIANT;
            else if ( dlen<0 )
                csq->type.type |= CSQ_INFRAME_DELETION;
            else
                csq->type.type |= CSQ_INFRAME_INSERTION;
        }
        else
        {
            for (i=0; i<hap->tref.l; i++) 
                if ( hap->tref.s[i] != hap->tseq.s[i] ) break;
            if ( i==hap->tref.l )
                csq->type.type |= CSQ_SYNONYMOUS_VARIANT;
            else if ( hap->tref.s[i] ==  '*' )
                csq->type.type |= CSQ_STOP_LOST;
            else if ( hap->tseq.s[i] ==  '*' )
                csq->type.type |= CSQ_STOP_GAINED;
            else
                csq->type.type |= CSQ_MISSENSE_VARIANT;
        }
    }
    if ( has_upstream_stop ) csq->type.type |= CSQ_UPSTREAM_STOP;
    csq->type.type &= ~rm_csq;

    if ( hap->stack[ibeg].node->type == HAP_SSS  )
    {
        node->csq_list[icsq].type.type   |= hap->stack[ibeg].node->csq & ~rm_csq;
        node->csq_list[icsq].type.ref     = hap->stack[ibeg].node->rec;
        node->csq_list[icsq].type.biotype = tr->type;
        csq_push(args, node->csq_list+icsq, hap->stack[ibeg].node->rec);
        return;
    }

    kstring_t str = node->csq_list[icsq].type.vstr;
    str.l = 0;

    // create the aa variant string
    int aa_rbeg = tr->strand==STRAND_FWD ? node2rbeg(ibeg)/3+1 : (hap->tr->nsref - 2*N_REF_PAD - node2rend(iend))/3+1;
    int aa_sbeg = tr->strand==STRAND_FWD ? node2sbeg(ibeg)/3+1 : (tlen - node2send(iend))/3+1;
    kputc_('|', &str);
    kputw(aa_rbeg, &str);
    kputs(hap->tref.s, &str);
    if ( !(csq->type.type & CSQ_SYNONYMOUS_VARIANT) )
    {
        kputc_('>', &str);
        kputw(aa_sbeg, &str);
        kputs(hap->tseq.s, &str);
    }
    kputc_('|', &str);

    // create the dna variant string and, in case of combined variants,
    // insert silent CSQ_PRINTED_UPSTREAM variants
    for (i=ibeg; i<=iend; i++)
    {
        if ( i>ibeg ) kputc_('+', &str);
        kputw(node2rpos(i)+1, &str);
        kputs(hap->stack[i].node->var, &str);
    }
    node->csq_list[icsq].type.vstr = str;
    csq_push(args, node->csq_list+icsq, hap->stack[ref_node].node->rec);

    for (i=ibeg; i<=iend; i++)
    {
        // csq are printed at one position only for combined variants, the rest is
        // silent and references the first
        if ( hap->stack[i].node->csq & ~CSQ_COMPOUND )
        {
            node->ncsq_list++;
            hts_expand0(csq_t,node->ncsq_list,node->mcsq_list,node->csq_list);
            csq_t *tmp_csq = &node->csq_list[node->ncsq_list - 1];
            tmp_csq->pos  = hap->stack[i].node->rec->pos;
            tmp_csq->type.trid    = tr->id;
            tmp_csq->type.gene    = tr->gene->name;
            tmp_csq->type.strand  = tr->strand;
            tmp_csq->type.type    = hap->stack[i].node->csq & ~CSQ_COMPOUND & ~rm_csq;
            tmp_csq->type.biotype = tr->type;
            tmp_csq->type.vstr.l  = 0;
            kputs(str.s,&tmp_csq->type.vstr);
            csq_push(args, tmp_csq, hap->stack[i].node->rec);
        }
        if ( i!=ref_node && (node->csq_list[icsq].type.type & CSQ_COMPOUND || !(hap->stack[i].node->csq & ~CSQ_COMPOUND)) )
        {
            node->ncsq_list++;
            hts_expand0(csq_t,node->ncsq_list,node->mcsq_list,node->csq_list);
            csq_t *tmp_csq = &node->csq_list[node->ncsq_list - 1];
            tmp_csq->pos  = hap->stack[i].node->rec->pos;
            tmp_csq->type.trid    = tr->id;
            tmp_csq->type.gene    = tr->gene->name;
            tmp_csq->type.strand  = tr->strand;
            tmp_csq->type.type    = CSQ_PRINTED_UPSTREAM | hap->stack[i].node->csq;
            tmp_csq->type.biotype = tr->type;
            tmp_csq->type.ref     = hap->stack[ref_node].node->rec;
            tmp_csq->type.vstr.l  = 0;
            csq_push(args, tmp_csq, hap->stack[i].node->rec);
        }
    }
}

void hap_finalize(args_t *args, hap_t *hap)
{
    tscript_t *tr = hap->tr;
    if ( !tr->sref )
        tscript_splice_ref(tr);

    kstring_t sref;
    sref.s = tr->sref;
    sref.l = tr->nsref;
    sref.m = sref.l;

    int istack = 0;
    hts_expand(hstack_t,1,hap->mstack,hap->stack);

    hap->sseq.l = 0;
    hap->tseq.l = 0;
    hap->stack[0].node = tr->root;
    hap->stack[0].ichild = -1;
    hap->stack[0].slen = 0;
    hap->stack[0].dlen = 0;

    while ( istack>=0 )
    {
        hstack_t *stack  = &hap->stack[istack];
        hap_node_t *node = hap->stack[istack].node;
        while ( ++hap->stack[istack].ichild < node->nchild )
        {
            if ( node->child[stack->ichild] ) break;
        }
        if ( stack->ichild == node->nchild ) { istack--; continue; }

        node = node->child[stack->ichild];

        istack++;
        hts_expand(hstack_t,istack+1,hap->mstack,hap->stack);
        stack = &hap->stack[istack-1];

        hap->stack[istack].node = node;
        hap->stack[istack].ichild = -1;

        hap->sseq.l = stack->slen;
        if ( node->type==HAP_CDS ) kputs(node->seq, &hap->sseq);
        hap->stack[istack].slen = hap->sseq.l;
        hap->stack[istack].dlen = hap->stack[istack-1].dlen + node->dlen;

        if ( !node->nend ) continue;    // not a leaf node

        // The spliced sequence has been built for the current haplotype and stored
        // in hap->sseq. Now we break it and output as independent parts
        
        kstring_t sseq;
        sseq.m = sref.m - 2*N_REF_PAD + hap->stack[istack].dlen;  // total length of the spliced query transcript
        hap->upstream_stop = 0;

        int i = 1, dlen = 0, ibeg, indel = 0;
        while ( i<istack && hap->stack[i].node->type == HAP_SSS ) i++;
        hap->sbeg = hap->stack[i].node->sbeg;

        if ( tr->strand==STRAND_FWD )
        {
            i = 0, ibeg = -1;
            while ( ++i <= istack )
            {
                if ( hap->stack[i].node->type == HAP_SSS )
                {
                    // start/stop/splice site overlap: don't know how to build the haplotypes correctly, skipping
                    hap_add_csq(args,hap,node,0,i,i,0,0);
                    continue;
                }
                dlen += hap->stack[i].node->dlen;
                if ( hap->stack[i].node->dlen ) indel = 1;
                if ( i<istack )
                {
                    if ( dlen%3 )   // frameshift
                    {
                        if ( ibeg==-1 ) ibeg = i;
                        continue;
                    }
                    int icur  = node2sbeg(i);
                    int inext = node2sbeg(i+1);
                    if ( icur/3 == inext/3 )    // in the same codon, can't be flushed yet
                    {
                        if ( ibeg==-1 ) ibeg = i;
                        continue;
                    }
                }
                if ( ibeg<0 ) ibeg = i;

                int ioff = node2soff(ibeg);
                int icur = node2sbeg(ibeg);
                int rbeg = node2rbeg(ibeg);
                int rend = node2rend(i);
                int fill = dlen%3;

                // alt
                if ( hap->sseq.l )
                {
                    sseq.l = hap->stack[i].slen - ioff;
                    sseq.s = hap->sseq.s + ioff;
                }
                else    // splice site overlap, see #1475227917
                    sseq.l = fill = 0;
                cds_translate(&sref, &sseq, icur,rbeg,rend, tr->strand, &hap->tseq, fill);

                // ref
                sseq.l = node2rend(i) - rbeg;
                sseq.s = sref.s + N_REF_PAD + rbeg;
                sseq.m = sref.m - 2*N_REF_PAD;
                cds_translate(&sref, &sseq, rbeg,rbeg,rend, tr->strand, &hap->tref, fill);
                sseq.m = sref.m - 2*N_REF_PAD + hap->stack[istack].dlen;

                hap_add_csq(args,hap,node,0, ibeg,i,dlen,indel);
                ibeg = -1;
                dlen = 0;
                indel = 0;
            }
        }
        else
        {
            i = istack + 1, ibeg = -1;
            while ( --i > 0 )
            {
                if ( hap->stack[i].node->type == HAP_SSS )
                {
                    hap_add_csq(args,hap,node,0,i,i,0,0);
                    continue;
                }
                dlen += hap->stack[i].node->dlen;
                if ( hap->stack[i].node->dlen ) indel = 1;
                if ( i>1 && hap->stack[i-1].node->type != HAP_SSS )
                {
                    if ( dlen%3 )
                    {
                        if ( ibeg==-1 ) ibeg = i;
                        continue;
                    }
                    int icur  = sseq.m - 1 - node2sbeg(i);
                    int inext = sseq.m - 1 - node2sbeg(i-1);
                    if ( icur/3 == inext/3 )
                    {
                        if ( ibeg==-1 ) ibeg = i;
                        continue;
                    }
                }
                if ( ibeg<0 ) ibeg = i;
                int ioff = node2soff(i);
                int icur = node2sbeg(i);
                int rbeg = node2rbeg(i);
                int rend = node2rend(ibeg);
                int fill = dlen%3;

                // alt
                if ( hap->sseq.l )
                {
                    sseq.l = hap->stack[ibeg].slen - ioff;
                    sseq.s = hap->sseq.s + ioff;
                }
                else    // splice site overlap, see #1475227917
                    sseq.l = fill = 0;
                cds_translate(&sref, &sseq, icur,rbeg,rend, tr->strand, &hap->tseq, fill);

                // ref
                sseq.l = node2rend(ibeg) - rbeg;
                sseq.s = sref.s + N_REF_PAD + rbeg;
                sseq.m = sref.m - 2*N_REF_PAD;
                cds_translate(&sref, &sseq, rbeg,rbeg,rend, tr->strand, &hap->tref, fill);
                sseq.m = sref.m - 2*N_REF_PAD + hap->stack[istack].dlen;

                hap_add_csq(args,hap,node,sseq.m, i,ibeg,dlen,indel);
                ibeg = -1;
                dlen = 0;
                indel = 0;
            }
        }
    }
}

static inline void csq_print_text(args_t *args, csq_t *csq, int ismpl, int ihap)
{
    if ( csq->type.type & CSQ_PRINTED_UPSTREAM ) return;

    char *smpl = ismpl >= 0 ? args->hdr->samples[ismpl] : "-";
    const char *chr = bcf_hdr_id2name(args->hdr,args->rid);

    fprintf(args->out,"CSQ\t%s\t", smpl);
    if ( ihap>0 )
        fprintf(args->out,"%d", ihap);
    else
        fprintf(args->out,"-");

    args->str.l = 0;
    kput_vcsq(args, &csq->type, &args->str);
    fprintf(args->out,"\t%s\t%d\t%s\n",chr,csq->pos+1,args->str.s);
}
static inline void hap_print_text(args_t *args, tscript_t *tr, int ismpl, int ihap, hap_node_t *node)
{
    if ( !node || !node->ncsq_list ) return;

    char *smpl = ismpl >= 0 ? args->hdr->samples[ismpl] : "-";
    const char *chr = bcf_hdr_id2name(args->hdr,args->rid);

    int i;
    for (i=0; i<node->ncsq_list; i++)
    {
        csq_t *csq = node->csq_list + i;
        if ( csq->type.type & CSQ_PRINTED_UPSTREAM ) continue;
        assert( csq->type.vstr.l );

        fprintf(args->out,"CSQ\t%s\t", smpl);
        if ( ihap>0 )
            fprintf(args->out,"%d", ihap);
        else
            fprintf(args->out,"-");

        args->str.l = 0;
        kput_vcsq(args, &csq->type, &args->str);
        fprintf(args->out,"\t%s\t%d\t%s\n",chr,csq->pos+1,args->str.s);
    }
}

static inline void hap_stage_vcf(args_t *args, tscript_t *tr, int ismpl, int ihap, hap_node_t *node)
{
    if ( !node || !node->ncsq_list || ismpl<0 ) return;

    int i;
    for (i=0; i<node->ncsq_list; i++)
    {
        csq_t *csq = node->csq_list + i;
        vrec_t *vrec = csq->vrec;
        int icsq = 2*csq->idx + ihap;
        if ( icsq >= args->ncsq_max ) // more than ncsq_max consequences, so can't fit it in FMT
        {
            int print_warning = 1;
            if ( args->quiet )
            {
                if ( args->quiet > 1 || args->ncsq_small_warned ) print_warning = 0;
                args->ncsq_small_warned = 1;
            }
            if ( print_warning )
            {
                fprintf(stderr,"Warning: --ncsq %d is too small to annotate %s at %s:%d with %d-th csq\n",
                        args->ncsq_max/2,args->hdr->samples[ismpl],bcf_hdr_id2name(args->hdr,args->rid),vrec->line->pos+1,csq->idx+1);
                if ( args->quiet ) fprintf(stderr,"(This warning is printed only once)\n");
            }
            break;
        }
        if ( vrec->nfmt < 1 + icsq/32 ) vrec->nfmt = 1 + icsq/32;
        vrec->smpl[ismpl*args->nfmt_bcsq + icsq/32] |= 1 << (icsq % 32);
    }
}

void hap_flush(args_t *args, uint32_t pos)
{
    int i,j;
    tr_heap_t *heap = args->active_tr;

    while ( heap->ndat && heap->dat[0]->end<=pos )
    {
        tscript_t *tr = heap->dat[0];
        khp_delete(trhp, heap);

        args->hap->tr = tr;
        if ( tr->root && tr->root->nchild ) // normal, non-localized calling
        {
            hap_finalize(args, args->hap);

            if ( args->output_type==FT_TAB_TEXT )   // plain text output, not a vcf
            {
                if ( args->phase==PHASE_DROP_GT )
                    hap_print_text(args, tr, -1,0, tr->hap[0]);
                else
                {
                    for (i=0; i<args->smpl->n; i++)
                    {
                        for (j=0; j<2; j++)
                            hap_print_text(args, tr, args->smpl->idx[i],j+1, tr->hap[i*2+j]);
                    }
                }
            }
            else if ( args->phase!=PHASE_DROP_GT )
            {
                for (i=0; i<args->smpl->n; i++)
                {
                    for (j=0; j<2; j++)
                        hap_stage_vcf(args, tr, args->smpl->idx[i],j, tr->hap[i*2+j]);
                }
            }
        }

        // mark the transcript for deletion. Cannot delete it immediately because
        // by-position VCF output will need them when flushed by vcf_buf_push
        args->nrm_tr++;
        hts_expand(tscript_t*,args->nrm_tr,args->mrm_tr,args->rm_tr);
        args->rm_tr[args->nrm_tr-1] = tr;
    }
}

#define SWAP(type_t, a, b) { type_t t = a; a = b; b = t; }

void vbuf_push(args_t *args, bcf1_t **rec_ptr)
{
    int i;

    assert(rec_ptr);
    bcf1_t *rec = *rec_ptr;

    // check for duplicate records
    i = args->vcf_rbuf.n ? rbuf_last(&args->vcf_rbuf) : -1;
    if ( i<0 || args->vcf_buf[i]->vrec[0]->line->pos!=rec->pos ) 
    {
        // vcf record with a new pos
        rbuf_expand0(&args->vcf_rbuf, vbuf_t*, args->vcf_rbuf.n+1, args->vcf_buf);
        i = rbuf_append(&args->vcf_rbuf);
        if ( !args->vcf_buf[i] ) args->vcf_buf[i] = (vbuf_t*) calloc(1,sizeof(vbuf_t));
        args->vcf_buf[i]->n = 0;
    }
    vbuf_t *vbuf = args->vcf_buf[i];
    vbuf->n++;
    hts_expand0(vrec_t*, vbuf->n, vbuf->m, vbuf->vrec);
    if ( !vbuf->vrec[vbuf->n - 1] )
        vbuf->vrec[vbuf->n - 1] = (vrec_t*) calloc(1,sizeof(vrec_t));

    vrec_t *vrec = vbuf->vrec[vbuf->n - 1];
    if ( args->phase!=PHASE_DROP_GT && args->smpl->n )
    {
        if ( !vrec->smpl ) vrec->smpl = (uint32_t*) calloc(args->hdr_nsmpl,sizeof(*vrec->smpl) * args->nfmt_bcsq);
        else memset(vrec->smpl,0,args->hdr_nsmpl*sizeof(*vrec->smpl) * args->nfmt_bcsq);
    }
    if ( !vrec->line ) vrec->line = bcf_init1();
    SWAP(bcf1_t*, (*rec_ptr), vrec->line);

    int ret;
    khint_t k = kh_put(pos2vbuf, args->pos2vbuf, (int)rec->pos, &ret);
    kh_val(args->pos2vbuf,k) = vbuf;
}

void vbuf_flush(args_t *args)
{
    if ( args->active_tr->ndat ) return; // cannot output buffered VCF lines (args.vbuf) until all active transcripts are gone

    int i,j;
    while ( (i=rbuf_shift(&args->vcf_rbuf))>=0 )
    {
        vbuf_t *vbuf = args->vcf_buf[i];
        for (i=0; i<vbuf->n; i++)
        {
            vrec_t *vrec = vbuf->vrec[i];
            if ( !args->out_fh ) // not a VCF output
            {
                vrec->nvcsq = 0;
                continue;
            }
            if ( !vrec->nvcsq )
            {
                bcf_write(args->out_fh, args->hdr, vrec->line);
                continue;
            }
            
            args->str.l = 0;
            kput_vcsq(args, &vrec->vcsq[0], &args->str);
            for (j=1; j<vrec->nvcsq; j++)
            {
                kputc_(',', &args->str);
                kput_vcsq(args, &vrec->vcsq[j], &args->str);
            }
            bcf_update_info_string(args->hdr, vrec->line, args->bcsq_tag, args->str.s);
            if ( args->hdr_nsmpl )
            {
                if ( vrec->nfmt < args->nfmt_bcsq )
                    for (j=1; j<args->hdr_nsmpl; j++) memcpy(vrec->smpl+j*vrec->nfmt, vrec->smpl+j*args->nfmt_bcsq, vrec->nfmt*sizeof(*vrec->smpl));
                bcf_update_format_int32(args->hdr, vrec->line, args->bcsq_tag, vrec->smpl, args->hdr_nsmpl*vrec->nfmt);
            }
            vrec->nvcsq = 0;
            bcf_write(args->out_fh, args->hdr, vrec->line);
        }
        if ( vbuf->n )
        {
            khint_t k = kh_get(pos2vbuf, args->pos2vbuf, vbuf->vrec[0]->line->pos);
            if ( k != kh_end(args->pos2vbuf) ) kh_del(pos2vbuf, args->pos2vbuf, k);
        }
        vbuf->n = 0;
    }

    for (i=0; i<args->nrm_tr; i++)
    {
        tscript_t *tr = args->rm_tr[i];
        if ( tr->root ) hap_destroy(tr->root);
        tr->root = NULL;
        free(tr->hap);
        free(tr->ref);
        free(tr->sref);
    }
    args->nrm_tr = 0;
    args->ncsq_buf = 0;
}

void tscript_init_ref(args_t *args, tscript_t *tr, const char *chr)
{
    int i, len;
    int pad_beg = tr->beg >= N_REF_PAD ? N_REF_PAD : tr->beg;

    tr->ref = faidx_fetch_seq(args->fai, chr, tr->beg - pad_beg, tr->end + N_REF_PAD, &len);
    if ( !tr->ref )
        error("faidx_fetch_seq failed %s:%d-%d\n", chr,tr->beg+1,tr->end+1);

    int pad_end = len - (tr->end - tr->beg + 1 + pad_beg);
    if ( pad_beg + pad_end != 2*N_REF_PAD )
    {
        char *ref = (char*) malloc(tr->end - tr->beg + 1 + 2*N_REF_PAD);
        for (i=0; i < N_REF_PAD - pad_beg; i++) ref[i] = 'N';
        memcpy(ref+i, tr->ref, len);
        for (i=0; i < N_REF_PAD - pad_end; i++) ref[i+len] = 'N';
        free(tr->ref);
        tr->ref = ref;
    }
}

static void sanity_check_ref(args_t *args, tscript_t *tr, bcf1_t *rec)
{
    char *ref = tr->ref + (rec->pos + N_REF_PAD >= tr->beg ? rec->pos - tr->beg + N_REF_PAD : 0);
    char *vcf = rec->d.allele[0] + (rec->pos + N_REF_PAD >= tr->beg ? 0 : tr->beg - N_REF_PAD - rec->pos);
    assert( vcf - rec->d.allele[0] < strlen(rec->d.allele[0]) );
    while ( *ref && *vcf )
    {
        if ( *ref!=*vcf && toupper(*ref)!=toupper(*vcf) ) 
            error("Error: the fasta reference does not match the VCF REF allele at %s:%d .. %s\n", bcf_seqname(args->hdr,rec),rec->pos+1,rec->d.allele[0]);
        ref++;
        vcf++;
    }
}

int test_cds_local(args_t *args, bcf1_t *rec)
{
    int i,j, ret = 0;
    const char *chr = bcf_seqname(args->hdr,rec);
    // note that the off-by-one extension of rlen is deliberate to account for insertions
    if ( !regidx_overlap(args->idx_cds,chr,rec->pos,rec->pos+rec->rlen, args->itr) ) return 0;

    // structures to fake the normal test_cds machinery
    hap_node_t root, node;
    root.type  = HAP_ROOT;
    kstring_t *tref = &args->hap->tref, *tseq = &args->hap->tseq;

    while ( regitr_overlap(args->itr) )
    {
        gf_cds_t *cds = regitr_payload(args->itr,gf_cds_t*);
        tscript_t *tr = cds->tr;
        if ( !GF_is_coding(tr->type) ) continue;
        ret = 1;

        if ( !tr->ref )
        {
            tscript_init_ref(args, tr, chr);
            tscript_splice_ref(tr);
            khp_insert(trhp, args->active_tr, &tr);     // only to clean the reference afterwards
        }

        sanity_check_ref(args, tr, rec);

        kstring_t sref;
        sref.s = tr->sref;
        sref.l = tr->nsref;
        sref.m = sref.l;

        for (i=1; i<rec->n_allele; i++)
        {
            if ( hap_init(args, &root, &node, cds, rec, i)!=0 ) continue;

            csq_t csq; 
            memset(&csq, 0, sizeof(csq_t));
            csq.pos          = rec->pos;
            csq.type.biotype = tr->type;
            csq.type.strand  = tr->strand;
            csq.type.trid    = tr->id;
            csq.type.gene    = tr->gene->name;

            int csq_type = node.csq;

            // code repetition: it would be nice to reuse the code from hap_add_csq, needs have refactoring though
            if ( node.type == HAP_SSS )
            {
                csq.type.type = csq_type;
                csq_stage(args, &csq, rec);
            }
            else
            {
                kstring_t sseq;
                sseq.m = sref.m - 2*N_REF_PAD + node.dlen;
                sseq.s = node.seq;
                int alen = sseq.l = strlen(sseq.s);
                int fill = node.dlen%3 && alen ? 1 : 0; // see #1475227917
                cds_translate(&sref, &sseq, node.sbeg,node.sbeg,node.sbeg+node.rlen, tr->strand, tseq, fill);

                sseq.m = sref.m - 2*N_REF_PAD;
                sseq.s = sref.s + N_REF_PAD + node.sbeg;
                sseq.l = node.rlen;
                cds_translate(&sref, &sseq, node.sbeg,node.sbeg,node.sbeg+node.rlen, tr->strand, tref, fill);

                // check for truncating stops
                for (j=0; j<tref->l; j++)
                    if ( tref->s[j]=='*' ) break;
                if ( j!=tref->l )
                {
                    tref->l = j+1;
                    tref->s[j+1] = 0;
                }
                for (j=0; j<tseq->l; j++)
                    if ( tseq->s[j]=='*' ) break;
                if ( j!=tseq->l )
                {
                    tseq->l = j+1;
                    tseq->s[j+1] = 0;
                }
                if ( csq_type & CSQ_STOP_LOST )
                {
                    if ( tref->s[tref->l-1]=='*' && tref->s[tref->l-1] == tseq->s[tseq->l-1] ) 
                    {
                        csq_type &= ~CSQ_STOP_LOST;
                        csq_type |= CSQ_STOP_RETAINED;
                    }
                    else if (tref->s[tref->l-1]!='*' )
                    {
                        // This is CDS 3' incomplete ENSG00000173376/synon.vcf, can also be missense
                        // We observe in real data a change to a stop, ENST00000528237/retained-stop-incomplete-cds.vcf
                        if ( tseq->s[tseq->l-1] == '*' )
                        {
                            csq_type &= ~CSQ_STOP_GAINED;
                            csq_type |= CSQ_STOP_RETAINED;
                        }
                        else
                            csq_type |= CSQ_INCOMPLETE_CDS;
                    }
                }
                if ( csq_type & CSQ_START_LOST && tref->s[0]!='M' )
                    csq_type &= ~CSQ_START_LOST;
                if ( node.dlen!=0 )
                {
                    if ( node.dlen%3 )
                        csq_type |= CSQ_FRAMESHIFT_VARIANT;
                    else if ( node.dlen<0 )
                        csq_type |= CSQ_INFRAME_DELETION;
                    else
                        csq_type |= CSQ_INFRAME_INSERTION;
                }
                else
                {
                    for (j=0; j<tref->l; j++) 
                        if ( tref->s[j] != tseq->s[j] ) break;
                    if ( j==tref->l )
                        csq_type |= CSQ_SYNONYMOUS_VARIANT;
                    else if ( tref->s[j] ==  '*' )
                        csq_type |= CSQ_STOP_LOST;
                    else if ( tseq->s[j] ==  '*' )
                        csq_type |= CSQ_STOP_GAINED;
                    else
                        csq_type |= CSQ_MISSENSE_VARIANT;
                }
                if ( csq_type & CSQ_COMPOUND )
                {
                    // create the aa variant string
                    kstring_t str = {0,0,0};
                    int aa_rbeg = tr->strand==STRAND_FWD ? node.sbeg/3+1 : (tr->nsref - 2*N_REF_PAD - node.sbeg - node.rlen)/3+1;
                    int aa_sbeg = tr->strand==STRAND_FWD ? node.sbeg/3+1 : (tr->nsref - 2*N_REF_PAD + node.dlen - node.sbeg - alen)/3+1;
                    kputc_('|', &str);
                    kputw(aa_rbeg, &str);
                    kputs(tref->s, &str);
                    if ( !(csq_type & CSQ_SYNONYMOUS_VARIANT) )
                    {
                        kputc_('>', &str);
                        kputw(aa_sbeg, &str);
                        kputs(tseq->s, &str);
                    }
                    kputc_('|', &str);
                    kputw(rec->pos+1, &str);
                    kputs(node.var, &str);
                    csq.type.vstr = str;
                    csq.type.type = csq_type & CSQ_COMPOUND;
                    csq_stage(args, &csq, rec);

                    // all this only to clean vstr when vrec is flushed
                    if ( !tr->root )
                        tr->root = (hap_node_t*) calloc(1,sizeof(hap_node_t));
                    tr->root->ncsq_list++;
                    hts_expand0(csq_t,tr->root->ncsq_list,tr->root->mcsq_list,tr->root->csq_list);
                    csq_t *rm_csq = tr->root->csq_list + tr->root->ncsq_list - 1;
                    rm_csq->type.vstr = str;            
                }
                if ( csq_type & ~CSQ_COMPOUND )
                {
                    csq.type.type = csq_type & ~CSQ_COMPOUND;
                    csq.type.vstr.l = 0;
                    csq_stage(args, &csq, rec);
                }
            }
            free(node.seq);
            free(node.var);
        }
    }
    return ret;
}

int test_cds(args_t *args, bcf1_t *rec)
{
    int i, ret = 0, hap_ret;
    const char *chr = bcf_seqname(args->hdr,rec);
    // note that the off-by-one extension of rlen is deliberate to account for insertions
    if ( !regidx_overlap(args->idx_cds,chr,rec->pos,rec->pos+rec->rlen, args->itr) ) return 0;
    while ( regitr_overlap(args->itr) )
    {
        gf_cds_t *cds = regitr_payload(args->itr,gf_cds_t*);
        tscript_t *tr = cds->tr;
        if ( !GF_is_coding(tr->type) ) continue;
        ret = 1;
        if ( !tr->root )
        {
            // initialize the transcript and its haplotype tree, fetch the reference sequence
            tscript_init_ref(args, tr, chr);

            tr->root = (hap_node_t*) calloc(1,sizeof(hap_node_t));
            tr->nhap = args->phase==PHASE_DROP_GT ? 1 : 2*args->smpl->n;     // maximum ploidy = diploid
            tr->hap  = (hap_node_t**) malloc(tr->nhap*sizeof(hap_node_t*));
            for (i=0; i<tr->nhap; i++) tr->hap[i] = NULL;
            tr->root->nend = tr->nhap;
            tr->root->type = HAP_ROOT;

            khp_insert(trhp, args->active_tr, &tr);
        }

        sanity_check_ref(args, tr, rec);

        if ( args->phase==PHASE_DROP_GT )
        {
            if ( rec->d.allele[1][0]=='<' || rec->d.allele[1][0]=='*' ) { continue; }
            hap_node_t *parent = tr->hap[0] ? tr->hap[0] : tr->root;
            hap_node_t *child  = (hap_node_t*)calloc(1,sizeof(hap_node_t));
            if ( (hap_ret=hap_init(args, parent, child, cds, rec, 1))!=0 )
            {
                // overlapping or intron variant, cannot apply
                if ( hap_ret==1 )
                {
                    if ( !args->quiet )
                        fprintf(stderr,"Warning: Skipping overlapping variants at %s:%d\t%s>%s\n", chr,rec->pos+1,rec->d.allele[0],rec->d.allele[1]);
                    if ( args->out ) 
                        fprintf(args->out,"LOG\tWarning: Skipping overlapping variants at %s:%d\t%s>%s\n", chr,rec->pos+1,rec->d.allele[0],rec->d.allele[1]);
                }
                else ret = 1;   // prevent reporting as intron in test_tscript
                free(child);
                continue;
            }
            parent->nend--;
            parent->nchild = 1;
            parent->mchild = 1;
            parent->child  = (hap_node_t**) malloc(sizeof(hap_node_t*));
            parent->child[0] = child;
            tr->hap[0] = child;
            tr->hap[0]->nend = 1;
            continue;
        }

        // apply the VCF variants and extend the haplotype tree
        int j, ismpl, ihap, ngts = bcf_get_genotypes(args->hdr, rec, &args->gt_arr, &args->mgt_arr);
        ngts /= bcf_hdr_nsamples(args->hdr);
        if ( ngts!=1 && ngts!=2 ) 
        {
            if ( !args->quiet )
                fprintf(stderr,"Warning: Skipping site with non-diploid/non-haploid genotypes at %s:%d\t%s>%s\n", chr,rec->pos+1,rec->d.allele[0],rec->d.allele[1]);
            if ( args->out ) 
                fprintf(args->out,"LOG\tWarning: Skipping site with non-diploid/non-haploid genotypes at %s:%d\t%s>%s\n", chr,rec->pos+1,rec->d.allele[0],rec->d.allele[1]);
            continue;
        }
        for (ismpl=0; ismpl<args->smpl->n; ismpl++)
        {
            int32_t *gt = args->gt_arr + args->smpl->idx[ismpl]*ngts;
            if ( gt[0]==bcf_gt_missing ) continue;

            if ( ngts>1 && gt[0]!=gt[1] && gt[1]!=bcf_gt_missing && gt[1]!=bcf_int32_vector_end )
            {
                if ( args->phase==PHASE_MERGE )
                {
                    if ( !bcf_gt_allele(gt[0]) ) gt[0] = gt[1];
                }
                if ( !bcf_gt_is_phased(gt[0]) && !bcf_gt_is_phased(gt[1]) )
                {
                    if ( args->phase==PHASE_REQUIRE )
                        error("Unphased genotype at %s:%d, sample %s. See the --phase option.\n", chr,rec->pos+1,args->hdr->samples[args->smpl->idx[ismpl]]);
                    if ( args->phase==PHASE_SKIP )
                        continue;
                    if ( args->phase==PHASE_NON_REF )
                    {
                        if ( !bcf_gt_allele(gt[0]) ) gt[0] = gt[1];
                        else if ( !bcf_gt_allele(gt[1]) ) gt[1] = gt[0];
                    }
                }
            }

            for (ihap=0; ihap<ngts; ihap++)
            {
                if ( gt[ihap]==bcf_gt_missing || gt[ihap]==bcf_int32_vector_end ) continue;

                i = 2*ismpl + ihap;

                int ial = bcf_gt_allele(gt[ihap]);
                if ( !ial ) continue;
                assert( ial < rec->n_allele );
                if ( rec->d.allele[ial][0]=='<' || rec->d.allele[ial][0]=='*' ) { continue; }

                hap_node_t *parent = tr->hap[i] ? tr->hap[i] : tr->root;
                if ( parent->cur_rec==rec && parent->cur_child[ial]>=0 )
                {
                    // this haplotype has been seen in another sample
                    tr->hap[i] = parent->child[ parent->cur_child[ial] ];
                    tr->hap[i]->nend++;
                    parent->nend--;
                    continue;
                }

                hap_node_t *child = (hap_node_t*)calloc(1,sizeof(hap_node_t));
                if ( (hap_ret=hap_init(args, parent, child, cds, rec, ial))!=0 )
                {
                    // overlapping or intron variant, cannot apply
                    if ( hap_ret==1 )
                    {
                        if ( !args->quiet )
                            fprintf(stderr,"Warning: Skipping overlapping variants at %s:%d, sample %s\t%s>%s\n",
                                    chr,rec->pos+1,args->hdr->samples[args->smpl->idx[ismpl]],rec->d.allele[0],rec->d.allele[ial]);
                        if ( args->out  )
                            fprintf(args->out,"LOG\tWarning: Skipping overlapping variants at %s:%d, sample %s\t%s>%s\n",
                                    chr,rec->pos+1,args->hdr->samples[args->smpl->idx[ismpl]],rec->d.allele[0],rec->d.allele[ial]);
                    }
                    free(child);
                    continue;
                }

                if ( parent->cur_rec!=rec )
                {
                    hts_expand(int,rec->n_allele,parent->mcur_child,parent->cur_child);
                    for (j=0; j<rec->n_allele; j++) parent->cur_child[j] = -1;
                    parent->cur_rec = rec;
                }

                j = parent->nchild++;
                hts_expand0(hap_node_t*,parent->nchild,parent->mchild,parent->child);
                parent->cur_child[ial] = j;
                parent->child[j] = child;
                tr->hap[i] = child;
                tr->hap[i]->nend++;
                parent->nend--;
            }
        }
    }
    return ret;
}

void csq_stage(args_t *args, csq_t *csq, bcf1_t *rec)
{
    // known issues: tab output leads to unsorted output. This is because
    // coding haplotypes are printed in one go and buffering is not used
    // with tab output. VCF output is OK though.
    if ( csq_push(args, csq, rec)!=0 ) return;    // the consequence already exists

    int i,j,ngt = 0;
    if ( args->phase!=PHASE_DROP_GT )
    {
        ngt = bcf_get_genotypes(args->hdr, rec, &args->gt_arr, &args->mgt_arr);
        if ( ngt>0 ) ngt /= bcf_hdr_nsamples(args->hdr);
    }
    if ( ngt<=0 )
    {
        if ( args->output_type==FT_TAB_TEXT )
            csq_print_text(args, csq, -1,0);
        return;
    }
    assert( ngt<=2 );

    if ( args->output_type==FT_TAB_TEXT )
    {
        for (i=0; i<args->smpl->n; i++)
        {
            int32_t *gt = args->gt_arr + args->smpl->idx[i]*ngt;
            for (j=0; j<ngt; j++)
            {
                if ( gt[j]==bcf_gt_missing || gt[j]==bcf_int32_vector_end || !bcf_gt_allele(gt[j]) ) continue;
                csq_print_text(args, csq, args->smpl->idx[i],j+1);
            }
        }
        return;
    }

    vrec_t *vrec = csq->vrec;
    for (i=0; i<args->smpl->n; i++)
    {
        int32_t *gt = args->gt_arr + args->smpl->idx[i]*ngt;
        for (j=0; j<ngt; j++)
        {
            if ( gt[j]==bcf_gt_missing || gt[j]==bcf_int32_vector_end || !bcf_gt_allele(gt[j]) ) continue;

            int icsq = 2*csq->idx + j;
            if ( icsq >= args->ncsq_max ) // more than ncsq_max consequences, so can't fit it in FMT
            {
                int ismpl = args->smpl->idx[i];
                int print_warning = 1;
                if ( args->quiet )
                {
                    if ( args->quiet > 1 || args->ncsq_small_warned ) print_warning = 0;
                    args->ncsq_small_warned = 1;
                }
                if ( print_warning )
                {
                    fprintf(stderr,"Warning: --ncsq %d is too small to annotate %s at %s:%d with %d-th csq\n",
                            args->ncsq_max/2,args->hdr->samples[ismpl],bcf_hdr_id2name(args->hdr,args->rid),vrec->line->pos+1,csq->idx+1);
                    if ( args->quiet ) fprintf(stderr,"(This warning is printed only once)\n");
                }
                break;
            }
            if ( vrec->nfmt < 1 + icsq/32 ) vrec->nfmt = 1 + icsq/32;
            vrec->smpl[i*args->nfmt_bcsq + icsq/32] |= 1 << (icsq % 32);
        }
    }
}
int test_utr(args_t *args, bcf1_t *rec)
{
    const char *chr = bcf_seqname(args->hdr,rec);
    // note that the off-by-one extension of rlen is deliberate to account for insertions
    if ( !regidx_overlap(args->idx_utr,chr,rec->pos,rec->pos+rec->rlen, args->itr) ) return 0;

    splice_t splice;
    splice_init(&splice, rec);

    int i, ret = 0;
    while ( regitr_overlap(args->itr) )
    {
        gf_utr_t *utr = regitr_payload(args->itr, gf_utr_t*);
        tscript_t *tr = splice.tr = utr->tr;
        for (i=1; i<rec->n_allele; i++)
        {
            if ( rec->d.allele[1][0]=='<' || rec->d.allele[1][0]=='*' ) { continue; }
            splice.vcf.alt = rec->d.allele[i];
            int splice_ret = splice_csq(args, &splice, utr->beg, utr->end);
            if ( splice_ret!=SPLICE_INSIDE && splice_ret!=SPLICE_OVERLAP ) continue;
            csq_t csq; 
            memset(&csq, 0, sizeof(csq_t));
            csq.pos          = rec->pos;
            csq.type.type    = utr->which==prime5 ? CSQ_UTR5 : CSQ_UTR3;
            csq.type.biotype = tr->type;
            csq.type.strand  = tr->strand;
            csq.type.trid    = tr->id;
            csq.type.gene    = tr->gene->name;
            csq_stage(args, &csq, rec);
            ret = 1;
        }
    }
    assert(!splice.kref.s);
    assert(!splice.kalt.s);
    return ret;
}
int test_splice(args_t *args, bcf1_t *rec)
{
    const char *chr = bcf_seqname(args->hdr,rec);
    if ( !regidx_overlap(args->idx_exon,chr,rec->pos,rec->pos + rec->rlen, args->itr) ) return 0;

    splice_t splice;
    splice_init(&splice, rec);
    splice.check_acceptor = splice.check_donor = 1;

    int i, ret = 0;
    while ( regitr_overlap(args->itr) )
    {
        gf_exon_t *exon = regitr_payload(args->itr, gf_exon_t*);
        splice.tr = exon->tr;
        if ( !splice.tr->ncds ) continue;  // not a coding transcript, no interest in splice sites

        splice.check_region_beg = splice.tr->beg==exon->beg ? 0 : 1;
        splice.check_region_end = splice.tr->end==exon->end ? 0 : 1;

        for (i=1; i<rec->n_allele; i++)
        {
            if ( rec->d.allele[1][0]=='<' || rec->d.allele[1][0]=='*' ) { continue; }
            splice.vcf.alt = rec->d.allele[i];
            splice_csq(args, &splice, exon->beg, exon->end);
            if ( splice.csq ) ret = 1;
        }
    }
    free(splice.kref.s);
    free(splice.kalt.s);
    return ret;
}
int test_tscript(args_t *args, bcf1_t *rec)
{
    const char *chr = bcf_seqname(args->hdr,rec);
    if ( !regidx_overlap(args->idx_tscript,chr,rec->pos,rec->pos+rec->rlen, args->itr) ) return 0;

    splice_t splice;
    splice_init(&splice, rec);

    int i, ret = 0;
    while ( regitr_overlap(args->itr) )
    {
        tscript_t *tr = splice.tr = regitr_payload(args->itr, tscript_t*);
        for (i=1; i<rec->n_allele; i++)
        {
            if ( rec->d.allele[1][0]=='<' || rec->d.allele[1][0]=='*' ) { continue; }
            splice.vcf.alt = rec->d.allele[i];
            int splice_ret = splice_csq(args, &splice, tr->beg, tr->end);
            if ( splice_ret!=SPLICE_INSIDE && splice_ret!=SPLICE_OVERLAP ) continue;    // SPLICE_OUTSIDE or SPLICE_REF
            csq_t csq; 
            memset(&csq, 0, sizeof(csq_t));
            csq.pos          = rec->pos;
            csq.type.type    = GF_is_coding(tr->type) ? CSQ_INTRON : CSQ_NON_CODING;
            csq.type.biotype = tr->type;
            csq.type.strand  = tr->strand;
            csq.type.trid    = tr->id;
            csq.type.gene    = tr->gene->name;
            csq_stage(args, &csq, rec);
            ret = 1;
        }
    }
    assert(!splice.kref.s);
    assert(!splice.kalt.s);
    return ret;
}

void process(args_t *args, bcf1_t **rec_ptr)
{
    if ( !rec_ptr )
    {
        hap_flush(args, REGIDX_MAX);
        vbuf_flush(args);
        return;
    }

    bcf1_t *rec = *rec_ptr;

    int call_csq = 1;
    if ( !rec->n_allele ) call_csq = 0;   // no alternate allele
    else if ( rec->n_allele==2 && (rec->d.allele[1][0]=='<' || rec->d.allele[1][0]=='*') ) call_csq = 0;     // gVCF, no alt allele
    else if ( args->filter )
    {
        call_csq = filter_test(args->filter, rec, NULL);
        if ( args->filter_logic==FLT_EXCLUDE ) call_csq = call_csq ? 0 : 1;
    }
    if ( !call_csq )
    {
        if ( !args->out_fh ) return;    // not a VCF output
        vbuf_push(args, rec_ptr);
        vbuf_flush(args);
        return;
    }

    if ( args->rid != rec->rid ) 
    {
        hap_flush(args, REGIDX_MAX);
        vbuf_flush(args);
    }
    args->rid = rec->rid;
    vbuf_push(args, rec_ptr);

    int hit = args->local_csq ? test_cds_local(args, rec) : test_cds(args, rec);
    hit += test_utr(args, rec);
    hit += test_splice(args, rec);
    if ( !hit ) test_tscript(args, rec);

    hap_flush(args, rec->pos-1);
    vbuf_flush(args);

    return;
}

static const char *usage(void)
{
    return 
        "\n"
        "About: Haplotype-aware consequence caller.\n"
        "Usage: bcftools csq [options] in.vcf\n"
        "\n"
        "Required options:\n"
        "   -f, --fasta-ref <file>          reference file in fasta format\n"
        "   -g, --gff-annot <file>          gff3 annotation file\n"
        "\n"
        "CSQ options:\n"
        "   -c, --custom-tag <string>       use this tag instead of the default BCSQ\n"
        "   -l, --local-csq                 localized predictions, consider only one VCF record at a time\n"
        "   -n, --ncsq <int>                maximum number of consequences to consider per site [16]\n"
        "   -p, --phase <a|m|r|R|s>         how to construct haplotypes and how to deal with unphased data: [r]\n"
        "                                     a: take GTs as is, create haplotypes regardless of phase (0/1 -> 0|1)\n"
        "                                     m: merge *all* GTs into a single haplotype (0/1 -> 1, 1/2 -> 1)\n"
        "                                     r: require phased GTs, throw an error on unphased het GTs\n"
        "                                     R: create non-reference haplotypes if possible (0/1 -> 1|1, 1/2 -> 1|2)\n"
        "                                     s: skip unphased GTs\n"
        "Options:\n"
        "   -e, --exclude <expr>            exclude sites for which the expression is true\n"
        "   -i, --include <expr>            select sites for which the expression is true\n"
        "   -o, --output <file>             write output to a file [standard output]\n"
        "   -O, --output-type <b|u|z|v|t>   b: compressed BCF, u: uncompressed BCF, z: compressed VCF\n"
        "                                   v: uncompressed VCF, t: plain tab-delimited text output [v]\n"
        "   -q, --quiet                     suppress warning messages. Can be given two times for even less messages\n"
        "   -r, --regions <region>          restrict to comma-separated list of regions\n"
        "   -R, --regions-file <file>       restrict to regions listed in a file\n"
        "   -s, --samples <-|list>          samples to include or \"-\" to apply all variants and ignore samples\n"
        "   -S, --samples-file <file>       samples to include\n"
        "   -t, --targets <region>          similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file <file>       similar to -R but streams rather than index-jumps\n"
        "\n"
        "Example:\n"
        "   bcftools csq -f hs37d5.fa -g Homo_sapiens.GRCh37.82.gff3.gz in.vcf\n"
        "\n"
        "   # GFF3 annotation files can be downloaded from Ensembl. e.g. for human:\n"
        "   ftp://ftp.ensembl.org/pub/current_gff3/homo_sapiens/\n"
        "   ftp://ftp.ensembl.org/pub/grch37/release-84/gff3/homo_sapiens/\n"
        "\n";
}

int main_csq(int argc, char *argv[])
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc = argc; args->argv = argv;
    args->output_type = FT_VCF;
    args->bcsq_tag = "BCSQ";
    args->ncsq_max = 2*16;

    static struct option loptions[] =
    {
        {"help",0,0,'h'},
        {"ncsq",1,0,'n'},
        {"custom-tag",1,0,'c'},
        {"local-csq",0,0,'l'},
        {"gff-annot",1,0,'g'},
        {"fasta-ref",1,0,'f'},
        {"include",1,0,'i'},
        {"exclude",1,0,'e'},
        {"output",1,0,'o'},
        {"output-type",1,NULL,'O'},
        {"phase",1,0,'p'},
        {"quiet",0,0,'q'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {0,0,0,0}
    };
    int c, targets_is_file = 0, regions_is_file = 0; 
    char *targets_list = NULL, *regions_list = NULL;
    while ((c = getopt_long(argc, argv, "?hr:R:t:T:i:e:f:o:O:g:s:S:p:qc:ln:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'l': args->local_csq = 1; break;
            case 'c': args->bcsq_tag = optarg; break;
            case 'q': args->quiet++; break;
            case 'p':
                switch (optarg[0]) 
                {
                    case 'a': args->phase = PHASE_AS_IS; break;
                    case 'm': args->phase = PHASE_MERGE; break;
                    case 'r': args->phase = PHASE_REQUIRE; break;
                    case 'R': args->phase = PHASE_NON_REF; break;
                    case 's': args->phase = PHASE_SKIP; break;
                    default: error("The -p code \"%s\" not recognised\n", optarg);
                }
                break;
            case 'f': args->fa_fname = optarg; break;
            case 'g': args->gff_fname = optarg; break;
            case 'n': 
                args->ncsq_max = 2 * atoi(optarg);
                if ( args->ncsq_max <=0 ) error("Expected positive integer with -n, got %s\n", optarg);
                break;
            case 'o': args->output_fname = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 't': args->output_type = FT_TAB_TEXT; break;
                          case 'b': args->output_type = FT_BCF_GZ; break;
                          case 'u': args->output_type = FT_BCF; break;
                          case 'z': args->output_type = FT_VCF_GZ; break;
                          case 'v': args->output_type = FT_VCF; break;
                          default: error("The output type \"%s\" not recognised\n", optarg);
                      }
                      break;
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'r': regions_list = optarg; break;
            case 'R': regions_list = optarg; regions_is_file = 1; break;
            case 's': args->sample_list = optarg; break;
            case 'S': args->sample_list = optarg; args->sample_is_file = 1; break;
            case 't': targets_list = optarg; break;
            case 'T': targets_list = optarg; targets_is_file = 1; break;
            case 'h':
            case '?': error("%s",usage());
            default: error("The option not recognised: %s\n\n", optarg); break;
        }
    }
    char *fname = NULL;
    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin
        else error("%s", usage());
    }
    else fname = argv[optind];
    if ( argc - optind>1 ) error("%s", usage());
    if ( !args->fa_fname ) error("Missing the --fa-ref option\n");
    if ( !args->gff_fname ) error("Missing the --gff option\n");
    args->sr = bcf_sr_init();
    if ( targets_list && bcf_sr_set_targets(args->sr, targets_list, targets_is_file, 0)<0 )
        error("Failed to read the targets: %s\n", targets_list);
    if ( regions_list && bcf_sr_set_regions(args->sr, regions_list, regions_is_file)<0 )
        error("Failed to read the regions: %s\n", regions_list);
    if ( !bcf_sr_add_reader(args->sr, fname) )
        error("Failed to open %s: %s\n", fname,bcf_sr_strerror(args->sr->errnum));
    args->hdr = bcf_sr_get_header(args->sr,0);

    init_data(args);
    while ( bcf_sr_next_line(args->sr) )
    {
        process(args, &args->sr->readers[0].buffer[0]);
    }
    process(args,NULL);

    destroy_data(args);
    bcf_sr_destroy(args->sr);
    free(args);

    return 0;
}

