/* The MIT License

   Copyright (c) 2016-2024 Genome Research Ltd.

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
        - dynamic N_REF_PAD
        - for stop-lost events (also in frameshifts) report the number of truncated aa's
        - memory could be greatly reduced by indexing gff (but it is quite compact already)
        - deletions that go beyond transcript boundaries are not checked at sequence level
            - alloc tscript->ref in hap_finalize, introduce fa_off_beg:16,fa_off_end:16
            - see test/csq/ENST00000573314/insertion-overlap.vcf #1476288882

    Read about transcript types here
        http://vega.sanger.ac.uk/info/about/gene_and_transcript_types.html
        https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
        https://www.gencodegenes.org/pages/biotypes.html

    List of supported biotypes
        antisense
        IG_C_gene
        IG_D_gene
        IG_J_gene
        IG_LV_gene
        IG_V_gene
        lincRNA
        lncRNA      .. generic term for 3prime_overlapping_ncRNA, antisense, bidirectional_promoter_lncRNA, lincRNA, macro_lncRNA, non_coding, processed_transcript, sense_intronic, sense_overlapping
        macro_lncRNA
        miRNA
        misc_RNA
        Mt_rRNA
        Mt_tRNA
        polymorphic_pseudogene
        processed_transcript
        protein_coding, mRNA
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

            B .. transcript line referencing A with supported biotype
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
        start_retained_variant  .. start codon retained by indel realignment
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
#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <inttypes.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/khash.h>
#include <htslib/khash_str2int.h>
#include <htslib/kseq.h>
#include <htslib/faidx.h>
#include <htslib/bgzf.h>
#include <errno.h>
#include <unistd.h>
#include <ctype.h>
#include <strings.h>
#include "bcftools.h"
#include "filter.h"
#include "regidx.h"
#include "kheap.h"
#include "smpl_ilist.h"
#include "rbuf.h"
#include "gff.h"

#ifndef __FUNCTION__
#  define __FUNCTION__ __func__
#endif

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define N_REF_PAD 10    // number of bases to avoid boundary effects

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
#define CSQ_ELONGATION          (1<<22)     // symbolic insertion
#define CSQ_START_RETAINED      (1<<23)

// Haplotype-aware consequences, printed in one vcf record only, the rest has a reference @12345
#define CSQ_COMPOUND (CSQ_SYNONYMOUS_VARIANT|CSQ_MISSENSE_VARIANT|CSQ_STOP_LOST|CSQ_STOP_GAINED| \
                      CSQ_INFRAME_DELETION|CSQ_INFRAME_INSERTION|CSQ_FRAMESHIFT_VARIANT| \
                      CSQ_START_LOST|CSQ_STOP_RETAINED|CSQ_INFRAME_ALTERING|CSQ_INCOMPLETE_CDS| \
                      CSQ_UPSTREAM_STOP|CSQ_START_RETAINED)
#define CSQ_START_STOP          (CSQ_STOP_LOST|CSQ_STOP_GAINED|CSQ_STOP_RETAINED|CSQ_START_LOST|CSQ_START_RETAINED)

#define CSQ_PRN_STRAND(csq)     ((csq)&CSQ_COMPOUND && !((csq)&(CSQ_SPLICE_ACCEPTOR|CSQ_SPLICE_DONOR|CSQ_SPLICE_REGION)))
#define CSQ_PRN_TSCRIPT         (~(CSQ_INTRON|CSQ_NON_CODING))
#define CSQ_PRN_NMD             (~(CSQ_INTRON|CSQ_NON_CODING))
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
    "coding_sequence",
    "feature_elongation",
    "start_retained"
};

/*
    Structures related to VCF output:

    vcsq_t
        information required to assemble consequence lines such as "inframe_deletion|XYZ|ENST01|+|5TY>5I|121ACG>A+124TA>T"

    vrec_t
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
    uint32_t vcf_ial;
    uint32_t biotype;   // one of GF_* types
    char *gene;         // gene name
    bcf1_t *ref;        // if type&CSQ_PRINTED_UPSTREAM, ref consequence "@1234"
    kstring_t vstr;     // variant string, eg 5TY>5I|121ACG>A+124TA>T
};
typedef struct
{
    bcf1_t *line;
    uint32_t *fmt_bm;   // bitmask of sample consequences with first/second haplotype interleaved
    uint32_t nfmt:4,    // the bitmask size (the number of integers per sample)
             nvcsq:28, mvcsq;
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
    uint32_t keep_until;    // the maximum transcript end position
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
    int vcf_ial;                // which VCF allele generated this node
    uint32_t nend;              // number of haplotypes ending in this node
    int *cur_child, mcur_child; // mapping from the allele to the currently active child
    csq_t *csq_list;            // list of haplotype's consequences, broken by position (each corresponds to a VCF record)
    int ncsq_list, mcsq_list;
};
#define TSCRIPT_AUX(x) ((tscript_t*)(x)->aux)
typedef struct
{
    char *ref;          // reference sequence, padded with N_REF_PAD bases on both ends
    char *sref;         // spliced reference sequence, padded with N_REF_PAD bases on both ends
    hap_node_t *root;   // root of the haplotype tree
    hap_node_t **hap;   // pointer to haplotype leaves, two for each sample
    int nhap, nsref;    // number of haplotypes and length of sref, including 2*N_REF_PAD
}
tscript_t;
static inline int cmp_tscript(gf_tscript_t **a, gf_tscript_t **b)
{
    return ( (*a)->end  < (*b)->end ) ? 1 : 0;
}
KHEAP_INIT(trhp, gf_tscript_t*, cmp_tscript)
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
    gf_tscript_t *tr;   // tr->ref: spliced transcript on ref strand
    kstring_t sseq;     // spliced haplotype sequence on ref strand
    kstring_t tseq;     // the variable part of translated haplotype transcript, coding strand
    kstring_t tref;     // the variable part of translated reference transcript, coding strand
    uint32_t sbeg;      // stack's sbeg, for cases first node's type is HAP_SSS
    int upstream_stop;
}
hap_t;

typedef struct _args_t
{
    // the main regidx lookups, from chr:beg-end to overlapping features and
    // index iterator
    gff_t *gff;
    regidx_t *idx_cds, *idx_utr, *idx_exon, *idx_tscript;
    regitr_t *itr;

    // text tab-delimited output (out) or vcf/bcf output (out_fh)
    FILE *out;
    htsFile *out_fh;
    char *index_fn;
    int write_index;
    char *dump_gff;

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
    int argc, output_type, clevel;
    int phase, verbosity, local_csq, record_cmd_line;
    int ncsq2_max, nfmt_bcsq;   // maximum number of csq per site that can be accessed from FORMAT/BCSQ (*2 and 1 bit skipped to avoid BCF missing values)
    int ncsq2_small_warned;
    int brief_predictions;
    int unify_chr_names;
    char *chr_name;
    int mchr_name;
    struct {
        int unknown_chr,unknown_tscript_biotype,unknown_strand,unknown_phase,duplicate_id;
        int unknown_cds_phase,incomplete_cds,wrong_phase,overlapping_cds;
    } warned;

    int rid;                    // current chromosome
    tr_heap_t *active_tr;       // heap of active transcripts for quick flushing
    hap_t *hap;                 // transcript haplotype recursion
    vbuf_t **vcf_buf;           // buffered VCF lines to annotate with CSQ and flush
    rbuf_t vcf_rbuf;            // round buffer indexes to vcf_buf
    kh_pos2vbuf_t *pos2vbuf;    // fast lookup of buffered lines by position
    gf_tscript_t **rm_tr;       // buffer of transcripts to clean
    int nrm_tr, mrm_tr;
    csq_t *csq_buf;             // pool of csq not managed by hap_node_t, i.e. non-CDS csqs
    int ncsq_buf, mcsq_buf;
    int force;                  // force run under various conditions. Currently only to skip out-of-phase transcripts
    int n_threads;              // extra compression/decompression threads

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

static inline int ncsq2_to_nfmt(int ncsq2)
{
    return 1 + (ncsq2 - 1) / 30;
}
static inline void icsq2_to_bit(int icsq2, int *ival, int *ibit)
{
    *ival = icsq2 / 30;
    *ibit = icsq2 % 30;
}

void init_data(args_t *args)
{
    args->nfmt_bcsq = ncsq2_to_nfmt(args->ncsq2_max);

    args->fai = fai_load(args->fa_fname);
    if ( !args->fai ) error("Failed to load the fai index: %s\n", args->fa_fname);

    args->gff = gff_init(args->gff_fname);
    gff_set(args->gff,verbosity,args->verbosity);
    gff_set(args->gff,strip_chr_names,args->unify_chr_names);
    gff_set(args->gff,force_out_of_phase,args->force);
    gff_set(args->gff,dump_fname,args->dump_gff);
    gff_parse(args->gff);
    args->idx_cds  = gff_get(args->gff,idx_cds);
    args->idx_utr  = gff_get(args->gff,idx_utr);
    args->idx_exon = gff_get(args->gff,idx_exon);
    args->idx_tscript = gff_get(args->gff,idx_tscript);
    args->itr = regitr_init(NULL);

    args->rid = -1;

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);

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
            if (bcf_hdr_set_samples(args->hdr,NULL,0) < 0)
                error_errno("[%s] Couldn't build sample filter", __func__);
        }
        args->phase = PHASE_DROP_GT;
    }
    else
        args->smpl = smpl_ilist_init(args->hdr, args->sample_list, args->sample_is_file, SMPL_STRICT);
    args->hdr_nsmpl = args->phase==PHASE_DROP_GT ? 0 : bcf_hdr_nsamples(args->hdr);

    if ( args->output_type==FT_TAB_TEXT )
    {
        args->out = args->output_fname ? fopen(args->output_fname,"w") : stdout;
        if ( !args->out ) error("Failed to write to %s: %s\n", !strcmp("-",args->output_fname)?"standard output":args->output_fname,strerror(errno));

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
        char wmode[8];
        set_wmode(wmode,args->output_type,args->output_fname,args->clevel);
        args->out_fh = hts_open(args->output_fname ? args->output_fname : "-", wmode);
        if ( args->out_fh == NULL ) error("[%s] Error: cannot write to %s: %s\n", __func__,args->output_fname? args->output_fname : "standard output", strerror(errno));
        if ( args->n_threads > 0)
            hts_set_opt(args->out_fh, HTS_OPT_THREAD_POOL, args->sr->p);
        if ( args->record_cmd_line ) bcf_hdr_append_version(args->hdr,args->argc,args->argv,"bcftools/csq");
        bcf_hdr_printf(args->hdr,"##INFO=<ID=%s,Number=.,Type=String,Description=\"%s consequence annotation from BCFtools/csq, see http://samtools.github.io/bcftools/howtos/csq-calling.html for details. Format: Consequence|gene|transcript|biotype|strand|amino_acid_change|dna_change\">",args->bcsq_tag, args->local_csq ? "Local" : "Haplotype-aware");
        if ( args->hdr_nsmpl )
            bcf_hdr_printf(args->hdr,"##FORMAT=<ID=%s,Number=.,Type=Integer,Description=\"Bitmask of indexes to INFO/BCSQ, with interleaved first/second haplotype. Use \\\"bcftools query -f'[%%CHROM\\t%%POS\\t%%SAMPLE\\t%%TBCSQ\\n]'\\\" to translate.\">",args->bcsq_tag);
        if ( bcf_hdr_write(args->out_fh, args->hdr)!=0 ) error("[%s] Error: cannot write the header to %s\n", __func__,args->output_fname?args->output_fname:"standard output");
        if ( init_index2(args->out_fh,args->hdr,args->output_fname,
                         &args->index_fn, args->write_index) < 0 )
            error("Error: failed to initialise index for %s\n",args->output_fname);
    }
    if ( args->verbosity > 0 ) fprintf(stderr,"Calling...\n");
}

void destroy_data(args_t *args)
{
    if ( args->ncsq2_small_warned )
        fprintf(stderr,
            "Note: Some samples had too many consequences to be represented in %d bytes. If you need to record them all,\n"
            "      the limit can be increased by running with `--ncsq %d`.\n",ncsq2_to_nfmt(args->ncsq2_max)/8,1+args->ncsq2_small_warned/2);

    regitr_destroy(args->itr);
    gff_destroy(args->gff);

    if ( args->filter )
        filter_destroy(args->filter);

    khp_destroy(trhp,args->active_tr);
    kh_destroy(pos2vbuf,args->pos2vbuf);
    if ( args->smpl ) smpl_ilist_destroy(args->smpl);
    int i,j,ret;
    if ( args->out_fh )
    {
        if ( args->write_index )
        {
            if ( bcf_idx_save(args->out_fh)<0 )
            {
                if ( hts_close(args->out_fh)!=0 ) error("Error: close failed .. %s\n", args->output_fname?args->output_fname:"stdout");
                error("Error: cannot write to index %s\n", args->index_fn);
            }
            free(args->index_fn);
        }
        ret = hts_close(args->out_fh);
    }
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
            free(vbuf->vrec[j]->fmt_bm);
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
    free(args->chr_name);
}

/*
    The splice_* functions are for consequences around splice sites: start,stop,splice_*
 */
#define SPLICE_VAR_REF 0   // ref: ACGT>ACGT, csq not applicable, skip completely
#define SPLICE_OUTSIDE 1   // splice acceptor or similar; csq set and is done, does not overlap the region
#define SPLICE_INSIDE  2   // overlaps coding region; csq can be set but coding prediction is needed
#define SPLICE_OVERLAP 3   // indel overlaps region boundary, csq set but could not determine csq
typedef struct
{
    gf_tscript_t *tr;
    struct {
        int32_t pos, rlen, alen, ial;
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
    splice->csq      = 0;
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
        kputsn(TSCRIPT_AUX(splice->tr)->ref + N_REF_PAD + rbeg - splice->tr->beg, splice->vcf.pos - rbeg, &splice->kref);
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
            kputsn(TSCRIPT_AUX(splice->tr)->ref + N_REF_PAD + end - splice->tr->beg, rlen - splice->kref.l, &splice->kref);
    }
#if XDBG
fprintf(stderr,"r3: %s\n",splice->kref.s);
#endif


    int aoff;
    if ( abeg < splice->vcf.pos )
    {
        assert( splice->tr->beg <= abeg );
        kputsn(TSCRIPT_AUX(splice->tr)->ref + N_REF_PAD + abeg - splice->tr->beg, splice->vcf.pos - abeg, &splice->kalt);
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
            kputsn(TSCRIPT_AUX(splice->tr)->ref + aoff + N_REF_PAD + end - splice->tr->beg, alen - splice->kalt.l, &splice->kalt);
    }
#if XDBG
fprintf(stderr,"a3: %s\n",splice->kalt.s);
fprintf(stderr," [%s]\n [%s]\n\n",splice->kref.s,splice->kalt.s);
#endif
}
void csq_stage(args_t *args, csq_t *csq, bcf1_t *rec);
static inline int csq_stage_utr(args_t *args, regitr_t *itr, bcf1_t *rec, uint32_t trid, uint32_t type, int ial)
{
    while ( regitr_overlap(itr) )
    {
        gf_utr_t *utr = regitr_payload(itr, gf_utr_t*);
        gf_tscript_t *tr = utr->tr;
        if ( tr->id != trid ) continue;
        csq_t csq;
        memset(&csq, 0, sizeof(csq_t));
        csq.pos          = rec->pos;
        csq.type.type    = (utr->which==prime5 ? CSQ_UTR5 : CSQ_UTR3) | type;
        csq.type.biotype = tr->type;
        csq.type.strand  = tr->strand;
        csq.type.trid    = tr->id;
        csq.type.vcf_ial = ial;
        csq.type.gene    = tr->gene->name;
        csq_stage(args, &csq, rec);
        return csq.type.type;
    }
    return 0;
}
static inline void csq_stage_splice(args_t *args, bcf1_t *rec, gf_tscript_t *tr, uint32_t type, int ial)
{
#if XDBG
fprintf(stderr,"csq_stage_splice %d: type=%d\n",(int)rec->pos+1,type);
#endif
    if ( !type ) return;
    csq_t csq;
    memset(&csq, 0, sizeof(csq_t));
    csq.pos          = rec->pos;
    csq.type.type    = type;
    csq.type.biotype = tr->type;
    csq.type.strand  = tr->strand;
    csq.type.trid    = tr->id;
    csq.type.vcf_ial = ial;
    csq.type.gene    = tr->gene->name;
    csq_stage(args, &csq, rec);
}
static inline const char *drop_chr_prefix(args_t *args, const char *chr)
{
    if ( !args->unify_chr_names ) return chr;
    if ( !strncasecmp("chr",chr,3) ) return chr+3;
    return chr;
}
static inline const char *add_chr_prefix(args_t *args, const char *chr)
{
    if ( !args->unify_chr_names ) return chr;
    int len = strlen(chr);
    hts_expand(char,len+4,args->mchr_name,args->chr_name);
    memcpy(args->chr_name,"chr",3);
    memcpy(args->chr_name+3,chr,len+1);
    return args->chr_name;
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
            const char *chr = drop_chr_prefix(args, bcf_seqname(args->hdr,splice->vcf.rec));
            if ( regidx_overlap(args->idx_utr,chr,splice->ref_beg+1,splice->ref_beg+1, itr) )     // adjacent utr
            {
                ret = csq_stage_utr(args, itr, splice->vcf.rec, splice->tr->id, splice->csq, splice->vcf.ial);
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
        csq_stage_splice(args, splice->vcf.rec, splice->tr, splice->csq, splice->vcf.ial);
        return SPLICE_OUTSIDE;
    }
    if ( splice->ref_end < ex_beg || (splice->ref_end == ex_beg && !splice->check_region_beg) )    // fully outside, before the exon
    {
        if ( splice->check_utr )
        {
            regitr_t *itr = regitr_init(NULL);
            const char *chr = drop_chr_prefix(args, bcf_seqname(args->hdr,splice->vcf.rec));
            if ( regidx_overlap(args->idx_utr,chr,splice->ref_end-1,splice->ref_end-1, itr) )     // adjacent utr
            {
                ret = csq_stage_utr(args, itr, splice->vcf.rec, splice->tr->id, splice->csq, splice->vcf.ial);
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
        csq_stage_splice(args, splice->vcf.rec, splice->tr, splice->csq, splice->vcf.ial);
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
    csq_stage_splice(args, splice->vcf.rec, splice->tr, splice->csq, splice->vcf.ial);
    return SPLICE_INSIDE;
}

int shifted_del_synonymous(args_t *args, splice_t *splice, uint32_t ex_beg, uint32_t ex_end)
{
    static int small_ref_padding_warned = 0;
    gf_tscript_t *tr = splice->tr;

    // We know the VCF record overlaps the exon, but does it overlap the start codon?
    if ( tr->strand==STRAND_REV && splice->vcf.pos + splice->vcf.rlen + 2 <= ex_end ) return 0;
    if ( tr->strand==STRAND_FWD && splice->vcf.pos >= ex_beg + 3 ) return 0;

#if XDBG
    fprintf(stderr,"shifted_del_synonymous: %d-%d  %s\n",ex_beg,ex_end, tr->strand==STRAND_FWD?"fwd":(tr->strand==STRAND_REV?"rev":"unk"));
    fprintf(stderr,"   %d  ..  %s > %s\n",splice->vcf.pos+1,splice->vcf.ref,splice->vcf.alt);
#endif

    // is there enough ref sequence for the extension? All coordinates are 0-based
    int ref_len = strlen(splice->vcf.ref);
    int alt_len = strlen(splice->vcf.alt);
    assert( ref_len > alt_len );
    int ndel = ref_len - alt_len;

    if ( tr->strand==STRAND_REV )
    {
        int32_t vcf_ref_end = splice->vcf.pos + ref_len - 1;  // end pos of the VCF REF allele
        int32_t tr_ref_end  = splice->tr->end + N_REF_PAD;    // the end pos of accessible cached ref seq
        if ( vcf_ref_end + ndel > tr_ref_end )
        {
            if ( !small_ref_padding_warned )
            {
                fprintf(stderr,"Warning: Could not verify synonymous start/stop at %s:%d due to small N_REF_PAD. (Improve me?)\n",bcf_seqname(args->hdr,splice->vcf.rec),splice->vcf.pos+1);
                small_ref_padding_warned = 1;
            }
            return 0;
        }

        char *ptr_vcf = splice->vcf.ref + alt_len;                         // the first deleted base in the VCF REF allele
        char *ptr_ref = TSCRIPT_AUX(splice->tr)->ref + N_REF_PAD + (vcf_ref_end + 1 - splice->tr->beg);  // the first ref base after the ndel bases deleted
#if XDBG
        fprintf(stderr,"vcf: %s\nref: %s\n",ptr_vcf,ptr_ref);
#endif
        int i = 0;
        while ( ptr_vcf[i] && ptr_vcf[i]==ptr_ref[i] ) i++;
        if ( ptr_vcf[i] ) return 0;       // the deleted sequence cannot be replaced
    }
    else if ( tr->strand==STRAND_FWD )
    {
        // STRAND_FWD
        int32_t vcf_block_beg = splice->vcf.pos + ref_len - 2*ndel;        // the position of the first base of the ref block that could potentially replace the deletion
        if ( vcf_block_beg < 0 ) return 0;

#if XDBG
        fprintf(stderr,"vcf_block_beg: %d\n",vcf_block_beg+1);
#endif

        if ( N_REF_PAD + vcf_block_beg < ex_beg )
        {
            if ( !small_ref_padding_warned )
            {
                fprintf(stderr,"Warning: Could not verify synonymous start/stop at %s:%d due to small N_REF_PAD. (Improve me?)\n",bcf_seqname(args->hdr,splice->vcf.rec),splice->vcf.pos+1);
                small_ref_padding_warned = 1;
            }
            return 0;
        }

        char *ptr_vcf = splice->vcf.ref + alt_len;                                      // the first deleted base in the VCF REF allele
        char *ptr_ref = TSCRIPT_AUX(splice->tr)->ref + N_REF_PAD + vcf_block_beg - splice->tr->beg;  // the replacement ref block
#if XDBG
        fprintf(stderr,"vcf: %s\nref: %s\n",ptr_vcf,ptr_ref);
#endif

        int i = 0;
        while ( ptr_vcf[i] && ptr_vcf[i]==ptr_ref[i] ) i++;
        if ( ptr_vcf[i] ) return 0;       // the deleted sequence cannot be replaced
    }

    return 1;
}

static inline int splice_csq_del(args_t *args, splice_t *splice, uint32_t ex_beg, uint32_t ex_end)
{
    if ( splice->check_start )
    {
        // check for synonymous start
        //      test/csq/ENST00000375992/incorrect-synon-del-not-start-lost.txt
        //      test/csq/ENST00000368801.2/start-lost.txt
        //      test/csq/ENST00000318249.2/synonymous-start-lost.txt
        int is_synonymous = shifted_del_synonymous(args, splice, ex_beg, ex_end);
        if ( is_synonymous )
        {
            splice->csq |= CSQ_START_RETAINED;
            return SPLICE_OVERLAP;
        }
    }

    // coordinates that matter for consequences, eg AC>ACG trimmed to C>CG
    splice->ref_beg = splice->vcf.pos + splice->tbeg - 1;                       // 1b before the deleted base
    splice->ref_end = splice->vcf.pos + splice->vcf.rlen - splice->tend - 1;    // the last deleted base

#if XDBG
fprintf(stderr,"splice_csq_del: %s>%s .. ex=%d,%d  beg,end=%d,%d  tbeg,tend=%d,%d  check_utr=%d start,stop,beg,end=%d,%d,%d,%d\n", splice->vcf.ref,splice->vcf.alt,ex_beg,ex_end,splice->ref_beg,splice->ref_end,splice->tbeg,splice->tend,splice->check_utr,splice->check_start,splice->check_stop,splice->check_region_beg,splice->check_region_end);
#endif

    if ( splice->ref_beg + 1 < ex_beg )     // the part before the exon; ref_beg is off by -1
    {
        if ( splice->check_region_beg )
        {
            int csq = 0;
            if ( splice->check_utr )
            {
                regitr_t *itr = regitr_init(NULL);
                const char *chr = drop_chr_prefix(args, bcf_seqname(args->hdr,splice->vcf.rec));
                if ( regidx_overlap(args->idx_utr,chr,splice->ref_beg,ex_beg-1, itr) )     // adjacent utr
                    csq = csq_stage_utr(args, itr, splice->vcf.rec, splice->tr->id, splice->csq, splice->vcf.ial);
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
                const char *chr = drop_chr_prefix(args, bcf_seqname(args->hdr,splice->vcf.rec));
                if ( regidx_overlap(args->idx_utr,chr,ex_end+1,splice->ref_end, itr) )     // adjacent utr
                    csq = csq_stage_utr(args, itr, splice->vcf.rec, splice->tr->id, splice->csq, splice->vcf.ial);
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
        csq_stage_splice(args, splice->vcf.rec, splice->tr, splice->csq, splice->vcf.ial);
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
            int ref_beg = splice->ref_beg + splice->kalt.l - 1;     // 0 for AAA>A, 1 for AAA>AC
            if ( ref_beg < splice->ref_end )
                splice->csq |= (splice->ref_end - ref_beg)%3 ? CSQ_FRAMESHIFT_VARIANT : CSQ_INFRAME_DELETION;
            return SPLICE_OVERLAP;
        }
    }
    csq_stage_splice(args, splice->vcf.rec, splice->tr, splice->csq, splice->vcf.ial);
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
                const char *chr = drop_chr_prefix(args, bcf_seqname(args->hdr,splice->vcf.rec));
                if ( regidx_overlap(args->idx_utr,chr,splice->ref_beg,ex_beg-1, itr) )     // adjacent utr
                    csq = csq_stage_utr(args, itr, splice->vcf.rec, splice->tr->id, splice->csq, splice->vcf.ial);
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
                const char *chr = drop_chr_prefix(args, bcf_seqname(args->hdr,splice->vcf.rec));
                if ( regidx_overlap(args->idx_utr,chr,ex_end+1,splice->ref_end, itr) )     // adjacent utr
                    csq = csq_stage_utr(args, itr, splice->vcf.rec, splice->tr->id, splice->csq, splice->vcf.ial);
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
        csq_stage_splice(args, splice->vcf.rec, splice->tr, splice->csq, splice->vcf.ial);
        return SPLICE_OUTSIDE;
    }

    if ( splice->ref_beg < ex_beg + 3 )
    {
        if ( splice->check_region_beg ) splice->csq |= CSQ_SPLICE_REGION;
        if ( splice->tr->strand==STRAND_FWD ) { if ( splice->check_start ) splice->csq |= CSQ_START_LOST; }
        else if ( splice->tr->strand==STRAND_REV ) { if ( splice->check_stop ) splice->csq |= CSQ_STOP_LOST; }
    }
    if ( splice->ref_end > ex_end - 3 )
    {
        if ( splice->check_region_end ) splice->csq |= CSQ_SPLICE_REGION;
        if ( splice->tr->strand==STRAND_REV ) { if ( splice->check_start ) splice->csq |= CSQ_START_LOST; }
        else if ( splice->tr->strand==STRAND_FWD ) { if ( splice->check_stop ) splice->csq |= CSQ_STOP_LOST; }
    }
    if ( splice->set_refalt )
    {
        splice->vcf.rlen -= splice->tbeg + splice->tend;
        splice->kref.l = 0; kputsn(splice->vcf.ref + splice->tbeg, splice->vcf.rlen, &splice->kref);
        splice->kalt.l = 0; kputsn(splice->vcf.alt + splice->tbeg, splice->vcf.rlen, &splice->kalt);
    }
    csq_stage_splice(args, splice->vcf.rec, splice->tr, splice->csq, splice->vcf.ial);
    return SPLICE_INSIDE;
}
static inline int splice_csq(args_t *args, splice_t *splice, uint32_t ex_beg, uint32_t ex_end)
{
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
    gf_tscript_t *tr = cds->tr;
    child->icds = cds->icds;     // index of cds in the tscript's list of exons
    child->vcf_ial = ial;

    splice_t splice;
    splice_init(&splice, rec);
    splice.tr = tr;
    splice.vcf.ial  = ial;
    splice.vcf.alt  = rec->d.allele[ial];
    splice.check_acceptor = splice.check_donor = splice.set_refalt = splice.check_utr = 1;
    if ( !(tr->trim & TRIM_5PRIME) )
    {
        if ( tr->strand==STRAND_FWD ) { if ( child->icds==0 ) splice.check_start = 1; }
        else if ( tr->strand==STRAND_REV ) { if ( child->icds==tr->ncds-1 ) splice.check_start = 1; }
    }
    if ( !(tr->trim & TRIM_3PRIME) )
    {
        if ( tr->strand==STRAND_FWD ) { if ( child->icds==tr->ncds-1 ) splice.check_stop = 1; }
        else if ( tr->strand==STRAND_REV ) { if ( child->icds==0 ) splice.check_stop = 1; }
    }
    if ( splice.check_start )   // do not check starts in incomplete CDS, defined as not starting with M
    {
        if ( tr->strand==STRAND_FWD ) { if ( dna2aa(TSCRIPT_AUX(tr)->ref+N_REF_PAD+cds->beg-tr->beg) != 'M' ) splice.check_start = 0; }
        else if ( tr->strand==STRAND_REV ) { if ( cdna2aa(TSCRIPT_AUX(tr)->ref+N_REF_PAD+cds->beg-tr->beg+cds->len-3) != 'M' ) splice.check_start = 0; }
    }
    if ( child->icds!=0 ) splice.check_region_beg = 1;
    if ( child->icds!=tr->ncds-1 ) splice.check_region_end = 1;

#if XDBG
fprintf(stderr,"\nhap_init: %d [%s][%s]   check start:%d,stop:%d\n",splice.vcf.pos+1,splice.vcf.ref,splice.vcf.alt,splice.check_start,splice.check_stop);
#endif
    int ret = splice_csq(args, &splice, cds->beg, cds->beg + cds->len - 1);
#if XDBG
fprintf(stderr,"cds splice_csq: %d [%s][%s] .. beg,end=%d %d, ret=%d, csq=%d\n\n",splice.vcf.pos+1,splice.kref.s,splice.kalt.s,splice.ref_beg+1,splice.ref_end+1,ret,splice.csq);
#endif

    if ( ret==SPLICE_VAR_REF ) return 2;  // not a variant, eg REF=CA ALT=CA
    if ( ret==SPLICE_OUTSIDE || ret==SPLICE_OVERLAP || splice.csq==CSQ_START_LOST )  // not a coding csq
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

    assert( parent->type!=HAP_SSS );
    if ( parent->type==HAP_CDS )
    {
        i = parent->icds;
        if ( i!=cds->icds )
        {
            // the variant is on a new exon, finish up the previous
            int len = tr->cds[i]->len - parent->rbeg - parent->rlen + tr->cds[i]->beg;
            if ( len > 0 )
                kputsn_(TSCRIPT_AUX(tr)->ref + N_REF_PAD + parent->rbeg + parent->rlen - tr->beg, len, &str);
        }

        // append any skipped non-variant exons
        while ( ++i < cds->icds )
            kputsn_(TSCRIPT_AUX(tr)->ref + N_REF_PAD + tr->cds[i]->beg - tr->beg, tr->cds[i]->len, &str);

        if ( parent->icds==child->icds )
        {
            int len = splice.ref_beg - parent->rbeg - parent->rlen;
            if ( len < 0 )   // overlapping variants
            {
                free(str.s);
                free(splice.kref.s);
                free(splice.kalt.s);
                return 1;
            }
            kputsn_(TSCRIPT_AUX(tr)->ref + N_REF_PAD + parent->rbeg + parent->rlen - tr->beg, len, &str);
        }
        else
            kputsn_(TSCRIPT_AUX(tr)->ref + N_REF_PAD + cds->beg - tr->beg, splice.ref_beg - cds->beg, &str);
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
fprintf(stderr,"\ntranslate: %d %d %d  fill=%d  seq.l=%d\n",sbeg,rbeg,rend,fill,(int)_seq->l);
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
    else if ( strand==STRAND_REV )
    {
        // right padding - number of bases to take from ref
        npad = (seq.m - (sbeg + seq.l)) % 3;
#if DBG>1
        fprintf(stderr,"    npad: %d\n",npad);
#endif
        if ( !(npad>=0 && sbeg+seq.l+npad<=seq.m) ) fprintf(stderr,"sbeg=%d  seq.l=%d seq.m=%d npad=%d\n",sbeg,(int)seq.l,(int)seq.m,npad);
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
        if(i==1)fprintf(stderr,"[0]  %c\n",tmp[2]);
        if(i==0)fprintf(stderr,"[0] %c%c\n",tmp[1],tmp[2]);
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
    else error("Should not happen: %d\n", strand);
    kputc_(0,tseq); tseq->l--;
#if DBG
 fprintf(stderr,"    tseq: %s\n", tseq->s);
#endif
}

void tscript_splice_ref(gf_tscript_t *tr)
{
    int i, len = 0;
    for (i=0; i<tr->ncds; i++)
        len += tr->cds[i]->len;

    TSCRIPT_AUX(tr)->nsref = len + 2*N_REF_PAD;
    TSCRIPT_AUX(tr)->sref  = (char*) malloc(len + 1 + 2*N_REF_PAD);
    len = 0;

    memcpy(TSCRIPT_AUX(tr)->sref, TSCRIPT_AUX(tr)->ref + tr->cds[0]->beg - tr->beg, N_REF_PAD);
    len += N_REF_PAD;

    for (i=0; i<tr->ncds; i++)
    {
        memcpy(TSCRIPT_AUX(tr)->sref + len, TSCRIPT_AUX(tr)->ref + N_REF_PAD + tr->cds[i]->beg - tr->beg, tr->cds[i]->len);
        len += tr->cds[i]->len;
    }
    memcpy(TSCRIPT_AUX(tr)->sref + len, TSCRIPT_AUX(tr)->ref + N_REF_PAD + tr->cds[tr->ncds-1]->beg - tr->beg, N_REF_PAD);
    len += N_REF_PAD;

    TSCRIPT_AUX(tr)->sref[len] = 0;
}

// returns: 0 if consequence was added, 1 if it already exists or could not be added
int csq_push(args_t *args, csq_t *csq, bcf1_t *rec)
{
#if XDBG
fprintf(stderr,"csq_push: %d .. %d\n",(int)rec->pos+1,csq->type.type);
#endif
    khint_t k = kh_get(pos2vbuf, args->pos2vbuf, (int)csq->pos);
    vbuf_t *vbuf = (k == kh_end(args->pos2vbuf)) ? NULL : kh_val(args->pos2vbuf, k);
    if ( !vbuf ) error("This should not happen. %s:%d  %s\n",bcf_seqname(args->hdr,rec),csq->pos+1,csq->type.vstr.s);

    int i;
    for (i=0; i<vbuf->n; i++)
        if ( vbuf->vrec[i]->line==rec ) break;
    if ( i==vbuf->n ) error("This should not happen.. %s:%d  %s\n", bcf_seqname(args->hdr,rec),csq->pos+1,csq->type.vstr.s);
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
            if ( csq->type.vcf_ial != vrec->vcsq[i].vcf_ial ) continue;
            if ( (csq->type.type&CSQ_UPSTREAM_STOP)^(vrec->vcsq[i].type&CSQ_UPSTREAM_STOP) ) continue;  // both must or mustn't have upstream_stop
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
    csq->idx  = vrec->nvcsq++;
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

    int has_csq = 0, i, n = sizeof(csq_strings)/sizeof(char*);
    for (i=1; i<n; i++)
        if ( csq_strings[i] && csq->type&(1<<i) ) { has_csq = 1; kputs(csq_strings[i],str); break; }
    i++;
    for (; i<n; i++)
        if ( csq_strings[i] && csq->type&(1<<i) ) { has_csq = 1; kputc_('&',str); kputs(csq_strings[i],str); }

    if ( (csq->biotype==GF_NMD) && (csq->type & CSQ_PRN_NMD) )
    {
        if ( has_csq ) kputc_('&',str); // just in case, this should always be true
        kputs("NMD_transcript",str);
    }

    kputc_('|', str);
    if ( csq->gene ) kputs(csq->gene , str);

    kputc_('|', str);
//    if ( csq->type & CSQ_PRN_TSCRIPT ) kputs(args->tscript_ids.str[csq->trid], str);
    if ( csq->type & CSQ_PRN_TSCRIPT ) kputs(gff_id2string(args->gff,transcript,csq->trid), str);

    kputc_('|', str);
    kputs(gf_type2gff_string(csq->biotype), str);

    if ( CSQ_PRN_STRAND(csq->type) || csq->vstr.l )
        kputs(csq->strand==STRAND_FWD ? "|+" : (csq->strand==STRAND_REV ? "|-" : "|."), str);

    if ( csq->vstr.l )
        kputs(csq->vstr.s, str);
}

void kprint_aa_prediction(args_t *args, int beg, kstring_t *aa, kstring_t *str)
{
    if ( !args->brief_predictions || (int)aa->l - args->brief_predictions < 3 )
        kputs(aa->s, str);
    else
    {
        int i, len = aa->l;
        if ( aa->s[len-1]=='*' ) len--;
        for (i=0; i<len && i<args->brief_predictions; i++) kputc(aa->s[i], str);
        kputs("..", str);
        kputw(beg+len, str);
    }
}

void hap_add_csq(args_t *args, hap_t *hap, hap_node_t *node, int tlen, int ibeg, int iend, int dlen, int indel)
{
    int i;
    gf_tscript_t *tr = hap->tr;
    assert( tr->strand==STRAND_FWD || tr->strand==STRAND_REV );
    int ref_node = tr->strand==STRAND_FWD ? ibeg : iend;
    int icsq = node->ncsq_list++;
    hts_expand0(csq_t,node->ncsq_list,node->mcsq_list,node->csq_list);
    csq_t *csq = &node->csq_list[icsq];
    csq->pos  = hap->stack[ref_node].node->rec->pos;
    csq->type.trid    = tr->id;
    csq->type.vcf_ial = node->vcf_ial;
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
            if ( hap->tref.s[hap->tref.l-1]!='*' && hap->tseq.s[hap->tseq.l-1]=='*' )
                csq->type.type |= CSQ_STOP_GAINED;
        }
        else
        {
            int aa_change = 0;
            for (i=0; i<hap->tref.l; i++)
            {
                if ( hap->tref.s[i] == hap->tseq.s[i] ) continue;
                aa_change = 1;
                if ( hap->tref.s[i] ==  '*' )
                    csq->type.type |= CSQ_STOP_LOST;
                else if ( hap->tseq.s[i] ==  '*' )
                    csq->type.type |= CSQ_STOP_GAINED;
                else
                    csq->type.type |= CSQ_MISSENSE_VARIANT;
            }
            if ( !aa_change )
                csq->type.type |= CSQ_SYNONYMOUS_VARIANT;
        }
    }
    // Check if compound inframe variants are real inframes, or if the stop codon occurs before the frameshift can be restored
    if ( ibeg!=iend && (csq->type.type & (CSQ_INFRAME_DELETION|CSQ_INFRAME_INSERTION|CSQ_INFRAME_ALTERING)) && hap->tseq.s[hap->tseq.l-1]=='*' )
    {
        rm_csq |= CSQ_INFRAME_DELETION | CSQ_INFRAME_INSERTION | CSQ_INFRAME_ALTERING;
        csq->type.type |= CSQ_FRAMESHIFT_VARIANT | CSQ_STOP_GAINED;
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
    int aa_rbeg = tr->strand==STRAND_FWD ? node2rbeg(ibeg)/3+1 : (TSCRIPT_AUX(hap->tr)->nsref - 2*N_REF_PAD - node2rend(iend))/3+1;
    int aa_sbeg = tr->strand==STRAND_FWD ? node2sbeg(ibeg)/3+1 : (tlen - node2send(iend))/3+1;
    kputc_('|', &str);
    kputw(aa_rbeg, &str);
    kprint_aa_prediction(args,aa_rbeg,&hap->tref,&str);
    if ( !(csq->type.type & CSQ_SYNONYMOUS_VARIANT) )
    {
        kputc_('>', &str);
        kputw(aa_sbeg, &str);
        kprint_aa_prediction(args,aa_sbeg,&hap->tseq,&str);
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
            //??tmp_csq->type.vcf_ial = node->vcf_ial;  .. this should not be needed for non-compound variants
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
            //??tmp_csq->type.vcf_ial = node->vcf_ial;  .. this should not be needed for non-compound variants
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
    gf_tscript_t *tr = hap->tr;
    if ( !TSCRIPT_AUX(tr)->sref )
        tscript_splice_ref(tr);

    kstring_t sref;
    sref.s = TSCRIPT_AUX(tr)->sref;
    sref.l = TSCRIPT_AUX(tr)->nsref;
    sref.m = sref.l;

    int istack = 0;
    hts_expand(hstack_t,1,hap->mstack,hap->stack);

    hap->sseq.l = 0;
    hap->tseq.l = 0;
    hap->stack[0].node = TSCRIPT_AUX(tr)->root;
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
        hap->sbeg = hap->stack[i].node->sbeg;
        assert( hap->stack[istack].node->type != HAP_SSS );
        if ( tr->strand==STRAND_FWD )
        {
            i = 0, ibeg = -1;
            while ( ++i <= istack )
            {
                assert( hap->stack[i].node->type != HAP_SSS );

                dlen += hap->stack[i].node->dlen;
                if ( hap->stack[i].node->dlen ) indel = 1;

                // This condition extends compound variants.
                if ( i<istack )
                {
                    if ( dlen%3 )   // frameshift
                    {
                        if ( ibeg==-1 ) ibeg = i;
                        continue;
                    }
                    // the last base of the current variant vs the first base of the next
                    // variant: are they in the same codon? (forward strand)
                    int icur  = node2sbeg(i);
                    int inext = node2sbeg(i+1);
                    if ( hap->stack[i].node->dlen > 0 ) icur += hap->stack[i].node->dlen;
                    else if ( hap->stack[i].node->dlen < 0 ) icur++;
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
        else if ( tr->strand==STRAND_REV )
        {
            i = istack + 1, ibeg = -1;
            while ( --i > 0 )
            {
                assert ( hap->stack[i].node->type != HAP_SSS );
                dlen += hap->stack[i].node->dlen;
                if ( hap->stack[i].node->dlen ) indel = 1;
                if ( i>1 )
                {
                    if ( dlen%3 )
                    {
                        if ( ibeg==-1 ) ibeg = i;
                        continue;
                    }
                    // the last base of the current variant vs the first base of the next
                    // variant: are they in the same codon? (reverse strand)
                    int icur  = sseq.m - 1 - node2sbeg(i);
                    int inext = sseq.m - 1 - node2sbeg(i-1);
                    if ( hap->stack[i].node->dlen > 0 ) icur += hap->stack[i].node->dlen - 1;
                    else if ( hap->stack[i].node->dlen < 0 ) icur -= hap->stack[i].node->dlen;
                    if ( hap->stack[i-1].node->dlen > 0 ) inext -= hap->stack[i-1].node->dlen;
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
static inline void hap_print_text(args_t *args, gf_tscript_t *tr, int ismpl, int ihap, hap_node_t *node)
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

static inline void hap_stage_vcf(args_t *args, gf_tscript_t *tr, int ismpl, int ihap, hap_node_t *node)
{
    if ( !node || !node->ncsq_list || ismpl<0 ) return;

    int i;
    for (i=0; i<node->ncsq_list; i++)
    {
        csq_t *csq = node->csq_list + i;
        vrec_t *vrec = csq->vrec;
        int icsq2 = 2*csq->idx + ihap;
        if ( icsq2 >= args->ncsq2_max ) // more than ncsq2_max consequences, so can't fit it in FMT
        {
            if ( args->verbosity && (!args->ncsq2_small_warned || args->verbosity > 1) )
            {
                fprintf(stderr,
                    "Warning: Too many consequences for sample %s at %s:%"PRId64", keeping the first %d and skipping the rest.\n",
                    args->hdr->samples[ismpl],bcf_hdr_id2name(args->hdr,args->rid),(int64_t) vrec->line->pos+1,csq->idx);
                if ( !args->ncsq2_small_warned )
                    fprintf(stderr,"         The limit can be increased by setting the --ncsq parameter. This warning is printed only once.\n");
            }
            if ( args->ncsq2_small_warned < icsq2 ) args->ncsq2_small_warned = icsq2;
            break;
        }
        int ival, ibit;
        icsq2_to_bit(icsq2, &ival,&ibit);
        if ( vrec->nfmt < 1 + ival ) vrec->nfmt = 1 + ival;
        vrec->fmt_bm[ismpl*args->nfmt_bcsq + ival] |= 1 << ibit;
    }
}

void hap_flush(args_t *args, uint32_t pos)
{
    int i,j;
    tr_heap_t *heap = args->active_tr;
    while ( heap->ndat && heap->dat[0]->end<=pos )
    {
        gf_tscript_t *tr = heap->dat[0];
        khp_delete(trhp, heap);
        args->hap->tr = tr;
        if ( TSCRIPT_AUX(tr)->root && TSCRIPT_AUX(tr)->root->nchild ) // normal, non-localized calling
        {
            hap_finalize(args, args->hap);

            if ( args->output_type==FT_TAB_TEXT )   // plain text output, not a vcf
            {
                if ( args->phase==PHASE_DROP_GT )
                    hap_print_text(args, tr, -1,0, TSCRIPT_AUX(tr)->hap[0]);
                else
                {
                    for (i=0; i<args->smpl->n; i++)
                    {
                        for (j=0; j<2; j++)
                            hap_print_text(args, tr, args->smpl->idx[i],j+1, TSCRIPT_AUX(tr)->hap[i*2+j]);
                    }
                }
            }
            else if ( args->phase!=PHASE_DROP_GT )
            {
                for (i=0; i<args->smpl->n; i++)
                {
                    for (j=0; j<2; j++)
                        hap_stage_vcf(args, tr, args->smpl->idx[i],j, TSCRIPT_AUX(tr)->hap[i*2+j]);
                }
            }
        }

        // mark the transcript for deletion. Cannot delete it immediately because
        // by-position VCF output will need them when flushed by vcf_buf_push
        args->nrm_tr++;
        hts_expand(gf_tscript_t*,args->nrm_tr,args->mrm_tr,args->rm_tr);
        args->rm_tr[args->nrm_tr-1] = tr;
    }
}

#define SWAP(type_t, a, b) { type_t t = a; a = b; b = t; }

vbuf_t *vbuf_push(args_t *args, bcf1_t **rec_ptr)
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
        args->vcf_buf[i]->keep_until = 0;
    }
    vbuf_t *vbuf = args->vcf_buf[i];
    vbuf->n++;
    hts_expand0(vrec_t*, vbuf->n, vbuf->m, vbuf->vrec);
    if ( !vbuf->vrec[vbuf->n - 1] )
        vbuf->vrec[vbuf->n - 1] = (vrec_t*) calloc(1,sizeof(vrec_t));

    vrec_t *vrec = vbuf->vrec[vbuf->n - 1];
    if ( args->phase!=PHASE_DROP_GT && args->smpl->n )
    {
        if ( !vrec->fmt_bm ) vrec->fmt_bm = (uint32_t*) calloc(args->hdr_nsmpl,sizeof(*vrec->fmt_bm) * args->nfmt_bcsq);
        else memset(vrec->fmt_bm,0,args->hdr_nsmpl*sizeof(*vrec->fmt_bm) * args->nfmt_bcsq);
    }
    if ( !vrec->line ) vrec->line = bcf_init1();
    SWAP(bcf1_t*, (*rec_ptr), vrec->line);

    int ret;
    khint_t k = kh_put(pos2vbuf, args->pos2vbuf, (int)rec->pos, &ret);
    kh_val(args->pos2vbuf,k) = vbuf;

    return vbuf;
}

void vbuf_flush(args_t *args, uint32_t pos)
{
    int i,j;
    while ( args->vcf_rbuf.n )
    {
        vbuf_t *vbuf;
        if ( !args->local_csq && args->active_tr->ndat )
        {
            // check if the first active transcript starts beyond the first buffered VCF record,
            // cannot output buffered VCF lines (args.vbuf) until the active transcripts are gone
            vbuf = args->vcf_buf[ args->vcf_rbuf.f ];
            if ( vbuf->keep_until > pos ) break;
            assert( vbuf->n );
        }

        i = rbuf_shift(&args->vcf_rbuf);
        assert( i>=0 );
        vbuf = args->vcf_buf[i];
        int pos = vbuf->n ? vbuf->vrec[0]->line->pos : -1;
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
                if ( bcf_write(args->out_fh, args->hdr, vrec->line)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname?args->output_fname:"standard output");
                int save_pos = vrec->line->pos;
                bcf_empty(vrec->line);
                vrec->line->pos = save_pos;  // this is necessary for compound variants
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
                    for (j=1; j<args->hdr_nsmpl; j++)
                        memmove(&vrec->fmt_bm[j*vrec->nfmt], &vrec->fmt_bm[j*args->nfmt_bcsq], vrec->nfmt*sizeof(*vrec->fmt_bm));
                bcf_update_format_int32(args->hdr, vrec->line, args->bcsq_tag, vrec->fmt_bm, args->hdr_nsmpl*vrec->nfmt);
            }
            vrec->nvcsq = 0;
            if ( bcf_write(args->out_fh, args->hdr, vrec->line)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname?args->output_fname:"standard output");
            int save_pos = vrec->line->pos;
            bcf_empty(vrec->line);
            vrec->line->pos = save_pos;
        }
        if ( pos!=-1 )
        {
            khint_t k = kh_get(pos2vbuf, args->pos2vbuf, pos);
            if ( k != kh_end(args->pos2vbuf) ) kh_del(pos2vbuf, args->pos2vbuf, k);
        }
        vbuf->n = 0;
    }
    if ( args->active_tr->ndat ) return;

    for (i=0; i<args->nrm_tr; i++)
    {
        gf_tscript_t *tr = args->rm_tr[i];
        tscript_t *aux = TSCRIPT_AUX(tr);
        if ( aux->root ) hap_destroy(aux->root);
        aux->root = NULL;
        free(aux->hap);
        free(aux->ref);
        free(aux->sref);
        free(aux);
        tr->aux = NULL;
    }
    args->nrm_tr = 0;
    args->ncsq_buf = 0;
}

void tscript_init_ref(args_t *args, gf_tscript_t *tr, const char *chr)
{
    int i, len;
    int pad_beg = tr->beg >= N_REF_PAD ? N_REF_PAD : tr->beg;

    const char *tmp_chr = chr;
    if ( !faidx_has_seq(args->fai,tmp_chr) )
    {
        tmp_chr = drop_chr_prefix(args,chr);
        if ( !faidx_has_seq(args->fai,tmp_chr) ) tmp_chr = add_chr_prefix(args,chr);
    }
    TSCRIPT_AUX(tr)->ref = faidx_fetch_seq(args->fai, tmp_chr, tr->beg - pad_beg, tr->end + N_REF_PAD, &len);
    if ( !TSCRIPT_AUX(tr)->ref )
        error("faidx_fetch_seq failed %s:%d-%d\n", chr,tr->beg+1,tr->end+1);

    int pad_end = len - (tr->end - tr->beg + 1 + pad_beg);
    if ( pad_beg + pad_end != 2*N_REF_PAD )
    {
        char *ref = (char*) malloc(tr->end - tr->beg + 1 + 2*N_REF_PAD + 1);
        for (i=0; i < N_REF_PAD - pad_beg; i++) ref[i] = 'N';
        memcpy(ref+i, TSCRIPT_AUX(tr)->ref, len);
        len += i;
        for (i=0; i < N_REF_PAD - pad_end; i++) ref[i+len] = 'N';
        ref[i+len] = 0;
        free(TSCRIPT_AUX(tr)->ref);
        TSCRIPT_AUX(tr)->ref = ref;
    }
}

static void sanity_check_ref(args_t *args, gf_tscript_t *tr, bcf1_t *rec)
{
    int vbeg = 0;
    int rbeg = rec->pos - tr->beg + N_REF_PAD;
    if ( rbeg < 0 ) { vbeg += abs(rbeg); rbeg = 0; }
    char *ref = TSCRIPT_AUX(tr)->ref + rbeg;
    char *vcf = rec->d.allele[0] + vbeg;
    assert( vcf - rec->d.allele[0] < strlen(rec->d.allele[0]) && ref - TSCRIPT_AUX(tr)->ref < tr->end - tr->beg + 2*N_REF_PAD );
    int i = 0;
    while ( ref[i] && vcf[i] )
    {
        if ( ref[i]!=vcf[i] && toupper(ref[i])!=toupper(vcf[i]) )
            error("Error: the fasta reference does not match the VCF REF allele at %s:%"PRId64" .. fasta=%c vcf=%c\n",
                    bcf_seqname(args->hdr,rec),(int64_t) rec->pos+vbeg+1,ref[i],vcf[i]);
        i++;
    }
}

int test_cds_local(args_t *args, bcf1_t *rec)
{
    int i,j, ret = 0;
    const char *chr = drop_chr_prefix(args, bcf_seqname(args->hdr,rec));
    // note that the off-by-one extension of rlen is deliberate to account for insertions
    if ( !regidx_overlap(args->idx_cds,chr,rec->pos,rec->pos+rec->rlen, args->itr) ) return 0;

    // structures to fake the normal test_cds machinery
    hap_node_t root, node;
    root.type  = HAP_ROOT;
    kstring_t *tref = &args->hap->tref, *tseq = &args->hap->tseq;

    while ( regitr_overlap(args->itr) )
    {
        gf_cds_t *cds = regitr_payload(args->itr,gf_cds_t*);
        gf_tscript_t *tr = cds->tr;
        if ( !GF_is_coding(tr->type) ) continue;
        ret = 1;

        if ( !TSCRIPT_AUX(tr) )
        {
            tr->aux = calloc(sizeof(tscript_t),1);
            tscript_init_ref(args, tr, chr);
            tscript_splice_ref(tr);
            khp_insert(trhp, args->active_tr, &tr);     // only to clean the reference afterwards
        }

        sanity_check_ref(args, tr, rec);

        kstring_t sref;
        sref.s = TSCRIPT_AUX(tr)->sref;
        sref.l = TSCRIPT_AUX(tr)->nsref;
        sref.m = sref.l;

        for (i=1; i<rec->n_allele; i++)
        {
            if ( rec->d.allele[i][0]=='<' || rec->d.allele[i][0]=='*' ) { continue; }
            if ( hap_init(args, &root, &node, cds, rec, i)!=0 ) continue;

            csq_t csq;
            memset(&csq, 0, sizeof(csq_t));
            csq.pos          = rec->pos;
            csq.type.biotype = tr->type;
            csq.type.strand  = tr->strand;
            csq.type.trid    = tr->id;
            csq.type.vcf_ial = i;
            csq.type.gene    = tr->gene->name;

            int csq_type = node.csq;

            // code repetition: it would be nice to reuse the code from hap_add_csq, needs refactoring though
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
                    if ( tref->s[tref->l-1]!='*' && tseq->s[tseq->l-1]=='*' )
                        csq_type |= CSQ_STOP_GAINED;
                }
                else
                {
                    int aa_change = 0;
                    for (j=0; j<tref->l; j++)
                    {
                        if ( tref->s[j] == tseq->s[j] ) continue;
                        aa_change = 1;
                        if ( tref->s[j] ==  '*' )
                            csq_type |= CSQ_STOP_LOST;
                        else if ( tseq->s[j] ==  '*' )
                            csq_type |= CSQ_STOP_GAINED;
                        else
                            csq_type |= CSQ_MISSENSE_VARIANT;
                    }
                    if ( !aa_change )
                        csq_type |= CSQ_SYNONYMOUS_VARIANT;
                }
                if ( csq_type & CSQ_COMPOUND )
                {
                    // create the aa variant string
                    kstring_t str = {0,0,0};
                    int aa_rbeg = tr->strand==STRAND_FWD ? node.sbeg/3+1 : (TSCRIPT_AUX(tr)->nsref - 2*N_REF_PAD - node.sbeg - node.rlen)/3+1;
                    int aa_sbeg = tr->strand==STRAND_FWD ? node.sbeg/3+1 : (TSCRIPT_AUX(tr)->nsref - 2*N_REF_PAD + node.dlen - node.sbeg - alen)/3+1;
                    kputc_('|', &str);
                    kputw(aa_rbeg, &str);
                    kprint_aa_prediction(args,aa_rbeg,tref,&str);
                    if ( !(csq_type & CSQ_SYNONYMOUS_VARIANT) )
                    {
                        kputc_('>', &str);
                        kputw(aa_sbeg, &str);
                        kprint_aa_prediction(args,aa_sbeg,tseq,&str);
                    }
                    kputc_('|', &str);
                    kputw(rec->pos+1, &str);
                    kputs(node.var, &str);
                    csq.type.vstr = str;
                    csq.type.type = csq_type & CSQ_COMPOUND;
                    csq_stage(args, &csq, rec);

                    // all this only to clean vstr when vrec is flushed
                    if ( !TSCRIPT_AUX(tr)->root )
                        TSCRIPT_AUX(tr)->root = (hap_node_t*) calloc(1,sizeof(hap_node_t));
                    TSCRIPT_AUX(tr)->root->ncsq_list++;
                    hts_expand0(csq_t,TSCRIPT_AUX(tr)->root->ncsq_list,TSCRIPT_AUX(tr)->root->mcsq_list,TSCRIPT_AUX(tr)->root->csq_list);
                    csq_t *rm_csq = TSCRIPT_AUX(tr)->root->csq_list + TSCRIPT_AUX(tr)->root->ncsq_list - 1;
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

int test_cds(args_t *args, bcf1_t *rec, vbuf_t *vbuf)
{
    static int overlaps_warned = 0, multiploid_warned = 0;

    int i, ret = 0, hap_ret;
    const char *chr = drop_chr_prefix(args, bcf_seqname(args->hdr,rec));
    // note that the off-by-one extension of rlen is deliberate to account for insertions
    if ( !regidx_overlap(args->idx_cds,chr,rec->pos,rec->pos+rec->rlen, args->itr) ) return 0;
    while ( regitr_overlap(args->itr) )
    {
        gf_cds_t *cds = regitr_payload(args->itr,gf_cds_t*);
        gf_tscript_t *tr = cds->tr;
        if ( !GF_is_coding(tr->type) ) continue;
        if ( vbuf->keep_until < tr->end ) vbuf->keep_until = tr->end;
        ret = 1;
        if ( !TSCRIPT_AUX(tr) )
        {
            // initialize the transcript and its haplotype tree, fetch the reference sequence
            tr->aux = calloc(sizeof(tscript_t),1);
            tscript_init_ref(args, tr, chr);

            TSCRIPT_AUX(tr)->root = (hap_node_t*) calloc(1,sizeof(hap_node_t));
            TSCRIPT_AUX(tr)->nhap = args->phase==PHASE_DROP_GT ? 1 : 2*args->smpl->n;     // maximum ploidy = diploid
            TSCRIPT_AUX(tr)->hap  = (hap_node_t**) malloc(TSCRIPT_AUX(tr)->nhap*sizeof(hap_node_t*));
            for (i=0; i<TSCRIPT_AUX(tr)->nhap; i++) TSCRIPT_AUX(tr)->hap[i] = NULL;
            TSCRIPT_AUX(tr)->root->nend = TSCRIPT_AUX(tr)->nhap;
            TSCRIPT_AUX(tr)->root->type = HAP_ROOT;

            khp_insert(trhp, args->active_tr, &tr);
        }

        sanity_check_ref(args, tr, rec);

        if ( args->phase==PHASE_DROP_GT )
        {
            if ( rec->d.allele[1][0]=='<' || rec->d.allele[1][0]=='*' ) { continue; }
            hap_node_t *parent = TSCRIPT_AUX(tr)->hap[0] ? TSCRIPT_AUX(tr)->hap[0] : TSCRIPT_AUX(tr)->root;
            hap_node_t *child  = (hap_node_t*)calloc(1,sizeof(hap_node_t));
            hap_ret = hap_init(args, parent, child, cds, rec, 1);
            if ( hap_ret!=0 )
            {
                // overlapping or intron variant, cannot apply
                if ( hap_ret==1 )
                {
                    if ( args->verbosity && (!overlaps_warned || args->verbosity > 1) )
                    {
                        fprintf(stderr,
                            "Warning: Skipping overlapping variants at %s:%"PRId64"\t%s>%s.\n",
                            chr,(int64_t) rec->pos+1,rec->d.allele[0],rec->d.allele[1]);
                        if ( !overlaps_warned )
                            fprintf(stderr,"         This message is printed only once, the verbosity can be increased with `--verbose 2`\n");
                        overlaps_warned = 1;
                    }
                    if ( args->out )
                        fprintf(args->out,"LOG\tWarning: Skipping overlapping variants at %s:%"PRId64"\t%s>%s\n", chr,(int64_t) rec->pos+1,rec->d.allele[0],rec->d.allele[1]);
                }
                else ret = 1;   // prevent reporting as intron in test_tscript
                hap_destroy(child);
                continue;
            }
            if ( child->type==HAP_SSS )
            {
                csq_t csq;
                memset(&csq, 0, sizeof(csq_t));
                csq.pos          = rec->pos;
                csq.type.biotype = tr->type;
                csq.type.strand  = tr->strand;
                csq.type.trid    = tr->id;
                csq.type.vcf_ial = 1;
                csq.type.gene    = tr->gene->name;
                csq.type.type = child->csq;
                csq_stage(args, &csq, rec);
                hap_destroy(child);
                ret = 1;
                continue;
            }
            parent->nend--;
            parent->nchild = 1;
            parent->mchild = 1;
            parent->child  = (hap_node_t**) malloc(sizeof(hap_node_t*));
            parent->child[0] = child;
            TSCRIPT_AUX(tr)->hap[0] = child;
            TSCRIPT_AUX(tr)->hap[0]->nend = 1;
            continue;
        }

        // apply the VCF variants and extend the haplotype tree
        int j, ismpl, ihap, ngts = bcf_get_genotypes(args->hdr, rec, &args->gt_arr, &args->mgt_arr);
        ngts /= bcf_hdr_nsamples(args->hdr);
        if ( ngts!=1 && ngts!=2 )
        {
            if ( args->verbosity && (!multiploid_warned || args->verbosity > 1) )
            {
                fprintf(stderr,
                    "Warning: Skipping site with non-diploid/non-haploid genotypes at %s:%"PRId64"\t%s>%s.\n",
                    chr,(int64_t) rec->pos+1,rec->d.allele[0],rec->d.allele[1]);
                if ( !multiploid_warned )
                    fprintf(stderr,"         This message is printed only once, the verbosity can be increased with `--verbose 2`\n");
                multiploid_warned = 1;
            }
            if ( args->out )
                fprintf(args->out,"LOG\tWarning: Skipping site with non-diploid/non-haploid genotypes at %s:%"PRId64"\t%s>%s\n", chr,(int64_t) rec->pos+1,rec->d.allele[0],rec->d.allele[1]);
            continue;
        }
        for (ismpl=0; ismpl<args->smpl->n; ismpl++)
        {
            int32_t *gt = args->gt_arr + args->smpl->idx[ismpl]*ngts;
            if ( gt[0]==bcf_gt_missing ) continue;

            if ( ngts>1 && gt[1]!=bcf_gt_missing && gt[1]!=bcf_int32_vector_end && bcf_gt_allele(gt[0])!=bcf_gt_allele(gt[1]) )
            {
                if ( args->phase==PHASE_MERGE )
                {
                    if ( !bcf_gt_allele(gt[0]) ) gt[0] = gt[1];
                }
                if ( !bcf_gt_is_phased(gt[0]) && !bcf_gt_is_phased(gt[1]) )
                {
                    if ( args->phase==PHASE_REQUIRE )
                        error("Unphased heterozygous genotype at %s:%"PRId64", sample %s. See the --phase option.\n", chr,(int64_t) rec->pos+1,args->hdr->samples[args->smpl->idx[ismpl]]);
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

                hap_node_t *parent = TSCRIPT_AUX(tr)->hap[i] ? TSCRIPT_AUX(tr)->hap[i] : TSCRIPT_AUX(tr)->root;
                if ( parent->cur_rec==rec && parent->cur_child[ial]>=0 )
                {
                    // this haplotype has been seen in another sample
                    TSCRIPT_AUX(tr)->hap[i] = parent->child[ parent->cur_child[ial] ];
                    TSCRIPT_AUX(tr)->hap[i]->nend++;
                    parent->nend--;
                    continue;
                }

                hap_node_t *child = (hap_node_t*)calloc(1,sizeof(hap_node_t));
                hap_ret = hap_init(args, parent, child, cds, rec, ial);
                if ( hap_ret!=0 )
                {
                    // overlapping or intron variant, cannot apply
                    if ( hap_ret==1 )
                    {
                        if ( args->verbosity && (!overlaps_warned || args->verbosity > 1) )
                        {
                            fprintf(stderr,
                                    "Warning: Skipping overlapping variants at %s:%"PRId64", sample %s\t%s>%s.\n",
                                    chr,(int64_t) rec->pos+1,args->hdr->samples[args->smpl->idx[ismpl]],rec->d.allele[0],rec->d.allele[ial]);
                            if ( !overlaps_warned )
                                fprintf(stderr,"         This message is printed only once, the verbosity can be increased with `--verbose 2`\n");
                            overlaps_warned = 1;
                        }
                        if ( args->out  )
                            fprintf(args->out,"LOG\tWarning: Skipping overlapping variants at %s:%"PRId64", sample %s\t%s>%s\n",
                                    chr,(int64_t) rec->pos+1,args->hdr->samples[args->smpl->idx[ismpl]],rec->d.allele[0],rec->d.allele[ial]);
                    }
                    hap_destroy(child);
                    continue;
                }
                if ( child->type==HAP_SSS )
                {
                    csq_t csq;
                    memset(&csq, 0, sizeof(csq_t));
                    csq.pos          = rec->pos;
                    csq.type.biotype = tr->type;
                    csq.type.strand  = tr->strand;
                    csq.type.trid    = tr->id;
                    csq.type.vcf_ial = ial;
                    csq.type.gene    = tr->gene->name;
                    csq.type.type = child->csq;
                    csq_stage(args, &csq, rec);
                    hap_destroy(child);
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
                TSCRIPT_AUX(tr)->hap[i] = child;
                TSCRIPT_AUX(tr)->hap[i]->nend++;
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
    if ( csq_push(args, csq, rec)!=0 && args->phase==PHASE_DROP_GT ) return;    // the consequence already exists

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
                if ( gt[j]==bcf_gt_missing || gt[j]==bcf_int32_vector_end ) continue;
                int ial = bcf_gt_allele(gt[j]);
                if ( !ial || ial!=csq->type.vcf_ial ) continue;
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
            if ( gt[j]==bcf_gt_missing || gt[j]==bcf_int32_vector_end ) continue;
            int ial = bcf_gt_allele(gt[j]);
            if ( !ial || ial!=csq->type.vcf_ial ) continue;

            int icsq2 = 2*csq->idx + j;
            if ( icsq2 >= args->ncsq2_max ) // more than ncsq_max consequences, so can't fit it in FMT
            {
                int ismpl = args->smpl->idx[i];
                if ( args->verbosity && (!args->ncsq2_small_warned || args->verbosity > 1) )
                {
                    fprintf(stderr,
                            "Warning: Too many consequences for sample %s at %s:%"PRId64", keeping the first %d and skipping the rest.\n",
                            args->hdr->samples[ismpl],bcf_hdr_id2name(args->hdr,args->rid),(int64_t) vrec->line->pos+1,icsq2+1);
                    if ( !args->ncsq2_small_warned )
                        fprintf(stderr,"         The limit can be increased by setting the --ncsq parameter. This warning is printed only once.\n");
                    args->ncsq2_small_warned = 1;
                }
                if ( args->ncsq2_small_warned < icsq2 ) args->ncsq2_small_warned = icsq2;
                break;
            }
            int ival, ibit;
            icsq2_to_bit(icsq2, &ival,&ibit);
            if ( vrec->nfmt < 1 + ival ) vrec->nfmt = 1 + ival;
            vrec->fmt_bm[i*args->nfmt_bcsq + ival] |= 1 << ibit;
        }
    }
}
int test_utr(args_t *args, bcf1_t *rec)
{
    const char *chr = drop_chr_prefix(args, bcf_seqname(args->hdr,rec));
    // note that the off-by-one extension of rlen is deliberate to account for insertions
    if ( !regidx_overlap(args->idx_utr,chr,rec->pos,rec->pos+rec->rlen, args->itr) ) return 0;

    splice_t splice;
    splice_init(&splice, rec);

    int i, ret = 0;
    while ( regitr_overlap(args->itr) )
    {
        gf_utr_t *utr = regitr_payload(args->itr, gf_utr_t*);
        gf_tscript_t *tr = splice.tr = utr->tr;
        for (i=1; i<rec->n_allele; i++)
        {
            if ( rec->d.allele[i][0]=='<' || rec->d.allele[i][0]=='*' ) { continue; }
            splice.vcf.alt = rec->d.allele[i];
            splice.csq     = 0;
            int splice_ret = splice_csq(args, &splice, utr->beg, utr->end);
            if ( splice_ret!=SPLICE_INSIDE && splice_ret!=SPLICE_OVERLAP ) continue;
            csq_t csq;
            memset(&csq, 0, sizeof(csq_t));
            csq.pos          = rec->pos;
            csq.type.type    = utr->which==prime5 ? CSQ_UTR5 : CSQ_UTR3;
            csq.type.biotype = tr->type;
            csq.type.strand  = tr->strand;
            csq.type.trid    = tr->id;
            csq.type.vcf_ial = i;
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
    const char *chr = drop_chr_prefix(args, bcf_seqname(args->hdr,rec));
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
            splice.csq     = 0;
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
    const char *chr = drop_chr_prefix(args, bcf_seqname(args->hdr,rec));
    if ( !regidx_overlap(args->idx_tscript,chr,rec->pos,rec->pos+rec->rlen, args->itr) ) return 0;

    splice_t splice;
    splice_init(&splice, rec);

    int i, ret = 0;
    while ( regitr_overlap(args->itr) )
    {
        gf_tscript_t *tr = splice.tr = regitr_payload(args->itr, gf_tscript_t*);
        for (i=1; i<rec->n_allele; i++)
        {
            if ( rec->d.allele[i][0]=='<' || rec->d.allele[i][0]=='*' ) { continue; }
            splice.vcf.alt = rec->d.allele[i];
            splice.csq     = 0;
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

void test_symbolic_alt(args_t *args, bcf1_t *rec)
{
    static int warned = 0;
    if ( args->verbosity && (!warned && args->verbosity > 0) )
    {
        fprintf(stderr,"Warning: The support for symbolic ALT insertions is experimental.\n");
        warned = 1;
    }

    const char *chr = drop_chr_prefix(args, bcf_seqname(args->hdr,rec));

    // only insertions atm
    int beg = rec->pos + 1;
    int end = beg;
    int csq_class = CSQ_ELONGATION;

    int hit = 0;
    if ( regidx_overlap(args->idx_cds,chr,beg,end, args->itr) )
    {
        while ( regitr_overlap(args->itr) )
        {
            csq_t csq;
            memset(&csq, 0, sizeof(csq_t));
            gf_cds_t *cds    = regitr_payload(args->itr,gf_cds_t*);
            gf_tscript_t *tr = cds->tr;
            csq.type.type    = (GF_is_coding(tr->type) ? CSQ_CODING_SEQUENCE : CSQ_NON_CODING) | csq_class;
            csq.pos          = rec->pos;
            csq.type.biotype = tr->type;
            csq.type.strand  = tr->strand;
            csq.type.trid    = tr->id;
            csq.type.gene    = tr->gene->name;
            csq_stage(args, &csq, rec);
            hit = 1;
        }
    }
    if ( regidx_overlap(args->idx_utr,chr,beg,end, args->itr) )
    {
        while ( regitr_overlap(args->itr) )
        {
            csq_t csq;
            memset(&csq, 0, sizeof(csq_t));
            gf_utr_t *utr    = regitr_payload(args->itr, gf_utr_t*);
            gf_tscript_t *tr = utr->tr;
            csq.type.type    = (utr->which==prime5 ? CSQ_UTR5 : CSQ_UTR3) | csq_class;
            csq.pos          = rec->pos;
            csq.type.biotype = tr->type;
            csq.type.strand  = tr->strand;
            csq.type.trid    = tr->id;
            csq.type.gene    = tr->gene->name;
            csq_stage(args, &csq, rec);
            hit = 1;
        }
    }
    if ( regidx_overlap(args->idx_exon,chr,beg,end, args->itr) )
    {
        splice_t splice;
        splice_init(&splice, rec);
        splice.check_acceptor = splice.check_donor = 1;

        while ( regitr_overlap(args->itr) )
        {
            gf_exon_t *exon = regitr_payload(args->itr, gf_exon_t*);
            splice.tr = exon->tr;
            if ( !splice.tr->ncds ) continue;  // not a coding transcript, no interest in splice sites
            splice.check_region_beg = splice.tr->beg==exon->beg ? 0 : 1;
            splice.check_region_end = splice.tr->end==exon->end ? 0 : 1;
            splice.vcf.alt = rec->d.allele[1];
            splice.csq     = csq_class;
            splice_csq(args, &splice, exon->beg, exon->end);
            if ( splice.csq ) hit = 1;
        }
    }
    if ( !hit && regidx_overlap(args->idx_tscript,chr,beg,end, args->itr) )
    {
        splice_t splice;
        splice_init(&splice, rec);

        while ( regitr_overlap(args->itr) )
        {
            csq_t csq;
            memset(&csq, 0, sizeof(csq_t));
            gf_tscript_t *tr = splice.tr = regitr_payload(args->itr, gf_tscript_t*);
            splice.vcf.alt = rec->d.allele[1];
            splice.csq     = csq_class;
            int splice_ret = splice_csq(args, &splice, tr->beg, tr->end);
            if ( splice_ret!=SPLICE_INSIDE && splice_ret!=SPLICE_OVERLAP ) continue;    // SPLICE_OUTSIDE or SPLICE_REF
            csq.type.type    = (GF_is_coding(tr->type) ? CSQ_INTRON : CSQ_NON_CODING) | csq_class;
            csq.pos          = rec->pos;
            csq.type.biotype = tr->type;
            csq.type.strand  = tr->strand;
            csq.type.trid    = tr->id;
            csq.type.gene    = tr->gene->name;
            csq_stage(args, &csq, rec);
        }
    }
}

void debug_print_buffers(args_t *args, int pos)
{
    int i,j;
    fprintf(stderr,"debug_print_buffers at %d\n", pos);
    fprintf(stderr,"vbufs:\n");
    for (i=0; i<args->vcf_rbuf.n; i++)
    {
        int k = rbuf_kth(&args->vcf_rbuf, i);
        vbuf_t *vbuf = args->vcf_buf[k];

        fprintf(stderr,"\tvbuf %d:\n", i);
        for (j=0; j<vbuf->n; j++)
        {
            vrec_t *vrec = vbuf->vrec[j];
            fprintf(stderr,"\t\t%"PRId64" .. nvcsq=%d\n", (int64_t) vrec->line->pos+1, vrec->nvcsq);
        }
    }
    fprintf(stderr,"pos2vbuf:");
    khint_t k;
    for (k = 0; k < kh_end(args->pos2vbuf); ++k)
        if (kh_exist(args->pos2vbuf, k)) fprintf(stderr," %d",1+(int)kh_key(args->pos2vbuf, k));
    fprintf(stderr,"\n");
    fprintf(stderr,"active_tr: %d\n", args->active_tr->ndat);
}

static void process(args_t *args, bcf1_t **rec_ptr)
{
    if ( !rec_ptr )
    {
        hap_flush(args, REGIDX_MAX);
        vbuf_flush(args, REGIDX_MAX);
        return;
    }

    bcf1_t *rec = *rec_ptr;
    static int32_t prev_rid = -1, prev_pos = -1;
    if ( prev_rid!=rec->rid )
    {
        prev_rid = rec->rid;
        prev_pos = rec->pos;

        // Common error is to use different naming conventions in the fasta and the VCF (e.g. X vs chrX).
        // Perform a simple sanity check (that does not catch much), the chromosome must be present in the
        // reference file
        if ( !faidx_has_seq(args->fai,bcf_seqname(args->hdr,rec)) )
        {
            if ( !faidx_has_seq(args->fai,drop_chr_prefix(args,bcf_seqname(args->hdr,rec))) && !faidx_has_seq(args->fai,add_chr_prefix(args,bcf_seqname(args->hdr,rec))) )
                error("Error: the chromosome \"%s\" is not present in %s\n",bcf_seqname(args->hdr,rec),args->fa_fname);
        }
    }
    if ( prev_pos > rec->pos )
        error("Error: The file is not sorted, %s:%d comes before %s:%"PRId64"\n",bcf_seqname(args->hdr,rec),prev_pos+1,bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);

    int call_csq = 1;
    if ( rec->n_allele < 2 ) call_csq = 0;   // no alternate allele
    else if ( rec->n_allele==2 && (rec->d.allele[1][0]=='*' || rec->d.allele[1][1]=='*') ) call_csq = 0;     // gVCF, not an alt allele
    else if ( rec->d.allele[1][0]=='<' )
    {
        if ( strncmp("<INS",rec->d.allele[1], 4) ) call_csq = 0;    // only <INS[:.*]> is supported at the moment
    }
    if ( call_csq && args->filter )
    {
        call_csq = filter_test(args->filter, rec, NULL);
        if ( args->filter_logic==FLT_EXCLUDE ) call_csq = call_csq ? 0 : 1;
    }
    if ( !call_csq )
    {
        if ( !args->out_fh ) return;    // not a VCF output
        vbuf_push(args, rec_ptr);
        hap_flush(args, rec->pos-1);
        vbuf_flush(args, rec->pos-1);
        return;
    }

    if ( args->rid != rec->rid )
    {
        hap_flush(args, REGIDX_MAX);
        vbuf_flush(args, REGIDX_MAX);
    }
    args->rid = rec->rid;
    vbuf_t *vbuf = vbuf_push(args, rec_ptr);

    if ( rec->d.allele[1][0]!='<' )
    {
        int hit = args->local_csq ? test_cds_local(args, rec) : test_cds(args, rec, vbuf);
        hit += test_utr(args, rec);
        hit += test_splice(args, rec);
        if ( !hit ) test_tscript(args, rec);
    }
    else
        test_symbolic_alt(args, rec);

    if ( rec->pos > 0 )
    {
        hap_flush(args, rec->pos-1);
        vbuf_flush(args, rec->pos-1);
    }

    return;
}

static const char *usage(void)
{
    return
        "\n"
        "About: Haplotype-aware consequence caller.\n"
        "Usage: bcftools csq [OPTIONS] in.vcf\n"
        "\n"
        "Required options:\n"
        "   -f, --fasta-ref FILE              Reference file in fasta format\n"
        "   -g, --gff-annot FILE              GFF3 annotation file\n"
        "\n"
        "CSQ options:\n"
        "   -B, --trim-protein-seq INT        Abbreviate protein-changing predictions to max INT aminoacids\n"
        "   -c, --custom-tag STRING           Use this tag instead of the default BCSQ\n"
        "   -l, --local-csq                   Localized predictions, consider only one VCF record at a time\n"
        "   -n, --ncsq INT                    Maximum number of per-haplotype consequences to consider for each site [15]\n"
        "   -p, --phase a|m|r|R|s             How to handle unphased heterozygous genotypes: [r]\n"
        "                                       a: take GTs as is, create haplotypes regardless of phase (0/1 -> 0|1)\n"
        "                                       m: merge *all* GTs into a single haplotype (0/1 -> 1, 1/2 -> 1)\n"
        "                                       r: require phased GTs, throw an error on unphased het GTs\n"
        "                                       R: create non-reference haplotypes if possible (0/1 -> 1|1, 1/2 -> 1|2)\n"
        "                                       s: skip unphased hets\n"
        "GFF options:\n"
        "       --dump-gff FILE.gz            Dump the parsed GFF file (for debugging purposes)\n"
        "       --force                       Run even if some sanity checks fail\n"
        "       --unify-chr-names 1|0         Automatically unify chromosome naming (e.g. chrX vs X) in GFF, fasta, and VCF [1]\n"
        "General options:\n"
        "   -e, --exclude EXPR                Exclude sites for which the expression is true\n"
        "   -i, --include EXPR                Select sites for which the expression is true\n"
        "       --no-version                  Do not append version and command line to the header\n"
        "   -o, --output FILE                 Write output to a file [standard output]\n"
        "   -O, --output-type b|u|z|v|t[0-9]  b: compressed BCF, u: uncompressed BCF, z: compressed VCF\n"
        "                                     v: uncompressed VCF, t: plain tab-delimited text output, 0-9: compression level [v]\n"
        "   -r, --regions REGION              Restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE           Restrict to regions listed in a file\n"
        "       --regions-overlap 0|1|2       Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n"
        "   -s, --samples -|LIST              Samples to include or \"-\" to apply all variants and ignore samples\n"
        "   -S, --samples-file FILE           Samples to include\n"
        "   -t, --targets REGION              Similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE           Similar to -R but streams rather than index-jumps\n"
        "       --targets-overlap 0|1|2       Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n"
        "       --threads INT                 Use multithreading with <int> worker threads [0]\n"
        "   -v, --verbose INT                 Verbosity level 0-2 [1]\n"
        "   -W, --write-index[=FMT]           Automatically index the output files [off]\n"
        "\n"
        "Example:\n"
        "   bcftools csq -f hs37d5.fa -g Homo_sapiens.GRCh37.87.gff3.gz in.vcf\n"
        "\n"
        "   # GFF3 annotation files can be downloaded from Ensembl. e.g. for human:\n"
        "   http://ftp.ensembl.org/pub/current_gff3/homo_sapiens/\n"
        "   http://ftp.ensembl.org/pub/grch37/current/gff3/homo_sapiens/\n"
        "\n";
}

int main_csq(int argc, char *argv[])
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc = argc; args->argv = argv;
    args->output_type = FT_VCF;
    args->bcsq_tag = "BCSQ";
    args->ncsq2_max = 2*(16-1);      // 1 bit is reserved for BCF missing values
    args->verbosity = 1;
    args->record_cmd_line = 1;
    args->clevel = -1;
    args->unify_chr_names = 1;

    static struct option loptions[] =
    {
        {"force",0,0,1},
        {"threads",required_argument,NULL,2},
        {"help",0,0,'h'},
        {"ncsq",1,0,'n'},
        {"brief-predictions",no_argument,0,'b'},
        {"trim-protein-seq",required_argument,0,'B'},
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
        {"verbose",1,0,'v'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"regions-overlap",required_argument,NULL,4},
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"targets-overlap",required_argument,NULL,5},
        {"no-version",no_argument,NULL,3},
        {"write-index",optional_argument,NULL,'W'},
        {"dump-gff",required_argument,NULL,7},
        {"unify-chr-names",required_argument,NULL,8},
        {0,0,0,0}
    };
    int c, targets_is_file = 0, regions_is_file = 0;
    int regions_overlap = 1;
    int targets_overlap = 0;
    char *targets_list = NULL, *regions_list = NULL, *tmp;
    while ((c = getopt_long(argc, argv, "?hr:R:t:T:i:e:f:o:O:g:s:S:p:qc:ln:bB:v:W::",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case  1 : args->force = 1; break;
            case  2 :
                args->n_threads = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --threads  %s\n", optarg);
                break;
            case  3 : args->record_cmd_line = 0; break;
            case 'b':
                    args->brief_predictions = 1;
                    fprintf(stderr,"Warning: The -b option will be removed in future versions. Please use -B 1 instead.\n");
                    break;
            case 'B':
                    args->brief_predictions = strtol(optarg,&tmp,10);
                    if ( *tmp || args->brief_predictions<1 ) error("Could not parse argument: --trim-protein-seq %s\n", optarg);
                    break;
            case 'l': args->local_csq = 1; break;
            case 'c': args->bcsq_tag = optarg; break;
            case 'q': error("Error: the -q option has been deprecated, use -v, --verbose instead.\n"); break;
            case 'v':
                args->verbosity = atoi(optarg);
                if ( args->verbosity<0 || args->verbosity>2 ) error("Error: expected integer 0-2 with -v, --verbose\n");
                break;
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
                args->ncsq2_max = 2 * atoi(optarg);
                if ( args->ncsq2_max <= 0 ) error("Expected positive integer with -n, got %s\n", optarg);
                break;
            case 'o': args->output_fname = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 't': args->output_type = FT_TAB_TEXT; break;
                          case 'b': args->output_type = FT_BCF_GZ; break;
                          case 'u': args->output_type = FT_BCF; break;
                          case 'z': args->output_type = FT_VCF_GZ; break;
                          case 'v': args->output_type = FT_VCF; break;
                          default:
                          {
                              args->clevel = strtol(optarg,&tmp,10);
                              if ( *tmp || args->clevel<0 || args->clevel>9 ) error("The output type \"%s\" not recognised\n", optarg);
                          }
                      }
                      if ( optarg[1] )
                      {
                          args->clevel = strtol(optarg+1,&tmp,10);
                          if ( *tmp || args->clevel<0 || args->clevel>9 ) error("Could not parse argument: --output-type %s\n", optarg+1);
                      }
                      break;
            case 'e':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'r': regions_list = optarg; break;
            case 'R': regions_list = optarg; regions_is_file = 1; break;
            case 's': args->sample_list = optarg; break;
            case 'S': args->sample_list = optarg; args->sample_is_file = 1; break;
            case 't': targets_list = optarg; break;
            case 'T': targets_list = optarg; targets_is_file = 1; break;
            case  4 :
                regions_overlap = parse_overlap_option(optarg);
                if ( regions_overlap < 0 ) error("Could not parse: --regions-overlap %s\n",optarg);
                break;
            case  5 :
                targets_overlap = parse_overlap_option(optarg);
                if ( targets_overlap < 0 ) error("Could not parse: --targets-overlap %s\n",optarg);
                break;
            case 'W':
                if (!(args->write_index = write_index_parse(optarg)))
                    error("Unsupported index format '%s'\n", optarg);
                break;
            case  7 : args->dump_gff = optarg; break;
            case  8 :
                if ( !strcmp(optarg,"0") ) args->unify_chr_names = 0;
                else if ( !strcmp(optarg,"1") ) args->unify_chr_names = 1;
                else error("Could not parse: --unify-chr-names %s\n",optarg);
                break;
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
    if ( !args->fa_fname ) error("Missing the --fasta-ref option\n");
    if ( !args->gff_fname ) error("Missing the --gff option\n");
    args->sr = bcf_sr_init();
    if ( targets_list )
    {
        bcf_sr_set_opt(args->sr,BCF_SR_TARGETS_OVERLAP,targets_overlap);
        if ( bcf_sr_set_targets(args->sr, targets_list, targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", targets_list);
    }
    if ( regions_list )
    {
        bcf_sr_set_opt(args->sr,BCF_SR_REGIONS_OVERLAP,regions_overlap);
        if ( bcf_sr_set_regions(args->sr, regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", regions_list);
    }
    if ( bcf_sr_set_threads(args->sr, args->n_threads)<0 ) error("Failed to create %d extra threads\n", args->n_threads);
    if ( !bcf_sr_add_reader(args->sr, fname) )
        error("Failed to read from %s: %s\n", !strcmp("-",fname)?"standard input":fname,bcf_sr_strerror(args->sr->errnum));
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

