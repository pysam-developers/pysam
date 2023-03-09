/*  vcfstats.c -- Produces stats which can be plotted using plot-vcfstats.

    Copyright (C) 2012-2023 Genome Research Ltd.

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
THE SOFTWARE.  */

/*
    Notes and known issues:
        - SN ts/tv calculation includes all non-ref alleles listed in ALT while per-sample ts/tv
        takes the first non-ref allele only, something to consider with many non-ref HETs.
*/
#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <htslib/faidx.h>
#include <inttypes.h>
#include "bcftools.h"
#include "filter.h"
#include "bin.h"
#include "dist.h"

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define HWE_STATS 1
#define QUAL_STATS 1
#define IRC_STATS 1
#define IRC_RLEN 10
#define NA_STRING "0"

typedef struct
{
    char *tag;
    float min, max;
    uint64_t *vals_ts, *vals_tv;
    void *val;
    int nbins, type, m_val, idx;
}
user_stats_t;

typedef struct
{
    int min, max, step, m_vals;
    uint64_t *vals;
}
idist_t;

typedef struct
{
    uint64_t n_snps, n_indels, n_mnps, n_others, n_mals, n_snp_mals, n_records, n_noalts;
    int *af_ts, *af_tv, *af_snps;   // first bin of af_* stats are singletons
    #if HWE_STATS
        int *af_hwe;
    #endif
    #if IRC_STATS
        int n_repeat[IRC_RLEN][4], n_repeat_na;    // number of indels which are repeat-consistent, repeat-inconsistent (dels and ins), and not applicable
        int *af_repeats[3];
    #endif
    int ts_alt1, tv_alt1;
    #if QUAL_STATS
        // Values are rounded to one significant digit and 1 is added (Q*10+1); missing and negative values go in the first bin
        // Only SNPs that are the 1st alternate allele are counted
        dist_t *qual_ts, *qual_tv, *qual_indels;
    #endif
    int *insertions, *deletions, m_indel;   // maximum indel length
    int in_frame, out_frame, na_frame, in_frame_alt1, out_frame_alt1, na_frame_alt1;
    int subst[15];
    int *smpl_hets, *smpl_homRR, *smpl_homAA, *smpl_ts, *smpl_tv, *smpl_indels, *smpl_ndp, *smpl_sngl;
    int *smpl_hapRef, *smpl_hapAlt, *smpl_missing;
    int *smpl_ins_hets, *smpl_del_hets, *smpl_ins_homs, *smpl_del_homs;
    int *smpl_frm_shifts; // not-applicable, in-frame, out-frame
    unsigned long int *smpl_dp;
    idist_t dp, dp_sites;
    int nusr;
    user_stats_t *usr;
    double *dvaf;   // distribution of the mean indel-allele frequency by length: -m_indel,-(m_indel-1),...-1,0,1,..,m_indel
    uint32_t *nvaf;
}
stats_t;

typedef struct
{
    uint64_t gt2gt[5][5];   // number of RR->RR, RR->RA, etc. matches/mismatches; see type2stats
    /*
        Pearson's R^2 is used for aggregate R^2
        y, yy .. sum of dosage and squared dosage in the query VCF (second file)
        x, xx .. sum of squared dosage in the truth VCF (first file)
        n     .. number of genotypes
     */
    double y, yy, x, xx, yx, n;
}
gtcmp_t;

typedef struct
{
    char *seq;
    int pos, cnt, len;
}
_idc1_t;
typedef struct
{
    faidx_t *ref;
    _idc1_t *dat;
    int ndat, mdat;
}
indel_ctx_t;

typedef struct
{
    // stats
    stats_t stats[3];
    int *tmp_iaf, ntmp_iaf, m_af, m_qual, naf_hwe, mtmp_frm;
    uint8_t *tmp_frm;
    int dp_min, dp_max, dp_step;
    gtcmp_t *smpl_gts_snps, *smpl_gts_indels;
    gtcmp_t *af_gts_snps, *af_gts_indels; // first bin of af_* stats are singletons
    bin_t *af_bins;
    float *farr;
    int mfarr;

    // indel context
    indel_ctx_t *indel_ctx;
    char *ref_fname;

    // user stats
    int nusr;
    user_stats_t *usr;

    // other
    bcf_srs_t *files;
    bcf_sr_regions_t *exons;
    char **argv, *exons_fname, *regions_list, *samples_list, *targets_list, *af_bins_list, *af_tag;
    int argc, verbose_sites, first_allele_only, samples_is_file;
    int split_by_id, nstats;

    filter_t *filter[2];
    char *filter_str;
    int filter_logic;   // include or exclude sites which match the filters? One of FLT_INCLUDE/FLT_EXCLUDE
    int n_threads;
}
args_t;

static int type2dosage[6], type2ploidy[6], type2stats[7];

static void idist_init(idist_t *d, int min, int max, int step)
{
    d->min = min; d->max = max; d->step = step;
    d->m_vals = 4 + (d->max - d->min)/d->step;
    d->vals = (uint64_t*) calloc(d->m_vals,sizeof(uint64_t));
}
static void idist_destroy(idist_t *d)
{
    if ( d->vals ) free(d->vals);
}
static inline uint64_t *idist(idist_t *d, int val)
{
    if ( val < d->min ) return &d->vals[0];
    if ( val > d->max ) return &d->vals[d->m_vals-1];
    return &d->vals[1 + (val - d->min) / d->step];
}
static inline int idist_i2bin(idist_t *d, int i)
{
    if ( i<=0 ) return d->min;
    if ( i>= d->m_vals ) return d->max;
    return i-1+d->min;
}

#define IC_DBG 0
#if IC_DBG
static void _indel_ctx_print1(_idc1_t *idc)
{
    int i;
    fprintf(stdout, "%d\t", idc->cnt);
    for (i=0; i<idc->len; i++)
        fputc(idc->seq[i], stdout);
    fputc('\n', stdout);
}
static void _indel_ctx_print(indel_ctx_t *ctx)
{
    int i;
    for (i=0; i<ctx->ndat; i++)
        _indel_ctx_print1(&ctx->dat[i]);
    fputc('\n',stdout);
}
#endif
static int _indel_ctx_lookup(indel_ctx_t *ctx, char *seq, int seq_len, int *hit)
{
    // binary search
    int min = 0, max = ctx->ndat - 1;
    while ( min<=max )
    {
        int i = (min+max)/2;
        int cmp = strncmp(seq, ctx->dat[i].seq, seq_len);
        if ( cmp<0 ) max = i - 1;
        else if ( cmp>0 ) min = i + 1;
        else
        {
            if ( seq_len==ctx->dat[i].len )
            {
                *hit = 1;
                return i;
            }
            else if ( seq_len<ctx->dat[i].len ) max = i - 1;
            else min = i + 1;
        }
    }
    *hit = 0;
    return max;
}
static void _indel_ctx_insert(indel_ctx_t *ctx, char *seq, int seq_len, int pos)
{
    int idat, hit, i;
    idat = _indel_ctx_lookup(ctx, seq, seq_len, &hit);
    if ( !hit )
    {
        if ( pos>0 ) return;
        idat++;
        ctx->ndat++;
        hts_expand(_idc1_t, ctx->ndat+1, ctx->mdat, ctx->dat);
        if ( idat<ctx->ndat && ctx->ndat>1 )
            memmove(&ctx->dat[idat+1], &ctx->dat[idat], (ctx->ndat - idat - 1)*sizeof(_idc1_t));
        ctx->dat[idat].len = seq_len;
        ctx->dat[idat].cnt = 1;
        ctx->dat[idat].pos = pos;
        ctx->dat[idat].seq = (char*) malloc(sizeof(char)*(seq_len+1));
        for (i=0; i<seq_len; i++) ctx->dat[idat].seq[i] = seq[i];
        ctx->dat[idat].seq[i] = 0;
        return;
    }
    if ( ctx->dat[idat].pos + seq_len == pos )
    {
        ctx->dat[idat].cnt++;
        ctx->dat[idat].pos = pos;
    }
}
indel_ctx_t *indel_ctx_init(char *fa_ref_fname)
{
    indel_ctx_t *ctx = (indel_ctx_t *) calloc(1,sizeof(indel_ctx_t));
    ctx->ref = fai_load(fa_ref_fname);
    if ( !ctx->ref )
    {
        free(ctx);
        return NULL;
    }
    return ctx;
}
void indel_ctx_destroy(indel_ctx_t *ctx)
{
    fai_destroy(ctx->ref);
    if ( ctx->mdat ) free(ctx->dat);
    free(ctx);
}
/**
 * indel_ctx_type() - determine indel context type
 * @ctx:
 * @chr:    chromosome name
 * @pos:    position of the first @ref base, 1-based
 * @ref:    reference allele
 * @alt:    alternate allele. Only first of multiple comma-separated alleles is
 *          considered
 * @nrep:   number of repeated elements (w)
 * @nlen:   length of a single repeat element (w)
 *
 * Returns the INDEL length, negative for deletions, positive for insertions
 */
int indel_ctx_type(indel_ctx_t *ctx, char *chr, int pos, char *ref, char *alt, int *nrep, int *nlen)
{
    const int win_size = 50;             // hard-wired for now
    const int rep_len  = IRC_RLEN;       // hard-wired for now

    int ref_len = strlen(ref);
    int alt_len = 0;
    while ( alt[alt_len] && alt[alt_len]!=',' ) alt_len++;

    int i, fai_ref_len;
    char *fai_ref = faidx_fetch_seq(ctx->ref, chr, pos-1, pos+win_size, &fai_ref_len);
    for (i=0; i<fai_ref_len; i++)
        if ( (int)fai_ref[i]>96 ) fai_ref[i] -= 32;

    // Sanity check: the reference sequence must match the REF allele
    for (i=0; i<fai_ref_len && i<ref_len; i++)
        if ( ref[i] != fai_ref[i] && ref[i] - 32 != fai_ref[i] && !iupac_consistent(fai_ref[i], ref[i]) )
            error("\nSanity check failed, the reference sequence differs: %s:%d+%d .. %c vs %c\n", chr, pos, i, ref[i],fai_ref[i]);

    // Count occurrences of all possible kmers
    ctx->ndat = 0;
    for (i=0; i<win_size; i++)
    {
        int k, kmax = rep_len <= i ? rep_len : i+1;
        for (k=0; k<kmax; k++)
            _indel_ctx_insert(ctx, &fai_ref[i-k+1], k+1, i-k);
    }

    #if IC_DBG
    fprintf(stdout,"ref: %s\n", ref);
    fprintf(stdout,"alt: %s\n", alt);
    fprintf(stdout,"ctx: %s\n", fai_ref);
    _indel_ctx_print(ctx);
    #endif

    int max_cnt = 0, max_len = 0;
    for (i=0; i<ctx->ndat; i++)
    {
        if ( max_cnt < ctx->dat[i].cnt || (max_cnt==ctx->dat[i].cnt && max_len < ctx->dat[i].len) )
        {
            max_cnt = ctx->dat[i].cnt;
            max_len = ctx->dat[i].len;
        }
        free(ctx->dat[i].seq);
    }
    free(fai_ref);

    *nrep = max_cnt;
    *nlen = max_len;
    return alt_len - ref_len;
}

static void add_user_stats(args_t *args, char *str)
{
    args->nusr++;
    args->usr = (user_stats_t*) realloc(args->usr,sizeof(user_stats_t)*args->nusr);
    user_stats_t *usr = &args->usr[args->nusr-1];
    memset(usr,0,sizeof(*usr));
    usr->min   = 0;
    usr->max   = 1;
    usr->nbins = 100;
    usr->idx   = 0;

    char *tmp = str;
    while ( *tmp && *tmp!=':' ) tmp++;

    // Tag with an index or just tag? (e.g. PV4[1] vs DP)
    if ( tmp > str && tmp[-1]==']' )
    {
        char *ptr = tmp;
        while ( ptr>str && *ptr!='[' ) ptr--;
        if ( *ptr=='[' )
        {
            char *ptr2;
            usr->idx = strtol(ptr+1, &ptr2, 10);
            if ( ptr+1==ptr2 || ptr2 != tmp-1 ) error("Could not parse the index in \"%s\" (ptr=%s;ptr2=%s(%p),tmp=%s(%p),idx=%d)\n", str,ptr,ptr2,ptr2,tmp,tmp,usr->idx);
            if ( usr->idx<0 ) error("Error: negative index is not allowed: \"%s\"\n", str);
            *ptr = 0;
        }
    }

    usr->tag = (char*)calloc(tmp-str+2,sizeof(char));
    memcpy(usr->tag,str,tmp-str);

    if ( *tmp )
    {
        char *ptr = ++tmp;
        usr->min = strtod(tmp, &ptr);
        if ( tmp==ptr ) error("Could not parse %s\n", str);
        tmp = ptr+1;
    }
    if ( *tmp )
    {
        char *ptr = tmp;
        usr->max = strtod(tmp, &ptr);
        if ( tmp==ptr ) error("Could not parse %s\n", str);
        tmp = ptr+1;
    }
    if ( *tmp )
    {
        char *ptr = tmp;
        usr->nbins = strtol(tmp, &ptr, 10);
        if ( tmp==ptr ) error("Could not parse %s\n", str);
        if ( usr->nbins<=0 ) error("Number of bins does not make sense (%d): %s.\n", usr->nbins, str);
    }
}
static void init_user_stats(args_t *args, bcf_hdr_t *hdr, stats_t *stats)
{
    stats->nusr = args->nusr;
    stats->usr = (user_stats_t*)malloc(sizeof(user_stats_t)*args->nusr);
    memcpy(stats->usr,args->usr,args->nusr*sizeof(user_stats_t));
    int i;
    for (i=0; i<stats->nusr; i++)
    {
        user_stats_t *usr = &stats->usr[i];
        usr->vals_ts = (uint64_t*)calloc(usr->nbins,sizeof(uint64_t));
        usr->vals_tv = (uint64_t*)calloc(usr->nbins,sizeof(uint64_t));
        int id = bcf_hdr_id2int(hdr,BCF_DT_ID,usr->tag);
        if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_INFO,id) ) error("The INFO tag \"%s\" is not defined in the header\n", usr->tag);
        usr->type = bcf_hdr_id2type(hdr,BCF_HL_INFO,id);
        if ( usr->type!=BCF_HT_REAL && usr->type!=BCF_HT_INT ) error("The INFO tag \"%s\" is not of Float or Integer type (%d)\n", usr->tag, usr->type);
    }
}
static void init_stats(args_t *args)
{
    int i;
    args->nstats = args->files->nreaders==1 ? 1 : 3;
    if ( args->split_by_id ) args->nstats = 2;

    if ( args->filter_str )
    {
        args->filter[0] = filter_init(bcf_sr_get_header(args->files,0), args->filter_str);
        if ( args->files->nreaders==2 )
            args->filter[1] = filter_init(bcf_sr_get_header(args->files,1), args->filter_str);
        args->files->max_unpack |= filter_max_unpack(args->filter[0]);
    }

    // AF corresponds to AC but is more robust to mixtures of haploid and diploid GTs
    if ( !args->af_bins_list )
    {
        args->m_af = 101;
        for (i=0; i<args->files->nreaders; i++)
            if ( bcf_hdr_nsamples(args->files->readers[i].header) + 1> args->m_af )
                args->m_af = bcf_hdr_nsamples(args->files->readers[i].header) + 1;
    }
    else
    {
        args->af_bins = bin_init(args->af_bins_list,0,1);

        // m_af is used also for other af arrays, where the first bin is for
        // singletons. However, since the last element is unused in af_bins
        // (n boundaries form n-1 intervals), the m_af count is good for both.
        args->m_af = bin_get_size(args->af_bins);
    }

    bcf_hdr_t *hdr = bcf_sr_get_header(args->files,0);
    if ( args->af_tag && !bcf_hdr_idinfo_exists(hdr,BCF_HL_INFO,bcf_hdr_id2int(hdr,BCF_DT_ID,args->af_tag)) )
        error("No such INFO tag: %s\n", args->af_tag);

    #if QUAL_STATS
        args->m_qual = 999;
    #endif
    #if HWE_STATS
        args->naf_hwe = 100;
    #endif

    if ( args->samples_list )
    {
        if ( !bcf_sr_set_samples(args->files,args->samples_list,args->samples_is_file) )
        {
            if ( !bcf_hdr_nsamples(args->files->readers[0].header) )
                error("No sample columns in %s\n", args->files->readers[0].fname);
            error("Unable to parse the samples: \"%s\"\n", args->samples_list);
        }
        args->af_gts_snps     = (gtcmp_t *) calloc(args->m_af,sizeof(gtcmp_t));
        args->af_gts_indels   = (gtcmp_t *) calloc(args->m_af,sizeof(gtcmp_t));
        args->smpl_gts_snps   = (gtcmp_t *) calloc(args->files->n_smpl,sizeof(gtcmp_t));
        args->smpl_gts_indels = (gtcmp_t *) calloc(args->files->n_smpl,sizeof(gtcmp_t));
    }
    for (i=0; i<args->nstats; i++)
    {
        stats_t *stats = &args->stats[i];
        stats->m_indel     = 60;
        stats->insertions  = (int*) calloc(stats->m_indel,sizeof(int));
        stats->deletions   = (int*) calloc(stats->m_indel,sizeof(int));
        stats->af_ts       = (int*) calloc(args->m_af,sizeof(int));
        stats->af_tv       = (int*) calloc(args->m_af,sizeof(int));
        stats->af_snps     = (int*) calloc(args->m_af,sizeof(int));
        int j;
        for (j=0; j<3; j++) stats->af_repeats[j] = (int*) calloc(args->m_af,sizeof(int));
        #if QUAL_STATS
            stats->qual_ts     = dist_init(5);
            stats->qual_tv     = dist_init(5);
            stats->qual_indels = dist_init(5);
        #endif
        if ( args->files->n_smpl )
        {
            stats->smpl_missing = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_hets   = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_homAA  = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_homRR  = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_hapRef = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_hapAlt = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_ins_hets = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_del_hets = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_ins_homs = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_del_homs = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_ts     = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_tv     = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_indels = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_dp     = (unsigned long int *) calloc(args->files->n_smpl,sizeof(unsigned long int));
            stats->smpl_ndp    = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_sngl   = (int *) calloc(args->files->n_smpl,sizeof(int));
            #if HWE_STATS
                stats->af_hwe  = (int*) calloc(args->m_af*args->naf_hwe,sizeof(int));
            #endif
            if ( args->exons_fname )
                stats->smpl_frm_shifts = (int*) calloc(args->files->n_smpl*3,sizeof(int));
            stats->nvaf = (uint32_t*) calloc(stats->m_indel*2+1,sizeof(*stats->nvaf));
            stats->dvaf = (double*) calloc(stats->m_indel*2+1,sizeof(*stats->dvaf));
        }
        idist_init(&stats->dp, args->dp_min,args->dp_max,args->dp_step);
        idist_init(&stats->dp_sites, args->dp_min,args->dp_max,args->dp_step);
        init_user_stats(args, i!=1 ? args->files->readers[0].header : args->files->readers[1].header, stats);
    }

    if ( args->exons_fname )
    {
        args->exons = bcf_sr_regions_init(args->exons_fname,1,0,1,2);
        if ( !args->exons )
            error("Error occurred while reading, was the file compressed with bgzip: %s?\n", args->exons_fname);
    }

    #if IRC_STATS
    if ( args->ref_fname )
        args->indel_ctx = indel_ctx_init(args->ref_fname);
    #endif

    type2dosage[GT_HOM_RR] = 0;
    type2dosage[GT_HET_RA] = 1;
    type2dosage[GT_HOM_AA] = 2;
    type2dosage[GT_HET_AA] = 2;
    type2dosage[GT_HAPL_R] = 0;
    type2dosage[GT_HAPL_A] = 1;

    type2ploidy[GT_HOM_RR] = 1;
    type2ploidy[GT_HET_RA] = 1;
    type2ploidy[GT_HOM_AA] = 1;
    type2ploidy[GT_HET_AA] = 1;
    type2ploidy[GT_HAPL_R] = -1;
    type2ploidy[GT_HAPL_A] = -1;

    type2stats[GT_HOM_RR] = 0;
    type2stats[GT_HET_RA] = 1;
    type2stats[GT_HOM_AA] = 2;
    type2stats[GT_HET_AA] = 3;
    type2stats[GT_HAPL_R] = 0;
    type2stats[GT_HAPL_A] = 2;
    type2stats[GT_UNKN]   = 4;

}
static void destroy_stats(args_t *args)
{
    int id, j;
    for (id=0; id<args->nstats; id++)
    {
        stats_t *stats = &args->stats[id];
        if (stats->af_ts) free(stats->af_ts);
        if (stats->af_tv) free(stats->af_tv);
        if (stats->af_snps) free(stats->af_snps);
        for (j=0; j<3; j++)
            if (stats->af_repeats[j]) free(stats->af_repeats[j]);
        #if QUAL_STATS
            if (stats->qual_ts) dist_destroy(stats->qual_ts);
            if (stats->qual_tv) dist_destroy(stats->qual_tv);
            if (stats->qual_indels) dist_destroy(stats->qual_indels);
        #endif
        #if HWE_STATS
            free(stats->af_hwe);
        #endif
        free(stats->insertions);
        free(stats->deletions);
        free(stats->smpl_missing);
        free(stats->smpl_hets);
        free(stats->smpl_homAA);
        free(stats->smpl_homRR);
        free(stats->smpl_hapRef);
        free(stats->smpl_hapAlt);
        free(stats->smpl_ins_homs);
        free(stats->smpl_del_homs);
        free(stats->smpl_ins_hets);
        free(stats->smpl_del_hets);
        free(stats->smpl_ts);
        free(stats->smpl_tv);
        free(stats->smpl_indels);
        free(stats->smpl_dp);
        free(stats->smpl_ndp);
        free(stats->smpl_sngl);
        idist_destroy(&stats->dp);
        idist_destroy(&stats->dp_sites);
        for (j=0; j<stats->nusr; j++)
        {
            free(stats->usr[j].vals_ts);
            free(stats->usr[j].vals_tv);
            free(stats->usr[j].val);
        }
        free(stats->usr);
        if ( args->exons ) free(stats->smpl_frm_shifts);
        free(stats->nvaf);
        free(stats->dvaf);
    }
    for (j=0; j<args->nusr; j++) free(args->usr[j].tag);
    if ( args->af_bins ) bin_destroy(args->af_bins);
    free(args->farr);
    free(args->usr);
    free(args->tmp_frm);
    free(args->tmp_iaf);
    if (args->exons) bcf_sr_regions_destroy(args->exons);
    free(args->af_gts_snps);
    free(args->af_gts_indels);
    free(args->smpl_gts_snps);
    free(args->smpl_gts_indels);
    if (args->indel_ctx) indel_ctx_destroy(args->indel_ctx);
    if (args->filter[0]) filter_destroy(args->filter[0]);
    if (args->filter[1]) filter_destroy(args->filter[1]);
}

static void init_iaf(args_t *args, bcf_sr_t *reader)
{
    bcf1_t *line = reader->buffer[0];
    hts_expand(int32_t,line->n_allele,args->ntmp_iaf,args->tmp_iaf);

    int i, ret;
    if ( args->af_tag )
    {
        ret = bcf_get_info_float(reader->header, line, args->af_tag, &args->farr, &args->mfarr);
        if ( ret<=0 || ret!=line->n_allele-1 )
        {
            // the AF tag is not present or wrong number of values, put in the singletons/unknown bin
            for (i=0; i<line->n_allele; i++) args->tmp_iaf[i] = 0;
            return;
        }
        args->tmp_iaf[0] = 0;
        for (i=1; i<line->n_allele; i++)
        {
            float af = args->farr[i-1];
            if ( af<0 ) af = 0;
            else if ( af>1 ) af = 1;
            int iaf = args->af_bins ? bin_get_idx(args->af_bins,af) : af*(args->m_af-2);
            args->tmp_iaf[i] = iaf + 1;     // the first tmp_iaf bin is reserved for singletons
        }
        return;
    }

    // tmp_iaf is first filled with AC counts in calc_ac and then transformed to
    //  an index to af_gts_snps
    ret = bcf_calc_ac(reader->header, line, args->tmp_iaf, args->samples_list ? BCF_UN_INFO|BCF_UN_FMT : BCF_UN_INFO);
    if ( !ret )
    {
        for (i=0; i<line->n_allele; i++) args->tmp_iaf[i] = 0;      // singletons/unknown bin
        return;
    }

    int an = 0;
    for (i=0; i<line->n_allele; i++)
        an += args->tmp_iaf[i];

    args->tmp_iaf[0] = 0;
    for (i=1; i<line->n_allele; i++)
    {
        if ( args->tmp_iaf[i]==1 )
            args->tmp_iaf[i] = 0;   // singletons into the first bin
        else if ( !an )
            args->tmp_iaf[i] = 1;   // no genotype at all, put to the AF=0 bin
        else
        {
            float af = (float) args->tmp_iaf[i] / an;
            if ( af<0 ) af = 0;
            else if ( af>1 ) af = 1;
            int iaf = args->af_bins ? bin_get_idx(args->af_bins,af) : af*(args->m_af-2);
            args->tmp_iaf[i] = iaf + 1;
        }
    }
}

static inline void do_mnp_stats(args_t *args, stats_t *stats, bcf_sr_t *reader)
{
    stats->n_mnps++;
}

static inline void do_other_stats(args_t *args, stats_t *stats, bcf_sr_t *reader)
{
    stats->n_others++;
}

static void do_indel_stats(args_t *args, stats_t *stats, bcf_sr_t *reader)
{
    stats->n_indels++;

    bcf1_t *line = reader->buffer[0];

    #if QUAL_STATS
        int iqual = (isnan(line->qual) || line->qual<0) ? 0 : 1 + (int)(line->qual*10);
        dist_insert(stats->qual_indels, iqual);
    #endif

    // Check if the indel is near an exon for the frameshift statistics
    int i, exon_overlap = 0;
    if ( args->exons )
    {
        if ( !bcf_sr_regions_overlap(args->exons, bcf_seqname(reader->header,line),line->pos,line->pos) ) exon_overlap = 1;
        hts_expand(uint8_t,line->n_allele,args->mtmp_frm,args->tmp_frm);
        for (i=0; i<line->n_allele; i++) args->tmp_frm[i] = 0;
    }

    for (i=1; i<line->n_allele; i++)
    {
        if ( args->first_allele_only && i>1 ) break;
        int is_indel = bcf_has_variant_type(line,i,VCF_INDEL);
        if (is_indel < 0) error("bcf_has_variant_type() failed.");
        if ( !is_indel ) continue;
        int len = bcf_variant_length(line, i);

        #if IRC_STATS
        // Indel repeat consistency
        if ( args->indel_ctx )
        {
            int nrep, nlen, ndel;
            ndel = indel_ctx_type(args->indel_ctx, (char*)reader->header->id[BCF_DT_CTG][line->rid].key, line->pos+1, line->d.allele[0], line->d.allele[i], &nrep, &nlen);
            if ( nlen<=1 || nrep<=1 )
            {
                // not a repeat or a single base repeat
                stats->n_repeat_na++;
                stats->af_repeats[2][ args->tmp_iaf[i] ]++;
            }
            else
            {
                if ( abs(ndel) % nlen )
                {
                    // the length of the inserted/deleted sequence is not consistent with the repeat element
                    stats->n_repeat[nlen-1][ndel<0 ? 1 : 3]++;
                    stats->af_repeats[1][ args->tmp_iaf[i] ]++;
                }
                else
                {
                    // the length consistent with the repeat
                    stats->n_repeat[nlen-1][ndel<0 ? 0 : 2]++;
                    stats->af_repeats[0][ args->tmp_iaf[i] ]++;
                }
            }
        }
        else
            stats->af_repeats[2][ args->tmp_iaf[i] ]++;
        #endif

        // Check the frameshifts
        int tlen = 0;
        if ( args->exons && exon_overlap )   // there is an exon
        {
            if ( len>0 )
            {
                // insertion
                if ( args->exons->start <= line->pos && args->exons->end > line->pos ) tlen = abs(len);
            }
            else if ( args->exons->start <= line->pos + abs(len) )
            {
                // deletion
                tlen = abs(len);
                if ( line->pos < args->exons->start )              // trim the beginning
                    tlen -= args->exons->start - line->pos + 1;
                if ( args->exons->end < line->pos + abs(len) )     // trim the end
                    tlen -= line->pos + abs(len) - args->exons->end;
            }
        }
        if ( tlen )     // there are some deleted/inserted bases in the exon
        {
            if ( tlen%3 ) { stats->out_frame++; args->tmp_frm[i] = 2; }
            else { stats->in_frame++; args->tmp_frm[i] = 1; }

            if ( i==1 )
            {
                if ( tlen%3 ) stats->out_frame_alt1++;
                else stats->in_frame_alt1++;
            }
        }
        else            // no exon affected
        {
            if ( i==1 ) stats->na_frame_alt1++;
            stats->na_frame++;
        }


        // Indel length distribution
        int *ptr = stats->insertions;
        if ( len<0 )
        {
            len *= -1;
            ptr = stats->deletions;
        }
        if ( --len >= stats->m_indel ) len = stats->m_indel-1;
        ptr[len]++;
    }
}

static void do_user_stats(stats_t *stats, bcf_sr_t *reader, int is_ts)
{
    int i, nval;
    for (i=0; i<stats->nusr; i++)
    {
        user_stats_t *usr = &stats->usr[i];
        uint64_t *vals = is_ts ? usr->vals_ts : usr->vals_tv;
        float val;
        if ( usr->type==BCF_HT_REAL )
        {
            if ( (nval=bcf_get_info_float(reader->header,reader->buffer[0],usr->tag,&usr->val,&usr->m_val))<=0 ) continue;
            if ( usr->idx >= nval ) continue;
            val = ((float*)usr->val)[usr->idx];
        }
        else
        {
            if ( (nval=bcf_get_info_int32(reader->header,reader->buffer[0],usr->tag,&usr->val,&usr->m_val))<=0 ) continue;
            if ( usr->idx >= nval ) continue;
            val = ((int32_t*)usr->val)[usr->idx];
        }
        int idx;
        if ( val<=usr->min ) idx = 0;
        else if ( val>=usr->max ) idx = usr->nbins - 1;
        else idx = (val - usr->min)/(usr->max - usr->min) * (usr->nbins-1);
        vals[idx]++;
    }
}

static void do_snp_stats(args_t *args, stats_t *stats, bcf_sr_t *reader)
{
    stats->n_snps++;

    bcf1_t *line = reader->buffer[0];
    int ref = bcf_acgt2int(*line->d.allele[0]);
    if ( ref<0 ) return;

    #if QUAL_STATS
        int iqual = (isnan(line->qual) || line->qual<0) ? 0 : 1 + (int)(line->qual*10);
    #endif

    int i;
    for (i=1; i<line->n_allele; i++)
    {
        if ( args->first_allele_only && i>1 ) break;
        if ( !(bcf_get_variant_type(line,i)&VCF_SNP) ) continue;
        int alt = bcf_acgt2int(*line->d.allele[i]);
        if ( alt<0 || ref==alt ) continue;
        stats->subst[ref<<2|alt]++;
        int iaf = args->tmp_iaf[i];
        stats->af_snps[iaf]++;
        if ( abs(ref-alt)==2 )
        {
            if (i==1)
            {
                stats->ts_alt1++;
                #if QUAL_STATS
                    dist_insert(stats->qual_ts,iqual);
                #endif
                do_user_stats(stats, reader, 1);
            }
            stats->af_ts[iaf]++;
        }
        else
        {
            if (i==1)
            {
                stats->tv_alt1++;
                #if QUAL_STATS
                    dist_insert(stats->qual_tv,iqual);
                #endif
                do_user_stats(stats, reader, 0);
            }
            stats->af_tv[iaf]++;
        }
    }
}

static inline void update_dvaf(stats_t *stats, bcf1_t *line, bcf_fmt_t *fmt, int ismpl, int ial, int jal)
{
    if ( !fmt ) return;

    float dvaf;
    #define BRANCH_INT(type_t,missing,vector_end) { \
        type_t *p = (type_t *) (fmt->p + fmt->size*ismpl); \
        if ( p[ial]==vector_end || p[jal]==vector_end ) return; \
        if ( p[ial]==missing || p[jal]==missing ) return; \
        if ( !p[ial] && !p[jal] ) return; \
        dvaf = (float)p[ial]/(p[ial]+p[jal]); \
    }
    switch (fmt->type) {
        case BCF_BT_INT8:  BRANCH_INT(int8_t,  bcf_int8_missing, bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH_INT(int16_t, bcf_int16_missing, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH_INT(int32_t, bcf_int32_missing, bcf_int32_vector_end); break;
        default: fprintf(stderr, "[E::%s] todo: %d\n", __func__, fmt->type); exit(1); break;
    }
    #undef BRANCH_INT

    int len = line->d.var[ial].n;
    if ( len < -stats->m_indel ) len = -stats->m_indel;
    else if ( len > stats->m_indel ) len = stats->m_indel;
    int bin = stats->m_indel + len;
    stats->nvaf[bin]++;
    stats->dvaf[bin] += dvaf;
}

static void do_sample_stats(args_t *args, stats_t *stats, bcf_sr_t *reader, int matched)
{
    bcf_srs_t *files = args->files;
    bcf1_t *line = reader->buffer[0];
    bcf_fmt_t *fmt_ptr;
    int nref_tot = 0, nhet_tot = 0, nalt_tot = 0;
    int line_type = bcf_get_variant_types(line);

    if ( (fmt_ptr = bcf_get_fmt(reader->header,reader->buffer[0],"GT")) )
    {
        bcf_fmt_t *ad_fmt_ptr = bcf_get_variant_types(line)&VCF_INDEL ? bcf_get_fmt(reader->header,reader->buffer[0],"AD") : NULL;

        int ref = bcf_acgt2int(*line->d.allele[0]);
        int is, n_nref = 0, i_nref = 0;
        for (is=0; is<args->files->n_smpl; is++)
        {
            int ial, jal;
            int gt = bcf_gt_type(fmt_ptr, reader->samples[is], &ial, &jal);
            if ( gt==GT_UNKN )
            {
                stats->smpl_missing[is]++;
                continue;
            }
            if ( gt==GT_HAPL_R || gt==GT_HAPL_A )
            {
                if ( line_type&VCF_INDEL && stats->smpl_frm_shifts )
                {
                    assert( ial<line->n_allele );
                    stats->smpl_frm_shifts[is*3 + args->tmp_frm[ial]]++;
                }
                if ( gt == GT_HAPL_R ) stats->smpl_hapRef[is]++;
                if ( gt == GT_HAPL_A ) stats->smpl_hapAlt[is]++;
                continue;
            }
            if ( gt != GT_HOM_RR ) { n_nref++; i_nref = is; }
            #if HWE_STATS
                switch (gt)
                {
                    case GT_HOM_RR: nref_tot++; break;
                    case GT_HET_RA: nhet_tot++; break;
                    case GT_HET_AA:
                    case GT_HOM_AA: nalt_tot++; break;
                }
            #endif
            int var_type = 0;
            if ( ial>0 ) var_type |= bcf_get_variant_type(line,ial);
            if ( jal>0 ) var_type |= bcf_get_variant_type(line,jal);
            if ( var_type&VCF_SNP || var_type==VCF_REF )  // count ALT=. as SNP
            {
                if ( gt == GT_HET_RA ) stats->smpl_hets[is]++;
                else if ( gt == GT_HET_AA ) stats->smpl_hets[is]++;
                else if ( gt == GT_HOM_RR ) stats->smpl_homRR[is]++;
                else if ( gt == GT_HOM_AA ) stats->smpl_homAA[is]++;
                if ( gt != GT_HOM_RR && line->d.var[ial].type&VCF_SNP ) // this is safe, bcf_get_variant_types has been already called
                {
                    int alt = bcf_acgt2int(*line->d.allele[ial]);
                    if ( alt<0 ) continue;
                    if ( abs(ref-alt)==2 )
                        stats->smpl_ts[is]++;
                    else
                        stats->smpl_tv[is]++;
                }
            }
            if ( var_type&VCF_INDEL )
            {
                if ( gt != GT_HOM_RR )
                {
                    stats->smpl_indels[is]++;

                    if ( gt==GT_HET_RA || gt==GT_HET_AA )
                    {
                        int is_ins = 0, is_del = 0;
                        if ( bcf_get_variant_type(line,ial)&VCF_INDEL )
                        {
                            if ( line->d.var[ial].n < 0 ) is_del = 1;
                            else is_ins = 1;
                            update_dvaf(stats,line,ad_fmt_ptr,is,ial,jal);
                        }
                        if ( bcf_get_variant_type(line,jal)&VCF_INDEL )
                        {
                            if ( line->d.var[jal].n < 0 ) is_del = 1;
                            else is_ins = 1;
                            update_dvaf(stats,line,ad_fmt_ptr,is,jal,ial);
                        }
                        // Note that alt-het genotypes with both ins and del allele are counted twice!!
                        if ( is_del ) stats->smpl_del_hets[is]++;
                        if ( is_ins ) stats->smpl_ins_hets[is]++;
                    }
                    else if ( gt==GT_HOM_AA )
                    {
                        if ( line->d.var[ial].n < 0 ) stats->smpl_del_homs[is]++;
                        else stats->smpl_ins_homs[is]++;
                    }
                }
                if ( stats->smpl_frm_shifts )
                {
                    assert( ial<line->n_allele && jal<line->n_allele );
                    stats->smpl_frm_shifts[is*3 + args->tmp_frm[ial]]++;
                    stats->smpl_frm_shifts[is*3 + args->tmp_frm[jal]]++;
                }
            }
        }
        if ( n_nref==1 ) stats->smpl_sngl[i_nref]++;
    }

    #if HWE_STATS
        if ( nhet_tot + nref_tot + nalt_tot )
        {
            float het_frac = (float)nhet_tot/(nhet_tot + nref_tot + nalt_tot);
            int idx = het_frac*(args->naf_hwe - 1);
//check me: what is this?
            if ( line->n_allele>1 ) idx += args->naf_hwe*args->tmp_iaf[1];
            stats->af_hwe[idx]++;
        }
    #endif

    if ( (fmt_ptr = bcf_get_fmt(reader->header,reader->buffer[0],"DP")) )
    {
        #define BRANCH_INT(type_t,missing,vector_end) { \
            int is; \
            for (is=0; is<args->files->n_smpl; is++) \
            { \
                type_t *p = (type_t *) (fmt_ptr->p + fmt_ptr->size*is); \
                if ( *p==vector_end ) continue; \
                if ( *p!=missing ) \
                { \
                    (*idist(&stats->dp, *p))++; \
                    stats->smpl_ndp[is]++; \
                    stats->smpl_dp[is] += *p; \
                } \
            } \
        }
        switch (fmt_ptr->type) {
            case BCF_BT_INT8:  BRANCH_INT(int8_t,  bcf_int8_missing, bcf_int8_vector_end); break;
            case BCF_BT_INT16: BRANCH_INT(int16_t, bcf_int16_missing, bcf_int16_vector_end); break;
            case BCF_BT_INT32: BRANCH_INT(int32_t, bcf_int32_missing, bcf_int32_vector_end); break;
            default: fprintf(stderr, "[E::%s] todo: %d\n", __func__, fmt_ptr->type); exit(1); break;
        }
        #undef BRANCH_INT
    }
    else if ( (fmt_ptr = bcf_get_fmt(reader->header,reader->buffer[0],"AD")) )
    {
        #define BRANCH_INT(type_t,missing,vector_end) { \
            int is,iv; \
            for (is=0; is<args->files->n_smpl; is++) \
            { \
                type_t *p = (type_t *) (fmt_ptr->p + fmt_ptr->size*is); \
                int dp = 0, has_value = 0; \
                for (iv=0; iv<fmt_ptr->n; iv++) \
                { \
                    if ( p[iv]==vector_end ) break; \
                    if ( p[iv]==missing ) continue; \
                    has_value = 1; \
                    dp += p[iv]; \
                } \
                if ( has_value ) \
                { \
                    (*idist(&stats->dp, dp))++; \
                    stats->smpl_ndp[is]++; \
                    stats->smpl_dp[is] += dp; \
                } \
            } \
        }
        switch (fmt_ptr->type) {
            case BCF_BT_INT8:  BRANCH_INT(int8_t,  bcf_int8_missing, bcf_int8_vector_end); break;
            case BCF_BT_INT16: BRANCH_INT(int16_t, bcf_int16_missing, bcf_int16_vector_end); break;
            case BCF_BT_INT32: BRANCH_INT(int32_t, bcf_int32_missing, bcf_int32_vector_end); break;
            default: fprintf(stderr, "[E::%s] todo: %d\n", __func__, fmt_ptr->type); exit(1); break;
        }
        #undef BRANCH_INT
    }

    if ( matched==3 )
    {
        int is;
        bcf_fmt_t *fmt0, *fmt1;
        fmt0 = bcf_get_fmt(files->readers[0].header,files->readers[0].buffer[0],"GT"); if ( !fmt0 ) return;
        fmt1 = bcf_get_fmt(files->readers[1].header,files->readers[1].buffer[0],"GT"); if ( !fmt1 ) return;

        // only the first ALT allele is considered
        if (args->ntmp_iaf <= 1) return; // Do not consider invariate sites
        int iaf = args->tmp_iaf[1];
        int line_type = bcf_get_variant_types(files->readers[0].buffer[0]);
        gtcmp_t *af_stats = line_type&VCF_SNP ? args->af_gts_snps : args->af_gts_indels;
        gtcmp_t *smpl_stats = line_type&VCF_SNP ? args->smpl_gts_snps : args->smpl_gts_indels;

        for (is=0; is<files->n_smpl; is++)
        {
            // Simplified comparison: only 0/0, 0/1, 1/1 is looked at as the identity of
            //  actual alleles can be enforced by running without the -c option.
            int gt0 = bcf_gt_type(fmt0, files->readers[0].samples[is], NULL, NULL);
            int gt1 = bcf_gt_type(fmt1, files->readers[1].samples[is], NULL, NULL);

            int idx0 = type2stats[gt0];
            int idx1 = type2stats[gt1];
            af_stats[iaf].gt2gt[idx0][idx1]++;
            smpl_stats[is].gt2gt[idx0][idx1]++;

            if ( gt0 == GT_UNKN || gt1 == GT_UNKN ) continue;
            if ( type2ploidy[gt0]*type2ploidy[gt1] == -1 ) continue;   // cannot compare diploid and haploid genotypes

            float y = type2dosage[gt0];
            float x = type2dosage[gt1];

            smpl_stats[is].yx += y*x;
            smpl_stats[is].x  += x;
            smpl_stats[is].xx += x*x;
            smpl_stats[is].y  += y;
            smpl_stats[is].yy += y*y;
            smpl_stats[is].n  += 1;

            af_stats[iaf].yx += y*x;
            af_stats[iaf].x  += x;
            af_stats[iaf].xx += x*x;
            af_stats[iaf].y  += y;
            af_stats[iaf].yy += y*y;
            af_stats[iaf].n  += 1;
        }

        if ( args->verbose_sites )
        {
            int nm = 0, nmm = 0, nrefm = 0;
            for (is=0; is<files->n_smpl; is++)
            {
                int gt = bcf_gt_type(fmt0, files->readers[0].samples[is], NULL, NULL);
                if ( gt == GT_UNKN ) continue;
                int gt2 = bcf_gt_type(fmt1, files->readers[1].samples[is], NULL, NULL);
                if ( gt2 == GT_UNKN ) continue;
                if ( gt != gt2 )
                {
                    nmm++;
                    bcf_sr_t *reader = &files->readers[0];
                    printf("DBG\t%s\t%"PRId64"\t%s\t%d\t%d\n",reader->header->id[BCF_DT_CTG][reader->buffer[0]->rid].key,(int64_t) reader->buffer[0]->pos+1,files->samples[is],gt,gt2);
                }
                else
                {
                    if ( gt!=GT_HOM_RR ) nrefm++;
                    nm++;
                }
            }
            float nrd = nrefm+nmm ? 100.*nmm/(nrefm+nmm) : 0;
            printf("PSD\t%s\t%"PRId64"\t%d\t%d\t%f\n", reader->header->id[BCF_DT_CTG][reader->buffer[0]->rid].key,(int64_t) reader->buffer[0]->pos+1,nm,nmm,nrd);
        }
    }
}

static void do_vcf_stats(args_t *args)
{
    bcf_srs_t *files = args->files;
    assert( sizeof(int)>files->nreaders );
    while ( bcf_sr_next_line(files) )
    {
        bcf_sr_t *reader = NULL;
        bcf1_t *line = NULL;
        int ret = 0, i, pass = 1;
        for (i=0; i<files->nreaders; i++)
        {
            if ( !bcf_sr_has_line(files,i) ) continue;
            if ( args->filter[i] )
            {
                int is_ok = filter_test(args->filter[i], bcf_sr_get_line(files,i), NULL);
                if ( args->filter_logic & FLT_EXCLUDE ) is_ok = is_ok ? 0 : 1;
                if ( !is_ok ) { pass = 0; break; }
            }
            ret |= 1<<i;
            if ( !reader )
            {
                reader = &files->readers[i];
                line = bcf_sr_get_line(files,i);
            }

        }
        if ( !pass ) continue;

        int line_type = bcf_get_variant_types(line);
        init_iaf(args, reader);

        stats_t *stats = &args->stats[ret-1];
        if ( args->split_by_id && line->d.id[0]=='.' && !line->d.id[1] )
            stats = &args->stats[1];

        stats->n_records++;

        if ( line_type==VCF_REF )
            stats->n_noalts++;
        if ( line_type&VCF_SNP )
            do_snp_stats(args, stats, reader);
        if ( line_type&VCF_INDEL )
            do_indel_stats(args, stats, reader);
        if ( line_type&VCF_MNP )
            do_mnp_stats(args, stats, reader);
        if ( line_type&VCF_OTHER )
            do_other_stats(args, stats, reader);

        if ( line->n_allele>2 )
        {
            stats->n_mals++;
            if ( line_type == VCF_SNP ) stats->n_snp_mals++;    // note: this will be fooled by C>C,T
        }

        if ( files->n_smpl )
            do_sample_stats(args, stats, reader, ret);

        if ( bcf_get_info_int32(reader->header,line,"DP",&args->tmp_iaf,&args->ntmp_iaf)==1 )
            (*idist(&stats->dp_sites, args->tmp_iaf[0]))++;
    }
}

static void print_header(args_t *args)
{
    int i;
    printf("# This file was produced by bcftools stats (%s+htslib-%s) and can be plotted using plot-vcfstats.\n", bcftools_version(),hts_version());
    printf("# The command line was:\tbcftools %s ", args->argv[0]);
    for (i=1; i<args->argc; i++)
        printf(" %s",args->argv[i]);
    printf("\n#\n");

    printf("# Definition of sets:\n# ID\t[2]id\t[3]tab-separated file names\n");
    if ( args->files->nreaders==1 )
    {
        const char *fname = strcmp("-",args->files->readers[0].fname) ? args->files->readers[0].fname : "<STDIN>";
        if ( args->split_by_id )
        {
            printf("ID\t0\t%s:known (sites with ID different from \".\")\n", fname);
            printf("ID\t1\t%s:novel (sites where ID column is \".\")\n", fname);
        }
        else
            printf("ID\t0\t%s\n", fname);
    }
    else
    {
        const char *fname0 = strcmp("-",args->files->readers[0].fname) ? args->files->readers[0].fname : "<STDIN>";
        const char *fname1 = strcmp("-",args->files->readers[1].fname) ? args->files->readers[1].fname : "<STDIN>";
        printf("ID\t0\t%s\n", fname0);
        printf("ID\t1\t%s\n", fname1);
        printf("ID\t2\t%s\t%s\n", fname0,fname1);

        if ( args->verbose_sites )
        {
            printf(
                    "# Verbose per-site discordance output.\n"
                    "# PSD\t[2]CHROM\t[3]POS\t[4]Number of matches\t[5]Number of mismatches\t[6]NRD\n");
            printf(
                    "# Verbose per-site and per-sample output. Genotype codes: %d:HomRefRef, %d:HomAltAlt, %d:HetAltRef, %d:HetAltAlt, %d:haploidRef, %d:haploidAlt\n"
                    "# DBG\t[2]CHROM\t[3]POS\t[4]Sample\t[5]GT in %s\t[6]GT in %s\n",
                    GT_HOM_RR, GT_HOM_AA, GT_HET_RA, GT_HET_AA, GT_HAPL_R, GT_HAPL_A, fname0,fname1);
        }
    }
}

#define T2S(x) type2stats[x]
static void print_stats(args_t *args)
{
    int i, j,k, id;
    printf("# SN, Summary numbers:\n");
    printf("#   number of records   .. number of data rows in the VCF\n");
    printf("#   number of no-ALTs   .. reference-only sites, ALT is either \".\" or identical to REF\n");
    printf("#   number of SNPs      .. number of rows with a SNP\n");
    printf("#   number of MNPs      .. number of rows with a MNP, such as CC>TT\n");
    printf("#   number of indels    .. number of rows with an indel\n");
    printf("#   number of others    .. number of rows with other type, for example a symbolic allele or\n");
    printf("#                          a complex substitution, such as ACT>TCGA\n");
    printf("#   number of multiallelic sites     .. number of rows with multiple alternate alleles\n");
    printf("#   number of multiallelic SNP sites .. number of rows with multiple alternate alleles, all SNPs\n");
    printf("# \n");
    printf("#   Note that rows containing multiple types will be counted multiple times, in each\n");
    printf("#   counter. For example, a row with a SNP and an indel increments both the SNP and\n");
    printf("#   the indel counter.\n");
    printf("# \n");
    printf("# SN\t[2]id\t[3]key\t[4]value\n");
    for (id=0; id<args->files->nreaders; id++)
        printf("SN\t%d\tnumber of samples:\t%d\n", id, bcf_hdr_nsamples(args->files->readers[id].header));
    for (id=0; id<args->nstats; id++)
    {
        stats_t *stats = &args->stats[id];
        printf("SN\t%d\tnumber of records:\t%"PRIu64"\n", id, stats->n_records);
        printf("SN\t%d\tnumber of no-ALTs:\t%"PRIu64"\n", id, stats->n_noalts);
        printf("SN\t%d\tnumber of SNPs:\t%"PRIu64"\n", id, stats->n_snps);
        printf("SN\t%d\tnumber of MNPs:\t%"PRIu64"\n", id, stats->n_mnps);
        printf("SN\t%d\tnumber of indels:\t%"PRIu64"\n", id, stats->n_indels);
        printf("SN\t%d\tnumber of others:\t%"PRIu64"\n", id, stats->n_others);
        printf("SN\t%d\tnumber of multiallelic sites:\t%"PRIu64"\n", id, stats->n_mals);
        printf("SN\t%d\tnumber of multiallelic SNP sites:\t%"PRIu64"\n", id, stats->n_snp_mals);
    }
    printf("# TSTV, transitions/transversions:\n# TSTV\t[2]id\t[3]ts\t[4]tv\t[5]ts/tv\t[6]ts (1st ALT)\t[7]tv (1st ALT)\t[8]ts/tv (1st ALT)\n");
    for (id=0; id<args->nstats; id++)
    {
        stats_t *stats = &args->stats[id];
        int ts=0,tv=0;
        for (i=0; i<args->m_af; i++) { ts += stats->af_ts[i]; tv += stats->af_tv[i];  }
        printf("TSTV\t%d\t%d\t%d\t%.2f\t%d\t%d\t%.2f\n", id,ts,tv,tv?(float)ts/tv:0, stats->ts_alt1,stats->tv_alt1,stats->tv_alt1?(float)stats->ts_alt1/stats->tv_alt1:0);
    }
    if ( args->exons_fname )
    {
        printf("# FS, Indel frameshifts:\n# FS\t[2]id\t[3]in-frame\t[4]out-frame\t[5]not applicable\t[6]out/(in+out) ratio\t[7]in-frame (1st ALT)\t[8]out-frame (1st ALT)\t[9]not applicable (1st ALT)\t[10]out/(in+out) ratio (1st ALT)\n");
        for (id=0; id<args->nstats; id++)
        {
            int in=args->stats[id].in_frame, out=args->stats[id].out_frame, na=args->stats[id].na_frame;
            int in1=args->stats[id].in_frame_alt1, out1=args->stats[id].out_frame_alt1, na1=args->stats[id].na_frame_alt1;
            printf("FS\t%d\t%d\t%d\t%d\t%.2f\t%d\t%d\t%d\t%.2f\n", id, in,out,na,out?(float)out/(in+out):0,in1,out1,na1,out1?(float)out1/(in1+out1):0);
        }
    }
    if ( args->indel_ctx )
    {
        printf("# ICS, Indel context summary:\n# ICS\t[2]id\t[3]repeat-consistent\t[4]repeat-inconsistent\t[5]not applicable\t[6]c/(c+i) ratio\n");
        for (id=0; id<args->nstats; id++)
        {
            int nc = 0, ni = 0, na = args->stats[id].n_repeat_na;
            for (i=0; i<IRC_RLEN; i++)
            {
                nc += args->stats[id].n_repeat[i][0] + args->stats[id].n_repeat[i][2];
                ni += args->stats[id].n_repeat[i][1] + args->stats[id].n_repeat[i][3];
            }
            printf("ICS\t%d\t%d\t%d\t%d\t%.4f\n", id, nc,ni,na,nc+ni ? (float)nc/(nc+ni) : 0.0);
        }
        printf("# ICL, Indel context by length:\n# ICL\t[2]id\t[3]length of repeat element\t[4]repeat-consistent deletions)\t[5]repeat-inconsistent deletions\t[6]consistent insertions\t[7]inconsistent insertions\t[8]c/(c+i) ratio\n");
        for (id=0; id<args->nstats; id++)
        {
            for (i=1; i<IRC_RLEN; i++)
            {
                int nc = args->stats[id].n_repeat[i][0]+args->stats[id].n_repeat[i][2], ni = args->stats[id].n_repeat[i][1]+args->stats[id].n_repeat[i][3];
                printf("ICL\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\n", id, i+1,
                    args->stats[id].n_repeat[i][0],args->stats[id].n_repeat[i][1],args->stats[id].n_repeat[i][2],args->stats[id].n_repeat[i][3],
                    nc+ni ? (float)nc/(nc+ni) : 0.0);
            }
        }
    }
    printf("# SiS, Singleton stats:\n# SiS\t[2]id\t[3]allele count\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable\n");
    for (id=0; id<args->nstats; id++)
    {
        stats_t *stats = &args->stats[id];
        printf("SiS\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", id,1,stats->af_snps[0],stats->af_ts[0],stats->af_tv[0],
            stats->af_repeats[0][0]+stats->af_repeats[1][0]+stats->af_repeats[2][0],stats->af_repeats[0][0],stats->af_repeats[1][0],stats->af_repeats[2][0]);
        // put the singletons stats into the first AF bin, note that not all of the stats is transferred (i.e. nrd mismatches)
        stats->af_snps[1]       += stats->af_snps[0];
        stats->af_ts[1]         += stats->af_ts[0];
        stats->af_tv[1]         += stats->af_tv[0];
        stats->af_repeats[0][1] += stats->af_repeats[0][0];
        stats->af_repeats[1][1] += stats->af_repeats[1][0];
        stats->af_repeats[2][1] += stats->af_repeats[2][0];
    }
    // move the singletons stats into the first AF bin, singleton stats was collected separately because of init_iaf
    if ( args->af_gts_snps )
    {
        args->af_gts_snps[1].y    += args->af_gts_snps[0].y;
        args->af_gts_snps[1].yy   += args->af_gts_snps[0].yy;
        args->af_gts_snps[1].xx   += args->af_gts_snps[0].xx;
        args->af_gts_snps[1].yx   += args->af_gts_snps[0].yx;
        args->af_gts_snps[1].n    += args->af_gts_snps[0].n;
    }
    if ( args->af_gts_indels )
    {
        args->af_gts_indels[1].y  += args->af_gts_indels[0].y;
        args->af_gts_indels[1].yy += args->af_gts_indels[0].yy;
        args->af_gts_indels[1].xx += args->af_gts_indels[0].xx;
        args->af_gts_indels[1].yx += args->af_gts_indels[0].yx;
        args->af_gts_indels[1].n  += args->af_gts_indels[0].n;
    }

    printf("# AF, Stats by non-reference allele frequency:\n# AF\t[2]id\t[3]allele frequency\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable\n");
    for (id=0; id<args->nstats; id++)
    {
        stats_t *stats = &args->stats[id];
        for (i=1; i<args->m_af; i++) // note that af[1] now contains also af[0], see SiS stats output above
        {
            if ( stats->af_snps[i]+stats->af_ts[i]+stats->af_tv[i]+stats->af_repeats[0][i]+stats->af_repeats[1][i]+stats->af_repeats[2][i] == 0  ) continue;
            double af = args->af_bins ? (bin_get_value(args->af_bins,i)+bin_get_value(args->af_bins,i-1))*0.5 : (double)(i-1)/(args->m_af-1);
            printf("AF\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", id,af,stats->af_snps[i],stats->af_ts[i],stats->af_tv[i],
                stats->af_repeats[0][i]+stats->af_repeats[1][i]+stats->af_repeats[2][i],stats->af_repeats[0][i],stats->af_repeats[1][i],stats->af_repeats[2][i]);
        }
    }
    #if QUAL_STATS
        printf("# QUAL, Stats by quality\n# QUAL\t[2]id\t[3]Quality\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t[6]number of transversions (1st ALT)\t[7]number of indels\n");
        for (id=0; id<args->nstats; id++)
        {
            stats_t *stats = &args->stats[id];
            int ndist_ts = dist_nbins(stats->qual_ts);
            int ndist_tv = dist_nbins(stats->qual_tv);
            int ndist_in = dist_nbins(stats->qual_indels);
            int ndist_max = ndist_ts;
            if ( ndist_max < ndist_tv ) ndist_max = ndist_tv;
            if ( ndist_max < ndist_in ) ndist_max = ndist_in;
            uint32_t beg, end;
            uint32_t nts, ntv, nin;
            for (i=0; i<ndist_max; i++)
            {
                nts = ntv = nin = 0;
                float qval = -1;
                if ( i < ndist_ts )
                {
                    nts = dist_get(stats->qual_ts, i, &beg, &end);
                    qval = beg>0 ? 0.1*(beg - 1) : -1;
                }
                if ( i < ndist_tv )
                {
                    ntv = dist_get(stats->qual_tv, i, &beg, &end);
                    if ( qval==-1 ) qval = beg > 0 ? 0.1*(beg - 1) : -1;
                }
                if ( i < ndist_in )
                {
                    nin = dist_get(stats->qual_indels, i, &beg, &end);
                    if ( qval==-1 ) qval = beg > 0 ? 0.1*(beg - 1) : -1;
                }
                if ( nts+ntv+nin==0 ) continue;

                printf("QUAL\t%d\t",id);
                if ( qval==-1 ) printf(".");
                else printf("%.1f",qval);
                printf("\t%d\t%d\t%d\t%d\n",nts+ntv,nts,ntv,nin);
            }
        }
    #endif
    for (i=0; i<args->nusr; i++)
    {
        printf("# USR:%s/%d\t[2]id\t[3]%s/%d\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t[6]number of transversions (1st ALT)\n",
            args->usr[i].tag,args->usr[i].idx,args->usr[i].tag,args->usr[i].idx);
        for (id=0; id<args->nstats; id++)
        {
            user_stats_t *usr = &args->stats[id].usr[i];
            int j;
            for (j=0; j<usr->nbins; j++)
            {
                if ( usr->vals_ts[j]+usr->vals_tv[j] == 0 ) continue;   // skip empty bins
                float val = usr->min + (usr->max - usr->min)*j/(usr->nbins-1);
                const char *fmt = usr->type==BCF_HT_REAL ? "USR:%s/%d\t%d\t%e\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n" : "USR:%s/%d\t%d\t%.0f\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n";
                printf(fmt,usr->tag,usr->idx,id,val,usr->vals_ts[j]+usr->vals_tv[j],usr->vals_ts[j],usr->vals_tv[j]);
            }
        }
    }
    printf("# IDD, InDel distribution:\n# IDD\t[2]id\t[3]length (deletions negative)\t[4]number of sites\t[5]number of genotypes\t[6]mean VAF\n");
    for (id=0; id<args->nstats; id++)
    {
        stats_t *stats = &args->stats[id];
        for (i=stats->m_indel-1; i>=0; i--)
        {
            if ( !stats->deletions[i] ) continue;
            // whops, differently organized arrow, dels are together with ins
            int bin = stats->m_indel - i - 1;
            printf("IDD\t%d\t%d\t%d\t", id,-i-1,stats->deletions[i]);
            if ( stats->nvaf && stats->nvaf[bin] )
                printf("%u\t%.2f",stats->nvaf[bin],stats->dvaf[bin]/stats->nvaf[bin]);
            else
                printf("0\t.");
            printf("\n");
        }
        for (i=0; i<stats->m_indel; i++)
        {
            if ( !stats->insertions[i] ) continue;
            int bin = stats->m_indel + i + 1;
            printf("IDD\t%d\t%d\t%d\t", id,i+1,stats->insertions[i]);
            if ( stats->nvaf && stats->nvaf[bin] )
                printf("%u\t%.2f",stats->nvaf[bin],stats->dvaf[bin]/stats->nvaf[bin]);
            else
                printf("0\t.");
            printf("\n");
        }
    }
    printf("# ST, Substitution types:\n# ST\t[2]id\t[3]type\t[4]count\n");
    for (id=0; id<args->nstats; id++)
    {
        int t;
        for (t=0; t<15; t++)
        {
            if ( t>>2 == (t&3) ) continue;
            printf("ST\t%d\t%c>%c\t%d\n", id, bcf_int2acgt(t>>2),bcf_int2acgt(t&3),args->stats[id].subst[t]);
        }
    }
    if ( args->files->nreaders>1 && args->files->n_smpl )
    {
        printf("SN\t%d\tnumber of samples:\t%d\n", 2, args->files->n_smpl);

        int x;
        for (x=0; x<2; x++)     // x=0: snps, x=1: indels
        {
            gtcmp_t *stats;
            if ( x==0 )
            {
                printf("# GCsAF, Genotype concordance by non-reference allele frequency (SNPs)\n# GCsAF\t[2]id\t[3]allele frequency\t[4]RR Hom matches\t[5]RA Het matches\t[6]AA Hom matches\t[7]RR Hom mismatches\t[8]RA Het mismatches\t[9]AA Hom mismatches\t[10]dosage r-squared\t[11]number of genotypes\n");
                stats = args->af_gts_snps;
            }
            else
            {
                printf("# GCiAF, Genotype concordance by non-reference allele frequency (indels)\n# GCiAF\t[2]id\t[3]allele frequency\t[4]RR Hom matches\t[5]RA Het matches\t[6]AA Hom matches\t[7]RR Hom mismatches\t[8]RA Het mismatches\t[9]AA Hom mismatches\t[10]dosage r-squared\t[11]number of genotypes\n");
                stats = args->af_gts_indels;
            }
            uint64_t nrd_m[4] = {0,0,0,0}, nrd_mm[4] = {0,0,0,0};   // across all bins
            for (i=0; i<args->m_af; i++)
            {
                int n = 0;
                uint64_t m[4] = {0,0,0,0}, mm[4] = {0,0,0,0};    // in i-th AF bin
                for (j=0; j<4; j++)     // rr, ra, aa hom, aa het, ./.
                    for (k=0; k<4; k++)
                    {
                        n += stats[i].gt2gt[j][k];
                        if ( j==k )
                        {
                            nrd_m[j] += stats[i].gt2gt[j][k];
                            m[j]     += stats[i].gt2gt[j][k];
                        }
                        else
                        {
                            nrd_mm[j] += stats[i].gt2gt[j][k];
                            mm[j]     += stats[i].gt2gt[j][k];
                        }
                    }
                if ( !i || !n ) continue;   // skip singleton stats and empty bins

                // Pearson's r2
                double r2 = 0;
                if ( stats[i].n )
                {
                    r2  = (stats[i].yx - stats[i].x*stats[i].y/stats[i].n);
                    r2 /= sqrt((stats[i].xx - stats[i].x*stats[i].x/stats[i].n) * (stats[i].yy - stats[i].y*stats[i].y/stats[i].n));
                    r2 *= r2;
                }
                double af = args->af_bins ? (bin_get_value(args->af_bins,i)+bin_get_value(args->af_bins,i-1))*0.5 : (double)(i-1)/(args->m_af-1);
                printf("GC%cAF\t2\t%f", x==0 ? 's' : 'i', af);
                printf("\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"", m[T2S(GT_HOM_RR)],m[T2S(GT_HET_RA)],m[T2S(GT_HOM_AA)]);
                printf("\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"", mm[T2S(GT_HOM_RR)],mm[T2S(GT_HET_RA)],mm[T2S(GT_HOM_AA)]);
                if ( stats[i].n && !isnan(r2) ) printf("\t%f", r2);
                else printf("\t"NA_STRING);
                printf("\t%.0f\n", stats[i].n);
            }

            if ( x==0 )
            {
                printf("# NRD and discordance is calculated as follows:\n");
                printf("#   m .. number of matches\n");
                printf("#   x .. number of mismatches\n");
                printf("#   NRD = 100 * (xRR + xRA + xAA) / (xRR + xRA + xAA + mRA + mAA)\n");
                printf("#   RR discordance = 100 * xRR / (xRR + mRR)\n");
                printf("#   RA discordance = 100 * xRA / (xRA + mRA)\n");
                printf("#   AA discordance = 100 * xAA / (xAA + mAA)\n");
                printf("# Non-Reference Discordance (NRD), SNPs\n# NRDs\t[2]id\t[3]NRD\t[4]Ref/Ref discordance\t[5]Ref/Alt discordance\t[6]Alt/Alt discordance\n");
            }
            else
                printf("# Non-Reference Discordance (NRD), indels\n# NRDi\t[2]id\t[3]NRD\t[4]Ref/Ref discordance\t[5]Ref/Alt discordance\t[6]Alt/Alt discordance\n");
            uint64_t m  = nrd_m[T2S(GT_HET_RA)] + nrd_m[T2S(GT_HOM_AA)] + nrd_m[T2S(GT_HET_AA)];
            uint64_t mm = nrd_mm[T2S(GT_HOM_RR)] + nrd_mm[T2S(GT_HET_RA)] + nrd_mm[T2S(GT_HOM_AA)] + nrd_mm[T2S(GT_HET_AA)];
            printf("NRD%c\t2\t%f\t%f\t%f\t%f\n", x==0 ? 's' : 'i',
                    m+mm ? mm*100.0/(m+mm) : 0,
                    nrd_m[T2S(GT_HOM_RR)]+nrd_mm[T2S(GT_HOM_RR)] ? nrd_mm[T2S(GT_HOM_RR)]*100.0/(nrd_m[T2S(GT_HOM_RR)]+nrd_mm[T2S(GT_HOM_RR)]) : 0,
                    nrd_m[T2S(GT_HET_RA)]+nrd_mm[T2S(GT_HET_RA)] ? nrd_mm[T2S(GT_HET_RA)]*100.0/(nrd_m[T2S(GT_HET_RA)]+nrd_mm[T2S(GT_HET_RA)]) : 0,
                    nrd_m[T2S(GT_HOM_AA)]+nrd_mm[T2S(GT_HOM_AA)] ? nrd_mm[T2S(GT_HOM_AA)]*100.0/(nrd_m[T2S(GT_HOM_AA)]+nrd_mm[T2S(GT_HOM_AA)]) : 0
                  );
        }

        for (x=0; x<2; x++) // x=0: snps, x=1: indels
        {
            gtcmp_t *stats;
            if ( x==0 )
            {
                printf("# GCsS, Genotype concordance by sample (SNPs)\n# GCsS\t[2]id\t[3]sample\t[4]non-reference discordance rate\t[5]RR Hom matches\t[6]RA Het matches\t[7]AA Hom matches\t[8]RR Hom mismatches\t[9]RA Het mismatches\t[10]AA Hom mismatches\t[11]dosage r-squared\n");
                stats = args->smpl_gts_snps;
            }
            else
            {
                printf("# GCiS, Genotype concordance by sample (indels)\n# GCiS\t[2]id\t[3]sample\t[4]non-reference discordance rate\t[5]RR Hom matches\t[6]RA Het matches\t[7]AA Hom matches\t[8]RR Hom mismatches\t[9]RA Het mismatches\t[10]AA Hom mismatches\t[11]dosage r-squared\n");
                stats = args->smpl_gts_indels;
            }
            for (i=0; i<args->files->n_smpl; i++)
            {
                uint64_t mm = 0, m = stats[i].gt2gt[T2S(GT_HET_RA)][T2S(GT_HET_RA)] + stats[i].gt2gt[T2S(GT_HOM_AA)][T2S(GT_HOM_AA)];
                for (j=0; j<3; j++)
                    for (k=0; k<3; k++)
                        if ( j!=k ) mm += stats[i].gt2gt[j][k];

                // Pearson's r2
                double r2 = 0;
                if ( stats[i].n )
                {
                    r2  = (stats[i].yx - stats[i].x*stats[i].y/stats[i].n);
                    r2 /= sqrt((stats[i].xx - stats[i].x*stats[i].x/stats[i].n) * (stats[i].yy - stats[i].y*stats[i].y/stats[i].n));
                    r2 *= r2;
                }
                printf("GC%cS\t2\t%s\t%.3f",  x==0 ? 's' : 'i', args->files->samples[i], m+mm ? mm*100.0/(m+mm) : 0);
                printf("\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"",
                    stats[i].gt2gt[T2S(GT_HOM_RR)][T2S(GT_HOM_RR)],
                    stats[i].gt2gt[T2S(GT_HET_RA)][T2S(GT_HET_RA)],
                    stats[i].gt2gt[T2S(GT_HOM_AA)][T2S(GT_HOM_AA)]);
                printf("\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"",
                    stats[i].gt2gt[T2S(GT_HOM_RR)][T2S(GT_HET_RA)] + stats[i].gt2gt[T2S(GT_HOM_RR)][T2S(GT_HOM_AA)],
                    stats[i].gt2gt[T2S(GT_HET_RA)][T2S(GT_HOM_RR)] + stats[i].gt2gt[T2S(GT_HET_RA)][T2S(GT_HOM_AA)],
                    stats[i].gt2gt[T2S(GT_HOM_AA)][T2S(GT_HOM_RR)] + stats[i].gt2gt[T2S(GT_HOM_AA)][T2S(GT_HET_RA)]);
                if ( stats[i].n && !isnan(r2) ) printf("\t%f\n", r2);
                else printf("\t"NA_STRING"\n");
            }
        }
        for (x=0; x<2; x++) // x=0: snps, x=1: indels
        {
                //printf("# GCiS, Genotype concordance by sample (indels)\n# GCiS\t[2]id\t[3]sample\t[4]non-reference discordance rate\t[5]RR Hom matches\t[6]RA Het matches\t[7]AA Hom matches\t[8]RR Hom mismatches\t[9]RA Het mismatches\t[10]AA Hom mismatches\t[11]dosage r-squared\n");

            gtcmp_t *stats;
            if ( x==0 )
            {
                printf("# GCTs, Genotype concordance table (SNPs)\n# GCTs");
                stats = args->smpl_gts_snps;
            }
            else
            {
                printf("# GCTi, Genotype concordance table (indels)\n# GCTi");
                stats = args->smpl_gts_indels;
            }
            i = 1;
            printf("\t[%d]sample", ++i);
            printf("\t[%d]RR Hom -> RR Hom", ++i);
            printf("\t[%d]RR Hom -> RA Het", ++i);
            printf("\t[%d]RR Hom -> AA Hom", ++i);
            printf("\t[%d]RR Hom -> AA Het", ++i);
            printf("\t[%d]RR Hom -> missing", ++i);
            printf("\t[%d]RA Het -> RR Hom", ++i);
            printf("\t[%d]RA Het -> RA Het", ++i);
            printf("\t[%d]RA Het -> AA Hom", ++i);
            printf("\t[%d]RA Het -> AA Het", ++i);
            printf("\t[%d]RA Het -> missing", ++i);
            printf("\t[%d]AA Hom -> RR Hom", ++i);
            printf("\t[%d]AA Hom -> RA Het", ++i);
            printf("\t[%d]AA Hom -> AA Hom", ++i);
            printf("\t[%d]AA Hom -> AA Het", ++i);
            printf("\t[%d]AA Hom -> missing", ++i);
            printf("\t[%d]AA Het -> RR Hom", ++i);
            printf("\t[%d]AA Het -> RA Het", ++i);
            printf("\t[%d]AA Het -> AA Hom", ++i);
            printf("\t[%d]AA Het -> AA Het", ++i);
            printf("\t[%d]AA Het -> missing", ++i);
            printf("\t[%d]missing -> RR Hom", ++i);
            printf("\t[%d]missing -> RA Het", ++i);
            printf("\t[%d]missing -> AA Hom", ++i);
            printf("\t[%d]missing -> AA Het", ++i);
            printf("\t[%d]missing -> missing\n", ++i);

            for (i=0; i<args->files->n_smpl; i++)
            {
                printf("GCT%c\t%s",  x==0 ? 's' : 'i', args->files->samples[i]);
                for (j=0; j<5; j++)
                    for (k=0; k<5; k++)
                        printf("\t%"PRIu64, stats[i].gt2gt[j][k]);
                printf("\n");
            }
        }
    }

    printf("# DP, Depth distribution\n# DP\t[2]id\t[3]bin\t[4]number of genotypes\t[5]fraction of genotypes (%%)\t[6]number of sites\t[7]fraction of sites (%%)\n");
    for (id=0; id<args->nstats; id++)
    {
        stats_t *stats = &args->stats[id];
        long unsigned int sum = 0, sum_sites = 0;
        for (i=0; i<stats->dp.m_vals; i++) { sum += stats->dp.vals[i]; sum_sites += stats->dp_sites.vals[i]; }
        for (i=0; i<stats->dp.m_vals; i++)
        {
            if ( stats->dp.vals[i]==0 && stats->dp_sites.vals[i]==0 ) continue;
            printf("DP\t%d\t", id);
            if ( i==0 ) printf("<%d", stats->dp.min);
            else if ( i+1==stats->dp.m_vals ) printf(">%d", stats->dp.max);
            else printf("%d", idist_i2bin(&stats->dp,i));
            printf("\t%"PRIu64"\t%f", stats->dp.vals[i], sum ? stats->dp.vals[i]*100./sum : 0);
            printf("\t%"PRIu64"\t%f\n", stats->dp_sites.vals[i], sum_sites ? stats->dp_sites.vals[i]*100./sum_sites : 0);
        }
    }

    if ( args->files->n_smpl )
    {
        printf("# PSC, Per-sample counts. Note that the ref/het/hom counts include only SNPs, for indels see PSI. The rest include both SNPs and indels.\n");
        printf("# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTransitions\t[8]nTransversions\t[9]nIndels\t[10]average depth\t[11]nSingletons"
            "\t[12]nHapRef\t[13]nHapAlt\t[14]nMissing\n");
        for (id=0; id<args->nstats; id++)
        {
            stats_t *stats = &args->stats[id];
            for (i=0; i<args->files->n_smpl; i++)
            {
                float dp = stats->smpl_ndp[i] ? stats->smpl_dp[i]/(float)stats->smpl_ndp[i] : 0;
                printf("PSC\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%d\t%d\t%d\t%d\n", id,args->files->samples[i],
                    stats->smpl_homRR[i], stats->smpl_homAA[i], stats->smpl_hets[i], stats->smpl_ts[i],
                    stats->smpl_tv[i], stats->smpl_indels[i],dp, stats->smpl_sngl[i], stats->smpl_hapRef[i],
                    stats->smpl_hapAlt[i], stats->smpl_missing[i]);
            }
        }

        printf("# PSI, Per-Sample Indels. Note that alt-het genotypes with both ins and del allele are counted twice, in both nInsHets and nDelHets.\n");
        printf("# PSI\t[2]id\t[3]sample\t[4]in-frame\t[5]out-frame\t[6]not applicable\t[7]out/(in+out) ratio\t[8]nInsHets\t[9]nDelHets\t[10]nInsAltHoms\t[11]nDelAltHoms\n");
        for (id=0; id<args->nstats; id++)
        {
            stats_t *stats = &args->stats[id];
            for (i=0; i<args->files->n_smpl; i++)
            {
                int na = 0, in = 0, out = 0;
                if ( args->exons )
                {
                    na  = stats->smpl_frm_shifts[i*3 + 0];
                    in  = stats->smpl_frm_shifts[i*3 + 1];
                    out = stats->smpl_frm_shifts[i*3 + 2];
                }
                printf("PSI\t%d\t%s\t%d\t%d\t%d\t%.2f\t%d\t%d\t%d\t%d\n", id,args->files->samples[i], in,out,na,in+out?1.0*out/(in+out):0,
                    stats->smpl_ins_hets[i],stats->smpl_del_hets[i],stats->smpl_ins_homs[i],stats->smpl_del_homs[i]);
            }
        }

        #ifdef HWE_STATS
        printf("# HWE\n# HWE\t[2]id\t[3]1st ALT allele frequency\t[4]Number of observations\t[5]25th percentile\t[6]median\t[7]75th percentile\n");
        for (id=0; id<args->nstats; id++)
        {
            stats_t *stats = &args->stats[id];
            for (i=0; i<args->naf_hwe; i++) stats->af_hwe[i+args->naf_hwe] += stats->af_hwe[i]; // singletons
            for (i=1; i<args->m_af; i++)
            {
                unsigned int sum_tot = 0, sum_tmp = 0;
                int j, *ptr = &stats->af_hwe[i*args->naf_hwe];
                for (j=0; j<args->naf_hwe; j++) sum_tot += ptr[j];
                if ( !sum_tot ) continue;

                double af = args->af_bins ? (bin_get_value(args->af_bins,i)+bin_get_value(args->af_bins,i-1))*0.5 : (double)(i-1)/(args->m_af-1);

                int nprn = 3;
                printf("HWE\t%d\t%f\t%d",id,af,sum_tot);
                for (j=0; j<args->naf_hwe; j++)
                {
                    sum_tmp += ptr[j];
                    float frac = (float)sum_tmp/sum_tot;
                    if ( frac >= 0.75 )
                    {
                        while (nprn>0) { printf("\t%f", (float)j/args->naf_hwe); nprn--; }
                        break;
                    }
                    if ( frac >= 0.5 )
                    {
                        while (nprn>1) { printf("\t%f", (float)j/args->naf_hwe); nprn--; }
                        continue;
                    }
                    if ( frac >= 0.25 )
                    {
                        while (nprn>2) { printf("\t%f", (float)j/args->naf_hwe); nprn--; }
                    }
                }
                assert(nprn==0);
                printf("\n");
            }
        }
        #endif
    }
}

static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Parses VCF or BCF and produces stats which can be plotted using plot-vcfstats.\n");
    fprintf(stderr, "         When two files are given, the program generates separate stats for intersection\n");
    fprintf(stderr, "         and the complements. By default only sites are compared, -s/-S must given to include\n");
    fprintf(stderr, "         also sample columns.\n");
    fprintf(stderr, "Usage:   bcftools stats [options] <A.vcf.gz> [<B.vcf.gz>]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "        --af-bins LIST               Allele frequency bins, a list (0.1,0.5,1) or a file (0.1\\n0.5\\n1)\n");
    fprintf(stderr, "        --af-tag STRING              Allele frequency tag to use, by default estimated from AN,AC or GT\n");
    fprintf(stderr, "    -1, --1st-allele-only            Include only 1st allele at multiallelic sites\n");
    fprintf(stderr, "    -c, --collapse STRING            Treat as identical records with <snps|indels|both|all|some|none>, see man page for details [none]\n");
    fprintf(stderr, "    -d, --depth INT,INT,INT          Depth distribution: min,max,bin size [0,500,1]\n");
    fprintf(stderr, "    -e, --exclude EXPR               Exclude sites for which the expression is true (see man page for details)\n");
    fprintf(stderr, "    -E, --exons FILE.gz              Tab-delimited file with exons for indel frameshifts (chr,beg,end; 1-based, inclusive, bgzip compressed)\n");
    fprintf(stderr, "    -f, --apply-filters LIST         Require at least one of the listed FILTER strings (e.g. \"PASS,.\")\n");
    fprintf(stderr, "    -F, --fasta-ref FILE             Faidx indexed reference sequence file to determine INDEL context\n");
    fprintf(stderr, "    -i, --include EXPR               Select sites for which the expression is true (see man page for details)\n");
    fprintf(stderr, "    -I, --split-by-ID                Collect stats for sites with ID separately (known vs novel)\n");
    fprintf(stderr, "    -r, --regions REGION             Restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file FILE          Restrict to regions listed in a file\n");
    fprintf(stderr, "        --regions-overlap 0|1|2      Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n");
    fprintf(stderr, "    -s, --samples LIST               List of samples for sample stats, \"-\" to include all samples\n");
    fprintf(stderr, "    -S, --samples-file FILE          File of samples to include\n");
    fprintf(stderr, "    -t, --targets REGION             Similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "    -T, --targets-file FILE          Similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "        --targets-overlap 0|1|2      Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n");
    fprintf(stderr, "    -u, --user-tstv TAG[:min:max:n]  Collect Ts/Tv stats for any tag using the given binning [0:1:100]\n");
    fprintf(stderr, "                                       A subfield can be selected as e.g. 'PV4[0]', here the first value of the PV4 tag\n");
    fprintf(stderr, "        --threads INT                Use multithreading with <int> worker threads [0]\n");
    fprintf(stderr, "    -v, --verbose                    Produce verbose per-site and per-sample output\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfstats(int argc, char *argv[])
{
    int c;
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->files  = bcf_sr_init();
    args->argc   = argc; args->argv = argv;
    args->dp_min = 0; args->dp_max = 500; args->dp_step = 1;
    int regions_is_file = 0, targets_is_file = 0;
    int regions_overlap = 1;
    int targets_overlap = 0;

    static struct option loptions[] =
    {
        {"af-bins",1,0,1},
        {"af-tag",1,0,2},
        {"1st-allele-only",0,0,'1'},
        {"include",1,0,'i'},
        {"exclude",1,0,'e'},
        {"help",0,0,'h'},
        {"collapse",1,0,'c'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"regions-overlap",required_argument,NULL,3},
        {"verbose",0,0,'v'},
        {"depth",1,0,'d'},
        {"apply-filters",1,0,'f'},
        {"exons",1,0,'E'},
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {"split-by-ID",0,0,'I'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"targets-overlap",required_argument,NULL,4},
        {"fasta-ref",1,0,'F'},
        {"user-tstv",1,0,'u'},
        {"threads",1,0,9},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "hc:r:R:e:s:S:d:i:t:T:F:f:1u:vIE:",loptions,NULL)) >= 0) {
        switch (c) {
            case  1 : args->af_bins_list = optarg; break;
            case  2 : args->af_tag = optarg; break;
            case 'u': add_user_stats(args,optarg); break;
            case '1': args->first_allele_only = 1; break;
            case 'F': args->ref_fname = optarg; break;
            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; targets_is_file = 1; break;
            case 'c':
                if ( !strcmp(optarg,"snps") ) args->files->collapse |= COLLAPSE_SNPS;
                else if ( !strcmp(optarg,"indels") ) args->files->collapse |= COLLAPSE_INDELS;
                else if ( !strcmp(optarg,"both") ) args->files->collapse |= COLLAPSE_SNPS | COLLAPSE_INDELS;
                else if ( !strcmp(optarg,"any") ) args->files->collapse |= COLLAPSE_ANY;
                else if ( !strcmp(optarg,"all") ) args->files->collapse |= COLLAPSE_ANY;
                else if ( !strcmp(optarg,"some") ) args->files->collapse |= COLLAPSE_SOME;
                else if ( !strcmp(optarg,"none") ) args->files->collapse = COLLAPSE_NONE;
                else error("The --collapse string \"%s\" not recognised.\n", optarg);
                break;
            case 'v': args->verbose_sites = 1; break;
            case 'd':
                if ( sscanf(optarg,"%d,%d,%d",&args->dp_min,&args->dp_max,&args->dp_step)!=3 )
                    error("Could not parse --depth %s\n", optarg);
                if ( args->dp_min<0 || args->dp_min >= args->dp_max || args->dp_step > args->dp_max - args->dp_min + 1 )
                    error("Is this a typo? --depth %s\n", optarg);
                break;
            case 'f': args->files->apply_filters = optarg; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case 'E': args->exons_fname = optarg; break;
            case 's': args->samples_list = optarg; break;
            case 'S': args->samples_list = optarg; args->samples_is_file = 1; break;
            case 'I': args->split_by_id = 1; break;
            case 'e':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case  3 :
                regions_overlap = parse_overlap_option(optarg);
                if ( regions_overlap < 0 ) error("Could not parse: --regions-overlap %s\n",optarg);
                break;
            case  4 :
                targets_overlap = parse_overlap_option(optarg);
                if ( targets_overlap < 0 ) error("Could not parse: --targets-overlap %s\n",optarg);
                break;
            case  9 : args->n_threads = strtol(optarg, 0, 0); break;
            case 'h':
            case '?': usage(); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    char *fname = NULL;
    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin
        else usage();
    }
    else fname = argv[optind];

    if ( argc-optind>2 ) usage();
    if ( argc-optind>1 )
    {
        args->files->require_index = 1;
        if ( args->split_by_id ) error("Only one file can be given with -i.\n");
    }
    if ( !args->samples_list ) args->files->max_unpack = BCF_UN_INFO;
    else args->files->max_unpack = BCF_UN_FMT;
    if ( args->targets_list )
    {
        bcf_sr_set_opt(args->files,BCF_SR_TARGETS_OVERLAP,targets_overlap);
        if ( bcf_sr_set_targets(args->files, args->targets_list, targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
    }
    if ( args->regions_list)
    {
        bcf_sr_set_opt(args->files,BCF_SR_REGIONS_OVERLAP,regions_overlap);
        if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    if ( args->n_threads && bcf_sr_set_threads(args->files, args->n_threads)<0)
        error("Failed to create threads\n");

    while (fname)
    {
        if ( !bcf_sr_add_reader(args->files, fname) )
            error("Failed to read from %s: %s\n", !strcmp("-",fname)?"standard input":fname,bcf_sr_strerror(args->files->errnum));
        fname = ++optind < argc ? argv[optind] : NULL;
    }

    init_stats(args);
    print_header(args);
    do_vcf_stats(args);
    print_stats(args);
    destroy_stats(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}

