/*  vcfstats.c -- Produces stats which can be plotted using plot-vcfstats.

    Copyright (C) 2012-2015 Genome Research Ltd.

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
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <htslib/faidx.h>
#include <inttypes.h>
#include "bcftools.h"
#include "filter.h"

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
    int nbins, type, m_val;
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
    double x;
    double x2;
    double y;
    double y2;
    double xy;
    double n;
}
smpl_r_t;

typedef struct
{
    int n_snps, n_indels, n_mnps, n_others, n_mals, n_snp_mals, n_records, n_noalts;
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
        int *qual_ts, *qual_tv, *qual_snps, *qual_indels;
    #endif
    int *insertions, *deletions, m_indel;   // maximum indel length
    int in_frame, out_frame, na_frame, in_frame_alt1, out_frame_alt1, na_frame_alt1;
    int subst[15];
    int *smpl_hets, *smpl_homRR, *smpl_homAA, *smpl_ts, *smpl_tv, *smpl_indels, *smpl_ndp, *smpl_sngl;
    int *smpl_indel_hets, *smpl_indel_homs;
    int *smpl_frm_shifts; // not-applicable, in-frame, out-frame
    unsigned long int *smpl_dp;
    idist_t dp, dp_sites;
    int nusr;
    user_stats_t *usr;
}
stats_t;

typedef struct
{
    uint64_t m[3], mm[3];        // number of hom, het and non-ref hom matches and mismatches
    float r2sum;
    uint32_t r2n;
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
    gtcmp_t *af_gts_snps, *af_gts_indels, *smpl_gts_snps, *smpl_gts_indels; // first bin of af_* stats are singletons

    // indel context
    indel_ctx_t *indel_ctx;
    char *ref_fname;

    // user stats
    int nusr;
    user_stats_t *usr;

    // other
    bcf_srs_t *files;
    bcf_sr_regions_t *exons;
    char **argv, *exons_fname, *regions_list, *samples_list, *targets_list;
    int argc, verbose_sites, first_allele_only, samples_is_file;
    int split_by_id, nstats;

    filter_t *filter[2];
    char *filter_str;
    int filter_logic;   // include or exclude sites which match the filters? One of FLT_INCLUDE/FLT_EXCLUDE

    // Per Sample r working data arrays of size equal to number of samples
    smpl_r_t* smpl_r_snps;
    smpl_r_t* smpl_r_indels;
}
args_t;

static int type2dosage[6], type2ploidy[6], type2stats[6];

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
        if ( ref[i] != fai_ref[i] && ref[i] - 32 != fai_ref[i] )
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
    usr->min  = 0;
    usr->max  = 1;
    usr->nbins = 100;

    char *tmp = str;
    while ( *tmp && *tmp!=':' ) tmp++;
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
        if ( usr->type!=BCF_HT_REAL && usr->type!=BCF_HT_INT ) error("The INFO tag \"%s\" is not of Float or Integer type (%d)\n", usr->type);
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
    }

    // AF corresponds to AC but is more robust for mixture of haploid and diploid GTs
    args->m_af = 101;
    for (i=0; i<args->files->nreaders; i++)
        if ( bcf_hdr_nsamples(args->files->readers[i].header) + 1> args->m_af )
            args->m_af = bcf_hdr_nsamples(args->files->readers[i].header) + 1;

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
        args->smpl_r_snps = (smpl_r_t*) calloc(args->files->n_smpl, sizeof(smpl_r_t));
        args->smpl_r_indels = (smpl_r_t*) calloc(args->files->n_smpl, sizeof(smpl_r_t));
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
            stats->qual_ts     = (int*) calloc(args->m_qual,sizeof(int));
            stats->qual_tv     = (int*) calloc(args->m_qual,sizeof(int));
            stats->qual_snps   = (int*) calloc(args->m_qual,sizeof(int));
            stats->qual_indels = (int*) calloc(args->m_qual,sizeof(int));
        #endif
        if ( args->files->n_smpl )
        {
            stats->smpl_hets   = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_homAA  = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_homRR  = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_indel_hets = (int *) calloc(args->files->n_smpl,sizeof(int));
            stats->smpl_indel_homs = (int *) calloc(args->files->n_smpl,sizeof(int));
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
    type2stats[GT_HET_AA] = 1;
    type2stats[GT_HAPL_R] = 0;
    type2stats[GT_HAPL_A] = 2;

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
            if (stats->qual_ts) free(stats->qual_ts);
            if (stats->qual_tv) free(stats->qual_tv);
            if (stats->qual_snps) free(stats->qual_snps);
            if (stats->qual_indels) free(stats->qual_indels);
        #endif
        #if HWE_STATS
            //if ( args->files->n_smpl ) free(stats->af_hwe);
            free(stats->af_hwe);
        #endif
        free(stats->insertions);
        free(stats->deletions);
        if (stats->smpl_hets) free(stats->smpl_hets);
        if (stats->smpl_homAA) free(stats->smpl_homAA);
        if (stats->smpl_homRR) free(stats->smpl_homRR);
        if (stats->smpl_indel_homs) free(stats->smpl_indel_homs);
        if (stats->smpl_indel_hets) free(stats->smpl_indel_hets);
        if (stats->smpl_ts) free(stats->smpl_ts);
        if (stats->smpl_tv) free(stats->smpl_tv);
        if (stats->smpl_indels) free(stats->smpl_indels);
        if (stats->smpl_dp) free(stats->smpl_dp);
        if (stats->smpl_ndp) free(stats->smpl_ndp);
        if (stats->smpl_sngl) free(stats->smpl_sngl);
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
    }
    for (j=0; j<args->nusr; j++) free(args->usr[j].tag);
    free(args->usr);
    free(args->tmp_frm);
    free(args->tmp_iaf);
    if (args->exons) bcf_sr_regions_destroy(args->exons);
    free(args->af_gts_snps);
    free(args->af_gts_indels);
    free(args->smpl_gts_snps);
    free(args->smpl_gts_indels);
    free(args->smpl_r_snps);
    free(args->smpl_r_indels);
    if (args->indel_ctx) indel_ctx_destroy(args->indel_ctx);
    if (args->filter[0]) filter_destroy(args->filter[0]);
    if (args->filter[1]) filter_destroy(args->filter[1]);
}

static void init_iaf(args_t *args, bcf_sr_t *reader)
{
    bcf1_t *line = reader->buffer[0];
    if ( args->ntmp_iaf < line->n_allele )
    {
        args->tmp_iaf = (int*)realloc(args->tmp_iaf, line->n_allele*sizeof(int));
        args->ntmp_iaf = line->n_allele;
    }
    // tmp_iaf is first filled with AC counts in calc_ac and then transformed to
    //  an index to af_gts_snps
    int i, ret = bcf_calc_ac(reader->header, line, args->tmp_iaf, args->samples_list ? BCF_UN_INFO|BCF_UN_FMT : BCF_UN_INFO);
    if ( ret )
    {
        int an=0;
        for (i=0; i<line->n_allele; i++)
            an += args->tmp_iaf[i];

        args->tmp_iaf[0] = 0;
        for (i=1; i<line->n_allele; i++)
        {
            if ( args->tmp_iaf[i]==1 )
                args->tmp_iaf[i] = 0; // singletons into the first bin
            else if ( !an )
                args->tmp_iaf[i] = 1;   // no genotype at all, put to the AF=0 bin
            else
                args->tmp_iaf[i] = 1 + args->tmp_iaf[i] * (args->m_af-2.0) / an;
        }
    }
    else
        for (i=0; i<line->n_allele; i++)
            args->tmp_iaf[i] = 0;

    // todo: otherwise use AF
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
        int iqual = line->qual >= args->m_qual || isnan(line->qual) ? args->m_qual - 1 : line->qual;
        stats->qual_indels[iqual]++;
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
        if ( bcf_get_variant_type(line,i)!=VCF_INDEL ) continue;
        int len = line->d.var[i].n;

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
    int i;
    for (i=0; i<stats->nusr; i++)
    {
        user_stats_t *usr = &stats->usr[i];
        uint64_t *vals = is_ts ? usr->vals_ts : usr->vals_tv;
        float val;
        if ( usr->type==BCF_HT_REAL )
        {
            if ( bcf_get_info_float(reader->header,reader->buffer[0],usr->tag,&usr->val,&usr->m_val)<=0 ) continue;
            val = ((float*)usr->val)[0];
        }
        else
        {
            if ( bcf_get_info_int32(reader->header,reader->buffer[0],usr->tag,&usr->val,&usr->m_val)<=0 ) continue;
            val = ((int32_t*)usr->val)[0];
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
        int iqual = line->qual >= args->m_qual || isnan(line->qual) ? args->m_qual - 1 : line->qual;
        stats->qual_snps[iqual]++;
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
                    stats->qual_ts[iqual]++;
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
                    stats->qual_tv[iqual]++;
                #endif
                do_user_stats(stats, reader, 0);
            }
            stats->af_tv[iaf]++;
        }
    }
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
        int ref = bcf_acgt2int(*line->d.allele[0]);
        int is, n_nref = 0, i_nref = 0;
        for (is=0; is<args->files->n_smpl; is++)
        {
            int ial, jal;
            int gt = bcf_gt_type(fmt_ptr, reader->samples[is], &ial, &jal);
            if ( gt==GT_UNKN ) continue;
            if ( gt==GT_HAPL_R || gt==GT_HAPL_A )
            {
                if ( line_type&VCF_INDEL && stats->smpl_frm_shifts )
                {
                    assert( ial<line->n_allele );
                    stats->smpl_frm_shifts[is*3 + args->tmp_frm[ial]]++;
                }
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
            if ( line_type&VCF_SNP || line_type==VCF_REF )  // count ALT=. as SNP
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
            if ( line_type&VCF_INDEL )
            {
                if ( gt != GT_HOM_RR )
                {
                    stats->smpl_indels[is]++;
                    if ( gt==GT_HET_RA || gt==GT_HET_AA ) stats->smpl_indel_hets[is]++;
                    else if ( gt==GT_HOM_AA ) stats->smpl_indel_homs[is]++;
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

    if ( matched==3 )
    {
        int is;
        bcf_fmt_t *fmt0, *fmt1;
        fmt0 = bcf_get_fmt(files->readers[0].header,files->readers[0].buffer[0],"GT"); if ( !fmt0 ) return;
        fmt1 = bcf_get_fmt(files->readers[1].header,files->readers[1].buffer[0],"GT"); if ( !fmt1 ) return;

        // only the first ALT allele is considered
        int iaf = line->n_allele>1 ? args->tmp_iaf[1] : 1;
        int line_type = bcf_get_variant_types(files->readers[0].buffer[0]);
        gtcmp_t *af_stats = line_type&VCF_SNP ? args->af_gts_snps : args->af_gts_indels;
        gtcmp_t *smpl_stats = line_type&VCF_SNP ? args->smpl_gts_snps : args->smpl_gts_indels;

        //
        // Calculates r squared
        // x is mean dosage of x at given site
        // x2 is mean squared dosage of x at given site
        // y is mean dosage of x at given site
        // y2 is mean squared dosage of x at given site
        // xy is mean dosage of x*y at given site
        // r2sum += (xy - x*y)^2 / ( (x2 - x^2) * (y2 - y^2) )
        // r2n is number of sites considered
        // output as r2sum/r2n for each AF bin
        int r2n = 0;
        float x = 0, y = 0, xy = 0, x2 = 0, y2 = 0;
        // Select smpl_r
        smpl_r_t *smpl_r = NULL;
        if (line_type&VCF_SNP)
        {
            smpl_r = args->smpl_r_snps;
        }
        else if (line_type&VCF_INDEL)
        {
            smpl_r = args->smpl_r_indels;
        }
        for (is=0; is<files->n_smpl; is++)
        {
            // Simplified comparison: only 0/0, 0/1, 1/1 is looked at as the identity of
            //  actual alleles can be enforced by running without the -c option.
            int gt0 = bcf_gt_type(fmt0, files->readers[0].samples[is], NULL, NULL);
            if ( gt0 == GT_UNKN ) continue;

            int gt1 = bcf_gt_type(fmt1, files->readers[1].samples[is], NULL, NULL);
            if ( gt1 == GT_UNKN ) continue;

            if ( type2ploidy[gt0]*type2ploidy[gt1] == -1 ) continue;   // cannot compare diploid and haploid genotypes

            int dsg0 = type2dosage[gt0];
            int dsg1 = type2dosage[gt1];
            x   += dsg0;
            x2  += dsg0*dsg0;
            y   += dsg1;
            y2  += dsg1*dsg1;
            xy  += dsg0*dsg1;
            r2n++;

            int idx = type2stats[gt0];
            if ( gt0==gt1 )
            {
                af_stats[iaf].m[idx]++;
                smpl_stats[is].m[idx]++;
            }
            else
            {
                af_stats[iaf].mm[idx]++;
                smpl_stats[is].mm[idx]++;
            }

            // Now do it across samples

            if (smpl_r) {
                smpl_r[is].xy += dsg0*dsg1;
                smpl_r[is].x += dsg0;
                smpl_r[is].x2 += dsg0*dsg0;
                smpl_r[is].y += dsg1;
                smpl_r[is].y2 += dsg1*dsg1;
                ++(smpl_r[is].n);
            }
        }

        if ( r2n )
        {
            x /= r2n; y /= r2n; x2 /= r2n; y2 /= r2n; xy /= r2n;
            float cov  = xy - x*y;
            float var2 = (x2 - x*x) * (y2 - y*y);
            if ( var2!=0 )
            {
                af_stats[iaf].r2sum += cov*cov/var2;
                af_stats[iaf].r2n++;
            }
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
                    printf("DBG\t%s\t%d\t%s\t%d\t%d\n",reader->header->id[BCF_DT_CTG][reader->buffer[0]->rid].key,reader->buffer[0]->pos+1,files->samples[is],gt,gt2);
                }
                else
                {
                    if ( gt!=GT_HOM_RR ) nrefm++;
                    nm++;
                }
            }
            float nrd = nrefm+nmm ? 100.*nmm/(nrefm+nmm) : 0;
            printf("PSD\t%s\t%d\t%d\t%d\t%f\n", reader->header->id[BCF_DT_CTG][reader->buffer[0]->rid].key,reader->buffer[0]->pos+1,nm,nmm,nrd);
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
            if ( line_type == VCF_SNP ) stats->n_snp_mals++;
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
    int i, id;
    printf("# SN, Summary numbers:\n# SN\t[2]id\t[3]key\t[4]value\n");
    for (id=0; id<args->files->nreaders; id++)
        printf("SN\t%d\tnumber of samples:\t%d\n", id, bcf_hdr_nsamples(args->files->readers[id].header));
    for (id=0; id<args->nstats; id++)
    {
        stats_t *stats = &args->stats[id];
        printf("SN\t%d\tnumber of records:\t%d\n", id, stats->n_records);
        printf("SN\t%d\tnumber of no-ALTs:\t%d\n", id, stats->n_noalts);
        printf("SN\t%d\tnumber of SNPs:\t%d\n", id, stats->n_snps);
        printf("SN\t%d\tnumber of MNPs:\t%d\n", id, stats->n_mnps);
        printf("SN\t%d\tnumber of indels:\t%d\n", id, stats->n_indels);
        printf("SN\t%d\tnumber of others:\t%d\n", id, stats->n_others);
        printf("SN\t%d\tnumber of multiallelic sites:\t%d\n", id, stats->n_mals);
        printf("SN\t%d\tnumber of multiallelic SNP sites:\t%d\n", id, stats->n_snp_mals);
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
    printf("# AF, Stats by non-reference allele frequency:\n# AF\t[2]id\t[3]allele frequency\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable\n");
    for (id=0; id<args->nstats; id++)
    {
        stats_t *stats = &args->stats[id];
        for (i=1; i<args->m_af; i++) // note that af[1] now contains also af[0], see SiS stats output above
        {
            if ( stats->af_snps[i]+stats->af_ts[i]+stats->af_tv[i]+stats->af_repeats[0][i]+stats->af_repeats[1][i]+stats->af_repeats[2][i] == 0  ) continue;
            printf("AF\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", id,100.*(i-1)/(args->m_af-1),stats->af_snps[i],stats->af_ts[i],stats->af_tv[i],
                stats->af_repeats[0][i]+stats->af_repeats[1][i]+stats->af_repeats[2][i],stats->af_repeats[0][i],stats->af_repeats[1][i],stats->af_repeats[2][i]);
        }
    }
    #if QUAL_STATS
        printf("# QUAL, Stats by quality:\n# QUAL\t[2]id\t[3]Quality\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t[6]number of transversions (1st ALT)\t[7]number of indels\n");
        for (id=0; id<args->nstats; id++)
        {
            stats_t *stats = &args->stats[id];
            for (i=0; i<args->m_qual; i++)
            {
                if ( stats->qual_snps[i]+stats->qual_ts[i]+stats->qual_tv[i]+stats->qual_indels[i] == 0  ) continue;
                printf("QUAL\t%d\t%d\t%d\t%d\t%d\t%d\n", id,i,stats->qual_snps[i],stats->qual_ts[i],stats->qual_tv[i],stats->qual_indels[i]);
            }
        }
    #endif
    for (i=0; i<args->nusr; i++)
    {
        printf("# USR:%s, Stats by %s:\n# USR:%s\t[2]id\t[3]%s\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t[6]number of transversions (1st ALT)\n",
            args->usr[i].tag,args->usr[i].tag,args->usr[i].tag,args->usr[i].tag);
        for (id=0; id<args->nstats; id++)
        {
            user_stats_t *usr = &args->stats[id].usr[i];
            int j;
            for (j=0; j<usr->nbins; j++)
            {
                if ( usr->vals_ts[j]+usr->vals_tv[j] == 0 ) continue;   // skip empty bins
                float val = usr->min + (usr->max - usr->min)*j/(usr->nbins-1);
                const char *fmt = usr->type==BCF_HT_REAL ? "USR:%s\t%d\t%e\t%d\t%d\t%d\n" : "USR:%s\t%d\t%.0f\t%d\t%d\t%d\n";
                printf(fmt,usr->tag,id,val,usr->vals_ts[j]+usr->vals_tv[j],usr->vals_ts[j],usr->vals_tv[j]);
            }
        }
    }
    printf("# IDD, InDel distribution:\n# IDD\t[2]id\t[3]length (deletions negative)\t[4]count\n");
    for (id=0; id<args->nstats; id++)
    {
        stats_t *stats = &args->stats[id];
        for (i=stats->m_indel-1; i>=0; i--)
            if ( stats->deletions[i] ) printf("IDD\t%d\t%d\t%d\n", id,-i-1,stats->deletions[i]);
        for (i=0; i<stats->m_indel; i++)
            if ( stats->insertions[i] ) printf("IDD\t%d\t%d\t%d\n", id,i+1,stats->insertions[i]);
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
        for (x=0; x<2; x++)
        {
            gtcmp_t *stats;
            if ( x==0 )
            {
                printf("# GCsAF, Genotype concordance by non-reference allele frequency (SNPs)\n# GCsAF\t[2]id\t[3]allele frequency\t[4]RR Hom matches\t[5]RA Het matches\t[6]AA Hom matches\t[7]RR Hom mismatches\t[8]RA Het mismatches\t[9]AA Hom mismatches\t[10]dosage r-squared\t[11]number of sites\n");
                stats = args->af_gts_snps;
            }
            else
            {
                printf("# GCiAF, Genotype concordance by non-reference allele frequency (indels)\n# GCiAF\t[2]id\t[3]allele frequency\t[4]RR Hom matches\t[5]RA Het matches\t[6]AA Hom matches\t[7]RR Hom mismatches\t[8]RA Het mismatches\t[9]AA Hom mismatches\t[10]dosage r-squared\t[11]number of sites\n");
                stats = args->af_gts_indels;
            }
            uint64_t nrd_m[3] = {0,0,0}, nrd_mm[3] = {0,0,0};
            for (i=0; i<args->m_af; i++)
            {
                int j, n = 0;
                for (j=0; j<3; j++)
                {
                    n += stats[i].m[j] + stats[i].mm[j];
                    nrd_m[j]  += stats[i].m[j];
                    nrd_mm[j] += stats[i].mm[j];
                }
                if ( !i || !n ) continue;   // skip singleton stats and empty bins
                printf("GC%cAF\t2\t%f", x==0 ? 's' : 'i', 100.*(i-1)/(args->m_af-1));
                printf("\t%"PRId64"\t%"PRId64"\t%"PRId64"", stats[i].m[T2S(GT_HOM_RR)],stats[i].m[T2S(GT_HET_RA)],stats[i].m[T2S(GT_HOM_AA)]);
                printf("\t%"PRId64"\t%"PRId64"\t%"PRId64"", stats[i].mm[T2S(GT_HOM_RR)],stats[i].mm[T2S(GT_HET_RA)],stats[i].mm[T2S(GT_HOM_AA)]);
                printf("\t%f\t%"PRId32"\n", stats[i].r2n ? stats[i].r2sum/stats[i].r2n : -1.0, stats[i].r2n);
            }

            if ( x==0 )
            {
                printf("# NRD and discordance is calculated as follows:\n");
                printf("#   m .. number of matches\n");
                printf("#   x .. number of mismatches\n");
                printf("#   NRD = (xRR + xRA + xAA) / (xRR + xRA + xAA + mRA + mAA)\n");
                printf("#   RR discordance = xRR / (xRR + mRR)\n");
                printf("#   RA discordance = xRA / (xRA + mRA)\n");
                printf("#   AA discordance = xAA / (xAA + mAA)\n");
                printf("# Non-Reference Discordance (NRD), SNPs\n# NRDs\t[2]id\t[3]NRD\t[4]Ref/Ref discordance\t[5]Ref/Alt discordance\t[6]Alt/Alt discordance\n");
            }
            else
                printf("# Non-Reference Discordance (NRD), indels\n# NRDi\t[2]id\t[3]NRD\t[4]Ref/Ref discordance\t[5]Ref/Alt discordance\t[6]Alt/Alt discordance\n");
            uint64_t m  = nrd_m[T2S(GT_HET_RA)] + nrd_m[T2S(GT_HOM_AA)];
            uint64_t mm = nrd_mm[T2S(GT_HOM_RR)] + nrd_mm[T2S(GT_HET_RA)] + nrd_mm[T2S(GT_HOM_AA)];
            printf("NRD%c\t2\t%f\t%f\t%f\t%f\n", x==0 ? 's' : 'i',
                    m+mm ? mm*100.0/(m+mm) : 0,
                    nrd_m[T2S(GT_HOM_RR)]+nrd_mm[T2S(GT_HOM_RR)] ? nrd_mm[T2S(GT_HOM_RR)]*100.0/(nrd_m[T2S(GT_HOM_RR)]+nrd_mm[T2S(GT_HOM_RR)]) : 0,
                    nrd_m[T2S(GT_HET_RA)]+nrd_mm[T2S(GT_HET_RA)] ? nrd_mm[T2S(GT_HET_RA)]*100.0/(nrd_m[T2S(GT_HET_RA)]+nrd_mm[T2S(GT_HET_RA)]) : 0,
                    nrd_m[T2S(GT_HOM_AA)]+nrd_mm[T2S(GT_HOM_AA)] ? nrd_mm[T2S(GT_HOM_AA)]*100.0/(nrd_m[T2S(GT_HOM_AA)]+nrd_mm[T2S(GT_HOM_AA)]) : 0
                  );
        }

        for (x=0; x<2; x++)
        {
            gtcmp_t *stats;
            smpl_r_t *smpl_r_array;
            if ( x==0 )
            {
                printf("# GCsS, Genotype concordance by sample (SNPs)\n# GCsS\t[2]id\t[3]sample\t[4]non-reference discordance rate\t[5]RR Hom matches\t[6]RA Het matches\t[7]AA Hom matches\t[8]RR Hom mismatches\t[9]RA Het mismatches\t[10]AA Hom mismatches\t[11]dosage r-squared\n");
                stats = args->smpl_gts_snps;
                smpl_r_array = args->smpl_r_snps;
            }
            else
            {
                printf("# GCiS, Genotype concordance by sample (indels)\n# GCiS\t[2]id\t[3]sample\t[4]non-reference discordance rate\t[5]RR Hom matches\t[6]RA Het matches\t[7]AA Hom matches\t[8]RR Hom mismatches\t[9]RA Het mismatches\t[10]AA Hom mismatches\t[11]dosage r-squared\n");
                stats = args->smpl_gts_indels;
                smpl_r_array = args->smpl_r_indels;
            }
            for (i=0; i<args->files->n_smpl; i++)
            {
                uint64_t m  = stats[i].m[T2S(GT_HET_RA)] + stats[i].m[T2S(GT_HOM_AA)];
                uint64_t mm = stats[i].mm[T2S(GT_HOM_RR)] + stats[i].mm[T2S(GT_HET_RA)] + stats[i].mm[T2S(GT_HOM_AA)];
                // Calculate r by formula 19.2 - Biostatistical Analysis 4th edition - Jerrold H. Zar
                smpl_r_t *smpl_r = smpl_r_array + i;
                double r = 0.0;
                if (smpl_r->n) {
                    double sum_crossprod = smpl_r->xy-(smpl_r->x*smpl_r->y)/smpl_r->n;//per 17.3 machine formula
                    double x2_xx = smpl_r->x2-(smpl_r->x*smpl_r->x)/smpl_r->n;
                    double y2_yy = smpl_r->y2-(smpl_r->y*smpl_r->y)/smpl_r->n;
                    r = (sum_crossprod)/sqrt(x2_xx*y2_yy);
                }
                printf("GC%cS\t2\t%s\t%.3f",  x==0 ? 's' : 'i', args->files->samples[i], m+mm ? mm*100.0/(m+mm) : 0);
                printf("\t%"PRId64"\t%"PRId64"\t%"PRId64"", stats[i].m[T2S(GT_HOM_RR)],stats[i].m[T2S(GT_HET_RA)],stats[i].m[T2S(GT_HOM_AA)]);
                printf("\t%"PRId64"\t%"PRId64"\t%"PRId64"", stats[i].mm[T2S(GT_HOM_RR)],stats[i].mm[T2S(GT_HET_RA)],stats[i].mm[T2S(GT_HOM_AA)]);
                if (smpl_r->n && !isnan(r)) printf("\t%f\n", r*r);
                else printf("\t"NA_STRING"\n");
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
            printf("\t%"PRId64"\t%f", stats->dp.vals[i], sum ? stats->dp.vals[i]*100./sum : 0);
            printf("\t%"PRId64"\t%f\n", stats->dp_sites.vals[i], sum_sites ? stats->dp_sites.vals[i]*100./sum_sites : 0);
        }
    }

    if ( args->files->n_smpl )
    {
        printf("# PSC, Per-sample counts\n# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTransitions\t[8]nTransversions\t[9]nIndels\t[10]average depth\t[11]nSingletons\n");
        for (id=0; id<args->nstats; id++)
        {
            stats_t *stats = &args->stats[id];
            for (i=0; i<args->files->n_smpl; i++)
            {
                float dp = stats->smpl_ndp[i] ? stats->smpl_dp[i]/(float)stats->smpl_ndp[i] : 0;
                printf("PSC\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%d\n", id,args->files->samples[i],
                    stats->smpl_homRR[i], stats->smpl_homAA[i], stats->smpl_hets[i], stats->smpl_ts[i],
                    stats->smpl_tv[i], stats->smpl_indels[i],dp, stats->smpl_sngl[i]);
            }
        }


        printf("# PSI, Per-Sample Indels\n# PSI\t[2]id\t[3]sample\t[4]in-frame\t[5]out-frame\t[6]not applicable\t[7]out/(in+out) ratio\t[8]nHets\t[9]nAA\n");
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
                int nhom = stats->smpl_indel_homs[i];
                int nhet = stats->smpl_indel_hets[i];
                printf("PSI\t%d\t%s\t%d\t%d\t%d\t%.2f\t%d\t%d\n", id,args->files->samples[i], in,out,na,in+out?1.0*out/(in+out):0,nhet,nhom);
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

                int nprn = 3;
                printf("HWE\t%d\t%f\t%d",id,100.*(i-1)/(args->m_af-1),sum_tot);
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
    fprintf(stderr, "    -1, --1st-allele-only              include only 1st allele at multiallelic sites\n");
    fprintf(stderr, "    -c, --collapse <string>            treat as identical records with <snps|indels|both|all|some|none>, see man page for details [none]\n");
    fprintf(stderr, "    -d, --depth <int,int,int>          depth distribution: min,max,bin size [0,500,1]\n");
    fprintf(stderr, "    -e, --exclude <expr>               exclude sites for which the expression is true (see man page for details)\n");
    fprintf(stderr, "    -E, --exons <file.gz>              tab-delimited file with exons for indel frameshifts (chr,from,to; 1-based, inclusive, bgzip compressed)\n");
    fprintf(stderr, "    -f, --apply-filters <list>         require at least one of the listed FILTER strings (e.g. \"PASS,.\")\n");
    fprintf(stderr, "    -F, --fasta-ref <file>             faidx indexed reference sequence file to determine INDEL context\n");
    fprintf(stderr, "    -i, --include <expr>               select sites for which the expression is true (see man page for details)\n");
    fprintf(stderr, "    -I, --split-by-ID                  collect stats for sites with ID separately (known vs novel)\n");
    fprintf(stderr, "    -r, --regions <region>             restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file <file>          restrict to regions listed in a file\n");
    fprintf(stderr, "    -s, --samples <list>               list of samples for sample stats, \"-\" to include all samples\n");
    fprintf(stderr, "    -S, --samples-file <file>          file of samples to include\n");
    fprintf(stderr, "    -t, --targets <region>             similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "    -T, --targets-file <file>          similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "    -u, --user-tstv <TAG[:min:max:n]>  collect Ts/Tv stats for any tag using the given binning [0:1:100]\n");
    fprintf(stderr, "    -v, --verbose                      produce verbose per-site and per-sample output\n");
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

    static struct option loptions[] =
    {
        {"1st-allele-only",0,0,'1'},
        {"include",1,0,'i'},
        {"exclude",1,0,'e'},
        {"help",0,0,'h'},
        {"collapse",1,0,'c'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"verbose",0,0,'v'},
        {"depth",1,0,'d'},
        {"apply-filters",1,0,'f'},
        {"exons",1,0,'E'},
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {"split-by-ID",0,0,'I'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"fasta-ref",1,0,'F'},
        {"user-tstv",1,0,'u'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "hc:r:R:e:s:S:d:i:t:T:F:f:1u:vIE:",loptions,NULL)) >= 0) {
        switch (c) {
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
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'h':
            case '?': usage();
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
    if ( args->targets_list && bcf_sr_set_targets(args->files, args->targets_list, targets_is_file, 0)<0 )
        error("Failed to read the targets: %s\n", args->targets_list);
    if ( args->regions_list && bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
        error("Failed to read the regions: %s\n", args->regions_list);
    while (fname)
    {
        if ( !bcf_sr_add_reader(args->files, fname) )
            error("Failed to open %s: %s\n", fname,bcf_sr_strerror(args->files->errnum));
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

