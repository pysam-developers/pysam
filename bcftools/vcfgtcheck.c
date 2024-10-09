/*  vcfgtcheck.c -- Check sample identity.

    Copyright (C) 2013-2024 Genome Research Ltd.

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

#include <stdio.h>
#include <stdarg.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <htslib/kbitset.h>
#include <htslib/hts_os.h>
#include <htslib/bgzf.h>
#include <inttypes.h>
#include <sys/time.h>
#include "bcftools.h"
#include "extsort.h"
//#include "hclust.h"
#include "filter.h"

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef struct
{
    int iqry, igt;
}
pair_t;

typedef struct
{
    bcf_srs_t *files;           // first reader is the query VCF - single sample normally or multi-sample for cross-check
    bcf_hdr_t *gt_hdr, *qry_hdr; // VCF with genotypes to compare against and the query VCF
    char *cwd, **argv, *gt_samples, *qry_samples, *regions, *targets, *qry_fname, *gt_fname, *pair_samples;
    char *output_fname;
    int argc, gt_samples_is_file, qry_samples_is_file, regions_is_file, targets_is_file, pair_samples_is_file, output_type;;
    int regions_overlap, targets_overlap, clevel;
    int qry_use_GT,gt_use_GT, nqry_smpl,ngt_smpl, *qry_smpl,*gt_smpl;
    int nused[2][2];
    double *pdiff, *qry_prob, *gt_prob;
    uint32_t *ndiff,*ncnt,*nmatch,ncmp, npairs;
    int32_t *qry_arr,*gt_arr, nqry_arr,ngt_arr;
    uint8_t *qry_dsg, *gt_dsg;
    pair_t *pairs;
    double *hwe_prob, dsg2prob[8][3], pl2prob[256];
    double min_inter_err, max_intra_err;
    int all_sites, hom_only, ntop, cross_check, calc_hwe_prob, sort_by_hwe, dry_run, gt_err;
    BGZF *out_fh;
    unsigned int nskip_no_match, nskip_not_ba, nskip_mono, nskip_no_data, nskip_dip_GT, nskip_dip_PL, nskip_filter;
    kstring_t kstr;

    // for --distinctive-sites
    double distinctive_sites;
    kbitset_t *kbs_diff;
    size_t diff_sites_size;
    extsort_t *es;
    char *es_tmp_prefix, *es_max_mem;

    // include or exclude sites which match the filters
    filter_t *qry_filter, *gt_filter;
    char *qry_filter_str, *gt_filter_str;
    int qry_filter_logic, gt_filter_logic;       // FLT_INCLUDE or FLT_EXCLUDE
}
args_t;

static void set_cwd(args_t *args)
{
    int i;
    char *buf;
    size_t nbuf = 500;
    args->cwd = (char*) malloc(sizeof(char)*nbuf);
    for (i=0; i<5; i++)
    {
        if ( (buf = getcwd(args->cwd, nbuf)) ) break;
        nbuf *= 2;
        args->cwd = (char*) realloc(args->cwd, sizeof(char)*nbuf);
    }
    assert(buf);
}
static void print_header(args_t *args)
{
    args->kstr.l = 0;
    ksprintf(&args->kstr, "# This file was produced by bcftools (%s+htslib-%s), the command line was:\n", bcftools_version(), hts_version());
    ksprintf(&args->kstr, "# \t bcftools %s ", args->argv[0]);
    int i;
    for (i=1; i<args->argc; i++)
        ksprintf(&args->kstr, " %s",args->argv[i]);
    ksprintf(&args->kstr, "\n# and the working directory was:\n");
    ksprintf(&args->kstr, "# \t %s\n#\n", args->cwd);
    if ( bgzf_write(args->out_fh, args->kstr.s, args->kstr.l)!=args->kstr.l ) error("Failed to write to %s\n", args->output_fname);
}

static int cmp_int(const void *_a, const void *_b)
{
    int a = *((int*)_a);
    int b = *((int*)_b);
    if ( a < b ) return -1;
    if ( a > b ) return 1;
    return 0;
}
static int cmp_pair(const void *_a, const void *_b)
{
    pair_t *a = (pair_t*)_a;
    pair_t *b = (pair_t*)_b;
    if ( a->iqry < b->iqry ) return -1;
    if ( a->iqry > b->iqry ) return 1;
    if ( a->igt < b->igt ) return -1;
    if ( a->igt > b->igt ) return 1;
    return 0;
}

typedef struct
{
    uint32_t ndiff,rid,pos,rand; // rand is to shuffle sites with the same ndiff from across all chromosomes
    unsigned long kbs_dat[1];
}
diff_sites_t;
#if DBG
static void diff_sites_debug_print(args_t *args, diff_sites_t *ds)
{
    int i;
    memcpy(args->kbs_diff->b,ds->kbs_dat,args->kbs_diff->n*sizeof(unsigned long));
    fprintf(stderr,"%s:%d\t%d\t",bcf_hdr_id2name(args->qry_hdr,ds->rid),ds->pos+1,ds->ndiff);
    for (i=0; i<args->npairs; i++) fprintf(stderr,"%d",kbs_exists(args->kbs_diff,i)?1:0);
    fprintf(stderr,"\n");
}
#endif
static int diff_sites_cmp(const void *aptr, const void *bptr)
{
    diff_sites_t *a = *((diff_sites_t**)aptr);
    diff_sites_t *b = *((diff_sites_t**)bptr);
    if ( a->ndiff < b->ndiff ) return 1;        // descending order
    if ( a->ndiff > b->ndiff ) return -1;
    if ( a->rand < b->rand ) return -1;
    if ( a->rand > b->rand ) return 1;
    return 0;
}
static void diff_sites_init(args_t *args)
{
    int nsites = args->distinctive_sites<=1 ? args->npairs*args->distinctive_sites : args->distinctive_sites;
    if ( nsites<=0 ) error("The value for --distinctive-sites was set too low: %d\n",nsites);
    if ( nsites > args->npairs )
    {
        fprintf(stderr,"Warning: The value for --distinctive-sites is bigger than is the number of pairs, all discordant sites be printed.\n");
        nsites = args->npairs;
        args->distinctive_sites = args->npairs + 1;
    }
    else
        args->distinctive_sites = nsites;
    args->kbs_diff = kbs_init(args->npairs);
    size_t n = (args->npairs + KBS_ELTBITS-1) / KBS_ELTBITS;
    assert( n==args->kbs_diff->n );
    args->diff_sites_size = sizeof(diff_sites_t) + (n-1)*sizeof(unsigned long);
    args->es = extsort_alloc();
    extsort_set_opt(args->es,size_t,DAT_SIZE,args->diff_sites_size);
    extsort_set_opt(args->es,const char*,TMP_PREFIX,args->es_tmp_prefix);
    extsort_set_opt(args->es,const char*,MAX_MEM,args->es_max_mem);
    extsort_set_opt(args->es,extsort_cmp_f,FUNC_CMP,diff_sites_cmp);
    extsort_init(args->es);
}
static void diff_sites_destroy(args_t *args)
{
    kbs_destroy(args->kbs_diff);
    extsort_destroy(args->es);
}
static inline void diff_sites_reset(args_t *args)
{
    kbs_clear(args->kbs_diff);
}
static inline void diff_sites_push(args_t *args, int ndiff, int rid, int pos)
{
    diff_sites_t *dat = (diff_sites_t*) malloc(args->diff_sites_size);
    memset(dat,0,sizeof(*dat)); // for debugging: prevent warnings about uninitialized memory coming from struct padding (not needed after rand added)
    dat->ndiff = ndiff;
    dat->rid  = rid;
    dat->pos  = pos;
    dat->rand = hts_lrand48();
    memcpy(dat->kbs_dat,args->kbs_diff->b,args->kbs_diff->n*sizeof(unsigned long));
    extsort_push(args->es,dat);
}
static inline int diff_sites_shift(args_t *args, int *ndiff, int *rid, int *pos)
{
    diff_sites_t *dat = (diff_sites_t*) extsort_shift(args->es);
    if ( !dat ) return 0;
    *ndiff = dat->ndiff;
    *rid   = dat->rid;
    *pos   = dat->pos;
    memcpy(args->kbs_diff->b,dat->kbs_dat,args->kbs_diff->n*sizeof(unsigned long));
    return 1;
}

static void init_samples(char *list, int list_is_file, int **smpl, int *nsmpl, bcf_hdr_t *hdr, char *vcf_fname)
{
    int i;
    if ( !strcmp(list,"-") )
    {
        *nsmpl = bcf_hdr_nsamples(hdr);
        *smpl  = (int*) malloc(sizeof(**smpl)*(*nsmpl));
        for (i=0; i<*nsmpl; i++) (*smpl)[i] = i;
        return;
    }

    char **tmp = hts_readlist(list, list_is_file, nsmpl);
    if ( !tmp || !*nsmpl ) error("Failed to parse %s\n", list);
    *smpl = (int*) malloc(sizeof(**smpl)*(*nsmpl));
    for (i=0; i<*nsmpl; i++)
    {
        int idx = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, tmp[i]);
        if ( idx<0 ) error("No such sample in %s: [%s]\n",vcf_fname,tmp[i]);
        (*smpl)[i] = idx;
        free(tmp[i]);
    }
    free(tmp);
    qsort(*smpl,*nsmpl,sizeof(**smpl),cmp_int);
    // check for duplicates
    for (i=1; i<*nsmpl; i++)
        if ( (*smpl)[i-1]==(*smpl)[i] )
            error("Error: the sample \"%s\" is listed twice in %s\n", hdr->samples[(*smpl)[i]],list);
}

static void init_data(args_t *args)
{
    hts_srand48(0);

    args->files = bcf_sr_init();
    if ( args->regions )
    {
        bcf_sr_set_opt(args->files,BCF_SR_REGIONS_OVERLAP,args->regions_overlap);
        if ( bcf_sr_set_regions(args->files, args->regions, args->regions_is_file)<0 ) error("Failed to read the regions: %s\n", args->regions);
    }
    if ( args->targets )
    {
        bcf_sr_set_opt(args->files,BCF_SR_TARGETS_OVERLAP,args->targets_overlap);
        if ( bcf_sr_set_targets(args->files, args->targets, args->targets_is_file, 0)<0 ) error("Failed to read the targets: %s\n", args->targets);
    }

    if ( args->gt_fname ) bcf_sr_set_opt(args->files, BCF_SR_REQUIRE_IDX);
    if ( !bcf_sr_add_reader(args->files,args->qry_fname) ) error("Failed to open %s: %s\n", args->qry_fname,bcf_sr_strerror(args->files->errnum));
    if ( args->gt_fname && !bcf_sr_add_reader(args->files, args->gt_fname) )
        error("Failed to read from %s: %s\n", !strcmp("-",args->gt_fname)?"standard input":args->gt_fname,bcf_sr_strerror(args->files->errnum));

    args->qry_hdr = bcf_sr_get_header(args->files,0);
    if ( !bcf_hdr_nsamples(args->qry_hdr) ) error("No samples in %s?\n", args->qry_fname);
    if ( args->gt_fname )
    {
        args->gt_hdr = bcf_sr_get_header(args->files,1);
        if ( !bcf_hdr_nsamples(args->gt_hdr) ) error("No samples in %s?\n", args->gt_fname);
    }

    if ( args->gt_hdr && args->gt_filter_str ) args->gt_filter = filter_init(args->gt_hdr, args->gt_filter_str);
    if ( args->qry_hdr && args->qry_filter_str ) args->qry_filter = filter_init(args->qry_hdr, args->qry_filter_str);

    // Determine whether GT or PL will be used
    if ( args->qry_use_GT==-1 ) // not set by -u, qry uses PL by default
    {
        if ( bcf_hdr_id2int(args->qry_hdr,BCF_DT_ID,"PL")>=0 )
            args->qry_use_GT = 0;
        else if ( bcf_hdr_id2int(args->qry_hdr,BCF_DT_ID,"GT")>=0 )
            args->qry_use_GT = 1;
        else
            error("[E::%s] Neither PL nor GT tag is present in the header of %s\n", __func__, args->qry_fname);
    }
    else if ( args->qry_use_GT==1 )
    {
        if ( bcf_hdr_id2int(args->qry_hdr,BCF_DT_ID,"GT")<0 )
            error("[E::%s] The GT tag is not present in the header of %s\n", __func__, args->qry_fname);
    }
    else if ( bcf_hdr_id2int(args->qry_hdr,BCF_DT_ID,"PL")<0 )
        error("[E::%s] The PL tag is not present in the header of %s\n", __func__, args->qry_fname);

    if ( args->gt_hdr )
    {
        if ( args->gt_use_GT==-1 ) // not set by -u, gt uses GT by default
        {
            if ( bcf_hdr_id2int(args->gt_hdr,BCF_DT_ID,"GT")>=0 )
                args->gt_use_GT = 1;
            else if ( bcf_hdr_id2int(args->gt_hdr,BCF_DT_ID,"PL")>=0 )
                args->gt_use_GT = 0;
            else
                error("[E::%s] Neither PL nor GT tag is present in the header of %s\n", __func__, args->gt_fname);
        }
        else if ( args->gt_use_GT==1 )
        {
            if ( bcf_hdr_id2int(args->gt_hdr,BCF_DT_ID,"GT")<0 )
                error("[E::%s] The GT tag is not present in the header of %s\n", __func__, args->gt_fname);
        }
        else if ( bcf_hdr_id2int(args->gt_hdr,BCF_DT_ID,"PL")<0 )
            error("[E::%s] The PL tag is not present in the header of %s\n", __func__, args->gt_fname);
    }
    else
        args->gt_use_GT = args->qry_use_GT;

    // Prepare samples
    int i,j;
    args->nqry_smpl = bcf_hdr_nsamples(args->qry_hdr);
    if ( args->qry_samples )
    {
        init_samples(args->qry_samples, args->qry_samples_is_file, &args->qry_smpl, &args->nqry_smpl, args->qry_hdr, args->qry_fname);
    }
    if ( args->gt_samples )
    {
        init_samples(args->gt_samples, args->gt_samples_is_file, &args->gt_smpl, &args->ngt_smpl,
            args->gt_hdr ? args->gt_hdr : args->qry_hdr,
            args->gt_fname ? args->gt_fname : args->qry_fname);
    }
    else if ( args->pair_samples )
    {
        int npairs;
        char **tmp = hts_readlist(args->pair_samples, args->pair_samples_is_file, &npairs);
        if ( !tmp || !npairs ) error("Failed to parse %s\n", args->pair_samples);
        if ( !args->pair_samples_is_file && npairs%2 ) error("Expected even number of comma-delimited samples with -p\n");
        args->npairs = args->pair_samples_is_file ? npairs : npairs/2;
        args->pairs  = (pair_t*) calloc(args->npairs,sizeof(*args->pairs));
        if ( !args->pair_samples_is_file )
        {
            for (i=0; i<args->npairs; i++)
            {
                args->pairs[i].iqry = bcf_hdr_id2int(args->qry_hdr, BCF_DT_SAMPLE, tmp[2*i]);
                args->pairs[i].igt  = bcf_hdr_id2int(args->gt_hdr?args->gt_hdr:args->qry_hdr, BCF_DT_SAMPLE, tmp[2*i+1]);
                if ( args->pairs[i].iqry < 0 ) error("No such sample in %s: [%s]\n",args->qry_fname,tmp[2*i]);
                if ( args->pairs[i].igt  < 0 ) error("No such sample in %s: [%s]\n",args->gt_fname?args->gt_fname:args->qry_fname,tmp[2*i+1]);
                free(tmp[2*i]);
                free(tmp[2*i+1]);
            }
        }
        else
        {
            for (i=0; i<args->npairs; i++)
            {
                char *ptr = tmp[i];
                while ( *ptr && !isspace(*ptr) ) ptr++;
                if ( !*ptr ) error("Could not parse %s: %s\n",args->pair_samples,tmp[i]);
                *ptr = 0;
                args->pairs[i].iqry = bcf_hdr_id2int(args->qry_hdr, BCF_DT_SAMPLE, tmp[i]);
                if ( args->pairs[i].iqry < 0 ) error("No such sample in %s: [%s]\n",args->qry_fname,tmp[i]);
                ptr++;
                while ( *ptr && isspace(*ptr) ) ptr++;
                args->pairs[i].igt = bcf_hdr_id2int(args->gt_hdr?args->gt_hdr:args->qry_hdr, BCF_DT_SAMPLE, ptr);
                if ( args->pairs[i].igt < 0 ) error("No such sample in %s: [%s]\n",args->gt_fname?args->gt_fname:args->qry_fname,ptr);
                free(tmp[i]);
            }
        }
        free(tmp);
        qsort(args->pairs,args->npairs,sizeof(*args->pairs),cmp_pair);
    }
    else if ( args->gt_hdr )
        args->ngt_smpl = bcf_hdr_nsamples(args->gt_hdr);
    if ( !args->ngt_smpl )
    {
        args->ngt_smpl = args->nqry_smpl;
        args->gt_smpl  = args->qry_smpl;
        args->cross_check = 1;
    }

    // The data arrays
    if ( !args->npairs ) args->npairs = args->cross_check ? args->nqry_smpl*(args->nqry_smpl+1)/2 : args->ngt_smpl*args->nqry_smpl;
    if ( !args->pair_samples )
    {
        args->qry_dsg = (uint8_t*) malloc(args->nqry_smpl);
        args->gt_dsg  = args->cross_check ? args->qry_dsg : (uint8_t*) malloc(args->ngt_smpl);
    }
    if ( args->gt_err )
    {
        args->pdiff = (double*) calloc(args->npairs,sizeof(*args->pdiff));      // log probability of pair samples being the same
        args->qry_prob = (double*) malloc(3*args->nqry_smpl*sizeof(*args->qry_prob));
        args->gt_prob  = args->cross_check ? args->qry_prob : (double*) malloc(3*args->ngt_smpl*sizeof(*args->gt_prob));

        // Convert genotypes to genotype likelihoods given by -E, the probability of reading one allele incorrectly. In this
        // simple model we have:
        //     - probability of reading an allele incorrectly, eg. 0 as 1 or 1 as 0
        //         P(0|1) = P(1|0) = e
        //     - probability of genotype G={00,01,11} being correct given observed dosage {0,1,2} and the
        //       genotyping error probability e:
        //          P(00|0) = 1       P(00|1) = e       P(00|2) = e^2
        //          P(01|0) = e       P(01|1) = 1       P(01|2) = e
        //          P(11|0) = e^2     P(11|1) = e       P(11|2) = 1
        //
        // dsg2prob: the first index is bitmask of 8 possible dsg combinations (only 1<<0,1<<2,1<<3 are set, accessing
        // anything else indicated an error, this is just to reuse gt_to_dsg(); the second index are the corresponding
        // probabilities of 0/0, 0/1, and 1/1 genotypes
        //
        for (i=0; i<8; i++)
            for (j=0; j<3; j++)
                args->dsg2prob[i][j] = HUGE_VAL;
        double eprob = pow(10,-0.1*args->gt_err);      // convert from phred score to probability
        args->dsg2prob[1][0] = 0;               // P(00|0) = 1
        args->dsg2prob[1][1] = -log(eprob);     // P(01|0) = e
        args->dsg2prob[1][2] = -2*log(eprob);   // P(11|0) = e^2
        args->dsg2prob[2][0] = -log(eprob);     // P(00|1) = e
        args->dsg2prob[2][1] = 0;               // P(01|1) = 1
        args->dsg2prob[2][2] = -log(eprob);     // P(11|1) = e
        args->dsg2prob[4][0] = -2*log(eprob);   // P(00|2) = e^2
        args->dsg2prob[4][1] = -log(eprob);     // P(01|2) = e
        args->dsg2prob[4][2] = 0;               // P(11|2) = 1

        // lookup table to avoid exponentiation
        for (i=0; i<256; i++) args->pl2prob[i] = pow(10,-0.1*i);
    }
    else
        args->ndiff = (uint32_t*) calloc(args->npairs,sizeof(*args->ndiff));    // number of differing genotypes for each pair of samples
    args->ncnt   = (uint32_t*) calloc(args->npairs,sizeof(*args->ncnt));         // number of comparisons performed (non-missing data)
    if ( !args->ncnt ) error("Error: failed to allocate %.1f Mb\n", args->npairs*sizeof(*args->ncnt)/1e6);
    if ( args->calc_hwe_prob )
    {
        // prob of the observed sequence of matches given site AFs and HWE
        args->hwe_prob = (double*) calloc(args->npairs,sizeof(*args->hwe_prob));
        if ( !args->hwe_prob ) error("Error: failed to allocate %.1f Mb. Run with --no-HWE-prob to save some memory.\n", args->npairs*sizeof(*args->hwe_prob)/1e6);
        args->nmatch = (uint32_t*) calloc(args->npairs,sizeof(*args->ncnt));         // number of matches, used only with calc_hwe_prob
        if ( !args->nmatch ) error("Error: failed to allocate %.1f Mb.\n", args->npairs*sizeof(*args->ncnt)/1e6);
    }

    if ( args->distinctive_sites ) diff_sites_init(args);

    args->out_fh = bgzf_open(args->output_fname, args->output_type&FT_GZ ? "wg" : "wu");
    if ( args->out_fh == NULL )
        error("[%s] Error: cannot write to %s: %s\n", __func__,args->output_fname ? args->output_fname : "standard output", strerror(errno));

    print_header(args);
}

static void destroy_data(args_t *args)
{
    free(args->kstr.s);
    if ( args->gt_filter ) filter_destroy(args->gt_filter);
    if ( args->qry_filter ) filter_destroy(args->qry_filter);
    if ( args->gt_dsg!=args->qry_dsg ) free(args->gt_dsg);
    free(args->qry_dsg);
    if ( args->gt_prob!=args->qry_prob ) free(args->gt_prob);
    free(args->qry_prob);
    free(args->es_max_mem);
    if ( bgzf_close(args->out_fh)!=0 )  error("Error: close failed .. %s\n", args->output_fname?args->output_fname:"stdout");
    if ( args->distinctive_sites ) diff_sites_destroy(args);
    free(args->hwe_prob);
    free(args->cwd);
    free(args->qry_arr);
    if ( args->gt_hdr ) free(args->gt_arr);
    free(args->pdiff);
    free(args->ndiff);
    free(args->ncnt);
    free(args->nmatch);
    free(args->qry_smpl);
    if ( args->gt_smpl!=args->qry_smpl ) free(args->gt_smpl);
    free(args->pairs);
    bcf_sr_destroy(args->files);
}

static inline uint8_t gt_to_dsg(int32_t *ptr)
{
    if ( bcf_gt_is_missing(ptr[0]) || bcf_gt_is_missing(ptr[1]) || ptr[1]==bcf_int32_vector_end ) return 0;
    uint8_t dsg = (bcf_gt_allele(ptr[0])?1:0) + (bcf_gt_allele(ptr[1])?1:0);
    return 1<<dsg;
}
static inline uint8_t pl_to_dsg(int32_t *ptr)
{
    if ( ptr[0]==bcf_int32_missing || ptr[1]==bcf_int32_missing || ptr[2]==bcf_int32_missing ) return 0;
    if ( ptr[1]==bcf_int32_vector_end || ptr[2]==bcf_int32_vector_end ) return 0;
    int min_pl = ptr[0]<ptr[1] ? (ptr[0]<ptr[2]?ptr[0]:ptr[2]) : (ptr[1]<ptr[2]?ptr[1]:ptr[2]);
    uint8_t dsg = 0;
    if ( ptr[0]==min_pl ) dsg |= 1;
    if ( ptr[1]==min_pl ) dsg |= 2;
    if ( ptr[2]==min_pl ) dsg |= 4;
    return dsg;
}
static inline uint8_t gt_to_prob(args_t *args, int32_t *ptr, double *prob)
{
    uint8_t dsg = gt_to_dsg(ptr);
    if ( dsg )
    {
        prob[0] = args->dsg2prob[dsg][0];
        prob[1] = args->dsg2prob[dsg][1];
        prob[2] = args->dsg2prob[dsg][2];
    }
    return dsg;
}
static inline uint8_t pl_to_prob(args_t *args, int32_t *ptr, double *prob)
{
    uint8_t dsg = pl_to_dsg(ptr);
    if ( dsg )
    {
        prob[0] = (ptr[0]>=0 && ptr[0]<255) ? args->pl2prob[ptr[0]] : args->pl2prob[255];
        prob[1] = (ptr[1]>=0 && ptr[1]<255) ? args->pl2prob[ptr[1]] : args->pl2prob[255];
        prob[2] = (ptr[2]>=0 && ptr[2]<255) ? args->pl2prob[ptr[2]] : args->pl2prob[255];
        double sum = prob[0] + prob[1] + prob[2];
        prob[0] /= sum;
        prob[1] /= sum;
        prob[2] /= sum;
        prob[0] = -log(prob[0]);
        prob[1] = -log(prob[1]);
        prob[2] = -log(prob[2]);
    }
    return dsg;
}
static int set_data(args_t *args, bcf_hdr_t *hdr, bcf1_t *rec, int32_t **arr, int32_t *narr, int *narr1, int *use_GT)
{
    static int warn_dip_GT = 1;
    static int warn_dip_PL = 1;
    int i;
    for (i=0; i<2; i++)
    {
        if ( *use_GT )
        {
            int ret = bcf_get_genotypes(hdr,rec,arr,narr);
            if ( ret < 0 )
            {
                if ( !i ) { *use_GT = 0; continue; }
                args->nskip_no_data++;
                return -1;
            }
            if ( ret != 2*bcf_hdr_nsamples(hdr) )
            {
                if ( warn_dip_GT )
                {
                    fprintf(stderr,"INFO: skipping %s:%"PRIhts_pos", only diploid FORMAT/GT fields supported. (This is printed only once.)\n", bcf_seqname(hdr,rec),rec->pos+1);
                    warn_dip_GT = 0;
                }
                args->nskip_dip_GT++;
                return -1;
            }
            *narr1 = 2;
            return 0;
        }

        int ret = bcf_get_format_int32(hdr,rec,"PL",arr,narr);
        if ( ret < 0 )
        {
            if ( !i ) { *use_GT = 1; continue; }
            args->nskip_no_data++;
            return -1;
        }
        if ( ret != 3*bcf_hdr_nsamples(hdr) )
        {
            if ( warn_dip_PL )
            {
                fprintf(stderr,"INFO: skipping %s:%"PRIhts_pos", only diploid FORMAT/PL fields supported. (This is printed only once.)\n", bcf_seqname(hdr,rec),rec->pos+1);
                warn_dip_PL = 0;
            }
            args->nskip_dip_PL++;
            return -1;
        }
        *narr1 = 3;
        return 0;
    }
    return -1;  // should never reach
}
static void process_line(args_t *args)
{
    int i,j,k, nqry1, ngt1, ret;

    bcf1_t *gt_rec = NULL, *qry_rec = bcf_sr_get_line(args->files,0);   // the query file
    if ( args->qry_filter )
    {
        int pass = filter_test(args->qry_filter, qry_rec, NULL);
        if ( args->qry_filter_logic==FLT_EXCLUDE ) pass = pass ? 0 : 1;
        if ( !pass ) { args->nskip_filter++; return; }
    }

    int qry_use_GT = args->qry_use_GT;
    int gt_use_GT  = args->gt_use_GT;

    ret = set_data(args, args->qry_hdr, qry_rec, &args->qry_arr, &args->nqry_arr, &nqry1, &qry_use_GT);
    if ( ret<0 ) return;

    if ( args->gt_hdr )
    {
        gt_rec = bcf_sr_get_line(args->files,1);
        if ( args->gt_filter )
        {
            int pass = filter_test(args->gt_filter, gt_rec, NULL);
            if ( args->gt_filter_logic==FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( !pass ) { args->nskip_filter++; return; }
        }
        ret = set_data(args, args->gt_hdr, gt_rec, &args->gt_arr, &args->ngt_arr, &ngt1, &gt_use_GT);
        if ( ret<0 ) return;
    }
    else
    {
        ngt1 = nqry1;
        args->gt_arr = args->qry_arr;
    }

    // stats: number of compared sites, and used tags
    args->ncmp++;
    args->nused[qry_use_GT][gt_use_GT]++;

    double hwe_dsg[8];
    if ( args->calc_hwe_prob )
    {
        int ac[2];
        if ( args->gt_hdr )
        {
            if ( bcf_calc_ac(args->gt_hdr, gt_rec, ac, BCF_UN_INFO|BCF_UN_FMT)!=1 ) error("todo: bcf_calc_ac() failed\n");
        }
        else if ( bcf_calc_ac(args->qry_hdr, qry_rec, ac, BCF_UN_INFO|BCF_UN_FMT)!=1 ) error("todo: bcf_calc_ac() failed\n");

        // Calculate HWE probability for each possible qry+gt dosage combination. The alternate allele dosage
        // values returned by gt_to_prob() below are 0,1,2,4 (0=missing, 1<<0, 1<<1, 1<<2). We consider only
        // biallelic sites, therefore we work with eight genotype combinations.
        //
        // The array hwe_dsg is accessed with hwe_dsg[qry_dsg & gt_dsg] and is constructed to account for PL uncertainty
        // when we encounter less informative PL, such as PL=0,0,10, where multiple dosage values are equally
        // likely. If we allowed complete uncertainty (PL=0,0,0), we'd have up to eight possible genotype
        // mask combinations: from e.g. 0=(gt_dsg=1<<0 & qry_dsg=1<<1) to 7=(gt_dsg=1<<0|1<<1|1<<2 & qry_dsg=1<<0|1<<1|1<<2).
        // Note the extreme case of 1|2|4 is skipped, see pl_to_dsg().
        //
        // When the dosage is uncertain, we take the minimum of their corresponding HWE value, for example
        //      hwe[0] = 0
        //      hwe[1] = (1-AF)**2
        //      hwe[2] = 2*AF*(1-AF)
        //      hwe[4] = AF**2
        //      hwe[3] = min{hwe[1],hwe[2]}

        double hwe[3];  // while hwe_dsg iterates over dsg bitmasks (0..7), hwe iterates over dsg (0,1,2)
        double af = ac[0]+ac[1] ? (double)ac[1]/(ac[0]+ac[1]) : 1e-6;
        hwe[0] = -log((1-af)*(1-af));
        hwe[1] = -log(2*af*(1-af));
        hwe[2] = -log(af*af);
        hwe_dsg[0] = 0;
        for (i=1; i<8; i++)
        {
            hwe_dsg[i] = HUGE_VAL;
            for (k=0; k<3; k++)
            {
                if ( ((1<<k)&i) && hwe_dsg[i] > hwe[k] ) hwe_dsg[i] = hwe[k];
            }
        }
    }

    // The sample pairs were given explicitly via -p/-P options
    if ( args->pairs )
    {
        if ( !args->gt_err )
        {
            int ndiff = 0;
            if ( args->kbs_diff ) diff_sites_reset(args);

            for (i=0; i<args->npairs; i++)
            {
                int32_t *ptr;
                uint8_t qry_dsg, gt_dsg;

                ptr = args->gt_arr + args->pairs[i].igt*ngt1;
                gt_dsg = gt_use_GT ? gt_to_dsg(ptr) : pl_to_dsg(ptr);
                if ( !gt_dsg ) continue;                        // missing value
                if ( args->hom_only && !(gt_dsg&5) ) continue;  // not a hom

                ptr = args->qry_arr + args->pairs[i].iqry*nqry1;
                qry_dsg = qry_use_GT ? gt_to_dsg(ptr) : pl_to_dsg(ptr);
                if ( !qry_dsg ) continue;                       // missing value

                int match = qry_dsg & gt_dsg;
                if ( !match )
                {
                    args->ndiff[i]++;
                    if ( args->kbs_diff ) { ndiff++; kbs_insert(args->kbs_diff, i); }
                }
                else if ( args->calc_hwe_prob )
                {
                    args->hwe_prob[i] += hwe_dsg[match];
                    args->nmatch[i]++;
                }
                args->ncnt[i]++;
            }

            if ( ndiff ) diff_sites_push(args, ndiff, qry_rec->rid, qry_rec->pos);
        }
        else    // gt_err set
        {
            for (i=0; i<args->npairs; i++)
            {
                int32_t *ptr;
                double qry_prob[3], gt_prob[3];
                uint8_t qry_dsg, gt_dsg;

                ptr = args->gt_arr + args->pairs[i].igt*ngt1;
                gt_dsg = gt_use_GT ? gt_to_prob(args,ptr,gt_prob) : pl_to_prob(args,ptr,gt_prob);
                if ( !gt_dsg ) continue;                        // missing value
                if ( args->hom_only && !(gt_dsg&5) ) continue;  // not a hom

                ptr = args->qry_arr + args->pairs[i].iqry*nqry1;
                qry_dsg = qry_use_GT ? gt_to_prob(args,ptr,qry_prob) : pl_to_prob(args,ptr,qry_prob);
                if ( !qry_dsg ) continue;                       // missing value

                double min = qry_prob[0] + gt_prob[0];
                qry_prob[1] += gt_prob[1];
                if ( min > qry_prob[1] ) min = qry_prob[1];
                qry_prob[2] += gt_prob[2];
                if ( min > qry_prob[2] ) min = qry_prob[2];
                args->pdiff[i] += min;

                if ( args->calc_hwe_prob )
                {
                    int match = qry_dsg & gt_dsg;
                    args->hwe_prob[i] += hwe_dsg[match];
                    if ( match ) args->nmatch[i]++;
                }
                args->ncnt[i]++;
            }
        }
        return;
    }

    int idx=0;
    if ( !args->gt_err )
    {
        for (i=0; i<args->nqry_smpl; i++)
        {
            int iqry = args->qry_smpl ? args->qry_smpl[i] : i;
            int32_t *ptr = args->qry_arr + nqry1*iqry;
            args->qry_dsg[i] = qry_use_GT ? gt_to_dsg(ptr) : pl_to_dsg(ptr);
        }
        if ( !args->cross_check )   // in this case gt_dsg points to qry_dsg
        {
            for (i=0; i<args->ngt_smpl; i++)
            {
                int igt = args->gt_smpl ? args->gt_smpl[i] : i;
                int32_t *ptr = args->gt_arr + ngt1*igt;
                args->gt_dsg[i] = gt_use_GT ? gt_to_dsg(ptr) : pl_to_dsg(ptr);
                if ( args->hom_only && !(args->gt_dsg[i]&5) ) args->gt_dsg[i] = 0;      // not a hom, set to a missing value
            }
        }
        for (i=0; i<args->nqry_smpl; i++)
        {
            int ngt = args->cross_check ? i : args->ngt_smpl;       // two files or a sub-diagonal cross-check mode?
            if ( !args->qry_dsg[i] ) { idx += ngt; continue; }      // missing value
            for (j=0; j<ngt; j++)
            {
                if ( !args->gt_dsg[j] ) { idx++; continue; }        // missing value
                int match = args->qry_dsg[i] & args->gt_dsg[j];
                if ( !match ) args->ndiff[idx]++;
                else if ( args->calc_hwe_prob )
                {
                    args->hwe_prob[idx] += hwe_dsg[match];
                    args->nmatch[idx]++;
                }
                args->ncnt[idx]++;
                idx++;
            }
        }
    }
    else    // gt_err set
    {
        for (i=0; i<args->nqry_smpl; i++)
        {
            int iqry = args->qry_smpl ? args->qry_smpl[i] : i;
            int32_t *ptr = args->qry_arr + nqry1*iqry;
            args->qry_dsg[i] = qry_use_GT ? gt_to_prob(args,ptr,args->qry_prob+i*3) : pl_to_prob(args,ptr,args->qry_prob+i*3);
        }
        if ( !args->cross_check )   // in this case gt_dsg points to qry_dsg
        {
            for (i=0; i<args->ngt_smpl; i++)
            {
                int igt = args->gt_smpl ? args->gt_smpl[i] : i;
                int32_t *ptr = args->gt_arr + ngt1*igt;
                args->gt_dsg[i] = gt_use_GT ? gt_to_prob(args,ptr,args->gt_prob+i*3) : pl_to_prob(args,ptr,args->gt_prob+i*3);
                if ( args->hom_only && !(args->gt_dsg[i]&5) ) args->gt_dsg[i] = 0;      // not a hom, set to a missing value
            }
        }
        for (i=0; i<args->nqry_smpl; i++)
        {
            int ngt = args->cross_check ? i : args->ngt_smpl;       // two files or a sub-diagonal cross-check mode?
            if ( !args->qry_dsg[i] ) { idx += ngt; continue; }      // missing value
            for (j=0; j<ngt; j++)
            {
                if ( !args->gt_dsg[j] ) { idx++; continue; }        // missing value

                double min = args->qry_prob[i*3] + args->gt_prob[j*3];
                if ( min > args->qry_prob[i*3+1] + args->gt_prob[j*3+1] ) min = args->qry_prob[i*3+1] + args->gt_prob[j*3+1];
                if ( min > args->qry_prob[i*3+2] + args->gt_prob[j*3+2] ) min = args->qry_prob[i*3+2] + args->gt_prob[j*3+2];
                args->pdiff[idx] += min;

                if ( args->calc_hwe_prob )
                {
                    int match = args->qry_dsg[i] & args->gt_dsg[j];
                    args->hwe_prob[idx] += hwe_dsg[match];
                    if ( match ) args->nmatch[idx]++;
                }
                args->ncnt[idx]++;
                idx++;
            }
        }
    }
}


typedef struct
{
    int ism, idx;
    double val;
}
idbl_t;
static int cmp_idbl(const void *_a, const void *_b)
{
    idbl_t *a = (idbl_t*)_a;
    idbl_t *b = (idbl_t*)_b;
    if ( a->val < b->val ) return -1;
    if ( a->val > b->val ) return 1;
    return 0;
}
static void report_distinctive_sites(args_t *args)
{
    extsort_sort(args->es);

    args->kstr.l = 0;
    ksprintf(&args->kstr,"# DS, distinctive sites:\n");
    ksprintf(&args->kstr,"#     - chromosome\n");
    ksprintf(&args->kstr,"#     - position\n");
    ksprintf(&args->kstr,"#     - cumulative number of pairs distinguished by this block\n");
    ksprintf(&args->kstr,"#     - block id\n");
    ksprintf(&args->kstr,"#DS\t[2]Chromosome\t[3]Position\t[4]Cumulative number of distinct pairs\t[5]Block id\n");

    kbitset_t *kbs_blk = kbs_init(args->npairs);
    kbitset_iter_t itr;
    int i,ndiff,rid,pos,ndiff_tot = 0, iblock = 0;
    int ndiff_min = args->distinctive_sites <= args->npairs ? args->distinctive_sites : args->npairs;
    while ( diff_sites_shift(args,&ndiff,&rid,&pos) )
    {
        int ndiff_new = 0, ndiff_dbg = 0;
        kbs_start(&itr);
        while ( (i=kbs_next(args->kbs_diff, &itr))>=0 )
        {
            ndiff_dbg++;
            if ( kbs_exists(kbs_blk,i) ) continue;   // already set
            kbs_insert(kbs_blk,i);
            ndiff_new++;
        }
        if ( ndiff_dbg!=ndiff ) error("Corrupted data, fixme: %d vs %d\n",ndiff_dbg,ndiff);
        if ( !ndiff_new ) continue;     // no new pair distinguished by this site
        ndiff_tot += ndiff_new;
        args->kstr.l = 0;
        ksprintf(&args->kstr,"DS\t%s\t%d\t%d\t%d\n",bcf_hdr_id2name(args->qry_hdr,rid),pos+1,ndiff_tot,iblock);
        if ( bgzf_write(args->out_fh, args->kstr.s, args->kstr.l)!=args->kstr.l ) error("Failed to write to %s\n", args->output_fname);
        if ( ndiff_tot < ndiff_min ) continue;   // fewer than the requested number of pairs can be distinguished at this point
        iblock++;
        ndiff_tot = 0;
        kbs_clear(kbs_blk);
    }
    kbs_destroy(kbs_blk);
}
static void report(args_t *args)
{
    args->kstr.l = 0;
    ksprintf(&args->kstr,"INFO\tsites-compared\t%u\n",args->ncmp);
    ksprintf(&args->kstr,"INFO\tsites-skipped-no-match\t%u\n",args->nskip_no_match);
    ksprintf(&args->kstr,"INFO\tsites-skipped-multiallelic\t%u\n",args->nskip_not_ba);
    ksprintf(&args->kstr,"INFO\tsites-skipped-monoallelic\t%u\n",args->nskip_mono);
    ksprintf(&args->kstr,"INFO\tsites-skipped-no-data\t%u\n",args->nskip_no_data);
    ksprintf(&args->kstr,"INFO\tsites-skipped-GT-not-diploid\t%u\n",args->nskip_dip_GT);
    ksprintf(&args->kstr,"INFO\tsites-skipped-PL-not-diploid\t%u\n",args->nskip_dip_PL);
    ksprintf(&args->kstr,"INFO\tsites-skipped-filtering-expression\t%u\n",args->nskip_filter);
    ksprintf(&args->kstr,"INFO\tsites-used-PL-vs-PL\t%u\n",args->nused[0][0]);
    ksprintf(&args->kstr,"INFO\tsites-used-PL-vs-GT\t%u\n",args->nused[0][1]);
    ksprintf(&args->kstr,"INFO\tsites-used-GT-vs-PL\t%u\n",args->nused[1][0]);
    ksprintf(&args->kstr,"INFO\tsites-used-GT-vs-GT\t%u\n",args->nused[1][1]);
    ksprintf(&args->kstr,"# DCv2, discordance version 2:\n");
    ksprintf(&args->kstr,"#     - Query sample\n");
    ksprintf(&args->kstr,"#     - Genotyped sample\n");
    ksprintf(&args->kstr,"#     - Discordance, given either as an abstract score or number of mismatches, see the options -E/-u\n"
                         "#       in man page for details. Note that samples with high missingness have fewer sites compared,\n"
                         "#       which results in lower overall discordance. Therefore it is advisable to use the average score\n"
                         "#       per site rather than the absolute value, i.e. divide the value by the number of sites compared\n"
                         "#       (smaller value = better match)\n");
    ksprintf(&args->kstr,"#     - Average negative log of HWE probability at matching sites, attempts to quantify the following\n"
                         "#       intuition: rare genotype matches are more informative than common genotype matches, hence two\n"
                         "#       samples with similar discordance can be further stratified by the HWE score (bigger value = better\n"
                         "#       match, the observed concordance was less likely to occur by chance)\n");
    ksprintf(&args->kstr,"#     - Number of sites compared for this pair of samples (bigger = more informative)\n");
    ksprintf(&args->kstr,"#     - Number of matching genotypes\n");
    ksprintf(&args->kstr,"#DCv2\t[2]Query Sample\t[3]Genotyped Sample\t[4]Discordance\t[5]Average -log P(HWE)\t[6]Number of sites compared\t[7]Number of matching genotypes\n");
    if ( bgzf_write(args->out_fh, args->kstr.s, args->kstr.l)!=args->kstr.l ) error("Failed to write to %s\n", args->output_fname);

    int trim = args->ntop;
    if ( !args->pairs )
    {
        if ( !args->ngt_smpl && args->nqry_smpl <= args->ntop ) trim = 0;
        if ( args->ngt_smpl && args->ngt_smpl <= args->ntop  ) trim = 0;
    }

    if ( args->pairs )
    {
        int i;
        for (i=0; i<args->npairs; i++)
        {
            args->kstr.l = 0;
            int iqry = args->pairs[i].iqry;
            int igt  = args->pairs[i].igt;
            if ( args->ndiff )
            {
                ksprintf(&args->kstr,"DCv2\t%s\t%s\t%u\t%e\t%u\t%u\n",
                        args->qry_hdr->samples[iqry],
                        args->gt_hdr?args->gt_hdr->samples[igt]:args->qry_hdr->samples[igt],
                        args->ndiff[i],
                        (args->calc_hwe_prob && args->nmatch[i]) ? args->hwe_prob[i]/args->nmatch[i] : 0,
                        args->ncnt[i],
                        args->nmatch[i]);
            }
            else
            {
                ksprintf(&args->kstr,"DCv2\t%s\t%s\t%e\t%e\t%u\t%u\n",
                        args->qry_hdr->samples[iqry],
                        args->gt_hdr?args->gt_hdr->samples[igt]:args->qry_hdr->samples[igt],
                        args->pdiff[i],
                        (args->calc_hwe_prob && args->nmatch[i]) ? args->hwe_prob[i]/args->nmatch[i] : 0,
                        args->ncnt[i],
                        args->nmatch[i]);
            }
            if ( bgzf_write(args->out_fh, args->kstr.s, args->kstr.l)!=args->kstr.l ) error("Failed to write to %s\n", args->output_fname);
        }
    }
    else if ( !trim )
    {
        int i,j,idx=0;
        for (i=0; i<args->nqry_smpl; i++)
        {
            int iqry = args->qry_smpl ? args->qry_smpl[i] : i;
            int ngt  = args->cross_check ? i : args->ngt_smpl;
            for (j=0; j<ngt; j++)
            {
                args->kstr.l = 0;
                int igt = args->gt_smpl ? args->gt_smpl[j] : j;
                if ( args->ndiff )
                {
                    ksprintf(&args->kstr,"DCv2\t%s\t%s\t%u\t%e\t%u\t%u\n",
                            args->qry_hdr->samples[iqry],
                            args->gt_hdr?args->gt_hdr->samples[igt]:args->qry_hdr->samples[igt],
                            args->ndiff[idx],
                            (args->calc_hwe_prob && args->nmatch[idx]) ? args->hwe_prob[idx]/args->nmatch[idx] : 0,
                            args->ncnt[idx],
                            args->calc_hwe_prob ? args->nmatch[idx] : 0);
                }
                else
                {
                    ksprintf(&args->kstr,"DCv2\t%s\t%s\t%e\t%e\t%u\t%u\n",
                            args->qry_hdr->samples[iqry],
                            args->gt_hdr?args->gt_hdr->samples[igt]:args->qry_hdr->samples[igt],
                            args->pdiff[idx],
                            (args->calc_hwe_prob && args->nmatch[idx]) ? args->hwe_prob[idx]/args->nmatch[idx] : 0,
                            args->ncnt[idx],
                            args->calc_hwe_prob ? args->nmatch[idx] : 0);
                }
                if ( bgzf_write(args->out_fh, args->kstr.s, args->kstr.l)!=args->kstr.l ) error("Failed to write to %s\n", args->output_fname);
                idx++;
            }
        }
    }
    else if ( !args->cross_check )
    {
        idbl_t *arr = (idbl_t*)malloc(sizeof(*arr)*args->ngt_smpl);
        int i,j;
        for (i=0; i<args->nqry_smpl; i++)
        {
            int idx  = i*args->ngt_smpl;
            for (j=0; j<args->ngt_smpl; j++)
            {
                if ( args->sort_by_hwe )
                    arr[j].val = args->nmatch[idx] ? -args->hwe_prob[idx]/args->nmatch[idx] : 0;  // -args->hwe_prob[idx];
                else if ( args->ndiff )
                    arr[j].val = args->ncnt[idx] ? (double)args->ndiff[idx]/args->ncnt[idx] : 0;
                else
                    arr[j].val = args->ncnt[idx] ? args->pdiff[idx]/args->ncnt[idx] : 0;
                arr[j].ism = j;
                arr[j].idx = idx;
                idx++;
            }
            qsort(arr, args->ngt_smpl, sizeof(*arr), cmp_idbl);
            int iqry = args->qry_smpl ? args->qry_smpl[i] : i;
            for (j=0; j<args->ntop; j++)
            {
                args->kstr.l = 0;
                int idx = arr[j].idx;
                int igt = args->gt_smpl ? args->gt_smpl[arr[j].ism] : arr[j].ism;
                if ( args->ndiff )
                {
                    ksprintf(&args->kstr,"DCv2\t%s\t%s\t%u\t%e\t%u\t%u\n",
                            args->qry_hdr->samples[iqry],
                            args->gt_hdr?args->gt_hdr->samples[igt]:args->qry_hdr->samples[igt],
                            args->ndiff[idx],
                            (args->calc_hwe_prob && args->nmatch[idx]) ? args->hwe_prob[idx]/args->nmatch[idx] : 0,
                            args->ncnt[idx],
                            args->calc_hwe_prob ? args->nmatch[idx] : 0);
                }
                else
                {
                    ksprintf(&args->kstr,"DCv2\t%s\t%s\t%e\t%e\t%u\t%u\n",
                            args->qry_hdr->samples[iqry],
                            args->gt_hdr?args->gt_hdr->samples[igt]:args->qry_hdr->samples[igt],
                            args->pdiff[idx],
                            (args->calc_hwe_prob && args->nmatch[idx]) ? args->hwe_prob[idx]/args->nmatch[idx] : 0,
                            args->ncnt[idx],
                            args->calc_hwe_prob ? args->nmatch[idx] : 0);
                }
                if ( bgzf_write(args->out_fh, args->kstr.s, args->kstr.l)!=args->kstr.l ) error("Failed to write to %s\n", args->output_fname);
            }
        }
        free(arr);
    }
    else
    {
        int narr = args->nqry_smpl-1;
        idbl_t *arr = (idbl_t*)malloc(sizeof(*arr)*narr);
        int i,j,k,idx;
        for (i=0; i<args->nqry_smpl; i++)
        {
            k = 0, idx = i*(i-1)/2;
            for (j=0; j<i; j++)
            {
                if ( args->sort_by_hwe )
                    arr[k].val = args->nmatch[idx] ? -args->hwe_prob[idx]/args->nmatch[idx] : 0;
                else if ( args->ndiff )
                    arr[k].val = args->ncnt[idx] ? (double)args->ndiff[idx]/args->ncnt[idx] : 0;
                else
                    arr[k].val = args->ncnt[idx] ? args->pdiff[idx]/args->ncnt[idx] : 0;
                arr[k].ism = j;
                arr[k].idx = idx;
                idx++;
                k++;
            }
            for (; j<narr; j++)
            {
                idx = j*(j+1)/2 + i;
                if ( args->sort_by_hwe )
                    arr[k].val = args->nmatch[idx] ? -args->hwe_prob[idx]/args->nmatch[idx] : 0;
                else if ( args->ndiff )
                    arr[k].val = args->ncnt[idx] ? (double)args->ndiff[idx]/args->ncnt[idx] : 0;
                else
                    arr[k].val = args->ncnt[idx] ? args->pdiff[idx]/args->ncnt[idx] : 0;
                arr[k].ism = j + 1;
                arr[k].idx = idx;
                k++;
            }
            qsort(arr, narr, sizeof(*arr), cmp_idbl);
            int iqry = args->qry_smpl ? args->qry_smpl[i] : i;
            for (j=0; j<args->ntop; j++)
            {
                args->kstr.l = 0;
                if ( i <= arr[j].ism ) continue;
                int idx = arr[j].idx;
                int igt = args->qry_smpl ? args->qry_smpl[arr[j].ism] : arr[j].ism;
                if ( args->ndiff )
                {
                    ksprintf(&args->kstr,"DCv2\t%s\t%s\t%u\t%e\t%u\t%u\n",
                            args->qry_hdr->samples[iqry],
                            args->qry_hdr->samples[igt],
                            args->ndiff[idx],
                            (args->calc_hwe_prob && args->nmatch[idx]) ? args->hwe_prob[idx]/args->nmatch[idx] : 0,
                            args->ncnt[idx],
                            args->calc_hwe_prob ? args->nmatch[idx] : 0);
                }
                else
                {
                    ksprintf(&args->kstr,"DCv2\t%s\t%s\t%e\t%e\t%u\t%u\n",
                            args->qry_hdr->samples[iqry],
                            args->qry_hdr->samples[igt],
                            args->pdiff[idx],
                            (args->calc_hwe_prob && args->nmatch[idx]) ? args->hwe_prob[idx]/args->nmatch[idx] : 0,
                            args->ncnt[idx],
                            args->calc_hwe_prob ? args->nmatch[idx] : 0);
                }
                if ( bgzf_write(args->out_fh, args->kstr.s, args->kstr.l)!=args->kstr.l ) error("Failed to write to %s\n", args->output_fname);
            }
        }
        free(arr);
    }
}

static int is_input_okay(args_t *args, int nmatch)
{
    int i;
    const char *msg;
    bcf_hdr_t *hdr;
    bcf1_t *rec;
    if ( args->gt_hdr && nmatch!=2 )
    {
        if ( args->nskip_no_match++ ) return 0;
        for (i=0; i<2; i++)
        {
            rec = bcf_sr_get_line(args->files,i);
            if ( rec ) break;
        }
        hdr = bcf_sr_get_header(args->files,i);
        fprintf(stderr,"INFO: skipping %s:%"PRIhts_pos", no record with matching POS+ALT. (This is printed only once.)\n",
                bcf_seqname(hdr,rec),rec->pos+1);
        return 0;
    }
    for (i=0; i<2; i++)
    {
        hdr = bcf_sr_get_header(args->files,i);
        rec = bcf_sr_get_line(args->files,i);
        if ( rec->n_allele>2 )
        {
            if ( args->nskip_not_ba++ ) return 0;
            msg = "not a biallelic site, run `bcftools norm -m -` first";
            goto not_okay;
        }
        if ( bcf_get_variant_types(rec)==VCF_REF )
        {
            if ( args->nskip_mono++ ) return 0;
            msg = "monoallelic site";
            goto not_okay;
        }
        if ( !args->gt_hdr ) break;
    }
    return 1;

not_okay:
    fprintf(stderr,"INFO: skipping %s:%"PRIhts_pos", %s. (This is printed only once.)\n",
        bcf_seqname(hdr,rec),rec->pos+1,msg);
    return 0;
}

static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Check sample identity. With no -g BCF given, multi-sample cross-check is performed.\n");
    fprintf(stderr, "Usage:   bcftools gtcheck [options] [-g <genotypes.vcf.gz>] <query.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    //fprintf(stderr, "    -a, --all-sites                  Output comparison for all sites\n");
    //fprintf(stderr, "    -c, --cluster MIN,MAX            Min inter- and max intra-sample error [0.23,-0.3]\n");
    fprintf(stderr, "        --distinctive-sites            Find sites that can distinguish between at least NUM sample pairs.\n");
    fprintf(stderr, "                  NUM[,MEM[,TMP]]          If the number is smaller or equal to 1, it is interpreted as the fraction of pairs.\n");
    fprintf(stderr, "                                           The optional MEM string sets the maximum memory used for in-memory sorting [500M]\n");
#ifdef _WIN32
    fprintf(stderr, "                                           and TMP is a prefix of temporary files used by external sorting [/bcftools.XXXXXX]\n");
#else
    fprintf(stderr, "                                           and TMP is a prefix of temporary files used by external sorting [/tmp/bcftools.XXXXXX]\n");
#endif
    fprintf(stderr, "        --dry-run                      Stop after first record to estimate required time\n");
    fprintf(stderr, "    -E, --error-probability INT        Phred-scaled probability of genotyping error, 0 for faster but less accurate results [40]\n");
    fprintf(stderr, "    -e, --exclude [qry|gt]:EXPR        Exclude sites for which the expression is true\n");
    fprintf(stderr, "    -g, --genotypes FILE               Genotypes to compare against\n");
    fprintf(stderr, "    -H, --homs-only                    Homozygous genotypes only, useful with low coverage data (requires -g)\n");
    fprintf(stderr, "    -i, --include [qry|gt]:EXPR        Include sites for which the expression is true\n");
    fprintf(stderr, "        --n-matches INT                Print only top INT matches for each sample (sorted by average score), 0 for unlimited.\n");
    fprintf(stderr, "                                           Use negative value to sort by HWE probability rather than by discordance [0]\n");
    fprintf(stderr, "        --no-HWE-prob                  Disable calculation of HWE probability\n");
    fprintf(stderr, "    -o, --output FILE                  Write output to a file [standard output]\n");
    fprintf(stderr, "    -O, --output-type t|z              t: plain tab-delimited text output, z: compressed [t]\n");
    fprintf(stderr, "    -p, --pairs LIST                   Comma-separated sample pairs to compare (qry,gt[,qry,gt..] with -g or qry,qry[,qry,qry..] w/o)\n");
    fprintf(stderr, "    -P, --pairs-file FILE              File with tab-delimited sample pairs to compare (qry,gt with -g or qry,qry w/o)\n");
    fprintf(stderr, "    -r, --regions REGION               Restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file FILE            Restrict to regions listed in a file\n");
    fprintf(stderr, "        --regions-overlap 0|1|2        Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n");
    fprintf(stderr, "    -s, --samples [qry|gt]:LIST        List of query or -g samples, \"-\" to select all samples (by default all samples are compared)\n");
    fprintf(stderr, "    -S, --samples-file [qry|gt]:FILE   File with the query or -g samples to compare\n");
    fprintf(stderr, "    -t, --targets REGION               Similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "    -T, --targets-file FILE            Similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "        --targets-overlap 0|1|2        Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n");
    fprintf(stderr, "    -u, --use TAG1[,TAG2]              Which tag to use in the query file (TAG1) and the -g file (TAG2) [PL,GT]\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "   # Check discordance of all samples from B against all samples in A\n");
    fprintf(stderr, "   bcftools gtcheck -g A.bcf B.bcf\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   # Limit comparisons to the given list of samples\n");
    fprintf(stderr, "   bcftools gtcheck -s gt:a1,a2,a3 -s qry:b1,b2 -g A.bcf B.bcf\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   # Compare only two pairs a1,b1 and a1,b2\n");
    fprintf(stderr, "   bcftools gtcheck -p a1,b1,a1,b2 -g A.bcf B.bcf\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfgtcheck(int argc, char *argv[])
{
    int c;
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv; set_cwd(args);
    args->qry_use_GT = -1;
    args->gt_use_GT  = -1;
    args->calc_hwe_prob = 1;
    args->gt_err = 40;
    args->regions_overlap = 1;
    args->targets_overlap = 0;
    args->output_fname = "-";

    // external sort for --distinctive-sites
#ifdef _WIN32
    args->es_tmp_prefix = NULL;
#else
    args->es_tmp_prefix = "/tmp/bcftools-gtcheck";
#endif
    args->es_max_mem = strdup("500M");

    // In simulated sample swaps the minimum error was 0.3 and maximum intra-sample error was 0.23
    //    - min_inter: pairs with smaller err value will be considered identical
    //    - max_intra: pairs with err value bigger than abs(max_intra_err) will be considered
    //                  different. If negative, the cutoff may be heuristically lowered
    args->min_inter_err =  0.23;
    args->max_intra_err = -0.3;

    static struct option loptions[] =
    {
        {"error-probability",1,0,'E'},  // note this used to be 'e', but can easily auto-detect to assure backward compatibility
        {"exclude",required_argument,0,'e'},
        {"include",required_argument,0,'i'},
        {"output",required_argument,0,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"use",1,0,'u'},
        {"cluster",1,0,'c'},
        {"GTs-only",1,0,'G'},
        {"all-sites",0,0,'a'},
        {"homs-only",0,0,'H'},
        {"help",0,0,'h'},
        {"genotypes",1,0,'g'},
        {"plot",1,0,'p'},
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {"n-matches",1,0,2},
        {"no-HWE-prob",0,0,3},
        {"target-sample",1,0,4},
        {"dry-run",0,0,5},
        {"distinctive-sites",1,0,6},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"regions-overlap",required_argument,NULL,7},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"targets-overlap",required_argument,NULL,8},
        {"pairs",1,0,'p'},
        {"pairs-file",1,0,'P'},
        {0,0,0,0}
    };
    char *tmp;
    while ((c = getopt_long(argc, argv, "hg:p:s:S:p:P:Hr:R:at:T:G:c:u:e:E:i:o:O:",loptions,NULL)) >= 0) {
        switch (c) {
            case 'o': args->output_fname = optarg; break;
            case 'O':
                switch (optarg[0]) {
                    case 't': args->output_type = FT_TAB_TEXT; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
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
                if ( !strncasecmp("gt:",optarg,3) )
                {
                    if ( args->gt_filter_str ) error("Error: only one -i or -e expression can be given for gt:, and they cannot be combined\n");
                    args->gt_filter_str = optarg;
                    args->gt_filter_logic |= FLT_EXCLUDE;
                }
                else if ( !strncasecmp("qry:",optarg,4) )
                {
                    if ( args->qry_filter_str ) error("Error: only one -i or -e expression can be given for qry:, and they cannot be combined\n");
                    args->qry_filter_str = optarg;
                    args->qry_filter_logic |= FLT_EXCLUDE;
                }
                else
                {
                    // this could be the old -e, --error-probability option
                    args->gt_err = strtol(optarg,&tmp,10);
                    if ( !tmp || *tmp )
                    {
                        // it is not
                        args->gt_filter_str  = optarg;
                        args->qry_filter_str = optarg;
                        args->gt_filter_logic  |= FLT_EXCLUDE;
                        args->qry_filter_logic |= FLT_EXCLUDE;
                    }
                    else
                    {
                        fprintf(stderr,"[warning] auto-detected the old format --error-probability option, please switch from -e to -E.\n");
                    }
                }
                break;
            case 'i':
                if ( !strncasecmp("gt:",optarg,3) )
                {
                    if ( args->gt_filter_str ) error("Error: only one -i or -e expression can be given for gt:, and they cannot be combined\n");
                    args->gt_filter_str = optarg;
                    args->gt_filter_logic |= FLT_INCLUDE;
                }
                else if ( !strncasecmp("qry:",optarg,4) )
                {
                    if ( args->qry_filter_str ) error("Error: only one -i or -e expression can be given for qry:, and they cannot be combined\n");
                    args->qry_filter_str = optarg;
                    args->qry_filter_logic |= FLT_INCLUDE;
                }
                else
                {
                    args->gt_filter_str  = optarg;
                    args->qry_filter_str = optarg;
                    args->gt_filter_logic  |= FLT_INCLUDE;
                    args->qry_filter_logic |= FLT_INCLUDE;
                }
                break;
            case 'E':
                args->gt_err = strtol(optarg,&tmp,10);
                if ( !tmp || *tmp ) error("Could not parse: --error-probability %s\n", optarg);
                break;
            case 'u':
                {
                    int i,nlist;
                    char **list = hts_readlist(optarg, 0, &nlist);
                    if ( !list || nlist<=0 || nlist>2 ) error("Failed to parse --use %s\n", optarg);
                    if ( !strcasecmp("GT",list[0]) ) args->qry_use_GT = 1;
                    else if ( !strcasecmp("PL",list[0]) ) args->qry_use_GT = 0;
                    else error("Failed to parse --use %s; only GT and PL are supported\n", optarg);
                    if ( nlist==2 )
                    {
                        if ( !strcasecmp("GT",list[1]) ) args->gt_use_GT = 1;
                        else if ( !strcasecmp("PL",list[1]) ) args->gt_use_GT = 0;
                        else error("Failed to parse --use %s; only GT and PL are supported\n", optarg);
                    }
                    else args->gt_use_GT = args->qry_use_GT;
                    for (i=0; i<nlist; i++) free(list[i]);
                    free(list);
                }
                break;
            case 2 :
                args->ntop = strtol(optarg,&tmp,10);
                if ( !tmp || *tmp ) error("Could not parse: --n-matches %s\n", optarg);
                if ( args->ntop < 0 )
                {
                    args->sort_by_hwe = 1;
                    args->ntop *= -1;
                }
                break;
            case 3 : args->calc_hwe_prob = 0; break;
            case 4 : error("The option -S, --target-sample has been deprecated\n"); break;
            case 5 : args->dry_run = 1; break;
            case 6 :
                args->distinctive_sites = strtod(optarg,&tmp);
                if ( *tmp )
                {
                    if ( *tmp!=',' ) error("Could not parse: --distinctive-sites %s\n", optarg);
                    tmp++;
                    free(args->es_max_mem);
                    args->es_max_mem = strdup(tmp);
                    while ( *tmp && *tmp!=',' ) tmp++;
                    if ( *tmp ) { *tmp = 0; args->es_tmp_prefix = tmp+1; }
                }
                args->gt_err = 0;
                break;
            case 'c':
                error("The -c option is to be implemented, please open an issue on github\n");
                args->min_inter_err = strtod(optarg,&tmp);
                if ( *tmp )
                {
                    if ( *tmp!=',') error("Could not parse: -c %s\n", optarg);
                    args->max_intra_err = strtod(tmp+1,&tmp);
                    if ( *tmp ) error("Could not parse: -c %s\n", optarg);
                }
                break;
            case 'G': error("The option -G, --GTs-only has been deprecated\n"); break;
            case 'a': args->all_sites = 1; error("The -a option is to be implemented, please open an issue on github\n"); break;
            case 'H': args->hom_only = 1; break;
            case 'g': args->gt_fname = optarg; break;
//            case 'p': args->plot = optarg; break;
            case 's':
                if ( !strncasecmp("gt:",optarg,3) ) args->gt_samples = optarg+3;
                else if ( !strncasecmp("qry:",optarg,4) ) args->qry_samples = optarg+4;
                else error("Which one? Query samples (qry:%s) or genotype samples (gt:%s)?\n",optarg,optarg);
                break;
            case 'S':
                if ( !strncasecmp("gt:",optarg,3) ) args->gt_samples = optarg+3, args->gt_samples_is_file = 1;
                else if ( !strncasecmp("qry:",optarg,4) ) args->qry_samples = optarg+4, args->qry_samples_is_file = 1;
                else error("Which one? Query samples (qry:%s) or genotype samples (gt:%s)?\n",optarg,optarg);
                break;
            case 'p': args->pair_samples = optarg; break;
            case 'P': args->pair_samples = optarg; args->pair_samples_is_file = 1; break;
            case 'r': args->regions = optarg; break;
            case 'R': args->regions = optarg; args->regions_is_file = 1; break;
            case 't': args->targets = optarg; break;
            case 'T': args->targets = optarg; args->targets_is_file = 1; break;
            case  7 :
                args->regions_overlap = parse_overlap_option(optarg);
                if ( args->regions_overlap < 0 ) error("Could not parse: --regions-overlap %s\n",optarg);
                break;
            case  8 :
                args->targets_overlap = parse_overlap_option(optarg);
                if ( args->targets_overlap < 0 ) error("Could not parse: --targets-overlap %s\n",optarg);
                break;
            case 'h':
            case '?': usage(); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->qry_fname = "-";  // reading from stdin
        else usage();   // no files given
    }
    else args->qry_fname = argv[optind];
    if ( argc>optind+1 ) error("Error: too many files given, run with -h for help\n");  // too many files given
    if ( args->pair_samples )
    {
        if ( args->gt_samples || args->qry_samples ) error("The -p/-P option cannot be combined with -s/-S\n");
        if ( args->ntop ) error("The --n-matches option cannot be combined with -p/-P\n");
    }
    if ( args->distinctive_sites && !args->pair_samples ) error("The experimental option --distinctive-sites requires -p/-P\n");
    if ( args->hom_only && !args->gt_fname ) error("The option --homs-only requires --genotypes\n");
    if ( args->distinctive_sites && args->gt_err ) error("The option --distinctive-sites cannot be combined with --error-probability\n");

    init_data(args);

    int ret;
    while ( (ret=bcf_sr_next_line(args->files)) )
    {
        if ( !is_input_okay(args,ret) ) continue;

        // time one record to give the user an estimate with very big files
        struct timeval t0, t1;
        if ( !args->ncmp )  gettimeofday(&t0, NULL);

        process_line(args);

        if ( args->ncmp==1 )
        {
            gettimeofday(&t1, NULL);
            double delta = (t1.tv_sec - t0.tv_sec) * 1e6 + (t1.tv_usec - t0.tv_usec);
            fprintf(stderr,"INFO:\tTime required to process one record .. %f seconds\n",delta/1e6);
            args->kstr.l = 0;
            ksprintf(&args->kstr,"INFO\tTime required to process one record .. %f seconds\n",delta/1e6);
            if ( bgzf_write(args->out_fh, args->kstr.s, args->kstr.l)!=args->kstr.l ) error("Failed to write to %s\n", args->output_fname);
            if ( args->dry_run ) break;
        }
    }
    if ( !args->dry_run )
    {
        report(args);
        if ( args->distinctive_sites ) report_distinctive_sites(args);
    }

    destroy_data(args);
    free(args);
    return 0;
}

