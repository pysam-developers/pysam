/*  plugins/setGT.c -- set gentoypes to given values

    Copyright (C) 2015-2017 Genome Research Ltd.

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/kfunc.h>
#include <inttypes.h>
#include <getopt.h>
#include <ctype.h>
#include "bcftools.h"
#include "filter.h"

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef int (*cmp_f)(double a, double b);

static int cmp_eq(double a, double b) { return a==b ? 1 : 0; }
static int cmp_le(double a, double b) { return a<=b ? 1 : 0; }
static int cmp_ge(double a, double b) { return a>=b ? 1 : 0; }
static int cmp_lt(double a, double b) { return a<b ? 1 : 0; }
static int cmp_gt(double a, double b) { return a>b ? 1 : 0; }

typedef struct
{
    bcf_hdr_t *in_hdr, *out_hdr;
    int32_t *gts, mgts, *iarr, miarr;
    int *arr, marr;
    uint64_t nchanged;
    int tgt_mask, new_mask, new_gt;
    filter_t *filter;
    char *filter_str;
    int filter_logic;
    const uint8_t *smpl_pass;
    double binom_val;
    char *binom_tag;
    cmp_f binom_cmp;
}
args_t;

args_t *args = NULL;

#define GT_MISSING   1
#define GT_PARTIAL  (1<<1)
#define GT_REF      (1<<2)
#define GT_MAJOR    (1<<3)
#define GT_PHASED   (1<<4)
#define GT_UNPHASED (1<<5)
#define GT_ALL      (1<<6)
#define GT_QUERY    (1<<7)
#define GT_BINOM    (1<<8)

const char *about(void)
{
    return "Set genotypes: partially missing to missing, missing to ref/major allele, etc.\n";
}

const char *usage(void)
{
    return 
        "About: Sets genotypes. The target genotypes can be specified as:\n"
        "           ./.  .. completely missing (\".\" or \"./.\", depending on ploidy)\n"
        "           ./x  .. partially missing (e.g., \"./0\" or \".|1\" but not \"./.\")\n"
        "           .    .. partially or completely missing\n"
        "           a    .. all genotypes\n"
        "           b    .. heterozygous genotypes failing two-tailed binomial test (example below)\n"
        "           q    .. select genotypes using -i/-e options\n"
        "       and the new genotype can be one of:\n"
        "           .    .. missing (\".\" or \"./.\", keeps ploidy)\n"
        "           0    .. reference allele\n"
        "           M    .. major allele\n"
        "           p    .. phased genotype\n"
        "           u    .. unphase genotype and sort by allele (1|0 becomes 0/1)\n"
        "Usage: bcftools +setGT [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -e, --exclude <expr>        Exclude a genotype if true (requires -t q)\n"
        "   -i, --include <expr>        include a genotype if true (requires -t q)\n"
        "   -n, --new-gt <type>         Genotypes to set, see above\n"
        "   -t, --target-gt <type>      Genotypes to change, see above\n"
        "\n"
        "Example:\n"
        "   # set missing genotypes (\"./.\") to phased ref genotypes (\"0|0\")\n"
        "   bcftools +setGT in.vcf -- -t . -n 0p\n"
        "\n"
        "   # set missing genotypes with DP>0 and GQ>20 to ref genotypes (\"0/0\")\n"
        "   bcftools +setGT in.vcf -- -t q -n 0 -i 'GT=\".\" && FMT/DP>0 && GQ>20'\n"
        "\n"
        "   # set partially missing genotypes to completely missing\n"
        "   bcftools +setGT in.vcf -- -t ./x -n .\n"
        "\n"
        "   # set heterozygous genotypes to 0/0 if binom.test(nAlt,nRef+nAlt,0.5)<1e-3\n"
        "   bcftools +setGT in.vcf -- -t \"b:AD<1e-3\" -n 0\n"  // todo: make -i/-e recognise something like is_het or gt="het" so that this can be generalized?
        "\n";
}

void parse_binom_expr(args_t *args, char *str)
{
    if ( str[1]!=':' ) goto err;

    char *beg = str+2;
    while ( *beg && isspace(*beg) ) beg++;
    if ( !*beg ) goto err;
    char *end = beg;
    while ( *end )
    {
        if ( isspace(*end) || *end=='<' || *end=='=' || *end=='>' ) break;
        end++;
    }
    if ( !*end ) goto err;
    args->binom_tag = (char*) calloc(1,end-beg+1);
    memcpy(args->binom_tag,beg,end-beg);
    int tag_id = bcf_hdr_id2int(args->in_hdr,BCF_DT_ID,args->binom_tag);
    if ( !bcf_hdr_idinfo_exists(args->in_hdr,BCF_HL_FMT,tag_id) ) error("The FORMAT tag \"%s\" is not present in the VCF\n", args->binom_tag);
    
    while ( *end && isspace(*end) ) end++;
    if ( !*end ) goto err;

    if ( !strncmp(end,"<=",2) ) { args->binom_cmp = cmp_le; beg = end+2; }
    else if ( !strncmp(end,">=",2) ) { args->binom_cmp = cmp_ge; beg = end+2; }
    else if ( !strncmp(end,"==",2) ) { args->binom_cmp = cmp_eq; beg = end+2; }
    else if ( !strncmp(end,"<",1) ) { args->binom_cmp = cmp_lt; beg = end+1; }
    else if ( !strncmp(end,">",1) ) { args->binom_cmp = cmp_gt; beg = end+1; }
    else if ( !strncmp(end,"=",1) ) { args->binom_cmp = cmp_eq; beg = end+1; }
    else goto err;

    while ( *beg && isspace(*beg) ) beg++;
    if ( !*beg ) goto err;

    args->binom_val = strtod(beg, &end);
    while ( *end && isspace(*end) ) end++;
    if ( *end ) goto err;

    args->tgt_mask |= GT_BINOM;
    return;

err:
    error(
        "Error parsing the expression: %s\n"
        "Expected TAG CMP VAL, where\n"
        "   TAG .. one of the format tags\n"
        "   CMP .. operator, one of <, <=, >, >=\n"
        "   VAL .. value\n"
        "For example:\n"
        "   bcftools +setGT in.vcf -- -t \"b:AD>1e-3\" -n 0\n"
        "\n", str
        );
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    args = (args_t*) calloc(1,sizeof(args_t));
    args->in_hdr  = in;
    args->out_hdr = out;

    int c;
    static struct option loptions[] =
    {
        {"include",required_argument,NULL,'i'},
        {"exclude",required_argument,NULL,'e'},
        {"new-gt",required_argument,NULL,'n'},
        {"target-gt",required_argument,NULL,'t'},
        {NULL,0,NULL,0}
    };
    while ((c = getopt_long(argc, argv, "?hn:t:i:e:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'i': args->filter_str = optarg; args->filter_logic = FLT_INCLUDE; break;
            case 'e': args->filter_str = optarg; args->filter_logic = FLT_EXCLUDE; break;
            case 'n': args->new_mask = bcf_gt_phased(0); 
                if ( strchr(optarg,'.') ) args->new_mask |= GT_MISSING;
                if ( strchr(optarg,'0') ) args->new_mask |= GT_REF;
                if ( strchr(optarg,'M') ) args->new_mask |= GT_MAJOR;
                if ( strchr(optarg,'p') ) args->new_mask |= GT_PHASED;
                if ( strchr(optarg,'u') ) args->new_mask |= GT_UNPHASED;
                if ( args->new_mask==0 ) error("Unknown parameter to --new-gt: %s\n", optarg);
                break;
            case 't':
                if ( !strcmp(optarg,".") ) args->tgt_mask |= GT_MISSING|GT_PARTIAL;
                if ( !strcmp(optarg,"./x") ) args->tgt_mask |= GT_PARTIAL;
                if ( !strcmp(optarg,"./.") ) args->tgt_mask |= GT_MISSING;
                if ( !strcmp(optarg,"a") ) args->tgt_mask |= GT_ALL;
                if ( !strcmp(optarg,"q") ) args->tgt_mask |= GT_QUERY;
                if ( !strcmp(optarg,"?") ) args->tgt_mask |= GT_QUERY;        // for backward compatibility
                if ( strchr(optarg,'b') ) parse_binom_expr(args, strchr(optarg,'b'));
                if ( args->tgt_mask==0 ) error("Unknown parameter to --target-gt: %s\n", optarg);
                break;
            case 'h':
            case '?':
            default: fprintf(stderr,"%s", usage()); exit(1); break;
        }
    }

    if ( !args->new_mask ) error("Expected -n option\n");
    if ( !args->tgt_mask ) error("Expected -t option\n");

    if ( args->new_mask & GT_MISSING ) args->new_gt = bcf_gt_missing;
    if ( args->new_mask & GT_REF ) args->new_gt = args->new_mask&GT_PHASED ? bcf_gt_phased(0) : bcf_gt_unphased(0);

    if ( args->filter_str  && !(args->tgt_mask&GT_QUERY) ) error("Expected -tq with -i/-e\n");
    if ( !args->filter_str && args->tgt_mask&GT_QUERY ) error("Expected -i/-e with -tq\n");
    if ( args->filter_str ) args->filter = filter_init(in,args->filter_str);

    return 0;
}

static inline int unphase_gt(int32_t *ptr, int ngts)
{
    int j, changed = 0;
    for (j=0; j<ngts; j++)
    {
        if ( ptr[j]==bcf_int32_vector_end ) break;
        if ( !bcf_gt_is_phased(ptr[j]) ) continue;
        ptr[j] = bcf_gt_unphased(bcf_gt_allele(ptr[j]));    // remove phasing
        changed++;
    }

    // insertion sort
    int k, l;
    for (k=1; k<j; k++)
    {
        int32_t x = ptr[k];
        l = k;
        while ( l>0 && ptr[l-1]>x )
        {
            ptr[l] = ptr[l-1];
            l--;
        }
        ptr[l] = x;
    }
    return changed;
}
static inline int set_gt(int32_t *ptr, int ngts, int gt)
{
    int j, changed = 0;
    for (j=0; j<ngts; j++)
    {
        if ( ptr[j]==bcf_int32_vector_end ) break;
        if ( ptr[j] != gt ) changed++;
        ptr[j] = gt;
    }
    return changed;
}

static inline double calc_binom(int na, int nb)
{
    int N = na + nb;
    if ( !N ) return 1;

    /*
        kfunc.h implements kf_betai, which is the regularized beta function I_x(a,b) = P(X<=a/(a+b))
    */
    double prob = 2 * kf_betai(na, nb+1, 0.5);
    if ( prob > 1 ) prob = 1;

    return prob;
}

bcf1_t *process(bcf1_t *rec)
{
    if ( !rec->n_sample ) return rec;

    int ngts = bcf_get_genotypes(args->in_hdr, rec, &args->gts, &args->mgts);
    ngts /= rec->n_sample;
    int i, j, changed = 0;

    int nbinom = 0;
    if ( args->tgt_mask & GT_BINOM )
    {
        nbinom = bcf_get_format_int32(args->in_hdr, rec, args->binom_tag, &args->iarr, &args->miarr);
        if ( nbinom<0 ) nbinom = 0;
        nbinom /= rec->n_sample;
    }
    
    // Calculating allele frequency for each allele and determining major allele
    // only do this if use_major is true
    int an = 0, maxAC = -1, majorAllele = -1;
    if ( args->new_mask & GT_MAJOR )
    {
        hts_expand(int,rec->n_allele,args->marr,args->arr);
        int ret = bcf_calc_ac(args->in_hdr,rec,args->arr,BCF_UN_FMT);
        if ( ret<= 0 )
            error("Could not calculate allele count at %s:%d\n", bcf_seqname(args->in_hdr,rec),rec->pos+1);

        for(i=0; i < rec->n_allele; ++i)
        {
            an += args->arr[i];
            if (args->arr[i] > maxAC)
            {
                maxAC = args->arr[i];
                majorAllele = i;
            }
        }

        // replacing new_gt by major allele
        args->new_gt = args->new_mask & GT_PHASED ?  bcf_gt_phased(majorAllele) : bcf_gt_unphased(majorAllele);
    }

    // replace gts
    if ( nbinom && ngts>=2 )    // only diploid genotypes are considered: higher ploidy ignored further, haploid here
    {
        if ( args->filter ) filter_test(args->filter,rec,&args->smpl_pass);
        for (i=0; i<rec->n_sample; i++)
        {
            if ( args->smpl_pass )
            {
                if ( !args->smpl_pass[i] && args->filter_logic==FLT_INCLUDE ) continue;
                if (  args->smpl_pass[i] && args->filter_logic==FLT_EXCLUDE ) continue;
            }
            int32_t *ptr = args->gts + i*ngts;
            if ( bcf_gt_is_missing(ptr[0]) || bcf_gt_is_missing(ptr[1]) || ptr[1]==bcf_int32_vector_end ) continue;
            if ( ptr[0]==ptr[1] ) continue; // a hom
            int ia = bcf_gt_allele(ptr[0]); 
            int ib = bcf_gt_allele(ptr[1]); 
            if ( ia>=nbinom || ib>=nbinom ) 
                error("The sample %s has incorrect number of %s fields at %s:%d\n",
                        args->in_hdr->samples[i],args->binom_tag,bcf_seqname(args->in_hdr,rec),rec->pos+1);

            double prob = calc_binom(args->iarr[i*nbinom+ia],args->iarr[i*nbinom+ib]);
            if ( !args->binom_cmp(prob,args->binom_val) ) continue;

            if ( args->new_mask&GT_UNPHASED )
                changed += unphase_gt(ptr, ngts);
            else
                changed += set_gt(ptr, ngts, args->new_gt);
        }
    }
    else if ( args->tgt_mask&GT_QUERY )
    {
        int pass_site = filter_test(args->filter,rec,&args->smpl_pass);
        if ( (pass_site && args->filter_logic==FLT_EXCLUDE) || (!pass_site && args->filter_logic==FLT_INCLUDE) ) return rec;
        for (i=0; i<rec->n_sample; i++)
        {
            if ( args->smpl_pass )
            {
                if ( !args->smpl_pass[i] && args->filter_logic==FLT_INCLUDE ) continue;
                if (  args->smpl_pass[i] && args->filter_logic==FLT_EXCLUDE ) continue;
            }

            if ( args->new_mask&GT_UNPHASED )
                changed += unphase_gt(args->gts + i*ngts, ngts);
            else
                changed += set_gt(args->gts + i*ngts, ngts, args->new_gt);
        }
    }
    else
    {
        for (i=0; i<rec->n_sample; i++)
        {
            int ploidy = 0, nmiss = 0;
            int32_t *ptr = args->gts + i*ngts;
            for (j=0; j<ngts; j++)
            {
                if ( ptr[j]==bcf_int32_vector_end ) break;
                ploidy++;
                if ( ptr[j]==bcf_gt_missing ) nmiss++;
            }

            int do_set = 0;
            if ( args->tgt_mask&GT_ALL ) do_set = 1;
            else if ( args->tgt_mask&GT_PARTIAL && nmiss ) do_set = 1;
            else if ( args->tgt_mask&GT_MISSING && ploidy==nmiss ) do_set = 1;

            if ( !do_set ) continue;

            if ( args->new_mask&GT_UNPHASED )
                changed += unphase_gt(ptr, ngts);
            else
                changed += set_gt(ptr, ngts, args->new_gt);
        }
    }
    args->nchanged += changed;
    if ( changed ) bcf_update_genotypes(args->out_hdr, rec, args->gts, ngts*rec->n_sample);
    return rec;
}

void destroy(void)
{
    fprintf(stderr,"Filled %"PRId64" alleles\n", args->nchanged);
    free(args->binom_tag);
    if ( args->filter ) filter_destroy(args->filter);
    free(args->arr);
    free(args->iarr);
    free(args->gts);
    free(args);
}


