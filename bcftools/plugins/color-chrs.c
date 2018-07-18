/* The MIT License

   Copyright (c) 2015 Genome Research Ltd.

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
    Trio haplotypes: mother (A,B), father (C,D), child (E,F)
    Modeling the following states:
        01|23|02 
        01|23|03
        01|23|12
        01|23|13
        01|23|20
        01|23|30
        01|23|21
        01|23|31
    with the likelihoods of two haplotypes A,B segments sharing an allele:
        P(01|A==B)  .. e (P of error)
        P(00|A==B)  .. 1-e
    and
        P(ab,cd,ef|E=A,F=C) = P(ea|E=A)*P(fc|F=C)


    Unrelated samples: (A,B) and (C,D)
    Modeling the states:
        xxxx .. A!=C,A!=D,B!=C,B!=D
        0x0x .. A=C,B!=D
        0xx0 .. A=D,B!=C
        x00x .. B=C,A!=D
        x0x0 .. B=D,A!=C
        0101 .. A=C,B=D
        0110 .. A=D,B=C
    with the likelihoods
        P(01|A!=B)  .. f*(1-f)
        P(00|A!=B)  .. (1-f)*(1-f)
        P(11|A!=B)  .. f*f

    Assuming 2x30 crossovers, P=2e-8.
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <errno.h>
#include "bcftools.h"
#include "HMM.h"

#define C_TRIO 1
#define C_UNRL 2

// states for unrelated samples
#define UNRL_xxxx  0
#define UNRL_0x0x  1
#define UNRL_0xx0  2
#define UNRL_x00x  3
#define UNRL_x0x0  4
#define UNRL_0101  5
#define UNRL_0110  6

// trio states
#define TRIO_AC 0
#define TRIO_AD 1
#define TRIO_BC 2
#define TRIO_BD 3
#define TRIO_CA 4
#define TRIO_DA 5
#define TRIO_CB 6
#define TRIO_DB 7

typedef struct _args_t
{
    bcf_hdr_t *hdr;
    hmm_t *hmm;
    double *eprob, *tprob, pij, pgt_err;
    uint32_t *sites;
    int32_t *gt_arr;
    int nsites, msites, ngt_arr, prev_rid;
    int mode, nstates, nhet_father, nhet_mother;
    int imother,ifather,ichild, isample,jsample;
    void (*set_observed_prob) (bcf1_t *rec);
    char *prefix;
    FILE *fp;
}
args_t;

static args_t args;

#define SW_MOTHER 1
#define SW_FATHER 2
static int hap_switch[8][8];

static void set_observed_prob_trio(bcf1_t *rec);
static void set_observed_prob_unrelated(bcf1_t *rec);
static void init_hmm_trio(args_t *args);
static void init_hmm_unrelated(args_t *args);


const char *about(void)
{
    return "Color shared chromosomal segments, requires phased GTs.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Color shared chromosomal segments, requires phased GTs. The output\n"
        "       can be visualized using the color-chrs.pl script.\n"
        "Usage: bcftools +color-chrs [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -p, --prefix <path>     output files prefix\n"
        "   -t, --trio <m,f,c>      names of mother, father and the child\n"
        "   -u, --unrelated <a,b>   names of two unrelated samples\n"
        "\n"
        "Example:\n"
        "   bcftools +color-chrs in.vcf --\n"
        "\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    char *trio_samples = NULL, *unrelated_samples = NULL;
    memset(&args,0,sizeof(args_t));
    args.prev_rid = -1;
    args.hdr = in;
    args.pij = 2e-8;
    args.pgt_err = 1e-9;

    static struct option loptions[] =
    {
        {"prefix",1,0,'p'},
        {"trio",1,0,'t'},
        {"unrelated",1,0,'u'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "?ht:u:p:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'p': args.prefix = optarg; break;
            case 't': trio_samples = optarg; break;
            case 'u': unrelated_samples = optarg; break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }
    if ( optind != argc ) error("%s",usage());
    if ( trio_samples && unrelated_samples ) error("Expected only one of the -t/-u options\n");
    if ( !trio_samples && !unrelated_samples ) error("Expected one of the -t/-u options\n");
    if ( !args.prefix ) error("Expected the -p option\n");

    int ret = bcf_hdr_set_samples(args.hdr, trio_samples ? trio_samples : unrelated_samples, 0);
    if ( ret<0 ) error("Could not parse samples: %s\n", trio_samples ? trio_samples : unrelated_samples);
    else if ( ret>0 ) error("%d-th sample not found: %s\n", ret,trio_samples ? trio_samples : unrelated_samples);

    if ( trio_samples )
    {
        int i,n = 0;
        char **list = hts_readlist(trio_samples, 0, &n);
        if ( n!=3 ) error("Expected three sample names with -t\n");
        args.imother = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[0]);
        args.ifather = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[1]);
        args.ichild  = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[2]);
        for (i=0; i<n; i++) free(list[i]);
        free(list);
        args.set_observed_prob = set_observed_prob_trio;
        args.mode = C_TRIO;
        init_hmm_trio(&args);
    }
    else
    {
        int i,n = 0;
        char **list = hts_readlist(unrelated_samples, 0, &n);
        if ( n!=2 ) error("Expected two sample names with -u\n");
        args.isample = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[0]);
        args.jsample = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[1]);
        for (i=0; i<n; i++) free(list[i]);
        free(list);
        args.set_observed_prob = set_observed_prob_unrelated;
        args.mode = C_UNRL;
        init_hmm_unrelated(&args);
    }
    return 1;
}

static void init_hmm_trio(args_t *args)
{
    int i,j;
    args->nstates = 8;
    args->tprob   = (double*) malloc(sizeof(double)*args->nstates*args->nstates);

    for (i=0; i<args->nstates; i++)
        for (j=0; j<args->nstates; j++) hap_switch[i][j] = 0;

    hap_switch[TRIO_AD][TRIO_AC] = SW_FATHER;
    hap_switch[TRIO_BC][TRIO_AC] = SW_MOTHER;
    hap_switch[TRIO_BD][TRIO_AC] = SW_MOTHER | SW_FATHER;
    hap_switch[TRIO_AC][TRIO_AD] = SW_FATHER;
    hap_switch[TRIO_BC][TRIO_AD] = SW_MOTHER | SW_FATHER;
    hap_switch[TRIO_BD][TRIO_AD] = SW_MOTHER;
    hap_switch[TRIO_AC][TRIO_BC] = SW_MOTHER;
    hap_switch[TRIO_AD][TRIO_BC] = SW_MOTHER | SW_FATHER;
    hap_switch[TRIO_BD][TRIO_BC] = SW_FATHER;
    hap_switch[TRIO_AC][TRIO_BD] = SW_MOTHER | SW_FATHER;
    hap_switch[TRIO_AD][TRIO_BD] = SW_MOTHER;
    hap_switch[TRIO_BC][TRIO_BD] = SW_FATHER;

    hap_switch[TRIO_DA][TRIO_CA] = SW_FATHER;
    hap_switch[TRIO_CB][TRIO_CA] = SW_MOTHER;
    hap_switch[TRIO_DB][TRIO_CA] = SW_MOTHER | SW_FATHER;
    hap_switch[TRIO_CA][TRIO_DA] = SW_FATHER;
    hap_switch[TRIO_CB][TRIO_DA] = SW_MOTHER | SW_FATHER;
    hap_switch[TRIO_DB][TRIO_DA] = SW_MOTHER;
    hap_switch[TRIO_CA][TRIO_CB] = SW_MOTHER;
    hap_switch[TRIO_DA][TRIO_CB] = SW_MOTHER | SW_FATHER;
    hap_switch[TRIO_DB][TRIO_CB] = SW_FATHER;
    hap_switch[TRIO_CA][TRIO_DB] = SW_MOTHER | SW_FATHER;
    hap_switch[TRIO_DA][TRIO_DB] = SW_MOTHER;
    hap_switch[TRIO_CB][TRIO_DB] = SW_FATHER;

    for (i=0; i<args->nstates; i++)
    {
        for (j=0; j<args->nstates; j++)
        {
            if ( !hap_switch[i][j] ) MAT(args->tprob,args->nstates,i,j) = 0;
            else
            {
                MAT(args->tprob,args->nstates,i,j) = 1;
                if ( hap_switch[i][j] & SW_MOTHER ) MAT(args->tprob,args->nstates,i,j) *= args->pij;
                if ( hap_switch[i][j] & SW_FATHER ) MAT(args->tprob,args->nstates,i,j) *= args->pij;
            }
        }
    }
    for (i=0; i<args->nstates; i++)
    {
        double sum = 0;
        for (j=0; j<args->nstates; j++)
        {
            if ( i!=j ) sum += MAT(args->tprob,args->nstates,i,j);
        }
        MAT(args->tprob,args->nstates,i,i) = 1 - sum;
    }

    #if 0
    for (i=0; i<args->nstates; i++)
    {
        for (j=0; j<args->nstates; j++)
            fprintf(stderr,"\t%d",hap_switch[j][i]);
        fprintf(stderr,"\n");
    }
    for (i=0; i<args->nstates; i++)
    {
        for (j=0; j<args->nstates; j++)
            fprintf(stderr,"\t%e",MAT(args->tprob,args->nstates,j,i));
        fprintf(stderr,"\n");
    }
    #endif

    args->hmm = hmm_init(args->nstates, args->tprob, 10000);
}
static void init_hmm_unrelated(args_t *args)
{
    int i,j;
    args->nstates = 7;
    args->tprob   = (double*) malloc(sizeof(double)*args->nstates*args->nstates);

    for (i=0; i<args->nstates; i++)
    {
        for (j=0; j<args->nstates; j++)
            MAT(args->tprob,args->nstates,i,j) = args->pij;
    }
    MAT(args->tprob,args->nstates,UNRL_0101,UNRL_xxxx) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,UNRL_0110,UNRL_xxxx) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,UNRL_x0x0,UNRL_0x0x) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,UNRL_0110,UNRL_0x0x) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,UNRL_x00x,UNRL_0xx0) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,UNRL_0101,UNRL_0xx0) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,UNRL_0101,UNRL_x00x) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,UNRL_0110,UNRL_x0x0) = args->pij*args->pij;
    MAT(args->tprob,args->nstates,UNRL_0110,UNRL_0101) = args->pij*args->pij;

    for (i=0; i<args->nstates; i++)
    {
        for (j=i+1; j<args->nstates; j++)
            MAT(args->tprob,args->nstates,i,j) = MAT(args->tprob,args->nstates,j,i);
    }
    for (i=0; i<args->nstates; i++)
    {
        double sum = 0;
        for (j=0; j<args->nstates; j++)
            if ( i!=j ) sum += MAT(args->tprob,args->nstates,i,j);
        MAT(args->tprob,args->nstates,i,i) = 1 - sum;
    }

    #if 0
    for (i=0; i<args->nstates; i++)
    {
        for (j=0; j<args->nstates; j++)
            fprintf(stderr,"\t%e",MAT(args->tprob,args->nstates,j,i));
        fprintf(stderr,"\n");
    }
    #endif

    args->hmm = hmm_init(args->nstates, args->tprob, 10000);
}
static inline double prob_shared(float af, int a, int b)
{
    return a==b ? 1-args.pgt_err : args.pgt_err;
}
static inline double prob_not_shared(float af, int a, int b)
{
    if ( a!=b ) return af*(1-af);
    else if ( a==0 ) return (1-af)*(1-af);
    else return af*af;
}
static void set_observed_prob_unrelated(bcf1_t *rec)
{
    float af = 0.5;  // alternate allele frequency

    int ngt = bcf_get_genotypes(args.hdr, rec, &args.gt_arr, &args.ngt_arr);
    if ( ngt<0 ) return;
    if ( ngt!=4 ) return;   // chrX

    int32_t a,b,c,d;
    a = args.gt_arr[2*args.isample];
    b = args.gt_arr[2*args.isample+1];
    c = args.gt_arr[2*args.jsample];
    d = args.gt_arr[2*args.jsample+1];
    if ( bcf_gt_is_missing(a) || bcf_gt_is_missing(b) ) return;
    if ( bcf_gt_is_missing(c) || bcf_gt_is_missing(d) ) return;
    if ( !bcf_gt_is_phased(a) && !bcf_gt_is_phased(b) ) return; // only the second allele should be set when phased
    if ( !bcf_gt_is_phased(c) && !bcf_gt_is_phased(d) ) return;
    a = bcf_gt_allele(a);
    b = bcf_gt_allele(b);
    c = bcf_gt_allele(c);
    d = bcf_gt_allele(d);

    int m = args.msites;
    args.nsites++;
    hts_expand(uint32_t,args.nsites,args.msites,args.sites);
    if ( m!=args.msites ) args.eprob = (double*) realloc(args.eprob, sizeof(double)*args.msites*args.nstates);

    args.sites[args.nsites-1] = rec->pos;
    double *prob = args.eprob + args.nstates*(args.nsites-1);
    prob[UNRL_xxxx] = prob_not_shared(af,a,c) * prob_not_shared(af,a,d) * prob_not_shared(af,b,c) * prob_not_shared(af,b,d);
    prob[UNRL_0x0x] = prob_shared(af,a,c) * prob_not_shared(af,b,d);
    prob[UNRL_0xx0] = prob_shared(af,a,d) * prob_not_shared(af,b,c);
    prob[UNRL_x00x] = prob_shared(af,b,c) * prob_not_shared(af,a,d);
    prob[UNRL_x0x0] = prob_shared(af,b,d) * prob_not_shared(af,a,c);
    prob[UNRL_0101] = prob_shared(af,a,c) * prob_shared(af,b,d);
    prob[UNRL_0110] = prob_shared(af,a,d) * prob_shared(af,b,c);

#if 0
    static int x = 0;
    if ( !x++)
    {
        printf("p(0==0) .. %f\n", prob_shared(af,0,0));
        printf("p(0!=0) .. %f\n", prob_not_shared(af,0,0));
        printf("p(0==1) .. %f\n", prob_shared(af,0,1));
        printf("p(0!=1) .. %f\n", prob_not_shared(af,0,1));
    }
    printf("%d|%d %d|%d  x:%f 11:%f 12:%f 21:%f 22:%f 11,22:%f 12,21:%f  %d\n", a,b,c,d,
            prob[UNRL_xxxx], prob[UNRL_0x0x], prob[UNRL_0xx0], prob[UNRL_x00x], prob[UNRL_x0x0], prob[UNRL_0101], prob[UNRL_0110], rec->pos+1);
#endif
}
static void set_observed_prob_trio(bcf1_t *rec)
{
    int ngt = bcf_get_genotypes(args.hdr, rec, &args.gt_arr, &args.ngt_arr);
    if ( ngt<0 ) return;
    if ( ngt!=6 ) return;   // chrX

    int32_t a,b,c,d,e,f;
    a = args.gt_arr[2*args.imother];
    b = args.gt_arr[2*args.imother+1];
    c = args.gt_arr[2*args.ifather];
    d = args.gt_arr[2*args.ifather+1];
    e = args.gt_arr[2*args.ichild];
    f = args.gt_arr[2*args.ichild+1];
    if ( bcf_gt_is_missing(a) || bcf_gt_is_missing(b) ) return;
    if ( bcf_gt_is_missing(c) || bcf_gt_is_missing(d) ) return;
    if ( bcf_gt_is_missing(e) || bcf_gt_is_missing(f) ) return;
    if ( !bcf_gt_is_phased(a) && !bcf_gt_is_phased(b) ) return; // only the second allele should be set when phased
    if ( !bcf_gt_is_phased(c) && !bcf_gt_is_phased(d) ) return;
    if ( !bcf_gt_is_phased(e) && !bcf_gt_is_phased(f) ) return;
    a = bcf_gt_allele(a);
    b = bcf_gt_allele(b);
    c = bcf_gt_allele(c);
    d = bcf_gt_allele(d);
    e = bcf_gt_allele(e);
    f = bcf_gt_allele(f);

    int mother = (1<<a) | (1<<b);
    int father = (1<<c) | (1<<d);
    int child  = (1<<e) | (1<<f);
    if ( !(mother&child) || !(father&child) )  return;      // Mendelian-inconsistent site, skip

    if ( a!=b ) args.nhet_mother++;
    if ( c!=d ) args.nhet_father++;

    int m = args.msites;
    args.nsites++;
    hts_expand(uint32_t,args.nsites,args.msites,args.sites);
    if ( m!=args.msites ) args.eprob = (double*) realloc(args.eprob, sizeof(double)*args.msites*args.nstates);

    args.sites[args.nsites-1] = rec->pos;
    double *prob = args.eprob + args.nstates*(args.nsites-1);
    prob[TRIO_AC] = prob_shared(0,e,a) * prob_shared(0,f,c);
    prob[TRIO_AD] = prob_shared(0,e,a) * prob_shared(0,f,d);
    prob[TRIO_BC] = prob_shared(0,e,b) * prob_shared(0,f,c);
    prob[TRIO_BD] = prob_shared(0,e,b) * prob_shared(0,f,d);
    prob[TRIO_CA] = prob_shared(0,e,c) * prob_shared(0,f,a);
    prob[TRIO_DA] = prob_shared(0,e,d) * prob_shared(0,f,a);
    prob[TRIO_CB] = prob_shared(0,e,c) * prob_shared(0,f,b);
    prob[TRIO_DB] = prob_shared(0,e,d) * prob_shared(0,f,b);
}

void flush_viterbi(args_t *args)
{
    const char *s1, *s2, *s3 = NULL;
    if ( args->mode==C_UNRL )
    {
        s1 = bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,args->isample);
        s2 = bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,args->jsample);
    }
    else if ( args->mode==C_TRIO )
    {
        s1 = bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,args->imother);
        s3 = bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,args->ifather);
        s2 = bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,args->ichild);
    }
    else abort();

    if ( !args->fp )
    {
        kstring_t str = {0,0,0};
        kputs(args->prefix, &str);
        kputs(".dat", &str);
        args->fp = fopen(str.s,"w");
        if ( !args->fp ) error("%s: %s\n", str.s,strerror(errno));
        free(str.s);
        fprintf(args->fp,"# SG, shared segment\t[2]Chromosome\t[3]Start\t[4]End\t[5]%s:1\t[6]%s:2\n",s2,s2);
        fprintf(args->fp,"# SW, number of switches\t[3]Sample\t[4]Chromosome\t[5]nHets\t[5]nSwitches\t[6]switch rate\n");
    }

    hmm_run_viterbi(args->hmm,args->nsites,args->eprob,args->sites);
    uint8_t *vpath = hmm_get_viterbi_path(args->hmm);
    int i, iprev = -1, prev_state = -1, nstates = hmm_get_nstates(args->hmm);
    int nswitch_mother = 0, nswitch_father = 0;
    for (i=0; i<args->nsites; i++)
    {
        int state = vpath[i*nstates];
        if ( state!=prev_state || i+1==args->nsites )
        {
            uint32_t start = iprev>=0 ? args->sites[iprev]+1 : 1, end = i>0 ? args->sites[i-1] : 1;
            const char *chr = bcf_hdr_id2name(args->hdr,args->prev_rid);
            if ( args->mode==C_UNRL )
            {
                switch (prev_state)
                {
                    case UNRL_0x0x:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:1\t-\n", chr,start,end,s1); break;
                    case UNRL_0xx0:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t-\t%s:1\n", chr,start,end,s1); break;
                    case UNRL_x00x:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:2\t-\n", chr,start,end,s1); break;
                    case UNRL_x0x0:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t-\t%s:2\n", chr,start,end,s1); break;
                    case UNRL_0101:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:1\t%s:2\n", chr,start,end,s1,s1); break;
                    case UNRL_0110:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:2\t%s:1\n", chr,start,end,s1,s1); break;
                }
            }
            else if ( args->mode==C_TRIO )
            {
                switch (prev_state)
                {
                    case TRIO_AC:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:1\t%s:1\n", chr,start,end,s1,s3); break;
                    case TRIO_AD:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:1\t%s:2\n", chr,start,end,s1,s3); break;
                    case TRIO_BC:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:2\t%s:1\n", chr,start,end,s1,s3); break;
                    case TRIO_BD:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:2\t%s:2\n", chr,start,end,s1,s3); break;
                    case TRIO_CA:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:1\t%s:1\n", chr,start,end,s3,s1); break;
                    case TRIO_DA:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:2\t%s:1\n", chr,start,end,s3,s1); break;
                    case TRIO_CB:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:1\t%s:2\n", chr,start,end,s3,s1); break;
                    case TRIO_DB:
                        fprintf(args->fp,"SG\t%s\t%d\t%d\t%s:2\t%s:2\n", chr,start,end,s3,s1); break;
                }
                if ( hap_switch[state][prev_state] & SW_MOTHER ) nswitch_mother++;
                if ( hap_switch[state][prev_state] & SW_FATHER ) nswitch_father++;
            }
            iprev = i-1;
        }
        prev_state = state;
    }
    float mrate = args->nhet_mother>1 ? (float)nswitch_mother/(args->nhet_mother-1) : 0;
    float frate = args->nhet_father>1 ? (float)nswitch_father/(args->nhet_father-1) : 0;
    fprintf(args->fp,"SW\t%s\t%s\t%d\t%d\t%f\n", s1,bcf_hdr_id2name(args->hdr,args->prev_rid),args->nhet_mother,nswitch_mother,mrate);
    fprintf(args->fp,"SW\t%s\t%s\t%d\t%d\t%f\n", s3,bcf_hdr_id2name(args->hdr,args->prev_rid),args->nhet_father,nswitch_father,frate);
    args->nsites = 0;
    args->nhet_father = args->nhet_mother = 0;
}
    
bcf1_t *process(bcf1_t *rec)
{
    if ( args.prev_rid==-1 ) args.prev_rid = rec->rid;
    if ( args.prev_rid!=rec->rid ) flush_viterbi(&args);
    args.prev_rid = rec->rid;
    args.set_observed_prob(rec);
    return NULL;
}

void destroy(void)
{
    flush_viterbi(&args);
    fclose(args.fp);

    free(args.gt_arr);
    free(args.tprob);
    free(args.sites);
    free(args.eprob);
    hmm_destroy(args.hmm);
}



