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

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/khash_str2int.h>
#include "bcftools.h"

typedef struct
{
    int father, mother, child;      // VCF sample index
    int prev, ipop;
    uint32_t err, nswitch, ntest;
}
trio_t;

typedef struct
{
    char *name;
    uint32_t err, nswitch, ntest, ntrio;
    float pswitch;
}
pop_t;

typedef struct
{
    int argc;
    char **argv;
    bcf_hdr_t *hdr;
    trio_t *trio;
    int ntrio, mtrio;
    int32_t *gt_arr;
    int npop;
    pop_t *pop;
    int mgt_arr, prev_rid;
}
args_t;

args_t args;

const char *about(void)
{
    return "Calculate phase switch rate in trio samples, children samples must have phased GTs.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Calculate phase switch rate in trio children.\n"
        "Usage: bcftools +trio-swich-rate [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -p, --ped <file>        PED file with optional 7th column to group\n"
        "                           results by population\n"
        "\n"
        "Example:\n"
        "   bcftools +trio-switch-rate file.bcf -- -p file.ped\n"
        "\n";
}

void parse_ped(args_t *args, char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Could not read: %s\n", fname);

    kstring_t str = {0,0,0};
    if ( hts_getline(fp, KS_SEP_LINE, &str) <= 0 ) error("Empty file: %s\n", fname);

    void *pop2i = khash_str2int_init();

    int moff = 0, *off = NULL;
    do
    {
        // familyID    sampleID paternalID maternalID sex   phenotype   population relationship   siblings   secondOrder   thirdOrder   children    comment
        // BB03    HG01884 HG01885 HG01956 2   0   ACB child   0   0   0   0
        int ncols = ksplit_core(str.s,0,&moff,&off);
        if ( ncols<4 ) error("Could not parse the ped file: %s\n", str.s);

        int father = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[2]]);
        if ( father<0 ) continue;
        int mother = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[3]]);
        if ( mother<0 ) continue;
        int child = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[1]]);
        if ( child<0 ) continue;

        args->ntrio++;
        hts_expand0(trio_t,args->ntrio,args->mtrio,args->trio);
        trio_t *trio = &args->trio[args->ntrio-1];
        trio->father = father;
        trio->mother = mother;
        trio->child  = child;

        if (ncols>6) {
            char *pop_name = &str.s[off[6]];
            if ( !khash_str2int_has_key(pop2i,pop_name) )
            {
                pop_name = strdup(&str.s[off[6]]);
                khash_str2int_set(pop2i,pop_name,args->npop);
                args->npop++;
                args->pop = (pop_t*) realloc(args->pop,args->npop*sizeof(*args->pop));
                memset(args->pop+args->npop-1,0,sizeof(*args->pop));
                args->pop[args->npop-1].name = pop_name;
            }
            khash_str2int_get(pop2i,pop_name,&trio->ipop);
            args->pop[trio->ipop].ntrio++;
        }
    } while ( hts_getline(fp, KS_SEP_LINE, &str)>=0 );

    khash_str2int_destroy(pop2i);
    free(str.s);
    free(off);
    hts_close(fp);
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    memset(&args,0,sizeof(args_t));
    args.argc   = argc; args.argv = argv;
    args.prev_rid = -1;
    args.hdr = in;
    char *ped_fname = NULL;
    static struct option loptions[] =
    {
        {"ped",required_argument,NULL,'p'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "?hp:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'p': ped_fname = optarg; break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }
    if ( !ped_fname ) error("Expected the -p option\n");
    parse_ped(&args, ped_fname);
    return 1;
}

typedef struct
{
    int a, b, phased;
}
gt_t;

int parse_genotype(gt_t *gt, int32_t *ptr);

inline int parse_genotype(gt_t *gt, int32_t *ptr)
{
    if ( ptr[0]==bcf_gt_missing ) return 0;
    if ( ptr[1]==bcf_gt_missing ) return 0;
    if ( ptr[1]==bcf_int32_vector_end ) return 0;
    gt->phased = bcf_gt_is_phased(ptr[1]) ? 1 : 0;
    gt->a = bcf_gt_allele(ptr[0]); if ( gt->a > 1 ) return 0;  // consider only the first two alleles at biallelic sites
    gt->b = bcf_gt_allele(ptr[1]); if ( gt->b > 1 ) return 0;
    return 1;
}

bcf1_t *process(bcf1_t *rec)
{
    int ngt = bcf_get_genotypes(args.hdr, rec, &args.gt_arr, &args.mgt_arr);
    if ( ngt<0 ) return NULL;
    ngt /= bcf_hdr_nsamples(args.hdr);
    if ( ngt!=2 ) return NULL;

    int i;
    if ( rec->rid!=args.prev_rid )
    {
        args.prev_rid = rec->rid;
        for (i=0; i<args.ntrio; i++) args.trio[i].prev = 0;
    }

    gt_t child, father, mother;
    for (i=0; i<args.ntrio; i++)
    {
        trio_t *trio = &args.trio[i];

        if ( !parse_genotype(&child, args.gt_arr + ngt*trio->child) ) continue;
        if ( !child.phased ) continue;
        if ( child.a+child.b != 1 ) continue;       // child is not a het

        if ( !parse_genotype(&father, args.gt_arr + ngt*trio->father) ) continue;
        if ( !parse_genotype(&mother, args.gt_arr + ngt*trio->mother) ) continue;
        if ( father.a+father.b == 1 && mother.a+mother.b == 1 ) continue;     // both parents are hets
        if ( father.a+father.b == mother.a+mother.b ) { trio->err++; continue; }    // mendelian error

        int test_phase = 0; 
        if ( father.a==father.b ) test_phase = 1 + (child.a==father.a);
        else if ( mother.a==mother.b ) test_phase = 1 + (child.b==mother.a);
        if ( trio->prev > 0 )
        {
            if ( trio->prev!=test_phase ) trio->nswitch++;
        }
        trio->ntest++;
        trio->prev = test_phase;
    }
    return NULL;
}

void destroy(void)
{
    int i;
    printf("# This file was produced by: bcftools +trio-switch-rate(%s+htslib-%s)\n", bcftools_version(),hts_version());
    printf("# The command line was:\tbcftools +trio-switch-rate %s", args.argv[0]);
    for (i=1; i<args.argc; i++) printf(" %s",args.argv[i]);
    printf("\n#\n");
    printf("# TRIO\t[2]Father\t[3]Mother\t[4]Child\t[5]nTested\t[6]nMendelian Errors\t[7]nSwitch\t[8]nSwitch (%%)\n");
    for (i=0; i<args.ntrio; i++)
    {
        trio_t *trio = &args.trio[i];
        printf("TRIO\t%s\t%s\t%s\t%d\t%d\t%d\t%.2f\n",
            bcf_hdr_int2id(args.hdr,BCF_DT_SAMPLE,trio->father),
            bcf_hdr_int2id(args.hdr,BCF_DT_SAMPLE,trio->mother),
            bcf_hdr_int2id(args.hdr,BCF_DT_SAMPLE,trio->child),
            trio->ntest, trio->err, trio->nswitch, trio->ntest ? trio->nswitch*100./trio->ntest : 0
        );
        if (args.npop) {
            pop_t *pop = &args.pop[trio->ipop];
            pop->err     += trio->err;
            pop->nswitch += trio->nswitch;
            pop->ntest   += trio->ntest;
            pop->pswitch += trio->ntest ? trio->nswitch*100./trio->ntest : 0;
        }
    }
    printf("# POP\tpopulation or other grouping defined by an optional 7-th column of the PED file\n");
    printf("# POP\t[2]Name\t[3]Number of trios\t[4]avgTested\t[5]avgMendelian Errors\t[6]avgSwitch\t[7]avgSwitch (%%)\n");
    for (i=0; i<args.npop; i++)
    {
        pop_t *pop = &args.pop[i];
        printf("POP\t%s\t%d\t%.0f\t%.0f\t%.0f\t%.2f\n", pop->name,pop->ntrio,
            (float)pop->ntest/pop->ntrio,(float)pop->err/pop->ntrio,(float)pop->nswitch/pop->ntrio,
            pop->pswitch/pop->ntrio);
    }
    for (i=0; i<args.npop; i++) free(args.pop[i].name);
    free(args.pop);
    free(args.trio);
    free(args.gt_arr);
}
