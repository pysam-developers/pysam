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
#include <inttypes.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include "bcftools.h"
#include "bin.h"

typedef struct
{
    char *af_tag;
    bcf_hdr_t *hdr;
    int32_t *gt, ngt, naf;
    float *af, list_min, list_max;
    bin_t *dev_bins, *prob_bins;
    uint64_t *dev_dist, *prob_dist;
}
args_t;

args_t *args;

const char *about(void)
{
    return "AF and GT probability distribution stats.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Collect AF deviation stats and GT probability distribution\n"
        "       given AF and assuming HWE\n"
        "Usage: bcftools +af-dist [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -d, --dev-bins <list>       AF deviation bins\n"
        "   -l, --list <min,max>        list genotypes from the given bin (for debugging)\n"
        "   -p, --prob-bins <list>      probability distribution bins\n"
        "   -t, --af-tag <tag>          VCF INFO tag to use [AF]\n"
        "\n"
        "Default binning:\n"
        "   -d: 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1\n"
        "   -p: 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1\n"
        "Example:\n"
        "   bcftools +af-tag file.bcf -- -t EUR_AF -p bins.txt\n"
        "\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    args = (args_t*) calloc(1,sizeof(args_t));
    char *dev_bins  = "0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1";
    char *prob_bins = "0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1";
    args->hdr = in;
    args->af_tag = "AF";
    args->list_min = -1;
    static struct option loptions[] =
    {
        {"list",required_argument,NULL,'l'},
        {"dev-bins",required_argument,NULL,'d'},
        {"prob-bins",required_argument,NULL,'p'},
        {"af-tag",required_argument,NULL,'t'},
        {NULL,0,NULL,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "?ht:d:p:l:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'l': 
            {
                char *a,*b;
                args->list_min = strtod(optarg,&a);
                if ( a==optarg || *a!=',' ) error("Could not parse: --list %s\n", optarg);
                args->list_max = strtod(a+1,&b);
                if ( a+1==b || *b ) error("Could not parse: --list %s\n", optarg);
                break;
            }
            case 'd': dev_bins = optarg; break;
            case 'p': prob_bins = optarg; break;
            case 't': args->af_tag = optarg; break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }

    args->dev_bins = bin_init(dev_bins,0,1);
    int nbins = bin_get_size(args->dev_bins);
    args->dev_dist = (uint64_t*)calloc(nbins,sizeof(*args->dev_dist));

    args->prob_bins = bin_init(prob_bins,0,1);
    nbins = bin_get_size(args->prob_bins);
    args->prob_dist = (uint64_t*)calloc(nbins,sizeof(*args->prob_dist));

    printf("# This file was produced by: bcftools +af-dist(%s+htslib-%s)\n", bcftools_version(),hts_version());
    printf("# The command line was:\tbcftools +af-dist %s", argv[0]);
    for (c=1; c<argc; c++) printf(" %s",argv[c]);
    printf("\n#\n");

    if ( args->list_min!=-1 )
        printf("# GT, genotypes with P(AF) in [%f,%f]; [2]Chromosome\t[3]Position[4]Sample\t[5]Genotype\t[6]AF-based probability\n",args->list_min,args->list_max);

    return 1;
}

bcf1_t *process(bcf1_t *rec)
{
    int naf = bcf_get_info_float(args->hdr,rec,args->af_tag,&args->af,&args->naf);
    if ( naf<=0 ) return NULL;
    float af = args->af[0];

    float pRA = 2*af*(1-af);
    float pAA = af*af;
    int iRA = bin_get_idx(args->prob_bins,pRA);
    int iAA = bin_get_idx(args->prob_bins,pAA);

    int list_RA = args->list_min==-1 || pRA < args->list_min || pRA > args->list_max ? 0 : 1;
    int list_AA = args->list_min==-1 || pAA < args->list_min || pAA > args->list_max ? 0 : 1;
    const char *chr = bcf_seqname(args->hdr,rec);

    int ngt = bcf_get_genotypes(args->hdr, rec, &args->gt, &args->ngt);
    int i, j, nsmpl = bcf_hdr_nsamples(args->hdr);
    int nals = 0, nalt = 0;
    ngt /= nsmpl;
    for (i=0; i<nsmpl; i++)
    {
        int32_t *ptr = args->gt + i*ngt;
        int dosage = 0;
        for (j=0; j<ngt; j++)
        {
            if ( bcf_gt_is_missing(ptr[j]) ) break;
            if ( ptr[j]==bcf_int32_vector_end ) break;
            if ( bcf_gt_allele(ptr[j])==1 ) dosage++;
        }
        if ( j!=ngt ) continue;

        nals += j;
        nalt += dosage;

        if ( dosage==1 )
        {
            args->prob_dist[iRA]++;
            if ( list_RA ) printf("GT\t%s\t%d\t%s\t1\t%f\n",chr,rec->pos+1,args->hdr->samples[i],pRA);
        }
        else if ( dosage==2 )
        {
            args->prob_dist[iAA]++;
            if ( list_AA ) printf("GT\t%s\t%d\t%s\t2\t%f\n",chr,rec->pos+1,args->hdr->samples[i],pAA);
        }
    }

    if ( nals && (nalt || af) )
    {
        float af_dev = fabs(af - (float)nalt/nals);
        int iAF = bin_get_idx(args->dev_bins,af_dev);
        args->dev_dist[iAF]++;
    }

    return NULL;
}

void destroy(void)
{
    printf("# PROB_DIST, genotype probability distribution, assumes HWE\n");
    int i, n;
    n = bin_get_size(args->prob_bins);
    for (i=0; i<n-1; i++)
    {
        float min = bin_get_value(args->prob_bins,i);
        float max = bin_get_value(args->prob_bins,i+1);
        printf("PROB_DIST\t%f\t%f\t%"PRId64"\n", min,max,args->prob_dist[i]);
    }
    printf("# DEV_DIST, distribution of AF deviation, based on %s and INFO/AN, AC calculated on the fly\n", args->af_tag);
    n = bin_get_size(args->dev_bins);
    for (i=0; i<n-1; i++)
    {
        float min = bin_get_value(args->dev_bins,i);
        float max = bin_get_value(args->dev_bins,i+1);
        printf("DEV_DIST\t%f\t%f\t%"PRId64"\n", min,max,args->dev_dist[i]);
    }
    bin_destroy(args->dev_bins);
    bin_destroy(args->prob_bins);
    free(args->dev_dist);
    free(args->prob_dist);
    free(args->gt);
    free(args->af);
    free(args);
}


