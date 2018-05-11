/*  plugins/GTsubset.c -- output only positions where the selected samples exclusively
                  share a genotype, i.e. all selected samples must have the same
                  genotype (including both alleles) and none of the unselected
                  samples can have the same genotype

    Copyright (C) 2016 Computational Biology of Infection Research,
                       Helmholtz Centre for Infection Research, Braunschweig,
                       Germany

    Author: David Laehnemann <david.laehnemann@hhu.de>

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
#include <getopt.h>
#include <math.h>
#include <unistd.h>
#include <inttypes.h>

#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>

#include "bcftools.h"

typedef struct _args_t
{
    bcf_hdr_t *hdr;     /*! VCF file header */
    int *gt_arr;        /*! temporary array, to store GTs of current line/record */
    int ngt_arr;        /*! hold the number of current GT array entries */
    int nsmp;           /*! number of samples, can be determined from header but is needed in multiple contexts */
    int n_sel_smps;     /*! number of selected samples who should exclusively share genotypes */
    int *selected_smps; /*! pointer to start of array containing 1 at indices corresponding to selected samples in header dict and 0 at others*/
}
args_t;

static args_t args;

const char *about(void)
{
    return "Output only sites where the requested samples all exclusively share a genotype (GT).\n";
}


const char *usage(void)
{
    return
        "\n"
        "About:   Output only sites where the requested samples all exclusively share a genotype (GT), i.e.\n"
        "         all selected samples must have the same GT, while non of the others can have it.\n"
        "Usage:   bcftools +GTsubset <multisample.bcf/.vcf.gz> [General Options] -- [Plugin Options] \n"
        "\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "  -s,--sample-list     comma-separated list of samples; only those sites where all of these\n"
        "                       samples exclusively share their genotype are given as output\n"
        "\n"
        "Example:\n"
        "   bcftools +GTsubset in.vcf -- -s SMP1,SMP2 \n"
        "\n";
}


int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    memset(&args,0,sizeof(args_t));

    int i;

    static struct option loptions[] =
    {
        {"help",            no_argument,       0,'h'},
        {"sample-list",     required_argument, 0,'s'},
        {0,0,0,0}
    };

    char **smps_strs = NULL;

    int c;
    while ((c = getopt_long(argc, argv, "?s:h",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 's': smps_strs = hts_readlist(optarg,0,&(args.n_sel_smps));
                      if ( args.n_sel_smps == 0 )
                      {
                          fprintf(stderr, "Sample specification not valid.\n");
                          error("%s", usage());
                      }
                      break;
            case 'h': usage(); break;
            case '?':
            default: error("%s", usage()); break;
        }
    }
    if ( optind != argc )  usage();  // too many files given

    args.hdr = bcf_hdr_dup(in);

    // Samples parsing from header and input option
    if ( !bcf_hdr_nsamples(args.hdr) )
    {
        error("No samples in input file.\n");
    }
    args.nsmp = bcf_hdr_nsamples(args.hdr);
    args.selected_smps = (int*) calloc(args.nsmp,sizeof(int));
    for ( i = 0; i < args.n_sel_smps; i++ )
    {
        int ind = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, smps_strs[i]);
        if ( ind == -1 )
        {
            error("Sample '%s' not in input vcf file.\n", smps_strs[i]);
        } else {
            args.selected_smps[ind] = 1;
        }
        free(smps_strs[i]);
    }
    free(smps_strs);

    /*
    fprintf(stderr, "Selected samples array:[");
    for (i=0;i<args.nsmp;i++)
    {
        fprintf(stderr, " %i", args.selected_smps[i]);
    }
    fprintf(stderr, " ]\n");
    */

    if ( bcf_hdr_id2int(args.hdr, BCF_DT_ID, "GT")<0 ) error("[E::%s] GT not present in the header\n", __func__);

    args.gt_arr = NULL;

    return 0;
}


/*
 * GT field (genotype) comparison function.
 */
bcf1_t *process(bcf1_t *rec)
{
    uint64_t i;
    bcf_unpack(rec, BCF_UN_FMT); // unpack the Format fields, including the GT field
    int gte_smp = 0; // number GT array entries per sample (should be 2, one entry per allele)
    args.ngt_arr = 0;        /*! hold the number of current GT array entries */
    if ( (gte_smp = bcf_get_genotypes(args.hdr, rec, &(args.gt_arr), &(args.ngt_arr) ) ) <= 0 )
    {
        error("GT not present at %s: %d\n", args.hdr->id[BCF_DT_CTG][rec->rid].key, rec->pos+1);
    }

    gte_smp /= args.nsmp; // divide total number of genotypes array entries (= args.ngt_arr) by number of samples

    // initialize with missing genotype
    int a1 = 0;
    int a2 = 0;

    // initialize with first selected sample genotype that is not missing
    int gt = -1;
    while ( (a1 == 0) || (a2 == 0) )
    {
       gt++;
       if (gt == args.nsmp) break;
       if (args.selected_smps[gt] == 0) continue;
       a1 = (args.gt_arr + gte_smp * gt)[0];
       if ( gte_smp == 2 ) a2 = (args.gt_arr + gte_smp * gt)[1];
       else if ( gte_smp == 1 ) a2 = bcf_int32_vector_end;
       else error("GTsubset does not support ploidy higher than 2.\n");
    }
//    fprintf(stderr, "a1: %i  a2: %i\n", a1, a2);

    // check all genotypes if they match (for included samples) or disagree (for samples not included)
    gt = 0;
    for ( i = 0; i < args.nsmp; i++ )
    {
        int *gt_ptr = args.gt_arr + gte_smp * i;

        int b1 = gt_ptr[0];
        int b2;
        if ( gte_smp == 2 ) // two entries available per sample, padded with missing values for haploid genotypes
        {
            b2 = gt_ptr[1];
        }
        else if (gte_smp == 1 ) // use vector end value for second entry, if only one is available
        {
            b2 = bcf_int32_vector_end;
        }
        else
        {
            error("GTsubset does not support ploidy higher than 2.\n");
        }

 //      fprintf(stderr, "b1: %i  b2: %i\n", b1, b2);
        /* missing genotypes are counted as always passing, as they neither
         * mismatch the initial selected genotype for a selected sample, nor
         * do they match the initial selected genotype for an excluded sample's
         * genotype */
        if ( (b1 == 0) || (b2 == 0) )
        {
            gt++;
//            fprintf(stderr, "missing => pass\n");
            continue;
        }
        else if ( args.selected_smps[i] == 1 )
        {
            if ( (b1 == a1) && (b2 == a2) )
            {
                gt++;
//                fprintf(stderr, "match => pass\n");
                continue;
            }
            else
            {
//                fprintf(stderr, "no match => fail\n");
                break;
            }
        }
        else if ( args.selected_smps[i] == 0 )
        {
            if ( (b1 != a1 ) || (b2 != a2) )
            {
                gt++;
 //               fprintf(stderr, "no match => pass\n");
                continue;
            }
            else
            {
//                fprintf(stderr, "match => fail\n");
                break;
            }
        }
    }
    if ( gt == args.nsmp )
    {
        return rec;
    }
    else
    {
        return NULL;
    }
}

void destroy(void)
{
    /* freeing up args */
    bcf_hdr_destroy(args.hdr);
    free(args.gt_arr);
    free(args.selected_smps);
}
