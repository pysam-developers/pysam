/*  plugins/GTisec.c -- collect genotype intersection counts of all possible
                   subsets of the present samples and output in banker's
                   sequence order (in this sequence, the number of contained
                   samples increases monotonically, a property that is e.g.
                   useful for programatically creating plotting files for the
                   R package VennDiagram or the plotting tool circos from the
                   counts, as in the command line tools bankers2VennDiagram and
                   bankers2circos at htpps://github.com/dlaehnemann/bankers2)

    Copyright (C) 2016 Computational Biology of Infection Research,
                       Helmholtz Centre for Infection Research, Braunschweig,
                       Germany

    Author: David Laehnemann <david.laehnemann@helmholtz-hzi.de>

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
#include <htslib/khash.h>
KHASH_MAP_INIT_INT(gts2smps, uint32_t)

#include "bcftools.h"

/*!
 * Flag definitions for args.flag
 */
#define MISSING        (1<<0)
#define VERBOSE        (1<<1)
#define SMPORDER       (1<<2)

typedef struct _args_t
{
    bcf_srs_t *file;    /*! multi-sample VCF file */
    bcf_hdr_t *hdr;     /*! VCF file header */
    FILE *out;          /*! output file pointer */
    int nsmp; /*! number of samples, can be determined from header but is needed in multiple contexts */
    int nsmpp2; /*! 2^(nsmp) (is needed multiple times) */
    int *gt_arr; /*! temporary array, to store GTs of current line/record */
    int ngt_arr; /*! hold the number of current GT array entries */
    uint32_t *bankers; /*! array to store banker's sequence for all possible sample subsets for
                                programmatic indexing into smp_is for output printing, e.g. for three
                                samples A, B and C this would be the following order:
                                [   C,   B,   A,  CB,  CA,  BA, CBA ]
                                [ 100, 010, 001, 110, 101, 011, 111 ]
                                */
    uint64_t *quick; /*! array to store n choose k lookup table of choose() function */
    uint8_t flag; /*! several flags, for positions see above*/
    uint64_t *missing_gts; /*! array to count missing genotypes of each sample */
    uint64_t *smp_is; /*! array to track all possible intersections between
                 samples, with each bit in the index integer belonging to one
                 sample. E.g. for three samples A, B and C, count would be in
                 the following order:
                 [   A,   B,  AB,   C,  AC,  BC, ABC ]
                 [ 001, 010, 011, 100, 101, 110, 111 ]
                 */
}
args_t;

static args_t args;
uint32_t compute_bankers(unsigned long a);

const char *about(void)
{
    return "Count genotype intersections across all possible sample subsets in a vcf file.\n";
}


const char *usage(void)
{
    return
        "\n"
        "About:   Count genotype intersections across all possible sample subsets in a vcf file.\n"
        "Usage:   bcftools +GTisec <multisample.bcf/.vcf.gz> [General Options] -- [Plugin Options] \n"
        "\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -m, --missing                   if set, include count of missing genotypes per sample in output\n"
        "   -v, --verbose                   if set, annotate count rows with corresponding sample subset lists\n"
        "   -H, --human-readable            if set, create human readable output; i.e. sort output by sample and\n"
        "                                   print each subset's intersection count once for each sample contained\n"
        "                                   in the subset; implies verbose output (-v)\n"
        "\n"
        "Example:\n"
        "   bcftools +GTisec in.vcf -- -v # for verbose output\n"
        "   bcftools +GTisec in.vcf -- -H # for human readable output\n"
        "\n";
}


int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    memset(&args,0,sizeof(args_t));
    args.flag = 0;

    static struct option loptions[] =
    {
        {"help",            no_argument,      0,'h'},
        {"missing",         no_argument,      0,'m'},
        {"verbose",         no_argument,      0,'v'},
        {"human-readable",  no_argument,      0,'H'},
        {0,0,0,0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "?mvHh",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'm': args.flag |= MISSING; break;
            case 'v': args.flag |= VERBOSE; break;
            case 'H': args.flag |= ( SMPORDER | VERBOSE ); break;
            case 'h': usage(); break;
            case '?':
            default: error("%s", usage()); break;
        }
    }
    if ( optind != argc )  usage();  // too many files given


    args.hdr = in;

    if ( !bcf_hdr_nsamples(args.hdr) )
    {
        error("No samples in input file.\n");
    }

    args.nsmp = bcf_hdr_nsamples(args.hdr);
    if ( args.nsmp > 32 ) error("Too many samples. A maximum of 32 is supported.\n");
    args.nsmpp2 = pow( 2, args.nsmp);
    args.bankers = (uint32_t*) calloc( args.nsmpp2, sizeof(uint32_t) );
    args.quick = (uint64_t*) calloc((args.nsmp * (args.nsmp + 1)) / 4, sizeof(unsigned long));
    if ( args.flag & MISSING ) args.missing_gts = (uint64_t*) calloc( args.nsmp, sizeof(uint64_t));
    args.smp_is = (uint64_t*) calloc( args.nsmpp2, sizeof(uint64_t));
    if ( bcf_hdr_id2int(args.hdr, BCF_DT_ID, "GT")<0 ) error("[E::%s] GT not present in the header\n", __func__);

    args.gt_arr = NULL;
    args.ngt_arr = 0;

    args.out = stdout;

    /*! Header printing */
    FILE *fp = args.out;
    fprintf(fp, "# This file was produced by bcftools +GTisec (%s+htslib-%s)\n", bcftools_version(), hts_version());
    fprintf(fp, "# The command line was:\tbcftools +GTisec %s ", argv[0]);
    int i;
    for (i=1; i < argc; i++)
    {
        fprintf(fp, " %s", argv[i]);
    }
    fprintf(fp,"\n");
    fprintf(fp,"# This file can be used as input to the subset plotting tools at:\n"
               "#   https://github.com/dlaehnemann/bankers2\n");
    fprintf(fp,"# Genotype intersections across samples:\n");
    fprintf(fp,"@SMPS");
    for (i = args.nsmp-1; i >= 0; i--)
    {
        fprintf(fp," %s", bcf_hdr_int2id(args.hdr, BCF_DT_SAMPLE, i));
    }
    fprintf(fp,"\n");
    if ( args.flag & MISSING )
    {
        if ( args.flag & SMPORDER )
        {
            fprintf(fp, "# The first line of each sample contains its count of missing genotypes, with a '-' appended\n"
                        "#   to the sample name.\n");
        }
        else
        {
            fprintf(fp, "# The first %i lines contain the counts for missing values of each sample in the order provided\n"
                        "#   in the SMPS-line above. Intersection counts only start afterwards.\n", args.nsmp);
        }
    }
    if ( args.flag & SMPORDER )
    {
        fprintf(fp, "# Human readable output (-H) was requested. Subset intersection counts are therefore sorted by\n"
                    "#   sample and repeated for each contained sample. For each sample, counts are in banker's \n"
                    "#   sequence order regarding all other samples.\n");
    }
    else
    {
        fprintf(fp, "# Subset intersection counts are in global banker's sequence order.\n");
        if ( args.nsmp > 2 )
        {
            fprintf(fp, "#   After exclusive sample counts in order of the SMPS-line, banker's sequence continues with:\n"
                        "#   %s,%s   %s,%s   ...\n",
                            bcf_hdr_int2id(in, BCF_DT_SAMPLE, args.nsmp-1 ),
                            bcf_hdr_int2id(in, BCF_DT_SAMPLE, args.nsmp-2 ),
                            bcf_hdr_int2id(in, BCF_DT_SAMPLE, args.nsmp-1 ),
                            bcf_hdr_int2id(in, BCF_DT_SAMPLE, args.nsmp-3 )
                            );
        }
    }
    if (args.flag & VERBOSE )
    {
        fprintf(fp,"# [1] Number of shared non-ref genotypes \t[2] Samples sharing non-ref genotype (GT)\n");
    }
    else
    {
        fprintf(fp,"# [1] Number of shared non-ref genotypes\n");
    }

    /* Compute banker's sequence for following printing by sample and
     * with increasing subset size.
     */
    uint32_t j;
    for ( j = 0; j < args.nsmpp2; j++ )
    {
        args.bankers[j] = compute_bankers(j);
    }

    return 1;
}


/* ADAPTED CODE FROM CORIN LAWSON (START)
 * https://github.com/au-phiware/bankers/blob/master/c/bankers.c
 * who implemented ideas of Eric Burnett:
 * http://www.thelowlyprogrammer.com/2010/04/indexing-and-enumerating-subsets-of.html
 */

/*
 * Compute the binomial coefficient of `n choose k'.
 * Use the fact that binom(n, k) = binom(n, n - k).
 * Use a lookup table (triangle, actually) for speed.
 * Otherwise it's dumb (heart) recursion.
 * Added relative to Corin Lawson:
 * * Passing in of sample number through pointer to args struct
 * * Make quick lookup table external to keep it persistent with clean allocation
 *   and freeing
 */
uint64_t choose(unsigned int n, unsigned int k) {
    if (n == 0)
        return 0;
    if (n == k || k == 0)
        return 1;
    if (k > n / 2)
        k = n - k;

    unsigned int i = (n * (n - 1)) / 4 + k - 1;
    if (args.quick[i] == 0)
        args.quick[i] = choose(n - 1, k - 1) + choose(n - 1, k);

    return args.quick[i];
}

/*
 * Returns the Banker's number at the specified position, a.
 * Derived from the recursive bit flip method.
 * Added relative to Corin Lawson:
 * * Uses same lookup table solution as choose function, just
 *   maintained externally to persist across separate function calls.
 * * Uses bitwise symmetry of banker's sequence to use bitwise inversion
 *   instead of recursive bit flip for second half of sequence.
 */
uint32_t compute_bankers(unsigned long a)
{
    if (a == 0)
        return 0;

    if ( args.bankers[a] == 0 )
    {
        if ( a >= (args.nsmpp2 / 2) )
            return args.bankers[a] = ( compute_bankers(args.nsmpp2 - (a+1)) ^ (args.nsmpp2 - 1) ); // use bitwise symmetry of bankers sequence
        unsigned int c = 0;
        uint32_t n = args.nsmp;
        uint64_t e = a, binom;
        binom = choose(n, c);
        do {
            e -= binom;
        } while ((binom = choose(n, ++c)) <= e);

        do {
            if (e == 0 || (binom = choose(n - 1, c - 1)) > e)
                c--, args.bankers[a] |= 1;
            else
                e -= binom;
        } while (--n && c && ((args.bankers[a] <<= 1) || 1));
        args.bankers[a] <<= n;
    }

    return args.bankers[a];
}

// ADAPTED CODE FROM CORIN LAWSON END


/*
 * GT field (genotype) comparison function.
 */
bcf1_t *process(bcf1_t *rec)
{
    uint64_t i;
    bcf_unpack(rec, BCF_UN_FMT); // unpack the Format fields, including the GT field
    int gte_smp = 0; // number GT array entries per sample (should be 2, one entry per allele)
    if ( (gte_smp = bcf_get_genotypes(args.hdr, rec, &(args.gt_arr), &(args.ngt_arr) ) ) <= 0 )
    {
        error("GT not present at %s: %d\n", args.hdr->id[BCF_DT_CTG][rec->rid].key, rec->pos+1);
    }

    gte_smp /= args.nsmp; // divide total number of genotypes array entries (= args.ngt_arr) by number of samples
    int ret;

    // stick all genotypes in a hash as keys and store up to 32 samples in a corresponding flag as its value
    khiter_t bucket;
    khash_t(gts2smps) *gts = kh_init(gts2smps); // create hash
    for ( i = 0; i < args.nsmp; i++ )
    {
        int *gt_ptr = args.gt_arr + gte_smp * i;

        if (bcf_gt_is_missing(gt_ptr[0]) || ( gte_smp == 2 && bcf_gt_is_missing(gt_ptr[1]) ) )
        {
            if ( args.flag & MISSING ) args.missing_gts[i]++; // count missing genotypes, if requested
            continue; // don't do anything else for missing genotypes, their "sharing" gives no info...
        }

        int a = bcf_gt_allele(gt_ptr[0]);
        int b;
        if ( gte_smp == 2 ) // two entries available per sample, padded with missing values for haploid genotypes
        {
            b = bcf_gt_allele(gt_ptr[1]);
        }
        else if (gte_smp == 1 ) // use missing value for second entry in hash key generation below, if only one is available
        {
            b = bcf_gt_allele(bcf_int32_vector_end);
        }
        else
        {
            error("gtisec does not support ploidy higher than 2.\n");
        }

        int idx = bcf_alleles2gt(a,b); // generate genotype specific hash key

        bucket = kh_get(gts2smps, gts, idx); // get the genotype's hash bucket

        if ( bucket == kh_end(gts) ) { // means that key does not exist
            bucket = kh_put(gts2smps, gts, idx, &ret); // create bucket with genotype index as key and return its iterator
            kh_val(gts, bucket) = 0; // initialize the bucket with all sample bits unset
        }
        kh_value(gts, bucket) |= (1<<i); // set the sample's bit to 1 in this genotype's bucket
    }

    // iterate over genotypes and for each genotype increment the appropriate smp_is entry
    for ( bucket = kh_begin(gts); bucket != kh_end(gts); ++bucket ) // iterate over all genotypes at this position
    {
        if ( kh_exist(gts, bucket) ) // for existing genotype buckets
        {
            uint32_t s = kh_val(gts, bucket); // get the 32 bit flag
            args.smp_is[s]++; // add to the corresponding subset
        }
    }
    kh_destroy(gts2smps, gts); // destroy hash

    return NULL;
}

void destroy(void)
{
    int32_t i;
    int s;

    FILE *fp = args.out;

    /* Printing to File */
    if ( args.flag & SMPORDER )
    {
        /* Iterate over samples, printing out all subsets including
         * the current sample, with the current sample first. This
         * includes multiple printouts of the same sample but makes
         * output more readable and is also needed for circos files
         * printing.
         */
        for ( s = args.nsmp-1; s >= 0; s--)
        {
            if ( args.flag & MISSING ) // if missing genotype counts are requested, print them to standard output
            {
                fprintf(fp, "%"PRIu64"\t%s-\n", args.missing_gts[s], bcf_hdr_int2id(args.hdr, BCF_DT_SAMPLE, s));
            }
            for ( i = 1; i < args.nsmpp2; i++ )
            {
                if ( (args.bankers[i]>>s) & 1 )
                {
                    fprintf(fp, "%"PRIu64"\t", args.smp_is[ args.bankers[i] ]); // print out count of genotypes shared by samples in current banker's sequence position
                    int j;
                    /* Print sample list */
                    fprintf(fp, "%s", bcf_hdr_int2id(args.hdr, BCF_DT_SAMPLE, s)); // print current sample first
                    for ( j = args.nsmp-1; j >= 0; j-- )
                    {
                        if ( (args.bankers[i] ^ (1<<s)) & (1<<j) ) // exclude current sample from printing again
                        {
                            fprintf(fp, ",%s", bcf_hdr_int2id(args.hdr, BCF_DT_SAMPLE, j ) ); // print out sample list, starting with our current major sample
                        }
                    }
                    fprintf(fp, "\n" );
                }
            }
        }
    }
    else if ( args.flag & VERBOSE )
    {
        if ( args.flag & MISSING ) // if missing genotype counts are requested, print them to standard output
        {
            for ( s = args.nsmp-1; s >= 0; s--)
            {
                fprintf(fp, "%"PRIu64"\t%s-\n", args.missing_gts[s], bcf_hdr_int2id(args.hdr, BCF_DT_SAMPLE, s));
            }
        }
        for ( i = 1; i < args.nsmpp2; i++ )
        {
            fprintf(fp, "%"PRIu64"\t", args.smp_is[ args.bankers[i] ]); // print out count of genotypes shared by samples in current banker's sequence position
            int j = 0;
            for ( s = args.nsmp-1; s >= 0; s--)
            {
               if ( (args.bankers[i]>>s) & 1 )
               {
                   fprintf(fp, "%s%s", j ? "," : "", bcf_hdr_int2id(args.hdr, BCF_DT_SAMPLE, s) ); // samples in specified order
                   j = 1;
               }
            }
            fprintf(fp, "\n" );
        }
    }
    else
    {
        if ( args.flag & MISSING ) // if missing genotype counts are requested, print them to standard output
        {
            for ( s = args.nsmp-1; s >= 0; s--)
            {
                fprintf(fp, "%"PRIu64"\n", args.missing_gts[s]);
            }
        }
        for ( i = 1; i < args.nsmpp2; i++ )
        {
            fprintf(fp, "%"PRIu64"\n", args.smp_is[ args.bankers[i] ]); // print out count of genotypes shared by samples in current banker's sequence position
        }
    }
    fclose(fp);

    /* freeing up args */
    free(args.gt_arr);
    free(args.bankers);
    free(args.quick);
    if (args.flag & MISSING) free(args.missing_gts);
    free(args.smp_is);
}
