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
#include <htslib/kfunc.h>
#include <inttypes.h>
#include "bcftools.h"
#include "convert.h"

typedef struct
{
    int smpl,ctrl;      // VCF sample index
    const char *smpl_name, *ctrl_name;
}
pair_t;

typedef struct
{
    bcf_hdr_t *hdr;
    pair_t *pair;
    int npair, mpair, min_dp, min_alt_dp;
    int32_t *ad_arr;
    int mad_arr;
    double th;
    convert_t *convert;
    kstring_t str;
    uint64_t nsite,ncmp;
}
args_t;

args_t args;

const char *about(void)
{
    return "Find positions with wildly varying ALT allele frequency (Fisher test on FMT/AD).\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Find positions with wildly varying ALT allele frequency (Fisher test on FMT/AD).\n"
        "Usage: bcftools +ad-bias [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -a, --min-alt-dp <int>      Minimum required alternate allele depth [1]\n"
        "   -d, --min-dp <int>          Minimum required depth [0]\n"
        "   -f, --format <string>       Optional tags to append to output (`bcftools query` style of format)\n"
        "   -s, --samples <file>        List of sample pairs, one tab-delimited pair per line\n"
        "   -t, --threshold <float>     Output only hits with p-value smaller than <float> [1e-3]\n"
        "\n"
        "Example:\n"
        "   bcftools +ad-bias file.bcf -- -t 1e-3 -s samples.txt\n"
        "\n";
}

void parse_samples(args_t *args, char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Could not read: %s\n", fname);

    kstring_t str = {0,0,0};
    if ( hts_getline(fp, KS_SEP_LINE, &str) <= 0 ) error("Empty file: %s\n", fname);

    int moff = 0, *off = NULL;
    do
    {
        // HPSI0513i-veqz_6    HPSI0513pf-veqz
        int ncols = ksplit_core(str.s,'\t',&moff,&off);
        if ( ncols<2 ) error("Could not parse the sample file: %s\n", str.s);

        int smpl = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[0]]);
        if ( smpl<0 ) continue;
        int ctrl = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[1]]);
        if ( ctrl<0 ) continue;

        args->npair++;
        hts_expand0(pair_t,args->npair,args->mpair,args->pair);
        pair_t *pair = &args->pair[args->npair-1];
        pair->ctrl = ctrl;
        pair->smpl = smpl;
        pair->smpl_name = bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,pair->smpl);
        pair->ctrl_name = bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,pair->ctrl);
    } while ( hts_getline(fp, KS_SEP_LINE, &str)>=0 );

    free(str.s);
    free(off);
    hts_close(fp);
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    memset(&args,0,sizeof(args_t));
    args.hdr = in;
    args.th  = 1e-3;
    args.min_alt_dp = 1;
    char *fname = NULL, *format = NULL;
    static struct option loptions[] =
    {
        {"min-dp",required_argument,NULL,'d'},
        {"min-alt-dp",required_argument,NULL,'a'},
        {"format",required_argument,NULL,'f'},
        {"samples",required_argument,NULL,'s'},
        {"threshold",required_argument,NULL,'t'},
        {NULL,0,NULL,0}
    };
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "?hs:t:f:d:a:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'a':
                args.min_alt_dp = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse: -a %s\n", optarg);
                break;
            case 'd':
                args.min_dp = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse: -d %s\n", optarg);
                break;
            case 't':
                args.th = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: -t %s\n", optarg);
                break;
            case 's': fname = optarg; break;
            case 'f': format = optarg; break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }
    if ( !fname ) error("Expected the -s option\n");
    parse_samples(&args, fname);
    if ( format ) args.convert = convert_init(args.hdr, NULL, 0, format);
    printf("# This file was produced by: bcftools +ad-bias(%s+htslib-%s)\n", bcftools_version(),hts_version());
    printf("# The command line was:\tbcftools +ad-bias %s", argv[0]);
    for (c=1; c<argc; c++) printf(" %s",argv[c]);
    printf("\n#\n");
    printf("# FT, Fisher Test\t[2]Sample\t[3]Control\t[4]Chrom\t[5]Pos\t[6]smpl.nREF\t[7]smpl.nALT\t[8]ctrl.nREF\t[9]ctrl.nALT\t[10]P-value");
    if ( format ) printf("\t[11-]User data: %s", format);
    printf("\n");
    return 1;
}

bcf1_t *process(bcf1_t *rec)
{
    int nad = bcf_get_format_int32(args.hdr, rec, "AD", &args.ad_arr, &args.mad_arr);
    if ( nad<0 ) return NULL;
    nad /= bcf_hdr_nsamples(args.hdr);
    
    if ( args.convert ) convert_line(args.convert, rec, &args.str);
    args.nsite++;

    int i;
    for (i=0; i<args.npair; i++)
    {
        pair_t *pair = &args.pair[i];
        int32_t *aptr = args.ad_arr + nad*pair->smpl;
        int32_t *bptr = args.ad_arr + nad*pair->ctrl;

        if ( aptr[0]==bcf_int32_missing ) continue;
        if ( bptr[0]==bcf_int32_missing ) continue;
        if ( aptr[0]+aptr[1] < args.min_dp ) continue;
        if ( bptr[0]+bptr[1] < args.min_dp ) continue;
        if ( aptr[1] < args.min_alt_dp && bptr[1] < args.min_alt_dp ) continue;

        args.ncmp++;

        int n11 = aptr[0], n12 = aptr[1];
        int n21 = bptr[0], n22 = bptr[1];
        double left, right, fisher;
        kt_fisher_exact(n11,n12,n21,n22, &left,&right,&fisher);
        if ( fisher >= args.th ) continue;

        printf("FT\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%e",
            pair->smpl_name,pair->ctrl_name,
            bcf_hdr_id2name(args.hdr,rec->rid), rec->pos+1,
            n11,n12,n21,n22, fisher
            );
        if ( args.convert ) printf("\t%s", args.str.s);
        printf("\n");
    }
    return NULL;
}

void destroy(void)
{
    printf("# SN, Summary Numbers\t[2]Number of Pairs\t[3]Number of Sites\t[4]Number of comparisons\t[5]P-value output threshold\n");
    printf("SN\t%d\t%"PRId64"\t%"PRId64"\t%e\n",args.npair,args.nsite,args.ncmp,args.th);
    if (args.convert) convert_destroy(args.convert);
    free(args.str.s);
    free(args.pair);
    free(args.ad_arr);
}
