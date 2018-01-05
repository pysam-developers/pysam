/* 
    Copyright (C) 2014 Genome Research Ltd.

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
#include <strings.h>
#include <getopt.h>
#include <stdarg.h>
#include <ctype.h>
#include <htslib/vcf.h>
#include <htslib/kseq.h>
#include "bcftools.h"
#include "ploidy.h"

static bcf_hdr_t *in_hdr = NULL, *out_hdr = NULL;
static int *sample2sex = NULL;
static int n_sample = 0, nsex = 0, *sex2ploidy = NULL;
static int32_t ngt_arr = 0, *gt_arr = NULL, *gt_arr2 = NULL, ngt_arr2 = 0;
static ploidy_t *ploidy = NULL;

const char *about(void)
{
    return "Fix ploidy.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Fix ploidy\n"
        "Usage: bcftools +fixploidy [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -p, --ploidy <file>   space/tab-delimited list of CHROM,FROM,TO,SEX,PLOIDY\n"
        "   -s, --sex <file>      list of samples, \"NAME SEX\"\n"
        "   -t, --tags <list>     VCF tags to fix [GT]\n"
        "\n"
        "Example:\n"
        "   # Default ploidy, if -p not given. Unlisted regions have ploidy 2\n"
        "   X 1 60000 M 1\n"
        "   X 2699521 154931043 M 1\n"
        "   Y 1 59373566 M 1\n"
        "   Y 1 59373566 F 0\n"
        "   MT 1 16569 M 1\n"
        "   MT 1 16569 F 1\n"
        "   \n"
        "   # Example of -s file, sex of unlisted samples is \"F\"\n"
        "   sampleName1 M\n"
        "   \n"
        "   bcftools +fixploidy in.vcf -- -s samples.txt\n"
        "\n";
}

void set_samples(char *fname, bcf_hdr_t *hdr, ploidy_t *ploidy, int *sample2sex)
{
    kstring_t tmp = {0,0,0};

    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Could not read: %s\n", fname);
    while ( hts_getline(fp, KS_SEP_LINE, &tmp) > 0 )
    {
        char *ss = tmp.s;
        while ( *ss && isspace(*ss) ) ss++;
        if ( !*ss ) error("Could not parse: %s\n", tmp.s);
        if ( *ss=='#' ) continue;
        char *se = ss;
        while ( *se && !isspace(*se) ) se++;
        char x = *se; *se = 0;

        int ismpl = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, ss);
        if ( ismpl < 0 ) { fprintf(stderr,"Warning: No such sample in the VCF: %s\n",ss); continue; }

        *se = x;
        ss = se+1;
        while ( *ss && isspace(*ss) ) ss++;
        if ( !*ss )  error("Could not parse: %s\n", tmp.s);
        se = ss;
        while ( *se && !isspace(*se) ) se++;
        if ( se==ss ) error("Could not parse: %s\n", tmp.s);

        sample2sex[ismpl] = ploidy_add_sex(ploidy, ss);
    }
    if ( hts_close(fp) ) error("Close failed: %s\n", fname);
    free(tmp.s);
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    int c;
    char *tags_str = "GT";
    char *ploidy_fname = NULL, *sex_fname = NULL;

    static struct option loptions[] =
    {
        {"ploidy",1,0,'p'},
        {"sex",1,0,'s'},
        {"tags",1,0,'t'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "?ht:s:p:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'p': ploidy_fname = optarg; break;
            case 's': sex_fname = optarg; break;
            case 't': tags_str = optarg; break;
            case 'h':
            case '?':
            default: error("%s", usage()); break;
        }
    }
    if ( strcasecmp("GT",tags_str) ) error("Only -t GT is currently supported, sorry\n");

    n_sample   = bcf_hdr_nsamples(in);
    sample2sex = (int*) calloc(n_sample,sizeof(int));
    in_hdr     = in;
    out_hdr    = out;

    if ( ploidy_fname )
        ploidy = ploidy_init(ploidy_fname, 2);
    else
    {
        ploidy = ploidy_init_string(
                "X 1 60000 M 1\n"
                "X 2699521 154931043 M 1\n"
                "Y 1 59373566 M 1\n"
                "Y 1 59373566 F 0\n"
                "MT 1 16569 M 1\n"
                "MT 1 16569 F 1\n", 2);
    }
    if ( !ploidy ) return -1;

    // add default sex in case it was not included
    int i, dflt_sex_id = ploidy_add_sex(ploidy, "F");
    for (i=0; i<n_sample; i++) sample2sex[i] = dflt_sex_id; // by default all are F
    if ( sex_fname ) set_samples(sex_fname, in, ploidy, sample2sex);
    nsex = ploidy_nsex(ploidy);
    sex2ploidy = (int*) malloc(sizeof(int)*nsex);

    return 0;
}


bcf1_t *process(bcf1_t *rec)
{
    int i,j, max_ploidy;

    int ngts = bcf_get_genotypes(in_hdr, rec, &gt_arr, &ngt_arr);
    if ( ngts<0 )
        return rec;     // GT field not present

    if ( ngts % n_sample )
        error("Error at %s:%d: wrong number of GT fields\n",bcf_seqname(in_hdr,rec),rec->pos+1);

    ploidy_query(ploidy, (char*)bcf_seqname(in_hdr,rec), rec->pos, sex2ploidy,NULL,&max_ploidy);

    ngts /= n_sample;
    if ( ngts < max_ploidy )
    {
        hts_expand(int32_t,max_ploidy*n_sample,ngt_arr2,gt_arr2);
        for (i=0; i<n_sample; i++)
        {
            int ploidy = sex2ploidy[ sample2sex[i] ];
            int32_t *src = &gt_arr[i*ngts];
            int32_t *dst = &gt_arr2[i*max_ploidy];
            j = 0;
            if ( !ploidy ) { dst[j] = bcf_gt_missing; j++; }
            else
                while ( j<ngts && j<ploidy && src[j]!=bcf_int32_vector_end ) { dst[j] = src[j]; j++; }
            assert( j );
            while ( j<ploidy ) { dst[j] = dst[j-1]; j++; } // expand "." to "./." and "0" to "0/0"
            while ( j<max_ploidy ) { dst[j] = bcf_int32_vector_end; j++; }
        }
        if ( bcf_update_genotypes(out_hdr,rec,gt_arr2,n_sample*max_ploidy) )
            error("Could not update GT field at %s:%d\n", bcf_seqname(in_hdr,rec),rec->pos+1);
    }
    else if ( ngts!=1 || max_ploidy!=1 )
    {
        for (i=0; i<n_sample; i++)
        {
            int ploidy = sex2ploidy[ sample2sex[i] ];
            int32_t *gts = &gt_arr[i*ngts];
            j = 0;
            if ( !ploidy ) { gts[j] = bcf_gt_missing; j++; }
            else 
                while ( j<ngts && j<ploidy && gts[j]!=bcf_int32_vector_end ) j++;
            assert( j );
            while ( j<ploidy ) { gts[j] = gts[j-1]; j++; } // expand "." to "./." and "0" to "0/0"
            while ( j<ngts ) { gts[j] = bcf_int32_vector_end; j++; }
        }
        if ( bcf_update_genotypes(out_hdr,rec,gt_arr,n_sample*ngts) )
            error("Could not update GT field at %s:%d\n", bcf_seqname(in_hdr,rec),rec->pos+1);
    }
    return rec;
}


void destroy(void)
{
    free(gt_arr);
    free(gt_arr2);
    free(sample2sex);
    free(sex2ploidy);
    ploidy_destroy(ploidy);
}


