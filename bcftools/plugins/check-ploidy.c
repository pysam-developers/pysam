/* 
    Copyright (C) 2017 Genome Research Ltd.

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
#include <stdint.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <htslib/kseq.h>
#include <inttypes.h>
#include <unistd.h>
#include "bcftools.h"

typedef struct
{
    char *sample;
    int beg,end,ploidy;
}
dat_t;

typedef struct
{
    int argc;
    char **argv;
    int rid, gt_id, ndat;
    dat_t *dat;
    bcf_hdr_t *hdr;
}
args_t;

static args_t *args;

const char *about(void)
{
    return "Check if ploidy of samples is consistent for all sites\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Check if ploidy of samples is consistent for all sites.\n"
        "Usage: bcftools +check-ploidy [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Example:\n"
        "   bcftools +check-ploidy file.bcf\n"
        "\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->hdr = in;
    args->ndat = bcf_hdr_nsamples(args->hdr);
    args->dat  = (dat_t*) calloc(args->ndat,sizeof(dat_t));
    int i;
    for (i=0; i<args->ndat; i++) args->dat[i].sample = args->hdr->samples[i];
    args->rid = -1;
    args->gt_id = bcf_hdr_id2int(args->hdr,BCF_DT_ID,"GT");
    if ( args->gt_id<0 ) error("Error: GT field is not present\n");
    printf("# [1]Sample\t[2]Chromosome\t[3]Region Start\t[4]Region End\t[5]Ploidy\n");
    return 1;
}

bcf1_t *process(bcf1_t *rec)
{
    int i;

    bcf_unpack(rec, BCF_UN_FMT);
    bcf_fmt_t *fmt_gt = NULL;
    for (i=0; i<rec->n_fmt; i++)
        if ( rec->d.fmt[i].id==args->gt_id ) { fmt_gt = &rec->d.fmt[i]; break; }
    if ( !fmt_gt ) return NULL;    // no GT tag

    if ( args->ndat != rec->n_sample ) 
        error("Incorrect number of samples at %s:%d .. found %d, expected %d\n",bcf_seqname(args->hdr,rec),rec->pos+1,rec->n_sample,args->ndat);

    if ( args->rid!=rec->rid && args->rid!=-1 )
    {
        for (i=0; i<args->ndat; i++)
        {
            dat_t *dat = &args->dat[i];
            if ( dat->ploidy!=0 ) printf("%s\t%s\t%d\t%d\t%d\n", dat->sample,bcf_seqname(args->hdr,rec),dat->beg+1,dat->end+1,dat->ploidy); 
            dat->ploidy = 0;
        }
    }
    args->rid = rec->rid;

    #define BRANCH_INT(type_t,vector_end) \
    { \
        for (i=0; i<rec->n_sample; i++) \
        { \
            type_t *p = (type_t*) (fmt_gt->p + i*fmt_gt->size); \
            int nal, missing = 0; \
            for (nal=0; nal<fmt_gt->n; nal++) \
            { \
                if ( p[nal]==vector_end ) break; /* smaller ploidy */ \
                if ( bcf_gt_is_missing(p[nal]) ) { missing=1; break; } /* missing allele */ \
            } \
            if ( !nal || missing ) continue; /* missing genotype */ \
            dat_t *dat = &args->dat[i]; \
            if ( dat->ploidy==nal ) \
            { \
                dat->end = rec->pos; \
                continue; \
            } \
            if ( dat->ploidy!=0 ) \
                printf("%s\t%s\t%d\t%d\t%d\n", dat->sample,bcf_seqname(args->hdr,rec),dat->beg+1,dat->end+1,dat->ploidy); \
            dat->ploidy = nal; \
            dat->beg = rec->pos; \
            dat->end = rec->pos; \
        } \
    }
    switch (fmt_gt->type) {
        case BCF_BT_INT8:  BRANCH_INT(int8_t,  bcf_int8_vector_end); break;
        case BCF_BT_INT16: BRANCH_INT(int16_t, bcf_int16_vector_end); break;
        case BCF_BT_INT32: BRANCH_INT(int32_t, bcf_int32_vector_end); break;
        default: error("The GT type is not recognised: %d at %s:%d\n",fmt_gt->type, bcf_seqname(args->hdr,rec),rec->pos+1); break;
    }
    #undef BRANCH_INT

    return NULL;
}

void destroy(void)
{
    int i;
    for (i=0; i<args->ndat; i++)
    {
        dat_t *dat = &args->dat[i];
        if ( dat->ploidy!=0 ) printf("%s\t%s\t%d\t%d\t%d\n", dat->sample,bcf_hdr_id2name(args->hdr,args->rid),dat->beg+1,dat->end+1,dat->ploidy); 
        dat->ploidy = 0;
    }
    free(args->dat);
    free(args);
}

