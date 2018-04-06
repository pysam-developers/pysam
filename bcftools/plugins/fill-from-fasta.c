/*  plugin/fill-from-fasta.c -- fill-from-fasta plugin.

    Copyright (C) 2016 Genome Research Ltd.

    Author: Shane McCarthy <sm15@sanger.ac.uk>

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
#include <strings.h>
#include <getopt.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include <htslib/kseq.h>
#include "filter.h"
#include "bcftools.h"

const char *about(void)
{
    return "Fill INFO or REF field based on values in a fasta file\n";
}

const char *usage(void)
{
    return 
"\n"
"About:   Fill INFO or REF field based on values in a fasta file.\n"
"         The fasta file must be indexed with samtools faidx.\n"
"Usage:   bcftools +fill-from-fasta [General Options] -- [Plugin Options]\n"
"\n"
"General options:\n"
"   run \"bcftools plugin\" for a list of common options\n"
"\n"
"Plugin options:\n"
"   -c, --column <str>          REF or INFO tag, e.g. AA for ancestral allele\n"
"   -f, --fasta <file>          fasta file\n"
"   -h, --header-lines <file>   optional file containing header lines to append\n"
"   -i, --include <expr>        annotate only records passing filter expression\n"
"   -e, --exclude <expr>        annotate only records failing filter expression\n"

"\n"
"Examples:\n"
"   # fill ancestral allele as INFO/AA for SNP records\n"
"   echo '##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral allele\">' > aa.hdr\n"
"   bcftools +fill-from-fasta in.vcf -- -c AA -f aa.fasta -h aa.hdr -i 'TYPE=\"snp\"'\n"
"\n"
"   # fix the REF allele in VCFs where REF=N or other\n"
"   bcftools +fill-from-fasta in.vcf -- -c REF -f reference.fasta\n"
"\n"
"   # select sites marked as P (PASS) in the 1000G Phase3 mask\n"
"   echo '##INFO=<ID=P3_MASK,Number=1,Type=String,Description=\"1000G Phase 3 mask\">' > mask.hdr\n"
"   bcftools +fill-from-fasta in.vcf -Ou -- -c P3_MASK -f 1000G_mask.fasta -h mask.hdr | bcftools view -i 'P3_MASK=\"P\"'\n"
"\n";
}

bcf_hdr_t *in_hdr = NULL, *out_hdr = NULL;
faidx_t *faidx;
int anno = 0;
char *column = NULL;

#define ANNO_REF 1
#define ANNO_STRING 2
#define ANNO_INT 3

filter_t *filter;
char *filter_str;
int filter_logic;   // one of FLT_INCLUDE/FLT_EXCLUDE (-i or -e)

#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    int c;
    char *ref_fname = NULL, *header_fname = NULL;
    static struct option loptions[] =
    {
        {"exclude",required_argument,NULL,'e'},
        {"include",required_argument,NULL,'i'},
        {"column",required_argument,NULL,'c'},
        {"fasta",required_argument,NULL,'f'},
        {"header-lines",required_argument,NULL,'h'},
        {NULL,0,NULL,0}
    };
    while ((c = getopt_long(argc, argv, "c:f:?h:i:e:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'e': filter_str = optarg; filter_logic |= FLT_EXCLUDE; break;
            case 'i': filter_str = optarg; filter_logic |= FLT_INCLUDE; break;
            case 'c': column = optarg; break;
            case 'f': ref_fname = optarg; break;
            case 'h': header_fname = optarg; break;
            case '?':
            default: fprintf(stderr,"%s", usage()); exit(1); break;
        }
    }
    in_hdr  = in;
    out_hdr = out;
    if ( filter_logic == (FLT_EXCLUDE|FLT_INCLUDE) ) { fprintf(stderr,"Only one of -i or -e can be given.\n"); return -1; }

    if ( !column )
    {
        fprintf(stderr,"--column option is required.\n");
        return -1;
    }
    if (header_fname)
    {
        htsFile *file = hts_open(header_fname, "rb");
        if ( !file ) { fprintf(stderr,"Error reading %s\n", header_fname); return -1; }
        kstring_t str = {0,0,0};
        while ( hts_getline(file, KS_SEP_LINE, &str) > 0 )
        {
            if ( bcf_hdr_append(out_hdr, str.s) ) { fprintf(stderr,"Could not parse %s: %s\n", header_fname, str.s); return -1; }
        }
        hts_close(file);
        free(str.s);
        bcf_hdr_sync(out_hdr);
    }
    if (!strcasecmp("REF", column)) anno = ANNO_REF;
    else {
        if ( !strncasecmp(column,"INFO/",5) ) column += 5;
        int hdr_id = bcf_hdr_id2int(out_hdr, BCF_DT_ID, column);
        if (hdr_id<0) { fprintf(stderr,"No header ID found for %s. Header lines can be added with the --header-lines option\n", column); return -1; }
        switch ( bcf_hdr_id2type(out_hdr,BCF_HL_INFO,hdr_id) )
        {
            case BCF_HT_INT:
                anno=ANNO_INT;
                break;
            case BCF_HT_STR:
                anno=ANNO_STRING;
                break;
            default:
                fprintf(stderr,"The type of %s not recognised (%d)\n", column, bcf_hdr_id2type(out_hdr,BCF_HL_INFO,hdr_id));
                return -1;
        }
    }
    if ( !ref_fname )
    {
        fprintf(stderr,"No fasta given.\n");
        return -1;
    }
    faidx = fai_load(ref_fname);
    if ( filter_str )
        filter = filter_init(in, filter_str);
    return 0;
}

bcf1_t *process(bcf1_t *rec)
{
    // filter determines if we will annotate the record
    // return record unchanged if filter applied
    if ( filter )
    {
        int ret = filter_test(filter, rec, NULL);
        if ( filter_logic==FLT_INCLUDE ) { if ( !ret ) return rec; }
        else if ( ret ) return rec;
    }

    int i;
    char *ref = rec->d.allele[0];
    int ref_len = strlen(ref);
    int fa_len;
    // could be sped up here by fetching the whole chromosome? could assume
    // sorted, but revert to this when non-sorted records found?
    char *fa = faidx_fetch_seq(faidx, bcf_seqname(in_hdr,rec), rec->pos, rec->pos+ref_len-1, &fa_len);
    if ( !fa ) error("faidx_fetch_seq failed at %s:%d\n", bcf_hdr_id2name(in_hdr,rec->rid), rec->pos+1);
    for (i=0; i<fa_len; i++)
        if ( (int)fa[i]>96 ) fa[i] -= 32;

    assert(ref_len == fa_len);
    if (anno==ANNO_REF)
        strncpy(rec->d.allele[0], fa, fa_len);
    else if (anno==ANNO_STRING)
        bcf_update_info_string(out_hdr, rec, column, fa);
    else if (anno==ANNO_INT && ref_len==1)
    {
        int val = atoi(&fa[0]);
        bcf_update_info_int32(out_hdr, rec, column, &val, 1);
    }
    free(fa);
    return rec;
}

void destroy(void)
{
    fai_destroy(faidx);
    if (filter) filter_destroy(filter);
}
