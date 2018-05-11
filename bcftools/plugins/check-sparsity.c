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
    int argc;
    char **argv, *fname, *region, **regs;
    int region_is_file, nregs, regs_free;
    int *smpl, nsmpl, *nsites, min_sites, gt_id;
    kstring_t tmps;
    bcf1_t *rec;
    tbx_t *tbx;
    hts_idx_t *idx;
    hts_itr_t *itr;
    htsFile *fp;
    bcf_hdr_t *hdr;
}
args_t;

const char *about(void)
{
    return "Print samples without genotypes in a region or chromosome\n";
}

static const char *usage_text(void)
{
    return 
        "\n"
        "About: Print samples without genotypess in a region (-r/-R) or chromosome (the default)\n"
        "\n"
        "Usage: bcftools +check-sparsity <file.vcf.gz> [Plugin Options]\n"
        "Plugin options:\n"
        "   -n, --n-markers <int>           minimum number of required markers [1]\n"
        "   -r, --regions <chr:beg-end>     restrict to comma-separated list of regions\n"
        "   -R, --regions-file <file>       restrict to regions listed in a file\n"
        "\n";
}

static void init_data(args_t *args)
{
    args->fp = hts_open(args->fname,"r");
    if ( !args->fp ) error("Could not read %s\n", args->fname);
    args->hdr = bcf_hdr_read(args->fp);
    if ( !args->hdr ) error("Could not read the header: %s\n", args->fname);

    args->rec = bcf_init1();
    args->gt_id = bcf_hdr_id2int(args->hdr,BCF_DT_ID,"GT");
    if ( args->gt_id<0 ) error("Error: GT field is not present\n");

    int i;
    args->nsmpl  = bcf_hdr_nsamples(args->hdr);
    args->nsites = (int*) calloc(args->nsmpl, sizeof(int));
    args->smpl   = (int*) malloc(sizeof(int)*args->nsmpl);
    for (i=0; i<args->nsmpl; i++) args->smpl[i] = i;

    if ( strcmp("-",args->fname) )  // not reading from stdin
    {
        if ( hts_get_format(args->fp)->format==vcf )
        {
            args->tbx = tbx_index_load(args->fname);
            if ( !args->tbx && args->region ) error("Could not load the VCF index, please drop the -r/-R option\n");
        }
        else if ( hts_get_format(args->fp)->format==bcf )
        {
            args->idx = bcf_index_load(args->fname);
            if ( !args->idx && args->region ) error("Could not load the BCF index, please drop the -r/-R option\n");
        }
    }
    else if ( args->region ) error("Cannot use index with this file, please drop the -r/-R option\n");

    if ( args->tbx || args->idx )
    {
        if ( args->region )
        {
            args->regs = hts_readlist(args->region, args->region_is_file, &args->nregs);
            if ( !args->regs ) error("Could not parse regions: %s\n", args->region);
            args->regs_free = 1;
        }
        else
            args->regs = (char**) (args->tbx ? tbx_seqnames(args->tbx, &args->nregs) : bcf_index_seqnames(args->idx, args->hdr, &args->nregs));
    }
}
static void destroy_data(args_t *args)
{
    int i;
    if ( args->regs_free )
        for (i=0; i<args->nregs; i++) free(args->regs[i]);
    free(args->regs);
    bcf_hdr_destroy(args->hdr);
    bcf_destroy(args->rec);
    free(args->tmps.s);
    free(args->smpl);
    free(args->nsites);
    if ( args->itr ) hts_itr_destroy(args->itr);
    if ( args->tbx ) tbx_destroy(args->tbx);
    if ( args->idx ) hts_idx_destroy(args->idx);
    hts_close(args->fp);
}

static void report(args_t *args, const char *reg)
{
    int i;
    for (i=0; i<args->nsmpl; i++)
        printf("%s\t%s\n", reg, args->hdr->samples[args->smpl[i]]);
    args->nsmpl = bcf_hdr_nsamples(args->hdr);
    for (i=0; i<args->nsmpl; i++) args->smpl[i] = i;
    memset(args->nsites, 0, sizeof(int)*args->nsmpl);
}
static void test_region(args_t *args, char *reg)
{
    if ( args->tbx )
    {
        args->itr = tbx_itr_querys(args->tbx,reg);
        if ( !args->itr ) return;
    }
    else if ( args->idx )
    {
        args->itr = bcf_itr_querys(args->idx,args->hdr,reg);
        if ( !args->itr ) return;
    }

    int ret,i, rid = -1, nread = 0;
    while (1)
    {
        if ( args->tbx )
        {
            if ( (ret=tbx_itr_next(args->fp, args->tbx, args->itr, &args->tmps)) < 0 ) break;  // no more lines
            ret = vcf_parse1(&args->tmps, args->hdr, args->rec);
            if ( ret<0 ) error("Could not parse the line: %s\n", args->tmps.s);
        }
        else if ( args->idx )
        {
            ret = bcf_itr_next(args->fp, args->itr, args->rec);
            if ( ret < -1 ) error("Could not parse a line from %s\n", reg);
            if ( ret < 0 ) break; // no more lines or an error
        }
        else
        {
            if ( args->fp->format.format==vcf )
            {
                if ( (ret=hts_getline(args->fp, KS_SEP_LINE, &args->tmps)) < 0 ) break;   // no more lines
                ret = vcf_parse1(&args->tmps, args->hdr, args->rec);
                if ( ret<0 ) error("Could not parse the line: %s\n", args->tmps.s);
            }
            else if ( args->fp->format.format==bcf )
            {
                ret = bcf_read1(args->fp, args->hdr, args->rec);
                if ( ret < -1 ) error("Could not parse %s\n", args->fname);
                if ( ret < 0 ) break; // no more lines or an error
            }
            if ( rid!=-1 && rid!=args->rec->rid )
            {
                report(args, bcf_hdr_id2name(args->hdr,rid));
                nread = 0;
            }
            rid = args->rec->rid;
        }

        bcf_unpack(args->rec, BCF_UN_FMT);
        bcf_fmt_t *fmt_gt = NULL;
        for (i=0; i<args->rec->n_fmt; i++)
            if ( args->rec->d.fmt[i].id==args->gt_id ) { fmt_gt = &args->rec->d.fmt[i]; break; }
        if ( !fmt_gt ) continue;        // no GT tag
        if ( fmt_gt->n==0 ) continue;   // empty?!
        if ( fmt_gt->type!=BCF_BT_INT8 ) error("TODO: the GT fmt_type is not int8!\n");

        // update the array of missing samples
        for (i=0; i<args->nsmpl; i++)
        {
            int8_t *ptr = (int8_t*) (fmt_gt->p + args->smpl[i]*fmt_gt->size);
            int ial = 0;
            for (ial=0; ial<fmt_gt->n; ial++)
                if ( ptr[ial]==bcf_gt_missing || ptr[ial]==bcf_int8_vector_end ) break;
            if ( ial==0 ) continue;     // missing
            if ( ++args->nsites[i] < args->min_sites ) continue;
            if ( i+1<args->nsmpl )
            {
                memmove(args->smpl+i, args->smpl+i+1, sizeof(int)*(args->nsmpl-i-1));
                memmove(args->nsites+i, args->nsites+i+1, sizeof(int)*(args->nsmpl-i-1));
            }
            args->nsmpl--;
            i--;
        }
        nread = 1;
        if ( !args->nsmpl ) break;
    }
    if ( nread ) report(args, rid==-1 ? reg : bcf_hdr_id2name(args->hdr,rid));

    tbx_itr_destroy(args->itr);
    args->itr = NULL;
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->min_sites = 1;
    static struct option loptions[] =
    {
        {"n-markers",required_argument,NULL,'n'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {NULL,0,NULL,0}
    };
    int c,i;
    char *tmp;
    while ((c = getopt_long(argc, argv, "vr:R:n:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'n': 
                args->min_sites = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse: -n %s\n", optarg);
                break;
            case 'R': args->region_is_file = 1; 
            case 'r': args->region = optarg; break; 
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }

    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else error(usage_text());
    }
    else args->fname = argv[optind];
    init_data(args);

    for (i=0; i<args->nregs; i++) test_region(args, args->regs[i]);
    if ( !args->nregs ) test_region(args, NULL);

    destroy_data(args);
    free(args);
    return 0;
}

