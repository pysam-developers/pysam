#include "bcftools.pysam.h"

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
/*
    Split VCF by sample(s)

*/

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <getopt.h>
#include <stdarg.h>
#include <unistd.h>
#include <stdint.h>
#include <errno.h>
#include <ctype.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include "bcftools.h"
#include "filter.h"

#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef struct
{
    htsFile **fh;
    filter_t *filter;
    char *filter_str;
    int filter_logic;   // one of FLT_INCLUDE/FLT_EXCLUDE (-i or -e)
    uint8_t *info_tags, *fmt_tags;
    int ninfo_tags, minfo_tags, nfmt_tags, mfmt_tags, keep_info, keep_fmt;
    int argc, region_is_file, target_is_file, output_type;
    char **argv, *region, *target, *fname, *output_dir, *keep_tags, **bnames, *samples_fname;
    bcf_hdr_t *hdr_in, *hdr_out;
    bcf_srs_t *sr;
}
args_t;

const char *about(void)
{
    return "Split VCF by sample creating single-sample VCFs\n";
}

static const char *usage_text(void)
{
    return 
        "\n"
        "About: Split VCF by sample, creating single-sample VCFs.\n"
        "\n"
        "Usage: bcftools +split [Options]\n"
        "Plugin options:\n"
        "   -e, --exclude EXPR              exclude sites for which the expression is true (applied on the outputs)\n"
        "   -i, --include EXPR              include only sites for which the expression is true (applied on the outputs)\n"
        "   -k, --keep-tags LIST            list of tags to keep. By default all tags are preserved\n"
        "   -o, --output DIR                write output to the directory DIR\n"
        "   -O, --output-type b|u|z|v       b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n"
        "   -r, --regions REGION            restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE         restrict to regions listed in a file\n"
        "   -S, --samples-file FILE         list of samples to keep with second (optional) column for basename of the new file\n"
        "   -t, --targets REGION            similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE         similar to -R but streams rather than index-jumps\n"
        "Examples:\n"
        "   # Split a VCF file\n"
        "   bcftools +split input.bcf -Ob -o dir\n"
        "\n"
        "   # Exclude sites with missing or hom-ref genotypes\n"
        "   bcftools +split input.bcf -Ob -o dir -i'GT=\"alt\"'\n"
        "\n"
        "   # Keep all INFO tags but only GT and PL in FORMAT\n"
        "   bcftools +split input.bcf -Ob -o dir -k INFO,FMT/GT,PL\n"
        "\n"
        "   # Keep all FORMAT tags but drop all INFO tags\n"
        "   bcftools +split input.bcf -Ob -o dir -k FMT\n"
        "\n";
}

void mkdir_p(const char *fmt, ...);

char **set_file_base_names(args_t *args)
{
    int i, nsmpl = bcf_hdr_nsamples(args->hdr_in);
    char **fnames = (char**) calloc(nsmpl,sizeof(char*));
    if ( args->samples_fname )
    {
        kstring_t str = {0,0,0};
        int nsamples = 0;
        char **samples = NULL;
        samples = hts_readlines(args->samples_fname, &nsamples);
        for (i=0; i<nsamples; i++)
        {
            str.l = 0;
            int escaped = 0;
            char *ptr = samples[i];
            while ( *ptr )
            {
                if ( *ptr=='\\' && !escaped ) { escaped = 1; ptr++; continue; }
                if ( isspace(*ptr) && !escaped ) break;
                kputc(*ptr, &str);
                escaped = 0;
                ptr++;
            }
            int idx = bcf_hdr_id2int(args->hdr_in, BCF_DT_SAMPLE, str.s);
            if ( idx<0 )
            {
                fprintf(bcftools_stderr,"Warning: The sample \"%s\" is not present in %s\n", str.s,args->fname);
                continue;
            }
            while ( *ptr && isspace(*ptr) ) ptr++;
            if ( !*ptr )
            {
                fnames[idx] = strdup(str.s);
                continue;
            }
            str.l = 0;
            while ( *ptr )
            {
                if ( *ptr=='\\' && !escaped ) { escaped = 1; ptr++; continue; }
                if ( isspace(*ptr) && !escaped ) break;
                kputc(*ptr, &str);
                escaped = 0;
                ptr++;
            }
            fnames[idx] = strdup(str.s);
        }
        for (i=0; i<nsamples; i++) free(samples[i]);
        free(samples);
        free(str.s);
    }
    else
    {
        for (i=0; i<nsmpl; i++)
            fnames[i] = strdup(args->hdr_in->samples[i]);
    }
    return fnames;
}

static void init_data(args_t *args)
{
    args->sr = bcf_sr_init();
    if ( args->region )
    {
        args->sr->require_index = 1;
        if ( bcf_sr_set_regions(args->sr, args->region, args->region_is_file)<0 ) error("Failed to read the regions: %s\n",args->region);
    }
    if ( args->target && bcf_sr_set_targets(args->sr, args->target, args->target_is_file, 0)<0 ) error("Failed to read the targets: %s\n",args->target);
    if ( !bcf_sr_add_reader(args->sr,args->fname) ) error("Error: %s\n", bcf_sr_strerror(args->sr->errnum));
    args->hdr_in  = bcf_sr_get_header(args->sr,0);
    args->hdr_out = bcf_hdr_dup(args->hdr_in);

    if ( args->filter_str )
        args->filter = filter_init(args->hdr_in, args->filter_str);

    mkdir_p("%s/",args->output_dir);

    int i, nsmpl = bcf_hdr_nsamples(args->hdr_in);
    if ( !nsmpl ) error("No samples to split: %s\n", args->fname);
    args->fh = (htsFile**)calloc(nsmpl,sizeof(*args->fh));
    args->bnames = set_file_base_names(args);
    kstring_t str = {0,0,0};
    for (i=0; i<nsmpl; i++)
    {
        if ( !args->bnames[i] ) continue;
        str.l = 0;
        kputs(args->output_dir, &str);
        if ( str.s[str.l-1] != '/' ) kputc('/', &str);
        int k, l = str.l;
        kputs(args->bnames[i], &str);
        for (k=l; k<str.l; k++) if ( isspace(str.s[k]) ) str.s[k] = '_';
        if ( args->output_type & FT_BCF ) kputs(".bcf", &str);
        else if ( args->output_type & FT_GZ ) kputs(".vcf.gz", &str);
        else kputs(".vcf", &str);
        args->fh[i] = hts_open(str.s, hts_bcf_wmode(args->output_type));
        if ( args->fh[i] == NULL ) error("Can't write to \"%s\": %s\n", str.s, strerror(errno));
        bcf_hdr_nsamples(args->hdr_out) = 1;
        args->hdr_out->samples[0] = args->bnames[i];
        bcf_hdr_write(args->fh[i], args->hdr_out);
    }
    free(str.s);

    // parse tags
    int is_info = 0, is_fmt = 0;
    char *beg = args->keep_tags;
    while ( beg && *beg )
    {
        if ( !strncasecmp("INFO/",beg,5) ) { is_info = 1; is_fmt = 0; beg += 5; }
        else if ( !strcasecmp("INFO",beg) ) { args->keep_info = 1; break; }
        else if ( !strncasecmp("INFO,",beg,5) ) { args->keep_info = 1; beg += 5; continue; }
        else if ( !strncasecmp("FMT/",beg,4) ) { is_info = 0; is_fmt = 1; beg += 4; }
        else if ( !strncasecmp("FORMAT/",beg,7) ) { is_info = 0; is_fmt = 1; beg += 7; }
        else if ( !strcasecmp("FMT",beg) ) { args->keep_fmt = 1; break; }
        else if ( !strcasecmp("FORMAT",beg) ) { args->keep_fmt = 1; break; }
        else if ( !strncasecmp("FMT,",beg,4) ) { args->keep_fmt = 1; beg += 4; continue; }
        else if ( !strncasecmp("FORMAT,",beg,7) ) { args->keep_fmt = 1; beg += 7; continue; }
        char *end = beg;
        while ( *end && *end!=',' ) end++;
        char tmp = *end; *end = 0;
        int id = bcf_hdr_id2int(args->hdr_in, BCF_DT_ID, beg);
        beg = tmp ? end + 1 : end;
        if ( is_info && bcf_hdr_idinfo_exists(args->hdr_in,BCF_HL_INFO,id) )
        {
            if ( id >= args->ninfo_tags ) args->ninfo_tags = id + 1;
            hts_expand0(uint8_t, args->ninfo_tags, args->minfo_tags, args->info_tags);
            args->info_tags[id] = 1;
        }
        if ( is_fmt && bcf_hdr_idinfo_exists(args->hdr_in,BCF_HL_FMT,id) )
        {
            if ( id >= args->nfmt_tags ) args->nfmt_tags = id + 1;
            hts_expand0(uint8_t, args->nfmt_tags, args->mfmt_tags, args->fmt_tags);
            args->fmt_tags[id] = 1;
        }
    }
    if ( !args->keep_info && !args->keep_fmt && !args->ninfo_tags && !args->nfmt_tags )
    {
        args->keep_info = args->keep_fmt = 1;
    }
}
static void destroy_data(args_t *args)
{
    free(args->info_tags);
    free(args->fmt_tags);
    if ( args->filter )
        filter_destroy(args->filter);
    int i, nsmpl = bcf_hdr_nsamples(args->hdr_in);
    for (i=0; i<nsmpl; i++)
    {
        if ( args->fh[i] && hts_close(args->fh[i])!=0 ) error("Error: close failed!\n");
        free(args->bnames[i]);
    }
    free(args->bnames);
    free(args->fh);
    bcf_sr_destroy(args->sr);
    bcf_hdr_destroy(args->hdr_out);
    free(args);
}

static bcf1_t *rec_set_info(args_t *args, bcf1_t *rec)
{
    bcf1_t *out = bcf_init1();
    out->rid  = rec->rid;
    out->pos  = rec->pos;
    out->rlen = rec->rlen;
    out->qual = rec->qual;
    out->n_allele = rec->n_allele;
    out->n_sample = 1;
    if ( args->keep_info )
    {
        out->n_info = rec->n_info;
        out->shared.m = out->shared.l = rec->shared.l;
        out->shared.s = (char*) malloc(out->shared.l);
        memcpy(out->shared.s,rec->shared.s,out->shared.l);
        return out;
    }

    // build the BCF record
    kstring_t tmp = {0,0,0};
    char *ptr = rec->shared.s;
    kputsn_(ptr, rec->unpack_size[0], &tmp); ptr += rec->unpack_size[0]; // ID
    kputsn_(ptr, rec->unpack_size[1], &tmp); ptr += rec->unpack_size[1]; // REF+ALT
    kputsn_(ptr, rec->unpack_size[2], &tmp);                             // FILTER
    if ( args->ninfo_tags )
    {
        int i;
        for (i=0; i<rec->n_info; i++)
        {
            bcf_info_t *info = &rec->d.info[i];
            int id = info->key;
            if ( !args->info_tags[id] ) continue;
            kputsn_(info->vptr - info->vptr_off, info->vptr_len + info->vptr_off, &tmp);
            out->n_info++;
        }
    }
    out->shared.m = tmp.m;
    out->shared.s = tmp.s;
    out->shared.l = tmp.l;
    out->unpacked = 0;
    return out;
}

static bcf1_t *rec_set_format(args_t *args, bcf1_t *src, int ismpl, bcf1_t *dst)
{
    dst->n_fmt = 0;
    kstring_t tmp = dst->indiv; tmp.l = 0;
    int i;
    for (i=0; i<src->n_fmt; i++)
    {
        bcf_fmt_t *fmt = &src->d.fmt[i];
        int id = fmt->id;
        if ( !args->keep_fmt && !args->fmt_tags[id] ) continue;

        bcf_enc_int1(&tmp, id);
        bcf_enc_size(&tmp, fmt->n, fmt->type);
        kputsn_(fmt->p + ismpl*fmt->size, fmt->size, &tmp);

        dst->n_fmt++;
    }
    dst->indiv = tmp;
    return dst;
}

static void process(args_t *args)
{
    bcf1_t *rec = bcf_sr_get_line(args->sr,0);
    bcf_unpack(rec, BCF_UN_ALL);

    int i, site_pass = 1;
    const uint8_t *smpl_pass = NULL;
    if ( args->filter )
    {
        site_pass = filter_test(args->filter, rec, &smpl_pass);
        if ( args->filter_logic & FLT_EXCLUDE ) site_pass = site_pass ? 0 : 1;
    }
    bcf1_t *out = NULL; 
    for (i=0; i<rec->n_sample; i++)
    {
        if ( !args->fh[i] ) continue;
        if ( !smpl_pass && !site_pass ) continue;
        if ( smpl_pass )
        {
            int pass = args->filter_logic & FLT_EXCLUDE ? ( smpl_pass[i] ? 0 : 1) : smpl_pass[i];
            if ( !pass ) continue;
        }
        if ( !out ) out = rec_set_info(args, rec);
        rec_set_format(args, rec, i, out);
        bcf_write(args->fh[i], args->hdr_out, out);
    }
    if ( out ) bcf_destroy(out);
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->output_type  = FT_VCF;
    static struct option loptions[] =
    {
        {"keep-tags",required_argument,NULL,'k'},
        {"exclude",required_argument,NULL,'e'},
        {"include",required_argument,NULL,'i'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"samples-file",required_argument,NULL,'S'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {NULL,0,NULL,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "vr:R:t:T:o:O:i:e:k:S:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'k': args->keep_tags = optarg; break;
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'T': args->target = optarg; args->target_is_file = 1; break;
            case 't': args->target = optarg; break; 
            case 'R': args->region = optarg; args->region_is_file = 1;  break;
            case 'S': args->samples_fname = optarg; break;
            case 'r': args->region = optarg; break; 
            case 'o': args->output_dir = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 'b': args->output_type = FT_BCF_GZ; break;
                          case 'u': args->output_type = FT_BCF; break;
                          case 'z': args->output_type = FT_VCF_GZ; break;
                          case 'v': args->output_type = FT_VCF; break;
                          default: error("The output type \"%s\" not recognised\n", optarg);
                      }
                      break;
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }
    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else { error("%s", usage_text()); }
    }
    else if ( optind+1!=argc ) error("%s", usage_text());
    else args->fname = argv[optind];

    if ( !args->output_dir ) error("Missing the -o option\n");
    if ( args->filter_logic == (FLT_EXCLUDE|FLT_INCLUDE) ) error("Only one of -i or -e can be given.\n");

    init_data(args);
    
    while ( bcf_sr_next_line(args->sr) ) process(args);

    destroy_data(args);
    return 0;
}


