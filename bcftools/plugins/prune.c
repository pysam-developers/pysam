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
    Prune sites by missingness, LD

    See calc_ld() in vcfbuf.c for the actual LD calculation

*/

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <getopt.h>
#include <stdarg.h>
#include <unistd.h>
#include <stdint.h>
#include <errno.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include "bcftools.h"
#include "vcfbuf.h"
#include "filter.h"

#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef struct
{
    filter_t *filter;
    char *filter_str, *af_tag;
    int filter_logic;   // one of FLT_INCLUDE/FLT_EXCLUDE (-i or -e)
    vcfbuf_t *vcfbuf;
    int argc, region_is_file, target_is_file, output_type, filter_r2_id, rand_missing, nsites, ld_win;
    char **argv, *region, *target, *fname, *output_fname, *info_pos, *info_r2, *filter_r2;
    htsFile *out_fh;
    bcf_hdr_t *hdr;
    bcf_srs_t *sr;
    double max_ld;
}
args_t;

const char *about(void)
{
    return "Prune sites by missingness, linkage disequilibrium\n";
}

static const char *usage_text(void)
{
    return 
        "\n"
        "About: Prune sites by missingness or linkage disequilibrium.\n"
        "\n"
        "Usage: bcftools +prune [Options]\n"
        "Plugin options:\n"
        "       --AF-tag STR                use this tag with -n to determine allele frequency\n"
        "   -a, --annotate-info STR         add INFO/STR_POS and STR_R2 annotation: an upstream site with the biggest r2 value\n"
        "   -e, --exclude EXPR              exclude sites for which the expression is true\n"
        "   -f, --set-filter STR            annotate FILTER column with STR instead of discarding the site\n"
        "   -i, --include EXPR              include only sites for which the expression is true\n"
        "   -l, --max-LD R2                 remove sites with r2 bigger than R2 within within the -w window\n"
        "   -n, --nsites-per-win N          keep at most N sites in the -w window, removing sites with small AF first\n"
        "   -o, --output FILE               write output to the FILE [standard output]\n"
        "   -O, --output-type b|u|z|v       b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n"
        "       --randomize-missing         replace missing data with randomly assigned genotype based on site's allele frequency\n"
        "   -r, --regions REGION            restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE         restrict to regions listed in a file\n"
        "   -t, --targets REGION            similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE         similar to -R but streams rather than index-jumps\n"
        "   -w, --window INT[bp|kb]         the window size of INT sites/bp/kb for the -n/-l options [100kb]\n"
        "Examples:\n"
        "   # Discard records with r2 bigger than 0.6 in a window of 1000 sites\n"
        "   bcftools +prune -l 0.6 -w 1000 input.bcf -Ob -o output.bcf\n"
        "\n"
        "   # Set FILTER (but do not discard) records with r2 bigger than 0.4 in the default window of 100kb\n"
        "   bcftools +prune -l 0.4 -f MAX_R2 input.bcf -Ob -o output.bcf\n"
        "\n"
        "   # Annotate INFO field of all records with maximum r2 in a window of 1000 sites\n"
        "   bcftools +prune -l 0.6 -w 1000 -f MAX_R2 input.bcf -Ob -o output.bcf\n"
        "\n"
        "   # Discard records with r2 bigger than 0.6, first removing records with more than 2% of genotypes missing\n"
        "   bcftools +prune -l 0.6 -e'F_MISSING>=0.02' input.bcf -Ob -o output.bcf\n"
        "\n";
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
    args->hdr = bcf_sr_get_header(args->sr,0);

    args->out_fh = hts_open(args->output_fname,hts_bcf_wmode(args->output_type));
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
    if ( args->filter_r2 )
    {
        bcf_hdr_printf(args->hdr,"##FILTER=<ID=%s,Description=\"A site with r2>%e upstream within %d%s\">",args->filter_r2,args->max_ld,
                args->ld_win < 0 ? -args->ld_win/1000 : args->ld_win,
                args->ld_win < 0 ? "kb" : " sites");
    }
    if ( args->info_r2 )
    {
        bcf_hdr_printf(args->hdr,"##INFO=<ID=%s,Number=1,Type=Integer,Description=\"A site with r2>%e upstream\">",args->info_pos,args->max_ld);
        bcf_hdr_printf(args->hdr,"##INFO=<ID=%s,Number=1,Type=Float,Description=\"A site with r2>%e upstream\">",args->info_r2,args->max_ld);
    }
    bcf_hdr_write(args->out_fh, args->hdr);
    if ( args->filter_r2 )
        args->filter_r2_id = bcf_hdr_id2int(args->hdr, BCF_DT_ID, args->filter_r2);

    args->vcfbuf = vcfbuf_init(args->hdr, args->ld_win);
    vcfbuf_set_opt(args->vcfbuf,double,VCFBUF_LD_MAX,args->max_ld);
    if ( args->nsites ) vcfbuf_set_opt(args->vcfbuf,int,VCFBUF_NSITES,args->nsites);
    if ( args->af_tag ) vcfbuf_set_opt(args->vcfbuf,char*,VCFBUF_AF_TAG,args->af_tag);
    if ( args->rand_missing ) vcfbuf_set_opt(args->vcfbuf,int,VCFBUF_RAND_MISSING,1);
    vcfbuf_set_opt(args->vcfbuf,int,VCFBUF_SKIP_FILTER,args->filter_r2 ? 1 : 0);

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);
}
static void destroy_data(args_t *args)
{
    if ( args->filter )
        filter_destroy(args->filter);
    hts_close(args->out_fh);
    vcfbuf_destroy(args->vcfbuf);
    bcf_sr_destroy(args->sr);
    free(args->info_pos);
    free(args->info_r2);
    free(args);
}
static void flush(args_t *args, int flush_all)
{
    bcf1_t *rec;
    while ( (rec = vcfbuf_flush(args->vcfbuf, flush_all)) )
        bcf_write1(args->out_fh, args->hdr, rec);
}
static void process(args_t *args)
{
    bcf1_t *rec = bcf_sr_get_line(args->sr,0);
    if ( args->filter )
    {
        int ret  = filter_test(args->filter, rec, NULL);
        if ( args->filter_logic==FLT_INCLUDE ) { if ( !ret ) return; }
        else if ( ret ) return;
    }
    bcf_sr_t *sr = bcf_sr_get_reader(args->sr, 0);
    if ( args->max_ld )
    {
        double ld_val;
        bcf1_t *ld_rec = vcfbuf_max_ld(args->vcfbuf, rec, &ld_val);
        if ( ld_rec && ld_val > args->max_ld )
        {
            if ( !args->filter_r2 ) return;
            bcf_add_filter(args->hdr, rec, args->filter_r2_id);
        }
        if ( ld_rec && args->info_r2 )
        {
            float tmp = ld_val;
            bcf_update_info_float(args->hdr, rec, args->info_r2, &tmp, 1);
            bcf_update_info_int32(args->hdr, rec, args->info_pos, &ld_rec->pos, 1);
        }
    }
    sr->buffer[0] = vcfbuf_push(args->vcfbuf, rec, 1);
    flush(args,0);
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->output_type  = FT_VCF;
    args->output_fname = "-";
    args->ld_win = -100e3;
    static struct option loptions[] =
    {
        {"randomize-missing",no_argument,NULL,1},
        {"AF-tag",required_argument,NULL,2},
        {"exclude",required_argument,NULL,'e'},
        {"include",required_argument,NULL,'i'},
        {"annotate-info",required_argument,NULL,'a'},
        {"set-filter",required_argument,NULL,'f'},
        {"max-LD",required_argument,NULL,'l'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"nsites-per-win",required_argument,NULL,'n'},
        {"window",required_argument,NULL,'w'},
        {NULL,0,NULL,0}
    };
    int c;
    char *tmp;
    while ((c = getopt_long(argc, argv, "vr:R:t:T:l:o:O:a:f:i:e:n:w:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case  1 : args->rand_missing = 1; break;
            case  2 : args->af_tag = optarg; break;
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'a': 
                {
                    int l = strlen(optarg);
                    args->info_pos = (char*)malloc(l+5);
                    args->info_r2  = (char*)malloc(l+5);
                    sprintf(args->info_pos,"%s_POS", optarg);
                    sprintf(args->info_r2,"%s_R2", optarg);
                }
                break; 
            case 'f': args->filter_r2 = optarg; break;
            case 'n': 
                args->nsites = strtod(optarg,&tmp);
                if ( tmp==optarg || *tmp ) error("Could not parse: --nsites-per-win %s\n", optarg);
                break;
            case 'l': 
                args->max_ld = strtod(optarg,&tmp);
                if ( tmp==optarg || *tmp ) error("Could not parse: --max-LD %s\n", optarg);
                break;
            case 'w': 
                args->ld_win = strtod(optarg,&tmp);
                if ( !*tmp ) break;
                if ( tmp==optarg ) error("Could not parse: --window %s\n", optarg);
                else if ( !strcasecmp("bp",tmp) ) args->ld_win *= -1;
                else if ( !strcasecmp("kb",tmp) ) args->ld_win *= -1000;
                else error("Could not parse: --window %s\n", optarg);
                break;
            case 'T': args->target_is_file = 1; 
            case 't': args->target = optarg; break; 
            case 'R': args->region_is_file = 1; 
            case 'r': args->region = optarg; break; 
            case 'o': args->output_fname = optarg; break;
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
    if ( args->filter_logic == (FLT_EXCLUDE|FLT_INCLUDE) ) error("Only one of -i or -e can be given.\n");
    if ( !args->max_ld && !args->nsites ) error("%sError: Expected --max-LD, --nsites-per-win or both\n\n", usage_text());

    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else { error(usage_text()); }
    }
    else if ( optind+1!=argc ) error(usage_text());
    else args->fname = argv[optind];

    init_data(args);
    
    while ( bcf_sr_next_line(args->sr) ) process(args);
    flush(args,1);

    destroy_data(args);
    return 0;
}


