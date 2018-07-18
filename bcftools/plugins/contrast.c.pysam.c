#include "bcftools.pysam.h"

/* The MIT License

   Copyright (c) 2018 Genome Research Ltd.

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
#include <errno.h>
#include <unistd.h>     // for isatty
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>
#include "bcftools.h"
#include "filter.h"


// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef struct
{
    int argc, filter_logic, regions_is_file, targets_is_file, output_type;
    char **argv, *output_fname, *fname, *regions, *targets, *filter_str;
    char *bg_samples_str, *novel_samples_str;
    int *bg_smpl, *novel_smpl, nbg_smpl, nnovel_smpl;
    filter_t *filter;
    bcf_srs_t *sr;
    bcf_hdr_t *hdr, *hdr_out;
    htsFile *out_fh;
    int32_t *gts;
    int mgts;
    uint32_t *bg_gts;
    int nbg_gts, mbg_gts, ntotal, nskipped, ntested, nnovel_al, nnovel_gt;
    kstring_t novel_als_smpl, novel_gts_smpl;
}
args_t;

args_t args;

const char *about(void)
{
    return "Find novel alleles and genotypes in two groups of samples.\n";
}

static const char *usage_text(void)
{
    return 
        "\n"
        "About: Finds novel alleles and genotypes in two groups of samples. Adds\n"
        "       an annotation which lists samples with a novel allele (INFO/NOVELAL)\n"
        "       or a novel genotype (INFO/NOVELGT)\n"
        "Usage: bcftools +contrast [Plugin Options]\n"
        "Plugin options:\n"
        "   -0, --bg-samples <list>     list of background samples\n"
        "   -1, --novel-samples <list>  list of samples where novel allele or genotype are expected\n"
        "   -e, --exclude EXPR          exclude sites and samples for which the expression is true\n"
        "   -i, --include EXPR          include sites and samples for which the expression is true\n"
        "   -o, --output FILE           output file name [bcftools_stdout]\n"
        "   -O, --output-type <b|u|z|v> b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n"
        "   -r, --regions REG           restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE     restrict to regions listed in a file\n"
        "   -t, --targets REG           similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE     similar to -R but streams rather than index-jumps\n"
        "\n"
        "Example:\n"
        "   # Test if any of the samples a,b is different from the samples c,d,e\n"
        "   bcftools +contrast -0 c,d,e -1 a,b file.bcf\n"
        "\n";
}

static void init_data(args_t *args)
{
    args->sr = bcf_sr_init();
    if ( args->regions )
    {
        args->sr->require_index = 1;
        if ( bcf_sr_set_regions(args->sr, args->regions, args->regions_is_file)<0 ) error("Failed to read the regions: %s\n",args->regions);
    }
    if ( args->targets && bcf_sr_set_targets(args->sr, args->targets, args->targets_is_file, 0)<0 ) error("Failed to read the targets: %s\n",args->targets);
    if ( !bcf_sr_add_reader(args->sr,args->fname) ) error("Error: %s\n", bcf_sr_strerror(args->sr->errnum));
    args->hdr = bcf_sr_get_header(args->sr,0);
    args->hdr_out = bcf_hdr_dup(args->hdr);
    bcf_hdr_append(args->hdr_out, "##INFO=<ID=NOVELAL,Number=.,Type=String,Description=\"List of samples with novel alleles\">");
    bcf_hdr_append(args->hdr_out, "##INFO=<ID=NOVELGT,Number=.,Type=String,Description=\"List of samples with novel genotypes. Note that only samples w/o a novel allele are listed.\">");

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);

    int i;
    char **smpl = hts_readlist(args->bg_samples_str, 0, &args->nbg_smpl);
    args->bg_smpl = (int*) malloc(sizeof(int)*args->nbg_smpl);
    for (i=0; i<args->nbg_smpl; i++)
    {
        args->bg_smpl[i] = bcf_hdr_id2int(args->hdr, BCF_DT_SAMPLE, smpl[i]);
        if ( args->bg_smpl[i]<0 ) error("The sample not present in the VCF: \"%s\"\n", smpl[i]);
        free(smpl[i]);
    }
    free(smpl);

    smpl = hts_readlist(args->novel_samples_str, 0, &args->nnovel_smpl);
    args->novel_smpl = (int*) malloc(sizeof(int)*args->nnovel_smpl);
    for (i=0; i<args->nnovel_smpl; i++)
    {
        args->novel_smpl[i] = bcf_hdr_id2int(args->hdr, BCF_DT_SAMPLE, smpl[i]);
        if ( args->novel_smpl[i]<0 ) error("The sample not present in the VCF: \"%s\"\n", smpl[i]);
        free(smpl[i]);
    }
    free(smpl);

    args->out_fh = hts_open(args->output_fname,hts_bcf_wmode(args->output_type));
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
    bcf_hdr_write(args->out_fh, args->hdr_out);
}
static void destroy_data(args_t *args)
{
    bcf_hdr_destroy(args->hdr_out);
    hts_close(args->out_fh);
    free(args->novel_als_smpl.s);
    free(args->novel_gts_smpl.s);
    free(args->gts);
    free(args->bg_gts);
    free(args->bg_smpl);
    free(args->novel_smpl);
    if ( args->filter ) filter_destroy(args->filter);
    bcf_sr_destroy(args->sr);
    free(args);
}
static inline int binary_search(uint32_t val, uint32_t *dat, int ndat)
{
    int i = -1, imin = 0, imax = ndat - 1;
    while ( imin<=imax )
    {
        i = (imin+imax)/2;
        if ( dat[i] < val ) imin = i + 1;
        else if ( dat[i] > val ) imax = i - 1;
        else return 1;
    }
    return 0;
}
static inline void binary_insert(uint32_t val, uint32_t **dat, int *ndat, int *mdat)
{
    int i = -1, imin = 0, imax = *ndat - 1;
    while ( imin<=imax )
    {
        i = (imin+imax)/2;
        if ( (*dat)[i] < val ) imin = i + 1;
        else if ( (*dat)[i] > val ) imax = i - 1;
        else return;
    }
    while ( i>=0 && (*dat)[i]>val ) i--;

    (*ndat)++;
    hts_expand(uint32_t, (*ndat), (*mdat), (*dat));

    if ( *ndat > 1 )
        memmove(*dat + i + 1, *dat + i, sizeof(uint32_t)*(*ndat - i - 1));

    (*dat)[i+1] = val;
}
static int process_record(args_t *args, bcf1_t *rec)
{
    args->ntotal++;

    static int warned = 0;
    int ngts = bcf_get_genotypes(args->hdr, rec, &args->gts, &args->mgts);
    ngts /= rec->n_sample;
    if ( ngts>2 ) error("todo: ploidy=%d\n", ngts);

    args->nbg_gts = 0;
    uint32_t bg_als = 0;
    int i,j;
    for (i=0; i<args->nbg_smpl; i++)
    {
        uint32_t gt  = 0;
        int32_t *ptr = args->gts + args->bg_smpl[i]*ngts;
        for (j=0; j<ngts; j++)
        {
            if ( ptr[j]==bcf_int32_vector_end ) break;
            if ( bcf_gt_is_missing(ptr[j]) ) continue; 
            int ial = bcf_gt_allele(ptr[j]);
            if ( ial > 31 )
            {
                if ( !warned )
                {
                    fprintf(bcftools_stderr,"Too many alleles (>32) at %s:%d, skipping. (todo?)\n", bcf_seqname(args->hdr,rec),rec->pos+1);
                    warned = 1;
                }
                args->nskipped++;
                return -1;
            }
            bg_als |= 1<<ial;
            gt |= 1<<ial;
        }
        binary_insert(gt, &args->bg_gts, &args->nbg_gts, &args->mbg_gts);
    }
    if ( !bg_als )
    {
        // all are missing
        args->nskipped++;
        return -1;
    }

    args->novel_als_smpl.l = 0;
    args->novel_gts_smpl.l = 0;

    int has_gt = 0;
    for (i=0; i<args->nnovel_smpl; i++)
    {
        int novel_al = 0;
        uint32_t gt  = 0;
        int32_t *ptr = args->gts + args->novel_smpl[i]*ngts;
        for (j=0; j<ngts; j++)
        {
            if ( ptr[j]==bcf_int32_vector_end ) break;
            if ( bcf_gt_is_missing(ptr[j]) ) continue; 
            int ial = bcf_gt_allele(ptr[j]);
            if ( ial > 31 )
            {
                if ( !warned )
                {
                    fprintf(bcftools_stderr,"Too many alleles (>32) at %s:%d, skipping. (todo?)\n", bcf_seqname(args->hdr,rec),rec->pos+1);
                    warned = 1;
                }
                args->nskipped++;
                return -1;
            }
            if ( !(bg_als & (1<<ial)) ) novel_al = 1; 
            gt |= 1<<ial;
        }
        if ( !gt ) continue;
        has_gt = 1;

        char *smpl = args->hdr->samples[ args->novel_smpl[i] ];
        if ( novel_al )
        {
            if ( args->novel_als_smpl.l ) kputc(',', &args->novel_als_smpl);
            kputs(smpl, &args->novel_als_smpl);
        }
        else if ( !binary_search(gt, args->bg_gts, args->nbg_gts) )
        {
            if ( args->novel_gts_smpl.l ) kputc(',', &args->novel_gts_smpl);
            kputs(smpl, &args->novel_gts_smpl);
        }
    }
    if ( !has_gt )
    {
        // all are missing
        args->nskipped++;
        return -1;
    }
    if ( args->novel_als_smpl.l ) 
    {
        bcf_update_info_string(args->hdr_out, rec, "NOVELAL", args->novel_als_smpl.s);
        args->nnovel_al++;
    }
    if ( args->novel_gts_smpl.l ) 
    {
        bcf_update_info_string(args->hdr_out, rec, "NOVELGT", args->novel_gts_smpl.s);
        args->nnovel_gt++;
    }
    args->ntested++;
    return 0;
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->output_fname = "-";
    static struct option loptions[] =
    {
        {"bg-samples",required_argument,0,'0'},
        {"novel-samples",required_argument,0,'1'},
        {"include",required_argument,0,'i'},
        {"exclude",required_argument,0,'e'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {NULL,0,NULL,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "O:o:i:e:r:R:t:T:0:1:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case '0': args->bg_samples_str = optarg; break;
            case '1': args->novel_samples_str = optarg; break;
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 't': args->targets = optarg; break;
            case 'T': args->targets = optarg; args->targets_is_file = 1; break;
            case 'r': args->regions = optarg; break;
            case 'R': args->regions = optarg; args->regions_is_file = 1; break;
            case 'o': args->output_fname = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 'b': args->output_type = FT_BCF_GZ; break;
                          case 'u': args->output_type = FT_BCF; break;
                          case 'z': args->output_type = FT_VCF_GZ; break;
                          case 'v': args->output_type = FT_VCF; break;
                          default: error("The output type \"%s\" not recognised\n", optarg);
                      };
                      break;
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }
    if ( optind==argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args->fname = "-";  // reading from stdin
        else { error("%s",usage_text()); }
    }
    else if ( optind+1!=argc ) error("%s",usage_text());
    else args->fname = argv[optind];

    init_data(args);

    while ( bcf_sr_next_line(args->sr) )
    {
        bcf1_t *rec = bcf_sr_get_line(args->sr,0);
        if ( args->filter )
        {
            int pass = filter_test(args->filter, rec, NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( !pass ) continue;
        }
        process_record(args, rec);
        bcf_write(args->out_fh, args->hdr_out, rec);
    }

    fprintf(bcftools_stderr,"Total/processed/skipped/novel_allele/novel_gt:\t%d\t%d\t%d\t%d\t%d\n", args->ntotal, args->ntested, args->nskipped, args->nnovel_al, args->nnovel_gt);
    destroy_data(args);

    return 0;
}
