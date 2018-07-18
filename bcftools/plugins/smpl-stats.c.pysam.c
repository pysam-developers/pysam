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
#include <unistd.h>     // for isatty
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include "bcftools.h"
#include "filter.h"


// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef struct
{
    uint32_t
        npass,          // number of genotypes passing the filter
        nnon_ref,       // number of non-reference genotypes
        nhomRR,
        nhomAA,
        nhemi,
        nhet,
        nSNV,
        nIndel,
        nmissing,
        nsingleton,     // het different from everyone else
        nts, ntv;       // number of transitions and transversions
}
stats_t;

typedef struct
{
    stats_t *stats, site_stats;
    filter_t *filter;
    char *expr;
}
flt_stats_t;

typedef struct
{
    int argc, filter_logic, regions_is_file, targets_is_file;
    int nflt_str;
    char *filter_str, **flt_str;
    char **argv, *output_fname, *fname, *regions, *targets;
    bcf_srs_t *sr;
    bcf_hdr_t *hdr;
    flt_stats_t *filters;
    int nfilters, nsmpl;
    int32_t *gt_arr, *ac;
    int mgt_arr, mac;
}
args_t;

args_t args;

const char *about(void)
{
    return "Calculate basic per-sample stats scanning over a range of thresholds simultaneously.\n";
}

static const char *usage_text(void)
{
    return 
        "\n"
        "About: Calculates basic per-sample stats. Use curly brackets to scan a range of values simultaneously\n"
        "Usage: bcftools +smpl-stats [Plugin Options]\n"
        "Plugin options:\n"
        "   -e, --exclude EXPR          exclude sites and samples for which the expression is true\n"
        "   -i, --include EXPR          include sites and samples for which the expression is true\n"
        "   -o, --output FILE           output file name [bcftools_stdout]\n"
        "   -r, --regions REG           restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE     restrict to regions listed in a file\n"
        "   -t, --targets REG           similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE     similar to -R but streams rather than index-jumps\n"
        "\n"
        "Example:\n"
        "   bcftools +smpl-stats -i 'GQ>{10,20,30,40,50}' file.bcf\n"
        "\n";
}

static void parse_filters(args_t *args)
{
    if ( !args->filter_str ) return;
    int mflt = 1;
    args->nflt_str = 1;
    args->flt_str  = (char**) malloc(sizeof(char*));
    args->flt_str[0] = strdup(args->filter_str);
    while (1)
    {
        int i, expanded = 0;
        for (i=args->nflt_str-1; i>=0; i--)
        {
            char *exp_beg = strchr(args->flt_str[i], '{');
            if ( !exp_beg ) continue;
            char *exp_end = strchr(exp_beg+1, '}');
            if ( !exp_end ) error("Could not parse the expression: %s\n", args->filter_str);
            char *beg = exp_beg+1, *mid = beg;
            while ( mid<exp_end )
            {
                while ( mid<exp_end && *mid!=',' ) mid++;
                kstring_t tmp = {0,0,0};
                kputsn(args->flt_str[i], exp_beg - args->flt_str[i], &tmp);
                kputsn(beg, mid - beg, &tmp);
                kputs(exp_end+1, &tmp);
                args->nflt_str++;
                hts_expand(char*, args->nflt_str, mflt, args->flt_str);
                args->flt_str[args->nflt_str-1] = tmp.s;
                beg = ++mid;
            }
            expanded = 1;
            free(args->flt_str[i]);
            memmove(&args->flt_str[i], &args->flt_str[i+1], (args->nflt_str-i-1)*sizeof(*args->flt_str));
            args->nflt_str--;
            args->flt_str[args->nflt_str] = NULL;
        }
        if ( !expanded ) break;
    }
    
    fprintf(bcftools_stderr,"Collecting data for %d filtering expressions\n", args->nflt_str);
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

    parse_filters(args);

    int i;
    if ( !args->nflt_str )
    {
        args->filters = (flt_stats_t*) calloc(1, sizeof(flt_stats_t));
        args->nfilters = 1;
        args->filters[0].expr = strdup("all");
    }
    else
    {
        args->nfilters = args->nflt_str;
        args->filters = (flt_stats_t*) calloc(args->nfilters, sizeof(flt_stats_t));
        for (i=0; i<args->nfilters; i++)
        {
            args->filters[i].filter = filter_init(args->hdr, args->flt_str[i]);
            args->filters[i].expr   = strdup(args->flt_str[i]);

            // replace tab's with spaces so that the output stays parsable
            char *tmp = args->filters[i].expr;
            while ( *tmp )
            { 
                if ( *tmp=='\t' ) *tmp = ' '; 
                tmp++; 
            }
        }
    }
    args->nsmpl = bcf_hdr_nsamples(args->hdr);
    for (i=0; i<args->nfilters; i++)
        args->filters[i].stats = (stats_t*) calloc(args->nsmpl,sizeof(stats_t));
}
static void destroy_data(args_t *args)
{
    int i;
    for (i=0; i<args->nfilters; i++)
    {
        if ( args->filters[i].filter ) filter_destroy(args->filters[i].filter);
        free(args->filters[i].stats);
        free(args->filters[i].expr);
    }
    free(args->filters);
    for (i=0; i<args->nflt_str; i++) free(args->flt_str[i]);
    free(args->flt_str);
    bcf_sr_destroy(args->sr);
    free(args->ac);
    free(args->gt_arr);
    free(args);
}
static void report_stats(args_t *args)
{
    int i = 0,j;
    FILE *fh = !args->output_fname || !strcmp("-",args->output_fname) ? bcftools_stdout : fopen(args->output_fname,"w");
    if ( !fh ) error("Could not open the file for writing: %s\n", args->output_fname);
    fprintf(fh,"# CMD line shows the command line used to generate this output\n");
    fprintf(fh,"# DEF lines define expressions for all tested thresholds\n");
    fprintf(fh,"# FLT* lines report numbers for every threshold and every sample:\n");
    fprintf(fh,"#   %d) filter id\n", ++i);
    fprintf(fh,"#   %d) sample\n", ++i);
    fprintf(fh,"#   %d) number of genotypes which pass the filter\n", ++i);
    fprintf(fh,"#   %d) number of non-reference genotypes\n", ++i);
    fprintf(fh,"#   %d) number of homozygous ref genotypes (0/0 or 0)\n", ++i);
    fprintf(fh,"#   %d) number of homozygous alt genotypes (1/1, 2/2, etc)\n", ++i);
    fprintf(fh,"#   %d) number of heterozygous genotypes (0/1, 1/2, etc)\n", ++i);
    fprintf(fh,"#   %d) number of hemizygous genotypes (0, 1, etc)\n", ++i);
    fprintf(fh,"#   %d) number of SNVs\n", ++i);
    fprintf(fh,"#   %d) number of indels\n", ++i);
    fprintf(fh,"#   %d) number of singletons\n", ++i);
    fprintf(fh,"#   %d) number of missing genotypes (./., ., ./0, etc)\n", ++i);
    fprintf(fh,"#   %d) number of transitions (genotypes such as \"1/2\" are counted twice)\n", ++i);
    fprintf(fh,"#   %d) number of transversions (genotypes such as \"1/2\" are counted twice)\n", ++i);
    fprintf(fh,"#   %d) overall ts/tv\n", ++i);
    i = 0;
    fprintf(fh,"# SITE* lines report numbers for every threshold and site:\n");
    fprintf(fh,"#   %d) filter id\n", ++i);
    fprintf(fh,"#   %d) number of sites which pass the filter\n", ++i);
    fprintf(fh,"#   %d) number of SNVs\n", ++i);
    fprintf(fh,"#   %d) number of indels\n", ++i);
    fprintf(fh,"#   %d) number of singletons\n", ++i);
    fprintf(fh,"#   %d) number of transitions (counted at most once at multiallelic sites)\n", ++i);
    fprintf(fh,"#   %d) number of transversions (counted at most once at multiallelic sites)\n", ++i);
    fprintf(fh,"#   %d) overall ts/tv\n", ++i);
    fprintf(fh, "CMD\t%s", args->argv[0]);
    for (i=1; i<args->argc; i++) fprintf(fh, " %s",args->argv[i]);
    fprintf(fh, "\n");
    for (i=0; i<args->nfilters; i++)
    {
        flt_stats_t *flt = &args->filters[i];
        fprintf(fh,"DEF\tFLT%d\t%s\n", i, flt->expr);
    }
    for (i=0; i<args->nfilters; i++)
    {
        flt_stats_t *flt = &args->filters[i];
        for (j=0; j<args->nsmpl; j++)
        {
            fprintf(fh,"FLT%d", i);
            fprintf(fh,"\t%s",args->hdr->samples[j]);
            stats_t *stats = &flt->stats[j];
            fprintf(fh,"\t%d", stats->npass);
            fprintf(fh,"\t%d", stats->nnon_ref);
            fprintf(fh,"\t%d", stats->nhomRR);
            fprintf(fh,"\t%d", stats->nhomAA);
            fprintf(fh,"\t%d", stats->nhet);
            fprintf(fh,"\t%d", stats->nhemi);
            fprintf(fh,"\t%d", stats->nSNV);
            fprintf(fh,"\t%d", stats->nIndel);
            fprintf(fh,"\t%d", stats->nsingleton);
            fprintf(fh,"\t%d", stats->nmissing);
            fprintf(fh,"\t%d", stats->nts);
            fprintf(fh,"\t%d", stats->ntv);
            fprintf(fh,"\t%.2f", stats->ntv ? (float)stats->nts/stats->ntv : INFINITY);
            fprintf(fh,"\n");
        }
        fprintf(fh,"SITE%d", i);
        stats_t *stats = &flt->site_stats;
        fprintf(fh,"\t%d", stats->npass);
        fprintf(fh,"\t%d", stats->nSNV);
        fprintf(fh,"\t%d", stats->nIndel);
        fprintf(fh,"\t%d", stats->nsingleton);
        fprintf(fh,"\t%d", stats->nts);
        fprintf(fh,"\t%d", stats->ntv);
        fprintf(fh,"\t%.2f", stats->ntv ? (float)stats->nts/stats->ntv : INFINITY);
        fprintf(fh,"\n");
    }
    if ( fclose(fh)!=0 ) error("Close failed: %s\n", (!args->output_fname || !strcmp("-",args->output_fname)) ? "bcftools_stdout" : args->output_fname);
}

static inline int parse_genotype(int32_t *arr, int ngt1, int idx, int als[2])
{
    int32_t *ptr = arr + ngt1 * idx;
    if ( bcf_gt_is_missing(ptr[0]) ) return -1;
    als[0] = bcf_gt_allele(ptr[0]);

    if ( ngt1==1 || ptr[1]==bcf_int32_vector_end ) { ptr[1] = ptr[0]; return -2; }

    if ( bcf_gt_is_missing(ptr[1]) ) return -1;
    als[1] = bcf_gt_allele(ptr[1]);

    return 0;
}

static void process_record(args_t *args, bcf1_t *rec, flt_stats_t *flt)
{
    int i,j;
    uint8_t *smpl_pass = NULL;

    // Find out which trios pass and if the site passes
    if ( flt->filter )
    {
        int pass_site = filter_test(flt->filter, rec, (const uint8_t**) &smpl_pass);
        if ( args->filter_logic & FLT_EXCLUDE )
        {
            if ( pass_site )
            {
                if ( !smpl_pass ) return;
                pass_site = 0;
                for (i=0; i<args->nsmpl; i++)
                {
                    if ( smpl_pass[i] ) smpl_pass[i] = 0;
                    else { smpl_pass[i] = 1; pass_site = 1; }
                }
                if ( !pass_site ) return;
            }
            else
                for (i=0; i<args->nsmpl; i++) smpl_pass[i] = 1;
        }
        else if ( !pass_site ) return;
    }

    // Find out the allele counts. Try to use INFO/AC, if not present, determine from the genotypes
    hts_expand(int, rec->n_allele, args->mac, args->ac);
    if ( !bcf_calc_ac(args->hdr, rec, args->ac, BCF_UN_INFO|BCF_UN_FMT) ) return;

    // Get the genotypes
    int ngt = bcf_get_genotypes(args->hdr, rec, &args->gt_arr, &args->mgt_arr);
    if ( ngt<0 ) return;
    int ngt1 = ngt / rec->n_sample;
    

    // For ts/tv: numeric code of the reference allele, -1 for insertions
    int ref = !rec->d.allele[0][1] ? bcf_acgt2int(*rec->d.allele[0]) : -1;

    int star_allele = -1;
    for (i=1; i<rec->n_allele; i++)
        if ( !rec->d.allele[i][1] && rec->d.allele[i][0]=='*' ) { star_allele = i; break; }

    // Run the stats
    int site_pass      = 0;
    int site_SNV       = 0;
    int site_Indel     = 0;
    int site_has_ts    = 0;
    int site_has_tv    = 0;
    int site_singleton = 0;
    for (i=0; i<args->nsmpl; i++)
    {
        if ( smpl_pass && !smpl_pass[i] ) continue;
        stats_t *stats = &flt->stats[i];

        // Determine the alternate allele and the genotypes, skip if any of the alleles is missing.
        int als[2];
        int ret = parse_genotype(args->gt_arr, ngt1, i, als);
        if ( ret==-1 ) { stats->nmissing++; continue; }   // missing allele
        if ( ret==-2 ) stats->nhemi++;
        else if ( als[0]!=als[1] ) stats->nhet++;
        else if ( als[0]==0 ) stats->nhomRR++;
        else stats->nhomAA++;

        stats->npass++;
        site_pass = 1;

        // Is there an alternate allele other than *?
        int has_nonref = 0;
        for (j=0; j<2; j++)
        {
            if ( als[j]==star_allele ) continue;
            if ( als[j]==0 ) continue;
            has_nonref = 1;
        }
        if ( !has_nonref ) continue; // only ref or * in this genotype
        
        stats->nnon_ref++;

        // Calculate ts/tv, count SNPs, indels. It does the right thing and handles also HetAA genotypes
        {
            int has_ts = 0, has_tv = 0, has_snv = 0, has_indel = 0;
            for (j=0; j<2; j++)
            {
                if ( als[j]==0 || als[j]==star_allele ) continue;
                if ( als[j] >= rec->n_allele )
                    error("The GT index is out of range at %s:%d in %s\n", bcf_seqname(args->hdr,rec),rec->pos+1,args->hdr->samples[j]);

                if ( args->ac[als[j]]==1 ) { stats->nsingleton++; site_singleton = 1; }

                int var_type = bcf_get_variant_type(rec, als[j]);
                if ( var_type==VCF_SNP || var_type==VCF_MNP )
                {
                    int k = 0;
                    while ( rec->d.allele[0][k] && rec->d.allele[als[j]][k] )
                    {
                        if ( rec->d.allele[0][k]==rec->d.allele[als[j]][k] ) { k++; continue; }

                        int alt = bcf_acgt2int(rec->d.allele[als[j]][k]);
                        if ( abs(ref-alt)==2 ) has_ts = 1;
                        else has_tv = 1;
                        has_snv = 1;

                        k++;
                    }
                }
                else if ( var_type==VCF_INDEL ) has_indel = 1;
            }
            if ( has_ts ) { stats->nts++; site_has_ts = 1; }
            if ( has_tv ) { stats->ntv++; site_has_tv = 1; }
            if ( has_snv ) { stats->nSNV++; site_SNV = 1; }
            if ( has_indel ) { stats->nIndel++; site_Indel = 1; }
        }
    }
    flt->site_stats.npass  += site_pass;
    flt->site_stats.nSNV   += site_SNV;
    flt->site_stats.nIndel += site_Indel;
    flt->site_stats.nts    += site_has_ts;
    flt->site_stats.ntv    += site_has_tv;
    flt->site_stats.nsingleton += site_singleton;
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->output_fname = "-";
    static struct option loptions[] =
    {
        {"include",required_argument,0,'i'},
        {"exclude",required_argument,0,'e'},
        {"output",required_argument,NULL,'o'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {NULL,0,NULL,0}
    };
    int c, i;
    while ((c = getopt_long(argc, argv, "o:s:i:e:r:R:t:T:",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 't': args->targets = optarg; break;
            case 'T': args->targets = optarg; args->targets_is_file = 1; break;
            case 'r': args->regions = optarg; break;
            case 'R': args->regions = optarg; args->regions_is_file = 1; break;
            case 'o': args->output_fname = optarg; break;
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
        for (i=0; i<args->nfilters; i++)
            process_record(args, rec, &args->filters[i]);
    }

    report_stats(args);
    destroy_data(args);

    return 0;
}
