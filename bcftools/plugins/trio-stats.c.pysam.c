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

#define iCHILD  0
#define iFATHER 1
#define iMOTHER 2

typedef struct
{
    int idx[3];     // VCF sample index for father, mother and child
    int pass;       // do all three pass the filters?
}
trio_t;

typedef struct
{
    uint32_t
        npass,          // number of genotypes passing the filter
        nnon_ref,       // number of non-reference genotypes
        nmendel_err,    // number of mendelian errors
        nnovel,         // a singleton allele, but observed only in the child. Counted as mendel_err as well.
        nsingleton,     // het mother or father different from everyone else
        ndoubleton,     // het mother+child or father+child different from everyone else
        nts, ntv;       // number of transitions and transversions
}
trio_stats_t;

typedef struct
{
    trio_stats_t *stats;
    filter_t *filter;
    char *expr;
}
flt_stats_t;

typedef struct
{
    int argc, filter_logic, regions_is_file, targets_is_file;
    int nflt_str;
    char *filter_str, **flt_str;
    char **argv, *ped_fname, *output_fname, *fname, *regions, *targets;
    bcf_srs_t *sr;
    bcf_hdr_t *hdr;
    trio_t *trio;
    int ntrio, mtrio;
    flt_stats_t *filters;
    int nfilters;
    int32_t *gt_arr, *ac, *ac_trio;
    int mgt_arr, mac, mac_trio;
}
args_t;

args_t args;

const char *about(void)
{
    return "Calculate transmission rate and other stats in trio children.\n";
}

static const char *usage_text(void)
{
    return 
        "\n"
        "About: Calculate transmission rate in trio children. Use curly brackets to scan\n"
        "       a range of values simultaneously\n"
        "Usage: bcftools +trio-stats [Plugin Options]\n"
        "Plugin options:\n"
        "   -e, --exclude EXPR          exclude sites and samples for which the expression is true\n"
        "   -i, --include EXPR          include sites and samples for which the expression is true\n"
        "   -o, --output FILE           output file name [bcftools_stdout]\n"
        "   -p, --ped FILE              PED file\n"
        "   -r, --regions REG           restrict to comma-separated list of regions\n"
        "   -R, --regions-file FILE     restrict to regions listed in a file\n"
        "   -t, --targets REG           similar to -r but streams rather than index-jumps\n"
        "   -T, --targets-file FILE     similar to -R but streams rather than index-jumps\n"
        "\n"
        "Example:\n"
        "   bcftools +trio-stats -p file.ped -i 'GQ>{10,20,30,40,50}' file.bcf\n"
        "\n";
}

static int cmp_trios(const void *_a, const void *_b)
{
    trio_t *a = (trio_t *) _a;
    trio_t *b = (trio_t *) _b;
    int i;
    int amin = a->idx[0];
    for (i=1; i<3; i++)
        if ( amin > a->idx[i] ) amin = a->idx[i];
    int bmin = b->idx[0];
    for (i=1; i<3; i++)
        if ( bmin > b->idx[i] ) bmin = b->idx[i];
    if ( amin < bmin ) return -1;
    if ( amin > bmin ) return 1;
    return 0;
}

static void parse_ped(args_t *args, char *fname)
{
    htsFile *fp = hts_open(fname, "r");
    if ( !fp ) error("Could not read: %s\n", fname);

    kstring_t str = {0,0,0};
    if ( hts_getline(fp, KS_SEP_LINE, &str) <= 0 ) error("Empty file: %s\n", fname);

    int moff = 0, *off = NULL;
    do
    {
        // familyID    sampleID paternalID maternalID sex   phenotype   population relationship   siblings   secondOrder   thirdOrder   children    comment
        // BB03    HG01884 HG01885 HG01956 2   0   ACB child   0   0   0   0
        int ncols = ksplit_core(str.s,0,&moff,&off);
        if ( ncols<4 ) error("Could not parse the ped file: %s\n", str.s);

        int father = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[2]]);
        if ( father<0 ) continue;
        int mother = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[3]]);
        if ( mother<0 ) continue;
        int child = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,&str.s[off[1]]);
        if ( child<0 ) continue;

        args->ntrio++;
        hts_expand0(trio_t,args->ntrio,args->mtrio,args->trio);
        trio_t *trio = &args->trio[args->ntrio-1];
        trio->idx[iFATHER] = father;
        trio->idx[iMOTHER] = mother;
        trio->idx[iCHILD]  = child;
    }
    while ( hts_getline(fp, KS_SEP_LINE, &str)>=0 );

    fprintf(bcftools_stderr,"Identified %d complete trios in the VCF file\n", args->ntrio);

    // sort the sample by index so that they are accessed more or less sequentially
    qsort(args->trio,args->ntrio,sizeof(trio_t),cmp_trios);
    
    free(str.s);
    free(off);
    hts_close(fp);
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

    parse_ped(args, args->ped_fname);
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
    for (i=0; i<args->nfilters; i++)
        args->filters[i].stats = (trio_stats_t*) calloc(args->ntrio,sizeof(trio_stats_t));
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
    free(args->trio);
    free(args->ac);
    free(args->ac_trio);
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
    fprintf(fh,"# FLT* lines report numbers for every threshold and every trio:\n");
    fprintf(fh,"#   %d) filter id\n", ++i);
    fprintf(fh,"#   %d) child\n", ++i);
    fprintf(fh,"#   %d) father\n", ++i);
    fprintf(fh,"#   %d) mother\n", ++i);
    fprintf(fh,"#   %d) number of valid trio genotypes (all trio members pass filters, all non-missing)\n", ++i);
    fprintf(fh,"#   %d) number of non-reference trio GTs (at least one trio member carries an alternate allele)\n", ++i);
    fprintf(fh,"#   %d) number of Mendelian errors\n", ++i);
    fprintf(fh,"#   %d) number of novel singleton alleles in the child (counted also as a Mendelian error)\n", ++i);
    fprintf(fh,"#   %d) number of untransmitted singletons, present only in one parent\n", ++i);
    fprintf(fh,"#   %d) number of transmitted singletons, present only in one parent and the child\n", ++i);
    fprintf(fh,"#   %d) number of transitions, all ALT alleles present in the trio are considered\n", ++i);
    fprintf(fh,"#   %d) number of transversions, all ALT alleles present in the trio are considered\n", ++i);
    fprintf(fh,"#   %d) overall ts/tv, all ALT alleles present in the trio are considered\n", ++i);
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
        for (j=0; j<args->ntrio; j++)
        {
            fprintf(fh,"FLT%d", i);
            fprintf(fh,"\t%s",args->hdr->samples[args->trio[j].idx[iCHILD]]);
            fprintf(fh,"\t%s",args->hdr->samples[args->trio[j].idx[iFATHER]]);
            fprintf(fh,"\t%s",args->hdr->samples[args->trio[j].idx[iMOTHER]]);
            trio_stats_t *stats = &flt->stats[j];
            fprintf(fh,"\t%d", stats->npass);
            fprintf(fh,"\t%d", stats->nnon_ref);
            fprintf(fh,"\t%d", stats->nmendel_err);
            fprintf(fh,"\t%d", stats->nnovel);
            fprintf(fh,"\t%d", stats->nsingleton);
            fprintf(fh,"\t%d", stats->ndoubleton);
            fprintf(fh,"\t%d", stats->nts);
            fprintf(fh,"\t%d", stats->ntv);
            fprintf(fh,"\t%.2f", stats->ntv ? (float)stats->nts/stats->ntv : INFINITY);
            fprintf(fh,"\n");
        }
    }
    if ( fclose(fh)!=0 ) error("Close failed: %s\n", (!args->output_fname || !strcmp("-",args->output_fname)) ? "bcftools_stdout" : args->output_fname);
}

static inline int parse_genotype(int32_t *arr, int ngt1, int idx, int als[2])
{
    int32_t *ptr = arr + ngt1 * idx;
    if ( bcf_gt_is_missing(ptr[0]) ) return -1;
    als[0] = bcf_gt_allele(ptr[0]);

    // treat haploid GTs as homozygous diploid
    if ( ngt1==1 || ptr[1]==bcf_int32_vector_end ) { als[1] = als[0]; return 0; }

    if ( bcf_gt_is_missing(ptr[1]) ) return -1;
    als[1] = bcf_gt_allele(ptr[1]);

    return 0;
}

static void process_record(args_t *args, bcf1_t *rec, flt_stats_t *flt)
{
    int i,j;

    // Find out which trios pass and if the site passes
    if ( flt->filter )
    {
        uint8_t *smpl_pass = NULL;
        int pass_site = filter_test(flt->filter, rec, (const uint8_t**) &smpl_pass);
        if ( args->filter_logic & FLT_EXCLUDE )
        {
            if ( pass_site )
            {
                if ( !smpl_pass ) return;
                pass_site = 0;
                for (i=0; i<args->ntrio; i++)
                {
                    int pass_trio = 1;
                    for (j=0; j<3; j++)
                    {
                        int idx = args->trio[i].idx[j];
                        if ( smpl_pass[idx] ) { pass_trio = 0; break; }
                    }
                    args->trio[i].pass = pass_trio;
                    if ( pass_trio ) pass_site = 1;
                }
                if ( !pass_site ) return;
            }
            else
                for (i=0; i<args->ntrio; i++) args->trio[i].pass = 1;
        }
        else if ( !pass_site ) return;
        else if ( smpl_pass )
        {
            pass_site = 0;
            for (i=0; i<args->ntrio; i++)
            {
                int pass_trio = 1;
                for (j=0; j<3; j++)
                {
                    int idx = args->trio[i].idx[j];
                    if ( !smpl_pass[idx] ) { pass_trio = 0; break; }
                }
                args->trio[i].pass = pass_trio;
                if ( pass_trio ) pass_site = 1;
            }
            if ( !pass_site ) return;
        }
        else
            for (i=0; i<args->ntrio; i++) args->trio[i].pass = 1;
    }

    // Find out the allele counts. Try to use INFO/AC, if not present, determine from the genotypes
    hts_expand(int, rec->n_allele, args->mac, args->ac);
    if ( !bcf_calc_ac(args->hdr, rec, args->ac, BCF_UN_INFO|BCF_UN_FMT) ) return;
    hts_expand(int, rec->n_allele, args->mac_trio, args->ac_trio);

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
    for (i=0; i<args->ntrio; i++)
    {
        if ( flt->filter && !args->trio[i].pass ) continue;
        trio_stats_t *stats = &flt->stats[i];

        // Determine the alternate allele and the genotypes, skip if any of the alleles is missing.
        // the order is: child, father, mother
        int als[6], *als_child = als, *als_father = als+2, *als_mother = als+4; 
        if ( parse_genotype(args->gt_arr, ngt1, args->trio[i].idx[iCHILD], als_child) < 0 ) continue;
        if ( parse_genotype(args->gt_arr, ngt1, args->trio[i].idx[iFATHER], als_father) < 0 ) continue;
        if ( parse_genotype(args->gt_arr, ngt1, args->trio[i].idx[iMOTHER], als_mother) < 0 ) continue;

        stats->npass++;

        // Has the trio an alternate allele other than *?
        int has_star_allele = 0, has_nonref = 0;
        memset(args->ac_trio,0,rec->n_allele*sizeof(*args->ac_trio));
        for (j=0; j<6; j++)
        {
            if ( als[j]==star_allele ) { has_star_allele = 1; continue; }
            if ( als[j]==0 ) continue;
            has_nonref = 1;
            args->ac_trio[ als[j] ]++;
        }
        if ( !has_nonref ) continue;   // only ref or * in this trio
        
        stats->nnon_ref++;

        // Calculate ts/tv. It does the right thing and handles also HetAA genotypes
        if ( ref != -1 )
        {
            int has_ts = 0, has_tv = 0;
            for (j=0; j<6; j++)
            {
                if ( als[j]==0 || als[j]==star_allele ) continue;
                if ( als[j] >= rec->n_allele )
                    error("The GT index is out of range at %s:%d in %s\n", bcf_seqname(args->hdr,rec),rec->pos+1,args->hdr->samples[args->trio[i].idx[j/2]]);
                if ( rec->d.allele[als[j]][1] ) continue;

                int alt = bcf_acgt2int(rec->d.allele[als[j]][0]);
                if ( abs(ref-alt)==2 ) has_ts = 1;
                else has_tv = 1;
            }
            if ( has_ts ) stats->nts++;
            if ( has_tv ) stats->ntv++;
        }

        // Skip some stats if the star allele is present as it was already checked at the primary record, we do not want to count the same
        // thing multiple times. There can be other alternate allele, but we ignore that for simplicity.
        if ( has_star_allele ) continue;

        // Detect mendelian errors
        int mendel_ok = (als_child[0]==als_father[0] || als_child[0]==als_father[1]) && (als_child[1]==als_mother[0] || als_child[1]==als_mother[1]) ? 1 : 0;
        if ( !mendel_ok ) mendel_ok = (als_child[1]==als_father[0] || als_child[1]==als_father[1]) && (als_child[0]==als_mother[0] || als_child[0]==als_mother[1]) ? 1 : 0;
        if ( !mendel_ok ) stats->nmendel_err++;

        // Is this a singleton, doubleton, neither?
        for (j=1; j<rec->n_allele; j++)
        {
            if ( args->ac_trio[j]==1 && args->ac[j]==1 )  // singleton (in parent) or novel (in child)
            {
                if ( als_child[0]==j || als_child[1]==j ) stats->nnovel++;
                else stats->nsingleton++;
            }
            else if ( args->ac_trio[j]==2 && args->ac[j]==2 )   // possibly a doubleton
            {
                if ( (als_child[0]==j || als_child[1]==j) && (als_child[0]!=j || als_child[1]!=j) ) stats->ndoubleton++;
            }
        }
    }
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
        {"ped",required_argument,NULL,'p'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {NULL,0,NULL,0}
    };
    int c, i;
    while ((c = getopt_long(argc, argv, "p:o:s:i:e:r:R:t:T:",loptions,NULL)) >= 0)
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
            case 'p': args->ped_fname = optarg; break;
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

    if ( !args->ped_fname ) error("Missing the -p, --ped option\n");

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
