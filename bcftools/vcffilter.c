/*  vcffilter.c -- Apply fixed-threshold filters.

    Copyright (C) 2013-2014 Genome Research Ltd.

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
THE SOFTWARE.  */

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include "bcftools.h"
#include "filter.h"
#include "rbuf.h"

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

// FILTER columns annotation: replace or add to existing FILTERs; set FILTER to PASS at good sites?
#define ANNOT_ADD   1
#define ANNOT_RESET 2

// Set genotypes of filtered samples
#define SET_GTS_MISSING 1
#define SET_GTS_REF 2

typedef struct _args_t
{
    filter_t *filter;
    char *filter_str;
    int filter_logic;   // include or exclude sites which match the filters? One of FLT_INCLUDE/FLT_EXCLUDE
    const uint8_t *smpl_pass;
    int set_gts;
    char *soft_filter;  // drop failed sites or annotate FILTER column?
    int annot_mode;     // add to existing FILTER annotation or replace? Otherwise reset FILTER to PASS or leave as it is?
    int flt_fail, flt_pass;     // BCF ids of fail and pass filters
    int snp_gap, indel_gap, IndelGap_id, SnpGap_id;
    int32_t ntmpi, *tmpi, ntmp_ac, *tmp_ac;
    rbuf_t rbuf;
    bcf1_t **rbuf_lines;

    bcf_srs_t *files;
    bcf_hdr_t *hdr;
    htsFile *out_fh;
    int output_type, n_threads;

    char **argv, *output_fname, *targets_list, *regions_list;
    int argc, record_cmd_line;
}
args_t;

static void init_data(args_t *args)
{
    args->out_fh = hts_open(args->output_fname,hts_bcf_wmode(args->output_type));
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
    if ( args->n_threads ) hts_set_threads(args->out_fh, args->n_threads);

    args->hdr = args->files->readers[0].header;
    args->flt_pass = bcf_hdr_id2int(args->hdr,BCF_DT_ID,"PASS"); assert( !args->flt_pass );  // sanity check: required by BCF spec

    // -i or -e: append FILTER line
    if ( args->soft_filter && args->filter_logic )
    {
        kstring_t flt_name = {0,0,0};
        if ( strcmp(args->soft_filter,"+") )
            kputs(args->soft_filter, &flt_name);
        else
        {
            // Make up a filter name
            int i = 0, id = -1;
            do
            {
                ksprintf(&flt_name,"Filter%d", ++i);
                id = bcf_hdr_id2int(args->hdr,BCF_DT_ID,flt_name.s);
            }
            while ( bcf_hdr_idinfo_exists(args->hdr,BCF_HL_FLT,id) );
        }
        // escape quotes
        kstring_t tmp = {0,0,0};
        char *t = args->filter_str;
        while ( *t )
        {
            if ( *t=='"' ) kputc('\\',&tmp);
            kputc(*t,&tmp);
            t++;
        }
        int ret = bcf_hdr_printf(args->hdr, "##FILTER=<ID=%s,Description=\"Set if %s: %s\">", flt_name.s,args->filter_logic & FLT_INCLUDE ? "not true" : "true", tmp.s);
        if ( ret!=0 )
            error("Failed to append header line: ##FILTER=<ID=%s,Description=\"Set if %s: %s\">\n", flt_name.s,args->filter_logic & FLT_INCLUDE ? "not true" : "true", tmp.s);
        args->flt_fail = bcf_hdr_id2int(args->hdr,BCF_DT_ID,flt_name.s); assert( args->flt_fail>=0 );
        free(flt_name.s);
        free(tmp.s);
    }

    if ( args->snp_gap || args->indel_gap )
    {
        if ( !args->filter_logic && args->soft_filter && strcmp(args->soft_filter,"+") )
        {
            kstring_t tmp = {0,0,0};
            if ( args->snp_gap ) kputs("\"SnpGap\"", &tmp);
            if ( args->indel_gap )
            {
                if ( tmp.s ) kputs(" and ", &tmp);
                kputs("\"IndelGap\"", &tmp);
            }
            if ( strncmp(tmp.s+1,args->soft_filter,tmp.l-2) )
                fprintf(stderr,"Warning: using %s filter name instead of \"%s\"\n", tmp.s,args->soft_filter);
            free(tmp.s);
        }

        rbuf_init(&args->rbuf, 64);
        args->rbuf_lines = (bcf1_t**) calloc(args->rbuf.m, sizeof(bcf1_t*));
        if ( args->snp_gap )
        {
            bcf_hdr_printf(args->hdr, "##FILTER=<ID=SnpGap,Description=\"SNP within %d bp of an indel\">", args->snp_gap);
            args->SnpGap_id = bcf_hdr_id2int(args->hdr, BCF_DT_ID, "SnpGap");
            assert( args->SnpGap_id>=0 );
        }
        if ( args->indel_gap )
        {
            bcf_hdr_printf(args->hdr, "##FILTER=<ID=IndelGap,Description=\"Indel within %d bp of an indel\">", args->indel_gap);
            args->IndelGap_id = bcf_hdr_id2int(args->hdr, BCF_DT_ID, "IndelGap");
            assert( args->IndelGap_id>=0 );
        }
    }

    if (args->record_cmd_line) bcf_hdr_append_version(args->hdr, args->argc, args->argv, "bcftools_filter");

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);
}

static void destroy_data(args_t *args)
{
    if ( args->rbuf_lines )
    {
        int i;
        for (i=0; i<args->rbuf.m; i++)
            if ( args->rbuf_lines[i] ) bcf_destroy1(args->rbuf_lines[i]);
        free(args->rbuf_lines);
    }
    if ( args->filter )
        filter_destroy(args->filter);
    free(args->tmpi);
    free(args->tmp_ac);
}

static void flush_buffer(args_t *args, int n)
{
    int i, j;
    for (i=0; i<n; i++)
    {
        int k = rbuf_shift(&args->rbuf);
        bcf1_t *rec = args->rbuf_lines[k];

        int pass = 1;
        if ( !args->soft_filter )
        {
            for (j=0; j<rec->d.n_flt; j++)
            {
                if ( args->indel_gap && rec->d.flt[j]==args->IndelGap_id ) { pass = 0; break; }
                if ( args->snp_gap && rec->d.flt[j]==args->SnpGap_id ) { pass = 0; break; }
            }
        }
        if ( pass ) bcf_write1(args->out_fh, args->hdr, rec);
    }
}

#define SWAP(type_t, a, b) { type_t t = a; a = b; b = t; }
static void buffered_filters(args_t *args, bcf1_t *line)
{
    /**
     *  The logic of SnpGap=3. The SNPs at positions 1 and 7 are filtered,
     *  positions 0 and 8 are not:
     *           0123456789
     *      ref  .G.GT..G..
     *      del  .A.G-..A..
     *  Here the positions 1 and 6 are filtered, 0 and 7 are not:
     *           0123-456789
     *      ref  .G.G-..G..
     *      ins  .A.GT..A..
     *
     *  The logic of IndelGap=2. The second indel is filtered:
     *           012345678901
     *      ref  .GT.GT..GT..
     *      del  .G-.G-..G-..
     *  And similarly here, the second is filtered:
     *           01 23 456 78
     *      ref  .A-.A-..A-..
     *      ins  .AT.AT..AT..
     */

    // To avoid additional data structure, we abuse bcf1_t's var and var_type records.
    const int SnpGap_set     = VCF_OTHER<<1;
    const int IndelGap_set   = VCF_OTHER<<2;
    const int IndelGap_flush = VCF_OTHER<<3;

    int var_type = 0, i;
    if ( line )
    {
        // Still on the same chromosome?
        int ilast = rbuf_last(&args->rbuf);
        if ( ilast>=0 && line->rid != args->rbuf_lines[ilast]->rid )
            flush_buffer(args, args->rbuf.n); // new chromosome, flush everything

        rbuf_expand0(&args->rbuf,bcf1_t*,args->rbuf.n,args->rbuf_lines);

        // Insert the new record in the buffer. The line would be overwritten in
        // the next bcf_sr_next_line call, therefore we need to swap it with an
        // unused one
        ilast = rbuf_append(&args->rbuf);
        if ( !args->rbuf_lines[ilast] ) args->rbuf_lines[ilast] = bcf_init1();
        SWAP(bcf1_t*, args->files->readers[0].buffer[0], args->rbuf_lines[ilast]);

        var_type = bcf_get_variant_types(line);

        // Find out the size of an indel. The indel boundaries are based on REF
        // (POS+1,POS+rlen-1). This is not entirely correct: mpileup likes to
        // output REF=CAGAGAGAGA, ALT=CAGAGAGAGAGA where REF=C,ALT=CGA could be
        // used. This filter is therefore more strict and may remove some valid
        // SNPs.
        int len = 1;
        if ( var_type & VCF_INDEL )
        {
            for (i=1; i<line->n_allele; i++)
                if ( len < 1-line->d.var[i].n ) len = 1-line->d.var[i].n;
        }

        // Set the REF allele's length to max deletion length or to 1 if a SNP or an insertion.
        line->d.var[0].n = len;
    }

    int k_flush = 1;
    if ( args->indel_gap )
    {
        k_flush = 0;
        // Find indels which are too close to each other
        int last_to = -1;
        for (i=-1; rbuf_next(&args->rbuf,&i); )
        {
            bcf1_t *rec  = args->rbuf_lines[i];
            int rec_from = rec->pos;
            if ( last_to!=-1 && last_to < rec_from ) break;

            k_flush++;
            if ( !(rec->d.var_type & VCF_INDEL) ) continue;

            rec->d.var_type |= IndelGap_set;
            last_to = args->indel_gap + rec->pos + rec->d.var[0].n - 1;
        }
        if ( i==args->rbuf.f && line && last_to!=-1 ) k_flush = 0;
        if ( k_flush || !line )
        {
            // Select the best indel from the cluster of k_flush indels
            int k = 0, max_ac = -1, imax_ac = -1;
            for (i=-1; rbuf_next(&args->rbuf,&i) && k<k_flush; )
            {
                k++;
                bcf1_t *rec  = args->rbuf_lines[i];
                if ( !(rec->d.var_type & IndelGap_set) ) continue;
                hts_expand(int, rec->n_allele, args->ntmpi, args->tmpi);
                int ret = bcf_calc_ac(args->hdr, rec, args->tmpi, BCF_UN_ALL);
                if ( imax_ac==-1 || (ret && max_ac < args->tmpi[1]) ) { max_ac = args->tmpi[1]; imax_ac = i; }
            }

            // Filter all but the best indel (with max AF or first if AF not available)
            k = 0;
            for (i=-1; rbuf_next(&args->rbuf,&i) && k<k_flush; )
            {
                k++;
                bcf1_t *rec = args->rbuf_lines[i];
                if ( !(rec->d.var_type & IndelGap_set) ) continue;
                rec->d.var_type |= IndelGap_flush;
                if ( i!=imax_ac ) bcf_add_filter(args->hdr, args->rbuf_lines[i], args->IndelGap_id);
            }
        }
    }

    if ( !line )
    {
        // Finished: flush everything
        flush_buffer(args, args->rbuf.n);
        return;
    }

    int j_flush = 1;
    if ( args->snp_gap )
    {
        j_flush = 0;
        int last_from = line->pos;
        for (i=-1; rbuf_next(&args->rbuf,&i); )
        {
            bcf1_t *rec = args->rbuf_lines[i];
            int rec_to  = rec->pos + rec->d.var[0].n - 1;   // last position affected by the variant
            if ( rec_to + args->snp_gap < last_from )
                j_flush++;
            else if ( (var_type & VCF_INDEL) && (rec->d.var_type & VCF_SNP) && !(rec->d.var_type & SnpGap_set) )
            {
                // this SNP has not been SnpGap-filtered yet
                rec->d.var_type |= SnpGap_set;
                bcf_add_filter(args->hdr, rec, args->SnpGap_id);
            }
            else if ( (var_type & VCF_SNP) && (rec->d.var_type & VCF_INDEL) )
            {
                // the line which we are adding is a SNP and needs to be filtered
                line->d.var_type |= SnpGap_set;
                bcf_add_filter(args->hdr, line, args->SnpGap_id);
                break;
            }
        }
    }
    flush_buffer(args, j_flush < k_flush ? j_flush : k_flush);
}

static void set_genotypes(args_t *args, bcf1_t *line, int pass_site)
{
    int i,j;
    if ( !bcf_hdr_nsamples(args->hdr) ) return;
    if ( args->smpl_pass )
    {
        int npass = 0;
        for (i=0; i<bcf_hdr_nsamples(args->hdr); i++) npass += args->smpl_pass[i];

        // return if all samples pass
        if ( npass==bcf_hdr_nsamples(args->hdr) && (args->filter_logic & FLT_INCLUDE) ) return;
        if ( npass==0 && (args->filter_logic & FLT_EXCLUDE) ) return;
    }
    else if ( pass_site ) return;

    int an = 0, has_an = bcf_get_info_int32(args->hdr, line, "AN", &args->tmp_ac, &args->ntmp_ac);
    if ( has_an==1 ) an = args->tmp_ac[0];
    else has_an = 0;

    int has_ac = bcf_get_info_int32(args->hdr, line, "AC", &args->tmp_ac, &args->ntmp_ac);
    has_ac = has_ac==line->n_allele-1 ? 1 : 0;

    int new_gt = 0, ngts = bcf_get_format_int32(args->hdr, line, "GT", &args->tmpi, &args->ntmpi);
    ngts /= bcf_hdr_nsamples(args->hdr);
    if ( args->set_gts==SET_GTS_MISSING ) new_gt = bcf_gt_missing;
    else if ( args->set_gts==SET_GTS_REF ) new_gt = bcf_gt_unphased(0);
    else error("todo: set_gts=%d\n", args->set_gts);
    for (i=0; i<bcf_hdr_nsamples(args->hdr); i++)
    {
        if ( args->smpl_pass )
        {
            int pass = args->smpl_pass[i];
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( pass ) continue;
        }
        int32_t *gts = args->tmpi + ngts*i;
        for (j=0; j<ngts; j++)
        {
            if ( gts[j]==bcf_int32_vector_end ) break;
            if ( args->set_gts==SET_GTS_MISSING && !bcf_gt_is_missing(gts[j]) )
            {
                int ial = bcf_gt_allele(gts[j]);
                if ( has_ac && ial>0 && ial<=line->n_allele ) args->tmp_ac[ ial-1 ]--;
                an--;
            }
            else if ( args->set_gts==SET_GTS_REF )
            {
                int ial = bcf_gt_allele(gts[j]);
                if ( bcf_gt_is_missing(gts[j]) ) an++;
                else if ( has_ac && ial>0 && ial<=line->n_allele ) args->tmp_ac[ ial-1 ]--;
            }
            gts[j] = new_gt;
        }
    }
    bcf_update_genotypes(args->hdr,line,args->tmpi,ngts*bcf_hdr_nsamples(args->hdr));
    if ( has_an ) bcf_update_info_int32(args->hdr,line,"AN",&an,1);
    if ( has_ac )  bcf_update_info_int32(args->hdr,line,"AC",args->tmp_ac,line->n_allele-1);
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Apply fixed-threshold filters.\n");
    fprintf(stderr, "Usage:   bcftools filter [options] <in.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -e, --exclude <expr>          exclude sites for which the expression is true (see man page for details)\n");
    fprintf(stderr, "    -g, --SnpGap <int>            filter SNPs within <int> base pairs of an indel\n");
    fprintf(stderr, "    -G, --IndelGap <int>          filter clusters of indels separated by <int> or fewer base pairs allowing only one to pass\n");
    fprintf(stderr, "    -i, --include <expr>          include only sites for which the expression is true (see man page for details\n");
    fprintf(stderr, "    -m, --mode [+x]               \"+\": do not replace but add to existing FILTER; \"x\": reset filters at sites which pass\n");
    fprintf(stderr, "        --no-version              do not append version and command line to the header\n");
    fprintf(stderr, "    -o, --output <file>           write output to a file [standard output]\n");
    fprintf(stderr, "    -O, --output-type <b|u|z|v>   b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "    -r, --regions <region>        restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file <file>     restrict to regions listed in a file\n");
    fprintf(stderr, "    -s, --soft-filter <string>    annotate FILTER column with <string> or unique filter name (\"Filter%%d\") made up by the program (\"+\")\n");
    fprintf(stderr, "    -S, --set-GTs <.|0>           set genotypes of failed samples to missing (.) or ref (0)\n");
    fprintf(stderr, "    -t, --targets <region>        similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "    -T, --targets-file <file>     similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "        --threads <int>           number of extra output compression threads [0]\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcffilter(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->output_fname = "-";
    args->output_type = FT_VCF;
    args->n_threads = 0;
    args->record_cmd_line = 1;
    int regions_is_file = 0, targets_is_file = 0;

    static struct option loptions[] =
    {
        {"set-GTs",required_argument,NULL,'S'},
        {"mode",required_argument,NULL,'m'},
        {"soft-filter",required_argument,NULL,'s'},
        {"exclude",required_argument,NULL,'e'},
        {"include",required_argument,NULL,'i'},
        {"targets",required_argument,NULL,'t'},
        {"targets-file",required_argument,NULL,'T'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"threads",required_argument,NULL,9},
        {"SnpGap",required_argument,NULL,'g'},
        {"IndelGap",required_argument,NULL,'G'},
        {"no-version",no_argument,NULL,8},
        {NULL,0,NULL,0}
    };
    char *tmp;
    while ((c = getopt_long(argc, argv, "e:i:t:T:r:R:h?s:m:o:O:g:G:S:",loptions,NULL)) >= 0) {
        switch (c) {
            case 'g': 
                args->snp_gap = strtol(optarg,&tmp,10); 
                if ( *tmp ) error("Could not parse argument: --SnpGap %s\n", optarg);
                break;
            case 'G':
                args->indel_gap = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --IndelGap %s\n", optarg);
                break;
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
            case 's': args->soft_filter = optarg; break;
            case 'm':
                if ( strchr(optarg,'x') ) args->annot_mode |= ANNOT_RESET;
                if ( strchr(optarg,'+') ) args->annot_mode |= ANNOT_ADD;
                break;
            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; targets_is_file = 1; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'S':
                if ( !strcmp(".",optarg) ) args->set_gts = SET_GTS_MISSING;
                else if ( !strcmp("0",optarg) ) args->set_gts = SET_GTS_REF;
                else error("The argument to -S not recognised: %s\n", optarg);
                break;
            case  9 : args->n_threads = strtol(optarg, 0, 0); break;
            case  8 : args->record_cmd_line = 0; break;
            case 'h':
            case '?': usage(args);
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( args->filter_logic == (FLT_EXCLUDE|FLT_INCLUDE) ) error("Only one of -i or -e can be given.\n");
    char *fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin
        else usage(args);
    }
    else fname = argv[optind];

    // read in the regions from the command line
    if ( args->regions_list )
    {
        args->files->require_index = 1;
        if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    else if ( optind+1 < argc )
    {
        int i;
        kstring_t tmp = {0,0,0};
        kputs(argv[optind+1],&tmp);
        for (i=optind+2; i<argc; i++) { kputc(',',&tmp); kputs(argv[i],&tmp); }
        args->files->require_index = 1;
        if ( bcf_sr_set_regions(args->files, tmp.s, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
        free(tmp.s);
    }
    if ( args->targets_list )
    {
        if ( bcf_sr_set_targets(args->files, args->targets_list,targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
    }
    if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open %s: %s\n", fname,bcf_sr_strerror(args->files->errnum));

    init_data(args);
    bcf_hdr_write(args->out_fh, args->hdr);
    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = bcf_sr_get_line(args->files, 0);
        int pass = 1;
        if ( args->filter )
        {
            pass = filter_test(args->filter, line, &args->smpl_pass);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
        }
        if ( args->soft_filter || args->set_gts || pass )
        {
            if ( pass )
            {
                bcf_unpack(line,BCF_UN_FLT);
                if ( args->annot_mode & ANNOT_RESET || !line->d.n_flt ) bcf_add_filter(args->hdr, line, args->flt_pass);
            }
            else if ( args->soft_filter )
            {
                if ( (args->annot_mode & ANNOT_ADD) ) bcf_add_filter(args->hdr, line, args->flt_fail);
                else bcf_update_filter(args->hdr, line, &args->flt_fail, 1);
            }
            if ( args->set_gts ) set_genotypes(args, line, pass);
            if ( !args->rbuf_lines )
                bcf_write1(args->out_fh, args->hdr, line);
            else
                buffered_filters(args, line);
        }
    }
    buffered_filters(args, NULL);

    hts_close(args->out_fh);
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}
