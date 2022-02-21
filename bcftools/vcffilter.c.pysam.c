#include "bcftools.pysam.h"

/*  vcffilter.c -- Apply fixed-threshold filters.

    Copyright (C) 2013-2022 Genome Research Ltd.

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
#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
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
#include "regidx.h"

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
    int snp_gap, snp_gap_type, indel_gap, IndelGap_id, SnpGap_id;
    char *snp_gap_str;
    int32_t ntmpi, *tmpi, ntmp_ac, *tmp_ac;
    rbuf_t rbuf;
    bcf1_t **rbuf_lines;

    bcf_srs_t *files;
    bcf_hdr_t *hdr;
    htsFile *out_fh;
    int output_type, n_threads, clevel;

    char **argv, *output_fname, *targets_list, *regions_list, *mask_list;
    int argc, record_cmd_line, mask_is_file, mask_overlap, mask_negate;
    regidx_t *mask;
}
args_t;

static void init_data(args_t *args)
{
    char wmode[8];
    set_wmode(wmode,args->output_type,args->output_fname,args->clevel);
    args->out_fh = hts_open(args->output_fname ? args->output_fname : "-", wmode);
    if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
    if ( args->n_threads ) hts_set_threads(args->out_fh, args->n_threads);

    if ( args->mask_list )
    {
        if ( args->mask_list[0]=='^' ) args->mask_negate = 1;
        if ( args->mask_is_file )
            args->mask = regidx_init(args->mask_negate?args->mask_list+1:args->mask_list,NULL,NULL,0,NULL);
        else
        {
            char *rmme = strdup(args->mask_negate?args->mask_list+1:args->mask_list), *tmp = rmme;
            while ( *tmp )
            {
                if ( *tmp==',' ) *tmp = '\n';
                tmp++;
            }
            args->mask = regidx_init_string(rmme, regidx_parse_reg, NULL, 0, NULL);
            free(rmme);
        }
        if ( !args->mask )
            error("Could not initialize the mask: %s\n",args->mask_list);
    }

    args->hdr = args->files->readers[0].header;
    args->flt_pass = bcf_hdr_id2int(args->hdr,BCF_DT_ID,"PASS"); assert( !args->flt_pass );  // sanity check: required by BCF spec

    if ( args->soft_filter && (args->filter_logic || args->mask_list) )
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

        kstring_t tmp = {0,0,0};
        if ( args->filter_logic )
        {
            // -i or -e: append FILTER line
            ksprintf(&tmp,"Set if %s: ",args->filter_logic & FLT_INCLUDE ? "not true" : "true");

            // escape quotes
            char *t = args->filter_str;
            while ( *t )
            {
                if ( *t=='"' ) kputc('\\',&tmp);
                kputc(*t,&tmp);
                t++;
            }
        }
        else if ( args->mask_list )
            ksprintf(&tmp,"Record masked by region");

        int ret = bcf_hdr_printf(args->hdr, "##FILTER=<ID=%s,Description=\"%s\">", flt_name.s,tmp.s);
        if ( ret!=0 )
            error("Failed to append header line: ##FILTER=<ID=%s,Description=\"%s\">\n", flt_name.s,tmp.s);
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
                fprintf(bcftools_stderr,"Warning: using %s filter name instead of \"%s\"\n", tmp.s,args->soft_filter);
            free(tmp.s);
        }

        rbuf_init(&args->rbuf, 64);
        args->rbuf_lines = (bcf1_t**) calloc(args->rbuf.m, sizeof(bcf1_t*));
        if ( args->snp_gap )
        {
            bcf_hdr_printf(args->hdr, "##FILTER=<ID=SnpGap,Description=\"SNP within %d bp of %s\">", args->snp_gap,args->snp_gap_str);
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
    if ( args->mask ) regidx_destroy(args->mask);
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
        if ( pass && bcf_write1(args->out_fh, args->hdr, rec)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
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
    const int SnpGap_set     = 1 << (8*sizeof(int)/2);
    const int IndelGap_set   = 1 << (8*sizeof(int)/2-1);
    const int IndelGap_flush = 1 << (8*sizeof(int)/2-2);

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
        // Set the REF allele's length to max deletion length or to 1 if a SNP or an insertion.
        line->d.var[0].n = line->rlen;
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
            int k = 0, max_ac = -1, imax_ac = -1, max_qual = -1, imax_qual = -1;
            for (i=-1; rbuf_next(&args->rbuf,&i) && k<k_flush; )
            {
                k++;
                bcf1_t *rec  = args->rbuf_lines[i];
                if ( !(rec->d.var_type & IndelGap_set) ) continue;
                hts_expand(int, rec->n_allele, args->ntmpi, args->tmpi);
                int ret = bcf_calc_ac(args->hdr, rec, args->tmpi, BCF_UN_ALL);
                if ( imax_ac==-1 || (ret && max_ac < args->tmpi[1]) ) { max_ac = args->tmpi[1]; imax_ac = i; }
                if ( imax_qual==-1 || max_qual < rec->qual ) { max_qual = rec->qual; imax_qual = i; }
            }

            // Filter all but the best indel (with the best QUAL, bigger AC, or take the first if neither QUAL nor AC are available)
            k = 0;
            for (i=-1; rbuf_next(&args->rbuf,&i) && k<k_flush; )
            {
                k++;
                bcf1_t *rec = args->rbuf_lines[i];
                if ( !(rec->d.var_type & IndelGap_set) ) continue;
                rec->d.var_type |= IndelGap_flush;

                int do_filter = 0;
                if ( max_qual>0 )
                {
                    if ( i!=imax_qual ) do_filter = 1;
                }
                else if ( i!=imax_ac ) do_filter = 1;
                if ( do_filter ) bcf_add_filter(args->hdr, args->rbuf_lines[i], args->IndelGap_id);
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
            else if ( (var_type & args->snp_gap_type) && (rec->d.var_type & VCF_SNP) && !(rec->d.var_type & SnpGap_set) )
            {
                // this SNP has not been SnpGap-filtered yet
                rec->d.var_type |= SnpGap_set;
                bcf_add_filter(args->hdr, rec, args->SnpGap_id);
            }
            else if ( (var_type & VCF_SNP) && (rec->d.var_type & args->snp_gap_type) )
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

static void _set_variant_boundaries(bcf1_t *rec, hts_pos_t *beg, hts_pos_t *end)
{
    hts_pos_t off;
    if ( rec->n_allele )
    {
        off = rec->rlen;
        bcf_unpack(rec, BCF_UN_STR);
        int i;
        for (i=1; i<rec->n_allele; i++)
        {
            // Make symbolic alleles start at POS, although this is not strictly true for
            // <DEL>,<INS> where POS should be the position BEFORE the deletion/insertion.
            // However, since arbitrary symbolic alleles can be defined by the user, we
            // will simplify the interpretation of --targets-overlap and --region-overlap.
            int j = 0;
            char *ref = rec->d.allele[0];
            char *alt = rec->d.allele[i];
            while ( ref[j] && alt[j] && ref[j]==alt[j] ) j++;
            if ( off > j ) off = j;
            if ( !off ) break;
        }
    }
    else
        off = 0;

    *beg = rec->pos + off;
    *end = rec->pos + rec->rlen - 1;
}

static void usage(args_t *args)
{
    fprintf(bcftools_stderr, "\n");
    fprintf(bcftools_stderr, "About:   Apply fixed-threshold filters.\n");
    fprintf(bcftools_stderr, "Usage:   bcftools filter [options] <in.vcf.gz>\n");
    fprintf(bcftools_stderr, "\n");
    fprintf(bcftools_stderr, "Options:\n");
    fprintf(bcftools_stderr, "    -e, --exclude EXPR             Exclude sites for which the expression is true (see man page for details)\n");
    fprintf(bcftools_stderr, "    -g, --SnpGap INT[:TYPE]        Filter SNPs within <int> base pairs of an indel (the default) or any combination of indel,mnp,bnd,other,overlap\n");
    fprintf(bcftools_stderr, "    -G, --IndelGap INT             Filter clusters of indels separated by <int> or fewer base pairs allowing only one to pass\n");
    fprintf(bcftools_stderr, "    -i, --include EXPR             Include only sites for which the expression is true (see man page for details\n");
    fprintf(bcftools_stderr, "        --mask [^]REGION           Soft filter regions, \"^\" to negate\n");
    fprintf(bcftools_stderr, "    -M, --mask-file [^]FILE        Soft filter regions listed in a file, \"^\" to negate\n");
    fprintf(bcftools_stderr, "        --mask-overlap 0|1|2       Mask if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n");
    fprintf(bcftools_stderr, "    -m, --mode [+x]                \"+\": do not replace but add to existing FILTER; \"x\": reset filters at sites which pass\n");
    fprintf(bcftools_stderr, "        --no-version               Do not append version and command line to the header\n");
    fprintf(bcftools_stderr, "    -o, --output FILE              Write output to a file [standard output]\n");
    fprintf(bcftools_stderr, "    -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n");
    fprintf(bcftools_stderr, "    -r, --regions REGION           Restrict to comma-separated list of regions\n");
    fprintf(bcftools_stderr, "    -R, --regions-file FILE        Restrict to regions listed in a file\n");
    fprintf(bcftools_stderr, "        --regions-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n");
    fprintf(bcftools_stderr, "    -s, --soft-filter STRING       Annotate FILTER column with <string> or unique filter name (\"Filter%%d\") made up by the program (\"+\")\n");
    fprintf(bcftools_stderr, "    -S, --set-GTs .|0              Set genotypes of failed samples to missing (.) or ref (0)\n");
    fprintf(bcftools_stderr, "    -t, --targets REGION           Similar to -r but streams rather than index-jumps\n");
    fprintf(bcftools_stderr, "    -T, --targets-file FILE        Similar to -R but streams rather than index-jumps\n");
    fprintf(bcftools_stderr, "        --targets-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n");
    fprintf(bcftools_stderr, "        --threads INT              Use multithreading with <int> worker threads [0]\n");
    fprintf(bcftools_stderr, "\n");
    bcftools_exit(1);
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
    args->clevel = -1;
    int regions_is_file = 0, targets_is_file = 0;
    int regions_overlap = 1;
    int targets_overlap = 0;
    args->mask_overlap  = 1;

    static struct option loptions[] =
    {
        {"set-GTs",required_argument,NULL,'S'},
        {"mode",required_argument,NULL,'m'},
        {"mask",required_argument,NULL,10},
        {"mask-file",required_argument,NULL,'M'},
        {"mask-overlap",required_argument,NULL,11},
        {"soft-filter",required_argument,NULL,'s'},
        {"exclude",required_argument,NULL,'e'},
        {"include",required_argument,NULL,'i'},
        {"targets",required_argument,NULL,'t'},
        {"targets-file",required_argument,NULL,'T'},
        {"targets-overlap",required_argument,NULL,4},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"regions-overlap",required_argument,NULL,3},
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
                if ( *tmp && *tmp!=':' ) error("Could not parse argument: --SnpGap %s\n", optarg);
                if ( *tmp==':' )
                {
                    args->snp_gap_str = tmp+1;
                    int i,n;
                    char **keys = hts_readlist(tmp+1,0,&n);
                    for(i=0; i<n; i++)
                    {
                        if ( !strcasecmp(keys[i],"indel") ) args->snp_gap_type |= VCF_INDEL;
                        else if ( !strcasecmp(keys[i],"mnp") ) args->snp_gap_type |= VCF_MNP;
                        else if ( !strcasecmp(keys[i],"bnd") ) args->snp_gap_type |= VCF_BND;
                        else if ( !strcasecmp(keys[i],"other") ) args->snp_gap_type |= VCF_OTHER;
                        else if ( !strcasecmp(keys[i],"overlap") ) args->snp_gap_type |= VCF_OVERLAP;
                        else error("Could not parse \"%s\" in \"--SnpGap %s\"\n", keys[i], optarg);
                        free(keys[i]);
                    }
                    if ( n ) free(keys);
                }
                else
                {
                    args->snp_gap_type = VCF_INDEL;
                    args->snp_gap_str = "indel";
                }
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
                    default:
                    {
                        args->clevel = strtol(optarg,&tmp,10);
                        if ( *tmp || args->clevel<0 || args->clevel>9 ) error("The output type \"%s\" not recognised\n", optarg);
                    }
                }
                if ( optarg[1] )
                {
                    args->clevel = strtol(optarg+1,&tmp,10);
                    if ( *tmp || args->clevel<0 || args->clevel>9 ) error("Could not parse argument: --compression-level %s\n", optarg+1);
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
            case 'e':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'S':
                if ( !strcmp(".",optarg) ) args->set_gts = SET_GTS_MISSING;
                else if ( !strcmp("0",optarg) ) args->set_gts = SET_GTS_REF;
                else error("The argument to -S not recognised: %s\n", optarg);
                break;
            case  9 : args->n_threads = strtol(optarg, 0, 0); break;
            case  8 : args->record_cmd_line = 0; break;
            case  3 :
                regions_overlap = parse_overlap_option(optarg);
                if ( regions_overlap < 0 ) error("Could not parse: --regions-overlap %s\n",optarg);
                break;
            case  4 :
                targets_overlap = parse_overlap_option(optarg);
                if ( targets_overlap < 0 ) error("Could not parse: --targets-overlap %s\n",optarg);
                break;
            case  10 : args->mask_list = optarg; break;
            case 'M' : args->mask_list = optarg; args->mask_is_file = 1; break;
            case  11 :
                if ( !strcasecmp(optarg,"0") ) args->mask_overlap = 0;
                else if ( !strcasecmp(optarg,"1") ) args->mask_overlap = 1;
                else if ( !strcasecmp(optarg,"2") ) args->mask_overlap = 2;
                else error("Could not parse: --mask-overlap %s\n",optarg);
                break;
            case 'h':
            case '?': usage(args); break;
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

    if ( args->mask_list && !args->soft_filter ) error("The option --soft-filter is required with --mask and --mask-file options\n");

    // read in the regions from the command line
    if ( args->regions_list )
    {
        args->files->require_index = 1;
        bcf_sr_set_opt(args->files,BCF_SR_REGIONS_OVERLAP,regions_overlap);
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
        bcf_sr_set_opt(args->files,BCF_SR_REGIONS_OVERLAP,regions_overlap);
        if ( bcf_sr_set_regions(args->files, tmp.s, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
        free(tmp.s);
    }
    if ( args->targets_list )
    {
        bcf_sr_set_opt(args->files,BCF_SR_TARGETS_OVERLAP,targets_overlap);
        if ( bcf_sr_set_targets(args->files, args->targets_list,targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
    }
    if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to read from %s: %s\n", !strcmp("-",fname)?"standard input":fname,bcf_sr_strerror(args->files->errnum));

    init_data(args);
    if ( bcf_hdr_write(args->out_fh, args->hdr)!=0 ) error("[%s] Error: cannot write the header to %s\n", __func__,args->output_fname);
    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = bcf_sr_get_line(args->files, 0);
        int pass = 1;
        if ( args->filter )
        {
            pass = filter_test(args->filter, line, &args->smpl_pass);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
        }
        if ( args->mask )
        {
            hts_pos_t beg, end;
            if ( args->mask_overlap==0 ) beg = end = line->pos;
            else if ( args->mask_overlap==1 ) beg = line->pos, end = line->pos + line->rlen - 1;
            else _set_variant_boundaries(line,&beg,&end);
            int mpass = regidx_overlap(args->mask,bcf_seqname(args->hdr,line),beg,end,NULL) ? 0 : 1;
            if ( args->mask_negate ) mpass = mpass ? 0 : 1;
            pass &= mpass;
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
            {
                if ( bcf_write1(args->out_fh, args->hdr, line)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
            }
            else
                buffered_filters(args, line);
        }
    }
    buffered_filters(args, NULL);

    if ( hts_close(args->out_fh)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->output_fname);
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}
