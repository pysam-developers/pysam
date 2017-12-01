#include "pysam.h"

/*  vcfconvert.c -- convert between VCF/BCF and related formats.

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
#include <strings.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <htslib/kseq.h>
#include "bcftools.h"
#include "filter.h"
#include "convert.h"
#include "tsv2vcf.h"

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef struct _args_t args_t;
struct _args_t
{
    faidx_t *ref;
    filter_t *filter;
    char *filter_str;
    int filter_logic;   // include or exclude sites which match the filters? One of FLT_INCLUDE/FLT_EXCLUDE
    convert_t *convert;
    bcf_srs_t *files;
    bcf_hdr_t *header;
    void (*convert_func)(struct _args_t *);
    struct {
        int total, skipped, hom_rr, het_ra, hom_aa, het_aa, missing; 
    } n;
    kstring_t str;
    int32_t *gts;
    float *flt;
    int rev_als, output_vcf_ids, hap2dip, output_chrom_first_col;
    int nsamples, *samples, sample_is_file, targets_is_file, regions_is_file, output_type;
    char **argv, *sample_list, *targets_list, *regions_list, *tag, *columns;
    char *outfname, *infname, *ref_fname, *sex_fname;
    int argc, n_threads, record_cmd_line;
};

static void destroy_data(args_t *args)
{
    if ( args->ref ) fai_destroy(args->ref);
    if ( args->convert) convert_destroy(args->convert);
    if ( args->filter ) filter_destroy(args->filter);
    free(args->samples);
    if ( args->files ) bcf_sr_destroy(args->files);
}

static void open_vcf(args_t *args, const char *format_str)
{
    args->files = bcf_sr_init();
    if ( args->n_threads && bcf_sr_set_threads(args->files, args->n_threads)!=0 )
        error("Could not initialize --threads %d\n", args->n_threads);

    if ( args->regions_list )
    {
        if ( bcf_sr_set_regions(args->files, args->regions_list, args->regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    if ( args->targets_list )
    {
        if ( bcf_sr_set_targets(args->files, args->targets_list, args->targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
    }
    if ( !bcf_sr_add_reader(args->files, args->infname) )
        error("Failed to open %s: %s\n", args->infname,bcf_sr_strerror(args->files->errnum));

    args->header = args->files->readers[0].header;

    if ( args->filter_str )
        args->filter = filter_init(args->header, args->filter_str);

    int i, nsamples = 0, *samples = NULL;
    if ( args->sample_list && strcmp("-",args->sample_list) )
    {
        for (i=0; i<args->files->nreaders; i++)
        {
            int ret = bcf_hdr_set_samples(args->files->readers[i].header,args->sample_list,args->sample_is_file);
            if ( ret<0 ) error("Error parsing the sample list\n");
            else if ( ret>0 ) error("Sample name mismatch: sample #%d not found in the header\n", ret);
        }

        if ( args->sample_list[0]!='^' )
        {
            // the sample ordering may be different if not negated
            int n;
            char **smpls = hts_readlist(args->sample_list, args->sample_is_file, &n);
            if ( !smpls ) error("Could not parse %s\n", args->sample_list);
            if ( n!=bcf_hdr_nsamples(args->files->readers[0].header) )
                error("The number of samples does not match, perhaps some are present multiple times?\n");
            nsamples = bcf_hdr_nsamples(args->files->readers[0].header);
            samples = (int*) malloc(sizeof(int)*nsamples);
            for (i=0; i<n; i++)
            {
                samples[i] = bcf_hdr_id2int(args->files->readers[0].header, BCF_DT_SAMPLE,smpls[i]);
                free(smpls[i]);
            }
            free(smpls);
        }
    }
    if ( format_str ) args->convert = convert_init(args->header, samples, nsamples, format_str);
    free(samples);
}

static int tsv_setter_chrom_pos_ref_alt(tsv_t *tsv, bcf1_t *rec, void *usr)
{
    args_t *args = (args_t*) usr;

    char tmp, *se = tsv->ss, *ss = tsv->ss;
    while ( se < tsv->se && *se!=':' ) se++;
    if ( *se!=':' ) error("Could not parse CHROM in CHROM:POS_REF_ALT id: %s\n", tsv->ss);
    tmp = *se; *se = 0;
    rec->rid = bcf_hdr_name2id(args->header,ss); 
    if ( rec->rid<0 ) error("Could not determine sequence name or multiple sequences present: %s\n", tsv->ss);
    *se = tmp;

    // POS
    rec->pos = strtol(se+1,&ss,10);
    if ( ss==se+1 ) error("Could not parse POS in CHROM:POS_REF_ALT: %s\n", tsv->ss);
    rec->pos--;

    // REF,ALT
    args->str.l = 0;
    se = ++ss;
    while ( se < tsv->se && *se!='_' ) se++; 
    if ( *se!='_' ) error("Could not parse REF in CHROM:POS_REF_ALT id: %s\n", tsv->ss);
    kputsn(ss,se-ss,&args->str);
    ss = ++se;
    while ( se < tsv->se && *se!='_' && isspace(*tsv->se) ) se++;
    if ( se < tsv->se && *se!='_' && isspace(*tsv->se) ) error("Could not parse ALT in CHROM:POS_REF_ALT id: %s\n", tsv->ss);
    kputc(',',&args->str);
    kputsn(ss,se-ss,&args->str);
    bcf_update_alleles_str(args->header, rec, args->str.s);

    // END - optional
    if (*se && *se=='_') {
        long end = strtol(se+1,&ss,10);
        if ( ss==se+1 ) error("Could not parse END in CHROM:POS_REF_ALT_END: %s\n", tsv->ss);
        bcf_update_info_int32(args->header, rec, "END", &end, 1);
    }

    return 0;
}
static int tsv_setter_verify_pos(tsv_t *tsv, bcf1_t *rec, void *usr)
{
    char *se;
    int pos = strtol(tsv->ss,&se,10);
    if ( tsv->ss==se ) error("Could not parse POS: %s\n", tsv->ss);
    if ( rec->pos != pos-1 ) error("POS mismatch: %s\n", tsv->ss);
    return 0;
}
static int tsv_setter_verify_ref_alt(tsv_t *tsv, bcf1_t *rec, void *usr)
{
    args_t *args = (args_t*) usr;
    args->rev_als = 0;
    char tmp = *tsv->se; *tsv->se = 0;
    if ( strcmp(tsv->ss,rec->d.allele[0]) )
    {
        if ( strcmp(tsv->ss,rec->d.allele[1]) ) { *tsv->se = tmp; error("REF/ALT mismatch: [%s][%s]\n", tsv->ss,rec->d.allele[1]); }
        args->rev_als = 1;
    }
    *tsv->se = tmp;
    while ( *tsv->se && isspace(*tsv->se) ) tsv->se++;
    tsv->ss = tsv->se;
    while ( *tsv->se && !isspace(*tsv->se) ) tsv->se++;
    tmp = *tsv->se; *tsv->se = 0;
    if ( !args->rev_als && strcmp(tsv->ss,rec->d.allele[1]) ) { *tsv->se = tmp; error("REF/ALT mismatch: [%s][%s]\n", tsv->ss,rec->d.allele[1]); }
    else if ( args->rev_als && strcmp(tsv->ss,rec->d.allele[0]) ) { *tsv->se = tmp; error("REF/ALT mismatch: [%s][%s]\n", tsv->ss,rec->d.allele[0]); }
    *tsv->se = tmp;
    return 0;
}
static int tsv_setter_gt_gp(tsv_t *tsv, bcf1_t *rec, void *usr)
{
    args_t *args = (args_t*) usr;
    int i, nsamples = bcf_hdr_nsamples(args->header);
    for (i=0; i<nsamples; i++)
    {
        float aa,ab,bb;
        aa = strtod(tsv->ss, &tsv->se);
        if ( tsv->ss==tsv->se ) { fprintf(pysam_stderr,"Could not parse first value of %d-th sample\n", i+1); return -1; }
        tsv->ss = tsv->se+1;
        ab = strtod(tsv->ss, &tsv->se);
        if ( tsv->ss==tsv->se ) { fprintf(pysam_stderr,"Could not parse second value of %d-th sample\n", i+1); return -1; }
        tsv->ss = tsv->se+1;
        bb = strtod(tsv->ss, &tsv->se);
        if ( tsv->ss==tsv->se ) { fprintf(pysam_stderr,"Could not parse third value of %d-th sample\n", i+1); return -1; }
        tsv->ss = tsv->se+1;

        if ( args->rev_als ) { float tmp = bb; bb = aa; aa = tmp; }
        args->flt[3*i+0] = aa;
        args->flt[3*i+1] = ab;
        args->flt[3*i+2] = bb;

        if ( aa >= ab )
        {
            if ( aa >= bb ) args->gts[2*i+0] = args->gts[2*i+1] = bcf_gt_unphased(0);
            else args->gts[2*i+0] = args->gts[2*i+1] = bcf_gt_unphased(1); 
        }
        else if ( ab >= bb ) 
        {
            args->gts[2*i+0] = bcf_gt_unphased(0);
            args->gts[2*i+1] = bcf_gt_unphased(1); 
        }
        else args->gts[2*i+0] = args->gts[2*i+1] = bcf_gt_unphased(1);
    }
    if ( *tsv->se ) error("Could not parse: %s\n", tsv->ss);
    if ( bcf_update_genotypes(args->header,rec,args->gts,nsamples*2) ) error("Could not update GT field\n");
    if ( bcf_update_format_float(args->header,rec,"GP",args->flt,nsamples*3) ) error("Could not update GP field\n");
    return 0;
}
static int tsv_setter_haps(tsv_t *tsv, bcf1_t *rec, void *usr)
{
    args_t *args = (args_t*) usr;
    int i, nsamples = bcf_hdr_nsamples(args->header);

    int32_t a0, a1;
    if ( args->rev_als ) { a0 = bcf_gt_phased(1); a1 = bcf_gt_phased(0); }
    else { a0 = bcf_gt_phased(0); a1 = bcf_gt_phased(1); }

    // up is short for "unphased"
    int nup = 0; 
    for (i=0; i<nsamples; i++)
    {
        char *ss = tsv->ss + 4*i + nup;
        int up = 0, all;

        for (all=0; all < 2; all++){
            // checking for premature ending
            if ( !ss[0] || !ss[1] || !ss[2] ||
                 (up && (!ss[3] || !ss[4]) ) )
            {
                fprintf(pysam_stderr,"Wrong number of fields at %d-th sample ([%c][%c][%c]). ",i+1,ss[0],ss[1],ss[2]);
                return -1;
            }

            switch(ss[all*2+up]){
            case '0':
                args->gts[2*i+all] = a0;
                break;
            case '1' :
                args->gts[2*i+all] = a1;
                break;
            case '?' :
                // there is no macro to express phased missing allele
                args->gts[2*i+all] = bcf_gt_phased(-1);
                break;
            case '-' :
                args->gts[2*i+all] = bcf_int32_vector_end;
                break;
            default :
                fprintf(pysam_stderr,"Could not parse: [%c][%s]\n", ss[all*2+up],tsv->ss);
                return -1; 
            }
            if( ss[all*2+up+1]=='*' ) up = up + 1;
        }
        
        if(up && up != 2)
        {
            fprintf(pysam_stderr,"Missing unphased marker '*': [%c][%s]", ss[2+up], tsv->ss);
            return -1;
        }

        // change alleles to unphased if the alleles are unphased
        if ( up )
        {
            args->gts[2*i] = bcf_gt_unphased(bcf_gt_allele(args->gts[2*i]));
            args->gts[2*i+1] = bcf_gt_unphased(bcf_gt_allele(args->gts[2*i+1]));
        }
        nup = nup + up;
    }
    if ( tsv->ss[(nsamples-1)*4+3+nup] )
    {
        fprintf(pysam_stderr,"nup: %d", nup);
        fprintf(pysam_stderr,"Wrong number of fields (%d-th column = [%c]). ", nsamples*2,tsv->ss[(nsamples-1)*4+nup]);
        return -1;
    }

    if ( bcf_update_genotypes(args->header,rec,args->gts,nsamples*2) ) error("Could not update GT field\n");
    return 0;
}
static void gensample_to_vcf(args_t *args)
{
    /*
     *  Inpute: IMPUTE2 output (indentation changed here for clarity): 
     *
     *      20:62116619_C_T 20:62116619     62116619 C T 0.969 0.031 0 ...
     *      ---             20:62116698_C_A 62116698 C A 1     0     0 ...
     *
     *  Second column is expected in the form of CHROM:POS_REF_ALT. We use second
     *  column because the first can be empty ("--") when filling sites from reference 
     *  panel.
     *
     *  Output: VCF with filled GT,GP
     *
     */
    kstring_t line = {0,0,0};

    char *gen_fname = NULL, *sample_fname = NULL;
    sample_fname = strchr(args->infname,',');
    if ( !sample_fname )
    {
        args->str.l = 0;
        ksprintf(&args->str,"%s.gen.gz", args->infname);
        gen_fname = strdup(args->str.s);
        args->str.l = 0;
        ksprintf(&args->str,"%s.samples", args->infname);
        sample_fname = strdup(args->str.s);
    }
    else
    {
        *sample_fname = 0;
        gen_fname = strdup(args->infname);
        sample_fname = strdup(sample_fname+1);
    }
    htsFile *gen_fh = hts_open(gen_fname, "r");
    if ( !gen_fh ) error("Could not read: %s\n", gen_fname);
    if ( hts_getline(gen_fh, KS_SEP_LINE, &line) <= 0 ) error("Empty file: %s\n", gen_fname);

    // Find out the chromosome name, sample names, init and print the VCF header
    args->str.l = 0;
    char *ss, *se = line.s;
    while ( *se && !isspace(*se) ) se++;
    if ( !*se ) error("Could not parse %s: %s\n", gen_fname,line.s);
    ss = se+1;
    se = strchr(ss,':');
    if ( !se ) error("Expected CHROM:POS_REF_ALT in second column of %s\n", gen_fname);
    kputsn(ss, se-ss, &args->str);

    tsv_t *tsv = tsv_init("-,CHROM_POS_REF_ALT,POS,REF_ALT,GT_GP");
    tsv_register(tsv, "CHROM_POS_REF_ALT", tsv_setter_chrom_pos_ref_alt, args);
    tsv_register(tsv, "POS", tsv_setter_verify_pos, NULL);
    tsv_register(tsv, "REF_ALT", tsv_setter_verify_ref_alt, args);
    tsv_register(tsv, "GT_GP", tsv_setter_gt_gp, args);

    args->header = bcf_hdr_init("w");
    bcf_hdr_append(args->header, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">");
    bcf_hdr_append(args->header, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(args->header, "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype Probabilities\">");
    bcf_hdr_printf(args->header, "##contig=<ID=%s,length=%d>", args->str.s,0x7fffffff);   // MAX_CSI_COOR
    if (args->record_cmd_line) bcf_hdr_append_version(args->header, args->argc, args->argv, "bcftools_convert");

    int i, nsamples;
    char **samples = hts_readlist(sample_fname, 1, &nsamples);
    if ( !samples ) error("Could not read %s\n", sample_fname);
    for (i=2; i<nsamples; i++)
    {
        se = samples[i]; while ( *se && !isspace(*se) ) se++;
        *se = 0;
        bcf_hdr_add_sample(args->header,samples[i]);
    }
    for (i=0; i<nsamples; i++) free(samples[i]);
    free(samples);

    htsFile *out_fh = hts_open(args->outfname,hts_bcf_wmode(args->output_type));
    if ( out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->outfname, strerror(errno));
    if ( args->n_threads ) hts_set_threads(out_fh, args->n_threads);
    bcf_hdr_write(out_fh,args->header);
    bcf1_t *rec = bcf_init();

    nsamples -= 2;
    args->gts = (int32_t *) malloc(sizeof(int32_t)*nsamples*2);
    args->flt = (float *) malloc(sizeof(float)*nsamples*3);

    do
    {
        bcf_clear(rec);
        args->n.total++;
        if ( !tsv_parse(tsv, rec, line.s) )
            bcf_write(out_fh, args->header, rec);
        else
            error("Error occurred while parsing: %s\n", line.s);
    }
    while ( hts_getline(gen_fh, KS_SEP_LINE, &line)>0 );

    if ( hts_close(out_fh) ) error("Close failed: %s\n", args->outfname);
    if ( hts_close(gen_fh) ) error("Close failed: %s\n", gen_fname);
    bcf_hdr_destroy(args->header);
    bcf_destroy(rec);
    free(sample_fname);
    free(gen_fname);
    free(args->str.s);
    free(line.s);
    free(args->gts);
    free(args->flt);
    tsv_destroy(tsv);

    fprintf(pysam_stderr,"Number of processed rows: \t%d\n", args->n.total);
}

static void haplegendsample_to_vcf(args_t *args)
{
    /*
     *  Convert from IMPUTE2 hap/legend/sample output files to VCF
     *
     *      hap:
     *          0 1 0 1
     *      legend:
     *          id position a0 a1
     *          1:186946386_G_T 186946386 G T
     *      sample:
     *          sample population group sex
     *          sample1 sample1 sample1 2
     *          sample2 sample2 sample2 2
     *
     *  Output: VCF with filled GT
     */
    kstring_t line = {0,0,0};

    char *hap_fname = NULL, *leg_fname = NULL, *sample_fname = NULL;
    sample_fname = strchr(args->infname,',');
    if ( !sample_fname )
    {
        args->str.l = 0;
        ksprintf(&args->str,"%s.hap.gz", args->infname);
        hap_fname = strdup(args->str.s);
        args->str.l = 0;
        ksprintf(&args->str,"%s.samples", args->infname);
        sample_fname = strdup(args->str.s);
        args->str.l = 0;
        ksprintf(&args->str,"%s.legend.gz", args->infname);
        leg_fname = strdup(args->str.s);
    }
    else
    {
        char *ss = sample_fname, *se = strchr(ss+1,',');
        if ( !se ) error("Could not parse hap/legend/sample file names: %s\n", args->infname);
        *ss = 0;
        *se = 0;
        hap_fname = strdup(args->infname);
        leg_fname = strdup(ss+1);
        sample_fname = strdup(se+1);
    }
    htsFile *hap_fh = hts_open(hap_fname, "r");
    if ( !hap_fh ) error("Could not read: %s\n", hap_fname);

    htsFile *leg_fh = hts_open(leg_fname,"r");
    if ( !leg_fh ) error("Could not read: %s\n", leg_fname);

    // Eat up first legend line, then determine chromosome name
    if ( hts_getline(leg_fh, KS_SEP_LINE, &line) <= 0 ) error("Empty file: %s\n", leg_fname);
    if ( hts_getline(leg_fh, KS_SEP_LINE, &line) <= 0 ) error("Empty file: %s\n", leg_fname);

    // Find out the chromosome name, sample names, init and print the VCF header
    args->str.l = 0;
    char *se = strchr(line.s,':');
    if ( !se ) error("Expected CHROM:POS_REF_ALT in first column of %s\n", leg_fname);
    kputsn(line.s, se-line.s, &args->str);

    tsv_t *leg_tsv = tsv_init("CHROM_POS_REF_ALT,POS,REF_ALT");
    tsv_register(leg_tsv, "CHROM_POS_REF_ALT", tsv_setter_chrom_pos_ref_alt, args);
    tsv_register(leg_tsv, "POS", tsv_setter_verify_pos, NULL);
    tsv_register(leg_tsv, "REF_ALT", tsv_setter_verify_ref_alt, args);

    tsv_t *hap_tsv = tsv_init("HAPS");
    tsv_register(hap_tsv, "HAPS", tsv_setter_haps, args);

    args->header = bcf_hdr_init("w");
    bcf_hdr_append(args->header, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">");
    bcf_hdr_append(args->header, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_printf(args->header, "##contig=<ID=%s,length=%d>", args->str.s,0x7fffffff);   // MAX_CSI_COOR
    if (args->record_cmd_line) bcf_hdr_append_version(args->header, args->argc, args->argv, "bcftools_convert");

    int i, nrows, nsamples;
    char **samples = hts_readlist(sample_fname, 1, &nrows);
    if ( !samples ) error("Could not read %s\n", sample_fname);
    nsamples = nrows - 1;

    // sample_fname should contain a header line, so need to ignore first row
    // returned from hts_readlist (i=1, and not i=0)
    for (i=1; i<nrows; i++)
    {
        se = samples[i]; while ( *se && !isspace(*se) ) se++;
        *se = 0;
        bcf_hdr_add_sample(args->header,samples[i]);
    }
    bcf_hdr_add_sample(args->header,NULL);
    for (i=0; i<nrows; i++) free(samples[i]);
    free(samples);

    htsFile *out_fh = hts_open(args->outfname,hts_bcf_wmode(args->output_type));
    if ( out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->outfname, strerror(errno));
    if ( args->n_threads ) hts_set_threads(out_fh, args->n_threads);
    bcf_hdr_write(out_fh,args->header);
    bcf1_t *rec = bcf_init();

    args->gts = (int32_t *) malloc(sizeof(int32_t)*nsamples*2);

    while (1)
    {
        bcf_clear(rec);
        args->n.total++;
        if ( tsv_parse(leg_tsv, rec, line.s) )
            error("Error occurred while parsing %s: %s\n", leg_fname,line.s);

        if ( hts_getline(hap_fh,  KS_SEP_LINE, &line)<=0 )
            error("Different number of records in %s and %s?\n", leg_fname,hap_fname);

        if ( tsv_parse(hap_tsv, rec, line.s) )
            error("Error occurred while parsing %s: %s\n", hap_fname,line.s);

        bcf_write(out_fh, args->header, rec);

        if ( hts_getline(leg_fh, KS_SEP_LINE, &line)<=0 )
        {
            if ( hts_getline(hap_fh, KS_SEP_LINE, &line)>0 )
                error("Different number of records in %s and %s?\n", leg_fname,hap_fname);
            break;
        }
    }

    if ( hts_close(out_fh) ) error("Close failed: %s\n", args->outfname);
    if ( hts_close(hap_fh) ) error("Close failed: %s\n", hap_fname);
    if ( hts_close(leg_fh) ) error("Close failed: %s\n", leg_fname);
    bcf_hdr_destroy(args->header);
    bcf_destroy(rec);
    free(sample_fname);
    free(hap_fname);
    free(leg_fname);
    free(args->str.s);
    free(line.s);
    free(args->gts);
    tsv_destroy(hap_tsv);
    tsv_destroy(leg_tsv);

    fprintf(pysam_stderr,"Number of processed rows: \t%d\n", args->n.total);
}

static void hapsample_to_vcf(args_t *args)
{
    /*
     *  Input: SHAPEIT output
     *
     *      20:19995888_A_G 20:19995888 19995888 A G 0 0 0 0 ...
     *
     *  First column is expected in the form of CHROM:POS_REF_ALT
     *
     *  Output: VCF with filled GT
     *
     */
    kstring_t line = {0,0,0};

    char *hap_fname = NULL, *sample_fname = NULL;
    sample_fname = strchr(args->infname,',');
    if ( !sample_fname )
    {
        args->str.l = 0;
        ksprintf(&args->str,"%s.hap.gz", args->infname);
        hap_fname = strdup(args->str.s);
        args->str.l = 0;
        ksprintf(&args->str,"%s.samples", args->infname);
        sample_fname = strdup(args->str.s);
    }
    else
    {
        *sample_fname = 0;
        hap_fname = strdup(args->infname);
        sample_fname = strdup(sample_fname+1);
    }
    htsFile *hap_fh = hts_open(hap_fname, "r");
    if ( !hap_fh ) error("Could not read: %s\n", hap_fname);
    if ( hts_getline(hap_fh, KS_SEP_LINE, &line) <= 0 ) error("Empty file: %s\n", hap_fname);

    // Find out the chromosome name, sample names, init and print the VCF header
    args->str.l = 0;
    char *se = strchr(line.s,':');
    if ( !se ) error("Expected CHROM:POS_REF_ALT in first column of %s\n", hap_fname);
    kputsn(line.s, se-line.s, &args->str);

    tsv_t *tsv = tsv_init("CHROM_POS_REF_ALT,-,POS,REF_ALT,HAPS");
    tsv_register(tsv, "CHROM_POS_REF_ALT", tsv_setter_chrom_pos_ref_alt, args);
    tsv_register(tsv, "POS", tsv_setter_verify_pos, NULL);
    tsv_register(tsv, "REF_ALT", tsv_setter_verify_ref_alt, args);
    tsv_register(tsv, "HAPS", tsv_setter_haps, args);

    args->header = bcf_hdr_init("w");
    bcf_hdr_append(args->header, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">");
    bcf_hdr_append(args->header, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_printf(args->header, "##contig=<ID=%s,length=%d>", args->str.s,0x7fffffff);   // MAX_CSI_COOR
    if (args->record_cmd_line) bcf_hdr_append_version(args->header, args->argc, args->argv, "bcftools_convert");

    int i, nsamples;
    char **samples = hts_readlist(sample_fname, 1, &nsamples);
    if ( !samples ) error("Could not read %s\n", sample_fname);
    for (i=2; i<nsamples; i++)
    {
        se = samples[i]; while ( *se && !isspace(*se) ) se++;
        *se = 0;
        bcf_hdr_add_sample(args->header,samples[i]);
    }
    bcf_hdr_add_sample(args->header,NULL);
    for (i=0; i<nsamples; i++) free(samples[i]);
    free(samples);

    htsFile *out_fh = hts_open(args->outfname,hts_bcf_wmode(args->output_type));
    if ( out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->outfname, strerror(errno));
    if ( args->n_threads ) hts_set_threads(out_fh, args->n_threads);
    bcf_hdr_write(out_fh,args->header);
    bcf1_t *rec = bcf_init();

    nsamples -= 2;
    args->gts = (int32_t *) malloc(sizeof(int32_t)*nsamples*2);

    do
    {
        bcf_clear(rec);
        args->n.total++;
        if ( !tsv_parse(tsv, rec, line.s) )
            bcf_write(out_fh, args->header, rec);
        else
            error("Error occurred while parsing: %s\n", line.s);
    }
    while ( hts_getline(hap_fh, KS_SEP_LINE, &line)>0 );

    if ( hts_close(out_fh) ) error("Close failed: %s\n", args->outfname);
    if ( hts_close(hap_fh) ) error("Close failed: %s\n", hap_fname);
    bcf_hdr_destroy(args->header);
    bcf_destroy(rec);
    free(sample_fname);
    free(hap_fname);
    free(args->str.s);
    free(line.s);
    free(args->gts);
    tsv_destroy(tsv);

    fprintf(pysam_stderr,"Number of processed rows: \t%d\n", args->n.total);
}

char *init_sample2sex(bcf_hdr_t *hdr, char *sex_fname)
{
    int i, nlines;
    char *sample2sex = (char*) calloc(bcf_hdr_nsamples(hdr),1);
    char **lines = hts_readlist(sex_fname, 1, &nlines);
    if ( !lines ) error("Could not read %s\n", sex_fname);
    for (i=0; i<nlines; i++)
    {
        char *se = lines[i]; while ( *se && !isspace(*se) ) se++;
        char tmp = *se;
        *se = 0;
        int id = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, lines[i]);
        *se = tmp;
        if ( id<0 ) continue;
        while ( *se && isspace(*se) ) se++;
        if ( *se=='M' ) sample2sex[id] = '1';
        else if ( *se=='F' ) sample2sex[id] = '2';
        else error("Could not parse %s: %s\n", sex_fname,lines[i]);
    }
    for (i=0; i<nlines; i++) free(lines[i]);
    free(lines);
    for (i=0; i<bcf_hdr_nsamples(hdr); i++) 
        if ( !sample2sex[i] ) error("Missing sex for sample %s in %s\n", bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i),sex_fname);
    return sample2sex;
}

static void vcf_to_gensample(args_t *args)
{
    kstring_t str = {0,0,0};

    // insert chrom as first column if needed
    if(args->output_chrom_first_col)
        kputs("%CHROM ", &str);
    else
        kputs("%CHROM:%POS\\_%REF\\_%FIRST_ALT ", &str);

    // insert rsid as second column if needed
    if(args->output_vcf_ids)
        kputs("%ID ", &str);
    else
        kputs("%CHROM:%POS\\_%REF\\_%FIRST_ALT ", &str);

    kputs("%POS %REF %FIRST_ALT", &str);
    if ( !args->tag || !strcmp(args->tag,"GT") ) kputs("%_GT_TO_PROB3",&str);
    else if ( !strcmp(args->tag,"PL") ) kputs("%_PL_TO_PROB3",&str);
    else if ( !strcmp(args->tag,"GP") ) kputs("%_GP_TO_PROB3",&str);
    else error("todo: --tag %s\n", args->tag);
    kputs("\n", &str);
    open_vcf(args,str.s);

    int ret, gen_compressed = 1, sample_compressed = 0;
    char *gen_fname = NULL, *sample_fname = NULL;
    str.l = 0;
    kputs(args->outfname,&str);
    int n_files = 0, i;
    char **files = hts_readlist(str.s, 0, &n_files);
    if ( n_files==1 )
    {
        int l = str.l;
        kputs(".samples",&str);
        sample_fname = strdup(str.s);
        str.l = l;
        kputs(".gen.gz",&str);
        gen_fname = strdup(str.s);
    }
    else if ( n_files==2 )
    {
        if (strlen(files[0]) && strcmp(files[0],".")!=0) gen_fname = strdup(files[0]);
        if (strlen(files[1]) && strcmp(files[1],".")!=0) sample_fname = strdup(files[1]);
    }
    else
    {
        error("Error parsing --gensample filenames: %s\n", args->outfname);
    }
    for (i=0; i<n_files; i++) free(files[i]);
    free(files);

    if ( gen_fname && (strlen(gen_fname)<3 || strcasecmp(".gz",gen_fname+strlen(gen_fname)-3)) ) gen_compressed = 0;
    if ( sample_fname && strlen(sample_fname)>3 && strcasecmp(".gz",sample_fname+strlen(sample_fname)-3)==0 ) sample_compressed = 0;

    if (gen_fname) fprintf(pysam_stderr, "Gen file: %s\n", gen_fname);
    if (sample_fname) fprintf(pysam_stderr, "Sample file: %s\n", sample_fname);

    // write samples file
    if (sample_fname) 
    {
        char *sample2sex = NULL;
        if ( args->sex_fname ) sample2sex = init_sample2sex(args->header,args->sex_fname);

        int i;
        BGZF *sout = bgzf_open(sample_fname, sample_compressed ? "wg" : "wu");
        str.l = 0;
        kputs(sample2sex ? "ID_1 ID_2 missing sex\n0 0 0 0\n" : "ID_1 ID_2 missing\n0 0 0\n", &str);
        ret = bgzf_write(sout, str.s, str.l);
        if ( ret != str.l ) error("Error writing %s: %s\n", sample_fname, strerror(errno));
        for (i=0; i<bcf_hdr_nsamples(args->header); i++)
        {
            str.l = 0;
            if ( sample2sex )
                ksprintf(&str, "%s %s 0 %c\n", args->header->samples[i],args->header->samples[i],sample2sex[i]);
            else
                ksprintf(&str, "%s %s 0\n", args->header->samples[i],args->header->samples[i]);
            ret = bgzf_write(sout, str.s, str.l);
            if ( ret != str.l ) error("Error writing %s: %s\n", sample_fname, strerror(errno));
        }
        if ( bgzf_close(sout)!=0 ) error("Error closing %s: %s\n", sample_fname, strerror(errno));
        free(sample_fname);
        free(sample2sex);
    }
    if (!gen_fname) {
        if ( str.m ) free(str.s);
        return;
    }

    int prev_rid = -1, prev_pos = -1;
    int no_alt = 0, non_biallelic = 0, filtered = 0, ndup = 0, nok = 0;
    BGZF *gout = bgzf_open(gen_fname, gen_compressed ? "wg" : "wu");
    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = bcf_sr_get_line(args->files,0);
        if ( args->filter )
        {
            int pass = filter_test(args->filter, line, NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( !pass ) { filtered++; continue; }
        }

        // ALT allele is required
        if ( line->n_allele<2 ) { no_alt++; continue; }

        // biallelic required
        if ( line->n_allele>2 ) {
            if (!non_biallelic)
                fprintf(pysam_stderr, "Warning: non-biallelic records are skipped. Consider splitting multi-allelic records into biallelic records using 'bcftools norm -m-'.\n");
            non_biallelic++;
            continue;
        }

        // skip duplicate lines, or otherwise shapeit complains
        if ( prev_rid==line->rid && prev_pos==line->pos ) { ndup++; continue; }
        prev_rid = line->rid;
        prev_pos = line->pos;

        str.l = 0;
        convert_line(args->convert, line, &str);
        if ( str.l )
        {
            int ret = bgzf_write(gout, str.s, str.l);
            if ( ret!= str.l ) error("Error writing %s: %s\n", gen_fname,strerror(errno));
            nok++;
        }
    }
    fprintf(pysam_stderr, "%d records written, %d skipped: %d/%d/%d/%d no-ALT/non-biallelic/filtered/duplicated\n", 
        nok, no_alt+non_biallelic+filtered+ndup, no_alt, non_biallelic, filtered, ndup);

    if ( str.m ) free(str.s);
    if ( bgzf_close(gout)!=0 ) error("Error closing %s: %s\n", gen_fname,strerror(errno));
    free(gen_fname);
}

static void vcf_to_haplegendsample(args_t *args)
{
    kstring_t str = {0,0,0};
    if ( args->hap2dip )
        kputs("%_GT_TO_HAP2\n", &str);
    else
        kputs("%_GT_TO_HAP\n", &str);
    open_vcf(args,str.s);

    int ret, hap_compressed = 1, legend_compressed = 1, sample_compressed = 0;
    char *hap_fname = NULL, *legend_fname = NULL, *sample_fname = NULL;
    str.l = 0;
    kputs(args->outfname,&str);
    int n_files = 0, i;
    char **files = hts_readlist(str.s, 0, &n_files);
    if ( n_files==1 )
    {
        int l = str.l;
        kputs(".samples",&str);
        sample_fname = strdup(str.s);
        str.l = l;
        kputs(".legend.gz",&str);
        legend_fname = strdup(str.s);
        str.l = l;
        kputs(".hap.gz",&str);
        hap_fname = strdup(str.s);
    }
    else if ( n_files==3 )
    {
        if (strlen(files[0]) && strcmp(files[0],".")!=0) hap_fname = strdup(files[0]);
        if (strlen(files[1]) && strcmp(files[1],".")!=0) legend_fname = strdup(files[1]);
        if (strlen(files[2]) && strcmp(files[2],".")!=0) sample_fname = strdup(files[2]);
    }
    else
    {
        error("Error parsing --hapslegendsample filenames: %s\n", args->outfname);
    }
    for (i=0; i<n_files; i++) free(files[i]);
    free(files);

    if ( hap_fname && (strlen(hap_fname)<3 || strcasecmp(".gz",hap_fname+strlen(hap_fname)-3)) ) hap_compressed = 0;
    if ( legend_fname && (strlen(legend_fname)<3 || strcasecmp(".gz",legend_fname+strlen(legend_fname)-3)) ) legend_compressed = 0;
    if ( sample_fname && strlen(sample_fname)>3 && strcasecmp(".gz",sample_fname+strlen(sample_fname)-3)==0 ) sample_compressed = 0;

    if (hap_fname) fprintf(pysam_stderr, "Hap file: %s\n", hap_fname);
    if (legend_fname) fprintf(pysam_stderr, "Legend file: %s\n", legend_fname);
    if (sample_fname) fprintf(pysam_stderr, "Sample file: %s\n", sample_fname);

    // write samples file
    if (sample_fname)
    {
        char *sample2sex = NULL;
        if ( args->sex_fname ) sample2sex = init_sample2sex(args->header,args->sex_fname);
        
        int i;
        BGZF *sout = bgzf_open(sample_fname, sample_compressed ? "wg" : "wu");
        str.l = 0;
        kputs("sample population group sex\n", &str);
        ret = bgzf_write(sout, str.s, str.l);
        if ( ret != str.l ) error("Error writing %s: %s\n", sample_fname, strerror(errno));
        for (i=0; i<bcf_hdr_nsamples(args->header); i++)
        {
            str.l = 0;
            ksprintf(&str, "%s %s %s %c\n", args->header->samples[i], args->header->samples[i], args->header->samples[i], sample2sex ? sample2sex[i] : '2');
            ret = bgzf_write(sout, str.s, str.l);
            if ( ret != str.l ) error("Error writing %s: %s\n", sample_fname, strerror(errno));
        }
        if ( bgzf_close(sout)!=0 ) error("Error closing %s: %s\n", sample_fname, strerror(errno));
        free(sample_fname);
        free(sample2sex);
    }
    if (!hap_fname && !legend_fname) {
        if ( str.m ) free(str.s);
        return;
    }

    // open haps and legend outputs
    BGZF *hout = hap_fname ? bgzf_open(hap_fname, hap_compressed ? "wg" : "wu") : NULL;
    if ( hap_compressed && args->n_threads ) bgzf_thread_pool(hout, args->files->p->pool, args->files->p->qsize);
    BGZF *lout = legend_fname ? bgzf_open(legend_fname, legend_compressed ? "wg" : "wu") : NULL;
    if (legend_fname) {
        str.l = 0;
        kputs("id position a0 a1\n", &str);
        ret = bgzf_write(lout, str.s, str.l);
        if ( ret != str.l ) error("Error writing %s: %s\n", legend_fname, strerror(errno));
    }

    int no_alt = 0, non_biallelic = 0, filtered = 0, nok = 0;
    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = bcf_sr_get_line(args->files,0);
        if ( args->filter )
        {
            int pass = filter_test(args->filter, line, NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( !pass ) { filtered++; continue; }
        }

        // ALT allele is required
        if ( line->n_allele<2 ) { no_alt++; continue; }
        // biallelic required
        if ( line->n_allele>2 ) {
            if (!non_biallelic)
                fprintf(pysam_stderr, "Warning: non-biallelic records are skipped. Consider splitting multi-allelic records into biallelic records using 'bcftools norm -m-'.\n");
            non_biallelic++;
            continue;
        }

        str.l = 0;
        convert_line(args->convert, line, &str);
        if ( !str.l ) continue;

        // write haps file
        if (hap_fname) {
            ret = bgzf_write(hout, str.s, str.l); // write hap file
            if ( ret != str.l ) error("Error writing %s: %s\n", hap_fname, strerror(errno));
        }
        if (legend_fname) {
            str.l = 0;
            if ( args->output_vcf_ids && (line->d.id[0]!='.' || line->d.id[1]!=0) )
                ksprintf(&str, "%s %d %s %s\n", line->d.id, line->pos+1, line->d.allele[0], line->d.allele[1]);
            else
                ksprintf(&str, "%s:%d_%s_%s %d %s %s\n", bcf_seqname(args->header, line), line->pos+1, line->d.allele[0], line->d.allele[1], line->pos+1, line->d.allele[0], line->d.allele[1]);

            // write legend file
            ret = bgzf_write(lout, str.s, str.l);
            if ( ret != str.l ) error("Error writing %s: %s\n", legend_fname, strerror(errno));
        }
        nok++;
    }
    fprintf(pysam_stderr, "%d records written, %d skipped: %d/%d/%d no-ALT/non-biallelic/filtered\n", nok,no_alt+non_biallelic+filtered, no_alt, non_biallelic, filtered);
    if ( str.m ) free(str.s);
    if ( hout && bgzf_close(hout)!=0 ) error("Error closing %s: %s\n", hap_fname, strerror(errno));
    if ( lout && bgzf_close(lout)!=0 ) error("Error closing %s: %s\n", legend_fname, strerror(errno));
    if (hap_fname) free(hap_fname);
    if (legend_fname) free(legend_fname);
}

static void vcf_to_hapsample(args_t *args)
{
    /*
     *  WTCCC style haplotypes file
     *  see https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#hapsample
     *
     *  These are essentially the haplotypes from the impute2 format with some
     *  legend info tacked on to the first 5 columns
     *
     */
    kstring_t str = {0,0,0};

    // print ID instead of CHROM:POS_REF_ALT1
    if ( args->output_vcf_ids )
        kputs("%CHROM %ID %POS %REF %FIRST_ALT ", &str);
    else
        kputs("%CHROM %CHROM:%POS\\_%REF\\_%FIRST_ALT %POS %REF %FIRST_ALT ", &str);
    
    if ( args->hap2dip )
        kputs("%_GT_TO_HAP2\n", &str);
    else
        kputs("%_GT_TO_HAP\n", &str);
    open_vcf(args,str.s);

    int ret, hap_compressed = 1, sample_compressed = 0;
    char *hap_fname = NULL, *sample_fname = NULL;
    str.l = 0;
    kputs(args->outfname,&str);
    int n_files = 0, i;
    char **files = hts_readlist(str.s, 0, &n_files);
    if ( n_files==1 )
    {
        int l = str.l;
        kputs(".sample",&str);
        sample_fname = strdup(str.s);
        str.l = l;
        kputs(".hap.gz",&str);
        hap_fname = strdup(str.s);
    }
    else if ( n_files==2 )
    {
        if (strlen(files[0]) && strcmp(files[0],".")!=0) hap_fname = strdup(files[0]);
        if (strlen(files[1]) && strcmp(files[1],".")!=0) sample_fname = strdup(files[1]);
    }
    else
    {
        error("Error parsing --hapsample filenames: %s\n", args->outfname);
    }
    for (i=0; i<n_files; i++) free(files[i]);
    free(files);

    if ( hap_fname && (strlen(hap_fname)<3 || strcasecmp(".gz",hap_fname+strlen(hap_fname)-3)) ) hap_compressed = 0;
    if ( sample_fname && strlen(sample_fname)>3 && strcasecmp(".gz",sample_fname+strlen(sample_fname)-3)==0 ) sample_compressed = 0;

    if (hap_fname) fprintf(pysam_stderr, "Hap file: %s\n", hap_fname);
    if (sample_fname) fprintf(pysam_stderr, "Sample file: %s\n", sample_fname);

    // write samples file
    if (sample_fname)
    {
        char *sample2sex = NULL;
        if ( args->sex_fname ) sample2sex = init_sample2sex(args->header,args->sex_fname);

        int i;
        BGZF *sout = bgzf_open(sample_fname, sample_compressed ? "wg" : "wu");
        str.l = 0;
        kputs(sample2sex ? "ID_1 ID_2 missing sex\n0 0 0 0\n" : "ID_1 ID_2 missing\n0 0 0\n", &str);
        ret = bgzf_write(sout, str.s, str.l);
        if ( ret != str.l ) error("Error writing %s: %s\n", sample_fname, strerror(errno));
        for (i=0; i<bcf_hdr_nsamples(args->header); i++)
        {
            str.l = 0;
            if ( sample2sex )
                ksprintf(&str, "%s %s 0 %c\n", args->header->samples[i],args->header->samples[i],sample2sex[i]);
            else
                ksprintf(&str, "%s %s 0\n", args->header->samples[i],args->header->samples[i]);
            ret = bgzf_write(sout, str.s, str.l);
            if ( ret != str.l ) error("Error writing %s: %s\n", sample_fname, strerror(errno));
        }
        if ( bgzf_close(sout)!=0 ) error("Error closing %s: %s\n", sample_fname, strerror(errno));
        free(sample_fname);
        free(sample2sex);
    }
    if (!hap_fname) {
        if ( str.m ) free(str.s);
        return;
    }

    // open haps output
    BGZF *hout = hap_fname ? bgzf_open(hap_fname, hap_compressed ? "wg" : "wu") : NULL;
    if ( hap_compressed && args->n_threads ) bgzf_thread_pool(hout, args->files->p->pool, args->files->p->qsize);

    int no_alt = 0, non_biallelic = 0, filtered = 0, nok = 0;
    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = bcf_sr_get_line(args->files,0);
        if ( args->filter )
        {
            int pass = filter_test(args->filter, line, NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( !pass ) { filtered++; continue; }
        }

        // ALT allele is required
        if ( line->n_allele<2 ) { no_alt++; continue; }
        // biallelic required
        if ( line->n_allele>2 ) {
            if (!non_biallelic)
                fprintf(pysam_stderr, "Warning: non-biallelic records are skipped. Consider splitting multi-allelic records into biallelic records using 'bcftools norm -m-'.\n");
            non_biallelic++;
            continue;
        }

        str.l = 0;
        convert_line(args->convert, line, &str);
        if ( !str.l ) continue;

        // write haps file
        if (hap_fname) {
            ret = bgzf_write(hout, str.s, str.l); // write hap file
            if ( ret != str.l ) error("Error writing %s: %s\n", hap_fname, strerror(errno));
        }
        nok++;
    }
    fprintf(pysam_stderr, "%d records written, %d skipped: %d/%d/%d no-ALT/non-biallelic/filtered\n", nok, no_alt+non_biallelic+filtered, no_alt, non_biallelic, filtered);
    if ( str.m ) free(str.s);
    if ( hout && bgzf_close(hout)!=0 ) error("Error closing %s: %s\n", hap_fname, strerror(errno));
    if (hap_fname) free(hap_fname);
}

static void bcf_hdr_set_chrs(bcf_hdr_t *hdr, faidx_t *fai)
{
    int i, n = faidx_nseq(fai);
    for (i=0; i<n; i++)
    {
        const char *seq = faidx_iseq(fai,i);
        int len = faidx_seq_len(fai, seq);
        bcf_hdr_printf(hdr, "##contig=<ID=%s,length=%d>", seq,len);
    }
}
static inline int acgt_to_5(char base)
{
    if ( base=='A' ) return 0;
    if ( base=='C' ) return 1;
    if ( base=='G' ) return 2;
    if ( base=='T' ) return 3;
    return 4;
}
static inline int tsv_setter_aa1(args_t *args, char *ss, char *se, int alleles[], int *nals, int ref, int32_t *gts)
{
    if ( se - ss > 2 ) return -1;   // currently only SNPs

    if ( ss[0]=='-' )
    {
        // missing GT
        gts[0] = bcf_gt_missing;
        gts[1] = bcf_int32_vector_end;
        args->n.missing++;
        return 0;
    }
    if ( ss[0]=='I' ) return -2;    // skip insertions/deletions for now
    if ( ss[0]=='D' ) return -2;

    int a0 = acgt_to_5(toupper(ss[0]));
    int a1 = ss[1] ? acgt_to_5(toupper(ss[1])) : a0;
    if ( alleles[a0]<0 ) alleles[a0] = (*nals)++;
    if ( alleles[a1]<0 ) alleles[a1] = (*nals)++;

    gts[0] = bcf_gt_unphased(alleles[a0]); 
    gts[1] = ss[1] ? bcf_gt_unphased(alleles[a1]) : bcf_int32_vector_end;

    if ( ref==a0 && ref==a1  ) args->n.hom_rr++;    // hom ref: RR
    else if ( ref==a0 ) args->n.het_ra++;           // het: RA
    else if ( ref==a1 ) args->n.het_ra++;           // het: AR
    else if ( a0==a1 ) args->n.hom_aa++;            // hom-alt: AA
    else args->n.het_aa++;                          // non-ref het: AA

    return 0;
}
static int tsv_setter_aa(tsv_t *tsv, bcf1_t *rec, void *usr)
{
    args_t *args = (args_t*) usr;

    int len;
    char *ref = faidx_fetch_seq(args->ref, (char*)bcf_hdr_id2name(args->header,rec->rid), rec->pos, rec->pos, &len);
    if ( !ref ) error("faidx_fetch_seq failed at %s:%d\n", bcf_hdr_id2name(args->header,rec->rid), rec->pos+1);

    int nals = 1, alleles[5] = { -1, -1, -1, -1, -1 };    // a,c,g,t,n
    ref[0] = toupper(ref[0]);
    int iref = acgt_to_5(ref[0]);
    alleles[iref] = 0;

    rec->n_sample = bcf_hdr_nsamples(args->header);

    int i, ret;
    for (i=0; i<rec->n_sample; i++)
    {
        if ( i>0 )
        {
            ret = tsv_next(tsv);
            if ( ret==-1 ) error("Too few columns for %d samples at %s:%d\n", rec->n_sample,bcf_hdr_id2name(args->header,rec->rid), rec->pos+1);
        }
        ret = tsv_setter_aa1(args, tsv->ss, tsv->se, alleles, &nals, iref, args->gts+i*2);
        if ( ret==-1 ) error("Error parsing the site %s:%d, expected two characters\n", bcf_hdr_id2name(args->header,rec->rid), rec->pos+1);
        if ( ret==-2 ) 
        {
            // something else than a SNP
            free(ref);
            return -1;
        }
    }

    args->str.l = 0;
    kputc(ref[0], &args->str);
    for (i=0; i<5; i++) 
    {
        if ( alleles[i]>0 )
        {
            kputc(',', &args->str);
            kputc("ACGTN"[i], &args->str);
        }
    }
    bcf_update_alleles_str(args->header, rec, args->str.s);
    if ( bcf_update_genotypes(args->header,rec,args->gts,rec->n_sample*2) ) error("Could not update the GT field\n");

    free(ref);
    return 0;
}

static void tsv_to_vcf(args_t *args)
{
    if ( !args->ref_fname ) error("--tsv2vcf requires the --fasta-ref option\n");
    if ( !args->sample_list ) error("--tsv2vcf requires the --samples option\n");

    args->ref = fai_load(args->ref_fname);
    if ( !args->ref ) error("Could not load the reference %s\n", args->ref_fname);

    args->header = bcf_hdr_init("w");
    bcf_hdr_set_chrs(args->header, args->ref);
    bcf_hdr_append(args->header, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    if (args->record_cmd_line) bcf_hdr_append_version(args->header, args->argc, args->argv, "bcftools_convert");

    int i, n;
    char **smpls = hts_readlist(args->sample_list, args->sample_is_file, &n);
    if ( !smpls ) error("Could not parse %s\n", args->sample_list);
    for (i=0; i<n; i++)
    {
        bcf_hdr_add_sample(args->header, smpls[i]);
        free(smpls[i]);
    }
    free(smpls);
    bcf_hdr_add_sample(args->header, NULL);
    args->gts = (int32_t *) malloc(sizeof(int32_t)*n*2);

    htsFile *out_fh = hts_open(args->outfname,hts_bcf_wmode(args->output_type));
    if ( out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->outfname, strerror(errno));
    if ( args->n_threads ) hts_set_threads(out_fh, args->n_threads);
    bcf_hdr_write(out_fh,args->header);

    tsv_t *tsv = tsv_init(args->columns ? args->columns : "ID,CHROM,POS,AA");
    if ( tsv_register(tsv, "CHROM", tsv_setter_chrom, args->header) < 0 ) error("Expected CHROM column\n");
    if ( tsv_register(tsv, "POS", tsv_setter_pos, NULL) < 0 ) error("Expected POS column\n");
    if ( tsv_register(tsv, "ID", tsv_setter_id, args->header) < 0 && !args->columns ) error("Expected ID column\n");
    if ( tsv_register(tsv, "AA", tsv_setter_aa, args) < 0 ) error("Expected AA column\n");

    bcf1_t *rec = bcf_init();
    bcf_float_set_missing(rec->qual);

    kstring_t line = {0,0,0};
    htsFile *in_fh = hts_open(args->infname, "r");
    if ( !in_fh ) error("Could not read: %s\n", args->infname);
    while ( hts_getline(in_fh, KS_SEP_LINE, &line) > 0 )
    {
        if ( line.s[0]=='#' ) continue;     // skip comments
        bcf_clear(rec);

        args->n.total++;
        if ( !tsv_parse(tsv, rec, line.s) )
            bcf_write(out_fh, args->header, rec);
        else
            args->n.skipped++;
    }
    if ( hts_close(in_fh) ) error("Close failed: %s\n", args->infname);
    free(line.s);

    bcf_hdr_destroy(args->header);
    hts_close(out_fh);
    tsv_destroy(tsv);
    bcf_destroy(rec);
    free(args->str.s);
    free(args->gts);

    fprintf(pysam_stderr,"Rows total: \t%d\n", args->n.total);
    fprintf(pysam_stderr,"Rows skipped: \t%d\n", args->n.skipped);
    fprintf(pysam_stderr,"Missing GTs: \t%d\n", args->n.missing);
    fprintf(pysam_stderr,"Hom RR: \t%d\n", args->n.hom_rr);
    fprintf(pysam_stderr,"Het RA: \t%d\n", args->n.het_ra);
    fprintf(pysam_stderr,"Hom AA: \t%d\n", args->n.hom_aa);
    fprintf(pysam_stderr,"Het AA: \t%d\n", args->n.het_aa);
}

static void vcf_to_vcf(args_t *args)
{
    open_vcf(args,NULL);
    htsFile *out_fh = hts_open(args->outfname,hts_bcf_wmode(args->output_type));
    if ( out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->outfname, strerror(errno));
    if ( args->n_threads ) hts_set_threads(out_fh, args->n_threads);

    bcf_hdr_t *hdr = bcf_sr_get_header(args->files,0);
    bcf_hdr_write(out_fh,hdr);

    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = bcf_sr_get_line(args->files,0);
        if ( args->filter )
        {
            int pass = filter_test(args->filter, line, NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( !pass ) continue;
        }
        bcf_write(out_fh,hdr,line);
    }
    hts_close(out_fh);
}

static void gvcf_to_vcf(args_t *args)
{
    if ( !args->ref_fname ) error("--gvcf2vcf requires the --fasta-ref option\n");

    args->ref = fai_load(args->ref_fname);
    if ( !args->ref ) error("Could not load the fai index for reference %s\n", args->ref_fname);

    open_vcf(args,NULL);
    htsFile *out_fh = hts_open(args->outfname,hts_bcf_wmode(args->output_type));
    if ( out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->outfname, strerror(errno));
    if ( args->n_threads ) hts_set_threads(out_fh, args->n_threads);

    bcf_hdr_t *hdr = bcf_sr_get_header(args->files,0);
    if (args->record_cmd_line) bcf_hdr_append_version(hdr, args->argc, args->argv, "bcftools_convert");
    bcf_hdr_write(out_fh,hdr);

    int32_t *itmp = NULL, nitmp = 0;

    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = bcf_sr_get_line(args->files,0);
        if ( args->filter )
        {
            int pass = filter_test(args->filter, line, NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( !pass ) continue;
        }

        if (!bcf_has_filter(hdr,line,"PASS"))
        {
            bcf_write(out_fh,hdr,line);
            continue;
        }

        // check if alleles compatible with being a gVCF record
        int i, gallele = -1;
        if (line->n_allele==1)
            gallele = 0; // illumina/bcftools-call gvcf (if INFO/END present)
        else
        {
            if ( line->d.allele[1][0]!='<' ) continue;
            for (i=1; i<line->n_allele; i++)
            {
                if ( line->d.allele[i][1]=='*' && line->d.allele[i][2]=='>' && line->d.allele[i][3]=='\0' ) { gallele = i; break; } // mpileup/spec compliant gVCF
                if ( line->d.allele[i][1]=='X' && line->d.allele[i][2]=='>' && line->d.allele[i][3]=='\0' ) { gallele = i; break; } // old mpileup gVCF
                if ( strcmp(line->d.allele[i],"<NON_REF>")==0 ) { gallele = i; break; }               // GATK gVCF
            }
        }

        // no gVCF compatible alleles
        if (gallele<0)
        {
            bcf_write(out_fh,hdr,line);
            continue;
        }

        int nend = bcf_get_info_int32(hdr,line,"END",&itmp,&nitmp);
        if ( nend!=1 )
        {
            // No INFO/END => not gVCF record
            bcf_write(out_fh,hdr,line);
            continue;
        }
        bcf_update_info_int32(hdr,line,"END",NULL,0);
        int pos, len;
        for (pos=line->pos; pos<itmp[0]; pos++)
        {
            line->pos = pos;
            char *ref = faidx_fetch_seq(args->ref, (char*)bcf_hdr_id2name(hdr,line->rid), line->pos, line->pos, &len);
            if ( !ref ) error("faidx_fetch_seq failed at %s:%d\n", bcf_hdr_id2name(hdr,line->rid), line->pos+1);
            strncpy(line->d.allele[0],ref,len);
            bcf_write(out_fh,hdr,line);
            free(ref);
        }
    }
    free(itmp);
    hts_close(out_fh);
}

static void usage(void)
{
    fprintf(pysam_stderr, "\n");
    fprintf(pysam_stderr, "About:   Converts VCF/BCF to other formats and back. See man page for file\n");
    fprintf(pysam_stderr, "         formats details. When specifying output files explicitly instead\n");
    fprintf(pysam_stderr, "         of with <prefix>, one can use '-' for pysam_stdout and '.' to suppress.\n");
    fprintf(pysam_stderr, "Usage:   bcftools convert [OPTIONS] <input_file>\n");
    fprintf(pysam_stderr, "\n");
    fprintf(pysam_stderr, "VCF input options:\n");
    fprintf(pysam_stderr, "   -e, --exclude <expr>        exclude sites for which the expression is true\n");
    fprintf(pysam_stderr, "   -i, --include <expr>        select sites for which the expression is true\n");
    fprintf(pysam_stderr, "   -r, --regions <region>      restrict to comma-separated list of regions\n");
    fprintf(pysam_stderr, "   -R, --regions-file <file>   restrict to regions listed in a file\n");
    fprintf(pysam_stderr, "   -s, --samples <list>        list of samples to include\n");
    fprintf(pysam_stderr, "   -S, --samples-file <file>   file of samples to include\n");
    fprintf(pysam_stderr, "   -t, --targets <region>      similar to -r but streams rather than index-jumps\n");
    fprintf(pysam_stderr, "   -T, --targets-file <file>   similar to -R but streams rather than index-jumps\n");
    fprintf(pysam_stderr, "\n");
    fprintf(pysam_stderr, "VCF output options:\n");
    fprintf(pysam_stderr, "       --no-version               do not append version and command line to the header\n");
    fprintf(pysam_stderr, "   -o, --output <file>            output file name [pysam_stdout]\n");
    fprintf(pysam_stderr, "   -O, --output-type <b|u|z|v>    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(pysam_stderr, "       --threads <int>            number of extra output compression threads [0]\n");
    fprintf(pysam_stderr, "\n");
    fprintf(pysam_stderr, "GEN/SAMPLE conversion (input/output from IMPUTE2):\n");
    fprintf(pysam_stderr, "   -G, --gensample2vcf <...>   <prefix>|<gen-file>,<sample-file>\n");
    fprintf(pysam_stderr, "   -g, --gensample <...>       <prefix>|<gen-file>,<sample-file>\n");
    fprintf(pysam_stderr, "       --tag <string>          tag to take values for .gen file: GT,PL,GL,GP [GT]\n");
    fprintf(pysam_stderr, "       --chrom                 output chromosome in first column instead of CHROM:POS_REF_ALT\n");
    fprintf(pysam_stderr, "       --sex <file>            output sex column in the sample-file, input format is: Sample\\t[MF]\n");
    fprintf(pysam_stderr, "       --vcf-ids               output VCF IDs in second column instead of CHROM:POS_REF_ALT\n");
    fprintf(pysam_stderr, "\n");
    fprintf(pysam_stderr, "gVCF conversion:\n");
    fprintf(pysam_stderr, "       --gvcf2vcf              expand gVCF reference blocks\n");
    fprintf(pysam_stderr, "   -f, --fasta-ref <file>      reference sequence in fasta format\n");
    fprintf(pysam_stderr, "\n");
    fprintf(pysam_stderr, "HAP/SAMPLE conversion (output from SHAPEIT):\n");
    fprintf(pysam_stderr, "       --hapsample2vcf <...>   <prefix>|<hap-file>,<sample-file>\n");
    fprintf(pysam_stderr, "       --hapsample <...>       <prefix>|<hap-file>,<sample-file>\n");
    fprintf(pysam_stderr, "       --haploid2diploid       convert haploid genotypes to diploid homozygotes\n");
    fprintf(pysam_stderr, "       --sex <file>            output sex column in the sample-file, input format is: Sample\\t[MF]\n");
    fprintf(pysam_stderr, "       --vcf-ids               output VCF IDs instead of CHROM:POS_REF_ALT\n");
    fprintf(pysam_stderr, "\n");
    fprintf(pysam_stderr, "HAP/LEGEND/SAMPLE conversion:\n");
    fprintf(pysam_stderr, "   -H, --haplegendsample2vcf <...>  <prefix>|<hap-file>,<legend-file>,<sample-file>\n");
    fprintf(pysam_stderr, "   -h, --haplegendsample <...>      <prefix>|<hap-file>,<legend-file>,<sample-file>\n");
    fprintf(pysam_stderr, "       --haploid2diploid            convert haploid genotypes to diploid homozygotes\n");
    fprintf(pysam_stderr, "       --sex <file>                 output sex column in the sample-file, input format is: Sample\\t[MF]\n");
    fprintf(pysam_stderr, "       --vcf-ids                    output VCF IDs instead of CHROM:POS_REF_ALT\n");
    fprintf(pysam_stderr, "\n");
    fprintf(pysam_stderr, "TSV conversion:\n");
    fprintf(pysam_stderr, "       --tsv2vcf <file>        \n");
    fprintf(pysam_stderr, "   -c, --columns <string>      columns of the input tsv file [ID,CHROM,POS,AA]\n");
    fprintf(pysam_stderr, "   -f, --fasta-ref <file>      reference sequence in fasta format\n");
    fprintf(pysam_stderr, "   -s, --samples <list>        list of sample names\n");
    fprintf(pysam_stderr, "   -S, --samples-file <file>   file of sample names\n");
    fprintf(pysam_stderr, "\n");
    // fprintf(pysam_stderr, "PLINK options:\n");
    // fprintf(pysam_stderr, "   -p, --plink <prefix>|<ped>,<map>,<fam>|<bed>,<bim>,<fam>|<tped>,<tfam>\n");
    // fprintf(pysam_stderr, "       --tped              make tped file instead\n");
    // fprintf(pysam_stderr, "       --bin               make binary bed/fam/bim files\n");
    // fprintf(pysam_stderr, "\n");
    // fprintf(pysam_stderr, "PBWT options:\n");
    // fprintf(pysam_stderr, "   -b, --pbwt          <prefix> or <pbwt>,<sites>,<sample>,<missing>\n");
    // fprintf(pysam_stderr, "\n");
    exit(1);
}

int main_vcfconvert(int argc, char *argv[])
{
    int c;
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    args->outfname = "-";
    args->output_type = FT_VCF;
    args->n_threads = 0;
    args->record_cmd_line = 1;

    static struct option loptions[] =
    {
        {"include",required_argument,NULL,'i'},
        {"exclude",required_argument,NULL,'e'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"threads",required_argument,NULL,9},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"targets",required_argument,NULL,'t'},
        {"targets-file",required_argument,NULL,'T'},
        {"samples",required_argument,NULL,'s'},
        {"samples-file",required_argument,NULL,'S'},
        {"sex",required_argument,NULL,11},
        {"gensample",required_argument,NULL,'g'},
        {"gensample2vcf",required_argument,NULL,'G'},
        {"tag",required_argument,NULL,1},
        {"chrom",no_argument,NULL,8},        
        {"tsv2vcf",required_argument,NULL,2},
        {"hapsample",required_argument,NULL,7},
        {"hapsample2vcf",required_argument,NULL,3},
        {"vcf-ids",no_argument,NULL,4},
        {"haploid2diploid",no_argument,NULL,5},
        {"gvcf2vcf",no_argument,NULL,6},
        {"haplegendsample",required_argument,NULL,'h'},
        {"haplegendsample2vcf",required_argument,NULL,'H'},
        {"columns",required_argument,NULL,'c'},
        {"fasta-ref",required_argument,NULL,'f'},
        {"no-version",no_argument,NULL,10},
        {NULL,0,NULL,0}
    };
    while ((c = getopt_long(argc, argv, "?h:r:R:s:S:t:T:i:e:g:G:o:O:c:f:H:",loptions,NULL)) >= 0) {
        switch (c) {
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; args->regions_is_file = 1; break;
            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; args->targets_is_file = 1; break;
            case 's': args->sample_list = optarg; break;
            case 'S': args->sample_list = optarg; args->sample_is_file = 1; break;
            case 'g': args->convert_func = vcf_to_gensample; args->outfname = optarg; break;
            case 'G': args->convert_func = gensample_to_vcf; args->infname = optarg; break;
            case  1 : args->tag = optarg; break;
            case  2 : args->convert_func = tsv_to_vcf; args->infname = optarg; break;
            case  3 : args->convert_func = hapsample_to_vcf; args->infname = optarg; break;
            case  4 : args->output_vcf_ids = 1; break;
            case  5 : args->hap2dip = 1; break;
            case  6 : args->convert_func = gvcf_to_vcf; break;
            case  7 : args->convert_func = vcf_to_hapsample; args->outfname = optarg; break;
            case  8 : args->output_chrom_first_col = 1; break;
            case 'H': args->convert_func = haplegendsample_to_vcf; args->infname = optarg; break;
            case 'f': args->ref_fname = optarg; break;
            case 'c': args->columns = optarg; break;
            case 'o': args->outfname = optarg; break;
            case 'O':
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                }
                break;
            case 'h': args->convert_func = vcf_to_haplegendsample; args->outfname = optarg; break;
            case  9 : args->n_threads = strtol(optarg, 0, 0); break;
            case 10 : args->record_cmd_line = 0; break;
            case 11 : args->sex_fname = optarg; break;
            case '?': usage();
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( !args->infname )
    {
        if ( optind>=argc )
        {
            if ( !isatty(fileno((FILE *)stdin)) ) args->infname = "-";
        }
        else args->infname = argv[optind];
    }
    if ( !args->infname ) usage();
    
    if ( args->convert_func ) args->convert_func(args);
    else vcf_to_vcf(args);

    destroy_data(args);
    free(args);
    return 0;
}
