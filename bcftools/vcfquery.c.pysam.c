#include "bcftools.pysam.h"

/*  vcfquery.c -- Extracts fields from VCF/BCF file.

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
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/khash_str2int.h>
#include <htslib/vcfutils.h>
#include "bcftools.h"
#include "filter.h"
#include "convert.h"
#include "smpl_ilist.h"


// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef struct
{
    filter_t *filter;
    char *filter_str;
    int filter_logic;   // include or exclude sites which match the filters? One of FLT_INCLUDE/FLT_EXCLUDE
    uint8_t *smpl_pass;
    convert_t *convert;
    bcf_srs_t *files;
    bcf_hdr_t *header;
    int sample_is_file;
    char **argv, *format_str, *sample_list, *targets_list, *regions_list, *vcf_list, *fn_out;
    int argc, list_columns, print_header, allow_undef_tags, force_samples;
    FILE *out;
}
args_t;

static void destroy_list(char **list, int n)
{
    int i;
    for (i=0; i<n; i++)
        free(list[i]);
    free(list);
}

static void init_data(args_t *args)
{
    args->header = args->files->readers[0].header;

    int i, nsamples = 0, *samples = NULL;
    if ( args->sample_list && strcmp("-",args->sample_list) )
    {
        for (i=0; i<args->files->nreaders; i++)
        {
            // This tells htslib to subset samples directly when reading. Also the header is modified to
            // include only the requested samples
            int ret = bcf_hdr_set_samples(args->files->readers[i].header,args->sample_list,args->sample_is_file);
            if ( ret<0 ) error("Error parsing the sample list\n");
            else if ( ret>0 && !args->force_samples )
                error("Error: sample #%d not found in the header, user --force-samples to proceed anyway\n", ret);
        }

        int flags = SMPL_REORDER;
        smpl_ilist_t *ilist = smpl_ilist_init(args->files->readers[0].header, args->sample_list, args->sample_is_file, flags);
        nsamples = ilist->n;
        samples = (int*) malloc(sizeof(int)*nsamples);
        for (i=0; i<ilist->n; i++)
            samples[i] = ilist->idx[i];
        smpl_ilist_destroy(ilist);
    }
    args->convert = convert_init(args->header, samples, nsamples, args->format_str);
    convert_set_option(args->convert, subset_samples, &args->smpl_pass);
    if ( args->allow_undef_tags ) convert_set_option(args->convert, allow_undef_tags, 1);
    free(samples);

    int max_unpack = convert_max_unpack(args->convert);
    if ( args->filter_str )
    {
        args->filter = filter_init(args->header, args->filter_str);
        max_unpack |= filter_max_unpack(args->filter);
    }
    args->files->max_unpack = max_unpack;
}

static void destroy_data(args_t *args)
{
    convert_destroy(args->convert);
    if ( args->filter )
        filter_destroy(args->filter);
}

static void query_vcf(args_t *args)
{
    kstring_t str = {0,0,0};

    if ( args->print_header )
    {
        convert_header(args->convert,&str);
        if ( fwrite(str.s, str.l, 1, args->out)!=1 ) error("[%s] Error: cannot write to %s\n", __func__,args->fn_out?args->fn_out:"standard output");
    }

    int i,max_convert_unpack = convert_max_unpack(args->convert);
    int max_filter_unpack = args->filter ? filter_max_unpack(args->filter) : 0;
    while ( bcf_sr_next_line(args->files) )
    {
        if ( !bcf_sr_has_line(args->files,0) ) continue;
        bcf1_t *line = args->files->readers[0].buffer[0];
        bcf_unpack(line, args->files->max_unpack);

        if ( args->filter )
        {
            int pass = filter_test(args->filter, line, (const uint8_t**) &args->smpl_pass);
            if ( args->filter_logic & FLT_EXCLUDE )
            {
                // This code addresses this problem:
                //  -i can include a site but exclude a sample
                //  -e exclude a site but include a sample

                if ( pass )
                {
                    if ( !args->smpl_pass ) continue;
                    if ( !(max_convert_unpack & BCF_UN_FMT) && !(max_filter_unpack & BCF_UN_FMT) ) continue;

                    pass = 0;
                    for (i=0; i<line->n_sample; i++)
                    {
                        if ( args->smpl_pass[i] ) args->smpl_pass[i] = 0;
                        else { args->smpl_pass[i] = 1; pass = 1; }
                    }
                    if ( !pass ) continue;
                }
                else if ( args->smpl_pass )
                    for (i=0; i<line->n_sample; i++) args->smpl_pass[i] = 1;
            }
            else if ( !pass ) continue;
        }

        str.l = 0;
        convert_line(args->convert, line, &str);
        if ( str.l && fwrite(str.s, str.l, 1, args->out)!=1 ) error("[%s] Error: cannot write to %s\n", __func__,args->fn_out?args->fn_out:"standard output");
    }
    if ( str.m ) free(str.s);
}

static void list_columns(args_t *args)
{
    int negate = 0;
    int i;
    bcf_sr_t *reader = &args->files->readers[0];
    void *has_sample = NULL;
    if ( args->sample_list )
    {
        if ( args->sample_list[0]=='^' ) negate = 1;
        has_sample = khash_str2int_init();
        int i, nsmpl;
        char **smpl = hts_readlist(negate ? args->sample_list+1 : args->sample_list, args->sample_is_file, &nsmpl);
        if ( !smpl ) error("Error: failed to read %s\n", negate ? args->sample_list+1 : args->sample_list);
        for (i=0; i<nsmpl; i++)
        {
            if ( bcf_hdr_id2int(reader->header,BCF_DT_SAMPLE,smpl[i])<0 && !args->force_samples )
                error("Error: sample #%d not found in the header, user --force-samples to proceed anyway\n", i+1);
            khash_str2int_inc(has_sample, smpl[i]);
        }
        free(smpl);
    }

    for (i=0; i<bcf_hdr_nsamples(reader->header); i++)
    {
        int skip = 0;
        if ( negate )
        {
            if ( khash_str2int_has_key(has_sample, reader->header->samples[i]) ) skip = 1;
        }
        else if ( has_sample && !khash_str2int_has_key(has_sample, reader->header->samples[i]) ) skip = 1;
        if ( skip ) continue;
        fprintf(bcftools_stdout, "%s\n", reader->header->samples[i]);
    }

    if ( has_sample )
        khash_str2int_destroy_free(has_sample);
}

static char **copy_header(bcf_hdr_t *hdr, char **src, int nsrc)
{
    char **dst = (char**) malloc(sizeof(char*)*nsrc);
    int i;
    for (i=0; i<nsrc; i++) dst[i] = strdup(src[i]);
    return dst;
}
static int compare_header(bcf_hdr_t *hdr, char **a, int na, char **b, int nb)
{
    if ( na!=nb ) return na-nb;
    int i;
    for (i=0; i<na; i++)
        if ( strcmp(a[i],b[i]) ) return 1;
    return 0;
}


static void usage(void)
{
    fprintf(bcftools_stderr, "\n");
    fprintf(bcftools_stderr, "About:   Extracts fields from VCF/BCF file and prints them in user-defined format\n");
    fprintf(bcftools_stderr, "Usage:   bcftools query [options] <A.vcf.gz> [<B.vcf.gz> [...]]\n");
    fprintf(bcftools_stderr, "\n");
    fprintf(bcftools_stderr, "Options:\n");
    fprintf(bcftools_stderr, "    -e, --exclude EXPR                Exclude sites for which the expression is true (see man page for details)\n");
    fprintf(bcftools_stderr, "        --force-samples               Only warn about unknown subset samples\n");
    fprintf(bcftools_stderr, "    -f, --format STRING               See man page for details\n");
    fprintf(bcftools_stderr, "    -H, --print-header                Print header\n");
    fprintf(bcftools_stderr, "    -i, --include EXPR                Select sites for which the expression is true (see man page for details)\n");
    fprintf(bcftools_stderr, "    -l, --list-samples                Print the list of samples and exit\n");
    fprintf(bcftools_stderr, "    -o, --output FILE                 Output file name [bcftools_stdout]\n");
    fprintf(bcftools_stderr, "    -r, --regions REGION              Restrict to comma-separated list of regions\n");
    fprintf(bcftools_stderr, "    -R, --regions-file FILE           Restrict to regions listed in a file\n");
    fprintf(bcftools_stderr, "        --regions-overlap 0|1|2       Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n");
    fprintf(bcftools_stderr, "    -s, --samples LIST                List of samples to include\n");
    fprintf(bcftools_stderr, "    -S, --samples-file FILE           File of samples to include\n");
    fprintf(bcftools_stderr, "    -t, --targets REGION              Similar to -r but streams rather than index-jumps\n");
    fprintf(bcftools_stderr, "    -T, --targets-file FILE           Similar to -R but streams rather than index-jumps\n");
    fprintf(bcftools_stderr, "        --targets-overlap 0|1|2       Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n");
    fprintf(bcftools_stderr, "    -u, --allow-undef-tags            Print \".\" for undefined tags\n");
    fprintf(bcftools_stderr, "    -v, --vcf-list FILE               Process multiple VCFs listed in the file\n");
    fprintf(bcftools_stderr, "\n");
    fprintf(bcftools_stderr, "Examples:\n");
    fprintf(bcftools_stderr, "\tbcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%SAMPLE=%%GT]\\n' file.vcf.gz\n");
    fprintf(bcftools_stderr, "\n");
    bcftools_exit(1);
}

int main_vcfquery(int argc, char *argv[])
{
    int c, collapse = 0;
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    int regions_is_file = 0, targets_is_file = 0;
    int regions_overlap = 1;
    int targets_overlap = 0;

    static struct option loptions[] =
    {
        {"help",0,0,'h'},
        {"list-samples",0,0,'l'},
        {"include",1,0,'i'},
        {"exclude",1,0,'e'},
        {"format",1,0,'f'},
        {"force-samples",0,0,3},
        {"output-file",1,0,'o'},
        {"output",1,0,'o'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"regions-overlap",required_argument,NULL,1},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
        {"targets-overlap",required_argument,NULL,2},
        {"annots",1,0,'a'},
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {"print-header",0,0,'H'},
        {"collapse",1,0,'c'},
        {"vcf-list",1,0,'v'},
        {"allow-undef-tags",0,0,'u'},
        {0,0,0,0}
    };
    while ((c = getopt_long(argc, argv, "hlr:R:f:a:s:S:Ht:T:c:v:i:e:o:u",loptions,NULL)) >= 0) {
        switch (c) {
            case 'o': args->fn_out = optarg; break;
            case 'f': args->format_str = strdup(optarg); break;
            case 'H': args->print_header = 1; break;
            case 'v': args->vcf_list = optarg; break;
            case 'c':
                error("The --collapse option is obsolete, pipe through `bcftools norm -c` instead.\n");
                break;
            case 'a':
                {
                    kstring_t str = {0,0,0};
                    kputs("%CHROM\t%POS\t%MASK\t%REF\t%ALT\t%", &str);
                    char *p = optarg;
                    while ( *p )
                    {
                        if ( *p==',' )
                            kputs("\t%", &str);
                        else
                            kputc(*p, &str);
                        p++;
                    }
                    kputc('\n', &str);
                    args->format_str = str.s;
                    break;
                }
            case 'e':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; targets_is_file = 1; break;
            case 'l': args->list_columns = 1; break;
            case 'u': args->allow_undef_tags = 1; break;
            case 's': args->sample_list = optarg; break;
            case 'S': args->sample_list = optarg; args->sample_is_file = 1; break;
            case  1 :
                regions_overlap = parse_overlap_option(optarg);
                if ( regions_overlap < 0 ) error("Could not parse: --regions-overlap %s\n",optarg);
                break;
            case  2 :
                targets_overlap = parse_overlap_option(optarg);
                if ( targets_overlap < 0 ) error("Could not parse: --targets-overlap %s\n",optarg);
                break;
            case  3 : args->force_samples = 1; break;
            case 'h':
            case '?': usage(); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    char *fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";
    }
    else fname = argv[optind];

    if ( args->list_columns )
    {
        if ( !fname ) error("Missing the VCF file name\n");
        args->files = bcf_sr_init();
        if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to read from %s: %s\n", !strcmp("-",fname)?"standard input":fname,bcf_sr_strerror(args->files->errnum));
        list_columns(args);
        bcf_sr_destroy(args->files);
        free(args);
        return 0;
    }

    if ( !args->format_str )
    {
        if ( argc==1 && !fname ) usage();
        error("Error: Missing the --format option\n");
    }
    args->out = args->fn_out ? fopen(args->fn_out, "w") : bcftools_stdout;
    if ( !args->out ) error("%s: %s\n", args->fn_out,strerror(errno));

    if ( !args->vcf_list )
    {
        if ( !fname ) usage();
        args->files = bcf_sr_init();
        if ( optind+1 < argc ) args->files->require_index = 1;
        if ( args->regions_list )
        {
            bcf_sr_set_opt(args->files,BCF_SR_REGIONS_OVERLAP,regions_overlap);
            if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
                error("Failed to read the regions: %s\n", args->regions_list);
        }
        if ( args->targets_list )
        {
            bcf_sr_set_opt(args->files,BCF_SR_TARGETS_OVERLAP,targets_overlap);
            if ( bcf_sr_set_targets(args->files, args->targets_list, targets_is_file, 0)<0 )
                error("Failed to read the targets: %s\n", args->targets_list);
        }
        while ( fname )
        {
            if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to read from %s: %s\n", !strcmp("-",fname)?"standard input":fname,bcf_sr_strerror(args->files->errnum));
            fname = ++optind < argc ? argv[optind] : NULL;
        }
        init_data(args);
        query_vcf(args);
        free(args->format_str);
        destroy_data(args);
        bcf_sr_destroy(args->files);
        if ( fclose(args->out)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->fn_out);
        free(args);
        return 0;
    }

    // multiple VCFs
    int i, k, nfiles, prev_nsamples = 0;
    char **fnames, **prev_samples = NULL;
    fnames = hts_readlist(args->vcf_list, 1, &nfiles);
    if ( !fnames ) error("Error: failed to read %s\n", args->vcf_list);
    if ( !nfiles ) error("No files in %s?\n", args->vcf_list);
    for (i=0; i<nfiles; i++)
    {
        args->files = bcf_sr_init();
        args->files->collapse = collapse;
        if ( args->regions_list && bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
        if ( optind < argc ) args->files->require_index = 1;
        if ( args->targets_list )
        {
            if ( bcf_sr_set_targets(args->files, args->targets_list,targets_is_file, 0)<0 )
                error("Failed to read the targets: %s\n", args->targets_list);
        }
        if ( !bcf_sr_add_reader(args->files, fnames[i]) ) error("Failed to open %s: %s\n", fnames[i],bcf_sr_strerror(args->files->errnum));
        for (k=optind; k<argc; k++)
            if ( !bcf_sr_add_reader(args->files, argv[k]) ) error("Failed to open %s: %s\n", argv[k],bcf_sr_strerror(args->files->errnum));
        init_data(args);
        if ( i==0 )
        {
            prev_samples = copy_header(args->header, args->files->readers[0].header->samples, bcf_hdr_nsamples(args->files->readers[0].header));
            prev_nsamples = bcf_hdr_nsamples(args->files->readers[0].header);
        }
        else
        {
            args->print_header = 0;
            if ( compare_header(args->header, args->files->readers[0].header->samples, bcf_hdr_nsamples(args->files->readers[0].header), prev_samples, prev_nsamples) )
                error("Different samples in %s and %s\n", fnames[i-1],fnames[i]);
        }
        query_vcf(args);
        destroy_data(args);
        bcf_sr_destroy(args->files);
    }
    if ( fclose(args->out)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->fn_out);;
    destroy_list(fnames, nfiles);
    destroy_list(prev_samples, prev_nsamples);
    free(args->format_str);
    free(args);
    return 0;
}


