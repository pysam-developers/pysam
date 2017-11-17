/*  vcfquery.c -- Extracts fields from VCF/BCF file.

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
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/khash_str2int.h>
#include <htslib/vcfutils.h>
#include "bcftools.h"
#include "filter.h"
#include "convert.h"


// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef struct
{
    filter_t *filter;
    char *filter_str;
    int filter_logic;   // include or exclude sites which match the filters? One of FLT_INCLUDE/FLT_EXCLUDE
    convert_t *convert;
    bcf_srs_t *files;
    bcf_hdr_t *header;
    int nsamples, *samples, sample_is_file;
    char **argv, *format_str, *sample_list, *targets_list, *regions_list, *vcf_list, *fn_out;
    int argc, list_columns, print_header, allow_undef_tags;
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
    args->convert = convert_init(args->header, samples, nsamples, args->format_str);
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
    free(args->samples);
}

static void query_vcf(args_t *args)
{
    kstring_t str = {0,0,0};

    if ( args->print_header )
    {
        convert_header(args->convert,&str);
        fwrite(str.s, str.l, 1, args->out);
    }

    while ( bcf_sr_next_line(args->files) )
    {
        if ( !bcf_sr_has_line(args->files,0) ) continue;
        bcf1_t *line = args->files->readers[0].buffer[0];
        bcf_unpack(line, args->files->max_unpack);

        if ( args->filter )
        {
            int pass = filter_test(args->filter, line, NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( !pass ) continue;
        }

        str.l = 0;
        convert_line(args->convert, line, &str);
        if ( str.l )
            fwrite(str.s, str.l, 1, args->out);
    }
    if ( str.m ) free(str.s);
}

static void list_columns(args_t *args)
{
    void *has_sample = NULL;
    if ( args->sample_list )
    {
        has_sample = khash_str2int_init();
        int i, nsmpl;
        char **smpl = hts_readlist(args->sample_list, args->sample_is_file, &nsmpl);
        for (i=0; i<nsmpl; i++) khash_str2int_inc(has_sample, smpl[i]);
        free(smpl);
    }

    int i;
    bcf_sr_t *reader = &args->files->readers[0];
    for (i=0; i<bcf_hdr_nsamples(reader->header); i++)
    {
        if ( has_sample && !khash_str2int_has_key(has_sample, reader->header->samples[i]) ) continue;
        printf("%s\n", reader->header->samples[i]);
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
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Extracts fields from VCF/BCF file and prints them in user-defined format\n");
    fprintf(stderr, "Usage:   bcftools query [options] <A.vcf.gz> [<B.vcf.gz> [...]]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -c, --collapse <string>           collapse lines with duplicate positions for <snps|indels|both|all|some|none>, see man page [none]\n");
    fprintf(stderr, "    -e, --exclude <expr>              exclude sites for which the expression is true (see man page for details)\n");
    fprintf(stderr, "    -f, --format <string>             see man page for details\n");
    fprintf(stderr, "    -H, --print-header                print header\n");
    fprintf(stderr, "    -i, --include <expr>              select sites for which the expression is true (see man page for details)\n");
    fprintf(stderr, "    -l, --list-samples                print the list of samples and exit\n");
    fprintf(stderr, "    -o, --output-file <file>          output file name [stdout]\n");
    fprintf(stderr, "    -r, --regions <region>            restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file <file>         restrict to regions listed in a file\n");
    fprintf(stderr, "    -s, --samples <list>              list of samples to include\n");
    fprintf(stderr, "    -S, --samples-file <file>         file of samples to include\n");
    fprintf(stderr, "    -t, --targets <region>            similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "    -T, --targets-file <file>         similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "    -u, --allow-undef-tags            print \".\" for undefined tags\n");
    fprintf(stderr, "    -v, --vcf-list <file>             process multiple VCFs listed in the file\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "\tbcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT[\\t%%SAMPLE=%%GT]\\n' file.vcf.gz\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfquery(int argc, char *argv[])
{
    int c, collapse = 0;
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;
    int regions_is_file = 0, targets_is_file = 0;

    static struct option loptions[] =
    {
        {"help",0,0,'h'},
        {"list-samples",0,0,'l'},
        {"include",1,0,'i'},
        {"exclude",1,0,'e'},
        {"format",1,0,'f'},
        {"output-file",1,0,'o'},
        {"regions",1,0,'r'},
        {"regions-file",1,0,'R'},
        {"targets",1,0,'t'},
        {"targets-file",1,0,'T'},
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
                if ( !strcmp(optarg,"snps") ) collapse |= COLLAPSE_SNPS;
                else if ( !strcmp(optarg,"indels") ) collapse |= COLLAPSE_INDELS;
                else if ( !strcmp(optarg,"both") ) collapse |= COLLAPSE_SNPS | COLLAPSE_INDELS;
                else if ( !strcmp(optarg,"any") ) collapse |= COLLAPSE_ANY;
                else if ( !strcmp(optarg,"all") ) collapse |= COLLAPSE_ANY;
                else if ( !strcmp(optarg,"some") ) collapse |= COLLAPSE_SOME;
                else error("The --collapse string \"%s\" not recognised.\n", optarg);
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
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; targets_is_file = 1; break;
            case 'l': args->list_columns = 1; break;
            case 'u': args->allow_undef_tags = 1; break;
            case 's': args->sample_list = optarg; break;
            case 'S': args->sample_list = optarg; args->sample_is_file = 1; break;
            case 'h':
            case '?': usage();
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
        if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open %s: %s\n", fname,bcf_sr_strerror(args->files->errnum));
        list_columns(args);
        bcf_sr_destroy(args->files);
        free(args);
        return 0;
    }

    if ( !args->format_str ) usage();
    args->out = args->fn_out ? fopen(args->fn_out, "w") : stdout;
    if ( !args->out ) error("%s: %s\n", args->fn_out,strerror(errno));

    if ( !args->vcf_list )
    {
        if ( !fname ) usage();
        args->files = bcf_sr_init();
        args->files->collapse = collapse;
        if ( optind+1 < argc ) args->files->require_index = 1;
        if ( args->regions_list && bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
        if ( args->targets_list )
        {
            if ( bcf_sr_set_targets(args->files, args->targets_list, targets_is_file, 0)<0 )
                error("Failed to read the targets: %s\n", args->targets_list);
        }
        while ( fname )
        {
            if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open %s: %s\n", fname,bcf_sr_strerror(args->files->errnum));
            fname = ++optind < argc ? argv[optind] : NULL;
        }
        init_data(args);
        query_vcf(args);
        free(args->format_str);
        destroy_data(args);
        bcf_sr_destroy(args->files);
        fclose(args->out);
        free(args);
        return 0;
    }

    // multiple VCFs
    int i, k, nfiles, prev_nsamples = 0;
    char **fnames, **prev_samples = NULL;
    fnames = hts_readlist(args->vcf_list, 1, &nfiles);
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
            prev_samples = copy_header(args->header, args->files->readers[0].header->samples, bcf_hdr_nsamples(args->files->readers[0].header));
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
    fclose(args->out);
    destroy_list(fnames, nfiles);
    destroy_list(prev_samples, prev_nsamples);
    free(args->format_str);
    free(args);
    return 0;
}


