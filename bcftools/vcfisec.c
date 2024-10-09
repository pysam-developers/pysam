/*  vcfisec.c -- Create intersections, unions and complements of VCF files.

    Copyright (C) 2012-2023 Genome Research Ltd.

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
#include <htslib/vcfutils.h>
#include <htslib/hts_os.h>
#include <htslib/hts_defs.h>
#include "bcftools.h"
#include "filter.h"

#define OP_PLUS 1
#define OP_MINUS 2
#define OP_EQUAL 3
#define OP_VENN 4
#define OP_COMPLEMENT 5
#define OP_EXACT 6

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef struct
{
    int isec_op, isec_n, *write, iwrite, nwrite, output_type, n_threads, clevel;
    int nflt, *flt_logic;
    filter_t **flt;
    char **flt_expr;
    bcf_srs_t *files;
    FILE *fh_log, *fh_sites;
    htsFile **fh_out;
    char **argv, *prefix, *output_fname, **fnames, *write_files, *targets_list, *regions_list;
    char *isec_exact, *file_list;
    int argc, record_cmd_line;
    char *index_fn;
    int write_index;
}
args_t;

/**
 *  mkdir_p() - create new directory for a file $fname
 *  @fname:   the file name to create the directory for, the part after last "/" is ignored
 */
void HTS_FORMAT(HTS_PRINTF_FMT, 1, 2)
mkdir_p(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(NULL, 0, fmt, ap) + 2;
    va_end(ap);

    char *tmp = (char*)malloc(n);
    if (!tmp) error("Couldn't allocate space for path: %s\n", strerror(errno));
    va_start(ap, fmt);
    vsnprintf(tmp, n, fmt, ap);
    va_end(ap);

    char *p = tmp+1;
    while (*p)
    {
        while (*p && *p!='/') p++;
        if ( !*p ) break;
        char ctmp = *p;
        *p = 0;
        int ret = mkdir(tmp,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if ( ret!=0 && errno!=EEXIST ) error("Error creating directory %s: %s\n", tmp,strerror(errno));
        *p = ctmp;
        while ( *p && *p=='/' ) p++;
    }
    free(tmp);
}

/**
 *  open_file() - open new file creating the file name using vsnprintf
 *  @fname:  if not NULL, on output will point to newly allocated fname string
 *  @mode:   if NULL, only the file name string will be created
 *  @fmt:    vsnprintf format and args
 *
 *  Returns open file descriptor or NULL if mode is NULL.
 */
FILE * HTS_FORMAT(HTS_PRINTF_FMT, 3, 4)
open_file(char **fname, const char *mode, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(NULL, 0, fmt, ap) + 2;
    va_end(ap);

    char *str = (char*)malloc(n);
    va_start(ap, fmt);
    vsnprintf(str, n, fmt, ap);
    va_end(ap);

    mkdir_p("%s", str);
    if ( !mode )
    {
        if ( !fname ) error("Uh: expected fname or mode\n");
        *fname = str;
        return NULL;
    }

    FILE *fp = fopen(str,mode);
    if ( fname ) *fname = str;
    else free(str);
    return fp;
}

void isec_vcf(args_t *args)
{
    bcf_srs_t *files = args->files;
    kstring_t str = {0,0,0};
    htsFile *out_fh = NULL;

    // When only one VCF is output, print VCF to stdout or -o file
    int out_std = 0;
    if ( args->nwrite==1 && !args->prefix ) out_std = 1;
    if ( args->targets_list && files->nreaders==1 ) out_std = 1;
    if ( out_std )
    {
        char wmode[8];
        set_wmode(wmode,args->output_type,args->output_fname,args->clevel);
        out_fh = hts_open(args->output_fname ? args->output_fname : "-", wmode);
        if ( out_fh == NULL ) error("Can't write to %s: %s\n", args->output_fname? args->output_fname : "standard output", strerror(errno));
        if ( args->n_threads ) hts_set_threads(out_fh, args->n_threads);
        if (args->record_cmd_line) bcf_hdr_append_version(files->readers[args->iwrite].header,args->argc,args->argv,"bcftools_isec");
        if ( bcf_hdr_write(out_fh, files->readers[args->iwrite].header)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname?args->output_fname:"standard output");
        if ( init_index2(out_fh,files->readers[args->iwrite].header,
                         args->output_fname,&args->index_fn,
                         args->write_index)<0 )
            error("Error: failed to initialise index for %s\n",
                  args->output_fname?args->output_fname:"standard output");
    }
    if ( !args->nwrite && !out_std && !args->prefix )
        fprintf(stderr,"Note: -w option not given, printing list of sites...\n");

    int n;
    while ( (n=bcf_sr_next_line(files)) )
    {
        bcf_sr_t *reader = NULL;
        bcf1_t *line = NULL;
        int i, ret = 0;
        for (i=0; i<files->nreaders; i++)
        {
            if ( !bcf_sr_has_line(files,i) ) continue;

            if ( args->nflt && args->flt[i] )
            {
                bcf1_t *rec = bcf_sr_get_line(files, i);
                int pass = filter_test(args->flt[i], rec, NULL);
                if ( args->flt_logic[i] & FLT_EXCLUDE ) pass = pass ? 0 : 1;
                if ( !pass )
                {
                    files->has_line[i] = 0;
                    n--;
                    continue;
                }
            }

            if ( !line )
            {
                line = files->readers[i].buffer[0];
                reader = &files->readers[i];
            }
            ret |= 1<<i;    // this may overflow for many files, but will be used only with two (OP_VENN)
        }
        if ( !line ) continue;  // the site has been filtered in all files

        switch (args->isec_op)
        {
            case OP_COMPLEMENT: if ( n!=1 || !bcf_sr_has_line(files,0) ) continue; break;
            case OP_EQUAL: if ( n != args->isec_n ) continue; break;
            case OP_PLUS: if ( n < args->isec_n ) continue; break;
            case OP_MINUS: if ( n > args->isec_n ) continue; break;
            case OP_EXACT:
                for (i=0; i<files->nreaders; i++)
                    if ( files->has_line[i] != args->isec_exact[i] ) break;
                if ( i<files->nreaders ) continue;
                break;
        }

        if ( out_std )
        {
            if ( bcf_sr_has_line(files,args->iwrite) && bcf_write1(out_fh, files->readers[args->iwrite].header, files->readers[args->iwrite].buffer[0])!=0 )
                error("[%s] Error: cannot write to %s\n", __func__, args->output_fname ? args->output_fname : "standard output");
            continue;
        }
        else if ( args->fh_sites )
        {
            str.l = 0;
            kputs(reader->header->id[BCF_DT_CTG][line->rid].key, &str); kputc('\t', &str);
            kputw(line->pos+1, &str); kputc('\t', &str);
            if (line->n_allele > 0) kputs(line->d.allele[0], &str);
            else kputc('.', &str);
            kputc('\t', &str);
            if (line->n_allele > 1) kputs(line->d.allele[1], &str);
            else kputc('.', &str);
            for (i=2; i<line->n_allele; i++)
            {
                kputc(',', &str);
                kputs(line->d.allele[i], &str);
            }
            kputc('\t', &str);
            for (i=0; i<files->nreaders; i++)
                kputc(bcf_sr_has_line(files,i)?'1':'0', &str);
            kputc('\n', &str);
            if ( fwrite(str.s,sizeof(char),str.l,args->fh_sites)!=str.l )
                error("[%s] Error: failed to write %d bytes to %s\n", __func__,(int)str.l,args->output_fname ? args->output_fname : "standard output");
        }

        if ( args->prefix )
        {
            if ( args->isec_op==OP_VENN && ret==3 )
            {
                if ( !args->nwrite || args->write[0] )
                {
                    if ( bcf_write1(args->fh_out[2], bcf_sr_get_header(files,0), bcf_sr_get_line(files,0))!=0 )
                         error("[%s] Error: cannot write\n", __func__);
                }
                if ( !args->nwrite || args->write[1] )
                {
                    if ( bcf_write1(args->fh_out[3], bcf_sr_get_header(files,1), bcf_sr_get_line(files,1))!=0 )
                        error("[%s] Error: cannot write\n", __func__);
                }
            }
            else
            {
                for (i=0; i<files->nreaders; i++)
                {
                    if ( !bcf_sr_has_line(files,i) ) continue;
                    if ( args->write && !args->write[i] ) continue;
                    if ( bcf_write1(args->fh_out[i], files->readers[i].header, files->readers[i].buffer[0])!=0 ) error("[%s] Error: cannot write\n", __func__);
                }
            }
        }
    }
    if ( str.s ) free(str.s);
    if ( out_fh )
    {
        if ( args->write_index )
        {
            if ( bcf_idx_save(out_fh)<0 )
            {
                if ( hts_close(out_fh)!=0 ) error("Error: close failed .. %s\n", args->output_fname?args->output_fname:"stdout");
                error("Error: cannot write to index %s\n", args->index_fn);
            }
            free(args->index_fn);
        }
        if ( hts_close(out_fh)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->output_fname? args->output_fname : "-");
    }
}

static void add_filter(args_t *args, char *expr, int logic)
{
    args->nflt++;
    args->flt_expr = (char**) realloc(args->flt_expr,sizeof(char*)*args->nflt);
    args->flt_logic = (int*) realloc(args->flt_logic,sizeof(int)*args->nflt);
    args->flt = (filter_t**) realloc(args->flt,sizeof(filter_t*)*args->nflt);
    if ( expr[0]=='-' && expr[1]==0 )
    {
        args->flt_expr[args->nflt-1] = NULL;
        args->flt[args->nflt-1] = NULL;
    }
    else
        args->flt_expr[args->nflt-1] = expr;
    args->flt_logic[args->nflt-1] = logic;
}

static void destroy_data(args_t *args);
static void init_data(args_t *args)
{
    int i;
    if ( args->nflt )
    {
        if ( args->nflt > 1 && args->nflt!=args->files->nreaders )
            error("Error: expected either one -i/-e option or as many as there are input files\n");
        if ( args->nflt < args->files->nreaders )
        {
            if ( !args->flt_expr[0] ) error("Error: useless use of -i/-e\n");
            args->nflt = args->files->nreaders;
            args->flt_expr = (char**) realloc(args->flt_expr,sizeof(char*)*args->nflt);
            args->flt_logic = (int*) realloc(args->flt_logic,sizeof(int)*args->nflt);
            args->flt = (filter_t**) realloc(args->flt,sizeof(filter_t*)*args->nflt);
            for (i=1; i<args->nflt; i++)
            {
                args->flt_expr[i]  = args->flt_expr[0];
                args->flt_logic[i] = args->flt_logic[0];
                args->flt[i] = filter_init(args->files->readers[i].header,args->flt_expr[i]);
            }
            args->flt[0] = filter_init(args->files->readers[0].header,args->flt_expr[0]);
        }
        else
        {
            for (i=0; i<args->files->nreaders; i++)
            {
                if ( !args->flt_expr[i] ) continue;
                args->flt[i] = filter_init(args->files->readers[i].header,args->flt_expr[i]);
            }
        }
    }

    if ( args->isec_op==OP_EXACT )
    {
        if ( strlen(args->isec_exact)!=args->files->nreaders )
            error("The number of files does not match the bitmask: %d vs %s\n", args->files->nreaders,args->isec_exact);
        for (i=0; i<args->files->nreaders; i++)
            if ( args->isec_exact[i]!='0' && args->isec_exact[i]!='1' ) error("Unexpected bitmask: %s\n",args->isec_exact);
        for (i=0; i<args->files->nreaders; i++)
            args->isec_exact[i] -= '0';
    }

    // Which files to write: parse the string passed with -w
    char *p = args->write_files;
    while (p && *p)
    {
        if ( !args->write ) args->write = (int*) calloc(args->files->nreaders,sizeof(int));
        if ( sscanf(p,"%d",&i)!=1 ) error("Could not parse --write %s\n", args->write_files);
        if ( i<=0 || i>args->files->nreaders ) error("The index is out of range: %d (-w %s)\n", i, args->write_files);
        args->write[i-1] = 1;
        args->iwrite = i-1;
        args->nwrite++;
        while (*p && *p!=',') p++;
        if ( *p==',' ) p++;
    }
    if ( args->nwrite>1 && !args->prefix ) error("Expected -p when multiple output files given: --write %s\n", args->write_files);
    if ( args->isec_op==OP_COMPLEMENT && args->nwrite )
    {
        if ( args->nwrite>1 ) error("Multiple files to -w make no sense with -C\n");
        if ( !args->write[0] ) error("Only -w1 makes sense with -C\n");
    }

    if ( args->prefix )
    {
        // Init output directory and create the readme file
        args->fh_log = open_file(NULL,"w","%s/README.txt", args->prefix);
        if ( !args->fh_log ) error("%s/README.txt: %s\n", args->prefix, strerror(errno));

        fprintf(args->fh_log,"This file was produced by vcfisec.\n");
        fprintf(args->fh_log,"The command line was:\tbcftools %s ", args->argv[0]);
        int i;
        for (i=1; i<args->argc; i++) fprintf(args->fh_log," %s",args->argv[i]);
        fprintf(args->fh_log,"\n\nUsing the following file names:\n");

        const char *suffix = "vcf";
        if ( args->output_type & FT_BCF ) suffix = "bcf";
        else if ( args->output_type & FT_GZ ) suffix = "vcf.gz";

        // Open output files and write the legend
        if ( args->isec_op==OP_VENN )
        {
            args->fh_out = (htsFile**) malloc(sizeof(htsFile*)*4);
            args->fnames = (char**) calloc(4,sizeof(char*));

            #define OPEN_FILE(i,j) { \
                open_file(&args->fnames[i], NULL, "%s/%04d.%s", args->prefix, i, suffix); \
                char wmode[8]; \
                set_wmode(wmode,args->output_type,args->fnames[i],args->clevel); \
                args->fh_out[i] = hts_open(args->fnames[i], wmode); \
                if ( !args->fh_out[i] ) error("Could not open %s\n", args->fnames[i]); \
                if ( args->n_threads ) hts_set_threads(args->fh_out[i], args->n_threads); \
                if (args->record_cmd_line) bcf_hdr_append_version(args->files->readers[j].header,args->argc,args->argv,"bcftools_isec"); \
                if ( bcf_hdr_write(args->fh_out[i], args->files->readers[j].header)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->fnames[i]); \
            }
            if ( !args->nwrite || args->write[0] )
            {
                OPEN_FILE(0,0);
                fprintf(args->fh_log,"%s\tfor records private to\t%s\n", args->fnames[0], args->files->readers[0].fname);
            }
            if ( !args->nwrite || args->write[1] )
            {
                OPEN_FILE(1,1);
                fprintf(args->fh_log,"%s\tfor records private to\t%s\n", args->fnames[1], args->files->readers[1].fname);
            }
            if ( !args->nwrite || args->write[0] )
            {
                OPEN_FILE(2,0);
                fprintf(args->fh_log,"%s\tfor records from %s shared by both\t%s %s\n", args->fnames[2], args->files->readers[0].fname, args->files->readers[0].fname, args->files->readers[1].fname);
            }
            if ( !args->nwrite || args->write[1] )
            {
                OPEN_FILE(3,1);
                fprintf(args->fh_log,"%s\tfor records from %s shared by both\t%s %s\n", args->fnames[3], args->files->readers[1].fname, args->files->readers[0].fname, args->files->readers[1].fname);
            }
        }
        else
        {
            // Init one output file for each reader
            args->fh_out = (htsFile**) calloc(args->files->nreaders, sizeof(htsFile*));
            args->fnames = (char**) calloc(args->files->nreaders, sizeof(char*));

            for (i=0; i<args->files->nreaders; i++)
            {
                if ( args->write && !args->write[i] ) continue;
                if ( args->isec_op==OP_COMPLEMENT && i>0 ) break;
                OPEN_FILE(i,i);
                fprintf(args->fh_log,"%s\tfor stripped\t%s\n", args->fnames[i], args->files->readers[i].fname);
            }
            #undef OPEN_FILE
        }
        args->fh_sites = open_file(NULL, "w", "%s/sites.txt", args->prefix);
        if ( !args->fh_sites ) error("%s/sites.txt: %s\n", args->prefix, strerror(errno));
    }
    else {
        if (args->output_fname) {
            args->fh_sites = fopen(args->output_fname, "w");
            if ( args->fh_sites == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
        }
        else
            args->fh_sites = stdout;
    }
}

static void destroy_data(args_t *args)
{
    int i;
    if ( args->nflt )
    {
        for (i=0; i<args->nflt; i++)
        {
            if ( !args->flt[i] ) continue;
            filter_destroy(args->flt[i]);
        }
        free(args->flt_expr);
        free(args->flt);
        free(args->flt_logic);
    }
    if ( args->prefix )
    {
        fclose(args->fh_log);
        int n = args->isec_op==OP_VENN ? 4 : args->files->nreaders;
        for (i=0; i<n; i++)
        {
            if ( !args->fnames[i] ) continue;
            if ( hts_close(args->fh_out[i])!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->fnames[i]);
            int is_tbi = !args->write_index 
                      || (args->write_index&127) == HTS_FMT_TBI;
            if ( args->output_type==FT_VCF_GZ && is_tbi )
            {
                tbx_conf_t conf = tbx_conf_vcf;
                tbx_index_build(args->fnames[i], -1, &conf);
            }
            else if ( args->output_type==FT_BCF_GZ || !is_tbi )
            {
                if ( bcf_index_build(args->fnames[i],14) ) error("Could not index %s\n", args->fnames[i]);
            }
            free(args->fnames[i]);
        }
        free(args->fh_out);
        free(args->fnames);
        if ( args->fh_sites ) fclose(args->fh_sites);
        if ( args->write ) free(args->write);
    }
}

static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Create intersections, unions and complements of VCF files.\n");
    fprintf(stderr, "Usage:   bcftools isec [options] <A.vcf.gz> <B.vcf.gz> [...]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -c, --collapse STRING          Treat as identical records with <snps|indels|both|all|some|none>, see man page for details [none]\n");
    fprintf(stderr, "    -C, --complement               Output positions present only in the first file but missing in the others\n");
    fprintf(stderr, "    -e, --exclude EXPR             Exclude sites for which the expression is true\n");
    fprintf(stderr, "    -f, --apply-filters LIST       Require at least one of the listed FILTER strings (e.g. \"PASS,.\")\n");
    fprintf(stderr, "    -i, --include EXPR             Include only sites for which the expression is true\n");
    fprintf(stderr, "    -l, --file-list FILE           Read the input file names from the file\n");
    fprintf(stderr, "        --no-version               Do not append version and command line to the header\n");
    fprintf(stderr, "    -n, --nfiles [+-=~]INT         Output positions present in this many (=), this many or more (+), this many or fewer (-), the exact (~) files\n");
    fprintf(stderr, "    -o, --output FILE              Write output to a file [standard output]\n");
    fprintf(stderr, "    -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n");
    fprintf(stderr, "    -p, --prefix DIR               If given, subset each of the input files accordingly, see also -w\n");
    fprintf(stderr, "    -r, --regions REGION           Restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file FILE        Restrict to regions listed in a file\n");
    fprintf(stderr, "        --regions-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n");
    fprintf(stderr, "    -t, --targets REGION           Similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "    -T, --targets-file FILE        Similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "        --targets-overlap 0|1|2    Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n");
    fprintf(stderr, "        --threads INT              Use multithreading with <int> worker threads [0]\n");
    fprintf(stderr, "    -w, --write LIST               List of files to write with -p given as 1-based indexes. By default, all files are written\n");
    fprintf(stderr, "    -W, --write-index[=FMT]        Automatically index the output files [off]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "   # Create intersection and complements of two sets saving the output in dir/*\n");
    fprintf(stderr, "   bcftools isec A.vcf.gz B.vcf.gz -p dir\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   # Filter sites in A and B (but not in C) and create intersection\n");
    fprintf(stderr, "   bcftools isec -e'MAF<0.01' -i'dbSNP=1' -e - A.vcf.gz B.vcf.gz C.vcf.gz -p dir\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   # Extract and write records from A shared by both A and B using exact allele match\n");
    fprintf(stderr, "   bcftools isec A.vcf.gz B.vcf.gz -p dir -n =2 -w 1\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   # Extract and write records from C found in A and C but not in B\n");
    fprintf(stderr, "   bcftools isec A.vcf.gz B.vcf.gz C.vcf.gz -p dir -n~101 -w 3\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "   # Extract records private to A or B comparing by position only\n");
    fprintf(stderr, "   bcftools isec A.vcf.gz B.vcf.gz -p dir -n -1 -c all\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfisec(int argc, char *argv[])
{
    int c;
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->files  = bcf_sr_init();
    args->argc   = argc; args->argv = argv;
    args->output_fname = NULL;
    args->output_type = FT_VCF;
    args->n_threads = 0;
    args->record_cmd_line = 1;
    args->clevel = -1;
    int targets_is_file = 0, regions_is_file = 0;
    int regions_overlap = 1;
    int targets_overlap = 0;

    static struct option loptions[] =
    {
        {"help",no_argument,NULL,'h'},
        {"exclude",required_argument,NULL,'e'},
        {"include",required_argument,NULL,'i'},
        {"collapse",required_argument,NULL,'c'},
        {"complement",no_argument,NULL,'C'},
        {"apply-filters",required_argument,NULL,'f'},
        {"file-list",required_argument,NULL,'l'},
        {"nfiles",required_argument,NULL,'n'},
        {"prefix",required_argument,NULL,'p'},
        {"write",required_argument,NULL,'w'},
        {"targets",required_argument,NULL,'t'},
        {"targets-file",required_argument,NULL,'T'},
        {"targets-overlap",required_argument,NULL,4},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"regions-overlap",required_argument,NULL,3},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"threads",required_argument,NULL,9},
        {"no-version",no_argument,NULL,8},
        {"write-index",optional_argument,NULL,'W'},
        {NULL,0,NULL,0}
    };
    char *tmp;
    while ((c = getopt_long(argc, argv, "hc:r:R:p:n:w:t:T:Cf:o:O:i:e:l:W::",loptions,NULL)) >= 0) {
        switch (c) {
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
            case 'c':
                if ( !strcmp(optarg,"snps") ) args->files->collapse |= COLLAPSE_SNPS;
                else if ( !strcmp(optarg,"indels") ) args->files->collapse |= COLLAPSE_INDELS;
                else if ( !strcmp(optarg,"both") ) args->files->collapse |= COLLAPSE_SNPS | COLLAPSE_INDELS;
                else if ( !strcmp(optarg,"any") ) args->files->collapse |= COLLAPSE_ANY;
                else if ( !strcmp(optarg,"all") ) args->files->collapse |= COLLAPSE_ANY;
                else if ( !strcmp(optarg,"some") ) args->files->collapse |= COLLAPSE_SOME;
                else if ( !strcmp(optarg,"none") ) args->files->collapse = COLLAPSE_NONE;
                else error("The --collapse string \"%s\" not recognised.\n", optarg);
                break;
            case 'f': args->files->apply_filters = optarg; break;
            case 'C':
                if ( args->isec_op!=0 && args->isec_op!=OP_COMPLEMENT ) error("Error: either -C or -n should be given, not both.\n");
                args->isec_op = OP_COMPLEMENT; break;
            case 'l': args->file_list = optarg; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; targets_is_file = 1; break;
            case 'p': args->prefix = optarg; break;
            case 'w':
                if ( args->write_files ) error("The option -w accepts a list of indices and can be given only once\n");
                args->write_files = optarg;
                break;
            case 'i': add_filter(args, optarg, FLT_INCLUDE); break;
            case 'e': add_filter(args, optarg, FLT_EXCLUDE); break;
            case 'n':
                {
                    if ( args->isec_op!=0 && args->isec_op==OP_COMPLEMENT ) error("Error: either -C or -n should be given, not both.\n");
                    if ( args->isec_op!=0 ) error("Error: -n should be given only once.\n");
                    char *p = optarg;
                    if ( *p=='-' ) { args->isec_op = OP_MINUS; p++; }
                    else if ( *p=='+' ) { args->isec_op = OP_PLUS; p++; }
                    else if ( *p=='=' ) { args->isec_op = OP_EQUAL; p++; }
                    else if ( *p=='~' ) { args->isec_op = OP_EXACT; p++; }
                    else if ( isdigit(*p) ) args->isec_op = OP_EQUAL;
                    else error("Could not parse --nfiles %s\n", optarg);
                    if ( args->isec_op == OP_EXACT ) args->isec_exact = p;
                    else if ( sscanf(p,"%d",&args->isec_n)!=1 ) error("Could not parse --nfiles %s\n", optarg);
                }
                break;
            case  3 :
                regions_overlap = parse_overlap_option(optarg);
                if ( regions_overlap < 0 ) error("Could not parse: --regions-overlap %s\n",optarg);
                break;
            case  4 :
                targets_overlap = parse_overlap_option(optarg);
                if ( targets_overlap < 0 ) error("Could not parse: --targets-overlap %s\n",optarg);
                break;
            case  9 : args->n_threads = strtol(optarg, 0, 0); break;
            case  8 : args->record_cmd_line = 0; break;
            case 'W':
                if (!(args->write_index = write_index_parse(optarg)))
                    error("Unsupported index format '%s'\n", optarg);
                break;
            case 'h':
            case '?': usage(); break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if ( argc-optind<1 && !args->file_list ) usage();   // no file given

    int nfiles = 0,i;
    char **files = NULL;
    if ( args->file_list )
    {
        files = hts_readlines(args->file_list, &nfiles);
        if ( !files ) error("Failed to read from %s\n", args->file_list);
    }
    if ( optind<argc )
    {
        int n = argc - optind;
        files = (char**)realloc(files,sizeof(*files)*(n+nfiles));
        for (i=nfiles; i>0; i--) files[n+i-1] = files[n+i-2];
        for (i=0; i<n; i++) files[i] = strdup(argv[optind+i]);
        nfiles += n;
    }

    if ( args->targets_list )
    {
        bcf_sr_set_opt(args->files,BCF_SR_TARGETS_OVERLAP,targets_overlap);
        if ( bcf_sr_set_targets(args->files, args->targets_list, targets_is_file,0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
    }
    if ( args->regions_list )
    {
        bcf_sr_set_opt(args->files,BCF_SR_REGIONS_OVERLAP,regions_overlap);
        if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    if ( nfiles==2 && !args->isec_op )
    {
        args->isec_op = OP_VENN;
        if ( !args->prefix ) error("Expected the -p option\n");
    }
    if ( !args->isec_op )
    {
        args->isec_op = OP_PLUS;
        args->isec_n  = 1;
    }
    args->files->require_index = 1;
    for (i=0; i<nfiles; i++)
    {
        if ( !bcf_sr_add_reader(args->files, files[i]) ) error("Failed to open %s: %s\n", files[i],bcf_sr_strerror(args->files->errnum));
        free(files[i]);
    }
    free(files);

    init_data(args);
    isec_vcf(args);
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}

