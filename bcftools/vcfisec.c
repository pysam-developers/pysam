/*  vcfisec.c -- Create intersections, unions and complements of VCF files.

    Copyright (C) 2012-2014 Genome Research Ltd.

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
    int isec_op, isec_n, *write, iwrite, nwrite, output_type, n_threads;
    int nflt, *flt_logic;
    filter_t **flt;
    char **flt_expr;
    bcf_srs_t *files;
    FILE *fh_log, *fh_sites;
    htsFile **fh_out;
    char **argv, *prefix, *output_fname, **fnames, *write_files, *targets_list, *regions_list;
    char *isec_exact;
    int argc, record_cmd_line;
}
args_t;

/**
 *  mkdir_p() - create new directory for a file $fname
 *  @fname:   the file name to create the directory for, the part after last "/" is ignored
 */
void mkdir_p(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(NULL, 0, fmt, ap) + 2;
    va_end(ap);

    char *path = (char*)malloc(n);
    va_start(ap, fmt);
    vsnprintf(path, n, fmt, ap);
    va_end(ap);

    char *tmp = strdup(path), *p = tmp+1;
    while (*p)
    {
        while (*p && *p!='/') p++;
        if ( !*p ) break;
        char ctmp = *p;
        *p = 0;
        int ret = mkdir(tmp,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if ( ret!=0 && errno!=EEXIST ) error("Error creating directory %s: %s\n", path,strerror(errno));
        *p = ctmp;
        while ( *p && *p=='/' ) p++;
    }
    free(tmp);
    free(path);
}

/**
 *  open_file() - open new file creating the file name using vsnprintf
 *  @fname:  if not NULL, on output will point to newly allocated fname string
 *  @mode:   if NULL, only the file name string will be created
 *  @fmt:    vsnprintf format and args
 *
 *  Returns open file descriptor or NULL if mode is NULL.
 */
FILE *open_file(char **fname, const char *mode, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    int n = vsnprintf(NULL, 0, fmt, ap) + 2;
    va_end(ap);

    char *str = (char*)malloc(n);
    va_start(ap, fmt);
    vsnprintf(str, n, fmt, ap);
    va_end(ap);

    mkdir_p(str);
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
        out_fh = hts_open(args->output_fname? args->output_fname : "-",hts_bcf_wmode(args->output_type));
        if ( out_fh == NULL ) error("Can't write to %s: %s\n", args->output_fname? args->output_fname : "standard output", strerror(errno));
        if ( args->n_threads ) hts_set_threads(out_fh, args->n_threads);
        if (args->record_cmd_line) bcf_hdr_append_version(files->readers[args->iwrite].header,args->argc,args->argv,"bcftools_isec");
        bcf_hdr_write(out_fh, files->readers[args->iwrite].header);
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
            if ( bcf_sr_has_line(files,args->iwrite) )
                bcf_write1(out_fh, files->readers[args->iwrite].header, files->readers[args->iwrite].buffer[0]);
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
            fwrite(str.s,sizeof(char),str.l,args->fh_sites);
        }

        if ( args->prefix )
        {
            if ( args->isec_op==OP_VENN && ret==3 )
            {
                if ( !args->nwrite || args->write[0] )
                    bcf_write1(args->fh_out[2], bcf_sr_get_header(files,0), bcf_sr_get_line(files,0));
                if ( !args->nwrite || args->write[1] )
                    bcf_write1(args->fh_out[3], bcf_sr_get_header(files,1), bcf_sr_get_line(files,1));
            }
            else
            {
                for (i=0; i<files->nreaders; i++)
                {
                    if ( !bcf_sr_has_line(files,i) ) continue;
                    if ( args->write && !args->write[i] ) continue;
                    bcf_write1(args->fh_out[i], files->readers[i].header, files->readers[i].buffer[0]);
                }
            }
        }
    }
    if ( str.s ) free(str.s);
    if ( out_fh ) hts_close(out_fh);
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
        if ( i<0 || i>args->files->nreaders ) error("The index is out of range: %d (%s)\n", i, args->write_files);
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
                args->fh_out[i] = hts_open(args->fnames[i], hts_bcf_wmode(args->output_type));  \
                if ( !args->fh_out[i] ) error("Could not open %s\n", args->fnames[i]); \
                if ( args->n_threads ) hts_set_threads(args->fh_out[i], args->n_threads); \
                if (args->record_cmd_line) bcf_hdr_append_version(args->files->readers[j].header,args->argc,args->argv,"bcftools_isec"); \
                bcf_hdr_write(args->fh_out[i], args->files->readers[j].header); \
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

            args->fh_sites = open_file(NULL, "w", "%s/sites.txt", args->prefix);
            if ( !args->fh_sites ) error("%s/sites.txt: %s\n", args->prefix, strerror(errno));
        }
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
            hts_close(args->fh_out[i]);
            if ( args->output_type==FT_VCF_GZ )
            {
                tbx_conf_t conf = tbx_conf_vcf;
                tbx_index_build(args->fnames[i], -1, &conf);
            }
            else if ( args->output_type==FT_BCF_GZ )
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
    fprintf(stderr, "    -c, --collapse <string>       treat as identical records with <snps|indels|both|all|some|none>, see man page for details [none]\n");
    fprintf(stderr, "    -C, --complement              output positions present only in the first file but missing in the others\n");
    fprintf(stderr, "    -e, --exclude <expr>          exclude sites for which the expression is true\n");
    fprintf(stderr, "    -f, --apply-filters <list>    require at least one of the listed FILTER strings (e.g. \"PASS,.\")\n");
    fprintf(stderr, "    -i, --include <expr>          include only sites for which the expression is true\n");
    fprintf(stderr, "        --no-version                  do not append version and command line to the header\n");
    fprintf(stderr, "    -n, --nfiles [+-=~]<int>      output positions present in this many (=), this many or more (+), this many or fewer (-), the exact (~) files\n");
    fprintf(stderr, "    -o, --output <file>           write output to a file [standard output]\n");
    fprintf(stderr, "    -O, --output-type <b|u|z|v>   b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "    -p, --prefix <dir>            if given, subset each of the input files accordingly, see also -w\n");
    fprintf(stderr, "    -r, --regions <region>        restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file <file>     restrict to regions listed in a file\n");
    fprintf(stderr, "    -t, --targets <region>        similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "    -T, --targets-file <file>     similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "        --threads <int>           number of extra output compression threads [0]\n");
    fprintf(stderr, "    -w, --write <list>            list of files to write with -p given as 1-based indexes. By default, all files are written\n");
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
    int targets_is_file = 0, regions_is_file = 0;

    static struct option loptions[] =
    {
        {"help",no_argument,NULL,'h'},
        {"exclude",required_argument,NULL,'e'},
        {"include",required_argument,NULL,'i'},
        {"collapse",required_argument,NULL,'c'},
        {"complement",no_argument,NULL,'C'},
        {"apply-filters",required_argument,NULL,'f'},
        {"nfiles",required_argument,NULL,'n'},
        {"prefix",required_argument,NULL,'p'},
        {"write",required_argument,NULL,'w'},
        {"targets",required_argument,NULL,'t'},
        {"targets-file",required_argument,NULL,'T'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"threads",required_argument,NULL,9},
        {"no-version",no_argument,NULL,8},
        {NULL,0,NULL,0}
    };
    while ((c = getopt_long(argc, argv, "hc:r:R:p:n:w:t:T:Cf:o:O:i:e:",loptions,NULL)) >= 0) {
        switch (c) {
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
            case 'C': args->isec_op = OP_COMPLEMENT; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;
            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; targets_is_file = 1; break;
            case 'p': args->prefix = optarg; break;
            case 'w': args->write_files = optarg; break;
            case 'i': add_filter(args, optarg, FLT_INCLUDE); break;
            case 'e': add_filter(args, optarg, FLT_EXCLUDE); break;
            case 'n':
                {
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
            case  9 : args->n_threads = strtol(optarg, 0, 0); break;
            case  8 : args->record_cmd_line = 0; break;
            case 'h':
            case '?': usage();
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if ( argc-optind<1 ) usage();   // no file given
    if ( args->targets_list && bcf_sr_set_targets(args->files, args->targets_list, targets_is_file,0)<0 )
        error("Failed to read the targets: %s\n", args->targets_list);
    if ( args->regions_list && bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
        error("Failed to read the regions: %s\n", args->regions_list);
    if ( argc-optind==2 && !args->isec_op )
    {
        args->isec_op = OP_VENN;
        if ( !args->prefix ) error("Expected the -p option\n");
    }
    if ( !args->targets_list )
    {
        if ( argc-optind<2  ) error("Expected multiple files or the --targets option\n");
        if ( !args->isec_op ) error("Expected two file names or one of the options --complement, --nfiles or --targets\n");
    }
    args->files->require_index = 1;
    while (optind<argc)
    {
        if ( !bcf_sr_add_reader(args->files, argv[optind]) ) error("Failed to open %s: %s\n", argv[optind],bcf_sr_strerror(args->files->errnum));
        optind++;
    }
    init_data(args);
    isec_vcf(args);
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}

