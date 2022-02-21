/*  vcfplugin.c -- plugin modules for operating on VCF/BCF files.

    Copyright (C) 2013-2021 Genome Research Ltd.

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

#include "config.h"
#include <stdio.h>
#include <strings.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kseq.h>
#include <htslib/khash_str2int.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif
#include "bcftools.h"
#include "vcmp.h"
#include "filter.h"

#ifdef ENABLE_BCF_PLUGINS

typedef struct _plugin_t plugin_t;

/**
 *   Plugin API:
 *   ----------
 *   const char *about(void)
 *      - short description used by 'bcftools plugin -lv'
 *
 *   const char *usage(void)
 *      - longer description used by 'bcftools +name -h'
 *
 *   int run(int argc, char **argv)
 *      - if implemented, the control is immediately handed over to the plugin,
 *      none of the init/process/destroy functions is called.  Return 0 on
 *      success or non-zero value on error.
 *
 *   int init(int argc, char **argv, bcf_hdr_t *in_hdr, bcf_hdr_t *out_hdr)
 *      - called once at startup, allows to initialize local variables.
 *      Return 1 to suppress normal VCF/BCF header output, -1 on critical
 *      errors, 0 otherwise.
 *
 *   bcf1_t *process(bcf1_t *rec)
 *      - called for each VCF record, return NULL for no output
 *
 *   void destroy(void)
 *      - called after all lines have been processed to clean up
 */
typedef void (*dl_version_f) (const char **, const char **);
typedef int (*dl_run_f) (int, char **);
typedef int (*dl_init_f) (int, char **, bcf_hdr_t *, bcf_hdr_t *);
typedef char* (*dl_about_f) (void);
typedef char* (*dl_usage_f) (void);
typedef bcf1_t* (*dl_process_f) (bcf1_t *);
typedef void (*dl_destroy_f) (void);

struct _plugin_t
{
    int argc;
    char *name, **argv;
    dl_version_f version;
    dl_run_f run;
    dl_init_f init;
    dl_about_f about;
    dl_usage_f usage;
    dl_process_f process;
    dl_destroy_f destroy;
    void *handle;
};


struct _args_t;

typedef struct _rm_tag_t
{
    char *key;
    int hdr_id;
    void (*handler)(struct _args_t *, bcf1_t *, struct _rm_tag_t *);
}
rm_tag_t;

typedef struct
{
    char **cols;
    int ncols, mcols;
    char **als;
    int nals, mals;
    kstring_t line;
    int rid, start, end;
}
annot_line_t;

typedef struct _annot_col_t
{
    int icol, replace;
    char *hdr_key;
    int (*setter)(struct _args_t *, bcf1_t *, struct _annot_col_t *, void*);
}
annot_col_t;

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

typedef struct _args_t
{
    bcf_srs_t *files;
    bcf_hdr_t *hdr, *hdr_out;
    htsFile *out_fh;
    int output_type, n_threads, clevel;

    filter_t *filter;
    char *filter_str;
    int filter_logic;   // include or exclude sites which match the filters? One of FLT_INCLUDE/FLT_EXCLUDE

    plugin_t plugin;
    int nplugin_paths;
    char **plugin_paths;

    char **argv, *output_fname, *regions_list, *targets_list;
    int argc, drop_header, verbose, record_cmd_line, plist_only;
}
args_t;

char *msprintf(const char *fmt, ...);

static void add_plugin_paths(args_t *args, const char *path)
{
    while (1)
    {
        size_t len = strcspn(path, HTS_PATH_SEPARATOR_STR);

        if ( len == 0 )
        {
#ifdef PLUGINPATH
            add_plugin_paths(args, PLUGINPATH);
#endif
        }
        else
        {
            char *dir = (char *) malloc(len + 1);
            strncpy(dir, path, len);
            dir[len] = '\0';

            struct stat st;
            if ( stat(dir, &st) == 0 )
            {
                args->plugin_paths = (char**) realloc(args->plugin_paths,sizeof(char*)*(args->nplugin_paths+1));
                args->plugin_paths[args->nplugin_paths] = dir;
                args->nplugin_paths++;
                if ( args->verbose > 1 && strcmp(".",dir) ) fprintf(stderr, "plugin directory %s .. ok\n", dir);
            }
            else
            {
                if ( args->verbose > 1 ) fprintf(stderr, "plugin directory %s .. %s\n", dir, strerror(errno));
                free(dir);
            }

        }

        path += len;
        if ( *path == HTS_PATH_SEPARATOR_CHAR ) path++;
        else break;
    }
}

static void init_plugin_paths(args_t *args)
{
    if ( args->nplugin_paths!=-1 ) return;

    args->nplugin_paths = 0;
    args->plugin_paths = NULL;

    char *path = getenv("BCFTOOLS_PLUGINS");
    add_plugin_paths(args, path ? path : "");
}

static void *dlopen_plugin(args_t *args, const char *fname)
{
    init_plugin_paths(args);

    void *handle;
    char *tmp;
    int is_absolute_path = 0;
#ifdef _WIN32
    // Windows accepts both forward slash (/) and backslash (\) as folder separator
    // and can have any path prefixed by the drive letter and a colon (:).
    if ( fname[0]=='/' || fname[0]=='\\') is_absolute_path = 1;
    else if ( fname[0] && fname[1]==':' && (fname[2]=='/' || fname[2]=='\\') ) is_absolute_path = 1;
#else
    if ( fname[0]=='/' ) is_absolute_path = 1;
#endif

    kstring_t err = {0,0,0};
    if ( !is_absolute_path )
    {
        int i;
        for (i=0; i<args->nplugin_paths; i++)
        {
            tmp = msprintf("%s/%s%s", args->plugin_paths[i], fname, PLUGIN_EXT);
#ifdef _WIN32
            handle = LoadLibraryA(tmp);
#else
            handle = dlopen(tmp, RTLD_NOW); // valgrind complains about unfreed memory, not our problem though
#endif
            if ( !handle )
#ifdef _WIN32
                ksprintf(&err,"LoadLibraryA   .. %lu\n", GetLastError());
#else
                ksprintf(&err,"%s:\n\tdlopen   .. %s\n", tmp,dlerror());
#endif
            else if ( args->verbose > 1 )
                fprintf(stderr,"%s:\n\tplugin open   .. ok\n", tmp);
            free(tmp);
            if ( handle ) return handle;
        }
    }

#ifdef _WIN32
    handle = LoadLibraryA(fname);
#else
    handle = dlopen(fname, RTLD_NOW);
#endif
    if ( !handle )
#ifdef _WIN32
        ksprintf(&err,"LoadLibraryA   .. %lu\n", GetLastError());
#else
        ksprintf(&err,"%s:\n\tdlopen   .. %s\n", fname,dlerror());
#endif
    else if ( args->verbose > 1 )
        fprintf(stderr,"%s:\n\tplugin open   .. ok\n", fname);

    if ( !handle && (!args->plist_only || args->verbose>1) )
        fprintf(stderr,"%s",err.s);
    free(err.s);

    return handle;
}

static void print_plugin_usage_hint(const char *name)
{
    if ( name )
        fprintf(stderr, "\nThe bcftools plugin \"%s\" was not found or is not functional", name);
    else
        fprintf(stderr, "\nNo functional bcftools plugins were found");
    if ( !getenv("BCFTOOLS_PLUGINS") )
    {
        fprintf(stderr,". The environment variable BCFTOOLS_PLUGINS is not set");
#ifdef PLUGINPATH
        fprintf(stderr,"\nand no usable plugins were found in %s", PLUGINPATH);
#endif
        fprintf(stderr,".\n\n");
    }
    else
    {
        fprintf(stderr,
                " in\n\tBCFTOOLS_PLUGINS=\"%s\".\n\n"
                "- Is the plugin path correct?\n\n"
                "- Run \"bcftools plugin -l\" or \"bcftools plugin -lvv\" for a list of available plugins.\n"
                "\n",
                getenv("BCFTOOLS_PLUGINS")
               );
    }
}

static int load_plugin(args_t *args, const char *fname, int exit_on_error, plugin_t *plugin)
{
    plugin->name = strdup(fname);

    plugin->handle = dlopen_plugin(args, fname);
    if ( !plugin->handle )
    {
        if ( exit_on_error )
        {
            print_plugin_usage_hint(fname);
            error("Could not load \"%s\".\n\n", fname);
        }
        return -1;
    }

#ifdef _WIN32
    plugin->init = (dl_init_f) GetProcAddress(plugin->handle, "init");
    if ( plugin->init && args->verbose > 1 ) fprintf(stderr,"\tinit     .. ok\n");

    plugin->run = (dl_run_f) GetProcAddress(plugin->handle, "run");
    if ( plugin->run && args->verbose > 1 ) fprintf(stderr,"\trun     .. ok\n");

    if ( !plugin->init && !plugin->run )
    {
        if ( exit_on_error ) error("Could not initialize %s, neither run or init found \n", plugin->name);
        else if ( args->verbose > 1 ) fprintf(stderr,"\tinit/run .. not found\n");
        return -1;
    }

    plugin->version = (dl_version_f) GetProcAddress(plugin->handle, "version");
    if ( !plugin->version )
    {
        if ( exit_on_error ) error("Could not initialize %s: version string not found\n", plugin->name);
        else if ( args->verbose > 1 ) fprintf(stderr,"\tversion  .. not found\n");
        return -1;
    }

    plugin->about = (dl_about_f) GetProcAddress(plugin->handle, "about");
    if ( !plugin->about )
    {
        if ( exit_on_error ) error("Could not initialize %s: about string not found\n", plugin->name);
        return -1;
    }

    plugin->usage = (dl_about_f) GetProcAddress(plugin->handle, "usage");
    if ( !plugin->usage )
        plugin->usage = plugin->about;

    if ( plugin->run ) return 0;

    plugin->process = (dl_process_f) GetProcAddress(plugin->handle, "process");
    if ( !plugin->process )
    {
        if ( exit_on_error ) error("Could not initialize %s: process method not found\n", plugin->name);
        return -1;
    }

    plugin->destroy = (dl_destroy_f) GetProcAddress(plugin->handle, "destroy");
    if ( !plugin->destroy )
    {
        if ( exit_on_error ) error("Could not initialize %s: destroy method not found\n", plugin->name);
        return -1;
    }
#else
    dlerror();
    plugin->init = (dl_init_f) dlsym(plugin->handle, "init");
    char *ret = dlerror();
    if ( ret )
        plugin->init = NULL;
    else
        if ( args->verbose > 1 ) fprintf(stderr,"\tinit     .. ok\n");

    plugin->run = (dl_run_f) dlsym(plugin->handle, "run");
    ret = dlerror();
    if ( ret )
        plugin->run = NULL;
    else
        if ( args->verbose > 1 ) fprintf(stderr,"\trun      .. ok\n");

    if ( !plugin->init && !plugin->run )
    {
        if ( exit_on_error ) error("Could not initialize %s, neither run or init found \n", plugin->name);
        else if ( args->verbose > 1 ) fprintf(stderr,"\tinit/run .. not found\n");
        return -1;
    }

    plugin->version = (dl_version_f) dlsym(plugin->handle, "version");
    ret = dlerror();
    if ( ret )
    {
        if ( exit_on_error ) error("Could not initialize %s, version string not found\n", plugin->name);
        else if ( args->verbose > 1 ) fprintf(stderr,"\tversion  .. not found\n");
        return -1;
    }

    plugin->about = (dl_about_f) dlsym(plugin->handle, "about");
    ret = dlerror();
    if ( ret )
    {
        if ( exit_on_error ) error("Could not initialize %s: %s\n", plugin->name, ret);
        return -1;
    }

    plugin->usage = (dl_about_f) dlsym(plugin->handle, "usage");
    ret = dlerror();
    if ( ret )
        plugin->usage = plugin->about;

    if ( plugin->run ) return 0;

    plugin->process = (dl_process_f) dlsym(plugin->handle, "process");
    ret = dlerror();
    if ( ret )
    {
        if ( exit_on_error ) error("Could not initialize %s: %s\n", plugin->name, ret);
        return -1;
    }

    plugin->destroy = (dl_destroy_f) dlsym(plugin->handle, "destroy");
    ret = dlerror();
    if ( ret )
    {
        if ( exit_on_error ) error("Could not initialize %s: %s\n", plugin->name, ret);
        return -1;
    }
#endif

    return 0;
}

static void check_version(args_t *args)
{
    static int warned_bcftools = 0, warned_htslib = 0;
    const char *bver, *hver;
    args->plugin.version(&bver, &hver);
    if ( strcmp(bver,bcftools_version()) && !warned_bcftools )
    {
        fprintf(stderr,"WARNING: bcftools version mismatch .. bcftools at %s, the plugin \"%s\" at %s\n", bcftools_version(),args->plugin.name,bver);
        warned_bcftools = 1;
    }
    if ( strcmp(hver,hts_version()) && !warned_htslib )
    {
        fprintf(stderr,"WARNING: htslib version mismatch .. bcftools at %s, the plugin \"%s\" at %s\n", hts_version(),args->plugin.name,hver);
        warned_htslib = 1;
    }
}

static void init_plugin(args_t *args)
{
    int ret = args->plugin.init(args->plugin.argc,args->plugin.argv,args->hdr,args->hdr_out);
    if ( ret<0 ) error("The plugin exited with an error.\n");
    check_version(args);
    args->drop_header += ret;
}

static int cmp_plugin_name(const void *p1, const void *p2)
{
    plugin_t *a = (plugin_t*) p1;
    plugin_t *b = (plugin_t*) p2;
    return strcmp(a->name,b->name);
}

// If args=NULL then returns the number of plugins available. Otherwise prints the
// plugins on stdout and returns 0 on success.
static int list_plugins(args_t *args)
{
    plugin_t *plugins = NULL;
    int nplugins = 0, mplugins = 0;

    int count_only = 0;
    args_t _args;
    if ( !args )
    {
        memset(&_args,0,sizeof(_args));
        args = &_args;
        args->nplugin_paths = -1;
        count_only = 1;
    }
    init_plugin_paths(args);

    kstring_t str = {0,0,0};
    int plugin_ext_len = strlen(PLUGIN_EXT);
    int i;
    for (i=0; i<args->nplugin_paths; i++)
    {
        DIR *dp = opendir(args->plugin_paths[i]);
        if ( dp==NULL ) continue;

        struct dirent *ep;
        while ( (ep=readdir(dp)) )
        {
            int len = strlen(ep->d_name);
            if ( strcasecmp(PLUGIN_EXT,ep->d_name+len-plugin_ext_len) ) continue;
            str.l = 0;
            ksprintf(&str,"%s/%s", args->plugin_paths[i],ep->d_name);
            hts_expand(plugin_t, nplugins+1, mplugins, plugins);
            if ( load_plugin(args, str.s, 0, &plugins[nplugins]) < 0 ) continue;
            nplugins++;
            str.l = 0;
            kputs(ep->d_name, &str);
            int l = str.l - 1;
            while ( l>=0 && str.s[l]!='.' ) l--;
            if ( l>=0 ) str.s[l] = 0;
            free(plugins[nplugins-1].name);
            plugins[nplugins-1].name = strdup(str.s);  // use a short name
        }
        closedir(dp);
    }
    if ( count_only )
    {
        free(str.s);
        return nplugins;
    }
    if ( nplugins )
    {
        qsort(plugins, nplugins, sizeof(plugins[0]), cmp_plugin_name);

        for (i=0; i<nplugins; i++)
        {
            if ( args->verbose )
                printf("\n-- %s --\n%s", plugins[i].name, plugins[i].about());
            else
                printf("%s\n", plugins[i].name);
        }
        if ( args->verbose ) printf("\n");
    }
    else
        print_plugin_usage_hint(NULL);
    free(str.s);
    return nplugins ? 0 : 1;
}
int count_plugins(void)
{
    return list_plugins(NULL);
}

static void init_data(args_t *args)
{
    args->hdr = args->files->readers[0].header;
    args->hdr_out = bcf_hdr_dup(args->hdr);

    init_plugin(args);

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);

    if (args->record_cmd_line) bcf_hdr_append_version(args->hdr_out, args->argc, args->argv, "bcftools_plugin");
    if ( !args->drop_header )
    {
        char wmode[8];
        set_wmode(wmode,args->output_type,args->output_fname,args->clevel);
        args->out_fh = hts_open(args->output_fname ? args->output_fname : "-", wmode);
        if ( args->out_fh == NULL ) error("Can't write to \"%s\": %s\n", args->output_fname, strerror(errno));
        if ( args->n_threads ) hts_set_threads(args->out_fh, args->n_threads);
        if ( bcf_hdr_write(args->out_fh, args->hdr_out)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
    }
}

static void destroy_data(args_t *args)
{
    free(args->plugin.name);
    if ( args->plugin.destroy ) args->plugin.destroy();
#ifdef _WIN32
    FreeLibrary(args->plugin.handle);
#else
    dlclose(args->plugin.handle);
#endif
    if ( args->hdr_out ) bcf_hdr_destroy(args->hdr_out);
    if ( args->nplugin_paths>0 )
    {
        int i;
        for (i=0; i<args->nplugin_paths; i++) free(args->plugin_paths[i]);
        free(args->plugin_paths);
    }
    if ( args->filter )
        filter_destroy(args->filter);
    if (args->out_fh && hts_close(args->out_fh)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->output_fname);
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Run user defined plugin\n");
    fprintf(stderr, "Usage:   bcftools plugin <name> [OPTIONS] <file> [-- PLUGIN_OPTIONS]\n");
    fprintf(stderr, "         bcftools +name [OPTIONS] <file>  [-- PLUGIN_OPTIONS]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "VCF input options:\n");
    fprintf(stderr, "   -e, --exclude EXPR             Exclude sites for which the expression is true\n");
    fprintf(stderr, "   -i, --include EXPR             Select sites for which the expression is true\n");
    fprintf(stderr, "   -r, --regions REGION           Restrict to comma-separated list of regions\n");
    fprintf(stderr, "   -R, --regions-file FILE        Restrict to regions listed in a file\n");
    fprintf(stderr, "        --regions-overlap 0|1|2   Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n");
    fprintf(stderr, "   -t, --targets REGION           Similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "   -T, --targets-file FILE        Similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "        --targets-overlap 0|1|2   Include if POS in the region (0), record overlaps (1), variant overlaps (2) [0]\n");
    fprintf(stderr, "VCF output options:\n");
    fprintf(stderr, "       --no-version               Do not append version and command line to the header\n");
    fprintf(stderr, "   -o, --output FILE              Write output to a file [standard output]\n");
    fprintf(stderr, "   -O, --output-type u|b|v|z[0-9] u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n");
    fprintf(stderr, "       --threads INTT             Use multithreading with <int> worker threads [0]\n");
    fprintf(stderr, "Plugin options:\n");
    fprintf(stderr, "   -h, --help                     List plugin's options\n");
    fprintf(stderr, "   -l, --list-plugins             List available plugins. See BCFTOOLS_PLUGINS environment variable and man page for details\n");
    fprintf(stderr, "   -v, --verbose                  Print verbose information, -vv increases verbosity\n");
    fprintf(stderr, "   -V, --version                  Print version string and exit\n");
    fprintf(stderr, "\n");
    exit(1);
}

static int is_verbose(int argc, char *argv[])
{
    int c, verbose = 0, opterr_ori = opterr;
    static struct option loptions[] =
    {
        {"verbose",no_argument,NULL,'v'},
        {NULL,0,NULL,0}
    };
    opterr = 0;
    while ((c = getopt_long(argc, argv, "-v",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'v': verbose++; break;
            case 1:
            default: break;
        }
    }
    opterr = opterr_ori;
    optind = 0;
    return verbose;
}
int main_plugin(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->output_fname = "-";
    args->output_type = FT_VCF;
    args->n_threads = 0;
    args->record_cmd_line = 1;
    args->nplugin_paths = -1;
    args->clevel = -1;
    int regions_is_file = 0, targets_is_file = 0, usage_only = 0, version_only = 0;
    int regions_overlap = 1;
    int targets_overlap = 0;

    if ( argc==1 ) usage(args);
    char *plugin_name = NULL;
    if ( argv[1][0]!='-' )
    {
        args->verbose = is_verbose(argc, argv);
        plugin_name = argv[1]; 
        argc--; 
        argv++; 
        load_plugin(args, plugin_name, 1, &args->plugin);
        if ( args->plugin.run )
        {
            check_version(args);
            int ret = args->plugin.run(argc, argv);
            destroy_data(args);
            free(args);
            return ret;
        }
    }

    static struct option loptions[] =
    {
        {"version",no_argument,NULL,'V'},
        {"verbose",no_argument,NULL,'v'},
        {"help",no_argument,NULL,'h'},
        {"list-plugins",no_argument,NULL,'l'},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"threads",required_argument,NULL,9},
        {"include",required_argument,NULL,'i'},
        {"exclude",required_argument,NULL,'e'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"regions-overlap",required_argument,NULL,1},
        {"targets",required_argument,NULL,'t'},
        {"targets-file",required_argument,NULL,'T'},
        {"targets-overlap",required_argument,NULL,2},
        {"no-version",no_argument,NULL,8},
        {NULL,0,NULL,0}
    };
    char *tmp;
    while ((c = getopt_long(argc, argv, "h?o:O:r:R:t:T:li:e:vV",loptions,NULL)) >= 0)
    {
        switch (c) {
            case 'V': version_only = 1; break;
            case 'v': args->verbose++; break;
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
                };
                if ( optarg[1] )
                {
                    args->clevel = strtol(optarg+1,&tmp,10);
                    if ( *tmp || args->clevel<0 || args->clevel>9 ) error("Could not parse argument: --compression-level %s\n", optarg+1);
                }
                break;
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
            case 'l': args->plist_only = 1; break;
            case  1 :
                regions_overlap = parse_overlap_option(optarg);
                if ( regions_overlap < 0 ) error("Could not parse: --regions-overlap %s\n",optarg);
                break;
            case  2 :
                targets_overlap = parse_overlap_option(optarg);
                if ( targets_overlap < 0 ) error("Could not parse: --targets-overlap %s\n",optarg);
                break;
            case  9 : args->n_threads = strtol(optarg, 0, 0); break;
            case  8 : args->record_cmd_line = 0; break;
            case '?':
            case 'h': usage_only = 1; break;
            default: error("Unknown argument: %s\n", optarg);
        }
    }
    if ( args->plist_only )  return list_plugins(args);
    if ( !plugin_name ) usage(args);

    if ( version_only )
    {
        const char *bver, *hver;
        args->plugin.version(&bver, &hver);
        printf("bcftools  %s using htslib %s\n", bcftools_version(), hts_version());
        printf("plugin at %s using htslib %s\n\n", bver, hver);
        return 0;
    }

    if ( usage_only )
    {
        if ( args->plugin.usage )
            fprintf(stderr,"%s",args->plugin.usage());
        else
            fprintf(stderr,"Usage: bcftools +%s [General Options] -- [Plugin Options]\n",plugin_name);
        return 0;
    }

    char *fname = NULL;
    if ( optind>=argc || (argv[optind][0]=='-' && argv[optind][1]) )
    {
        args->plugin.argc = argc - optind + 1;
        args->plugin.argv = argv + optind - 1;

        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin
        else if ( optind>=argc ) usage(args);
        else
        {
            optind = 1;
            init_plugin(args);
        }
    }
    else
    {
        fname = argv[optind];
        args->plugin.argc = argc - optind;
        args->plugin.argv = argv + optind;
    }
    optind = 0;

    args->files = bcf_sr_init();
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
        args->files->collapse |= COLLAPSE_SOME;
    }
    if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to read from %s: %s\n", !strcmp("-",fname)?"standard input":fname,bcf_sr_strerror(args->files->errnum));

    init_data(args);
    while ( bcf_sr_next_line(args->files) )
    {
        bcf1_t *line = bcf_sr_get_line(args->files,0);
        if ( args->filter )
        {
            int pass = filter_test(args->filter, line, NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) pass = pass ? 0 : 1;
            if ( !pass ) continue;
        }
        line = args->plugin.process(line);
        if ( line )
        {
            if ( line->errcode ) error("[E::main_plugin] Unchecked error (%d), exiting\n",line->errcode);
            if ( bcf_write1(args->out_fh, args->hdr_out, line)!=0 ) error("[%s] Error: cannot write to %s\n", __func__,args->output_fname);
        }
    }
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return 0;
}

#else /* ENABLE_BCF_PLUGINS */

int main_plugin(int argc, char *argv[])
{
    fprintf(stderr, "bcftools plugins are disabled.  To use them, you will need to rebuild\n"
	    "bcftools from the source distribution with plugins enabled.\n");
    return 1;
}

#endif /* ENABLE_BCF_PLUGINS */
