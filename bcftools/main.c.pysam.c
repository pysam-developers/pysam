#include "bcftools.pysam.h"

/*  main.c -- main bcftools command front-end.

    Copyright (C) 2012-2021 Genome Research Ltd.

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
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <htslib/hts.h>
#include "version.h"
#include "bcftools.h"

int main_tabix(int argc, char *argv[]);
int main_vcfindex(int argc, char *argv[]);
int main_vcfstats(int argc, char *argv[]);
int main_vcfisec(int argc, char *argv[]);
int main_vcfmerge(int argc, char *argv[]);
int main_vcfquery(int argc, char *argv[]);
int main_vcffilter(int argc, char *argv[]);
int main_vcfsom(int argc, char *argv[]);
int main_vcfnorm(int argc, char *argv[]);
int main_vcfgtcheck(int argc, char *argv[]);
int main_vcfview(int argc, char *argv[]);
int main_vcfhead(int argc, char *argv[]);
int main_vcfcall(int argc, char *argv[]);
int main_vcfannotate(int argc, char *argv[]);
int main_vcfroh(int argc, char *argv[]);
int main_vcfconcat(int argc, char *argv[]);
int main_reheader(int argc, char *argv[]);
int main_vcfconvert(int argc, char *argv[]);
int main_vcfcnv(int argc, char *argv[]);
#if USE_GPL
int main_polysomy(int argc, char *argv[]);
#endif
#ifdef ENABLE_BCF_PLUGINS
int main_plugin(int argc, char *argv[]);
int count_plugins(void);
#endif
int main_consensus(int argc, char *argv[]);
int main_csq(int argc, char *argv[]);
int main_mpileup(int argc, char *argv[]);
int main_sort(int argc, char *argv[]);

typedef struct
{
    int (*func)(int, char*[]);
    const char *alias, *help;
}
cmd_t;

static cmd_t cmds[] =
{
    { .func  = NULL,
      .alias = "Indexing",
      .help  = NULL
    },
    { .func = main_vcfindex,
      .alias = "index",
      .help = "index VCF/BCF files"
    },
    { .func = main_tabix,
      .alias = "tabix",
      .help = "-tabix for BGZF'd BED, GFF, SAM, VCF and more" // do not advertise; only keep here for testing
    },

    { .func  = NULL,
      .alias = "VCF/BCF manipulation",
      .help  = NULL
    },

    { .func  = main_vcfannotate,
      .alias = "annotate",
      .help  = "annotate and edit VCF/BCF files",
    },
    { .func  = main_vcfconcat,
      .alias = "concat",
      .help  = "concatenate VCF/BCF files from the same set of samples"
    },
    { .func  = main_vcfconvert,
      .alias = "convert",
      .help  = "convert VCF/BCF files to different formats and back"
    },
    { .func  = main_vcfhead,
      .alias = "head",
      .help  = "view VCF/BCF file headers"
    },
    { .func  = main_vcfisec,
      .alias = "isec",
      .help  = "intersections of VCF/BCF files"
    },
    { .func  = main_vcfmerge,
      .alias = "merge",
      .help  = "merge VCF/BCF files files from non-overlapping sample sets"
    },
    { .func  = main_vcfnorm,
      .alias = "norm",
      .help  = "left-align and normalize indels"
    },
#ifdef ENABLE_BCF_PLUGINS
    { .func  = main_plugin,
      .alias = "plugin",
      .help  = "user-defined plugins"
    },
#endif
    { .func  = main_vcfquery,
      .alias = "query",
      .help  = "transform VCF/BCF into user-defined formats"
    },
    { .func  = main_reheader,
      .alias = "reheader",
      .help  = "modify VCF/BCF header, change sample names"
    },
    { .func  = main_sort,
      .alias = "sort",
      .help  = "sort VCF/BCF file"
    },
    { .func  = main_vcfview,
      .alias = "view",
      .help  = "VCF/BCF conversion, view, subset and filter VCF/BCF files"
    },

    { .func  = NULL,
      .alias = "VCF/BCF analysis",
      .help  = NULL
    },

    { .func  = main_vcfcall,
      .alias = "call",
      .help  = "SNP/indel calling"
    },
    { .func  = main_consensus,
      .alias = "consensus",
      .help  = "create consensus sequence by applying VCF variants"
    },
    { .func  = main_vcfcnv,
      .alias = "cnv",
      .help  = "HMM CNV calling"
    },
    { .func  = main_csq,
      .alias = "csq",
      .help  = "call variation consequences"
    },
    { .func  = main_vcffilter,
      .alias = "filter",
      .help  = "filter VCF/BCF files using fixed thresholds"
    },
    { .func  = main_vcfgtcheck,
      .alias = "gtcheck",
      .help  = "check sample concordance, detect sample swaps and contamination"
    },
    { .func  = main_mpileup,
        .alias = "mpileup",
        .help  = "multi-way pileup producing genotype likelihoods"
    },
#if USE_GPL
    { .func  = main_polysomy,
      .alias = "polysomy",
      .help  = "detect number of chromosomal copies",
    },
#endif
    { .func  = main_vcfroh,
      .alias = "roh",
      .help  = "identify runs of autozygosity (HMM)",
    },
    { .func  = main_vcfstats,
      .alias = "stats",
      .help  = "produce VCF/BCF stats"
    },

    { .func  = main_vcfsom,
      .alias = "som",
      .help  = "-filter using Self-Organized Maps (experimental)"   // do not advertise

    },
    { .func  = NULL,
      .alias = NULL,
      .help  = NULL
    }
};

char *bcftools_version(void)
{
    return BCFTOOLS_VERSION;
}

static void usage(FILE *fp)
{
    fprintf(fp, "\n");
    fprintf(fp, "Program: bcftools (Tools for variant calling and manipulating VCFs and BCFs)\n");
#if USE_GPL
    fprintf(fp, "License: GNU GPLv3+, due to use of the GNU Scientific Library\n");
#endif
    fprintf(fp, "Version: %s (using htslib %s)\n", bcftools_version(), hts_version());
    fprintf(fp, "\n");
    fprintf(fp, "Usage:   bcftools [--version|--version-only] [--help] <command> <argument>\n");
    fprintf(fp, "\n");
    fprintf(fp, "Commands:\n");

    int i = 0;
    const char *sep = NULL;
    while (cmds[i].alias)
    {
        if ( !cmds[i].func ) sep = cmds[i].alias;
        if ( sep )
        {
            fprintf(fp, "\n -- %s\n", sep);
            sep = NULL;
        }
        if ( cmds[i].func && cmds[i].help[0]!='-' ) fprintf(fp, "    %-12s %s\n", cmds[i].alias, cmds[i].help);
        i++;
    }
#if ENABLE_BCF_PLUGINS
    fprintf(fp,"\n -- Plugins (collection of programs for calling, file manipulation & analysis)\n");
    int nplugins = count_plugins();
    if ( nplugins )
        fprintf(fp,"    %d plugins available, run \"bcftools plugin -lv\" to see a complete list\n", nplugins);
    else
        fprintf(fp,"    0 plugins available, run \"bcftools plugin -l\" for help\n");
#endif
    fprintf(fp,"\n");
    fprintf(fp,
            " Most commands accept VCF, bgzipped VCF, and BCF with the file type detected\n"
            " automatically even when streaming from a pipe. Indexed VCF and BCF will work\n"
            " in all situations. Un-indexed VCF and BCF and streams will work in most but\n"
            " not all situations.\n");
    fprintf(fp,"\n");
}

// This is a tricky one, but on Windows the filename wildcard expansion is done by
// the application and not by the shell, as traditionally it never had a "shell".
// Even now, DOS and Powershell do not do this expansion (but bash does).
//
// This means that Mingw/Msys implements code before main() that takes e.g. "*" and
// expands it up to a list of matching filenames.  This in turn breaks things like
// specifying "*" as a region (all the unmapped reads).  We take a hard line here -
// filename expansion is the task of the shell, not our application!
#ifdef _WIN32
int _CRT_glob = 0;
#endif

int bcftools_main(int argc, char *argv[])
{
    if (argc < 2) { usage(bcftools_stderr); return 1; }

    if (strcmp(argv[1], "version") == 0 || strcmp(argv[1], "--version") == 0 || strcmp(argv[1], "-v") == 0) {
        fprintf(bcftools_stdout, "bcftools %s\nUsing htslib %s\nCopyright (C) 2023 Genome Research Ltd.\n", bcftools_version(), hts_version());
#if USE_GPL
        fprintf(bcftools_stdout, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n");
#else
        fprintf(bcftools_stdout, "License Expat: The MIT/Expat license\n");
#endif
        fprintf(bcftools_stdout, "This is free software: you are free to change and redistribute it.\nThere is NO WARRANTY, to the extent permitted by law.\n");
        return 0;
    }
    else if (strcmp(argv[1], "--version-only") == 0) {
        fprintf(bcftools_stdout, "%s+htslib-%s\n", bcftools_version(), hts_version());
        return 0;
    }
    else if (strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) {
        if (argc == 2) { usage(bcftools_stdout); return 0; }
        // Otherwise change "bcftools help COMMAND [...]" to "bcftools COMMAND";
        // main_xyz() functions by convention display the subcommand's usage
        // when invoked without any arguments.
        argv++;
        argc = 2;
    }
    else if ( argv[1][0]=='+' )
    {
        // "bcftools plugin name" can be run as "bcftools +name"
        argv[1]++;
        argv[0] = "plugin";
        argv--;
        argc++;
    }

    int i = 0;
    while (cmds[i].alias)
    {
        if (cmds[i].func && strcmp(argv[1],cmds[i].alias)==0)
        {
            return cmds[i].func(argc-1,argv+1);
        }
        i++;
    }
    fprintf(bcftools_stderr, "[E::%s] unrecognized command '%s'\n", __func__, argv[1]);
    return 1;
}

