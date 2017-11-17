/*  main.c -- main bcftools command front-end.

    Copyright (C) 2012-2016 Genome Research Ltd.

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
int main_plugin(int argc, char *argv[]);
int main_consensus(int argc, char *argv[]);
int main_csq(int argc, char *argv[]);
int bam_mpileup(int argc, char *argv[]);
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
    { .func  = main_plugin,
      .alias = "plugin",
#ifdef ENABLE_BCF_PLUGINS
      .help  = "user-defined plugins"
#else
      /* Do not advertise when plugins disabled. */
      .help  = "-user-defined plugins"
#endif
    },
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
    { .func  = bam_mpileup,
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
    fprintf(fp,"\n");
    fprintf(fp,
            " Most commands accept VCF, bgzipped VCF, and BCF with the file type detected\n"
            " automatically even when streaming from a pipe. Indexed VCF and BCF will work\n"
            " in all situations. Un-indexed VCF and BCF and streams will work in most but\n"
            " not all situations.\n");
    fprintf(fp,"\n");
}

int main(int argc, char *argv[])
{
    if (argc < 2) { usage(stderr); return 1; }

    if (strcmp(argv[1], "version") == 0 || strcmp(argv[1], "--version") == 0 || strcmp(argv[1], "-v") == 0) {
        printf("bcftools %s\nUsing htslib %s\nCopyright (C) 2016 Genome Research Ltd.\n", bcftools_version(), hts_version());
#if USE_GPL
        printf("License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n");
#else
        printf("License Expat: The MIT/Expat license\n");
#endif
        printf("This is free software: you are free to change and redistribute it.\nThere is NO WARRANTY, to the extent permitted by law.\n");
        return 0;
    }
    else if (strcmp(argv[1], "--version-only") == 0) {
        printf("%s+htslib-%s\n", bcftools_version(), hts_version());
        return 0;
    }
    else if (strcmp(argv[1], "help") == 0 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) {
        if (argc == 2) { usage(stdout); return 0; }
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
    fprintf(stderr, "[E::%s] unrecognized command '%s'\n", __func__, argv[1]);
    return 1;
}

