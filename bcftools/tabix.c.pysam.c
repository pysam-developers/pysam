#include "bcftools.pysam.h"

/*  tabix.c -- tabix subcommand.

    Copyright (C) 2012 Broad Institute.
    Copyright (C) 2013, 2016 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

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
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <htslib/bgzf.h>
#include <htslib/tbx.h>

int main_tabix(int argc, char *argv[])
{
    int c, min_shift = -1, is_force = 0, is_all = 0, detect = 1;
    tbx_conf_t conf = tbx_conf_gff;
    while ((c = getopt(argc, argv, "0fap:s:b:e:S:c:m:")) >= 0)
        if (c == '0') conf.preset |= TBX_UCSC;
        else if (c == 'f') is_force = 1;
        else if (c == 'a') is_all = 1;
        else if (c == 'm') min_shift = atoi(optarg);
        else if (c == 's') conf.sc = atoi(optarg);
        else if (c == 'b') conf.bc = atoi(optarg);
        else if (c == 'e') conf.ec = atoi(optarg);
        else if (c == 'c') conf.meta_char = *optarg;
        else if (c == 'S') conf.line_skip = atoi(optarg);
        else if (c == 'p') {
            if (strcmp(optarg, "gff") == 0) conf = tbx_conf_gff;
            else if (strcmp(optarg, "bed") == 0) conf = tbx_conf_bed;
            else if (strcmp(optarg, "sam") == 0) conf = tbx_conf_sam;
            else if (strcmp(optarg, "vcf") == 0) conf = tbx_conf_vcf;
            else {
                fprintf(bcftools_stderr, "The type '%s' not recognised\n", optarg);
                return 1;
            detect = 0;
            }

        }
    if (optind == argc) {
        fprintf(bcftools_stderr, "\nUsage: bcftools tabix [options] <in.gz> [reg1 [...]]\n\n");
        fprintf(bcftools_stderr, "Options: -p STR    preset: gff, bed, sam or vcf [gff]\n");
        fprintf(bcftools_stderr, "         -s INT    column number for sequence names (suppressed by -p) [1]\n");
        fprintf(bcftools_stderr, "         -b INT    column number for region start [4]\n");
        fprintf(bcftools_stderr, "         -e INT    column number for region end (if no end, set INT to -b) [5]\n");
        fprintf(bcftools_stderr, "         -0        specify coordinates are zero-based\n");
        fprintf(bcftools_stderr, "         -S INT    skip first INT lines [0]\n");
        fprintf(bcftools_stderr, "         -c CHAR   skip lines starting with CHAR [null]\n");
        fprintf(bcftools_stderr, "         -a        print all records\n");
        fprintf(bcftools_stderr, "         -f        force to overwrite existing index\n");
        fprintf(bcftools_stderr, "         -m INT    set the minimal interval size to 1<<INT; 0 for the old tabix index [0]\n");
        fprintf(bcftools_stderr, "\n");
        return 1;
    }
    if (is_all) { // read without random access
        kstring_t s;
        BGZF *fp;
        s.l = s.m = 0; s.s = 0;
        fp = bgzf_open(argv[optind], "r");
        while (bgzf_getline(fp, '\n', &s) >= 0) fputs(s.s, bcftools_stdout) & fputc('\n', bcftools_stdout);
        bgzf_close(fp);
        free(s.s);
    } else if (optind + 2 > argc) { // create index
        if ( detect )
        {
            // auto-detect file type by file name
            int l = strlen(argv[optind]);
            int strcasecmp(const char *s1, const char *s2);
            if (l>=7 && strcasecmp(argv[optind]+l-7, ".gff.gz") == 0) conf = tbx_conf_gff;
            else if (l>=7 && strcasecmp(argv[optind]+l-7, ".bed.gz") == 0) conf = tbx_conf_bed;
            else if (l>=7 && strcasecmp(argv[optind]+l-7, ".sam.gz") == 0) conf = tbx_conf_sam;
            else if (l>=7 && strcasecmp(argv[optind]+l-7, ".vcf.gz") == 0) conf = tbx_conf_vcf;
        }

        if (!is_force) {
            char *fn;
            FILE *fp;
            fn = (char*)malloc(strlen(argv[optind]) + 5);
            strcat(strcpy(fn, argv[optind]), min_shift <= 0? ".tbi" : ".csi");
            if ((fp = fopen(fn, "rb")) != 0) {
                fclose(fp);
                free(fn);
                fprintf(bcftools_stderr, "[E::%s] the index file exists; use option '-f' to overwrite\n", __func__);
                return 1;
            }
            free(fn);
        }
        if ( tbx_index_build(argv[optind], min_shift, &conf) )
        {
            fprintf(bcftools_stderr,"tbx_index_build failed: Is the file bgzip-compressed? Was wrong -p [type] option used?\n");
            return 1;
        }
    } else { // read with random access
        tbx_t *tbx;
        BGZF *fp;
        kstring_t s;
        int i;
        if ((tbx = tbx_index_load(argv[optind])) == 0) return 1;
        if ((fp = bgzf_open(argv[optind], "r")) == 0) return 1;
        s.s = 0; s.l = s.m = 0;
        for (i = optind + 1; i < argc; ++i) {
            hts_itr_t *itr;
            if ((itr = tbx_itr_querys(tbx, argv[i])) == 0) continue;
            while (tbx_bgzf_itr_next(fp, tbx, itr, &s) >= 0) fputs(s.s, bcftools_stdout) & fputc('\n', bcftools_stdout);
            tbx_itr_destroy(itr);
        }
        free(s.s);
        bgzf_close(fp);
        tbx_destroy(tbx);
    }
    return 0;
}
