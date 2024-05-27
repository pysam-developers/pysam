/*  vcfhead.c -- view VCF/BCF file headers.

    Copyright (C) 2021 University of Glasgow.
    Copyright (C) 2023 Genome Research Ltd.

    Author: John Marshall <jmarshall@hey.com>

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

#include <getopt.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <htslib/kstring.h>
#include <htslib/vcf.h>

#include "bcftools.h"

int main_vcfhead(int argc, char *argv[])
{
    static const char usage[] =
"\n"
"About: Displays VCF/BCF headers and optionally the first few variant records\n"
"Usage: bcftools head [OPTION]... [FILE]\n"
"\n"
"Options:\n"
"  -h, --headers INT    Display INT header lines [all]\n"
"  -n, --records INT    Display INT variant record lines [none]\n"
"  -s, --samples INT    Display INT records starting with the #CHROM header line [none]\n"
"\n";

    static const struct option loptions[] = {
        { "headers", required_argument, NULL, 'h' },
        { "records", required_argument, NULL, 'n' },
        { "samples", required_argument, NULL, 's' },
        { NULL, 0, NULL, 0 }
    };

    int all_headers = 1;
    int samples = 0;
    uint64_t nheaders = 0;
    uint64_t nrecords = 0;

    int c, nargs;
    while ((c = getopt_long(argc, argv, "h:n:s:", loptions, NULL)) >= 0)
        switch (c) {
        case 'h': all_headers = 0; nheaders = strtoull(optarg, NULL, 0); break;
        case 'n': nrecords = strtoull(optarg, NULL, 0); break;
        case 's': nrecords = strtoull(optarg, NULL, 0); samples = 1; break;
        default:
            fputs(usage, stderr);
            return EXIT_FAILURE;
        }

    if ( samples && all_headers ) all_headers = 0;

    nargs = argc - optind;
    if (nargs == 0 && isatty(STDIN_FILENO)) {
        fputs(usage, stdout);
        return EXIT_SUCCESS;
    }
    else if (nargs > 1) {
        fputs(usage, stderr);
        return EXIT_FAILURE;
    }

    const char *fname = (nargs == 1)? argv[optind] : "-";
    vcfFile *fp = bcf_open(fname, "r");
    if (fp == NULL) {
        if (strcmp(fname, "-") != 0)
            error_errno("[%s] Can't open \"%s\"", __func__, fname);
        else
            error_errno("[%s] Can't open standard input", __func__);
    }

    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (hdr == NULL) {
        bcf_close(fp);
        if (strcmp(fname, "-") != 0)
            error("[%s] Can't read headers from \"%s\"\n", __func__, fname);
        else
            error("[%s] Can't read headers\n", __func__);
    }

    kstring_t str = KS_INITIALIZE;

    if (all_headers) {
        bcf_hdr_format(hdr, 0, &str);
        fputs(ks_str(&str), stdout);
    }
    else if (nheaders > 0 || samples ) {
        bcf_hdr_format(hdr, 0, &str);
        char *lim = str.s;
        uint64_t n;
        int samples_printed = 0;
        for (n = 0; n < nheaders; n++) {
            if ( samples && !strncmp(lim,"#CHROM\t",7) ) samples_printed = 1;
            lim = strchr(lim, '\n');
            if (lim) lim++;
            else break;
        }
        if ( nheaders )
        {
            char tmp;
            if (lim) { tmp = *lim; *lim = '\0'; }
            fputs(ks_str(&str), stdout);
            if (lim) *lim = tmp;
        }
        if ( lim && samples && !samples_printed )
        {
            while ( lim && *lim )
            {
                if ( !strncmp(lim,"#CHROM\t",7) ) { fputs(lim, stdout); break; }
                lim = strchr(lim, '\n');
                if (lim) lim++;
                else break;
            }
        }
    }

    if (nrecords > 0) {
        bcf1_t *rec = bcf_init();
        uint64_t n;
        for (n = 0; n < nrecords && bcf_read(fp, hdr, rec) >= 0; n++) {
            ks_clear(&str);
            if (vcf_format(hdr, rec, &str) >= 0)
                fputs(ks_str(&str), stdout);
            else
                fprintf(stderr, "[%s] Record #%"PRIu64 " is invalid\n", __func__, n+1);
        }
        bcf_destroy(rec);
    }

    ks_free(&str);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);

    return EXIT_SUCCESS;
}
