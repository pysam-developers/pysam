/*  version.c -- report version numbers for plugins.

    Copyright (C) 2014-2021 Genome Research Ltd.

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <errno.h>
#include <htslib/hts.h>
#include "bcftools.h"
#include "version.h"

void version(const char **bcftools_version, const char **htslib_version)
{
    *bcftools_version = BCFTOOLS_VERSION;
    *htslib_version = hts_version();
}

void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    exit(-1);
}

void error_errno(const char *format, ...)
{
    va_list ap;
    int e = errno;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    if (e) {
        fprintf(stderr, ": %s\n", strerror(e));
    } else {
        fprintf(stderr, "\n");
    }
    exit(-1);
}

const char *hts_bcf_wmode(int file_type)
{
    if ( file_type == FT_BCF ) return "wbu";    // uncompressed BCF
    if ( file_type & FT_BCF ) return "wb";      // compressed BCF
    if ( file_type & FT_GZ ) return "wz";       // compressed VCF
    return "w";                                 // uncompressed VCF
}

const char *hts_bcf_wmode2(int file_type, const char *fname)
{
    if ( !fname ) return hts_bcf_wmode(file_type);
    int len = strlen(fname);
    if ( len >= 4 && !strcasecmp(".bcf",fname+len-4) ) return hts_bcf_wmode(FT_BCF|FT_GZ);
    if ( len >= 4 && !strcasecmp(".vcf",fname+len-4) ) return hts_bcf_wmode(FT_VCF);
    if ( len >= 7 && !strcasecmp(".vcf.gz",fname+len-7) ) return hts_bcf_wmode(FT_VCF|FT_GZ);
    if ( len >= 8 && !strcasecmp(".vcf.bgz",fname+len-8) ) return hts_bcf_wmode(FT_VCF|FT_GZ);
    return hts_bcf_wmode(file_type);
}

void set_wmode(char dst[8], int file_type, const char *fname, int clevel)
{
    const char *ret = NULL;
    int len = fname ? strlen(fname) : 0;
    if ( len >= 4 && !strcasecmp(".bcf",fname+len-4) ) ret = hts_bcf_wmode(FT_BCF|FT_GZ);
    else if ( len >= 4 && !strcasecmp(".vcf",fname+len-4) ) ret = hts_bcf_wmode(FT_VCF);
    else if ( len >= 7 && !strcasecmp(".vcf.gz",fname+len-7) ) ret = hts_bcf_wmode(FT_VCF|FT_GZ);
    else if ( len >= 8 && !strcasecmp(".vcf.bgz",fname+len-8) ) ret = hts_bcf_wmode(FT_VCF|FT_GZ);
    else ret = hts_bcf_wmode(file_type);
    if ( clevel>=0 && clevel<=9 )
    {
        if ( strchr(ret,'v') || strchr(ret,'u') ) error("Error: compression level (%d) cannot be set on uncompressed streams (%s)\n",clevel,fname);
        len = strlen(ret);
        if ( len>6 ) error("Fixme: %s\n", ret);
        sprintf(dst, "%s%d", ret, clevel);
    }
    else
        strcpy(dst, ret);
}

int parse_overlap_option(const char *arg)
{
    if ( strcasecmp(arg, "pos") == 0 || strcmp(arg, "0") == 0 ) return 0;
    else if ( strcasecmp(arg, "record") == 0 || strcmp(arg, "1") == 0 ) return 1;
    else if ( strcasecmp(arg, "variant") == 0 || strcmp(arg, "2") == 0 ) return 2;
    else return -1;
}
