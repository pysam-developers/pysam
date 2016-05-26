/*  bcftools.h -- utility function declarations.

    Copyright (C) 2013 Genome Research Ltd.

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

#ifndef BCFTOOLS_H
#define BCFTOOLS_H

#include <stdarg.h>
#include <htslib/hts_defs.h>
#include <htslib/vcf.h>
#include <math.h>

#define FT_GZ 1
#define FT_VCF 2
#define FT_VCF_GZ (FT_GZ|FT_VCF)
#define FT_BCF (1<<2)
#define FT_BCF_GZ (FT_GZ|FT_BCF)
#define FT_STDIN (1<<3)

char *bcftools_version(void);
void error(const char *format, ...) HTS_NORETURN;
void bcf_hdr_append_version(bcf_hdr_t *hdr, int argc, char **argv, const char *cmd);
const char *hts_bcf_wmode(int file_type);

void *smalloc(size_t size);     // safe malloc

static inline char gt2iupac(char a, char b)
{
    static const char iupac[4][4] = { {'A','M','R','W'},{'M','C','S','Y'},{'R','S','G','K'},{'W','Y','K','T'} };
    if ( a>='a' ) a -= 'a' - 'A';
    if ( b>='a' ) b -= 'a' - 'A';
    if ( a=='A' ) a = 0;
    else if ( a=='C' ) a = 1;
    else if ( a=='G' ) a = 2;
    else if ( a=='T' ) a = 3;
    else return 'N';
    if ( b=='A' ) b = 0;
    else if ( b=='C' ) b = 1;
    else if ( b=='G' ) b = 2;
    else if ( b=='T' ) b = 3;
    else return 'N';
    return iupac[(int)a][(int)b];
}

static inline double phred_score(double prob)
{
    if ( prob==0 ) return 99;
    prob = -4.3429*log(prob);
    return prob>99 ? 99 : prob;
}

#endif
