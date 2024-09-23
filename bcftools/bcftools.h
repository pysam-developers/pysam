/*  bcftools.h -- utility function declarations.

    Copyright (C) 2013-2024 Genome Research Ltd.

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
#include <htslib/synced_bcf_reader.h>
#include <htslib/kfunc.h>
#include <math.h>

#define FT_TAB_TEXT 0       // custom tab-delimited text file
#define FT_GZ 1
#define FT_VCF 2
#define FT_VCF_GZ (FT_GZ|FT_VCF)
#define FT_BCF (1<<2)
#define FT_BCF_GZ (FT_GZ|FT_BCF)
#define FT_STDIN (1<<3)

char *bcftools_version(void);

/// Report an error and exit -1
void error(const char *format, ...) HTS_NORETURN HTS_FORMAT(HTS_PRINTF_FMT, 1, 2);

/// Report an error and exit -1.  If errno != 0, appends strerror(errno).
//  Note: unlike error() above, the message should not end with "\n" as a
//  newline will be added by the function.
void error_errno(const char *format, ...) HTS_NORETURN HTS_FORMAT(HTS_PRINTF_FMT, 1, 2);

// For on the fly index creation with --write-index
int init_index2(htsFile *fh, bcf_hdr_t *hdr, const char *fname, char **idx_fname, int idx_fmt);
int init_index(htsFile *fh, bcf_hdr_t *hdr, const char *fname, char **idx_fname);

// Used to set args->write_index in CLI.
// It will be true if set correctly.
// Note due to HTS_FMT_CSI being zero we have to use an additional bit.
int write_index_parse(char *arg);

void bcf_hdr_append_version(bcf_hdr_t *hdr, int argc, char **argv, const char *cmd);
const char *hts_bcf_wmode(int file_type);
const char *hts_bcf_wmode2(int file_type, const char *fname);
void set_wmode(char dst[8], int file_type, const char *fname, int compression_level);  // clevel: 0-9 with or zb type, -1 unset
char *init_tmp_prefix(const char *prefix);
int read_AF(bcf_sr_regions_t *tgt, bcf1_t *line, double *alt_freq);
int parse_overlap_option(const char *arg);

// Default sort order: chr,pos,alleles
int cmp_bcf_pos(const void *aptr, const void *bptr);
int cmp_bcf_pos_ref_alt(const void *aptr, const void *bptr);

static inline int iupac2bitmask(char iupac)
{
    const int A = 1;
    const int C = 2;
    const int G = 4;
    const int T = 8;
    if ( iupac >= 97 ) iupac -= 32;
    if ( iupac == 'A' ) return A;
    if ( iupac == 'C' ) return C;
    if ( iupac == 'G' ) return G;
    if ( iupac == 'T' ) return T;
    if ( iupac == 'M' ) return A|C;
    if ( iupac == 'R' ) return A|G;
    if ( iupac == 'W' ) return A|T;
    if ( iupac == 'S' ) return C|G;
    if ( iupac == 'Y' ) return C|T;
    if ( iupac == 'K' ) return G|T;
    if ( iupac == 'V' ) return A|C|G;
    if ( iupac == 'H' ) return A|C|T;
    if ( iupac == 'D' ) return A|G|T;
    if ( iupac == 'B' ) return C|G|T;
    if ( iupac == 'N' ) return A|C|G|T;
    return -1;
}
static inline char bitmask2iupac(int bitmask)
{
    const char iupac[16] = {'.','A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N'};
    if ( bitmask <= 0 || bitmask > 15 ) return 0;
    return iupac[bitmask];
}

static inline int iupac_consistent(char iupac, char nt)
{
    static const char iupac_mask[90] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,14,2,
        13,0,0,4,11,0,0,12,0,3,15,0,0,0,5,6,8,0,7,9,0,10
    };
    if ( iupac > 89 ) return 0;
    if ( nt > 90 ) nt -=  32;  // lowercase
    if ( nt=='A' ) nt = 1;
    else if ( nt=='C' ) nt = 2;
    else if ( nt=='G' ) nt = 4;
    else if ( nt=='T' ) nt = 8;
    return iupac_mask[(int)iupac] & nt ? 1 : 0;
}

static inline char nt_to_upper(char nt)
{
    if ( nt < 97 ) return nt;
    return nt - 32;
}

static inline double phred_score(double prob)
{
    if ( prob==0 ) return 99;
    prob = -4.3429*log(prob);
    return prob>99 ? 99 : prob;
}

static inline double calc_binom_two_sided(int na, int nb, double aprob)
{
    if ( !na && !nb ) return -1;
    if ( na==nb ) return 1;

    // kfunc.h implements kf_betai, which is the regularized beta function  P(X<=k/N;p) = I_{1-p}(N-k,k+1)

    double prob = na > nb ? 2 * kf_betai(na, nb+1, aprob) : 2 * kf_betai(nb, na+1, aprob);

    if ( prob > 1 ) prob = 1;   // this can happen, machine precision error, eg. kf_betai(1,0,0.5)
    return prob;
}
static inline double calc_binom_one_sided(int na, int nb, double aprob, int ge)
{
    return ge ? kf_betai(na, nb + 1, aprob) : kf_betai(nb, na + 1, 1 - aprob);
}

static const uint64_t bcf_double_missing    = 0x7ff0000000000001;
static const uint64_t bcf_double_vector_end = 0x7ff0000000000002;
static inline void bcf_double_set(double *ptr, uint64_t value)
{
    union { uint64_t i; double d; } u;
    u.i = value;
    *ptr = u.d;
}
static inline int bcf_double_test(double d, uint64_t value)
{
    union { uint64_t i; double d; } u;
    u.d = d;
    return u.i==value ? 1 : 0;
}
#define bcf_double_set_vector_end(x) bcf_double_set(&(x),bcf_double_vector_end)
#define bcf_double_set_missing(x)    bcf_double_set(&(x),bcf_double_missing)
#define bcf_double_is_vector_end(x)  bcf_double_test((x),bcf_double_vector_end)
#define bcf_double_is_missing(x)     bcf_double_test((x),bcf_double_missing)
#define bcf_double_is_missing_or_vector_end(x)     (bcf_double_test((x),bcf_double_missing) || bcf_double_test((x),bcf_double_vector_end))

static inline int get_unseen_allele(bcf1_t *line)
{
    int i;
    for (i=1; i<line->n_allele; i++)
    {
        if ( !strcmp(line->d.allele[i],"<*>") ) return i;
        if ( !strcmp(line->d.allele[i],"<NON_REF>") ) return i;
        if ( !strcmp(line->d.allele[i],"<X>") ) return i;
    }
    return 0;
}

#endif
