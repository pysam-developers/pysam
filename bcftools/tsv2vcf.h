/*  tsv2vcf.h -- convert from whitespace-separated fields to VCF

    Copyright (C) 2014 Genome Research Ltd.

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
    THE SOFTWARE.
*/

#ifndef __TSV2VCF_H__
#define __TSV2VCF_H__

#include <htslib/vcf.h>

typedef struct _tsv_t tsv_t;
typedef int (*tsv_setter_t)(tsv_t *, bcf1_t *, void *);

typedef struct
{
    char *name;
    tsv_setter_t setter;
    void *usr;
}
tsv_col_t;

struct _tsv_t
{
    int ncols, icol;
    tsv_col_t *cols;
    char *se, *ss;
};

tsv_t *tsv_init(const char *str);
void tsv_destroy(tsv_t *tsv);
int tsv_register(tsv_t *tsv, const char *id, tsv_setter_t setter, void *usr);

/**
 *  tsv_parse() - parse tsv line and fill VCF record
 *  Returns 0 on success or -1 on parse error
 */
int tsv_parse(tsv_t *tsv, bcf1_t *rec, char *str);

/**
 *  tstv_next() - position ss,se to next field; first pass with ss=se=str
 *  Returns 0 on success, or -1 if no more fields
 */
static inline int tsv_next(tsv_t *tsv)
{
    if ( !*tsv->se ) return -1;
    if ( tsv->ss==tsv->se )
    {
        while ( *tsv->se && !isspace(*tsv->se) ) tsv->se++;
        return 0;
    }
    while ( *tsv->se && isspace(*tsv->se) ) tsv->se++;
    tsv->ss = tsv->se;
    while ( *tsv->se && !isspace(*tsv->se) ) tsv->se++;
    return 0;
}

/**
 *  The setters return 0 on success or negative value if the line is to be skipped.
 */
int tsv_setter_chrom(tsv_t *tsv, bcf1_t *rec, void *usr);
int tsv_setter_pos(tsv_t *tsv, bcf1_t *rec, void *usr);
int tsv_setter_id(tsv_t *tsv, bcf1_t *rec, void *usr);

#endif

