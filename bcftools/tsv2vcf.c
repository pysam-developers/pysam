/*  tsv2vcf.c -- convert from whitespace-separated fields to VCF

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
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
*/

#include <ctype.h>
#include <strings.h>
#include "tsv2vcf.h"

tsv_t *tsv_init(const char *str)
{
    tsv_t *tsv = (tsv_t *) calloc(1,sizeof(tsv_t));
    kstring_t tmp = {0,0,0};
    const char *ss = str, *se = ss;
    tsv->ncols = 0;
    while ( *ss )
    {
        if ( *se && *se!=',' ) { se++; continue; }
        tsv->ncols++;
        tsv->cols = (tsv_col_t*) realloc(tsv->cols,sizeof(tsv_col_t)*tsv->ncols);
        tsv->cols[tsv->ncols-1].name   = NULL;
        tsv->cols[tsv->ncols-1].setter = NULL;
        tmp.l = 0;
        kputsn(ss, se-ss, &tmp);
        if ( strcasecmp("-",tmp.s) )
            tsv->cols[tsv->ncols-1].name = strdup(tmp.s);
        if ( !*se ) break;
        ss = ++se;
    }
    free(tmp.s);
    return tsv;
}

void tsv_destroy(tsv_t *tsv)
{
    int i;
    for (i=0; i<tsv->ncols; i++) free(tsv->cols[i].name);
    free(tsv->cols);
    free(tsv);
}

int tsv_register(tsv_t *tsv, const char *id, tsv_setter_t setter, void *usr)
{
    int i;
    for (i=0; i<tsv->ncols; i++)
    {
        if ( !tsv->cols[i].name || strcasecmp(tsv->cols[i].name,id) ) continue;
        tsv->cols[i].setter = setter;
        tsv->cols[i].usr    = usr;
        return 0;
    }
    return -1;
}

int tsv_parse(tsv_t *tsv, bcf1_t *rec, char *str)
{
    int status = 0;
    tsv->icol = 0;
    tsv->ss = tsv->se = str;
    while ( *tsv->ss && tsv->icol < tsv->ncols )
    {
        while ( *tsv->se && !isspace(*tsv->se) ) tsv->se++;
        if ( tsv->cols[tsv->icol].setter )
        {
            int ret = tsv->cols[tsv->icol].setter(tsv,rec,tsv->cols[tsv->icol].usr);
            if ( ret<0 ) return -1;
            status++;
        }
        while ( *tsv->se && isspace(*tsv->se) ) tsv->se++;
        tsv->ss = tsv->se;
        tsv->icol++;
    }
    return status ? 0 : -1;
}

int tsv_setter_chrom(tsv_t *tsv, bcf1_t *rec, void *usr)
{
    char tmp = *tsv->se;
    *tsv->se = 0;
    rec->rid = bcf_hdr_name2id((bcf_hdr_t*)usr, tsv->ss);
    *tsv->se = tmp;
    return rec->rid==-1 ? -1 : 0;
}

int tsv_setter_pos(tsv_t *tsv, bcf1_t *rec, void *usr)
{
    char *endptr;
    rec->pos = strtol(tsv->ss, &endptr, 10) - 1;
    if ( tsv->ss==endptr ) return -1;
    return 0;
}

int tsv_setter_id(tsv_t *tsv, bcf1_t *rec, void *usr)
{
    char tmp = *tsv->se;
    *tsv->se = 0;
    bcf_update_id((bcf_hdr_t*)usr, rec, tsv->ss);
    *tsv->se = tmp;
    return 0;
}

int tsv_setter_ref_alt(tsv_t *tsv, bcf1_t *rec, void *usr)
{
    bcf_hdr_t *hdr = (bcf_hdr_t*)usr;
    char *sb = tsv->ss;
    while ( *sb && !isspace(*sb) ) sb++;
    if ( !*sb ) return -1;
    char tmp = *sb;
    *sb = ',';
    bcf_update_alleles_str(hdr, rec, tsv->ss);
    *sb = tmp;
    return 0;
}


