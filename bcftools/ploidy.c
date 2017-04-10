/* 
    Copyright (C) 2014-2016 Genome Research Ltd.

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

#include <htslib/khash_str2int.h>
#include <htslib/kseq.h>
#include <htslib/hts.h>
#include "bcftools.h"
#include "ploidy.h"

struct _ploidy_t
{
    int nsex, msex;     // number of genders, m:number of allocated elements in id2sex array
    int dflt, min, max; // ploidy: default, min and max (only explicitly listed)
    int *sex2dflt;
    regidx_t *idx;
    regitr_t *itr;
    void *sex2id;
    char **id2sex;
    kstring_t tmp_str;
};

typedef struct
{
    int sex, ploidy;
}
sex_ploidy_t;


regidx_t *ploidy_regions(ploidy_t *ploidy)
{
    return ploidy->idx;
}

int ploidy_parse(const char *line, char **chr_beg, char **chr_end, uint32_t *beg, uint32_t *end, void *payload, void *usr)
{
    int i, ret;
    ploidy_t *ploidy = (ploidy_t*) usr;
    void *sex2id = ploidy->sex2id;

    // Check for special case of default ploidy "* * * <sex> <ploidy>"
    int default_ploidy_def = 0;

    char *ss = (char*) line;
    while ( *ss && isspace(*ss) ) ss++;
    if ( ss[0]=='*' && (!ss[1] || isspace(ss[1])) )
        default_ploidy_def = 1; // definition of default ploidy, chr="*"
    else
    {
        // Fill CHR,FROM,TO
        ret = regidx_parse_tab(line,chr_beg,chr_end,beg,end,NULL,NULL);
        if ( ret!=0 ) return ret;
    }

    // Skip the fields already parsed by regidx_parse_tab
    ss = (char*) line;
    while ( *ss && isspace(*ss) ) ss++;
    for (i=0; i<3; i++)
    {
        while ( *ss && !isspace(*ss) ) ss++;
        if ( !*ss ) return -2;  // wrong number of fields
        while ( *ss && isspace(*ss) ) ss++;
    }
    if ( !*ss ) return -2;

    // Parse the payload
    char *se = ss;
    while ( *se && !isspace(*se) ) se++;
    if ( !*se || se==ss ) error("Could not parse: %s\n", line);
    ploidy->tmp_str.l = 0;
    kputsn(ss,se-ss,&ploidy->tmp_str);

    sex_ploidy_t *sp = (sex_ploidy_t*) payload;
    if ( khash_str2int_get(sex2id, ploidy->tmp_str.s, &sp->sex) != 0 )
    {
        ploidy->nsex++;
        hts_expand0(char*,ploidy->nsex,ploidy->msex,ploidy->id2sex);
        ploidy->id2sex[ploidy->nsex-1] = strdup(ploidy->tmp_str.s);
        sp->sex = khash_str2int_inc(ploidy->sex2id, ploidy->id2sex[ploidy->nsex-1]);
        ploidy->sex2dflt = (int*) realloc(ploidy->sex2dflt,sizeof(int)*ploidy->nsex);
        ploidy->sex2dflt[ploidy->nsex-1] = -1;
    }

    ss = se;
    while ( *se && isspace(*se) ) se++;
    if ( !*se ) error("Could not parse: %s\n", line);
    sp->ploidy = strtol(ss,&se,10);
    if ( ss==se ) error("Could not parse: %s\n", line);
    if ( ploidy->min<0 || sp->ploidy < ploidy->min ) ploidy->min = sp->ploidy;
    if ( ploidy->max<0 || sp->ploidy > ploidy->max ) ploidy->max = sp->ploidy;

    // Special case, chr="*" stands for a default value
    if ( default_ploidy_def )
    {
        ploidy->sex2dflt[ploidy->nsex-1] = sp->ploidy;
        return -1;
    }

    return 0;
}

static void _set_defaults(ploidy_t *ploidy, int dflt)
{
    int i;
    if ( khash_str2int_get(ploidy->sex2id, "*", &i) == 0 ) dflt = ploidy->sex2dflt[i];
    for (i=0; i<ploidy->nsex; i++)
        if ( ploidy->sex2dflt[i]==-1 ) ploidy->sex2dflt[i] = dflt;

    ploidy->dflt = dflt;
    if ( ploidy->min<0 || dflt < ploidy->min ) ploidy->min = dflt;
    if ( ploidy->max<0 || dflt > ploidy->max ) ploidy->max = dflt;
}

ploidy_t *ploidy_init(const char *fname, int dflt)
{
    ploidy_t *pld = (ploidy_t*) calloc(1,sizeof(ploidy_t));
    if ( !pld ) return NULL;

    pld->min = pld->max = -1;
    pld->sex2id = khash_str2int_init();
    pld->idx = regidx_init(fname,ploidy_parse,NULL,sizeof(sex_ploidy_t),pld);
    if ( !pld->idx )
    {
        ploidy_destroy(pld);
        return NULL;
    }
    pld->itr = regitr_init(pld->idx);
    _set_defaults(pld,dflt);
    return pld;
}

ploidy_t *ploidy_init_string(const char *str, int dflt)
{
    ploidy_t *pld = (ploidy_t*) calloc(1,sizeof(ploidy_t));
    if ( !pld ) return NULL;

    pld->min = pld->max = -1;
    pld->sex2id = khash_str2int_init();
    pld->idx = regidx_init(NULL,ploidy_parse,NULL,sizeof(sex_ploidy_t),pld);
    pld->itr = regitr_init(pld->idx);

    kstring_t tmp = {0,0,0};
    const char *ss = str;
    while ( *ss )
    {
        while ( *ss && isspace(*ss) ) ss++;
        const char *se = ss;
        while ( *se && *se!='\r' && *se!='\n' ) se++;
        tmp.l = 0;
        kputsn(ss, se-ss, &tmp);
        regidx_insert(pld->idx,tmp.s);
        while ( *se && isspace(*se) ) se++;
        ss = se;
    }
    free(tmp.s);

    _set_defaults(pld,dflt);
    return pld;
}

void ploidy_destroy(ploidy_t *ploidy)
{
    if ( ploidy->sex2id ) khash_str2int_destroy_free(ploidy->sex2id);
    if ( ploidy->itr ) regitr_destroy(ploidy->itr);
    if ( ploidy->idx ) regidx_destroy(ploidy->idx);
    free(ploidy->id2sex);
    free(ploidy->tmp_str.s);
    free(ploidy->sex2dflt);
    free(ploidy);
}

int ploidy_query(ploidy_t *ploidy, char *seq, int pos, int *sex2ploidy, int *min, int *max)
{
    int i, ret = regidx_overlap(ploidy->idx, seq,pos,pos, ploidy->itr);

    if ( !sex2ploidy && !min && !max ) return ret;

    if ( !ret )
    {
        // no overlap
        if ( min ) *min = ploidy->dflt;
        if ( max ) *max = ploidy->dflt;
        if ( sex2ploidy )
            for (i=0; i<ploidy->nsex; i++) sex2ploidy[i] = ploidy->sex2dflt[i];
        return 0;
    }

    int _min = INT_MAX, _max = -1;
    if ( sex2ploidy ) for (i=0; i<ploidy->nsex; i++) sex2ploidy[i] = ploidy->dflt;

    while ( regitr_overlap(ploidy->itr) )
    {
        int sex = regitr_payload(ploidy->itr,sex_ploidy_t).sex;
        int pld = regitr_payload(ploidy->itr,sex_ploidy_t).ploidy;
        if ( pld!=ploidy->dflt ) 
        {
            if ( sex2ploidy ) sex2ploidy[ sex ] = pld;
            if ( _min > pld ) _min = pld;
            if ( _max < pld ) _max = pld;
        }
    }
    if ( _max==-1 ) _max = _min = ploidy->dflt;
    if ( max ) *max = _max;
    if ( min ) *min = _min;

    return 1;
}

int ploidy_nsex(ploidy_t *ploidy)
{
    return ploidy->nsex;
}

char *ploidy_id2sex(ploidy_t *ploidy, int id)
{
    if ( id<0 || id>=ploidy->nsex ) return NULL;
    return ploidy->id2sex[id];
}

int ploidy_sex2id(ploidy_t *ploidy, char *sex)
{
    int id;
    if ( khash_str2int_get(ploidy->sex2id,sex,&id)!=0 ) return -1;
    return id;
}

int ploidy_add_sex(ploidy_t *ploidy, const char *sex)
{
    int id;
    if ( khash_str2int_get(ploidy->sex2id, sex, &id)==0 ) return id;
    ploidy->nsex++;
    hts_expand0(char*,ploidy->nsex,ploidy->msex,ploidy->id2sex);
    ploidy->id2sex[ploidy->nsex-1] = strdup(sex);
    ploidy->sex2dflt = (int*) realloc(ploidy->sex2dflt,sizeof(int)*ploidy->nsex);
    ploidy->sex2dflt[ploidy->nsex-1] = ploidy->dflt;
    return khash_str2int_inc(ploidy->sex2id, ploidy->id2sex[ploidy->nsex-1]);
}

int ploidy_max(ploidy_t *ploidy)
{
    return ploidy->dflt > ploidy->max ? ploidy->dflt : ploidy->max;
}

int ploidy_min(ploidy_t *ploidy)
{
    return ploidy->dflt < ploidy->min ? ploidy->dflt : ploidy->min;
}

