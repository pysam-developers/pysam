#include "bcftools.pysam.h"

/* 
    Copyright (C) 2016-2021 Genome Research Ltd.

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
#include "bcftools.h"
#include "smpl_ilist.h"

void smpl_ilist_destroy(smpl_ilist_t *smpl)
{
    int i;
    if ( smpl->pair )
    {
        for (i=0; i<smpl->n; i++) free(smpl->pair[i]);
        free(smpl->pair);
    }
    free(smpl->idx);
    free(smpl);
}

static inline int is_space_or_escaped(const char *min, const char *str)
{
    if ( !isspace(*str) ) return 0;
    int n = 0;
    while ( --str>=min && *str=='\\' ) n++;
    return n%2 ? 0 : 1;
}
smpl_ilist_t *smpl_ilist_init(bcf_hdr_t *hdr, char *sample_list, int is_file, int flags)
{
    smpl_ilist_t *smpl = (smpl_ilist_t*) calloc(1,sizeof(smpl_ilist_t));

    int i;
    if ( !sample_list )
    {
        smpl->n = bcf_hdr_nsamples(hdr);
        smpl->idx = (int*) malloc(sizeof(int)*smpl->n);
        for (i=0; i<smpl->n; i++) smpl->idx[i] = i;
        return smpl;
    }

    int negate = sample_list[0]=='^' ? 1 : 0;
    int nlist;
    char **list = hts_readlist(negate?sample_list+1:sample_list, is_file, &nlist);
    if ( !list ) error("Could not parse %s\n", sample_list);

    if ( negate && (flags&SMPL_REORDER) ) flags &= ~SMPL_REORDER;

    // preserve the VCF order unless flags&SMPL_REORDER is set
    int j = 0;
    int *tmp = (int*)calloc(bcf_hdr_nsamples(hdr),sizeof(int));
    char **pair = NULL;
    for (i=0; i<nlist; i++)
    {
        char *smpl1 = list[i];
        char *smpl2 = NULL;

        char *ptr = list[i];
        while ( *ptr && !is_space_or_escaped(list[i], ptr) ) ptr++;
        if ( *ptr )
        {
            *ptr = 0;
            smpl2 = ptr+1;
        }

        char *smpl_name = flags&SMPL_PAIR2 && smpl2 ? smpl2 : smpl1;
        int idx = bcf_hdr_id2int(hdr, BCF_DT_SAMPLE, smpl_name);
        if ( idx<0 )
        {
            if ( !(flags&SMPL_STRICT) )
            {
                if ( flags&SMPL_VERBOSE ) fprintf(bcftools_stderr,"No such sample: \"%s\"\n",smpl_name);
                continue;
            }
            error("No such sample: \"%s\"\n", smpl_name);
        }

        if ( flags & SMPL_REORDER )
        {
            tmp[j++] = idx;
        }
        else
        {
            tmp[idx] = 1;
            if ( smpl2 )
            {
                if ( !pair ) pair = (char**)calloc(bcf_hdr_nsamples(hdr),sizeof(char*));
                if ( flags&SMPL_PAIR2 ) pair[idx] = strdup(smpl1);
                else if ( flags&SMPL_PAIR1 ) pair[idx] = strdup(smpl2);
            }
        }
        smpl->n++;
    }

    if ( flags & SMPL_REORDER )
    {
        smpl->idx = tmp;
        for (i=0; i<nlist; i++) free(list[i]);
        free(list);
        return smpl;
    }

    if ( negate ) smpl->n = bcf_hdr_nsamples(hdr) - smpl->n;
    smpl->idx = (int*) malloc(sizeof(int)*smpl->n);

    j = 0;
    if ( !negate )
    {
        if ( pair ) smpl->pair = (char**) calloc(smpl->n,sizeof(char*));
        for (i=0; i<bcf_hdr_nsamples(hdr); i++)
        {
            if ( !tmp[i] ) continue;
            smpl->idx[j] = i;
            if ( pair && pair[i] ) smpl->pair[j] = pair[i];
            j++;
        }
    }
    else
    {
        for (i=0; i<bcf_hdr_nsamples(hdr); i++)
            if ( !tmp[i] ) smpl->idx[j++] = i;
    }

    free(tmp);
    free(pair);
    for (i=0; i<nlist; i++) free(list[i]);
    free(list);

    return smpl;
}

smpl_ilist_t *smpl_ilist_map(bcf_hdr_t *hdr_a, bcf_hdr_t *hdr_b, int flags)
{
    if ( flags&SMPL_STRICT && bcf_hdr_nsamples(hdr_a)!=bcf_hdr_nsamples(hdr_b) )
        error("Different number of samples: %d vs %d\n", bcf_hdr_nsamples(hdr_a),bcf_hdr_nsamples(hdr_b));

    smpl_ilist_t *smpl = (smpl_ilist_t*) calloc(1,sizeof(smpl_ilist_t));

    int i;
    smpl->n = bcf_hdr_nsamples(hdr_a);
    smpl->idx = (int*) malloc(sizeof(int)*smpl->n);
    for (i=0; i<smpl->n; i++)
    {
        const char *name = bcf_hdr_int2id(hdr_a, BCF_DT_SAMPLE, i);
        smpl->idx[i] = bcf_hdr_id2int(hdr_b, BCF_DT_SAMPLE, name);
        if ( flags&SMPL_STRICT && smpl->idx[i]<0 ) 
            error("The sample %s is not present in the second file\n", name);
    }
    return smpl;
}

