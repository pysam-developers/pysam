#include "bcftools.pysam.h"

/*  gvcf.c -- support for gVCF files.

    Copyright (C) 2014-2015 Genome Research Ltd.

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

#include "gvcf.h"
#include "bcftools.h"

struct _gvcf_t
{
    int *dp_range, ndp_range;   // per-sample DP ranges
    int prev_range;             // 0 if not in a block
    int32_t *dp, mdp, *pl, mpl, npl;
    int32_t *tmp, mtmp, *gts, ngts,mgts, nqsum,mqsum;
    float *qsum;
    int32_t rid, start, end, min_dp;
    kstring_t als;
    bcf1_t *line;
};

void gvcf_update_header(gvcf_t *gvcf, bcf_hdr_t *hdr)
{
    bcf_hdr_append(hdr,"##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">");
    bcf_hdr_append(hdr,"##INFO=<ID=MinDP,Number=1,Type=Integer,Description=\"Minimum per-sample depth in this gVCF block\">");
}

gvcf_t *gvcf_init(const char *dp_ranges)
{
    gvcf_t *gvcf = (gvcf_t*) calloc(1,sizeof(gvcf_t));
    gvcf->line = bcf_init();

    int n = 1;
    const char *ss = dp_ranges;
    while ( *ss )
    {
        if ( *ss==',' ) n++;
        ss++;
    }
    gvcf->ndp_range = n;
    gvcf->dp_range  = (int*) malloc(sizeof(int)*gvcf->ndp_range);

    n  = 0;
    ss = dp_ranges;
    while ( *ss )
    {
        char *se = (char*) ss;
        gvcf->dp_range[n++] = strtol(ss,&se,10);
        if ( se==ss ) return NULL;
        if ( *se==',' && se[1] ) { ss = se+1; continue; }
        else if ( !*se ) break;
        return NULL;
    }
    return gvcf;
}

void gvcf_destroy(gvcf_t *gvcf)
{
    free(gvcf->dp_range);
    free(gvcf->dp);
    free(gvcf->pl);
    free(gvcf->tmp);
    free(gvcf->qsum);
    free(gvcf->gts);
    free(gvcf->als.s);
    if ( gvcf->line ) bcf_destroy(gvcf->line);
    free(gvcf);
}

bcf1_t *gvcf_write(gvcf_t *gvcf, htsFile *fh, bcf_hdr_t *hdr, bcf1_t *rec, int is_ref)
{
    int i, ret, nsmpl = bcf_hdr_nsamples(hdr);
    int can_collapse = is_ref ? 1 : 0;
    int32_t dp_range = 0, min_dp = 0;

    // No record and nothing to flush?
    if ( !rec && !gvcf->prev_range ) return NULL;

    // Flush gVCF block if there are no more records, chr changed, a gap
    // encountered, or other conditions not met (block broken by a non-ref or DP too low).
    int needs_flush = can_collapse ? 0 : 1;


    // Can the record be included in a gVCF block? That is, is this a ref-only site?
    if ( rec && can_collapse )
    {
        bcf_unpack(rec, BCF_UN_ALL);

        // per-sample depth
        ret = bcf_get_format_int32(hdr, rec, "DP", &gvcf->tmp, &gvcf->mtmp);
        if ( ret==nsmpl )
        {
            min_dp = gvcf->tmp[0];
            for (i=1; i<nsmpl; i++)
                if ( min_dp > gvcf->tmp[i] ) min_dp = gvcf->tmp[i];

            for (i=0; i<gvcf->ndp_range; i++)
                if ( min_dp < gvcf->dp_range[i] ) break;

            dp_range = i;
            if ( !dp_range )
            {
                // leave the record unchanged, DP is too small. Alternatively, return NULL here
                // to skip these sites
                needs_flush  = 1;
                can_collapse = 0;
            }
        }
        else
            needs_flush = 1;       // DP field not present
    }

    if ( gvcf->prev_range && gvcf->prev_range!=dp_range ) needs_flush = 1;
    if ( !rec || gvcf->rid!=rec->rid || rec->pos > gvcf->end+1 ) needs_flush = 1;

    // If prev_range is set, something can be flushed
    if ( gvcf->prev_range && needs_flush )
    {
        // mpileup can output two records with the same position, SNP and
        // indel. Make sure the end position does not include the non-variant
        // SNP position just before the indel.
        if ( rec && rec->rid==gvcf->rid && rec->pos==gvcf->end ) gvcf->end--;

        gvcf->end++;    // from 0-based to 1-based coordinate

        bcf_clear1(gvcf->line);
        gvcf->line->rid  = gvcf->rid;
        gvcf->line->pos  = gvcf->start;
        gvcf->line->rlen = gvcf->end - gvcf->start;
        bcf_update_alleles_str(hdr, gvcf->line, gvcf->als.s);
        if ( gvcf->start+1 < gvcf->end )    // create gVCF record only if it spans at least two sites
            bcf_update_info_int32(hdr, gvcf->line, "END", &gvcf->end, 1);
        bcf_update_info_int32(hdr, gvcf->line, "MinDP", &gvcf->min_dp, 1);
        if ( gvcf->nqsum>0 )
            bcf_update_info_float(hdr, gvcf->line, "QS", gvcf->qsum, gvcf->nqsum);
        if ( gvcf->ngts )
            bcf_update_genotypes(hdr,gvcf->line,gvcf->gts,gvcf->ngts);
        if ( gvcf->npl>0 )
            bcf_update_format_int32(hdr, gvcf->line, "PL", gvcf->pl, gvcf->npl);
        bcf_update_format_int32(hdr, gvcf->line, "DP", gvcf->dp, nsmpl);
        if ( bcf_write1(fh, hdr, gvcf->line)!=0 ) error("[%s] Error: failed to write the record\n", __func__);
        gvcf->prev_range = 0;
        gvcf->rid  = -1;
        gvcf->npl  = 0;
        gvcf->nqsum = 0;
        gvcf->ngts  = 0;

        if ( !rec ) return NULL;     // just flushing the buffer, this was last record
    }

    if ( can_collapse )
    {
        if ( !gvcf->prev_range )
        {
            hts_expand(int32_t,nsmpl,gvcf->mdp,gvcf->dp);
            memcpy(gvcf->dp,gvcf->tmp,nsmpl*sizeof(int32_t));   // tmp still contains DP from rec
            gvcf->npl = bcf_get_format_int32(hdr, rec, "PL", &gvcf->pl, &gvcf->mpl);

            gvcf->nqsum = bcf_get_info_float(hdr,rec,"QS",&gvcf->qsum,&gvcf->mqsum);
            gvcf->ngts  = bcf_get_genotypes(hdr,rec,&gvcf->gts,&gvcf->mgts);

            gvcf->rid    = rec->rid;
            gvcf->start  = rec->pos;
            gvcf->als.l = 0;
            kputs(rec->d.allele[0],&gvcf->als);
            for (i=1; i<rec->n_allele; i++)
            {
                kputc(',',&gvcf->als);
                kputs(rec->d.allele[i],&gvcf->als);
            }
            gvcf->min_dp = min_dp;
        }
        else
        {
            if ( gvcf->min_dp > min_dp ) gvcf->min_dp = min_dp;
            for (i=0; i<nsmpl; i++)
                if ( gvcf->dp[i] > gvcf->tmp[i] ) gvcf->dp[i] = gvcf->tmp[i];
            ret = bcf_get_format_int32(hdr, rec, "PL", &gvcf->tmp, &gvcf->mtmp);
            if ( ret>=0 )
            {
                if ( ret!=nsmpl*3 ) error("Unexpected number of PL fields\n");
                for (i=0; i<nsmpl; i++)
                {
                    if ( gvcf->pl[3*i+1] > gvcf->tmp[3*i+1] )
                    {
                        gvcf->pl[3*i+1] = gvcf->tmp[3*i+1];
                        gvcf->pl[3*i+2] = gvcf->tmp[3*i+2];
                    }
                    else if ( gvcf->pl[3*i+1]==gvcf->tmp[3*i+1] && gvcf->pl[3*i+2] > gvcf->tmp[3*i+2] )
                        gvcf->pl[3*i+2] = gvcf->tmp[3*i+2];
                }
            }
            else
                gvcf->npl = 0;
        }
        gvcf->prev_range = dp_range;
        if ( bcf_get_info_int32(hdr,rec,"END",&gvcf->tmp,&gvcf->mtmp)==1 )
            gvcf->end = gvcf->tmp[0] - 1;   // from 1-based to 0-based
        else
            gvcf->end = rec->pos;
        return NULL;
    }

    if ( is_ref && min_dp )
        bcf_update_info_int32(hdr, rec, "MinDP", &min_dp, 1);

    return rec;
}

