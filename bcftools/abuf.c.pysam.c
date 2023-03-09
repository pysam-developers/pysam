#include "bcftools.pysam.h"

/* The MIT License

   Copyright (c) 2021-2023 Genome Research Ltd.

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

#include <assert.h>
#include <strings.h>
#include <htslib/vcf.h>
#include <ctype.h>
#include "bcftools.h"
#include "abuf.h"
#include "rbuf.h"

typedef enum
{
    M_FIRST, M_SUM
}
merge_rule_t;

typedef struct
{
    kstring_t ref, alt;
    int ial;        // the index of the original ALT allele, 1-based
    int beg, end;   // 0-based inclusive offsets to ref,alt
}
atom_t;

typedef struct
{
    bcf1_t *rec;
    int nori, nout;     // number of ALTs in the input, and VCF rows on output
    uint8_t *tbl;       // nori columns, nout rows; indicates allele contribution to output rows, see "The atomization works as follows" below
    uint8_t *overlaps;  // is the star allele needed for this variant?
    atom_t **atoms;
    int matoms, mtbl, moverlaps;
    char *info_tag;
}
split_t;

struct _abuf_t
{
    abuf_opt_t mode;
    split_t split;
    atom_t *atoms;
    int natoms, matoms;
    const bcf_hdr_t *hdr;
    bcf_hdr_t *out_hdr;
    bcf1_t **vcf;       // dimensions stored in rbuf
    rbuf_t rbuf;

    kstring_t tmps;
    void *tmp, *tmp2;
    int32_t *gt, *tmpi;
    int ngt, mgt, ntmpi, mtmpi, mtmp, mtmp2;
    int star_allele;
};

abuf_t *abuf_init(const bcf_hdr_t *hdr, abuf_opt_t mode)
{
    if ( mode!=SPLIT ) error("todo\n");
    abuf_t *buf = (abuf_t*) calloc(1,sizeof(abuf_t));
    buf->hdr  = hdr;
    buf->out_hdr = (bcf_hdr_t*) hdr;
    buf->mode = mode;
    buf->star_allele = 1;
    rbuf_init(&buf->rbuf, 0);
    return buf;
}

void abuf_destroy(abuf_t *buf)
{
    int i;
    for (i=0; i<buf->matoms; i++)
    {
        free(buf->atoms[i].ref.s);
        free(buf->atoms[i].alt.s);
    }
    free(buf->atoms);
    free(buf->split.atoms);
    free(buf->split.overlaps);
    free(buf->split.tbl);
    for (i=0; i<buf->rbuf.m; i++)
        if ( buf->vcf[i] ) bcf_destroy(buf->vcf[i]);
    free(buf->vcf);
    free(buf->gt);
    free(buf->tmpi);
    free(buf->tmp);
    free(buf->tmp2);
    free(buf->tmps.s);
    free(buf);
}

void abuf_set(abuf_t *buf, abuf_opt_t key, void *value)
{
    if ( key==BCF_HDR ) { buf->out_hdr = *((bcf_hdr_t**)value); return; }
    if ( key==INFO_TAG )
    {
        buf->split.info_tag = *((char**)value);
        bcf_hdr_printf(buf->out_hdr,"##INFO=<ID=%s,Number=1,Type=String,Description=\"Original variant. Format: CHR|POS|REF|ALT|USED_ALT_IDX\">",buf->split.info_tag);
        return;
    }
    if ( key==STAR_ALLELE ) { buf->star_allele = *((int*)value); return; }
}

/*
    Split alleles into primitivs, e.g.
        CC>TT  becomes  C>T,C>T
        GCGT>GTGA  becomes C>T,T>A

    There is no sequence alignment, just trimming and hungry matching
    from left side.
*/
static void _atomize_allele(abuf_t *buf, bcf1_t *rec, int ial)
{
    // Trim identical sequence from right
    char *ref = rec->d.allele[0];
    char *alt = rec->d.allele[ial];
    int rlen = strlen(ref);
    int alen = strlen(alt);
    while ( rlen>1 && alen>1 && ref[rlen-1]==alt[alen-1] ) rlen--, alen--;
    int Mlen = rlen > alen ? rlen : alen;

    atom_t *atom = NULL;
    int i;
    for (i=0; i<Mlen; i++)
    {
        char refb = i<rlen ? ref[i] : '-';
        char altb = i<alen ? alt[i] : '-';
        if ( refb!=altb )
        {
            if ( refb=='-' || altb=='-' )
            {
                assert(atom);
                if ( altb!='-' ) kputc(altb, &atom->alt);
                if ( refb!='-' ) { kputc(refb, &atom->ref); atom->end++; }
                continue;
            }
            buf->natoms++;
            hts_expand0(atom_t,buf->natoms,buf->matoms,buf->atoms);
            atom = &buf->atoms[buf->natoms-1];
            atom->ref.l = 0;
            atom->alt.l = 0;
            kputc(refb, &atom->ref);
            kputc(altb, &atom->alt);
            atom->beg = atom->end = i;
            atom->ial = ial;

            if ( rlen!=alen && (i+1>=rlen || i+1>=alen) )   // the next base is an indel combined with SNV, e.g. C>GGG?
            {
                buf->natoms++;
                hts_expand0(atom_t,buf->natoms,buf->matoms,buf->atoms);
                atom = &buf->atoms[buf->natoms-1];
                atom->ref.l = 0;
                atom->alt.l = 0;
                kputc(refb, &atom->ref);
                kputc(refb, &atom->alt);
                atom->beg = atom->end = i;
                atom->ial = ial;
            }
            continue;
        }
        if ( i+1>=rlen || i+1>=alen )   // is the next base an indel?
        {
            buf->natoms++;
            hts_expand0(atom_t,buf->natoms,buf->matoms,buf->atoms);
            atom = &buf->atoms[buf->natoms-1];
            atom->ref.l = 0;
            atom->alt.l = 0;
            kputc(refb, &atom->ref);
            kputc(altb, &atom->alt);
            atom->beg = atom->end = i;
            atom->ial = ial;
        }
    }
}
static int _atoms_inconsistent(const atom_t *a, const atom_t *b)
{
    if ( a->beg < b->beg ) return -1;
    if ( a->beg > b->beg ) return 1;
    int rcmp = strcasecmp(a->ref.s,b->ref.s);
    if ( rcmp ) return rcmp;
    return strcasecmp(a->alt.s,b->alt.s);
}
/*
    For reproducibility of tests on different platforms, we need to guarantee the same order of identical
    atoms originating from different source ALTs.  Even though they are consistent, different values can be
    picked for VCF annotations as currently the values from the one that comes first are used.
*/
static int _cmp_atoms(const void *aptr, const void *bptr)
{
    const atom_t *a = (const atom_t*) aptr;
    const atom_t *b = (const atom_t*) bptr;
    int rcmp = _atoms_inconsistent(a,b);
    if ( rcmp ) return rcmp;
    if ( a->ial < b->ial ) return -1;
    if ( a->ial > b->ial ) return 1;
    return 0;
}
static void _split_table_init(abuf_t *buf, bcf1_t *rec, int natoms)
{
    buf->split.rec  = rec;
    buf->split.nori = rec->n_allele - 1;
    buf->split.nout = 0;
    hts_expand(uint8_t,buf->split.nori*natoms,buf->split.mtbl,buf->split.tbl);
    hts_expand(atom_t*,natoms,buf->split.matoms,buf->split.atoms);
    hts_expand(uint8_t,natoms,buf->split.moverlaps,buf->split.overlaps);
    memset(buf->split.overlaps,0,sizeof(*buf->split.overlaps)*natoms);
}
static void _split_table_new(abuf_t *buf, atom_t *atom)
{
    int i, iout = buf->split.nout++;
    buf->split.atoms[iout] = atom;
    uint8_t *ptr = buf->split.tbl + iout*buf->split.nori;
    for (i=0; i<buf->split.nori; i++) ptr[i] = 0;
    ptr[atom->ial-1] = 1;
}
static void _split_table_overlap(abuf_t *buf, int iout, atom_t *atom)
{
    uint8_t *ptr = buf->split.tbl + iout*buf->split.nori;
    ptr[atom->ial-1] = _atoms_inconsistent(atom,buf->split.atoms[iout]) ? 2 : 1;
    buf->split.overlaps[iout] = 1;
}
#if 0
static void _split_table_print(abuf_t *buf)
{
    int i,j;
    for (i=0; i<buf->split.nout; i++)
    {
        atom_t *atom = buf->split.atoms[i];
        uint8_t *ptr = buf->split.tbl + i*buf->split.nori;
        fprintf(bcftools_stderr,"%d\t%s\t%s",(int)buf->split.rec->pos+1+atom->beg,atom->ref.s,atom->alt.s);
        for (j=0; j<buf->split.nori; j++) fprintf(bcftools_stderr,"\t%d",(int)ptr[j]);
        fprintf(bcftools_stderr,"\n");
    }
}
static void _split_table_print_atoms(abuf_t *buf)
{
    int i;
    for (i=0; i<buf->natoms; i++)
    {
        atom_t *atom = &buf->atoms[i];
        fprintf(bcftools_stderr,"atom%d %p: ialt=%d %s>%s %d-%d\n",i,atom,atom->ial,atom->ref.s,atom->alt.s,atom->beg,atom->end);
    }
}
#endif
static inline uint8_t _has_star_allele(abuf_t *buf, int iout)
{
    if ( !buf->star_allele ) return 0;
    return buf->split.overlaps[iout];
}
static inline int _split_table_get_ial(abuf_t *buf, int irow, int ial)
{
    if ( !ial ) return ial;
    return buf->split.tbl[irow*buf->split.nori + ial - 1];
}
static void _split_table_set_chrom_qual(abuf_t *buf)
{
    int iout,j;
    bcf1_t *rec = buf->split.rec;
    for (iout=0; iout<buf->split.nout; iout++)
    {
        rbuf_expand0(&buf->rbuf, bcf1_t*, buf->rbuf.n+1, buf->vcf);
        j = rbuf_append(&buf->rbuf);
        if ( !buf->vcf[j] ) buf->vcf[j] = bcf_init1();
        bcf1_t *out = buf->vcf[j];
        bcf_clear1(out);

        atom_t *atom = buf->split.atoms[iout];
        out->rid = rec->rid;
        out->pos = rec->pos + atom->beg;
        bcf_update_id(buf->out_hdr, out, rec->d.id);

        const char *als[3];
        als[0] = atom->ref.s;
        als[1] = atom->alt.s;
        als[2] = "*";
        int nals = _has_star_allele(buf,iout) ? 3 : 2;
        bcf_update_alleles(buf->out_hdr, out, als, nals);

        if ( bcf_float_is_missing(rec->qual) )
            bcf_float_set_missing(out->qual);
        else
            out->qual = rec->qual;

        bcf_update_filter(buf->out_hdr, out, rec->d.flt, rec->d.n_flt);
    }
}
int copy_string_field(char *src, int isrc, int src_len, kstring_t *dst, int idst);
static void _split_table_set_info(abuf_t *buf, bcf_info_t *info, merge_rule_t mode)
{
    const char *tag = bcf_hdr_int2id(buf->hdr,BCF_DT_ID,info->key);
    int type = bcf_hdr_id2type(buf->hdr,BCF_HL_INFO,info->key);
    int len  = bcf_hdr_id2length(buf->hdr,BCF_HL_INFO,info->key);
    if ( len==BCF_VL_G ) return;                                                // todo: Number=G INFO tags
    if ( type==BCF_HT_LONG ) return;                                            // todo: 64bit integers

    bcf1_t *rec = buf->split.rec;
    int mtmp = ( type==BCF_HT_INT || type==BCF_HT_REAL ) ? buf->mtmp/4 : buf->mtmp;
    int nval = bcf_get_info_values(buf->hdr,rec,tag,&buf->tmp,&mtmp,type);
    if ( type==BCF_HT_INT || type==BCF_HT_REAL ) buf->mtmp = mtmp*4;

    // Check for incorrect number of values. Note this check does not consider all values missing
    // and will remove annotations that don't pass.
    if ( type==BCF_HT_INT || type==BCF_HT_REAL )
    {
        if ( (len==BCF_VL_A && nval != rec->n_allele - 1) || (len==BCF_VL_R && nval != rec->n_allele) ) return;
    }

    if ( buf->mtmp2 < buf->mtmp )
    {
        buf->tmp2  = realloc(buf->tmp2, buf->mtmp);
        if ( !buf->tmp2 ) error("Failed to alloc %d bytes\n", buf->mtmp);
        buf->mtmp2 = buf->mtmp;
    }

    const int num_size = 4;
    assert( num_size==sizeof(int32_t) && num_size==sizeof(float) );
    int32_t missing = bcf_int32_missing;
    void *missing_ptr = (void*)&missing;
    if ( type==BCF_HT_REAL ) bcf_float_set_missing(*((float*)missing_ptr));
    int32_t vector_end = bcf_int32_vector_end;
    void *vector_end_ptr = (void*)&vector_end;
    if ( type==BCF_HT_REAL ) bcf_float_set_vector_end(*((float*)vector_end_ptr));

    int iout,i;
    for (iout=0; iout<buf->split.nout; iout++)
    {
        bcf1_t *out = buf->vcf[rbuf_kth(&buf->rbuf,iout)];
        int star_allele = _has_star_allele(buf,iout);
        int ret = 0;
        if ( len==BCF_VL_FIXED || len==BCF_VL_VAR )
            ret = bcf_update_info(buf->out_hdr, out, tag, type==BCF_HT_FLAG ? NULL : buf->tmp, nval, type);
        else if ( len==BCF_VL_A && type!=BCF_HT_STR )
        {
            int iori = buf->split.atoms[iout]->ial - 1;
            assert( iori<nval );
            if ( !memcmp(vector_end_ptr,buf->tmp+num_size*iori,num_size) )
                memcpy(buf->tmp2,missing_ptr,num_size);
            else
                memcpy(buf->tmp2,buf->tmp+num_size*iori,num_size);
            if ( star_allele )
                memcpy(buf->tmp2+num_size,missing_ptr,num_size);
            ret = bcf_update_info(buf->out_hdr, out, tag, buf->tmp2, 1 + star_allele, type);
        }
        else if ( len==BCF_VL_A && type==BCF_HT_STR )
        {
            int iori = buf->split.atoms[iout]->ial - 1;
            kstring_t dst;
            dst.l = 0; dst.m = buf->mtmp2; dst.s = (char*)buf->tmp2;
            kputc('.',&dst);
            if ( star_allele ) kputs(",.",&dst);
            copy_string_field(buf->tmp, iori, nval, &dst, 0);
            if ( star_allele ) copy_string_field(".", 0, 1, &dst, 1);
            buf->mtmp2 = dst.m;
            buf->tmp2  = dst.s;
            ret = bcf_update_info(buf->out_hdr, out, tag, buf->tmp2, dst.l, type);
        }
        else if ( len==BCF_VL_R && type!=BCF_HT_STR )
        {
            memcpy(buf->tmp2,buf->tmp,num_size);   // REF contributes to all records
            int iori = buf->split.atoms[iout]->ial;
            assert( iori<nval && iori<=buf->split.nori );
            if ( !memcmp(vector_end_ptr,buf->tmp+num_size*iori,num_size) )
                memcpy(buf->tmp2+num_size,missing_ptr,num_size);
            else
                memcpy(buf->tmp2+num_size,buf->tmp+num_size*iori,num_size);
            if ( type==BCF_HT_INT && mode==M_SUM )
            {
                uint8_t *tbl = buf->split.tbl + iout*buf->split.nori;
                for (i=iori; i<buf->split.nori; i++)
                {
                    if ( tbl[i]==1 ) ((int32_t*)buf->tmp2)[1] += ((int32_t*)buf->tmp)[i+1];
                }
            }
            if ( star_allele )
                memcpy(buf->tmp2+2*num_size,missing_ptr,num_size);
            ret = bcf_update_info(buf->out_hdr, out, tag, buf->tmp2, 2 + star_allele, type);
        }
        else if ( len==BCF_VL_R && type==BCF_HT_STR )
        {
            int iori = buf->split.atoms[iout]->ial - 1;
            kstring_t dst;
            dst.l = 0; dst.m = buf->mtmp2; dst.s = (char*)buf->tmp2;
            kputs(".,.",&dst);
            if ( star_allele ) kputs(",.",&dst);
            copy_string_field(buf->tmp, 0, nval, &dst, 0);
            copy_string_field(buf->tmp, iori+1, nval, &dst, 1);
            if ( star_allele ) copy_string_field(".", 0, 1, &dst, 2);
            buf->mtmp2 = dst.m;
            buf->tmp2  = dst.s;
            ret = bcf_update_info(buf->out_hdr, out, tag, buf->tmp2, dst.l, type);
        }
        if ( ret!=0 ) error("An error occurred while updating INFO/%s\n",tag);
    }
}
static void _split_table_set_history(abuf_t *buf)
{
    int i,j;
    bcf1_t *rec = buf->split.rec;
    buf->tmps.l = 0;
    ksprintf(&buf->tmps,"%s|%"PRIhts_pos"|%s|",bcf_seqname(buf->hdr,rec),rec->pos+1,rec->d.allele[0]);
    for (i=1; i<rec->n_allele; i++)
    {
        kputs(rec->d.allele[i],&buf->tmps);
        if ( i+1<rec->n_allele ) kputc(',',&buf->tmps);
        else kputc(',',&buf->tmps);
    }
    int len = buf->tmps.l;
    buf->tmps.s[buf->tmps.l-1] = '|';

    for (i=0; i<buf->split.nout; i++)
    {
        buf->tmps.l = len;
        bcf1_t *out = buf->vcf[rbuf_kth(&buf->rbuf,i)];
        uint8_t *ptr = buf->split.tbl + i*buf->split.nori;
        for (j=0; j<buf->split.nori; j++)
        {
            if ( ptr[j]!=1 ) continue;
            kputw(j+1,&buf->tmps);
            kputc(',',&buf->tmps);
        }
        buf->tmps.s[--buf->tmps.l] = 0;
        if ( (bcf_update_info_string(buf->out_hdr, out, buf->split.info_tag, buf->tmps.s))!=0 )
            error("An error occurred while updating INFO/%s\n",buf->split.info_tag);
    }
}
static void _split_table_set_gt(abuf_t *buf)
{
    int nsmpl = bcf_hdr_nsamples(buf->hdr);
    if ( !nsmpl ) return;

    bcf1_t *rec = buf->split.rec;
    buf->ngt = bcf_get_genotypes(buf->hdr, rec, &buf->gt, &buf->mgt);
    if ( buf->ngt<=0 ) return;
    else
        hts_expand(int32_t,buf->ngt,buf->mtmpi,buf->tmpi);

    int iout,i,j;
    for (iout=0; iout<buf->split.nout; iout++)
    {
        bcf1_t *out = buf->vcf[rbuf_kth(&buf->rbuf,iout)];
        int star_allele = _has_star_allele(buf,iout);
        int max_ploidy = buf->ngt/nsmpl;
        int32_t *src = buf->gt, *dst = buf->tmpi;
        for (i=0; i<nsmpl; i++)
        {
            for (j=0; j<max_ploidy; j++)
            {
                if ( src[j]==bcf_int32_vector_end || bcf_gt_is_missing(src[j]) )
                {
                    dst[j] = src[j];
                    continue;
                }
                int iori = bcf_gt_allele(src[j]);
                if ( iori<0 || iori>=rec->n_allele )
                    error("Out-of-bounds genotypes at %s:%"PRIhts_pos"\n",bcf_seqname(buf->hdr,rec),rec->pos+1);
                int ial = _split_table_get_ial(buf,iout,iori);
                if ( ial==2 && !star_allele )
                {
                    dst[j] = bcf_gt_missing;
                    if ( bcf_gt_is_phased(src[j]) ) dst[j] |= 1;
                }
                else
                    dst[j] = bcf_gt_is_phased(src[j]) ? bcf_gt_phased(ial) : bcf_gt_unphased(ial);
            }
            src += max_ploidy;
            dst += max_ploidy;
        }
        bcf_update_genotypes(buf->out_hdr,out,buf->tmpi,buf->ngt);
    }
}
static void _split_table_set_format(abuf_t *buf, bcf_fmt_t *fmt, merge_rule_t mode)
{
    int nsmpl = bcf_hdr_nsamples(buf->hdr);
    if ( !nsmpl ) return;

    const char *tag = bcf_hdr_int2id(buf->hdr,BCF_DT_ID,fmt->id);
    if ( tag[0]=='G' && tag[1]=='T' && !tag[2] )        // FORMAT/GT
    {
        _split_table_set_gt(buf);
        return;
    }

    int type = bcf_hdr_id2type(buf->hdr,BCF_HL_FMT,fmt->id);
    int len  = bcf_hdr_id2length(buf->hdr,BCF_HL_FMT,fmt->id);
    if ( type==BCF_HT_STR && len==BCF_VL_G ) return;                            // possible todo: Number=G for strings
    if ( type==BCF_HT_LONG ) return;                                            // todo: 64bit integers

    const int num_size = 4;
    assert( num_size==sizeof(int32_t) && num_size==sizeof(float) );
    int32_t missing = bcf_int32_missing;
    void *missing_ptr = (void*)&missing;
    if ( type==BCF_HT_REAL ) bcf_float_set_missing(*((float*)missing_ptr));
    int32_t vector_end = bcf_int32_vector_end;
    void *vector_end_ptr = (void*)&vector_end;
    if ( type==BCF_HT_REAL ) bcf_float_set_vector_end(*((float*)vector_end_ptr));

    bcf1_t *rec = buf->split.rec;
    int mtmp = ( type==BCF_HT_INT || type==BCF_HT_REAL ) ? buf->mtmp/num_size : buf->mtmp;  // number of items
    int nval = bcf_get_format_values(buf->hdr,rec,tag,&buf->tmp,&mtmp,type);
    if ( type==BCF_HT_INT || type==BCF_HT_REAL ) buf->mtmp = mtmp*num_size;                 // number of bytes

    if ( type==BCF_HT_INT || type==BCF_HT_REAL )
    {
        if ( len==BCF_VL_G && nval!=nsmpl*rec->n_allele && nval!=nsmpl*rec->n_allele*(rec->n_allele+1)/2 ) return;      // not haploid nor diploid

        // Check for incorrect number of values. Note this check does not consider all values missing
        // and will remove annotations that don't pass.
        if ( (len==BCF_VL_A && nval != nsmpl*(rec->n_allele - 1)) || (len==BCF_VL_R && nval != nsmpl*rec->n_allele) ) return;
    }

    // Increase buffer size to accommodate star allele
    int nval1 = nval / nsmpl;
    mtmp = buf->mtmp;
    if ( type==BCF_HT_INT || type==BCF_HT_REAL )
    {
        if ( (len==BCF_VL_A || len==BCF_VL_R) && mtmp < num_size*nsmpl*(nval1+1) ) mtmp = num_size*nsmpl*(nval1+1); // +1 for the possibility of the star allele
        else if ( len==BCF_VL_G && mtmp < num_size*nsmpl*(nval1+3) ) mtmp = num_size*nsmpl*(nval1+3);
    }
    else if ( type==BCF_HT_STR )
    {
        if ( (len==BCF_VL_A || len==BCF_VL_R) && mtmp < nsmpl*(nval1+2) ) mtmp = nsmpl*(nval1+2); // +2 for the possibility of the star allele, ",."
        else if ( len==BCF_VL_G && mtmp < nsmpl*(nval1+6) ) mtmp = nsmpl*(nval1+6);
    }

    if ( buf->mtmp2 < mtmp )
    {
        buf->tmp2  = realloc(buf->tmp2, mtmp);
        if ( !buf->tmp2 ) error("Failed to alloc %d bytes\n", mtmp);
        buf->mtmp2 = mtmp;
    }

    int iout, i, j;
    for (iout=0; iout<buf->split.nout; iout++)
    {
        int star_allele = _has_star_allele(buf,iout);
        bcf1_t *out = buf->vcf[rbuf_kth(&buf->rbuf,iout)];
        int ret = 0;
        if ( len==BCF_VL_FIXED || len==BCF_VL_VAR )
            ret = bcf_update_format(buf->out_hdr, out, tag, buf->tmp, nval, type);
        else if ( len==BCF_VL_A && type!=BCF_HT_STR )
        {
            int iori = buf->split.atoms[iout]->ial - 1;
            assert( iori<nval );
            for (i=0; i<nsmpl; i++)
            {
                void *src = buf->tmp  + nval1*num_size*i;
                void *dst = buf->tmp2 + num_size*i*(star_allele+1);
                if ( !memcmp(vector_end_ptr,src+iori*num_size,num_size) )
                    memcpy(dst,missing_ptr,num_size);
                else
                    memcpy(dst,src+iori*num_size,num_size);
                if ( star_allele )
                    memcpy(dst+num_size,missing_ptr,num_size);
            }
            ret = bcf_update_format(buf->out_hdr, out, tag, buf->tmp2, nsmpl*(star_allele+1), type);
        }
        else if ( (len==BCF_VL_A || len==BCF_VL_R) && type==BCF_HT_STR )
        {
            int ioff = len==BCF_VL_R ? 1 : 0;
            int iori = buf->split.atoms[iout]->ial - 1;
            int nval1_dst = star_allele ? nval1 + 2 : nval1;
            memset(buf->tmp2,0,nval1_dst*nsmpl);
            for (i=0; i<nsmpl; i++)
            {
                kstring_t dst;
                dst.l = 0; dst.m = nval1_dst; dst.s = (char*)buf->tmp2 + nval1_dst*i;
                kputc_('.',&dst);
                if ( star_allele ) kputsn_(",.",2,&dst);
                if ( len==BCF_VL_R )
                {
                    kputsn_(",.",2,&dst);
                    copy_string_field(buf->tmp+nval1*i, 0, nval1, &dst, 0);
                }
                copy_string_field(buf->tmp+nval1*i, iori+ioff, nval1, &dst, 0+ioff);
                if ( star_allele ) copy_string_field(".", 0, 1, &dst, 1+ioff);
            }
            ret = bcf_update_format(buf->out_hdr, out, tag, buf->tmp2, nval1_dst*nsmpl, type);
        }
        else if ( len==BCF_VL_R && type!=BCF_HT_STR )
        {
            int iori = buf->split.atoms[iout]->ial;
            assert( iori<=nval );
            for (i=0; i<nsmpl; i++)
            {
                void *src = buf->tmp  + nval1*num_size*i;
                void *dst = buf->tmp2 + num_size*i*(star_allele+2);
                memcpy(dst,src,num_size);
                memcpy(dst+num_size,src+iori*num_size,num_size);
                if ( type==BCF_HT_INT && mode==M_SUM )
                {
                    uint8_t *tbl = buf->split.tbl + iout*buf->split.nori;
                    for (j=iori; j<buf->split.nori; j++)
                        if ( tbl[j]==1 ) ((int32_t*)dst)[1] += ((int32_t*)src)[j+1];
                }
                if ( star_allele )
                    memcpy(dst+num_size*2,missing_ptr,num_size);
            }
            ret = bcf_update_format(buf->out_hdr, out, tag, buf->tmp2, nsmpl*(star_allele+2), type);
        }
        else if ( len==BCF_VL_G && type!=BCF_HT_STR )
        {
            int iori = buf->split.atoms[iout]->ial;
            int i01  = bcf_alleles2gt(0,iori);
            int i11  = bcf_alleles2gt(iori,iori);
            assert( iori<nval );
            #define BRANCH(type_t, is_missing, is_vector_end, set_missing, set_vector_end) { \
                for (i=0; i<nsmpl; i++) \
                { \
                    type_t *src = (type_t*)buf->tmp + i*nval1; \
                    type_t *dst = (type_t*)buf->tmp2 + i*3*(1+star_allele); \
                    int n=0; /* determine ploidy of this genotype */ \
                    while ( n<nval1 && !(is_vector_end) ) { n++; src++; } \
                    src = (type_t*)buf->tmp + i*nval1; \
                    memcpy(dst++,src,sizeof(type)); \
                    int nmiss = 0, nend = 0; \
                    if ( n==rec->n_allele ) /* haploid */ \
                    { \
                        memcpy(dst++,src+iori,sizeof(type)); \
                        if ( star_allele ) { nmiss = 1; nend = 3; } \
                        else nend = 1; \
                    } \
                    else if ( n==nval1 ) \
                    { \
                        memcpy(dst++,src+i01,sizeof(type)); \
                        memcpy(dst++,src+i11,sizeof(type)); \
                        if ( star_allele ) nmiss = 3; \
                    } \
                    else if ( n==1 && is_missing ) \
                    { \
                        if ( star_allele ) nend = 5; \
                        else nend = 2; \
                    } \
                    else  \
                        error("Incorrect number of values at %s:%"PRIhts_pos" .. tag=FORMAT/%s Number=G nAlleles=%d nValues=%d, %d-th sample\n", \
                                bcf_seqname(buf->hdr,rec),rec->pos+1,tag,rec->n_allele,n,i+1); \
                    for (j=0; j<nmiss; j++) { set_missing; dst++; } \
                    for (j=0; j<nend; j++) { set_vector_end; dst++; } \
                } \
            }
            switch (type)
            {
                case BCF_HT_INT:  BRANCH(int32_t, *src==bcf_int32_missing, *src==bcf_int32_vector_end, *dst=bcf_int32_missing, *dst=bcf_int32_vector_end); break;
                case BCF_HT_REAL: BRANCH(float, bcf_float_is_missing(*src), bcf_float_is_vector_end(*src), bcf_float_set_missing(*dst), bcf_float_set_vector_end(*dst)); break;
                default: error("Unexpected case: %d\n", type);
            }
            #undef BRANCH
            ret = bcf_update_format(buf->out_hdr, out, tag, buf->tmp2, 3*(1+star_allele)*nsmpl, type);
        }
        if ( ret!=0 ) error("An error occurred while updating FORMAT/%s\n",tag);
    }
}
static inline int _is_acgtn(char *seq)
{
    while ( *seq )
    {
        char c = toupper(*seq);
        if ( c!='A' && c!='C' && c!='G' && c!='T' && c!='N' ) return 0;
        seq++;
    }
    return 1;
}
/*
    The atomization works as follows:
    - Atomize each alternate allele separately by leaving out sequence identical to the reference. No
      alignment is performed, just greedy trimming of the end, then from left. This operation returns
      a list of atoms (atom_t) which carry fragments of REF,ALT and their positions as 0-based offsets
      to the original REF allele
    - Sort atoms by POS, REF and ALT. Each unique atom (POS+REF+ALT) forms a new VCF record, each
      with a single ALT.
    - For each new VCF record determine how to translate the original allele index (iori) to this new
      record:
        - 1: the original allele matches the atom
        - 0: the original allele does not overlap this atom or the overlapping part matches the REF
             allele
        - 2 (or equivalently "."): there is a mismatch between the original allele and the atom
      The mapping is encoded in a table with columns corresponding to the original ALTs and rows
      to the new POS+ALTs (atoms). The table is initialized to 0, then we set 1's for matching
      atoms and 2's for overlapping mismatching atoms.

    Note that different ALT alleles can result in the same atom (the same output line) and this code
    does not know how to reconcile possibly conflicting VCF annotations. This could be improved
    and merge logic provided, similarly to `merge -l`. For example, the allelic depths (AD) should
    be summed for the same atomized output allele. However, this level of complexity is not addressed
    in this initial draft. Higher priority for now is to provide the inverse "join" operation.

    Update 2021-04-09:
        Tags QS,AD are now automatically incremented as they should be, for both INFO and FORMAT.
        Note that the code will fail on missing values (todo) and it needs to be generalized and
        made customizable.
*/
void _abuf_split(abuf_t *buf, bcf1_t *rec)
{
    int i,j;
    if ( rec->n_allele < 2 )
    {
        rbuf_expand0(&buf->rbuf, bcf1_t*, buf->rbuf.n+1, buf->vcf);
        int j = rbuf_append(&buf->rbuf);
        if ( buf->vcf[j] ) bcf_destroy(buf->vcf[j]);
        buf->vcf[j] = bcf_dup(rec);
        return;
    }
    for (i=0; i<rec->n_allele; i++)
    {
        if ( _is_acgtn(rec->d.allele[i]) ) continue;
        rbuf_expand0(&buf->rbuf, bcf1_t*, buf->rbuf.n+1, buf->vcf);
        int j = rbuf_append(&buf->rbuf);
        if ( buf->vcf[j] ) bcf_destroy(buf->vcf[j]);
        buf->vcf[j] = bcf_dup(rec);
        return;
    }

    buf->natoms = 0;
    for (i=1; i<rec->n_allele; i++) _atomize_allele(buf,rec,i);
    qsort(buf->atoms,buf->natoms,sizeof(*buf->atoms),_cmp_atoms);
    _split_table_init(buf,rec,buf->natoms);
    for (i=0; i<buf->natoms; i++)
    {
        if ( i && !_atoms_inconsistent(&buf->atoms[i-1],&buf->atoms[i]) ) continue;
        _split_table_new(buf, &buf->atoms[i]);  // add a new unique output atom
    }
    for (i=0; i<buf->natoms; i++)
    {
        // Looping over sorted list of all atoms with possible duplicates from different source ALT alleles
        atom_t *atom = &buf->atoms[i];
        for (j=0; j<buf->split.nout; j++)
        {
            atom_t *out = buf->split.atoms[j];
            if ( atom == out ) continue;            // table already set to 1
            if ( atom->beg > out->end ) continue;   // cannot overlap this output atom
            if ( atom->end < out->beg ) break;      // this atom is ahead of all subsequent output records
            _split_table_overlap(buf, j, atom);
        }
    }
    // _split_table_print(buf);
    // _split_table_print_atoms(buf);
    assert( !buf->rbuf.n ); // all records should be flushed first in the SPLIT mode

    // Create the output records, transferring all annotations:
    // CHROM-QUAL
    _split_table_set_chrom_qual(buf);

    // INFO
    for (i=0; i<rec->n_info; i++)
    {
        // this implementation of merging rules is temporary: generalize and made customizable through the API
        merge_rule_t mode = M_FIRST;
        const char *tag = bcf_hdr_int2id(buf->hdr,BCF_DT_ID,rec->d.info[i].key);
        if ( !strcmp(tag,"QS") || !strcmp(tag,"AD") ) mode = M_SUM;

        _split_table_set_info(buf, &rec->d.info[i], mode);
    }

    // Set INFO tag showing the original record
    if ( buf->split.info_tag )
        _split_table_set_history(buf);

    // FORMAT
    for (i=0; i<rec->n_fmt; i++)
    {
        // this implementation of merging rules is temporary: generalize and made customizable through the API
        merge_rule_t mode = M_FIRST;
        const char *tag = bcf_hdr_int2id(buf->hdr,BCF_DT_ID,rec->d.fmt[i].id);
        if ( !strcmp(tag,"QS") || !strcmp(tag,"AD") ) mode = M_SUM;

        _split_table_set_format(buf, &rec->d.fmt[i], mode);
    }

    // Check that at least one FORMAT field was added, if not, the number of samples must be set manually
    for (i=0; i<buf->split.nout; i++)
    {
        bcf1_t *out = buf->vcf[rbuf_kth(&buf->rbuf,i)];
        if ( !out->n_sample ) out->n_sample = rec->n_sample;
    }
}

void abuf_push(abuf_t *buf, bcf1_t *rec)
{
    bcf_unpack(rec, BCF_UN_ALL);
    if ( buf->mode==SPLIT ) _abuf_split(buf,rec);
}

bcf1_t *abuf_flush(abuf_t *buf, int flush_all)
{
    int i;

    if ( buf->rbuf.n==0 ) return NULL;
    if ( flush_all ) goto ret;

ret:
    i = rbuf_shift(&buf->rbuf);
    return buf->vcf[i];
}

