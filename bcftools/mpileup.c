/*  mpileup.c -- mpileup subcommand. Previously bam_plcmd.c from samtools

    Copyright (C) 2008-2024 Genome Research Ltd.
    Portions copyright (C) 2009-2012 Broad Institute.

    Author: Heng Li <lh3@sanger.ac.uk>

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <limits.h>
#include <inttypes.h>
#include <errno.h>
#include <sys/stat.h>
#include <getopt.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/khash_str2int.h>
#include <htslib/hts_os.h>
#include <assert.h>
#include "regidx.h"
#include "bcftools.h"
#include "bam2bcf.h"
#include "bam_sample.h"
#include "gvcf.h"

#define MPLP_BCF        1
#define MPLP_VCF        (1<<1)
#define MPLP_NO_COMP    (1<<2)
#define MPLP_NO_ORPHAN  (1<<3)
#define MPLP_REALN      (1<<4)
#define MPLP_NO_INDEL   (1<<5)
#define MPLP_REDO_BAQ   (1<<6)
#define MPLP_ILLUMINA13 (1<<7)
#define MPLP_IGNORE_RG  (1<<8)
#define MPLP_PRINT_POS  (1<<9)
#define MPLP_PRINT_MAPQ (1<<10)
#define MPLP_PER_SAMPLE (1<<11)
#define MPLP_SMART_OVERLAPS (1<<12)
#define MPLP_REALN_PARTIAL  (1<<13)

typedef struct _mplp_aux_t mplp_aux_t;
typedef struct _mplp_pileup_t mplp_pileup_t;

// Data shared by all bam files
typedef struct {
    int min_mq, flag, min_baseQ, max_baseQ, delta_baseQ, capQ_thres, max_depth,
        max_indel_depth, max_read_len, ambig_reads;
    uint32_t fmt_flag;
    int rflag_skip_any_unset, rflag_skip_all_unset, rflag_skip_any_set, rflag_skip_all_set, output_type;
    int openQ, extQ, tandemQ, min_support, indel_win_size; // for indels
    int seqQ_offset;
    double min_frac; // for indels
    double indel_bias, poly_mqual;
    double del_bias; // compensate for diff deletion vs insertion error rates
    double vs_ref;
    char *reg_fname, *pl_list, *fai_fname, *output_fname;
    int reg_is_file, record_cmd_line, n_threads, clevel;
    faidx_t *fai;
    regidx_t *bed, *reg;    // bed: skipping regions, reg: index-jump to regions
    regitr_t *bed_itr, *reg_itr;
    int bed_logic;          // 1: include region, 0: exclude region
    gvcf_t *gvcf;

    // auxiliary structures for calling
    bcf_callaux_t *bca;
    bcf_callret1_t *bcr;
    bcf_call_t bc;
    bam_mplp_t iter;
    mplp_aux_t **mplp_data;
    int nfiles;
    char **files;
    mplp_pileup_t *gplp;
    int *n_plp;
    const bam_pileup1_t **plp;
    bam_smpl_t *bsmpl;
    kstring_t buf;
    bcf1_t *bcf_rec;
    htsFile *bcf_fp;
    bcf_hdr_t *bcf_hdr;
    int indels_v20;
    int edlib;
    int argc;
    char **argv;
    int write_index;
    char *index_fn;
} mplp_conf_t;

typedef struct {
    char *ref[2];
    int ref_id[2];
    int ref_len[2];
} mplp_ref_t;

#define MPLP_REF_INIT {{NULL,NULL},{-1,-1},{0,0}}

// Data specific to each bam file
struct _mplp_aux_t {
    samFile *fp;
    hts_itr_t *iter;
    bam_hdr_t *h;
    mplp_ref_t *ref;
    const mplp_conf_t *conf;
    int bam_id;
    hts_idx_t *idx;     // maintained only with more than one -r regions
};

// Data passed to htslib/mpileup
struct _mplp_pileup_t {
    int n;
    int *n_plp, *m_plp;
    bam_pileup1_t **plp;
};

static int mplp_get_ref(mplp_aux_t *ma, int tid,  char **ref, int *ref_len) {
    mplp_ref_t *r = ma->ref;

    //printf("get ref %d {%d/%p, %d/%p}\n", tid, r->ref_id[0], r->ref[0], r->ref_id[1], r->ref[1]);

    if (!r || !ma->conf->fai) {
        *ref = NULL;
        return 0;
    }

    // Do we need to reference count this so multiple mplp_aux_t can
    // track which references are in use?
    // For now we just cache the last two. Sufficient?
    if (tid == r->ref_id[0]) {
        *ref = r->ref[0];
        *ref_len = r->ref_len[0];
        return 1;
    }
    if (tid == r->ref_id[1]) {
        // Last, swap over
        int tmp;
        tmp = r->ref_id[0];  r->ref_id[0]  = r->ref_id[1];  r->ref_id[1]  = tmp;
        tmp = r->ref_len[0]; r->ref_len[0] = r->ref_len[1]; r->ref_len[1] = tmp;

        char *tc;
        tc = r->ref[0]; r->ref[0] = r->ref[1]; r->ref[1] = tc;
        *ref = r->ref[0];
        *ref_len = r->ref_len[0];
        return 1;
    }

    // New, so migrate to old and load new
    free(r->ref[1]);
    r->ref[1]     = r->ref[0];
    r->ref_id[1]  = r->ref_id[0];
    r->ref_len[1] = r->ref_len[0];

    r->ref_id[0] = tid;
    r->ref[0] = faidx_fetch_seq(ma->conf->fai,
                                ma->h->target_name[r->ref_id[0]],
                                0,
                                INT_MAX,
                                &r->ref_len[0]);

    if (!r->ref[0]) {
        r->ref[0] = NULL;
        r->ref_id[0] = -1;
        r->ref_len[0] = 0;
        *ref = NULL;
        return 0;
    }

    *ref = r->ref[0];
    *ref_len = r->ref_len[0];
    return 1;
}

static int mplp_func(void *data, bam1_t *b)
{
    char *ref;
    mplp_aux_t *ma = (mplp_aux_t*)data;
    int ret, ref_len;
    while (1)
    {
        int has_ref;
        ret = ma->iter? sam_itr_next(ma->fp, ma->iter, b) : sam_read1(ma->fp, ma->h, b);
        if (ret < 0) break;
        // The 'B' cigar operation is not part of the specification, considering as obsolete.
        //  bam_remove_B(b);
        if (b->core.tid < 0 || (b->core.flag&BAM_FUNMAP)) continue; // exclude unmapped reads
        if (ma->conf->rflag_skip_any_unset && (ma->conf->rflag_skip_any_unset&b->core.flag)!=ma->conf->rflag_skip_any_unset) continue;
        if (ma->conf->rflag_skip_all_set && (ma->conf->rflag_skip_all_set&b->core.flag)==ma->conf->rflag_skip_all_set) continue;
        if (ma->conf->rflag_skip_all_unset && !(ma->conf->rflag_skip_all_unset&b->core.flag)) continue;
        if (ma->conf->rflag_skip_any_set && ma->conf->rflag_skip_any_set&b->core.flag) continue;
        if (ma->conf->bed)
        {
            // test overlap
            regitr_t *itr = ma->conf->bed_itr;
            int beg = b->core.pos, end = bam_endpos(b)-1;
            int overlap = regidx_overlap(ma->conf->bed, ma->h->target_name[b->core.tid],beg,end, itr);
            if ( !ma->conf->bed_logic && !overlap )
            {
                // exclude only reads which are fully contained in the region
                while ( regitr_overlap(itr) )
                {
                    if ( beg < itr->beg ) { overlap = 1; break; }
                    if ( end > itr->end ) { overlap = 1; break; }
                }
            }
            if ( !overlap ) continue;
        }
        if ( bam_smpl_get_sample_id(ma->conf->bsmpl,ma->bam_id,b)<0 ) continue;
        if (ma->conf->flag & MPLP_ILLUMINA13) {
            int i;
            uint8_t *qual = bam_get_qual(b);
            for (i = 0; i < b->core.l_qseq; ++i)
                qual[i] = qual[i] > 31? qual[i] - 31 : 0;
        }

        if (ma->conf->fai && b->core.tid >= 0) {
            has_ref = mplp_get_ref(ma, b->core.tid, &ref, &ref_len);
            if (has_ref && ref_len <= b->core.pos) { // exclude reads outside of the reference sequence
                fprintf(stderr,"[%s] Skipping because %"PRId64" is outside of %d [ref:%d]\n",
                        __func__, (int64_t) b->core.pos, ref_len, b->core.tid);
                continue;
            }
        } else {
            has_ref = 0;
        }

        // Allow sufficient room for bam_aux_append of ZQ tag without
        // a realloc and consequent breakage of pileup's cached pointers.
        if (has_ref && (ma->conf->flag &MPLP_REALN) && !bam_aux_get(b, "ZQ")) {
            // Doing sam_prob_realn later is problematic as it adds to
            // the tag list (ZQ or BQ), which causes a realloc of b->data.
            // This happens after pileup has built a hash table on the
            // read name.  It's a deficiency in pileup IMO.

            // We could implement a new sam_prob_realn that returns ZQ
            // somewhere else and cache it ourselves (pileup clientdata),
            // but for now we simply use a workaround.
            //
            // We create a fake tag of the correct length, which we remove
            // just prior calling sam_prob_realn so we can guarantee there is
            // room. (We can't just make room now as bam_copy1 removes it
            // again).
            if (b->core.l_qseq > 500) {
                uint8_t *ZQ = malloc((uint32_t)b->core.l_qseq+1);
                memset(ZQ, '@', b->core.l_qseq);
                ZQ[b->core.l_qseq] = 0;
                bam_aux_append(b, "_Q", 'Z', b->core.l_qseq+1, ZQ);
                free(ZQ);
            } else {
                static uint8_t ZQ[501] =
                    "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
                    "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
                    "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
                    "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
                    "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
                    "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
                    "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
                    "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
                    "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
                    "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@";
                ZQ[b->core.l_qseq] = 0;
                bam_aux_append(b, "_Q", 'Z', b->core.l_qseq+1, ZQ);
                ZQ[b->core.l_qseq] = '@';
            }
        }

        if (has_ref && ma->conf->capQ_thres > 10) {
            int q = sam_cap_mapq(b, ref, ref_len, ma->conf->capQ_thres);
            if (q < 0) continue;    // skip
            else if (b->core.qual > q) b->core.qual = q;
        }
        if (b->core.qual < ma->conf->min_mq) continue;
        else if ((ma->conf->flag&MPLP_NO_ORPHAN) && (b->core.flag&BAM_FPAIRED) && !(b->core.flag&BAM_FPROPER_PAIR)) continue;

        return ret;
    };
    return ret;
}

// Called once per new bam added to the pileup.
// We cache sample information here so we don't have to keep recomputing this
// on each and every pileup column. If FMT/SCR annotation is requested, a flag
// is set to indicate the presence of a soft clip.
static int pileup_constructor(void *data, const bam1_t *b, bam_pileup_cd *cd)
{
    cd->p = calloc(1,sizeof(plp_cd_t));

    PLP_NM(cd) = PLP_NM_UNSET;

    mplp_aux_t *ma = (mplp_aux_t *)data;
    int n = bam_smpl_get_sample_id(ma->conf->bsmpl, ma->bam_id, (bam1_t *)b);
    PLP_SET_SAMPLE_ID(cd, n);

    // Whether read has a soft-clip is used in mplp_realn's heuristics.
    // TODO: consider whether clip length is beneficial to use?
    int i;
    for (i=0; i<b->core.n_cigar; i++) {
        int cig = bam_get_cigar(b)[i] & BAM_CIGAR_MASK;
        if (cig == BAM_CSOFT_CLIP) {
            PLP_SET_SOFT_CLIP(cd);
            break;
        }
    }

    if (ma->conf->flag & MPLP_REALN) {
        int i;
        // int tot_ins = 0;
        // int p = 0;
        uint32_t *cigar = bam_get_cigar(b);
        for (i=0; i<b->core.n_cigar; i++) {
            int cig = cigar[i] & BAM_CIGAR_MASK;
            // if (bam_cigar_type(cig) & 2)
            //     p += cigar[i] >> BAM_CIGAR_SHIFT;
            if (cig == BAM_CINS || cig == BAM_CDEL || cig == BAM_CREF_SKIP) {
                // tot_ins += cigar[i] >> BAM_CIGAR_SHIFT;
                // Possible further optimsation, check tot_ins==1 later
                // (and remove break) so we can detect single bp indels.
                // We may want to focus BAQ on more complex regions only.
                PLP_SET_INDEL(cd);
                break;
            }

            // TODO: proper p->cd struct and have cd->i as a size rather
            // than a flag.

            // Then aggregate together the sizes and if just 1 size for all
            // reads or 2 sizes for approx 50/50 split in all reads, then
            // treat this as a well-aligned variant and don't run BAQ.
        }
    }

    return 0;
}
static int pileup_destructor(void *data, const bam1_t *b, bam_pileup_cd *cd)
{
    free(cd->p);
    return 0;
}

static void group_smpl(mplp_pileup_t *m, bam_smpl_t *bsmpl, int n, int *n_plp, const bam_pileup1_t **plp)
{
    int i, j;
    memset(m->n_plp, 0, m->n * sizeof(int));
    for (i = 0; i < n; ++i) // iterate over all bams
    {
        for (j = 0; j < n_plp[i]; ++j)  // iterate over all reads available at this position
        {
            const bam_pileup1_t *p = plp[i] + j;
            int id = PLP_SAMPLE_ID(&(p->cd));
            if (m->n_plp[id] == m->m_plp[id])
            {
                m->m_plp[id] = m->m_plp[id]? m->m_plp[id]<<1 : 8;
                m->plp[id] = (bam_pileup1_t*) realloc(m->plp[id], sizeof(bam_pileup1_t) * m->m_plp[id]);
            }
            m->plp[id][m->n_plp[id]++] = *p;
        }
    }
}

static void flush_bcf_records(mplp_conf_t *conf, htsFile *fp, bcf_hdr_t *hdr, bcf1_t *rec)
{
    if ( !conf->gvcf )
    {
        if ( rec && bcf_write1(fp, hdr, rec)!=0 ) error("[%s] Error: failed to write the record to %s\n", __func__,conf->output_fname?conf->output_fname:"standard output");
        return;
    }

    if ( !rec )
    {
        gvcf_write(conf->gvcf, fp, hdr, NULL, 0);
        return;
    }

    int is_ref = 0;
    if ( rec->n_allele==1 ) is_ref = 1;
    else if ( rec->n_allele==2 )
    {
        // second allele is mpileup's X, not a variant
        if ( rec->d.allele[1][0]=='<' && rec->d.allele[1][1]=='*' && rec->d.allele[1][2]=='>' ) is_ref = 1;
    }
    rec = gvcf_write(conf->gvcf, fp, hdr, rec, is_ref);
    if ( rec && bcf_write1(fp,hdr,rec)!=0 ) error("[%s] Error: failed to write the record to %s\n", __func__,conf->output_fname?conf->output_fname:"standard output");
}

/*
 * Loops for an indel at this position.
 *
 * Only reads that overlap an indel loci get realigned.  This considerably
 * reduces the cost of running BAQ while keeping the main benefits.
 *
 * TODO: also consider only realigning reads that don't span the indel
 * by more than a certain amount either-side.  Ie focus BAQ only on reads
 * ending adjacent to the indel, where the alignment is most likely to
 * be wrong.  (2nd TODO: do this based on sequence context; STRs bad, unique
 * data good.)
 *
 * NB: this may sadly realign after we've already used the data.  Hmm...
 */
static void mplp_realn(int n, int *n_plp, const bam_pileup1_t **plp,
                       int flag, int max_read_len,
                       char *ref, int ref_len, int pos) {
    int i, j, has_indel = 0, has_clip = 0, nt = 0;
    int min_indel = INT_MAX, max_indel = INT_MIN;

    // Is an indel present.
    // NB: don't bother even checking if very long as almost guaranteed
    // to have indel (and likely soft-clips too).
    for (i = 0; i < n; i++) { // iterate over bams
        nt += n_plp[i];
        for (j = 0; j < n_plp[i]; j++) { // iterate over reads
            bam_pileup1_t *p = (bam_pileup1_t *)plp[i] + j;
            has_indel += (PLP_HAS_INDEL(&p->cd) || p->indel) ? 1 : 0;
            // Has_clip is almost always true for very long reads
            // (eg PacBio CCS), but these rarely matter as the clip
            // is likely a long way from this indel.
            has_clip  += (PLP_HAS_SOFT_CLIP(&p->cd))         ? 1 : 0;
            if (max_indel < p->indel)
                max_indel = p->indel;
            if (min_indel > p->indel)
                min_indel = p->indel;
        }
    }

    if (flag & MPLP_REALN_PARTIAL) {
        if (has_indel == 0 ||
            (has_clip < 0.2*nt && max_indel == min_indel &&
             (has_indel < 0.1*nt /*|| has_indel > 0.9*nt*/ || has_indel == 1)))
            return;
    }

    // Realign
    for (i = 0; i < n; i++) { // iterate over bams
        for (j = 0; j < n_plp[i]; j++) { // iterate over reads
            const bam_pileup1_t *p = plp[i] + j;
            bam1_t *b = p->b;

            // Avoid doing multiple times.
            //
            // Note we cannot modify p->cd.i here with a PLP_SET macro
            // because the cd item is held by mpileup in an lbnode_t
            // struct and copied over to the pileup struct for each
            // iteration, essentially making p->cd.i read only.
            //
            // We could use our own structure (p->cd.p), allocated during
            // the constructor, but for simplicity we play dirty and
            // abuse an unused flag bit instead.
            if ( PLP_IS_REALN(&(p->cd)) ) continue;
            PLP_SET_REALN(&(p->cd));

            if (b->core.l_qseq > max_read_len)
                continue;

            // Check p->cigar_ind and see what cigar elements are before
            // and after.  How close is this location to the end of the
            // read?  Only realign if we don't span by more than X bases.
            //
            // Again, best only done on deeper data as BAQ helps
            // disproportionately more on shallow data sets.
            //
            // This rescues some of the false negatives that are caused by
            // systematic reduction in quality due to sample vs ref alignment.

// At deep coverage we skip realigning more reads as we have sufficient depth.
// This rescues for false negatives.  At shallow depth we pay for this with
// more FP so are more stringent on spanning size.
#define REALN_DIST (40+10*(nt<40)+10*(nt<20))
            uint32_t *cig = bam_get_cigar(b);
            int ncig = b->core.n_cigar;

            // Don't realign reads where indel is in middle?
            // On long read data we don't care about soft-clips at the ends.
            // For short read data, we always calc BAQ on these as they're
            // a common source of false positives.
            if ((flag & MPLP_REALN_PARTIAL) && nt > 15 && ncig > 1) {
                // Left & right cigar op match.
                int lr = b->core.l_qseq > 500;
                int lm = 0, rm = 0, k, nm = 0;
                for (k = 0; k < ncig; k++) {
                    int cop = bam_cigar_op(cig[k]);
                    if (lr && (cop == BAM_CHARD_CLIP || cop == BAM_CSOFT_CLIP))
                        continue;

                    if (cop == BAM_CMATCH || cop == BAM_CDIFF ||
                        cop == BAM_CEQUAL) {
                        lm += bam_cigar_oplen(cig[k]);
                        nm++;
                    } else {
                        break;
                    }
                }

                // if everything is a match (or sequence (mis)match) then move on
                // because we don't have an indel in the middle
                if (nm != ncig) {
                    for (k = ncig-1; k >= 0; k--) {
                        int cop = bam_cigar_op(cig[k]);
                        if (lr && (cop == BAM_CHARD_CLIP || cop == BAM_CSOFT_CLIP))
                            continue;

                        if (cop == BAM_CMATCH || cop == BAM_CDIFF ||
                            cop == BAM_CEQUAL)
                            rm += bam_cigar_oplen(cig[k]);
                        else
                            break;
                    }

                    if (lm >= REALN_DIST*4 && rm >= REALN_DIST*4)
                        continue;

                    if (lm >= REALN_DIST && rm >= REALN_DIST &&
                        has_clip < (0.15+0.05*(nt>20))*nt)
                        continue;
                }
            }

            if (b->core.l_qseq > 500) {
                // don't do BAQ on long-read data if it's going to
                // cause us to have a large band-with and costly in CPU
                int rl = bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
                if (abs(rl - b->core.l_qseq) * b->core.l_qseq >= 500000)
                    continue;
            }

            // Fudge: make room for ZQ tag.
            uint8_t *_Q = bam_aux_get(b, "_Q");
            if (_Q) bam_aux_del(b, _Q);
            sam_prob_realn(b, ref, ref_len, (flag & MPLP_REDO_BAQ) ? 7 : 3);
        }
    }

    return;
}

static int mpileup_reg(mplp_conf_t *conf, uint32_t beg, uint32_t end)
{
    bam_hdr_t *hdr = conf->mplp_data[0]->h; // header of first file in input list

    int ret, i, tid, pos, ref_len;
    char *ref;

    while ( (ret=bam_mplp_auto(conf->iter, &tid, &pos, conf->n_plp, conf->plp)) > 0)
    {
        if ( pos<beg || pos>end ) continue;
        if ( conf->bed && tid >= 0 )
        {
            int overlap = regidx_overlap(conf->bed, hdr->target_name[tid], pos, pos, NULL);
            if ( !conf->bed_logic ) overlap = overlap ? 0 : 1;
            if ( !overlap ) continue;
        }
        int has_ref = mplp_get_ref(conf->mplp_data[0], tid, &ref, &ref_len);
        if (has_ref && (conf->flag & MPLP_REALN))
            mplp_realn(conf->nfiles, conf->n_plp, conf->plp, conf->flag, conf->max_read_len, ref, ref_len, pos);

        int total_depth, _ref0, ref16;
        for (i = total_depth = 0; i < conf->nfiles; ++i) total_depth += conf->n_plp[i];
        group_smpl(conf->gplp, conf->bsmpl, conf->nfiles, conf->n_plp, conf->plp);
        _ref0 = (ref && pos < ref_len)? ref[pos] : 'N';
        ref16 = seq_nt16_table[_ref0];
        bcf_callaux_clean(conf->bca, &conf->bc);
        for (i = 0; i < conf->gplp->n; ++i)
            bcf_call_glfgen(conf->gplp->n_plp[i], conf->gplp->plp[i], ref16, conf->bca, conf->bcr + i);
        conf->bc.tid = tid; conf->bc.pos = pos;
        bcf_call_combine(conf->gplp->n, conf->bcr, conf->bca, ref16, &conf->bc);
        bcf_clear1(conf->bcf_rec);
        bcf_call2bcf(&conf->bc, conf->bcf_rec, conf->bcr, conf->fmt_flag, conf->bca, 0);
        flush_bcf_records(conf, conf->bcf_fp, conf->bcf_hdr, conf->bcf_rec);

        // call indels; todo: subsampling with total_depth>max_indel_depth instead of ignoring?
        // check me: rghash in bcf_call_gap_prep() should have no effect, reads mplp_func already excludes them
        if (!(conf->flag&MPLP_NO_INDEL) && total_depth < conf->max_indel_depth)
        {
            bcf_callaux_clean(conf->bca, &conf->bc);
            conf->bca->chr = tid>=0 ? hdr->target_name[tid] : NULL;
            int iret;
            if (conf->edlib)
                iret = bcf_edlib_gap_prep(conf->gplp->n, conf->gplp->n_plp, conf->gplp->plp, pos, conf->bca, ref, ref_len);
            else if ( conf->indels_v20 )
                iret = bcf_iaux_gap_prep(conf->gplp->n, conf->gplp->n_plp, conf->gplp->plp, pos, conf->bca, ref);
            else
                iret = bcf_call_gap_prep(conf->gplp->n, conf->gplp->n_plp, conf->gplp->plp, pos, conf->bca, ref);
            if ( iret>=0 )
            {
                for (i = 0; i < conf->gplp->n; ++i)
                    bcf_call_glfgen(conf->gplp->n_plp[i], conf->gplp->plp[i], -1, conf->bca, conf->bcr + i);
                if (bcf_call_combine(conf->gplp->n, conf->bcr, conf->bca, -1, &conf->bc) >= 0)
                {
                    bcf_clear1(conf->bcf_rec);
                    bcf_call2bcf(&conf->bc, conf->bcf_rec, conf->bcr, conf->fmt_flag, conf->bca, ref);
                    flush_bcf_records(conf, conf->bcf_fp, conf->bcf_hdr, conf->bcf_rec);
                }
            }
        }
    }
    return ret;
}

static int mpileup(mplp_conf_t *conf)
{
    if (conf->nfiles == 0) {
        fprintf(stderr,"[%s] no input file/data given\n", __func__);
        exit(EXIT_FAILURE);
    }

    mplp_ref_t mp_ref = MPLP_REF_INIT;
    conf->gplp = (mplp_pileup_t *) calloc(1,sizeof(mplp_pileup_t));
    conf->mplp_data = (mplp_aux_t**) calloc(conf->nfiles, sizeof(mplp_aux_t*));
    conf->plp = (const bam_pileup1_t**) calloc(conf->nfiles, sizeof(bam_pileup1_t*));
    conf->n_plp = (int*) calloc(conf->nfiles, sizeof(int));

    // Allow to run mpileup on multiple regions in one go. This comes at cost: the bai index
    // must be kept in the memory for the whole time which can be a problem with many bams.
    // Therefore if none or only one region is requested, we initialize the bam iterator as
    // before and free the index. Only when multiple regions are queried, we keep the index.
    int nregs = 0;
    if ( conf->reg_fname )
    {
        if ( conf->reg_is_file )
        {
            conf->reg = regidx_init(conf->reg_fname,NULL,NULL,0,NULL);
            if ( !conf->reg ) {
                fprintf(stderr,"Could not parse the regions: %s\n", conf->reg_fname);
                exit(EXIT_FAILURE);
            }
        }
        else
        {
            conf->reg = regidx_init(NULL,regidx_parse_reg,NULL,sizeof(char*),NULL);
            if ( regidx_insert_list(conf->reg,conf->reg_fname,',') !=0 ) {
                fprintf(stderr,"Could not parse the regions: %s\n", conf->reg_fname);
                exit(EXIT_FAILURE);
            }
        }
        nregs = regidx_nregs(conf->reg);
        if ( nregs )
        {
            // the regions list can be empty, see #2250
            conf->reg_itr = regitr_init(conf->reg);
            regitr_loop(conf->reg_itr);   // region iterator now positioned at the first region
        }
    }

    // read the header of each file in the list and initialize data
    // beware: mpileup has always assumed that tid's are consistent in the headers, add sanity check at least!
    bam_hdr_t *hdr = NULL;      // header of first file in input list
    int i;
    for (i = 0; i < conf->nfiles; ++i) {
        bam_hdr_t *h_tmp;
        conf->mplp_data[i] = (mplp_aux_t*) calloc(1, sizeof(mplp_aux_t));
        conf->mplp_data[i]->fp = sam_open(conf->files[i], "rb");
        if ( !conf->mplp_data[i]->fp )
        {
            fprintf(stderr, "[%s] failed to open %s: %s\n", __func__, conf->files[i], strerror(errno));
            exit(EXIT_FAILURE);
        }
        if (hts_set_opt(conf->mplp_data[i]->fp, CRAM_OPT_DECODE_MD, 1)) {
            fprintf(stderr, "Failed to set CRAM_OPT_DECODE_MD value\n");
            exit(EXIT_FAILURE);
        }
        if (conf->fai_fname && hts_set_fai_filename(conf->mplp_data[i]->fp, conf->fai_fname) != 0) {
            fprintf(stderr, "[%s] failed to process %s: %s\n",
                    __func__, conf->fai_fname, strerror(errno));
            exit(EXIT_FAILURE);
        }
        conf->mplp_data[i]->conf = conf;
        conf->mplp_data[i]->ref = &mp_ref;
        h_tmp = sam_hdr_read(conf->mplp_data[i]->fp);
        if ( !h_tmp ) {
            fprintf(stderr,"[%s] fail to read the header of %s\n", __func__, conf->files[i]);
            exit(EXIT_FAILURE);
        }
        conf->mplp_data[i]->h = i ? hdr : h_tmp; // for j==0, "h" has not been set yet
        conf->mplp_data[i]->bam_id = bam_smpl_add_bam(conf->bsmpl,h_tmp->text,conf->files[i]);
        if ( conf->mplp_data[i]->bam_id<0 )
        {
            // no usable readgroups in this bam, it can be skipped
            sam_close(conf->mplp_data[i]->fp);
            free(conf->mplp_data[i]);
            bam_hdr_destroy(h_tmp);
            free(conf->files[i]);
            if ( i+1<conf->nfiles ) memmove(&conf->files[i],&conf->files[i+1],sizeof(*conf->files)*(conf->nfiles-i-1));
            conf->nfiles--;
            i--;
            continue;
        }
        if (conf->reg && nregs) {
            hts_idx_t *idx = sam_index_load(conf->mplp_data[i]->fp, conf->files[i]);
            if (idx == NULL) {
                fprintf(stderr, "[%s] fail to load index for %s\n", __func__, conf->files[i]);
                exit(EXIT_FAILURE);
            }
            conf->buf.l = 0;
            ksprintf(&conf->buf,"%s:%u-%u",conf->reg_itr->seq,conf->reg_itr->beg+1,conf->reg_itr->end+1);
            conf->mplp_data[i]->iter = sam_itr_querys(idx, conf->mplp_data[i]->h, conf->buf.s);
            if ( !conf->mplp_data[i]->iter )
            {
                conf->mplp_data[i]->iter = sam_itr_querys(idx, conf->mplp_data[i]->h, conf->reg_itr->seq);
                if ( conf->mplp_data[i]->iter ) {
                    fprintf(stderr,"[E::%s] fail to parse region '%s'\n", __func__, conf->buf.s);
                    exit(EXIT_FAILURE);
                }
                fprintf(stderr,"[E::%s] the sequence \"%s\" not found: %s\n",__func__,conf->reg_itr->seq,conf->files[i]);
                exit(EXIT_FAILURE);
            }
            if ( nregs==1 ) // no need to keep the index in memory
               hts_idx_destroy(idx);
            else
                conf->mplp_data[i]->idx = idx;
        }

        if ( !hdr ) hdr = h_tmp; /* save the header of first file in list */
        else {
            // FIXME: check consistency between h and h_tmp
            bam_hdr_destroy(h_tmp);

            // we store only the first file's header; it's (alleged to be)
            // compatible with the i-th file's target_name lookup needs
            conf->mplp_data[i]->h = hdr;
        }
    }
    if ( !hdr ) {
        fprintf(stderr, "[%s] failed to find a file header with usable read groups\n", __func__);
        exit(EXIT_FAILURE);
    }
    // allocate data storage proportionate to number of samples being studied sm->n
    bam_smpl_get_samples(conf->bsmpl, &conf->gplp->n);
    conf->gplp->n_plp = (int*) calloc(conf->gplp->n, sizeof(int));
    conf->gplp->m_plp = (int*) calloc(conf->gplp->n, sizeof(int));
    conf->gplp->plp = (bam_pileup1_t**) calloc(conf->gplp->n, sizeof(bam_pileup1_t*));

    fprintf(stderr, "[%s] %d samples in %d input files\n", __func__, conf->gplp->n, conf->nfiles);
    // write the VCF header
    char wmode[8];
    set_wmode(wmode,conf->output_type,conf->output_fname,conf->clevel);
    conf->bcf_fp = hts_open(conf->output_fname ? conf->output_fname : "-", wmode);
    if (conf->bcf_fp == NULL) {
        fprintf(stderr, "[%s] failed to write to %s: %s\n", __func__, conf->output_fname? conf->output_fname : "standard output", strerror(errno));
        exit(EXIT_FAILURE);
    }
    if ( conf->n_threads ) hts_set_threads(conf->bcf_fp, conf->n_threads);

    // BCF header creation
    conf->bcf_hdr = bcf_hdr_init("w");
    conf->buf.l = 0;

    if (conf->record_cmd_line)
    {
        ksprintf(&conf->buf, "##bcftoolsVersion=%s+htslib-%s\n",bcftools_version(),hts_version());
        bcf_hdr_append(conf->bcf_hdr, conf->buf.s);

        conf->buf.l = 0;
        ksprintf(&conf->buf, "##bcftoolsCommand=mpileup");
        for (i=1; i<conf->argc; i++) ksprintf(&conf->buf, " %s", conf->argv[i]);
        kputc('\n', &conf->buf);
        bcf_hdr_append(conf->bcf_hdr, conf->buf.s);
    }

    if (conf->fai_fname)
    {
        conf->buf.l = 0;
        ksprintf(&conf->buf, "##reference=file://%s\n", conf->fai_fname);
        bcf_hdr_append(conf->bcf_hdr, conf->buf.s);
    }

    // Translate BAM @SQ tags to BCF ##contig tags
    // todo: use/write new BAM header manipulation routines, fill also UR, M5
    for (i=0; i<hdr->n_targets; i++)
    {
        conf->buf.l = 0;
        ksprintf(&conf->buf, "##contig=<ID=%s,length=%d>", hdr->target_name[i], hdr->target_len[i]);
        bcf_hdr_append(conf->bcf_hdr, conf->buf.s);
    }
    conf->buf.l = 0;

    bcf_hdr_append(conf->bcf_hdr,"##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">");
    if ( conf->fmt_flag&B2B_INFO_IDV )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=IDV,Number=1,Type=Integer,Description=\"Maximum number of raw reads supporting an indel\">");
    if ( conf->fmt_flag&B2B_INFO_IMF )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=IMF,Number=1,Type=Float,Description=\"Maximum fraction of raw reads supporting an indel\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">");
    if ( conf->fmt_flag&B2B_INFO_VDB )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=VDB,Number=1,Type=Float,Description=\"Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)\",Version=\"3\">");

    if ( conf->fmt_flag&B2B_INFO_RPBZ )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=RPBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Read Position Bias (closer to 0 is better)\">");
    if ( conf->fmt_flag&B2B_INFO_MQBZ )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=MQBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Mapping Quality Bias (closer to 0 is better)\">");
    if ( conf->fmt_flag&B2B_INFO_BQBZ )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=BQBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Base Quality Bias (closer to 0 is better)\">");
    if ( conf->fmt_flag&B2B_INFO_MQSBZ )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=MQSBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Mapping Quality vs Strand Bias (closer to 0 is better)\">");
    if ( conf->fmt_flag&B2B_INFO_MIN_PL_SUM )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=MIN_PL_SUM,Number=1,Type=Integer,Description=\"Sum of min PLs across all samples before normalization (experimental)\">");
    if ( conf->fmt_flag&B2B_INFO_NM )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=NM,Number=2,Type=Float,Description=\"Average number of mismatches in ref and alt reads (approximate, experimental, make me localized?)\">");
    if ( conf->fmt_flag&B2B_INFO_NMBZ )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=NMBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Number of Mismatches within supporting reads (closer to 0 is better; approximate, experimental, make me localized?)\">");
    if ( conf->fmt_flag&B2B_FMT_NMBZ )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=NMBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Number of Mismatches within supporting reads (closer to 0 is better)\">");
    if ( conf->fmt_flag&B2B_INFO_SCBZ )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=SCBZ,Number=1,Type=Float,Description=\"Mann-Whitney U-z test of Soft-Clip Length Bias (closer to 0 is better)\">");
    if ( conf->fmt_flag&B2B_INFO_FS )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=FS,Number=1,Type=Float,Description=\"Fisher's exact test P-value to detect strand bias\">");
    if ( conf->fmt_flag&B2B_INFO_SGB )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=SGB,Number=1,Type=Float,Description=\"Segregation based metric, http://samtools.github.io/bcftools/rd-SegBias.pdf\">");
    if ( conf->fmt_flag&B2B_INFO_MQ0F )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=MQ0F,Number=1,Type=Float,Description=\"Fraction of MQ0 reads (smaller is better)\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=I16,Number=16,Type=Float,Description=\"Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h\">");
    bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=QS,Number=R,Type=Float,Description=\"Auxiliary tag used for calling\">");
    bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">");
    if ( conf->fmt_flag&B2B_FMT_DP )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of high-quality bases\">");
    if ( conf->fmt_flag&B2B_FMT_DV )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"Number of high-quality non-reference bases\">");
    if ( conf->fmt_flag&B2B_FMT_DPR )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=DPR,Number=R,Type=Integer,Description=\"Number of high-quality bases observed for each allele\">");
    if ( conf->fmt_flag&B2B_INFO_DPR )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=DPR,Number=R,Type=Integer,Description=\"Number of high-quality bases observed for each allele\">");
    if ( conf->fmt_flag&B2B_FMT_DP4 )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=DP4,Number=4,Type=Integer,Description=\"Number of high-quality ref-fwd, ref-reverse, alt-fwd and alt-reverse bases\">");
    if ( conf->fmt_flag&B2B_FMT_SP )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=SP,Number=1,Type=Integer,Description=\"Phred-scaled strand bias P-value\">");
    if ( conf->fmt_flag&B2B_FMT_AD )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths (high-quality bases)\">");
    if ( conf->fmt_flag&B2B_FMT_ADF )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=ADF,Number=R,Type=Integer,Description=\"Allelic depths on the forward strand (high-quality bases)\">");
    if ( conf->fmt_flag&B2B_FMT_ADR )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=ADR,Number=R,Type=Integer,Description=\"Allelic depths on the reverse strand (high-quality bases)\">");
    if ( conf->fmt_flag&B2B_FMT_QS )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=QS,Number=R,Type=Integer,Description=\"Phred-score allele quality sum used by `call -mG` and `+trio-dnm`\">");
    if ( conf->fmt_flag&B2B_INFO_AD )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total allelic depths (high-quality bases)\">");
    if ( conf->fmt_flag&B2B_INFO_ADF )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=ADF,Number=R,Type=Integer,Description=\"Total allelic depths on the forward strand (high-quality bases)\">");
    if ( conf->fmt_flag&B2B_INFO_SCR )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=SCR,Number=1,Type=Integer,Description=\"Number of soft-clipped reads (at high-quality bases)\">");
    if ( conf->fmt_flag&B2B_FMT_SCR )
        bcf_hdr_append(conf->bcf_hdr,"##FORMAT=<ID=SCR,Number=1,Type=Integer,Description=\"Per-sample number of soft-clipped reads (at high-quality bases)\">");
    if ( conf->fmt_flag&B2B_INFO_ADR )
        bcf_hdr_append(conf->bcf_hdr,"##INFO=<ID=ADR,Number=R,Type=Integer,Description=\"Total allelic depths on the reverse strand (high-quality bases)\">");
    if ( conf->gvcf )
        gvcf_update_header(conf->gvcf, conf->bcf_hdr);

    int nsmpl;
    const char **smpl = bam_smpl_get_samples(conf->bsmpl, &nsmpl);
    for (i=0; i<nsmpl; i++)
        bcf_hdr_add_sample(conf->bcf_hdr, smpl[i]);
    if ( bcf_hdr_write(conf->bcf_fp, conf->bcf_hdr)!=0 ) error("[%s] Error: failed to write the header to %s\n",__func__,conf->output_fname?conf->output_fname:"standard output");
    if ( init_index2(conf->bcf_fp,conf->bcf_hdr,conf->output_fname,
                     &conf->index_fn, conf->write_index) < 0 )
        error("Error: failed to initialise index for %s\n",conf->output_fname);

    conf->bca = bcf_call_init(-1., conf->min_baseQ, conf->max_baseQ,
                              conf->delta_baseQ);
    conf->bcr = (bcf_callret1_t*) calloc(nsmpl, sizeof(bcf_callret1_t));
    conf->bca->openQ = conf->openQ, conf->bca->extQ = conf->extQ, conf->bca->tandemQ = conf->tandemQ;
    conf->bca->indel_bias = conf->indel_bias;
    conf->bca->del_bias = conf->del_bias;
    conf->bca->min_frac = conf->min_frac;
    conf->bca->min_support = conf->min_support;
    conf->bca->per_sample_flt = conf->flag & MPLP_PER_SAMPLE;
    conf->bca->fmt_flag = conf->fmt_flag;
    conf->bca->ambig_reads = conf->ambig_reads;
    conf->bca->indel_win_size = conf->indel_win_size;
    conf->bca->indels_v20 = conf->indels_v20;
    conf->bca->edlib = conf->edlib;
    conf->bca->seqQ_offset = conf->seqQ_offset;
    conf->bca->poly_mqual = conf->poly_mqual;
    conf->bca->vs_ref = conf->vs_ref;

    conf->bc.bcf_hdr = conf->bcf_hdr;
    conf->bc.n  = nsmpl;
    conf->bc.PL = (int32_t*) malloc(15 * nsmpl * sizeof(*conf->bc.PL));
    conf->bc.QS = (int32_t*) malloc(nsmpl*sizeof(*conf->bc.QS)*B2B_MAX_ALLELES);
    for (i=0; i<nsmpl; i++)
        conf->bcr[i].QS = conf->bc.QS + i*B2B_MAX_ALLELES;
    if (conf->fmt_flag)
    {
        assert( sizeof(float)==sizeof(int32_t) );
        conf->bc.DP4 = (int32_t*) malloc(nsmpl * sizeof(int32_t) * 4);
        conf->bc.fmt_arr = (uint8_t*) malloc(nsmpl * sizeof(float)); // all fmt_flag fields, float and int32
        if ( conf->fmt_flag&(B2B_INFO_DPR|B2B_FMT_DPR|B2B_INFO_AD|B2B_INFO_ADF|B2B_INFO_ADR|B2B_FMT_AD|B2B_FMT_ADF|B2B_FMT_ADR) )
        {
            // first B2B_MAX_ALLELES fields for total numbers, the rest per-sample
            conf->bc.ADR = (int32_t*) malloc((nsmpl+1)*B2B_MAX_ALLELES*sizeof(int32_t));
            conf->bc.ADF = (int32_t*) malloc((nsmpl+1)*B2B_MAX_ALLELES*sizeof(int32_t));
            for (i=0; i<nsmpl; i++)
            {
                conf->bcr[i].ADR = conf->bc.ADR + (i+1)*B2B_MAX_ALLELES;
                conf->bcr[i].ADF = conf->bc.ADF + (i+1)*B2B_MAX_ALLELES;
            }
        }
        if ( conf->fmt_flag&(B2B_INFO_SCR|B2B_FMT_SCR) )
            conf->bc.SCR = (int32_t*) malloc((nsmpl+1)*sizeof(*conf->bc.SCR));
    }
    int nnmbz = (conf->fmt_flag&B2B_FMT_NMBZ) ? nsmpl + 1 : 1;
    conf->bc.ref_nm = (int32_t*) malloc(sizeof(*conf->bc.ref_nm) * nnmbz * B2B_N_NM);
    conf->bc.alt_nm = (int32_t*) malloc(sizeof(*conf->bc.alt_nm) * nnmbz * B2B_N_NM);
    conf->bc.mwu_nm = (float*) malloc((nsmpl+1)*sizeof(*conf->bc.mwu_nm));
    conf->bca->ref_nm = conf->bc.ref_nm;     // this is just to make the arrays available in bcf_call_glfgen()
    conf->bca->alt_nm = conf->bc.alt_nm;
    if ( conf->fmt_flag&B2B_FMT_NMBZ )
    {
        for (i=0; i<nsmpl; i++) conf->bcr[i].ref_nm = conf->bc.ref_nm + (i+1)*B2B_N_NM;
        for (i=0; i<nsmpl; i++) conf->bcr[i].alt_nm = conf->bc.alt_nm + (i+1)*B2B_N_NM;
    }

    // init mpileup
    conf->iter = bam_mplp_init(conf->nfiles, mplp_func, (void**)conf->mplp_data);
    if ( conf->flag & MPLP_SMART_OVERLAPS ) bam_mplp_init_overlaps(conf->iter);
    fprintf(stderr, "[%s] maximum number of reads per input file set to -d %d\n",  __func__, conf->max_depth);
    if ( (double)conf->max_depth * conf->nfiles > 1<<20)
        fprintf(stderr, "Warning: Potential memory hog, up to %.0fM reads in the pileup!\n", (double)conf->max_depth*conf->nfiles);
    if ( (double)conf->max_depth * conf->nfiles / nsmpl < 250 )
        fprintf(stderr, "Note: The maximum per-sample depth with -d %d is %.1fx\n", conf->max_depth,(double)conf->max_depth * conf->nfiles / nsmpl);
    bam_mplp_set_maxcnt(conf->iter, conf->max_depth);
    conf->max_indel_depth = conf->max_indel_depth * nsmpl;
    conf->bcf_rec = bcf_init1();
    bam_mplp_constructor(conf->iter, pileup_constructor);
    bam_mplp_destructor(conf->iter, pileup_destructor);


    // Run mpileup for multiple regions
    int ret = 0;
    if ( nregs )
    {
        int ireg = 0;
        do
        {
            // first region is already positioned
            if ( ireg++ > 0 )
            {
                conf->buf.l = 0;
                ksprintf(&conf->buf,"%s:%u-%u",conf->reg_itr->seq,conf->reg_itr->beg+1,conf->reg_itr->end+1);

                for (i=0; i<conf->nfiles; i++)
                {
                    hts_itr_destroy(conf->mplp_data[i]->iter);
                    conf->mplp_data[i]->iter = sam_itr_querys(conf->mplp_data[i]->idx, conf->mplp_data[i]->h, conf->buf.s);
                    if ( !conf->mplp_data[i]->iter )
                    {
                        conf->mplp_data[i]->iter = sam_itr_querys(conf->mplp_data[i]->idx, conf->mplp_data[i]->h, conf->reg_itr->seq);
                        if ( conf->mplp_data[i]->iter ) {
                            fprintf(stderr,"[E::%s] fail to parse region '%s'\n", __func__, conf->buf.s);
                            exit(EXIT_FAILURE);
                        }
                        fprintf(stderr,"[E::%s] the sequence \"%s\" not found: %s\n",__func__,conf->reg_itr->seq,conf->files[i]);
                        exit(EXIT_FAILURE);
                    }
                    bam_mplp_reset(conf->iter);
                }
            }
            ret = mpileup_reg(conf,conf->reg_itr->beg,conf->reg_itr->end);
            if ( ret<0 ) break;
        }
        while ( regitr_loop(conf->reg_itr) );
    }
    else if ( !conf->reg )
        ret = mpileup_reg(conf,0,UINT32_MAX);
    if ( ret<0 )
    {
        fprintf(stderr, "[%s] failed to read from input file\n", __func__);
        exit(EXIT_FAILURE);
    }

    flush_bcf_records(conf, conf->bcf_fp, conf->bcf_hdr, NULL);

    // clean up
    free(conf->bc.tmp.s);
    bcf_destroy1(conf->bcf_rec);
    if (conf->bcf_fp)
    {
        if ( conf->write_index )
        {
            if ( bcf_idx_save(conf->bcf_fp)<0 )
            {
                if ( hts_close(conf->bcf_fp)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,conf->output_fname);
                error("Error: cannot write to index %s\n",conf->index_fn);
            }
            free(conf->index_fn);
        }
        if ( hts_close(conf->bcf_fp)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,conf->output_fname);
        bcf_hdr_destroy(conf->bcf_hdr);
        bcf_call_destroy(conf->bca);
        free(conf->bc.PL);
        free(conf->bc.DP4);
        free(conf->bc.ADR);
        free(conf->bc.ADF);
        free(conf->bc.SCR);
        free(conf->bc.QS);
        free(conf->bc.ref_nm);
        free(conf->bc.alt_nm);
        free(conf->bc.fmt_arr);
        free(conf->bc.mwu_nm);
        free(conf->bcr);
    }
    if ( conf->gvcf ) gvcf_destroy(conf->gvcf);
    free(conf->buf.s);
    for (i = 0; i < conf->gplp->n; ++i) free(conf->gplp->plp[i]);
    free(conf->gplp->plp); free(conf->gplp->n_plp); free(conf->gplp->m_plp); free(conf->gplp);
    bam_mplp_destroy(conf->iter);
    bam_hdr_destroy(hdr);
    for (i = 0; i < conf->nfiles; ++i) {
        if ( nregs>1 ) hts_idx_destroy(conf->mplp_data[i]->idx);
        sam_close(conf->mplp_data[i]->fp);
        if ( conf->mplp_data[i]->iter) hts_itr_destroy(conf->mplp_data[i]->iter);
        free(conf->mplp_data[i]);
    }
    if ( conf->reg_itr ) regitr_destroy(conf->reg_itr);
    free(conf->mplp_data); free(conf->plp); free(conf->n_plp);
    free(mp_ref.ref[0]);
    free(mp_ref.ref[1]);
    return 0;
}

static int is_url(const char *s)
{
    static const char uri_scheme_chars[] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+.-";
    return s[strspn(s, uri_scheme_chars)] == ':';
}

#define MAX_PATH_LEN 1024
int read_file_list(const char *file_list,int *n,char **argv[])
{
    char buf[MAX_PATH_LEN];
    int len, nfiles = 0;
    char **files = NULL;
    struct stat sb;

    *n = 0;
    *argv = NULL;

    FILE *fh = fopen(file_list,"r");
    if ( !fh )
    {
        fprintf(stderr,"%s: %s\n", file_list,strerror(errno));
        return 1;
    }

    files = (char**) calloc(nfiles,sizeof(char*));
    nfiles = 0;
    while ( fgets(buf,MAX_PATH_LEN,fh) )
    {
        // allow empty lines and trailing spaces
        len = strlen(buf);
        while ( len>0 && isspace(buf[len-1]) ) len--;
        if ( !len ) continue;

        // check sanity of the file list
        buf[len] = 0;
        if (! (is_url(buf) || stat(buf, &sb) == 0))
        {
            // no such file, check if it is safe to print its name
            int i, safe_to_print = 1;
            for (i=0; i<len; i++)
                if (!isprint(buf[i])) { safe_to_print = 0; break; }
            if ( safe_to_print )
                fprintf(stderr,"The file list \"%s\" appears broken, could not locate: %s\n", file_list,buf);
            else
                fprintf(stderr,"Does the file \"%s\" really contain a list of files and do all exist?\n", file_list);
            return 1;
        }

        nfiles++;
        files = (char**) realloc(files,nfiles*sizeof(char*));
        files[nfiles-1] = strdup(buf);
    }
    if ( fclose(fh)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,file_list);
    if ( !nfiles )
    {
        fprintf(stderr,"No files read from %s\n", file_list);
        return 1;
    }
    *argv = files;
    *n    = nfiles;
    return 0;
}
#undef MAX_PATH_LEN

#define SET_FMT_FLAG(str,bit,msg) \
    if (!strcasecmp(tag,str) || !strcasecmp(tag,"FMT/"str) || !strcasecmp(tag,"FORMAT/"str)) \
    { \
        if ( *msg ) fprintf(stderr,"%s",msg); \
        if ( exclude ) \
            *flag &= ~bit; \
        else \
            *flag |= bit; \
        free(tags[i]); \
        continue; \
    }
#define SET_INFO_FLAG(str,bit,msg) if (!strcasecmp(tag,"INFO/"str)) \
    { \
        if ( exclude ) \
            *flag &= ~bit; \
        else \
            *flag |= bit; \
        free(tags[i]); \
        continue; \
    }

void parse_format_flag(uint32_t *flag, const char *str)
{
    int i, n_tags;
    char **tags = hts_readlist(str, 0, &n_tags);
    for(i=0; i<n_tags; i++)
    {
        int exclude = tags[i][0]=='-' ? 1 : 0;
        char *tag = exclude ? tags[i]+1 : tags[i];
        SET_FMT_FLAG("AD", B2B_FMT_AD, "");
        SET_FMT_FLAG("ADF", B2B_FMT_ADF, "");
        SET_FMT_FLAG("ADR", B2B_FMT_ADR, "");
        SET_FMT_FLAG("DP", B2B_FMT_DP, "");
        SET_FMT_FLAG("DP4", B2B_FMT_DP4, "[warning] tag DP4 functional, but deprecated. Please switch to `ADF` and `ADR` in future.\n");
        SET_FMT_FLAG("DPR", B2B_FMT_DPR, "[warning] tag DPR functional, but deprecated. Please switch to `AD` in future.\n");
        SET_FMT_FLAG("DV", B2B_FMT_DV, "[warning] tag DV functional, but deprecated. Please switch to `AD` in future.\n");
        SET_FMT_FLAG("NMBZ", B2B_FMT_NMBZ, "");
        SET_FMT_FLAG("QS", B2B_FMT_QS, "");
        SET_FMT_FLAG("SP", B2B_FMT_SP, "");
        SET_FMT_FLAG("SCR", B2B_FMT_SCR, "");
        SET_INFO_FLAG("DPR", B2B_INFO_DPR, "[warning] tag INFO/DPR functional, but deprecated. Please switch to `INFO/AD` in future.\n");
        SET_INFO_FLAG("AD", B2B_INFO_AD, "");
        SET_INFO_FLAG("ADF", B2B_INFO_ADF, "");
        SET_INFO_FLAG("ADR", B2B_INFO_ADR, "");
        SET_INFO_FLAG("BQBZ", B2B_INFO_BQBZ, "");
        SET_INFO_FLAG("FS", B2B_INFO_FS, "");
        SET_INFO_FLAG("IDV", B2B_INFO_IDV, "");
        SET_INFO_FLAG("IMF", B2B_INFO_IMF, "");
        SET_INFO_FLAG("MIN_PL_SUM", B2B_INFO_MIN_PL_SUM, "");
        SET_INFO_FLAG("MQ0F", B2B_INFO_MQ0F, "");
        SET_INFO_FLAG("MQBZ", B2B_INFO_MQBZ, "");
        SET_INFO_FLAG("NM", B2B_INFO_NM, "");
        SET_INFO_FLAG("NMBZ", B2B_INFO_NMBZ, "");
        SET_INFO_FLAG("RPBZ", B2B_INFO_RPBZ, "");
        SET_INFO_FLAG("SCBZ", B2B_INFO_SCBZ, "");
        SET_INFO_FLAG("SCR", B2B_INFO_SCR, "");
        SET_INFO_FLAG("SGB", B2B_INFO_SGB, "");
        SET_INFO_FLAG("VDB", B2B_INFO_VDB, "");
        fprintf(stderr,"Could not parse tag \"%s\" in \"%s\"\n", tag, str);
        exit(EXIT_FAILURE);
    }
    if (n_tags) free(tags);
}

// todo: make it possible to turn off some annotations or change the defaults,
//      specifically RPB, VDB, MWU, SGB tests. It would be good to do some
//      benchmarking first to see if it's worth it.
static void list_annotations(FILE *fp)
{
    fprintf(fp,
        "Annotations added by default are in this list prefixed with \"*\". To suppress their output, run with\n"
        "e.g. \"-a -FORMAT/AD\".\n"
        "\n"
        "FORMAT annotation tags available (\"FORMAT/\" prefix is optional):\n"
        "\n"
        "* FORMAT/AD   .. Allelic depth (Number=R,Type=Integer)\n"
        "  FORMAT/ADF  .. Allelic depths on the forward strand (Number=R,Type=Integer)\n"
        "  FORMAT/ADR  .. Allelic depths on the reverse strand (Number=R,Type=Integer)\n"
        "  FORMAT/DP   .. Number of high-quality bases (Number=1,Type=Integer)\n"
        "  FORMAT/NMBZ .. Mann-Whitney U-z test of Number of Mismatches within supporting reads (Number=1,Type=Float)\n"
        "  FORMAT/QS   .. Allele phred-score quality sum for use with `call -mG` and +trio-dnm (Number=R,Type=Integer)\n"
        "  FORMAT/SP   .. Phred-scaled strand bias P-value (Number=1,Type=Integer)\n"
        "  FORMAT/SCR  .. Number of soft-clipped reads (Number=1,Type=Integer)\n"
        "\n"
        "INFO annotation tags available:\n"
        "\n"
        "  INFO/AD    .. Total allelic depth (Number=R,Type=Integer)\n"
        "  INFO/ADF   .. Total allelic depths on the forward strand (Number=R,Type=Integer)\n"
        "  INFO/ADR   .. Total allelic depths on the reverse strand (Number=R,Type=Integer)\n"
        "* INFO/BQBZ  .. Mann-Whitney U test of Base Quality Bias (Number=1,Type=Float)\n"
        "  INFO/FS    .. Fisher's exact test P-value to detect strand bias (Number=1,Type=Float)\n"
        "* INFO/IDV   .. Maximum number of raw reads supporting an indel (Number=1,Type=Integer)\n"
        "* INFO/IMF   .. Maximum fraction of raw reads supporting an indel (Number=1,Type=Float)\n"
        "  INFO/MIN_PL_SUM\n"
        "             .. Sum of min PL across all samples before normalization, experimental (Number=1,Type=Integer)\n"
        "* INFO/MQ0F  .. Fraction of reads with zero mapping quality (Number=1,Type=Float)\n"
        "* INFO/MQBZ  .. Mann-Whitney U test of Mapping Quality Bias (Number=1,Type=Float)\n"
        "* INFO/MQSBZ .. Mann-Whitney U-z test of Mapping Quality vs Strand Bias (Number=1,Type=Float)\n"
        "  INFO/NM    .. Approximate average number of mismatches in ref and alt reads, experimental (Number=2,Type=Float)\n"
        "  INFO/NMBZ  .. Mann-Whitney U-z test of Number of Mismatches within supporting reads (Number=1,Type=Float)\n"
        "* INFO/RPBZ  .. Mann-Whitney U test of Read Position Bias (Number=1,Type=Float)\n"
        "* INFO/SCBZ  .. Mann-Whitney U-z test of Soft-Clip Length Bias (Number=1,Type=Float)\n"
        "  INFO/SCR   .. Number of soft-clipped reads (Number=1,Type=Integer)\n"
        "* INFO/SGB   .. Segregation based metric, http://samtools.github.io/bcftools/rd-SegBias.pdf (Number=1,Type=Float)\n"
        "* INFO/VDB   .. Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (Number=1,Type=Float)\n"
        "\n");
}

static void print_usage(FILE *fp, const mplp_conf_t *mplp)
{
    char *tmp_skip_all_set = bam_flag2str(mplp->rflag_skip_all_set);
    char *tmp_skip_any_unset = bam_flag2str(mplp->rflag_skip_any_unset);
    char *tmp_skip_all_unset = bam_flag2str(mplp->rflag_skip_all_unset);
    char *tmp_skip_any_set = bam_flag2str(mplp->rflag_skip_any_set);

    // Display usage information, formatted for the standard 80 columns.
    // (The unusual string formatting here aids the readability of this
    // source code in 80 columns, to the extent that's possible.)

    fprintf(fp,
        "\n"
        "Usage: bcftools mpileup [options] in1.bam [in2.bam [...]]\n"
        "\n"
        "Input options:\n"
        "  -6, --illumina1.3+      Quality is in the Illumina-1.3+ encoding\n"
        "  -A, --count-orphans     Include anomalous read pairs, with flag PAIRED but not PROPER_PAIR set\n"
        "  -b, --bam-list FILE     List of input BAM filenames, one per line\n"
        "  -B, --no-BAQ            Disable BAQ (per-Base Alignment Quality)\n"
        "  -C, --adjust-MQ INT     Adjust mapping quality [0]\n"
        "  -D, --full-BAQ          Apply BAQ everywhere, not just in problematic regions\n"
        "  -d, --max-depth INT     Max raw per-file depth; avoids excessive memory usage [%d]\n", mplp->max_depth);
            fprintf(fp,
        "  -E, --redo-BAQ          Recalculate BAQ on the fly, ignore existing BQs\n"
        "  -f, --fasta-ref FILE    Faidx indexed reference sequence file\n"
        "      --no-reference      Do not require fasta reference file\n"
        "  -G, --read-groups FILE  Select or exclude read groups listed in the file\n"
        "  -q, --min-MQ INT        Skip alignments with mapQ smaller than INT [%d]\n", mplp->min_mq);
    fprintf(fp,
        "  -Q, --min-BQ INT        Skip bases with baseQ/BAQ smaller than INT [%d]\n", mplp->min_baseQ);
    fprintf(fp,
        "      --max-BQ INT        Limit baseQ/BAQ to no more than INT [%d]\n", mplp->max_baseQ);
    fprintf(fp,
        "      --delta-BQ INT      Use neighbour_qual + INT if less than qual [%d]\n", mplp->delta_baseQ);
    fprintf(fp,
        "  -r, --regions REG[,...] Comma separated list of regions in which pileup is generated\n"
        "  -R, --regions-file FILE Restrict to regions listed in a file\n"
        "      --ignore-RG         Ignore RG tags (one BAM = one sample)\n"
        "  --ls, --skip-all-set STR|INT  Skip reads with all of the bits set []\n");
    fprintf(fp,
        "  --ns, --skip-any-set STR|INT  Skip reads with any of the bits set [%s]\n", tmp_skip_any_set);
    fprintf(fp,
        "  --lu, --skip-all-unset STR|INT  Skip reads with all of the bits unset []\n"
        "  --nu, --skip-any-unset STR|INT  Skip reads with any of the bits unset []\n");
    fprintf(fp,
        "  -s, --samples LIST      Comma separated list of samples to include\n"
        "  -S, --samples-file FILE File of samples to include\n"
        "  -t, --targets REG[,...] Similar to -r but streams rather than index-jumps\n"
        "  -T, --targets-file FILE Similar to -R but streams rather than index-jumps\n"
        "  -x, --ignore-overlaps   Disable read-pair overlap detection\n"
        "      --seed INT          Random number seed used for sampling deep regions [0]\n"
        "\n"
        "Output options:\n"
        "  -a, --annotate LIST     Optional tags to output; '\\?' to list available tags []\n"
        "  -g, --gvcf INT[,...]    Group non-variant sites into gVCF blocks according\n"
        "                          To minimum per-sample DP\n"
        "      --no-version        Do not append version and command line to the header\n"
        "  -o, --output FILE       Write output to FILE [standard output]\n"
        "  -O, --output-type TYPE  'b' compressed BCF; 'u' uncompressed BCF;\n"
        "                          'z' compressed VCF; 'v' uncompressed VCF; 0-9 compression level [v]\n"
        "      --threads INT       Use multithreading with INT worker threads [0]\n"
        "  -W, --write-index[=FMT] Automatically index the output files [off]\n"
        "\n"
        "SNP/INDEL genotype likelihoods options:\n"
        "  -X, --config STR        Specify platform profile (use \"-X list\" for details)\n"
        "  -e, --ext-prob INT      Phred-scaled gap extension seq error probability [%d]\n", mplp->extQ);
    fprintf(fp,
        "  -F, --gap-frac FLOAT    Minimum fraction of gapped reads [%g]\n", mplp->min_frac);
    fprintf(fp,
        "  -h, --tandem-qual INT   Coefficient for homopolymer errors [%d]\n", mplp->tandemQ);
    fprintf(fp,
        "  -I, --skip-indels       Do not perform indel calling\n"
        "  -L, --max-idepth INT    Maximum per-file depth for INDEL calling [%d]\n", mplp->max_indel_depth);
    fprintf(fp,
        "  -m, --min-ireads INT    Minimum number gapped reads for indel candidates [%d]\n", mplp->min_support);
    fprintf(fp,
        "  -M, --max-read-len INT  Maximum length of read to pass to BAQ algorithm [%d]\n", mplp->max_read_len);
    fprintf(fp,
        "  -o, --open-prob INT     Phred-scaled gap open seq error probability [%d]\n", mplp->openQ);
    fprintf(fp,
        "  -p, --per-sample-mF     Apply -m and -F per-sample for increased sensitivity\n"
        "  -P, --platforms STR     Comma separated list of platforms for indels [all]\n"
        "  --ar, --ambig-reads STR   What to do with ambiguous indel reads: drop,incAD,incAD0 [drop]\n");
    fprintf(fp,
        "      --indel-bias FLOAT  Raise to favour recall over precision [%.2f]\n", mplp->indel_bias);
    fprintf(fp,
        "      --del-bias FLOAT    Relative likelihood of insertion to deletion [%.2f]\n", mplp->del_bias);
    fprintf(fp,
        "      --score-vs-ref FLOAT\n"
        "                          Ratio of score vs ref (1) or 2nd-best allele (0) [%.2f]\n", mplp->vs_ref);
    fprintf(fp,
        "      --indel-size INT    Approximate maximum indel size considered [%d]\n", mplp->indel_win_size);
    fprintf(fp,
        "      --indels-2.0        New EXPERIMENTAL indel calling model (diploid reference consensus)\n"
        "      --indels-cns        New EXPERIMENTAL indel calling model with edlib\n"
        "      --seqq-offset       Indel-cns tuning for indel seq-qual scores [120]\n"
        "      --no-indels-cns     Disable CNS mode, to use after a -X profile\n"
        "      --poly-mqual        (Edlib mode) Use minimum quality within homopolymers\n");
    fprintf(fp,"\n");
    fprintf(fp,
            "Notes: Assuming diploid individuals.\n\n"
            "Example:\n"
            "   # See also http://samtools.github.io/bcftools/howtos/variant-calling.html\n"
            "   bcftools mpileup -Ou -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf\n"
            "\n");

    free(tmp_skip_all_set);
    free(tmp_skip_any_unset);
    free(tmp_skip_all_unset);
    free(tmp_skip_any_set);
}

static void print_profiles(void) {
    printf(
"Configuration profiles activated with -X, --config:\n\n"
"1.12\n"
"    -Q13 -h100 -m1 -F0.002\n\n"
"bgi, bgi-1.20\n"
"    --indels-cns -B --indel-size 80 -F0.1 --indel-bias 0.9 --seqq-offset 120\n\n"
"illumina-1.18\n"
"    --indel-size 110\n\n"
"illumina\n"
"illumina-1.20\n"
"    --indels-cns --indel-size 110\n\n"
"ont\n"
"    -B -Q5 --max-BQ 30 -I\n\n"
"ont-sup, ont-sup-1.20\n"
"    --indels-cns -B -Q1 --max-BQ 35 -F0.2 -o15 -e1 -h110 --delta-BQ 99\\\n"
"    --del-bias 0.4 --indel-bias 0.7 --poly-mqual --seqq-offset 130\\\n"
"    --indel-size 80\n\n"
"pacbio-ccs-1.18\n"
"    -D -Q5 --max-BQ 50 -F0.1 -o25 -e1 --delta-BQ 10 \\\n"
"    -M99999 --indel-size 110\n\n"
"pacbio-ccs, pacbio-ccs-1.20\n"
"    --indels-cns -B -Q5 --max-BQ 50 -F0.1 -o25 -e1 -h300 --delta-BQ 10 \\\n"
"    --del-bias 0.4 --poly-mqual --indel-bias 0.9 --seqq-offset 118\\\n"
"    --indel-size 80 --score-vs-ref 0.7\n\n"
"ultima, ultima-1.20\n"
"    --indels-cns -B -Q1 --max-BQ 30 -F0.15 -o20 -e10 -h250 --delta-BQ 10 \\\n"
"    --del-bias 0.3 --indel-bias 0.7 --poly-mqual --seqq-offset 140 \\\n"
"    --indel-size 80 --score-vs-ref 0.3\n\n"
"\n");
}

int main_mpileup(int argc, char *argv[])
{
    int c, i, ret = 1;
    const char *file_list = NULL;
    char **fn = NULL;
    int nfiles = 0, use_orphan = 0, noref = 0;
    mplp_conf_t mplp;
    memset(&mplp, 0, sizeof(mplp_conf_t));
    mplp.min_baseQ = 1;
    mplp.max_baseQ = 60;
    mplp.delta_baseQ = 30;
    mplp.capQ_thres = 0;
    mplp.max_depth = 250; mplp.max_indel_depth = 250;
    mplp.openQ = 40; mplp.extQ = 20; mplp.tandemQ = 500;
    mplp.min_frac = 0.05; mplp.indel_bias = 1.0; mplp.min_support = 2;
    mplp.vs_ref = 0;
    mplp.flag = MPLP_NO_ORPHAN | MPLP_REALN | MPLP_REALN_PARTIAL
              | MPLP_SMART_OVERLAPS;
    mplp.argc = argc; mplp.argv = argv;
    mplp.rflag_skip_any_set = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
    mplp.output_fname = NULL;
    mplp.output_type = FT_VCF;
    mplp.record_cmd_line = 1;
    mplp.n_threads = 0;
    mplp.bsmpl = bam_smpl_init();
    // the default to be changed in future, see also parse_format_flag()
    mplp.fmt_flag = B2B_INFO_BQBZ|B2B_INFO_IDV|B2B_INFO_IMF|B2B_INFO_MQ0F|B2B_INFO_MQBZ|B2B_INFO_MQSBZ|B2B_INFO_RPBZ|B2B_INFO_SCBZ|B2B_INFO_SGB|B2B_INFO_VDB|B2B_FMT_AD;
    mplp.max_read_len = 500;
    mplp.ambig_reads = B2B_DROP;
    mplp.indel_win_size = 110;
    mplp.poly_mqual = 0;
    mplp.seqQ_offset = 120;
    mplp.clevel = -1;
    mplp.del_bias = 0; // even insertion and deletion likelhoods.
    hts_srand48(0);

    static const struct option lopts[] =
    {
        {"nu", required_argument, NULL, 16},
        {"lu", required_argument, NULL, 17},
        {"rf", required_argument, NULL, 17},   // old --rf, --incl-flags = --lu, --skip-all-unset
        {"ns", required_argument, NULL, 18},
        {"ff", required_argument, NULL, 18},   // old --ff, --excl-flags = --ns, --skip-any-set
        {"ls", required_argument, NULL, 19},
        {"skip-any-unset", required_argument, NULL, 16},
        {"skip-all-unset", required_argument, NULL, 17},
        {"skip-any-set", required_argument, NULL, 18},
        {"skip-all-set", required_argument, NULL, 19},
        {"output", required_argument, NULL, 3},
        {"open-prob", required_argument, NULL, 4},
        {"ignore-RG", no_argument, NULL, 5},
        {"ignore-rg", no_argument, NULL, 5},
        {"gvcf", required_argument, NULL, 'g'},
        {"no-reference", no_argument, NULL, 7},
        {"no-version", no_argument, NULL, 8},
        {"threads",required_argument,NULL,9},
        {"illumina1.3+", no_argument, NULL, '6'},
        {"count-orphans", no_argument, NULL, 'A'},
        {"bam-list", required_argument, NULL, 'b'},
        {"no-BAQ", no_argument, NULL, 'B'},
        {"no-baq", no_argument, NULL, 'B'},
        {"full-BAQ", no_argument, NULL, 'D'},
        {"full-baq", no_argument, NULL, 'D'},
        {"adjust-MQ", required_argument, NULL, 'C'},
        {"adjust-mq", required_argument, NULL, 'C'},
        {"max-depth", required_argument, NULL, 'd'},
        {"redo-BAQ", no_argument, NULL, 'E'},
        {"redo-baq", no_argument, NULL, 'E'},
        {"fasta-ref", required_argument, NULL, 'f'},
        {"read-groups", required_argument, NULL, 'G'},
        {"region", required_argument, NULL, 'r'},
        {"regions", required_argument, NULL, 'r'},
        {"regions-file", required_argument, NULL, 'R'},
        {"targets", required_argument, NULL, 't'},
        {"targets-file", required_argument, NULL, 'T'},
        {"min-MQ", required_argument, NULL, 'q'},
        {"min-mq", required_argument, NULL, 'q'},
        {"min-BQ", required_argument, NULL, 'Q'},
        {"min-bq", required_argument, NULL, 'Q'},
        {"max-bq", required_argument, NULL, 11},
        {"max-BQ", required_argument, NULL, 11},
        {"delta-BQ", required_argument, NULL, 12},
        {"ignore-overlaps", no_argument, NULL, 'x'},
        {"output-type", required_argument, NULL, 'O'},
        {"samples", required_argument, NULL, 's'},
        {"samples-file", required_argument, NULL, 'S'},
        {"annotate", required_argument, NULL, 'a'},
        {"ext-prob", required_argument, NULL, 'e'},
        {"gap-frac", required_argument, NULL, 'F'},
        {"indel-bias", required_argument, NULL, 10},
        {"indel-size", required_argument, NULL, 15},
        {"indels-2.0", no_argument, NULL, 20},
        {"indels-cns", no_argument, NULL, 22},
        {"no-indels-cns", no_argument, NULL, 25},
        {"tandem-qual", required_argument, NULL, 'h'},
        {"skip-indels", no_argument, NULL, 'I'},
        {"max-idepth", required_argument, NULL, 'L'},
        {"min-ireads", required_argument, NULL, 'm'},
        {"per-sample-mF", no_argument, NULL, 'p'},
        {"per-sample-mf", no_argument, NULL, 'p'},
        {"platforms", required_argument, NULL, 'P'},
        {"max-read-len", required_argument, NULL, 'M'},
        {"config", required_argument, NULL, 'X'},
        {"seed", required_argument, NULL, 13},
        {"ambig-reads", required_argument, NULL, 14},
        {"ar", required_argument, NULL, 14},
        {"write-index",optional_argument,NULL,'W'},
        {"del-bias", required_argument, NULL, 23},
        {"poly-mqual", no_argument, NULL, 24},
        {"no-poly-mqual", no_argument, NULL, 26},
        {"score-vs-ref",required_argument, NULL, 27},
        {"seqq-offset", required_argument, NULL, 28},
        {NULL, 0, NULL, 0}
    };
    while ((c = getopt_long(argc, argv, "Ag:f:r:R:q:Q:C:BDd:L:b:P:po:e:h:Im:F:EG:6O:xa:s:S:t:T:M:X:UW::",lopts,NULL)) >= 0) {
        switch (c) {
        case 'x': mplp.flag &= ~MPLP_SMART_OVERLAPS; break;
        case  16 :
            mplp.rflag_skip_any_unset = bam_str2flag(optarg);
            if ( mplp.rflag_skip_any_unset <0 ) {
                fprintf(stderr,"Could not parse --nf %s\n", optarg);
                goto err;
            }
            break;
        case  17 :
            mplp.rflag_skip_all_unset = bam_str2flag(optarg);
            if ( mplp.rflag_skip_all_unset<0 ) {
                fprintf(stderr,"Could not parse --if %s\n", optarg);
                goto err;
            }
            break;
        case  18 :
            mplp.rflag_skip_any_set = bam_str2flag(optarg);
            if ( mplp.rflag_skip_any_set <0 ) {
                fprintf(stderr,"Could not parse --ef %s\n", optarg);
                goto err;
            }
            break;
        case  19 :
            mplp.rflag_skip_all_set = bam_str2flag(optarg);
            if ( mplp.rflag_skip_all_set <0 ) {
                fprintf(stderr,"Could not parse --df %s\n", optarg);
                goto err;
            }
            break;
        case  3 : mplp.output_fname = optarg; break;
        case  4 : mplp.openQ = atoi(optarg); break;
        case  5 : bam_smpl_ignore_readgroups(mplp.bsmpl); break;
        case 'g':
            mplp.gvcf = gvcf_init(optarg);
            if ( !mplp.gvcf ) error("Could not parse: --gvcf %s\n", optarg);
            break;
        case 'f':
            mplp.fai = fai_load(optarg);
            if (mplp.fai == NULL)
                goto err;
            mplp.fai_fname = optarg;
            break;
        case  7 : noref = 1; break;
        case  8 : mplp.record_cmd_line = 0; break;
        case  9 : mplp.n_threads = strtol(optarg, 0, 0); break;
        case 'd': mplp.max_depth = atoi(optarg); break;
        case 'r': mplp.reg_fname = strdup(optarg); break;
        case 'R': mplp.reg_fname = strdup(optarg); mplp.reg_is_file = 1; break;
        case 't':
                  // In the original version the whole BAM was streamed which is inefficient
                  //  with few BED intervals and big BAMs. Todo: devise a heuristic to determine
                  //  best strategy, that is streaming or jumping.
                  if ( optarg[0]=='^' ) optarg++;
                  else mplp.bed_logic = 1;
                  mplp.bed = regidx_init(NULL,regidx_parse_reg,NULL,0,NULL);
                  mplp.bed_itr = regitr_init(mplp.bed);
                  if ( regidx_insert_list(mplp.bed,optarg,',') !=0 )
                  {
                      fprintf(stderr,"Could not parse the targets: %s\n", optarg);
                      exit(EXIT_FAILURE);
                  }
                  break;
        case 'T':
                  if ( optarg[0]=='^' ) optarg++;
                  else mplp.bed_logic = 1;
                  mplp.bed = regidx_init(optarg,NULL,NULL,0,NULL);
                  if (!mplp.bed) {
                      fprintf(stderr, "bcftools mpileup: Could not read file \"%s\"", optarg);
                      goto err;
                  }
                  break;
        case 'P': mplp.pl_list = strdup(optarg); break;
        case 'p': mplp.flag |= MPLP_PER_SAMPLE; break;
        case 'B': mplp.flag &= ~MPLP_REALN; break;
        case 'D': mplp.flag &= ~MPLP_REALN_PARTIAL; break;
        case 'I': mplp.flag |= MPLP_NO_INDEL; break;
        case 'E': mplp.flag |= MPLP_REDO_BAQ; break;
        case '6': mplp.flag |= MPLP_ILLUMINA13; break;
        case 's': if ( bam_smpl_add_samples(mplp.bsmpl,optarg,0)<0 ) error("Could not read samples: %s\n",optarg); break;
        case 'S': if ( bam_smpl_add_samples(mplp.bsmpl,optarg,1)<0 ) error("Could not read samples: %s\n",optarg); break;
        case 'O':
            switch (optarg[0]) {
                case 'b': mplp.output_type = FT_BCF_GZ; break;
                case 'u': mplp.output_type = FT_BCF; break;
                case 'z': mplp.output_type = FT_VCF_GZ; break;
                case 'v': mplp.output_type = FT_VCF; break;
                default:
                {
                    char *tmp;
                    mplp.clevel = strtol(optarg,&tmp,10);
                    if ( *tmp || mplp.clevel<0 || mplp.clevel>9 ) error("The output type \"%s\" not recognised\n", optarg);
                }
            }
            if ( optarg[1] )
            {
                char *tmp;
                mplp.clevel = strtol(optarg+1,&tmp,10);
                if ( *tmp || mplp.clevel<0 || mplp.clevel>9 ) error("Could not parse argument: --output-type %s\n", optarg+1);
            }
            break;
        case 'C': mplp.capQ_thres = atoi(optarg); break;
        case 'q': mplp.min_mq = atoi(optarg); break;
        case 'Q': mplp.min_baseQ = atoi(optarg); break;
        case  11: mplp.max_baseQ = atoi(optarg); break;
        case  12: mplp.delta_baseQ = atoi(optarg); break;
        case 'b': file_list = optarg; break;
        case 'o': {
                char *end;
                long value = strtol(optarg, &end, 10);
                // Distinguish between -o INT and -o FILE (a bit of a hack!)
                if (*end == '\0') mplp.openQ = value;
                else mplp.output_fname = optarg;
            }
            break;
        case 'e': mplp.extQ = atoi(optarg); break;
        case 'h': mplp.tandemQ = atoi(optarg); break;
        case 10: // --indel-bias (inverted so higher => more indels called)
            if (atof(optarg) < 1e-2)
                mplp.indel_bias = 1/1e2;
            else
                mplp.indel_bias = 1/atof(optarg);
            break;
        case 27:
            mplp.vs_ref = atof(optarg);
            //if (mplp.vs_ref < 0) mplp.vs_ref = 0;
            if (mplp.vs_ref > 1) mplp.vs_ref = 1;
            break;
        case  15: {
                char *tmp;
                mplp.indel_win_size = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --indel-size %s\n", optarg);
                if ( mplp.indel_win_size < 20 )
                {
                    mplp.indel_win_size = 20;
                    fprintf(stderr,"Warning: running with --indel-size %d, the requested value is too small\n",mplp.indel_win_size);
                }
            }
            break;
        case  20: mplp.indels_v20 = 1; mplp.edlib = 0; break;
        case 'W':
            if (!(mplp.write_index = write_index_parse(optarg)))
                error("Unsupported index format '%s'\n", optarg);
            break;
        case  22: mplp.edlib = 1; mplp.indels_v20 = 0; break;
        case  25: mplp.edlib = 0; break;
        case  28:
            mplp.seqQ_offset = atoi(optarg);
            if (mplp.seqQ_offset < 100)
                mplp.seqQ_offset = 100;
            if (mplp.seqQ_offset > 200)
                mplp.seqQ_offset = 200;
            break;
        case  23: mplp.del_bias = atof(optarg); break;
        case  24: mplp.poly_mqual = 1; break;
        case  26: mplp.poly_mqual = 0; break;
        case 'A': use_orphan = 1; break;
        case 'F': mplp.min_frac = atof(optarg); break;
        case 'm': mplp.min_support = atoi(optarg); break;
        case 'L': mplp.max_indel_depth = atoi(optarg); break;
        case 'G': bam_smpl_add_readgroups(mplp.bsmpl, optarg, 1); break;
        case 'a':
            if (optarg[0]=='?') {
                list_annotations(stderr);
                goto err;
            }
            parse_format_flag(&mplp.fmt_flag,optarg);
        break;
        case 'M': mplp.max_read_len = atoi(optarg); break;
        case 'X':
            if (strcasecmp(optarg, "pacbio-ccs-1.18") == 0) {
                mplp.min_frac = 0.1;
                mplp.min_baseQ = 5;
                mplp.max_baseQ = 50;
                mplp.delta_baseQ = 10;
                mplp.openQ = 25;
                mplp.extQ = 1;
                mplp.flag |= MPLP_REALN_PARTIAL;
                mplp.max_read_len = 99999;

            } else if (strcasecmp(optarg, "pacbio-ccs") == 0 ||
                strcasecmp(optarg, "pacbio-ccs-1.20") == 0) {
                mplp.min_frac = 0.1;
                mplp.min_baseQ = 5;
                mplp.max_baseQ = 50;
                mplp.delta_baseQ = 10;
                mplp.tandemQ = 300;
                mplp.openQ = 25;
                mplp.extQ = 1;
                mplp.flag &= ~MPLP_REALN;
                mplp.del_bias = 0.4;
                mplp.indel_bias = 1/.9;
                mplp.seqQ_offset = 118;
                mplp.poly_mqual = 1;
                mplp.edlib = 1;
                mplp.vs_ref = 0.7;
                mplp.indel_win_size = 80;

            } else if (strcasecmp(optarg, "ont") == 0) {
                fprintf(stderr, "With old ONT data may be beneficial to also run bcftools call with "
                        "a higher -P, eg -P0.01 or -P 0.1\n");
                mplp.min_baseQ = 5;
                mplp.max_baseQ = 30;
                mplp.flag &= ~MPLP_REALN;
                mplp.flag |= MPLP_NO_INDEL;

            } else if (strcasecmp(optarg, "ont-sup") == 0 ||
                       strcasecmp(optarg, "ont-sup-1.20") == 0) {
                mplp.min_frac = 0.2;
                mplp.min_baseQ = 1;
                mplp.max_baseQ = 35;
                mplp.delta_baseQ = 99;
                mplp.openQ = 15;
                mplp.extQ = 1;
                mplp.flag &= ~MPLP_REALN;
                mplp.max_read_len = 9999999;
                mplp.del_bias = 0.4;
                mplp.poly_mqual = 1;
                mplp.edlib = 1;
                // If we increase -h then we can increase bias denominator too
                mplp.tandemQ = 110;
                mplp.indel_bias = 1/0.7;
                mplp.seqQ_offset = 130;
                mplp.indel_win_size = 80;

            } else if (strcasecmp(optarg, "ultima") == 0 ||
                       strcasecmp(optarg, "ultima-1.20") == 0) {
                mplp.min_frac = 0.15;
                mplp.min_baseQ = 1;
                mplp.max_baseQ = 30;
                mplp.delta_baseQ = 10;
                mplp.openQ = 20;
                mplp.extQ = 10;
                mplp.tandemQ = 250;
                mplp.flag &= ~MPLP_REALN;
                mplp.del_bias = 0.3;
                mplp.poly_mqual = 1;
                mplp.edlib = 1;
                mplp.indel_bias = 1/0.7;
                mplp.seqQ_offset = 140;
                mplp.vs_ref = 0.3;
                mplp.indel_win_size = 80;

            } else if (strcasecmp(optarg, "1.12") == 0) {
                // 1.12 and earlier
                mplp.min_frac = 0.002;
                mplp.min_support = 1;
                mplp.min_baseQ = 13;
                mplp.tandemQ = 100;
                mplp.flag &= ~MPLP_REALN_PARTIAL;
                mplp.flag |= MPLP_REALN;

            } else if (strcasecmp(optarg, "illumina-1.18") == 0) {
                mplp.indel_win_size = 110;
                mplp.flag |= MPLP_REALN_PARTIAL;

            } else if (strcasecmp(optarg, "illumina") == 0 ||
                       strcasecmp(optarg, "illumina-1.20") == 0) {
                mplp.edlib = 1;
                mplp.indel_win_size = 110;
                mplp.flag |= MPLP_REALN_PARTIAL;
                mplp.indel_bias = 1;
                mplp.seqQ_offset = 125;
                //mplp.indel_win_size = 80; TEST?

            } else if (strcasecmp(optarg, "bgi") == 0 ||
                       strcasecmp(optarg, "bgi-1.20") == 0) {
                mplp.min_frac = 0.1;
                mplp.edlib = 1;
                mplp.indel_bias = 1;
                mplp.seqQ_offset = 120;
                mplp.flag |= MPLP_REALN_PARTIAL;
                mplp.indel_win_size = 80;

            } else if (strcasecmp(optarg, "list") == 0 ||
                       strcasecmp(optarg, "help") == 0) {
                print_profiles();
                goto err;
            } else {
                fprintf(stderr, "Unknown configuration name '%s'\n"
                        "Please use '-X list' to show available choices.\n",
                        optarg);
                goto err;
            }
            break;
        case 13: hts_srand48(atoi(optarg)); break;
        case 14:
            if ( !strcasecmp(optarg,"drop") ) mplp.ambig_reads = B2B_DROP;
            else if ( !strcasecmp(optarg,"incAD") ) mplp.ambig_reads = B2B_INC_AD;
            else if ( !strcasecmp(optarg,"incAD0") ) mplp.ambig_reads = B2B_INC_AD0;
            else error("The option to --ambig-reads not recognised: %s\n",optarg);
            break;
        default:
            fprintf(stderr,"Invalid option: '%c'\n", c);
            goto err;
        }
    }

    if ( mplp.gvcf && !(mplp.fmt_flag&B2B_FMT_DP) )
    {
        fprintf(stderr,"[warning] The -a DP option is required with --gvcf, switching on.\n");
        mplp.fmt_flag |= B2B_FMT_DP;
    }
    if ( mplp.flag&(MPLP_BCF|MPLP_VCF|MPLP_NO_COMP) )
    {
        if ( mplp.flag&MPLP_VCF )
        {
            if ( mplp.flag&MPLP_NO_COMP ) mplp.output_type = FT_VCF;
            else mplp.output_type = FT_VCF_GZ;
        }
        else if ( mplp.flag&MPLP_BCF )
        {
            if ( mplp.flag&MPLP_NO_COMP ) mplp.output_type = FT_BCF;
            else mplp.output_type = FT_BCF_GZ;
        }
    }
    if ( !(mplp.flag&MPLP_REALN) && mplp.flag&MPLP_REDO_BAQ )
    {
        fprintf(stderr,"Error: The -B option cannot be combined with -E\n");
        goto err;
    }
    if (use_orphan) mplp.flag &= ~MPLP_NO_ORPHAN;
    if (argc == 1)
    {
        print_usage(stderr, &mplp);
        goto err;
    }
    if (!mplp.fai && !noref) {
        fprintf(stderr,"Error: mpileup requires the --fasta-ref option by default; use --no-reference to run without a fasta reference\n");
        goto err;
    }

    if (file_list)
    {
        if ( read_file_list(file_list,&nfiles,&fn) )
            goto err;
        mplp.files  = fn;
        mplp.nfiles = nfiles;
    }
    else
    {
        mplp.nfiles = argc - optind;
        mplp.files  = (char**) malloc(mplp.nfiles*sizeof(char*));
        for (i=0; i<mplp.nfiles; i++) mplp.files[i] = strdup(argv[optind+i]);
    }
    ret = mpileup(&mplp);

    for (i=0; i<mplp.nfiles; i++) free(mplp.files[i]);
    free(mplp.files);
    free(mplp.reg_fname); free(mplp.pl_list);
    if (mplp.fai) fai_destroy(mplp.fai);
    if (mplp.bed) regidx_destroy(mplp.bed);
    if (mplp.bed_itr) regitr_destroy(mplp.bed_itr);
    if (mplp.reg) regidx_destroy(mplp.reg);

 err:
    bam_smpl_destroy(mplp.bsmpl);

    return ret;
}
