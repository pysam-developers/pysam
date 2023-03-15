#include "bcftools.pysam.h"

/*  bam2bcf_iaux.c -- modified indel caller

    Copyright (C) 2022 Genome Research Ltd.

    Author: pd3@sanger, jkb

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
    DEALINGS IN THE SOFTWARE
*/

#include <assert.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/khash_str2int.h>
#include "bcftools.h"
#include "bam2bcf.h"
#include "read_consensus.h"
#include "cigar_state.h"

#include <htslib/ksort.h>
KSORT_INIT_STATIC_GENERIC(uint32_t)

#ifndef DEBUG_ALN
#define DEBUG_ALN 0
#endif

#define MAX_TYPES 64

typedef struct
{
    int pos;    // current position
    char *chr;  // current chromosome
    int nsmpl;  // number of samples
    int *nplp;              // per-sample number of reads
    bam_pileup1_t **plp;    // per-sample reads
    bcf_callaux_t *bca;     // auxiliary bam2bcf structure
    const char *ref;        // reference genome (ASCII)
    uint32_t *uitmp;        // temporary unsigned int array
    char *inscns;           // insertions consensus "ACGTN"[itype*max_ins_len+i]
    int muitmp, minscns;    // size of uitmp, inscns
    int iref_type, ntypes, types[MAX_TYPES];   // indel types
    int max_ins_len;        // largest insertion
    int left, right;        // consensus sequence boundaries, 0-based fa ref coordinates
    read_cns_t *rcns;       // read consensus
    cns_seq_t *cns_seq;     // array of consensus sequences
    int *cns_pos;           // array of relative pos indexes within cns_seq sequences
    uint8_t *ref_seq, *qry_seq; // reference and query sequence to align
    int nref_seq, nqry_seq;     // the allocated size of ref_seq and qry_seq
    uint8_t *qual;
    int nqual;
    int *read_scores,           // read scores for each indel type [ntypes*iread+itype]
        mread_scores,
        ref_qual[MAX_TYPES],    // refseq quality at pos for each indel type in the context of homopolymer runs
        sum_qual[MAX_TYPES];    // qual contributions to each indel type from all reads
}
indel_aux_t;

#if DEBUG_ALN
static void debug_print_types(indel_aux_t *iaux)
{
    int i,j;
    fprintf(bcftools_stderr,"types at %s:%d ntypes=%d... ",iaux->chr,iaux->pos+1,iaux->ntypes);
    for (i=0; i<iaux->ntypes; i++)
    {
        fprintf(bcftools_stderr," type%d=",i);
        if ( iaux->types[i]<=0 )
        {
            if ( i==iaux->iref_type ) fprintf(bcftools_stderr,"%d(ref)",iaux->types[i]);
            else fprintf(bcftools_stderr,"%d",iaux->types[i]);
            continue;
        }
        char *cns = &iaux->inscns[i*iaux->max_ins_len];
        for (j=0; j<iaux->types[i]; j++) fprintf(bcftools_stderr,"%c","ACGTN"[(int)cns[j]]);
    }
    fprintf(bcftools_stderr,"\n");
}
#else
#define debug_print_types(iaux)
#endif

void bcf_iaux_destroy(bcf_callaux_t *bca)
{
    if ( !bca->iaux ) return;
    indel_aux_t *iaux = (indel_aux_t*)bca->iaux;
    free(iaux->uitmp);
    free(iaux->inscns);
    free(iaux->ref_seq);
    free(iaux->qry_seq);
    free(iaux->qual);
    free(iaux->read_scores);
    rcns_destroy(iaux->rcns);
    free(iaux);
}

static void iaux_init_sequence_context(indel_aux_t *iaux)
{
    // Calculate left and right boundary. The array types is sorted in ascending order, the first
    // element is the largest deletion (if a deletion present)
    iaux->left  = iaux->pos > iaux->bca->indel_win_size ? iaux->pos - iaux->bca->indel_win_size : 0;
    iaux->right = iaux->pos + iaux->bca->indel_win_size;
    if ( iaux->types[0] < 0 ) iaux->right -= iaux->types[0];    // extend by the largest deletion length

    // In case the alignments stand out the reference
    int i;
    for (i=iaux->pos; i<iaux->right; i++)
        if ( !iaux->ref[i] ) break;
    iaux->right = i;

    // Sequence quality in the context of homopolymers for each indel type
    int l_run = bcf_cgp_l_run(iaux->ref, iaux->pos);    // The length of the homopolymer run around the current position
    for (i=0; i<iaux->ntypes; i++)
    {
        int l = iaux->types[i];

        // This is the original est_seqQ() code. FIXME: check if the inserted sequence is consistent with the homopolymer run
        int q  = iaux->bca->openQ + iaux->bca->extQ * (abs(l) - 1);
        int qh = l_run >= 3? (int)(iaux->bca->tandemQ * (double)abs(l) / l_run + .499) : 1000;
        if ( q > qh ) q = qh;

        iaux->ref_qual[i] = q < 255 ? q : 255;
    }

    // Determine the indel region, this makes the difference between e.g. T>TA vs TA>TAA
    iaux->bca->indelreg = 0;
    for (i=0; i<iaux->ntypes; i++)
    {
        if ( !iaux->types[i] ) continue;
        int ireg;
        if ( iaux->types[i] > 0 )
            ireg = est_indelreg(iaux->pos, iaux->ref, iaux->types[i], &iaux->inscns[i*iaux->max_ins_len]);
        else
            ireg = est_indelreg(iaux->pos, iaux->ref, -iaux->types[i], 0);
        if ( ireg > iaux->bca->indelreg ) iaux->bca->indelreg = ireg;
    }
}

static int iaux_init_scores(indel_aux_t *iaux, int ismpl)
{
    int n = iaux->nplp[ismpl] * iaux->ntypes;
    if ( iaux->mread_scores < n )
    {
        int *tmp = (int*) realloc(iaux->read_scores,n*sizeof(int));
        if ( !tmp ) return -1;
        iaux->mread_scores = n;
        iaux->read_scores  = tmp;
    }
    memset(iaux->read_scores,0,n);
    return 0;
}

static int _have_indel_reads(indel_aux_t *iaux)
{
    int i,j;
    for (i=0; i<iaux->nsmpl; i++)
    {
        for (j=0; j<iaux->nplp[i]; j++)
            if ( iaux->plp[i][j].indel ) return 1;
    }
    return 0;
}

// For insertions only their sizes were collected so far. Now go through the reads and
// create consensus sequence for each insert, therefore note that there can be only one
// sequence per insertion length
static int iaux_init_ins_types(indel_aux_t *iaux)
{
    if ( !iaux->max_ins_len ) return 0;

    uint32_t *aux;
    int naux = 5 * iaux->ntypes * iaux->max_ins_len;
    if ( iaux->muitmp < naux )
    {
        aux = (uint32_t*) realloc(iaux->uitmp,naux*sizeof(*aux));
        if ( !aux ) return -1;
        iaux->uitmp  = aux;
        iaux->muitmp = naux;
    }
    else aux = iaux->uitmp;
    memset(aux,0,naux*sizeof(*aux));

    // count the number of occurrences of each base at each position for each type of insertion
    int t,s,i,j;
    for (t=0; t<iaux->ntypes; t++)
    {
        if ( iaux->types[t] <= 0) continue;
        for (s=0; s<iaux->nsmpl; s++)
        {
            for (i=0; i<iaux->nplp[s]; i++)
            {
                bam_pileup1_t *plp = iaux->plp[s] + i;
                if ( plp->indel != iaux->types[t] ) continue;
                uint8_t *seq = bam_get_seq(plp->b);
                for (j=0; j<plp->indel; j++)
                {
                    int c = seq_nt16_int[bam_seqi(seq, plp->qpos+j+1)];
                    assert(c<5);
                    aux[5*(t*iaux->max_ins_len+j) + c]++;
                }
            }
        }
    }

    char *cns;
    int ncns = iaux->ntypes * iaux->max_ins_len;
    if ( iaux->minscns < ncns )
    {
        cns = (char*) realloc(iaux->inscns,naux*sizeof(*aux));
        if ( !cns ) return -1;
        iaux->inscns  = cns;
        iaux->minscns = ncns;
    }
    else cns = iaux->inscns;
    memset(aux,0,ncns*sizeof(*cns));

    // use the majority rule to construct the consensus
    for (t=0; t<iaux->ntypes; t++)
    {
        for (i=0; i<iaux->types[t]; i++)    // this naturally includes only insertions
        {
            uint32_t *tmp = &aux[5*(t*iaux->max_ins_len+i)], max = tmp[0], max_j = 0;
            for (j=1; j<5; j++)
                if ( max < tmp[j] ) max = tmp[j], max_j = j;
            cns[t*iaux->max_ins_len + i] = max ? max_j : 4;
            if ( max_j==4 ) { iaux->types[t] = 0; break; } // discard insertions which contain N's
        }
    }
    return 0;
}

#define MINUS_CONST 0x10000000
static int iaux_init_types(indel_aux_t *iaux)
{
    if ( !_have_indel_reads(iaux) ) return 0;

    iaux->bca->max_support = 0;
    memset(iaux->sum_qual,0,MAX_TYPES*sizeof(*iaux->sum_qual));

    int i,j, nreads = 0;
    for (i=0; i<iaux->nsmpl; i++) nreads += iaux->nplp[i];

    uint32_t *aux;
    if ( iaux->muitmp < nreads+1 )
    {
        aux = (uint32_t*) realloc(iaux->uitmp,(nreads+1)*sizeof(*iaux->uitmp));
        if ( !aux ) return -1;
        iaux->uitmp  = aux;
        iaux->muitmp = nreads+1;
    }
    else aux = iaux->uitmp;
    memset(aux,0,(nreads+1)*sizeof(*aux));

    int naux = 0, indel_support_ok = 0, n_alt = 0, n_tot = 0;
    int max_rd_len = 0;   // max sequence length that includes ref+del bases

    // Fill out aux[] array with all the non-zero indel sizes. This is an unsorted list with as many
    // entries as there are reads
    aux[naux++] = MINUS_CONST;  // zero indel is always a type (REF)
    for (i=0; i<iaux->nsmpl; i++)
    {
        int nalt = naux, ntot = 0;  // per sample values
        for (j=0; j<iaux->nplp[i]; j++)
        {
            const bam_pileup1_t *plp = iaux->plp[i] + j;
            ntot++;
            if ( plp->indel ) aux[naux++] = MINUS_CONST + plp->indel;
            if ( !PLP_QLEN(&plp->cd) ) PLP_QLEN(&plp->cd) = bam_cigar2qlen(plp->b->core.n_cigar, bam_get_cigar(plp->b));
            if ( PLP_QLEN(&plp->cd) > max_rd_len ) max_rd_len = PLP_QLEN(&plp->cd);
        }
        nalt = naux - nalt;
        if ( iaux->bca->per_sample_flt )
        {
            double frac = (double)nalt/naux;
            if ( nalt >= iaux->bca->min_support && frac >= iaux->bca->min_frac ) indel_support_ok = 1;
            if ( nalt > iaux->bca->max_support && frac > 0 ) iaux->bca->max_support = nalt, iaux->bca->max_frac = frac;
        }
        else
        {
            n_alt += nalt;
            n_tot += ntot;
        }
    }

    // Check if the minimum required number of indel reads has been observed
    if ( !iaux->bca->per_sample_flt && n_alt >= iaux->bca->min_support && (double)n_alt/n_tot >= iaux->bca->min_frac ) indel_support_ok = 1;
    if ( naux==1 || !indel_support_ok ) return 0;

    // To prevent long stretches of N's to be mistaken for indels (sometimes thousands of bases), check the number of N's in the
    // sequence and skip places where half or more reference bases in the sequence that follows pos are Ns
    int nN = 0, i_end = iaux->pos + (iaux->bca->indel_win_size < max_rd_len ? iaux->bca->indel_win_size : max_rd_len);
    for (i=iaux->pos; i<i_end && iaux->ref[i]; i++)
        if ( iaux->ref[i] == 'N' ) nN++;
    if ( 2*nN > i - iaux->pos ) return -1;

    // Sort aux[] and dedup indel types
    int n_types = 1;
    ks_introsort(uint32_t, naux, aux);
    for (i=1; i<naux; i++)
        if ( aux[i] != aux[i-1] ) n_types++;

    if ( n_types >= MAX_TYPES )
    {
        static int warned = 0;
        if ( !warned )
        {
            fprintf(bcftools_stderr, "Warning: excessive number of INDEL alleles at %s:%d, skipping. (This warning is printed only once)\n",iaux->chr,iaux->pos+1);
            warned = 1;
        }
        return -1;
    }

    // Fill out the types[] array detailing the size of insertion or deletion.
    iaux->ntypes = 0;
    iaux->max_ins_len = 0;
    for (i=0; i<naux; i++)
    {
        int isize = (int32_t)(aux[i] - MINUS_CONST);
        for (j=i+1; j<naux; j++)
            if ( aux[j] != aux[i] ) break;

        // Only include the REF type and types with sufficient support. Note that the position
        // already passed, this is just to reduce the number of indel types. The check is
        // permissive, the thresholds min_support and min_frac are not enforced in per-sample mode
        int is_ok = 0;
        if ( !isize )
        {
            is_ok = 1;
            iaux->iref_type = iaux->ntypes;
        }
        else
        {
            if ( j-i >= iaux->bca->min_support ) is_ok = 1;
            // What is the best way to handle the -pmF options:
            //  - consider only sites where a single indel type passes the -mF threshold, as opposed to all indel types cumulatively
            //  - once a site passes, include all indel types in the evaluation, as opposed to considering only the strong candidates
            // In this implementation sites are selected by counting reads from all indel types cumulatively and all indel types
            // are considered.
            // Uncomment the following condition to consider only strong indel candidates once the site has been selected
            //      if ( !iaux->bca->per_sample_flt && (double)(j-i) / n_tot < iaux->bca->min_frac ) is_ok = 0;
        }
        if ( is_ok )
        {
            iaux->types[iaux->ntypes++] = isize;
            if ( isize > 0 && isize > iaux->max_ins_len ) iaux->max_ins_len = isize;
        }
        i = j-1;
    }
    if ( iaux->ntypes <= 1 ) return 0;

    // Init insertion types, including their sequence
    if ( iaux_init_ins_types(iaux) < 0 ) return -1;

    iaux_init_sequence_context(iaux);

    return iaux->ntypes;
}
#undef MINUS_CONST

static int iaux_set_consensus(indel_aux_t *iaux, int ismpl)
{
    if ( !iaux->rcns )
        iaux->rcns = rcns_init(iaux->pos, iaux->left, iaux->right);
    else
        rcns_reset(iaux->rcns, iaux->pos, iaux->left, iaux->right);

    rcns_set_reads(iaux->rcns, iaux->plp[ismpl], iaux->nplp[ismpl]);

    iaux->cns_seq = rcns_get_consensus(iaux->rcns, iaux->ref + iaux->left);

// todo:
//  rcns should also collect localized number of mismatches as a substitute
//  for uninformative MQ. This would not affect calling but would help with
//  filtering

    return 0;
}

#if 0
// Finds the smallest index in the seq_pos array holding value equal to pos, or if there is no
// such value, the largest index with value smaller than pos. Starts at initial guess ioff.
// This could use a binary search but the assumption is that the initial guess is indel-size close
// to the actuall coordinate.
//
// TODO: remove this function and seq_pos from cns creation as it seems unnecessary
static int find_ref_offset(hts_pos_t pos, hts_pos_t *seq_pos, int nseq_pos, int ioff)
{
    if ( ioff<0 ) ioff = 0;
    else if ( ioff >= nseq_pos ) ioff = nseq_pos - 1;
    if ( seq_pos[ioff] < pos )
    {
        while ( ioff+1 < nseq_pos && seq_pos[ioff] < pos ) ioff++;
        if ( seq_pos[ioff] > pos ) ioff--;
        return ioff;
    }
    while ( ioff > 0 && seq_pos[ioff-1] >= pos ) ioff--;
    return ioff;
}
#endif

static int iaux_align_read(indel_aux_t *iaux, bam1_t *bam, uint8_t *ref_seq, int nref_seq)
{
    if ( bam->core.flag & BAM_FUNMAP ) return 1;   // skip unmapped reads

    // Trim both ref and qry to the window of interest
    hts_pos_t ref_beg = iaux->left;     // fa ref coordinates
    hts_pos_t ref_end = iaux->right < ref_beg + nref_seq ? iaux->right : ref_beg + nref_seq - 1;

    cigar_state_t cigar;
    cstate_init(&cigar,bam);
    int qry_off1, qry_off2, ref_off1, ref_off2;
    if ( ref_beg > bam->core.pos )
    {
        // the read needs trimming from left
        qry_off1 = cstate_seek_fwd(&cigar, &ref_beg, 1);
        ref_off1 = ref_beg - iaux->left;

        if ( ref_beg + (bam->core.l_qseq - qry_off1) > ref_end )
        {
            // the read needs trimming from right
            qry_off2 = ref_end - ref_beg + qry_off1;
            ref_off2 = ref_end - iaux->left;
        }
        else
        {
            // the ref template needs trimming from right
            qry_off2 = bam->core.l_qseq - 1;
            ref_off2 = ref_off1 + qry_off2 - qry_off1;
        }
    }
    else
    {
        // the ref template needs trimming from left
        qry_off1 = 0;
        ref_off1 = bam->core.pos - ref_beg;

        if ( bam->core.pos + bam->core.l_qseq - 1 > ref_end )
        {
            // the read needs trimming from right
            ref_off2 = ref_end - iaux->left;
            qry_off2 = ref_off2 - ref_off1;
        }
        else
        {
            // the ref template needs trimming from right
            qry_off2 = bam->core.l_qseq - 1;
            ref_off2 = ref_off1 + qry_off2 - qry_off1;
        }
    }
//fprintf(bcftools_stderr,"xtrim: %s ..  left,right=%d,%d  rbeg,end=%d,%d  qpos=%d  qlen=%d  qoff=%d,%d  roff=%d,%d rlen=%d\n",bam_get_qname(bam),iaux->left,iaux->right,(int)ref_beg,(int)ref_end,(int)bam->core.pos,bam->core.l_qseq, qry_off1,qry_off2,ref_off1,ref_off2,nref_seq);

    assert( qry_off1<=qry_off2 );
    assert( qry_off1>=0 && qry_off1<bam->core.l_qseq );
    assert( qry_off2>=0 && qry_off2<bam->core.l_qseq );

    assert( ref_off1<=ref_off2 );
    assert( ref_off1>=0 && ref_off1<nref_seq );
    assert( ref_off2>=0 && ref_off2<nref_seq );

    // prepare query sequence
    int i, qlen = qry_off2 - qry_off1 + 1, rlen = ref_off2 - ref_off1 + 1;
    if ( iaux->nqry_seq < qlen )
    {
        uint8_t *tmp = (uint8_t*) realloc(iaux->qry_seq, qlen);
        if ( !tmp ) return -1;  // critical error
        iaux->qry_seq  = tmp;
        iaux->nqry_seq = qlen;
    }
    uint8_t *seq = bam_get_seq(bam);
    for (i=qry_off1; i<=qry_off2; i++) iaux->qry_seq[i-qry_off1] = seq_nt16_int[bam_seqi(seq,i)];

    // prepare qualities, either BQ or BAQ qualities (ZQ)
    if ( iaux->nqual < qlen )
    {
        uint8_t *tmp = (uint8_t*) realloc(iaux->qual, qlen);
        if ( !tmp ) return -1;  // critical error
        iaux->qual  = tmp;
        iaux->nqual = qlen;
    }
    uint8_t *qual = iaux->qual;
    const uint8_t *qq = bam_get_qual(bam);
    const uint8_t *bq = (uint8_t*)bam_aux_get(bam, "ZQ");
    if ( bq ) bq++; // skip type
    for (i=qry_off1; i<=qry_off2; i++)
    {
        int j = i - qry_off1;
        qual[j] = bq ? qq[i] + (bq[i] - 64) : qq[i];
        if ( qual[j] > 30 ) qual[j] = 30;
        if ( qual[j] < 7 ) qual[j] = 7;
    }

// Illumina
probaln_par_t apf = { 1e-4, 1e-2, 10 };

    // align
    int score = probaln_glocal(ref_seq + ref_off1, rlen, iaux->qry_seq, qlen, qual, &apf, 0, 0);
    int adj_score = (int)(100. * score / qlen + .499) * iaux->bca->indel_bias;

#if DEBUG_ALN
    fprintf(bcftools_stderr,"aln: %d/%d\t%s\n\tref:  ",score,adj_score,bam_get_qname(bam));
    for (i=0; i<rlen; i++) fprintf(bcftools_stderr,"%c","ACGTN"[(int)ref_seq[ref_off1 + i]]);
    fprintf(bcftools_stderr,"\n\tqry:  ");
    for (i=0; i<qlen; i++) fprintf(bcftools_stderr,"%c","ACGTN"[(int)iaux->qry_seq[i]]);
    fprintf(bcftools_stderr,"\n\tqual: ");
    for (i=0; i<qlen; i++) fprintf(bcftools_stderr,"%c",(char)(qual[i]+64));
    fprintf(bcftools_stderr,"\n\ttrim: qry_len=%d qry_off=%d,%d   ref_len=%d ref_off=%d,%d ref_beg,end=%d,%d\n",qlen,qry_off1,qry_off2,rlen,ref_off1,ref_off2,(int)ref_beg,(int)ref_end);
#endif

    if ( adj_score > 255 ) adj_score = 255;
    return score<<8 | adj_score;
}

// Score all reads for this sample and indel type using the up to two consensus sequence templates.
// On output sets iaux->read_scores[iread*ntypes+itype] = (raw_score<<8 | length_adjusted_score)
static int iaux_score_reads(indel_aux_t *iaux, int ismpl, int itype)
{
    int i;
    cns_seq_t *cns = iaux->cns_seq;
    while ( cns->nseq )
    {
        // Resize buffers if necessary
        int ref_len = cns->nseq + iaux->types[itype];
        if ( iaux->nref_seq < ref_len )
        {
            uint8_t *ref_buf = (uint8_t*) realloc(iaux->ref_seq,sizeof(uint8_t)*ref_len);
            if ( !ref_buf ) return -1;
            iaux->ref_seq  = ref_buf;
            iaux->nref_seq = ref_len;
        }

        // Apply the indel and create the template ref sequence...
        memcpy(iaux->ref_seq,cns->seq,(cns->ipos+1)*sizeof(*iaux->ref_seq));
        if ( iaux->types[itype] < 0 )   // deletion
            memcpy(iaux->ref_seq + cns->ipos + 1, cns->seq + cns->ipos + 1 - iaux->types[itype], (cns->nseq - cns->ipos - 1 + iaux->types[itype])*sizeof(*iaux->ref_seq));
        else
        {
            char *ins = &iaux->inscns[itype*iaux->max_ins_len];
            for (i=0; i<iaux->types[itype]; i++) iaux->ref_seq[cns->ipos+1+i] = ins[i];
            memcpy(iaux->ref_seq + cns->ipos + 1 + iaux->types[itype], cns->seq + 1 + cns->ipos, (cns->nseq - cns->ipos - 1)*sizeof(*iaux->ref_seq));
        }

#if DEBUG_ALN
    fprintf(bcftools_stderr,"template %d, type %d, sample %d: ",cns==iaux->cns_seq?0:1,itype,ismpl);
    for (i=0; i<ref_len; i++) fprintf(bcftools_stderr,"%c","ACGTN"[(int)iaux->ref_seq[i]]);
    fprintf(bcftools_stderr,"\n");
#endif

        // Align and score reads
        for (i=0; i<iaux->nplp[ismpl]; i++)
        {
            const bam_pileup1_t *plp = iaux->plp[ismpl] + i;
            int aln_score = iaux_align_read(iaux, plp->b, iaux->ref_seq, ref_len);
            int *score = &iaux->read_scores[i*iaux->ntypes+itype];
            if ( cns==iaux->cns_seq || *score > aln_score ) *score = aln_score;
        }
        cns++;
    }
    return 0;
}

// Determines indel quality for each read and populates 22 bits of pileup aux field with
// three integers as follows
//      plp->aux = indel_type << 16 | seqQ << 8 | indelQ
static int iaux_eval_scored_reads(indel_aux_t *iaux, int ismpl)
{
    int i,j;
    for (i=0; i<iaux->nplp[ismpl]; i++)
    {
        bam_pileup1_t *plp = iaux->plp[ismpl] + i;

        // Find the best indel type and the ref type, their scores difference is the indel quality
        int *score = &iaux->read_scores[i*iaux->ntypes];
        int alt_score = INT_MAX, alt_j = 0;
        for (j=0; j<iaux->iref_type; j++)
            if ( alt_score > score[j] ) alt_score = score[j], alt_j = j;
        for (j=iaux->iref_type+1; j<iaux->ntypes; j++)
            if ( alt_score > score[j] ) alt_score = score[j], alt_j = j;
        int ref_score = score[iaux->iref_type];
        int sc0, sc1, j0;
        if ( alt_score < ref_score ) sc0 = alt_score, sc1 = ref_score, j0 = alt_j;
        else sc0 = ref_score, sc1 = alt_score, j0 = iaux->iref_type;

        int indelQ = (sc1>>8) - (sc0>>8);           // low=bad, high=good
        int seqQ   = iaux->ref_qual[alt_j];

        // Reduce indelQ. High length-normalized alignment scores (i.e. bad alignments)
        // lower the quality more (e.g. gnuplot> plot [0:111] (1-x/111.)*255)
        int len_normQ  = sc0 & 0xff;    // length-normalized score of the best match (ref or alt)
        int adj_indelQ;                 // final indelQ used in calling
        if ( len_normQ > 111 )
        {
            // In the original code reads matching badly to any indel type or reference had indelQ set to 0
            // here and thus would be effectively removed from calling. This leads to problems when there are
            // many soft clipped reads and a few good matching indel reads (see noisy-softclips.bam in
            // mpileup-tests). Only the few good quality indel reads would become visible to the caller and
            // the indel would be called with high quality. Here we change the logic to make the badly matching
            // reads low quality reference reads. The threshold was set to make the test case still be called
            // as an indel, but with very low quality.
            //
            // Original code:
            //  adj_indelQ = 0;
            //
            adj_indelQ = 12;
            j0 = iaux->iref_type;
        }
        else
            adj_indelQ = (int)((1. - len_normQ/111.) * indelQ + .499);

#if DEBUG_ALN
        // Prints the selected indel type (itype); adjusted indelQ which will be used if bigger than seqQ;
        //  raw indelQ; length-normalized indelQ and sequence context quality; ref and best alt indel type
        //  and their raw and length-normalized scores
        fprintf(bcftools_stderr,"itype=%d adj_indelQ=%d\trawQ=%d\tlen_normQ=%d\tseqQ=%d\tref:%d=%d/%d alt:%d=%d/%d)\t%s\n",
            j0,adj_indelQ,indelQ,len_normQ,seqQ,iaux->iref_type,ref_score>>8,ref_score&0xff,alt_j,alt_score>>8,alt_score&0xff,bam_get_qname(plp->b));
#endif

        if ( adj_indelQ > seqQ ) adj_indelQ = seqQ;     // seqQ already capped at 255
        plp->aux = j0<<16 | seqQ<<8 | adj_indelQ;       // use 22 bits in total
        iaux->sum_qual[j0] += adj_indelQ;
    }
    return 0;
}

// Find the best indel types, include the ref type plus maximum three alternate indel alleles.
static int iaux_eval_best_indels(indel_aux_t *iaux)
{
    bcf_callaux_t *bca = iaux->bca;
    bca->maxins = iaux->max_ins_len;
    bca->inscns = (char*) realloc(bca->inscns, bca->maxins * 4);
    if ( bca->maxins && !bca->inscns ) return -1;

    // insertion sort, descending, high-quality indels come first
    int i,j,t, tmp, *sumq = iaux->sum_qual, ntypes = iaux->ntypes;
    for (t=0; t<ntypes; t++) sumq[t] = sumq[t]<<6 | t;
    for (t=1; t<ntypes; t++)
        for (j=t; j>0 && sumq[j] > sumq[j-1]; j--)
            tmp = sumq[j], sumq[j] = sumq[j-1], sumq[j-1] = tmp;
    for (t=0; t<ntypes; t++)  // look for the reference type
        if ( (sumq[t]&0x3f)==iaux->iref_type ) break;
    if ( t )
    {
        // move the reference type to the first
        tmp = sumq[t];
        for (; t>0; t--) sumq[t] = sumq[t-1];
        sumq[0] = tmp;
    }

    // Initialize bca's structures and create a mapping between old and new types
    int old2new_type[MAX_TYPES];
    for (t=0; t<iaux->ntypes; t++)
    {
        int itype = sumq[t] & 0x3f;
        old2new_type[itype] = t;
        if ( t>=4 ) continue;
        bca->indel_types[t] = iaux->types[itype];
        if ( bca->indel_types[t] <= 0 ) continue;
        memcpy(&bca->inscns[t*bca->maxins], &iaux->inscns[itype*iaux->max_ins_len], bca->maxins);
    }

    // Update indel type in plp->aux for all reads
    int ismpl, n_alt = 0;
    for (ismpl=0; ismpl<iaux->nsmpl; ismpl++)
    {
        for (i=0; i<iaux->nplp[ismpl]; i++)
        {
            bam_pileup1_t *plp = iaux->plp[ismpl] + i;
            int itype_old = (plp->aux >> 16) & 0x3f;
            int itype_new = old2new_type[itype_old];
            plp->aux = itype_new<<16 | (itype_new>=4 ? 0 : (plp->aux & 0xffff));
            if ( itype_new>0 ) n_alt++;
        }
    }
    return n_alt;
}

/*
    notes:
    - n .. number of samples
    - the routine sets bam_pileup1_t.aux (27 bits) of each read as follows:
        - 5: unused
        - 6: the call; index to bcf_callaux_t.indel_types   .. (aux>>16)&0x3f
        - 8: estimated sequence quality                     .. (aux>>8)&0xff
        - 8: indel quality                                  .. aux&0xff
 */
int bcf_iaux_gap_prep(int n, int *n_plp, bam_pileup1_t **plp, int pos, bcf_callaux_t *bca, const char *ref)
{
assert(!(ref == 0 || bca == 0));    // can this ever happen? when?
    if (ref == 0 || bca == 0) return -1;

    if ( !bca->iaux ) bca->iaux = calloc(1,sizeof(indel_aux_t));
    indel_aux_t *iaux = bca->iaux;
    iaux->nsmpl = n;
    iaux->nplp  = n_plp;
    iaux->plp   = plp;
    iaux->bca   = bca;
    iaux->ref   = ref;
    iaux->pos   = pos;
    iaux->chr   = bca->chr;

    // Check if there is an indel at this position and if yes, find all indel types and determine
    // window boundaries. todo: We want this information cached so that for long reads we don't keep
    // redoing the whole analysis again and again
    int ntypes = iaux_init_types(iaux);
    if ( ntypes<=0 ) return -1;

    debug_print_types(iaux);

    // Create two template consensus sequences for each sample (assuming max diploid organism).
    // Then apply each indel type on top of the templates, realign every read and remember score
    int i,j;
    for (i=0; i<iaux->nsmpl; i++)
    {
        iaux_set_consensus(iaux, i);
        iaux_init_scores(iaux, i);
        for (j=0; j<ntypes; j++) iaux_score_reads(iaux, i, j);
        iaux_eval_scored_reads(iaux, i);
    }
    int nalt = iaux_eval_best_indels(iaux);
    return nalt > 0 ? 0 : -1;
}

