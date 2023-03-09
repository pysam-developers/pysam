#include "bcftools.pysam.h"

/*  read_consensus.c -- create and maintain consensus of reads

    Copyright (C) 2022 Genome Research Ltd.

    Author: pd3@sanger

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

#include <assert.h>
#include <math.h>
#include "bcftools.h"
#include "read_consensus.h"
#include "cigar_state.h"
#include "kheap.h"


// Frequency arrays for each variant type
#define NI 10 // number of alternative insertion sequences at one position in a single sample
typedef struct
{
    char *nt16_seq[NI];
    int len[NI];
    int freq[NI];
}
ins_freq_t;

typedef struct
{
    int len[NI];
    int freq[NI];
}
del_freq_t;

#define BF_DEL 5
typedef struct
{
    int base[6];    // frequencies of A,C,G,T,N,deletion
}
base_freq_t;


// Candidate variants for each interesting position to build consensus haplotypes
enum variant_type { snv, ins, del, done };
typedef struct
{
    enum variant_type vtype;
    hts_pos_t pos;      // variant position (reference sequence coordinates), indels follow VCF convention
    int idx;            // temporary 0-based index to rcns.cvar
    int which,          // base/ins/del in rcns.[base|ins|del]_freq array
        depth;          // coverage at the position
    float af, af_dev;   // variant allele frequency (just for debugging printout) and absolute af deviation from 0.5
}
candidate_var_t;
static inline int cvar_not_preferred(candidate_var_t *a, candidate_var_t *b)
{
    if ( a->af_dev == b->af_dev ) return a->depth < b->depth ? 1 : 0;
    return a->af_dev > b->af_dev ? 1 : 0;
}
KHEAP_INIT(cvh, candidate_var_t, cvar_not_preferred);
typedef khp_cvh_t cvar_heap_t;

#define MAX_NCVAR 8         // This results in alloc() of 2^MAX_NCVAR possible haplotypes
#define NHAP (1<<MAX_NCVAR) // The number of possible haplotypes
struct _read_cns_t
{
    hts_pos_t pos, beg, end;    // current position and window boundaries (0-based, inclusive, ref seq coordinates)
    int band,                   // maximum absolute deviation from the diagonal, used for BAQ alignment
        max_del;                // maximum deletion lentgth starting at the tested position
    base_freq_t *base_freq;     // frequency of each variant type: base, ins, del
    ins_freq_t *ins_freq;
    del_freq_t *del_freq;
    char *stmp;             // temporary array
    int mstmp, mfreq;       // allocated size of stmp and *_freq arrays
    cvar_heap_t *cv_heap;   // heap to maintain the top MAX_NCVAR variants
    int ncvar;              // cvar and cv_heap size
    candidate_var_t cvar[MAX_NCVAR];    // candidate variants, sorted by position and type
    int hap_freq[NHAP];     // haplotype frequencies
    bam_pileup1_t *plp;     // reads to construct consensus from
    int nplp;               // number of reads in the pileup
    int cns_hap[2], ncns;   // the top two consensus haplotypes and the number of haplotypes to use
    int mcns;               // the allocated size of cns.seq and cns.pos buffers
    cns_seq_t cns[3];       // the consensus sequences to fill
};

void rcns_destroy(read_cns_t *rcns)
{
    if ( !rcns ) return;
    int i,j;
    for (i=0; i<rcns->mfreq; i++)
    {
        ins_freq_t *ifrq = &rcns->ins_freq[i];
        for (j=0; j<NI && ifrq->nt16_seq[j]; j++) free(ifrq->nt16_seq[j]);
    }
    for (i=0; i<2; i++)
        free(rcns->cns[i].seq);
    free(rcns->ins_freq);
    free(rcns->del_freq);
    free(rcns->base_freq);
    free(rcns->stmp);
    khp_destroy(cvh,rcns->cv_heap);
    free(rcns);
}
static int init_arrays(read_cns_t *rcns)
{
    int i,j,n = rcns->end - rcns->beg + 1;
    if ( n > rcns->mfreq )
    {
        ins_freq_t *ifrq = (ins_freq_t*) realloc(rcns->ins_freq,sizeof(*rcns->ins_freq)*n);
        if ( !ifrq ) return -1;
        rcns->ins_freq = ifrq;
        memset(ifrq+rcns->mfreq,0,sizeof(*rcns->ins_freq)*(n-rcns->mfreq));

        del_freq_t *dfrq = (del_freq_t*) realloc(rcns->del_freq,sizeof(*rcns->del_freq)*n);
        if ( !dfrq ) return -1;
        rcns->del_freq = dfrq;
        memset(dfrq+rcns->mfreq,0,sizeof(*rcns->del_freq)*(n-rcns->mfreq));

        base_freq_t *bfrq = (base_freq_t*) realloc(rcns->base_freq,sizeof(*rcns->base_freq)*n);
        if ( !bfrq ) return -1;
        rcns->base_freq = bfrq;
        memset(bfrq+rcns->mfreq,0,sizeof(*rcns->base_freq)*(n-rcns->mfreq));

        rcns->mfreq = n;
    }
    memset(rcns->base_freq,0,sizeof(*rcns->base_freq)*n);
    memset(rcns->del_freq,0,sizeof(*rcns->del_freq)*n);
    for (i=0; i<n; i++)
    {
        ins_freq_t *ifrq = &rcns->ins_freq[i];
        for (j=0; j<NI && ifrq->nt16_seq[j]; j++) free(ifrq->nt16_seq[j]);
    }
    memset(rcns->ins_freq,0,sizeof(*rcns->ins_freq)*n);
    return 0;
}
int rcns_reset(read_cns_t *rcns, hts_pos_t pos, hts_pos_t beg, hts_pos_t end)
{
    rcns->band = 0;
    rcns->pos  = pos;
    rcns->beg  = beg;
    rcns->end  = end;
    int i;
    for (i=0; i<2; i++) rcns->cns[i].nseq = rcns->cns[i].ipos = 0;
    // this should not be necessary if the caller did run all steps
    while (rcns->cv_heap->ndat) khp_delete(cvh, rcns->cv_heap);
    return init_arrays(rcns);
}

static inline void add_base(read_cns_t *rcns, int ref_pos, int nt16)
{
    int i = ref_pos - rcns->beg;
    rcns->base_freq[i].base[seq_nt16_int[nt16]]++;
}
static void add_ins(read_cns_t *rcns, int ref_pos, int seq_pos, uint8_t *raw_seq, int len)
{
    int i = ref_pos - rcns->beg;
    ins_freq_t *ifrq = &rcns->ins_freq[i];
    char *str;
    if ( rcns->mstmp < len )
    {
        str = realloc(rcns->stmp,len*sizeof(*str));
        if ( !str ) return;
        rcns->mstmp = len;
        rcns->stmp  = str;
    }
    else
        str = rcns->stmp;
    for (i=0; i<len; i++) str[i] = bam_seqi(raw_seq,i+seq_pos);

    for (i=0; i<NI && ifrq->nt16_seq[i]; i++)
        if ( ifrq->len[i]==len && !memcmp(ifrq->nt16_seq[i],str,len) ) break;

    if ( i>=NI ) return;    // too many choices, typically homopolymers in long reads; discard

    if ( !ifrq->nt16_seq[i] )      // new insertion
    {
        if ( !(ifrq->nt16_seq[i]=malloc(len)) ) return;
        memcpy(ifrq->nt16_seq[i], str, len);
        ifrq->len[i] = len;
    }
    ifrq->freq[i]++;
}
static void add_del(read_cns_t *rcns, int ref_pos, int len)
{
    int i = ref_pos - rcns->beg;
    int j,n = rcns->end - rcns->beg + 1;
    if ( i + len + 1 < n ) n = i + len + 1;
    for (j=i+1; j<n; j++)
        rcns->base_freq[j].base[BF_DEL]++;

    del_freq_t *dfrq = &rcns->del_freq[i];
    for (i=0; i<NI && dfrq->len[i]; i++)
        if ( dfrq->len[i]==len ) break;

    if ( i>=NI ) return;    // too many choices, typically homopolymers in long reads; discard

    if ( !dfrq->len[i] ) dfrq->len[i] = len;    // new deletion
    dfrq->freq[i]++;
}

read_cns_t *rcns_init(hts_pos_t pos, hts_pos_t beg, hts_pos_t end)
{
    read_cns_t *rcns = (read_cns_t*) calloc(1,sizeof(read_cns_t));
    rcns->pos  = pos;
    rcns->beg  = beg;
    rcns->end  = end;
    rcns->cv_heap = khp_init(cvh);
    if ( init_arrays(rcns)!=0 )
    {
        rcns_destroy(rcns);
        return NULL;
    }
    return rcns;
}

int rcns_set_reads(read_cns_t *rcns, bam_pileup1_t *plp, int nplp)
{
    // save the reads for phasing, this can be called multiple times
    rcns->plp  = plp;
    rcns->nplp = nplp;

    // fill consensus arrays
    int i,j,k, local_band_max = 0;  // maximum absolute deviation from diagonal
    for (i=0; i<nplp; i++) // for each read...
    {
        const bam_pileup1_t *p = plp + i;
        bam1_t *b = p->b;
        int x = b->core.pos;  // ref coordinate
        int y = 0;            // seq coordinate
        uint32_t *cigar = bam_get_cigar(b);
        uint8_t *seq = bam_get_seq(b);
        int local_band = 0; // current deviation from diagonal
        for (k = 0; k < b->core.n_cigar; ++k)
        {
            int op  = cigar[k] &  BAM_CIGAR_MASK;
            int len = cigar[k] >> BAM_CIGAR_SHIFT;
            if ( op==BAM_CSOFT_CLIP ) y += len;
            else if ( op==BAM_CMATCH || op==BAM_CEQUAL || op==BAM_CDIFF )
            {
                if ( x<rcns->end && x+len>rcns->beg )
                {
                    int j_beg = rcns->beg > x ? rcns->beg - x : 0;  // how many bases to skip in the ref and qry
                    int j_end = rcns->end < x + len - 1 ? rcns->end - x : len - 1;
                    x += j_beg;
                    y += j_beg;
                    for (j=j_beg; j<=j_end; j++, x++, y++) add_base(rcns,x,bam_seqi(seq,y));
                }
                else
                {
                    x += len;
                    y += len;
                }
            }
            else if ( op==BAM_CINS )
            {
                if ( x>rcns->beg && x<rcns->end )
                {
                    local_band += p->indel;
                    add_ins(rcns,x-1,y,seq,len);    // x-1: one base before as in VCF
                }
                y += len;
            }
            else if ( op==BAM_CDEL )
            {
                if ( x>rcns->beg && x+len-1<=rcns->end )
                {
                    local_band += -p->indel;
                    add_del(rcns,x-1,len);          // x-1: one base before as in VCF
                }
                x += len;
            }
            else if ( op==BAM_CHARD_CLIP ) continue;
            else error("rcns_set_reads todo: unknown cigar operator %d\n",op);
            if ( local_band_max < local_band ) local_band_max = local_band;
        }

        // Track the biggest deviation +/- from diagonal, used in BAQ alignment step.
        if ( rcns->band < local_band_max ) rcns->band = local_band_max;
    }

    return 0;
}

#if DEBUG_RCNS
static void debug_print_base_freqs(read_cns_t *rcns, const char *ref)
{
    int i,j,k,n = rcns->end - rcns->beg + 1;
    fprintf(bcftools_stderr,"beg,end,pos=%d %d %d\n",(int)rcns->beg,(int)rcns->end,(int)rcns->pos);
    base_freq_t *bfreq = rcns->base_freq;
    ins_freq_t *ifreq  = rcns->ins_freq;
    del_freq_t *dfreq  = rcns->del_freq;
    for (i=0; i<n && ref[i]; i++)
    {
        fprintf(bcftools_stderr,"%"PRIhts_pos" %c\t",rcns->beg+i+1,ref[i]);
        for (j=0; j<6; j++)
            fprintf(bcftools_stderr,"\t%d%s",bfreq[i].base[j],ref[i]=="ACGTNi"[j]?"*":"");
        fprintf(bcftools_stderr,"\t");
        for (j=0; j<NI && dfreq[i].len[j]; j++)
            fprintf(bcftools_stderr," -%d:%d",dfreq[i].len[j],dfreq[i].freq[j]);
        fprintf(bcftools_stderr,"\t");
        for (j=0; j<NI && ifreq[i].len[j]; j++)
        {
            fprintf(bcftools_stderr," +");
            for (k=0; k<ifreq[i].len[j]; k++) fprintf(bcftools_stderr,"%c",seq_nt16_str[(int)ifreq[i].nt16_seq[j][k]]);
            fprintf(bcftools_stderr,":%d",ifreq[i].freq[j]);
        }
        fprintf(bcftools_stderr,"\n");
    }
}
static const char *vtype2string(enum variant_type vtype)
{
    if ( vtype==snv ) return "snv";
    if ( vtype==ins ) return "ins";
    if ( vtype==del ) return "del";
    return "???";
}
static void debug_print_candidate_variants(read_cns_t *rcns)
{
    int i;
    fprintf(bcftools_stderr,"Candidate variants:\n");
    for (i=0; i<rcns->ncvar; i++)
    {
        candidate_var_t *var = &rcns->cvar[i];
        fprintf(bcftools_stderr,"\tvar%d  pos=%"PRIhts_pos" idx=%d vtype=%s which=%d depth=%d af=%f af_dev=%f\n",
            i,var->pos+1,var->idx,vtype2string(var->vtype),var->which,var->depth,var->af,var->af_dev);
    }
}
static void debug_print_haplotype_frequency_spectrum(read_cns_t *rcns)
{
    int i,j;
    fprintf(bcftools_stderr,"Haplotype frequencies (bits from left correspond to var0,1,..):\n");
    for (i=0; i<NHAP; i++)
    {
        if ( !rcns->hap_freq[i] ) continue;
        fprintf(bcftools_stderr,"\t%d: ",i);
        for (j=0; j<rcns->ncvar; j++)
            fprintf(bcftools_stderr,"%d", i&(1<<j) ? 1 : 0);
        fprintf(bcftools_stderr,"\t%d\n", rcns->hap_freq[i]);
    }
}
static void debug_print_consensus(read_cns_t *rcns, const char *ref)
{
    int i,j,n = rcns->end - rcns->beg + 1;
    fprintf(bcftools_stderr,"ref:        ");
    for (i=0; i<n && ref[i]; i++) fprintf(bcftools_stderr,"%c",ref[i]);
    fprintf(bcftools_stderr,"\n");
    for (i=0; i<2; i++)
    {
        if ( !rcns->cns[i].nseq ) break;
        fprintf(bcftools_stderr,"Consensus%d: ",i);
        for (j=0; j<=rcns->cns[i].ipos; j++)
            fprintf(bcftools_stderr,"%c","ACGTN"[(int)rcns->cns[i].seq[j]]);
        fprintf(bcftools_stderr,"#");
        for (; j<rcns->cns[i].nseq; j++)
            fprintf(bcftools_stderr,"%c","ACGTN"[(int)rcns->cns[i].seq[j]]);
        fprintf(bcftools_stderr,"\n");
    }
}
#else
#define debug_print_base_freqs(rcns,ref)
#define debug_print_candidate_variants(rcns)
#define debug_print_haplotype_frequency_spectrum(rcns)
#define debug_print_consensus(rcns,ref)
#endif

static int cvar_pos_cmp(const void *aptr, const void *bptr)
{
    candidate_var_t *a = (candidate_var_t*)aptr;
    candidate_var_t *b = (candidate_var_t*)bptr;
    if ( a->pos < b->pos ) return -1;
    if ( a->pos > b->pos ) return 1;
    if ( a->vtype < b->vtype ) return -1;
    if ( a->vtype > b->vtype ) return 1;
    if ( a->which < b->which ) return -1;
    if ( a->which > b->which ) return 1;
    return 0;
}
static void register_variant(read_cns_t *rcns, enum variant_type vtype, int cns_pos, int which, int depth, float freq)
{
    cvar_heap_t *cv_heap = rcns->cv_heap;
    if ( vtype==done )
    {
        rcns->ncvar = 0;
        while (cv_heap->ndat)
        {
            rcns->cvar[rcns->ncvar++] = cv_heap->dat[0];
            khp_delete(cvh,cv_heap);
        }
        // sort the variants by pos,type,which to make determination of haplotypes from reads faster
        if ( rcns->ncvar )
            qsort(rcns->cvar, rcns->ncvar, sizeof(*rcns->cvar), cvar_pos_cmp);
        return;
    }

    candidate_var_t var;
    var.pos    = cns_pos + rcns->beg;
    var.which  = which;
    var.vtype  = vtype;
    var.depth  = depth;
    var.af_dev = fabs(0.5-freq);
    var.af     = freq;

    int free_slot;

    // keep the number of variants small, maximum MAX_NCVAR
    if ( rcns->ncvar==MAX_NCVAR )
    {
        if ( cvar_not_preferred(&var,&cv_heap->dat[0]) ) return;  // no need to add, the new variant is worse than the heap's worst one
        free_slot = cv_heap->dat[0].idx;
        khp_delete(cvh,cv_heap);
    }
    else
        free_slot = rcns->ncvar++;
    var.idx =  free_slot;
    rcns->cvar[free_slot] = var;
    khp_insert(cvh,cv_heap,&var);
}

// Identify candidate variant positions. (Note that homozygous variants are not considered
// as those will be added trivially by taking the consensus base.) The detection limit is
// for now hard-wired. This has only indirect effect on sensitivity, will just not contribute
// to the consensus template when realigning.
static int select_candidate_variants(read_cns_t *rcns, const char *ref)
{
    const float af_th = 0.1;
    int i,j, n = rcns->end - rcns->beg + 1;
    int max_ins_len = 0;    // maximum total length of all insertions applied to allocate big enough buffers
    base_freq_t *bfreq = rcns->base_freq;
    ins_freq_t *ifreq  = rcns->ins_freq;
    del_freq_t *dfreq  = rcns->del_freq;
    for (i=0; i<n && ref[i]; i++)
    {
        for (j=0; j<NI && ifreq[i].len[j]; j++) max_ins_len += ifreq[i].len[j];

        if ( i==rcns->pos - rcns->beg ) continue;   // creating consensus from everything but the variants at the current position

        int dp = 0;
        for (j=0; j<4; j++) dp += bfreq[i].base[j];
        for (j=0; j<NI && dfreq[i].len[j]; j++) dp += dfreq[i].freq[j];
        for (j=0; j<NI && ifreq[i].len[j]; j++) dp += ifreq[i].freq[j];
        float af = 0;   // allele frequency
        for (j=0; j<4; j++)
        {
            if ( !bfreq[i].base[j] || ref[i]=="ACGTN"[j] ) continue;   // ref base or no coverage
            af = (float)bfreq[i].base[j]/dp;
            if ( af>af_th && af<(1-af_th) ) register_variant(rcns,snv,i,j,dp,af);
        }
        for (j=0; j<NI && dfreq[i].len[j]; j++)
        {
            af = (float)dfreq[i].freq[j]/dp;
            if ( af>af_th && af<(1-af_th) ) register_variant(rcns,del,i,j,dp,af);
        }
        for (j=0; j<NI && ifreq[i].len[j]; j++)
        {
            af = (float)ifreq[i].freq[j]/dp;
            if ( af>af_th && af<(1-af_th) ) register_variant(rcns,ins,i,j,dp,af);
        }
    }
    register_variant(rcns,done,0,0,0,0);    // finalize

    // Reallocate buffers
    if ( rcns->mcns < n + max_ins_len )
    {
        n += max_ins_len;
        for (i=0; i<2; i++)
        {
            char *seq = (char*) realloc(rcns->cns[i].seq,sizeof(char)*n);
            if ( !seq ) return -1;
            rcns->cns[i].seq = seq;
        }
        rcns->mcns = n;
    }

    // Find the longest deletion at the query position
    i = rcns->pos - rcns->beg;
    rcns->max_del = 0;
    for (j=0; j<NI && j<dfreq[i].len[j]; j++)
    {
        if ( rcns->max_del < dfreq[i].len[j] ) rcns->max_del = dfreq[i].len[j];
    }

    return 0;
}
static int create_haplotype_frequency_spectrum(read_cns_t *rcns)
{
    memset(rcns->hap_freq,0,sizeof(rcns->hap_freq));

    int i;
    for (i=0; i<rcns->nplp; i++) // for each read...
    {
        const bam_pileup1_t *p = rcns->plp + i;
        cigar_state_t cigar;
        cstate_init(&cigar,p->b);

        int j,k,hap = 0;
        for (j=0; j<rcns->ncvar; j++)
        {
            candidate_var_t *cvar = &rcns->cvar[j];
            if ( cvar->vtype==snv )
            {
                int iseq = cstate_seek_op_fwd(&cigar, cvar->pos, BAM_CMATCH, NULL);
                if ( iseq==-2 ) break;
                if ( iseq==-1 ) continue;
                int nt16 = bam_seqi(cigar.seq, iseq);
                if ( seq_nt16_int[nt16]==cvar->which ) hap |= 1<<j;
            }
            else if ( cvar->vtype==ins )
            {
                int len;
                ins_freq_t *ifrq = &rcns->ins_freq[cvar->pos - rcns->beg];
                int iseq = cstate_seek_op_fwd(&cigar, cvar->pos+1, BAM_CINS, &len);
                if ( iseq==-2 ) break;
                if ( iseq==-1 ) continue;
                if ( len!=ifrq->len[cvar->which] ) continue;
                for (k=0; k<ifrq->len[cvar->which]; k++)
                    if ( bam_seqi(cigar.seq,iseq+k)!=ifrq->nt16_seq[cvar->which][k] ) break;
                if ( k==ifrq->len[cvar->which] ) hap |= 1<<j;
            }
            else if ( cvar->vtype==del )
            {
                int len;
                del_freq_t *dfrq = &rcns->del_freq[cvar->pos - rcns->beg];
                int ret = cstate_seek_op_fwd(&cigar, cvar->pos+1, BAM_CDEL, &len);
                if ( ret==-2 ) break;
                if ( ret==-1 ) continue;
                if ( len!=dfrq->len[cvar->which] ) continue;
                hap |= 1<<j;
            }
        }
        rcns->hap_freq[hap]++;
    }
    return 0;
}

typedef struct
{
    int haplotype, count;
}
ii_t;

static int ii_cmp(const void *a, const void *b)
{
    if ( ((ii_t*)a)->count > ((ii_t*)b)->count ) return -1;
    if ( ((ii_t*)a)->count < ((ii_t*)b)->count ) return 1;
    return 0;
}

// Select two most common haplotypes trying to account for 1bp errors. Haplotypes
// are represented as 8-bit numbers, each bit corresponds to one candidate variant.
static int correct_haplotype_errors(read_cns_t *rcns)
{
    int i,j, tot = 0;
    ii_t freq[NHAP];
    for (i=0; i<NHAP; i++)
    {
        freq[i].haplotype = i;
        freq[i].count = rcns->hap_freq[i];
        tot += rcns->hap_freq[i];
    }
    qsort(freq, NHAP, sizeof(ii_t), ii_cmp);    // sort haplotypes in descending order
    for (i=NHAP-1; i>1; i--)
    {
        if ( !freq[i].count ) continue;
        if ( freq[1].count > tot - freq[0].count - freq[1].count ) break;   // the top2 hapotypes cannot change anymore

        // Find a similar haplotype with the highest frequency. Assuming errors go in 0->1
        // direction only and considering one error only.
        int count = freq[i].count, max_hap = 0;
        for (j=0; j<MAX_NCVAR; j++)
        {
            if ( !(freq[i].haplotype & (1U<<j)) ) continue; // j-th bit not set in this haplotype
            int hap = freq[i].haplotype ^ (1U<<j);          // toggle j-th bit
            assert( hap>=0 && hap<NHAP );
            if ( count < rcns->hap_freq[hap] ) count = rcns->hap_freq[hap], max_hap = hap;
        }
        if ( count == freq[i].count ) continue;

        // Update frequency and sort the two modified elements
        count = freq[i].count;
        freq[i].count = 0;
        rcns->hap_freq[freq[i].haplotype] = 0;
        rcns->hap_freq[max_hap] += count;
        for (j=i+1; j<NHAP; j++)
        {
            if ( !freq[j].count ) break;
            ii_t tmp = freq[j-1]; freq[j-1] = freq[j]; freq[j] = tmp;
        }
        for (j=i-1; j>=0; j--)
        {
            if ( freq[j].haplotype==max_hap ) freq[j].count += count;   // update the best matching haplotype
            if ( freq[j].count < freq[j+1].count )
            {
                ii_t tmp = freq[j]; freq[j] = freq[j+1]; freq[j+1] = tmp;
            }
        }
    }

    // Use only one consensus if the next best haplotype is populated by less than 10% of reads
    rcns->ncns = ((float)freq[1].count / (freq[0].count + freq[1].count) < 0.1) ? 1 : 2;

    // Remove unused candidate variants from the top two haplotypes
    int hap0 = freq[0].haplotype;
    int hap1 = rcns->ncns==2 ? freq[1].haplotype : 0;
    rcns->cns_hap[0] = 0;
    rcns->cns_hap[1] = 0;
    for (i=0,j=0; i<MAX_NCVAR; i++)
    {
        if ( !((hap0|hap1) & (1U<<i)) ) continue;   // unused candidate variant, skip
        if ( i!=j ) rcns->cvar[j] = rcns->cvar[i];
        if ( hap0 & (1U<<i) ) rcns->cns_hap[0] |= 1U<<j;
        if ( hap1 & (1U<<i) ) rcns->cns_hap[1] |= 1U<<j;
        j++;
    }
    rcns->ncvar = j;

#if DEBUG_RCNS
    // This only matters for debugging print
    memset(rcns->hap_freq,0,NHAP*sizeof(*rcns->hap_freq));
    rcns->hap_freq[rcns->cns_hap[1]] = freq[1].count;   // NB: the order matters when ncns==1
    rcns->hap_freq[rcns->cns_hap[0]] = freq[0].count;
#endif

    return 0;
}


// Check how frequent are insertions adjacent to the j-th position. Note that reads with an
// insertion usually increment also bfreq counts at this position, but not necessarily so,
// therefore the counts are approximate
static inline void apply_consensus_insertion(read_cns_t *rcns, cns_seq_t *cns, int j, int ivar)
{
    // Only apply consensus insertions that are not being tested by bam2bcf_iaux, i.e. not at the current pos
    hts_pos_t ref_pos = rcns->beg + j;
    if ( rcns->pos == ref_pos ) return;

    // Only apply when there is no insertion at this position registered as a variant
    while ( ivar < rcns->ncvar && rcns->cvar[ivar].pos == ref_pos )
    {
        if ( rcns->cvar[ivar].vtype == ins ) return;
        ivar++;
    }

    base_freq_t *bfreq = rcns->base_freq;
    ins_freq_t *ifreq  = rcns->ins_freq;
    int k, nreads = 0;
    for (k=0; k<BF_DEL; k++) nreads += bfreq[j].base[k];
    int max_freq = 0, kmax = 0;
    for (k=0; k<NI && ifreq[j].len[k]; k++)
        if ( max_freq < ifreq[j].freq[k] ) max_freq = ifreq[j].freq[k], kmax = k;

    // Include consensus insertion only if it has more than half of the reads
    if ( nreads > max_freq*2 ) return;

    int len = ifreq[j].len[kmax];
    char *seq = ifreq[j].nt16_seq[kmax];
    for (k=0; k<len; k++)
        cns->seq[cns->nseq++] = seq_nt16_int[(int)seq[k]];
}

// For each position of the realignment window apply either the candidate variants
// from ith haplotype or decide on the base/ins/del by majority vote
static void create_consensus(read_cns_t *rcns, const char *ref, int ith)
{
    int n = rcns->end - rcns->beg + 1;
    cns_seq_t *cns = &rcns->cns[ith];
    base_freq_t *bfreq = rcns->base_freq;
    ins_freq_t *ifreq  = rcns->ins_freq;
    del_freq_t *dfreq  = rcns->del_freq;
    hts_pos_t prev_pos = 0;
    int j,k, ivar = 0;
    for (j=0; j<n; j++)
    {
        hts_pos_t ref_pos = rcns->beg + j;
        if ( rcns->pos == ref_pos ) cns->ipos = cns->nseq;

        while ( ivar < rcns->ncvar && rcns->cvar[ivar].pos < ref_pos ) ivar++;

        if ( ivar >= rcns->ncvar || rcns->cvar[ivar].pos != ref_pos )
        {
            // This position is not recognised as a het variant so take the most frequent base, including
            // a deletion if that is most frequent. However, for deleted bases make sure they are not part
            // of the deletion that is being tested at this positions
            int max_freq = 0, kmax = seq_nt16_int[seq_nt16_table[(int)ref[j]]];
            int nk = ( ref_pos < rcns->pos || ref_pos > rcns->pos + rcns->max_del ) ? BF_DEL+1 : BF_DEL;
            for (k=0; k<nk; k++)
                if ( max_freq < bfreq[j].base[k] ) max_freq = bfreq[j].base[k], kmax = k;

            if ( kmax!=BF_DEL )  // the most frequent base can be a deletion
            {
                prev_pos = ref_pos;
                cns->seq[cns->nseq++] = kmax;
            }
            // Only apply consensus insertions that are not being tested by bam2bcf_iaux, i.e. not at the current pos
            apply_consensus_insertion(rcns, cns, j, ivar);
            continue;
        }
        int which = rcns->cvar[ivar].which;
        if ( !(rcns->cns_hap[ith] & (1U<<ivar)) )
        {
            // This position has a heterozygous variant but not in this haplotype. Take the
            // most frequent base different from the ivar-th variant
            int max_freq = 0, kmax = seq_nt16_int[seq_nt16_table[(int)ref[j]]];
            for (k=0; k<6; k++)
            {
                if ( rcns->cvar[ivar].vtype==snv && rcns->cvar[ivar].which==k ) continue;
                if ( max_freq < bfreq[j].base[k] ) max_freq = bfreq[j].base[k], kmax = k;
            }
            if ( kmax!=BF_DEL && (!cns->nseq || prev_pos != ref_pos) )
            {
                prev_pos = ref_pos;
                cns->seq[cns->nseq++] = kmax;
            }
            apply_consensus_insertion(rcns, cns, j, ivar);
            continue;
        }
        if ( rcns->cvar[ivar].vtype == snv )
        {
            prev_pos = ref_pos;
            cns->seq[cns->nseq++] = which;
            apply_consensus_insertion(rcns, cns, j, ivar);
            continue;
        }

        // There can be multiple variants at this position, for example snv+ins. SNVs come first
        // thanks to cvar_pos_cmp(), make sure the base has not been added already.
        if ( !cns->nseq || prev_pos != ref_pos )
        {
            int max_freq = 0, kmax = seq_nt16_int[seq_nt16_table[(int)ref[j]]];
            for (k=0; k<6; k++)
            {
                if ( rcns->cvar[ivar].vtype==snv && rcns->cvar[ivar].which==k ) continue;
                if ( max_freq < bfreq[j].base[k] ) max_freq = bfreq[j].base[k], kmax = k;
            }
            if ( kmax!=BF_DEL )
            {
                prev_pos = ref_pos;
                cns->seq[cns->nseq++] = kmax;
            }
        }
        if ( rcns->cvar[ivar].vtype == ins )
        {
            int len = ifreq[j].len[which];
            char *seq = ifreq[j].nt16_seq[which];
            for (k=0; k<len; k++)
            {
                prev_pos = ref_pos;
                cns->seq[cns->nseq++] = seq_nt16_int[(int)seq[k]];
            }
        }
        else if ( rcns->cvar[ivar].vtype == del ) j += dfreq[j].len[which];
    }
}

// The algorithm:
//  1. Identify heterozygous variant positions
//  2. Sort variants by abs(variant_allele_freq-0.5) in descending order
//  3. Take the top sorted variants (up to 8 to fit in uint8_t) and count the number of
//      corresponding reads to create frequency spectrum
//  4. Correct errors, collapse to the requested number of haplotypes (consensus sequences)
//      using majority vote for the distribution tail
cns_seq_t *rcns_get_consensus(read_cns_t *rcns, const char *ref)
{
    debug_print_base_freqs(rcns, ref);

    select_candidate_variants(rcns, ref);
    debug_print_candidate_variants(rcns);

    if ( rcns->ncvar )
    {
        create_haplotype_frequency_spectrum(rcns);
        debug_print_haplotype_frequency_spectrum(rcns);

        correct_haplotype_errors(rcns);
        debug_print_candidate_variants(rcns);
        debug_print_haplotype_frequency_spectrum(rcns);
    }
    else
    {
        rcns->cns_hap[0] = 0;
        rcns->ncns = 1;
    }

    // create consensus
    int i;
    for (i=0; i<rcns->ncns; i++) create_consensus(rcns,ref,i);
    debug_print_consensus(rcns,ref);

    return rcns->cns;
}
