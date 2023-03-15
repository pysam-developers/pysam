#include "bcftools.pysam.h"

/* The MIT License

   Copyright (c) 2014-2023 Genome Research Ltd.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <unistd.h>
#include <ctype.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kseq.h>
#include <htslib/bgzf.h>
#include "regidx.h"
#include "bcftools.h"
#include "rbuf.h"
#include "filter.h"
#include "smpl_ilist.h"

// Logic of the filters: include or exclude sites which match the filters?
#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define PICK_REF   1
#define PICK_ALT   2
#define PICK_LONG  4
#define PICK_SHORT 8
#define PICK_IUPAC 16

#define TO_UPPER 0
#define TO_LOWER 1

typedef struct
{
    int num;                // number of ungapped blocks in this chain
    int *block_lengths;     // length of the ungapped blocks in this chain
    int *ref_gaps;          // length of the gaps on the reference sequence between blocks
    int *alt_gaps;          // length of the gaps on the alternative sequence between blocks
    int ori_pos;
    int ref_last_block_ori; // start position on the reference sequence of the following ungapped block (0-based)
    int alt_last_block_ori; // start position on the alternative sequence of the following ungapped block (0-based)
}
chain_t;

#define MASK_LC 1
#define MASK_UC 2
#define MASK_SKIP(x) (((x)->with!=MASK_LC && (x)->with!=MASK_UC) ? 1 : 0)
typedef struct
{
    char *fname, with;
    regidx_t *idx;
    regitr_t *itr;
}
mask_t;

typedef struct
{
    kstring_t fa_buf;   // buffered reference sequence
    int fa_ori_pos;     // start position of the fa_buffer (wrt original sequence)
    int fa_frz_pos;     // protected position to avoid conflicting variants (last pos for SNPs/ins)
    int fa_mod_off;     // position difference of fa_frz_pos in the ori and modified sequence (ins positive)
    int fa_frz_mod;     // the fa_buf offset of the protected fa_frz_pos position, includes the modified sequence
    int fa_end_pos;     // region's end position in the original sequence
    int fa_length;      // region's length in the original sequence (in case end_pos not provided in the FASTA header)
    int fa_case;        // output upper case or lower case: TO_UPPER|TO_LOWER
    int fa_src_pos;     // last genomic coordinate read from the input fasta (0-based)
    char prev_base;     // this is only to validate the REF allele in the VCF - the modified fa_buf cannot be used for inserts following deletions, see 600#issuecomment-383186778
    int prev_base_pos;  // the position of prev_base
    int prev_is_insert;

    rbuf_t vcf_rbuf;
    bcf1_t **vcf_buf;
    int nvcf_buf, rid;
    char *chr, *chr_prefix;

    mask_t *mask;
    int nmask;

    int chain_id;       // chain_id, to provide a unique ID to each chain in the chain output
    chain_t *chain;     // chain structure to store the sequence of ungapped blocks between the ref and alt sequences
                        // Note that the chain is re-initialised for each chromosome/seq_region

    filter_t *filter;
    char *filter_str;
    int filter_logic;   // include or exclude sites which match the filters? One of FLT_INCLUDE/FLT_EXCLUDE

    bcf_srs_t *files;
    bcf_hdr_t *hdr;
    FILE *fp_out;
    FILE *fp_chain;
    char **argv;
    int argc, output_iupac, iupac_GTs, haplotype, allele, isample, napplied;
    uint8_t *iupac_bitmask, *iupac_als;
    int miupac_bitmask, miupac_als;
    char *fname, *ref_fname, *sample, *sample_fname, *output_fname, *mask_fname, *chain_fname, missing_allele, absent_allele;
    char mark_del, mark_ins, mark_snv;
    smpl_ilist_t *smpl;
}
args_t;

static void destroy_chain(chain_t *chain)
{
    if ( !chain ) return;
    free(chain->ref_gaps);
    free(chain->alt_gaps);
    free(chain->block_lengths);
    free(chain);
}
static chain_t* init_chain(chain_t *chain, int ref_ori_pos)
{
    if ( chain ) destroy_chain(chain);
    chain = (chain_t*) calloc(1,sizeof(chain_t));
    chain->num = 0;
    chain->block_lengths = NULL;
    chain->ref_gaps = NULL;
    chain->alt_gaps = NULL;
    chain->ori_pos = ref_ori_pos;
    chain->ref_last_block_ori = ref_ori_pos;
    chain->alt_last_block_ori = ref_ori_pos;
    return chain;
}

static void print_chain(args_t *args)
{
    /*
        Example chain format (see: https://genome.ucsc.edu/goldenPath/help/chain.html):
        chain 1 500 + 480 500 1 501 + 480 501 1
        12 3 1
        1 0 3
        484

        chain line is:
        - chain
        - score (sum of the length of ungapped block in this case)
        - ref_seqname (from the fasta header, parsed by htslib)
        - ref_seqlength (from the fasta header)
        - ref_strand (+ or -; always + for bcf-consensus)
        - ref_start (as defined in the fasta header)
        - ref_end (as defined in the fasta header)
        - alt_seqname (same as ref_seqname as bcf-consensus only considers SNPs and indels)
        - alt_seqlength (adjusted to match the length of the alt sequence)
        - alt_strand (+ or -; always + for bcf-consensus)
        - alt_start (same as ref_start, as no edits are recorded/applied before that position)
        - alt_end (adjusted to match the length of the alt sequence)
        - chain_num (just an auto-increment id)

        the other (sorted) lines are:
        - length of the ungapped alignment block
        - gap on the ref sequence between this and the next block (all but the last line)
        - gap on the alt sequence between this and the next block (all but the last line)
    */
    chain_t *chain = args->chain;
    int n = chain->num;
    int ref_end_pos = args->fa_length + chain->ori_pos;
    int last_block_size = ref_end_pos - chain->ref_last_block_ori;
    int alt_end_pos = chain->alt_last_block_ori + last_block_size;
    int score = 0;
    for (n=0; n<chain->num; n++) {
        score += chain->block_lengths[n];
    }
    score += last_block_size;
    fprintf(args->fp_chain, "chain %d %s %d + %d %d %s %d + %d %d %d\n", score, args->chr, ref_end_pos, chain->ori_pos, ref_end_pos, args->chr, alt_end_pos, chain->ori_pos, alt_end_pos, ++args->chain_id);
    for (n=0; n<chain->num; n++) {
        fprintf(args->fp_chain, "%d %d %d\n", chain->block_lengths[n], chain->ref_gaps[n], chain->alt_gaps[n]);
    }
    fprintf(args->fp_chain, "%d\n\n", last_block_size);
}

static void push_chain_gap(chain_t *chain, int ref_start, int ref_len, int alt_start, int alt_len)
{
    // fprintf(bcftools_stderr, "push_chain_gap(chain=%p, ref_start=%d, ref_len=%d, alt_start=%d, alt_len=%d)\n", chain, ref_start, ref_len, alt_start, alt_len);
    int num = chain->num;

    if (num && ref_start <= chain->ref_last_block_ori) {
        // In case this variant is back-to-back with the previous one
        chain->ref_last_block_ori = ref_start + ref_len;
        chain->alt_last_block_ori = alt_start + alt_len;
        chain->ref_gaps[num-1] += ref_len;
        chain->alt_gaps[num-1] += alt_len;

    } else {
        // Extend the ungapped blocks, store the gap length
        chain->block_lengths = (int*) realloc(chain->block_lengths, (num + 1) * sizeof(int));
        chain->ref_gaps = (int*) realloc(chain->ref_gaps, (num + 1) * sizeof(int));
        chain->alt_gaps = (int*) realloc(chain->alt_gaps, (num + 1) * sizeof(int));
        chain->block_lengths[num] = ref_start - chain->ref_last_block_ori;
        chain->ref_gaps[num] = ref_len;
        chain->alt_gaps[num] = alt_len;
        // Update the start positions of the next block
        chain->ref_last_block_ori = ref_start + ref_len;
        chain->alt_last_block_ori = alt_start + alt_len;
        // Increment the number of ungapped blocks
        chain->num++;
    }
}

static void init_data(args_t *args)
{
    args->files = bcf_sr_init();
    args->files->require_index = 1;
    if ( !bcf_sr_add_reader(args->files,args->fname) ) error("Failed to read from %s: %s\n", !strcmp("-",args->fname)?"standard input":args->fname, bcf_sr_strerror(args->files->errnum));
    args->hdr = args->files->readers[0].header;
    args->isample = -1;
    if ( !args->sample )
        args->smpl = smpl_ilist_init(args->hdr,NULL,0,SMPL_NONE|SMPL_VERBOSE);
    else if ( args->sample && strcmp("-",args->sample) )
    {
        args->smpl = smpl_ilist_init(args->hdr,args->sample,0,SMPL_NONE|SMPL_VERBOSE);
        if ( args->smpl && !args->smpl->n ) error("No matching sample found\n");
    }
    else if ( args->sample_fname )
    {
        args->smpl = smpl_ilist_init(args->hdr,args->sample_fname,1,SMPL_NONE|SMPL_VERBOSE);
        if ( args->smpl && !args->smpl->n ) error("No matching sample found\n");
    }
    if ( args->smpl )
    {
        if ( args->haplotype || args->allele )
        {
            if ( args->smpl->n > 1 ) error("Too many samples, only one can be used with -H\n");
            args->isample = args->smpl->idx[0];
        }
        else
            args->iupac_GTs = 1;
    }
    int i;
    for (i=0; i<args->nmask; i++)
    {
        mask_t *mask = &args->mask[i];
        mask->idx = regidx_init(mask->fname,NULL,NULL,0,NULL);
        if ( !mask->idx ) error("Failed to initialize mask regions\n");
        mask->itr = regitr_init(mask->idx);
    }
    // In case we want to store the chains
    if ( args->chain_fname )
    {
        args->fp_chain = fopen(args->chain_fname,"w");
        if ( ! args->fp_chain ) error("Failed to create %s: %s\n", args->chain_fname, strerror(errno));
        args->chain_id = 0;
    }
    rbuf_init(&args->vcf_rbuf, 100);
    args->vcf_buf = (bcf1_t**) calloc(args->vcf_rbuf.m, sizeof(bcf1_t*));
    if ( args->output_fname ) {
        args->fp_out = fopen(args->output_fname,"w");
        if ( ! args->fp_out ) error("Failed to create %s: %s\n", args->output_fname, strerror(errno));
    }
    else args->fp_out = bcftools_stdout;
    if ( args->isample<0 && !args->iupac_GTs ) fprintf(bcftools_stderr,"Note: the --samples option not given, applying all records regardless of the genotype\n");
    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);
    args->rid = -1;
}
static void add_mask(args_t *args, char *fname)
{
    args->nmask++;
    args->mask = (mask_t*)realloc(args->mask,args->nmask*sizeof(*args->mask));
    mask_t *mask = &args->mask[args->nmask-1];
    mask->fname = fname;
    mask->with  = 'N';
}
static void add_mask_with(args_t *args, char *with)
{
    if ( !args->nmask ) error("The --mask-with option must follow --mask\n");
    mask_t *mask = &args->mask[args->nmask-1];
    if ( !strcasecmp(with,"uc") ) mask->with = MASK_UC;
    else if ( !strcasecmp(with,"lc") ) mask->with = MASK_LC;
    else if ( strlen(with)!=1 ) error("Expected \"lc\", \"uc\", or a single character with the --mask-with option\n");
    else mask->with = *with;
}
static void destroy_data(args_t *args)
{
    free(args->iupac_als);
    free(args->iupac_bitmask);
    if (args->filter) filter_destroy(args->filter);
    if ( args->smpl ) smpl_ilist_destroy(args->smpl);
    bcf_sr_destroy(args->files);
    int i;
    for (i=0; i<args->vcf_rbuf.m; i++)
        if ( args->vcf_buf[i] ) bcf_destroy1(args->vcf_buf[i]);
    free(args->vcf_buf);
    free(args->fa_buf.s);
    free(args->chr);
    for (i=0; i<args->nmask; i++)
    {
        mask_t *mask = &args->mask[i];
        regidx_destroy(mask->idx);
        regitr_destroy(mask->itr);
    }
    free(args->mask);
    if ( args->chain_fname )
        if ( fclose(args->fp_chain) ) error("Close failed: %s\n", args->chain_fname);
    if ( fclose(args->fp_out) ) error("Close failed: %s\n", args->output_fname);
    destroy_chain(args->chain);
}

static void init_region(args_t *args, char *line)
{
    char *ss, *se = line;
    while ( *se && !isspace(*se) && *se!=':' ) se++;
    int from = 0, to = 0;
    char tmp = 0, *tmp_ptr = NULL;
    if ( *se )
    {
        tmp = *se; *se = 0; tmp_ptr = se;
        ss = ++se;
        from = strtol(ss,&se,10);
        if ( ss==se || !*se || *se!='-' ) from = 0;
        else
        {
            from--;
            ss = ++se;
            to = strtol(ss,&se,10);
            if ( ss==se || (*se && !isspace(*se)) ) { from = 0; to = 0; }
            else to--;
        }
    }
    free(args->chr);
    args->chr = strdup(line);
    args->rid = bcf_hdr_name2id(args->hdr,line);
    if ( args->rid<0 ) fprintf(bcftools_stderr,"Warning: Sequence \"%s\" not in %s\n", line,args->fname);
    args->prev_base_pos = -1;
    args->fa_buf.l  = 0;
    args->fa_length = 0;
    args->fa_end_pos = to;
    args->fa_ori_pos = from;
    args->fa_src_pos = from;
    args->fa_mod_off = 0;
    args->fa_frz_pos = -1;
    args->fa_frz_mod = -1;
    args->fa_case    = -1;
    args->vcf_rbuf.n = 0;
    bcf_sr_seek(args->files,line,args->fa_ori_pos);
    if ( tmp_ptr ) *tmp_ptr = tmp;
    fprintf(args->fp_out,">%s%s\n",args->chr_prefix?args->chr_prefix:"",line);
    if ( args->chain_fname )
        args->chain = init_chain(args->chain, args->fa_ori_pos);
}

static bcf1_t **next_vcf_line(args_t *args)
{
    if ( args->vcf_rbuf.n )
    {
        int i = rbuf_shift(&args->vcf_rbuf);
        return &args->vcf_buf[i];
    }
    while ( bcf_sr_next_line(args->files) )
    {
        if ( args->filter )
        {
            int is_ok = filter_test(args->filter, bcf_sr_get_line(args->files,0), NULL);
            if ( args->filter_logic & FLT_EXCLUDE ) is_ok = is_ok ? 0 : 1;
            if ( !is_ok ) continue;
        }
        return &args->files->readers[0].buffer[0];
    }
    return NULL;
}
static void unread_vcf_line(args_t *args, bcf1_t **rec_ptr)
{
    bcf1_t *rec = *rec_ptr;
    if ( args->vcf_rbuf.n >= args->vcf_rbuf.m )
        error("FIXME: too many overlapping records near %s:%"PRId64"\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);

    // Insert the new record in the buffer. The line would be overwritten in
    // the next bcf_sr_next_line call, therefore we need to swap it with an
    // unused one
    int i = rbuf_append(&args->vcf_rbuf);
    if ( !args->vcf_buf[i] ) args->vcf_buf[i] = bcf_init1();
    bcf1_t *tmp = rec; *rec_ptr = args->vcf_buf[i]; args->vcf_buf[i] = tmp;
}
static void flush_fa_buffer(args_t *args, int len)
{
    if ( !args->fa_buf.l ) return;
    int nwr = 0;
    while ( nwr + 60 <= args->fa_buf.l )
    {
        if ( fwrite(args->fa_buf.s+nwr,1,60,args->fp_out) != 60 ) error("Could not write: %s\n", args->output_fname);
        if ( fwrite("\n",1,1,args->fp_out) != 1 ) error("Could not write: %s\n", args->output_fname);
        nwr += 60;
    }
    if ( nwr )
        args->fa_ori_pos += nwr;

    args->fa_frz_mod -= nwr;

    if ( len )
    {
        // not finished on this chr yet and the buffer cannot be emptied completely
        if ( nwr && nwr < args->fa_buf.l )
            memmove(args->fa_buf.s,args->fa_buf.s+nwr,args->fa_buf.l-nwr);
        args->fa_buf.l -= nwr;
        return;
    }

    // empty the whole buffer
    if ( nwr == args->fa_buf.l ) { args->fa_buf.l = 0; return; }

    if ( fwrite(args->fa_buf.s+nwr,1,args->fa_buf.l - nwr,args->fp_out) != args->fa_buf.l - nwr ) error("Could not write: %s\n", args->output_fname);
    if ( fwrite("\n",1,1,args->fp_out) != 1 ) error("Could not write: %s\n", args->output_fname);

    args->fa_ori_pos += args->fa_buf.l - nwr - args->fa_mod_off;
    args->fa_mod_off = 0;
    args->fa_buf.l = 0;
}
static void apply_absent(args_t *args, hts_pos_t pos)
{
    if ( !args->fa_buf.l || pos <= args->fa_frz_pos + 1 || pos <= args->fa_ori_pos ) return;

    int ie = pos && pos - args->fa_ori_pos + args->fa_mod_off < args->fa_buf.l ? pos - args->fa_ori_pos + args->fa_mod_off : args->fa_buf.l;
    int ib = args->fa_frz_mod < 0 ? 0 : args->fa_frz_mod;
    int i;
    for (i=ib; i<ie; i++)
        args->fa_buf.s[i] = args->absent_allele;
}
static void freeze_ref(args_t *args, bcf1_t *rec)
{
    if ( args->fa_frz_pos >= rec->pos + rec->rlen - 1 ) return;
    args->fa_frz_pos = rec->pos + rec->rlen - 1;
    args->fa_frz_mod = rec->pos - args->fa_ori_pos + args->fa_mod_off + rec->rlen;
}
static char *mark_del(char *ref, int rlen, char *alt, int mark)
{
    char *out = malloc(rlen+1);
    int i;
    if ( alt )
    {
        int nalt = strlen(alt);
        for (i=0; i<nalt; i++) out[i] = alt[i];
    }
    else    // symbolic <DEL>
    {
        int nref = strlen(ref);
        for (i=0; i<nref; i++) out[i] = ref[i];
    }
    for (; i<rlen; i++) out[i] = mark;
    out[rlen] = 0;
    return out;
}
static void mark_ins(char *ref, char *alt, char mark)
{
    int i, nref = strlen(ref), nalt = strlen(alt);
    if ( mark=='l' )
        for (i=nref; i<nalt; i++) alt[i] = tolower(alt[i]);
    else
        for (i=nref; i<nalt; i++) alt[i] = toupper(alt[i]);
}
static void mark_snv(char *ref, char *alt, char mark)
{
    int i, nref = strlen(ref), nalt = strlen(alt);
    int n = nref < nalt ? nref : nalt;
    if ( mark=='l' )
    {
        for (i=0; i<n; i++)
            if ( tolower(ref[i])!=tolower(alt[i]) ) alt[i] = tolower(alt[i]);
    }
    else
    {
        for (i=0; i<n; i++)
            if ( tolower(ref[i])!=tolower(alt[i]) ) alt[i] = toupper(alt[i]);
    }
}
static void iupac_init(args_t *args, bcf1_t *rec)
{
    int i;
    hts_resize(uint8_t, rec->n_allele, &args->miupac_als, &args->iupac_als, 0);
    for (i=0; i<args->miupac_als; i++) args->iupac_als[i] = 0;
}
static int iupac_add_gt(args_t *args, bcf1_t *rec, uint8_t *gt, int ngt)
{
    int i, is_set = 0;
    for (i=0; i<ngt; i++)
    {
        if ( bcf_gt_is_missing(gt[i]) ) continue;
        if ( gt[i]==(uint8_t)bcf_int8_vector_end ) break;
        int ial = bcf_gt_allele(gt[i]);
        if ( ial >= rec->n_allele ) error("Invalid VCF, too few ALT alleles at %s:%"PRId64"\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
        args->iupac_als[ial] = 1;
        is_set = 1;
    }
    return is_set;
}
static int iupac_set_allele(args_t *args, bcf1_t *rec)
{
    int i,j, max_len = 0, alt_len = 0, ialt = -1, fallback_alt = -1;
    for (i=0; i<rec->n_allele; i++)
    {
        if ( !args->iupac_als[i] ) continue;
        if ( fallback_alt <=0 ) fallback_alt = i;
        int l = strlen(rec->d.allele[i]);
        for (j=0; j<l; j++)
            if ( iupac2bitmask(rec->d.allele[i][j]) < 0 ) break;
        if ( j<l ) continue; // symbolic allele, breakpoint or invalid character in the allele
        if ( l>max_len )
        {
            hts_resize(uint8_t, l, &args->miupac_bitmask, &args->iupac_bitmask, HTS_RESIZE_CLEAR);
            for (j=max_len; j<l; j++) args->iupac_bitmask[j] = 0;
            max_len = l;
        }
        if ( i>0 && l>alt_len )
        {
            alt_len = l;
            ialt = i;
        }
        for (j=0; j<l; j++)
            args->iupac_bitmask[j] |= iupac2bitmask(rec->d.allele[i][j]);
    }
    if ( alt_len > 0 )
    {
        for (j=0; j<alt_len; j++) rec->d.allele[ialt][j] = bitmask2iupac(args->iupac_bitmask[j]);
        return ialt;
    }
    if ( fallback_alt >= 0 ) return fallback_alt;
    return ialt;
}
static void apply_variant(args_t *args, bcf1_t *rec)
{
    static int warned_haplotype = 0;

    if ( args->absent_allele ) apply_absent(args, rec->pos);
    if ( rec->n_allele==1 && !args->missing_allele && !args->absent_allele ) { return; }

    int i,j;
    if ( args->mask )
    {
        char *chr = (char*)bcf_hdr_id2name(args->hdr,args->rid);
        int start = rec->pos;
        int end   = rec->pos + rec->rlen - 1;
        for (i=0; i<args->nmask; i++)
        {
            mask_t *mask = &args->mask[i];
            if ( MASK_SKIP(mask) && regidx_overlap(mask->idx, chr,start,end,NULL) ) return;
        }
    }

    int ialt = 1;    // the alternate allele
    if ( args->iupac_GTs )
    {
        bcf_unpack(rec, BCF_UN_FMT);
        bcf_fmt_t *fmt = bcf_get_fmt(args->hdr, rec, "GT");
        if ( !fmt ) return;
        if ( fmt->type!=BCF_BT_INT8 )
            error("Todo: GT field represented with BCF_BT_INT8, too many alleles at %s:%"PRId64"?\n",bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
        ialt = -1;
        int is_set = 0;
        iupac_init(args,rec);
        for (i=0; i<args->smpl->n; i++)
        {
            uint8_t *ptr = fmt->p + fmt->size*args->smpl->idx[i];
            is_set += iupac_add_gt(args, rec, ptr, fmt->n);
        }
        if ( !is_set && !args->missing_allele ) return;
        if ( is_set ) ialt = iupac_set_allele(args, rec);
    }
    else if ( args->isample >= 0 )
    {
        bcf_unpack(rec, BCF_UN_FMT);
        bcf_fmt_t *fmt = bcf_get_fmt(args->hdr, rec, "GT");
        if ( !fmt ) return;

        if ( fmt->type!=BCF_BT_INT8 )
            error("Todo: GT field represented with BCF_BT_INT8, too many alleles at %s:%"PRId64"?\n",bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
        uint8_t *ptr = fmt->p + fmt->size*args->isample;

        enum { use_hap, use_iupac, pick_one } action = use_hap;
        if ( args->allele==PICK_IUPAC )
        {
            if ( !args->haplotype ) action = use_iupac;
            if ( !bcf_gt_is_phased(ptr[0]) && !bcf_gt_is_phased(ptr[fmt->n-1]) ) action = use_iupac;
        }
        else if ( args->output_iupac ) action = use_iupac;
        else if ( !args->haplotype ) action = pick_one;

        if ( action==use_hap )
        {
            if ( args->haplotype > fmt->n )
            {
                if ( bcf_gt_is_missing(ptr[fmt->n-1]) || bcf_gt_is_missing(ptr[0]) )
                {
                    if ( !args->missing_allele ) return;
                    ialt = -1;
                }
                else
                {
                    if ( !warned_haplotype )
                    {
                        fprintf(bcftools_stderr, "Can't apply %d-th haplotype at %s:%"PRId64". (This warning is printed only once.)\n", args->haplotype,bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
                        warned_haplotype = 1;
                    }
                    return;
                }
            }
            else
            {
                ialt = (int8_t)ptr[args->haplotype-1];
                if ( bcf_gt_is_missing(ialt) || ialt==bcf_int8_vector_end )
                {
                    if ( !args->missing_allele ) return;
                    ialt = -1;
                }
                else
                    ialt = bcf_gt_allele(ialt);
            }
        }
        else if ( action==use_iupac )
        {
            ialt = -1;
            iupac_init(args,rec);
            int is_set = iupac_add_gt(args, rec, ptr, fmt->n);
            if ( !is_set && !args->missing_allele ) return;
            if ( is_set ) ialt = iupac_set_allele(args, rec);
        }
        else
        {
            int is_hom = 1;
            for (i=0; i<fmt->n; i++)
            {
                if ( bcf_gt_is_missing(ptr[i]) )
                {
                    if ( !args->missing_allele ) return;  // ignore missing or half-missing genotypes
                    ialt = -1;
                    break;
                }
                if ( ptr[i]==(uint8_t)bcf_int8_vector_end ) break;
                ialt = bcf_gt_allele(ptr[i]);
                if ( i>0 && ialt!=bcf_gt_allele(ptr[i-1]) ) { is_hom = 0; break; }
            }
            if ( !is_hom )
            {
                int prev_len = 0, jalt;
                for (i=0; i<fmt->n; i++)
                {
                    if ( ptr[i]==(uint8_t)bcf_int8_vector_end ) break;
                    jalt = bcf_gt_allele(ptr[i]);
                    if ( rec->n_allele <= jalt ) error("Broken VCF, too few alts at %s:%"PRId64"\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
                    if ( args->allele & (PICK_LONG|PICK_SHORT) )
                    {
                        int len = jalt==0 ? rec->rlen : strlen(rec->d.allele[jalt]);
                        if ( i==0 ) ialt = jalt, prev_len = len;
                        else if ( len == prev_len )
                        {
                            if ( args->allele & PICK_REF && jalt==0 ) ialt = jalt, prev_len = len;
                            else if ( args->allele & PICK_ALT && ialt==0 ) ialt = jalt, prev_len = len;
                        }
                        else if ( args->allele & PICK_LONG && len > prev_len ) ialt = jalt, prev_len = len;
                        else if ( args->allele & PICK_SHORT && len < prev_len ) ialt = jalt, prev_len = len;
                    }
                    else
                    {
                        if ( args->allele & PICK_REF && jalt==0 ) ialt = jalt;
                        else if ( args->allele & PICK_ALT && ialt==0 ) ialt = jalt;
                    }
                }
            }
        }
        if ( !ialt )
        {
            // ref allele
            if ( args->absent_allele ) freeze_ref(args,rec);
            return;
        }
        if ( rec->n_allele <= ialt ) error("Broken VCF, too few alts at %s:%"PRId64"\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
    }
    else if ( args->output_iupac && rec->n_allele>1 )
    {
        int ialt, alen = 0, mlen = 0;
        for (i=0; i<rec->n_allele; i++)
        {
            int l = strlen(rec->d.allele[i]);
            for (j=0; j<l; j++)
                if ( iupac2bitmask(rec->d.allele[i][j]) < 0 ) break;
            if ( j<l ) continue;    // symbolic allele, breakpoint or invalid character in the allele

            if ( l > mlen )
            {
                hts_expand(uint8_t,l,args->miupac_bitmask,args->iupac_bitmask);
                for (j=mlen; j<l; j++) args->iupac_bitmask[j] = 0;
                mlen = l;
            }
            if ( i>0 && l>alen )
            {
                alen = l;
                ialt = i;
            }
            for (j=0; j<l; j++)
                args->iupac_bitmask[j] |= iupac2bitmask(rec->d.allele[i][j]);
        }
        if ( alen > 0 )
            for (j=0; j<alen; j++) rec->d.allele[ialt][j] = bitmask2iupac(args->iupac_bitmask[j]);
        else
            ialt = 1;
    }

    if ( rec->n_allele==1 && ialt!=-1 )
    {
        // non-missing reference
        if ( args->absent_allele ) freeze_ref(args,rec);
        return;
    }
    if ( ialt==-1 )
    {
        char alleles[4];
        alleles[0] = rec->d.allele[0][0];
        alleles[1] = ',';
        alleles[2] = args->missing_allele;
        alleles[3] = 0;
        bcf_update_alleles_str(args->hdr, rec, alleles);
        ialt = 1;
    }

    // For some variant types POS+REF refer to the base *before* the event; in such case set trim_beg
    int trim_beg = 0;
    int var_type = bcf_get_variant_type(rec,ialt);
    int var_len  = rec->d.var[ialt].n;
    if ( var_type & VCF_INDEL )
    {
        // normally indel starts one base after, but not if the first base of the fa reference is deleted
        if ( rec->d.allele[0][0] == rec->d.allele[ialt][0] )
            trim_beg = 1;
        else
            trim_beg = 0;
    }
    else if ( (var_type & VCF_OTHER) && !strcasecmp(rec->d.allele[ialt],"<DEL>") )
    {
        trim_beg = 1;
        var_len  = 1 - rec->rlen;
    }
    else if ( (var_type & VCF_OTHER) && !strncasecmp(rec->d.allele[ialt],"<INS",4) ) trim_beg = 1;

    // Overlapping variant?
    if ( rec->pos <= args->fa_frz_pos )
    {
        // Can be still OK iff this is an insertion (and which does not follow another insertion, see #888).
        // This still may not be enough for more complicated cases with multiple duplicate positions
        // and other types in between. In such case let the user normalize the VCF and remove duplicates.

        int overlap = 0;
        if ( rec->pos < args->fa_frz_pos || !trim_beg || var_len==0 || args->prev_is_insert ) overlap = 1;

        if ( overlap )
        {
            fprintf(bcftools_stderr,"The site %s:%"PRId64" overlaps with another variant, skipping...\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
            return;
        }

    }

    char *alt_allele = rec->d.allele[ialt];
    int rmme_alt = 0;

    int len_diff = 0, alen = 0;
    int idx = rec->pos - args->fa_ori_pos + args->fa_mod_off;
    if ( idx<0 )
    {
        fprintf(bcftools_stderr,"Warning: ignoring overlapping variant starting at %s:%"PRId64"\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
        return;
    }
    if ( rec->rlen > args->fa_buf.l - idx )
    {
        rec->rlen = args->fa_buf.l - idx;
        if ( alt_allele[0]!='<' )
        {
            alen = strlen(alt_allele);
            if ( alen > rec->rlen )
            {
                fprintf(bcftools_stderr,"Warning: trimming variant \"%s\" starting at %s:%"PRId64"\n", alt_allele,bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
                alt_allele[rec->rlen] = 0;
            }
        }
    }
    if ( idx>=args->fa_buf.l )
        error("FIXME: %s:%"PRId64" .. idx=%d, ori_pos=%d, len=%"PRIu64", off=%d\n",bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1,idx,args->fa_ori_pos,(uint64_t)args->fa_buf.l,args->fa_mod_off);

    // sanity check the reference base
    if ( alt_allele[0]=='<' )
    {
        // TODO: symbolic deletions probably need more work above with PICK_SHORT|PICK_LONG

        if ( strcasecmp(alt_allele,"<DEL>") && strcasecmp(alt_allele,"<*>") && strcasecmp(alt_allele,"<NON_REF>") )
            error("Symbolic alleles other than <DEL>, <*> or <NON_REF> are currently not supported, e.g. \"%s\" at %s:%"PRId64".\n"
                  "Please use filtering expressions to exclude such sites, for example by running with: -e 'ALT~\"<.*>\"'\n",
                alt_allele,bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
        if ( !strcasecmp(alt_allele,"<DEL>") )
        {
            static int multibase_ref_del_warned = 0;
            if ( rec->d.allele[0][1]!=0 && !multibase_ref_del_warned )
            {
                fprintf(bcftools_stderr,
                    "Warning: one REF base is expected with <DEL>, assuming the actual deletion starts at POS+1 at %s:%"PRId64".\n"
                    "         (This warning is printed only once.)\n", bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1);
                multibase_ref_del_warned = 1;
            }
            if ( args->mark_del )   // insert dashes instead of delete sequence
            {
                alt_allele = mark_del(rec->d.allele[0], rec->rlen, NULL, args->mark_del);
                alen = rec->rlen;
                len_diff = 0;
                rmme_alt = 1;
            }
            else
            {
                len_diff = 1-rec->rlen;
                alt_allele = rec->d.allele[0];     // according to VCF spec, the first REF base must precede the event
                alen = 1;
            }
        }
        else
        {
            // <*>  or <NON_REF> .. gVCF, evidence for the reference allele throughout the whole block
            freeze_ref(args,rec);
            return;
        }
    }
    else if ( strncasecmp(rec->d.allele[0],args->fa_buf.s+idx,rec->rlen) )
    {
        // This is hacky, handle a special case: if SNP or an insert follows a deletion (AAC>A, C>CAA),
        // the reference base in fa_buf is lost and the check fails. We do not keep a buffer
        // with the original sequence as it should not be necessary, we should encounter max
        // one base overlap

        int fail = 1;
        if ( args->prev_base_pos==rec->pos && toupper(rec->d.allele[0][0])==toupper(args->prev_base) )
        {
            if ( rec->rlen==1 ) fail = 0;
            else if ( !strncasecmp(rec->d.allele[0]+1,args->fa_buf.s+idx+1,rec->rlen-1) ) fail = 0;
        }

        if ( fail )
        {
            char tmp = 0;
            if ( args->fa_buf.l - idx > rec->rlen )
            {
                tmp = args->fa_buf.s[idx+rec->rlen];
                args->fa_buf.s[idx+rec->rlen] = 0;
            }
            error(
                    "The fasta sequence does not match the REF allele at %s:%"PRId64":\n"
                    "   REF .vcf: [%s]\n"
                    "   ALT .vcf: [%s]\n"
                    "   REF .fa : [%s]%c%s\n",
                    bcf_seqname(args->hdr,rec),(int64_t) rec->pos+1, rec->d.allele[0], alt_allele, args->fa_buf.s+idx,
                    tmp?tmp:' ',tmp?args->fa_buf.s+idx+rec->rlen+1:""
                 );
        }
        alen = strlen(alt_allele);
        len_diff = alen - rec->rlen;

        if ( args->mark_del && len_diff<0 )
        {
            alt_allele = mark_del(rec->d.allele[0], rec->rlen, alt_allele, args->mark_del);
            alen = rec->rlen;
            len_diff = 0;
            rmme_alt = 1;
        }
    }
    else
    {
        alen = strlen(alt_allele);
        len_diff = alen - rec->rlen;

        if ( args->mark_del && len_diff<0 )
        {
            alt_allele = mark_del(rec->d.allele[0], rec->rlen, alt_allele, args->mark_del);
            alen = rec->rlen;
            len_diff = 0;
            rmme_alt = 1;
        }
    }

    args->fa_case = toupper(args->fa_buf.s[idx])==args->fa_buf.s[idx] ? TO_UPPER : TO_LOWER;
    if ( args->fa_case==TO_UPPER )
        for (i=0; i<alen; i++) alt_allele[i] = toupper(alt_allele[i]);
    else
        for (i=0; i<alen; i++) alt_allele[i] = tolower(alt_allele[i]);

    if ( args->mark_ins && len_diff>0 )
        mark_ins(rec->d.allele[0], alt_allele, args->mark_ins);
    if ( args->mark_snv )
        mark_snv(rec->d.allele[0], alt_allele, args->mark_snv);

    if ( len_diff <= 0 )
    {
        // deletion or same size event
        assert( args->fa_buf.l >= idx+rec->rlen );
        args->prev_base = args->fa_buf.s[idx+rec->rlen-1];
        args->prev_base_pos = rec->pos + rec->rlen - 1;
        args->prev_is_insert = 0;
        args->fa_frz_mod = idx + alen;

        for (i=trim_beg; i<alen; i++)
            args->fa_buf.s[idx+i] = alt_allele[i];

        if ( len_diff )
            memmove(args->fa_buf.s+idx+alen,args->fa_buf.s+idx+rec->rlen,args->fa_buf.l-idx-rec->rlen);
    }
    else
    {
        args->prev_is_insert = 1;
        args->prev_base_pos = rec->pos;

        // insertion
        ks_resize(&args->fa_buf, args->fa_buf.l + len_diff);
        memmove(args->fa_buf.s + idx + rec->rlen + len_diff, args->fa_buf.s + idx + rec->rlen, args->fa_buf.l - idx - rec->rlen);

        // This can get tricky, make sure the bases unchanged by the insertion do not overwrite preceeding variants.
        // For example, here we want to get TAA:
        //      POS REF ALT
        //      1   C   T
        //      1   C   CAA
        int ibeg = 0;
        while ( ibeg<alen && rec->d.allele[0][ibeg]==alt_allele[ibeg] && rec->pos + ibeg <= args->prev_base_pos  ) ibeg++;
        for (i=ibeg; i<alen; i++)
            args->fa_buf.s[idx+i] = alt_allele[i];

        args->fa_frz_mod = idx + alen - ibeg + 1;
    }
    if (args->chain && len_diff != 0)
    {
        // If first nucleotide of both REF and ALT are the same... (indels typically include the nucleotide before the variant)
        if ( strncasecmp(rec->d.allele[0],alt_allele,1) == 0)
        {
            // ...extend the block by 1 bp: start is 1 bp further and alleles are 1 bp shorter
            push_chain_gap(args->chain, rec->pos + 1, rec->rlen - 1, rec->pos + 1 + args->fa_mod_off, alen - 1);
        }
        else
        {
            // otherwise, just the coordinates of the variant as given
            push_chain_gap(args->chain, rec->pos, rec->rlen, rec->pos + args->fa_mod_off, alen);
        }
    }
    args->fa_buf.l += len_diff;
    args->fa_mod_off += len_diff;
    args->fa_frz_pos  = rec->pos + rec->rlen - 1;
    args->napplied++;
    if ( rmme_alt ) free(alt_allele);
}


static void mask_region(args_t *args, char *seq, int len)
{
    int start = args->fa_src_pos - len;
    int end   = args->fa_src_pos;
    int i;

    for (i=0; i<args->nmask; i++)
    {
        mask_t *mask = &args->mask[i];
        if ( !regidx_overlap(mask->idx, args->chr,start,end, mask->itr) ) continue;

        int idx_start, idx_end, j;
        while ( regitr_overlap(mask->itr) )
        {
            idx_start = mask->itr->beg - start;
            idx_end   = mask->itr->end - start;
            if ( idx_start < 0 ) idx_start = 0;
            if ( idx_end >= len ) idx_end = len - 1;
            if ( mask->with==MASK_UC )
                for (j=idx_start; j<=idx_end; j++) seq[j] = toupper(seq[j]);
            else if ( mask->with==MASK_LC )
                for (j=idx_start; j<=idx_end; j++) seq[j] = tolower(seq[j]);
            else
                for (j=idx_start; j<=idx_end; j++) seq[j] = mask->with;
        }
    }
}

static void consensus(args_t *args)
{
    BGZF *fasta = bgzf_open(args->ref_fname, "r");
    if ( !fasta ) error("Error reading %s\n", args->ref_fname);
    kstring_t str = {0,0,0};
    while ( bgzf_getline(fasta, '\n', &str) > 0 )
    {
        if ( str.s[0]=='>' )
        {
            // new sequence encountered
            if ( args->chain ) print_chain(args);

            // apply all cached variants and variants that might have been missed because of short fasta (see test/consensus.9.*)
            bcf1_t **rec_ptr = NULL;
            while ( args->rid>=0 && (rec_ptr = next_vcf_line(args)) )
            {
                bcf1_t *rec = *rec_ptr;
                if ( rec->rid!=args->rid || ( args->fa_end_pos && rec->pos > args->fa_end_pos ) ) break;
                apply_variant(args, rec);
            }
            if ( args->absent_allele )
            {
                int pos = 0;
                if ( args->vcf_rbuf.n && args->vcf_buf[args->vcf_rbuf.f]->rid==args->rid )
                    pos = args->vcf_buf[args->vcf_rbuf.f]->pos;
                apply_absent(args, pos);
            }
            flush_fa_buffer(args, 0);
            init_region(args, str.s+1);
            continue;
        }
        args->fa_length  += str.l;
        args->fa_src_pos += str.l;

        // determine if uppercase or lowercase is used in this fasta file
        if ( args->fa_case==-1 ) args->fa_case = toupper(str.s[0])==str.s[0] ? 1 : 0;

        if ( args->mask ) mask_region(args, str.s, str.l);
        kputs(str.s, &args->fa_buf);

        bcf1_t **rec_ptr = NULL;
        while ( args->rid>=0 && (rec_ptr = next_vcf_line(args)) )
        {
            bcf1_t *rec = *rec_ptr;

            // still the same chr and the same region? if not, fasta buf can be flushed
            if ( rec->rid!=args->rid || ( args->fa_end_pos && rec->pos > args->fa_end_pos ) )
            {
                // save the vcf record until next time and flush
                unread_vcf_line(args, rec_ptr);
                rec_ptr = NULL;
                break;
            }

            // is the vcf record well beyond cached fasta buffer? if yes, the buf can be flushed
            if ( args->fa_ori_pos + args->fa_buf.l - args->fa_mod_off <= rec->pos )
            {
                unread_vcf_line(args, rec_ptr);
                rec_ptr = NULL;
                break;
            }

            // is the cached fasta buffer full enough? if not, read more fasta, no flushing
            if ( args->fa_ori_pos + args->fa_buf.l - args->fa_mod_off < rec->pos + rec->rlen )
            {
                unread_vcf_line(args, rec_ptr);
                break;
            }
            apply_variant(args, rec);
        }
        if ( !rec_ptr )
        {
            if ( args->absent_allele ) apply_absent(args, args->fa_ori_pos - args->fa_mod_off + args->fa_buf.l);
            flush_fa_buffer(args, 60);
        }
    }
    bcf1_t **rec_ptr = NULL;
    while ( args->rid>=0 && (rec_ptr = next_vcf_line(args)) )
    {
        bcf1_t *rec = *rec_ptr;
        if ( rec->rid!=args->rid ) break;
        if ( args->fa_end_pos && rec->pos > args->fa_end_pos ) break;
        if ( args->fa_ori_pos + args->fa_buf.l - args->fa_mod_off <= rec->pos ) break;
        apply_variant(args, rec);
    }
    if (args->chain) print_chain(args);
    if ( args->absent_allele ) apply_absent(args, HTS_POS_MAX);
    flush_fa_buffer(args, 0);
    bgzf_close(fasta);
    free(str.s);
    fprintf(bcftools_stderr,"Applied %d variants\n", args->napplied);
}

static void usage(args_t *args)
{
    fprintf(bcftools_stderr, "\n");
    fprintf(bcftools_stderr, "About: Create consensus sequence by applying VCF variants to a reference fasta\n");
    fprintf(bcftools_stderr, "       file. By default, the program will apply all ALT variants. Using the\n");
    fprintf(bcftools_stderr, "       --samples (and, optionally, --haplotype) option will apply genotype\n");
    fprintf(bcftools_stderr, "       (or haplotype) calls from FORMAT/GT. The program ignores allelic depth\n");
    fprintf(bcftools_stderr, "       information, such as INFO/AD or FORMAT/AD.\n");
    fprintf(bcftools_stderr, "Usage: bcftools consensus [OPTIONS] <file.vcf.gz>\n");
    fprintf(bcftools_stderr, "Options:\n");
    fprintf(bcftools_stderr, "    -c, --chain FILE               Write a chain file for liftover\n");
    fprintf(bcftools_stderr, "    -a, --absent CHAR              Replace positions absent from VCF with CHAR\n");
    fprintf(bcftools_stderr, "    -e, --exclude EXPR             Exclude sites for which the expression is true (see man page for details)\n");
    fprintf(bcftools_stderr, "    -f, --fasta-ref FILE           Reference sequence in fasta format\n");
    fprintf(bcftools_stderr, "    -H, --haplotype WHICH          Choose which allele to use from the FORMAT/GT field, note\n");
    fprintf(bcftools_stderr, "                                   the codes are case-insensitive:\n");
    fprintf(bcftools_stderr, "                                       1: first allele from GT, regardless of phasing\n");
    fprintf(bcftools_stderr, "                                       2: second allele from GT, regardless of phasing\n");
    fprintf(bcftools_stderr, "                                       R: REF allele in het genotypes\n");
    fprintf(bcftools_stderr, "                                       A: ALT allele\n");
    fprintf(bcftools_stderr, "                                       I: IUPAC code for all genotypes\n");
    fprintf(bcftools_stderr, "                                       LR,LA: longer allele and REF/ALT if equal length\n");
    fprintf(bcftools_stderr, "                                       SR,SA: shorter allele and REF/ALT if equal length\n");
    fprintf(bcftools_stderr, "                                       1pIu,2pIu: first/second allele for phased and IUPAC code for unphased GTs\n");
    fprintf(bcftools_stderr, "    -i, --include EXPR             Select sites for which the expression is true (see man page for details)\n");
    fprintf(bcftools_stderr, "    -I, --iupac-codes              Output IUPAC codes based on FORMAT/GT, use -s/-S to subset samples\n");
    fprintf(bcftools_stderr, "        --mark-del CHAR            Instead of removing sequence, insert CHAR for deletions\n");
    fprintf(bcftools_stderr, "        --mark-ins uc|lc           Highlight insertions in uppercase (uc) or lowercase (lc), leaving the rest as is\n");
    fprintf(bcftools_stderr, "        --mark-snv uc|lc           Highlight substitutions in uppercase (uc) or lowercase (lc), leaving the rest as is\n");
    fprintf(bcftools_stderr, "    -m, --mask FILE                Replace regions according to the next --mask-with option. The default is --mask-with N\n");
    fprintf(bcftools_stderr, "        --mask-with CHAR|uc|lc     Replace with CHAR (skips overlapping variants); change to uppercase (uc) or lowercase (lc)\n");
    fprintf(bcftools_stderr, "    -M, --missing CHAR             Output CHAR instead of skipping a missing genotype \"./.\"\n");
    fprintf(bcftools_stderr, "    -o, --output FILE              Write output to a file [standard output]\n");
    fprintf(bcftools_stderr, "    -p, --prefix STRING            Prefix to add to output sequence names\n");
    fprintf(bcftools_stderr, "    -s, --samples LIST             Comma-separated list of samples to include, \"-\" to ignore samples and use REF,ALT\n");
    fprintf(bcftools_stderr, "    -S, --samples-file FILE        File of samples to include\n");
    fprintf(bcftools_stderr, "Examples:\n");
    fprintf(bcftools_stderr, "   # Get the consensus for one region. The fasta header lines are then expected\n");
    fprintf(bcftools_stderr, "   # in the form \">chr:from-to\".\n");
    fprintf(bcftools_stderr, "   samtools faidx ref.fa 8:11870-11890 | bcftools consensus in.vcf.gz > out.fa\n");
    fprintf(bcftools_stderr, "\n");
    fprintf(bcftools_stderr, "   # See also http://samtools.github.io/bcftools/howtos/consensus-sequence.html\n");
    fprintf(bcftools_stderr, "\n");
    bcftools_exit(1);
}

int main_consensus(int argc, char *argv[])
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;

    static struct option loptions[] =
    {
        {"mark-del",required_argument,NULL,1},
        {"mark-ins",required_argument,NULL,2},
        {"mark-snv",required_argument,NULL,3},
        {"mask-with",1,0,4},
        {"exclude",required_argument,NULL,'e'},
        {"include",required_argument,NULL,'i'},
        {"sample",1,0,'s'},
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {"iupac-codes",0,0,'I'},
        {"haplotype",1,0,'H'},
        {"output",1,0,'o'},
        {"fasta-ref",1,0,'f'},
        {"mask",1,0,'m'},
        {"missing",1,0,'M'},
        {"absent",1,0,'a'},
        {"chain",1,0,'c'},
        {"prefix",required_argument,0,'p'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "h?s:S:1Ii:e:H:f:o:m:c:M:p:a:",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case  1 : args->mark_del = optarg[0]; break;
            case  2 :
                if ( !strcasecmp(optarg,"uc") ) args->mark_ins = 'u';
                else if ( !strcasecmp(optarg,"lc") ) args->mark_ins = 'l';
                else error("The argument is not recognised: --mark-ins %s\n",optarg);
                break;
            case  3 :
                if ( !strcasecmp(optarg,"uc") ) args->mark_snv = 'u';
                else if ( !strcasecmp(optarg,"lc") ) args->mark_snv = 'l';
                else error("The argument is not recognised: --mark-snv %s\n",optarg);
                break;
            case 'p': args->chr_prefix = optarg; break;
            case 's': args->sample = optarg; break;
            case 'S': args->sample_fname = optarg; break;
            case 'o': args->output_fname = optarg; break;
            case 'I': args->output_iupac = 1; break;
            case 'e':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i':
                if ( args->filter_str ) error("Error: only one -i or -e expression can be given, and they cannot be combined\n");
                args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;
            case 'f': args->ref_fname = optarg; break;
            case 'm': add_mask(args,optarg); break;
            case  4 : add_mask_with(args,optarg); break;
            case 'a':
                args->absent_allele = optarg[0];
                if ( optarg[1]!=0 ) error("Expected single character with -a, got \"%s\"\n", optarg);
                break;
            case 'M':
                args->missing_allele = optarg[0];
                if ( optarg[1]!=0 ) error("Expected single character with -M, got \"%s\"\n", optarg);
                break;
            case 'c': args->chain_fname = optarg; break;
            case 'H':
                if ( !strcasecmp(optarg,"R") ) args->allele |= PICK_REF;
                else if ( !strcasecmp(optarg,"A") ) args->allele |= PICK_ALT;
                else if ( !strcasecmp(optarg,"L") ) args->allele |= PICK_LONG|PICK_REF;
                else if ( !strcasecmp(optarg,"S") ) args->allele |= PICK_SHORT|PICK_REF;
                else if ( !strcasecmp(optarg,"LR") ) args->allele |= PICK_LONG|PICK_REF;
                else if ( !strcasecmp(optarg,"LA") ) args->allele |= PICK_LONG|PICK_ALT;
                else if ( !strcasecmp(optarg,"SR") ) args->allele |= PICK_SHORT|PICK_REF;
                else if ( !strcasecmp(optarg,"SA") ) args->allele |= PICK_SHORT|PICK_ALT;
                else if ( !strcasecmp(optarg,"I") ) args->allele |= PICK_IUPAC;
                else if ( !strcasecmp(optarg,"1pIu") ) args->allele |= PICK_IUPAC, args->haplotype = 1;
                else if ( !strcasecmp(optarg,"2pIu") ) args->allele |= PICK_IUPAC, args->haplotype = 2;
                else
                {
                    char *tmp;
                    args->haplotype = strtol(optarg, &tmp, 10);
                    if ( tmp==optarg || *tmp ) error("Error: Could not parse --haplotype %s, expected numeric argument\n", optarg);
                    if ( args->haplotype <=0 ) error("Error: Expected positive integer with --haplotype\n");
                }
                break;
            default: usage(args); break;
        }
    }
    if ( optind>=argc ) usage(args);
    args->fname = argv[optind];

    if ( !args->ref_fname && !isatty(fileno((FILE *)stdin)) ) args->ref_fname = "-";
    if ( !args->ref_fname ) usage(args);

    init_data(args);
    consensus(args);
    destroy_data(args);
    free(args);

    return 0;
}


