/* The MIT License

   Copyright (c) 2014 Genome Research Ltd.

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
#include <errno.h>
#include <getopt.h>
#include <unistd.h>
#include <ctype.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/kseq.h>
#include <htslib/regidx.h>
#include "bcftools.h"
#include "rbuf.h"

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


typedef struct
{
    kstring_t fa_buf;   // buffered reference sequence
    int fa_ori_pos;     // start position of the fa_buffer (wrt original sequence)
    int fa_frz_pos;     // protected position to avoid conflicting variants (last pos for SNPs/ins)
    int fa_mod_off;     // position difference of fa_frz_pos in the ori and modified sequence (ins positive)
    int fa_end_pos;     // region's end position in the original sequence
    int fa_length;      // region's length in the original sequence (in case end_pos not provided in the FASTA header)
    int fa_case;        // output upper case or lower case?
    int fa_src_pos;     // last genomic coordinate read from the input fasta (0-based)

    rbuf_t vcf_rbuf;
    bcf1_t **vcf_buf;
    int nvcf_buf, rid;

    regidx_t *mask;

    int chain_id;       // chain_id, to provide a unique ID to each chain in the chain output
    chain_t *chain;     // chain structure to store the sequence of ungapped blocks between the ref and alt sequences
                        // Note that the chain is re-initialised for each chromosome/seq_region

    bcf_srs_t *files;
    bcf_hdr_t *hdr;
    FILE *fp_out;
    FILE *fp_chain;
    char **argv;
    int argc, output_iupac, haplotype, isample;
    char *fname, *ref_fname, *sample, *output_fname, *mask_fname, *chain_fname;
}
args_t;

static chain_t* init_chain(chain_t *chain, int ref_ori_pos)
{
//     fprintf(stderr, "init_chain(*chain, ref_ori_pos=%d)\n", ref_ori_pos);
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

static void destroy_chain(args_t *args)
{
    chain_t *chain = args->chain;
    free(chain->ref_gaps);
    free(chain->alt_gaps);
    free(chain->block_lengths);
    free(chain);
    chain = NULL;
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
    fprintf(args->fp_chain, "chain %d %s %d + %d %d %s %d + %d %d %d\n", score, bcf_hdr_id2name(args->hdr,args->rid), ref_end_pos, chain->ori_pos, ref_end_pos, bcf_hdr_id2name(args->hdr,args->rid), alt_end_pos, chain->ori_pos, alt_end_pos, ++args->chain_id);
    for (n=0; n<chain->num; n++) {
        fprintf(args->fp_chain, "%d %d %d\n", chain->block_lengths[n], chain->ref_gaps[n], chain->alt_gaps[n]);
    }
    fprintf(args->fp_chain, "%d\n\n", last_block_size);
}

static void push_chain_gap(chain_t *chain, int ref_start, int ref_len, int alt_start, int alt_len)
{
//     fprintf(stderr, "push_chain_gap(*chain, ref_start=%d, ref_len=%d, alt_start=%d, alt_len=%d)\n", ref_start, ref_len, alt_start, alt_len);
    int num = chain->num;

    if (ref_start <= chain->ref_last_block_ori) {
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
    if ( !bcf_sr_add_reader(args->files,args->fname) ) error("Failed to open %s: %s\n", args->fname, bcf_sr_strerror(args->files->errnum));
    args->hdr = args->files->readers[0].header;
    args->isample = -1;
    if ( args->sample )
    {
        args->isample = bcf_hdr_id2int(args->hdr,BCF_DT_SAMPLE,args->sample);
        if ( args->isample<0 ) error("No such sample: %s\n", args->sample);
    }
    if ( args->haplotype && args->isample<0 )
    {
        if ( bcf_hdr_nsamples(args->hdr) > 1 ) error("The --sample option is expected with --haplotype\n");
        args->isample = 0;
    }
    if ( args->mask_fname )
    {
        args->mask = regidx_init(args->mask_fname,NULL,NULL,0,NULL);
        if ( !args->mask ) error("Failed to initialize mask regions\n");
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
    else args->fp_out = stdout;
}

static void destroy_data(args_t *args)
{
    bcf_sr_destroy(args->files);
    int i;
    for (i=0; i<args->vcf_rbuf.m; i++)
        if ( args->vcf_buf[i] ) bcf_destroy1(args->vcf_buf[i]);
    free(args->vcf_buf);
    free(args->fa_buf.s);
    if ( args->mask ) regidx_destroy(args->mask);
    if ( args->chain_fname )
        if ( fclose(args->fp_chain) ) error("Close failed: %s\n", args->chain_fname);
    if ( fclose(args->fp_out) ) error("Close failed: %s\n", args->output_fname);
}

static void init_region(args_t *args, char *line)
{
    char *ss, *se = line;
    while ( *se && !isspace(*se) && *se!=':' ) se++;
    int from = 0, to = 0;
    char tmp, *tmp_ptr = NULL;
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
    args->rid = bcf_hdr_name2id(args->hdr,line);
    if ( args->rid<0 ) fprintf(stderr,"Warning: Sequence \"%s\" not in %s\n", line,args->fname);
    args->fa_buf.l = 0;
    args->fa_length = 0;
    args->fa_end_pos = to;
    args->fa_ori_pos = from;
    args->fa_src_pos = from;
    args->fa_mod_off = 0;
    args->fa_frz_pos = -1;
    args->fa_case    = -1;
    args->vcf_rbuf.n = 0;
    bcf_sr_seek(args->files,line,args->fa_ori_pos);
    if ( tmp_ptr ) *tmp_ptr = tmp;
    fprintf(args->fp_out,">%s\n",line);
    if (args->chain_fname )
    {
        args->chain = init_chain(args->chain, args->fa_ori_pos);
    } else {
        args->chain = NULL;
    }
}

static bcf1_t **next_vcf_line(args_t *args)
{
    if ( args->vcf_rbuf.n )
    {
        int i = rbuf_shift(&args->vcf_rbuf);
        return &args->vcf_buf[i];
    }
    else if ( bcf_sr_next_line(args->files) )
        return &args->files->readers[0].buffer[0];

    return NULL;
}
static void unread_vcf_line(args_t *args, bcf1_t **rec_ptr)
{
    bcf1_t *rec = *rec_ptr;
    if ( args->vcf_rbuf.n >= args->vcf_rbuf.m )
        error("FIXME: too many overlapping records near %s:%d\n", bcf_seqname(args->hdr,rec),rec->pos+1);

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
static void apply_variant(args_t *args, bcf1_t *rec)
{
    if ( rec->n_allele==1 ) return;

    if ( rec->pos <= args->fa_frz_pos )
    {
        fprintf(stderr,"The site %s:%d overlaps with another variant, skipping...\n", bcf_seqname(args->hdr,rec),rec->pos+1);
        return;
    }
    if ( args->mask )
    {
        char *chr = (char*)bcf_hdr_id2name(args->hdr,args->rid);
        int start = rec->pos;
        int end   = rec->pos + rec->rlen - 1;
        if ( regidx_overlap(args->mask, chr,start,end,NULL) ) return;
    }

    int i, ialt = 1;
    if ( args->isample >= 0 )
    {
        bcf_fmt_t *fmt = bcf_get_fmt(args->hdr, rec, "GT");
        if ( !fmt ) return;
        if ( args->haplotype )
        {
            if ( args->haplotype > fmt->n ) error("Can't apply %d-th haplotype at %s:%d\n", args->haplotype,bcf_seqname(args->hdr,rec),rec->pos+1);
            uint8_t *ignore, *ptr = fmt->p + fmt->size*args->isample + args->haplotype - 1;
            ialt = bcf_dec_int1(ptr, fmt->type, &ignore);
            if ( bcf_gt_is_missing(ialt) || ialt==bcf_int32_vector_end ) return;
            ialt = bcf_gt_allele(ialt);
        }
        else if ( args->output_iupac ) 
        {
            uint8_t *ignore, *ptr = fmt->p + fmt->size*args->isample;
            ialt = bcf_dec_int1(ptr, fmt->type, &ignore);
            if ( bcf_gt_is_missing(ialt) || ialt==bcf_int32_vector_end ) return;
            ialt = bcf_gt_allele(ialt);

            int jalt;
            if ( fmt->n>1 )
            {
                ptr = fmt->p + fmt->size*args->isample + 1;
                jalt = bcf_dec_int1(ptr, fmt->type, &ignore);
                if ( bcf_gt_is_missing(jalt) || jalt==bcf_int32_vector_end ) jalt = ialt;
                else jalt = bcf_gt_allele(jalt);
            }
            else jalt = ialt;
            if ( rec->n_allele <= ialt || rec->n_allele <= jalt ) error("Broken VCF, too few alts at %s:%d\n", bcf_seqname(args->hdr,rec),rec->pos+1);
            if ( ialt!=jalt && !rec->d.allele[ialt][1] && !rec->d.allele[jalt][1] ) // is this a het snp?
            {
                char ial = rec->d.allele[ialt][0];
                char jal = rec->d.allele[jalt][0];
                rec->d.allele[ialt][0] = gt2iupac(ial,jal);
            }
        }
        else
        {
            for (i=0; i<fmt->n; i++)
            {
                uint8_t *ignore, *ptr = fmt->p + fmt->size*args->isample + i;
                ialt = bcf_dec_int1(ptr, fmt->type, &ignore);
                if ( bcf_gt_is_missing(ialt) || ialt==bcf_int32_vector_end ) return;
                ialt = bcf_gt_allele(ialt);
                if ( ialt ) break;
            }
        }
        if ( !ialt ) return;  // ref allele
        if ( rec->n_allele <= ialt ) error("Broken VCF, too few alts at %s:%d\n", bcf_seqname(args->hdr,rec),rec->pos+1);
    }
    else if ( args->output_iupac && !rec->d.allele[0][1] && !rec->d.allele[1][1] )
    {
        char ial = rec->d.allele[0][0];
        char jal = rec->d.allele[1][0];
        rec->d.allele[1][0] = gt2iupac(ial,jal);
    }

    int idx = rec->pos - args->fa_ori_pos + args->fa_mod_off;
    if ( idx<0 || idx>=args->fa_buf.l ) 
        error("FIXME: %s:%d .. idx=%d, ori_pos=%d, len=%d, off=%d\n",bcf_seqname(args->hdr,rec),rec->pos+1,idx,args->fa_ori_pos,args->fa_buf.l,args->fa_mod_off);

    // sanity check the reference base
    int len_diff = 0, alen = 0;
    if ( rec->d.allele[ialt][0]=='<' )
    {
        if ( strcasecmp(rec->d.allele[ialt], "<DEL>") )
            error("Symbolic alleles other than <DEL> are currently not supported: %s at %s:%d\n",rec->d.allele[ialt],bcf_seqname(args->hdr,rec),rec->pos+1);
        assert( rec->d.allele[0][1]==0 );           // todo: for now expecting strlen(REF) = 1
        len_diff = 1-rec->rlen;
        rec->d.allele[ialt] = rec->d.allele[0];     // according to VCF spec, REF must precede the event
        alen = strlen(rec->d.allele[ialt]);
    }
    else if ( strncasecmp(rec->d.allele[0],args->fa_buf.s+idx,rec->rlen) )
    {
        // fprintf(stderr,"%d .. [%s], idx=%d ori=%d off=%d\n",args->fa_ori_pos,args->fa_buf.s,idx,args->fa_ori_pos,args->fa_mod_off);
        char tmp = 0;
        if ( args->fa_buf.l - idx > rec->rlen ) 
        { 
            tmp = args->fa_buf.s[idx+rec->rlen];
            args->fa_buf.s[idx+rec->rlen] = 0;
        }
        error(
            "The fasta sequence does not match the REF allele at %s:%d:\n"
            "   .vcf: [%s]\n" 
            "   .vcf: [%s] <- (ALT)\n" 
            "   .fa:  [%s]%c%s\n",
            bcf_seqname(args->hdr,rec),rec->pos+1, rec->d.allele[0], rec->d.allele[ialt], args->fa_buf.s+idx, 
            tmp?tmp:' ',tmp?args->fa_buf.s+idx+rec->rlen+1:""
            );
    }
    else
    {
        alen = strlen(rec->d.allele[ialt]);
        len_diff = alen - rec->rlen;
    }

    if ( args->fa_case )
        for (i=0; i<alen; i++) rec->d.allele[ialt][i] = toupper(rec->d.allele[ialt][i]);
    else
        for (i=0; i<alen; i++) rec->d.allele[ialt][i] = tolower(rec->d.allele[ialt][i]);

    if ( len_diff <= 0 )
    {
        // deletion or same size event
        for (i=0; i<alen; i++)
            args->fa_buf.s[idx+i] = rec->d.allele[ialt][i];
        if ( len_diff )
            memmove(args->fa_buf.s+idx+alen,args->fa_buf.s+idx+rec->rlen,args->fa_buf.l-idx-rec->rlen);
    }
    else
    {
        // insertion
        ks_resize(&args->fa_buf, args->fa_buf.l + len_diff);
        memmove(args->fa_buf.s + idx + rec->rlen + len_diff, args->fa_buf.s + idx + rec->rlen, args->fa_buf.l - idx - rec->rlen);
        for (i=0; i<alen; i++)
            args->fa_buf.s[idx+i] = rec->d.allele[ialt][i];
    }
    if (args->chain && len_diff != 0)
    {
        // If first nucleotide of both REF and ALT are the same... (indels typically include the nucleotide before the variant)
        if ( strncasecmp(rec->d.allele[0],rec->d.allele[ialt],1) == 0)
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
}


static void mask_region(args_t *args, char *seq, int len)
{
    char *chr = (char*)bcf_hdr_id2name(args->hdr,args->rid);
    int start = args->fa_src_pos - len;
    int end   = args->fa_src_pos;

    regitr_t itr;
    if ( !regidx_overlap(args->mask, chr,start,end, &itr) ) return;

    int idx_start, idx_end, i;
    while ( REGITR_OVERLAP(itr,start,end) )
    {
        idx_start = REGITR_START(itr) - start;
        idx_end   = REGITR_END(itr) - start;
        if ( idx_start < 0 ) idx_start = 0;
        if ( idx_end >= len ) idx_end = len - 1;
        for (i=idx_start; i<=idx_end; i++) seq[i] = 'N';
        itr.i++;
    }
}

static void consensus(args_t *args)
{
    htsFile *fasta = hts_open(args->ref_fname, "rb");
    if ( !fasta ) error("Error reading %s\n", args->ref_fname);
    kstring_t str = {0,0,0};
    while ( hts_getline(fasta, KS_SEP_LINE, &str) > 0 )
    {
        if ( str.s[0]=='>' )
        {
            // new sequence encountered, apply all chached variants
            while ( args->vcf_rbuf.n )
            {
                if (args->chain) {
                    print_chain(args);
                    destroy_chain(args);
                }
                bcf1_t *rec = args->vcf_buf[args->vcf_rbuf.f];
                if ( rec->rid!=args->rid || ( args->fa_end_pos && rec->pos > args->fa_end_pos ) ) break;
                int i = rbuf_shift(&args->vcf_rbuf);
                apply_variant(args, args->vcf_buf[i]);
            }
            flush_fa_buffer(args, 0);
            init_region(args, str.s+1);
            continue;
        }
        args->fa_length  += str.l;
        args->fa_src_pos += str.l;

        // determine if uppercase or lowercase is used in this fasta file
        if ( args->fa_case==-1 ) args->fa_case = toupper(str.s[0])==str.s[0] ? 1 : 0;

        if ( args->mask && args->rid>=0) mask_region(args, str.s, str.l);
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
        if ( !rec_ptr ) flush_fa_buffer(args, 60);
    }
    if (args->chain) {
        print_chain(args);
        destroy_chain(args);
    }
    flush_fa_buffer(args, 0);
    hts_close(fasta);
    free(str.s);
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Create consensus sequence by applying VCF variants to a reference\n");
    fprintf(stderr, "         fasta file.\n");
    fprintf(stderr, "Usage:   bcftools consensus [OPTIONS] <file.vcf>\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -f, --fasta-ref <file>     reference sequence in fasta format\n");
    fprintf(stderr, "    -H, --haplotype <1|2>      apply variants for the given haplotype\n");
    fprintf(stderr, "    -i, --iupac-codes          output variants in the form of IUPAC ambiguity codes\n");
    fprintf(stderr, "    -m, --mask <file>          replace regions with N\n");
    fprintf(stderr, "    -o, --output <file>        write output to a file [standard output]\n");
    fprintf(stderr, "    -c, --chain <file>         write a chain file for liftover\n");
    fprintf(stderr, "    -s, --sample <name>        apply variants of the given sample\n");
    fprintf(stderr, "Examples:\n");
    fprintf(stderr, "   # Get the consensus for one region. The fasta header lines are then expected\n");
    fprintf(stderr, "   # in the form \">chr:from-to\".\n");
    fprintf(stderr, "   samtools faidx ref.fa 8:11870-11890 | bcftools consensus in.vcf.gz > out.fa\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_consensus(int argc, char *argv[])
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    args->argc   = argc; args->argv = argv;

    static struct option loptions[] = 
    {
        {"sample",1,0,'s'},
        {"iupac-codes",0,0,'i'},
        {"haplotype",1,0,'H'},
        {"output",1,0,'o'},
        {"fasta-ref",1,0,'f'},
        {"mask",1,0,'m'},
        {"chain",1,0,'c'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "h?s:1iH:f:o:m:c:",loptions,NULL)) >= 0) 
    {
        switch (c) 
        {
            case 's': args->sample = optarg; break;
            case 'o': args->output_fname = optarg; break;
            case 'i': args->output_iupac = 1; break;
            case 'f': args->ref_fname = optarg; break;
            case 'm': args->mask_fname = optarg; break;
            case 'c': args->chain_fname = optarg; break;
            case 'H': 
                args->haplotype = optarg[0] - '0'; 
                if ( args->haplotype <=0 ) error("Expected positive integer with --haplotype\n");
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


