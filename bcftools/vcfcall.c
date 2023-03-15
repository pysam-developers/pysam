/*  vcfcall.c -- SNP/indel variant calling from VCF/BCF.

    Copyright (C) 2013-2022 Genome Research Ltd.

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
THE SOFTWARE.  */

#include <stdarg.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <errno.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <htslib/vcf.h>
#include <time.h>
#include <stdarg.h>
#include <htslib/kfunc.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/khash_str2int.h>
#include <ctype.h>
#include "bcftools.h"
#include "call.h"
#include "prob1.h"
#include "ploidy.h"
#include "gvcf.h"
#include "regidx.h"
#include "vcfbuf.h"

void error(const char *format, ...);

#define CF_NO_GENO      1
#define CF_INS_MISSED   (1<<1)
#define CF_CCALL        (1<<2)
//                      (1<<3)
//                      (1<<4)
//                      (1<<5)
#define CF_ACGT_ONLY    (1<<6)
#define CF_QCALL        (1<<7)
#define CF_ADJLD        (1<<8)
#define CF_NO_INDEL     (1<<9)
#define CF_ANNO_MAX     (1<<10)
#define CF_MCALL        (1<<11)
#define CF_PAIRCALL     (1<<12)
#define CF_QCNT         (1<<13)
#define CF_INDEL_ONLY   (1<<14)

typedef struct
{
    tgt_als_t *als;
    int nmatch_als, ibuf;
}
rec_tgt_t;

typedef struct
{
    int flag;   // combination of CF_* flags above
    int output_type, n_threads, record_cmd_line, clevel;
    htsFile *bcf_in, *out_fh;
    char *bcf_fname, *output_fname;
    char **samples;             // for subsampling and ploidy
    int nsamples, *samples_map; // mapping from output sample names to original VCF
    char *regions, *targets;    // regions to process
    int regions_is_file, targets_is_file, regions_overlap;
    regidx_t *tgt_idx;
    regitr_t *tgt_itr, *tgt_itr_prev, *tgt_itr_tmp;
    vcfbuf_t *vcfbuf;

    char *samples_fname;
    int samples_is_file;
    int *sample2sex;    // mapping for ploidy. If negative, interpreted as -1*ploidy
    int *sex2ploidy, *sex2ploidy_prev, nsex;
    ploidy_t *ploidy;
    gvcf_t *gvcf;

    bcf1_t *missed_line;
    call_t aux;     // parameters and temporary data
    kstring_t str;

    int argc;
    char **argv;

    //  int flag, prior_type, n1, n_sub, *sublist, n_perm;
    //  uint32_t *trio_aux;
    //  char *prior_file, **subsam;
    //  uint8_t *ploidy;
    //  double theta, pref, indel_frac, min_smpl_frac, min_lrt;
    // Permutation tests
    //  int n_perm, *seeds;
    //  double min_perm_p;
    //  void *bed;
}
args_t;

static char **add_sample(void *name2idx, char **lines, int *nlines, int *mlines, char *name, char sex, int *ith)
{
    int ret = khash_str2int_get(name2idx, name, ith);
    if ( ret==0 ) return lines;

    hts_expand(char*,(*nlines+1),*mlines,lines);
    int len = strlen(name);
    lines[*nlines] = (char*) malloc(len+3);
    memcpy(lines[*nlines],name,len);
    lines[*nlines][len]   = ' ';
    lines[*nlines][len+1] = sex;
    lines[*nlines][len+2] = 0;
    *ith = *nlines;
    (*nlines)++;
    khash_str2int_set(name2idx, strdup(name), *ith);
    return lines;
}

typedef struct
{
    const char *alias, *about, *ploidy;
}
ploidy_predef_t;

static ploidy_predef_t ploidy_predefs[] =
{
    { .alias  = "GRCh37",
      .about  = "Human Genome reference assembly GRCh37 / hg19",
      .ploidy =
          "X 1 60000 M 1\n"
          "X 2699521 154931043 M 1\n"
          "Y 1 59373566 M 1\n"
          "Y 1 59373566 F 0\n"
          "MT 1 16569 M 1\n"
          "MT 1 16569 F 1\n"
          "chrX 1 60000 M 1\n"
          "chrX 2699521 154931043 M 1\n"
          "chrY 1 59373566 M 1\n"
          "chrY 1 59373566 F 0\n"
          "chrM 1 16569 M 1\n"
          "chrM 1 16569 F 1\n"
          "*  * *     M 2\n"
          "*  * *     F 2\n"
    },
    { .alias  = "GRCh38",
      .about  = "Human Genome reference assembly GRCh38 / hg38",
      .ploidy =
          "X 1 9999 M 1\n"
          "X 2781480 155701381 M 1\n"
          "Y 1 57227415 M 1\n"
          "Y 1 57227415 F 0\n"
          "MT 1 16569 M 1\n"
          "MT 1 16569 F 1\n"
          "chrX 1 9999 M 1\n"
          "chrX 2781480 155701381 M 1\n"
          "chrY 1 57227415 M 1\n"
          "chrY 1 57227415 F 0\n"
          "chrM 1 16569 M 1\n"
          "chrM 1 16569 F 1\n"
          "*  * *     M 2\n"
          "*  * *     F 2\n"
    },
    { .alias  = "X",
      .about  = "Treat male samples as haploid and female as diploid regardless of the chromosome name",
      .ploidy =
          "*  * *     M 1\n"
          "*  * *     F 2\n"
    },
    { .alias  = "Y",
      .about  = "Treat male samples as haploid and female as no-copy, regardless of the chromosome name",
      .ploidy =
          "*  * *     M 1\n"
          "*  * *     F 0\n"
    },
    { .alias  = "1",
      .about  = "Treat all samples as haploid",
      .ploidy =
          "*  * *     * 1\n"
    },
    { .alias  = "2",
      .about  = "Treat all samples as diploid",
      .ploidy =
          "*  * *     * 2\n"
    },
    {
        .alias  = NULL,
        .about  = NULL,
        .ploidy = NULL,
    }
};

// only 5 columns are required and the first is ignored:
//  ignored,sample,father(or 0),mother(or 0),sex(1=M,2=F)
static char **parse_ped_samples(call_t *call, char **vals, int nvals, int *nsmpl)
{
    int i, j, mlines = 0, nlines = 0;
    kstring_t str = {0,0,0}, fam_str = {0,0,0};
    void *name2idx = khash_str2int_init();
    char **lines = NULL;
    for (i=0; i<nvals; i++)
    {
        str.l = 0;
        kputs(vals[i], &str);
        char *col_ends[5], *tmp = str.s;
        j = 0;
        while ( *tmp && j<5 )
        {
            if ( isspace(*tmp) )
            {
                *tmp = 0;
                ++tmp;
                while ( isspace(*tmp) ) tmp++;  // allow multiple spaces
                col_ends[j] = tmp-1;
                j++;
                continue;
            }
            tmp++;
        }
        if ( j<4 ) break;

        char sex;
        if ( col_ends[3][1]=='1' ) sex = 'M';
        else if ( col_ends[3][1]=='2' ) sex = 'F';
        else break;

        lines = add_sample(name2idx, lines, &nlines, &mlines, col_ends[0]+1, sex, &j);
        if ( strcmp(col_ends[1]+1,"0") && strcmp(col_ends[2]+1,"0") )   // father and mother
        {
            call->nfams++;
            hts_expand(family_t, call->nfams, call->mfams, call->fams);
            family_t *fam = &call->fams[call->nfams-1];
            fam_str.l = 0;
            ksprintf(&fam_str,"father=%s, mother=%s, child=%s", col_ends[1]+1,col_ends[2]+1,col_ends[0]+1);
            fam->name = strdup(fam_str.s);

            if ( !khash_str2int_has_key(name2idx, col_ends[1]+1) )
                lines = add_sample(name2idx, lines, &nlines, &mlines, col_ends[1]+1, 'M', &fam->sample[FATHER]);
            if ( !khash_str2int_has_key(name2idx, col_ends[2]+1) )
                lines = add_sample(name2idx, lines, &nlines, &mlines, col_ends[2]+1, 'F', &fam->sample[MOTHER]);

            khash_str2int_get(name2idx, col_ends[0]+1, &fam->sample[CHILD]);
            khash_str2int_get(name2idx, col_ends[1]+1, &fam->sample[FATHER]);
            khash_str2int_get(name2idx, col_ends[2]+1, &fam->sample[MOTHER]);
        }
    }
    free(str.s);
    free(fam_str.s);
    khash_str2int_destroy_free(name2idx);

    if ( i!=nvals ) return NULL; // not a ped file

    *nsmpl = nlines;
    return lines;
}


/*
 *  Reads sample names and their ploidy (optional) from a file.
 *  Alternatively, if no such file exists, the file name is interpreted
 *  as a comma-separated list of samples. When ploidy is not present,
 *  the default ploidy 2 is assumed.
 */
static void set_samples(args_t *args, const char *fn, int is_file)
{
    int i, nlines;
    char **lines = hts_readlist(fn, is_file, &nlines);
    if ( !lines ) error("Could not read the file: %s\n", fn);

    int nsmpls;
    char **smpls = parse_ped_samples(&args->aux, lines, nlines, &nsmpls);
    if ( smpls )
    {
        for (i=0; i<nlines; i++) free(lines[i]);
        free(lines);
        lines = smpls;
        nlines = nsmpls;
    }
    else if ( is_file )
        fprintf(stderr,"Note: could not parse as PED: %s\n",fn);

    args->samples_map = (int*) malloc(sizeof(int)*bcf_hdr_nsamples(args->aux.hdr)); // for subsetting
    args->sample2sex  = (int*) malloc(sizeof(int)*bcf_hdr_nsamples(args->aux.hdr));
    int dflt_sex_id = ploidy_nsex(args->ploidy) - 1;
    for (i=0; i<bcf_hdr_nsamples(args->aux.hdr); i++) args->sample2sex[i] = dflt_sex_id;

    int *old2new = (int*) malloc(sizeof(int)*bcf_hdr_nsamples(args->aux.hdr));
    for (i=0; i<bcf_hdr_nsamples(args->aux.hdr); i++) old2new[i] = -1;

    int nsmpl = 0, map_needed = 0;
    for (i=0; i<nlines; i++)
    {
        char *ss = lines[i];
        while ( *ss && isspace(*ss) ) ss++;
        if ( !*ss ) error("Could not parse: %s\n", lines[i]);
        if ( *ss=='#' ) continue;
        char *se = ss;
        while ( *se && !isspace(*se) ) se++;
        char x = *se, *xptr = se; *se = 0;

        int ismpl = bcf_hdr_id2int(args->aux.hdr, BCF_DT_SAMPLE, ss);
        if ( ismpl < 0 ) { fprintf(stderr,"Warning: No such sample in the VCF: %s\n",ss); continue; }
        if ( old2new[ismpl] != -1 ) { fprintf(stderr,"Warning: The sample is listed multiple times: %s\n",ss); continue; }

        ss = se+(x != '\0');
        while ( *ss && isspace(*ss) ) ss++;
        if ( !*ss ) ss = "2";   // default ploidy
        se = ss;
        while ( *se && !isspace(*se) ) se++;
        if ( se==ss ) { *xptr = x; error("Could not parse: \"%s\"\n", lines[i]); }

        if ( ss[1]==0 && (ss[0]=='0' || ss[0]=='1' || ss[0]=='2') )
            args->sample2sex[nsmpl] = -1*(ss[0]-'0');
        else
            args->sample2sex[nsmpl] = ploidy_add_sex(args->ploidy, ss);

        if ( ismpl!=nsmpl ) map_needed = 1;
        args->samples_map[nsmpl] = ismpl;
        old2new[ismpl] = nsmpl;
        nsmpl++;
    }

    for (i=0; i<args->aux.nfams; i++)
    {
        int j, nmiss = 0;
        family_t *fam = &args->aux.fams[i];
        for (j=0; j<3; j++)
        {
            fam->sample[i] = old2new[fam->sample[i]];
            if ( fam->sample[i]<0 ) nmiss++;
        }
        assert( nmiss==0 || nmiss==3 );
    }
    free(old2new);

    if ( !map_needed ) { free(args->samples_map); args->samples_map = NULL; }

    args->nsamples = nsmpl;
    args->samples = lines;
}

static void init_missed_line(args_t *args)
{
    int i;
    for (i=0; i<bcf_hdr_nsamples(args->aux.hdr); i++)
    {
        args->aux.gts[i*2]   = bcf_gt_missing;
        args->aux.gts[i*2+1] = bcf_int32_vector_end;
    }
    args->missed_line = bcf_init1();
    bcf_update_genotypes(args->aux.hdr, args->missed_line, args->aux.gts, 2*bcf_hdr_nsamples(args->aux.hdr));
    bcf_float_set_missing(args->missed_line->qual);
}

static int tgt_parse(const char *line, char **chr_beg, char **chr_end, uint32_t *beg, uint32_t *end, void *payload, void *usr)
{
    char *ss = (char*) line;
    while ( *ss && isspace(*ss) ) ss++;
    if ( !*ss ) { fprintf(stderr,"Could not parse the line: %s\n", line); return -2; }
    if ( *ss=='#' ) return -1;  // skip comments

    char *se = ss;
    while ( *se && !isspace(*se) ) se++;

    *chr_beg = ss;
    *chr_end = se-1;

    if ( !*se ) { fprintf(stderr,"Could not parse the line: %s\n", line); return -2; }

    ss = se+1;
    *beg = strtod(ss, &se);
    if ( ss==se ) { fprintf(stderr,"Could not parse tab line: %s\n", line); return -2; }
    if ( *beg==0 ) { fprintf(stderr,"Could not parse tab line, expected 1-based coordinate: %s\n", line); return -2; }
    (*beg)--;
    *end = *beg;

    if ( !usr ) return 0; // allele information not required

    ss = se+1;
    tgt_als_t *als = (tgt_als_t*)payload;
    als->used   = 0;
    als->n      = 0;
    als->allele = NULL;
    while ( *ss )
    {
        se = ss;
        while ( *se && *se!=',' ) se++;
        als->n++;
        als->allele = (char**)realloc(als->allele,als->n*sizeof(*als->allele));
        als->allele[als->n-1] = (char*)malloc(se-ss+1);
        memcpy(als->allele[als->n-1],ss,se-ss);
        als->allele[als->n-1][se-ss] = 0;
        ss = se+1;
        if ( !*se ) break;
    }
    return 0;
}
static void tgt_free(void *payload)
{
    tgt_als_t *als = (tgt_als_t*)payload;
    int i;
    for (i=0; i<als->n; i++) free(als->allele[i]);
    free(als->allele);
}
static void tgt_flush_region(args_t *args, char *chr, uint32_t beg, uint32_t end)
{
    if ( !regidx_overlap(args->tgt_idx, chr,beg,end,args->tgt_itr_tmp) ) return;
    while ( regitr_overlap(args->tgt_itr_tmp) )
    {
        if ( args->tgt_itr_tmp->beg < beg ) continue;

        tgt_als_t *tgt_als = &regitr_payload(args->tgt_itr_tmp,tgt_als_t);
        if ( tgt_als->used ) continue;

        args->missed_line->rid  = bcf_hdr_name2id(args->aux.hdr,chr);
        args->missed_line->pos  = args->tgt_itr_tmp->beg;
        bcf_unpack(args->missed_line,BCF_UN_ALL);
        bcf_update_alleles(args->aux.hdr, args->missed_line, (const char**)tgt_als->allele, tgt_als->n);
        tgt_als->used = 1;
        if ( bcf_write1(args->out_fh, args->aux.hdr, args->missed_line)!=0 ) error("[%s] Error: failed to write to %s\n", __func__,args->output_fname);
    }
}
static void tgt_flush(args_t *args, bcf1_t *rec)
{
    if ( rec )
    {
        char *chr = (char*)bcf_seqname(args->aux.hdr,rec);

        if ( !args->tgt_itr_prev )                  // first record
            tgt_flush_region(args,chr,0,rec->pos-1);

        else if ( strcmp(chr,args->tgt_itr_prev->seq) )  // first record on a new chromosome
        {
            tgt_flush_region(args,args->tgt_itr_prev->seq,args->tgt_itr_prev->beg+1,REGIDX_MAX);
            tgt_flush_region(args,chr,0,rec->pos-1);
        }
        else                                        // another record on the same chromosome
            tgt_flush_region(args,args->tgt_itr_prev->seq,args->tgt_itr_prev->beg,rec->pos-1);
    }
    else
    {
        // flush everything
        if ( args->tgt_itr_prev )
            tgt_flush_region(args,args->tgt_itr_prev->seq,args->tgt_itr_prev->beg,REGIDX_MAX);

        int i, nchr = 0;
        char **chr = regidx_seq_names(args->tgt_idx, &nchr);
        for (i=0; i<nchr; i++)
            tgt_flush_region(args,chr[i],0,REGIDX_MAX);
    }
}
inline static int is_indel(int nals, char **als)
{
    // This is mpileup output, we can make some assumption:
    //  - no MNPs
    //  - "<*>" is not present at indels sites and there are no other symbolic alleles than <*>
    if ( als[1][0]=='<' ) return 0;

    int i;
    for (i=0; i<nals; i++)
    {
        if ( als[i][0]=='<' ) continue;
        if ( als[i][1] ) return 1;
    }
    return 0;
}
bcf1_t *next_line(args_t *args)
{
    bcf1_t *rec = NULL;
    if ( !args->vcfbuf )
    {
        while ( bcf_sr_next_line(args->aux.srs) )
        {
            rec = args->aux.srs->readers[0].buffer[0];
            if ( args->aux.srs->errnum || rec->errcode ) error("Error: could not parse the input VCF\n");
            if ( args->tgt_idx )
            {
                if ( !regidx_overlap(args->tgt_idx, bcf_seqname(args->aux.hdr,rec),rec->pos,rec->pos,args->tgt_itr) ) continue;

                // For backward compatibility: require the exact position, not an interval overlap
                int pos_match = 0;
                while ( regitr_overlap(args->tgt_itr) )
                {
                    if ( args->tgt_itr->beg != rec->pos ) continue;
                    pos_match = 1;
                    break;
                }
                if ( !pos_match ) continue;
            }
            if ( args->samples_map ) bcf_subset(args->aux.hdr, rec, args->nsamples, args->samples_map);
            bcf_unpack(rec, BCF_UN_STR);
            return rec;
        }
        return NULL;
    }

    // If we are here,-C alleles was given and vcfbuf and tgt_idx are set

    // Fill the buffer with duplicate lines
    int vcfbuf_full = 1;
    int nbuf = vcfbuf_nsites(args->vcfbuf);
    bcf1_t *rec0 = NULL, *recN = NULL;
    if ( nbuf==0 ) vcfbuf_full = 0;
    else if ( nbuf==1 )
    {
        vcfbuf_full = 0;
        rec0 = vcfbuf_peek(args->vcfbuf, 0);
    }
    else
    {
        rec0 = vcfbuf_peek(args->vcfbuf, 0);
        recN = vcfbuf_peek(args->vcfbuf, nbuf-1);
        if ( rec0->rid == recN->rid && rec0->pos == recN->pos ) vcfbuf_full = 0;
    }
    if ( !vcfbuf_full )
    {
        while ( bcf_sr_next_line(args->aux.srs) )
        {
            rec = args->aux.srs->readers[0].buffer[0];
            if ( args->aux.srs->errnum || rec->errcode ) error("Error: could not parse the input VCF\n");
            if ( !regidx_overlap(args->tgt_idx, bcf_seqname(args->aux.hdr,rec),rec->pos,rec->pos,args->tgt_itr) ) continue;
            // as above: require the exact position, not an interval overlap
            int exact_match = 0;
            while ( regitr_overlap(args->tgt_itr) )
            {
                if ( args->tgt_itr->beg != rec->pos ) continue;
                exact_match = 1;
                break;
            }
            if ( !exact_match ) continue;

            if ( args->samples_map ) bcf_subset(args->aux.hdr, rec, args->nsamples, args->samples_map);
            bcf_unpack(rec, BCF_UN_STR);
            if ( !rec0 ) rec0 = rec;
            recN = rec;
            args->aux.srs->readers[0].buffer[0] = vcfbuf_push(args->vcfbuf, rec);
            if ( rec0->rid!=recN->rid || rec0->pos!=recN->pos ) break;
        }
    }

    nbuf = vcfbuf_nsites(args->vcfbuf);
    int n, i,j;
    for (n=nbuf; n>1; n--)
    {
        recN = vcfbuf_peek(args->vcfbuf, n-1);
        if ( rec0->rid==recN->rid && rec0->pos==recN->pos ) break;
    }
    if ( n==0 )
    {
        assert( !nbuf );
        return NULL;
    }

    // Find the VCF and tab record with the best matching combination of alleles, prioritize
    // records of the same type (snp vs indel)
    rec_tgt_t rec_tgt;
    memset(&rec_tgt,0,sizeof(rec_tgt));
    regidx_overlap(args->tgt_idx, bcf_seqname(args->aux.hdr,rec0),rec0->pos,rec0->pos,args->tgt_itr);
    regitr_t *tmp_itr = regitr_init(args->tgt_idx);
    for (i=0; i<n; i++)
    {
        regitr_copy(tmp_itr, args->tgt_itr);
        rec = vcfbuf_peek(args->vcfbuf, i);
        int rec_indel = is_indel(rec->n_allele, rec->d.allele) ? 1 : -1;
        while ( regitr_overlap(tmp_itr) )
        {
            if ( tmp_itr->beg != rec->pos ) continue;
            tgt_als_t *als = &regitr_payload(tmp_itr,tgt_als_t);
            if ( als->used ) continue;
            int nmatch_als = 0;
            vcmp_t *vcmp = vcmp_init();
            int ret = vcmp_set_ref(vcmp, rec->d.allele[0], als->allele[0]);
            if ( ret==0 )
            {
                nmatch_als++;
                if ( rec->n_allele > 1 && als->n > 1 )
                {
                    for (j=1; j<als->n; j++)
                    {
                        if ( vcmp_find_allele(vcmp, rec->d.allele+1, rec->n_allele-1, als->allele[j])>=0 ) nmatch_als++;
                    }
                }
            }
            int als_indel = is_indel(als->n, als->allele) ? 1 : -1;
            nmatch_als *= rec_indel*als_indel;
            if ( nmatch_als > rec_tgt.nmatch_als || !rec_tgt.als )
            {
                rec_tgt.nmatch_als = nmatch_als;
                rec_tgt.als  = als;
                rec_tgt.ibuf = i;
            }
            vcmp_destroy(vcmp);
        }
    }
    regitr_destroy(tmp_itr);

    args->aux.tgt_als = rec_tgt.als;
    if ( rec_tgt.als ) rec_tgt.als->used = 1;

    rec = vcfbuf_remove(args->vcfbuf, rec_tgt.ibuf);
    return rec;
}

static void init_data(args_t *args)
{
    args->aux.srs = bcf_sr_init();

    // Open files for input and output, initialize structures
    if ( args->targets )
    {
        args->tgt_idx = regidx_init(args->targets, tgt_parse, args->aux.flag&CALL_CONSTR_ALLELES ? tgt_free : (regidx_free_f) NULL, sizeof(tgt_als_t), args->aux.flag&CALL_CONSTR_ALLELES ? args : NULL);
        args->tgt_itr = regitr_init(args->tgt_idx);
        args->tgt_itr_tmp = regitr_init(args->tgt_idx);
    }

    if ( args->regions )
    {
        bcf_sr_set_opt(args->aux.srs,BCF_SR_REGIONS_OVERLAP,args->regions_overlap);
        if ( bcf_sr_set_regions(args->aux.srs, args->regions, args->regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions);
    }

    if ( !bcf_sr_add_reader(args->aux.srs, args->bcf_fname) )
        error("Failed to read from %s: %s\n", !strcmp("-",args->bcf_fname)?"standard input":args->bcf_fname,bcf_sr_strerror(args->aux.srs->errnum));
    args->aux.hdr = bcf_sr_get_header(args->aux.srs,0);

    int i;
    if ( args->samples_fname )
    {
        set_samples(args, args->samples_fname, args->samples_is_file);
        if ( args->aux.flag&CALL_CONSTR_TRIO )
        {
            if ( 3*args->aux.nfams!=args->nsamples ) error("Expected only trios in %s, sorry!\n", args->samples_fname);
            fprintf(stderr,"Detected %d samples in %d trio families\n", args->nsamples,args->aux.nfams);
        }
    }
    if ( args->ploidy  )
    {
        args->nsex = ploidy_nsex(args->ploidy);
        args->sex2ploidy = (int*) calloc(args->nsex,sizeof(int));
        args->sex2ploidy_prev = (int*) calloc(args->nsex,sizeof(int));
        if ( !args->nsamples )
        {
            args->nsamples = bcf_hdr_nsamples(args->aux.hdr);
            args->sample2sex = (int*) malloc(sizeof(int)*args->nsamples);
            for (i=0; i<args->nsamples; i++) args->sample2sex[i] = args->nsex - 1;
        }
    }
    if ( args->nsamples )
    {
        args->aux.ploidy = (uint8_t*) malloc(args->nsamples);
        for (i=0; i<args->nsamples; i++) args->aux.ploidy[i] = ploidy_max(args->ploidy);
        for (i=0; i<args->nsex; i++) args->sex2ploidy_prev[i] = ploidy_max(args->ploidy);
        for (i=0; i<args->nsamples; i++)
            if ( args->sample2sex[i] >= args->nsex ) args->sample2sex[i] = args->nsex - 1;
    }

    if ( args->gvcf )
    {
        int id = bcf_hdr_id2int(args->aux.hdr,BCF_DT_ID,"DP");
        if ( id<0 || !bcf_hdr_idinfo_exists(args->aux.hdr,BCF_HL_FMT,id) ) error("--gvcf output mode requires FORMAT/DP tag, which is not present in the input header\n");
        gvcf_update_header(args->gvcf, args->aux.hdr);
    }

    if ( args->samples_map )
    {
        args->aux.hdr = bcf_hdr_subset(bcf_sr_get_header(args->aux.srs,0), args->nsamples, args->samples, args->samples_map);
        if ( !args->aux.hdr ) error("Error occurred while subsetting samples\n");
        for (i=0; i<args->nsamples; i++)
            if ( args->samples_map[i]<0 ) error("No such sample: %s\n", args->samples[i]);
        if ( !bcf_hdr_nsamples(args->aux.hdr) ) error("No matching sample found\n");
    }
    else
    {
        args->aux.hdr = bcf_hdr_dup(bcf_sr_get_header(args->aux.srs,0));
        if ( args->samples )
        {
            for (i=0; i<args->nsamples; i++)
                if ( bcf_hdr_id2int(args->aux.hdr,BCF_DT_SAMPLE,args->samples[i])<0 )
                    error("No such sample: %s\n", args->samples[i]);
        }
    }

    if ( args->aux.flag & CALL_CONSTR_ALLELES )
        args->vcfbuf = vcfbuf_init(args->aux.hdr, 0);

    char wmode[8];
    set_wmode(wmode,args->output_type,args->output_fname,args->clevel);
    args->out_fh = hts_open(args->output_fname ? args->output_fname : "-", wmode);
    if ( args->out_fh == NULL ) error("Error: cannot write to \"%s\": %s\n", args->output_fname, strerror(errno));
    if ( args->n_threads ) hts_set_threads(args->out_fh, args->n_threads);

    if ( args->flag & CF_QCALL )
        return;

    if ( args->flag & CF_MCALL )
        mcall_init(&args->aux);

    if ( args->flag & CF_CCALL )
        ccall_init(&args->aux);

    bcf_hdr_remove(args->aux.hdr, BCF_HL_INFO, "QS");
    bcf_hdr_remove(args->aux.hdr, BCF_HL_INFO, "I16");

    if (args->record_cmd_line) bcf_hdr_append_version(args->aux.hdr, args->argc, args->argv, "bcftools_call");
    if ( bcf_hdr_write(args->out_fh, args->aux.hdr)!=0 ) error("[%s] Error: cannot write the header to %s\n", __func__,args->output_fname);

    if ( args->flag&CF_INS_MISSED ) init_missed_line(args);
}

static void destroy_data(args_t *args)
{
    if ( args->vcfbuf ) vcfbuf_destroy(args->vcfbuf);
    if ( args->tgt_idx )
    {
        regidx_destroy(args->tgt_idx);
        regitr_destroy(args->tgt_itr);
        regitr_destroy(args->tgt_itr_tmp);
        if ( args->tgt_itr_prev ) regitr_destroy(args->tgt_itr_prev);
    }
    if ( args->flag & CF_CCALL ) ccall_destroy(&args->aux);
    else if ( args->flag & CF_MCALL ) mcall_destroy(&args->aux);
    else if ( args->flag & CF_QCALL ) qcall_destroy(&args->aux);
    int i;
    if ( args->samples )
    {
        for (i=0; i<args->nsamples; i++) free(args->samples[i]);
    }
    if ( args->aux.fams )
    {
        for (i=0; i<args->aux.nfams; i++) free(args->aux.fams[i].name);
        free(args->aux.fams);
    }
    if ( args->missed_line ) bcf_destroy(args->missed_line);
    ploidy_destroy(args->ploidy);
    free(args->sex2ploidy);
    free(args->sex2ploidy_prev);
    free(args->samples);
    free(args->samples_map);
    free(args->sample2sex);
    free(args->aux.ploidy);
    free(args->str.s);
    if ( args->gvcf ) gvcf_destroy(args->gvcf);
    bcf_hdr_destroy(args->aux.hdr);
    if ( hts_close(args->out_fh)!=0 ) error("[%s] Error: close failed .. %s\n", __func__,args->output_fname);
    bcf_sr_destroy(args->aux.srs);
}

void parse_novel_rate(args_t *args, const char *str)
{
    if ( sscanf(str,"%le,%le,%le",&args->aux.trio_Pm_SNPs,&args->aux.trio_Pm_del,&args->aux.trio_Pm_ins)==3 )  // explicit for all
    {
        args->aux.trio_Pm_SNPs = 1 - args->aux.trio_Pm_SNPs;
        args->aux.trio_Pm_del  = 1 - args->aux.trio_Pm_del;
        args->aux.trio_Pm_ins  = 1 - args->aux.trio_Pm_ins;
    }
    else if ( sscanf(str,"%le,%le",&args->aux.trio_Pm_SNPs,&args->aux.trio_Pm_del)==2 )   // dynamic for indels
    {
        args->aux.trio_Pm_SNPs = 1 - args->aux.trio_Pm_SNPs;
        args->aux.trio_Pm_ins  = -1;    // negative value for dynamic calculation
    }
    else if ( sscanf(str,"%le",&args->aux.trio_Pm_SNPs)==1 )  // same for all
    {
        args->aux.trio_Pm_SNPs = 1 - args->aux.trio_Pm_SNPs;
        args->aux.trio_Pm_del  = -1;
        args->aux.trio_Pm_ins  = -1;
    }
    else error("Could not parse --novel-rate %s\n", str);
}

static void list_annotations(FILE *fp)
{
    fprintf(fp,
        "\n"
        "Optional INFO annotations available with -m (\"INFO/\" prefix is optional):\n"
        "  INFO/PV4   .. P-values for strand bias, baseQ bias, mapQ bias and tail distance bias (Number=4,Type=Float)\n"
        "\n"
        "Optional FORMAT annotations available with -m (\"FORMAT/\" prefix is optional):\n"
        "  FORMAT/GQ  .. Phred-scaled genotype quality (Number=1,Type=Integer)\n"
        "  FORMAT/GP  .. Phred-scaled genotype posterior probabilities (Number=G,Type=Float)\n"
        "\n");
}

static int parse_output_tags(const char *str)
{
    int flag = 0;
    const char *ss = str;
    while ( *ss )
    {
        const char *se = ss;
        while ( *se && *se!=',' ) se++;
        if ( !strncasecmp(ss,"GQ",se-ss) || !strncasecmp(ss,"FORMAT/GQ",se-ss) || !strncasecmp(ss,"FMT/GQ",se-ss)  ) flag |= CALL_FMT_GQ;
        else if ( !strncasecmp(ss,"GP",se-ss) || !strncasecmp(ss,"FORMAT/GP",se-ss) || !strncasecmp(ss,"FMT/GP",se-ss) ) flag |= CALL_FMT_GP;
        else if ( !strncasecmp(ss,"PV4",se-ss) || !strncasecmp(ss,"INFO/PV4",se-ss) ) flag |= CALL_FMT_PV4;
        else
        {
            fprintf(stderr,"Could not parse \"%s\"\n", str);
            exit(1);
        }
        if ( !*se ) break;
        ss = se + 1;
    }
    return flag;
}

static void set_ploidy(args_t *args, bcf1_t *rec)
{
    ploidy_query(args->ploidy,(char*)bcf_seqname(args->aux.hdr,rec),rec->pos,args->sex2ploidy,NULL,NULL);

    int i;
    for (i=0; i<args->nsex; i++)
        if ( args->sex2ploidy[i]!=args->sex2ploidy_prev[i] ) break;

    if ( i==args->nsex ) return;    // ploidy same as previously

    for (i=0; i<args->nsamples; i++)
    {
        if ( args->sample2sex[i]<0 )
            args->aux.ploidy[i] = -1*args->sample2sex[i];
        else
            args->aux.ploidy[i] = args->sex2ploidy[args->sample2sex[i]];
    }
    int *tmp = args->sex2ploidy; args->sex2ploidy = args->sex2ploidy_prev; args->sex2ploidy_prev = tmp;
}

ploidy_t *init_ploidy(char *alias)
{
    const ploidy_predef_t *pld = ploidy_predefs;

    int detailed = 0, len = strlen(alias);
    if ( alias[len-1]=='?' ) { detailed = 1; alias[len-1] = 0; }

    while ( pld->alias && strcasecmp(alias,pld->alias) ) pld++;

    if ( !pld->alias )
    {
        fprintf(stderr,"\nPRE-DEFINED PLOIDY FILES\n\n");
        fprintf(stderr," * Columns are: CHROM,FROM,TO,SEX,PLOIDY\n");
        fprintf(stderr," * Coordinates are 1-based inclusive.\n");
        fprintf(stderr," * A '*' means any value not otherwise defined.\n\n");
        pld = ploidy_predefs;
        while ( pld->alias )
        {
            fprintf(stderr,"%s\n   .. %s\n\n", pld->alias,pld->about);
            if ( detailed )
                fprintf(stderr,"%s\n", pld->ploidy);
            pld++;
        }
        fprintf(stderr,"Run as --ploidy <alias> (e.g. --ploidy GRCh37).\n");
        fprintf(stderr,"To see the detailed ploidy definition, append a question mark (e.g. --ploidy GRCh37?).\n");
        fprintf(stderr,"\n");
        exit(-1);
    }
    else if ( detailed )
    {
        fprintf(stderr,"%s", pld->ploidy);
        exit(-1);
    }
    return ploidy_init_string(pld->ploidy,2);
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   SNP/indel variant calling from VCF/BCF. To be used in conjunction with bcftools mpileup.\n");
    fprintf(stderr, "         This command replaces the former \"bcftools view\" caller. Some of the original\n");
    fprintf(stderr, "         functionality has been temporarily lost in the process of transition to htslib,\n");
    fprintf(stderr, "         but will be added back on popular demand. The original calling model can be\n");
    fprintf(stderr, "         invoked with the -c option.\n");
    fprintf(stderr, "Usage:   bcftools call [options] <in.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "File format options:\n");
    fprintf(stderr, "       --no-version                Do not append version and command line to the header\n");
    fprintf(stderr, "   -o, --output FILE               Write output to a file [standard output]\n");
    fprintf(stderr, "   -O, --output-type b|u|z|v       Output type: 'b' compressed BCF; 'u' uncompressed BCF; 'z' compressed VCF; 'v' uncompressed VCF [v]\n");
    fprintf(stderr, "   -O, --output-type u|b|v|z[0-9]  u/b: un/compressed BCF, v/z: un/compressed VCF, 0-9: compression level [v]\n");
    fprintf(stderr, "       --ploidy ASSEMBLY[?]        Predefined ploidy, 'list' to print available settings, append '?' for details [2]\n");
    fprintf(stderr, "       --ploidy-file FILE          Space/tab-delimited list of CHROM,FROM,TO,SEX,PLOIDY\n");
    fprintf(stderr, "   -r, --regions REGION            Restrict to comma-separated list of regions\n");
    fprintf(stderr, "   -R, --regions-file FILE         Restrict to regions listed in a file\n");
    fprintf(stderr, "       --regions-overlap 0|1|2     Include if POS in the region (0), record overlaps (1), variant overlaps (2) [1]\n");
    fprintf(stderr, "   -s, --samples LIST              List of samples to include [all samples]\n");
    fprintf(stderr, "   -S, --samples-file FILE         PED file or a file with an optional column with sex (see man page for details) [all samples]\n");
    fprintf(stderr, "   -t, --targets REGION            Similar to -r but streams rather than index-jumps\n");
    fprintf(stderr, "   -T, --targets-file FILE         Similar to -R but streams rather than index-jumps\n");
    fprintf(stderr, "       --threads INT               Use multithreading with INT worker threads [0]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Input/output options:\n");
    fprintf(stderr, "   -A, --keep-alts                 Keep all possible alternate alleles at variant sites\n");
    fprintf(stderr, "   -a, --annotate LIST             Optional tags to output (lowercase allowed); '?' to list available tags\n");
    fprintf(stderr, "   -F, --prior-freqs AN,AC         Use prior allele frequencies, determined from these pre-filled tags\n");
    fprintf(stderr, "   -G, --group-samples FILE|-      Group samples by population (file with \"sample\\tgroup\") or \"-\" for single-sample calling.\n");
    fprintf(stderr, "                                   This requires FORMAT/QS or other Number=R,Type=Integer tag such as FORMAT/AD\n");
    fprintf(stderr, "       --group-samples-tag TAG     The tag to use with -G, by default FORMAT/QS and FORMAT/AD are checked automatically\n");
    fprintf(stderr, "   -g, --gvcf INT,[...]            Group non-variant sites into gVCF blocks by minimum per-sample DP\n");
    fprintf(stderr, "   -i, --insert-missed             Output also sites missed by mpileup but present in -T\n");
    fprintf(stderr, "   -M, --keep-masked-ref           Keep sites with masked reference allele (REF=N)\n");
    fprintf(stderr, "   -V, --skip-variants TYPE        Skip indels/snps\n");
    fprintf(stderr, "   -v, --variants-only             Output variant sites only\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Consensus/variant calling options:\n");
    fprintf(stderr, "   -c, --consensus-caller          The original calling method (conflicts with -m)\n");
    fprintf(stderr, "   -C, --constrain STR             One of: alleles, trio (see manual)\n");
    fprintf(stderr, "   -m, --multiallelic-caller       Alternative model for multiallelic and rare-variant calling (conflicts with -c)\n");
    fprintf(stderr, "   -n, --novel-rate FLOAT,[...]    Likelihood of novel mutation for constrained trio calling, see man page for details [1e-8,1e-9,1e-9]\n");
    fprintf(stderr, "   -p, --pval-threshold FLOAT      Variant if P(ref|D)<FLOAT with -c [0.5]\n");
    fprintf(stderr, "   -P, --prior FLOAT               Mutation rate (use bigger for greater sensitivity), use with -m [1.1e-3]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Example:\n");
    fprintf(stderr, "   # See also http://samtools.github.io/bcftools/howtos/variant-calling.html\n");
    fprintf(stderr, "   bcftools mpileup -Ou -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf\n");

    // todo (and more)
    // fprintf(stderr, "\nContrast calling and association test options:\n");
    // fprintf(stderr, "       -1 INT    number of group-1 samples [0]\n");
    // fprintf(stderr, "       -C FLOAT  posterior constrast for LRT<FLOAT and P(ref|D)<0.5 [%g]\n", args->aux.min_lrt);
    // fprintf(stderr, "       -U INT    number of permutations for association testing (effective with -1) [0]\n");
    // fprintf(stderr, "       -X FLOAT  only perform permutations for P(chi^2)<FLOAT [%g]\n", args->aux.min_perm_p);
    fprintf(stderr, "\n");
    exit(-1);
}

int main_vcfcall(int argc, char *argv[])
{
    char *ploidy_fname = NULL, *ploidy = NULL;
    args_t args;
    memset(&args, 0, sizeof(args_t));
    args.argc = argc; args.argv = argv;
    args.aux.prior_type = -1;
    args.aux.indel_frac = -1;
    args.aux.theta      = 1.1e-3;
    args.aux.pref       = 0.5;
    args.aux.min_perm_p = 0.01;
    args.aux.min_lrt    = 1;
    args.flag           = CF_ACGT_ONLY;
    args.output_fname   = "-";
    args.output_type    = FT_VCF;
    args.n_threads = 0;
    args.record_cmd_line = 1;
    args.aux.trio_Pm_SNPs = 1 - 1e-8;
    args.aux.trio_Pm_ins  = args.aux.trio_Pm_del  = 1 - 1e-9;
    args.regions_overlap = 1;
    args.clevel = -1;

    int c;
    static struct option loptions[] =
    {
        {"help",no_argument,NULL,'h'},
        {"format-fields",required_argument,NULL,'f'},
        {"annotate",required_argument,NULL,'a'},
        {"prior-freqs",required_argument,NULL,'F'},
        {"gvcf",required_argument,NULL,'g'},
        {"group-samples",required_argument,NULL,'G'},
        {"group-samples-tag",required_argument,NULL,3},
        {"output",required_argument,NULL,'o'},
        {"output-type",required_argument,NULL,'O'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"regions-overlap",required_argument,NULL,4},
        {"samples",required_argument,NULL,'s'},
        {"samples-file",required_argument,NULL,'S'},
        {"targets",required_argument,NULL,'t'},
        {"targets-file",required_argument,NULL,'T'},
        {"threads",required_argument,NULL,9},
        {"keep-alts",no_argument,NULL,'A'},
        {"insert-missed",no_argument,NULL,'i'},
        {"skip-Ns",no_argument,NULL,'N'},            // now the new default
        {"keep-masked-refs",no_argument,NULL,'M'},
        {"skip-variants",required_argument,NULL,'V'},
        {"variants-only",no_argument,NULL,'v'},
        {"consensus-caller",no_argument,NULL,'c'},
        {"constrain",required_argument,NULL,'C'},
        {"multiallelic-caller",no_argument,NULL,'m'},
        {"pval-threshold",required_argument,NULL,'p'},
        {"prior",required_argument,NULL,'P'},
        {"novel-rate",required_argument,NULL,'n'},
        {"ploidy",required_argument,NULL,1},
        {"ploidy-file",required_argument,NULL,2},
        {"chromosome-X",no_argument,NULL,'X'},
        {"chromosome-Y",no_argument,NULL,'Y'},
        {"no-version",no_argument,NULL,8},
        {NULL,0,NULL,0}
    };

    char *tmp = NULL;
    while ((c = getopt_long(argc, argv, "h?o:O:r:R:s:S:t:T:ANMV:vcmp:C:n:P:f:a:ig:XYF:G:", loptions, NULL)) >= 0)
    {
        switch (c)
        {
            case  2 : ploidy_fname = optarg; break;
            case  1 : ploidy = optarg; break;
            case 'X': ploidy = "X"; fprintf(stderr,"Warning: -X will be deprecated, please use --ploidy instead.\n"); break;
            case 'Y': ploidy = "Y"; fprintf(stderr,"Warning: -Y will be deprecated, please use --ploidy instead.\n"); break;
            case 'G': args.aux.sample_groups = optarg; break;
            case  3 : args.aux.sample_groups_tag = optarg; break;
            case 'f': fprintf(stderr,"Warning: -f, --format-fields will be deprecated, please use -a, --annotate instead.\n");
            case 'a':
                      if (optarg[0]=='?') { list_annotations(stderr); return 1; }
                      args.aux.output_tags |= parse_output_tags(optarg);
                      break;
            case 'M': args.flag &= ~CF_ACGT_ONLY; break;     // keep sites where REF is N
            case 'N': args.flag |= CF_ACGT_ONLY; break;      // omit sites where first base in REF is N (the new default)
            case 'A': args.aux.flag |= CALL_KEEPALT; break;
            case 'c': args.flag |= CF_CCALL; break;          // the original EM based calling method
            case 'i': args.flag |= CF_INS_MISSED; break;
            case 'v': args.aux.flag |= CALL_VARONLY; break;
            case 'F':
                args.aux.prior_AN = optarg;
                args.aux.prior_AC = strchr(optarg,',');
                if ( !args.aux.prior_AC ) error("Expected two tags with -F (e.g. AN,AC), got \"%s\"\n",optarg);
                *args.aux.prior_AC = 0;
                args.aux.prior_AC++;
                break;
            case 'g':
                args.gvcf = gvcf_init(optarg);
                if ( !args.gvcf ) error("Could not parse: --gvcf %s\n", optarg);
                break;
            case 'o': args.output_fname = optarg; break;
            case 'O':
                      switch (optarg[0]) {
                          case 'b': args.output_type = FT_BCF_GZ; break;
                          case 'u': args.output_type = FT_BCF; break;
                          case 'z': args.output_type = FT_VCF_GZ; break;
                          case 'v': args.output_type = FT_VCF; break;
                          default:
                          {
                              args.clevel = strtol(optarg,&tmp,10);
                              if ( *tmp || args.clevel<0 || args.clevel>9 ) error("The output type \"%s\" not recognised\n", optarg);
                          }
                      }
                      if ( optarg[1] )
                      {
                          args.clevel = strtol(optarg+1,&tmp,10);
                          if ( *tmp || args.clevel<0 || args.clevel>9 ) error("Could not parse argument: --compression-level %s\n", optarg+1);
                      }
                      break;
            case 'C':
                      if ( !strcasecmp(optarg,"alleles") ) args.aux.flag |= CALL_CONSTR_ALLELES;
                      else if ( !strcasecmp(optarg,"trio") ) args.aux.flag |= CALL_CONSTR_TRIO;
                      else error("Unknown argument to -C: \"%s\"\n", optarg);
                      break;
            case 'V':
                      if ( !strcasecmp(optarg,"snps") ) args.flag |= CF_INDEL_ONLY;
                      else if ( !strcasecmp(optarg,"indels") ) args.flag |= CF_NO_INDEL;
                      else error("Unknown skip category \"%s\" (-S argument must be \"snps\" or \"indels\")\n", optarg);
                      break;
            case 'm': args.flag |= CF_MCALL; break;         // multiallelic calling method
            case 'p':
                args.aux.pref = strtod(optarg,&tmp);
                if ( *tmp ) error("Could not parse: --pval-threshold %s\n", optarg);
                break;
            case 'P': args.aux.theta = strtod(optarg,&tmp);
                      if ( *tmp ) error("Could not parse, expected float argument: -P %s\n", optarg);
                      break;
            case 'n': parse_novel_rate(&args,optarg); break;
            case 'r': args.regions = optarg; break;
            case 'R': args.regions = optarg; args.regions_is_file = 1; break;
            case 't': args.targets = optarg; break;
            case 'T': args.targets = optarg; args.targets_is_file = 1; break;
            case 's': args.samples_fname = optarg; break;
            case 'S': args.samples_fname = optarg; args.samples_is_file = 1; break;
            case  9 : args.n_threads = strtol(optarg, 0, 0); break;
            case  8 : args.record_cmd_line = 0; break;
            case  4 :
                args.regions_overlap = parse_overlap_option(optarg);
                if ( args.regions_overlap < 0 ) error("Could not parse: --regions-overlap %s\n",optarg);
                break;
            default: usage(&args);
        }
    }
    // Sanity check options and initialize
    if ( ploidy_fname ) args.ploidy = ploidy_init(ploidy_fname, 2);
    else if ( ploidy ) args.ploidy = init_ploidy(ploidy);

    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) args.bcf_fname = "-";  // reading from stdin
        else usage(&args);
    }
    else args.bcf_fname = argv[optind++];

    if ( !ploidy_fname && !ploidy )
    {
        if ( !args.samples_is_file ) fprintf(stderr,"Note: none of --samples-file, --ploidy or --ploidy-file given, assuming all sites are diploid\n");
        args.ploidy = ploidy_init_string("* * * 0 0\n* * * 1 1\n* * * 2 2\n",2);
    }

    if ( !args.ploidy ) error("Could not initialize ploidy\n");
    if ( (args.flag & CF_CCALL ? 1 : 0) + (args.flag & CF_MCALL ? 1 : 0) + (args.flag & CF_QCALL ? 1 : 0) > 1 ) error("Only one of -c or -m options can be given\n");
    if ( !(args.flag & CF_CCALL) && !(args.flag & CF_MCALL) && !(args.flag & CF_QCALL) ) error("Expected -c or -m option\n");
    if ( (args.flag & CF_CCALL ? 1: 0) && args.gvcf ) error("gvcf -g option not functional with -c calling mode yet\n");
    if ( args.aux.n_perm && args.aux.ngrp1_samples<=0 ) error("Expected -1 with -U\n");    // not sure about this, please fix
    if ( args.aux.flag & CALL_CONSTR_ALLELES )
    {
        if ( !args.targets ) error("Expected -t or -T with \"-C alleles\"\n");
        if ( !(args.flag & CF_MCALL) ) error("The \"-C alleles\" mode requires -m\n");
    }
    if ( args.flag & CF_INS_MISSED && !(args.aux.flag&CALL_CONSTR_ALLELES) ) error("The -i option requires -C alleles\n");
    if ( args.aux.flag&CALL_VARONLY && args.gvcf ) error("The two options cannot be combined: --variants-only and --gvcf\n");
    if ( args.aux.sample_groups && !(args.flag & CF_MCALL) ) error("The -G feature is supported only with the -m calling mode\n");
    init_data(&args);

    bcf1_t *bcf_rec;
    while ( (bcf_rec = next_line(&args)) )
    {
        // Skip duplicate positions with all matching `-C alleles -T` used up
        if ( args.aux.flag&CALL_CONSTR_ALLELES && !args.aux.tgt_als ) continue;

        // Skip unwanted sites
        int i, is_indel = bcf_is_snp(bcf_rec) ? 0 : 1;
        if ( (args.flag & CF_INDEL_ONLY) && !is_indel ) continue;
        if ( (args.flag & CF_NO_INDEL) && is_indel ) continue;
        if ( (args.flag & CF_ACGT_ONLY) && (bcf_rec->d.allele[0][0]=='N' || bcf_rec->d.allele[0][0]=='n') ) continue;   // REF[0] is 'N'

        // Which allele is symbolic? All SNPs should have it, but not indels
        args.aux.unseen = 0;
        for (i=1; i<bcf_rec->n_allele; i++)
        {
            if ( bcf_rec->d.allele[i][0]=='X' ) { args.aux.unseen = i; break; }  // old X
            if ( bcf_rec->d.allele[i][0]=='<' )
            {
                if ( bcf_rec->d.allele[i][1]=='X' && bcf_rec->d.allele[i][2]=='>' ) { args.aux.unseen = i; break; } // old <X>
                if ( bcf_rec->d.allele[i][1]=='*' && bcf_rec->d.allele[i][2]=='>' ) { args.aux.unseen = i; break; } // new <*>
            }
        }
        int is_ref = (bcf_rec->n_allele==1 || (bcf_rec->n_allele==2 && args.aux.unseen>0)) ? 1 : 0;

        if ( is_ref && args.aux.flag&CALL_VARONLY )
            continue;

        bcf_unpack(bcf_rec, BCF_UN_ALL);
        if ( args.nsex ) set_ploidy(&args, bcf_rec);

        // Various output modes: QCall output (todo)
        if ( args.flag & CF_QCALL )
        {
            qcall(&args.aux, bcf_rec);
            continue;
        }

        if ( args.flag & CF_INS_MISSED )
        {
            tgt_flush(&args,bcf_rec);
            if ( !args.tgt_itr_prev ) args.tgt_itr_prev = regitr_init(args.tgt_idx);
            regitr_copy(args.tgt_itr_prev, args.tgt_itr);
        }

        // Calling modes which output VCFs
        int ret;
        if ( args.flag & CF_MCALL )
            ret = mcall(&args.aux, bcf_rec);
        else
            ret = ccall(&args.aux, bcf_rec);
        if ( ret==-1 ) error("Something is wrong\n");
        else if ( ret==-2 ) continue;   // skip the site

        // Normal output
        if ( (args.aux.flag & CALL_VARONLY) && ret==0 && !args.gvcf ) continue;     // not a variant
        if ( args.gvcf )
            bcf_rec = gvcf_write(args.gvcf, args.out_fh, args.aux.hdr, bcf_rec, ret==1?1:0);
        if ( bcf_rec && bcf_write1(args.out_fh, args.aux.hdr, bcf_rec)!=0 ) error("[%s] Error: failed to write to %s\n", __func__,args.output_fname);
    }
    if ( args.gvcf ) gvcf_write(args.gvcf, args.out_fh, args.aux.hdr, NULL, 0);
    if ( args.flag & CF_INS_MISSED ) tgt_flush(&args,NULL);
    destroy_data(&args);
    return 0;
}

