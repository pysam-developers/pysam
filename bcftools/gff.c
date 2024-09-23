/* The MIT License

   Copyright (c) 2023 Genome Research Ltd.

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

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <inttypes.h>
#include <string.h>
#include <strings.h>
#include <htslib/hts.h>
#include <htslib/khash.h>
#include <htslib/khash_str2int.h>
#include <htslib/kseq.h>
#include <htslib/bgzf.h>
#include <errno.h>
#include "bcftools.h"
#include "regidx.h"
#include "gff.h"

/*
    Helper structures, only for initialization

    ftr_t
        temporary list of all exons, CDS, UTRs
*/
KHASH_MAP_INIT_INT(int2tscript, gf_tscript_t*)
KHASH_MAP_INIT_INT(int2gene, gf_gene_t*)
typedef struct
{
    int type;           // GF_CDS, GF_EXON, GF_5UTR, GF_3UTR
    uint32_t beg;
    uint32_t end;
    uint32_t trid;
    uint32_t strand:2;  // STRAND_{REV,FWD,UNK}
    uint32_t phase:2;   // 0, 1, 2, or 3 for unknown
    uint32_t iseq:28;
}
ftr_t;

/*
    Mapping from GFF ID string (such as ENST00000450305 or Zm00001d027230_P001)
    to integer id.  To keep the memory requirements low, the original version
    relied on IDs in the form of a string prefix and a numerical id.  However,
    it turns out that this assumption is not valid for some ensembl GFFs, see
    for example Zea_mays.AGPv4.36.gff3.gz
 */
typedef struct
{
    void *str2id;       // khash_str2int
    int nstr, mstr;
    char **str;         // numeric id to string
}
id_tbl_t;

typedef struct
{
    // all exons, CDS, UTRs
    ftr_t *ftr;
    int nftr, mftr;

    // mapping from gene id to gf_gene_t
    kh_int2gene_t *gid2gene;

    // mapping from transcript id to tscript, for quick CDS anchoring
    kh_int2tscript_t *id2tr;

    // sequences
    void *seq2int;  // str2int hash
    char **seq;
    int nseq, mseq;

    // ignored biotypes
    void *ignored_biotypes;

    id_tbl_t gene_ids;   // temporary table for mapping between gene id (eg. Zm00001d027245) and a numeric idx

    // pointers to the current partially processed line
    char *id, *id_end, *parent, *parent_end, *biotype, *biotype_end,
         *chr, *chr_end, *name, *name_end, *type, *type_end;
}
aux_t;

struct gff_t_
{
    const char *fname, *dump_fname;

    // the main regidx lookups, from chr:beg-end to overlapping features and
    // index iterator
    regidx_t *idx_cds, *idx_utr, *idx_exon, *idx_tscript;

    // temporary structures, deleted after initializtion
    aux_t init;

    // mapping between transcript id (eg. Zm00001d027245_T001) and a numeric idx
    id_tbl_t tscript_ids;

    int strip_chr_names, verbosity;
    int force;      // force run under various conditions. Currently only to skip out-of-phase transcripts

    struct {
        int unknown_chr,unknown_tscript_biotype,unknown_strand,unknown_phase,duplicate_id;
        int unknown_cds_phase,incomplete_cds,wrong_phase,overlapping_cds;
    } warned;
};

static const char *gf_strings_noncoding[] =
{
    "MT_rRNA", "MT_tRNA", "lincRNA", "miRNA", "misc_RNA", "rRNA", "snRNA", "snoRNA", "processed_transcript",
    "antisense", "macro_lncRNA", "ribozyme", "sRNA", "scRNA", "scaRNA", "sense_intronic", "sense_overlapping",
    "pseudogene", "processed_pseudogene", "artifact", "IG_pseudogene", "IG_C_pseudogene", "IG_J_pseudogene",
    "IG_V_pseudogene", "TR_V_pseudogene", "TR_J_pseudogene", "MT_tRNA_pseudogene", "misc_RNA_pseudogene",
    "miRNA_pseudogene", "ribozyme", "retained_intron", "retrotransposed", "Trna_pseudogene", "transcribed_processed_pseudogene",
    "transcribed_unprocessed_pseudogene", "transcribed_unitary_pseudogene",    "translated_unprocessed_pseudogene",
    "translated_processed_pseudogene", "known_ncRNA", "unitary_pseudogene", "unprocessed_pseudogene",
    "LRG_gene", "3_prime_overlapping_ncRNA", "disrupted_domain", "vaultRNA", "bidirectional_promoter_lncRNA", "ambiguous_orf",
    "lncRNA"
};
static const char *gf_strings_coding[] = { "protein_coding", "polymorphic_pseudogene", "IG_C", "IG_D", "IG_J", "IG_LV", "IG_V", "TR_C", "TR_D", "TR_J", "TR_V", "NMD", "non_stop_decay"};
static const char *gf_strings_special[] = { "CDS", "exon", "3_prime_UTR", "5_prime_UTR" };

int gff_set(gff_t *gff, gff_opt_t key, ...)
{
    va_list args;
    switch (key)
    {
        case dump_fname:
            va_start(args, key);
            gff->dump_fname = va_arg(args,char*);
            va_end(args);
            return 0;

        case force_out_of_phase:
            va_start(args, key);
            gff->force = va_arg(args,int);
            va_end(args);
            return 0;

        case strip_chr_names:
            va_start(args, key);
            gff->strip_chr_names = va_arg(args,int);
            va_end(args);
            return 0;

        case verbosity:
            va_start(args, key);
            gff->verbosity = va_arg(args,int);
            va_end(args);
            return 0;

        default:
            error("The key %d is not supported with gff_set\n",key);
    }
    return 0;
}

void *gff_get(gff_t *gff, gff_opt_t key)
{
    switch (key)
    {
        case idx_cds: return gff->idx_cds;
        case idx_utr: return gff->idx_utr;
        case idx_exon: return gff->idx_exon;
        case idx_tscript: return gff->idx_tscript;
        default:
            error("The key %d is not supported with gff_get\n",key);
    }
    return NULL;
}

const char *gff_id2string(gff_t *gff, id_type_t type, int id)    // currently only transcript ids
{
    return gff->tscript_ids.str[id];
}

const char *gf_type2gff_string(int type)
{
    if ( !GF_is_coding(type) )
    {
        if ( type < (1<<GF_coding_bit) ) return gf_strings_noncoding[type-1];
        type &= (1<<(GF_coding_bit+1)) - 1;
        return gf_strings_special[type - 1];
    }
    type &= (1<<GF_coding_bit) - 1;
    return gf_strings_coding[type - 1];
}

/*
    gff parsing functions
*/
static inline int feature_set_seq(gff_t *gff, char *chr_beg, char *chr_end)
{
    aux_t *aux = &gff->init;
    char tmp = chr_end[1];
    chr_end[1] = 0;
    int iseq;
    if ( khash_str2int_get(aux->seq2int, chr_beg, &iseq)!=0 )
    {
        char *new_chr = strdup(chr_beg);
        hts_expand(char*, aux->nseq+1, aux->mseq, aux->seq);
        aux->seq[aux->nseq] = new_chr;
        iseq = khash_str2int_inc(aux->seq2int, aux->seq[aux->nseq]);
        aux->nseq++;
        assert( aux->nseq < 1<<29 );  // see gf_gene_t.iseq and ftr_t.iseq
    }
    chr_end[1] = tmp;
    return iseq;
}
static inline char *gff_skip(const char *line, char *ss)
{
    while ( *ss && *ss!='\t' ) ss++;
    if ( !*ss ) error("[%s:%d %s] Could not parse the line: %s\n",__FILE__,__LINE__,__FUNCTION__,line);
    return ss+1;
}
static inline void gff_parse_chr(gff_t *gff, const char *line, char **chr_beg, char **chr_end)
{
    char *se = (char*) line;
    while ( *se && *se!='\t' ) se++;
    if ( !*se ) error("[%s:%d %s] Could not parse the line: %s\n",__FILE__,__LINE__,__FUNCTION__,line);
    if ( gff->strip_chr_names && !strncasecmp("chr",line,3) ) line += 3;
    *chr_beg = (char*) line;
    *chr_end = se-1;
}
static inline char *gff_parse_beg_end(const char *line, char *ss, uint32_t *beg, uint32_t *end)
{
    char *se = ss;
    *beg = strtol(ss, &se, 10) - 1;
    if ( ss==se ) error("[%s:%d %s] Could not parse the line:\n\t%s\n\t%s\n",__FILE__,__LINE__,__FUNCTION__,line,ss);
    ss = se+1;
    *end = strtol(ss, &se, 10) - 1;
    if ( ss==se ) error("[%s:%d %s] Could not parse the line: %s\n",__FILE__,__LINE__,__FUNCTION__,line);
    return se+1;
}
static void gff_id_init(id_tbl_t *tbl)
{
    memset(tbl, 0, sizeof(*tbl));
    tbl->str2id = khash_str2int_init();
}
static void gff_id_destroy(id_tbl_t *tbl)
{
    khash_str2int_destroy_free(tbl->str2id);
    free(tbl->str);
}
static inline int gff_id_register(id_tbl_t *tbl, char *beg, char *end, uint32_t *id_ptr)
{
    char tmp = end[1];
    end[1] = 0;
    int id;
    if ( khash_str2int_get(tbl->str2id, beg, &id) < 0 )
    {
        id = tbl->nstr++;
        hts_expand(char*, tbl->nstr, tbl->mstr, tbl->str);
        tbl->str[id] = strdup(beg);
        khash_str2int_set(tbl->str2id, tbl->str[id], id);
    }
    end[1] = tmp;
    *id_ptr = id;
    return 0;
}
static inline int gff_parse_biotype(char *line)
{
    if ( !line ) return -1;
    switch (*line)
    {
        case 'p':
            if ( !strncmp(line,"protein_coding",14) ) return GF_PROTEIN_CODING;
            else if ( !strncmp(line,"pseudogene",10) ) return GF_PSEUDOGENE;
            else if ( !strncmp(line,"processed_transcript",20) ) return GF_PROCESSED_TRANSCRIPT;
            else if ( !strncmp(line,"processed_pseudogene",20) ) return GF_PROCESSED_PSEUDOGENE;
            else if ( !strncmp(line,"polymorphic_pseudogene",22) ) return GF_POLYMORPHIC_PSEUDOGENE;
            break;
        case 'a':
            if ( !strncmp(line,"artifact",8) ) return GF_ARTIFACT;
            else if ( !strncmp(line,"antisense",9) ) return GF_ANTISENSE;
            else if ( !strncmp(line,"ambiguous_orf",13) ) return GF_AMBIGUOUS_ORF;
            break;
        case 'I':
            if ( !strncmp(line,"IG_pseudogene",13) ) return GF_IG_PSEUDOGENE;
            else if ( !strncmp(line,"IG_C_pseudogene",15) ) return GF_IG_C_PSEUDOGENE;
            else if ( !strncmp(line,"IG_J_pseudogene",15) ) return GF_IG_J_PSEUDOGENE;
            else if ( !strncmp(line,"IG_V_pseudogene",15) ) return GF_IG_V_PSEUDOGENE;
            else if ( !strncmp(line,"IG_C",4) ) return GF_IG_C;
            else if ( !strncmp(line,"IG_D",4) ) return GF_IG_D;
            else if ( !strncmp(line,"IG_J",4) ) return GF_IG_J;
            else if ( !strncmp(line,"IG_V",4) ) return GF_IG_V;
            else if ( !strncmp(line,"IG_LV",5) ) return GF_IG_LV;
            break;
        case 'T':
            if ( !strncmp(line,"TR_V_pseudogene",15) ) return GF_TR_V_PSEUDOGENE;
            else if ( !strncmp(line,"TR_J_pseudogene",15) ) return GF_TR_J_PSEUDOGENE;
            else if ( !strncmp(line,"TR_C",4) ) return GF_TR_C;
            else if ( !strncmp(line,"TR_D",4) ) return GF_TR_D;
            else if ( !strncmp(line,"TR_J",4) ) return GF_TR_J;
            else if ( !strncmp(line,"TR_V",4) ) return GF_TR_V;
            break;
        case 'M':
            if ( !strncmp(line,"Mt_tRNA_pseudogene",18) ) return GF_MT_tRNA_PSEUDOGENE;
            else if ( !strncasecmp(line,"Mt_tRNA",7) ) return GF_MT_tRNA;
            else if ( !strncasecmp(line,"Mt_rRNA",7) ) return GF_MT_tRNA;
            else if ( !strncasecmp(line,"MRNA",4) ) return GF_PROTEIN_CODING;
            break;
        case 'l':
            if ( !strncmp(line,"lincRNA",7) ) return GF_lincRNA;
            if ( !strncmp(line,"lncRNA",7) ) return GF_lncRNA;
            break;
        case 'm':
            if ( !strncmp(line,"macro_lncRNA",12) ) return GF_macro_lncRNA;
            else if ( !strncmp(line,"misc_RNA_pseudogene",19) ) return GF_misc_RNA_PSEUDOGENE;
            else if ( !strncmp(line,"miRNA_pseudogene",16) ) return GF_miRNA_PSEUDOGENE;
            else if ( !strncmp(line,"miRNA",5) ) return GF_miRNA;
            else if ( !strncmp(line,"misc_RNA",8) ) return GF_MISC_RNA;
            else if ( !strncasecmp(line,"mRNA",4) ) return GF_PROTEIN_CODING;
            break;
        case 'r':
            if ( !strncmp(line,"rRNA",4) ) return GF_rRNA;
            else if ( !strncmp(line,"ribozyme",8) ) return GF_RIBOZYME;
            else if ( !strncmp(line,"retained_intron",15) ) return GF_RETAINED_INTRON;
            else if ( !strncmp(line,"retrotransposed",15) ) return GF_RETROTRANSPOSED;
            break;
        case 's':
            if ( !strncmp(line,"snRNA",5) ) return GF_snRNA;
            else if ( !strncmp(line,"sRNA",4) ) return GF_sRNA;
            else if ( !strncmp(line,"scRNA",5) ) return GF_scRNA;
            else if ( !strncmp(line,"scaRNA",6) ) return GF_scaRNA;
            else if ( !strncmp(line,"snoRNA",6) ) return GF_snoRNA;
            else if ( !strncmp(line,"sense_intronic",14) ) return GF_SENSE_INTRONIC;
            else if ( !strncmp(line,"sense_overlapping",17) ) return GF_SENSE_OVERLAPPING;
            break;
        case 't':
            if ( !strncmp(line,"tRNA_pseudogene",15) ) return GF_tRNA_PSEUDOGENE;
            else if ( !strncmp(line,"transcribed_processed_pseudogene",32) ) return GF_TRANSCRIBED_PROCESSED_PSEUDOGENE;
            else if ( !strncmp(line,"transcribed_unprocessed_pseudogene",34) ) return GF_TRANSCRIBED_UNPROCESSED_PSEUDOGENE;
            else if ( !strncmp(line,"transcribed_unitary_pseudogene",30) ) return GF_TRANSCRIBED_UNITARY_PSEUDOGENE;
            else if ( !strncmp(line,"translated_unprocessed_pseudogene",33) ) return GF_TRANSLATED_UNPROCESSED_PSEUDOGENE;
            else if ( !strncmp(line,"translated_processed_pseudogene",31) ) return GF_TRANSLATED_PROCESSED_PSEUDOGENE;
            break;
        case 'n':
            if ( !strncmp(line,"nonsense_mediated_decay",23) ) return GF_NMD;
            else if ( !strncmp(line,"non_stop_decay",14) ) return GF_NON_STOP_DECAY;
            break;
        case 'N':
            if ( !strncmp(line,"NMD",3) ) return GF_NMD;
            break;
        case 'k':
            if ( !strncmp(line,"known_ncrna",11) ) return GF_KNOWN_NCRNA;
            break;
        case 'u':
            if ( !strncmp(line,"unitary_pseudogene",18) ) return GF_UNITARY_PSEUDOGENE;
            else if ( !strncmp(line,"unprocessed_pseudogene",22) ) return GF_UNPROCESSED_PSEUDOGENE;
            break;
        case 'L':
            if ( !strncmp(line,"LRG_gene",8) ) return GF_LRG_GENE;
            break;
        case '3':
            if ( !strncasecmp(line,"3prime_overlapping_ncRNA",24) ) return GF_3PRIME_OVERLAPPING_ncRNA;
            else if ( !strncasecmp(line,"3_prime_overlapping_ncRNA",25) ) return GF_3PRIME_OVERLAPPING_ncRNA;
            break;
        case 'd':
            if ( !strncmp(line,"disrupted_domain",16) ) return GF_DISRUPTED_DOMAIN;
            break;
        case 'v':
            if ( !strncmp(line,"vaultRNA",8) ) return GF_vaultRNA;
            break;
        case 'b':
            if ( !strncmp(line,"bidirectional_promoter_lncRNA",29) ) return GF_BIDIRECTIONAL_PROMOTER_lncRNA;
            break;
    }
    return 0;
}
static inline int gff_ignored_biotype(gff_t *gff, char *ss, char *se)
{
    if ( !ss ) return 0;

    char tmp = se[1];
    se[1] = 0;

    char *key = ss;
    int n = 0;
    if ( khash_str2int_get(gff->init.ignored_biotypes, ss, &n)!=0 ) key = strdup(ss);
    khash_str2int_set(gff->init.ignored_biotypes, key, n+1);

    se[1] = tmp;
    return 1;
}
static gf_gene_t *gene_init(aux_t *aux, uint32_t gene_id)
{
    khint_t k = kh_get(int2gene, aux->gid2gene, (int)gene_id);
    gf_gene_t *gene = (k == kh_end(aux->gid2gene)) ? NULL : kh_val(aux->gid2gene, k);
    if ( !gene )
    {
        gene = (gf_gene_t*) calloc(1,sizeof(gf_gene_t));
        int ret;
        k = kh_put(int2gene, aux->gid2gene, (int)gene_id, &ret);
        kh_val(aux->gid2gene,k) = gene;
    }
    return gene;
}
static void gff_parse_transcript(gff_t *gff, const char *line, ftr_t *ftr)
{
    aux_t *aux = &gff->init;

    ftr->type = gff_parse_biotype(aux->biotype);
    if ( ftr->type <= 0 )
    {
        char tmp = aux->type_end[1];
        aux->type_end[1] = 0;
        ftr->type = gff_parse_biotype(aux->type);
        aux->type_end[1] = tmp;
    }
    if ( ftr->type <= 0 )
    {
        if ( !gff_ignored_biotype(gff,aux->biotype,aux->biotype_end) )
        {
            if ( gff->verbosity > 0 )
            {
                if ( !gff->warned.unknown_tscript_biotype || gff->verbosity > 1 )
                    fprintf(stderr,"Warning: Ignoring transcript with unknown biotype .. %s\n", line);
                gff->warned.unknown_tscript_biotype++;
            }
        }
        return;
    }

    if ( !aux->id )
        error("[%s:%d %s] Could not parse the line, neither \"ID=transcript:\" nor \"ID=\" substring is present: %s\n",__FILE__,__LINE__,__FUNCTION__,line);
    if ( !aux->parent )
        error("[%s:%d %s] Could not parse the line, neither \"Parent=gene:\" nor \"Parent=\" substring is present: %s\n",__FILE__,__LINE__,__FUNCTION__,line);

    uint32_t trid,gene_id;
    gff_id_register(&gff->tscript_ids, aux->id, aux->id_end, &trid);
    gff_id_register(&aux->gene_ids, aux->parent, aux->parent_end, &gene_id);

    gf_tscript_t *tr = (gf_tscript_t*) calloc(1,sizeof(gf_tscript_t));
    tr->id     = trid;
    tr->strand = ftr->strand;
    tr->gene   = gene_init(aux, gene_id);
    tr->type   = ftr->type;
    tr->beg    = ftr->beg;
    tr->end    = ftr->end;

    khint_t k;
    int ret;
    k = kh_put(int2tscript, aux->id2tr, (int)trid, &ret);
    kh_val(aux->id2tr,k) = tr;
}
// register exon, CDS, UTR
static void gff_parse_exon(gff_t *gff, const char *line, ftr_t *ftr)
{
    aux_t *aux = &gff->init;
    if ( !aux->parent )
        error("[%s:%d %s] Could not parse the line, neither \"Parent=transcript:\" nor \"Parent=\" substring found: %s\n",__FILE__,__LINE__,__FUNCTION__,line);

    // associate with transcript id
    gff_id_register(&gff->tscript_ids, aux->parent, aux->parent_end, &ftr->trid);

    if ( ftr->strand==STRAND_UNK && gff->verbosity > 0 )
    {
        if ( !gff->warned.unknown_strand || gff->verbosity > 1 )
            fprintf(stderr,"Warning: Ignoring GFF feature with unknown strand .. %s\n",line);
        gff->warned.unknown_strand++;
    }
    if ( ftr->phase==CDS_PHASE_UNKN && gff->verbosity > 0 )
    {
        if ( !gff->warned.unknown_phase|| gff->verbosity > 1 )
            fprintf(stderr,"Warning: Ignoring GFF feature with unknown phase .. %s\n",line);
        gff->warned.unknown_phase++;
    }
    ftr->iseq = feature_set_seq(gff, aux->chr,aux->chr_end);
}
static void gff_parse_gene(gff_t *gff, const char *line, ftr_t *ftr)
{
    aux_t *aux = &gff->init;
    if ( !aux->id ) return;

    uint32_t gene_id;
    gff_id_register(&aux->gene_ids, aux->id, aux->id_end, &gene_id);

    gf_gene_t *gene = gene_init(aux, gene_id);
    if ( gene->name )
    {
        if ( !gff->warned.duplicate_id || gff->verbosity > 1 )
            fprintf(stderr,"Warning: The GFF contains features with duplicate id .. %s\n",line);
        gff->warned.duplicate_id++;
        return;
    }

    gene->iseq   = feature_set_seq(gff, aux->chr,aux->chr_end);
    gene->beg    = ftr->beg;
    gene->end    = ftr->end;
    gene->strand = ftr->strand;
    gene->id     = gene_id;

    if ( aux->name )
    {
        gene->name = (char*) malloc(aux->name_end - aux->name + 2);
        memcpy(gene->name,aux->name,aux->name_end - aux->name + 1);
        gene->name[aux->name_end - aux->name + 1] = 0;
    }
    else
        gene->name = strdup(aux->gene_ids.str[gene_id]); // Name=<GeneName> field is not present, use the gene ID instead
}

// Returns 0 for exons,CDS,UTRs to indicate these need to be pruned later and regidx built on them,
// or -1 to indicate the structure needs not be saved (either because of an error or because saved
// as transcript or gene.)
static int gff_parse_line(gff_t *gff, char *line, ftr_t *ftr)
{
    // - skip empty lines and commented lines
    // - columns
    //      1.      chr
    //      2.      <skip>
    //      3.      CDS, transcript, gene, ...
    //      4-5.    beg,end
    //      6.      <skip>
    //      7.      strand
    //      8.      phase
    //      9.      Parent=transcript:ENST(\d+);ID=...;biotype=... etc

    char *ss = line;
    if ( !*ss ) return -1;      // skip blank lines
    if ( *ss=='#' ) return -1;  // skip comments

    aux_t *aux = &gff->init;
    gff_parse_chr(gff, line, &aux->chr, &aux->chr_end);
    ss = gff_skip(line, aux->chr_end + 2);

    // 3rd column: is this a CDS, transcript, gene, etc.. The parsing order by frequency in Homo_sapiens.GRCh37.87.gff3
    int is_gene_line = 0;
    ftr->type = 0;
    aux->type = ss;
    if ( !strncmp("exon\t",ss,5) ) { ftr->type = GF_EXON; ss += 5; }
    else if ( !strncmp("CDS\t",ss,4) ) { ftr->type = GF_CDS; ss += 4; }
    else if ( !strncmp("three_prime_UTR\t",ss,16) ) { ftr->type = GF_UTR3; ss += 16; }
    else if ( !strncmp("five_prime_UTR\t",ss,15) ) { ftr->type = GF_UTR5; ss += 15; }
    else if ( !strncmp("biological_region\t",ss,18) ) { return -1; }    // skip
    else if ( !strncmp("gene\t",ss,5) ) { is_gene_line = 1; ss += 5; }
    else ss = gff_skip(line, ss);
    aux->type_end = ss - 1;

    // 4-5th columns: beg,end
    ss = gff_parse_beg_end(line, ss, &ftr->beg,&ftr->end);

    // 6th column: skip
    ss = gff_skip(line, ss);

    // 7th column: strand
    ftr->strand = -1;
    if ( *ss == '+' ) ftr->strand = STRAND_FWD;
    else if ( *ss == '-' ) ftr->strand = STRAND_REV;
    else ftr->strand = STRAND_UNK;
    ss += 2;

    // 8th column: phase (codon offset)
    ftr->phase = CDS_PHASE_UNKN;
    if ( *ss == '0' ) ftr->phase = 0;
    else if ( *ss == '1' ) ftr->phase = 1;
    else if ( *ss == '2' ) ftr->phase = 2;
    else if ( *ss == '.' ) ftr->phase = CDS_PHASE_UNKN;     // exons and even CDS in some GFFs do not have phase
    ss += 2;

    // 9th column: id, parent, name, biotype
    aux->name = NULL, aux->id = NULL, aux->parent = NULL, aux->biotype = NULL;
    while ( *ss )
    {
        char *es = ss;
        while ( *es && *es!=';' ) es++;
        if ( !strncmp(ss,"ID=",3) )
        {
            ss += 3;
            aux->id_end = es - 1;
            aux->id = ss;
            if ( !strncmp(ss,"gene:",5) ) { aux->id += 5; is_gene_line = 1; }
            else if ( !strncmp(ss,"transcript:",11) ) aux->id += 11;
        }
        else if ( !strncmp(ss,"Name=",5) ) { aux->name = ss + 5; aux->name_end = es - 1; }
        else if ( !strncmp(ss,"Parent=",7) )
        {
            ss += 7;
            aux->parent_end = es - 1;
            aux->parent = ss;
            if ( !strncmp(ss,"gene:",5) ) aux->parent += 5;
            else if ( !strncmp(ss,"transcript:",11) ) aux->parent += 11;
        }
        else if ( !strncmp(ss,"biotype=",8) ) { aux->biotype = ss + 8; aux->biotype_end = es - 1; }
        else if ( !strncmp(ss,"gene_biotype=",13) ) { aux->biotype = ss + 13; aux->biotype_end = es - 1; }
        if ( !*es ) break;
        ss = es + 1;
    }

    if ( is_gene_line || !aux->parent )
    {
        gff_parse_gene(gff, line, ftr);
        return -1;
    }

    if ( ftr->type )
    {
        gff_parse_exon(gff, line, ftr);
        return 0;
    }

    gff_parse_transcript(gff, line, ftr);
    return -1;
}

static int cmp_cds_ptr(const void *a, const void *b)
{
    // comparison function for qsort of transcripts's CDS
    if ( (*((gf_cds_t**)a))->beg < (*((gf_cds_t**)b))->beg ) return -1;
    if ( (*((gf_cds_t**)a))->beg > (*((gf_cds_t**)b))->beg ) return 1;
    return 0;
}

static inline void chr_beg_end(aux_t *aux, int iseq, char **chr_beg, char **chr_end)
{
    *chr_beg = *chr_end = aux->seq[iseq];
    while ( (*chr_end)[1] ) (*chr_end)++;
}
static gf_tscript_t *tscript_init(aux_t *aux, uint32_t trid)
{
    khint_t k = kh_get(int2tscript, aux->id2tr, (int)trid);
    gf_tscript_t *tr = (k == kh_end(aux->id2tr)) ? NULL : kh_val(aux->id2tr, k);
    assert( tr );
    return tr;
}
static void register_cds(gff_t *gff, ftr_t *ftr)
{
    // Make the CDS searchable via idx_cds. Note we do not malloc tr->cds just yet.
    //  ftr is the result of parsing a gff CDS line
    aux_t *aux = &gff->init;

    gf_tscript_t *tr = tscript_init(aux, ftr->trid);
    if ( tr->strand != ftr->strand ) error("Conflicting strand in transcript %"PRIu32" .. %d vs %d\n",ftr->trid,tr->strand,ftr->strand);

    gf_cds_t *cds = (gf_cds_t*) malloc(sizeof(gf_cds_t));
    cds->tr    = tr;
    cds->beg   = ftr->beg;
    cds->len   = ftr->end - ftr->beg + 1;
    cds->icds  = 0;     // to keep valgrind on mac happy
    cds->phase = ftr->phase;

    hts_expand(gf_cds_t*,tr->ncds+1,tr->mcds,tr->cds);
    tr->cds[tr->ncds++] = cds;
}
static void register_utr(gff_t *gff, ftr_t *ftr)
{
    aux_t *aux = &gff->init;
    gf_utr_t *utr = (gf_utr_t*) malloc(sizeof(gf_utr_t));
    utr->which = ftr->type==GF_UTR3 ? prime3 : prime5;
    utr->beg   = ftr->beg;
    utr->end   = ftr->end;
    utr->tr    = tscript_init(aux, ftr->trid);

    char *chr_beg, *chr_end;
    chr_beg_end(&gff->init, utr->tr->gene->iseq, &chr_beg, &chr_end);
    regidx_push(gff->idx_utr, chr_beg,chr_end, utr->beg,utr->end, &utr);
}
static void register_exon(gff_t *gff, ftr_t *ftr)
{
    aux_t *aux = &gff->init;
    gf_exon_t *exon = (gf_exon_t*) malloc(sizeof(gf_exon_t));
    exon->beg = ftr->beg;
    exon->end = ftr->end;
    exon->tr  = tscript_init(aux, ftr->trid);

    char *chr_beg, *chr_end;
    chr_beg_end(&gff->init, exon->tr->gene->iseq, &chr_beg, &chr_end);
    regidx_push(gff->idx_exon, chr_beg,chr_end, exon->beg - N_SPLICE_REGION_INTRON, exon->end + N_SPLICE_REGION_INTRON, &exon);
}

static void tscript_init_cds(gff_t *gff)
{
    aux_t *aux = &gff->init;

    // Sort CDS in all transcripts, set offsets, check their phase, length, create index (idx_cds)
    khint_t k;
    for (k=0; k<kh_end(aux->id2tr); k++)
    {
        if ( !kh_exist(aux->id2tr, k) ) continue;
        gf_tscript_t *tr = (gf_tscript_t*) kh_val(aux->id2tr, k);

        // position-to-tscript lookup
        char *chr_beg, *chr_end;
        chr_beg_end(aux, tr->gene->iseq, &chr_beg, &chr_end);
        regidx_push(gff->idx_tscript, chr_beg, chr_end, tr->beg, tr->end, &tr);

        if ( !tr->ncds ) continue;      // transcript with no CDS

        // sort CDs
        qsort(tr->cds, tr->ncds, sizeof(gf_cds_t*), cmp_cds_ptr);

        // trim non-coding start
        int i, len = 0;
        if ( tr->strand==STRAND_FWD )
        {
            if ( tr->cds[0]->phase != CDS_PHASE_UNKN )
            {
                if ( tr->cds[0]->phase ) tr->trim |= TRIM_5PRIME;
                tr->cds[0]->beg += tr->cds[0]->phase;
                tr->cds[0]->len -= tr->cds[0]->phase;
                tr->cds[0]->phase = 0;
            }

            // sanity check phase; the phase number in gff tells us how many bases to skip in this
            // feature to reach the first base of the next codon
            int tscript_ok = 1;
            for (i=0; i<tr->ncds; i++)
            {
                if ( tr->cds[i]->phase == CDS_PHASE_UNKN )
                {
                    if ( gff->verbosity > 0 )
                    {
                        if ( !gff->warned.unknown_cds_phase || gff->verbosity > 1 )
                            fprintf(stderr,"Warning: CDS with unknown phase, could not verify reading frame in transcript %s\n",gff->tscript_ids.str[tr->id]);
                        gff->warned.unknown_cds_phase++;
                    }
                    len += tr->cds[i]->len;
                    continue;
                }
                int phase = tr->cds[i]->phase ? 3 - tr->cds[i]->phase : 0;
                if ( phase!=len%3 )
                {
                    if ( !gff->force )
                        error("Error: GFF3 assumption failed for transcript %s, CDS=%"PRIu32": phase!=len%%3 (phase=%d, len=%d). Use the --force option to proceed anyway (at your own risk).\n",
                                gff->tscript_ids.str[tr->id],tr->cds[i]->beg+1,phase,len);
                    if ( gff->verbosity > 0 )
                    {
                        if ( !gff->warned.wrong_phase || gff->verbosity > 1 )
                            fprintf(stderr,"Warning: The GFF has inconsistent phase column in transcript %s, skipping. CDS pos=%"PRIu32": phase!=len%%3 (phase=%d, len=%d)\n",
                                    gff->tscript_ids.str[tr->id],tr->cds[i]->beg+1,phase,len);
                        gff->warned.wrong_phase++;
                    }
                    tscript_ok = 0;
                    break;
                }
                len += tr->cds[i]->len;
            }
            if ( !tscript_ok ) continue;    // skip this transcript
        }
        else if ( tr->strand==STRAND_REV )
        {
            if ( tr->cds[tr->ncds-1]->phase != CDS_PHASE_UNKN )
            {
                // Check that the phase is not bigger than CDS length. Curiously, this can really happen,
                // see Mus_musculus.GRCm38.85.gff3.gz, transcript:ENSMUST00000163141.
                // This also fixes phase of 5' incomplete CDS, see test/csq/ENST00000520868/ENST00000520868.gff
                // todo: the same for the fwd strand
                i = tr->ncds - 1;
                int phase = tr->cds[i]->phase;
                if ( phase ) tr->trim |= TRIM_5PRIME;
                while ( i>=0 && phase > tr->cds[i]->len )
                {
                    phase -= tr->cds[i]->len;
                    tr->cds[i]->phase = 0;
                    tr->cds[i]->len   = 0;
                    i--;
                }
                if ( gff->verbosity > 0 && tr->cds[i]->phase )
                {
                    if ( !gff->warned.incomplete_cds || gff->verbosity > 1 )
                        fprintf(stderr,"Note: truncated transcript %s with incomplete CDS (this is very common)\n",gff->tscript_ids.str[tr->id]);
                    gff->warned.incomplete_cds++;
                }
                tr->cds[i]->len  -= tr->cds[i]->phase;
                tr->cds[i]->phase = 0;
            }

            // sanity check phase
            int tscript_ok = 1;
            for (i=tr->ncds-1; i>=0; i--)
            {
                if ( tr->cds[i]->phase == CDS_PHASE_UNKN )
                {
                    if ( gff->verbosity > 0 )
                    {
                        if ( !gff->warned.unknown_cds_phase || gff->verbosity > 1 )
                            fprintf(stderr,"Warning: CDS with unknown phase, could not verify reading frame in transcript %s\n",gff->tscript_ids.str[tr->id]);
                        gff->warned.unknown_cds_phase++;
                    }
                    len += tr->cds[i]->len;
                    continue;
                }
                int phase = tr->cds[i]->phase ? 3 - tr->cds[i]->phase : 0;
                if ( phase!=len%3 )
                {
                    if ( !gff->force )
                        error("Error: GFF3 assumption failed for transcript %s, CDS=%"PRIu32": phase!=len%%3 (phase=%d, len=%d). Use the --force option to proceed anyway (at your own risk).\n",
                                gff->tscript_ids.str[tr->id],tr->cds[i]->beg+1,phase,len);
                    if ( gff->verbosity > 0 )
                    {
                        if ( !gff->warned.wrong_phase || gff->verbosity > 1 )
                            fprintf(stderr,"Warning: The GFF has inconsistent phase column in transcript %s, skipping. CDS pos=%"PRIu32": phase!=len%%3 (phase=%d, len=%d)\n",
                                    gff->tscript_ids.str[tr->id],tr->cds[i]->beg+1,phase,len);
                        gff->warned.wrong_phase++;
                    }
                    tscript_ok = 0;
                    break;
                }
                len += tr->cds[i]->len;
            }
            if ( !tscript_ok ) continue;    // skip this transcript
        }
        else
            continue;   // unknown strand

        // set len. At the same check that CDS within a transcript do not overlap
        len = 0;
        for (i=0; i<tr->ncds; i++)
        {
            tr->cds[i]->icds = i;
            len += tr->cds[i]->len;
            if ( !i ) continue;

            gf_cds_t *a = tr->cds[i-1];
            gf_cds_t *b = tr->cds[i];
            if ( a->beg + a->len - 1 >= b->beg )
            {
                if ( gff->verbosity > 0 )
                {
                    if ( !gff->warned.overlapping_cds || gff->verbosity > 1 )
                        fprintf(stderr,"Warning: GFF contains overlapping CDS %s, %"PRIu32"-%"PRIu32" and %"PRIu32"-%"PRIu32" (ribosomal slippage?)\n",
                                gff->tscript_ids.str[tr->id], a->beg+1,a->beg+a->len, b->beg+1,b->beg+b->len);
                    gff->warned.overlapping_cds++;
                }
            }
        }

        if ( len%3 != 0 )
        {
            // There are 13k transcripts with incomplete 3' CDS. See for example ENST00000524289
            //  http://sep2015.archive.ensembl.org/Homo_sapiens/Transcript/Sequence_cDNA?db=core;g=ENSG00000155868;r=5:157138846-157159019;t=ENST00000524289
            // Also, the incomplete CDS can be too short (1 or 2bp), so it is not enough to trim the last one.

            if ( gff->verbosity > 0 )
            {
                if ( !gff->warned.incomplete_cds || gff->verbosity > 1 )
                    fprintf(stderr,"Note: truncated transcript %s with incomplete CDS (this is very common)\n",gff->tscript_ids.str[tr->id]);
                gff->warned.incomplete_cds++;
            }

            tr->trim |= TRIM_3PRIME;
            if ( tr->strand==STRAND_FWD )
            {
                i = tr->ncds - 1;
                while ( i>=0 && len%3 )
                {
                    int dlen = tr->cds[i]->len >= len%3 ? len%3 : tr->cds[i]->len;
                    tr->cds[i]->len -= dlen;
                    len -= dlen;
                    i--;
                }
            }
            else if ( tr->strand==STRAND_REV )
            {
                i = 0;
                while ( i<tr->ncds && len%3 )
                {
                    int dlen = tr->cds[i]->len >= len%3 ? len%3 : tr->cds[i]->len;
                    tr->cds[i]->len -= dlen;
                    tr->cds[i]->beg += dlen;
                    len -= dlen;
                    i++;
                }
            }
        }

        // set CDS offsets and insert into regidx
        len=0;
        for (i=0; i<tr->ncds; i++)
        {
            tr->cds[i]->pos = len;
            len += tr->cds[i]->len;
            regidx_push(gff->idx_cds, chr_beg,chr_end, tr->cds[i]->beg,tr->cds[i]->beg+tr->cds[i]->len-1, &tr->cds[i]);
        }
    }
}

static void regidx_free_gf(void *payload) { free(*((gf_cds_t**)payload)); }
static void regidx_free_tscript(void *payload) { gf_tscript_t *tr = *((gf_tscript_t**)payload); free(tr->cds); free(tr); }

static int gff_dump(gff_t *gff, const char *fname)
{
    BGZF *out = bgzf_open(fname,"wg");
    if ( !out ) error("Failed to open %s: %s\n", fname, strerror(errno));

    kstring_t str = {0,0,0};

    khint_t k;
    for (k=0; k<kh_end(gff->init.gid2gene); k++)
    {
        if ( !kh_exist(gff->init.gid2gene, k) ) continue;
        gf_gene_t *gene = (gf_gene_t*) kh_val(gff->init.gid2gene, k);
        char *gene_id = gff->init.gene_ids.str[gene->id];
        str.l = 0;
        ksprintf(&str,"%s\t.\tgene\t%"PRIu32"\t%"PRIu32"\t.\t%c\t.\tID=%s;Name=%s;used=%d\n",gff->init.seq[gene->iseq],gene->beg+1,gene->end+1,gene->strand==STRAND_FWD?'+':(gene->strand==STRAND_REV?'-':'.'),gene_id,gene->name,gene->used);
        if ( bgzf_write(out, str.s, str.l) != str.l ) error("Error writing %s: %s\n", fname, strerror(errno));
    }

    regitr_t *itr = regitr_init(gff->idx_tscript);
    while ( regitr_loop(itr) )
    {
        gf_tscript_t *tr = regitr_payload(itr, gf_tscript_t*);
        char *gene_id =  gff->init.gene_ids.str[tr->gene->id];
        const char *type = tr->type==GF_PROTEIN_CODING ? "mRNA" : gf_type2gff_string(tr->type);
        str.l = 0;
        ksprintf(&str,"%s\t.\t%s\t%"PRIu32"\t%"PRIu32"\t.\t%c\t.\tID=%s;Parent=%s;biotype=%s;used=%d\n",itr->seq,type,itr->beg+1,itr->end+1,tr->strand==STRAND_FWD?'+':(tr->strand==STRAND_REV?'-':'.'),gff->tscript_ids.str[tr->id],gene_id,gf_type2gff_string(tr->type),tr->used);
        if ( bgzf_write(out, str.s, str.l) != str.l ) error("Error writing %s: %s\n", fname, strerror(errno));
    }
    regitr_destroy(itr);

    itr = regitr_init(gff->idx_cds);
    while ( regitr_loop(itr) )
    {
        gf_cds_t *cds = regitr_payload(itr,gf_cds_t*);
        gf_tscript_t *tr = cds->tr;
        str.l = 0;
        ksprintf(&str,"%s\t.\tCDS\t%"PRIu32"\t%"PRIu32"\t.\t%c\t%c\tParent=%s\n",itr->seq,cds->beg+1,cds->beg+cds->len,tr->strand==STRAND_FWD?'+':(tr->strand==STRAND_REV?'-':'.'),cds->phase==3?'.':cds->phase+(int)'0',gff->tscript_ids.str[tr->id]);
        if ( bgzf_write(out, str.s, str.l) != str.l ) error("Error writing %s: %s\n", fname, strerror(errno));
    }
    regitr_destroy(itr);

    itr = regitr_init(gff->idx_utr);
    while ( regitr_loop(itr) )
    {
        gf_utr_t *utr = regitr_payload(itr,gf_utr_t*);
        gf_tscript_t *tr = utr->tr;
        str.l = 0;
        ksprintf(&str,"%s\t.\t%s_prime_UTR\t%"PRIu32"\t%"PRIu32"\t.\t%c\t.\tParent=%s\n",itr->seq,utr->which==prime3?"three":"five",utr->beg+1,utr->end+1,tr->strand==STRAND_FWD?'+':(tr->strand==STRAND_REV?'-':'.'),gff->tscript_ids.str[tr->id]);
        if ( bgzf_write(out, str.s, str.l) != str.l ) error("Error writing %s: %s\n", fname, strerror(errno));
    }
    regitr_destroy(itr);

    itr = regitr_init(gff->idx_exon);
    while ( regitr_loop(itr) )
    {
        gf_exon_t *exon = regitr_payload(itr,gf_exon_t*);
        gf_tscript_t *tr = exon->tr;
        str.l = 0;
        ksprintf(&str,"%s\t.\texon\t%"PRIu32"\t%"PRIu32"\t.\t%c\t.\tParent=%s\n",itr->seq,exon->beg+1,exon->end+1,tr->strand==STRAND_FWD?'+':(tr->strand==STRAND_REV?'-':'.'),gff->tscript_ids.str[tr->id]);
        if ( bgzf_write(out, str.s, str.l) != str.l ) error("Error writing %s: %s\n", fname, strerror(errno));
    }
    regitr_destroy(itr);

    if ( bgzf_close(out)!=0 ) error("Error: close failed .. %s\n", fname);
    free(str.s);

    return 0;
}

int gff_parse(gff_t *gff)
{
    if ( gff->verbosity > 0 ) fprintf(stderr,"Parsing %s ...\n", gff->fname);

    aux_t *aux = &gff->init;
    aux->seq2int   = khash_str2int_init();   // chrom's numeric id
    aux->gid2gene  = kh_init(int2gene);      // gene id to gf_gene_t, for idx_gene
    aux->id2tr     = kh_init(int2tscript);   // transcript id to tscript_t
    gff->idx_tscript = regidx_init(NULL, NULL, regidx_free_tscript, sizeof(gf_tscript_t*), NULL);
    aux->ignored_biotypes = khash_str2int_init();
    gff_id_init(&aux->gene_ids);
    gff_id_init(&gff->tscript_ids);

    // parse gff
    kstring_t str = {0,0,0};
    htsFile *fp = hts_open(gff->fname,"r");
    if ( !fp ) error("Failed to read %s\n", gff->fname);
    while ( hts_getline(fp, KS_SEP_LINE, &str) > 0 )
    {
        hts_expand(ftr_t, aux->nftr+1, aux->mftr, aux->ftr);
        int ret = gff_parse_line(gff, str.s, aux->ftr + aux->nftr);
        if ( !ret ) aux->nftr++;
    }
    free(str.s);
    if ( hts_close(fp)!=0 ) error("Close failed: %s\n", gff->fname);


    // process gff information: connect CDS and exons to transcripts
    gff->idx_cds  = regidx_init(NULL, NULL, regidx_free_gf, sizeof(gf_cds_t*), NULL);
    gff->idx_utr  = regidx_init(NULL, NULL, regidx_free_gf, sizeof(gf_utr_t*), NULL);
    gff->idx_exon = regidx_init(NULL, NULL, regidx_free_gf, sizeof(gf_exon_t*), NULL);

    int i;
    for (i=0; i<aux->nftr; i++)
    {
        ftr_t *ftr = &aux->ftr[i];

        // check whether to keep this feature: is there a mapping trid -> gene_id -> gene?
        khint_t k = kh_get(int2tscript, aux->id2tr, (int)ftr->trid);
        if ( k==kh_end(aux->id2tr) ) continue;       // no corresponding transcript registered, must be an unsupported biotype

        gf_tscript_t *tr = kh_val(aux->id2tr,k);
        tr->used = 1;
        tr->gene->used = 1;

        // populate regidx by category:
        //      ftr->type   .. GF_CDS, GF_EXON, GF_UTR3, GF_UTR5
        //      gene->type  .. GF_PROTEIN_CODING, GF_MT_rRNA, GF_IG_C, ...
        if ( ftr->type==GF_CDS ) register_cds(gff, ftr);
        else if ( ftr->type==GF_EXON ) register_exon(gff, ftr);
        else if ( ftr->type==GF_UTR5 ) register_utr(gff, ftr);
        else if ( ftr->type==GF_UTR3 ) register_utr(gff, ftr);
        else
            error("something: %s\t%"PRIu32"\t%"PRIu32"\t%s\t%s\n", aux->seq[ftr->iseq],ftr->beg+1,ftr->end+1,gff->tscript_ids.str[ftr->trid],gf_type2gff_string(ftr->type));
    }
    tscript_init_cds(gff);

    if ( gff->verbosity > 0 )
    {
        fprintf(stderr,"Indexed %d transcripts, %d exons, %d CDSs, %d UTRs\n",
                regidx_nregs(gff->idx_tscript),
                regidx_nregs(gff->idx_exon),
                regidx_nregs(gff->idx_cds),
                regidx_nregs(gff->idx_utr));
    }

    if ( gff->verbosity > 0 && khash_str2int_size(aux->ignored_biotypes) )
    {
        khash_t(str2int) *ign = (khash_t(str2int)*)aux->ignored_biotypes;
        fprintf(stderr,"Ignored the following biotypes:\n");
        for (i = kh_begin(ign); i < kh_end(ign); i++)
        {
            if ( !kh_exist(ign,i)) continue;
            const char *biotype = kh_key(ign,i);
            if ( !strcmp(biotype,"TCE") ) biotype = "TCE (\"To be Experimentally Confirmed\")";
            fprintf(stderr,"\t%dx\t.. %s\n", kh_value(ign,i), biotype);
        }
    }
    khash_str2int_destroy_free(aux->ignored_biotypes);

    // warned about unprinted warnings
    if ( gff->verbosity > 0 )
    {
        int nwarn = 0;
        #define INC_NWARN(X) if (gff->warned.X) nwarn += gff->verbosity > 1 ? 0 : gff->warned.X - 1;
        INC_NWARN(unknown_chr);
        INC_NWARN(unknown_tscript_biotype);
        INC_NWARN(unknown_strand);
        INC_NWARN(unknown_phase);
        INC_NWARN(duplicate_id);
        INC_NWARN(unknown_cds_phase);
        INC_NWARN(incomplete_cds);
        INC_NWARN(wrong_phase);
        INC_NWARN(overlapping_cds);
        if ( nwarn > 0 )
            fprintf(stderr,"Warning: %d warnings were suppressed, increase verbosity to see them all\n",nwarn);
    }

    if ( gff->dump_fname ) gff_dump(gff, gff->dump_fname);

    if (  !regidx_nregs(gff->idx_tscript) )
        error("Error: No usable transcripts found, likely a failure to parse a non-standard GFF file. Please check if the misc/gff2gff\n"
              "       or misc/gff2gff.py script can fix the problem (both do different things). See also the man page for the description\n"
              "       of the expected format http://samtools.github.io/bcftools/bcftools-man.html#csq\n");

    free(aux->seq);
    free(aux->ftr);
    khash_str2int_destroy_free(aux->seq2int);
    // keeping only to destroy the genes at the end: kh_destroy(int2gene,aux->gid2gene);
    kh_destroy(int2tscript,aux->id2tr);
    gff_id_destroy(&aux->gene_ids);

    return 0;
}

gff_t *gff_init(const char *fname)
{
    gff_t *gff = calloc(sizeof(gff_t),1);
    gff->fname = fname;
    return gff;
}
void gff_destroy(gff_t *gff)
{
    khint_t k;
    if ( gff->init.gid2gene )
    {
        for (k=0; k<kh_end(gff->init.gid2gene); k++)
        {
            if ( !kh_exist(gff->init.gid2gene, k) ) continue;
            gf_gene_t *gene = (gf_gene_t*) kh_val(gff->init.gid2gene, k);
            free(gene->name);
            free(gene);
        }
        kh_destroy(int2gene,gff->init.gid2gene);
    }

    regidx_destroy(gff->idx_cds);
    regidx_destroy(gff->idx_utr);
    regidx_destroy(gff->idx_exon);
    regidx_destroy(gff->idx_tscript);

    gff_id_destroy(&gff->tscript_ids);
    free(gff);
}

