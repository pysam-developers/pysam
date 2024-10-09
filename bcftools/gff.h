/* The MIT License

   Copyright (c) 2023-2024 Genome Research Ltd.

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
/*
    GFF parsing code refactored from csq.c

    Things that would be nice to have
        - dynamic N_REF_PAD
        - for stop-lost events (also in frameshifts) report the number of truncated aa's
        - memory could be greatly reduced by indexing gff (but it is quite compact already)
        - deletions that go beyond transcript boundaries are not checked at sequence level
            - alloc tscript->ref in hap_finalize, introduce fa_off_beg:16,fa_off_end:16
            - see test/csq/ENST00000573314/insertion-overlap.vcf #1476288882

    Read about transcript types here
        http://vega.sanger.ac.uk/info/about/gene_and_transcript_types.html
        https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
        https://www.gencodegenes.org/pages/biotypes.html

    List of supported biotypes
        antisense
        IG_C_gene
        IG_D_gene
        IG_J_gene
        IG_LV_gene
        IG_V_gene
        lincRNA
        lncRNA      .. generic term for 3prime_overlapping_ncRNA, antisense, bidirectional_promoter_lncRNA, lincRNA, macro_lncRNA, non_coding, processed_transcript, sense_intronic, sense_overlapping
        macro_lncRNA
        miRNA
        misc_RNA
        Mt_rRNA
        Mt_tRNA
        polymorphic_pseudogene
        processed_transcript
        protein_coding, mRNA
        ribozyme
        rRNA
        sRNA
        scRNA
        scaRNA
        sense_intronic
        sense_overlapping
        snRNA
        snoRNA
        TR_C_gene
        TR_D_gene
        TR_J_gene
        TR_V_gene

    The gff parsing logic
        We collect features such by combining gff lines A,B,C as follows:
            A .. gene line with a supported biotype
                    A.ID=~/^gene:/

            B .. transcript line referencing A with supported biotype
                    B.ID=~/^transcript:/ && B.Parent=~/^gene:A.ID/

            C .. corresponding CDS, exon, and UTR lines:
                    C[3] in {"CDS","exon","three_prime_UTR","five_prime_UTR"} && C.Parent=~/^transcript:B.ID/

        For coding biotypes ("protein_coding" or "polymorphic_pseudogene") the
        complete chain link C -> B -> A is required. For the rest, link B -> A suffices.


    The supported consequence types, sorted by impact:
        splice_acceptor_variant .. end region of an intron changed (2bp at the 3' end of an intron)
        splice_donor_variant    .. start region of an intron changed (2bp at the 5' end of an intron)
        stop_gained             .. DNA sequence variant resulting in a stop codon
        frameshift_variant      .. number of inserted/deleted bases not a multiple of three, disrupted translational frame
        stop_lost               .. elongated transcript, stop codon changed
        start_lost              .. the first codon changed
        inframe_altering        .. combination of indels leading to unchanged reading frame and length
        inframe_insertion       .. inserted coding sequence, unchanged reading frame
        inframe_deletion        .. deleted coding sequence, unchanged reading frame
        missense_variant        .. amino acid (aa) change, unchanged length
        splice_region_variant   .. change within 1-3 bases of the exon or 3-8 bases of the intron
        synonymous_variant      .. DNA sequence variant resulting in no amino acid change
        stop_retained_variant   .. different stop codon
        start_retained_variant  .. start codon retained by indel realignment
        non_coding_variant      .. variant in non-coding sequence, such as RNA gene
        5_prime_UTR_variant
        3_prime_UTR_variant
        intron_variant          .. reported only if none of the above
        intergenic_variant      .. reported only if none of the above


    The annotation algorithm.
        The algorithm checks if the variant falls in a region of a supported type. The
        search is performed in the following order, until a match is found:
            1. idx_cds(gf_cds_t) - lookup CDS by position, create haplotypes, call consequences
            2. idx_utr(gf_utr_t) - check UTR hits
            3. idx_exon(gf_exon_t) - check for splice variants
            4. idx_tscript(tscript_t) - check for intronic variants, RNAs, etc.

        These regidx indexes are created by parsing a gff3 file as follows:
            1.  create the array "ftr" of all UTR, CDS, exons. This will be
            processed later and pruned based on transcript types we want to keep.
            In the same go, create the hash "id2tr" of transcripts to keep
            (based on biotype) which maps from transcript_id to a transcript. At
            the same time also build the hash "gid2gene" which maps from gene_id to
            gf_gene_t pointer.

            2.  build "idx_cds", "idx_tscript", "idx_utr" and "idx_exon" indexes.
            Use only features from "ftr" which are present in "id2tr".

            3.  clean data that won't be needed anymore: ftr, id2tr, gid2gene.

    Data structures.
        idx_cds, idx_utr, idx_exon, idx_tscript:
            as described above, regidx structures for fast lookup of exons/transcripts
            overlapping a region, the payload is a pointer to tscript.cds
*/

#ifndef GFF_H__
#define GFF_H__

#include <stdint.h>

#ifndef __FUNCTION__
#  define __FUNCTION__ __func__
#endif

// Definition of splice_region, splice_acceptor and splice_donor
#define N_SPLICE_DONOR         2
#define N_SPLICE_REGION_EXON   3
#define N_SPLICE_REGION_INTRON 8

#define STRAND_REV 0
#define STRAND_FWD 1
#define STRAND_UNK 2

#define TRIM_NONE   0
#define TRIM_5PRIME 1
#define TRIM_3PRIME 2


// GFF line types
#define GFF_UNKN_LINE    0
#define GFF_TSCRIPT_LINE 1
#define GFF_GENE_LINE    2


/*
    Genomic features, for fast lookup by position to overlapping features
*/
#define GF_coding_bit 6
#define GF_is_coding(x) ((x) & (1<<GF_coding_bit))
#define GF_MT_rRNA                       1                      // non-coding: 1, 2, ...
#define GF_MT_tRNA                       2
#define GF_lincRNA                       3
#define GF_miRNA                         4
#define GF_MISC_RNA                      5
#define GF_rRNA                          6
#define GF_snRNA                         7
#define GF_snoRNA                        8
#define GF_PROCESSED_TRANSCRIPT          9
#define GF_ANTISENSE                    10
#define GF_macro_lncRNA                 11
#define GF_ribozyme                     12
#define GF_sRNA                         13
#define GF_scRNA                        14
#define GF_scaRNA                       15
#define GF_SENSE_INTRONIC               16
#define GF_SENSE_OVERLAPPING            17
#define GF_PSEUDOGENE                   18
#define GF_PROCESSED_PSEUDOGENE         19
#define GF_ARTIFACT                     20
#define GF_IG_PSEUDOGENE                21
#define GF_IG_C_PSEUDOGENE              22
#define GF_IG_J_PSEUDOGENE              23
#define GF_IG_V_PSEUDOGENE              24
#define GF_TR_V_PSEUDOGENE              25
#define GF_TR_J_PSEUDOGENE              26
#define GF_MT_tRNA_PSEUDOGENE           27
#define GF_misc_RNA_PSEUDOGENE          28
#define GF_miRNA_PSEUDOGENE             29
#define GF_RIBOZYME                     30
#define GF_RETAINED_INTRON              31
#define GF_RETROTRANSPOSED              32
#define GF_tRNA_PSEUDOGENE              33
#define GF_TRANSCRIBED_PROCESSED_PSEUDOGENE     34
#define GF_TRANSCRIBED_UNPROCESSED_PSEUDOGENE   35
#define GF_TRANSCRIBED_UNITARY_PSEUDOGENE       36
#define GF_TRANSLATED_UNPROCESSED_PSEUDOGENE    37
#define GF_TRANSLATED_PROCESSED_PSEUDOGENE      38
#define GF_KNOWN_NCRNA                          39
#define GF_UNITARY_PSEUDOGENE                   40
#define GF_UNPROCESSED_PSEUDOGENE               41
#define GF_LRG_GENE                             42
#define GF_3PRIME_OVERLAPPING_ncRNA             43
#define GF_DISRUPTED_DOMAIN                     44
#define GF_vaultRNA                             45
#define GF_BIDIRECTIONAL_PROMOTER_lncRNA        46
#define GF_AMBIGUOUS_ORF                        47
#define GF_lncRNA                               48
#define GF_PROTEIN_CODING               (1|(1<<GF_coding_bit))  // coding: 65, 66, ...
#define GF_POLYMORPHIC_PSEUDOGENE       (2|(1<<GF_coding_bit))
#define GF_IG_C                         (3|(1<<GF_coding_bit))
#define GF_IG_D                         (4|(1<<GF_coding_bit))
#define GF_IG_J                         (5|(1<<GF_coding_bit))
#define GF_IG_LV                        (6|(1<<GF_coding_bit))
#define GF_IG_V                         (7|(1<<GF_coding_bit))
#define GF_TR_C                         (8|(1<<GF_coding_bit))
#define GF_TR_D                         (9|(1<<GF_coding_bit))
#define GF_TR_J                        (10|(1<<GF_coding_bit))
#define GF_TR_V                        (11|(1<<GF_coding_bit))
#define GF_NMD                         (12|(1<<GF_coding_bit))
#define GF_NON_STOP_DECAY              (13|(1<<GF_coding_bit))
#define GF_CDS      ((1<<(GF_coding_bit+1))+1)                  // special types: 129, 130, ...
#define GF_EXON     ((1<<(GF_coding_bit+1))+2)
#define GF_UTR3     ((1<<(GF_coding_bit+1))+3)
#define GF_UTR5     ((1<<(GF_coding_bit+1))+4)
// GF_MAX = (1<<30)-1, see hap_node_t

#define CDS_PHASE_UNKN 3
typedef struct gf_tscript_t_ gf_tscript_t;
typedef struct
{
    gf_tscript_t *tr;   // transcript
    uint32_t beg;       // the start coordinate of the CDS (on the reference strand, 0-based)
    uint32_t pos;       // 0-based index of the first exon base within the transcript (only to
                        //  update hap_node_t.sbeg in hap_init, could be calculated on the fly)
    uint32_t len;       // exon length
    uint32_t icds:30,   // exon index within the transcript
             phase:2;   // offset of the CDS: 0,1,2 or 3 for unknown
}
gf_cds_t;
typedef struct
{
    char *name;                     // human readable name, e.g. ORF45
    uint32_t iseq;
    uint32_t id,beg,end,strand:31,  // used only by --dump-gff
             used:1;                // does it have any exons, CDS, UTR?
}
gf_gene_t;
typedef struct
{
    uint32_t beg,end;
    gf_tscript_t *tr;
}
gf_exon_t;
typedef enum { prime3, prime5 } utr_t;
typedef struct
{
    utr_t which;
    uint32_t beg,end;
    gf_tscript_t *tr;
}
gf_utr_t;
struct gf_tscript_t_
{
    uint32_t id;        // transcript id
    uint32_t beg,end;   // transcript's beg and end coordinate (ref strand, 0-based, inclusive)
    uint32_t strand:2,  // STRAND_REV,FWD,UNK
             used:1,    // does it have any exons, UTRs, CDS?
             ncds:29,   // number of exons
             mcds;
    gf_cds_t **cds;     // ordered list of exons
    uint32_t trim:2,    // complete, 5' or 3' trimmed, see TRIM_* types
             type:30;   // one of GF_* types
    gf_gene_t *gene;
    void *aux;          // auxiliary user data
};

typedef enum
{
    // write options
    verbosity,          // int, 0-2
    strip_chr_names,    // int, 0 to leave as is, 1 to strip 'chr' prefix
    force_out_of_phase, // int, 1 to proceed even CDS exon out of expected phase
    dump_fname,         // const char*, dump the parsed GFF into this file, for debugging purposes

    // read options
    idx_cds,
    idx_utr,
    idx_exon,
    idx_tscript,
}
gff_opt_t;

typedef enum { transcript } id_type_t;  // for gff_id2str

typedef struct gff_t_ gff_t;

gff_t *gff_init(const char *fname);
int gff_parse(gff_t *gff);
void gff_destroy(gff_t *gff);

int gff_set(gff_t *gff, gff_opt_t key, ...);   // returns 0 on success
void *gff_get(gff_t *gff, gff_opt_t key);
const char *gff_id2string(gff_t *gff, id_type_t type, int id);
const char *gf_type2gff_string(int type);

#endif
