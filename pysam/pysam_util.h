#ifndef PYSAM_UTIL_H
#define PYSAM_UTIL_H

// This file contains some of the definitions from bam_index.c

#define BAM_MIN_CHUNK_GAP 32768
// 1<<14 is the size of minimum bin.
#define BAM_LIDX_SHIFT    14
// =(8^6-1)/7+1
#define MAX_BIN 37450

typedef struct
{
  uint64_t u, v;
} pair64_t;

int is_overlap(uint32_t beg, uint32_t end, const bam1_t *b);

int pysam_bam_plbuf_push(const bam1_t *b, bam_plbuf_t *buf, int cont);

// accessor functions - necessary as bam_plbuf_t is hidden
// among the implementation
int pysam_get_pos( const bam_plbuf_t *buf);
int pysam_get_tid( const bam_plbuf_t *buf);
bam_pileup1_t * pysam_get_pileup( const bam_plbuf_t *buf);

int pysam_bam_fetch_init(bamFile fp,
			 const bam_index_t *idx,
			 int tid, int beg, int end,
			 pair64_t ** offp);

int pysam_bam_fetch_is_overlap(uint32_t beg, uint32_t end, const bam1_t *b);

int pysam_dispatch(int argc, char *argv[] );
  


#endif
