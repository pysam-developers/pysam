#include <ctype.h>
#include <assert.h>
#include "khash.h"
#include "ksort.h"
#include "sam.h"
#include "hts.h"
#include "knetfile.h"
#include "htslib_util.h"
#include <stdio.h>

// for sequence parsing
#include "kseq.h"

#ifndef inline
#define inline __inline
#endif

// #######################################################
// utility routines to avoid using callbacks in bam_fetch
// taken from bam_index.c
// The order of the following declarations is important.
// #######################################################
#define BAM_MAX_BIN 37450 // =(8^6-1)/7+1

// initialize hashes
typedef struct
{
  uint64_t u, v;
} pair64_t;

#define pair64_lt(a,b) ((a).u < (b).u)

KSORT_INIT(myoff, pair64_t, pair64_lt);

typedef struct {
	uint32_t m, n;
	pair64_t *list;
} bam_binlist_t;

typedef struct {
	int32_t n, m;
	uint64_t *offset;
} bam_lidx_t;


// initialize hashes ('i' and 's' are idenditifiers)
KHASH_MAP_INIT_INT(i, bam_binlist_t);
KHASH_MAP_INIT_STR(s, int)

/*
struct __bam_index_t
{
  int32_t n;
  uint64_t n_no_coor; // unmapped reads without coordinate
  khash_t(i) **index;
  bam_lidx_t *index2;
};

typedef struct __linkbuf_t {
	bam1_t b;
	uint32_t beg, end;
	struct __linkbuf_t *next;
} lbnode_t;

typedef struct {
	int cnt, n, max;
	lbnode_t **buf;
} mempool_t;

struct __bam_plbuf_t {
	mempool_t *mp;
	lbnode_t *head, *tail, *dummy;
	bam_pileup_f func;
	void *func_data;
	int32_t tid, pos, max_tid, max_pos;
	int max_pu, is_eof;
	bam_pileup1_t *pu;
	int flag_mask;
};
  
// the following code has been taken from bam_plbuf_push
// and modified such that instead of a function call
// the function returns and will continue (if cont is true).
// from where it left off.

// returns
// 1: if buf is full and can be emitted
// 0: if b has been added
// -1: if there was an error
int pysam_pileup_next(const bam1_t *b,
		      bam_plbuf_t *buf,
		      bam_pileup1_t ** plp,
		      int * tid,
		      int * pos,
		      int * n_plp )
{
  *plp = bam_plp_next(buf->iter, tid, pos, n_plp);
  if (plp == NULL) return 0;
  return 1;
}

typedef struct __bmc_aux_t {
	int max;
	uint32_t *info;
	uint16_t *info16;
	errmod_t *em;
} bmc_aux_t;
*/

// Return number of mapped reads on tid.
// If tid < 0, return mapped reads without a coordinate (0)
uint32_t pysam_get_mapped(const hts_idx_t *idx, const int tid)
{
  // return no values if index data not present
  if (idx==NULL) return 0;

  // TODO
  /* if (tid >= 0) */
  /*   { */
  /*     khint_t k; */
  /*     khash_t(i) *h = idx->index[tid]; */
  /*     k = kh_get(i, h, BAM_MAX_BIN); */

  /*     if (k != kh_end(h)) */
  /* 	return kh_val(h, k).list[1].u; */
  /*     else */
  /* 	return 0; */
  /*   } */

  return 0;
}

uint32_t pysam_get_unmapped( const hts_idx_t *idx, const int tid)
{
  // TODO
  return 0;
  
  /* if (tid >= 0) */
  /*   { */
  /*     khint_t k; */
  /*     khash_t(i) *h = idx->index[tid]; */
  /*     k = kh_get(i, h, BAM_MAX_BIN); */

  /*     if (k != kh_end(h)) */
  /* 	return kh_val(h, k).list[1].v; */
  /*     else */
  /* 	return 0; */
  /*   } */

  /* return idx->n_no_coor; */
}

// taken from samtools/bam_import.c
static inline uint8_t *alloc_data(bam1_t *b, size_t size)
{
  if (b->m_data < size)
    {
      b->m_data = size;
      kroundup32(b->m_data);
      b->data = (uint8_t*)realloc(b->data, b->m_data);
    }
  return b->data;
}

// update the variable length data within a bam1_t entry.
// Adds *nbytes_new* - *nbytes_old* into the variable length data of *src* at *pos*.
// Data within the bam1_t entry is moved so that it is
// consistent with the data field lengths.
bam1_t * pysam_bam_update(bam1_t * b,
			  const size_t nbytes_old,
			  const size_t nbytes_new, 
			  uint8_t * field_start)
{
  int d = nbytes_new - nbytes_old;
  int new_size;
  size_t nbytes_before;

  // no change
  if (d == 0)
    return b;

  // new size of total data
  new_size = d + b->l_data;
  
  // fields before field in data
  nbytes_before = field_start - b->data;

  if (b->l_data != 0)
    {
      assert(nbytes_before >= 0);
      assert(nbytes_before <= b->l_data);
    }

  // increase memory if required
  if (d > 0)
    {
      alloc_data(b, new_size);
      field_start = b->data + nbytes_before;
    }
  
  // move data after field to new location
  memmove(field_start + nbytes_new,
	  field_start + nbytes_old,
	  b->l_data - (nbytes_before + nbytes_old));

  // adjust l_data
  b->l_data = new_size;
      
  return b;
}

// translate a nucleotide character to binary code
unsigned char pysam_translate_sequence(const unsigned char s)
{
  return seq_nt16_table[s];
}

// Auxiliary functions for B support
void bam_aux_appendB(bam1_t *b,
		     const char tag[2],
		     char type,
		     char subtype,
		     int len,
		     uint8_t *data)
{

  int ori_len;
  int l_data;

  // check that type is 'B'
  if('B' != type) return;

  ori_len = b->l_data;

  l_data = len * aux_type2size(subtype);
  // infer the data length from the sub-type
  b->l_data += 8 + l_data;

  // htslib: obsolete?
  // b->l_aux += 8 + l_data;
  if (b->m_data < b->l_data) 
    {
      b->m_data = b->l_data;
      kroundup32(b->m_data);
      b->data = (uint8_t*)realloc(b->data, b->m_data);
    }

  b->data[ori_len] = tag[0];
  b->data[ori_len + 1] = tag[1];
  // tag
  b->data[ori_len + 2] = type;
  // type
  b->data[ori_len + 3] = subtype;
  // subtype
  (*(int32_t*)(b->data + ori_len + 4)) = len;
  // size
  memcpy(b->data + ori_len + 8, data, l_data);
  // data
}

int aux_type2size(uint8_t type)
{
	switch (type) {
	case 'A': case 'c': case 'C':
		return 1;
	case 's': case 'S':
		return 2;
	case 'i': case 'I': case 'f':
		return 4;
	case 'd':
		return 8;
	case 'Z': case 'H': case 'B':
		return type;
	default:
		return 0;
	}
}



