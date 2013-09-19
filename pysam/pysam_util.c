#include <ctype.h>
#include <assert.h>
#include "bam.h"
#include "khash.h"
#include "ksort.h"
#include "bam_endian.h"
#include "knetfile.h"
#include "pysam_util.h"
#include "errmod.h" // for pysam_dump 

// for sequence parsing
#include "kseq.h"

#ifndef inline
#define inline __inline
#endif

// Definition of pysamerr
#include "stdio.h"
FILE * pysamerr = NULL;

FILE * pysam_set_stderr(int fd)
{
  if (pysamerr != NULL)
    fclose(pysamerr);
  pysamerr = fdopen(fd, "w");
  return pysamerr;
}

void pysam_unset_stderr()
{
  if (pysamerr != NULL)
    fclose(pysamerr);
  pysamerr = fopen("/dev/null", "w");
}

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
  
static mempool_t *mp_init()
{
	mempool_t *mp;
	mp = (mempool_t*)calloc(1, sizeof(mempool_t));
	return mp;
}
static void mp_destroy(mempool_t *mp)
{
	int k;
	for (k = 0; k < mp->n; ++k) {
		free(mp->buf[k]->b.data);
		free(mp->buf[k]);
	}
	free(mp->buf);
	free(mp);
}
static inline lbnode_t *mp_alloc(mempool_t *mp)
{
	++mp->cnt;
	if (mp->n == 0) return (lbnode_t*)calloc(1, sizeof(lbnode_t));
	else return mp->buf[--mp->n];
}
static inline void mp_free(mempool_t *mp, lbnode_t *p)
{
	--mp->cnt; p->next = 0; // clear lbnode_t::next here
	if (mp->n == mp->max) {
		mp->max = mp->max? mp->max<<1 : 256;
		mp->buf = (lbnode_t**)realloc(mp->buf, sizeof(lbnode_t*) * mp->max);
	}
	mp->buf[mp->n++] = p;
}

static inline int resolve_cigar(bam_pileup1_t *p, uint32_t pos)
{
	unsigned k;
	bam1_t *b = p->b;
	bam1_core_t *c = &b->core;
	uint32_t x = c->pos, y = 0;
	int ret = 1, is_restart = 1;

	if (c->flag&BAM_FUNMAP) return 0; // unmapped read
	assert(x <= pos); // otherwise a bug
	p->qpos = -1; p->indel = 0; p->is_del = p->is_head = p->is_tail = 0;
	for (k = 0; k < c->n_cigar; ++k) {
		int op = bam1_cigar(b)[k] & BAM_CIGAR_MASK; // operation
		int l = bam1_cigar(b)[k] >> BAM_CIGAR_SHIFT; // length
		if (op == BAM_CMATCH) { // NOTE: this assumes the first and the last operation MUST BE a match or a clip
			if (x + l > pos) { // overlap with pos
				p->indel = p->is_del = 0;
				p->qpos = y + (pos - x);
				if (x == pos && is_restart) p->is_head = 1;
				if (x + l - 1 == pos) { // come to the end of a match
					if (k < c->n_cigar - 1) { // there are additional operation(s)
						uint32_t cigar = bam1_cigar(b)[k+1]; // next CIGAR
						int op_next = cigar&BAM_CIGAR_MASK; // next CIGAR operation
						if (op_next == BAM_CDEL) p->indel = -(int32_t)(cigar>>BAM_CIGAR_SHIFT); // del
						else if (op_next == BAM_CINS) p->indel = cigar>>BAM_CIGAR_SHIFT; // ins
						if (op_next == BAM_CSOFT_CLIP || op_next == BAM_CREF_SKIP || op_next == BAM_CHARD_CLIP)
							p->is_tail = 1; // tail
					} else p->is_tail = 1; // this is the last operation; set tail
				}
			}
			x += l; y += l;
		} else if (op == BAM_CDEL) { // then set ->is_del
			if (x + l > pos) {
				p->indel = 0; p->is_del = 1;
				p->qpos = y + (pos - x);
			}
			x += l;
		} else if (op == BAM_CREF_SKIP) x += l;
		else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) y += l;
		is_restart = (op == BAM_CREF_SKIP || op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP);
		if (x > pos) {
			if (op == BAM_CREF_SKIP) ret = 0; // then do not put it into pileup at all
			break;
		}
	}
	assert(x > pos); // otherwise a bug
	return ret;

}
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

// Return number of mapped reads on tid.
// If tid < 0, return mapped reads without a coordinate (0)
uint32_t pysam_get_mapped( const bam_index_t *idx, const int tid )
{
  // return no values if index data not present
  if (idx==NULL) return 0;
  
  if (tid >= 0)
    {
      khint_t k;
      khash_t(i) *h = idx->index[tid];
      k = kh_get(i, h, BAM_MAX_BIN);

      if (k != kh_end(h))
	return kh_val(h, k).list[1].u;
      else
	return 0;
    }

  return 0;
}

uint32_t pysam_get_unmapped( const bam_index_t *idx, const int tid )
{

  if (tid >= 0)
    {
      khint_t k;
      khash_t(i) *h = idx->index[tid];
      k = kh_get(i, h, BAM_MAX_BIN);

      if (k != kh_end(h))
	return kh_val(h, k).list[1].v;
      else
	return 0;
    }

  return idx->n_no_coor;
}

/* uint32_t pysam_glf_depth( glf1_t * g ) */
/* { */
/*   return g->depth; */
/* } */


/* void pysam_dump_glf( glf1_t * g, bam_maqcns_t * c ) */
/* { */
/*   int x = 0; */
/*   fprintf(stderr, */
/* 	  "glf: ref_base=%i, max_mapQ=%i, min_lk=%i, depth=%i", */
/* 	  g->ref_base, */
/* 	  g->max_mapQ, */
/* 	  g->min_lk, */
/* 	  g->depth ); */

/*   for (x = 0; x < 10; ++x)  */
/*     fprintf(stderr, ", lk%x=%i, ", x, g->lk[x]); */

/*   fprintf(stderr, */
/* 	  "maqcns: het_rate=%f, theta=%f, n_hap=%i, cap_mapQ=%i, errmod=%i, min_baseQ=%i, eta=%f, q_r=%f, aux_max=%i", */
/* 	  c->het_rate, */
/* 	  c->theta, */
/* 	  c->n_hap, */
/* 	  c->cap_mapQ, */
/* 	  c->errmod, */
/* 	  c->min_baseQ, */
/* 	  c->eta, */
/* 	  c->q_r, */
/* 	  c->aux->max); */
  
/*   for (x = 0; x < c->aux->max; ++x) */
/*     { */
/*       fprintf(stderr, ", info-%i=%i ", x, c->aux->info[x]); */
/*       if (c->aux->info[x] == 0) break; */
/*     } */
  
/*   for (x = 0; x < c->aux->max; ++x) */
/*     { */
/*       fprintf(stderr, ", info16-%i=%i ", x, c->aux->info16[x]); */
/*       if (c->aux->info16[x] == 0) break; */
/*     } */
/* } */


  

// pysam dispatch function to emulate the samtools
// command line within python.
// taken from the main function in bamtk.c
// added code to reset getopt
int bam_taf2baf(int argc, char *argv[]);
int bam_mpileup(int argc, char *argv[]);
int bam_merge(int argc, char *argv[]);
int bam_index(int argc, char *argv[]);
int bam_sort(int argc, char *argv[]);
int bam_tview_main(int argc, char *argv[]);
int bam_mating(int argc, char *argv[]);
int bam_rmdup(int argc, char *argv[]);
int bam_flagstat(int argc, char *argv[]);
int bam_fillmd(int argc, char *argv[]);
int bam_idxstats(int argc, char *argv[]);
int main_samview(int argc, char *argv[]);
int main_import(int argc, char *argv[]);
int main_reheader(int argc, char *argv[]);
int main_cut_target(int argc, char *argv[]);
int main_phase(int argc, char *argv[]);
int main_cat(int argc, char *argv[]);
int main_depth(int argc, char *argv[]);
int main_bam2fq(int argc, char *argv[]);
int main_pad2unpad(int argc, char *argv[]);
int main_bedcov(int argc, char *argv[]);
int main_bamshuf(int argc, char *argv[]);

int faidx_main(int argc, char *argv[]);


int pysam_dispatch(int argc, char *argv[] )
{
  extern int optind;
#ifdef _WIN32
  setmode(fileno(stdout), O_BINARY);
  setmode(fileno(stdin),  O_BINARY);
#ifdef _USE_KNETFILE
  knet_win32_init();
#endif
#endif

  // reset getopt
  optind = 1;

  if (argc < 2) return 1;
  int retval = 0;
  
  if (strcmp(argv[1], "view") == 0) retval = main_samview(argc-1, argv+1);
  else if (strcmp(argv[1], "import") == 0) retval = main_import(argc-1, argv+1);
  else if (strcmp(argv[1], "mpileup") == 0) retval = bam_mpileup(argc-1, argv+1);
  else if (strcmp(argv[1], "merge") == 0) retval = bam_merge(argc-1, argv+1);
  else if (strcmp(argv[1], "sort") == 0) retval = bam_sort(argc-1, argv+1);
  else if (strcmp(argv[1], "index") == 0) retval = bam_index(argc-1, argv+1);
  else if (strcmp(argv[1], "faidx") == 0) retval = faidx_main(argc-1, argv+1);
  else if (strcmp(argv[1], "idxstats") == 0) retval = bam_idxstats(argc-1, argv+1);
  else if (strcmp(argv[1], "fixmate") == 0) retval = bam_mating(argc-1, argv+1);
  else if (strcmp(argv[1], "rmdup") == 0) retval = bam_rmdup(argc-1, argv+1);
  else if (strcmp(argv[1], "flagstat") == 0) retval = bam_flagstat(argc-1, argv+1);
  else if (strcmp(argv[1], "calmd") == 0) retval = bam_fillmd(argc-1, argv+1);
  else if (strcmp(argv[1], "fillmd") == 0) retval = bam_fillmd(argc-1, argv+1);
  else if (strcmp(argv[1], "reheader") == 0) retval = main_reheader(argc-1, argv+1);
  else if (strcmp(argv[1], "cat") == 0) retval = main_cat(argc-1, argv+1);
  else if (strcmp(argv[1], "targetcut") == 0) retval = main_cut_target(argc-1, argv+1);
  else if (strcmp(argv[1], "phase") == 0) retval = main_phase(argc-1, argv+1);
  else if (strcmp(argv[1], "depth") == 0) retval = main_depth(argc-1, argv+1);
  else if (strcmp(argv[1], "bam2fq") == 0) retval = main_bam2fq(argc-1, argv+1);
  else if (strcmp(argv[1], "pad2unpad") == 0) retval = main_pad2unpad(argc-1, argv+1);
  else if (strcmp(argv[1], "depad") == 0) retval = main_pad2unpad(argc-1, argv+1);
  else if (strcmp(argv[1], "bedcov") == 0) retval = main_bedcov(argc-1, argv+1);
  else if (strcmp(argv[1], "bamshuf") == 0) retval = main_bamshuf(argc-1, argv+1);
  
#if _CURSES_LIB != 0
  else if (strcmp(argv[1], "tview") == 0) retval = bam_tview_main(argc-1, argv+1);
#endif
  else 
    {
      fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
      return 1;
    }
  fflush( stdout );
  
  return retval;
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
bam1_t * pysam_bam_update( bam1_t * b,
			   const size_t nbytes_old,
			   const size_t nbytes_new, 
			   uint8_t * pos )
{
  int d = nbytes_new-nbytes_old;
  int new_size;
  size_t offset;

  // no change
  if (d == 0) return b;

  new_size = d + b->data_len;
  offset = pos - b->data;

  //printf("d=%i, old=%i, new=%i, old_size=%i, new_size=%i\n",
  // d, nbytes_old, nbytes_new, b->data_len, new_size);
  
  // increase memory if required
  if (d > 0)
    {
      alloc_data( b, new_size );
      pos = b->data + offset;
    }
  
  if (b->data_len != 0)
    {
      if (offset < 0 || offset > b->data_len)
	fprintf(stderr, "[pysam_bam_insert] illegal offset: '%i'\n", (int)offset);
    }
  
  // printf("dest=%p, src=%p, n=%i\n", pos+nbytes_new, pos + nbytes_old, b->data_len - (offset+nbytes_old));
  memmove( pos + nbytes_new,
	   pos + nbytes_old,
	   b->data_len - (offset + nbytes_old));
    
  b->data_len = new_size;
      
  return b;
}

// translate a nucleotide character to binary code
unsigned char pysam_translate_sequence( const unsigned char s )
{
  return bam_nt16_table[s];
}


void bam_init_header_hash(bam_header_t *header);

// translate a reference string *s* to a tid
// code taken from bam_parse_region
int pysam_reference2tid( bam_header_t *header, const char * s )
{
  
  khiter_t iter;
  khash_t(s) *h;
  
  bam_init_header_hash(header);
  h = (khash_t(s)*)header->hash;

  iter = kh_get(s, h, s); /* get the ref_id */
  if (iter == kh_end(h)) { // name not found
    return -1;
  }

  return kh_value(h, iter);
}

// Auxiliary functions for B support
void bam_aux_appendB(bam1_t *b, const char tag[2], char type, char subtype, int len, uint8_t *data)
{

  int ori_len;

  int data_len;

  // check that type is 'B'
  if('B' != type) return;

  ori_len = b->data_len;

  data_len = len * bam_aux_type2size(subtype);
  // infer the data length from the sub-type
  b->data_len += 8 + data_len;

  b->l_aux += 8 + data_len;

  if (b->m_data < b->data_len) 
    {

      b->m_data = b->data_len;

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
  memcpy(b->data + ori_len + 8, data, data_len);
  // data
}

/*
// return size of auxiliary type
int bam_aux_type2size(int x)
{
  if (x == 'C' || x == 'c' || x == 'A') return 1;
  else if (x == 'S' || x == 's') return 2;
  else if (x == 'I' || x == 'i' || x == 'f') return 4;
  else return 0;
}
*/





