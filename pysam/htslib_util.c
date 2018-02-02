#include <ctype.h>
#include <assert.h>
#include "htslib/khash.h"
#include "htslib/ksort.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/knetfile.h"
#include "htslib/kseq.h"
#include "htslib_util.h"
#include <stdio.h>

#ifndef inline
#define inline __inline
#endif

// set htslib verbosity level
extern int hts_verbose;
int hts_set_verbosity(int verbosity)
{
  int old_verbosity = hts_verbose;
  hts_verbose = verbosity;
  return old_verbosity;
}

int hts_get_verbosity(void)
{
  return hts_verbose;
}


// taken from samtools/bam_import.c
static inline uint8_t * alloc_data(bam1_t *b, size_t size)
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
// Return NULL on error (memory allocation)
bam1_t * pysam_bam_update(bam1_t * b,
			  const size_t nbytes_old,
			  const size_t nbytes_new, 
			  uint8_t * field_start)
{
  int d = nbytes_new - nbytes_old;
  int new_size;
  size_t nbytes_before;
  uint8_t * retval = NULL;
    
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
      retval = alloc_data(b, new_size);
      if (retval == NULL)
	return NULL;
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
