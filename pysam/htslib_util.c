#include <assert.h>
#include "htslib/khash.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib_util.h"

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

// ASCII table of base complement characters
const char pysam_seq_comp_table[256] = {
    0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
    0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
    0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27, 0x28, 0x29, 0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f,
    0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f,
    0x40, 'T',  'V',  'G',  'H',  0x45, 0x46, 'C',  'D',  0x49, 0x4a, 'M',  0x4c, 'K',  'N',  0x4f,
    0x50, 0x51, 'Y',  'S',  'A',  'A',  'B',  'W',  0x58, 'R',  0x5a, 0x5b, 0x5c, 0x5d, 0x5e, 0x5f,
    0x60, 't',  'v',  'g',  'h',  0x65, 0x66, 'c',  'd',  0x69, 0x6a, 'm',  0x6c, 'k',  'n',  0x6f,
    0x70, 0x71, 'y',  's',  'a',  'a',  'b',  'w',  0x78, 'r',  0x7a, 0x7b, 0x7c, 0x7d, 0x7e, 0x7f,
    0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89, 0x8a, 0x8b, 0x8c, 0x8d, 0x8e, 0x8f,
    0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98, 0x99, 0x9a, 0x9b, 0x9c, 0x9d, 0x9e, 0x9f,
    0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7, 0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf,
    0xb0, 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7, 0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf,
    0xc0, 0xc1, 0xc2, 0xc3, 0xc4, 0xc5, 0xc6, 0xc7, 0xc8, 0xc9, 0xca, 0xcb, 0xcc, 0xcd, 0xce, 0xcf,
    0xd0, 0xd1, 0xd2, 0xd3, 0xd4, 0xd5, 0xd6, 0xd7, 0xd8, 0xd9, 0xda, 0xdb, 0xdc, 0xdd, 0xde, 0xdf,
    0xe0, 0xe1, 0xe2, 0xe3, 0xe4, 0xe5, 0xe6, 0xe7, 0xe8, 0xe9, 0xea, 0xeb, 0xec, 0xed, 0xee, 0xef,
    0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xff,
};
