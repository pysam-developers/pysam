#ifndef HTSLIB_UTIL_H
#define HTSLIB_UTIL_H

#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/khash.h"

int hts_useek(htsFile *fp, long uoffset, int where);
long hts_utell(htsFile *fp);

int hts_set_verbosity(int verbosity);
int hts_get_verbosity();


KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

KHASH_DECLARE(s2i, kh_cstr_t, int64_t)
typedef khash_t(s2i) s2i_t;
		     
//////////////////////////////////////////////////////////////////
/*! set pysam standard error to point to file descriptor

  Setting the stderr will close the previous stderr.
 */
// FILE * pysam_set_stderr(int fd);

//////////////////////////////////////////////////////////////////
/*! set pysam standard error to /dev/null.
  
  Unsetting the stderr will close the previous stderr.
 */
// void pysam_unset_stderr();

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// various helper functions
//

/*!
  @abstract Update the variable length data within a bam1_t entry

  Old data is deleted and the data within b are re-arranged to 
  make place for new data.
  
  @discussion Returns b

  @param  b           bam1_t data
  @param  nbytes_old  size of old data
  @param  nbytes_new  size of new data
  @param  pos         position of data
*/
bam1_t * pysam_bam_update(bam1_t * b,
			  const size_t nbytes_old,
			  const size_t nbytes_new,
			  uint8_t * pos);

// translate a nucleotide character to binary code
unsigned char pysam_translate_sequence(const unsigned char s);

// return byte size of type
int aux_type2size(uint8_t type);


//-------------------------------------------------------
// Wrapping accessor macros in sam.h
static inline int pysam_bam_is_rev(bam1_t * b) {
  return bam_is_rev(b);};

static inline int pysam_bam_is_mrev(bam1_t * b) {
  return bam_is_mrev(b);}

static inline char * pysam_bam_get_qname(bam1_t * b) {
  return bam_get_qname(b);}

static inline uint32_t * pysam_bam_get_cigar(bam1_t * b) {
  return bam_get_cigar(b);}

static inline uint8_t * pysam_bam_get_seq(bam1_t * b) {
  return bam_get_seq(b);}

static inline uint8_t * pysam_bam_get_qual(bam1_t * b) {
  return bam_get_qual(b);}

static inline uint8_t * pysam_bam_get_aux(bam1_t * b) {
  return bam_get_aux(b);}

static inline int pysam_bam_get_l_aux(bam1_t * b) {
  return bam_get_l_aux(b); }

static inline char pysam_bam_seqi(uint8_t * s, int i) {
  return bam_seqi(s,i);}

// Wrapping bit field access in bam1_core_t
// bit fields not supported in cython and due
// to endian-ness it is not clear which part
// of the bit-field is in the higher or lower bytes.
static inline uint16_t pysam_get_bin(bam1_t * b) {
  return b->core.bin;}

static inline uint8_t pysam_get_qual(bam1_t * b) {
  return b->core.qual;}

static inline uint8_t pysam_get_l_qname(bam1_t * b) {
  return b->core.l_qname;}

static inline uint16_t pysam_get_flag(bam1_t * b) {
  return b->core.flag;}

static inline uint16_t pysam_get_n_cigar(bam1_t * b) {
  return b->core.n_cigar;}

static inline void pysam_set_bin(bam1_t * b, uint16_t v) {
  b->core.bin=v;}

static inline void pysam_set_qual(bam1_t * b, uint8_t v) {
  b->core.qual=v;}

static inline void pysam_set_l_qname(bam1_t * b, uint8_t v) {
  b->core.l_qname=v;}

static inline void pysam_set_flag(bam1_t * b, uint16_t v) {
  b->core.flag=v;}

static inline void pysam_set_n_cigar(bam1_t * b, uint16_t v) {
  b->core.n_cigar=v;}

static inline void pysam_update_flag(bam1_t * b, uint16_t v, uint16_t flag) {
  if (v)
    b->core.flag |= flag;
  else
    b->core.flag &= ~flag;
}

  
#endif
