#ifndef PYSAM_UTIL_H
#define PYSAM_UTIL_H

#include "kseq.h"

// #######################################################
// fastq parsing
KSEQ_INIT(gzFile, gzread)

//////////////////////////////////////////////////////////////////
/*! set pysam standard error to point to file descriptor

  Setting the stderr will close the previous stderr.
 */
FILE * pysam_set_stderr( int fd );

//////////////////////////////////////////////////////////////////
/*! set pysam standard error to /dev/null.
  
  Unsetting the stderr will close the previous stderr.
 */
void pysam_unset_stderr();

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
// various helper functions
//
// fill pileup buffer for next position.

int pysam_pileup_next(const bam1_t *b,
		      bam_plbuf_t *buf,
		      bam_pileup1_t ** plp,
		      int * tid,
		      int * pos,
		      int * n_plp);

int pysam_dispatch(int argc, char *argv[] );

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
bam1_t * pysam_bam_update( bam1_t * b,
			   const size_t nbytes_old,
			   const size_t nbytes_new,
			   uint8_t * pos );

// translate a nucleotide character to binary code
unsigned char pysam_translate_sequence( const unsigned char s );

// defined in bam_import.c
extern unsigned char bam_nt16_table[256];

// translate a reference string *s* to a tid
int pysam_reference2tid( bam_header_t *header, const char * s );

// return number of mapped reads for tid
uint32_t pysam_get_mapped( const bam_index_t *idx, const int tid );

// return number of unmapped reads for tid
uint32_t pysam_get_unmapped( const bam_index_t *idx, const int tid );

// debugging functions
/* #include "glf.h" */
/* uint32_t pysam_glf_depth( glf1_t * g); */

/* #include "bam_maqcns.h" */
/* void pysam_dump_glf( glf1_t * g, bam_maqcns_t * c ); */


// return size of auxilliary type
// int bam_aux_type2size(int x);

#endif
