#ifndef PYSAM_UTIL_H
#define PYSAM_UTIL_H

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


#endif
