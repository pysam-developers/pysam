#include "samfile_util.h"
#include "htslib/sam.h"

// taken from bam_md.c
// replace bam1_{qual,seq,cigar} with bam_get_{qual,seq,cigar}
// bam1_seqi -> bam_seqi
// bam_nt16_table -> seq_nt16_table

#include <math.h>
#include <string.h>
#include <stdlib.h>

char bam_nt16_nt4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

 

