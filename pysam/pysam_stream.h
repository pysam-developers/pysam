#ifndef PYSAM_STREAM_H
#define PYSAM_STREAM_H

#include "htslib/kseq.h"

// #######################################################
// fastq parsing
KSEQ_INIT(gzFile, gzread)

//KSTREAM_INIT( gzFile, gzread, 16384)

#endif
