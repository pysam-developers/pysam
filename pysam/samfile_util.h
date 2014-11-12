#ifndef SAMFILE_UTIL_H
#define SAMFILE_UTIL_H

#include "htslib/sam.h"

int bam_cap_mapQ(bam1_t *b, char *ref, int thres);
int bam_prob_realn(bam1_t *b, const char *ref);

#endif

