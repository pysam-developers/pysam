/*

   bam_md doesn't have a header file but we want to use its internal methods,
   so we explicitly expose them here, as it makes wrapping easier

*/

#ifndef BAM_MD_H
#define BAM_MD_H

#include <htslib/sam.h>

#ifdef __cplusplus
extern "C" {
#endif
    void bam_fillmd1_core(bam1_t *b, char *ref, int ref_len, int flag, int max_nm);
#ifdef __cplusplus
}
#endif

#endif // BAM_MD_H
