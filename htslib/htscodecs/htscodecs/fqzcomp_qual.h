/*
 * Copyright (c) 2011-2013, 2018-2019 Genome Research Ltd.
 * Author(s): James Bonfield
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
 *       Institute nor the names of its contributors may be used to endorse
 *       or promote products derived from this software without specific
 *       prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
 * LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef FQZ_COMP_QUAL_H
#define FQZ_COMP_QUAL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>

/* Bit flags, deliberately mirroring BAM ones */
#define FQZ_FREVERSE 16
#define FQZ_FREAD2 128

/* Current FQZ format version */
#define FQZ_VERS 5

#define FQZ_MAX_STRAT 3

/*
 * Minimal per-record information taken from a cram slice.
 *
 * To compress we need to know the junction from one quality string to
 * the next (len), whether it is first/second read and whether it is
 * reverse complemented (flags).
 */
typedef struct {
    int num_records;
    uint32_t *len;    // of size num_records
    uint32_t *flags;  // of size num_records
} fqz_slice;


// Global flags
static const int GFLAG_MULTI_PARAM = 1;
static const int GFLAG_HAVE_STAB   = 2;
static const int GFLAG_DO_REV      = 4;

// Param flags
// Add PFLAG_HAVE_DMAP and a dmap[] for delta incr?
static const int PFLAG_DO_DEDUP    = 2;
static const int PFLAG_DO_LEN      = 4;
static const int PFLAG_DO_SEL      = 8;
static const int PFLAG_HAVE_QMAP   = 16;
static const int PFLAG_HAVE_PTAB   = 32;
static const int PFLAG_HAVE_DTAB   = 64;
static const int PFLAG_HAVE_QTAB   = 128;

/*
 * FQZ parameters.  These may be simply passed in as NULL to fqz_compress
 * and it'll automatically choose, but if we wish to have complete control
 * then this (long) struct contains all the details.
 *
 * TODO: document all this!
 */

// A single parameter block
typedef struct {
    // Starting context value
    uint16_t context;

    // flags
    unsigned int pflags;
    unsigned int do_sel, do_dedup, store_qmap, fixed_len;
    unsigned char use_qtab, use_dtab, use_ptab;

    // context bits and locations
    unsigned int qbits, qloc;
    unsigned int pbits, ploc;
    unsigned int dbits, dloc;
    unsigned int sbits, sloc;

    // models
    int max_sym, nsym, max_sel;

    // tables / maps
    unsigned int qmap[256];
    unsigned int qtab[256];
    unsigned int ptab[1024];
    unsigned int dtab[256];

    // Not stored paramters, but computed as part of encoder
    // parameterisation.
    int qshift;
    int pshift;
    int dshift;
    int sshift;
    unsigned int qmask; // (1<<qbits)-1
    int do_r2, do_qa;
} fqz_param;

// The global params, which is a collection of parameter blocks plus
// a few pieces of meta-data.
typedef struct {
    int vers;               // Format version; Set to FQZ_VERS
    unsigned int gflags;    // global param flags
    int nparam;             // Number of fqz_param blocks
    int max_sel;            // Number of selector values
    unsigned int stab[256]; // Selector to parameter no. table

    int max_sym;            // max symbol value across all sub-params

    fqz_param *p;           // 1 or more parameter blocks
} fqz_gparams;


/** Compress a block of quality values.
 *
 * @param vers          The CRAM version number (<<8) plus fqz strategy (0-3)
 * @param s             Length and flag data CRAM per-record
 * @param in            Buffer of concatenated quality values (no separator)
 * @param in_size       Size of in buffer
 * @param out_size      Size of returned output
 * @param strat         FQZ compression strategy (0 to FQZ_MAX_STRAT)
 * @param gp            Optional fqzcomp paramters (may be NULL).
 *
 * @return              The compressed quality buffer on success,
 *                      NULL on failure.
 */
char *fqz_compress(int vers, fqz_slice *s, char *in, size_t in_size,
                   size_t *out_size, int strat, fqz_gparams *gp);

/** Decompress a block of quality values.
 *
 * @param in            Buffer of compressed quality values
 * @param in_size       Size of in buffer
 * @param out_size      Size of returned output
 * @param lengths       Optional array filled out with record lengths.
 *                      May be NULL.  If not, preallocate it to correct size.
 *
 * @return              The uncompressed concatenated qualities on success,
 *                      NULL on failure.
 */
char *fqz_decompress(char *in, size_t in_size, size_t *out_size,
                     int *lengths, int nlengths);

/** A utlity function to analyse a quality buffer to gather statistical
 *  information.  This is written into qhist and pm.  This function is only
 *  useful if you intend on passing your own fqz_gparams block to
 *  fqz_compress.
 */
void fqz_qual_stats(fqz_slice *s,
                    unsigned char *in, size_t in_size,
                    fqz_param *pm,
                    uint32_t qhist[256],
                    int one_param);

#ifdef __cplusplus
}
#endif

#endif
