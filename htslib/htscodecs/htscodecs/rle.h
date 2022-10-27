/*
 * Copyright (c) 2019 Genome Research Ltd.
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

#ifndef HTS_RLE_H
#define HTS_RLE_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Performs run length encoding of a byte stream, turning it into a
 * list of lengths and a list of literals.
 *
 * The method used is a bit different to traditional run length
 * encoding.  It always outputs run-lengths for symbols in the
 * 'rle_syms' list (even if that length is +0 more), and never outputs
 * lengths for symbols not in that list.
 *
 * "run" should be preallocated to be large enough;
 * e.g at least data_len bytes long as a worse case.
 * "rle_syms" should be allocated to be at least 256 bytes.
 *
 * If *rle_nsyms is zero this function will survey the input data
 * first to choose symbols automatically, writing back to rle_syms and
 * rle_nsyms.
 *
 * The "out" buffer may be passed in as NULL in which case it is
 * allocated and returned (and is up to the caller to free).
 * Otherwise if specified as non-NULL it will be written to, but
 * it is up to the caller to ensure the buffer size is large enough.
 * A worst case scenario is 2*data_len.
 *
 * Returns the literal buffer on success with new length in out_len,
 *         also fills out run buffer and run_len,  and potentially
 *         updates rle_syms / rle_nsyms too.
 * Returns NULL of failure
 */
uint8_t *hts_rle_encode(uint8_t *data, uint64_t data_len,
                        uint8_t *run,  uint64_t *run_len,
                        uint8_t *rle_syms, int *rle_nsyms,
                        uint8_t *out, uint64_t *out_len);

/*
 * Expands a run lengthed data steam from a pair of literal and
 * run-length buffers.
 *
 * On input *out_len holds the length of the supplied out
 * buffer.  On exit, it holds the used portion of this buffer.
 *
 * Returns uncompressed data (out) on success,
 *         NULL on failure.
 */
uint8_t *hts_rle_decode(uint8_t *lit, uint64_t lit_len,
                        uint8_t *run, uint64_t run_len,
                        uint8_t *rle_syms, int rle_nsyms,
                        uint8_t *out, uint64_t *out_len);

// TODO: Add rle scanning func to compute rle_syms.

#ifdef __cplusplus
}
#endif

#endif /* HTS_RLE_H */
