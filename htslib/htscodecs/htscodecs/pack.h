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

#ifndef HTS_PACK_H
#define HTS_PACK_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Packs multiple symbols into a single byte if the total alphabet of symbols
 * used is <= 16.  Each new symbol takes up 1, 2, 4 or 8 bits, or 0 if the
 * alphabet used is 1 (constant).
 *
 * If successful, out_meta/out_meta_len are set to hold the mapping table
 * to be used during decompression.
 *
 * Returns the packed buffer on success with new length in out_len,
 *         NULL of failure
 */
uint8_t *hts_pack(uint8_t *data, int64_t len,
                  uint8_t *out_meta, int *out_meta_len, uint64_t *out_len);

/*
 * Unpacks the meta-data portions of the hts_pack algorithm.
 * This consists of the count of symbols and their values.
 *
 * The "map" array is filled out with the used symbols.
 * "nsym" is set to contain the number of symbols per byte;
 * 0, 1, 2, 4 or 8.
 *
 * Returns number of bytes of data[] consumed on success,
 *         zero on failure.
 */
uint8_t hts_unpack_meta(uint8_t *data, uint32_t data_len,
                        uint64_t udata_len, uint8_t *map, int *nsym);

/*
 * Unpacks a packed data steam (given the unpacked meta-data).
 *
 * "map" is the pack map, mapping 0->n to the expanded symbols.
 * The "out" buffer must be preallocated by the caller to be the correct
 * size.  For error checking purposes, out_len is set to the size of
 * this buffer.
 *
 * Returns uncompressed data (out) on success,
 *         NULL on failure.
 */
uint8_t *hts_unpack(uint8_t *data, int64_t len, uint8_t *out, uint64_t out_len, int nsym, uint8_t *map);

#ifdef __cplusplus
}
#endif

#endif /* HTS_PACK_H */
