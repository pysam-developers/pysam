/*
 * Copyright (c) 2017, 2019 Genome Research Ltd.
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
 *    Institute nor the names of its contributors may be used to endorse
 *    or promote products derived from this software without specific
 *    prior written permission.
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

#ifndef _TOKENISE_NAME3_H_
#define _TOKENISE_NAME3_H_

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Converts a line or \0 separated block of reading names to a compressed buffer.
 * The code can only encode whole lines and will not attempt a partial line.
 * Use the "last_start_p" return value to identify the partial line start
 * offset, for continuation purposes.
 *
 * Returns a malloced buffer holding compressed data of size *out_len,
 *         or NULL on failure
 */
uint8_t *tok3_encode_names(char *blk, int len, int level, int use_arith,
                           int *out_len, int *last_start_p);

/*
 * Decodes a compressed block of read names into \0 separated names.
 * The size of the data returned (malloced) is in *out_len.
 *
 * Returns NULL on failure.
 */
uint8_t *tok3_decode_names(uint8_t *in, uint32_t sz, uint32_t *out_len);

#ifdef __cplusplus
}
#endif

#endif /* _TOKENISE_NAME3_H_ */
