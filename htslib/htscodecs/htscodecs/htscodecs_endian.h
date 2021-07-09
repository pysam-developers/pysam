/*
 * Copyright (c) 2020 Genome Research Ltd.
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

#ifndef HTSCODECS_ENDIAN_H
#define HTSCODECS_ENDIAN_H

// Endianness checking.

// Sets HTSCODECS_ENDIAN_KNOWN if system type detected and either
// HTSCODECS_LITTLE_ENDIAN or HTSCODECS_BIG_ENDIAN.

/*
 * In general our preferred route is to write code in an endian agnostic
 * fashion, but our data formats are natively little endian.  Therefore
 * in time critical code it's sometimes best to exploit that.
 *
 * Therefore we'll optimise code along the lines of:
 *
 * #ifdef HTSCODECS_LITTLE_ENDIAN
 *     // do something little endian specific
 * #else
 *     // do something in an endian agnostic fashion
 * #endif
 *
 * This means our code works even if we fail to recognise the
 * specific machine type.
 */

#if (defined(__i386__)      \
 ||  defined(__i386)        \
 ||  defined(__amd64__)     \
 ||  defined(__amd64)       \
 ||  defined(__x86_64__)    \
 ||  defined(__x86_64)      \
 ||  defined(__i686__)      \
 ||  defined(__i686))       \
 || (defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__) \
 ||  defined(__LITTLE_ENDIAN__)                                            \
 ||  defined(__ARMEL__)     \
 ||  defined(__THUMBEL__)   \
 ||  defined(__AARCH64EL__) \
 ||  defined(_MIPSEL)       \
 ||  defined(__MIPSEL)      \
 ||  defined(__MIPSEL__) 
    // Little endian
#   define HTSCODECS_LITTLE_ENDIAN
#   define HTSCODECS_ENDIAN_KNOWN
#elif (defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__) \
   ||  defined(__BIG_ENDIAN__) \
   || defined(__ARMEB__)       \
   || defined(__THUMBEB__)     \
   || defined(__AAARCHEB__)    \
   || defined(_MIPSEB)         \
   || defined(__MIPSEB)        \
   || defined(__MIPSEB__)
    // Big endian
#   define HTSCODECS_BIG_ENDIAN
#   define HTSCODECS_ENDIAN_KNOWN
#else
//    Unknown - code will need to check HTSCODES_ENDIAN_KNOWN and do endian agnostic
#endif

#endif /* HTSCODECS_ENDIAN_H */
