/*
 * Copyright (c) 2021-2024 Genome Research Ltd.
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

#ifndef HTSCODECS_H
#define HTSCODECS_H

/*
 * Version X.Y.Z encoded as XYYYZZ.
 * We mainly increment X and Y.  Z *may* get bumped in between official
 * releases in order to distinguish untagged github checkouts from
 * official release tarballs.
 *
 * Note currently this needs manually editing as it isn't automatically
 * updated by autoconf.
 */
#define HTSCODECS_VERSION 100601

/*
 * A const string form of the HTSCODECS_VERSION define.
 * NB: This is obtained from the auto-generated version.h, so
 * we can include release number and git hash.
 */
const char *htscodecs_version(void);

#endif /* HTSCODECS_H */
