/*
 * Copyright (c) 2022 Genome Research Ltd.
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

#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>

#include "utils.h"

#ifndef NO_THREADS
#include <pthread.h>
#endif

//#define TLS_DEBUG

#ifndef NO_THREADS
/*
 * Thread local storage per thread in the pool.
 *
 * We have some large memory blocks for rANS which we cannot store on the
 * stack due to various system limitations.  Allocaitng them can be
 * expensive as some OSes use mmap and will pass the pages back to the OS
 * on each free.  This unfortunately then means zeroing the pages out again
 * on each new malloc, plus additional switching into the kernel.
 *
 * Instead where available, we use pthread_once to allocate a small arena
 * of memory buffers and we continually reuse these same buffers.  We don't
 * need to memset it (calloc equivalent) either as we're sure that any
 * leakage of data is simply an earlier set of precomputed frequency
 * lookups, and not something more sinister such as an encryption key.
 *
 * If we don't have pthreads, then we have to fall back to the slow
 * traditional calloc instead.
 */

#define MAX_TLS_BUFS 10
typedef struct {
    void   *bufs[MAX_TLS_BUFS];
    size_t sizes[MAX_TLS_BUFS];
    int     used[MAX_TLS_BUFS];
} tls_pool;

static pthread_once_t rans_once = PTHREAD_ONCE_INIT;
static pthread_key_t rans_key;

/*
 * Frees all local storage for this thread.
 * Note: this isn't a function to free a specific allocated item.
 */
static void htscodecs_tls_free_all(void *ptr) {
    tls_pool *tls = (tls_pool *)ptr;
    if (!tls)
        return;

    int i;
    for (i = 0; i < MAX_TLS_BUFS; i++) {
#ifdef TLS_DEBUG
        if (tls->bufs[i])
            fprintf(stderr, "Free %ld = %p\n", tls->sizes[i], tls->bufs[i]);
#endif
        if (tls->used[i]) {
            fprintf(stderr, "Closing thread while TLS data is in use\n");
        }
        free(tls->bufs[i]);
    }

    free(tls);
}

static void htscodecs_tls_init(void) {
    pthread_key_create(&rans_key, htscodecs_tls_free_all);
}

/*
 * Allocates size bytes from the global Thread Local Storage pool.
 * This is shared by all subsequent calls within this thread.
 *
 * An simpler alternative could be possible where we have a fixed number
 * of types of alloc, say 5, and specify the correct label when allocating.
 * Eg histogram, name_context, fqzcomp, rans.  We can have multiple types
 * in use in different stack frames (such name_context + hist + rans), but
 * the number is very limited.  That then paves the way to simply check and
 * realloc without needing to keep track of use status or overflowing
 * the maximum number permitted.
 */
void *htscodecs_tls_alloc(size_t size) {
    int i;

    int err = pthread_once(&rans_once, htscodecs_tls_init);
    if (err != 0) {
        fprintf(stderr, "Initialising TLS data failed: pthread_once: %s\n",
                strerror(err));
        return NULL;
    }

    // Initialise tls_pool on first usage
    tls_pool *tls = pthread_getspecific(rans_key);
    if (!tls) {
        if (!(tls = calloc(1, sizeof(*tls))))
            return NULL;
        pthread_setspecific(rans_key, tls);
    }

    // Query pool for size
    int avail = -1;
    for (i = 0; i < MAX_TLS_BUFS; i++) {
        if (!tls->used[i]) {
            if (size <= tls->sizes[i]) {
                tls->used[i] = 1;
#ifdef TLS_DEBUG
                fprintf(stderr, "Reuse %d: %ld/%ld = %p\n",
                        i, size, tls->sizes[i], tls->bufs[i]);
#endif
                return tls->bufs[i];
            } else if (avail == -1) {
                avail = i;
            }
        }
    }

    if (i == MAX_TLS_BUFS && avail == -1) {
        // Shouldn't happen given our very limited use of this function
        fprintf(stderr, "Error: out of rans_tls_alloc slots\n");
        return NULL;
    }

    if (tls->bufs[avail])
        free(tls->bufs[avail]);
    if (!(tls->bufs[avail] = calloc(1, size)))
        return NULL;
#ifdef TLS_DEBUG
    fprintf(stderr, "Alloc %d: %ld = %p\n", avail, size, tls->bufs[avail]);
#endif
    tls->sizes[avail] = size;
    tls->used[avail] = 1;

    return tls->bufs[avail];
}

void *htscodecs_tls_calloc(size_t nmemb, size_t size) {
#ifdef TLS_DEBUG
    fprintf(stderr, "htscodecs_tls_calloc(%ld)\n", nmemb*size);
#endif
    void *ptr = htscodecs_tls_alloc(nmemb * size);
    if (ptr)
        memset(ptr, 0, nmemb * size);
    return ptr;
}

void htscodecs_tls_free(void *ptr) {
    if (!ptr)
        return;

    tls_pool *tls = pthread_getspecific(rans_key);

    int i;
    for (i = 0; i < MAX_TLS_BUFS; i++) {
        if (tls->bufs[i] == ptr)
            break;
    }
#ifdef TLS_DEBUG
    fprintf(stderr, "Fake free %d size %ld ptr %p\n",
            i, tls->sizes[i], tls->bufs[i]);
#endif
    if (i == MAX_TLS_BUFS) {
        fprintf(stderr, "Attempt to htscodecs_tls_free a buffer not allocated"
                " with htscodecs_tls_alloc\n");
        return;
    }
    if (!tls->used[i]) {
        fprintf(stderr, "Attempt to htscodecs_tls_free a buffer twice\n");
        return;
    }
    tls->used[i] = 0;
}

#else
/*
 * Calloc/free equivalents instead.
 *
 * We use calloc instead of malloc as a sufficiently malformed set of input
 * frequencies may not sum to the expected total frequency size, leaving
 * some elements uninitialised.  It's unlikely, but potentially a crafty
 * attacker could somehow exploit this to pull out parts of this allocated
 * buffer and leak them into the decompressed data stream, potentially
 * compromising previous buffers such as encryption keys.  (Although
 * frankly any well-written crypto library should be zeroing such memory
 * before freeing it to ensure it's never visible to a subsequent malloc.)
 */
void *htscodecs_tls_alloc(size_t size) {
    return calloc(1, size);
}

void *htscodecs_tls_calloc(size_t nmemb, size_t size) {
    return calloc(nmemb, size);
}

void htscodecs_tls_free(void *ptr) {
    free(ptr);
}
#endif
