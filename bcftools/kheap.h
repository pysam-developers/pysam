/* The MIT License

   Copyright (C) 2016 Genome Research Ltd.

   Author: Petr Danecek <pd3@sanger.ac.uk>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:
   
   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.
   
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

 */
/*
    Usage example:

        #include "kheap.h"

        // First we prepare the user data to store, in this example it is a
        // struct with a single element "key", and a comparator function
        // "is_smaller".  In this example the comparator defines a min heap (as
        // opposed to a max heap).
        typedef struct
        {
            uint32_t key; 
        } 
        data_t;
        static inline int is_smaller(data_t *a, data_t *b)
        {
            return a->key < b->key ? 1 : 0; 
        }
        data_t data[3] = { {3}, {2}, {1} };


        // Heap declaration, "mh" is an arbitrary string.  The typedef is not
        // required, it is just a convenience shortcut so that we can use
        // "heap_t" instead of the generic "khp_mh_t" automatically created by
        // the KHEAP_INIT macro.
        KHEAP_INIT(mh, data_t, is_smaller)
        typedef khp_mh_t heap_t;

        // Initialize the heap, insert the test data, then retrieve them back,
        // sorted. Multiple heaps with the same name "mh" can be created and
        // used simultaneously, as long as they all use the same data type
        // "data_t".
        heap_t *heap = khp_init(mh);

        // When inserting a new element, the heap stores a copy of the memory
        // area pointed to by the third argument.
        for (int i=0; i<3; i++)
            khp_insert(mh, heap, &data[i]);

        while (heap->ndat)
        {
            printf("%d\n", heap->dat[0].pos);
            khp_delete(mh, heap);
        }

        // Clean up
        khp_destroy(mh, heap);

*/

#ifndef __KHEAP_H__
#define __KHEAP_H__

#include <stdlib.h>

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#ifndef kh_inline
#ifdef _MSC_VER
#define kh_inline __inline
#else
#define kh_inline inline
#endif
#endif /* kh_inline */

#ifndef klib_unused
#if (defined __clang__ && __clang_major__ >= 3) || (defined __GNUC__ && __GNUC__ >= 3)
#define klib_unused __attribute__ ((__unused__))
#else
#define klib_unused
#endif
#endif /* klib_unused */


#define __KHEAP_TYPE(name, kheap_t) \
    typedef struct {                \
        int ndat, mdat;             \
        kheap_t *dat;               \
        kheap_t tmp;                \
    } khp_##name##_t;

#define khp_parent(i) (((i)-1)/2)
#define khp_lchild(i) (2*(i)+1)
#define khp_rchild(i) (2*(i)+2)
#define khp_swap(hp,i,j) {               \
        ((hp)->tmp)    = ((hp)->dat[i]); \
        ((hp)->dat[i]) = ((hp)->dat[j]); \
        ((hp)->dat[j]) = ((hp)->tmp);    \
    }

#define __KHEAP_IMPL(name, SCOPE, kheap_t, __cmp)                       \
    SCOPE khp_##name##_t *khp_init_##name(void)                         \
    {                                                                   \
        return (khp_##name##_t*)calloc(1, sizeof(khp_##name##_t));      \
    }                                                                   \
    SCOPE void khp_destroy_##name(khp_##name##_t *heap)                 \
    {                                                                   \
        if (heap) free(heap->dat);                                      \
        free(heap);                                                     \
    }                                                                   \
    SCOPE int khp_insert_##name(khp_##name##_t *heap, kheap_t *dat)     \
    {                                                                   \
        heap->ndat++;                                                   \
        if ( heap->ndat > heap->mdat )                                  \
        {                                                               \
            heap->mdat = heap->ndat;                                    \
            kroundup32(heap->mdat);                                     \
            heap->dat = (kheap_t*)realloc(heap->dat, heap->mdat*sizeof(kheap_t));         \
            memset(heap->dat + heap->ndat, 0, (heap->mdat - heap->ndat)*sizeof(kheap_t)); \
        }                                                               \
        int i = heap->ndat - 1;                                         \
        while ( i && __cmp(dat,&heap->dat[khp_parent(i)]) )             \
        {                                                               \
            heap->dat[i] = heap->dat[khp_parent(i)];                    \
            i = khp_parent(i);                                          \
        }                                                               \
        heap->dat[i] = *dat;                                            \
        return i;                                                       \
    }                                                                   \
    SCOPE void khp_heapify_##name(khp_##name##_t *heap, int i)          \
    {                                                                   \
/*todo: loop instead of a recursive function? */ \
        int extreme = khp_lchild(i) < heap->ndat && __cmp(&heap->dat[khp_lchild(i)],&heap->dat[i]) ? khp_lchild(i) : i;     \
        if ( khp_rchild(i) < heap->ndat && __cmp(&heap->dat[khp_rchild(i)],&heap->dat[extreme]) ) extreme = khp_rchild(i);  \
        if ( extreme != i )                                             \
        {                                                               \
            khp_swap(heap,i,extreme);                                   \
            khp_heapify_##name(heap,extreme);                           \
        }                                                               \
    }                                                                   \
    SCOPE void khp_delete_##name(khp_##name##_t *heap)                  \
    {                                                                   \
        if ( !heap || !heap->ndat ) return;                             \
        heap->dat[0] = heap->dat[--heap->ndat];                         \
        khp_heapify_##name(heap, 0);                                    \
    }                                                                   \

#define KHEAP_INIT(name, kheap_t, __cmp)            \
    __KHEAP_TYPE(name, kheap_t)                     \
    __KHEAP_IMPL(name, static kh_inline klib_unused, kheap_t, __cmp)

#define khp_init(name) khp_init_##name()
#define khp_destroy(name, heap) khp_destroy_##name(heap)
#define khp_insert(name, heap, dat) khp_insert_##name(heap, dat)
#define khp_delete(name, heap) khp_delete_##name(heap)

#endif
