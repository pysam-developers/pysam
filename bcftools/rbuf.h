/*  rbuf.h -- round buffers.

    Copyright (C) 2013-2014, 2017 Genome Research Ltd.

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
THE SOFTWARE.  */

#ifndef __RBUF_H__
#define __RBUF_H__

#include <string.h>

typedef struct
{
    int m,n,f;    // m: allocated size, n: number of elements in the buffer, f: first element
}
rbuf_t;

/**
 *  rbuf_init() - initialize round buffer
 *  @rbuf:  the rbuf_t holder
 *  @size:  the maximum number of elements
 *
 */
static inline void rbuf_init(rbuf_t *rbuf, int size)
{
    rbuf->m = size; rbuf->n = rbuf->f = 0;
}
/**
 *  rbuf_kth() - get index of the k-th element of the round buffer
 *  @rbuf:  the rbuf_t holder
 *  @k:     0-based index. If negative, return k-th element from the end, 1-based
 */
static inline int rbuf_kth(rbuf_t *rbuf, int k)
{
    if ( k >= rbuf->n ) return -1;
    if ( k < 0 )
    {
        k = rbuf->n + k; 
        if ( k < 0 ) return -1;
    }
    int i = k + rbuf->f;
    if ( i >= rbuf->m ) i -= rbuf->m;
    return i;
}
/**
 *  rbuf_last() - get index of the last element of the round buffer
 *  @rbuf:  the rbuf_t holder
 */
#define rbuf_last(rbuf) rbuf_kth(rbuf, -1)

/**
 *  rbuf_l2ridx() - get 0-based rbuf index which corresponds to i-th linear index
 *  @rbuf:  the rbuf_t holder
 *  @idx:   0-based linear index
 *
 *  Returns 0-based circular index or -1 if out of bounds
 */
static inline int rbuf_l2ridx(rbuf_t *rbuf, int idx)
{
    if ( idx < 0 || idx >= rbuf->n ) return -1;
    if ( idx >= rbuf->f )
    {
        int i = idx - rbuf->f;
        if ( i >= rbuf->n ) return -1;
        return i;
    }
    int i = rbuf->m - rbuf->f + idx;
    if ( i >= rbuf->n ) return -1;
    return i;
}

/**
 *  rbuf_next() - get index of the next element in the round buffer
 *  @rbuf:  the rbuf_t holder
 *  @i:     pointer to the last rbuf index. Set to -1 before the first call.
 *
 *  Sets i to the next position in the buffer. The return value indicates if
 *  the position points to a valid element (1) or if there are no more elements
 *  after *i (0). When the end is reached, *i is set to the first element in the
 *  buffer.
 */
static inline int rbuf_next(rbuf_t *rbuf, int *i)
{
    if ( !rbuf->n ) return 0;
    if ( *i==-1 ) { *i = rbuf->f; return 1; }
    int n = (rbuf->f <= *i) ? *i - rbuf->f + 1 : *i + rbuf->m - rbuf->f + 1;
    if ( ++(*i) >= rbuf->m ) *i = 0;
    if ( n < rbuf->n ) return 1;
    *i = rbuf->f;
    return 0;
}
/**
 *  rbuf_prev() - get index of the previous element in the round buffer
 *  @rbuf:  the rbuf_t holder
 *  @i:     pointer to the last rbuf index. Set to -1 before the first call.
 *
 *  Sets i to the previous position in the buffer. The return value indicates if
 *  the position points to a valid element (1) or if there are no more elements
 *  before *i (0).
 */
static inline int rbuf_prev(rbuf_t *rbuf, int *i)
{
    if ( !rbuf->n || *i==rbuf->f ) return 0;
    if ( *i==-1 )
    {
        *i = rbuf_last(rbuf);
        return 1;
    }
    if ( --(*i) < 0 ) *i = rbuf->m - 1;
    return 1;
}
/**
 *  rbuf_prepend() - register new element at the start of the round buffer
 *  @rbuf:  the rbuf_t holder
 *
 *  Returns index of the newly inserted element.
 */
static inline int rbuf_prepend(rbuf_t *rbuf)
{
    if ( rbuf->n < rbuf->m ) rbuf->n++;

    rbuf->f = rbuf->f > 0 ? rbuf->f - 1 : rbuf->m - 1;
    return rbuf->f;
}
/**
 *  rbuf_append() - register new element at the end of the round buffer
 *  @rbuf:  the rbuf_t holder
 *
 *  Returns index of the newly inserted element.
 */
static inline int rbuf_append(rbuf_t *rbuf)
{
    if ( rbuf->n < rbuf->m )
    {
        rbuf->n++;
        int i = rbuf->f + rbuf->n;
        return i <= rbuf->m ? i - 1 : i - rbuf->m - 1;
    }

    rbuf->f++;
    if ( rbuf->f >= rbuf->m )
    {
        rbuf->f = 0;
        return rbuf->m - 1;
    }
    return rbuf->f - 1;
}
/**
 *  rbuf_shift() - removes first element from the buffer
 *  @rbuf:  the rbuf_t holder
 *
 *  Returns index of the removed element.
 */
static inline int rbuf_shift(rbuf_t *rbuf)
{
    if ( !rbuf->n ) return -1;
    int ret = rbuf->f;
    rbuf->f++;
    if ( rbuf->f >= rbuf->m ) rbuf->f = 0;
    rbuf->n--;
    return ret;
}
/**
 *  rbuf_shift_n() - removes first n elements from the buffer
 *  @rbuf:  the rbuf_t holder
 *  @n:     number of elements to remove
 */
static inline void rbuf_shift_n(rbuf_t *rbuf, int n)
{
    if ( n >= rbuf->n )
    {
        rbuf->n = rbuf->f = 0;
        return;
    }
    rbuf->n -= n;
    rbuf->f += n;
    if ( rbuf->f >= rbuf->m ) rbuf->f -= rbuf->m;
}

/**
 *  rbuf_expand0() - expand round buffer and set the newly allocated elements to 0
 *  @rbuf:      the rbuf holder
 *  @type_t:    data type
 *  @n:         requested number of elements
 *  @data:      data array to be realloced
 *
 *  Note: The new array is linearized and leaves the rbuf.f offset untouched,
 *  thus the size of the new buffer is determined by the current position.
 */
#define rbuf_expand0(rbuf,type_t,n,data) \
{ \
    if ( n > (rbuf)->m ) \
    { \
        int m = n + (rbuf)->f; \
        m--, m|=m>>1, m|=m>>2, m|=m>>4, m|=m>>8, m|=m>>16, m++; /* kroundup32 */ \
        data = (type_t*) realloc(data, sizeof(type_t)*m); \
        type_t *ptr = data; \
        memset(ptr+(rbuf)->m,0,sizeof(type_t)*(m-(rbuf)->m)); \
        if ( (rbuf)->f ) \
        { \
            memcpy(ptr+(rbuf)->m,ptr,sizeof(type_t)*(rbuf)->f); \
            memset(ptr,0,sizeof(type_t)*(rbuf)->f); \
        } \
        (rbuf)->m = m; \
    } \
}

/**
 *  rbuf_remove_kth() - remove k-th rbuf element (0-based) and memmove the data block
 *  @rbuf:      the rbuf holder
 *  @type_t:    data type
 *  @k:         k-th element to remove
 *  @data:      data array to be modified
 */
#define rbuf_remove_kth(rbuf, type_t, kth, data) \
{ \
    int k = rbuf_kth(rbuf, kth); \
    if ( k < (rbuf)->f )    /* shrink from back */ \
    { \
        int l = rbuf_kth(rbuf, -1); \
        if ( k < l ) \
        { \
            type_t tmp = (data)[k]; \
            memmove(data+k, data+k+1, (l - k)*sizeof(type_t)); \
            (data)[l] = tmp; \
        } \
        (rbuf)->n--; \
    } \
    else                    /* shrink from front */ \
    { \
        if ( k > (rbuf)->f ) \
        { \
            type_t tmp = (data)[k]; \
            memmove(&data[(rbuf)->f+1], &data[(rbuf)->f], (k - (rbuf)->f)*sizeof(type_t)); \
            (data)[(rbuf)->f] = tmp; \
        } \
        (rbuf)->f++; \
        (rbuf)->n--; \
        if ( (rbuf)->f == (rbuf)->m ) (rbuf)->f = 0; \
    } \
}

#endif
