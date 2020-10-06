/*  ext-sort.h -- sort on disk

   Copyright (C) 2020 Genome Research Ltd.

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

#ifndef __EXTSORT_H__
#define __EXTSORT_H__

//todo: return status to all functions

typedef struct _extsort_t extsort_t;

typedef int (*extsort_cmp_f) (const void *aptr, const void *bptr);

// Modes of operation
typedef enum
{
    DAT_SIZE,       // size_t        .. assuming constant size records for now
    TMP_PREFIX,     // const char*   .. prefix of temporary files, XXXXXX will be appended
    MAX_MEM,        // const char*   .. maximum memory to use, e.g. 100MB
    FUNC_CMP,       // extsort_cmp_f .. sort function
}
extsort_opt_t;

#define extsort_set_opt(es,type,key,value) { type tmp = value; extsort_set(es, key, (void*)&tmp); }

extsort_t *extsort_alloc(void);
void extsort_set(extsort_t *es, extsort_opt_t key, void *value);
void extsort_init(extsort_t *es);
void extsort_push(extsort_t *es, void *dat);    // dat will be freed by extsort later
void extsort_sort(extsort_t *es);
void *extsort_shift(extsort_t *es);
void extsort_destroy(extsort_t *es);

#endif
