/*  test/test-rbuf.c -- rbuf_t test harness.

    Copyright (C) 2014 Genome Research Ltd.

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <stdlib.h>
#include "rbuf.h"

void debug_print(rbuf_t *rbuf, int *dat)
{
    int i;
    for (i=-1; rbuf_next(rbuf, &i); ) printf(" %2d", i);
    printf("\n");
    for (i=-1; rbuf_next(rbuf, &i); ) printf(" %2d", dat[i]);
    printf("\n");
}

int main(int argc, char **argv)
{
    int i, j, *dat = (int*)calloc(10,sizeof(int));
    rbuf_t rbuf;
    rbuf_init(&rbuf,10);

    rbuf.f = 5; // force wrapping
    for (i=0; i<9; i++)
    {
        j = rbuf_append(&rbuf);
        dat[j] = i+1;
    }
    printf("Inserted 1-9 starting at offset 5:\n");
    debug_print(&rbuf, dat);

    i = rbuf_kth(&rbuf, 3);
    printf("4th is %d\n", dat[i]);

    printf("Deleting 1-2:\n");
    rbuf_shift_n(&rbuf, 2);
    debug_print(&rbuf, dat);

    printf("Prepending 0-8:\n");
    for (i=0; i<9; i++)
    {
        j = rbuf_prepend(&rbuf);
        dat[j] = i;
    }
    debug_print(&rbuf, dat);

    printf("Expanding:\n");
    rbuf_expand0(&rbuf,int,rbuf.n+1,dat);
    debug_print(&rbuf, dat);

    free(dat);
    return 0;
}

