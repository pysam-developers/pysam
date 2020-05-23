/* 
    Copyright (C) 2019 Genome Research Ltd.
    
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
    DEALINGS IN THE SOFTWARE.
*/

#ifndef __COLS_H__
#define __COLS_H__

#include <stdlib.h>

typedef struct
{
    int n,m;
    char **off, *rmme;
}
cols_t;

/*
    cols_split() can be called repeatedly to split new strings, memory is allocated
    and deallocated automatically
*/
cols_t *cols_split(const char *line, cols_t *cols, char delim);

/* 
    Although cols_append() can be combined with cols_split(), it is much slower and
    the string must exist throughout the life of cols unless initialized with cols_split().
*/
void cols_append(cols_t *cols, char *str);
void cols_clear(cols_t *cols);
void cols_destroy(cols_t *cols);

#endif
