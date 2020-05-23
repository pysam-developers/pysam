#include "bcftools.pysam.h"

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

#include <string.h>
#include "cols.h"

cols_t *cols_split(const char *line, cols_t *cols, char delim)
{
    if ( !cols ) cols = (cols_t*) calloc(1,sizeof(cols_t));
    if ( cols->rmme ) free(cols->rmme);
    cols->n = 0;
    cols->rmme = strdup(line);
    char *ss = cols->rmme;
    while (1)
    {
        char *se = ss;
        while ( *se && *se!=delim ) se++;
        char tmp = *se;
        *se = 0;
        cols->n++;
        if ( cols->n > cols->m )
        {
            cols->m += 10;
            cols->off = (char**) realloc(cols->off, sizeof(*cols->off)*cols->m);
        }
        cols->off[ cols->n - 1 ] = ss;
        if ( !tmp ) break;
        ss = se + 1;
    }
    return cols;
}

void cols_append(cols_t *cols, char *str)
{
    if ( cols->rmme )
    {
        size_t str_len = strlen(str);
        size_t lst_len = strlen(cols->off[ cols->n - 1 ]);
        size_t tot_len = 2 + str_len + lst_len + (cols->off[ cols->n - 1 ] - cols->rmme);

        cols_t *tmp_cols = (cols_t*)calloc(1,sizeof(cols_t));
        tmp_cols->rmme = (char*) calloc(tot_len,1);
        tmp_cols->off  = (char**) calloc(cols->n+1,sizeof(*tmp_cols->off));

        char *ptr = tmp_cols->rmme;
        int i;
        for (i=0; i<cols->n; i++)
        {
            size_t len = strlen(cols->off[i]);
            memcpy(ptr, cols->off[i], len);
            tmp_cols->off[i] = ptr;
            ptr += len + 1;
        }
        memcpy(ptr, str, str_len);
        tmp_cols->off[i] = ptr;

        free(cols->off);
        free(cols->rmme);
        cols->rmme = tmp_cols->rmme;
        cols->off  = tmp_cols->off;
        cols->n    = cols->n+1;
        cols->m    = cols->n;
        free(tmp_cols);
        return;
    }
    cols->n++;
    if ( cols->n > cols->m )
    {
        cols->m++;
        cols->off = (char**) realloc(cols->off,sizeof(*cols->off)*cols->m);
    }
    cols->off[cols->n-1] = str;
}
void cols_clear(cols_t *cols)
{
    if ( !cols ) return;
    free(cols->rmme);
    free(cols->off);
    cols->rmme = NULL;
    cols->off  = NULL;
}
void cols_destroy(cols_t *cols)
{
    if ( !cols ) return;
    cols_clear(cols);
    free(cols);
}

