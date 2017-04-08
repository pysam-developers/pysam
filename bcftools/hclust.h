/* The MIT License

   Copyright (c) 2016 Genome Research Ltd.

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
    Simple hierarchical clustering
*/

#ifndef __HCLUST_H__
#define __HCLUST_H__

#include <stdio.h>

typedef struct _hclust_t hclust_t;

typedef struct
{
    float dist;
    int nmemb, *memb;
}
cluster_t;

#define PDIST(mat,a,b) (mat)[((a)>(b)?((a)*((a)-1)/2+(b)):((b)*((b)-1)/2+(a)))]

/*
 *  hclust_init() - init and run clustering
 *  @n:     number of elements
 *  @pdist: pairwise distances. The array will be modified by hclust and
 *          must exist until hclust_destroy() is called
 */
hclust_t *hclust_init(int n, float *pdist);
void hclust_destroy(hclust_t *clust);

/*
 *  hclust_create_list() - returns a list of clusters
 *  @min_inter_dist: minimum inter-cluster distance. If smaller, elements are considered
 *                   homogenous, belonging to the same cluster.
 *  @max_intra_dist: maximum intra-cluster distance allowed. If smaller than 0,
 *                   the threshold can be heuristically lowered, otherwise considered
 *                   a fixed cutoff. The pointer will be filled to the cutoff actually used.
 */
cluster_t *hclust_create_list(hclust_t *clust, float min_inter_dist, float *max_intra_dist, int *nclust);
void hclust_destroy_list(cluster_t *clust, int nclust);

/* 
 *  Access debugging data used in the decision making process.  Note that this
 *  must be called immediately after hclust_create_list because other calls,
 *  such as hclust_create_dot(), invalidate the temporary data structures.
 */
char **hclust_explain(hclust_t *clust, int *nlines);

char *hclust_create_dot(hclust_t *clust, char **labels, float th);

#endif

