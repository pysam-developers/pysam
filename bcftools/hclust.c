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

#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <stdlib.h>
#include <assert.h>
#include "bcftools.h"
#include "hclust.h"

typedef struct _node_t
{
    struct _node_t *akid, *bkid, *next, *prev, *parent;
    int id, idx;    // id: unique node id; idx: current index to pdist
    float value;    // max pairwise dist of elements within the node
}
node_t;

struct _hclust_t
{
    int ndat, nclust;       // ndat: number of elements (pdist matrix size); nclust: current number of clusters
    float *pdist;           // pairwise cluster distances, diagonal matrix accessed via the PDIST macro
    node_t *first, *last;   // clusters are maintained in a double-linked list
    node_t **rmme;          // convenience array to remove all allocated nodes at the end
    int nrmme;
    kstring_t str;          // (for debugging) pointer to str.s is returned by create_dot()
    char **dbg;             // (for debugging) created by create_list() via set_threshold() and returned by explain()
    int ndbg, mdbg;
};

node_t *append_node(hclust_t *clust, int idx)
{
    node_t *node = (node_t*) calloc(1,sizeof(node_t));

    clust->nclust++;
    node->id  = clust->nrmme;
    node->idx = idx;
    if ( !clust->first )
    {
        clust->first = node; 
        clust->last  = node; 
    }
    else
    {
        node->prev = clust->last;
        clust->last->next = node; 
        clust->last = node; 
    }
    
    if ( clust->nrmme >= clust->ndat*2 ) error("hclust fixme: %d vs %d\n",clust->nrmme,clust->ndat);
    clust->rmme[clust->nrmme++] = node;

    return node;
}
void remove_node(hclust_t *clust, node_t *node)
{
    if ( node==clust->first ) clust->first = node->next;
    if ( node==clust->last ) clust->last = node->prev;
    if ( node->next ) node->next->prev = node->prev;
    if ( node->prev ) node->prev->next = node->next;
    clust->nclust--;
}

#if DEBUG
void hclust_debug(hclust_t *clust)
{
    int i;
    fprintf(stderr,"nrmme=%d  nclust=%d\n", clust->nrmme,clust->nclust);
    for (i=0; i<clust->nrmme; i++)
    {
        node_t *node = clust->rmme[i];
        int akid  = node->akid ? node->akid->id : -1;
        int bkid  = node->bkid ? node->bkid->id : -1;
        int akidx = node->akid ? node->akid->idx : -1;
        int bkidx = node->bkid ? node->bkid->idx : -1;
        fprintf(stderr,"\t%d\t%d\t%f\t%d %d\t%d %d\n",node->id,node->idx,node->value,akid,bkid,akidx,bkidx);
    }

    int j;
    for (i=1; i<clust->ndat; i++)
    {
        int active = 0;
        node_t *node = clust->first;
        while (node)
        {
            if ( node->idx==i ) { active=1; break; }
            node = node->next;
        }
        fprintf(stderr,"%2d%c ",i,active?'*':' ');
        for (j=0; j<i; j++)
        {
            if ( PDIST(clust->pdist,i,j)==9 )
                fprintf(stderr,"  -----  ");
            else
                fprintf(stderr," %f", PDIST(clust->pdist,i,j));
        }
        fprintf(stderr,"\n");
    }
    for (j=0; j<clust->ndat-1; j++) fprintf(stderr,"  %6d ",j); fprintf(stderr,"\n");
}
#endif

hclust_t *hclust_init(int n, float *pdist)
{
    hclust_t *clust = (hclust_t*) calloc(1,sizeof(hclust_t));
    clust->ndat  = n;
    clust->pdist = pdist;
    clust->rmme = (node_t**) calloc(n*2,sizeof(node_t*));

    // init clusters
    int i;
    for (i=0; i<clust->ndat; i++) append_node(clust,i);

    // build the tree
    while ( clust->nclust>1 )
    {
        // find two clusters with minimum distance
        float min_value = HUGE_VAL;
        node_t *iclust = clust->first->next;
        node_t *min_iclust = NULL, *min_jclust = NULL;
        while ( iclust )
        {
            node_t *jclust = clust->first;
            while ( jclust!=iclust )
            {
                float value = PDIST(clust->pdist,iclust->idx,jclust->idx);
                if ( value < min_value ) 
                { 
                    min_value  = value;
                    min_iclust = iclust;
                    min_jclust = jclust; 
                }
                jclust = jclust->next;
            }
            iclust = iclust->next;
        }
        assert( min_iclust && min_jclust ); // pdist contains inf or nan, fix the caller
        remove_node(clust,min_iclust);
        remove_node(clust,min_jclust);

        // update the pairwise distances. We keep the matrix and as we are moving up the
        // tree, we use fewer columns/rows as the number of clusters decreases: we reuse
        // i-th and leave j-th unused. Inter-cluster distance is defined as maximum distance
        // between pairwise distances of elements within the cluster.
        iclust = clust->first;
        while ( iclust )
        {
            if ( PDIST(clust->pdist,iclust->idx,min_iclust->idx) < PDIST(clust->pdist,iclust->idx,min_jclust->idx) )
                PDIST(clust->pdist,iclust->idx,min_iclust->idx) = PDIST(clust->pdist,iclust->idx,min_jclust->idx);
            iclust = iclust->next;
        }

        node_t *node = append_node(clust,min_iclust->idx);
        node->akid  = min_iclust;
        node->bkid  = min_jclust;
        node->value = min_value;
        node->akid->parent = node;
        node->bkid->parent = node;
    }

    return clust;
}
void hclust_destroy(hclust_t *clust)
{
    int i;
    for (i=0; i<clust->nrmme; i++) free(clust->rmme[i]);
    free(clust->rmme);
    free(clust->dbg);
    free(clust->str.s);
    free(clust);
}

char *hclust_create_dot(hclust_t *clust, char **labels, float th)
{
    clust->str.l = 0;
    ksprintf(&clust->str,"digraph myGraph {");

    int i;
    for (i=0; i<clust->nrmme; i++)
    {
        node_t *node = clust->rmme[i];
        if ( node->value )
            ksprintf(&clust->str,"\"%d\" [label=\"%f\"];", node->id,node->value);
        else
            ksprintf(&clust->str,"\"%d\" [label=\"%s\"];", node->id,labels[node->idx]);
    }
    for (i=0; i<clust->nrmme; i++)
    {
        node_t *node = clust->rmme[i];
        if ( node->akid )
        {
            if ( node->value >= th && node->akid && node->akid->value < th )
                ksprintf(&clust->str,"\"%d\" -> \"%d\" [color=\"#D43F3A\" penwidth=3];", node->id,node->akid->id);
            else
                ksprintf(&clust->str,"\"%d\" -> \"%d\";", node->id,node->akid->id);
        }

        if ( node->bkid )
        {
            if ( node->value >= th && node->bkid && node->bkid->value < th )
                ksprintf(&clust->str,"\"%d\" -> \"%d\" [color=\"#D43F3A\" penwidth=3];", node->id,node->bkid->id);
            else
                ksprintf(&clust->str,"\"%d\" -> \"%d\";", node->id,node->bkid->id);
        }
    }
    ksprintf(&clust->str,"};");
    return clust->str.s;
}
char **hclust_explain(hclust_t *clust, int *nlines)
{
    clust->ndbg = 0;
    char *beg = clust->str.s;
    while ( *beg )
    {
        char *end = beg;
        while ( *end && *end!='\n' ) end++;
        clust->ndbg++;
        hts_expand(char*,clust->ndbg,clust->mdbg,clust->dbg);
        clust->dbg[clust->ndbg-1] = beg;
        if ( !*end ) break;
        *end = 0;
        beg = end + 1;
    }

    *nlines = clust->ndbg;
    return clust->dbg;
}

cluster_t *append_cluster(node_t *node, cluster_t *cluster, int *nclust, node_t **stack)
{
    (*nclust)++;
    cluster = (cluster_t*) realloc(cluster,sizeof(cluster_t)*(*nclust));
    cluster_t *clust = &cluster[*nclust-1];
    clust->nmemb = 0;
    clust->memb  = NULL;
    clust->dist  = node->value;

    int nstack = 1;
    stack[0] = node;

    while ( nstack )
    {
        node_t *node = stack[--nstack];
        node_t *akid = node->akid;
        node_t *bkid = node->bkid;
        if ( node->akid )
        {
            stack[nstack++] = akid;
            stack[nstack++] = bkid;
        }
        else    
        {
            clust->nmemb++;
            clust->memb = (int*) realloc(clust->memb,sizeof(int)*clust->nmemb);
            clust->memb[clust->nmemb-1] = node->id;
        }
    }
    return cluster;
}

int cmp_nodes(const void *a, const void *b)
{
    const node_t *an = *((const node_t**) a);
    const node_t *bn = *((const node_t**) b);
    if ( an->value < bn->value ) return -1;
    if ( an->value > bn->value ) return 1;
    return 0;
}

float calc_dev(node_t **dat, int n)
{
    float avg = 0, dev = 0;
    int i;
    for (i=0; i<n; i++) avg += dat[i]->value;
    avg /= n;
    for (i=0; i<n; i++) dev += (dat[i]->value - avg)*(dat[i]->value - avg);
    return sqrt(dev/n);
}

/*
    Heuristics to determine clustering cutoff: sort nodes by distance and
    split into two groups by minimizing the standard deviation.
    This works best when two elements from a single different sample are
    included in the mix.
        - min_inter_dist .. smaller values are always considered identical
        - max_intra_dist .. larger values are always considered different
 */
float hclust_set_threshold(hclust_t *clust, float min_inter_dist, float max_intra_dist)
{
    node_t **dat = clust->rmme + clust->ndat;
    int i, ndat = clust->nrmme - clust->ndat;
 
    qsort(dat, ndat, sizeof(dat), cmp_nodes);

    clust->str.l = 0;
    float th, min_dev = HUGE_VAL;
    int imin = -1;
    for (i=0; i<ndat; i++)
    {
        float dev = 0;
        if ( i>0 ) dev += calc_dev(dat,i);
        if ( i+1<ndat ) dev += calc_dev(dat+i,ndat-i);
        th  = dat[i]->value;
        ksprintf(&clust->str,"DEV\t%f\t%f\n",th,dev);
        if ( min_dev > dev && th >= min_inter_dist ) { min_dev = dev; imin = i; }
    }
    if ( max_intra_dist > 0 )
        th = max_intra_dist;  // use fixed cutoff, the above was only for debugging output
    else
    {
        // dynamic cutoff
        max_intra_dist = fabs(max_intra_dist);
        th = imin==-1 ? max_intra_dist : dat[imin]->value;
        if ( th > max_intra_dist ) th = max_intra_dist;
    }
    ksprintf(&clust->str,"TH\t%f\n", th);
    ksprintf(&clust->str,"MAX_DIST\t%f\n", dat[ndat-1]->value);
    ksprintf(&clust->str,"MIN_INTER\t%f\n", min_inter_dist);
    ksprintf(&clust->str,"MAX_INTRA\t%f\n", max_intra_dist);
    return th;
} 

cluster_t *hclust_create_list(hclust_t *clust, float min_inter_dist, float *max_intra_dist, int *nclust)
{
    float cutoff = *max_intra_dist = hclust_set_threshold(clust, min_inter_dist, *max_intra_dist);

    node_t **stack = (node_t**) malloc(sizeof(node_t*)*clust->ndat);
    node_t **tmp = (node_t**) malloc(sizeof(node_t*)*clust->ndat);
    stack[0] = clust->first;
    int nstack = 1;
    
    cluster_t *cluster = NULL;
    int ncluster = 0;

    if ( stack[0]->value < cutoff )
    {
        // all values are within the limits - create a single cluster
        cluster = append_cluster(stack[0], cluster, &ncluster, tmp);
        nstack = 0;
    }

    while ( nstack )
    {
        node_t *node = stack[--nstack];
        node_t *akid = node->akid;
        node_t *bkid = node->bkid;
        if ( !akid )
        {
            cluster = append_cluster(node, cluster, &ncluster, tmp);
            continue;
        }

        if ( node->value >= cutoff && akid->value < cutoff )
            cluster = append_cluster(akid, cluster, &ncluster, tmp);
        else    
            stack[nstack++] = akid;

        if ( node->value >= cutoff && bkid->value < cutoff )
            cluster = append_cluster(bkid, cluster, &ncluster, tmp);
        else    
            stack[nstack++] = bkid;
    }

    free(tmp);
    free(stack);

    *nclust = ncluster;
    return cluster;
}

void hclust_destroy_list(cluster_t *clust, int nclust)
{
    int i;
    for (i=0; i<nclust; i++) free(clust[i].memb);
    free(clust);
}


