/* The MIT License

   Copyright (c) 2015 Genome Research Ltd.

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <errno.h>
#include "bcftools.h"

#define MODE_COUNT     1
#define MODE_LIST_GOOD 2
#define MODE_LIST_BAD  4
#define MODE_DELETE    8

typedef struct
{
    int nok, nbad;
    int imother,ifather,ichild;
}
trio_t;

typedef struct _args_t
{
    bcf_hdr_t *hdr;
    int32_t *gt_arr;
    int mode;
    int ngt_arr, nrec;
    trio_t *trios;
    int ntrios;
}
args_t;

static args_t args;

const char *about(void)
{
    return "Count Mendelian consistent / inconsistent genotypes.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Count Mendelian consistent / inconsistent genotypes.\n"
        "Usage: bcftools +mendelian [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        "Plugin options:\n"
        "   -c, --count             count the number of consistent sites\n"
        "   -d, --delete            delete inconsistent genotypes (set to \"./.\")\n"
        "   -l, --list [+x]         list consistent (+) or inconsistent (x) sites\n"
        "   -t, --trio <m,f,c>      names of mother, father and the child\n"
        "   -T, --trio-file <file>  list of trios, one per line\n"
        "\n"
        "Example:\n"
        "   bcftools +mendelian in.vcf -- -t Mother,Father,Child -c\n"
        "\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    char *trio_samples = NULL, *trio_file = NULL;
    memset(&args,0,sizeof(args_t));
    args.hdr  = in;
    args.mode = 0;

    static struct option loptions[] =
    {
        {"trio",1,0,'t'},
        {"trio-file",1,0,'T'},
        {"delete",0,0,'d'},
        {"list",1,0,'l'},
        {"count",0,0,'c'},
        {0,0,0,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "?ht:T:l:cd",loptions,NULL)) >= 0)
    {
        switch (c) 
        {
            case 'd': args.mode |= MODE_DELETE; break;
            case 'c': args.mode |= MODE_COUNT; break;
            case 'l': 
                if ( !strcmp("+",optarg) ) args.mode |= MODE_LIST_GOOD; 
                else if ( !strcmp("x",optarg) ) args.mode |= MODE_LIST_BAD; 
                else error("The argument not recognised: --list %s\n", optarg);
                break;
            case 't': trio_samples = optarg; break;
            case 'T': trio_file = optarg; break;
            case 'h':
            case '?':
            default: error(usage()); break;
        }
    }
    if ( optind != argc ) error(usage());
    if ( !trio_samples && !trio_file ) error("Expected the -t/T option\n");
    if ( !args.mode ) error("Expected one of the -c, -d or -l options\n");
    if ( args.mode&MODE_DELETE && !(args.mode&(MODE_LIST_GOOD|MODE_LIST_BAD)) ) args.mode |= MODE_LIST_GOOD|MODE_LIST_BAD;

    int i, n = 0;
    char **list;
    if ( trio_samples )
    {
        args.ntrios = 1;
        args.trios = (trio_t*) calloc(1,sizeof(trio_t));
        list = hts_readlist(trio_samples, 0, &n);
        if ( n!=3 ) error("Expected three sample names with -t\n");
        args.trios[0].imother = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[0]);
        args.trios[0].ifather = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[1]);
        args.trios[0].ichild  = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, list[2]);
        for (i=0; i<n; i++) free(list[i]);
        free(list);
    }
    if ( trio_file )
    {
        list = hts_readlist(trio_file, 1, &n);
        args.ntrios = n;
        args.trios = (trio_t*) calloc(n,sizeof(trio_t));
        for (i=0; i<n; i++)
        {
            char *ss = list[i], *se;
            se = strchr(ss, ',');
            if ( !se ) error("Could not parse %s: %s\n",trio_file, ss);
            *se = 0;
            args.trios[i].imother = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, ss);
            if ( args.trios[i].imother<0 ) error("No such sample: \"%s\"\n", ss);
            ss = ++se; 
            se = strchr(ss, ',');
            if ( !se ) error("Could not parse %s\n",trio_file);
            *se = 0;
            args.trios[i].ifather = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, ss);
            if ( args.trios[i].ifather<0 ) error("No such sample: \"%s\"\n", ss);
            ss = ++se; 
            if ( *ss=='\0' ) error("Could not parse %s\n",trio_file);
            args.trios[i].ichild = bcf_hdr_id2int(args.hdr, BCF_DT_SAMPLE, ss);
            if ( args.trios[i].ichild<0 ) error("No such sample: \"%s\"\n", ss);
            free(list[i]);
        }
        free(list);
    }
    return args.mode&(MODE_LIST_GOOD|MODE_LIST_BAD) ? 0 : 1;
}

bcf1_t *process(bcf1_t *rec)
{
    bcf1_t *dflt = args.mode&MODE_LIST_GOOD ? rec : NULL;
    args.nrec++;

    int ngt = bcf_get_genotypes(args.hdr, rec, &args.gt_arr, &args.ngt_arr);
    if ( ngt<0 ) return dflt;
    if ( ngt!=2*bcf_hdr_nsamples(args.hdr) && ngt!=bcf_hdr_nsamples(args.hdr) ) return dflt;
    ngt /= bcf_hdr_nsamples(args.hdr);

    int i, has_bad = 0, needs_update = 0;
    for (i=0; i<args.ntrios; i++)
    {
        int32_t a,b,c,d,e,f;
        trio_t *trio = &args.trios[i];

        a = args.gt_arr[ngt*trio->imother];
        b = ngt==2 ? args.gt_arr[ngt*trio->imother+1] : bcf_int32_vector_end;
        c = args.gt_arr[ngt*trio->ifather];
        d = ngt==2 ? args.gt_arr[ngt*trio->ifather+1] : bcf_int32_vector_end;
        e = args.gt_arr[ngt*trio->ichild];
        f = ngt==2 ? args.gt_arr[ngt*trio->ichild+1] : bcf_int32_vector_end;

        // missing allele in the child: missing data or daugther's chrY
        if ( bcf_gt_is_missing(e) || bcf_gt_is_missing(f) ) continue;

        // missing data in father
        if ( bcf_gt_is_missing(c) || bcf_gt_is_missing(d) ) continue; 

        uint64_t mother,father,child1,child2;

        // sex chr - father is haploid
        if ( d==bcf_int32_vector_end )
        {
            if ( bcf_gt_is_missing(a) && (bcf_gt_is_missing(b) || b==bcf_int32_vector_end) )
            {
                // either the child does not match the father or the child is diploid
                if ( bcf_gt_allele(c)!=bcf_gt_allele(e) || f!=bcf_int32_vector_end )
                {
                    trio->nbad++;
                    has_bad = 1;
                    if ( args.mode&MODE_DELETE )
                    {   
                        args.gt_arr[ngt*trio->imother] = bcf_gt_missing;
                        if ( b!=bcf_int32_vector_end ) args.gt_arr[ngt*trio->imother+1] = bcf_gt_missing; // should be always true 
                        args.gt_arr[ngt*trio->ifather] = bcf_gt_missing;
                        args.gt_arr[ngt*trio->ichild]  = bcf_gt_missing;
                        if ( ngt==2 ) args.gt_arr[ngt*trio->ichild+1] = bcf_int32_vector_end;
                        needs_update = 1;
                    }
                }
                else
                    trio->nok++;
                continue;
            }
            // son's chrX
            if ( f==bcf_int32_vector_end )
            {
                // chrY - no data in mother
                if ( bcf_gt_allele(e)!=bcf_gt_allele(a) && bcf_gt_allele(e)!=bcf_gt_allele(b) )
                {
                    trio->nbad++;
                    has_bad = 1;
                    if ( args.mode&MODE_DELETE )
                    {
                        args.gt_arr[ngt*trio->imother] = bcf_gt_missing;
                        if ( ngt==2 ) args.gt_arr[ngt*trio->imother+1] = bcf_gt_missing;
                        args.gt_arr[ngt*trio->ifather] = bcf_gt_missing;
                        if ( d!=bcf_int32_vector_end ) args.gt_arr[ngt*trio->ifather+1] = bcf_gt_missing;
                        args.gt_arr[ngt*trio->ichild]  = bcf_gt_missing;
                        needs_update = 1;
                    }
                }
                else
                    trio->nok++;
                continue;
            }

            // skip genotypes which do not fit in a 64bit bitmask
            if ( bcf_gt_allele(a) > 63 || bcf_gt_allele(b) > 63 || bcf_gt_allele(c) > 63 ) continue; 
            if ( bcf_gt_allele(e) > 63 || bcf_gt_allele(f) > 63 ) continue; 

            // daughter's chrX
            mother = (1<<bcf_gt_allele(a)) | (1<<bcf_gt_allele(b));
            father = 1<<bcf_gt_allele(c);
            child1 = 1<<bcf_gt_allele(e);
            child2 = 1<<bcf_gt_allele(f);
        }
        else
        {
            mother = (1<<bcf_gt_allele(a)) | (1<<bcf_gt_allele(b));
            father = (1<<bcf_gt_allele(c)) | (1<<bcf_gt_allele(d));
            child1 = 1<<bcf_gt_allele(e);
            child2 = 1<<bcf_gt_allele(f);
        }

        if ( (mother&child1 && father&child2) || (mother&child2 && father&child1) )
        {
            trio->nok++;
        }
        else
        {
            trio->nbad++;
            has_bad = 1;
            if ( args.mode&MODE_DELETE )
            {
                args.gt_arr[ngt*trio->imother] = bcf_gt_missing;
                if ( b!=bcf_int32_vector_end ) args.gt_arr[ngt*trio->imother+1] = bcf_gt_missing; // should be always true 
                args.gt_arr[ngt*trio->ifather] = bcf_gt_missing;
                if ( d!=bcf_int32_vector_end ) args.gt_arr[ngt*trio->ifather+1] = bcf_gt_missing;
                args.gt_arr[ngt*trio->ichild] = bcf_gt_missing;
                if ( f!=bcf_int32_vector_end ) args.gt_arr[ngt*trio->ichild+1]  = bcf_gt_missing;
                needs_update = 1;
            }
        }
    }

    if ( needs_update && bcf_update_genotypes(args.hdr,rec,args.gt_arr,ngt*bcf_hdr_nsamples(args.hdr)) )
        error("Could not update GT field at %s:%d\n", bcf_seqname(args.hdr,rec),rec->pos+1);

    if ( args.mode&MODE_DELETE ) return rec;
    if ( args.mode&MODE_LIST_GOOD ) return has_bad ? NULL : rec;
    if ( args.mode&MODE_LIST_BAD ) return has_bad ? rec : NULL;

    return NULL;
}

void destroy(void)
{
    int i;
    fprintf(stderr,"# [1]nOK\t[2]nBad\t[3]nSkipped\t[4]Trio\n");
    for (i=0; i<args.ntrios; i++)
    {
        trio_t *trio = &args.trios[i];
        fprintf(stderr,"%d\t%d\t%d\t%s,%s,%s\n", 
            trio->nok,trio->nbad,args.nrec-(trio->nok+trio->nbad),
            bcf_hdr_int2id(args.hdr, BCF_DT_SAMPLE, trio->imother),
            bcf_hdr_int2id(args.hdr, BCF_DT_SAMPLE, trio->ifather),
            bcf_hdr_int2id(args.hdr, BCF_DT_SAMPLE, trio->ichild)
            );
    }
    free(args.gt_arr);
    free(args.trios);
}



