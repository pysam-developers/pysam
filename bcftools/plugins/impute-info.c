/*  plugins/impute-info.c -- adds info metrics to a VCF file.

    Copyright (C) 2015-2016 Genome Research Ltd.

    Author: Shane McCarthy <sm15@sanger.ac.uk>

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


/*

Marchini & Howie, Nature Genetics, 11 (July 2010)

Let G_ij in {0,1,2} be the genotype of the ith individual at the
jth SNP in a study cohort of N samples. Let

    p_ijk = P(G_ij = k | H,G)

be the probability that the genotype at the jth SNP of the ith
individual is k.

Let the expected allele dosage for the genotype at the jth SNP
of the ith individual be
    
    e_ij = p_ij1 + 2 * p_ij2

and define

    f_ij = p_ij1 + 4 * p_ij2

Let theta_j denote the (unknown) population allele frequency of the jth SNP
with:

    theta_j = SUM[i=1..N] e_ij / 2 * N

The IMPUTE2 information measure is then:

    if theta_j in (0,1):
        I(theta_j) = 1 - SUM[i=1..N](f_ij - e_ij^2) / 2 * N * theta_j * (1 - theta_j)
    else:
        I(theta_j) = 1

*/

#include <stdio.h>
#include <stdlib.h>
#include <htslib/vcf.h>
#include <math.h>
#include <getopt.h>

const char *about(void)
{
    return "Add imputation information metrics to the INFO field based on selected FORMAT tags.\n";
}

const char *usage(void)
{
    return 
        "\n"
        "About: Add imputation information metrics to the INFO field based\n"
        "       on selected FORMAT tags. Only the IMPUTE2 INFO metric from\n"
        "       FORMAT/GP tags is currently available.\n"
        "Usage: bcftools +impute-info [General Options] -- [Plugin Options]\n"
        "Options:\n"
        "   run \"bcftools plugin\" for a list of common options\n"
        "\n"
        // "Plugin options:\n"
        // "   -i, --info <list>   information metrics to add [INFO]\n" // [BEAGLE_R2,MACH_R2]
        // "   -t, --tags <tag>    VCF tags to determine the information from [GP]\n"
        // "\n"
        "Example:\n"
        "   bcftools +impute-info in.vcf\n"
        "\n";
}

bcf_hdr_t *in_hdr = NULL, *out_hdr = NULL;
int gp_type = BCF_HT_REAL;
uint8_t *buf = NULL;
int nbuf = 0;   // NB: number of elements, not bytes
int nrec = 0, nskip_gp = 0, nskip_dip = 0;

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    in_hdr  = in;
    out_hdr = out;
    bcf_hdr_append(out_hdr,"##INFO=<ID=INFO,Number=1,Type=Float,Description=\"IMPUTE2 info score\">");
    return 0;
}

bcf1_t *process(bcf1_t *rec)
{
    int nval = 0, i, j, nret = bcf_get_format_values(in_hdr,rec,"GP",(void**)&buf,&nbuf,gp_type);
    if ( nret<0 )
    {
        if (!nskip_gp) fprintf(stderr, "[impute-info.c] Warning: info tag not added to sites without GP tag\n");
        nskip_gp++;
        return rec; // require FORMAT/GP tag, return site unchanged
    }

    nret /= rec->n_sample;
    if ( nret != 3 )
    {
        if (!nskip_dip) fprintf(stderr, "[impute-info.c] Warning: info tag not added to sites that are not biallelic diploid\n");
        nskip_dip++;
        return rec; // require biallelic diploid, return site unchanged
    }

    double esum = 0, e2sum = 0, fsum = 0;
    #define BRANCH(type_t,is_missing,is_vector_end) \
    { \
        type_t *ptr = (type_t*) buf; \
        for (i=0; i<rec->n_sample; i++) \
        { \
            double vals[3] = {0,0,0}; \
            for (j=0; j<nret; j++) \
            { \
                if ( is_missing || is_vector_end ) break; \
                vals[j] = ptr[j]; \
            } \
            double norm = vals[0]+vals[1]+vals[2]; \
            if ( norm ) for (j=0; j<3; j++) vals[j] /= norm; \
            esum  += vals[1] + 2*vals[2]; \
            e2sum += (vals[1] + 2*vals[2]) * (vals[1] + 2*vals[2]); \
            fsum  += vals[1] + 4*vals[2]; \
            ptr   += nret; \
            nval++; \
        } \
    }
    switch (gp_type)
    {
        case BCF_HT_INT:  BRANCH(int32_t,ptr[j]==bcf_int32_missing,ptr[j]==bcf_int32_vector_end); break;
        case BCF_HT_REAL: BRANCH(float,bcf_float_is_missing(ptr[j]),bcf_float_is_vector_end(ptr[j])); break;
    }
    #undef BRANCH

    double theta = esum / (2 * (double)nval);
    float info  = (theta>0 && theta<1) ? (float)(1 - (fsum - e2sum) / (2 * (double)nval * theta * (1.0 - theta))) : 1;

    bcf_update_info_float(out_hdr, rec, "INFO", &info, 1);
    nrec++;
    return rec;
}


void destroy(void)
{
    fprintf(stderr,"Lines total/info-added/unchanged-no-tag/unchanged-not-biallelic-diploid:\t%d/%d/%d/%d\n", nrec+nskip_gp+nskip_dip, nrec, nskip_gp, nskip_dip);
    free(buf);
}
