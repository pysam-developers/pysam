/*  vcfview.c -- VCF/BCF conversion, view, subset and filter VCF/BCF files.

    Copyright (C) 2013-2014 Genome Research Ltd.

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
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.  */

#include <stdio.h>
#include <strings.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include "bcftools.h"
#include "filter.h"
#include "htslib/khash_str2int.h"

#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

#define ALLELE_NONREF 1
#define ALLELE_MINOR 2
#define ALLELE_ALT1 3
#define ALLELE_MAJOR 4
#define ALLELE_NONMAJOR 5

#define GT_NEED_HOM 1
#define GT_NEED_HET 2
#define GT_NO_HOM   3
#define GT_NO_HET   4
#define GT_NEED_MISSING 5
#define GT_NO_MISSING 6

typedef struct _args_t
{
    filter_t *filter;
    char *filter_str;
    int filter_logic;   // one of FLT_INCLUDE/FLT_EXCLUDE (-i or -e)

    bcf_srs_t *files;
    bcf_hdr_t *hdr, *hnull, *hsub; // original header, sites-only header, subset header
    char **argv, *format, *sample_names, *subset_fname, *targets_list, *regions_list;
    int argc, clevel, n_threads, output_type, print_header, update_info, header_only, n_samples, *imap, calc_ac;
    int trim_alts, sites_only, known, novel, min_alleles, max_alleles, private_vars, uncalled, phased;
    int min_ac, min_ac_type, max_ac, max_ac_type, min_af_type, max_af_type, gt_type;
    int *ac, mac;
    float min_af, max_af;
    char *fn_ref, *fn_out, **samples;
    int sample_is_file, force_samples;
    char *include_types, *exclude_types;
    int include, exclude;
    int record_cmd_line;
    htsFile *out;
}
args_t;

static void init_data(args_t *args)
{
    int i;
    args->hdr = args->files->readers[0].header;

    if (args->calc_ac && args->update_info)
    {
        bcf_hdr_append(args->hdr,"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">");
        bcf_hdr_append(args->hdr,"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
    }
    if (args->record_cmd_line) bcf_hdr_append_version(args->hdr, args->argc, args->argv, "bcftools_view");
    else bcf_hdr_sync(args->hdr);

    // setup sample data
    if (args->sample_names)
    {
        void *hdr_samples = khash_str2int_init();
        for (i=0; i<bcf_hdr_nsamples(args->hdr); i++)
            khash_str2int_inc(hdr_samples, bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,i));

        void *exclude = (args->sample_names[0]=='^') ? khash_str2int_init() : NULL;
        int nsmpl;
        char **smpl = NULL;
        args->samples = NULL; args->n_samples = 0;
        smpl = hts_readlist(exclude ? &args->sample_names[1] : args->sample_names, args->sample_is_file, &nsmpl);
        if ( !smpl )
        {
            error("Could not read the list: \"%s\"\n", exclude ? &args->sample_names[1] : args->sample_names);
        }

        if ( exclude )
        {
            for (i=0; i<nsmpl; i++) {
                if (!khash_str2int_has_key(hdr_samples,smpl[i])) {
                    if (args->force_samples) {
                        fprintf(stderr, "Warn: exclude called for sample that does not exist in header: \"%s\"... skipping\n", smpl[i]);
                    } else {
                        error("Error: exclude called for sample that does not exist in header: \"%s\". Use \"--force-samples\" to ignore this error.\n", smpl[i]);
                    }
                }
                khash_str2int_inc(exclude, smpl[i]);
            }

            for (i=0; i<bcf_hdr_nsamples(args->hdr); i++)
            {
                if ( exclude && khash_str2int_has_key(exclude,bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,i))  ) continue;
                args->samples = (char**) realloc(args->samples, (args->n_samples+1)*sizeof(const char*));
                args->samples[args->n_samples++] = strdup(bcf_hdr_int2id(args->hdr,BCF_DT_SAMPLE,i));
            }
            khash_str2int_destroy(exclude);
        }
        else
        {
            for (i=0; i<nsmpl; i++) {
                if (!khash_str2int_has_key(hdr_samples,smpl[i])) {
                    if (args->force_samples) {
                        fprintf(stderr, "Warn: subset called for sample that does not exist in header: \"%s\"... skipping\n", smpl[i]);
                        continue;
                    } else {
                        error("Error: subset called for sample that does not exist in header: \"%s\". Use \"--force-samples\" to ignore this error.\n", smpl[i]);
                    }
                }
                args->samples = (char**) realloc(args->samples, (args->n_samples+1)*sizeof(const char*));
                args->samples[args->n_samples++] = strdup(smpl[i]);
            }
        }
        for (i=0; i<nsmpl; i++) free(smpl[i]);
        free(smpl);
        khash_str2int_destroy(hdr_samples);
        if (args->n_samples == 0) {
            fprintf(stderr, "Warn: subsetting has removed all samples\n");
            args->sites_only = 1;
        }
    }

    if (args->n_samples)
        args->imap = (int*)malloc(args->n_samples * sizeof(int));

    // determine variant types to include/exclude
    if (args->include_types || args->exclude_types) {
        if (args->include_types && args->exclude_types) {
            fprintf(stderr, "Error: only supply one of --include-types, --exclude-types options\n");
            exit(1);
        }
        char **type_list = 0;
        int m = 0, n = 0;
        const char *q, *p;
        for (q = p = args->include_types ? args->include_types : args->exclude_types;; ++p) {
            if (*p == ',' || *p == 0) {
                if (m == n) {
                    m = m? m<<1 : 16;
                    type_list = (char**)realloc(type_list, m * sizeof(char*));
                }
                type_list[n] = (char*)calloc(p - q + 1, 1);
                strncpy(type_list[n++], q, p - q);
                q = p + 1;
                if (*p == 0) break;
            }
        }
        type_list = (char**)realloc(type_list, n * sizeof(char*));

        if (args->include_types) {
            args->include = 0;
            for (i = 0; i < n; ++i) {
                if (strcmp(type_list[i], "snps") == 0) args->include |= VCF_SNP<<1;
                else if (strcmp(type_list[i], "indels") == 0) args->include |= VCF_INDEL<<1;
                else if (strcmp(type_list[i], "mnps") == 0) args->include |= VCF_MNP<<1;
                else if (strcmp(type_list[i], "other") == 0) args->include |= VCF_OTHER<<1;
                else if (strcmp(type_list[i], "ref") == 0) args->include |= VCF_OTHER<<1;
                else if (strcmp(type_list[i], "bnd") == 0) args->include |= VCF_BND<<1;
                else {
                    fprintf(stderr, "[E::%s] unknown type\n", type_list[i]);
                    fprintf(stderr, "Accepted types are snps, indels, mnps, other\n");
                    exit(1);
                }
            }
        }
        if (args->exclude_types) {
            args->exclude = 0;
            for (i = 0; i < n; ++i) {
                if (strcmp(type_list[i], "snps") == 0) args->exclude |= VCF_SNP<<1;
                else if (strcmp(type_list[i], "indels") == 0) args->exclude |= VCF_INDEL<<1;
                else if (strcmp(type_list[i], "mnps") == 0) args->exclude |= VCF_MNP<<1;
                else if (strcmp(type_list[i], "other") == 0) args->exclude |= VCF_OTHER<<1;
                else if (strcmp(type_list[i], "ref") == 0) args->exclude |= VCF_OTHER<<1;
                else if (strcmp(type_list[i], "bnd") == 0) args->exclude |= VCF_BND<<1;
                else {
                    fprintf(stderr, "[E::%s] unknown type\n", type_list[i]);
                    fprintf(stderr, "Accepted types are snps, indels, mnps, other\n");
                    exit(1);
                }
            }
        }
        for (i = 0; i < n; ++i)
            free(type_list[i]);
        free(type_list);
    }

    // setup output
    char modew[8];
    strcpy(modew, "w");
    if (args->clevel >= 0 && args->clevel <= 9) sprintf(modew + 1, "%d", args->clevel);
    if (args->output_type==FT_BCF) strcat(modew, "bu");         // uncompressed BCF
    else if (args->output_type & FT_BCF) strcat(modew, "b");    // compressed BCF
    else if (args->output_type & FT_GZ) strcat(modew,"z");      // compressed VCF
    args->out = hts_open(args->fn_out ? args->fn_out : "-", modew);
    if ( !args->out ) error("%s: %s\n", args->fn_out,strerror(errno));
    if ( args->n_threads > 0)
        hts_set_opt(args->out, HTS_OPT_THREAD_POOL, args->files->p);

    // headers: hdr=full header, hsub=subset header, hnull=sites only header
    if (args->sites_only){
        args->hnull = bcf_hdr_subset(args->hdr, 0, 0, 0);
        bcf_hdr_remove(args->hnull, BCF_HL_FMT, NULL);
    }
    if (args->n_samples > 0)
    {
        args->hsub = bcf_hdr_subset(args->hdr, args->n_samples, args->samples, args->imap);
        if ( !args->hsub ) error("Error occurred while subsetting samples\n");
        if ( args->n_samples != bcf_hdr_nsamples(args->hsub) )
        {
            int i;
            for (i=0; i<args->n_samples; i++)
                if ( args->imap[i]<0 ) error("Error: No such sample: \"%s\"\n", args->samples[i]);
        }
    }

    if ( args->filter_str )
        args->filter = filter_init(args->hdr, args->filter_str);
}

static void destroy_data(args_t *args)
{
    int i;
    if ( args->imap ) {
        for (i = 0; i < args->n_samples; ++i)
            free(args->samples[i]);
        free(args->samples);
        free(args->imap);
    }
    if (args->hnull) bcf_hdr_destroy(args->hnull);
    if (args->hsub) bcf_hdr_destroy(args->hsub);
    if ( args->filter )
        filter_destroy(args->filter);
    free(args->ac);
}

// true if all samples are phased.
// haploid genotypes are considered phased
// ./. => not phased, .|. => phased
int bcf_all_phased(const bcf_hdr_t *header, bcf1_t *line)
{
    bcf_unpack(line, BCF_UN_FMT);
    bcf_fmt_t *fmt_ptr = bcf_get_fmt(header, line, "GT");
    int all_phased = 1;
    if ( fmt_ptr )
    {
        int i, isample;
        for (isample=0; isample<line->n_sample; isample++)
        {
            int sample_phased = 0;
            #define BRANCH_INT(type_t,vector_end) { \
                type_t *p = (type_t*) (fmt_ptr->p + isample*fmt_ptr->size); \
                for (i=0; i<fmt_ptr->n; i++) \
                { \
                    if (fmt_ptr->n == 1 || (p[i] == vector_end && i == 1)) { sample_phased = 1; break; } /* haploid phased by definition */ \
                    if ( p[i] == vector_end ) { break; }; /* smaller ploidy */ \
                    if ( bcf_gt_is_missing(p[i]) ) continue; /* missing allele */ \
                    if ((p[i])&1) { \
                        sample_phased = 1; \
                        break; \
                    } \
                } \
            }
            switch (fmt_ptr->type) {
                case BCF_BT_INT8:  BRANCH_INT(int8_t,  bcf_int8_vector_end); break;
                case BCF_BT_INT16: BRANCH_INT(int16_t, bcf_int16_vector_end); break;
                case BCF_BT_INT32: BRANCH_INT(int32_t, bcf_int32_vector_end); break;
                default: fprintf(stderr, "[E::%s] todo: fmt_type %d\n", __func__, fmt_ptr->type); exit(1); break;
            }
            #undef BRANCH_INT
            if (!sample_phased) {
                all_phased = 0;
                break;
            }
        }
    }
    return all_phased;
}

int subset_vcf(args_t *args, bcf1_t *line)
{
    if ( args->min_alleles && line->n_allele < args->min_alleles ) return 0; // min alleles
    if ( args->max_alleles && line->n_allele > args->max_alleles ) return 0; // max alleles
    if (args->novel || args->known)
    {
        if ( args->novel && (line->d.id[0]!='.' || line->d.id[1]!=0) ) return 0; // skip sites which are known, ID != '.'
        if ( args->known && line->d.id[0]=='.' && line->d.id[1]==0 ) return 0;  // skip sites which are novel, ID == '.'
    }

    if (args->include || args->exclude)
    {
        int line_type = bcf_get_variant_types(line);
        if ( args->include && !((line_type<<1) & args->include) ) return 0; // include only given variant types
        if ( args->exclude &&   (line_type<<1) & args->exclude  ) return 0; // exclude given variant types
    }

    if ( args->filter )
    {
        int ret = filter_test(args->filter, line, NULL);
        if ( args->filter_logic==FLT_INCLUDE ) { if ( !ret ) return 0; }
        else if ( ret ) return 0;
    }

    hts_expand(int, line->n_allele, args->mac, args->ac);
    int i, an = 0, non_ref_ac = 0;
    if (args->calc_ac) {
        bcf_calc_ac(args->hdr, line, args->ac, BCF_UN_INFO|BCF_UN_FMT); // get original AC and AN values from INFO field if available, otherwise calculate
        for (i=1; i<line->n_allele; i++)
            non_ref_ac += args->ac[i];
        for (i=0; i<line->n_allele; i++)
            an += args->ac[i];
    }

    if (args->n_samples)
    {
        int non_ref_ac_sub = 0, *ac_sub = (int*) calloc(line->n_allele,sizeof(int));
        bcf_subset(args->hdr, line, args->n_samples, args->imap);
        if (args->calc_ac) {
            bcf_calc_ac(args->hsub, line, ac_sub, BCF_UN_FMT); // recalculate AC and AN
            an = 0;
            for (i=0; i<line->n_allele; i++) {
                args->ac[i] = ac_sub[i];
                an += ac_sub[i];
            }
            for (i=1; i<line->n_allele; i++)
                non_ref_ac_sub += ac_sub[i];
            if (args->private_vars) {
                if (args->private_vars == FLT_INCLUDE && !(non_ref_ac_sub > 0 && non_ref_ac == non_ref_ac_sub)) { free(ac_sub); return 0; } // select private sites
                if (args->private_vars == FLT_EXCLUDE && non_ref_ac_sub > 0 && non_ref_ac == non_ref_ac_sub) { free(ac_sub); return 0; } // exclude private sites
            }
            non_ref_ac = non_ref_ac_sub;
        }
        free(ac_sub);
    }

    bcf_fmt_t *gt_fmt;
    if ( args->gt_type && (gt_fmt=bcf_get_fmt(args->hdr,line,"GT")) )
    {
        int nhet = 0, nhom = 0, nmiss = 0;
        for (i=0; i<bcf_hdr_nsamples(args->hdr); i++)
        {
            int type = bcf_gt_type(gt_fmt,i,NULL,NULL);
            if ( type==GT_HET_RA || type==GT_HET_AA )
            {
                if ( args->gt_type==GT_NO_HET ) return 0;
                nhet = 1;
            }
            else if ( type==GT_UNKN )
            {
                if ( args->gt_type==GT_NO_MISSING ) return 0;
                nmiss = 1;
            }
            else
            {
                if ( args->gt_type==GT_NO_HOM ) return 0;
                nhom = 1;
            }
        }
        if ( args->gt_type==GT_NEED_HOM && !nhom ) return 0;
        else if ( args->gt_type==GT_NEED_HET && !nhet ) return 0;
        else if ( args->gt_type==GT_NEED_MISSING && !nmiss ) return 0;
    }

    int minor_ac = 0;
    int major_ac = 0;
    if ( args->calc_ac )
    {
        minor_ac = args->ac[0];
        major_ac = args->ac[0];
        for (i=1; i<line->n_allele; i++){
            if (args->ac[i] < minor_ac) { minor_ac = args->ac[i]; }
            if (args->ac[i] > major_ac) { major_ac = args->ac[i]; }
        }
    }

    if (args->min_ac!=-1)
    {
        if (args->min_ac_type == ALLELE_NONREF && args->min_ac>non_ref_ac) return 0; // min AC
        else if (args->min_ac_type == ALLELE_MINOR && args->min_ac>minor_ac) return 0; // min minor AC
        else if (args->min_ac_type == ALLELE_ALT1 && args->min_ac>args->ac[1]) return 0; // min 1st alternate AC
        else if (args->min_ac_type == ALLELE_MAJOR && args->min_ac > major_ac) return 0; // min major AC
        else if (args->min_ac_type == ALLELE_NONMAJOR && args->min_ac > an-major_ac) return 0; // min non-major AC
    }
    if (args->max_ac!=-1)
    {
        if (args->max_ac_type == ALLELE_NONREF && args->max_ac<non_ref_ac) return 0; // max AC
        else if (args->max_ac_type == ALLELE_MINOR && args->max_ac<minor_ac) return 0; // max minor AC
        else if (args->max_ac_type == ALLELE_ALT1 && args->max_ac<args->ac[1]) return 0; // max 1st alternate AC
        else if (args->max_ac_type == ALLELE_MAJOR && args->max_ac < major_ac) return 0; // max major AC
        else if (args->max_ac_type == ALLELE_NONMAJOR && args->max_ac < an-major_ac) return 0; // max non-major AC
    }
    if (args->min_af!=-1)
    {
        if (an == 0) return 0; // freq not defined, skip site
        if (args->min_af_type == ALLELE_NONREF && args->min_af>non_ref_ac/(double)an) return 0; // min AF
        else if (args->min_af_type == ALLELE_MINOR && args->min_af>minor_ac/(double)an) return 0; // min minor AF
        else if (args->min_af_type == ALLELE_ALT1 && args->min_af>args->ac[1]/(double)an) return 0; // min 1st alternate AF
        else if (args->min_af_type == ALLELE_MAJOR && args->min_af > major_ac/(double)an) return 0; // min major AF
        else if (args->min_af_type == ALLELE_NONMAJOR && args->min_af > (an-major_ac)/(double)an) return 0; // min non-major AF
    }
    if (args->max_af!=-1)
    {
        if (an == 0) return 0; // freq not defined, skip site
        if (args->max_af_type == ALLELE_NONREF && args->max_af<non_ref_ac/(double)an) return 0; // max AF
        else if (args->max_af_type == ALLELE_MINOR && args->max_af<minor_ac/(double)an) return 0; // max minor AF
        else if (args->max_af_type == ALLELE_ALT1 && args->max_af<args->ac[1]/(double)an) return 0; // max 1st alternate AF
        else if (args->max_af_type == ALLELE_MAJOR && args->max_af < major_ac/(double)an) return 0; // max major AF
        else if (args->max_af_type == ALLELE_NONMAJOR && args->max_af < (an-major_ac)/(double)an) return 0; // max non-major AF
    }
    if (args->uncalled) {
        if (args->uncalled == FLT_INCLUDE && an > 0) return 0; // select uncalled
        if (args->uncalled == FLT_EXCLUDE && an == 0) return 0; // skip if uncalled
    }
    if (args->calc_ac && args->update_info) {
        bcf_update_info_int32(args->hdr, line, "AC", &args->ac[1], line->n_allele-1);
        bcf_update_info_int32(args->hdr, line, "AN", &an, 1);
    }
    if (args->trim_alts)
    {
        int ret = bcf_trim_alleles(args->hsub ? args->hsub : args->hdr, line);
        if ( ret<0 ) error("Error: Could not trim alleles at %s:%d\n", bcf_seqname(args->hsub ? args->hsub : args->hdr, line), line->pos+1);
    }
    if (args->phased) {
        int phased = bcf_all_phased(args->hdr, line);
        if (args->phased == FLT_INCLUDE && !phased) { return 0; } // skip unphased
        if (args->phased == FLT_EXCLUDE && phased) { return 0; } // skip phased
    }
    if (args->sites_only) bcf_subset(args->hsub ? args->hsub : args->hdr, line, 0, 0);
    return 1;
}

void set_allele_type (int *atype, char *atype_string)
{
    *atype = ALLELE_NONREF;
    if (strcmp(atype_string, "minor") == 0) {
        *atype = ALLELE_MINOR;
    }
    else if (strcmp(atype_string, "alt1") == 0) {
        *atype = ALLELE_ALT1;
    }
    else if (strcmp(atype_string, "nref") == 0) {
        *atype = ALLELE_NONREF;
    }
    else if (strcmp(atype_string, "major") == 0) {
        *atype = ALLELE_MAJOR;
    }
    else if (strcmp(atype_string, "nonmajor") == 0) {
        *atype = ALLELE_NONMAJOR;
    }
    else {
        error("Error: allele type not recognised. Expected one of nref|alt1|minor|major|nonmajor, got \"%s\".\n", atype_string);
    }
}

static void usage(args_t *args)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   VCF/BCF conversion, view, subset and filter VCF/BCF files.\n");
    fprintf(stderr, "Usage:   bcftools view [options] <in.vcf.gz> [region1 [...]]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output options:\n");
    fprintf(stderr, "    -G,   --drop-genotypes              drop individual genotype information (after subsetting if -s option set)\n");
    fprintf(stderr, "    -h/H, --header-only/--no-header     print the header only/suppress the header in VCF output\n");
    fprintf(stderr, "    -l,   --compression-level [0-9]     compression level: 0 uncompressed, 1 best speed, 9 best compression [%d]\n", args->clevel);
    fprintf(stderr, "          --no-version                  do not append version and command line to the header\n");
    fprintf(stderr, "    -o,   --output-file <file>          output file name [stdout]\n");
    fprintf(stderr, "    -O,   --output-type <b|u|z|v>       b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "    -r, --regions <region>              restrict to comma-separated list of regions\n");
    fprintf(stderr, "    -R, --regions-file <file>           restrict to regions listed in a file\n");
    fprintf(stderr, "    -t, --targets [^]<region>           similar to -r but streams rather than index-jumps. Exclude regions with \"^\" prefix\n");
    fprintf(stderr, "    -T, --targets-file [^]<file>        similar to -R but streams rather than index-jumps. Exclude regions with \"^\" prefix\n");
    fprintf(stderr, "        --threads <int>                 number of extra (de)compression threads [0]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Subset options:\n");
    fprintf(stderr, "    -a, --trim-alt-alleles        trim alternate alleles not seen in the subset\n");
    fprintf(stderr, "    -I, --no-update               do not (re)calculate INFO fields for the subset (currently INFO/AC and INFO/AN)\n");
    fprintf(stderr, "    -s, --samples [^]<list>       comma separated list of samples to include (or exclude with \"^\" prefix)\n");
    fprintf(stderr, "    -S, --samples-file [^]<file>  file of samples to include (or exclude with \"^\" prefix)\n");
    fprintf(stderr, "        --force-samples           only warn about unknown subset samples\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Filter options:\n");
    fprintf(stderr, "    -c/C, --min-ac/--max-ac <int>[:<type>]      minimum/maximum count for non-reference (nref), 1st alternate (alt1), least frequent\n");
    fprintf(stderr, "                                                   (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]\n");
    fprintf(stderr, "    -f,   --apply-filters <list>                require at least one of the listed FILTER strings (e.g. \"PASS,.\")\n");
    fprintf(stderr, "    -g,   --genotype [^]<hom|het|miss>          require one or more hom/het/missing genotype or, if prefixed with \"^\", exclude sites with hom/het/missing genotypes\n");
    fprintf(stderr, "    -i/e, --include/--exclude <expr>            select/exclude sites for which the expression is true (see man page for details)\n");
    fprintf(stderr, "    -k/n, --known/--novel                       select known/novel sites only (ID is not/is '.')\n");
    fprintf(stderr, "    -m/M, --min-alleles/--max-alleles <int>     minimum/maximum number of alleles listed in REF and ALT (e.g. -m2 -M2 for biallelic sites)\n");
    fprintf(stderr, "    -p/P, --phased/--exclude-phased             select/exclude sites where all samples are phased\n");
    fprintf(stderr, "    -q/Q, --min-af/--max-af <float>[:<type>]    minimum/maximum frequency for non-reference (nref), 1st alternate (alt1), least frequent\n");
    fprintf(stderr, "                                                   (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]\n");
    fprintf(stderr, "    -u/U, --uncalled/--exclude-uncalled         select/exclude sites without a called genotype\n");
    fprintf(stderr, "    -v/V, --types/--exclude-types <list>        select/exclude comma-separated list of variant types: snps,indels,mnps,ref,bnd,other [null]\n");
    fprintf(stderr, "    -x/X, --private/--exclude-private           select/exclude sites where the non-reference alleles are exclusive (private) to the subset samples\n");
    fprintf(stderr, "\n");
    exit(1);
}

int main_vcfview(int argc, char *argv[])
{
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    args->argc    = argc; args->argv = argv;
    args->files   = bcf_sr_init();
    args->clevel  = -1;
    args->print_header = 1;
    args->update_info = 1;
    args->output_type = FT_VCF;
    args->n_threads = 0;
    args->record_cmd_line = 1;
    args->min_ac = args->max_ac = args->min_af = args->max_af = -1;
    int targets_is_file = 0, regions_is_file = 0;

    static struct option loptions[] =
    {
        {"genotype",required_argument,NULL,'g'},
        {"compression-level",required_argument,NULL,'l'},
        {"threads",required_argument,NULL,9},
        {"header-only",no_argument,NULL,'h'},
        {"no-header",no_argument,NULL,'H'},
        {"exclude",required_argument,NULL,'e'},
        {"include",required_argument,NULL,'i'},
        {"trim-alt-alleles",no_argument,NULL,'a'},
        {"no-update",no_argument,NULL,'I'},
        {"drop-genotypes",no_argument,NULL,'G'},
        {"private",no_argument,NULL,'x'},
        {"exclude-private",no_argument,NULL,'X'},
        {"uncalled",no_argument,NULL,'u'},
        {"exclude-uncalled",no_argument,NULL,'U'},
        {"apply-filters",required_argument,NULL,'f'},
        {"known",no_argument,NULL,'k'},
        {"novel",no_argument,NULL,'n'},
        {"min-alleles",required_argument,NULL,'m'},
        {"max-alleles",required_argument,NULL,'M'},
        {"samples",required_argument,NULL,'s'},
        {"samples-file",required_argument,NULL,'S'},
        {"force-samples",no_argument,NULL,1},
        {"output-type",required_argument,NULL,'O'},
        {"output-file",required_argument,NULL,'o'},
        {"types",required_argument,NULL,'v'},
        {"exclude-types",required_argument,NULL,'V'},
        {"targets",required_argument,NULL,'t'},
        {"targets-file",required_argument,NULL,'T'},
        {"regions",required_argument,NULL,'r'},
        {"regions-file",required_argument,NULL,'R'},
        {"min-ac",required_argument,NULL,'c'},
        {"max-ac",required_argument,NULL,'C'},
        {"min-af",required_argument,NULL,'q'},
        {"max-af",required_argument,NULL,'Q'},
        {"phased",no_argument,NULL,'p'},
        {"exclude-phased",no_argument,NULL,'P'},
        {"no-version",no_argument,NULL,8},
        {NULL,0,NULL,0}
    };
    char *tmp;
    while ((c = getopt_long(argc, argv, "l:t:T:r:R:o:O:s:S:Gf:knv:V:m:M:auUhHc:C:Ii:e:xXpPq:Q:g:",loptions,NULL)) >= 0)
    {
        char allele_type[8] = "nref";
        switch (c)
        {
            case 'O':
                switch (optarg[0]) {
                    case 'b': args->output_type = FT_BCF_GZ; break;
                    case 'u': args->output_type = FT_BCF; break;
                    case 'z': args->output_type = FT_VCF_GZ; break;
                    case 'v': args->output_type = FT_VCF; break;
                    default: error("The output type \"%s\" not recognised\n", optarg);
                };
                break;
            case 'l':
                args->clevel = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --compression-level %s\n", optarg);
                args->output_type |= FT_GZ; 
                break;
            case 'o': args->fn_out = optarg; break;
            case 'H': args->print_header = 0; break;
            case 'h': args->header_only = 1; break;

            case 't': args->targets_list = optarg; break;
            case 'T': args->targets_list = optarg; targets_is_file = 1; break;
            case 'r': args->regions_list = optarg; break;
            case 'R': args->regions_list = optarg; regions_is_file = 1; break;

            case 's': args->sample_names = optarg; break;
            case 'S': args->sample_names = optarg; args->sample_is_file = 1; break;
            case  1 : args->force_samples = 1; break;
            case 'a': args->trim_alts = 1; args->calc_ac = 1; break;
            case 'I': args->update_info = 0; break;
            case 'G': args->sites_only = 1; break;

            case 'f': args->files->apply_filters = optarg; break;
            case 'k': args->known = 1; break;
            case 'n': args->novel = 1; break;
            case 'm':
                args->min_alleles = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --min-alleles %s\n", optarg);
                break;
            case 'M': 
                args->max_alleles = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --max-alleles %s\n", optarg);
                break;
            case 'v': args->include_types = optarg; break;
            case 'V': args->exclude_types = optarg; break;
            case 'e': args->filter_str = optarg; args->filter_logic |= FLT_EXCLUDE; break;
            case 'i': args->filter_str = optarg; args->filter_logic |= FLT_INCLUDE; break;

            case 'c':
            {
                args->min_ac_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%d:%s",&args->min_ac, allele_type)!=2 && sscanf(optarg,"%d",&args->min_ac)!=1 )
                    error("Error: Could not parse --min-ac %s\n", optarg);
                set_allele_type(&args->min_ac_type, allele_type);
                args->calc_ac = 1;
                break;
            }
            case 'C':
            {
                args->max_ac_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%d:%s",&args->max_ac, allele_type)!=2 && sscanf(optarg,"%d",&args->max_ac)!=1 )
                    error("Error: Could not parse --max-ac %s\n", optarg);
                set_allele_type(&args->max_ac_type, allele_type);
                args->calc_ac = 1;
                break;
            }
            case 'q':
            {
                args->min_af_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%f:%s",&args->min_af, allele_type)!=2 && sscanf(optarg,"%f",&args->min_af)!=1 )
                    error("Error: Could not parse --min_af %s\n", optarg);
                set_allele_type(&args->min_af_type, allele_type);
                args->calc_ac = 1;
                break;
            }
            case 'Q':
            {
                args->max_af_type = ALLELE_NONREF;
                if ( sscanf(optarg,"%f:%s",&args->max_af, allele_type)!=2 && sscanf(optarg,"%f",&args->max_af)!=1 )
                    error("Error: Could not parse --min_af %s\n", optarg);
                set_allele_type(&args->max_af_type, allele_type);
                args->calc_ac = 1;
                break;
            }

            case 'x': args->private_vars |= FLT_INCLUDE; args->calc_ac = 1; break;
            case 'X': args->private_vars |= FLT_EXCLUDE; args->calc_ac = 1; break;
            case 'u': args->uncalled |= FLT_INCLUDE; args->calc_ac = 1; break;
            case 'U': args->uncalled |= FLT_EXCLUDE; args->calc_ac = 1; break;
            case 'p': args->phased |= FLT_INCLUDE; break; // phased
            case 'P': args->phased |= FLT_EXCLUDE; break; // exclude-phased
            case 'g':
            {
                if ( !strcasecmp(optarg,"hom") ) args->gt_type = GT_NEED_HOM;
                else if ( !strcasecmp(optarg,"het") ) args->gt_type = GT_NEED_HET;
                else if ( !strcasecmp(optarg,"miss") ) args->gt_type = GT_NEED_MISSING;
                else if ( !strcasecmp(optarg,"^hom") ) args->gt_type = GT_NO_HOM;
                else if ( !strcasecmp(optarg,"^het") ) args->gt_type = GT_NO_HET;
                else if ( !strcasecmp(optarg,"^miss") ) args->gt_type = GT_NO_MISSING;
                else error("The argument to -g not recognised. Expected one of hom/het/miss/^hom/^het/^miss, got \"%s\".\n", optarg);
                break;
            }
            case  9 : args->n_threads = strtol(optarg, 0, 0); break;
            case  8 : args->record_cmd_line = 0; break;
            case '?': usage(args);
            default: error("Unknown argument: %s\n", optarg);
        }
    }

    if ( args->filter_logic == (FLT_EXCLUDE|FLT_INCLUDE) ) error("Only one of -i or -e can be given.\n");
    if ( args->private_vars > FLT_EXCLUDE ) error("Only one of -x or -X can be given.\n");
    if ( args->uncalled > FLT_EXCLUDE ) error("Only one of -u or -U can be given.\n");
    if ( args->phased > FLT_EXCLUDE ) error("Only one of -p or -P can be given.\n");

    if ( args->sample_names && args->update_info) args->calc_ac = 1;

    char *fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin
        else usage(args);
    }
    else fname = argv[optind];

    // read in the regions from the command line
    if ( args->regions_list )
    {
        if ( bcf_sr_set_regions(args->files, args->regions_list, regions_is_file)<0 )
            error("Failed to read the regions: %s\n", args->regions_list);
    }
    else if ( optind+1 < argc )
    {
        int i;
        kstring_t tmp = {0,0,0};
        kputs(argv[optind+1],&tmp);
        for (i=optind+2; i<argc; i++) { kputc(',',&tmp); kputs(argv[i],&tmp); }
        if ( bcf_sr_set_regions(args->files, tmp.s, 0)<0 )
            error("Failed to read the regions: %s\n", tmp.s);
        free(tmp.s);
    }
    if ( args->targets_list )
    {
        if ( bcf_sr_set_targets(args->files, args->targets_list, targets_is_file, 0)<0 )
            error("Failed to read the targets: %s\n", args->targets_list);
    }

    if ( bcf_sr_set_threads(args->files, args->n_threads)<0 ) error("Failed to create threads\n");
    if ( !bcf_sr_add_reader(args->files, fname) ) error("Failed to open %s: %s\n", fname,bcf_sr_strerror(args->files->errnum));

    init_data(args);
    bcf_hdr_t *out_hdr = args->hnull ? args->hnull : (args->hsub ? args->hsub : args->hdr);
    if (args->print_header)
        bcf_hdr_write(args->out, out_hdr);
    else if ( args->output_type & FT_BCF )
        error("BCF output requires header, cannot proceed with -H\n");

    int ret = 0;
    if (!args->header_only)
    {
        while ( bcf_sr_next_line(args->files) )
        {
            bcf1_t *line = args->files->readers[0].buffer[0];
            if ( line->errcode && out_hdr!=args->hdr ) error("Undefined tags in the header, cannot proceed in the sample subset mode.\n");
            if ( subset_vcf(args, line) )
                bcf_write1(args->out, out_hdr, line);
        }
        ret = args->files->errnum;
        if ( ret ) fprintf(stderr,"Error: %s\n", bcf_sr_strerror(args->files->errnum));
    }
    hts_close(args->out);
    destroy_data(args);
    bcf_sr_destroy(args->files);
    free(args);
    return ret;
}
