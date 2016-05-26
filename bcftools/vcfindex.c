
/*  vcfindex.c -- Index bgzip compressed VCF/BCF files for random access.

    Copyright (C) 2014-2016 Genome Research Ltd.

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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include <sys/stat.h>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>
#include "bcftools.h"

#define BCF_LIDX_SHIFT    14

static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Index bgzip compressed VCF/BCF files for random access.\n");
    fprintf(stderr, "Usage:   bcftools index [options] <in.bcf>|<in.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Indexing options:\n");
    fprintf(stderr, "    -c, --csi            generate CSI-format index for VCF/BCF files [default]\n");
    fprintf(stderr, "    -f, --force          overwrite index if it already exists\n");
    fprintf(stderr, "    -m, --min-shift INT  set minimal interval size for CSI indices to 2^INT [14]\n");
    fprintf(stderr, "    -t, --tbi            generate TBI-format index for VCF files\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Stats options:\n");
    fprintf(stderr, "    -n, --nrecords       print number of records based on existing index file\n");
    fprintf(stderr, "    -s, --stats   print per contig stats based on existing index file\n");
    fprintf(stderr, "\n");
    exit(1);
}

int vcf_index_stats(char *fname, int stats)
{
    char *fn_out = NULL;
    FILE *out;
    out = fn_out ? fopen(fn_out, "w") : stdout;

    const char **seq;
    int i, nseq;
    tbx_t *tbx = NULL;
    hts_idx_t *idx = NULL;

    htsFile *fp = hts_open(fname,"r");
    if ( !fp ) { fprintf(stderr,"Could not read %s\n", fname); return 1; }
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if ( !hdr ) { fprintf(stderr,"Could not read the header: %s\n", fname); return 1; }

    if ( hts_get_format(fp)->format==vcf )
    {
        tbx = tbx_index_load(fname);
        if ( !tbx ) { fprintf(stderr,"Could not load TBI index: %s\n", fname); return 1; }
    }
    else if ( hts_get_format(fp)->format==bcf )
    {
        idx = bcf_index_load(fname);
        if ( !idx ) { fprintf(stderr,"Could not load CSI index: %s\n", fname); return 1; }
    }
    else
    {
        fprintf(stderr,"Could not detect the file type as VCF or BCF: %s\n", fname);
        return 1;
    }

    seq = tbx ? tbx_seqnames(tbx, &nseq) : bcf_index_seqnames(idx, hdr, &nseq);
    uint64_t sum = 0;
    for (i=0; i<nseq; i++)
    {
        uint64_t records, v;
        hts_idx_get_stat(tbx ? tbx->idx : idx, i, &records, &v);
        sum+=records;
        if (stats&2 || !records) continue;
        bcf_hrec_t *hrec = bcf_hdr_get_hrec(hdr, BCF_HL_CTG, "ID", seq[i], NULL);
        int hkey = hrec ? bcf_hrec_find_key(hrec, "length") : -1;
        fprintf(out,"%s\t%s\t%" PRIu64 "\n", seq[i], hkey<0?".":hrec->vals[hkey], records);
    }
    if (!sum)
    {
        // No counts found.
        // Is this because index version has no stored count data, or no records?
        bcf1_t *rec = bcf_init1();
        if (bcf_read1(fp, hdr, rec) >= 0)
        {
            fprintf(stderr,"%s index of %s does not contain any count metadata. Please re-index with a newer version of bcftools or tabix.\n", tbx ? "TBI" : "CSI", fname);
            return 1;
        }
        bcf_destroy1(rec);
    }
    if (stats&2) fprintf(out, "%" PRIu64 "\n", sum);
    free(seq);
    fclose(out);
    hts_close(fp);
    bcf_hdr_destroy(hdr);
    if (tbx)
        tbx_destroy(tbx);
    if (idx)
        hts_idx_destroy(idx);
    return 0;
}

int main_vcfindex(int argc, char *argv[])
{
    int c, force = 0, tbi = 0, stats = 0;
    int min_shift = BCF_LIDX_SHIFT;

    static struct option loptions[] =
    {
        {"csi",no_argument,NULL,'c'},
        {"tbi",no_argument,NULL,'t'},
        {"force",no_argument,NULL,'f'},
        {"min-shift",required_argument,NULL,'m'},
        {"stats",no_argument,NULL,'s'},
        {"nrecords",no_argument,NULL,'n'},
        {NULL, 0, NULL, 0}
    };

    char *tmp;
    while ((c = getopt_long(argc, argv, "ctfm:sn", loptions, NULL)) >= 0)
    {
        switch (c)
        {
            case 'c': tbi = 0; break;
            case 't': tbi = 1; min_shift = 0; break;
            case 'f': force = 1; break;
            case 'm': 
                min_shift = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --min-shift %s\n", optarg);
                break;
            case 's': stats |= 1; break;
            case 'n': stats |= 2; break;
            default: usage();
        }
    }
    if ( optind==argc ) usage();
    if (stats>2)
    {
        fprintf(stderr, "[E::%s] expected only one of --stats or --nrecords options\n", __func__);
        return 1;
    }
    if (tbi && min_shift>0)
    {
        fprintf(stderr, "[E::%s] min-shift option only expected for CSI indices \n", __func__);
        return 1;
    }
    if (min_shift < 0 || min_shift > 30)
    {
        fprintf(stderr, "[E::%s] expected min_shift in range [0,30] (%d)\n", __func__, min_shift);
        return 1;
    }

    char *fname = argv[optind];
    if (stats) return vcf_index_stats(fname, stats);

    htsFile *fp = hts_open(fname,"r"); 
    if ( !fp ) error("Failed to read %s\n", fname);
    htsFormat type = *hts_get_format(fp);
    hts_close(fp);

    if ( (type.format!=bcf && type.format!=vcf) || type.compression!=bgzf )
    {
        fprintf(stderr, "[E::%s] unknown filetype; expected bgzip compressed VCF or BCF\n", __func__);
        if ( type.compression!=bgzf )
            fprintf(stderr, "[E::%s] was the VCF/BCF compressed with bgzip?\n", __func__);
        return 1;
    }
    if (tbi && type.format==bcf)
    {
        fprintf(stderr, "[Warning] TBI-index does not work for BCF files. Generating CSI instead.\n");
        tbi = 0; min_shift = BCF_LIDX_SHIFT;
    }
    if (min_shift == 0 && type.format==bcf)
    {
        fprintf(stderr, "[E::%s] Require min_shift>0 for BCF files.\n", __func__);
        return 1;
    }
    if (!tbi && type.format==vcf && min_shift == 0)
    {
        fprintf(stderr, "[Warning] min-shift set to 0 for VCF file. Generating TBI file.\n");
        tbi = 1;
    }

    if (!force)
    {
        // Before complaining about existing index, check if the VCF file isn't newer.
        char *idx_fname = (char*)alloca(strlen(fname) + 5);
        strcat(strcpy(idx_fname, fname), tbi ? ".tbi" : ".csi");
        struct stat stat_tbi, stat_file;
        if ( stat(idx_fname, &stat_tbi)==0 )
        {
            stat(fname, &stat_file);
            if ( stat_file.st_mtime <= stat_tbi.st_mtime )
            {
                fprintf(stderr,"[E::%s] the index file exists. Please use '-f' to overwrite.\n", __func__);
                return 1;
            }
        }
    }

    if (type.format==bcf)
    {
        if ( bcf_index_build(fname, min_shift) != 0 )
        {
            fprintf(stderr,"[E::%s] bcf_index_build failed for %s\n", __func__, fname);
            return 1;
        }
    }
    else
    {
        if ( tbx_index_build(fname, min_shift, &tbx_conf_vcf) != 0 )
        {
            fprintf(stderr,"[E::%s] tbx_index_build failed for %s\n", __func__, fname);
            return 1;
        }
    }
    return 0;
}
