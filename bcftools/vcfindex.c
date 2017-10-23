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
#include <htslib/kstring.h>
#include <htslib/bgzf.h>
#include "bcftools.h"

#define BCF_LIDX_SHIFT    14

static void usage(void)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   Index bgzip compressed VCF/BCF files for random access.\n");
    fprintf(stderr, "Usage:   bcftools index [options] <in.bcf>|<in.vcf.gz>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Indexing options:\n");
    fprintf(stderr, "    -c, --csi                generate CSI-format index for VCF/BCF files [default]\n");
    fprintf(stderr, "    -f, --force              overwrite index if it already exists\n");
    fprintf(stderr, "    -m, --min-shift INT      set minimal interval size for CSI indices to 2^INT [14]\n");
    fprintf(stderr, "    -o, --output-file FILE   optional output index file name\n");
    fprintf(stderr, "    -t, --tbi                generate TBI-format index for VCF files\n");
    fprintf(stderr, "        --threads            sets the number of threads [0]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Stats options:\n");
    fprintf(stderr, "    -n, --nrecords       print number of records based on existing index file\n");
    fprintf(stderr, "    -s, --stats          print per contig stats based on existing index file\n");
    fprintf(stderr, "\n");
    exit(1);
}

int vcf_index_stats(char *fname, int stats)
{
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
        if ( !tbx ) { fprintf(stderr,"Could not load index for VCF: %s\n", fname); return 1; }
    }
    else if ( hts_get_format(fp)->format==bcf )
    {
        idx = bcf_index_load(fname);
        if ( !idx ) { fprintf(stderr,"Could not load index for BCF file: %s\n", fname); return 1; }
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
        printf("%s\t%s\t%" PRIu64 "\n", seq[i], hkey<0?".":hrec->vals[hkey], records);
    }
    if (!sum)
    {
        // No counts found.
        // Is this because index version has no stored count data, or no records?
        bcf1_t *rec = bcf_init1();
        if (bcf_read1(fp, hdr, rec) >= 0)
        {
            fprintf(stderr,"index of %s does not contain any count metadata. Please re-index with a newer version of bcftools or tabix.\n", fname);
            return 1;
        }
        bcf_destroy1(rec);
    }
    if (stats&2) printf("%" PRIu64 "\n", sum);
    free(seq);
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
    int c, force = 0, tbi = 0, stats = 0, n_threads = 0;
    int min_shift = BCF_LIDX_SHIFT;
    char *outfn = NULL;

    static struct option loptions[] =
    {
        {"csi",no_argument,NULL,'c'},
        {"tbi",no_argument,NULL,'t'},
        {"force",no_argument,NULL,'f'},
        {"min-shift",required_argument,NULL,'m'},
        {"stats",no_argument,NULL,'s'},
        {"nrecords",no_argument,NULL,'n'},
        {"threads",required_argument,NULL,9},
        {"output-file",required_argument,NULL,'o'},
        {NULL, 0, NULL, 0}
    };

    char *tmp;
    while ((c = getopt_long(argc, argv, "ctfm:sno:", loptions, NULL)) >= 0)
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
            case 9:
                n_threads = strtol(optarg,&tmp,10);
                if ( *tmp ) error("Could not parse argument: --threads %s\n", optarg);
                break;
            case 'o': outfn = optarg; break;
            default: usage();
        }
    }
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

    char *fname = NULL;
    if ( optind>=argc )
    {
        if ( !isatty(fileno((FILE *)stdin)) ) fname = "-";  // reading from stdin
        else usage();
    }
    else fname = argv[optind];
    if (stats) return vcf_index_stats(fname, stats);

    kstring_t idx_fname = {0,0,0};
    if (outfn)
        kputs(outfn,&idx_fname);
    else
    {
        if (!strcmp(fname, "-")) { fprintf(stderr, "[E::%s] must specify an output path for index file when reading VCF/BCF from stdin\n", __func__); return 1; }
        ksprintf(&idx_fname, "%s.%s", fname, tbi ? "tbi" : "csi");
    }
    if (!force)
    {
        // Before complaining about existing index, check if the VCF file isn't newer.
        struct stat stat_tbi, stat_file;
        if ( stat(idx_fname.s, &stat_tbi)==0 )
        {
            stat(fname, &stat_file);
            if ( stat_file.st_mtime <= stat_tbi.st_mtime )
            {
                fprintf(stderr,"[E::%s] the index file exists. Please use '-f' to overwrite %s\n", __func__, idx_fname.s);
                free(idx_fname.s);
                return 1;
            }
        }

        // check for truncated files, allow only with -f
        BGZF *fp = bgzf_open(fname, "r");
        if ( !fp ) error("index: failed to open %s\n", fname);
        if ( bgzf_check_EOF(fp)!=1 ) error("index: the input is probably truncated, use -f to index anyway: %s\n", fname);
        if ( bgzf_close(fp)!=0 ) error("index: close failed: %s\n", fname);
    }

    int ret = bcf_index_build3(fname, idx_fname.s, min_shift, n_threads);
    free(idx_fname.s);
    if (ret != 0) {
        if (ret == -2)
            error("index: failed to open \"%s\"\n", fname);
        else if (ret == -3)
            error("index: \"%s\" is in a format that cannot be usefully indexed\n", fname);
        else
            error("index: failed to create index for \"%s\"\n", fname);
    }
    return 0;
}
