#include "samtools.pysam.h"

/*  bam_quickcheck.c -- quickcheck subcommand.

    Copyright (C) 2015 Genome Research Ltd.

    Author: Joshua C. Randall <jcrandall@alum.mit.edu>

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

#include <config.h>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

static void usage_quickcheck(FILE *write_to)
{
    fprintf(write_to,
"Usage: samtools quickcheck [options] <input> [...]\n"
"Options:\n"
"  -v              verbose output (repeat for more verbosity)\n"
"\n"
"Notes:\n"
"\n"
"1. In order to use this command effectively, you should check its exit status;\n"
"   without any -v options it will NOT print any output, even when some files\n"
"   fail the check. One way to use quickcheck might be as a check that all\n"
"   BAM files in a directory are okay:\n"
"\n"
"\tsamtools quickcheck *.bam && echo 'all ok' \\\n"
"\t   || echo 'fail!'\n"
"\n"
"   To also determine which files have failed, use the -v option:\n"
"\n"
"\tsamtools quickcheck -v *.bam > bad_bams.fofn \\\n"
"\t   && echo 'all ok' \\\n"
"\t   || echo 'some files failed check, see bad_bams.fofn'\n"
    );
}

int main_quickcheck(int argc, char** argv)
{
    int verbose = 0;
    hts_verbose = 0;

    const char* optstring = "v";
    int opt;
    while ((opt = getopt(argc, argv, optstring)) != -1) {
        switch (opt) {
        case 'v':
            verbose++;
            break;
        default:
            usage_quickcheck(samtools_stderr);
            return 1;
        }
    }

    argc -= optind;
    argv += optind;

    if (argc < 1) {
        usage_quickcheck(samtools_stdout);
        return 1;
    }

    if (verbose >= 2) {
        fprintf(samtools_stderr, "verbosity set to %d\n", verbose);
    }

    if (verbose >= 4) {
        hts_verbose = 3;
    }

    int ret = 0;
    int i;

    for (i = 0; i < argc; i++) {
        char* fn = argv[i];
        int file_state = 0;

        if (verbose >= 3) fprintf(samtools_stderr, "checking %s\n", fn);

        // attempt to open
        htsFile *hts_fp = hts_open(fn, "r");
        if (hts_fp == NULL) {
            if (verbose >= 2) fprintf(samtools_stderr, "%s could not be opened for reading.\n", fn);
            file_state |= 2;
        }
        else {
            if (verbose >= 3) fprintf(samtools_stderr, "opened %s\n", fn);
            // make sure we have sequence data
            const htsFormat *fmt = hts_get_format(hts_fp);
            if (fmt->category != sequence_data ) {
                if (verbose >= 2) fprintf(samtools_stderr, "%s was not identified as sequence data.\n", fn);
                file_state |= 4;
            }
            else {
                if (verbose >= 3) fprintf(samtools_stderr, "%s is sequence data\n", fn);
                // check header
                bam_hdr_t *header = sam_hdr_read(hts_fp);
                if (header == NULL) {
                    if (verbose >= 2) fprintf(samtools_stderr, "%s caused an error whilst reading its header.\n", fn);
                    file_state |= 8;
                } else {
                    if (header->n_targets <= 0) {
                        if (verbose >= 2) fprintf(samtools_stderr, "%s had no targets in header.\n", fn);
                        file_state |= 8;
                    }
                    else {
                        if (verbose >= 3) fprintf(samtools_stderr, "%s has %d targets in header.\n", fn, header->n_targets);
                    }
                    bam_hdr_destroy(header);
                }
            }
            // check EOF on formats that support this
            int ret;
            if ((ret = hts_check_EOF(hts_fp)) < 0) {
                if (verbose >= 2) fprintf(samtools_stderr, "%s caused an error whilst checking for EOF block.\n", fn);
                file_state |= 16;
            }
            else {
                switch (ret) {
                    case 0:
                        if (verbose >= 2) fprintf(samtools_stderr, "%s was missing EOF block when one should be present.\n", fn);
                        file_state |= 16;
                        break;
                    case 1:
                        if (verbose >= 3) fprintf(samtools_stderr, "%s has good EOF block.\n", fn);
                        break;
                    case 2:
                        if (verbose >= 3) fprintf(samtools_stderr, "%s cannot be checked for EOF block as it is not seekable.\n", fn);
                        break;
                    case 3:
                        if (verbose >= 3) fprintf(samtools_stderr, "%s cannot be checked for EOF block because its filetype does not contain one.\n", fn);
                        break;
                }
            }

            if (hts_close(hts_fp) < 0) {
                file_state |= 32;
                if (verbose >= 2) fprintf(samtools_stderr, "%s did not close cleanly.\n", fn);
            }
        }

        if (file_state > 0 && verbose >= 1) {
            fprintf(samtools_stdout, "%s\n", fn);
        }
        ret |= file_state;
    }

    return ret;
}
