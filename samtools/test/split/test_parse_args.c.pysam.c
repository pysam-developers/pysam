#include "pysam.h"

/*  test/split/test_parse_args.c -- split test cases.

    Copyright (C) 2014 Genome Research Ltd.

    Author: Martin O. Pollard <mp15@sanger.ac.uk>

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

#include "../../bam_split.c"
#include "../test.h"
#include <stdlib.h>
#include <unistd.h>

void setup_test_1(int* argc, char*** argv)
{
    *argc = 1;
    *argv = (char**)calloc(sizeof(char*), 1);
    (*argv)[0] = strdup("prog_name");
}

bool check_test_1(const parsed_opts_t* opts) {
    if ( opts->merged_input_name != NULL
        || opts->unaccounted_header_name != NULL
        || opts->unaccounted_name != NULL
        || strcmp(opts->output_format_string,"%*_%#.%.")
        || opts->verbose == true )
        return false;
    return true;
}

void setup_test_2(int* argc, char*** argv)
{
    *argc = 2;
    *argv = (char**)calloc(sizeof(char*), 2);
    (*argv)[0] = strdup("prog_name");
    (*argv)[1] = strdup("merged.bam");
}

bool check_test_2(const parsed_opts_t* opts) {
    if ( opts->merged_input_name == NULL
        || strcmp(opts->merged_input_name, "merged.bam")
        || opts->unaccounted_header_name != NULL
        || opts->unaccounted_name != NULL
        || strcmp(opts->output_format_string,"%*_%#.%.")
        || opts->verbose == true )
        return false;
    return true;
}

int samtools_test_parse_args_main(int argc, char**argv)
{
    // test state
    const int NUM_TESTS = 2;
    int verbose = 0;
    int success = 0;
    int failure = 0;

    int getopt_char;
    while ((getopt_char = getopt(argc, argv, "v")) != -1) {
        switch (getopt_char) {
            case 'v':
                ++verbose;
                break;
            default:
                fprintf(pysam_stdout, 
                       "usage: test_parse_args [-v]\n\n"
                       " -v verbose output\n"
                       );
                break;
        }
    }

    // Setup pysam_stdout and pysam_stderr redirect
    kstring_t res_pysam_stdout = { 0, 0, NULL };
    kstring_t res_pysam_stderr = { 0, 0, NULL };
    FILE* orig_pysam_stdout = fdopen(dup(STDOUT_FILENO), "a"); // Save pysam_stderr
    FILE* orig_pysam_stderr = fdopen(dup(STDERR_FILENO), "a"); // Save pysam_stderr
    char* tempfname_pysam_stdout = (optind < argc)? argv[optind] : "test_parse_args.tmp.o";
    char* tempfname_pysam_stderr = (optind < argc)? argv[optind] : "test_parse_args.tmp.e";
    FILE* check_pysam_stdout = NULL;
    FILE* check_pysam_stderr = NULL;

    // Cleanup getopt
    optind = 1;

    // setup
    if (verbose) fprintf(orig_pysam_stdout,"BEGIN test 1\n");  // test eliminating a tag that isn't there
    int argc_1;
    char** argv_1;
    setup_test_1(&argc_1, &argv_1);
    if (verbose > 1) {
        fprintf(orig_pysam_stdout, "argc: %d\n", argc_1);
    }
    if (verbose) fprintf(orig_pysam_stdout,"RUN test 1\n");

    // test
    xfreopen(tempfname_pysam_stdout, "w", pysam_stdout); // Redirect pysam_stdout to pipe
    xfreopen(tempfname_pysam_stderr, "w", pysam_stderr); // Redirect pysam_stderr to pipe
    parsed_opts_t* result_1 = parse_args(argc_1, argv_1);
    fclose(pysam_stdout);
    fclose(pysam_stderr);

    if (verbose) fprintf(orig_pysam_stdout, "END RUN test 1\n");
    if (verbose > 1) {
        fprintf(orig_pysam_stdout, "argc: %d\n", argc_1);
    }

    // check result
    res_pysam_stdout.l = res_pysam_stderr.l = 0;
    check_pysam_stdout = fopen(tempfname_pysam_stdout, "r");
    check_pysam_stderr = fopen(tempfname_pysam_stderr, "r");
    if ( !result_1
        && kgetline(&res_pysam_stdout, (kgets_func *)fgets, check_pysam_stdout) >= 0
        && !feof(check_pysam_stdout)
        && res_pysam_stdout.l > 0
        && kgetline(&res_pysam_stderr, (kgets_func *)fgets, check_pysam_stderr) < 0
        && (feof(check_pysam_stderr) || res_pysam_stderr.l == 0)) {
        ++success;
    } else {
        ++failure;
        if (verbose) fprintf(orig_pysam_stdout, "FAIL test 1\n");
    }
    fclose(check_pysam_stderr);
    fclose(check_pysam_stdout);

    // teardown
    cleanup_opts(result_1);
    int i = 0;
    for (i = 0; i < argc_1; ++i) {
        free(argv_1[i]);
    }
    free(argv_1);
    if (verbose) fprintf(orig_pysam_stdout, "END test 1\n");

    // Cleanup getopt
    optind = 1;

    if (verbose) fprintf(orig_pysam_stdout, "BEGIN test 2\n");  // test eliminating a tag that is there
    int argc_2;
    char** argv_2;
    setup_test_2(&argc_2, &argv_2);
    if (verbose > 1) {
        fprintf(orig_pysam_stdout, "argc: %d\n", argc_2);
    }
    if (verbose) fprintf(orig_pysam_stdout, "RUN test 2\n");

    // test
    xfreopen(tempfname_pysam_stdout, "w", pysam_stdout); // Redirect pysam_stdout to pipe
    xfreopen(tempfname_pysam_stderr, "w", pysam_stderr); // Redirect pysam_stderr to pipe
    parsed_opts_t* result_2 = parse_args(argc_2, argv_2);
    fclose(pysam_stdout);
    fclose(pysam_stderr);

    if (verbose) fprintf(orig_pysam_stdout, "END RUN test 2\n");
    if (verbose > 1) {
        fprintf(orig_pysam_stdout, "argc: %d\n", argc_2);
    }

    // check result
    res_pysam_stdout.l = res_pysam_stderr.l = 0;
    check_pysam_stdout = fopen(tempfname_pysam_stdout, "r");
    check_pysam_stderr = fopen(tempfname_pysam_stderr, "r");
    if ( result_2
        && check_test_2(result_2)
        && kgetline(&res_pysam_stdout, (kgets_func *)fgets, check_pysam_stdout) < 0
        && (feof(check_pysam_stdout) || res_pysam_stdout.l == 0)
        && kgetline(&res_pysam_stderr, (kgets_func *)fgets, check_pysam_stderr) < 0
        && (feof(check_pysam_stderr) || res_pysam_stderr.l == 0)) {
        ++success;
    } else {
        ++failure;
        if (verbose) fprintf(orig_pysam_stdout, "FAIL test 2\n");
    }
    fclose(check_pysam_stdout);
    fclose(check_pysam_stderr);

    // teardown
    cleanup_opts(result_2);
    int j = 0;
    for (j = 0; j < argc_2; ++j) {
        free(argv_2[j]);
    }
    free(argv_2);

    if (verbose) fprintf(orig_pysam_stdout, "END test 2\n");


    // Cleanup
    free(res_pysam_stdout.s);
    free(res_pysam_stderr.s);
    remove(tempfname_pysam_stdout);
    remove(tempfname_pysam_stderr);
    fclose(orig_pysam_stdout);
    if (failure > 0)
        fprintf(orig_pysam_stderr, "%d failures %d successes\n", failure, success);
    fclose(orig_pysam_stderr);

    return (success == NUM_TESTS)? EXIT_SUCCESS : EXIT_FAILURE;
}
