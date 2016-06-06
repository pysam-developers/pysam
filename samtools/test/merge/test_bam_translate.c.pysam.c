#include "pysam.h"

/*  test/merge/test_bam_translate.c -- header merging test harness.

    Copyright (C) 2013, 2014 Genome Research Ltd.

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

#include "../../bam_sort.c"
#include "../test.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

void dump_read(bam1_t* b) {
    fprintf(pysam_stdout, "->core.tid:(%d)\n", b->core.tid);
    fprintf(pysam_stdout, "->core.pos:(%d)\n", b->core.pos);
    fprintf(pysam_stdout, "->core.bin:(%d)\n", b->core.bin);
    fprintf(pysam_stdout, "->core.qual:(%d)\n", b->core.qual);
    fprintf(pysam_stdout, "->core.l_qname:(%d)\n", b->core.l_qname);
    fprintf(pysam_stdout, "->core.flag:(%d)\n", b->core.flag);
    fprintf(pysam_stdout, "->core.n_cigar:(%d)\n", b->core.n_cigar);
    fprintf(pysam_stdout, "->core.l_qseq:(%d)\n", b->core.l_qseq);
    fprintf(pysam_stdout, "->core.mtid:(%d)\n", b->core.mtid);
    fprintf(pysam_stdout, "->core.mpos:(%d)\n", b->core.mpos);
    fprintf(pysam_stdout, "->core.isize:(%d)\n", b->core.isize);
    if (b->data) {
        fprintf(pysam_stdout, "->data:");
        int i;
        for (i = 0; i < b->l_data; ++i) {
            fprintf(pysam_stdout, "%x ", b->data[i]);
        }
        fprintf(pysam_stdout, "\n");
    }
    if (b->core.l_qname) {
        fprintf(pysam_stdout, "qname: %s\n",bam_get_qname(b));
    }
    if (b->core.l_qseq) {
        fprintf(pysam_stdout, "qseq:");
        int i;
        for (i = 0; i < b->core.l_qseq; ++i) {
            fprintf(pysam_stdout, "%c",seq_nt16_str[seq_nt16_table[bam_seqi(bam_get_seq(b),i)]]);
        }
        fprintf(pysam_stdout, "\n");
        fprintf(pysam_stdout, "qual:");
        for (i = 0; i < b->core.l_qseq; ++i) {
            fprintf(pysam_stdout, "%c",bam_get_qual(b)[i]);
        }
        fprintf(pysam_stdout, "\n");

    }

    if (bam_get_l_aux(b)) {
        int i = 0;
        uint8_t* aux = bam_get_aux(b);

        while (i < bam_get_l_aux(b)) {
            fprintf(pysam_stdout, "%.2s:%c:",aux+i,*(aux+i+2));
            i += 2;
            switch (*(aux+i)) {
                case 'Z':
                    while (*(aux+1+i) != '\0') { putc(*(aux+1+i), pysam_stdout); ++i; }
                    break;
            }
            putc('\n',pysam_stdout);
            ++i;++i;
        }
    }
    fprintf(pysam_stdout, "\n");
}

void trans_tbl_test_init(trans_tbl_t* tbl, int32_t n_targets)
{
    tbl->tid_trans = (int*)calloc(n_targets, sizeof(int32_t));
    tbl->rg_trans = kh_init(c2c);
    tbl->pg_trans = kh_init(c2c);
}

void setup_test_1(bam1_t** b_in, trans_tbl_t* tbl) {
    bam1_t* b;

    b = bam_init1();
    trans_tbl_test_init(tbl, 4);

    tbl->tid_trans[0] = 5;
    tbl->tid_trans[1] = 6;
    tbl->tid_trans[2] = 7;
    tbl->tid_trans[3] = 8;


    b->core.tid = 0;
    b->core.pos = 1334;
    b->core.bin = 0;
    b->core.qual = 10;
    b->core.l_qname = 10;
    b->core.flag = 0;
    b->core.n_cigar = 1;
    b->core.l_qseq = 10;
    b->core.mtid = -1;
    b->core.mpos = 0;
    b->core.isize = -1;
    size_t data_len = 10 + 4 + 5 + 10 + 0;
    b->data = (uint8_t*)malloc(data_len);
    memcpy(b->data,
                     "123456789\0" // q_name
                     "\x00\x00\x00\xA0" // cigar
                     "\x00\x00\x00\x00\x00" // qseq
                     "\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF" // qual
                     "" // aux
           , data_len
                     );
    b->m_data = b->l_data = data_len;

    *b_in = b;
}

void setup_test_2(bam1_t** b_in, trans_tbl_t* tbl) {
    bam1_t* b;

    b = bam_init1();
    trans_tbl_test_init(tbl, 4);

    tbl->tid_trans[0] = 5;
    tbl->tid_trans[1] = 6;
    tbl->tid_trans[2] = 7;
    tbl->tid_trans[3] = 8;
    int in_there = 0;
    khiter_t iter = kh_put(c2c, tbl->rg_trans, strdup("hello"), &in_there);
    kh_value(tbl->rg_trans, iter) = strdup("goodbye");

    b->core.tid = 0;
    b->core.pos = 1334;
    b->core.bin = 0;
    b->core.qual = 10;
    b->core.l_qname = 10;
    b->core.flag = 0;
    b->core.n_cigar = 1;
    b->core.l_qseq = 10;
    b->core.mtid = -1;
    b->core.mpos = 0;
    b->core.isize = -1;
    size_t data_len = 10 + 4 + 5 + 10 + 9;
    b->data = (uint8_t*)malloc(data_len);
    memcpy(b->data,
           "123456789\0" // q_name
           "\x00\x00\x00\xA0" // cigar
           "\x00\x00\x00\x00\x00" // qseq
           "\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF" // qual
           "RGZhello\0" // aux
           , data_len
           );
    b->m_data = b->l_data = data_len;

    *b_in = b;
}

void setup_test_3(bam1_t** b_in, trans_tbl_t* tbl) {
    bam1_t* b;

    b = bam_init1();
    trans_tbl_test_init(tbl, 4);

    tbl->tid_trans[0] = 5;
    tbl->tid_trans[1] = 6;
    tbl->tid_trans[2] = 7;
    tbl->tid_trans[3] = 8;
    int in_there = 0;
    khiter_t iter = kh_put(c2c, tbl->pg_trans, strdup("hello"), &in_there);
    kh_value(tbl->pg_trans,iter) = strdup("goodbye");


    b->core.tid = 0;
    b->core.pos = 1334;
    b->core.bin = 0;
    b->core.qual = 10;
    b->core.l_qname = 10;
    b->core.flag = 0;
    b->core.n_cigar = 1;
    b->core.l_qseq = 10;
    b->core.mtid = -1;
    b->core.mpos = 0;
    b->core.isize = -1;
    size_t data_len = 10 + 4 + 5 + 10 + 9;
    b->data = (uint8_t*)malloc(data_len);
    memcpy(b->data,
           "123456789\0" // q_name
           "\x00\x00\x00\xA0" // cigar
           "\x00\x00\x00\x00\x00" // qseq
           "\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF" // qual
           "PGZhello\0" // aux
           , data_len
           );
    b->m_data = b->l_data = data_len;

    *b_in = b;
}

void setup_test_4(bam1_t** b_in, trans_tbl_t* tbl) {
    bam1_t* b;

    b = bam_init1();
    trans_tbl_test_init(tbl, 4);

    tbl->tid_trans[0] = 5;
    tbl->tid_trans[1] = 6;
    tbl->tid_trans[2] = 7;
    tbl->tid_trans[3] = 8;

    b->core.tid = 0;
    b->core.pos = 1334;
    b->core.bin = 0;
    b->core.qual = 10;
    b->core.l_qname = 10;
    b->core.flag = 0;
    b->core.n_cigar = 1;
    b->core.l_qseq = 10;
    b->core.mtid = -1;
    b->core.mpos = 0;
    b->core.isize = -1;
    size_t data_len = 10 + 4 + 5 + 10 + 12;
    b->data = (uint8_t*)malloc(data_len);
    memcpy(b->data,
           "123456789\0" // q_name
           "\x00\x00\x00\xA0" // cigar
           "\x00\x00\x00\x00\x00" // qseq
           "\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF" // qual
           "RGZrg4hello\0" // aux
           , data_len
           );
    b->m_data = b->l_data = data_len;

    *b_in = b;
}

void setup_test_5(bam1_t** b_in, trans_tbl_t* tbl) {
    bam1_t* b;

    b = bam_init1();
    trans_tbl_test_init(tbl, 4);

    tbl->tid_trans[0] = 5;
    tbl->tid_trans[1] = 6;
    tbl->tid_trans[2] = 7;
    tbl->tid_trans[3] = 8;


    b->core.tid = 0;
    b->core.pos = 1334;
    b->core.bin = 0;
    b->core.qual = 10;
    b->core.l_qname = 10;
    b->core.flag = 0;
    b->core.n_cigar = 1;
    b->core.l_qseq = 10;
    b->core.mtid = -1;
    b->core.mpos = 0;
    b->core.isize = -1;
    size_t data_len = 10 + 4 + 5 + 10 + 12;
    b->data = (uint8_t*)malloc(data_len);
    memcpy(b->data,
           "123456789\0" // q_name
           "\x00\x00\x00\xA0" // cigar
           "\x00\x00\x00\x00\x00" // qseq
           "\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF" // qual
           "PGZpg5hello\0" // aux
           , data_len
           );
    b->m_data = b->l_data = data_len;

    *b_in = b;
}

void setup_test_6(bam1_t** b_in, trans_tbl_t* tbl) {
    bam1_t* b;

    b = bam_init1();
    trans_tbl_test_init(tbl, 4);

    tbl->tid_trans[0] = 5;
    tbl->tid_trans[1] = 6;
    tbl->tid_trans[2] = 7;
    tbl->tid_trans[3] = 8;
    int in_there = 0;
    khiter_t iter_rg = kh_put(c2c, tbl->rg_trans, strdup("hello"), &in_there);
    kh_value(tbl->rg_trans, iter_rg) = strdup("goodbye");
    khiter_t iter_pg = kh_put(c2c, tbl->pg_trans, strdup("quail"), &in_there);
    kh_value(tbl->pg_trans, iter_pg) = strdup("bird");


    b->core.tid = 0;
    b->core.pos = 1334;
    b->core.bin = 0;
    b->core.qual = 10;
    b->core.l_qname = 10;
    b->core.flag = 0;
    b->core.n_cigar = 1;
    b->core.l_qseq = 10;
    b->core.mtid = -1;
    b->core.mpos = 0;
    b->core.isize = -1;
    size_t data_len = 10 + 4 + 5 + 10 + 18;
    b->data = (uint8_t*)malloc(data_len);
    memcpy(b->data,
           "123456789\0" // q_name
           "\x00\x00\x00\xA0" // cigar
           "\x00\x00\x00\x00\x00" // qseq
           "\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF\xFF" // qual
           "RGZhello\0PGZquail\0" // aux
           , data_len
           );
    b->m_data = b->l_data = data_len;

    *b_in = b;
}


int samtools_test_bam_translate_main(int argc, char**argv)
{
    // test state
    const int NUM_TESTS = 6;
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
                break;
        }
    }

    bam1_t* b;

    // Setup pysam_stderr redirect
    kstring_t res = { 0, 0, NULL };
    FILE* orig_pysam_stderr = fdopen(dup(STDERR_FILENO), "a"); // Save pysam_stderr
    char* tempfname = (optind < argc)? argv[optind] : "test_bam_translate.tmp";
    FILE* check = NULL;

    // setup
    if (verbose) fprintf(pysam_stdout, "BEGIN test 1\n");  // TID test
    trans_tbl_t tbl1;
    setup_test_1(&b,&tbl1);
    if (verbose > 1) {
        fprintf(pysam_stdout, "b\n");
        dump_read(b);
    }
    if (verbose) fprintf(pysam_stdout, "RUN test 1\n");

    // test
    xfreopen(tempfname, "w", pysam_stderr); // Redirect pysam_stderr to pipe
    bam_translate(b, &tbl1);
    fclose(pysam_stderr);

    if (verbose) fprintf(pysam_stdout, "END RUN test 1\n");
    if (verbose > 1) {
        fprintf(pysam_stdout, "b\n");
        dump_read(b);
    }

    // check result
    check = fopen(tempfname, "r");
    res.l = 0;
    if (kgetline(&res, (kgets_func *)fgets, check) < 0 &&
        (feof(check) || res.l == 0) ) {
        ++success;
    } else {
        ++failure;
        if (verbose) fprintf(pysam_stdout, "FAIL test 1\n");
    }
    fclose(check);

    // teardown
    bam_destroy1(b);
    trans_tbl_destroy(&tbl1);
    if (verbose) fprintf(pysam_stdout, "END test 1\n");

    // setup
    if (verbose) fprintf(pysam_stdout, "BEGIN test 2\n");  // RG exists and translate test
    trans_tbl_t tbl2;
    setup_test_2(&b,&tbl2);
    if (verbose > 1) {
        fprintf(pysam_stdout, "b\n");
        dump_read(b);
    }
    if (verbose) fprintf(pysam_stdout, "RUN test 2\n");

    // test
    xfreopen(tempfname, "w", pysam_stderr); // Redirect pysam_stderr to pipe
    bam_translate(b, &tbl2);
    fclose(pysam_stderr);

    if (verbose) fprintf(pysam_stdout, "END RUN test 2\n");
    if (verbose > 1) {
        fprintf(pysam_stdout, "b\n");
        dump_read(b);
    }

    // check result
    check = fopen(tempfname, "r");
    res.l = 0;
    if (kgetline(&res, (kgets_func *)fgets, check) < 0 &&
        (feof(check) || res.l == 0) ) {
        ++success;
    } else {
        ++failure;
        if (verbose) fprintf(pysam_stdout, "FAIL test 2\n");
    }
    fclose(check);

    // teardown
    bam_destroy1(b);
    trans_tbl_destroy(&tbl2);
    if (verbose) fprintf(pysam_stdout, "END test 2\n");

    if (verbose) fprintf(pysam_stdout, "BEGIN test 3\n");  // PG exists and translate  test
    // setup
    trans_tbl_t tbl3;
    setup_test_3(&b,&tbl3);
    if (verbose > 1) {
        fprintf(pysam_stdout, "b\n");
        dump_read(b);
    }
    if (verbose) fprintf(pysam_stdout, "RUN test 3\n");

    // test
    xfreopen(tempfname, "w", pysam_stderr); // Redirect pysam_stderr to pipe
    bam_translate(b, &tbl3);
    fclose(pysam_stderr);

    if (verbose) fprintf(pysam_stdout, "END RUN test 3\n");
    if (verbose > 1) {
        fprintf(pysam_stdout, "b\n");
        dump_read(b);
    }

    // check result
    check = fopen(tempfname, "r");
    res.l = 0;
    if (kgetline(&res, (kgets_func *)fgets, check) < 0 &&
        (feof(check) || res.l == 0)) {
        ++success;
    } else {
        ++failure;
        if (verbose) fprintf(pysam_stdout, "FAIL test 3\n");
    }
    fclose(check);

    // teardown
    bam_destroy1(b);
    trans_tbl_destroy(&tbl3);
    if (verbose) fprintf(pysam_stdout, "END test 3\n");

    if (verbose) fprintf(pysam_stdout, "BEGIN test 4\n");  // RG test non-existent
    // setup
    trans_tbl_t tbl4;
    setup_test_4(&b,&tbl4);
    if (verbose > 1) {
        fprintf(pysam_stdout, "b\n");
        dump_read(b);
    }
    if (verbose) fprintf(pysam_stdout, "RUN test 4\n");

    // test
    xfreopen(tempfname, "w", pysam_stderr); // Redirect pysam_stderr to pipe
    bam_translate(b, &tbl4);
    fclose(pysam_stderr);

    if (verbose) fprintf(pysam_stdout, "END RUN test 4\n");
    if (verbose > 1) {
        fprintf(pysam_stdout, "b\n");
        dump_read(b);
    }
    // check result
    check = fopen(tempfname, "r");
    res.l = 0;
    if (kgetline(&res, (kgets_func *)fgets, check) >= 0 &&
        strcmp("[bam_translate] RG tag \"rg4hello\" on read \"123456789\" encountered with no corresponding entry in header, tag lost. Unknown tags are only reported once per input file for each tag ID.",res.s) == 0) {
        ++success;
    } else {
        ++failure;
        if (verbose) fprintf(pysam_stdout, "FAIL test 4\n");
    }
    fclose(check);

    // teardown
    bam_destroy1(b);
    trans_tbl_destroy(&tbl4);
    if (verbose) fprintf(pysam_stdout, "END test 4\n");

    if (verbose) fprintf(pysam_stdout, "BEGIN test 5\n");  // PG test non-existent
    // setup
    trans_tbl_t tbl5;
    setup_test_5(&b,&tbl5);
    if (verbose > 1) {
        fprintf(pysam_stdout, "b\n");
        dump_read(b);
        fprintf(pysam_stdout, "RUN test 5\n");
    }
    // test
    xfreopen(tempfname, "w", pysam_stderr); // Redirect pysam_stderr to pipe
    bam_translate(b, &tbl5);
    fclose(pysam_stderr);

    if (verbose) fprintf(pysam_stdout, "END RUN test 5\n");
    if (verbose > 1) {
        fprintf(pysam_stdout, "b\n");
        dump_read(b);
    }

    // check result
    check = fopen(tempfname, "r");
    res.l = 0;
    if (kgetline(&res, (kgets_func *)fgets, check) >= 0 &&
        strcmp("[bam_translate] PG tag \"pg5hello\" on read \"123456789\" encountered with no corresponding entry in header, tag lost. Unknown tags are only reported once per input file for each tag ID.",res.s) == 0) {
        ++success;
    } else {
        ++failure;
        if (verbose) fprintf(pysam_stdout, "FAIL test 5\n");
    }
    fclose(check);

    // teardown
    bam_destroy1(b);
    trans_tbl_destroy(&tbl5);
    if (verbose) fprintf(pysam_stdout, "END test 5\n");

    if (verbose) fprintf(pysam_stdout, "BEGIN test 6\n");  // RG and PG exists and translate test
    // setup
    trans_tbl_t tbl6;
    setup_test_6(&b,&tbl6);
    if (verbose > 1) {
        fprintf(pysam_stdout, "b\n");
        dump_read(b);
    }
    if (verbose) fprintf(pysam_stdout, "RUN test 6\n");

    // test
    xfreopen(tempfname, "w", pysam_stderr); // Redirect pysam_stderr to pipe
    bam_translate(b, &tbl6);
    fclose(pysam_stderr);

    if (verbose) fprintf(pysam_stdout, "END RUN test 6\n");
    if (verbose > 1) {
        fprintf(pysam_stdout, "b\n");
        dump_read(b);
    }

    // check result
    check = fopen(tempfname, "r");
    res.l = 0;
    if (kgetline(&res, (kgets_func *)fgets, check) < 0 &&
        (feof(check) || res.l == 0) ) {
        ++success;
    } else {
        ++failure;
        if (verbose) fprintf(pysam_stdout, "FAIL test 6\n");
    }
    fclose(check);

    // teardown
    bam_destroy1(b);
    trans_tbl_destroy(&tbl6);
    if (verbose) fprintf(pysam_stdout, "END test 6\n");

    // Cleanup
    free(res.s);
    remove(tempfname);
    if (failure > 0)
        fprintf(orig_pysam_stderr, "%d failures %d successes\n", failure, success);
    fclose(orig_pysam_stderr);

    return (success == NUM_TESTS)? EXIT_SUCCESS : EXIT_FAILURE;
}
