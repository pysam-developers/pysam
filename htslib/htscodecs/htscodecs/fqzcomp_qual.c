/*
 * Copyright (c) 2011-2013, 2018-2022 Genome Research Ltd.
 * Author(s): James Bonfield
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
 *       Institute nor the names of its contributors may be used to endorse
 *       or promote products derived from this software without specific
 *       prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
 * LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// We use generic maps to turn 0-M into 0-N where N <= M
// before adding these into the context.  These are used
// for positions, running-diffs and quality values.
//
// This can be used as a simple divisor, eg pos/24 to get
// 2 bits of positional data for each quarter along a 100bp
// read, or it can be tailored for specific such as noting
// the first 5 cycles are poor, then we have stability and
// a gradual drop off in the last 20 or so.  Perhaps we then
// map pos 0-4=0, 5-79=1, 80-89=2, 90-99=3.
//
// We don't need to specify how many bits of data we are
// using (2 in the above example), as that is just implicit
// in the values in the map.  Specify not to use a map simply
// disables that context type (our map is essentially 0-M -> 0).

// Example of command line usage:
//
// f=~/scratch/data/q4
// cc -Wall -DTEST_MAIN -O3 -g fqzcomp_qual2.c -lm
// ./a.out $f > /tmp/_ && ./a.out -d < /tmp/_ > /tmp/__ && cmp /tmp/__ $f

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include <inttypes.h>
#include <sys/types.h>

#include "fqzcomp_qual.h"
#include "varint.h"
#include "utils.h"

#define CTX_BITS 16
#define CTX_SIZE (1<<CTX_BITS)

#ifndef MIN
#  define MIN(a,b) ((a)<(b)?(a):(b))
#  define MAX(a,b) ((a)>(b)?(a):(b))
#endif

#define QMAX 256
#define QBITS 12
#define QSIZE (1<<QBITS)

#define NSYM 2
#include "c_simple_model.h"

#undef NSYM
#define NSYM QMAX
//#include "c_escape_model.h"
#include "c_simple_model.h"
//#include "c_cdf_model.h"
//#include "c_cdf16_model.h"

// An array of 0,0,0, 1,1,1,1, 3, 5,5
// is turned into a run-length of 3x0, 4x1, 0x2, 1x4, 0x4, 2x5,
// which then becomes 3 4 0 1 0 2.
//
// NB: size > 255 therefore means we need to repeatedly read to find
// the actual run length.
// Alternatively we could bit-encode instead of byte encode, eg BETA.
static int store_array(unsigned char *out, unsigned int *array, int size) {
    unsigned char tmp[2048];

    int i, j, k;
    for (i = j = k = 0; i < size; j++) {
        int run_len = i;
        while (i < size && array[i] == j)
            i++;
        run_len = i-run_len;

        int r;
        do {
            r = MIN(255, run_len);
            tmp[k++] = r;
            run_len -= r;
        } while (r == 255);
    }
    while (i < size)
        tmp[k++] = 0, j++;

    // RLE on out.
    //    1 2 3 3 3 3 3 4 4    5
    // => 1 2 3 3 +3... 4 4 +0 5
    int last = -1;
    for (i = j = 0; j < k; i++) {
        out[i] = tmp[j++];
        if (out[i] == last) {
            int n = j;
            while (j < k && tmp[j] == last)
                j++;
            out[++i] = j-n;
        } else {
            last = out[i];
        }
    }
    k = i;

//    fprintf(stderr, "Store_array %d => %d {", size, k);
//    for (i = 0; i < k; i++)
//      fprintf(stderr, "%d,", out[i]);
//    fprintf(stderr, "}\n");
    return k;
}

static int read_array(unsigned char *in, size_t in_size, unsigned int *array, int size) {
    unsigned char R[1024];
    int i, j, z, last = -1, nb = 0;

    size = MIN(1024, size);

    // Remove level one of run-len encoding
    for (i = j = z = 0; z < size && i < in_size; i++) {
        int run = in[i];
        R[j++] = run;
        z += run;
        if (run == last) {
            if (i+1 >= in_size)
                return -1;
            int copy = in[++i];
            z += run * copy;
            while (copy-- && z <= size && j < 1024)
                R[j++] = run;
        }
        if (j >= 1024)
            return -1;
        last = run;
    }
    nb = i;

    // Now expand inner level of run-length encoding
    int R_max = j;
    for (i = j = z = 0; j < size; i++) {
        int run_len = 0;
        int run_part;
        if (z >= R_max)
            return -1;
        do {
            run_part = R[z++];
            run_len += run_part;
        } while (run_part == 255 && z < R_max);
        if (run_part == 255)
            return -1;

        while (run_len && j < size)
            run_len--, array[j++] = i;
    }

    return nb;
}

// FIXME: how to auto-tune these rather than trial and error?
// r2 = READ2
// qa = qual avg (0, 2, 4)
static int strat_opts[][12] = {
//   qb  qs pb ps db ds ql sl pl  dl  r2 qa
    {10, 5, 4,-1, 2, 1, 0, 14, 10, 14, 0,-1}, // basic options (level < 7)
    {8,  5, 7, 0, 0, 0, 0, 14, 8,  14, 1,-1}, // e.g. HiSeq 2000
    {12, 6, 2, 0, 2, 3, 0, 9,  12, 14, 0, 0}, // e.g. MiSeq
    {12, 6, 0, 0, 0, 0, 0, 12, 0,  0,  0, 0}, // e.g. IonTorrent; adaptive O1
    {0,  0, 0, 0, 0, 0, 0, 0,  0,  0,  0, 0}, // custom
};
static int nstrats = sizeof(strat_opts) / sizeof(*strat_opts);

#ifdef HAVE_BUILTIN_PREFETCH
static inline void mm_prefetch(void *x) {
    __builtin_prefetch(x);
}
#else
static inline void mm_prefetch(void *x) {
    // Fetch and discard is quite close to a genuine prefetch
    *(volatile int *)x;
}
#endif

typedef struct {
    unsigned int qctx;  // quality sub-context
    unsigned int p;     // pos (bytes remaining)
    unsigned int delta; // delta running total
    unsigned int prevq; // previous quality
    unsigned int s;     // selector
    unsigned int qtot, qlen;
    unsigned int first_len;
    unsigned int last_len;
    ssize_t rec;
    unsigned int ctx;
} fqz_state;

static void dump_table(unsigned int *tab, int size, char *name) {
    int i, last = -99, run = 0;
    fprintf(stderr, "\t%s\t{", name);
    for (i = 0; i < size; i++) {
        if (tab[i] == last) {
            run++;
        } else if (run == 1 && tab[i] == last+1) {
            int first = last;
            do {
                last = tab[i];
                i++;
            } while (i < size && tab[i] == last+1);
            i--;

            // Want 0,1,2,3,3,3 as 0..2 3x3, not 0..3 3x2
            if (tab[i] == tab[i+1])
                i--;
            if (tab[i] != first)
                fprintf(stderr, "..%d", tab[i]);
            run = 1;
            last = -99;
        } else {
            if (run > 1)
                fprintf(stderr, " x %d%s%d", run, i?", ":"", tab[i]);
            else
                fprintf(stderr, "%s%d", i?", ":"", tab[i]);
            run = 1;
            last = tab[i];
        }
    }
    if (run > 1)
        fprintf(stderr, " x %d", run);
    fprintf(stderr, "}\n");
}

static void dump_map(unsigned int *map, int size, char *name) {
    int i, c = 0;
    fprintf(stderr, "\t%s\t{", name);
    for (i = 0; i < size; i++)
        if (map[i] != INT_MAX)
            fprintf(stderr, "%s%d=%d", c++?", ":"", i, map[i]);
    fprintf(stderr, "}\n");
}

#pragma GCC diagnostic ignored "-Wunused-function"
static void dump_params(fqz_gparams *gp) {
    fprintf(stderr, "Global params = {\n");
    fprintf(stderr, "\tvers\t%d\n", gp->vers);
    fprintf(stderr, "\tgflags\t0x%02x\n", gp->gflags);
    fprintf(stderr, "\tnparam\t%d\n", gp->nparam);
    fprintf(stderr, "\tmax_sel\t%d\n", gp->max_sel);
    fprintf(stderr, "\tmax_sym\t%d\n", gp->max_sym);
    if (gp->gflags & GFLAG_HAVE_STAB)
        dump_table(gp->stab, 256, "stab");
    fprintf(stderr, "}\n");

    int i;
    for (i = 0; i < gp->nparam; i++) {
        fqz_param *pm = &gp->p[i];
        fprintf(stderr, "\nParam[%d] = {\n", i);
        fprintf(stderr, "\tcontext\t0x%04x\n", pm->context);
        fprintf(stderr, "\tpflags\t0x%02x\n",  pm->pflags);
        fprintf(stderr, "\tmax_sym\t%d\n",  pm->max_sym);
        fprintf(stderr, "\tqbits\t%d\n",   pm->qbits);
        fprintf(stderr, "\tqshift\t%d\n",  pm->qshift);
        fprintf(stderr, "\tqloc\t%d\n",    pm->qloc);
        fprintf(stderr, "\tsloc\t%d\n",    pm->sloc);
        fprintf(stderr, "\tploc\t%d\n",    pm->ploc);
        fprintf(stderr, "\tdloc\t%d\n",    pm->dloc);

        if (pm->pflags & PFLAG_HAVE_QMAP)
            dump_map(pm->qmap, 256, "qmap");

        if (pm->pflags & PFLAG_HAVE_QTAB)
            dump_table(pm->qtab, 256, "qtab");
        if (pm->pflags & PFLAG_HAVE_PTAB)
            dump_table(pm->ptab, 1024, "ptab");
        if (pm->pflags & PFLAG_HAVE_DTAB)
            dump_table(pm->dtab, 256, "dtab");
        fprintf(stderr, "}\n");
    }
}

typedef struct {
    SIMPLE_MODEL(QMAX,_) *qual;
    SIMPLE_MODEL(256,_)   len[4];
    SIMPLE_MODEL(2,_)     revcomp;
    SIMPLE_MODEL(256,_)   sel;
    SIMPLE_MODEL(2,_)     dup;
} fqz_model;

static int fqz_create_models(fqz_model *m, fqz_gparams *gp) {
    int i;

    if (!(m->qual = htscodecs_tls_alloc(sizeof(*m->qual) * CTX_SIZE)))
        return -1;

    for (i = 0; i < CTX_SIZE; i++)
        SIMPLE_MODEL(QMAX,_init)(&m->qual[i], gp->max_sym+1);

    for (i = 0; i < 4; i++)
        SIMPLE_MODEL(256,_init)(&m->len[i],256);

    SIMPLE_MODEL(2,_init)(&m->revcomp,2);
    SIMPLE_MODEL(2,_init)(&m->dup,2);
    if (gp->max_sel > 0)
        SIMPLE_MODEL(256,_init)(&m->sel, gp->max_sel+1);

    return 0;
}

static void fqz_destroy_models(fqz_model *m) {
    htscodecs_tls_free(m->qual);
}

static inline unsigned int fqz_update_ctx(fqz_param *pm, fqz_state *state, int q) {
    unsigned int last = 0; // pm->context
    state->qctx = (state->qctx << pm->qshift) + pm->qtab[q];
    last += (state->qctx & pm->qmask) << pm->qloc;

    // The final shifts have been factored into the tables already.
    last += pm->ptab[MIN(1023, state->p)];      // << pm->ploc
    last += pm->dtab[MIN(255,  state->delta)];  // << pm->dloc
    last += state->s << pm->sloc;

    // On the fly average is slow work.
    // However it can be slightly better than using a selector bit
    // as it's something we can compute on the fly and thus doesn't
    // consume output bits for storing the selector itself.
    //
    // Q4 (novaseq.bam)
    // qtot+=q*q -DQ1=8.84 -DQ2=8.51 -DQ3=7.70; 7203598 (-0.7%)
    // qtot+=q   -DQ1=2.96 -DQ2=2.85 -DQ3=2.69; 7207315
    // vs old delta;                            7255614 (default params)
    // vs 2 bit selector (no delta)             7203006 (-x 0x8261000e80)
    // vs 2 bit selector (no delta)             7199153 (-x 0x7270000e70) -0.8%
    // vs 2 bit selector (no delta)             7219668 (-x 0xa243000ea0)
    //{
    //  double qa = state->qtot / (state->qlen+.01);
    //  //fprintf(stderr, "%f\n", qa);
    //  int x = 0;
    //  if (qa>=Q1) x=3;
    //  else if (qa>=Q2) x=2;
    //  else if (qa>=Q3) x=1;
    //  else x=0;
    //  last += x << pm->dloc; // tmp reuse of delta pos
    //  state->qtot += q*q;
    //  state->qlen++;
    //}

    // Only update delta after 1st base.
    state->delta += (state->prevq != q);
    state->prevq = q;

    state->p--;

    return last & (CTX_SIZE-1);
}

// Build quality stats for qhist and set nsym, do_dedup and do_sel params.
// One_param is -1 to gather stats on all data, or >= 0 to gather data
// on one specific selector parameter.  Used only in TEST_MAIN via
// fqz_manual_parameters at the moment.
void fqz_qual_stats(fqz_slice *s,
                    unsigned char *in, size_t in_size,
                    fqz_param *pm,
                    uint32_t qhist[256],
                    int one_param) {
#define NP 32
    uint32_t qhistb[NP][256] = {{0}};  // both
    uint32_t qhist1[NP][256] = {{0}};  // READ1 only
    uint32_t qhist2[NP][256] = {{0}};  // READ2 only
    uint64_t t1[NP] = {0};             // Count for READ1
    uint64_t t2[NP] = {0};             // COUNT for READ2
    uint32_t avg[2560] = {0};          // Avg qual *and later* avg-to-selector map.

    int dir = 0;
    int last_len = 0;
    int do_dedup = 0;
    size_t rec;
    size_t i, j;
    int num_rec = 0;

    // See what info we've been given.
    // Do we have READ1 / READ2?
    // Do we have selector hidden in the top bits of flag?
    int max_sel = 0;
    int has_r2 = 0;
    for (rec = 0; rec < s->num_records; rec++) {
        if (one_param >= 0 && (s->flags[rec] >> 16) != one_param)
            continue;
        num_rec++;
        if (max_sel < (s->flags[rec] >> 16))
            max_sel = (s->flags[rec] >> 16);
        if (s->flags[rec] & FQZ_FREAD2)
            has_r2 = 1;
    }

    // Dedup detection and histogram stats gathering
    int *avg_qual = calloc((s->num_records+1), sizeof(int));
    if (!avg_qual)
        return;

    rec = i = j = 0;
    while (i < in_size) {
        if (one_param >= 0 && (s->flags[rec] >> 16) != one_param) {
            avg_qual[rec] = 0;
            i += s->len[rec++];
            continue;
        }
        if (rec < s->num_records) {
            j = s->len[rec];
            dir = s->flags[rec] & FQZ_FREAD2 ? 1 : 0;
            if (i > 0 && j == last_len
                && !memcmp(in+i-last_len, in+i, j))
                do_dedup++; // cache which records are dup?
        } else {
            j = in_size - i;
            dir = 0;
        }
        last_len = j;

        uint32_t (*qh)[256] = dir ? qhist2 : qhist1;
        uint64_t *th        = dir ? t2     : t1;

        uint32_t tot = 0;
        for (; i < in_size && j > 0; i++, j--) {
            tot += in[i];
            qhist[in[i]]++;
            qhistb[j & (NP-1)][in[i]]++;
            qh[j & (NP-1)][in[i]]++;
            th[j & (NP-1)]++;
        }
        tot = last_len ? (tot*10.0)/last_len+.5 : 0;

        avg_qual[rec] = tot;
        avg[MIN(2559, tot)]++;

        rec++;
    }
    pm->do_dedup = ((rec+1)/(do_dedup+1) < 500);

    last_len = 0;

    // Unique symbol count
    for (i = pm->max_sym = pm->nsym = 0; i < 256; i++) {
        if (qhist[i])
            pm->max_sym = i, pm->nsym++;
    }


    // Auto tune: does average quality helps us?
    if (pm->do_qa != 0) {
        // Histogram of average qual in avg[]
        // NB: we convert avg[] from count to selector index

        // Few symbols means high compression which means
        // selector bits become more significant fraction.
        // Reduce selector bits by skewing the distribution
        // to not be even binning.
        double qf0 = pm->nsym > 8 ? 0.2 : 0.05;
        double qf1 = pm->nsym > 8 ? 0.5 : 0.22;
        double qf2 = pm->nsym > 8 ? 0.8 : 0.60;

        int total = 0;
        i = 0;
        while (i < 2560) {
            total += avg[i];
            if (total > qf0 * num_rec) {
                //fprintf(stderr, "Q1=%d\n", (int)i);
                break;
            }
            avg[i++] = 0;
        }
        while (i < 2560) {
            total += avg[i];
            if (total > qf1 * num_rec) {
                //fprintf(stderr, "Q2=%d\n", (int)i);
                break;
            }
            avg[i++] = 1;
        }
        while (i < 2560) {
            total += avg[i];
            if (total > qf2 * num_rec) {
                //fprintf(stderr, "Q3=%d\n", (int)i);
                break;
            }
            avg[i++] = 2;
        }
        while (i < 2560)
            avg[i++] = 3;

        // Compute simple entropy of merged signal vs split signal.
        i = 0;
        rec = 0;

        int qbin4[4][NP][256] = {{{0}}};
        int qbin2[2][NP][256] = {{{0}}};
        int qbin1   [NP][256] = {{0}};
        int qcnt4[4][NP] = {{0}};
        int qcnt2[4][NP] = {{0}};
        int qcnt1   [NP] = {0};
        while (i < in_size) {
            if (one_param >= 0 && (s->flags[rec] >> 16) != one_param) {
                i += s->len[rec++];
                continue;
            }
            if ((rec & 7) && rec < s->num_records) {
                // subsample for speed
                i += s->len[rec++];
                continue;
            }
            if (rec < s->num_records)
                j = s->len[rec];
            else
                j = in_size - i;
            last_len = j;

            uint32_t tot = avg_qual[rec];
            int qb4 = avg[MIN(2559, tot)];
            int qb2 = qb4/2;

            for (; i < in_size && j > 0; i++, j--) {
                int x = j & (NP-1);
                qbin4[qb4][x][in[i]]++;  qcnt4[qb4][x]++;
                qbin2[qb2][x][in[i]]++;  qcnt2[qb2][x]++;
                qbin1     [x][in[i]]++;  qcnt1     [x]++;
            }
            rec++;
        }

        double e1 = 0, e2 = 0, e4 = 0;
        for (j = 0; j < NP; j++) {
            for (i = 0; i < 256; i++) {
                if (qbin1   [j][i]) e1 += qbin1   [j][i] * fast_log(qbin1   [j][i] / (double)qcnt1   [j]);
                if (qbin2[0][j][i]) e2 += qbin2[0][j][i] * fast_log(qbin2[0][j][i] / (double)qcnt2[0][j]);
                if (qbin2[1][j][i]) e2 += qbin2[1][j][i] * fast_log(qbin2[1][j][i] / (double)qcnt2[1][j]);
                if (qbin4[0][j][i]) e4 += qbin4[0][j][i] * fast_log(qbin4[0][j][i] / (double)qcnt4[0][j]);
                if (qbin4[1][j][i]) e4 += qbin4[1][j][i] * fast_log(qbin4[1][j][i] / (double)qcnt4[1][j]);
                if (qbin4[2][j][i]) e4 += qbin4[2][j][i] * fast_log(qbin4[2][j][i] / (double)qcnt4[2][j]);
                if (qbin4[3][j][i]) e4 += qbin4[3][j][i] * fast_log(qbin4[3][j][i] / (double)qcnt4[3][j]);
            }
        }
        e1 /= -log(2)/8;
        e2 /= -log(2)/8;
        e4 /= -log(2)/8;
        //fprintf(stderr, "E1=%f E2=%f E4=%f %f\n", e1, e2+s->num_records/8, e4+s->num_records/4, (e4+s->num_records/4)/(e2+s->num_records/8));

        // Note by using the selector we're robbing bits from elsewhere in
        // the context, which may reduce compression better.
        // We don't know how much by, so this is basically a guess!
        // For now we just say need 5% saving here.
        double qm = pm->do_qa > 0 ? 1 : 0.98;
        if ((pm->do_qa == -1 || pm->do_qa >= 4) &&
            e4 + s->num_records/4 < e2*qm + s->num_records/8 &&
            e4 + s->num_records/4 < e1*qm) {
            //fprintf(stderr, "do q4\n");
            for (i = 0; i < s->num_records; i++) {
                //fprintf(stderr, "%d -> %d -> %d, %d\n", (int)i, avg_qual[i], avg[MIN(2559, avg_qual[i])], s->flags[i]>>16);
                s->flags[i] |= avg[MIN(2559, avg_qual[i])] <<16;
            }
            pm->do_sel = 1;
            max_sel = 3;
        } else if ((pm->do_qa == -1 || pm->do_qa >= 2) && e2 + s->num_records/8 < e1*qm) {
            //fprintf(stderr, "do q2\n");
            for (i = 0; i < s->num_records; i++)
                s->flags[i] |= (avg[MIN(2559, avg_qual[i])]>>1) <<16;
            pm->do_sel = 1;
            max_sel = 1;
        }

        if (pm->do_qa == -1) {
            // assume qual, pos, delta in that order.
            if (pm->pbits > 0 && pm->dbits > 0) {
                // 1 from pos/delta
                pm->sloc = pm->dloc-1;
                pm->pbits--;
                pm->dbits--;
                pm->dloc++;
            } else if (pm->dbits >= 2) {
                // 2 from delta
                pm->sloc = pm->dloc;
                pm->dbits -= 2;
                pm->dloc += 2;
            } else if (pm->qbits >= 2) {
                pm->qbits -= 2;
                pm->ploc -= 2;
                pm->sloc = 16-2 - pm->do_r2;
                if (pm->qbits == 6 && pm->qshift == 5)
                    pm->qbits--;
            }
            pm->do_qa = 4;
        }
    }

    // Auto tune: does splitting up READ1 and READ2 help us?
    if (has_r2 || pm->do_r2) { // FIXME: && but debug for now
        double e1 = 0, e2 = 0; // entropy sum

        for (j = 0; j < NP; j++) {
            if (!t1[j] || !t2[j]) continue;
            for (i = 0; i < 256; i++) {
                if (!qhistb[j][i]) continue;
                e1 -= (qhistb[j][i])*log(qhistb[j][i] / (double)(t1[j]+t2[j]));
                if (qhist1[j][i])
                    e2 -= qhist1[j][i] * log(qhist1[j][i] / (double)t1[j]);
                if (qhist2[j][i])
                    e2 -= qhist2[j][i] * log(qhist2[j][i] / (double)t2[j]);
            }
        }
        e1 /= log(2)*8; // bytes
        e2 /= log(2)*8;

        //fprintf(stderr, "read1/2 entropy merge %f split %f\n", e1, e2);

        // Note by using the selector we're robbing bits from elsewhere in
        // the context, which may reduce compression better.
        // We don't know how much by, so this is basically a guess!
        // For now we just say need 5% saving here.
        double qm = pm->do_r2 > 0 ? 1 : 0.95;
        if (e2 + (8+s->num_records/8) < e1*qm) {
            for (rec = 0; rec < s->num_records; rec++) {
                if (one_param >= 0 && (s->flags[rec] >> 16) != one_param)
                    continue;
                int sel = s->flags[rec] >> 16;
                s->flags[rec] =  (s->flags[rec] & 0xffff)
                    | ((s->flags[rec] & FQZ_FREAD2)
                       ? ((sel*2)+1) << 16
                       : ((sel*2)+0) << 16);
                if (max_sel < (s->flags[rec]>>16))
                    max_sel = (s->flags[rec]>>16);
            }
        }
    }

    // We provided explicit selector data or auto-tuned it
    if (max_sel > 0) {
        pm->do_sel = 1;
        pm->max_sel = max_sel;
    }

    free(avg_qual);
}

static inline
int fqz_store_parameters1(fqz_param *pm, unsigned char *comp) {
    int comp_idx = 0, i, j;

    // Starting context
    comp[comp_idx++] = pm->context;
    comp[comp_idx++] = pm->context >> 8;

    comp[comp_idx++] = pm->pflags;
    comp[comp_idx++] = pm->max_sym;

    comp[comp_idx++] = (pm->qbits<<4)|pm->qshift;
    comp[comp_idx++] = (pm->qloc<<4)|pm->sloc;
    comp[comp_idx++] = (pm->ploc<<4)|pm->dloc;

    if (pm->store_qmap) {
        for (i = j = 0; i < 256; i++)
            if (pm->qmap[i] != INT_MAX)
                comp[comp_idx++] = i;
    }

    if (pm->qbits && pm->use_qtab)
        // custom qtab
        comp_idx += store_array(comp+comp_idx, pm->qtab, 256);

    if (pm->pbits && pm->use_ptab)
        // custom ptab
        comp_idx += store_array(comp+comp_idx, pm->ptab, 1024);

    if (pm->dbits && pm->use_dtab)
        // custom dtab
        comp_idx += store_array(comp+comp_idx, pm->dtab, 256);

    return comp_idx;
}

static 
int fqz_store_parameters(fqz_gparams *gp, unsigned char *comp) {
    int comp_idx = 0;
    comp[comp_idx++] = gp->vers; // Format number

    comp[comp_idx++] = gp->gflags;

    if (gp->gflags & GFLAG_MULTI_PARAM)
        comp[comp_idx++] = gp->nparam;

    if (gp->gflags & GFLAG_HAVE_STAB) {
        comp[comp_idx++] = gp->max_sel;
        comp_idx += store_array(comp+comp_idx, gp->stab, 256);
    }

    int i;
    for (i = 0; i < gp->nparam; i++)
        comp_idx += fqz_store_parameters1(&gp->p[i], comp+comp_idx);

    //fprintf(stderr, "Encoded %d bytes of param\n", comp_idx);
    return comp_idx;
}

// Choose a set of parameters based on quality statistics and
// some predefined options (selected via "strat").
static inline
int fqz_pick_parameters(fqz_gparams *gp,
                        int vers,
                        int strat,
                        fqz_slice *s,
                        unsigned char *in,
                        size_t in_size) {
    //approx sqrt(delta), must be sequential
    int dsqr[] = {
        0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
    };
    uint32_t qhist[256] = {0};

    if (strat >= nstrats) strat = nstrats-1;

    // Start with 1 set of parameters.
    // FIXME: add support for multiple params later.
    memset(gp, 0, sizeof(*gp));
    gp->vers = FQZ_VERS;

    if (!(gp->p = calloc(1, sizeof(fqz_param))))
        return -1;
    gp->nparam = 1;
    gp->max_sel = 0;

    if (vers == 3) // V3.0 doesn't store qual in original orientation
        gp->gflags |= GFLAG_DO_REV;

    fqz_param *pm = gp->p;

    // Programmed strategies, which we then amend based on our
    // statistical analysis of the quality stream.
    pm->qbits  = strat_opts[strat][0];
    pm->qshift = strat_opts[strat][1];
    pm->pbits  = strat_opts[strat][2];
    pm->pshift = strat_opts[strat][3];
    pm->dbits  = strat_opts[strat][4];
    pm->dshift = strat_opts[strat][5];
    pm->qloc   = strat_opts[strat][6];
    pm->sloc   = strat_opts[strat][7];
    pm->ploc   = strat_opts[strat][8];
    pm->dloc   = strat_opts[strat][9];

    // Params for controlling behaviour here.
    pm->do_r2 = strat_opts[strat][10];
    pm->do_qa = strat_opts[strat][11];

    // Validity check input lengths and buffer size
    size_t tlen = 0, i;
    for (i = 0; i < s->num_records; i++) {
        if (tlen + s->len[i] > in_size)
            // Oversized buffer
            s->len[i] = in_size - tlen;
        tlen += s->len[i];
    }
    if (s->num_records > 0 && tlen < in_size)
        // Undersized buffer
        s->len[s->num_records-1] += in_size - tlen;

    // Quality metrics, for all recs
    fqz_qual_stats(s, in, in_size, pm, qhist, -1);

    pm->store_qmap = (pm->nsym <= 8 && pm->nsym*2 < pm->max_sym);

    // Check for fixed length.
    uint32_t first_len = s->len[0];
    for (i = 1; i < s->num_records; i++) {
        if (s->len[i] != first_len)
            break;
    }
    pm->fixed_len = (i == s->num_records);
    pm->use_qtab = 0; // unused by current encoder

    if (strat >= nstrats-1)
        goto manually_set; // used in TEST_MAIN for debugging

    if (pm->pshift < 0)
        pm->pshift = MAX(0, log((double)s->len[0]/(1<<pm->pbits))/log(2)+.5);

    if (pm->nsym <= 4) {
        // NovaSeq
        pm->qshift = 2; // qmax 64, although we can store up to 256 if needed
        if (in_size < 5000000) {
            pm->pbits =2;
            pm->pshift=5;
        }
    } else if (pm->nsym <= 8) {
        // HiSeqX
        pm->qbits =MIN(pm->qbits,9);
        pm->qshift=3;
        if (in_size < 5000000)
            pm->qbits =6;
    }

    if (in_size < 300000) {
        pm->qbits=pm->qshift;
        pm->dbits=2;
    }

 manually_set:
//    fprintf(stderr, "-x 0x%x%x%x%x%x%x%x%x%x%x%x%x\n",
//          pm->qbits, pm->qshift,
//          pm->pbits, pm->pshift,
//          pm->dbits, pm->dshift,
//          pm->qloc, pm->sloc, pm->ploc, pm->dloc,
//          pm->do_r2, pm->do_qa);

    for (i = 0; i < sizeof(dsqr)/sizeof(*dsqr); i++)
        if (dsqr[i] > (1<<pm->dbits)-1)
            dsqr[i] = (1<<pm->dbits)-1;

    if (pm->store_qmap) {
        int j;
        for (i = j = 0; i < 256; i++)
            if (qhist[i])
                pm->qmap[i] = j++;
            else
                pm->qmap[i] = INT_MAX;
        pm->max_sym = pm->nsym;
    } else {
        pm->nsym = 255;
        for (i = 0; i < 256; i++)
            pm->qmap[i] = i;
    }
    if (gp->max_sym < pm->max_sym)
        gp->max_sym = pm->max_sym;

    // Produce ptab from pshift.
    if (pm->qbits) {
        for (i = 0; i < 256; i++) {
            pm->qtab[i] = i; // 1:1

            // Alternative mappings:
            //qtab[i] = i > 30 ? MIN(max_sym,i)-15 : i/2;  // eg for 9827 BAM
        }

    }
    pm->qmask = (1<<pm->qbits)-1;

    if (pm->pbits) {
        for (i = 0; i < 1024; i++)
            pm->ptab[i] = MIN((1<<pm->pbits)-1, i>>pm->pshift);

        // Alternatively via analysis of quality distributions we
        // may select a bunch of positions that are special and
        // have a non-uniform ptab[].
        // Manual experimentation on a NovaSeq run saved 2.8% here.
    }

    if (pm->dbits) {
        for (i = 0; i < 256; i++)
            pm->dtab[i] = dsqr[MIN(sizeof(dsqr)/sizeof(*dsqr)-1, i>>pm->dshift)];
    }

    pm->use_ptab = (pm->pbits > 0);
    pm->use_dtab = (pm->dbits > 0);

    pm->pflags =
        (pm->use_qtab   ?PFLAG_HAVE_QTAB :0)|
        (pm->use_dtab   ?PFLAG_HAVE_DTAB :0)|
        (pm->use_ptab   ?PFLAG_HAVE_PTAB :0)|
        (pm->do_sel     ?PFLAG_DO_SEL    :0)|
        (pm->fixed_len  ?PFLAG_DO_LEN    :0)|
        (pm->do_dedup   ?PFLAG_DO_DEDUP  :0)|
        (pm->store_qmap ?PFLAG_HAVE_QMAP :0);

    gp->max_sel = 0;
    if (pm->do_sel) {
        // 2 selectors values, but 1 parameter block.
        // We'll use the sloc instead to encode the selector bits into
        // the context.
        gp->max_sel = 1; // indicator to check recs
        gp->gflags |= GFLAG_HAVE_STAB;
        // NB: stab is already all zero
    }

    if (gp->max_sel && s->num_records) {
        int max = 0;
        for (i = 0; i < s->num_records; i++) {
            if (max < (s->flags[i] >> 16))
                max = (s->flags[i] >> 16);
        }
        gp->max_sel = max;
    }

    return 0;
}

static void fqz_free_parameters(fqz_gparams *gp) {
    if (gp && gp->p) free(gp->p);
}

static int compress_new_read(fqz_slice *s,
                             fqz_state *state,
                             fqz_gparams *gp,
                             fqz_param *pm,
                             fqz_model *model,
                             RangeCoder *rc,
                             unsigned char *in,
                             size_t *in_i, // in[in_i],
                             unsigned int *last) {
    ssize_t rec = state->rec;
    size_t i = *in_i;
    if (pm->do_sel || (gp->gflags & GFLAG_MULTI_PARAM)) {
        state->s = rec < s->num_records
            ? s->flags[rec] >> 16 // reuse spare bits
            : 0;
        SIMPLE_MODEL(256,_encodeSymbol)(&model->sel, rc, state->s);
    } else {
        state->s = 0;
    }
    int x = (gp->gflags & GFLAG_HAVE_STAB) ? gp->stab[state->s] : state->s;
    pm = &gp->p[x];

    int len = s->len[rec];
    if (!pm->fixed_len || state->first_len) {
        SIMPLE_MODEL(256,_encodeSymbol)(&model->len[0], rc, (len>> 0) & 0xff);
        SIMPLE_MODEL(256,_encodeSymbol)(&model->len[1], rc, (len>> 8) & 0xff);
        SIMPLE_MODEL(256,_encodeSymbol)(&model->len[2], rc, (len>>16) & 0xff);
        SIMPLE_MODEL(256,_encodeSymbol)(&model->len[3], rc, (len>>24) & 0xff);
        state->first_len = 0;
    }

    if (gp->gflags & GFLAG_DO_REV) {
        // no need to reverse complement for V4.0 as the core format
        // already has this feature.
        if (s->flags[rec] & FQZ_FREVERSE)
            SIMPLE_MODEL(2,_encodeSymbol)(&model->revcomp, rc, 1);
        else
            SIMPLE_MODEL(2,_encodeSymbol)(&model->revcomp, rc, 0);
    }

    state->rec++;

    state->qtot = 0;
    state->qlen = 0;

    state->p = len;
    state->delta = 0;
    state->qctx = 0;
    state->prevq = 0;

    *last = pm->context;

    if (pm->do_dedup) {
        // Possible dup of previous read?
        if (i && len == state->last_len &&
            !memcmp(in+i-state->last_len, in+i, len)) {
            SIMPLE_MODEL(2,_encodeSymbol)(&model->dup, rc, 1);
            i += len-1;
            state->p = 0;
            *in_i = i;
            return 1; // is a dup
        } else {
            SIMPLE_MODEL(2,_encodeSymbol)(&model->dup, rc, 0);
        }

        state->last_len = len;
    }

    *in_i = i;

    return 0; // not dup
}

static
unsigned char *compress_block_fqz2f(int vers,
                                    int strat,
                                    fqz_slice *s,
                                    unsigned char *in,
                                    size_t in_size,
                                    size_t *out_size,
                                    fqz_gparams *gp) {
    fqz_gparams local_gp;
    int free_params = 0;

    unsigned int last = 0;
    size_t i, j;
    ssize_t rec = 0;

    int comp_idx = 0;
    RangeCoder rc;

    // Pick and store params
    if (!gp) {
        gp = &local_gp;
        if (fqz_pick_parameters(gp, vers, strat, s, in, in_size) < 0)
            return NULL;
        free_params = 1;
    }

    // Worst case scenario assuming random input data and no way to compress
    // is NBytes*growth for some small growth factor (arith_dynamic uses 1.05),
    // plus fixed overheads for the header / params.  Growth can be high
    // here as we're modelling things and pathological cases may trigger a
    // bad probability model.
    //
    // Per read is 4-byte len if not fixed length (but less if avg smaller)
    //             up to 1 byte for selection state (log2(max_sel) bits)
    //             1-bit for reverse flag
    //             1-bit for dup-last flag (but then no quals)
    // Per qual is 1-byte (assuming QMAX==256)
    //
    // Header size is total guess, as depends on params, but it's almost
    // always tiny, so a few K extra should be sufficient.
    //
    // => Total of (s->num_records*4.25 + in_size)*growth + hdr
    int sel_bits = 0, sel = gp->max_sel;
    while (sel) {
        sel_bits++;
        sel >>= 1;
    }
    double len_sz = gp->p[0].fixed_len ? 0.25 : 4.25;
    len_sz += sel_bits / 8.0;
    size_t comp_sz = (s->num_records*len_sz + in_size)*1.1 + 10000;

    unsigned char *comp = (unsigned char *)malloc(comp_sz);
    unsigned char *compe = comp + (size_t)comp_sz;
    if (!comp)
        return NULL;

    //dump_params(gp);
    comp_idx = var_put_u32(comp, compe, in_size);
    comp_idx += fqz_store_parameters(gp, comp+comp_idx);

    fqz_param *pm;

    // Optimise tables to remove shifts in loop (NB: cannot do this in next vers)
    for (j = 0; j < gp->nparam; j++) {
        pm = &gp->p[j];

        for (i = 0; i < 1024; i++)
            pm->ptab[i] <<= pm->ploc;

        for (i = 0; i < 256; i++)
            pm->dtab[i] <<= pm->dloc;
    }

    // Create models and initialise range coder
    fqz_model model;
    if (fqz_create_models(&model, gp) < 0)
        return NULL;

    RC_SetOutput(&rc, (char *)comp+comp_idx);
    RC_SetOutputEnd(&rc, (char *)comp+comp_sz);
    RC_StartEncode(&rc);

    // For CRAM3.1, reverse upfront if needed
    pm = &gp->p[0];
    if (gp->gflags & GFLAG_DO_REV) {
        i = rec = j = 0;
        while (i < in_size) {
            int len = rec < s->num_records-1
                ? s->len[rec] : in_size - i;

            if (s->flags[rec] & FQZ_FREVERSE) {
                // Reverse complement sequence - note: modifies buffer
                int I,J;
                unsigned char *cp = in+i;
                for (I = 0, J = len-1; I < J; I++, J--) {
                    unsigned char c;
                    c = cp[I];
                    cp[I] = cp[J];
                    cp[J] = c;
                }
            }

            i += len;
            rec++;
        }
        rec = 0;
    }

    fqz_state state = {0};
    pm = &gp->p[0];
    state.p = 0;
    state.first_len = 1;
    state.last_len = 0;
    state.rec = rec;

    for (i = 0; i < in_size; i++) {
        if (state.p == 0) {
            if (state.rec >= s->num_records || s->len[state.rec] <= 0) {
                free(comp);
                comp = NULL;
                goto err;
            }

            if (compress_new_read(s, &state, gp, pm, &model, &rc,
                                  in, &i, /*&rec,*/ &last))
                continue;
        }

#if 0
        //                        fqz_qual_stats imp.
        // q40 6.876  6.852       5.96
        // q4  6.566              5.07
        // _Q  1.383              1.11
        unsigned char q = in[i];
        unsigned char qm = pm->qmap[q];

        SIMPLE_MODEL(QMAX,_encodeSymbol)(&model.qual[last], &rc, qm);
        last = fqz_update_ctx(pm, &state, qm);
#else
        //     gcc    clang            gcc+fqz_qual_stats imp.
        // q40 5.033  5.026     -27%   4.137 -38%
        // q4  5.595            -15%   4.011 -36%
        // _Q  1.225            -11%   0.956
        int j = -1;

        while (state.p >= 4 && i+j+4 < in_size) {
            int l1 = last, l2, l3, l4;
            // Model has symbols sorted by frequency, so most common are at
            // start.  So while model is approx 1Kb, the first cache line is
            // a big win.
            mm_prefetch(&model.qual[l1]);
            unsigned char qm1 = pm->qmap[in[i + ++j]];
            last = fqz_update_ctx(pm, &state, qm1); l2 = last;

            mm_prefetch(&model.qual[l2]);
            unsigned char qm2 = pm->qmap[in[i + ++j]];
            last = fqz_update_ctx(pm, &state, qm2); l3 = last;

            mm_prefetch(&model.qual[l3]);
            unsigned char qm3 = pm->qmap[in[i + ++j]];
            last = fqz_update_ctx(pm, &state, qm3); l4 = last;

            mm_prefetch(&model.qual[l4]);
            unsigned char qm4 = pm->qmap[in[i + ++j]];
            last = fqz_update_ctx(pm, &state, qm4);

            SIMPLE_MODEL(QMAX,_encodeSymbol)(&model.qual[l1], &rc, qm1);
            SIMPLE_MODEL(QMAX,_encodeSymbol)(&model.qual[l2], &rc, qm2);
            SIMPLE_MODEL(QMAX,_encodeSymbol)(&model.qual[l3], &rc, qm3);
            SIMPLE_MODEL(QMAX,_encodeSymbol)(&model.qual[l4], &rc, qm4);
        }

        while (state.p > 0) {
            int l2 = last;
            mm_prefetch(&model.qual[last]);
            unsigned char qm = pm->qmap[in[i + ++j]];
            last = fqz_update_ctx(pm, &state, qm);
            SIMPLE_MODEL(QMAX,_encodeSymbol)(&model.qual[l2], &rc, qm);
        }
        i += j;
#endif
    }

    if (RC_FinishEncode(&rc) < 0) {
        free(comp);
        comp = NULL;
        *out_size = 0;
        goto err;
    }

    // For CRAM3.1, undo our earlier reversal step
    rec = state.rec;
    if (gp->gflags & GFLAG_DO_REV) {
        i = rec = j = 0;
        while (i < in_size) {
            int len = rec < s->num_records-1
                ? s->len[rec]
                : in_size - i;

            if (s->flags[rec] & FQZ_FREVERSE) {
                // Reverse complement sequence - note: modifies buffer
                int I,J;
                unsigned char *cp = in+i;
                for (I = 0, J = len-1; I < J; I++, J--) {
                    unsigned char c;
                    c = cp[I];
                    cp[I] = cp[J];
                    cp[J] = c;
                }
            }

            i += len;
            rec++;
        }
    }

    // Clear selector abuse of flags
    for (rec = 0; rec < s->num_records; rec++)
        s->flags[rec] &= 0xffff;

    *out_size = comp_idx + RC_OutSize(&rc);
    //fprintf(stderr, "%d -> %d\n", (int)in_size, (int)*out_size);

 err:
    fqz_destroy_models(&model);
    if (free_params)
        fqz_free_parameters(gp);

    return comp;
}

// Read fqz paramaters.
//
// FIXME: pass in and check in_size.
//
// Returns number of bytes read on success,
//         -1 on failure.
static inline
int fqz_read_parameters1(fqz_param *pm, unsigned char *in, size_t in_size) {
    int in_idx = 0;
    size_t i;

    if (in_size < 7)
        return -1;

    // Starting context
    pm->context = in[in_idx] + (in[in_idx+1]<<8);
    in_idx += 2;

    // Bit flags
    pm->pflags     = in[in_idx++];
    pm->use_qtab   = pm->pflags & PFLAG_HAVE_QTAB;
    pm->use_dtab   = pm->pflags & PFLAG_HAVE_DTAB;
    pm->use_ptab   = pm->pflags & PFLAG_HAVE_PTAB;
    pm->do_sel     = pm->pflags & PFLAG_DO_SEL;
    pm->fixed_len  = pm->pflags & PFLAG_DO_LEN;
    pm->do_dedup   = pm->pflags & PFLAG_DO_DEDUP;
    pm->store_qmap = pm->pflags & PFLAG_HAVE_QMAP;
    pm->max_sym    = in[in_idx++];

    // Sub-context sizes and locations
    pm->qbits      = in[in_idx]>>4;
    pm->qmask      = (1<<pm->qbits)-1;
    pm->qshift     = in[in_idx++]&15;
    pm->qloc       = in[in_idx]>>4;
    pm->sloc       = in[in_idx++]&15;
    pm->ploc       = in[in_idx]>>4;
    pm->dloc       = in[in_idx++]&15;

    // Maps and tables
    if (pm->store_qmap) {
        for (i = 0; i < 256; i++) pm->qmap[i] = INT_MAX; // so dump_map works
        if (in_idx + pm->max_sym > in_size)
            return -1;
        for (i = 0; i < pm->max_sym; i++)
            pm->qmap[i] = in[in_idx++];
    } else {
        for (i = 0; i < 256; i++)
            pm->qmap[i] = i;
    }

    if (pm->qbits) {
        if (pm->use_qtab) {
            int used = read_array(in+in_idx, in_size-in_idx, pm->qtab, 256);
            if (used < 0)
                return -1;
            in_idx += used;
        } else {
            for (i = 0; i < 256; i++)
                pm->qtab[i] = i;
        }
    }

    if (pm->use_ptab) {
        int used = read_array(in+in_idx, in_size-in_idx, pm->ptab, 1024);
        if (used < 0)
            return -1;
        in_idx += used;
    } else {
        for (i = 0; i < 1024; i++)
            pm->ptab[i] = 0;
    }

    if (pm->use_dtab) {
        int used = read_array(in+in_idx, in_size-in_idx, pm->dtab, 256);
        if (used < 0)
            return -1;
        in_idx += used;
    } else {
        for (i = 0; i < 256; i++)
            pm->dtab[i] = 0;
    }

    return in_idx;
}

static
int fqz_read_parameters(fqz_gparams *gp, unsigned char *in, size_t in_size) {
    int in_idx = 0;
    int i;

    if (in_size < 10)
        return -1;

    // Format version
    gp->vers = in[in_idx++];
    if (gp->vers != FQZ_VERS)
        return -1;

    // Global glags
    gp->gflags = in[in_idx++];

    // Number of param blocks and param selector details
    gp->nparam = (gp->gflags & GFLAG_MULTI_PARAM) ? in[in_idx++] : 1;
    if (gp->nparam <= 0)
        return -1;
    gp->max_sel = gp->nparam > 1 ? gp->nparam : 0;

    if (gp->gflags & GFLAG_HAVE_STAB) {
        gp->max_sel = in[in_idx++];
        int used = read_array(in+in_idx, in_size-in_idx, gp->stab, 256);
        if (used < 0)
            goto err;
        in_idx += used;
    } else {
        for (i = 0; i < gp->nparam; i++)
            gp->stab[i] = i;
        for (; i < 256; i++)
            gp->stab[i] = gp->nparam-1;
    }

    // Load the individual parameter locks
    if (!(gp->p = malloc(gp->nparam * sizeof(*gp->p))))
        return -1;

    gp->max_sym = 0;
    for (i = 0; i < gp->nparam; i++) {
        int e = fqz_read_parameters1(&gp->p[i], in + in_idx, in_size-in_idx);
        if (e < 0)
            goto err;
        if (gp->p[i].do_sel && gp->max_sel == 0)
            goto err; // Inconsistent
        in_idx += e;

        if (gp->max_sym < gp->p[i].max_sym)
            gp->max_sym = gp->p[i].max_sym;
    }

    //fprintf(stderr, "Decoded %d bytes of param\n", in_idx);
    return in_idx;

 err:
    fqz_free_parameters(gp);
    gp->nparam = 0;
    return -1;
}

// Handles the state.p==0 section of uncompress_block_fqz2f
static int decompress_new_read(fqz_slice *s,
                               fqz_state *state,
                               fqz_gparams *gp,
                               fqz_param *pm,
                               fqz_model *model,
                               RangeCoder *rc,
                               unsigned char *in, ssize_t *in_i, // in[in_i],
                               unsigned char *uncomp, size_t *out_size,
                               int *rev, char *rev_a, int *len_a,
                               int *lengths, int nlengths) {
    size_t i = *in_i;
    ssize_t rec = state->rec;

    if (pm->do_sel) {
        state->s = SIMPLE_MODEL(256,_decodeSymbol)(&model->sel, rc);
    } else {
        state->s = 0;
    }

    int x = (gp->gflags & GFLAG_HAVE_STAB)
        ? gp->stab[MIN(255, state->s)]
        : state->s;
    if (x >= gp->nparam)
        return -1;
    pm = &gp->p[x];

    unsigned int len = state->last_len;
    if (!pm->fixed_len || state->first_len) {
        len  = SIMPLE_MODEL(256,_decodeSymbol)(&model->len[0], rc);
        len |= SIMPLE_MODEL(256,_decodeSymbol)(&model->len[1], rc)<<8;
        len |= SIMPLE_MODEL(256,_decodeSymbol)(&model->len[2], rc)<<16;
        len |= ((unsigned)SIMPLE_MODEL(256,_decodeSymbol)(&model->len[3], rc))<<24;
        state->first_len = 0;
        state->last_len = len;
    }
    if (len > *out_size-i || len <= 0)
        return -1;

    if (lengths && rec < nlengths)
        lengths[rec] = len;

    if (gp->gflags & GFLAG_DO_REV) {
        *rev = SIMPLE_MODEL(2,_decodeSymbol)(&model->revcomp, rc);
        rev_a[rec] = *rev;
        len_a[rec] = len;
    }

    if (pm->do_dedup) {
        if (SIMPLE_MODEL(2,_decodeSymbol)(&model->dup, rc)) {
            // Dup of last line
            if (len > i)
                return -1;
            memcpy(uncomp+i, uncomp+i-len, len);
            i += len;
            state->p = 0;
            state->rec++;
            *in_i = i;
            return 1; // dup => continue
        }
    }

    state->rec++;
    state->p = len;
    state->delta = 0;
    state->prevq = 0;
    state->qctx = 0;
    state->ctx = pm->context;

    *in_i = i;

    return 0;
}


static
unsigned char *uncompress_block_fqz2f(fqz_slice *s,
                                      unsigned char *in,
                                      size_t in_size,
                                      size_t *out_size,
                                      int *lengths,
                                      int nlengths) {
    fqz_gparams gp;
    fqz_param *pm;
    char *rev_a = NULL;
    int *len_a = NULL;
    memset(&gp, 0, sizeof(gp));

    uint32_t len;
    ssize_t i, rec = 0, in_idx;
    in_idx = var_get_u32(in, in+in_size, &len);
    *out_size = len;

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (len > 100000)
        return NULL;
#endif

    unsigned char *uncomp = NULL;
    RangeCoder rc;
    unsigned int last = 0;

    // Decode parameter blocks
    if ((i = fqz_read_parameters(&gp, in+in_idx, in_size-in_idx)) < 0)
        return NULL;
    //dump_params(&gp);
    in_idx += i;

    // Optimisations to remove shifts from main loop
    for (i = 0; i < gp.nparam; i++) {
        int j;
        pm = &gp.p[i];
        for (j = 0; j < 1024; j++)
            pm->ptab[j] <<= pm->ploc;
        for (j = 0; j < 256; j++)
            pm->dtab[j] <<= pm->dloc;
    }

    // Initialise models and entropy coder
    fqz_model model;
    if (fqz_create_models(&model, &gp) < 0)
        return NULL;

    RC_SetInput(&rc, (char *)in+in_idx, (char *)in+in_size);
    RC_StartDecode(&rc);


    // Allocate buffers
    uncomp = (unsigned char *)malloc(*out_size);
    if (!uncomp)
        goto err;

    int nrec = 1000;
    rev_a = malloc(nrec);
    len_a = malloc(nrec * sizeof(int));
    if (!rev_a || !len_a)
        goto err;

    // Main decode loop
    fqz_state state;
    state.delta = 0;
    state.prevq = 0;
    state.qctx = 0;
    state.p = 0;
    state.s = 0;
    state.first_len = 1;
    state.last_len = 0;
    state.rec = 0;
    state.ctx = last;

    int rev = 0;
    int x = 0;
    pm = &gp.p[x];
    for (i = 0; i < len; ) {
        if (state.rec >= nrec) {
            nrec *= 2;
            rev_a = realloc(rev_a, nrec);
            len_a = realloc(len_a, nrec*sizeof(int));
            if (!rev_a || !len_a)
                goto err;
        }

        if (state.p == 0) {
            int r = decompress_new_read(s, &state, &gp, pm, &model, &rc,
                                        in, &i, uncomp, out_size,
                                        &rev, rev_a, len_a,
                                        lengths, nlengths);
            if (r < 0)
                goto err;
            if (r > 0)
                continue;
            last = state.ctx;
        }

        // Decode and update context
        do {
            unsigned char Q = SIMPLE_MODEL(QMAX,_decodeSymbol)
                (&model.qual[last], &rc);

            last = fqz_update_ctx(pm, &state, Q);
            uncomp[i++] = pm->qmap[Q];
        } while (state.p != 0 && i < len);
    }

    rec = state.rec;
    if (rec >= nrec) {
        nrec *= 2;
        rev_a = realloc(rev_a, nrec);
        len_a = realloc(len_a, nrec*sizeof(int));
        if (!rev_a || !len_a)
            goto err;
    }
    rev_a[rec] = rev;
    len_a[rec] = len;

    if (gp.gflags & GFLAG_DO_REV) {
        for (i = rec = 0; i < len && rec < nrec; i += len_a[rec++]) {
            if (!rev_a[rec])
                continue;

            int I, J;
            unsigned char *cp = uncomp+i;
            for (I = 0, J = len_a[rec]-1; I < J; I++, J--) {
                unsigned char c;
                c = cp[I];
                cp[I] = cp[J];
                cp[J] = c;
            }
        }
    }

    if (RC_FinishDecode(&rc) < 0)
        goto err;

    fqz_destroy_models(&model);
    free(rev_a);
    free(len_a);
    fqz_free_parameters(&gp);

#ifdef TEST_MAIN
    s->num_records = rec;
#endif

    return uncomp;

 err:
    fqz_destroy_models(&model);
    free(rev_a);
    free(len_a);
    fqz_free_parameters(&gp);
    free(uncomp);

    return NULL;
}

char *fqz_compress(int vers, fqz_slice *s, char *in, size_t uncomp_size,
                   size_t *comp_size, int strat, fqz_gparams *gp) {
    if (uncomp_size > INT_MAX) {
        *comp_size = 0;
        return NULL;
    }

    return (char *)compress_block_fqz2f(vers, strat, s, (unsigned char *)in,
                                        uncomp_size, comp_size, gp);
}

char *fqz_decompress(char *in, size_t comp_size, size_t *uncomp_size,
                     int *lengths, int nlengths) {
    return (char *)uncompress_block_fqz2f(NULL, (unsigned char *)in,
                                          comp_size, uncomp_size, lengths, nlengths);
}
