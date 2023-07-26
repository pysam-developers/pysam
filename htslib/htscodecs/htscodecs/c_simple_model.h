/*
 * Copyright (c) 2012, 2018-2019 Genome Research Ltd.
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

#include <stdint.h>
#include "c_range_coder.h"

/*
 *--------------------------------------------------------------------------
 * A simple frequency model.
 *
 * Define NSYM to be an integer value before including this file.
 * It will then generate types and functions specific to that
 * maximum number of symbols.
 *
 * This keeps a list of symbols and their frequencies, approximately
 * sorted by symbol frequency. We allow for a single symbol to periodically
 * move up the list when emitted, effectively doing a single step of
 * bubble sort periodically. This means it's largely the same complexity
 * irrespective of alphabet size.
 * It's more efficient on strongly biased distributions than random data.
 *
 * There is no escape symbol, so the model is tailored to relatively
 * stationary samples (although we do have occasional normalisation to
 * avoid frequency counters getting too high).
 *--------------------------------------------------------------------------
 */

//-----------------------------------------------------------------------------
// Bits we want included once only - constants, types, etc
#ifndef C_SIMPLE_MODEL_H
#define C_SIMPLE_MODEL_H

#define MAX_FREQ (1<<16)-17
#define PASTE3(a,b,c) a##b##c
#define SIMPLE_MODEL(a,b) PASTE3(SIMPLE_MODEL,a,b)
#define STEP 16
typedef struct {
    uint16_t Freq;
    uint16_t Symbol;
} SymFreqs;
#endif /* C_SIMPLE_MODEL_H */


//-----------------------------------------------------------------------------
// Bits we regenerate for each NSYM value.

typedef struct {
    uint32_t TotFreq;  // Total frequency

    // Array of Symbols approximately sorted by Freq. 
    SymFreqs sentinel, F[NSYM+1], terminal;
} SIMPLE_MODEL(NSYM,_);


static inline void SIMPLE_MODEL(NSYM,_init)(SIMPLE_MODEL(NSYM,_) *m, int max_sym) {
    int i;
    
    for (i=0; i<max_sym; i++) {
        m->F[i].Symbol = i;
        m->F[i].Freq   = 1;
    }
    for (; i<NSYM; i++) {
        m->F[i].Symbol = i;
        m->F[i].Freq   = 0;
    }

    m->TotFreq         = max_sym;
    m->sentinel.Symbol = 0;
    m->sentinel.Freq   = MAX_FREQ; // Always first; simplifies sorting.
    m->terminal.Symbol = 0;
    m->terminal.Freq   = MAX_FREQ;
    m->F[NSYM].Freq    = 0; // terminates normalize() loop. See below.
}


static inline void SIMPLE_MODEL(NSYM,_normalize)(SIMPLE_MODEL(NSYM,_) *m) {
    SymFreqs *s;

    /* Faster than F[i].Freq for 0 <= i < NSYM */
    m->TotFreq=0;
    for (s = m->F; s->Freq; s++) {
        s->Freq -= s->Freq>>1;
        m->TotFreq += s->Freq;
    }
}

static inline void SIMPLE_MODEL(NSYM,_encodeSymbol)(SIMPLE_MODEL(NSYM,_) *m,
                                                    RangeCoder *rc, uint16_t sym) {
    SymFreqs *s = m->F;
    uint32_t AccFreq  = 0;

    while (s->Symbol != sym)
        AccFreq += s++->Freq;

    RC_Encode(rc, AccFreq, s->Freq, m->TotFreq);
    s->Freq    += STEP;
    m->TotFreq += STEP;

    if (m->TotFreq > MAX_FREQ)
        SIMPLE_MODEL(NSYM,_normalize)(m);

    /* Keep approx sorted */
    if (s[0].Freq > s[-1].Freq) {
        SymFreqs t = s[0];
        s[0] = s[-1];
        s[-1] = t;
    }
}

static inline uint16_t SIMPLE_MODEL(NSYM,_decodeSymbol)(SIMPLE_MODEL(NSYM,_) *m, RangeCoder *rc) {
    SymFreqs* s = m->F;
    uint32_t freq = RC_GetFreq(rc, m->TotFreq);
    uint32_t AccFreq;

    if (freq > MAX_FREQ)
        return 0; // error

    for (AccFreq = 0; (AccFreq += s->Freq) <= freq; s++)
        ;
    if (s - m->F > NSYM)
        return 0; // error

    AccFreq -= s->Freq;

    RC_Decode(rc, AccFreq, s->Freq, m->TotFreq);
    s->Freq    += STEP;
    m->TotFreq += STEP;

    if (m->TotFreq > MAX_FREQ)
        SIMPLE_MODEL(NSYM,_normalize)(m);

    /* Keep approx sorted */
    if (s[0].Freq > s[-1].Freq) {
        SymFreqs t = s[0];
        s[0] = s[-1];
        s[-1] = t;
        return t.Symbol;
    }

    return s->Symbol;
}
