/* The MIT License

   Copyright (c) 2014-2016 Genome Research Ltd.

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

#ifndef __HMM_H__
#define __HMM_H__

#define MAT(matrix,ndim,i,j) (matrix)[(ndim)*(i)+(j)]       // P(i|j), that is, transition j->i

typedef struct _hmm_t hmm_t;

typedef void (*set_tprob_f) (hmm_t *hmm, uint32_t prev_pos, uint32_t pos, void *data, double *tprob);

/**
 *   hmm_init() - initialize HMM
 *   @nstates:  number of states
 *   @tprob:    transition probabilities matrix (nstates x nstates), for elements ordering
 *              see the MAT macro above.
 *   @ntprob:   number of precalculated tprob matrices or 0 for constant probs, independent
 *              of distance
 */
hmm_t *hmm_init(int nstates, double *tprob, int ntprob);
void hmm_set_tprob(hmm_t *hmm, double *tprob, int ntprob);

#define HMM_VIT 1
#define HMM_FWD 2
#define HMM_BWD 4

/**
 *   hmm_init_states() - initial state probabilities
 *   @probs:  initial state probabilities or NULL to reset to default
 *
 *   If uncalled, all states are initialized with the same likelihood
 */
void hmm_init_states(hmm_t *hmm, double *probs);

/**
 *   hmm_snapshot() - take the model's snapshot, intended for sliding HMM
 *   @snapshot: NULL or snapshot returned by previous hmm_snapshot() call, must be free()-ed by the caller
 *   @pos:      take the snapshot at this position
 *
 *   If both restore() and snapshot() are needed, restore() must be called first.
 */
void *hmm_snapshot(hmm_t *hmm, void *snapshot, uint32_t pos);

/**
 *   hmm_restore() - restore model's snapshot, intended for sliding HMM
 *   @snapshot: snapshot returned by hmm_snapshot() call or NULL to reset
 *   @isite:    take the snapshot at i-th step
 *
 *   If both restore() and snapshot() are needed, restore() must be called first.
 */
void hmm_restore(hmm_t *hmm, void *snapshot);
void hmm_reset(hmm_t *hmm, void *snapshot);

/**
 *   hmm_get_tprob() - return the array of transition matrices, precalculated
 *      to ntprob positions. The first matrix is the initial tprob matrix
 *      set by hmm_init() or hmm_set_tprob()
 */
double *hmm_get_tprob(hmm_t *hmm);
int hmm_get_nstates(hmm_t *hmm);

/**
 *   hmm_set_tprob_func() - custom setter of transition probabilities
 */
void hmm_set_tprob_func(hmm_t *hmm, set_tprob_f set_tprob, void *data);

/**
 *   hmm_run_viterbi() - run Viterbi algorithm
 *   @nsites:   number of sites 
 *   @eprob:    emission probabilities for each site and state (nsites x nstates)
 *   @sites:    list of positions
 *
 *   When done, hmm->vpath[] contains the calculated Viterbi path. The states
 *   are indexed starting from 0, a state at i-th site can be accessed as
 *   vpath[nstates*i].
 */
void hmm_run_viterbi(hmm_t *hmm, int nsites, double *eprob, uint32_t *sites);

/**
 *   hmm_get_viterbi_path() - the viterbi path: state at ith site is the
 *      (nstates*isite)-th element
 */
uint8_t *hmm_get_viterbi_path(hmm_t *hmm);

/**
 *   hmm_run_fwd_bwd() - run the forward-backward algorithm
 *   @nsites:   number of sites 
 *   @eprob:    emission probabilities for each site and state (nsites x nstates)
 *   @sites:    list of positions
 */
void hmm_run_fwd_bwd(hmm_t *hmm, int nsites, double *eprob, uint32_t *sites);

/**
 *   hmm_get_fwd_bwd_prob() - the probability of i-th state at j-th site can
 *      be accessed as fwd_bwd[j*nstates+i].
 */
double *hmm_get_fwd_bwd_prob(hmm_t *hmm);

/**
 *   hmm_run_baum_welch() - run one iteration of Baum-Welch algorithm
 *   @nsites:   number of sites 
 *   @eprob:    emission probabilities for each site and state (nsites x nstates)
 *   @sites:    list of positions
 *
 *   Same as hmm_run_fwd_bwd, in addition a pointer to a matrix with the new
 *   transition probabilities is returned. In this version, emission
 *   probabilities are not updated.
 */
double *hmm_run_baum_welch(hmm_t *hmm, int nsites, double *eprob, uint32_t *sites);

void hmm_destroy(hmm_t *hmm);

#endif

