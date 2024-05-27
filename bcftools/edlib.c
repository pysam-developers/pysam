/*
 * A cut down C translated of the C++ edlib.cpp file.
 * Taken from edlib v0.1.0-166-g931be2b
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "edlib.h"

typedef uint64_t Word;
static const int WORD_SIZE = 64; // Size of Word in bits
static const Word WORD_1 = (Word)1;
static const Word HIGH_BIT_MASK = 1LL << 63;  // 100..00
//#define MAX_UCHAR 255
#define MAX_UCHAR 7 // better cache usage for our data

#ifndef MAX
#define MAX(a,b) ((a)>(b)?(a):(b))
#endif

typedef struct Block {
    Word P;  // Pvin
    Word M;  // Mvin
    int score; // score of last cell in block;
} Block;


/**
 * Defines equality relation on alphabet characters.
 * By default each character is always equal only to itself, but you can also provide additional equalities.
 */
typedef struct EqualityDefinition {
    bool matrix[MAX_UCHAR + 1][MAX_UCHAR + 1];
} EqualityDefinition;

static EqualityDefinition *
CreateEqualityDefinition(const char *alphabet, int alphabet_size,
			 const EdlibEqualityPair* additionalEqualities,
			 const int additionalEqualitiesLength) {
    EqualityDefinition *ed = malloc(sizeof(*ed));

    for (size_t i = 0; i < alphabet_size; i++) {
	for (size_t j = 0; j < alphabet_size; j++) {
	    ed->matrix[i][j] = (i == j);
	}
    }
    if (additionalEqualities != NULL) {
	for (int i = 0; i < additionalEqualitiesLength; i++) {
	    const char *firstTransformed = strchr(alphabet, additionalEqualities[i].first);
	    const char *secondTransformed = strchr(alphabet, additionalEqualities[i].second);
	    if (firstTransformed && alphabet_size) {
		ed->matrix[firstTransformed - alphabet][secondTransformed - alphabet] =
		ed->matrix[secondTransformed - alphabet][firstTransformed - alphabet]
		    = true;
	    }
	}
    }

    return ed;
}

/**
 * @param a  Element from transformed sequence.
 * @param b  Element from transformed sequence.
 * @return True if a and b are defined as equal, false otherwise.
 */
static inline const /* attribute pure or const? */
bool equalityDefinition_areEqual(const EqualityDefinition *ed, unsigned char a, unsigned char b) {
    return ed->matrix[a][b];
}

static int myersCalcEditDistanceSemiGlobal(const Word* Peq, int W, int maxNumBlocks,
                                           int queryLength,
                                           const unsigned char* target, int targetLength,
                                           int k, EdlibAlignMode mode,
                                           int* bestScore_, int** positions_, int* numPositions_);

static char *transformSequences(const char* queryOriginal, int queryLength,
				const char* targetOriginal, int targetLength,
				unsigned char** queryTransformed,
				unsigned char** targetTransformed,
                                int *alphabet_size);

static inline int ceilDiv(int x, int y);

static inline unsigned char* createReverseCopy(const unsigned char* seq, int length);

static inline Word* buildPeq(const int alphabetLength,
                             const unsigned char* query,
                             const int queryLength,
                             const EqualityDefinition* equalityDefinition);


/**
 * Main edlib method.
 */
EdlibAlignResult edlibAlign(const char* const queryOriginal, const int queryLength,
			    const char* const targetOriginal, const int targetLength,
			    const EdlibAlignConfig config) {
    EdlibAlignResult result;
    result.status = EDLIB_STATUS_OK;
    result.editDistance = -1;
    result.endLocations = result.startLocations = NULL;
    result.numLocations = 0;
    result.alignment = NULL;
    result.alignmentLength = 0;
    result.alphabetLength = 0;

    /*------------ TRANSFORM SEQUENCES AND RECOGNIZE ALPHABET -----------*/
    unsigned char* query, * target;
    int alphabet_size;
    char *alphabet = transformSequences(queryOriginal, queryLength, targetOriginal, targetLength,
                                         &query, &target, &alphabet_size);
    result.alphabetLength = alphabet_size;
    /*-------------------------------------------------------*/

    // Handle special situation when at least one of the sequences has length 0.
    if (queryLength == 0 || targetLength == 0) {
        if (config.mode == EDLIB_MODE_NW) {
            result.editDistance = MAX(queryLength, targetLength);
            result.endLocations = malloc(sizeof(int) * 1);
            result.endLocations[0] = targetLength - 1;
            result.numLocations = 1;
        } else if (config.mode == EDLIB_MODE_SHW || config.mode == EDLIB_MODE_HW) {
            result.editDistance = queryLength;
            result.endLocations = malloc(sizeof(int) * 1);
            result.endLocations[0] = -1;
            result.numLocations = 1;
        } else {
            result.status = EDLIB_STATUS_ERROR;
        }

        free(query);
        free(target);
        free(alphabet);
        return result;
    }

    /*--------------------- INITIALIZATION ------------------*/
    int maxNumBlocks = ceilDiv(queryLength, WORD_SIZE); // bmax in Myers
    int W = maxNumBlocks * WORD_SIZE - queryLength; // number of redundant cells in last level blocks
    EqualityDefinition *equalityDefinition =
	CreateEqualityDefinition(alphabet, alphabet_size, config.additionalEqualities, config.additionalEqualitiesLength);
    Word* Peq = buildPeq(alphabet_size, query, queryLength, equalityDefinition);
    /*-------------------------------------------------------*/

    /*------------------ MAIN CALCULATION -------------------*/
    // TODO: Store alignment data only after k is determined? That could make things faster.
//    int positionNW; // Used only when mode is NW.
//    AlignmentData* alignData = NULL;
    bool dynamicK = false;
    int k = config.k;
    if (k < 0) { // If valid k is not given, auto-adjust k until solution is found.
        dynamicK = true;
        k = WORD_SIZE; // Gives better results than smaller k.
    }

    do {
        if (config.mode == EDLIB_MODE_HW || config.mode == EDLIB_MODE_SHW) {
            myersCalcEditDistanceSemiGlobal(Peq, W, maxNumBlocks,
                                            queryLength, target, targetLength,
                                            k, config.mode, &(result.editDistance),
                                            &(result.endLocations), &(result.numLocations));
        } else {  // mode == EDLIB_MODE_NW
//            myersCalcEditDistanceNW(Peq, W, maxNumBlocks,
//                                    queryLength, target, targetLength,
//                                    k, &(result.editDistance), &positionNW,
//                                    false, &alignData, -1);
        }
        k *= 2;
    } while(dynamicK && result.editDistance == -1);

    if (result.editDistance >= 0) {  // If there is solution.
        // If NW mode, set end location explicitly.
        if (config.mode == EDLIB_MODE_NW) {
            result.endLocations = malloc(sizeof(int) * 1);
            result.endLocations[0] = targetLength - 1;
            result.numLocations = 1;
        }

        // Find starting locations.
        if (config.task == EDLIB_TASK_LOC || config.task == EDLIB_TASK_PATH) {
            result.startLocations = malloc(result.numLocations * sizeof(int));
            if (config.mode == EDLIB_MODE_HW) {  // If HW, I need to calculate start locations.
                const unsigned char* rTarget = createReverseCopy(target, targetLength);
                const unsigned char* rQuery  = createReverseCopy(query, queryLength);
                // Peq for reversed query.
                Word* rPeq = buildPeq(alphabet_size, rQuery, queryLength, equalityDefinition);
                for (int i = 0; i < result.numLocations; i++) {
                    int endLocation = result.endLocations[i];
                    if (endLocation == -1) {
                        // NOTE: Sometimes one of optimal solutions is that query starts before target, like this:
                        //                       AAGG <- target
                        //                   CCTT     <- query
                        //   It will never be only optimal solution and it does not happen often, however it is
                        //   possible and in that case end location will be -1. What should we do with that?
                        //   Should we just skip reporting such end location, although it is a solution?
                        //   If we do report it, what is the start location? -4? -1? Nothing?
                        // TODO: Figure this out. This has to do in general with how we think about start
                        //   and end locations.
                        //   Also, we have alignment later relying on this locations to limit the space of it's
                        //   search -> how can it do it right if these locations are negative or incorrect?
                        result.startLocations[i] = 0;  // I put 0 for now, but it does not make much sense.
                    } else {
                        int bestScoreSHW, numPositionsSHW;
                        int* positionsSHW;
                        myersCalcEditDistanceSemiGlobal(
                                rPeq, W, maxNumBlocks,
                                queryLength, rTarget + targetLength - endLocation - 1, endLocation + 1,
                                result.editDistance, EDLIB_MODE_SHW,
                                &bestScoreSHW, &positionsSHW, &numPositionsSHW);
                        // Taking last location as start ensures that alignment will not start with insertions
                        // if it can start with mismatches instead.
                        result.startLocations[i] = endLocation - positionsSHW[numPositionsSHW - 1];
                        free(positionsSHW);
                    }
                }
                free((void *)rTarget);
                free((void *)rQuery);
                free(rPeq);
            } else {  // If mode is SHW or NW
                for (int i = 0; i < result.numLocations; i++) {
                    result.startLocations[i] = 0;
                }
            }
        }
    }
    /*-------------------------------------------------------*/

    //--- Free memory ---//
    free(Peq);
    free(query);
    free(target);
    free(alphabet);
    free(equalityDefinition);
//    DestroyAlignmentData(alignData);
    //-------------------//

    return result;
}

/**
 * Build Peq table for given query and alphabet.
 * Peq is table of dimensions alphabetLength+1 x maxNumBlocks.
 * Bit i of Peq[s * maxNumBlocks + b] is 1 if i-th symbol from block b of query equals symbol s, otherwise it is 0.
 * NOTICE: free returned array with free()!
 */
static inline Word* buildPeq(const int alphabetLength,
                             const unsigned char* const query,
                             const int queryLength,
                             const EqualityDefinition* equalityDefinition) {
    int maxNumBlocks = ceilDiv(queryLength, WORD_SIZE);
    // table of dimensions alphabetLength+1 x maxNumBlocks. Last symbol is wildcard.
    Word* Peq = malloc((alphabetLength + 1) * maxNumBlocks * sizeof(*Peq));

    // Build Peq (1 is match, 0 is mismatch). NOTE: last column is wildcard(symbol that matches anything) with just 1s
    // Optimised Peq building avoiding branching.
    for (int symbol = 0; symbol < alphabetLength; symbol++) {
        for (int b = 0; b < maxNumBlocks; b++) {
            Word PeqW = 0;
            for (int r = (b+1) * WORD_SIZE - 1; r >= b * WORD_SIZE; r--) {
                PeqW = (PeqW<<1)
                     + (r >= queryLength
                        || equalityDefinition_areEqual(equalityDefinition,
                                                       query[r], symbol)); 
            }
            Peq[symbol * maxNumBlocks + b] = PeqW;
        }
    }
    {
        int symbol = alphabetLength;
        for (int b = 0; b < maxNumBlocks; b++) {
            // Last symbol is wildcard, so it is all 1s
            Peq[symbol * maxNumBlocks + b] = (Word)-1;
        }
    }

    return Peq;
}


/**
 * Returns new sequence that is reverse of given sequence.
 * Free returned array with free()
 */
static inline unsigned char* createReverseCopy(const unsigned char* const seq, const int length) {
    unsigned char* rSeq = malloc(length);
    for (int i = 0; i < length; i++) {
        rSeq[i] = seq[length - i - 1];
    }
    return rSeq;
}

/**
 * Corresponds to Advance_Block function from Myers.
 * Calculates one word(block), which is part of a column.
 * Highest bit of word (one most to the left) is most bottom cell of block from column.
 * Pv[i] and Mv[i] define vin of cell[i]: vin = cell[i] - cell[i-1].
 * @param [in] Pv  Bitset, Pv[i] == 1 if vin is +1, otherwise Pv[i] == 0.
 * @param [in] Mv  Bitset, Mv[i] == 1 if vin is -1, otherwise Mv[i] == 0.
 * @param [in] Eq  Bitset, Eq[i] == 1 if match, 0 if mismatch.
 * @param [in] hin  Will be +1, 0 or -1.
 * @param [out] PvOut  Bitset, PvOut[i] == 1 if vout is +1, otherwise PvOut[i] == 0.
 * @param [out] MvOut  Bitset, MvOut[i] == 1 if vout is -1, otherwise MvOut[i] == 0.
 * @param [out] hout  Will be +1, 0 or -1.
 */
static inline int calculateBlock(Word Pv, Word Mv, Word Eq, const int hin,
                                 Word *PvOut, Word *MvOut) {
    // hin can be 1, -1 or 0.
    // 1  -> 00...01
    // 0  -> 00...00
    // -1 -> 11...11 (2-complement)

    Word hinIsNeg = (Word)(hin >> 2) & WORD_1; // 00...001 if hin is -1, 00...000 if 0 or 1

    Word Xv = Eq | Mv;
    // This is instruction below written using 'if': if (hin < 0) Eq |= (Word)1;
    Eq |= hinIsNeg;
    Word Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;

    Word Ph = Mv | ~(Xh | Pv);
    Word Mh = Pv & Xh;

    int hout = 0;
    // This is instruction below written using 'if': if (Ph & HIGH_BIT_MASK) hout = 1;
    hout = (Ph & HIGH_BIT_MASK) >> (WORD_SIZE - 1);
    // This is instruction below written using 'if': if (Mh & HIGH_BIT_MASK) hout = -1;
    hout -= (Mh & HIGH_BIT_MASK) >> (WORD_SIZE - 1);

    Ph <<= 1;
    Mh <<= 1;

    // This is instruction below written using 'if': if (hin < 0) Mh |= (Word)1;
    Mh |= hinIsNeg;
    // This is instruction below written using 'if': if (hin > 0) Ph |= (Word)1;
    Ph |= (Word)((hin + 1) >> 1);

    *PvOut = Mh | ~(Xv | Ph);
    *MvOut = Ph & Xv;

    return hout;
}

/**
 * Does ceiling division x / y.
 * Note: x and y must be non-negative and x + y must not overflow.
 */
static inline int ceilDiv(const int x, const int y) {
    return x % y ? x / y + 1 : x / y;
}

static inline int min(const int x, const int y) {
    return x < y ? x : y;
}


/**
 * @param [in] block
 * @return Values of cells in block, starting with bottom cell in block.
 */
static inline int *getBlockCellValues(const Block block) {
    int *scores = malloc(WORD_SIZE * sizeof(*scores));
    int score = block.score;
    Word mask = HIGH_BIT_MASK;
    for (int i = 0; i < WORD_SIZE - 1; i++) {
        scores[i] = score;
        if (block.P & mask) score--;
        if (block.M & mask) score++;
        mask >>= 1;
    }
    scores[WORD_SIZE - 1] = score;
    return scores;
}

/**
 * @param [in] block
 * @param [in] k
 * @return True if all cells in block have value larger than k, otherwise false.
 */
static inline bool allBlockCellsLarger(const Block block, const int k) {
    int *scores = getBlockCellValues(block);
    for (int i = 0; i < WORD_SIZE; i++) {
        if (scores[i] <= k) {
            free(scores);
            return false;
        }
    }

    free(scores);
    return true;
}


/**
 * Uses Myers' bit-vector algorithm to find edit distance for one of semi-global alignment methods.
 * @param [in] Peq  Query profile.
 * @param [in] W  Size of padding in last block.
 *                TODO: Calculate this directly from query, instead of passing it.
 * @param [in] maxNumBlocks  Number of blocks needed to cover the whole query.
 *                           TODO: Calculate this directly from query, instead of passing it.
 * @param [in] queryLength
 * @param [in] target
 * @param [in] targetLength
 * @param [in] k
 * @param [in] mode  EDLIB_MODE_HW or EDLIB_MODE_SHW
 * @param [out] bestScore_  Edit distance.
 * @param [out] positions_  Array of 0-indexed positions in target at which best score was found.
                            Make sure to free this array with free().
 * @param [out] numPositions_  Number of positions in the positions_ array.
 * @return Status.
 */
static int myersCalcEditDistanceSemiGlobal(
        const Word* const Peq, const int W, const int maxNumBlocks,
        const int queryLength,
        const unsigned char* const target, const int targetLength,
        int k, const EdlibAlignMode mode,
        int* const bestScore_, int** const positions_, int* const numPositions_) {
    *positions_ = NULL;
    *numPositions_ = 0;

    // firstBlock is 0-based index of first block in Ukkonen band.
    // lastBlock is 0-based index of last block in Ukkonen band.
    int firstBlock = 0;
    int lastBlock = min(ceilDiv(k + 1, WORD_SIZE), maxNumBlocks) - 1; // y in Myers
    Block *bl; // Current block

    Block* blocks = malloc(maxNumBlocks * sizeof(*blocks));

    // For HW, solution will never be larger then queryLength.
    if (mode == EDLIB_MODE_HW) {
        k = min(queryLength, k);
    }

    // Each STRONG_REDUCE_NUM column is reduced in more expensive way.
    // This gives speed up of about 2 times for small k.
    const int STRONG_REDUCE_NUM = 2048;

    // Initialize P, M and score
    bl = blocks;
    for (int b = 0; b <= lastBlock; b++) {
        bl->score = (b + 1) * WORD_SIZE;
        bl->P = (Word)(-1); // All 1s
        bl->M = (Word)(0);
        bl++;
    }

    int bestScore = -1;
#define MAX_POS 100  // maximum number of positions returned.
    int positions[MAX_POS];
    int npositions = 0;
    const int startHout = mode == EDLIB_MODE_HW ? 0 : 1; // If 0 then gap before query is not penalized;
    const unsigned char* targetChar = target;
    for (int c = 0; c < targetLength; c++) { // for each column
        const Word* Peq_c = Peq + (*targetChar) * maxNumBlocks;

        //----------------------- Calculate column -------------------------//
        int hout = startHout;
        bl = blocks + firstBlock;
        Peq_c += firstBlock;
        for (int b = firstBlock; b <= lastBlock; b++) {
            hout = calculateBlock(bl->P, bl->M, *Peq_c, hout, &bl->P, &bl->M);
            bl->score += hout;
            bl++; Peq_c++;
        }
        bl--; Peq_c--;
        //------------------------------------------------------------------//

        //---------- Adjust number of blocks according to Ukkonen ----------//
        if ((lastBlock < maxNumBlocks - 1) && (bl->score - hout <= k) // bl is pointing to last block
            && ((*(Peq_c + 1) & WORD_1) || hout < 0)) { // Peq_c is pointing to last block
            // If score of left block is not too big, calculate one more block
            lastBlock++; bl++; Peq_c++;
            bl->P = (Word)(-1); // All 1s
            bl->M = (Word)(0);
            bl->score = (bl - 1)->score - hout + WORD_SIZE + calculateBlock(bl->P, bl->M, *Peq_c, hout, &bl->P, &bl->M);
        } else {
            while (lastBlock >= firstBlock && bl->score >= k + WORD_SIZE) {
                lastBlock--; bl--; Peq_c--;
            }
        }

        // Every some columns, do some expensive but also more efficient block reducing.
        // This is important!
        //
        // Reduce the band by decreasing last block if possible.
        if (c % STRONG_REDUCE_NUM == 0) {
            while (lastBlock >= 0 && lastBlock >= firstBlock && allBlockCellsLarger(*bl, k)) {
                lastBlock--; bl--; Peq_c--;
            }
        }
        // For HW, even if all cells are > k, there still may be solution in next
        // column because starting conditions at upper boundary are 0.
        // That means that first block is always candidate for solution,
        // and we can never end calculation before last column.
        if (mode == EDLIB_MODE_HW && lastBlock == -1) {
            lastBlock++; bl++; Peq_c++;
        }

        // Reduce band by increasing first block if possible. Not applicable to HW.
        if (mode != EDLIB_MODE_HW) {
            while (firstBlock <= lastBlock && blocks[firstBlock].score >= k + WORD_SIZE) {
                firstBlock++;
            }
            if (c % STRONG_REDUCE_NUM == 0) { // Do strong reduction every some blocks
                while (firstBlock <= lastBlock && allBlockCellsLarger(blocks[firstBlock], k)) {
                    firstBlock++;
                }
            }
        }

        // If band stops to exist finish
        if (lastBlock < firstBlock) {
            *bestScore_ = bestScore;
            if (bestScore != -1) {
                *positions_ = malloc(npositions * sizeof(int));
                *numPositions_ = npositions;
                memcpy(*positions_, positions, npositions * sizeof(int));
            }
            free(blocks);
            return EDLIB_STATUS_OK;
        }
        //------------------------------------------------------------------//

        //------------------------- Update best score ----------------------//
        if (lastBlock == maxNumBlocks - 1) {
            int colScore = bl->score;
            if (colScore <= k) { // Scores > k dont have correct values (so we cannot use them), but are certainly > k.
                // NOTE: Score that I find in column c is actually score from column c-W
                if (bestScore == -1 || colScore <= bestScore) {
                    if (colScore != bestScore) {
			npositions = 0;
                        bestScore = colScore;
                        // Change k so we will look only for equal or better
                        // scores then the best found so far.
                        k = bestScore;
                    }
		    if (npositions < MAX_POS)
			positions[npositions++] = c - W;
                }
            }
        }
        //------------------------------------------------------------------//

        targetChar++;
    }


    // Obtain results for last W columns from last column.
    if (lastBlock == maxNumBlocks - 1) {
	int *blockScores = getBlockCellValues(*bl);
        for (int i = 0; i < W; i++) {
            int colScore = blockScores[i + 1];
            if (colScore <= k && (bestScore == -1 || colScore <= bestScore)) {
                if (colScore != bestScore) {
                    npositions = 0;
                    k = bestScore = colScore;
                }
		if (npositions < MAX_POS)
		    positions[npositions++] = targetLength - W + i;
            }
        }
        free(blockScores);
    }

    *bestScore_ = bestScore;
    if (bestScore != -1) {
        *positions_ = malloc(npositions * sizeof(int));
        *numPositions_ = npositions;
        memcpy(*positions_, positions, npositions * sizeof(int));
    }

    free(blocks);
    return EDLIB_STATUS_OK;
}


/**
 * Takes char query and char target, recognizes alphabet and transforms them into unsigned char sequences
 * where elements in sequences are not any more letters of alphabet, but their index in alphabet.
 * Most of internal edlib functions expect such transformed sequences.
 * This function will allocate queryTransformed and targetTransformed, so make sure to free them when done.
 * Example:
 *   Original sequences: "ACT" and "CGT".
 *   Alphabet would be recognized as "ACTG". Alphabet length = 4.
 *   Transformed sequences: [0, 1, 2] and [1, 3, 2].
 * @param [in] queryOriginal
 * @param [in] queryLength
 * @param [in] targetOriginal
 * @param [in] targetLength
 * @param [out] queryTransformed  It will contain values in range [0, alphabet length - 1].
 * @param [out] targetTransformed  It will contain values in range [0, alphabet length - 1].
 * @return  Alphabet as a string of unique characters, where index of each character is its value in transformed
 *          sequences.
 */
static char *transformSequences(const char* const queryOriginal, const int queryLength,
				const char* const targetOriginal, const int targetLength,
				unsigned char** const queryTransformed,
				unsigned char** const targetTransformed,
                                int *alphabet_size) {
    // Alphabet is constructed from letters that are present in sequences.
    // Each letter is assigned an ordinal number, starting from 0 up to alphabetLength - 1,
    // and new query and target are created in which letters are replaced with their ordinal numbers.
    // This query and target are used in all the calculations later.
    *queryTransformed  = malloc(sizeof(unsigned char) * queryLength);
    *targetTransformed = malloc(sizeof(unsigned char) * targetLength);

    char *alphabet = malloc(MAX_UCHAR+1), *alphabet_cp = alphabet;

    // Alphabet information, it is constructed on fly while transforming sequences.
    // letterIdx[c] is index of letter c in alphabet.
    unsigned char letterIdx[MAX_UCHAR + 1];
    bool inAlphabet[MAX_UCHAR + 1]; // inAlphabet[c] is true if c is in alphabet
    for (int i = 0; i < MAX_UCHAR + 1; i++) inAlphabet[i] = false;

    for (int i = 0; i < queryLength; i++) {
        unsigned char c = queryOriginal[i];
        if (!inAlphabet[c]) {
            inAlphabet[c] = true;
            letterIdx[c] = alphabet_cp - alphabet;
            *alphabet_cp++ = queryOriginal[i];
        }
        (*queryTransformed)[i] = letterIdx[c];
    }
    for (int i = 0; i < targetLength; i++) {
        unsigned char c = targetOriginal[i];
        if (!inAlphabet[c]) {
            inAlphabet[c] = true;
            letterIdx[c] = alphabet_cp - alphabet;
            *alphabet_cp++ = targetOriginal[i];
        }
        (*targetTransformed)[i] = letterIdx[c];
    }

    *alphabet_size = alphabet_cp - alphabet;
    return alphabet;
}


EdlibAlignConfig edlibNewAlignConfig(int k, EdlibAlignMode mode, EdlibAlignTask task,
				     const EdlibEqualityPair* additionalEqualities,
				     int additionalEqualitiesLength) {
    EdlibAlignConfig config;
    config.k = k;
    config.mode = mode;
    config.task = task;
    config.additionalEqualities = additionalEqualities;
    config.additionalEqualitiesLength = additionalEqualitiesLength;
    return config;
}

EdlibAlignConfig edlibDefaultAlignConfig(void) {
    return edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);
}

void edlibFreeAlignResult(EdlibAlignResult result) {
    if (result.endLocations) free(result.endLocations);
    if (result.startLocations) free(result.startLocations);
    if (result.alignment) free(result.alignment);
}
