// Copyright Eugene Shelwien.
// Release into public domain.

// Modifications by James Bonfield (2019)


/*
 * Note it is up to the calling code to ensure that no overruns on input and
 * output buffers occur.
 *
 * Call the input() and output() functions to set and query the current
 * buffer locations.
 *

 */

#ifndef C_RANGER_CODER_H
#define C_RANGER_CODER_H

#define  DO(n)     int _;for (_=0; _<n; _++)
#define  TOP       (1<<24)
#define  Thres (unsigned)255*TOP

typedef unsigned char uc;

typedef struct {
    uint32_t low, code, range;
    uint32_t FFNum;  // Number of consecutive FFs
    uint32_t Cache;  // Top 8-bits of low ready to emit
    uint32_t Carry;  // Flag to indicate if we emit Cache or Cache+1
    uc *in_buf;
    uc *out_buf;
    uc *in_end;
} RangeCoder;

static inline void RC_SetInput(RangeCoder *rc, char *in, char *in_end) {
    rc->out_buf = rc->in_buf = (uc *)in;
    rc->in_end = (uc *)in_end;
}
static inline void RC_SetOutput(RangeCoder *rc, char *out) { rc->in_buf = rc->out_buf = (uc *)out; }
static inline char *RC_GetInput(RangeCoder *rc) { return (char *)rc->in_buf; }
static inline char *RC_GetOutput(RangeCoder *rc) { return (char *)rc->out_buf; }
static inline size_t RC_OutSize(RangeCoder *rc) { return rc->out_buf - rc->in_buf; }
static inline size_t RC_InSize(RangeCoder *rc) { return rc->in_buf - rc->out_buf; }

static inline void RC_StartEncode(RangeCoder *rc)
{ 
    rc->range = 0xFFFFFFFF;
    rc->low   = 0;
    rc->FFNum = 0;
    rc->Carry = 0;
    rc->Cache = 0;
    rc->code  = 0;
}

static inline void RC_StartDecode(RangeCoder *rc)
{ 
    rc->range = 0xFFFFFFFF;
    rc->low   = 0;
    rc->FFNum = 0;
    rc->Carry = 0;
    rc->Cache = 0;
    rc->code  = 0;
    if (rc->in_buf+5 > rc->in_end) {
        rc->in_buf = rc->in_end; // prevent decode
        return;
    }
    DO(5) rc->code = (rc->code<<8) | *rc->in_buf++;
}

static inline void RC_ShiftLow(RangeCoder *rc) {
    if (rc->low < Thres || rc->Carry) {
        *rc->out_buf++ = rc->Cache + rc->Carry;

        // Flush any stored FFs
        while (rc->FFNum) {
            *rc->out_buf++ = rc->Carry-1; // (Carry-1)&255;
            rc->FFNum--;
        }

        // Take copy of top byte ready for next flush
        rc->Cache = rc->low >> 24;
        rc->Carry = 0;
    } else {
        // Low if FFxx xxxx.  Bump FF count and shift in as before
        rc->FFNum++;
    }
    rc->low = rc->low<<8;
}

static inline void RC_FinishEncode(RangeCoder *rc) 
{ 
    DO(5) RC_ShiftLow(rc);
}

static inline void RC_FinishDecode(RangeCoder *rc) {}

static inline void RC_Encode (RangeCoder *rc, uint32_t cumFreq, uint32_t freq, uint32_t totFreq) 
{
    uint32_t tmp = rc->low;
    rc->low  += cumFreq * (rc->range/= totFreq);
    rc->range*= freq;

    rc->Carry += rc->low<tmp; // Overflow

    while (rc->range < TOP) {
        rc->range <<= 8;
        RC_ShiftLow(rc);
    }
}

static inline uint32_t RC_GetFreq (RangeCoder *rc, uint32_t totFreq) {
    //return rc->code/(rc->range/=totFreq);
    return (totFreq && rc->range >= totFreq) ? rc->code/(rc->range/=totFreq) : 0;
}

static inline void RC_Decode (RangeCoder *rc, uint32_t cumFreq, uint32_t freq, uint32_t totFreq) 
{
    rc->code -= cumFreq * rc->range;
    rc->range *= freq;
    while (rc->range < TOP) {
        if (rc->in_buf >= rc->in_end)
            return; // FIXME: could signal error, instead of caller just generating nonsense
        rc->code = (rc->code<<8) + *rc->in_buf++;
        rc->range <<= 8;
    }
}

#endif /* C_RANGER_CODER_H */
