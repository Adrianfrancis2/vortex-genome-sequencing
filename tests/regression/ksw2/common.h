#ifndef _COMMON_H_
#define _COMMON_H_

#include <stdint.h>

#define KSW2_NEG_INF (-0x40000000)

#define KSW2_MODE_DIAG   0 // one alignment per warp, anti-diagonal wavefront
#define KSW2_MODE_THREAD 1 // one alignment per thread, row by row

// Sequences are bytes 0..3 (ACGT). Scratch offsets count int32 words.
typedef struct {
  uint32_t q_off, t_off;
  uint32_t qlen, tlen;
  uint32_t scratch_off;
} pair_t;

typedef struct {
  int32_t score;        // global score, H at (tlen-1, qlen-1)
  int32_t max;          // best extension score, 0 if no cell is positive
  int32_t max_q, max_t; // its cell, -1 when max is 0
} result_t;

typedef struct {
  uint64_t pairs_addr;
  uint64_t seq_addr;
  uint64_t scratch_addr;
  uint64_t result_addr;
  uint32_t num_pairs;
  uint32_t mode;
  int32_t  w;           // band width, < 0 for none
  int32_t  match, mismatch, gapo, gape;
} kernel_arg_t;

// Scratch words per pair. Diag: H for anti-diagonals r-2, r-1, r and E, F for
// r-1, r (7 vectors of tlen), plus (h, t, q) per lane for the max reduction.
// Thread: H and E for the previous and current row (4 vectors of qlen + 1).
static inline uint32_t ksw2_scratch_words(uint32_t qlen, uint32_t tlen, uint32_t num_lanes) {
  uint32_t diag = 7 * tlen + 3 * num_lanes;
  uint32_t row  = 4 * (qlen + 1);
  return diag > row ? diag : row;
}

#endif
