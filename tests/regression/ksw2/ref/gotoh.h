#ifndef GOTOH_H
#define GOTOH_H

#include <stdint.h>

typedef struct {
	int32_t score;          // global score: H at (tlen-1, qlen-1)
	int32_t max;            // best extension score, 0 if no cell is positive
	int32_t max_q, max_t;   // its cell, -1 if max == 0
	int32_t n_max;          // number of cells scoring max (tie count)
} gotoh_res_t;

// Affine-gap extension from (0,0) with KSW2's ksw_extz recurrence and
// boundaries; a gap of length l costs q + l*e. w < 0 disables the band.
// Ties for max go to the smallest anti-diagonal t+q, then the smallest t.
void gotoh_extz(int qlen, const uint8_t *query, int tlen, const uint8_t *target,
                int8_t match, int8_t mismatch, int8_t q, int8_t e, int w,
                gotoh_res_t *res);

#endif
