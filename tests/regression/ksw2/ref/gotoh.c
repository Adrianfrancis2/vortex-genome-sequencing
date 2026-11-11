#include <stdlib.h>
#include "gotoh.h"

#define NEG_INF (-0x40000000)

static inline int32_t max2(int32_t a, int32_t b) { return a > b ? a : b; }

// Rows are target positions, columns query positions; index 0 is the
// boundary, so cell (I,J) is target[I-1] vs query[J-1].
void gotoh_extz(int qlen, const uint8_t *query, int tlen, const uint8_t *target,
                int8_t match, int8_t mismatch, int8_t q, int8_t e, int w,
                gotoh_res_t *res)
{
	const int32_t qe = q + e;
	if (w < 0) w = tlen > qlen ? tlen : qlen;

	int32_t *Hp = malloc((qlen + 1) * sizeof(int32_t));
	int32_t *Ep = malloc((qlen + 1) * sizeof(int32_t));
	int32_t *Hc = malloc((qlen + 1) * sizeof(int32_t));
	int32_t *Ec = malloc((qlen + 1) * sizeof(int32_t));

	// KSW2 opens the top boundary only up to column w and the left one
	// only up to row w+1; outside the band every value is -inf.
	Hp[0] = 0;
	Ep[0] = NEG_INF;
	for (int J = 1; J <= qlen; ++J) {
		Hp[J] = J <= w ? -(q + e * J) : NEG_INF;
		Ep[J] = NEG_INF;
	}

	int32_t max = 0, max_q = -1, max_t = -1, n_max = 0;

	for (int I = 1; I <= tlen; ++I) {
		Hc[0] = I - 1 <= w ? -(q + e * I) : NEG_INF;
		Ec[0] = NEG_INF;
		int32_t f = NEG_INF;
		for (int J = 1; J <= qlen; ++J) {
			if (abs(I - J) > w) {
				Hc[J] = Ec[J] = f = NEG_INF;
				continue;
			}
			int32_t s = target[I - 1] == query[J - 1] ? match : -mismatch;
			int32_t ecur = max2(Hp[J] - qe, Ep[J] - e);
			int32_t fcur = max2(Hc[J - 1] - qe, f - e);
			int32_t h = max2(Hp[J - 1] + s, max2(ecur, fcur));
			Hc[J] = h;
			Ec[J] = ecur;
			f = fcur;

			int32_t t = I - 1, qq = J - 1;
			if (h > max) {
				max = h, max_t = t, max_q = qq, n_max = 1;
			} else if (h == max && max > 0) {
				++n_max;
				int r = t + qq, best_r = max_t + max_q;
				if (r < best_r || (r == best_r && t < max_t))
					max_t = t, max_q = qq;
			}
		}
		int32_t *tmp;
		tmp = Hp, Hp = Hc, Hc = tmp;
		tmp = Ep, Ep = Ec, Ec = tmp;
	}

	res->score = abs(tlen - qlen) <= w ? Hp[qlen] : NEG_INF;
	res->max = max;
	res->max_q = max_q;
	res->max_t = max_t;
	res->n_max = n_max;

	free(Hp); free(Ep); free(Hc); free(Ec);
}
