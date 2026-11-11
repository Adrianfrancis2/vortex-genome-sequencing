// Runs KSW2 (ksw_extz2_sse and the scalar ksw_extz) and our gotoh_extz on
// every pair, checks that they agree, and writes the golden file:
//   id gscore max max_q max_t
// Positions come from gotoh_extz, whose tie rule is fixed; KSW2's positions
// are compared only where the max cell is unique, since its two versions
// break ties in different orders. Pairs with |tlen - qlen| > w are skipped
// under a band.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "ksw2.h"
#include "gotoh.h"

#define MATCH    2
#define MISMATCH 4
#define GAPO     4
#define GAPE     2

static int encode(char *s, int n)
{
	for (int i = 0; i < n; ++i) {
		switch (s[i]) {
		case 'A': s[i] = 0; break;
		case 'C': s[i] = 1; break;
		case 'G': s[i] = 2; break;
		case 'T': s[i] = 3; break;
		default: return -1;
		}
	}
	return 0;
}

int main(int argc, char *argv[])
{
	const char *out_path = 0;
	int w = -1, c;
	while ((c = getopt(argc, argv, "o:w:")) >= 0) {
		if (c == 'o') out_path = optarg;
		else if (c == 'w') w = atoi(optarg);
	}
	if (optind >= argc) {
		fprintf(stderr, "Usage: ksw2_ref [-w band] [-o golden.txt] pairs.txt\n");
		return 2;
	}
	FILE *in = fopen(argv[optind], "r");
	if (!in) { perror(argv[optind]); return 2; }
	FILE *out = out_path ? fopen(out_path, "w") : 0;
	if (out_path && !out) { perror(out_path); return 2; }

	// 4 bases plus KSW2's wildcard row/column, which scores 0.
	int8_t mat[25];
	for (int i = 0; i < 5; ++i)
		for (int j = 0; j < 5; ++j)
			mat[i * 5 + j] = (i == 4 || j == 4) ? 0 : (i == j ? MATCH : -MISMATCH);

	long n_skipped = 0, n_sse_band = 0;
	long n_pairs = 0, n_sse_scalar = 0, n_gotoh = 0, n_pos = 0, n_ties = 0, n_cells = 0;
	char *line = 0;
	size_t cap = 0;
	while (getline(&line, &cap, in) > 0) {
		char set[64];
		int id, qlen, tlen, off_q, off_t;
		if (sscanf(line, "%d %63s %d %d %n%*s %n", &id, set, &qlen, &tlen, &off_q, &off_t) < 4) {
			fprintf(stderr, "bad line: %s", line);
			return 2;
		}
		char *qs = line + off_q, *ts = line + off_t;
		ts[tlen] = 0;
		qs[qlen] = 0;
		if (encode(qs, qlen) || encode(ts, tlen)) {
			fprintf(stderr, "pair %d: non-ACGT base\n", id);
			return 2;
		}

		// ksw_extz writes past its row buffer when a whole target row lies
		// outside the band, which happens exactly when tlen - qlen > w.
		if (w >= 0 && abs(tlen - qlen) > w) {
			++n_skipped;
			if (out)
				fprintf(out, "%d skipped\n", id);
			continue;
		}

		ksw_extz_t ez_sse, ez_sc;
		memset(&ez_sse, 0, sizeof(ez_sse));
		memset(&ez_sc, 0, sizeof(ez_sc));
		ksw_extz2_sse(0, qlen, (uint8_t *)qs, tlen, (uint8_t *)ts, 5, mat, GAPO, GAPE, w, -1, 0, KSW_EZ_SCORE_ONLY, &ez_sse);
		ksw_extz(0, qlen, (uint8_t *)qs, tlen, (uint8_t *)ts, 5, mat, GAPO, GAPE, w, -1, KSW_EZ_SCORE_ONLY, &ez_sc);
		free(ez_sse.cigar);
		free(ez_sc.cigar);

		gotoh_res_t g;
		gotoh_extz(qlen, (uint8_t *)qs, tlen, (uint8_t *)ts, MATCH, MISMATCH, GAPO, GAPE, w, &g);

		++n_pairs;
		n_cells += (long)qlen * tlen;
		// With a band, SSE seeds uncomputed neighbors with a zero difference
		// instead of -inf, so paths just outside the band can leak in; the
		// scalar version's strict |t - q| <= w band is the reference there.
		int sse_differs = ez_sse.score != ez_sc.score || (int)ez_sse.max != (int)ez_sc.max;
		if (sse_differs && w < 0) {
			++n_sse_scalar;
			fprintf(stderr, "pair %d (%s): SSE score/max %d/%d vs scalar %d/%d\n",
			        id, set, ez_sse.score, (int)ez_sse.max, ez_sc.score, (int)ez_sc.max);
		} else if (sse_differs) {
			++n_sse_band;
		}
		if (g.score != ez_sc.score || g.max != (int)ez_sc.max) {
			++n_gotoh;
			fprintf(stderr, "pair %d (%s): gotoh score/max %d/%d vs KSW2 %d/%d\n",
			        id, set, g.score, g.max, ez_sc.score, (int)ez_sc.max);
		}
		if (g.n_max > 1) {
			++n_ties;
		} else if (g.max_q != ez_sc.max_q || g.max_t != ez_sc.max_t
		           || (!sse_differs && (g.max_q != ez_sse.max_q || g.max_t != ez_sse.max_t))) {
			++n_pos;
			fprintf(stderr, "pair %d (%s): unique max at q=%d t=%d, scalar %d,%d, SSE %d,%d\n",
			        id, set, g.max_q, g.max_t, ez_sc.max_q, ez_sc.max_t, ez_sse.max_q, ez_sse.max_t);
		}
		if (out)
			fprintf(out, "%d %d %d %d %d\n", id, g.score, g.max, g.max_q, g.max_t);
	}
	free(line);
	fclose(in);
	if (out) fclose(out);

	long fails = n_sse_scalar + n_gotoh + n_pos;
	fprintf(stderr, "w=%d pairs=%ld (skipped %ld) cells=%ld | SSE!=scalar %ld (band leak, info only: %ld) | "
	        "gotoh!=KSW2 %ld | unique-max position mismatches %ld | tied-max pairs %ld | %s\n",
	        w, n_pairs, n_skipped, n_cells, n_sse_scalar, n_sse_band, n_gotoh, n_pos, n_ties,
	        fails ? "FAIL" : "PASS");
	return fails ? 1 : 0;
}
