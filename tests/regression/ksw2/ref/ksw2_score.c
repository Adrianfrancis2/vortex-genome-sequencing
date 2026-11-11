// Global alignment score with ksw_extz2_sse, in the format align_check.py
// expects: reads "query target" pairs ("-" is empty) from the file named on
// the command line and prints "score seconds" per pair, in order.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ksw2.h"

#define MATCH    2
#define MISMATCH 4
#define AMBIG    1
#define GAPO     4
#define GAPE     2

static int encode(const char *s, uint8_t *out)
{
	int n = 0;
	for (; s[n]; ++n) {
		switch (s[n]) {
		case 'A': out[n] = 0; break;
		case 'C': out[n] = 1; break;
		case 'G': out[n] = 2; break;
		case 'T': out[n] = 3; break;
		case 'N': out[n] = 4; break;
		default: return -1;
		}
	}
	return n;
}

static double now(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return ts.tv_sec + ts.tv_nsec * 1e-9;
}

int main(int argc, char *argv[])
{
	if (argc != 2) {
		fprintf(stderr, "Usage: ksw2_score pairs.txt\n");
		return 2;
	}
	FILE *in = fopen(argv[1], "r");
	if (!in) { perror(argv[1]); return 2; }

	int8_t mat[25];
	for (int i = 0; i < 5; ++i)
		for (int j = 0; j < 5; ++j)
			mat[i * 5 + j] = (i == 4 || j == 4) ? -AMBIG : (i == j ? MATCH : -MISMATCH);

	char *line = 0;
	size_t cap = 0;
	while (getline(&line, &cap, in) > 0) {
		char *qs = strtok(line, " \t\r\n");
		char *ts = strtok(0, " \t\r\n");
		if (!qs || !ts) continue;
		if (!strcmp(qs, "-")) qs = "";
		if (!strcmp(ts, "-")) ts = "";

		uint8_t *q = malloc(strlen(qs) + 1), *t = malloc(strlen(ts) + 1);
		int qlen = encode(qs, q), tlen = encode(ts, t);
		if (qlen < 0 || tlen < 0) {
			fprintf(stderr, "unexpected base in pair %s %s\n", qs, ts);
			return 2;
		}

		ksw_extz_t ez;
		memset(&ez, 0, sizeof(ez));
		double t0 = now();
		ksw_extz2_sse(0, qlen, q, tlen, t, 5, mat, GAPO, GAPE, -1, -1, 0, KSW_EZ_SCORE_ONLY, &ez);
		double t1 = now();

		printf("%d %.9f\n", ez.score, t1 - t0);
		free(ez.cigar);
		free(q);
		free(t);
	}
	free(line);
	fclose(in);
	return 0;
}
