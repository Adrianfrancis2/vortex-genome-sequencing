# Genome alignment checker

`align_check.py` compiles a global-alignment program, runs it on the test
cases, and reports for every case whether the score is right and how fast the
alignment ran (GCUPS), followed by the total number of cases passed.

```bash
python3 tests/genome/align_check.py program.c
```

KSW2 (`git submodule update --init third_party/ksw2` first):

```bash
python3 tests/genome/align_check.py tests/regression/ksw2/ref/ksw2_score.c \
    third_party/ksw2/ksw2_extz2_sse.c --cflags "-O3 -march=native -Ithird_party/ksw2"
```

KSW2 has no score for an empty sequence (it returns `KSW_NEG_INF`), so the
cases with an empty query or target fail for it.

Every `.c` file given is compiled together with `gcc` (`--cc` to change it)
and `--cflags` (default `-O3 -march=native`).

## What the program must do

**Input**: the checker runs `<program> <pairs file>`. Each line of the file is
one pair of sequences, separated by a space:

```
ACGTACGT ACGTACGT
- ACG
ACNT ACGT
```

- Sequences use the uppercase letters `A C G T N`; `-` stands for an empty sequence.
- Lengths go up to about 300 bases.

**Scoring**: fixed, built into the program:

| match | mismatch | pair involving N | gap open | gap extend |
|---|---|---|---|---|
| +2 | −4 | −1 | 4 | 2 |

A gap of length k costs `4 + 2k`. The alignment is global (end to end).

**Output**: one line per pair, **in the same order as the input**:

```
16 0.000000412
-10 0.000000045
5 0.000000088
```

- The first number is the alignment score.
- The second is the time spent aligning that one pair, in seconds, measured
  inside the program around the alignment call only (for example with
  `clock_gettime(CLOCK_MONOTONIC, ...)`), not reading or printing. The checker
  divides the pair's cells (`len(query) * len(target)`) by it to get GCUPS. If
  it is left out, the score is still checked and GCUPS shows `n/a`.
- Print nothing else on standard output (use standard error for debugging),
  and exit with status 0.

## Output

One row per case (`--failures-only` prints only the failures):

```
result id                              expected       got    cells     GCUPS
PASS   c07_identical                         16        16       64    0.1174
```

then each suite's pass count and overall GCUPS (all cells divided by all
alignment time), and finally the total number of cases passed. The exit status
is 0 only if every case passes.

GCUPS for a tiny pair is mostly fixed per-call overhead, not throughput. The
overall GCUPS of the `random` suite is the meaningful number.

## Suites

Each suite is one run of the program on one pairs file.

| Suite | Cases | Where the expected score comes from |
|---|---|---|
| `hand` | 22 hand-picked cases in `cases/hand.tsv`, one rule of the problem each | Worked out by hand, confirmed by trying every alignment |
| `exhaustive` | Every pair of sequences up to `--max-len` (default 3): 7,225 pairs | Trying every alignment |
| `random` | `--random-count` random pairs up to `--random-max-len` bases, with mutations and occasional `N` | The three-matrix dynamic program, which `--self-test` checks against the definition |

`--suite hand|exhaustive|random|all` picks which ones run; the default is all three.

`--self-test` checks the checker itself: the values in `cases/hand.tsv` and the
dynamic program must agree with trying every alignment.
