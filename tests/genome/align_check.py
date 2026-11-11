#!/usr/bin/env python3
"""Compile a global-alignment program, run it on the test cases, and check every score.

    python3 align_check.py program.c [more.c ...] [--cflags "-Ipath"]

The program is run as `program PAIRS_FILE`, where each line of the file is
"QUERY TARGET", and prints one "<score> <seconds>" line per pair in the same
order; <seconds> is the time spent on that alignment itself (see README.md).
"""

import argparse
import itertools
import os
import random
import shlex
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
HAND_CASES = os.path.join(HERE, "cases", "hand.tsv")
BASES = "ACGT"
NEG = float("-inf")

# The scoring the program under test must use; cases/hand.tsv is computed with it.
MATCH, MISMATCH, AMBIG, GAP_OPEN, GAP_EXT = 2, 4, 1, 4, 2


class Case:
    def __init__(self, cid, query, target, expected):
        self.id = cid
        self.query = query
        self.target = target
        self.expected = expected

    @property
    def cells(self):
        return len(self.query) * len(self.target)


def pair_score(x, y):
    if x == "N" or y == "N":
        return -AMBIG
    return MATCH if x == y else -MISMATCH


def gotoh(query, target):
    """Score by the three-matrix dynamic program."""
    n, m = len(query), len(target)
    H = [[NEG] * (m + 1) for _ in range(n + 1)]
    E = [[NEG] * (m + 1) for _ in range(n + 1)]
    F = [[NEG] * (m + 1) for _ in range(n + 1)]
    H[0][0] = 0
    for j in range(1, m + 1):
        E[0][j] = H[0][j] = -(GAP_OPEN + j * GAP_EXT)
    for i in range(1, n + 1):
        F[i][0] = H[i][0] = -(GAP_OPEN + i * GAP_EXT)
    for i in range(1, n + 1):
        qi = query[i - 1]
        for j in range(1, m + 1):
            E[i][j] = max(E[i][j - 1] - GAP_EXT, H[i][j - 1] - GAP_OPEN - GAP_EXT)
            F[i][j] = max(F[i - 1][j] - GAP_EXT, H[i - 1][j] - GAP_OPEN - GAP_EXT)
            H[i][j] = max(H[i - 1][j - 1] + pair_score(qi, target[j - 1]), E[i][j], F[i][j])
    return int(H[n][m])


_alignment_cache = {}


def alignments(n, m):
    """Every alignment of lengths n and m, as strings of M (pair), I (query base vs gap), D (gap vs target base)."""
    key = (n, m)
    if key not in _alignment_cache:
        out = []

        def walk(i, j, ops):
            if i == n and j == m:
                out.append("".join(ops))
                return
            if i < n and j < m:
                walk(i + 1, j + 1, ops + ["M"])
            if i < n:
                walk(i + 1, j, ops + ["I"])
            if j < m:
                walk(i, j + 1, ops + ["D"])

        walk(0, 0, [])
        _alignment_cache[key] = out
    return _alignment_cache[key]


def by_definition(query, target):
    """Score by trying every alignment literally; only feasible for short sequences."""
    best = NEG
    for ops in alignments(len(query), len(target)):
        i = j = 0
        total = 0
        prev = ""
        for op in ops:
            if op == "M":
                total += pair_score(query[i], target[j])
                i += 1
                j += 1
            else:
                total -= GAP_EXT if op == prev else GAP_OPEN + GAP_EXT
                if op == "I":
                    i += 1
                else:
                    j += 1
            prev = op
        best = max(best, total)
    return int(best)


def load_hand_cases():
    cases = []
    with open(HAND_CASES) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) != 4:
                sys.exit(f"{HAND_CASES}: expected 4 columns, got {len(cols)}: {line}")
            cid, q, t, exp = cols
            cases.append(Case(cid, "" if q == "-" else q, "" if t == "-" else t, int(exp)))
    return cases


def exhaustive_cases(max_len):
    seqs = [""]
    for length in range(1, max_len + 1):
        seqs += ["".join(p) for p in itertools.product(BASES, repeat=length)]
    return [Case(f"x_{q or '-'}_{t or '-'}", q, t, by_definition(q, t))
            for q, t in itertools.product(seqs, seqs)]


def mutate(seq, rng, rate):
    out = []
    for base in seq:
        r = rng.random()
        if r < rate / 3:
            out.append(rng.choice([b for b in BASES if b != base]))
        elif r < 2 * rate / 3:
            continue
        elif r < rate:
            out.append(base)
            out.extend(rng.choice(BASES) for _ in range(rng.randint(1, 6)))
        else:
            out.append(base)
    return "".join(out)


def random_cases(count, max_len, seed):
    rng = random.Random(seed)
    cases = []
    for k in range(count):
        length = rng.randint(1, max_len)
        q = "".join(rng.choice(BASES) for _ in range(length))
        if rng.random() < 0.1:
            t = "".join(rng.choice(BASES) for _ in range(rng.randint(1, max_len)))
        else:
            t = mutate(q, rng, rng.choice([0.01, 0.05, 0.1, 0.2]))
        if rng.random() < 0.05:
            pos = rng.randrange(len(q))
            q = q[:pos] + "N" + q[pos + 1:]
        cases.append(Case(f"r{k:05d}", q, t, gotoh(q, t)))
    return cases


def compile_program(sources, cc, cflags, workdir):
    for src in sources:
        if not os.path.isfile(src):
            sys.exit(f"no such file: {src}")
    exe = os.path.join(workdir, "aligner")
    cmd = [cc, *shlex.split(cflags), "-o", exe, *sources, "-lm"]
    print("compiling:", " ".join(shlex.quote(c) for c in cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stdout + proc.stderr)
        sys.exit("compilation failed")
    return exe


def run_program(exe, cases, timeout, workdir):
    """Run the program once on all pairs; return one (score, seconds, error) per case, in order."""
    path = os.path.join(workdir, "pairs.txt")
    with open(path, "w") as f:
        for c in cases:
            f.write(f"{c.query or '-'} {c.target or '-'}\n")
    try:
        proc = subprocess.run([exe, path], capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        sys.exit(f"program timed out after {timeout}s")
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        sys.exit(f"program exited with status {proc.returncode}")
    lines = [line for line in proc.stdout.splitlines() if line.strip()]
    if len(lines) != len(cases):
        print(f"warning: program printed {len(lines)} lines for {len(cases)} pairs; "
              f"answers are matched to pairs by line position")
    results = []
    for k in range(len(cases)):
        if k >= len(lines):
            results.append((None, None, "no output line"))
            continue
        cols = lines[k].split()
        try:
            score = int(cols[0])
            seconds = float(cols[1]) if len(cols) > 1 else None
            results.append((score, seconds, None))
        except (IndexError, ValueError):
            results.append((None, None, f"unparseable output {lines[k][:200]!r}"))
    return results


def gcups(cells, seconds):
    if not cells or not seconds or seconds <= 0:
        return None
    return cells / seconds / 1e9


def fmt_gcups(value):
    return "n/a" if value is None else f"{value:.4f}"


def run_suite(name, exe, cases, timeout, failures_only, workdir):
    print(f"\n== {name}: {len(cases)} cases ==")
    results = run_program(exe, cases, timeout, workdir)
    print(f"{'result':6s} {'id':30s} {'expected':>9s} {'got':>9s} {'cells':>8s} {'GCUPS':>9s}")
    passed = 0
    timed_cells = 0
    timed_seconds = 0.0
    for c, (score, seconds, error) in zip(cases, results):
        ok = score == c.expected
        passed += ok
        if seconds is not None and c.cells:
            timed_cells += c.cells
            timed_seconds += seconds
        if failures_only and ok:
            continue
        got = "error" if score is None else str(score)
        print(f"{'PASS' if ok else 'FAIL':6s} {c.id:30s} {c.expected:9d} {got:>9s} {c.cells:8d} "
              f"{fmt_gcups(gcups(c.cells, seconds)):>9s}" + (f"  ({error})" if error else ""))
    print(f"{name}: {passed}/{len(cases)} passed, "
          f"overall GCUPS {fmt_gcups(gcups(timed_cells, timed_seconds))} "
          f"({timed_cells} cells in {timed_seconds:.6f}s)")
    return passed, timed_cells, timed_seconds


def self_test(max_len):
    """Check that the hand-computed values and the dynamic program both agree with the definition."""
    ok = True
    for c in load_hand_cases():
        d = by_definition(c.query, c.target) if len(c.query) + len(c.target) <= 20 else None
        g = gotoh(c.query, c.target)
        if g != c.expected or (d is not None and d != c.expected):
            print(f"hand case {c.id}: file={c.expected} dp={g} definition={d}")
            ok = False
    count = 0
    for c in exhaustive_cases(max_len):
        count += 1
        if gotoh(c.query, c.target) != c.expected:
            print(f"dp disagrees with definition on {c.id}")
            ok = False
    print(f"self-test: hand cases and {count} exhaustive pairs {'agree' if ok else 'DISAGREE'}")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sources", nargs="*", help="the program's .c file, plus any other .c files it needs")
    ap.add_argument("--cc", default="gcc")
    ap.add_argument("--cflags", default="-O3 -march=native",
                    help='compiler flags, e.g. "-O3 -march=native -I../../third_party/ksw2"')
    ap.add_argument("--suite", choices=["hand", "exhaustive", "random", "all"], default="all")
    ap.add_argument("--max-len", type=int, default=3,
                    help="exhaustive suite: every pair up to this length")
    ap.add_argument("--random-count", type=int, default=200)
    ap.add_argument("--random-max-len", type=int, default=150)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--timeout", type=float, default=600, help="seconds allowed per program run")
    ap.add_argument("--failures-only", action="store_true", help="print only the failing cases")
    ap.add_argument("--self-test", action="store_true", help="check this script's own expected values and exit")
    args = ap.parse_args()

    if args.self_test:
        sys.exit(0 if self_test(args.max_len) else 1)
    if not args.sources:
        ap.error("give the program's .c file")

    suites = ["hand", "exhaustive", "random"] if args.suite == "all" else [args.suite]
    workdir = tempfile.mkdtemp(prefix="align_check_")
    try:
        exe = compile_program(args.sources, args.cc, args.cflags, workdir)
        totals = []
        for name in suites:
            if name == "hand":
                cases = load_hand_cases()
            elif name == "exhaustive":
                cases = exhaustive_cases(args.max_len)
            else:
                cases = random_cases(args.random_count, args.random_max_len, args.seed)
            totals.append((name, len(cases),
                           *run_suite(name, exe, cases, args.timeout, args.failures_only, workdir)))
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    print("\n== summary ==")
    for name, count, passed, cells, seconds in totals:
        print(f"{name:12s} {passed:6d}/{count:<6d} passed   GCUPS {fmt_gcups(gcups(cells, seconds))}")
    total = sum(t[1] for t in totals)
    passed = sum(t[2] for t in totals)
    print(f"{'TOTAL':12s} {passed:6d}/{total:<6d} passed")
    sys.exit(0 if passed == total else 1)


if __name__ == "__main__":
    main()
