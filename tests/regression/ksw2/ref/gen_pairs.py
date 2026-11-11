#!/usr/bin/env python3
"""Seeded generator of query/target pairs for the KSW2 reference.

Writes one pair per line: id set qlen tlen query target
"""
import argparse
import random

BASES = "ACGT"


def rand_seq(rng, n):
    return "".join(rng.choice(BASES) for _ in range(n))


def geom_len(rng, mean):
    n = 1
    while rng.random() > 1.0 / mean:
        n += 1
    return n


def mutate(rng, seq, div):
    """Substitutions take 3/4 of the divergence, indels 1/4 (split evenly)."""
    p_sub, p_ins, p_del = 0.75 * div, 0.125 * div, 0.125 * div
    out = []
    i = 0
    while i < len(seq):
        r = rng.random()
        if r < p_sub:
            out.append(rng.choice([b for b in BASES if b != seq[i]]))
            i += 1
        elif r < p_sub + p_ins:
            out.append(rand_seq(rng, geom_len(rng, 2.0)))
        elif r < p_sub + p_ins + p_del:
            i += geom_len(rng, 2.0)
        else:
            out.append(seq[i])
            i += 1
    res = "".join(out)
    return res if res else seq[:1]


def edge_pairs(rng):
    base = rand_seq(rng, 64)
    pairs = [
        ("A", "A"),
        ("A", "C"),
        ("ACGTACGT", "ACGTACGT"),
        ("AAAAAAAA", "CCCCCCCC"),
        ("ACGT", "ACGTACGT"),
        ("ACGTACGT", "ACGT"),
        ("AAAA", "AAAAAAAA"),
        ("AAAAAAAA", "AAAA"),
        ("ACAC", "ACACACAC"),
        ("ACACACAC", "ACAC"),
        (base[:20] + base[40:], base),
        (base, base[:20] + base[40:]),
        (base[:32] + rand_seq(rng, 12) + base[32:], base),
    ]
    for _ in range(32):
        pairs.append((rand_seq(rng, rng.randint(1, 8)), rand_seq(rng, rng.randint(1, 8))))
    return pairs


def synth_pairs(rng, lengths, divs, per_cell):
    for n in lengths:
        for d in divs:
            for _ in range(per_cell):
                target = rand_seq(rng, n)
                yield mutate(rng, target, d), target, f"synth_L{n}_d{int(d * 100)}"


def phix_pairs(rng, ref, read_lens, per_len, div, pad):
    for n in read_lens:
        for _ in range(per_len):
            pos = rng.randint(0, len(ref) - n - pad)
            read = mutate(rng, ref[pos:pos + n], div)
            yield read, ref[pos:pos + n + pad], f"phix_L{n}"


def read_fasta(path):
    with open(path) as f:
        return "".join(l.strip().upper() for l in f if not l.startswith(">"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=254)
    ap.add_argument("--phix", default="data/phix174.fa")
    ap.add_argument("--per-cell", type=int, default=16)
    ap.add_argument("-o", "--out", default="data/pairs.txt")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    rows = [(q, t, "edge") for q, t in edge_pairs(rng)]
    rows += list(synth_pairs(rng, [64, 128, 256, 512, 1024], [0.01, 0.05, 0.15], args.per_cell))
    rows += list(phix_pairs(rng, read_fasta(args.phix), [150, 250], args.per_cell, 0.01, 16))

    with open(args.out, "w") as f:
        for i, (q, t, s) in enumerate(rows):
            f.write(f"{i} {s} {len(q)} {len(t)} {q} {t}\n")


if __name__ == "__main__":
    main()
