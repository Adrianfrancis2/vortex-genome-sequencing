#!/bin/bash
# Regenerates data/pairs.txt and data/golden.txt, then checks the banded
# recurrence against KSW2 at the band widths M5 sweeps.
set -e
cd "$(dirname "$0")"

make -s
python3 gen_pairs.py --out data/pairs.txt
./ksw2_ref -o data/golden.txt data/pairs.txt
for w in 16 32 64 128; do
	./ksw2_ref -w "$w" data/pairs.txt
done
