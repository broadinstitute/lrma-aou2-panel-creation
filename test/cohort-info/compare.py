#!/usr/bin/env python3
"""compare.py MERGED.tsv REFERENCE.tsv : CHROM POS ID REF ALT AF INFO per record, merged cohort vs one pop over all samples.

Pooling recovers each batch's sums from AF (6 decimals) and INFO (3 decimals), so pooled INFO is within 5e-4 of the exact
value, and both files print INFO to 3 decimals: |delta INFO| <= 1.5e-3. AF agrees to the 6th decimal. When the reference AF
prints as 0 or 1 its INFO is 1 by convention, while the pooled p of batches that each print a tiny nonzero AF need not be, so
INFO is not compared there.
"""
import sys

INFO_TOL, AF_TOL = 1.6e-3, 2e-6
m, r = (open(p).read().splitlines() for p in sys.argv[1:3])
if len(m) != len(r):
    sys.exit(f'FAIL record count {len(m)} vs {len(r)}')
bad, worst, n_cmp = 0, 0.0, 0
for a, b in zip(m, r):
    a, b = a.split('\t'), b.split('\t')
    if a[:5] != b[:5]:
        sys.exit(f'FAIL site {a[:5]} vs {b[:5]}')
    af_m, info_m, af_r, info_r = float(a[5]), float(a[6]), float(b[5]), float(b[6])
    if abs(af_m - af_r) > AF_TOL:
        bad += 1
        print('AF', a[:5], af_m, af_r)
    if 1e-6 <= af_r <= 1 - 1e-6:
        n_cmp += 1
        worst = max(worst, abs(info_m - info_r))
        if abs(info_m - info_r) > INFO_TOL:
            bad += 1
            print('INFO', a[:5], info_m, info_r)
print(f'{len(m)} records, INFO compared on {n_cmp}, max |delta INFO| {worst:.4f}, {bad} out of tolerance')
sys.exit(1 if bad else 0)
