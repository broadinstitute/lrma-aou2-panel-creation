#!/usr/bin/env python3
"""Synthetic bubble posteriors for test/cohort-info/run.sh. Writes into OUTDIR:

  main.vcf             GLIMPSE2-like ligated posteriors (one record per bubble path, GT:GP, INFO with RAF/AF/INFO)
  main.batch-NNN.vcf   the same, one file per batch of samples
  sites.vcf            the panel's bubble sites, ID= lists the atomic variants on each path
  ids.vcf              the panel's atomic (biallelic ID) variants, ID= and AF=

Allele frequencies run from absent to common so the AF 0, rare and multi-path cases all occur.
"""
import os
import random
import sys

N_SAMPLES, N_BATCHES, N_BUBBLES, CHROM_LEN = 48, 4, 3000, 2_000_000
FREQS = [0.0, 0.0, 0.002, 0.01, 0.03, 0.1, 0.25, 0.5]
BASES = 'ACGT'


def hap_probs(carried, rng):
    """GLIMPSE-like haplotype posterior: confident most of the time, sometimes wrong-ish."""
    c = rng.choice([0.999, 0.99, 0.95, 0.8, 0.6])
    return c if carried else (1 - c) * rng.random()


def gp_of(p0, p1):
    g = [(1 - p0) * (1 - p1), p0 * (1 - p1) + (1 - p0) * p1, p0 * p1]
    v = [round(x * 1000) for x in g]
    v[v.index(max(v))] += 1000 - sum(v)
    return ','.join('%.3f' % (x / 1000) for x in v)


def main():
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    rng = random.Random(20261115)
    samples = ['S%03d' % (i + 1) for i in range(N_SAMPLES)]
    atomic, records = [], []      # atomic: (pos, ref, alt, aid, raf); records: (pos, ref, alt, [aid], [cells], info)
    used = set()
    for k in range(N_BUBBLES):
        pos = 1000 + 600 * k
        n_paths = rng.choice([1, 1, 1, 2, 2, 3, 4])
        n_atomic = n_paths + rng.choice([0, 0, 1])
        aids = []
        for j in range(n_atomic):
            while True:
                apos = pos + rng.randrange(0, 40)
                ref = rng.choice(BASES)
                alt = ''.join(rng.choice(BASES) for _ in range(60)) if rng.random() < 0.1 else rng.choice([b for b in BASES if b != ref])
                if (apos, ref, alt) not in used:
                    used.add((apos, ref, alt))
                    break
            aid = 'a%d_%d' % (k, j)
            atomic.append((apos, ref, alt, aid, '%.6f' % rng.random()))
            aids.append(aid)
        paths = [[aids[i]] for i in range(n_paths)]
        for extra in aids[n_paths:]:
            rng.choice(paths).append(extra)
        if n_paths > 1 and rng.random() < 0.3:      # overlapping paths: a variant on several paths
            paths[1].append(paths[0][0])
        freqs = [rng.choice(FREQS) for _ in paths]
        # each haplotype carries at most one path
        haps = []
        for _ in range(N_SAMPLES):
            h = []
            for _ in range(2):
                r, which = rng.random(), None
                for i, f in enumerate(freqs):
                    if r < f:
                        which = i
                        break
                    r -= f
                h.append(which)
            haps.append(h)
        for i, path in enumerate(paths):
            cells = []
            for h in haps:
                p0, p1 = hap_probs(h[0] == i, rng), hap_probs(h[1] == i, rng)
                cells.append('%d|%d:%s' % (p0 > 0.5, p1 > 0.5, gp_of(p0, p1)))
            info = 'RAF=%.6f;AF=%.6f;INFO=%.3f' % (rng.random(), freqs[i], rng.random())
            records.append((pos, 'A', 'C' * (i + 1), path, cells, info))

    contig = '##contig=<ID=chr1,length=%d>' % CHROM_LEN
    head = ['##fileformat=VCFv4.2', contig,
            '##INFO=<ID=RAF,Number=A,Type=Float,Description="panel AF">',
            '##INFO=<ID=AF,Number=A,Type=Float,Description="AF">',
            '##INFO=<ID=INFO,Number=A,Type=Float,Description="INFO">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">',
            '##FORMAT=<ID=DS,Number=A,Type=Float,Description="dosage">',
            '##FORMAT=<ID=GP,Number=G,Type=Float,Description="genotype probs">']
    per = N_SAMPLES // N_BATCHES
    groups = [('main', range(N_SAMPLES))] + [('main.batch-%03d' % (b + 1), range(b * per, (b + 1) * per)) for b in range(N_BATCHES)]
    for name, cols in groups:
        with open(os.path.join(out, name + '.vcf'), 'w') as f:
            f.write('\n'.join(head) + '\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples[c] for c in cols) + '\n')
            for pos, ref, alt, _, cells, info in records:
                f.write('\t'.join(['chr1', str(pos), '.', ref, alt, '.', '.', info, 'GT:GP'] + [cells[c] for c in cols]) + '\n')
    with open(os.path.join(out, 'sites.vcf'), 'w') as f:
        f.write('##fileformat=VCFv4.2\n%s\n##INFO=<ID=ID,Number=.,Type=String,Description="atomic IDs">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' % contig)
        for pos, ref, alt, path, _, _ in records:
            f.write('\t'.join(['chr1', str(pos), '.', ref, alt, '.', '.', 'ID=' + ':'.join(path)]) + '\n')
    with open(os.path.join(out, 'ids.vcf'), 'w') as f:
        f.write('##fileformat=VCFv4.2\n%s\n##INFO=<ID=ID,Number=1,Type=String,Description="atomic ID">\n##INFO=<ID=AF,Number=A,Type=Float,Description="panel AF">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' % contig)
        for apos, ref, alt, aid, raf in sorted(atomic):
            f.write('\t'.join(['chr1', str(apos), 'orig_' + aid, ref, alt, '.', '.', 'ID=%s;AF=%s' % (aid, raf)]) + '\n')


if __name__ == '__main__':
    main()
