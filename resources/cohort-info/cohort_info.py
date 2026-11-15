#!/usr/bin/env python3
"""Pool per-batch AF/INFO into one cohort INFO, then write it onto popped BCFs.

    cohort_info.py pool --n-samples 250,250,... --sites sites.tsv.gz --report report.txt --out cohort.tsv.gz \
        info.batch-001.tsv.gz info.batch-002.tsv.gz ...
    bcftools view popped.batch-NNN.bcf | cohort_info.py apply --table cohort.tsv.gz > popped.batch-NNN.cohort.vcf

paste-vcfs copies INFO from the first input of each paste. The hierarchical paste pastes groups of batches, so `apply` is
run on the first batch of every level-0 group; every level then carries the cohort values.

Each info.tsv.gz has one line per record of the popped file, in file order: AF and INFO from pop-glimpse2 run on that
batch's bubble posteriors, then the INFO of the old deliverable (used only for the report). sites.tsv.gz is CHROM, POS, ID,
REF, ALT, INFO/ID, RAF, N_PATHS, N_PATHS_TOTAL from the re-pop of the first batch; RAF and the path counts come from the
panel and are the same in every batch. The caller has checked that all batches have the same site list.

Per batch b with n_b samples, pop-glimpse2 writes AF_b = DS_b / 2n_b and INFO_b = 1 - VAR_b / (2n_b AF_b (1 - AF_b)), so
DS_b and VAR_b are recovered from the two printed values. The cohort values are

    p    = sum DS_b / sum 2n_b
    INFO = 1 - sum VAR_b / (sum 2n_b p (1 - p)),   1 when p is 0 or 1 (GLIMPSE2 convention).

INFO_b is printed to 3 decimals. Because sum VAR_b <= sum 2n_b AF_b(1-AF_b) <= sum 2n_b p(1-p), the pooled INFO is off by at
most 5e-4 from the value of the unrounded sums.
"""
import argparse
import gzip
import subprocess
import sys
from itertools import zip_longest

TAGS = [
    '##INFO=<ID=ID,Number=1,Type=String,Description="Atomic variant ID">',
    '##INFO=<ID=RAF,Number=A,Type=Float,Description="ALT allele frequency in the reference panel">',
    '##INFO=<ID=AF,Number=A,Type=Float,Description="ALT allele frequency in the cohort, DS_SUM / N_HAPS">',
    '##INFO=<ID=INFO,Number=A,Type=Float,Description="Cohort IMPUTE INFO, 1 - VAR_SUM / (N_HAPS AF (1 - AF)); 1 when AF is 0 or 1">',
    '##INFO=<ID=N_PATHS,Number=1,Type=Integer,Description="Number of bubble paths that contain this variant">',
    '##INFO=<ID=N_PATHS_TOTAL,Number=1,Type=Integer,Description="Number of total non-reference paths in the bubble that contains this variant">',
    '##INFO=<ID=DS_SUM,Number=1,Type=Float,Description="Sum over cohort samples of the expected ALT dosage">',
    '##INFO=<ID=VAR_SUM,Number=1,Type=Float,Description="Sum over cohort samples of the dosage variance">',
    '##INFO=<ID=N_HAPS,Number=1,Type=Integer,Description="Number of cohort haplotypes">',
]

# report bins for INFO and for the change in INFO
INFO_EDGES = [0.0, 0.1, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.99, 1.0]
DELTA_EDGES = [-1.0, -0.5, -0.1, -0.01, 0.01, 0.1, 0.5, 1.0]


def zcat(path):
    return subprocess.Popen(['gzip', '-dc', path], stdout=subprocess.PIPE, text=True, bufsize=1 << 20).stdout


def num(s):
    return None if s in ('.', '') else float(s)


def bin_of(x, edges):
    for i in range(1, len(edges) - 1):
        if x < edges[i]:
            return i - 1
    return len(edges) - 2


def pool(a):
    n = [int(x) for x in a.n_samples.split(',')]
    if len(n) != len(a.info_tsvs):
        sys.exit(f'ERROR: {len(n)} sample counts for {len(a.info_tsvs)} batches')
    n_haps = 2 * sum(n)
    sites = zcat(a.sites)
    infos = [zcat(p) for p in a.info_tsvs]
    out = gzip.open(a.out, 'wt', compresslevel=1)
    out.write(f'#cohort_info\tbatches={len(n)}\tsamples={sum(n)}\n')
    n_rec = 0
    hist = {k: [0] * (len(INFO_EDGES) - 1) for k in ('cohort', 'cohort_multipath', 'cohort_sv')}
    delta = [0] * (len(DELTA_EDGES) - 1)            # cohort INFO - old INFO of the first batch
    # per batch and record: [rows, old INFO == 0, old INFO == 0 and re-popped INFO > 0]
    zeroed = {'all': [0, 0, 0], 'multipath': [0, 0, 0], 'sv': [0, 0, 0]}
    for site, *rows in zip_longest(sites, *infos):
        if site is None or any(r is None for r in rows):
            sys.exit(f'ERROR: inputs end at different records after {n_rec}')
        s = site.rstrip('\n').split('\t')
        iid, raf, n_paths, n_paths_total = s[5:9]
        multipath = n_paths not in ('.', '1')
        sv = s[4].startswith('<') or abs(len(s[3]) - len(s[4])) >= 50

        ds = var = 0.0
        first_old = None
        for i, (nb, r) in enumerate(zip(n, rows)):
            af_b, info_b, old_b = r.rstrip('\n').split('\t')
            af_b = float(af_b)
            info_b = float(info_b)
            ds += 2 * nb * af_b
            var += (1.0 - info_b) * 2 * nb * af_b * (1.0 - af_b)
            old_b = num(old_b)
            if i == 0:
                first_old = old_b
            for k, on in (('all', True), ('multipath', multipath), ('sv', sv)):
                if on:
                    zeroed[k][0] += 1
                    zeroed[k][1] += old_b == 0.0
                    zeroed[k][2] += old_b == 0.0 and info_b > 0.0
        p = ds / n_haps
        info = 1.0 if p <= 0.0 or p >= 1.0 else min(1.0, max(0.0, 1.0 - var / (n_haps * p * (1.0 - p))))

        new = [f'ID={iid}']
        if raf != '.':
            new.append(f'RAF={raf}')
        new.append(f'AF={p:.6g};INFO={info:.3f}')
        if n_paths != '.':
            new.append(f'N_PATHS={n_paths}')
        if n_paths_total != '.':
            new.append(f'N_PATHS_TOTAL={n_paths_total}')
        new.append(f'DS_SUM={ds:.4f};VAR_SUM={var:.4f};N_HAPS={n_haps}')
        out.write('\t'.join(s[:5] + [';'.join(new)]) + '\n')

        n_rec += 1
        hist['cohort'][bin_of(info, INFO_EDGES)] += 1
        if multipath:
            hist['cohort_multipath'][bin_of(info, INFO_EDGES)] += 1
        if sv:
            hist['cohort_sv'][bin_of(info, INFO_EDGES)] += 1
        if first_old is not None:
            delta[bin_of(info - first_old, DELTA_EDGES)] += 1
    out.close()

    with open(a.report, 'w') as r:
        r.write(f'records\t{n_rec}\nbatches\t{len(n)}\nsamples\t{sum(n)}\n')
        for k, h in hist.items():
            for i, c in enumerate(h):
                r.write(f'info_hist\t{k}\t[{INFO_EDGES[i]},{INFO_EDGES[i + 1]})\t{c}\n')
        for i, c in enumerate(delta):
            r.write(f'delta_vs_old_first_batch_info\t[{DELTA_EDGES[i]},{DELTA_EDGES[i + 1]})\t{c}\n')
        for k, (rows_, z, back) in zeroed.items():
            r.write(f'old_batch_info_zero\t{k}\tbatch_rows={rows_}\told_zero={z}\tnow_positive={back}\n')


def apply(a):
    table = gzip.open(a.table, 'rt')
    head = table.readline().rstrip('\n').split('\t')
    if head[0] != '#cohort_info':
        sys.exit('ERROR: not a cohort_info table')
    out, vcf = sys.stdout, sys.stdin
    for line in vcf:
        if line.startswith('##INFO='):
            continue
        if line.startswith('#CHROM'):
            out.write('##cohort_info=<Description="AF and INFO pooled over %s %s, by cohort_info.py">\n' % (head[1], head[2]))
            out.write('\n'.join(TAGS) + '\n')
            out.write(line)
            break
        out.write(line)
    n_rec = 0
    for rec, t in zip_longest(vcf, table):
        if rec is None or t is None:
            sys.exit(f'ERROR: VCF and table end at different records after {n_rec}')
        f = rec.split('\t', 8)
        t = t.rstrip('\n').split('\t')
        if f[:5] != t[:5]:
            sys.exit(f'ERROR: record {f[:5]} does not match the cohort table {t[:5]}')
        f[7] = t[5]
        out.write('\t'.join(f))
        n_rec += 1


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest='cmd', required=True)
    p = sub.add_parser('pool')
    p.add_argument('--n-samples', required=True, help='comma-separated sample count per batch, in input order')
    p.add_argument('--sites', required=True)
    p.add_argument('--report', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('info_tsvs', nargs='+')
    q = sub.add_parser('apply')
    q.add_argument('--table', required=True)
    a = ap.parse_args()
    pool(a) if a.cmd == 'pool' else apply(a)


if __name__ == '__main__':
    main()
