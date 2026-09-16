#!/usr/bin/env python3
"""pop-glimpse2-joint-opt: INFO/AF/RAF of popped records.

Runs the binary on a 3-path synthetic bubble (tests/bubble.vcf, sites.vcf, id.split.vcf) and checks
against an independent re-implementation of the tool's own posterior model (pseudo-haplotype ->
odds normalisation -> per-atomic sum) that:
  - AF and INFO are the IMPUTE statistics of the emitted posterior (not a path's);
  - RAF is the atomic variant's panel frequency (variant a sits on two paths with RAF 0.30 and
    0.20; the panel says 0.50);
  - INFO_BUBBLE_FIRST reproduces the pre-fix INFO and INFO_BUBBLE_MAX the max over carrying paths;
  - GT/DS/GP are unchanged by the annotation change (fixed expected strings).
Usage: test_info_post.py <path to pop-glimpse2-joint-opt binary>
"""
import math, os, subprocess, sys
here = os.path.dirname(os.path.abspath(__file__))
binary = sys.argv[1]
out = subprocess.run([binary, f"{here}/id.split.vcf", f"{here}/sites.vcf"], stdin=open(f"{here}/bubble.vcf"), capture_output=True, text=True, check=True).stdout
recs = {}
for line in out.splitlines():
    if line.startswith('#'): continue
    f = line.split('\t'); info = dict(kv.split('=', 1) for kv in f[7].split(';'))
    recs[info['ID']] = (info, f[9:])
# independent expectation
paths = [("0|0",1,0,0),("1|0",0.1,0.9,0),("0|0",0.6,0.4,0),("1|1",0,0.05,0.95)], [("0|1",0.2,0.8,0),("0|0",1,0,0),("0|0",0.9,0.1,0),("0|0",1,0,0)], [("0|0",1,0,0),("0|0",1,0,0),("0|1",0.5,0.5,0),("0|0",1,0,0)]
ids = [["a"], ["a","b"], ["c"]]; path_info = [0.9, 0.4, 0.7]; panel_raf = {"a": 0.5, "b": 0.2, "c": 0.05}
def hap(gt, gp1, gp2):
    p = (gp2+gp1, gp2) if gt.startswith("1|0") else (gp2, gp2+gp1) if gt.startswith("0|1") else (gp2+gp1/2, gp2+gp1/2)
    return tuple(min(max(x, 1e-5), 1-1e-5) for x in p)
exp = {}
for s in range(4):
    hp = [hap(*[paths[a][s][0], paths[a][s][2], paths[a][s][3]]) for a in range(3)]
    z0 = 1+sum(h[0]/(1-h[0]) for h in hp); z1 = 1+sum(h[1]/(1-h[1]) for h in hp)
    for a in range(3):
        n0 = (hp[a][0]/(1-hp[a][0]))/z0; n1 = (hp[a][1]/(1-hp[a][1]))/z1
        for v in ids[a]:
            d = exp.setdefault(v, [[0, 0] for _ in range(4)]); d[s][0] += n0; d[s][1] += n1
fails = 0
def check(cond, msg):
    global fails
    if not cond: fails += 1; print("FAIL:", msg)
for v, d in exp.items():
    info, samples = recs[v]
    # the tool zeroes haplotype mass below 2e-4 (the 1e-5 clamp floor surviving normalisation)
    z = lambda x: 0.0 if x < 2e-4 else x
    sd = sum(z(x) for q in d for x in q); sv = sum(z(x)*(1-z(x)) for q in d for x in q); af = sd/8; e_info = 1 - sv/(8*af*(1-af))
    check(abs(float(info['AF']) - af) < 1e-6, f"{v}: AF {info['AF']} != {af:.6f}")
    check(abs(float(info['INFO']) - e_info) < 5e-4, f"{v}: INFO {info['INFO']} != {e_info:.3f}")
    check(abs(float(info['RAF']) - panel_raf[v]) < 1e-9, f"{v}: RAF {info['RAF']} != panel {panel_raf[v]}")
    carrying = [a for a in range(3) if v in ids[a]]
    check(abs(float(info['INFO_BUBBLE_MAX']) - max(path_info[a] for a in carrying)) < 5e-4, f"{v}: INFO_BUBBLE_MAX")
    check(abs(float(info['INFO_BUBBLE_FIRST']) - path_info[carrying[0]]) < 5e-4, f"{v}: INFO_BUBBLE_FIRST")
    check(int(info['N_PATHS']) == len(carrying), f"{v}: N_PATHS")
    check(abs(float(info['DS_SUM']) - sd) < 1e-5 and abs(float(info['VAR_SUM']) - sv) < 1e-5 and info['N_HAPS'] == '8', f"{v}: raw sums")
    # genotype fields exactly as the pre-fix tool wrote them (regression guard)
    expected_gt = {"a": ["0|1:0.8:0.2,0.8,0", "1|0:0.9:0.1,0.9,0", "0|0:0.364:0.666,0.303,0.031", "1|1:1.95:0.001,0.049,0.95"],
                   "b": ["0|1:0.8:0.2,0.8,0", "0|0:0:1,0,0", "0|0:0.063:0.938,0.061,0.001", "0|0:0:1,0,0"],
                   "c": ["0|0:0:1,0,0", "0|0:0:1,0,0", "0|0:0.434:0.566,0.434,0", "0|0:0:1,0,0"]}
    check(samples == expected_gt[v], f"{v}: GT/DS/GP changed: {samples}")
mono = recs.get("m")
check(mono is not None and mono[0]['AF'] == '0.000000' and mono[0]['INFO'] == '.' and float(mono[0]['DS_SUM']) == 0.0, f"monomorphic site: {mono[0] if mono else None}")
hdr = [l for l in out.splitlines() if l.startswith('##INFO=<ID=')]
for tag in ("ID", "RAF", "AF", "INFO", "INFO_BUBBLE_MAX", "INFO_BUBBLE_FIRST", "N_PATHS", "DS_SUM", "VAR_SUM", "N_HAPS"):
    check(sum(1 for l in hdr if l.startswith(f'##INFO=<ID={tag},')) == 1, f"header: {tag} declared {sum(1 for l in hdr if l.startswith(f'##INFO=<ID={tag},'))} times")
print("pop-glimpse2 INFO test:", "FAIL" if fails else "ok (3 variants, 4 samples, 10 header tags)")
sys.exit(1 if fails else 0)
