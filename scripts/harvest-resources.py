#!/usr/bin/env python3
"""
Turn [RESOURCE] lines emitted by GLIMPSE2FromPreprocessedPLsJoint into the numbers nobody in
this pipeline has had yet: what each task actually used, versus what it was given.

Every memory figure in that WDL is a bound or an empirical cap. Phase is sized from an upper
bound known to overshoot (n_states is frequently well below Kpbwt); ligate's 32 GiB has no
model behind it at all, three having been proposed and refuted. This script exists so the next
batch settles both, without a dedicated profiling run.

Usage:
    gsutil cat 'gs://<bucket>/**/stdout' | scripts/harvest-resources.py
    scripts/harvest-resources.py cromwell-logs/*.log
    grep -rh '\\[RESOURCE\\]' logs/ | scripts/harvest-resources.py

Reads stdin when given no files. Ignores everything that is not a [RESOURCE] line, so it is
safe to pipe raw logs at it.
"""

import fileinput
import re
import statistics
import sys

GIB = 1024 ** 3


def parse(lines):
    """Extract k=v pairs from [RESOURCE] lines. Unknown keys are kept; order is irrelevant."""
    out = []
    for line in lines:
        if "[RESOURCE]" not in line:
            continue
        rec = dict(re.findall(r"(\w+)=([^\s]+)", line.split("[RESOURCE]", 1)[1]))
        if rec.get("task"):
            out.append(rec)
    return out


def num(rec, key):
    """Numeric field, or None for missing/NA -- instrumentation is best-effort by design."""
    v = rec.get(key)
    if v is None or v == "NA":
        return None
    try:
        return float(v)
    except ValueError:
        return None


def summarise(task, recs):
    peaks, ratios, walls = [], [], []
    for r in recs:
        peak, req = num(r, "peak_rss_gib"), num(r, "requested_mem_gib")
        w = num(r, "wall_s")
        if w is not None:
            walls.append(w)
        if peak is None:
            continue
        peaks.append(peak)
        if req:
            ratios.append(peak / req)

    killed = [r for r in recs if r.get("rc") == "137"]
    print(f"\n=== {task}  ({len(recs)} tasks, {len(peaks)} with peak RSS) ===")
    if killed:
        print(f"  *** {len(killed)} task(s) exited rc=137 (SIGKILL/OOM). These are the rows that")
        print(f"      matter most -- peak here is what the kernel killed the task at:")
        for r in killed[:5]:
            print(f"      requested {r.get('requested_mem_gib','?')} GiB, limit {r.get('limit_gib','?')},"
                  f" peak {r.get('peak_rss_gib','?')} GiB, region {r.get('region','-')}")
    if not peaks:
        print("  no peak RSS captured -- cgroup unreadable in this environment")
    else:
        req_set = sorted({int(x) for x in (num(r, "requested_mem_gib") for r in recs) if x})
        print(f"  requested GiB   : {req_set if len(req_set) < 12 else f'{req_set[0]}..{req_set[-1]}'}")
        print(f"  peak RSS GiB    : min {min(peaks):.2f}  median {statistics.median(peaks):.2f}  max {max(peaks):.2f}")
        if ratios:
            worst = max(ratios)
            print(f"  used / requested: median {statistics.median(ratios):.0%}  worst {worst:.0%}")
            if worst > 1.0:
                print(f"  *** OVER REQUEST on {sum(1 for x in ratios if x > 1.0)} task(s):"
                      " the request is too small and only cgroup slack prevented an OOM")
            elif worst < 0.5:
                print(f"  *** OVER-PROVISIONED: worst case used {worst:.0%} of its request."
                      " There is room to size down.")
    if walls:
        print(f"  wall seconds    : median {statistics.median(walls):.0f}  max {max(walls):.0f}")

    disk = [(num(r, "disk_used_gb"), num(r, "disk_total_gb")) for r in recs]
    disk = [(u, t) for u, t in disk if u and t]
    if disk:
        worst_u, worst_t = max(disk, key=lambda p: p[0] / p[1])
        print(f"  disk high-water : {worst_u:.0f} / {worst_t:.0f} GB ({worst_u / worst_t:.0%} of the smallest margin)")


def check_phase_model(recs):
    """
    Phase asks for 8 + 4 * threads * Kpbwt * L / 1e9 GB, i.e. an assumed 4 bytes per
    site per thread per state. Recover that coefficient from what actually happened.
    """
    pts = []
    for r in recs:
        peak, L = num(r, "peak_rss_gib"), num(r, "n_variants")
        t, kp = num(r, "threads"), num(r, "kpbwt")
        if None in (peak, L, t, kp) or L * t * kp == 0:
            continue
        pts.append((((peak * GIB) - 8 * 1e9) / (L * t * kp), L, peak))
    if not pts:
        return
    coeffs = sorted(c for c, _, _ in pts)
    print("\n=== phase memory model ===")
    print(f"  assumed coefficient : 4.00 bytes per site per thread per state (upper bound)")
    print(f"  observed            : median {statistics.median(coeffs):.2f}  max {max(coeffs):.2f}"
          f"   (n={len(coeffs)})")
    if max(coeffs) > 4.0:
        print("  *** the bound is NOT a bound -- some shard exceeded it. Do not size down.")
    else:
        print(f"  the bound holds; sizing on {max(coeffs):.2f} instead of 4.00 would cut the"
              f" matrix term by {(1 - max(coeffs) / 4.0):.0%}")


def main():
    files = sys.argv[1:]
    recs = parse(fileinput.input(files=files or ("-",)))
    if not recs:
        print("no [RESOURCE] lines found", file=sys.stderr)
        return 1

    for task in sorted({r["task"] for r in recs}):
        rows = [r for r in recs if r["task"] == task]
        summarise(task, rows)
        if task == "GLIMPSE2Phase":
            check_phase_model(rows)

    print("\nLigate is the one to look at first: 32 GiB is a cap with no surviving model,"
          "\nand a single peak from a seam that used to die at 12 GiB settles it.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
