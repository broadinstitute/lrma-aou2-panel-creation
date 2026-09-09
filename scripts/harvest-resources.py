#!/usr/bin/env python3
"""
Summarise the [RESOURCE] lines emitted by GLIMPSE2FromPreprocessedPLsJoint: what each task
actually used versus what it was given, and whether the phase memory bound holds.

    gsutil cat 'gs://<bucket>/**/stderr' | scripts/harvest-resources.py
    scripts/harvest-resources.py logs/*.log

Reads stdin when given no files, and ignores non-[RESOURCE] lines.
"""

import fileinput
import os
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
        # disk_used_gb is sampled once at task exit, so a transient mid-run peak (localization,
        # sort spill, checkpoint) can be higher. Report it as a lower bound, not a high-water mark.
        print(f"  disk at exit    : {worst_u:.0f} / {worst_t:.0f} GB ({worst_u / worst_t:.0%} of total;"
              f" exit-time sample = LOWER bound, do not size disk down on this alone)")


DEFAULT_WDL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "wdl",
                           "methods", "imputation", "GLIMPSE2FromPreprocessedPLsJoint.wdl")


def wdl_model(path=DEFAULT_WDL):
    """
    Recover (constant, bytes-per-site-per-thread-per-state) from the WDL itself rather than
    hardcoding them. These drifted once already -- the WDL was refitted from 4.0 to 8.0 while
    this script still asserted 4.0, which would have reported a bound breach on every shard.
    """
    try:
        m = re.search(r"Int computed_mem_gb\s*=\s*(\d+(?:\.\d+)?)\s*\+\s*"
                      r"ceil\(\(\(\((\d+(?:\.\d+)?)\s*\*\s*phase_threads\)",
                      open(path).read())
        if m:
            return float(m.group(1)), float(m.group(2))
    except OSError:
        pass
    return None, None


def check_phase_model(recs):
    """
    Phase asks for CONST + COEFF * threads * Kpbwt * L / 1e9, both read from the WDL.
    Recover the coefficient actually needed from what the shards did.
    """
    const, coeff = wdl_model()
    pts = []
    # Only fit on cleanly-successful (rc=0) shards. Classify the rest by cause, because they
    # mean different things for sizing:
    #   * rc=137 (SIGKILL, ~OOM): a memory-bound BREACH. peak_rss is censored (a lower bound on
    #     demand), so it is excluded from the fit AND flagged -- do not size memory down.
    #   * any other non-zero rc (e.g. rc=2): failed for a NON-memory reason. Also unreliable for
    #     the fit, so excluded, but it is NOT evidence about memory and must not read as a breach
    #     (that would discourage valid downsizing or motivate needless memory increases).
    oom = 0
    other_fail = 0
    for r in recs:
        rc = r.get("rc")
        if rc == "137":
            oom += 1
            continue
        if rc is not None and rc != "0":
            other_fail += 1
            continue
        peak, L = num(r, "peak_rss_gib"), num(r, "n_variants")
        t, kp = num(r, "threads"), num(r, "kpbwt")
        if None in (peak, L, t, kp) or L * t * kp == 0:
            continue
        # Both terms in GiB-derived bytes. The WDL computes CONST + matrix/1e9 and hands the
        # result to Cromwell as GiB, so its fixed term is CONST GiB while its matrix term is
        # decimal -- a ~7% conservative bias baked into the request. Subtracting CONST * GIB
        # keeps the recovered coefficient from inheriting that skew as a false bound breach.
        pts.append((((peak * GIB) - (const or 2) * GIB) / (L * t * kp), L, peak))
    print("\n=== phase memory model ===")
    if oom:
        print(f"  *** {oom} shard(s) OOM-killed (rc=137): the memory bound was BREACHED. Their peak")
        print(f"      RSS is censored (a lower bound), excluded from the fit -- which is therefore")
        print(f"      itself a LOWER bound. Do not size memory down.")
    if other_fail:
        print(f"  ({other_fail} shard(s) failed at a non-zero rc other than 137 -- a non-memory")
        print(f"   failure; excluded from the fit, NOT counted as a memory breach.)")
    if not pts:
        print("  no successful (rc=0) shards with peak RSS -- nothing to fit")
        return
    coeffs = sorted(c for c, _, _ in pts)
    if coeff is None:
        print("  could not read the sizing constants from the WDL; skipping comparison")
        return
    print(f"  WDL assumes         : {coeff:.2f} bytes per site per thread per state,"
          f" plus {const:.0f} GiB")
    print(f"  observed (rc=0 only): median {statistics.median(coeffs):.2f}  max {max(coeffs):.2f}"
          f"   (n={len(coeffs)})")
    if max(coeffs) > coeff:
        print("  *** the assumption is NOT an upper bound -- a successful shard exceeded it."
              " Do not size down.")
    elif oom:
        print("  the surviving shards fit under the assumption, but the OOM kill(s) above mean the"
              " memory bound did not hold in practice. Do not size down.")
    else:
        print(f"  holds with margin; the worst shard needed {max(coeffs):.2f},"
              f" {(1 - max(coeffs) / coeff):.0%} below what is requested")


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

    return 0


if __name__ == "__main__":
    sys.exit(main())
