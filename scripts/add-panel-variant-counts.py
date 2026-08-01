#!/usr/bin/env python3
"""
Backfill n_variants into an existing chunked_panel.json.

GLIMPSE2ChunkAndSplitPanel emits the field now, but a panel chunked before that change has a
JSON without it, and GLIMPSE2FromPreprocessedPLsJoint requires it. Rechunking the panel to get
one number per shard would mean regenerating hundreds of GB of .bin files; this reads the
counts off the sites-only BCFs the panel already ships and rewrites the JSON in place.

    scripts/add-panel-variant-counts.py \\
        --chunked-panel  aou_lr_phase2_v1.chunked_panel.json \\
        --panel-resources aou_lr_phase2_v1.pop_glimpse2_panel_resources.json \\
        --out            aou_lr_phase2_v1.chunked_panel.with_counts.json

Counts each shard's input (buffered) region with --regions-overlap 0, which is how GLIMPSE2
selects sites: a record counts when its POS falls inside, regardless of REF span. This
reproduces GLIMPSE2's reported L exactly (verified on chr22 s0 = 657,784, chr20 s4 = 965,039,
chr7 s12 = 1,346,888).

Needs bcftools and read access to the panel. Set GCS_OAUTH_TOKEN for gs:// paths.
"""

import argparse
import concurrent.futures as cf
import json
import re
import subprocess
import sys


def count_region(vcf, region):
    """
    Panel variants whose POS falls inside `region`.

    Streams rather than capturing: a dense shard is ~1M records, and with --jobs 8 buffering
    every one of them before counting holds hundreds of MB for a number.
    """
    proc = subprocess.Popen(
        ["bcftools", "view", "--no-version", "--threads", "1", "-H",
         "-r", region, "--regions-overlap", "0", vcf],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    n = sum(1 for _ in proc.stdout)
    proc.stdout.close()
    err = proc.stderr.read()
    proc.stderr.close()
    if proc.wait() != 0:
        raise RuntimeError(f"bcftools failed on {region}: {err.strip().splitlines()[-1:]}")
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunked-panel", required=True)
    ap.add_argument("--panel-resources", required=True,
                    help="pop_glimpse2_panel_resources.json, for the sites-only BCF paths")
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()

    panel = json.load(open(a.chunked_panel))
    res = json.load(open(a.panel_resources))

    for chrom in sorted(panel, key=lambda c: int(re.sub(r"\D", "", c) or 0)):
        entry = panel[chrom]
        if "n_variants" in entry:
            print(f"  {chrom:<6} already has n_variants, skipping", file=sys.stderr)
            continue
        sites = res[chrom]["panel_bubble_split_sites_only_vcf"]
        regions = entry["input_regions"]
        with cf.ThreadPoolExecutor(max_workers=a.jobs) as pool:
            counts = list(pool.map(lambda r: count_region(sites, r), regions))
        if any(c == 0 for c in counts):
            raise SystemExit(f"{chrom}: a region counted zero variants -- refusing to write")
        entry["n_variants"] = counts
        print(f"  {chrom:<6} {len(counts):>3} shards   L {min(counts):>9,} .. {max(counts):>9,}",
              file=sys.stderr)

    missing = [c for c, e in panel.items() if len(e["n_variants"]) != len(e["input_regions"])]
    if missing:
        raise SystemExit(f"length mismatch on {missing} -- refusing to write")

    with open(a.out, "w") as fh:
        json.dump(panel, fh, indent=2)
    print(f"wrote {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
