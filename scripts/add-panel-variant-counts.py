#!/usr/bin/env python3
"""
Write region-keyed variant counts into an existing chunked_panel.json.

GLIMPSE2ChunkAndSplitPanel emits the field now, but a panel chunked before that change has a
JSON without them, and GLIMPSE2FromPreprocessedPLsJoint requires them. Rechunking the panel to
get one number per shard would mean regenerating hundreds of GB of .bin files; this reads the
counts off the sites-only BCFs the panel already ships and writes a new JSON.

Counts are keyed by the exact input-region string. Parallel arrays are deliberately unsupported:
a production array was shuffled while retaining the correct length and values, silently assigning
chr7 shard 12's 1,346,888-site binary a 420,889-site memory estimate.

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
import os
import re
import subprocess
import sys
import tempfile


def count_region(vcf, region):
    """
    Panel variants whose POS falls inside `region`.

    Streams rather than capturing: a dense shard is ~1M records, and with --jobs 8 buffering
    every one of them before counting holds hundreds of MB for a number.
    """
    # Drain stderr to disk while stdout is streamed. Reading stdout before a piped stderr can
    # deadlock if bcftools fills the stderr pipe and blocks before closing stdout.
    with tempfile.TemporaryFile(mode="w+") as stderr:
        proc = subprocess.Popen(
            ["bcftools", "view", "--no-version", "--threads", "1", "-H",
             "-r", region, "--regions-overlap", "0", vcf],
            stdout=subprocess.PIPE, stderr=stderr, text=True)
        n = sum(1 for _ in proc.stdout)
        proc.stdout.close()
        returncode = proc.wait()
        if returncode != 0:
            stderr.seek(0)
            detail = stderr.read().strip().splitlines()[-1:]
            raise RuntimeError(f"bcftools failed on {region}: {detail}")
    return n


def validate_panel_entry(chrom, entry):
    required = ("input_regions", "output_regions", "panel_split_chunk_bins")
    missing = [field for field in required if field not in entry]
    if missing:
        raise SystemExit(f"{chrom}: missing required fields {missing}")

    regions = entry["input_regions"]
    lengths = {field: len(entry[field]) for field in required}
    if not regions:
        raise SystemExit(f"{chrom}: no input regions -- refusing to count")
    if len(set(lengths.values())) != 1:
        raise SystemExit(f"{chrom}: shard arrays have different lengths {lengths}")

    duplicates = sorted(region for region in set(regions) if regions.count(region) > 1)
    if duplicates:
        raise SystemExit(f"{chrom}: duplicate input regions {duplicates}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunked-panel", required=True)
    ap.add_argument("--panel-resources", required=True,
                    help="pop_glimpse2_panel_resources.json, for the sites-only BCF paths")
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()
    if a.jobs < 1:
        ap.error("--jobs must be at least 1")

    with open(a.chunked_panel) as fh:
        panel = json.load(fh)
    with open(a.panel_resources) as fh:
        res = json.load(fh)

    for chrom in sorted(panel, key=lambda c: int(re.sub(r"\D", "", c) or 0)):
        entry = panel[chrom]
        validate_panel_entry(chrom, entry)
        if chrom not in res:
            raise SystemExit(f"{chrom}: missing from panel resources")
        sites = res[chrom]["panel_bubble_split_sites_only_vcf"]
        regions = entry["input_regions"]
        with cf.ThreadPoolExecutor(max_workers=a.jobs) as pool:
            counts = list(pool.map(lambda r: count_region(sites, r), regions))
        if any(c == 0 for c in counts):
            raise SystemExit(f"{chrom}: a region counted zero variants -- refusing to write")
        entry.pop("n_variants", None)
        entry["n_variants_by_region"] = dict(zip(regions, counts, strict=True))
        print(f"  {chrom:<6} {len(counts):>3} shards   L {min(counts):>9,} .. {max(counts):>9,}",
              file=sys.stderr)

    invalid = [c for c, e in panel.items()
               if len(e["n_variants_by_region"]) != len(e["input_regions"])
               or set(e["n_variants_by_region"]) != set(e["input_regions"])
               or any(n <= 0 for n in e["n_variants_by_region"].values())]
    if invalid:
        raise SystemExit(f"region/count mismatch on {invalid} -- refusing to write")

    # Replace atomically so an interrupted write cannot corrupt the only finished manifest.
    out = os.path.abspath(a.out)
    with tempfile.NamedTemporaryFile(
            mode="w", dir=os.path.dirname(out), prefix=f".{os.path.basename(out)}.",
            suffix=".tmp", delete=False) as fh:
        temporary_out = fh.name
        json.dump(panel, fh, indent=2)
        fh.write("\n")
    os.replace(temporary_out, out)
    print(f"wrote {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
