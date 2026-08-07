#!/usr/bin/env python3
"""
Build the sample manifest for the 1kGP leaveout validation run.

The leaveout panel holds out samples that remain in the full panel. Those that are also in the
1kGP 3202 have a public DRAGEN 4.4.7 gVCF on AWS Open Data, so they can be imputed against the
leaveout panel and scored against the full panel's assembly-based genotypes without moving
anything between Terra and VWB. This intersects the panel's leaveout list against what is
actually published and writes the files the staging and preprocessing workflows read.

Run it twice. The first pass picks the samples and writes the staging workflow's inputs; the
second, after staging has run, turns that workflow's outputs into the FOFNs PreprocessPLsGVCF
reads. Staged paths are Cromwell execution paths and cannot be guessed, which is why the second
pass needs the outputs JSON.

    scripts/build-1kgp-validation-manifest.py \\
        --leaveout-samples leaveout.samples.txt \\
        --out-dir manifests/ --prefix 1kgp-validation.chr22 --regions chr22

    scripts/build-1kgp-validation-manifest.py \\
        --leaveout-samples leaveout.samples.txt \\
        --out-dir manifests/ --prefix 1kgp-validation.chr22 --regions chr22 \\
        --staged-outputs-json staging.outputs.json

Take the leaveout list from the panels themselves rather than from notes, so the manifest cannot
drift from what was actually built:

    bcftools query -l <full panel>.bcf     | sort > all.txt
    bcftools query -l <leaveout panel>.bcf | sort > kept.txt
    comm -23 all.txt kept.txt > leaveout.samples.txt

The available-sample list is cached in resources/1kgp-dragen-v4-4-7.samples.txt so this runs
without network access. Regenerate it with --refresh-available (needs the aws CLI; the bucket
is public, so no credentials).

Emits, into --out-dir:
    <prefix>.samples.txt        one sample per line, the intersection
    <prefix>.staging.json       inputs for StageDragen1kGPGvcfs
    <prefix>.unavailable.txt    leaveout samples with no published DRAGEN gVCF
    <prefix>.gvcfs.fofn         staged gVCF paths, parallel to samples.txt   (second pass)
    <prefix>.gvcf_idxs.fofn     staged index paths, parallel to samples.txt  (second pass)
"""

import argparse
import json
import os
import subprocess
import sys

# The published analysis directory. Pinned rather than globbed: the bucket also carries
# per-run directories (<sample>_dragen-germline-v4-4-7-<uuid>) holding the same data, and
# matching those would double every sample.
S3_BUCKET = "s3://1000genomes-dragen-v4-4-7"
S3_INDIVIDUALS = f"{S3_BUCKET}/data/individuals/hg38-alt_masked.cnv.graph.hla.methyl_cg.rna-11-r5.0-2"

DEFAULT_AVAILABLE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "resources", "1kgp-dragen-v4-4-7.samples.txt")


def refresh_available(path):
    """Re-list the bucket's sample directories. ~6400 entries, two per sample."""
    proc = subprocess.run(
        ["aws", "s3", "ls", "--no-sign-request", f"{S3_INDIVIDUALS}/"],
        capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"aws s3 ls failed: {proc.stderr.strip()}")

    samples = set()
    for line in proc.stdout.splitlines():
        # "                           PRE HG00096/" -- keep the bare IDs, drop the run dirs.
        parts = line.split()
        if len(parts) == 2 and parts[0] == "PRE":
            name = parts[1].rstrip("/")
            if "_dragen-germline" not in name:
                samples.add(name)

    with open(path, "w") as f:
        for s in sorted(samples):
            f.write(s + "\n")
    return samples


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--leaveout-samples", required=True,
                    help="samples held out of the leaveout panel, one per line")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--prefix", default="1kgp-validation.chr22",
                    help="also the staging workflow's output_prefix; carry the chromosome in it")
    ap.add_argument("--regions", default="chr22", help="comma-separated, passed to bcftools -r")
    ap.add_argument("--staged-outputs-json",
                    help="StageDragen1kGPGvcfs outputs JSON. Given this, the FOFNs are built "
                         "from the paths the workflow actually wrote. Omit on the first pass.")
    ap.add_argument("--available", default=DEFAULT_AVAILABLE)
    ap.add_argument("--refresh-available", action="store_true")
    a = ap.parse_args()

    if a.refresh_available:
        available = refresh_available(a.available)
        print(f"refreshed {a.available}: {len(available)} samples", file=sys.stderr)
    else:
        with open(a.available) as f:
            available = {ln.strip() for ln in f if ln.strip()}

    with open(a.leaveout_samples) as f:
        leaveout = [ln.strip() for ln in f if ln.strip()]

    # Order is the leaveout file's, not sorted: PreprocessPLsGVCF pairs sample_names_file
    # against the gVCF arrays positionally, so all three files must stay in lockstep.
    keep = [s for s in leaveout if s in available]
    drop = [s for s in leaveout if s not in available]

    print(f"leaveout samples:   {len(leaveout)}")
    print(f"  with DRAGEN gVCF: {len(keep)}")
    print(f"  unavailable:      {len(drop)}"
          + (f"  ({', '.join(drop[:8])}{', ...' if len(drop) > 8 else ''})" if drop else ""))

    if not keep:
        raise SystemExit("no samples selected; nothing to write")

    os.makedirs(a.out_dir, exist_ok=True)
    base = os.path.join(a.out_dir, a.prefix)
    written = ["samples.txt", "unavailable.txt", "staging.json"]

    with open(f"{base}.samples.txt", "w") as f:
        f.writelines(s + "\n" for s in keep)
    with open(f"{base}.unavailable.txt", "w") as f:
        f.writelines(s + "\n" for s in drop)

    with open(f"{base}.staging.json", "w") as f:
        json.dump({
            "StageDragen1kGPGvcfs.sample_ids": keep,
            "StageDragen1kGPGvcfs.regions": a.regions.split(","),
            "StageDragen1kGPGvcfs.output_prefix": a.prefix,
        }, f, indent=2)
        f.write("\n")

    if a.staged_outputs_json:
        # Cromwell delocalizes to a run-specific execution path, so the staged locations cannot
        # be predicted from the bucket alone -- they have to be read back. Match on the sample
        # ID in the basename rather than trusting array order, because a partially-rerun
        # scatter can come back reordered.
        outs = json.load(open(a.staged_outputs_json))
        outs = outs.get("outputs", outs)

        def index_by_sample(paths, what):
            found = {}
            for path in paths:
                name = os.path.basename(path)
                hits = [s for s in keep if f".{s}." in name]
                if len(hits) != 1:
                    raise SystemExit(f"cannot pin a single sample to {what} {name}: {hits}")
                found[hits[0]] = path
            return found

        gvcf_by_sample = index_by_sample(outs["StageDragen1kGPGvcfs.sliced_gvcfs"], "gvcf")
        idx_by_sample = index_by_sample(outs["StageDragen1kGPGvcfs.sliced_gvcf_idxs"], "index")

        missing = [s for s in keep if s not in gvcf_by_sample or s not in idx_by_sample]
        if missing:
            raise SystemExit(f"staging produced nothing for {len(missing)} samples: {missing[:8]}")

        # Same order as samples.txt: PreprocessPLsGVCF pairs the three files positionally.
        with open(f"{base}.gvcfs.fofn", "w") as f:
            f.writelines(gvcf_by_sample[s] + "\n" for s in keep)
        with open(f"{base}.gvcf_idxs.fofn", "w") as f:
            f.writelines(idx_by_sample[s] + "\n" for s in keep)
        written += ["gvcfs.fofn", "gvcf_idxs.fofn"]
    else:
        print("no --staged-outputs-json: wrote the staging inputs only. Rerun with it once "
              "StageDragen1kGPGvcfs finishes to get the FOFNs.")

    print(f"wrote {base}.{{{','.join(written)}}}")


if __name__ == "__main__":
    main()
