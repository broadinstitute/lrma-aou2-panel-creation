#!/usr/bin/env bash
# CI test for CohortInfoFromBubbles on synthetic data (48 samples, 4 batches of 12, one chromosome).
#   1. the pop that made the deliverables (joint-opt before 79a6f70), the info-fix joint-opt and lrma-sv-imputation-utils main
#      write identical GT:DS:GP, which RepopBatch relies on
#   2. pooling the per-batch AF/INFO and pasting the batches gives the same file as popping all 48 samples at once: GT:DS:GP
#      identical, AF and INFO within the rounding bound (compare.py)
#   3. with CROMWELL_JAR set, the WDL itself (re-pop and pooling) runs and its base file matches step 2
# Usage: UTILS=path/to/lrma-sv-imputation-utils bash test/cohort-info/run.sh
set -euo pipefail
REPO=$(cd "$(dirname "$0")/../.." && pwd)
UTILS=${UTILS:?set UTILS to a checkout of lrma-sv-imputation-utils}
W=$(mktemp -d); echo "work dir $W"
cd "$W"

say() { echo "== $*"; }
gtmd5() { bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT:%DS:%GP]\n' "$1" | md5sum | cut -d' ' -f1; }

say "build pop-glimpse2 (utils main, info-fix joint-opt, pre-fix joint-opt) and paste-vcfs"
cargo build -q --release --locked --manifest-path "$UTILS/resources/pop-glimpse2/Cargo.toml"
cargo build -q --release --locked --manifest-path "$UTILS/resources/paste-vcfs/Cargo.toml"
cp "$UTILS/resources/pop-glimpse2/target/release/pop-glimpse2" pop.sam
cp "$UTILS/resources/paste-vcfs/target/release/paste-vcfs" paste-vcfs
cargo build -q --release --locked --manifest-path "$REPO/resources/pop-glimpse2/Cargo.toml" --bin pop-glimpse2-joint-opt
cp "$REPO/resources/pop-glimpse2/target/release/pop-glimpse2-joint-opt" pop.infofix
git -C "$REPO" worktree add -q --detach "$W/prefix" 79a6f70^
cargo build -q --release --locked --manifest-path "$W/prefix/resources/pop-glimpse2/Cargo.toml" --bin pop-glimpse2-joint-opt
cp "$W/prefix/resources/pop-glimpse2/target/release/pop-glimpse2-joint-opt" pop.old
git -C "$REPO" worktree remove --force "$W/prefix"
if [ -n "${OUT_BIN:-}" ]; then mkdir -p "$OUT_BIN" && cp pop.sam "$OUT_BIN/pop-glimpse2" && cp paste-vcfs "$OUT_BIN/" && (cd "$OUT_BIN" && sha256sum pop-glimpse2 paste-vcfs > SHA256SUMS); fi

say "fixture"
python3 "$REPO/test/cohort-info/make_fixture.py" fx
bcftools view fx/sites.vcf -Oz -o sites.vcf.gz --write-index=tbi
bcftools view fx/ids.vcf -Oz -o ids.vcf.gz --write-index=tbi
REGIONS="chr1:1-700000 chr1:700001-1400000 chr1:1400001-2000000"

# pop $1 over bubble BCF $2 region by region, as production did, into $3
pop_regions() {
    local bin=$1 in=$2 out=$3 k=0 list=()
    for r in $REGIONS; do
        k=$((k+1))
        bcftools view -r "$r" --regions-overlap 0 sites.vcf.gz -Oz -o s.$k.vcf.gz
        bcftools view -r "$r" --regions-overlap 0 "$in" | "./$bin" ids.vcf.gz s.$k.vcf.gz 2> pop.log > raw.vcf
        bcftools sort -Ob -o "$out.$k.bcf" raw.vcf || {
            echo "FAIL sort of $bin output on $r"; grep '^##' raw.vcf | grep -v '^##contig'
            grep -v '^#' raw.vcf | head -3 | cut -c1-400; exit 1; }
        list+=("$out.$k.bcf")
    done
    bcftools concat --naive "${list[@]}" -Ob -o "$out" && bcftools index -f "$out" && rm "${list[@]}"
}

BATCHES="batch-001 batch-002 batch-003 batch-004"
for b in $BATCHES; do
    bcftools view fx/main.$b.vcf -Ob -o bubble.$b.bcf --write-index
    pop_regions pop.old bubble.$b.bcf popped.$b.bcf          # stands in for the deliverable
done
bcftools view fx/main.vcf -Ob -o bubble.all.bcf --write-index
pop_regions pop.sam bubble.all.bcf reference.bcf

say "1. GT:DS:GP identical across the three pops"
for b in $BATCHES; do
    pop_regions pop.sam bubble.$b.bcf repop.$b.bcf
    pop_regions pop.infofix bubble.$b.bcf infofix.$b.bcf
    o=$(gtmd5 popped.$b.bcf); s=$(gtmd5 repop.$b.bcf); i=$(gtmd5 infofix.$b.bcf)
    [ "$o" = "$s" ] && [ "$o" = "$i" ] || { echo "FAIL $b: old $o, utils $s, info-fix $i"; exit 1; }
    echo "$b: $(bcftools index -n popped.$b.bcf) records, GT:DS:GP md5 equal"
done

say "2. pool, paste, compare with one pop over all samples"
infos=()
for b in $BATCHES; do
    paste <(bcftools query -f '%INFO/AF\t%INFO/INFO\n' repop.$b.bcf) <(bcftools query -f '%INFO/INFO\n' popped.$b.bcf) | gzip > info.$b.tsv.gz
    infos+=(info.$b.tsv.gz)
done
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/ID\t%INFO/RAF\t%INFO/N_PATHS\t%INFO/N_PATHS_TOTAL\n' repop.batch-001.bcf | gzip > sites.batch-001.tsv.gz
bcftools view popped.batch-001.bcf | \
    python3 "$REPO/resources/cohort-info/cohort_info.py" --n-samples 12,12,12,12 --sites sites.batch-001.tsv.gz --report report.txt "${infos[@]}" | \
    bcftools view -Ob -o base.bcf --write-index
cat report.txt
tags=ID,RAF,AF,INFO,N_PATHS,N_PATHS_TOTAL,DS_SUM,VAR_SUM,N_HAPS
./paste-vcfs --info $tags --format GT,DS,GP -o merged.bcf base.bcf popped.batch-002.bcf popped.batch-003.bcf popped.batch-004.bcf
bcftools index -f merged.bcf
[ "$(bcftools query -l merged.bcf | md5sum)" = "$(bcftools query -l reference.bcf | md5sum)" ] || { echo "FAIL sample order"; exit 1; }
[ "$(gtmd5 merged.bcf)" = "$(gtmd5 reference.bcf)" ] || { echo "FAIL merged GT:DS:GP differ from one pop over all samples"; exit 1; }
q='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/INFO\n'
bcftools query -f "$q" merged.bcf > merged.tsv; bcftools query -f "$q" reference.bcf > reference.tsv
python3 "$REPO/test/cohort-info/compare.py" merged.tsv reference.tsv
[ "$(bcftools query -f '%INFO/N_HAPS\n' merged.bcf | sort -u)" = 96 ] || { echo "FAIL N_HAPS"; exit 1; }

if [ -n "${CROMWELL_JAR:-}" ]; then
    say "3. the WDL: re-pop and pooling (run_paste=false) on the same fixture"
    : > batches.tsv
    for b in $BATCHES; do printf '%s\t%s\t%s\t%s\t%s\n' $b $W/bubble.$b.bcf $W/bubble.$b.bcf.csi $W/popped.$b.bcf $W/popped.$b.bcf.csi >> batches.tsv; done
    python3 - "$W" "$REGIONS" > panel.json <<'EOF'
import json, sys
w, regions = sys.argv[1], sys.argv[2].split()
print(json.dumps({"chr1": {"panel_bubble_split_sites_only_vcf": f"{w}/sites.vcf.gz", "panel_bubble_split_sites_only_vcf_idx": f"{w}/sites.vcf.gz.tbi",
                           "panel_id_split_vcf_gz": f"{w}/ids.vcf.gz", "panel_id_split_vcf_gz_tbi": f"{w}/ids.vcf.gz.tbi", "pop_regions": regions}}))
EOF
    DOCKER_JSON=""
    [ -z "${DOCKER:-}" ] || DOCKER_JSON=", \"CohortInfoFromBubbles.docker\": \"$DOCKER\""
    cat > inputs.json <<EOF
{"CohortInfoFromBubbles.chromosome": "chr1", "CohortInfoFromBubbles.batches_tsv": "$W/batches.tsv",
 "CohortInfoFromBubbles.pop_glimpse2_panel_resources_json": "$W/panel.json",
 "CohortInfoFromBubbles.pop_glimpse2_binary": "$W/pop.sam", "CohortInfoFromBubbles.cohort_info_script": "$REPO/resources/cohort-info/cohort_info.py",
 "CohortInfoFromBubbles.paste_vcfs_binary": "$W/paste-vcfs", "CohortInfoFromBubbles.paste_regions": ["chr1"],
 "CohortInfoFromBubbles.expected_samples": 12, "CohortInfoFromBubbles.run_paste": false,
 "CohortInfoFromBubbles.RepopBatch.cpu": 2, "CohortInfoFromBubbles.RepopBatch.mem_gb": 4$DOCKER_JSON}
EOF
    java -jar "$WOMTOOL_JAR" validate "$REPO/wdl/methods/imputation/CohortInfoFromBubbles.wdl" -i inputs.json
    java -jar "$CROMWELL_JAR" run "$REPO/wdl/methods/imputation/CohortInfoFromBubbles.wdl" -i inputs.json -m meta.json > cromwell.log 2>&1 \
        || { tail -100 cromwell.log; find cromwell-executions -name stderr -exec tail -20 {} +; exit 1; }
    wbase=$(python3 -c 'import json;print(json.load(open("meta.json"))["outputs"]["CohortInfoFromBubbles.base_cohort_bcf"])')
    [ "$(bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\n' "$wbase" | md5sum)" = "$(bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO\n' base.bcf | md5sum)" ] \
        || { echo "FAIL WDL base INFO differs from the direct run"; exit 1; }
    [ "$(gtmd5 "$wbase")" = "$(gtmd5 popped.batch-001.bcf)" ] || { echo "FAIL WDL base genotypes"; exit 1; }
    find cromwell-executions -name report.txt -path '*RepopBatch*' -exec cat {} +
    echo "WDL base file matches the direct run"
fi
say "PASS"
