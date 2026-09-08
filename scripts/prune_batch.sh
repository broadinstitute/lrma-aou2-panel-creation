#!/usr/bin/env bash
# prune_batch.sh <BATCH>   (e.g. prune_batch.sh 003)
#
# After a batch's 22 chromosomes finish, secure the final popped BCFs (+indexes)
# to gs://<bucket>/deliverables/batch-<BATCH>/, VERIFY them, and only then delete
# that batch's Cromwell execution intermediates.
#
# Correctness guards (addressing review of the first version):
#   #1 batch identity is compared as an integer, never by unanchored substring —
#      batch 3 never matches batch-30/300, batch 1 never matches batch-100.
#   #2 index copies MUST succeed, and verification requires 250 samples AND
#      records>0 AND exactly the expected single contig AND an identical sample
#      set across all 22 chromosomes — a header-only or wrong-contig BCF fails.
#   #3 Cromwell queries are fully paginated and cross-checked against
#      totalResultsCount, so no workflow is silently excluded.
#   #4 idempotent/resumable: already-valid deliverables are reused (not recopied
#      from possibly-deleted execution dirs); any delete failure is propagated
#      (the script never prints DONE on a partial cleanup).
set -euo pipefail

BATCH_RAW="${1:?usage: prune_batch.sh <BATCH>  e.g. 003}"
[[ "$BATCH_RAW" =~ ^[0-9]+$ ]] || { echo "ERROR: batch must be numeric, got '$BATCH_RAW'"; exit 2; }
BATCH_N=$((10#$BATCH_RAW))                 # integer form for exact comparisons
BATCH=$(printf '%03d' "$BATCH_N")          # canonical zero-padded label
CROMWELL="${CROMWELL_URL:-http://localhost:8000}"
BUCKET="${PANEL_BUCKET:-gs://longreadsphase2imputation}"
EXPECT_SAMPLES="${EXPECT_SAMPLES:-250}"
EXEC_ROOT="$BUCKET/workflows/cromwell-executions/GLIMPSE2FromPreprocessedPLsJoint"
DELIV="$BUCKET/deliverables/batch-$BATCH"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
log(){ echo "[prune batch-$BATCH] $*"; }
die(){ echo "[prune batch-$BATCH] ERROR: $*" >&2; exit 1; }

# ---- #3 paginated Cromwell query: emit id<TAB>status<TAB>chr for the batch ----
query_workflows() {
  python3 - "$CROMWELL" "$BATCH" <<'PY'
import json,sys,urllib.request
base,batch=sys.argv[1],sys.argv[2]
rows=[]; page=1; size=100; total=None
while True:
    url=f"{base}/api/workflows/v1/query?label=batch:batch-{batch}&pageSize={size}&page={page}&additionalQueryResultFields=labels"
    try:
        d=json.load(urllib.request.urlopen(url,timeout=30))
    except Exception as e:
        sys.stderr.write(f"query failed page {page}: {e}\n"); sys.exit(3)
    r=d.get('results',[]); total=d.get('totalResultsCount',total)
    for w in r:
        rows.append((w['id'],w['status'],(w.get('labels') or {}).get('chromosome','?')))
    if len(r)<size: break
    page+=1
# #3 fail closed if pagination did not cover the reported total
if isinstance(total,int) and len(rows)!=total:
    sys.stderr.write(f"pagination mismatch: fetched {len(rows)} != total {total}\n"); sys.exit(3)
for i,s,c in rows: print(f"{i}\t{s}\t{c}")
PY
}

# ---- verify one deliverable BCF: samples, records, single expected contig ----
verify_bcf() {  # $1=gs bcf  $2=expected chrom (e.g. chr7)  -> prints sample-set md5, or fails
  local f="$1" chrom="$2"
  gcloud storage ls "$f.csi" >/dev/null 2>&1 || { echo "no .csi index" >&2; return 1; }   # #2
  local samples records contigs
  samples=$(bcftools query -l "$f" 2>/dev/null | wc -l | tr -d ' ')
  [ "$samples" = "$EXPECT_SAMPLES" ] || { echo "samples=$samples" >&2; return 1; }
  records=$(bcftools index -n "$f" 2>/dev/null || echo 0)
  [ "${records:-0}" -gt 0 ] || { echo "records=0 (header-only)" >&2; return 1; }            # #2
  contigs=$(bcftools index --stats "$f" 2>/dev/null | awk '$3>0{print $1}' | sort -u)
  [ "$contigs" = "$chrom" ] || { echo "contigs='$contigs' expected '$chrom'" >&2; return 1; } # #2
  bcftools query -l "$f" 2>/dev/null | sort | md5
}

log "checking Cromwell for this batch (paginated)"
query_workflows > "$WORK/wf.tsv" || die "Cromwell query failed"
nwf=$(wc -l < "$WORK/wf.tsv" | tr -d ' ')
[ "$nwf" -gt 0 ] || { log "no workflows labelled batch:batch-$BATCH — nothing to prune"; exit 0; }

# refuse unless every workflow is terminal
nonterm=$(awk -F'\t' '$2!="Succeeded" && $2!="Failed" && $2!="Aborted"' "$WORK/wf.tsv" | wc -l | tr -d ' ')
[ "$nonterm" = "0" ] || die "$nonterm workflow(s) still active — batch not complete"

# ---- #4 idempotency: are deliverables already complete & valid? ----
# Plain indexed array (bash 3.2 compatible; `declare -A` needs bash 4+, which the
# deploy environment does not guarantee). One md5 per verified chromosome; the loop
# runs chr1..22 exactly once each, so 22 entries means 22 distinct chromosomes.
SETMD5=()
secure_and_verify_chr() {  # $1=chrom
  local chrom="$1" existing bcf wid md5
  existing=$(gcloud storage ls "$DELIV/*.batch-$BATCH.$chrom.*popped.bcf" 2>/dev/null | head -1 || true)
  if [ -n "$existing" ] && md5=$(verify_bcf "$existing" "$chrom" 2>/dev/null); then
    SETMD5+=("$md5"); return 0                            # already secured & valid -> reuse (#4)
  fi
  # otherwise fetch the popped BCF path from a Succeeded workflow's outputs
  wid=$(awk -F'\t' -v c="$chrom" '$3==c && $2=="Succeeded"{print $1; exit}' "$WORK/wf.tsv")
  [ -n "$wid" ] || die "no Succeeded workflow for $chrom and no valid existing deliverable"
  bcf=$(python3 - "$CROMWELL" "$wid" <<'PY'
import json,sys,urllib.request
d=json.load(urllib.request.urlopen(f"{sys.argv[1]}/api/workflows/v1/{sys.argv[2]}/metadata?expandSubWorkflows=false",timeout=30))
o=d.get('outputs',{})
print(next((v for v in o.values() if isinstance(v,str) and v.endswith('popped.bcf')),''))
PY
)
  [ -n "$bcf" ] || die "no popped BCF output for $chrom ($wid)"
  gcloud storage cp "$bcf"     "$DELIV/" --quiet || die "copy failed: $bcf"
  gcloud storage cp "$bcf.csi" "$DELIV/" --quiet || die "index copy failed: $bcf.csi"  # #2 no '|| true'
  md5=$(verify_bcf "$bcf" "$chrom") || die "verification failed for $chrom"
  SETMD5+=("$md5")
}

log "securing + verifying 22 deliverables"
for c in $(seq 1 22); do secure_and_verify_chr "chr$c"; done

# ---- #2 cross-check: all 22 chromosomes share one identical sample set ----
uniq_sets=$(printf '%s\n' "${SETMD5[@]}" | sort -u | wc -l | tr -d ' ')
[ "${#SETMD5[@]}" = "22" ] || die "only ${#SETMD5[@]}/22 chromosomes verified"
[ "$uniq_sets" = "1" ] || die "sample sets differ across chromosomes ($uniq_sets distinct) — refusing to delete"
log "verified: 22/22 @ $EXPECT_SAMPLES samples, records>0, single-contig, one consistent sample set ✓"

# ---- delete: exact UUID dirs from Cromwell; propagate any failure (#4) ----
log "deleting execution intermediates"
fails=0
awk -F'\t' '{print $1}' "$WORK/wf.tsv" | sort -u | while read -r u; do
  [ -n "$u" ] || continue
  if gcloud storage rm -r "$EXEC_ROOT/$u" --quiet 2>/dev/null; then log "  del $u"
  else log "  FAILED to delete $u"; echo x >> "$WORK/fails"; fi
done
# #1 managed-backend outputs: match the batch number EXACTLY (integer), never substring
gcloud storage ls "$BUCKET/glimpse2_phase/" 2>/dev/null | while read -r p; do
  b=$(printf '%s\n' "$p" | grep -oE 'batch[_-]?0*[0-9]+' | grep -oE '[0-9]+$' | head -1 || true)
  [ -n "$b" ] && [ "$((10#$b))" -eq "$BATCH_N" ] || continue
  if gcloud storage rm -r "$p" --quiet 2>/dev/null; then log "  del managed $p"
  else log "  FAILED to delete $p"; echo x >> "$WORK/fails"; fi
done
[ -s "$WORK/fails" ] && die "$(wc -l < "$WORK/fails" | tr -d ' ') deletion(s) failed — cleanup incomplete, NOT reporting done"

log "DONE — batch-$BATCH pruned; 22 verified deliverables kept in $DELIV"
