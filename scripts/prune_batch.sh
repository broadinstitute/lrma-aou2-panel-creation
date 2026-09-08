#!/usr/bin/env bash
# prune_batch.sh <BATCH>  (e.g. prune_batch.sh 003)
#
# Permanently caps per-batch storage: after a batch's 22 chromosomes finish,
# secure the final popped BCFs to gs://<bucket>/deliverables/batch-<BATCH>/,
# VERIFY each has the expected sample count, and only THEN delete that batch's
# Cromwell execution intermediates (~250-300 GB/batch of throwaway shards/logs).
#
# Safety design (this script deletes cloud data, so it is defensive by default):
#   * Refuses unless EVERY workflow for the batch is in a terminal state
#     (never touches a Running/Submitted batch).
#   * Secures + verifies all 22 deliverables BEFORE any delete; aborts if any
#     is missing or has the wrong sample count.
#   * Deletes only the exact execution UUIDs Cromwell reports for this batch
#     (plus managed-backend glimpse2_phase/*batch-<BATCH>* prefixes) — no globbing
#     of other batches, no "delete everything except" logic.
#   * Idempotent: re-running after a prune is a no-op.
set -euo pipefail

BATCH="${1:?usage: prune_batch.sh <BATCH>  e.g. 003}"
CROMWELL="${CROMWELL_URL:-http://localhost:8000}"
BUCKET="${PANEL_BUCKET:-gs://longreadsphase2imputation}"
EXPECT_SAMPLES="${EXPECT_SAMPLES:-250}"
EXEC_ROOT="$BUCKET/workflows/cromwell-executions/GLIMPSE2FromPreprocessedPLsJoint"
DELIV="$BUCKET/deliverables/batch-$BATCH"
CHRS=$(seq 1 22)
log(){ echo "[prune batch-$BATCH] $*"; }

# --- 1. gather this batch's workflows from Cromwell (authoritative) ---
meta=$(curl -s --max-time 30 "$CROMWELL/api/workflows/v1/query?label=batch:batch-$BATCH&pageSize=200&additionalQueryResultFields=labels") \
  || { log "ERROR: cannot reach Cromwell at $CROMWELL"; exit 1; }
python3 - "$meta" <<'PY' > /tmp/prune_$BATCH.json
import json,sys
r=json.loads(sys.argv[1]).get('results',[])
out={'workflows':[{'id':w['id'],'status':w['status'],
       'chr':(w.get('labels') or {}).get('chromosome','?')} for w in r]}
json.dump(out,open(1,'w'))
PY
nwf=$(python3 -c "import json;print(len(json.load(open('/tmp/prune_$BATCH.json'))['workflows']))")
[ "$nwf" -gt 0 ] || { log "no workflows labelled batch:batch-$BATCH in Cromwell — nothing to prune"; exit 0; }

# --- 2. refuse if any workflow is non-terminal ---
nonterm=$(python3 -c "
import json
w=json.load(open('/tmp/prune_$BATCH.json'))['workflows']
print(sum(1 for x in w if x['status'] not in ('Succeeded','Failed','Aborted')))")
[ "$nonterm" = "0" ] || { log "ABORT: $nonterm workflow(s) still active — batch not complete"; exit 1; }

# --- 3. secure + verify deliverables (one popped BCF per chromosome) ---
log "securing deliverables -> $DELIV"
for c in $CHRS; do
  chrom="chr$c"
  # authoritative output path from a Succeeded workflow's metadata
  wid=$(python3 -c "
import json
w=json.load(open('/tmp/prune_$BATCH.json'))['workflows']
print(next((x['id'] for x in w if x['chr']=='$chrom' and x['status']=='Succeeded'),''))")
  [ -n "$wid" ] || { log "ABORT: no Succeeded workflow for $chrom"; exit 1; }
  bcf=$(curl -s --max-time 30 "$CROMWELL/api/workflows/v1/$wid/metadata?expandSubWorkflows=false" \
        | python3 -c "import json,sys;o=json.load(sys.stdin).get('outputs',{});print(next((v for v in o.values() if isinstance(v,str) and v.endswith('popped.bcf')),''))")
  [ -n "$bcf" ] || { log "ABORT: no popped BCF output for $chrom ($wid)"; exit 1; }
  gcloud storage cp "$bcf" "$DELIV/" --quiet
  gcloud storage cp "$bcf.csi" "$DELIV/" --quiet 2>/dev/null || true
done

log "verifying deliverables ($EXPECT_SAMPLES samples each)"
bad=0
for f in $(gcloud storage ls "$DELIV/*.popped.bcf" 2>/dev/null); do
  n=$(bcftools query -l "$f" 2>/dev/null | wc -l | tr -d ' ')
  [ "$n" = "$EXPECT_SAMPLES" ] || { log "  BAD: $(basename "$f") has $n samples"; bad=$((bad+1)); }
done
ndeliv=$(gcloud storage ls "$DELIV/*.popped.bcf" 2>/dev/null | grep -oE 'chr[0-9]+' | sort -u | wc -l | tr -d ' ')
{ [ "$ndeliv" = "22" ] && [ "$bad" = "0" ]; } || { log "ABORT: deliverables incomplete/invalid ($ndeliv/22 present, $bad bad) — NOT deleting"; exit 1; }
log "deliverables verified: 22/22 @ $EXPECT_SAMPLES samples ✓"

# --- 4. delete this batch's execution UUIDs (exact list) + managed prefixes ---
log "deleting execution intermediates"
python3 -c "
import json
w=json.load(open('/tmp/prune_$BATCH.json'))['workflows']
print('\n'.join(sorted(set(x['id'] for x in w))))" > /tmp/prune_uuids_$BATCH.txt
while read -r u; do
  [ -n "$u" ] && gcloud storage rm -r "$EXEC_ROOT/$u" --quiet 2>/dev/null && log "  del $u"
done < /tmp/prune_uuids_$BATCH.txt
# managed-backend outputs for this batch, if any
for p in $(gcloud storage ls "$BUCKET/glimpse2_phase/" 2>/dev/null | grep -iE "batch[_-]?0*$BATCH|batch$BATCH"); do
  gcloud storage rm -r "$p" --quiet 2>/dev/null && log "  del managed $p"
done

log "DONE — batch-$BATCH pruned; 22 deliverables kept in $DELIV"
