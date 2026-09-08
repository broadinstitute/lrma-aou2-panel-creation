#!/usr/bin/env bash
# prune_completed_batches.sh
#
# Driver for prune_batch.sh: finds every batch Cromwell knows about whose
# workflows are ALL terminal, and prunes each (secure->verify->delete).
# Idempotent and safe to run on a schedule or before each new batch launch —
# it silently skips batches that are still running or already pruned.
set -euo pipefail
CROMWELL="${CROMWELL_URL:-http://localhost:8000}"
HERE="$(cd "$(dirname "$0")" && pwd)"

# distinct batch numbers currently in Cromwell, and whether each is fully terminal
curl -s --max-time 30 "$CROMWELL/api/workflows/v1/query?pageSize=1000&additionalQueryResultFields=labels" \
 | python3 -c "
import json,sys,collections
r=json.load(sys.stdin).get('results',[])
st=collections.defaultdict(list)
for w in r:
    b=(w.get('labels') or {}).get('batch')
    if b and b.startswith('batch-'): st[b].append(w['status'])
TERM={'Succeeded','Failed','Aborted'}
for b in sorted(st):
    done=all(s in TERM for s in st[b])
    print(f\"{b.split('-')[-1]} {'READY' if done else 'RUNNING'} {len(st[b])}\")
" | while read -r batch state n; do
  if [ "$state" = "READY" ]; then
    echo "=== pruning batch-$batch ($n workflows, all terminal) ==="
    "$HERE/prune_batch.sh" "$batch" || echo "  (batch-$batch prune skipped/failed — left intact)"
  else
    echo "=== skip batch-$batch ($n workflows, still running) ==="
  fi
done
