#!/usr/bin/env bash
# prune_completed_batches.sh
#
# Driver for prune_batch.sh: finds every batch Cromwell knows about whose
# workflows are ALL terminal and prunes each. Idempotent and safe to run
# repeatedly; skips batches that are still running or already pruned.
#
# #3 The batch enumeration is fully paginated and cross-checked against
#    totalResultsCount, so no batch is dropped once the run exceeds one page.
set -euo pipefail
CROMWELL="${CROMWELL_URL:-http://localhost:8000}"
HERE="$(cd "$(dirname "$0")" && pwd)"

# paginated: emit "<batchnum> READY|RUNNING <count>" per batch
python3 - "$CROMWELL" <<'PY' > /tmp/prune_driver_batches.txt
import json,sys,urllib.request,collections
base=sys.argv[1]; rows=[]; page=1; size=100; total=None
while True:
    url=f"{base}/api/workflows/v1/query?pageSize={size}&page={page}&additionalQueryResultFields=labels"
    try:
        d=json.load(urllib.request.urlopen(url,timeout=30))
    except Exception as e:
        sys.stderr.write(f"query failed page {page}: {e}\n"); sys.exit(3)
    r=d.get('results',[]); total=d.get('totalResultsCount',total)
    rows.extend(r)
    if len(r)<size: break
    page+=1
if isinstance(total,int) and len(rows)!=total:
    sys.stderr.write(f"pagination mismatch: fetched {len(rows)} != total {total}\n"); sys.exit(3)
TERM={'Succeeded','Failed','Aborted'}
st=collections.defaultdict(list)
for w in rows:
    b=(w.get('labels') or {}).get('batch')
    if b and b.startswith('batch-'): st[b].append(w['status'])
for b in sorted(st):
    ready=all(s in TERM for s in st[b])
    print(f"{b.split('-')[-1]} {'READY' if ready else 'RUNNING'} {len(st[b])}")
PY

[ -s /tmp/prune_driver_batches.txt ] || { echo "no batches found in Cromwell"; exit 0; }
while read -r batch state n; do
  if [ "$state" = "READY" ]; then
    echo "=== pruning batch-$batch ($n workflows, all terminal) ==="
    "$HERE/prune_batch.sh" "$batch" || echo "  (batch-$batch prune failed — left intact, will retry next run)"
  else
    echo "=== skip batch-$batch ($n workflows, still running) ==="
  fi
done < /tmp/prune_driver_batches.txt
