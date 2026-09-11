#!/usr/bin/env bash
# prune_completed_batches.sh
#
# Driver for prune_batch.sh: finds every batch whose workflows are ALL terminal and
# prunes each. Idempotent and safe to run repeatedly; skips batches that are still
# running or already pruned.
#
# Backend (PRUNE_BACKEND): "managed" (default) enumerates VWB cloud-Cromwell jobs via
# `wb workflow job list` (display names bNNN-chrK); "local" queries the deprecated
# laptop Cromwell with the fully paginated, count-cross-checked query (#3).
set -euo pipefail
CROMWELL="${CROMWELL_URL:-http://localhost:8000}"
BACKEND="${PRUNE_BACKEND:-managed}"
HERE="$(cd "$(dirname "$0")" && pwd)"
# `wb` is a Java CLI; cron's environment has neither JAVA_HOME nor the brew JDK on PATH.
if [ -d /opt/homebrew/opt/openjdk/libexec/openjdk.jdk/Contents/Home ]; then
  export JAVA_HOME=/opt/homebrew/opt/openjdk/libexec/openjdk.jdk/Contents/Home
  export PATH=/opt/homebrew/opt/openjdk/bin:$PATH
fi
LIST="$(mktemp)"; WBJSON="$(mktemp)"; trap 'rm -f "$LIST" "$WBJSON"' EXIT

# emit "<batchnum> READY|RUNNING <count>" per batch
if [ "$BACKEND" = managed ]; then
  # (JSON via a file: `python3 -` already uses stdin for the program text)
  wb workflow job list --limit=1000 --format=JSON > "$WBJSON" 2>/dev/null || true
  python3 - "$WBJSON" <<'PY' > "$LIST"
import json,sys,re,collections
try: jobs=json.load(open(sys.argv[1]))
except Exception as e: sys.stderr.write(f"wb job list unparsable: {e}\n"); sys.exit(3)
if len(jobs)>=1000: sys.stderr.write("wb job list hit its limit; refusing (possible truncation)\n"); sys.exit(3)
TERM={'COMPLETED','FAILED','CANCELLED'}
st=collections.defaultdict(list); pat=re.compile(r'^b(\d+)-chr\d+$')
for j in jobs:
    m=pat.match(j.get('displayName') or '')
    if m: st[int(m.group(1))].append(j['status'])
for b in sorted(st):
    print(f"{b:03d} {'READY' if all(s in TERM for s in st[b]) else 'RUNNING'} {len(st[b])}")
PY
else
  python3 - "$CROMWELL" <<'PY' > "$LIST"
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
fi

[ -s "$LIST" ] || { echo "no batches found on $BACKEND backend"; exit 0; }
while read -r batch state n; do
  if [ "$state" = "READY" ]; then
    echo "=== pruning batch-$batch ($n workflows, all terminal) ==="
    "$HERE/prune_batch.sh" "$batch" || echo "  (batch-$batch prune failed — left intact, will retry next run)"
  else
    echo "=== skip batch-$batch ($n workflows, still running) ==="
  fi
done < "$LIST"
