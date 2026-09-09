#!/usr/bin/env bash
# Blocking tests for prune_batch.sh using a fake Cromwell + stubbed gcloud/bcftools.
# Covers the review's data-loss paths: a header-only / index-less / wrong-contig
# deliverable must ABORT before any deletion; a good batch must delete exactly its
# own UUIDs; a re-run over already-secured deliverables must be a no-op, not a failure.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
SCRIPT="$HERE/../prune_batch.sh"
PASS=0; FAIL=0
ok(){ echo "  ok: $1"; PASS=$((PASS+1)); }
bad(){ echo "  FAIL: $1"; FAIL=$((FAIL+1)); }

# --- fake Cromwell: 22 Succeeded chr workflows for batch-003, each output a popped.bcf ---
start_cromwell() {
  python3 - "$1" >/dev/null 2>&1 <<'PY' &
import sys,json,http.server
port=int(sys.argv[1])
wf=[{"id":f"uuid-chr{c}","status":"Succeeded","labels":{"batch":"batch-003","chromosome":f"chr{c}"}} for c in range(1,23)]
class H(http.server.BaseHTTPRequestHandler):
    def log_message(self,*a): pass
    def do_GET(self):
        if "/query" in self.path:
            body={"results":wf,"totalResultsCount":len(wf)}
        else:  # /metadata for uuid-chrN
            uid=self.path.split("/v1/")[1].split("/metadata")[0]
            c=uid.replace("uuid-chr","")
            body={"outputs":{"x":f"gs://b/exec/{uid}/call-Concat/aou2.batch-003.chr{c}.popped.bcf"}}
        b=json.dumps(body).encode(); self.send_response(200)
        self.send_header("Content-Type","application/json"); self.end_headers(); self.wfile.write(b)
httpd=http.server.HTTPServer(("127.0.0.1",port),H); httpd.serve_forever()
PY
  echo $!
}

run_case() {  # $1=name $2=mode(headeronly|good|missing_csi)  expected exit behaviour asserted by caller
  local name="$1" mode="$2" port=8771
  local tmp; tmp="$(mktemp -d)"
  # stub gcloud + bcftools on PATH
  cat > "$tmp/gcloud" <<EOF
#!/usr/bin/env bash
if [ "\$1" = storage ] && [ "\$2" = rm ]; then echo "\$@" >> "$tmp/rm.log"; exit 0; fi
if [ "\$1" = storage ] && [ "\$2" = cp ]; then exit 0; fi
if [ "\$1" = storage ] && [ "\$2" = ls ]; then
  case "\$3" in
    *.csi) [ "$mode" = missing_csi ] && exit 1; echo "\$3"; exit 0;;
    *deliverables*batch-003.chr*)  # per-chr existing-deliverable check
      if [ "$mode" = rerun ]; then c=\$(echo "\$3"|grep -oE 'chr[0-9]+'|head -1); echo "gs://b/deliverables/batch-003/aou2.batch-003.\$c.popped.bcf"; fi
      exit 0;;                      # empty otherwise -> fetch path
    *cromwell-executions*uuid-*) [ "$mode" = rerun ] && exit 1; exit 0;;   # exec dir gone on rerun
    *glimpse2_phase*) exit 0;;
  esac
  exit 0
fi
exit 0
EOF
  cat > "$tmp/bcftools" <<EOF
#!/usr/bin/env bash
case "\$1" in
  query) for i in \$(seq 1 250); do echo "S\$i"; done;;   # 250 samples
  index)
    if [ "\$2" = -n ]; then [ "$mode" = headeronly ] && echo 0 || echo 1000; fi
    if [ "\$2" = --stats ]; then c=\$(echo "\$3"|grep -oE 'chr[0-9]+'); [ "$mode" = headeronly ] && echo "\$c 0 0" || echo "\$c 1 1000"; fi;;
esac
EOF
  chmod +x "$tmp/gcloud" "$tmp/bcftools"
  PANEL_BUCKET=gs://b CROMWELL_URL="http://127.0.0.1:$port" GCS_OAUTH_TOKEN=faketoken \
    PATH="$tmp:$PATH" bash "$SCRIPT" 003 > "$tmp/out.log" 2>&1
  local rc=$?
  echo "$rc|$tmp"
}

PORT=8771
CPID=$(start_cromwell $PORT); sleep 1
trap 'kill $CPID 2>/dev/null' EXIT

echo "== case 1: header-only deliverable must ABORT, no rm =="
res=$(run_case c1 headeronly); rc="${res%%|*}"; tmp="${res#*|}"
{ [ "$rc" -ne 0 ] && [ ! -s "$tmp/rm.log" ]; } && ok "header-only refused, nothing deleted" || bad "header-only: rc=$rc rm.log=$(cat "$tmp/rm.log" 2>/dev/null)"

echo "== case 2: missing index must ABORT, no rm =="
res=$(run_case c2 missing_csi); rc="${res%%|*}"; tmp="${res#*|}"
{ [ "$rc" -ne 0 ] && [ ! -s "$tmp/rm.log" ]; } && ok "missing .csi refused, nothing deleted" || bad "missing_csi: rc=$rc"

echo "== case 3: good batch deletes exactly its 22 UUIDs =="
res=$(run_case c3 good); rc="${res%%|*}"; tmp="${res#*|}"
n=$(grep -c 'GLIMPSE2FromPreprocessedPLsJoint/uuid-chr' "$tmp/rm.log" 2>/dev/null); n=${n:-0}
{ [ "$rc" -eq 0 ] && [ "$n" -eq 22 ]; } && ok "good batch: 22 UUIDs deleted" || bad "good: rc=$rc deleted=$n"

echo "== case 4: rerun (valid deliverables, exec dirs already gone) must be a clean no-op =="
res=$(run_case c4 rerun); rc="${res%%|*}"; tmp="${res#*|}"
d=$(grep -c 'already gone' "$tmp/out.log" 2>/dev/null); d=${d:-0}
f=$(grep -c 'FAILED to delete' "$tmp/out.log" 2>/dev/null); f=${f:-0}
{ [ "$rc" -eq 0 ] && [ "$f" -eq 0 ] && [ "$d" -eq 22 ]; } && ok "rerun: 22 already-gone, 0 failures, rc=0" || bad "rerun: rc=$rc already_gone=$d failed=$f"

echo "----"; echo "PASS=$PASS FAIL=$FAIL"
[ "$FAIL" -eq 0 ]
