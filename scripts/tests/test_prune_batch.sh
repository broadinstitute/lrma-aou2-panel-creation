#!/usr/bin/env bash
# Blocking tests for prune_batch.sh (+ driver) with stubbed gcloud/bcftools/wb and a fake
# local Cromwell. Covers the data-loss paths: a header-only / index-less deliverable must
# ABORT before any deletion; a good batch must delete exactly its own tree; a re-run over
# already-secured deliverables must be a no-op; and on the managed backend a FAILED or
# still-RUNNING chromosome, or an ambiguous popped BCF, must refuse to delete anything.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
SCRIPT="$HERE/../prune_batch.sh"
DRIVER="$HERE/../prune_completed_batches.sh"
PASS=0; FAIL=0
ok(){ echo "  ok: $1"; PASS=$((PASS+1)); }
bad(){ echo "  FAIL: $1"; FAIL=$((FAIL+1)); }

# --- fake local Cromwell: 22 Succeeded chr workflows for batch-003, each output a popped.bcf ---
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

# --- fake `wb workflow job list` JSON for batch 003 (+ a batch-030 decoy and a cancelled dup) ---
write_wb_json() {  # $1=file $2=mode
  python3 - "$1" "$2" <<'PY'
import json,sys
out,mode=sys.argv[1],sys.argv[2]
def job(name,st,rid):
    root="gs://b/glimpse2_phase/batch003/"+name.split('-')[1]+"/"+name
    return {"runId":rid,"displayName":name,"status":st,"engineAttributes":{"gcpCromwell":{"jesGcsRoot":root}}}
jobs=[job(f"b003-chr{c}","COMPLETED",f"run-chr{c}") for c in range(1,23)]
jobs.append(job("b003-chr1","CANCELLED","run-chr1-dup"))          # cancelled duplicate must be ignored
jobs.append(job("b030-chr1","RUNNING","run-b030"))                 # other batch must not leak in
if mode=="mfailed":  jobs[4]["status"]="FAILED"
if mode=="mrunning": jobs[4]["status"]="RUNNING"
json.dump(jobs,open(out,"w"))
PY
}

run_case() {  # $1=name $2=mode  -> "rc|tmpdir"
  local name="$1" mode="$2" port=8771 backend=local
  case "$mode" in m*) backend=managed;; esac
  local tmp; tmp="$(mktemp -d)"
  write_wb_json "$tmp/wb.json" "$mode"
  cat > "$tmp/wb" <<EOF
#!/usr/bin/env bash
cat "$tmp/wb.json"
EOF
  # stub gcloud + bcftools on PATH
  cat > "$tmp/gcloud" <<EOF
#!/usr/bin/env bash
if [ "\$1" = storage ] && [ "\$2" = rm ]; then echo "\$@" >> "$tmp/rm.log"; exit 0; fi
if [ "\$1" = storage ] && [ "\$2" = cp ]; then exit 0; fi
if [ "\$1" = storage ] && [ "\$2" = ls ]; then
  case "\$3" in
    *.csi) [ "$mode" = missing_csi ] && exit 1; echo "\$3"; exit 0;;
    *deliverables*batch-003.chr*)  # per-chr existing-deliverable check
      if [ "$mode" = rerun ] || [ "$mode" = mrerun ]; then c=\$(echo "\$3"|grep -oE 'chr[0-9]+'|head -1); echo "gs://b/deliverables/batch-003/aou2.batch-003.\$c.popped.bcf"; fi
      exit 0;;                      # empty otherwise -> fetch path
    *cromwell-executions*uuid-*) [ "$mode" = rerun ] && exit 1; exit 0;;   # exec dir gone on rerun
    */GLIMPSE2FromPreprocessedPLsJoint/run-chr*/'**'/*)   # managed popped-BCF glob: shard decoys + the real one
      r=\$(echo "\$3" | sed 's|/\*\*/.*||'); c=\$(echo "\$3" | grep -oE 'run-chr[0-9]+' | sed 's/run-//')
      f="aou2.batch-003.\$c.glimpse2.popped.bcf"
      echo "\$r/call-PopAndMarginalizeCollisions/shard-0/\$f"
      echo "\$r/call-PopAndMarginalizeCollisions/shard-1/\$f"
      echo "\$r/call-ConcatPopAndMarginalizeCollisions/ConcatVcfs/sub-1/call-ConcatVcfs/\$f"
      [ "$mode" = mambig ] && echo "\$r/call-ConcatPopAndMarginalizeCollisions/ConcatVcfs/sub-2/call-ConcatVcfs/\$f"
      exit 0;;
    gs://b/glimpse2_phase/)          # top level: this batch's tree + a batch-030 decoy (exact-integer match)
      [ "$mode" = mrerun ] || echo "gs://b/glimpse2_phase/batch003/"
      echo "gs://b/glimpse2_phase/batch030/"; exit 0;;
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
  chmod +x "$tmp/gcloud" "$tmp/bcftools" "$tmp/wb"
  PANEL_BUCKET=gs://b CROMWELL_URL="http://127.0.0.1:$port" GCS_OAUTH_TOKEN=faketoken PRUNE_BACKEND=$backend \
    PATH="$tmp:$PATH" bash "$SCRIPT" 003 > "$tmp/out.log" 2>&1
  local rc=$?
  echo "$rc|$tmp"
}

PORT=8771
CPID=$(start_cromwell $PORT); sleep 1
trap 'kill $CPID 2>/dev/null' EXIT

echo "== local 1: header-only deliverable must ABORT, no rm =="
res=$(run_case c1 headeronly); rc="${res%%|*}"; tmp="${res#*|}"
{ [ "$rc" -ne 0 ] && [ ! -s "$tmp/rm.log" ]; } && ok "header-only refused, nothing deleted" || bad "header-only: rc=$rc rm.log=$(cat "$tmp/rm.log" 2>/dev/null)"

echo "== local 2: missing index must ABORT, no rm =="
res=$(run_case c2 missing_csi); rc="${res%%|*}"; tmp="${res#*|}"
{ [ "$rc" -ne 0 ] && [ ! -s "$tmp/rm.log" ]; } && ok "missing .csi refused, nothing deleted" || bad "missing_csi: rc=$rc"

echo "== local 3: good batch deletes exactly its 22 UUIDs =="
res=$(run_case c3 good); rc="${res%%|*}"; tmp="${res#*|}"
n=$(grep -c 'GLIMPSE2FromPreprocessedPLsJoint/uuid-chr' "$tmp/rm.log" 2>/dev/null); n=${n:-0}
{ [ "$rc" -eq 0 ] && [ "$n" -eq 22 ]; } && ok "good batch: 22 UUIDs deleted" || bad "good: rc=$rc deleted=$n"

echo "== local 4: rerun (valid deliverables, exec dirs already gone) must be a clean no-op =="
res=$(run_case c4 rerun); rc="${res%%|*}"; tmp="${res#*|}"
d=$(grep -c 'already gone' "$tmp/out.log" 2>/dev/null); d=${d:-0}
f=$(grep -c 'FAILED to delete' "$tmp/out.log" 2>/dev/null); f=${f:-0}
{ [ "$rc" -eq 0 ] && [ "$f" -eq 0 ] && [ "$d" -eq 22 ]; } && ok "rerun: 22 already-gone, 0 failures, rc=0" || bad "rerun: rc=$rc already_gone=$d failed=$f"

echo "== managed 5: 22 COMPLETED -> secure from ConcatVcfs dir (not shards), delete batch003 tree only =="
res=$(run_case m5 mgood); rc="${res%%|*}"; tmp="${res#*|}"
t=$(grep -c 'glimpse2_phase/batch003/' "$tmp/rm.log" 2>/dev/null); t=${t:-0}
x=$(grep -c 'batch030' "$tmp/rm.log" 2>/dev/null); x=${x:-0}
u=$(grep -c 'cromwell-executions' "$tmp/rm.log" 2>/dev/null); u=${u:-0}
{ [ "$rc" -eq 0 ] && [ "$t" -eq 1 ] && [ "$x" -eq 0 ] && [ "$u" -eq 0 ] && grep -q 'DONE' "$tmp/out.log"; } \
  && ok "managed good: batch003 tree deleted once, batch030 untouched, no exec-root deletes" \
  || bad "managed good: rc=$rc tree=$t decoy=$x exec=$u  $(tail -2 "$tmp/out.log")"

echo "== managed 6: one chromosome FAILED (no COMPLETED run) must ABORT, no rm =="
res=$(run_case m6 mfailed); rc="${res%%|*}"; tmp="${res#*|}"
{ [ "$rc" -ne 0 ] && [ ! -s "$tmp/rm.log" ]; } && ok "managed failed chr refused, nothing deleted" || bad "managed failed: rc=$rc rm=$(cat "$tmp/rm.log" 2>/dev/null)"

echo "== managed 7: one chromosome still RUNNING must ABORT, no rm =="
res=$(run_case m7 mrunning); rc="${res%%|*}"; tmp="${res#*|}"
{ [ "$rc" -ne 0 ] && [ ! -s "$tmp/rm.log" ] && grep -q 'still active' "$tmp/out.log"; } && ok "managed running chr refused, nothing deleted" || bad "managed running: rc=$rc"

echo "== managed 8: ambiguous popped BCF (2 candidates) must ABORT, no rm =="
res=$(run_case m8 mambig); rc="${res%%|*}"; tmp="${res#*|}"
{ [ "$rc" -ne 0 ] && [ ! -s "$tmp/rm.log" ] && grep -q 'expected exactly 1' "$tmp/out.log"; } && ok "ambiguous popped BCF refused" || bad "managed ambig: rc=$rc $(tail -1 "$tmp/out.log")"

echo "== managed 9: rerun (deliverables valid, tree already gone) must be a clean no-op =="
res=$(run_case m9 mrerun); rc="${res%%|*}"; tmp="${res#*|}"
{ [ "$rc" -eq 0 ] && [ ! -s "$tmp/rm.log" ] && grep -q 'DONE' "$tmp/out.log"; } && ok "managed rerun: rc=0, nothing to delete" || bad "managed rerun: rc=$rc $(tail -1 "$tmp/out.log")"

echo "== driver 10: managed discovery marks 003 READY and 030 RUNNING =="
tmp="$(mktemp -d)"; write_wb_json "$tmp/wb.json" mgood
printf '#!/usr/bin/env bash\ncat "%s/wb.json"\n' "$tmp" > "$tmp/wb"; chmod +x "$tmp/wb"
printf '#!/usr/bin/env bash\necho "stub prune $*"\n' > "$tmp/prune_batch.sh"; chmod +x "$tmp/prune_batch.sh"
cp "$DRIVER" "$tmp/driver.sh"
out=$(PRUNE_BACKEND=managed PATH="$tmp:$PATH" bash "$tmp/driver.sh" 2>&1)
{ echo "$out" | grep -q 'pruning batch-003 (23 workflows' && echo "$out" | grep -q 'skip batch-030' && echo "$out" | grep -q 'stub prune 003'; } \
  && ok "driver: 003 pruned, 030 skipped" || bad "driver: $out"

echo "----"; echo "PASS=$PASS FAIL=$FAIL"
[ "$FAIL" -eq 0 ]
