#!/usr/bin/env bash
# Runs the real command block of wdl/methods/imputation/SecureAndPruneBatch.wdl against stubbed
# gcloud / bcftools / md5sum / date, for the success path and every refusal path. Asserts that
# nothing is deleted and no PRUNED.ok is written unless all guards pass.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
WDL="$HERE/../../wdl/methods/imputation/SecureAndPruneBatch.wdl"
PASS=0; FAIL=0
ok(){ echo "  ok: $1"; PASS=$((PASS+1)); }
bad(){ echo "  FAIL: $1"; FAIL=$((FAIL+1)); }

STUBS="$(mktemp -d)"
cat > "$STUBS/gcloud" <<'EOF'
#!/usr/bin/env bash
set -u
echo "$*" >> "$STATE/calls.log"
nomatch() { echo "ERROR: (gcloud.storage.ls) One or more URLs matched no objects." >&2; exit 1; }
R=gs://b/glimpse2_phase/batch005; D=gs://b/deliverables/batch-005
if [ "$1" = storage ] && [ "$2" = ls ]; then
  shift 2; L=0; if [ "$1" = -l ]; then L=1; shift; fi; url="$1"; TS="2026-09-11T21:27:51Z"
  emit() { if [ $L = 1 ]; then printf '  %s  %s  %s\n' 2000000000 "$TS" "$1"; else printf '%s\n' "$1"; fi; }
  case "$url" in
    "$R/**")
      [ "$MODE" = listerr ] && { echo "ERROR: (gcloud.storage.ls) HTTPError 503: backend unavailable" >&2; exit 1; }
      [ "$MODE" = rerun ] && nomatch
      if [ -f "$STATE/deleted" ] && [ "$MODE" != stillthere ]; then nomatch; fi
      [ "$MODE" = recent ] && TS=$(python3 -c 'import datetime as d;print((d.datetime.now(d.timezone.utc)-d.timedelta(minutes=30)).strftime("%Y-%m-%dT%H:%M:%SZ"))')
      for k in $(seq 1 22); do
        if [ "$MODE" = missingchr ] && [ "$k" = 7 ]; then continue; fi
        base="aou2.batch-005.chr$k.prod-x.glimpse2.popped.bcf"; run="$R/b005/GLIMPSE2FromPreprocessedPLsJoint/run$k"
        emit "$run/call-ConcatPopAndMarginalizeCollisions/ConcatVcfs/sub$k/call-ConcatVcfs/$base"
        emit "$run/call-ConcatPopAndMarginalizeCollisions/ConcatVcfs/sub$k/call-ConcatVcfs/$base.csi"
        emit "$run/call-PopAndMarginalizeCollisions/shard-0/$base"
        emit "$run/call-GLIMPSE2Ligate/aou2.batch-005.chr$k.prod-x.glimpse2.bubble.bcf"
        if [ "$MODE" = dupsrc ] && [ "$k" = 3 ]; then emit "$run/call-ConcatPopAndMarginalizeCollisions/ConcatVcfs/subB/call-ConcatVcfs/$base"; fi
        if [ "$MODE" = duplisting ] && [ "$k" = 3 ]; then emit "$run/call-ConcatPopAndMarginalizeCollisions/ConcatVcfs/sub$k/call-ConcatVcfs/$base"; fi
      done
      [ $L = 1 ] && echo "TOTAL: 88 objects, 1 bytes (1 B)"
      exit 0;;
    "$D/*.glimpse2.popped.bcf")
      [ "$MODE" = delivlisterr ] && { echo "ERROR: (gcloud.storage.ls) HTTPError 500: internal" >&2; exit 1; }
      if [ "$MODE" = rerun ]; then for k in $(seq 1 22); do echo "$D/aou2.batch-005.chr$k.prod-x.glimpse2.popped.bcf"; done; exit 0; fi
      if [ "$MODE" = otherdeliv ]; then echo "$D/aou2.batch-005.chr2.prod-oldsha.glimpse2.popped.bcf"; exit 0; fi
      if [ "$MODE" = partialdeliv ] || [ "$MODE" = staledeliv ]; then echo "$D/aou2.batch-005.chr4.prod-x.glimpse2.popped.bcf"; exit 0; fi
      nomatch;;
    "$D/PRUNED.ok") [ -f "$STATE/marker" ] && { echo "$url"; exit 0; }; nomatch;;
    *) grep -qxF "$url" "$STATE/copied" 2>/dev/null && { echo "$url"; exit 0; }
       if [ "$MODE" = rerun ] || [ "$MODE" = partialdeliv ] || [ "$MODE" = staledeliv ]; then case "$url" in "$D"/*) echo "$url"; exit 0;; esac; fi
       nomatch;;
  esac
fi
if [ "$1" = storage ] && [ "$2" = cp ]; then
  shift 2; args=(); for a in "$@"; do [ "$a" = --quiet ] || args+=("$a"); done
  n=${#args[@]}; dest="${args[$((n-1))]}"
  if [ "${args[0]}" = - ]; then cat > "$STATE/marker"; exit 0; fi
  for ((i=0;i<n-1;i++)); do s="${args[$i]}"; b=$(basename "$s")
    case "$dest" in gs://*) echo "${dest%/}/$b" >> "$STATE/copied";; *) mkdir -p "$dest"; echo x > "$dest/$b";; esac
  done; exit 0
fi
if [ "$1" = storage ] && [ "$2" = objects ] && [ "$3" = describe ]; then
  url="$4"
  case "$url" in "$D"/*)
    [ "$MODE" = crcfail ] && { echo "ERROR: describe failed" >&2; exit 1; }
    if [ "$MODE" = crcmismatch ]; then case "$url" in *.bcf) echo "crc-different"; exit 0;; esac; fi
    if [ "$MODE" = staledeliv ] && ! grep -qxF "$url" "$STATE/copied" 2>/dev/null; then echo "crc-stale"; exit 0; fi;;
  esac
  echo "crc-$(basename "$url")"; exit 0
fi
if [ "$1" = storage ] && [ "$2" = rm ]; then touch "$STATE/deleted"; exit 0; fi
echo "unhandled gcloud call: $*" >&2; exit 9
EOF
cat > "$STUBS/bcftools" <<'EOF'
#!/usr/bin/env bash
case "$1" in
  query) n=250; [ "$MODE" = fewsamples ] && n=249; for i in $(seq 1 $n); do echo "S$i"; done;;
  index)
    if [ "$2" = -n ]; then if [ "$MODE" = headeronly ]; then echo 0; else echo 1000; fi
    elif [ "$2" = --stats ]; then c=$(basename "$3" | grep -oE 'chr[0-9]+' | head -1)
      if [ "$MODE" = headeronly ]; then echo "$c 0 0"; elif [ "$MODE" = wrongcontig ]; then echo "chr99 1 1000"; else echo "$c 1 1000"; fi; fi;;
esac
EOF
cat > "$STUBS/md5sum" <<'EOF'
#!/usr/bin/env bash
python3 -c 'import hashlib,sys;print(hashlib.md5(sys.stdin.buffer.read()).hexdigest()+"  -")'
EOF
cat > "$STUBS/date" <<'EOF'
#!/usr/bin/env bash
if [ "${2:-}" = -d ]; then
  python3 - "$3" "$4" <<'PY'
import sys,re,datetime as d
m=re.match(r'(\d+) hours ago',sys.argv[1]); t=d.datetime.now(d.timezone.utc)-d.timedelta(hours=int(m.group(1)))
print(t.strftime(sys.argv[2].lstrip('+')))
PY
else exec /bin/date "$@"; fi
EOF
chmod +x "$STUBS"/*

render() {  # $1=out script  $2=delete_run_folder
  python3 - "$WDL" "$1" "$2" <<'PY'
import sys,re
w=open(sys.argv[1]).read(); cmd=w.split('command <<<',1)[1].split('>>>',1)[0]
for k,v in {'batch':'5','bucket':'gs://b','expected_samples':'250','quiet_hours':'6','delete_run_folder':sys.argv[3]}.items():
    cmd=cmd.replace('~{'+k+'}',v)
assert '~{' not in cmd, re.findall(r'~\{[^}]*\}',cmd)
open(sys.argv[2],'w').write(cmd)
PY
}

run_case() {  # $1=mode $2=delete flag  -> sets RC, ST (state dir)
  ST="$(mktemp -d)"; : > "$ST/calls.log"; : > "$ST/copied"
  render "$ST/cmd.sh" "${2:-true}"
  ( cd "$ST" && MODE="$1" STATE="$ST" PATH="$STUBS:$PATH" bash cmd.sh > out.log 2>&1 ); RC=$?
}
nocopy(){ [ ! -s "$ST/copied" ]; }; nodelete(){ [ ! -f "$ST/deleted" ]; }; nomarker(){ [ ! -f "$ST/marker" ]; }
has(){ grep -q "$1" "$ST/out.log"; }

echo "== 1 good: secure 22, delete, confirm gone, then marker =="
run_case good
rmline=$(grep -n 'storage rm' "$ST/calls.log" | head -1 | cut -d: -f1); mkline=$(grep -n 'storage cp - gs://b/deliverables/batch-005/PRUNED.ok' "$ST/calls.log" | cut -d: -f1)
{ [ $RC = 0 ] && [ -f "$ST/deleted" ] && [ -f "$ST/marker" ] && [ "$(wc -l < "$ST/copied" | tr -d ' ')" = 44 ] && [ -n "$rmline" ] && [ -n "$mkline" ] && [ "$rmline" -lt "$mkline" ] && has 'DONE'; } \
  && ok "good: 44 objects copied, folder deleted before marker" || bad "good rc=$RC $(tail -3 "$ST/out.log")"

for m in listerr:3 delivlisterr:3; do mode=${m%%:*}; want=${m##*:}
  echo "== listing error ($mode) must abort with exit $want: no copy, no delete, no marker =="
  run_case $mode; { [ $RC = $want ] && nocopy && nodelete && nomarker; } && ok "$mode: aborted, nothing changed" || bad "$mode rc=$RC $(tail -2 "$ST/out.log")"
done

echo "== recent write within quiet window must refuse =="
run_case recent; { [ $RC = 1 ] && nocopy && nodelete && nomarker && has 'possible active run'; } && ok "recent: refused, nothing changed" || bad "recent rc=$RC $(tail -2 "$ST/out.log")"

for mode in missingchr dupsrc headeronly fewsamples wrongcontig otherdeliv; do
  echo "== $mode must refuse: no delete, no marker =="
  run_case $mode; { [ $RC = 1 ] && nodelete && nomarker; } && ok "$mode: refused" || bad "$mode rc=$RC $(tail -2 "$ST/out.log")"
done

echo "== duplicate listing of the same object is not a duplicate output =="
run_case duplisting; { [ $RC = 0 ] && [ -f "$ST/marker" ]; } && ok "duplisting: treated as one output" || bad "duplisting rc=$RC $(tail -2 "$ST/out.log")"

echo "== crc32c lookup failure after copy must abort, never compare empty hashes =="
run_case crcfail; { [ $RC != 0 ] && nodelete && nomarker; } && ok "crcfail: aborted" || bad "crcfail rc=$RC $(tail -2 "$ST/out.log")"
run_case crcmismatch; { [ $RC = 1 ] && nodelete && nomarker && has 'crc32c mismatch'; } && ok "crcmismatch: refused" || bad "crcmismatch rc=$RC $(tail -2 "$ST/out.log")"

echo "== objects remain after delete: no marker =="
run_case stillthere; { [ $RC = 1 ] && [ -f "$ST/deleted" ] && nomarker && has 'still has'; } && ok "stillthere: no marker written" || bad "stillthere rc=$RC $(tail -2 "$ST/out.log")"

echo "== rerun: folder already gone, 22 valid deliverables -> verify, marker, no delete call =="
run_case rerun; { [ $RC = 0 ] && nocopy && ! grep -q 'storage rm' "$ST/calls.log" && [ -f "$ST/marker" ] && has 'run folder absent'; } && ok "rerun: verified existing, marker written" || bad "rerun rc=$RC $(tail -2 "$ST/out.log")"

echo "== identical earlier copy for chr4 is kept, the other 21 are copied =="
run_case partialdeliv; { [ $RC = 0 ] && [ -f "$ST/marker" ] && [ "$(wc -l < "$ST/copied" | tr -d ' ')" = 42 ] && has 'chr4: deliverable already present and identical'; } && ok "partialdeliv: identical copy kept, 21 copied" || bad "partialdeliv rc=$RC copied=$(wc -l < "$ST/copied") $(tail -2 "$ST/out.log")"

echo "== stale earlier copy for chr4 (crc differs) is recopied =="
run_case staledeliv; { [ $RC = 0 ] && [ -f "$ST/marker" ] && [ "$(wc -l < "$ST/copied" | tr -d ' ')" = 44 ] && has 'chr4: secured'; } && ok "staledeliv: stale copy replaced, 22 secured" || bad "staledeliv rc=$RC copied=$(wc -l < "$ST/copied") $(tail -2 "$ST/out.log")"

echo "== delete_run_folder=false: secure only, no delete, no marker =="
run_case good false; { [ $RC = 0 ] && nodelete && nomarker && [ "$(wc -l < "$ST/copied" | tr -d ' ')" = 44 ]; } && ok "nodelete: secured, folder kept, no marker" || bad "nodelete rc=$RC $(tail -2 "$ST/out.log")"

echo "----"; echo "PASS=$PASS FAIL=$FAIL"
[ "$FAIL" -eq 0 ]
