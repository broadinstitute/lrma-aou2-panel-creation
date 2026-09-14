#!/bin/bash
# Gated cloud cleanup of finished batches, submitted from the laptop, executed in the cloud:
#   phase 1: secure-prune-batch-v2 on batch 005 (already cleaned; nothing can be deleted)
#   phase 2: v2 on batch 006, then independent storage verification
#   phase 3: v2 on the other 22 finished batches, each verified
# START_PHASE=N skips phases before N once they are proven. Job files go to glimpse2_prune/.
# Stops at the first failure. Every wb call goes through wbsafe. Logs: ~/imputation/cleanup_logs
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; W=$HERE/wbsafe; V=$HERE/verify_pruned_batch.sh
L=~/imputation/cleanup_logs; BK=gs://longreadsphase2imputation
export JAVA_HOME=/opt/homebrew/opt/openjdk/libexec/openjdk.jdk/Contents/Home PATH=/opt/homebrew/opt/openjdk/bin:$PATH CLOUDSDK_CORE_DISABLE_PROMPTS=1
say(){ echo "[$(date +%H:%M:%S)] $*"; }
stop(){ say "STOP: $*"; exit 1; }
logged_in(){ timeout 120 wb auth status 2>/dev/null | grep -qx 'LOGGED IN'; }

say "phase 0: waiting for Workbench login"
for i in $(seq 1 90); do logged_in && break; sleep 20; done
logged_in || stop "no Workbench login after about 45 minutes"
say "logged in"

submit(){  # $1=batch  $2=display name  $3=output path -> prints runId
  local resp
  resp=$($W workflow job run --workflow=secure-prune-batch-v2 --job-id="$2" \
      --inputs="{\"SecureAndPruneBatch.batch\": $1}" --output-bucket-id=Long_Reads_Phase_2_Imputation \
      --output-path="glimpse2_prune/$3" --format=JSON 2>&1) || { echo "$resp" | tail -3 >&2; return 1; }
  python3 -c 'import json,sys;print(json.loads(sys.stdin.read())["runId"])' <<<"$resp"
}
wait_all(){  # runIds...; writes "runId STATUS" lines to $L/wait_state.txt
  local deadline=$(( $(date +%s) + 7200 )) id st pending
  while :; do
    logged_in || stop "Workbench login lost while waiting"
    pending=0; : > $L/wait_state.txt
    for id in "$@"; do
      st=$($W workflow job describe --job-id="$id" 2>/dev/null | awk -F': *' '/^Status:/{print $2}')
      echo "$id ${st:-UNKNOWN}" >> $L/wait_state.txt
      case "$st" in COMPLETED|FAILED|CANCELLED) ;; *) pending=$((pending+1));; esac
    done
    [ "$pending" = 0 ] && return 0
    [ "$(date +%s)" -lt "$deadline" ] || return 1
    say "  $pending of $# still running"; sleep 60
  done
}
report(){  # $1=runId -> prints the job's report.txt, read from the output path Workbench records
  local p
  # the path itself contains "gs://", so strip the label instead of splitting on colons
  p=$($W workflow job describe --job-id="$1" 2>/dev/null | sed -n 's/^Output bucket path: *//p' | head -1)
  case "$p" in gs://*) ;; *) echo "no valid output path recorded for $1: '$p'" >&2; return 1;; esac
  gcloud storage cat "$p/call-SecureAndPrune/report.txt"
}

START_PHASE="${START_PHASE:-1}"
if [ "$START_PHASE" -le 1 ]; then
say "phase 1: v2 proof on batch 005 (run folder already gone, nothing to delete)"
id=$(submit 5 prune-v2-rerun-b005 batch005-v2rerun) || stop "submit failed for 005"
say "  submitted $id"; wait_all "$id" || stop "005 proof did not finish in 2 h"
st=$(awk '{print $2}' $L/wait_state.txt); [ "$st" = COMPLETED ] || stop "005 proof ended $st"
rep=$(report "$id"); echo "$rep" | grep -vE 'existing deliverable verified'
[[ "$rep" = *DONE* ]] && [ "$(grep -c 'existing deliverable verified' <<<"$rep")" = 22 ] || stop "005 proof report incomplete"
$V 5 || stop "005 verification failed after proof"
fi

if [ "$START_PHASE" -le 2 ]; then
say "phase 2: v2 real cleanup on batch 006"
gcloud storage ls "$BK/glimpse2_phase/batch006/**" >/dev/null 2>&1 || stop "batch006 run folder unexpectedly absent"
id=$(submit 6 prune-v2-b006 batch006) || stop "submit failed for 006"
say "  submitted $id"; wait_all "$id" || stop "006 did not finish in 2 h"
st=$(awk '{print $2}' $L/wait_state.txt); [ "$st" = COMPLETED ] || stop "006 ended $st"
rep=$(report "$id"); echo "$rep" | grep -vE ': secured |already present and identical'
[[ "$rep" = *"verified 22/22"* && "$rep" = *"deleted gs://"* && "$rep" = *DONE* ]] || stop "006 report incomplete"
$V 6 || stop "006 verification failed"
fi

say "phase 3: remaining 22 batches"
ids=(); : > $L/phase3_ids.txt
for n in 7 8 9 10 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 27 29 30; do
  B=$(printf '%03d' $n)
  if $V $n >/dev/null 2>&1; then say "  batch-$B already verified clean, skipping"; continue; fi
  gcloud storage ls "$BK/glimpse2_phase/batch$B/**" >/dev/null 2>&1 || stop "batch-$B has no run folder and is not verified clean; not submitting"
  id=$(submit $n prune-v2-b$B batch$B) || stop "submit failed for batch-$B"
  echo "$n $id" >> $L/phase3_ids.txt; ids+=("$id"); say "  submitted batch-$B $id"
done
[ "${#ids[@]}" -gt 0 ] && { wait_all "${ids[@]}" || stop "phase 3 jobs did not all finish in 2 h"; }
say "  job statuses: $(awk '{print $2}' $L/wait_state.txt | sort | uniq -c | tr '\n' ' ')"
okc=0; bad=""
for n in 7 8 9 10 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 27 29 30; do
  if $V $n; then okc=$((okc+1)); else bad="$bad $n"; fi
done
say "phase 3 result: verified clean $okc of 22${bad:+; NOT clean:$bad}"
[ -z "$bad" ] || exit 1
say "ALL DONE: batches 001-030 cleaned up and verified"
