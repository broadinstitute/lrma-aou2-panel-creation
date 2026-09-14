#!/bin/bash
# verify_pruned_batch.sh N: storage-only check that batch N was cleaned up correctly.
export CLOUDSDK_CORE_DISABLE_PROMPTS=1
B=$(printf '%03d' "$1"); BK=gs://longreadsphase2imputation; D=$BK/deliverables/batch-$B
out=$(gcloud storage ls "$D/" 2>&1); rc=$?
[ $rc = 0 ] || { echo "batch-$B: FAIL listing deliverables: $(echo "$out" | tail -1)"; exit 1; }
nb=$(grep -c 'popped\.bcf$' <<<"$out"); nc=$(grep -c 'popped\.bcf\.csi$' <<<"$out"); nm=$(grep -c 'PRUNED\.ok$' <<<"$out"); no=$(grep -vcE 'popped\.bcf(\.csi)?$|PRUNED\.ok$' <<<"$out")
run=$(gcloud storage ls "$BK/glimpse2_phase/batch$B/**" 2>&1 | tail -1)
case "$run" in *"matched no objects"*) gone=yes;; *) gone="NO: $run";; esac
export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token 2>/dev/null)
spot=""
for c in chr1 chr22; do f=$(grep -E "\.batch-$B\.$c\.[^/]*popped\.bcf$" <<<"$out" | head -1)
  spot="$spot $c:$(bcftools query -l "$f" 2>/dev/null | wc -l | tr -d ' ')s/$(bcftools index -n "$f" 2>/dev/null)r"; done
ok=FAIL; [ "$nb" = 22 ] && [ "$nc" = 22 ] && [ "$nm" = 1 ] && [ "$no" = 0 ] && [ "$gone" = yes ] && [[ "$spot" = *"chr1:250s"* ]] && [[ "$spot" = *"chr22:250s"* ]] && ok=OK
echo "batch-$B: $ok  bcf=$nb csi=$nc marker=$nm other=$no run_folder_gone=$gone spot:$spot"
[ "$ok" = OK ]
