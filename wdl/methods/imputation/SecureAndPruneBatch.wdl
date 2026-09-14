version 1.0

# SecureAndPruneBatch
#
# Cloud-side cleanup for one finished aou2_50k imputation batch. The laptop only submits it;
# nothing recurring runs locally. The task VM acts as the workspace service account.
#
# Submit only after confirming in Workbench that every job for the batch is terminal and each of
# the 22 chromosomes has a COMPLETED run. The task adds its own guards on top of that:
#
#   0. A storage listing that fails for any reason other than "matched no objects" aborts the
#      task, so a listing error is never mistaken for an absent folder or file.
#   1. Refuses if anything under <bucket>/glimpse2_phase/batchNNN/ was written in the last
#      quiet_hours hours, which would indicate a run or rerun still in progress.
#   2. Finds exactly one final popped BCF per chromosome, from the ConcatVcfs call directory.
#      Per-region PopAndMarginalizeCollisions shards share the basename and are not matched.
#   3. Downloads each BCF + .csi and verifies it: expected sample count, records > 0, exactly the
#      expected single contig, and one identical sample set across all 22 chromosomes.
#   4. Server-side copies BCF + .csi to <bucket>/deliverables/batch-NNN/ and requires both
#      destination crc32c values to equal the source. Refuses if a different deliverable already
#      exists for that chromosome.
#   5. Deletes the run folder, confirms by listing that nothing remains, and only then writes
#      deliverables/batch-NNN/PRUNED.ok.
#
# Every failure exits non-zero before the delete, or after it but before the marker. Rerunning is
# safe: a chromosome whose run output is already gone is accepted only if its existing
# deliverable passes the same verification. Tests: scripts/tests/test_secure_and_prune_batch.sh

workflow SecureAndPruneBatch {
    input {
        Int batch
        String bucket = "gs://longreadsphase2imputation"
        Int expected_samples = 250
        Int quiet_hours = 6
        Boolean delete_run_folder = true
        String docker = "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.23"
    }

    call SecureAndPrune {
        input:
            batch = batch,
            bucket = bucket,
            expected_samples = expected_samples,
            quiet_hours = quiet_hours,
            delete_run_folder = delete_run_folder,
            docker = docker
    }

    output {
        File report = SecureAndPrune.report
    }
}

task SecureAndPrune {
    input {
        Int batch
        String bucket
        Int expected_samples
        Int quiet_hours
        Boolean delete_run_folder
        String docker
    }

    command <<<
        set -euo pipefail
        B=$(printf '%03d' ~{batch})
        RUN="~{bucket}/glimpse2_phase/batch$B"
        DELIV="~{bucket}/deliverables/batch-$B"
        REPORT="$PWD/report.txt"
        : > "$REPORT"
        log() { echo "[batch-$B] $*" | tee -a "$REPORT"; }
        die() { log "REFUSING: $*"; exit 1; }

        for t in gcloud bcftools md5sum awk date; do
            command -v "$t" >/dev/null || { echo "missing tool in image: $t" >&2; exit 2; }
        done

        # Lists URL into FILE. Returns 0 if objects were listed, 1 if the URL matched no objects.
        # Any other failure aborts the whole task (exit inside a function, not a subshell).
        list_or_none() {  # $1=url  $2=out file  $3=optional -l
            local err
            err=$(mktemp)
            if gcloud storage ls ${3:-} "$1" > "$2" 2> "$err"; then rm -f "$err"; return 0; fi
            if grep -q 'matched no objects' "$err"; then : > "$2"; rm -f "$err"; return 1; fi
            log "ABORT: storage listing failed for $1: $(head -c 300 "$err")"
            exit 3
        }

        # crc32c of one object. Always assign it to a variable: a failed lookup then aborts under
        # set -e instead of comparing two empty strings as equal.
        crc() {
            local h
            h=$(gcloud storage objects describe "$1" --format='value(crc32c_hash)') || { echo "cannot read crc32c for $1" >&2; return 3; }
            [ -n "$h" ] || { echo "empty crc32c for $1" >&2; return 3; }
            printf '%s\n' "$h"
        }

        # Verifies a local BCF; prints its sample-set md5. Returns non-zero on any failure.
        verify_local() {  # $1=local bcf  $2=chrom
            local f="$1" c="$2" n r ctg
            [ -s "$f" ] && [ -s "$f.csi" ] || { echo "$c: missing bcf or .csi" >&2; return 1; }
            n=$(bcftools query -l "$f" | wc -l | tr -d ' ') || return 1
            [ "$n" = "~{expected_samples}" ] || { echo "$c: samples=$n" >&2; return 1; }
            r=$(bcftools index -n "$f") || return 1
            [ "${r:-0}" -gt 0 ] || { echo "$c: records=0 (header-only)" >&2; return 1; }
            ctg=$(bcftools index --stats "$f" | awk '$3>0{print $1}' | sort -u | tr '\n' ' ' | sed 's/ $//') || return 1
            [ "$ctg" = "$c" ] || { echo "$c: contigs='$ctg'" >&2; return 1; }
            bcftools query -l "$f" | sort | md5sum | cut -d' ' -f1
        }

        # --- 1. one snapshot of the run folder, with timestamps
        run_present=0
        list_or_none "$RUN/**" run_objects.txt -l || run_present=$?
        if [ "$run_present" = 0 ]; then
            newest=$(awk '$2 ~ /^[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9:]*Z$/ {print $2}' run_objects.txt | sort | tail -1)
            [ -n "$newest" ] || die "could not read object timestamps from the run folder listing"
            cutoff=$(date -u -d "~{quiet_hours} hours ago" +%Y-%m-%dT%H:%M:%SZ)
            log "run folder: $(awk '$NF ~ /^gs:/' run_objects.txt | wc -l | tr -d ' ') objects, newest write $newest, cutoff $cutoff"
            [[ "$newest" < "$cutoff" ]] || die "run folder written within ~{quiet_hours} h, possible active run"
        else
            log "run folder absent"
        fi

        # --- 2. candidate outputs and existing deliverables
        { awk '{print $NF}' run_objects.txt | grep -E '/call-ConcatVcfs/[^/]*\.glimpse2\.popped\.bcf$' || true; } | sort -u > sources.txt
        list_or_none "$DELIV/*.glimpse2.popped.bcf" deliverables.txt || true
        log "final popped BCFs in run folder: $(wc -l < sources.txt | tr -d ' '), existing deliverables: $(wc -l < deliverables.txt | tr -d ' ')"

        # --- 3/4. verify and secure each chromosome
        : > sample_md5s.txt
        mkdir -p work
        for k in $(seq 1 22); do
            c="chr$k"
            src=$(grep -E "\.batch-$B\.$c\.[^/]*\.glimpse2\.popped\.bcf$" sources.txt || true)
            have=$(grep -E "\.batch-$B\.$c\.[^/]*\.glimpse2\.popped\.bcf$" deliverables.txt || true)
            nsrc=$(printf '%s' "$src" | grep -c . || true)
            nhave=$(printf '%s' "$have" | grep -c . || true)
            [ "$nsrc" -le 1 ] || die "$c: $nsrc distinct final popped BCFs in the run folder"
            [ "$nhave" -le 1 ] || die "$c: $nhave deliverables already present"
            rm -rf work/*

            if [ "$nsrc" = 1 ]; then
                dest="$DELIV/$(basename "$src")"
                [ "$nhave" = 0 ] || [ "$have" = "$dest" ] || die "$c: existing deliverable $(basename "$have") differs from run output $(basename "$src")"
                gcloud storage cp "$src" "$src.csi" work/ --quiet
                md5=$(verify_local "work/$(basename "$src")" "$c") || die "$c: run output failed verification"
                s_bcf=$(crc "$src")
                s_csi=$(crc "$src.csi")
                d_bcf=""
                d_csi=""
                if [ "$nhave" = 1 ]; then
                    d_bcf=$(crc "$dest")
                    if list_or_none "$dest.csi" /dev/null; then d_csi=$(crc "$dest.csi"); fi
                fi
                if [ "$d_bcf" = "$s_bcf" ] && [ "$d_csi" = "$s_csi" ]; then
                    log "$c: deliverable already present and identical"
                else
                    gcloud storage cp "$src" "$src.csi" "$DELIV/" --quiet
                    d_bcf=$(crc "$dest")
                    d_csi=$(crc "$dest.csi")
                    [ "$d_bcf" = "$s_bcf" ] || die "$c: crc32c mismatch after copy (bcf)"
                    [ "$d_csi" = "$s_csi" ] || die "$c: crc32c mismatch after copy (csi)"
                    log "$c: secured $(basename "$src")"
                fi
            elif [ "$nhave" = 1 ]; then
                gcloud storage cp "$have" "$have.csi" work/ --quiet
                md5=$(verify_local "work/$(basename "$have")" "$c") || die "$c: no run output, and the existing deliverable failed verification"
                log "$c: no run output, existing deliverable verified"
            else
                die "$c: no run output and no deliverable"
            fi
            printf '%s\n' "$md5" >> sample_md5s.txt
        done
        rm -rf work

        n=$(wc -l < sample_md5s.txt | tr -d ' ')
        u=$(sort -u sample_md5s.txt | wc -l | tr -d ' ')
        [ "$n" = 22 ] || die "only $n/22 chromosomes verified"
        [ "$u" = 1 ] || die "sample sets differ across chromosomes ($u distinct)"
        log "verified 22/22 @ ~{expected_samples} samples, records>0, single contig, one sample set"

        # --- 5. delete, confirm gone, then marker
        if [ "~{delete_run_folder}" != "true" ]; then
            log "delete_run_folder=false: run folder kept, no PRUNED.ok written"
            log "DONE"
            exit 0
        fi
        if [ "$run_present" = 0 ]; then
            gcloud storage rm -r "$RUN/" --quiet
            log "deleted $RUN/"
        fi
        still=0
        list_or_none "$RUN/**" remaining.txt || still=$?
        [ "$still" = 1 ] || die "run folder still has $(wc -l < remaining.txt | tr -d ' ') objects after delete; no marker written"
        printf 'pruned %s by SecureAndPruneBatch (cloud); 22/22 verified at %s samples\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "~{expected_samples}" \
            | gcloud storage cp - "$DELIV/PRUNED.ok" --quiet
        list_or_none "$DELIV/PRUNED.ok" /dev/null || die "marker not visible after write"
        log "wrote $DELIV/PRUNED.ok"
        log "DONE"
    >>>

    output {
        File report = "report.txt"
    }

    runtime {
        docker: docker
        cpu: 2
        memory: "4 GiB"
        disks: "local-disk 30 HDD"
        preemptible: 3
        maxRetries: 1
    }
}
