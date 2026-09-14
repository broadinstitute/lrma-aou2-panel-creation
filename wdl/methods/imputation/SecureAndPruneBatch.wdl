version 1.0

# SecureAndPruneBatch
#
# Cloud-side replacement for the laptop cron (scripts/prune_batch.sh). Submit it once for a
# batch whose 22 per-chromosome GLIMPSE2FromPreprocessedPLsJoint runs have all COMPLETED.
# Nothing runs on a local machine; the task VM uses the workspace service account.
#
# For batch NNN it:
#   1. finds exactly one final popped BCF per chromosome under
#      <bucket>/glimpse2_phase/batchNNN/ (the ConcatVcfs call directory; the per-region
#      PopAndMarginalizeCollisions shards write the same basename and are not matched)
#   2. downloads each BCF + .csi and verifies it locally: expected sample count, records > 0,
#      exactly the expected single contig, and one identical sample set across all 22
#   3. server-side copies BCF + .csi to <bucket>/deliverables/batch-NNN/ and requires the
#      destination crc32c to equal the source
#   4. writes deliverables/batch-NNN/PRUNED.ok, then deletes <bucket>/glimpse2_phase/batchNNN/
#      (bubble posteriors, stray GLIMPSE2 checkpoints, logs)
#
# Any failure exits non-zero before step 4, so nothing is deleted unless all 22 chromosomes
# are secured and verified. Idempotent: a chromosome whose run folder is already gone is
# accepted only if its existing deliverable passes the same verification.

workflow SecureAndPruneBatch {
    input {
        Int batch
        String bucket = "gs://longreadsphase2imputation"
        Int expected_samples = 250
        Boolean delete_run_folder = true
        String docker = "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.23"
    }

    call SecureAndPrune {
        input:
            batch = batch,
            bucket = bucket,
            expected_samples = expected_samples,
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
        Boolean delete_run_folder
        String docker
    }

    command <<<
        set -euo pipefail
        B=$(printf '%03d' ~{batch})
        RUN="~{bucket}/glimpse2_phase/batch$B"
        DELIV="~{bucket}/deliverables/batch-$B"
        REPORT="$PWD/report.txt"
        log() { echo "[batch-$B] $*" | tee -a "$REPORT"; }

        for t in gcloud bcftools md5sum; do
            command -v "$t" >/dev/null || { echo "missing tool in image: $t" >&2; exit 2; }
        done

        # --- sources: one final popped BCF per chromosome (sort -u: ** can list an object twice)
        { gcloud storage ls "$RUN/**/call-ConcatVcfs/*.glimpse2.popped.bcf" 2>/dev/null || true; } \
            | sort -u > sources.txt
        log "run folder: $RUN  popped BCFs found: $(wc -l < sources.txt)"

        crc() { gcloud storage objects describe "$1" --format='value(crc32c_hash)'; }

        verify_local() {  # $1=local bcf  $2=chrom  -> prints sample-set md5
            local f="$1" c="$2" n r ctg
            [ -s "$f.csi" ] || { echo "no .csi for $f" >&2; return 1; }
            n=$(bcftools query -l "$f" | wc -l | tr -d ' ')
            [ "$n" = "~{expected_samples}" ] || { echo "$c samples=$n" >&2; return 1; }
            r=$(bcftools index -n "$f")
            [ "${r:-0}" -gt 0 ] || { echo "$c records=0 (header-only)" >&2; return 1; }
            ctg=$(bcftools index --stats "$f" | awk '$3>0{print $1}' | sort -u | tr '\n' ' ' | sed 's/ $//')
            [ "$ctg" = "$c" ] || { echo "$c contigs='$ctg'" >&2; return 1; }
            bcftools query -l "$f" | sort | md5sum | cut -d' ' -f1
        }

        : > sample_md5s.txt
        mkdir -p work
        for k in $(seq 1 22); do
            c="chr$k"
            src=$(grep -E "\.batch-$B\.$c\.[^/]*\.glimpse2\.popped\.bcf$" sources.txt || true)
            nsrc=$(printf '%s\n' "$src" | grep -c . || true)
            existing=$({ gcloud storage ls "$DELIV/*.batch-$B.$c.*.glimpse2.popped.bcf" 2>/dev/null || true; } | head -1)

            if [ "$nsrc" -gt 1 ]; then
                log "$c: $nsrc distinct popped BCFs in run folder, refusing"; printf '%s\n' "$src" >&2; exit 1
            elif [ "$nsrc" -eq 1 ]; then
                gcloud storage cp "$src" "$src.csi" work/ --quiet
                md5=$(verify_local "work/$(basename "$src")" "$c") || { log "$c: source failed verification"; exit 1; }
                dest="$DELIV/$(basename "$src")"
                if [ -n "$existing" ] && [ "$existing" = "$dest" ] \
                   && [ "$(crc "$src")" = "$(crc "$dest")" ] && [ "$(crc "$src.csi")" = "$(crc "$dest.csi" 2>/dev/null || echo none)" ]; then
                    log "$c: deliverable already present and identical"
                else
                    gcloud storage cp "$src" "$src.csi" "$DELIV/" --quiet
                    [ "$(crc "$src")" = "$(crc "$dest")" ] || { log "$c: crc32c mismatch after copy (bcf)"; exit 1; }
                    [ "$(crc "$src.csi")" = "$(crc "$dest.csi")" ] || { log "$c: crc32c mismatch after copy (csi)"; exit 1; }
                    log "$c: secured $(basename "$src")"
                fi
            elif [ -n "$existing" ]; then
                gcloud storage cp "$existing" "$existing.csi" work/ --quiet
                md5=$(verify_local "work/$(basename "$existing")" "$c") || { log "$c: existing deliverable failed verification, no source to recopy"; exit 1; }
                log "$c: no run output, existing deliverable verified"
            else
                log "$c: no run output and no deliverable, refusing"; exit 1
            fi
            echo "$md5" >> sample_md5s.txt
            rm -rf work/*
        done

        n=$(wc -l < sample_md5s.txt | tr -d ' ')
        u=$(sort -u sample_md5s.txt | wc -l | tr -d ' ')
        [ "$n" = 22 ] || { log "only $n/22 chromosomes verified, refusing"; exit 1; }
        [ "$u" = 1 ] || { log "sample sets differ across chromosomes ($u distinct), refusing"; exit 1; }
        log "verified 22/22 @ ~{expected_samples} samples, records>0, single contig, one sample set"

        printf 'pruned %s by SecureAndPruneBatch (cloud)\n' "$(date -u +%FT%TZ)" | gcloud storage cp - "$DELIV/PRUNED.ok" --quiet

        if [ "~{delete_run_folder}" = "true" ]; then
            if gcloud storage ls "$RUN/" >/dev/null 2>&1; then
                gcloud storage rm -r "$RUN/" --quiet
                log "deleted $RUN/"
            else
                log "run folder already gone"
            fi
        else
            log "delete_run_folder=false, run folder kept"
        fi
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
