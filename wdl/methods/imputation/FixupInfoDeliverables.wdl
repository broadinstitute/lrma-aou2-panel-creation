version 1.0

# FixupInfoDeliverables
#
# Backfill for batches imputed before GLIMPSE2FromPreprocessedPLsJoint recomputed INFO/AF after
# popping. pop-glimpse2 rebuilds each popped allele's GP but copies AF and INFO unchanged from the
# split record GLIMPSE2 scored, so those fields do not describe the genotypes in the deliverable.
# This runs Samuel Lee's fixup-info (resources/fixup-info) on each chromosome's popped BCF in
# deliverables/batch-NNN/ with the file's own sample count as a single minibatch.
#
# Per chromosome it refuses to publish unless:
#   - the sample count equals expected_samples and the record count is unchanged,
#   - no record lost its panel ID (fixup-info writes "." for a site absent from the panel),
#   - an md5 over CHROM, POS, ID, REF, ALT and every sample's GT:DS:GP is identical before and
#     after, i.e. only INFO changed.
# Publish copies the fixed BCF and index over the deliverable names server-side, requires matching
# crc32c, and writes deliverables/batch-NNN/INFO_FIXED.ok only after all 22 are replaced.
# Rerunning over a fixed batch reproduces identical INFO and is harmless.

workflow FixupInfoDeliverables {
    input {
        Int batch
        String bucket = "gs://longreadsphase2imputation"
        File fixup_info_binary
        File pop_glimpse2_panel_resources_json
        Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11",
                                     "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]
        Int expected_samples = 250
        Boolean publish = true
        # The popping image (Ubuntu 24.04, glibc 2.39, bcftools 1.23.1). fixup-info needs glibc >= 2.29.
        String fixup_docker = "us.gcr.io/broad-dsde-methods/slee/lrma-aou2-panel-creation-rust@sha256:0f25c4091c49d8eb0c3d8bcdb45e7680a093e0a691e74ccec2d3743758a7d22c"
        String gcloud_docker = "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.23"
    }

    call ListDeliverables { input: batch = batch, bucket = bucket, docker = gcloud_docker }
    Map[String, String] bcf_by_chrom = read_map(ListDeliverables.bcf_by_chrom)
    Map[String, Map[String, String]] panel = read_json(pop_glimpse2_panel_resources_json)

    scatter (c in chromosomes) {
        call FixupChromosome {
            input:
                chrom = c,
                popped_bcf = bcf_by_chrom[c],
                popped_bcf_csi = bcf_by_chrom[c] + ".csi",
                panel_id_split_vcf_gz = panel[c]["panel_id_split_vcf_gz"],
                fixup_info_binary = fixup_info_binary,
                expected_samples = expected_samples,
                docker = fixup_docker
        }
        String fixed_bcf_path = FixupChromosome.fixed_bcf
        String dest_path = bcf_by_chrom[c]
    }

    if (publish) {
        call Publish {
            input:
                batch = batch,
                bucket = bucket,
                fixed_bcfs = fixed_bcf_path,
                dest_bcfs = dest_path,
                n_expected = length(chromosomes),
                docker = gcloud_docker
        }
    }

    output {
        Array[File] reports = FixupChromosome.report
        File? publish_report = Publish.report
    }
}

task ListDeliverables {
    input {
        Int batch
        String bucket
        String docker
    }
    command <<<
        set -euo pipefail
        B=$(printf '%03d' ~{batch})
        gcloud storage ls "~{bucket}/deliverables/batch-$B/*.glimpse2.popped.bcf" > bcfs.txt
        # chromosome<TAB>path, exactly one per chromosome
        awk '{ if (match($0, /\.chr[0-9XY]+\./)) print substr($0, RSTART+1, RLENGTH-2) "\t" $0 }' bcfs.txt | sort -k1,1 > bcf_by_chrom.tsv
        dups=$(cut -f1 bcf_by_chrom.tsv | uniq -d | awk 'END{print NR}')
        [ "$dups" = 0 ] || { echo "ERROR: more than one popped BCF for a chromosome" >&2; cat bcf_by_chrom.tsv >&2; exit 1; }
        cat bcf_by_chrom.tsv
    >>>
    output { File bcf_by_chrom = "bcf_by_chrom.tsv" }
    runtime {
        docker: docker
        cpu: 1
        memory: "2 GiB"
        disks: "local-disk 10 HDD"
        preemptible: 3
        maxRetries: 1
    }
}

task FixupChromosome {
    input {
        String chrom
        File popped_bcf
        File popped_bcf_csi
        File panel_id_split_vcf_gz
        File fixup_info_binary
        Int expected_samples
        String docker
    }
    Int disk_gb = 20 + ceil(3 * size(popped_bcf, "GB") + size(panel_id_split_vcf_gz, "GB"))
    String out = basename(popped_bcf)

    command <<<
        set -euo pipefail
        cp "~{fixup_info_binary}" fixup-info && chmod +x fixup-info
        IN="~{popped_bcf}"
        n=$(bcftools query -l "$IN" | awk 'END{print NR}')
        [ "$n" = "~{expected_samples}" ] || { echo "ERROR: ~{chrom} has $n samples" >&2; exit 1; }

        bcftools view "$IN" | ./fixup-info ~{panel_id_split_vcf_gz} "$n" batched_info.tsv.gz | bcftools view -Ob -o "~{out}"
        bcftools index -f "~{out}"

        rin=$(bcftools index -n "$IN"); rout=$(bcftools index -n "~{out}")
        [ "$rin" = "$rout" ] || { echo "ERROR: ~{chrom} records $rin -> $rout" >&2; exit 1; }
        miss=$(bcftools query -i 'INFO/ID="."' -f '.\n' "~{out}" | awk 'END{print NR}')
        [ "$miss" = 0 ] || { echo "ERROR: ~{chrom} $miss records without a panel ID" >&2; exit 1; }
        m_in=$(bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT:%DS:%GP]\n' "$IN" | md5sum | cut -d' ' -f1)
        m_out=$(bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT:%DS:%GP]\n' "~{out}" | md5sum | cut -d' ' -f1)
        [ "$m_in" = "$m_out" ] || { echo "ERROR: ~{chrom} sites or genotypes changed" >&2; exit 1; }

        changed=$(paste <(bcftools query -f '%INFO/INFO\n' "$IN") <(bcftools query -f '%INFO/INFO\n' "~{out}") | awk '{d=$1-$2; if (d<0) d=-d; if (d>0.1) c++} END{print c+0}')
        printf '%s samples=%s records=%s ids_missing=0 genotypes_md5_equal=yes info_changed_gt_0.1=%s\n' \
            "~{chrom}" "$n" "$rout" "$changed" | tee report.txt
    >>>
    output {
        File fixed_bcf = "~{out}"
        File fixed_bcf_csi = "~{out}.csi"
        File report = "report.txt"
    }
    runtime {
        docker: docker
        cpu: 2
        memory: "16 GiB"
        disks: "local-disk " + disk_gb + " SSD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
    }
}

task Publish {
    input {
        Int batch
        String bucket
        Array[String] fixed_bcfs
        Array[String] dest_bcfs
        Int n_expected
        String docker
    }
    command <<<
        set -euo pipefail
        B=$(printf '%03d' ~{batch})
        crc() { local h; h=$(gcloud storage objects describe "$1" --format='value(crc32c_hash)') || return 3; [ -n "$h" ] || return 3; printf '%s\n' "$h"; }
        paste ~{write_lines(fixed_bcfs)} ~{write_lines(dest_bcfs)} > pairs.tsv
        [ "$(awk 'END{print NR}' pairs.tsv)" = "~{n_expected}" ] || { echo "ERROR: pair count" >&2; exit 1; }
        : > report.txt
        while IFS=$'\t' read -r SRC DST; do
            for ext in "" ".csi"; do
                gcloud storage cp "$SRC$ext" "$DST$ext" --quiet
                s=$(crc "$SRC$ext"); d=$(crc "$DST$ext")
                [ "$s" = "$d" ] || { echo "ERROR: crc mismatch $DST$ext" >&2; exit 1; }
            done
            echo "replaced $DST" | tee -a report.txt
        done < pairs.tsv
        if [ "~{n_expected}" = 22 ]; then
            printf 'INFO/AF recomputed from popped GP by fixup-info %s; 22/22 verified\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
                | gcloud storage cp - "~{bucket}/deliverables/batch-$B/INFO_FIXED.ok" --quiet
            echo "wrote INFO_FIXED.ok" | tee -a report.txt
        fi
    >>>
    output { File report = "report.txt" }
    runtime {
        docker: docker
        cpu: 1
        memory: "2 GiB"
        disks: "local-disk 10 HDD"
        preemptible: 3
        maxRetries: 1
    }
}
