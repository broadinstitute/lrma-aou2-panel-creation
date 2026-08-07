version 1.0

# Stage public 1kGP DRAGEN gVCFs into GCS, cut to the regions being validated.
#
# The AWS Open Data copy is one 4.3 GB gVCF per sample. Validating a single chromosome needs
# roughly a fiftieth of that, and the bucket ships a .tbi, so this streams the slice straight
# off the public HTTPS endpoint instead of copying the whole file down. For chr22 across the
# leaveout set that is ~70 MB per sample rather than 4.3 GB -- the difference between ~20 GB
# and ~1.2 TB of egress, localization time and GCS residency.
#
# TransferAWSOpenData in this repo does the whole-file copy and is the fallback if the backend
# turns out to block outbound HTTPS to s3.amazonaws.com but permits the aws CLI. That workflow
# is also the thing to reach for if a future evaluation wants whole genomes.
#
# Requires an htslib built with libcurl. Verified against
#   .../hg38-alt_masked.cnv.graph.hla.methyl_cg.rna-11-r5.0-2/<sample>/<sample>.hard-filtered.gvcf.gz
# which is DRAGEN 4.4.7, GATK-style reference blocks (END/GQ/MIN_DP) and plain FORMAT/PL --
# no LAA/LPL local alleles -- so extract-bubble-PLs reads it through its PL path unchanged.

workflow StageDragen1kGPGvcfs {
    input {
        Array[String] sample_ids
        Array[String] regions              # e.g. ["chr22"]; passed to bcftools -r
        String output_prefix

        # The published analysis directory. Pinned because the bucket also holds per-run
        # directories (<sample>_dragen-germline-v4-4-7-<uuid>) with the same data under a
        # different name.
        String s3_https_prefix = "https://1000genomes-dragen-v4-4-7.s3.amazonaws.com/data/individuals/hg38-alt_masked.cnv.graph.hla.methyl_cg.rna-11-r5.0-2"
    }

    scatter (sample_id in sample_ids) {
        call SliceRemoteGvcf {
            input:
                sample_id = sample_id,
                regions = regions,
                s3_https_prefix = s3_https_prefix,
                output_prefix = output_prefix
        }
    }

    output {
        Array[File] sliced_gvcfs = SliceRemoteGvcf.sliced_gvcf
        Array[File] sliced_gvcf_idxs = SliceRemoteGvcf.sliced_gvcf_idx
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    String? disk_type
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

task SliceRemoteGvcf {
    input {
        String sample_id
        Array[String] regions
        String s3_https_prefix
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    String gvcf_url = s3_https_prefix + "/" + sample_id + "/" + sample_id + ".hard-filtered.gvcf.gz"
    String out_name = output_prefix + "." + sample_id + ".gvcf.gz"

    command <<<
        set -euxo pipefail

        # htslib reads the .tbi from the same prefix on its own, so only the region list is
        # transferred. --no-version keeps the header byte-identical across reruns, which
        # matters because a rerun after preemption must not change the staged file.
        bcftools view \
            --no-version \
            -r ~{sep="," regions} \
            -Oz -o ~{out_name} \
            "~{gvcf_url}"

        bcftools index -t ~{out_name}

        # A silently empty slice means the region list and the gVCF's contig naming disagree
        # (this data is hg38 with chr-prefixed contigs). Fail here rather than let a sample of
        # nothing flow into preprocessing, where it would look like a genuinely uncovered
        # sample.
        n=$(bcftools index -n ~{out_name})
        echo "records in slice: ${n}"
        if [ "${n}" -eq 0 ]; then
            echo "ERROR: no records for ~{sample_id} in ~{sep="," regions}" >&2
            exit 1
        fi
    >>>

    output {
        File sliced_gvcf = "~{out_name}"
        File sliced_gvcf_idx = "~{out_name}.tbi"
    }

    #########################
    # A chr22 slice lands at ~70 MB and the task holds nothing else, so the disk is boot plus
    # slack. Streaming is network-bound, not CPU-bound; extra cores buy nothing here.
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            32,
        boot_disk_gb:       10,
        disk_type:          "HDD",
        preemptible_tries:  3,
        max_retries:        2,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.23"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type])
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
