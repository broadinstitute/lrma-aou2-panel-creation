version 1.0

workflow SelectSampleVcf {

    input {
        File joint_vcf
        File joint_vcf_idx
        String sample_name
        String output_tag
        String gcs_output_dir
        String? region
    }

    call SelectSampleVcf { input:
        vcf = joint_vcf,
        vcf_idx = joint_vcf_idx,
        sample_name = sample_name,
        output_tag = output_tag,
        gcs_output_dir = gcs_output_dir,
        region = region
    }

    output {
        String vcf_path = SelectSampleVcf.vcf_path
        String vcf_idx_path = SelectSampleVcf.vcf_idx_path
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Boolean? use_ssd
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

task SelectSampleVcf {

    input {
        File vcf
        File vcf_idx
        String sample_name
        String output_tag
        String gcs_output_dir
        String? region
        Int view_verbosity = 8
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        bcftools view --no-version "~{vcf}##idx##~{vcf_idx}" \
            ~{"--regions" + region} \
            --regions-overlap 0 \
            -s ~{sample_name} \
            --verbosity ~{view_verbosity} \
            -Ob -o ~{sample_name}.~{output_tag}.bcf
        bcftools index ~{sample_name}.~{output_tag}.bcf

        gcloud storage cp ~{sample_name}.~{output_tag}.bcf* ~{gcs_output_dir}/
    >>>

    output {
        String vcf_path = "~{gcs_output_dir}/~{sample_name}.~{output_tag}.bcf"
        String vcf_idx_path = "~{gcs_output_dir}/~{sample_name}.~{output_tag}.bcf.csi"
    }
    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            25,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsde-methods/slee/lr-bcftools-patched-gcloud/lr-bcftools-patched-gcloud:1.23"      # see https://github.com/broadinstitute/bcftools-patched
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
