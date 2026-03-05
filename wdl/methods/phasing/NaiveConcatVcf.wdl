version 1.0

workflow NaiveConcatVcf {

    input {
        Array[Array[File]] region_by_sample_vcf
        Array[String] sample_list
        Int memory
        String outputprefix
    }

    Array[Array[File]] sample_by_region_vcf_gzs = transpose(region_by_sample_vcf)

    scatter (i in range(length(sample_by_region_vcf_gzs))) {
        call BcftoolsConcatNaive {
            input:
                vcfs = sample_by_region_vcf_gzs[i],
                output_prefix = sample_list[i]
        }
    }

    output {
        Array[File] sample_concat_vcf = BcftoolsConcatNaive.concatenated_vcf
        Array[File] sample_concat_vcf_tbi = BcftoolsConcatNaive.concatenated_vcf_tbi
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

struct DataTypeParameters {
    Int num_shards
    String map_preset
}



task BcftoolsConcatNaive {

    input {
        Array[File] vcfs
        Array[File]? vcf_tbis
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        # Index all input VCF files if no precomputed tbis provided.
        if ! ~{defined(vcf_tbis)}; then
            for vcf in ~{sep=" " vcfs}; do
                bcftools index "$vcf"
            done
        fi

        bcftools concat \
            ~{sep=" " vcfs} \
            --naive-force \
            --no-version | \
            bcftools sort -Oz -o ~{output_prefix}.vcf.gz

        bcftools index -t ~{output_prefix}.vcf.gz

    >>>

    output {
        File concatenated_vcf = "~{output_prefix}.vcf.gz"
        File concatenated_vcf_tbi = "~{output_prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:            "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
