version 1.0

workflow ShapeitLigate {

    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String output_prefix
    }


    call LigateVcfs { input:
        vcfs = vcfs,
        vcf_idxs = vcf_idxs,
        output_prefix = output_prefix + ".shapeit4.ligated"
    }

    output {
        File ligated_vcf = LigateVcfs.ligated_vcf
        File ligated_vcf_idx = LigateVcfs.ligated_vcf_idx
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

task LigateVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 4 * ceil(size(vcfs, "GiB"))

    command <<<
        set -euxo pipefail

        ligate_static --input ~{write_lines(vcfs)} --threads $(nproc) --output ~{output_prefix}.ligate.bcf
        bcftools +fill-tags --no-version --threads $(nproc) ~{output_prefix}.ligate.bcf \
            -Ob -o ~{output_prefix}.bcf -- -t AF,AC,AN
        bcftools index ~{output_prefix}.bcf
    >>>

    output {
        File ligated_vcf = "~{output_prefix}.bcf"
        File ligated_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "SSD",
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsde-methods/slee/shapeit5:v1"
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
