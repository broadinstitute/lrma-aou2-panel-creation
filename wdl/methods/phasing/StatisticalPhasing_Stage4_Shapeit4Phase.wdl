version 1.0

workflow Shapeit4Phase {

    input {
        File vcf
        File vcf_idx
        File genetic_maps_tsv
        String chromosome
        String region
        String output_prefix

        String extra_args = "--thread $(nproc) --use-PS 0.0001 --pbwt-depth 4 --mcmc-iterations 4b,1p,1b,1p,4m"
    }

    Map[String, String] genetic_maps_dict = read_map(genetic_maps_tsv)

    call Shapeit4 { input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        genetic_map = genetic_maps_dict[chromosome],
        region = region,
        output_prefix = output_prefix + ".shapeit4",
        extra_args = extra_args
    }

    output {
        File phased_vcf = Shapeit4.phased_vcf
        File phased_vcf_idx = Shapeit4.phased_vcf_idx
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

task Shapeit4 {
    input {
        File vcf
        File vcf_idx
        File genetic_map
        String region
        String output_prefix
        String extra_args = "--thread $(nproc) --use-PS 0.0001 --pbwt-depth 4 --mcmc-iterations 4b,1p,1b,1p,4m"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 4 * ceil(size(vcf, "GiB"))
    
    command <<<
        set -euxo pipefail

        shapeit4.2 --input ~{vcf} \
                --map ~{genetic_map} \
                --region ~{region} \
                --sequencing \
                --output ~{output_prefix}.bcf \
                ~{extra_args}
        bcftools index ~{output_prefix}.bcf
    >>>

    output{
        File phased_vcf = "~{output_prefix}.bcf"
        File phased_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             16,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/shapeit4:v1"
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
