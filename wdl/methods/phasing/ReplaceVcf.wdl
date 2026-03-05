version 1.0


workflow StatisticalPhasing {

    input {
        File original_vcf
        File original_vcf_tbi
        File new_vcf
        File new_vcf_tbi
        Array[String] samplelist
        String region
        String output_prefix

    }
    call ReplaceSampleVcf { 
        input: 
            original_vcf= original_vcf,
            original_vcf_tbi= original_vcf_tbi,
            new_vcf= new_vcf,
            new_vcf_tbi= new_vcf_tbi,
            region = region,
            samplelist= samplelist,
            output_prefix= output_prefix
    }

    output {
        File combined_vcf = ReplaceSampleVcf.combined_vcf
        File combined_tbi = ReplaceSampleVcf.combined_tbi
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


task ReplaceSampleVcf {

    meta {
        description: "Replace a list of sample VCF file to a new vcf"
    }

    input {
        File original_vcf
        File original_vcf_tbi
        File new_vcf
        File new_vcf_tbi
        String region
        Array[String] samplelist
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([original_vcf, new_vcf], "GB")) + 100

    command <<<
        set -euxo pipefail

        bcftools view -s ^~{sep="," samplelist} -r ~{region} -Oz -o original_tmp.vcf.gz ~{original_vcf}
        bcftools index original_tmp.vcf.gz

        bcftools view -s ~{sep="," samplelist} -r ~{region} -Oz -o new_tmp.vcf.gz ~{new_vcf}
        bcftools index new_tmp.vcf.gz

        bcftools merge original_tmp.vcf.gz new_tmp.vcf.gz -Oz -o ~{output_prefix}.combined.vcf.gz
        bcftools index -t ~{output_prefix}.combined.vcf.gz
    >>>

    output {
        File combined_vcf = "~{output_prefix}.combined.vcf.gz"
        File combined_tbi = "~{output_prefix}.combined.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.2"
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
