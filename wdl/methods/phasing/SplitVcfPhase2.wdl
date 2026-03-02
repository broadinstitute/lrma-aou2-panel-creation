version 1.0

workflow SplitCohortVcf {

    input {
        File joint_vcf
        File joint_vcf_tbi
        String locus
        String output_prefix
    }


    call SubsetAndSplitVcf{ input:
        vcf_gz = joint_vcf,
        vcf_gz_tbi = joint_vcf_tbi,
        locus = locus,
        output_prefix = output_prefix
    }
    

    output {
        Array[String] split_vcf_paths = SubsetAndSplitVcf.split_vcf_paths
        Float number_of_calls = SubsetAndSplitVcf.number_of_calls
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


task SubsetAndSplitVcf {

    input {
        File vcf_gz
        File vcf_gz_tbi
        String locus
        String output_prefix
        String gcs_output
        Int view_verbosity = 3
        RuntimeAttr? runtime_attr_override
    }

    meta {
        description: "Subset a VCF file to a given locus and split"
    }

    parameter_meta {
        vcf_gz: {
            description: "VCF file to be subsetted",
            localization_optional: true
        }
        vcf_gz_tbi: {
            description: "Tabix index for the VCF file",
            localization_optional: true
        }
        locus: "Locus to be subsetted"
        runtime_attr_override: "Override default runtime attributes"
    }

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        mkdir output
        bcftools view --no-version ~{vcf_gz} --regions ~{locus} --regions-overlap 0 --verbosity ~{vcf_gz} -Ou | \
            bcftools +split -Ob -o output
        
        cd output
        for bcf in $(find . -name "*.bcf"); do
            bcf_basename=$(basename "$bcf")
            mv "$bcf" "$bcf_basename.~{output_prefix}.bcf"
        done

        bcftools view -H "$bcf_basename.~{output_prefix}.bcf" | wc -l > number.txt

        cd -

        gcloud storage cp output/*.bcf "~{gcs_output}/~{output_prefix}"
        gsutil ls "~{gcs_output}/~{output_prefix}/*bcf" > output_vcf.txt

    >>>

    output {
        Array[String] split_vcf_paths = read_lines("output_vcf.txt")
        Float number_of_calls = read_float("output/number.txt")

    }
    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            50,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.23"
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
