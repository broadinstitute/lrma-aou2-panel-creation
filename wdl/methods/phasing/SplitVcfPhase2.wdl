version 1.0

workflow SplitCohortVcf {

    input {
        File joint_vcf
        File joint_vcf_tbi
        String locus
        String output_prefix
    }


    call SubsetandSplitVcf{ input:
        vcf_gz = joint_vcf,
        vcf_gz_tbi = joint_vcf_tbi,
        locus = locus,
        output_prefix = output_prefix
    }
    

    output {
        Array[File] splitted_vcf = SubsetandSplitVcf.splitted_vcf
        Array[File] splitted_vcf_tbi = SubsetandSplitVcf.splitted_vcf_tbi
        Float number_of_call = SubsetandSplitVcf.number_of_call
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


task SubsetandSplitVcf {

    input {
        File vcf_gz
        File vcf_gz_tbi
        String locus
        String output_prefix
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
        bcftools view --no-version ~{vcf_gz} --regions ~{locus} -Oz -o ~{output_prefix}.vcf.gz
        bcftools index -t ~{output_prefix}.vcf.gz

        bcftools +split -Oz -o output ~{output_prefix}.vcf.gz
        
        cd output
        for vcf in $(find . -name "*.vcf.gz"); do
            vcf_basename=$(basename "$vcf")
            bcftools view "$vcf" -Oz -o "$vcf_basename.~{output_prefix}.~{locus}.vcf.gz"
            bcftools index -t "$vcf_basename.~{output_prefix}.~{locus}.vcf.gz"
            rm "$vcf" 
        done


        bcftools view -H "$vcf_basename.~{output_prefix}.~{locus}.vcf.gz" | wc -l > number.txt

        cd -

    >>>

    output {
        Array[File] splitted_vcf = glob("output/*.vcf.gz")
        Array[File] splitted_vcf_tbi = glob("output/*.vcf.gz.tbi")
        Float number_of_call = read_float("coverage.txt")

    }
    ###################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            10,
        boot_disk_gb:       10,
        use_ssd:            false,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.2"
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