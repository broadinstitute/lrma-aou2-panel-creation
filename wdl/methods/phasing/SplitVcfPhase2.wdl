version 1.0

workflow SplitCohortVcf {

    input {
        File joint_vcf
        File joint_vcf_tbi
        String locus 
        Int memory
    }


    call SubsetandSplitVcf{ input:
        vcf_gz = joint_vcf,
        vcf_gz_tbi = joint_vcf_tbi,
        locus = locus,
        memory = memory
    }
    

    output {
        Array[File] splitted_vcf = SubsetandSplitVcf.splitted_vcf
        Array[File] splitted_vcf_tbi = SubsetandSplitVcf.splitted_vcf_tbi
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


task SubsetandSplitVcf {

    input {
        File vcf_gz
        File vcf_gz_tbi
        String locus
        Int memory
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

        bcftools view --no-version ~{vcf_gz} --regions ~{locus} | \
             bcftools +split -Oz -o output 
        
        cd output
        for vcf in $(find . -name "*.vcf.gz"); do
            vcf_basename=$(basename "$vcf")
            bcftools view "$vcf" -Oz -o "$vcf_basename.~{locus}.vcf.gz"
            bcftools index -t "$vcf_basename.~{locus}.vcf.gz"
            rm "$vcf" 
        done
        cd -

    >>>

    output {
        Array[File] splitted_vcf = glob("output/*.vcf.gz")
        Array[File] splitted_vcf_tbi = glob("output/*.vcf.gz.tbi")

    }
    ###################
    runtime {
        cpu: 1
        memory: memory + " GiB"
        disks: "local-disk 500 LOCAL"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           1
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }
}