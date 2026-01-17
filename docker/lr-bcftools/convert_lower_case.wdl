version 1.0

workflow convert_lower_case {
    input {
        String sample_id
        File vcf_gz
        File vcf_gz_tbi
    }
    parameter_meta {
    }
    
    call convert {
        input:
            vcf = vcf_gz,
            samplename = sample_id
            
    }
    
    output {
    	File converted_svs = convert.subset_vcf
    	File converted_svs_tbi = convert.subset_tbi
        }
}


task convert {

    meta {
        description: ""
    }

    parameter_meta {
        vcf: "VCF file to be cleaned"
    }

    input {
        File vcf
        String samplename
    }

    Int disk_size = 2*ceil(size([vcf], "GB")) + 1
    String docker_dir = "/truvari_intrasample"
    String work_dir = "/cromwell_root/truvari_intrasample"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cp ~{docker_dir}/convert_lower_case.py ~{work_dir}/convert_lower_case.py
        cd ~{work_dir}

        python convert_lower_case.py -i ~{vcf} -o ~{samplename}_sv_cleaned.vcf
        bgzip ~{samplename}_sv_cleaned.vcf ~{samplename}_sv_cleaned.vcf.gz
        tabix -p vcf ~{samplename}_sv_cleaned.vcf.gz
    >>>

    output {
        File subset_vcf = "~{work_dir}/~{samplename}_sv_cleaned.vcf.gz"
        File subset_tbi = "~{work_dir}/~{samplename}_sv_cleaned.vcf.gz.tbi"
    }
    ###################
    runtime {
        cpu: 2
        memory:  "4 GiB"
        disks: "local-disk 50 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     3
        max_retries:           2
        docker:"hangsuunc/cleanvcf:v1"
    }

}
