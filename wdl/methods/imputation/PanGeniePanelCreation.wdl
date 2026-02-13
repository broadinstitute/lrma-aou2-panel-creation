version 1.0

struct RuntimeAttributes {
    Int? cpu
    Int? command_mem_gb
    Int? additional_mem_gb
    Int? disk_size_gb
    Int? boot_disk_size_gb
    Boolean? use_ssd
    Int? preemptible
    Int? max_retries
}

workflow PanGeniePanelCreation {
    input {
        File phased_vcf
        File reference_fasta
        File prepare_vcf_script
        File add_ids_script
        File merge_vcfs_script
        Float frac_missing = 0.2
        String output_prefix

        String docker
    }

    call PanGeniePanelCreation {
        input:
            phased_vcf = phased_vcf,
            reference_fasta = reference_fasta,
            prepare_vcf_script = prepare_vcf_script,
            add_ids_script = add_ids_script,
            merge_vcfs_script = merge_vcfs_script,
            frac_missing = frac_missing,
            output_prefix = output_prefix,
            docker = docker
    }

    output {
        File panel_vcf = PanGeniePanelCreation.panel_vcf
        File panel_vcf_idx = PanGeniePanelCreation.panel_vcf_idx
    }
}

task PanGeniePanelCreation {
    input {
        File phased_vcf
        File reference_fasta
        String output_prefix

        File prepare_vcf_script
        File add_ids_script
        File merge_vcfs_script
        Float frac_missing

        String docker

        RuntimeAttributes runtime_attributes = {}
    }

    command {
        set -euxo pipefail

        bcftools stats ~{phased_vcf} > ~{output_prefix}.stats.txt 

        # validate variants against reference, run PanGenie prepare-vcf and add-ids scripts, and split to biallelic
        bcftools norm --no-version --check-ref e --fasta-ref ~{reference_fasta} ~{phased_vcf} | \
            pypy ~{prepare_vcf_script} --missing ~{frac_missing} | \
            pypy ~{add_ids_script} | \
            bcftools norm --no-version -m-any -Ov -o prepare.id.split.vcf

        # run PanGenie merge script
        pypy -m pip install pyfaidx
        bcftools view --no-version -h prepare.id.split.vcf > header.txt
        pypy ~{merge_vcfs_script} merge \
            -vcf prepare.id.split.vcf \
            -header header.txt \
            -r ~{reference_fasta} \
            -ploidy 2 | \
            bcftools view --no-version -Ob -o ~{output_prefix}.prepare.id.split.mergehap.bcf
        bcftools index ~{output_prefix}.prepare.id.split.mergehap.bcf

        bcftools stats ~{output_prefix}.prepare.id.split.mergehap.bcf > ~{output_prefix}.prepare.id.split.mergehap.stats.txt
    }

    runtime {
        docker: docker
        cpu: select_first([runtime_attributes.cpu, 1])
        memory: select_first([runtime_attributes.command_mem_gb, 6]) + select_first([runtime_attributes.additional_mem_gb, 1]) + " GB"
        disks: "local-disk " + select_first([runtime_attributes.disk_size_gb, 500]) + if select_first([runtime_attributes.use_ssd, false]) then " SSD" else " HDD"
        bootDiskSizeGb: select_first([runtime_attributes.boot_disk_size_gb, 15])
        preemptible: select_first([runtime_attributes.preemptible, 2])
        maxRetries: select_first([runtime_attributes.max_retries, 1])
    }

    output {
        File input_stats = "~{output_prefix}.stats.txt"
        File panel_stats = "~{output_prefix}.prepare.id.split.mergehap.stats.txt"
        File panel_vcf = "~{output_prefix}.prepare.id.split.mergehap.bcf"
        File panel_vcf_idx = "~{output_prefix}.prepare.id.split.mergehap.bcf.csi"
    }
}
