version 1.0

workflow PanGeniePanelCreation {
    input {
        File phased_vcf
        File phased_vcf_idx
        File reference_fasta
        String region
        File prepare_vcf_script
        File add_ids_script
        File merge_vcfs_script
        String output_prefix

        Float frac_missing = 0.2
    }

    call PanGeniePanelCreation {
        input:
            phased_vcf = phased_vcf,
            phased_vcf_idx = phased_vcf_idx,
            reference_fasta = reference_fasta,
            region = region,
            prepare_vcf_script = prepare_vcf_script,
            add_ids_script = add_ids_script,
            merge_vcfs_script = merge_vcfs_script,
            frac_missing = frac_missing,
            output_prefix = output_prefix
    }

    output {
        File panel_vcf = PanGeniePanelCreation.panel_vcf
        File panel_vcf_idx = PanGeniePanelCreation.panel_vcf_idx
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

task PanGeniePanelCreation {
    input {
        File phased_vcf
        File phased_vcf_idx
        File reference_fasta
        String region
        String output_prefix

        File prepare_vcf_script
        File add_ids_script
        File merge_vcfs_script
        Float frac_missing

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 3 * ceil(size(phased_vcf, "GiB")) + ceil(size(reference_fasta, "GiB")) 

    command <<<
        set -euxo pipefail

        bcftools stats -r ~{region} --regions-overlap 0 ~{phased_vcf} > ~{output_prefix}.stats.txt

        # validate variants against reference, run PanGenie prepare-vcf and add-ids scripts, and split to biallelic
        bcftools norm --no-version -r ~{region} --regions-overlap 0 --check-ref e --fasta-ref ~{reference_fasta} ~{phased_vcf} | \
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
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             6,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsde-methods/slee/pangenie-panel-creation:v1"
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

    output {
        File input_stats = "~{output_prefix}.stats.txt"
        File panel_stats = "~{output_prefix}.prepare.id.split.mergehap.stats.txt"
        File panel_vcf = "~{output_prefix}.prepare.id.split.mergehap.bcf"
        File panel_vcf_idx = "~{output_prefix}.prepare.id.split.mergehap.bcf.csi"
    }
}
