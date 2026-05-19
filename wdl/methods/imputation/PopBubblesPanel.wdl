version 1.0

import "../ConcatVcfs.wdl" as ConcatVcfs

workflow PopBubblesPanel {
    input {
        Array[File] panel_vcfs               # per-chromosome, multiallelic bubble
        Array[File] panel_vcf_idxs
        Array[File] panel_id_split_vcfs      # per-chromosome
        Array[File] panel_id_split_vcf_idxs
        Array[String]+ chromosomes

        Array[String] output_prefixes

        File pop_python_script
    }

    scatter (i in range(length(panel_vcfs)))
    {
        scatter (j in range(length(chromosomes))) {
            call PopBubblesPanel as ChromosomePopBubblesPanel {
                input:
                    panel_vcf = panel_vcfs[j],
                    panel_vcf_idx = panel_vcf_idxs[j],
                    panel_id_split_vcf = panel_id_split_vcfs[j],
                    panel_id_split_vcf_idx = panel_id_split_vcf_idxs[j],
                    pop_python_script = pop_python_script,
                    chromosome = chromosomes[j],
                    output_prefix = output_prefixes[i] + "." + chromosomes[j]
            }
        }

        # concat across chromosomes
        call ConcatVcfs.ConcatVcfs as ConcatVcfs {
            input:
                vcfs = ChromosomePopBubblesPanel.popped_vcf,
                vcf_idxs = ChromosomePopBubblesPanel.popped_vcf_idx,
                output_prefix = output_prefixes[i]
        }
    }

    output {
        Array[File] popped_vcfs = ConcatVcfs.concatenated_vcf
        Array[File] popped_vcf_idxs = ConcatVcfs.concatenated_vcf_idx
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

task PopBubblesPanel {
    input {
        File panel_vcf
        File panel_vcf_idx
        File panel_id_split_vcf
        File panel_id_split_vcf_idx
        File pop_python_script
        String chromosome
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 3 * ceil(size(panel_vcf, "GB"))

    command <<<
        set -euox pipefail

        bcftools view ~{panel_vcf} | \
            pypy ~{pop_python_script} ~{panel_id_split_vcf} | \
            bcftools view -W=csi -Ob -o ~{output_prefix}.popped.bcf
    >>>

    output {
        File popped_vcf = "~{output_prefix}.popped.bcf"
        File popped_vcf_idx = "~{output_prefix}.popped.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             6,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsde-methods/slee/lrma-aou2-panel-creation-pypy:v1"
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
