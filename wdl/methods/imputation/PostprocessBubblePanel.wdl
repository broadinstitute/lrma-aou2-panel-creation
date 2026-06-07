version 1.0

import "../ConcatVcfs.wdl" as ConcatVcfs

workflow PostprocessBubblePanel {
    input {
        File panel_bubble_vcf
        File panel_bubble_vcf_idx
        File panel_id_split_vcf
        File panel_id_split_vcf_idx
        File reference_fasta_fai
        Array[String] regions
        Array[String]? leave_out_samples

        String output_prefix

        File pop_python_script
    }

    scatter (i in range(length(regions))) {
        call PopBubblesPanel { input:
            panel_bubble_vcf = panel_bubble_vcf,
            panel_bubble_vcf_idx = panel_bubble_vcf_idx,
            panel_id_split_vcf = panel_id_split_vcf,
            panel_id_split_vcf_idx = panel_id_split_vcf_idx,
            pop_python_script = pop_python_script,
            region = regions[i],
            output_prefix = output_prefix + ".region-" + i
        }

        call SplitBubblesPanel { input:
            panel_bubble_vcf = panel_bubble_vcf,
            panel_bubble_vcf_idx = panel_bubble_vcf_idx,
            reference_fasta_fai = reference_fasta_fai,
            region = regions[i],
            output_prefix = output_prefix + ".bubble.split.region-" + i,
            leave_out_samples = []
        }
    }

    call ConcatVcfs.ConcatVcfs as ConcatIdSplit { input:
        vcfs = PopBubblesPanel.id_split_vcf_gz,
        vcf_idxs = PopBubblesPanel.id_split_vcf_gz_tbi,
        output_prefix = output_prefix + ".id.split",
        do_bcf = false
    }

    call ConcatVcfs.ConcatVcfs as ConcatPopped { input:
        vcfs = PopBubblesPanel.popped_vcf,
        vcf_idxs = PopBubblesPanel.popped_vcf_idx,
        output_prefix = output_prefix + ".popped"
    }

    call ConcatVcfs.ConcatVcfs as ConcatBubbleSplit { input:
        vcfs = SplitBubblesPanel.split_bubbles_vcf,
        vcf_idxs = SplitBubblesPanel.split_bubbles_vcf_idx,
        output_prefix = output_prefix + ".bubble.split"
    }

    call ConcatVcfs.ConcatVcfs as ConcatBubbleSplitSitesOnly { input:
        vcfs = SplitBubblesPanel.split_bubbles_sites_only_vcf,
        vcf_idxs = SplitBubblesPanel.split_bubbles_sites_only_vcf_idx,
        output_prefix = output_prefix + ".bubble.split.sites"
    }

    if (defined(leave_out_samples)) {
        scatter (i in range(length(regions))) {
            call SplitBubblesPanel as SplitBubblesPanelLeaveOut { input:
                panel_bubble_vcf = panel_bubble_vcf,
                panel_bubble_vcf_idx = panel_bubble_vcf_idx,
                reference_fasta_fai = reference_fasta_fai,
                region = regions[i],
                output_prefix = output_prefix + ".bubble.split.leaveout.region-" + i,
                leave_out_samples = leave_out_samples
            }
        }

        call ConcatVcfs.ConcatVcfs as ConcatBubbleSplitLeaveOut { input:
            vcfs = SplitBubblesPanelLeaveOut.split_bubbles_vcf,
            vcf_idxs = SplitBubblesPanelLeaveOut.split_bubbles_vcf_idx,
            output_prefix = output_prefix + ".bubble.split.leaveout"
        }
        
        call ConcatVcfs.ConcatVcfs as ConcatBubbleSplitSitesOnlyLeaveOut { input:
            vcfs = SplitBubblesPanelLeaveOut.split_bubbles_sites_only_vcf,
            vcf_idxs = SplitBubblesPanelLeaveOut.split_bubbles_sites_only_vcf_idx,
            output_prefix = output_prefix + ".bubble.split.sites.leaveout"
        }
    }

    output {
        File panel_id_split_vcf_gz = ConcatIdSplit.concatenated_vcf
        File panel_id_split_vcf_gz_tbi = ConcatIdSplit.concatenated_vcf_idx
        File panel_popped_vcf = ConcatPopped.concatenated_vcf
        File panel_popped_vcf_idx = ConcatPopped.concatenated_vcf_idx
        File panel_bubble_split_vcf = ConcatBubbleSplit.concatenated_vcf
        File panel_bubble_split_vcf_idx = ConcatBubbleSplit.concatenated_vcf_idx
        File panel_bubble_split_sites_only_vcf = ConcatBubbleSplitSitesOnly.concatenated_vcf
        File panel_bubble_split_sites_only_vcf_idx = ConcatBubbleSplitSitesOnly.concatenated_vcf_idx
        File? panel_bubble_split_leaveout_vcf = ConcatBubbleSplitLeaveOut.concatenated_vcf
        File? panel_bubble_split_leaveout_vcf_idx = ConcatBubbleSplitLeaveOut.concatenated_vcf_idx
        File? panel_bubble_split_sites_only_leaveout_vcf = ConcatBubbleSplitSitesOnlyLeaveOut.concatenated_vcf
        File? panel_bubble_split_sites_only_leaveout_vcf_idx = ConcatBubbleSplitSitesOnlyLeaveOut.concatenated_vcf_idx
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
        File panel_bubble_vcf
        File panel_bubble_vcf_idx
        File panel_id_split_vcf
        File panel_id_split_vcf_idx
        File pop_python_script
        String region
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size([panel_bubble_vcf, panel_id_split_vcf], "GB"))

    command <<<
        set -euox pipefail

        bcftools view -r ~{region} --regions-overlap 0 ~{panel_id_split_vcf} -G -W=tbi -Oz -o ~{output_prefix}.id.split.vcf.gz
        bcftools view -r ~{region} --regions-overlap 0 ~{panel_bubble_vcf} | \
            pypy ~{pop_python_script} ~{output_prefix}.id.split.vcf.gz | \
            bcftools sort -W=csi -Ob -o ~{output_prefix}.popped.bcf
    >>>

    output {
        File id_split_vcf_gz = "~{output_prefix}.id.split.vcf.gz"
        File id_split_vcf_gz_tbi = "~{output_prefix}.id.split.vcf.gz.tbi"
        File popped_vcf = "~{output_prefix}.popped.bcf"
        File popped_vcf_idx = "~{output_prefix}.popped.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
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

task SplitBubblesPanel {
    input {
        File panel_bubble_vcf
        File panel_bubble_vcf_idx
        File reference_fasta_fai
        String region
        String output_prefix

        Array[String]? leave_out_samples

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size(panel_bubble_vcf, "GB"))

    Array[String] leave_out_samples_array = select_first([leave_out_samples, []])
    File leave_out_samples_list = write_lines(leave_out_samples_array)

    command <<<
        set -euox pipefail

        bcftools view -r ~{region} --regions-overlap 0 ~{panel_bubble_vcf} -Ou |
            bcftools norm -m-any -Ou -N |
            bcftools +setGT -Ou -- -t . -n 0p |
            bcftools +fill-AN-AC -Ou |
            bcftools reheader -f ~{reference_fasta_fai} |
            ~{if length(leave_out_samples_array) > 0 then "bcftools view -S ^" + leave_out_samples_list + " --force-samples -Ou |" else ""}
            bcftools sort -W=csi -Ob -o ~{output_prefix}.bcf

        bcftools view -G ~{output_prefix}.bcf -W=csi -Ob -o ~{output_prefix}.sites.bcf
    >>>

    output {
        File split_bubbles_vcf = "~{output_prefix}.bcf"
        File split_bubbles_vcf_idx = "~{output_prefix}.bcf.csi"
        File split_bubbles_sites_only_vcf = "~{output_prefix}.sites.bcf"
        File split_bubbles_sites_only_vcf_idx = "~{output_prefix}.sites.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             6,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
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
