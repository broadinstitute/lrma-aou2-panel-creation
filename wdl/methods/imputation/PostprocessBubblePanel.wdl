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

        Int? reduce_length_threshold
        Float? reduce_af_threshold
    }

    call ConvertToVcfGz as ConvertToVcfGzIdSplit { input:
        vcf = panel_id_split_vcf,
        vcf_idx = panel_id_split_vcf_idx,
        output_prefix = output_prefix + ".id.split"
    }

    if (defined(reduce_length_threshold) && defined(reduce_af_threshold)) {
        call ExtractSVIds { input:
            panel_id_split_vcf = panel_id_split_vcf,
            output_prefix = output_prefix + ".id.split.sv"
        }

        scatter (i in range(length(regions))) {
            call ReduceBubblePanel { input:
                panel_bubble_vcf = panel_bubble_vcf,
                panel_id_split_sv_vcf_gz = ExtractSVIds.panel_id_split_sv_vcf_gz,
                length_threshold = select_first([reduce_length_threshold]),
                af_threshold = select_first([reduce_af_threshold]),
                region = regions[i],
                output_prefix = output_prefix + ".reduced.bubble.region-" + i
            }
        }

        call ConcatVcfs.ConcatVcfs as ConcatReducedBubble { input:
            vcfs = ReduceBubblePanel.reduced_vcf,
            vcf_idxs = ReduceBubblePanel.reduced_vcf_idx,
            output_prefix = output_prefix + ".reduced.bubble"
        }
    }

    File panel_bubble_vcf_ = select_first([ConcatReducedBubble.concatenated_vcf, panel_bubble_vcf])
    File panel_bubble_vcf_idx_ = select_first([ConcatReducedBubble.concatenated_vcf_idx, panel_bubble_vcf_idx])
    String output_prefix_ = if defined(ConcatReducedBubble.concatenated_vcf_idx) then output_prefix + ".reduced" else output_prefix

    scatter (i in range(length(regions))) {
        call PopBubblesPanel { input:
            panel_bubble_vcf = panel_bubble_vcf_,
            panel_bubble_vcf_idx = panel_bubble_vcf_idx_,
            panel_id_split_vcf_gz = ConvertToVcfGzIdSplit.vcf_gz,
            panel_id_split_vcf_gz_tbi = ConvertToVcfGzIdSplit.vcf_gz_tbi,
            pop_python_script = pop_python_script,
            region = regions[i],
            output_prefix = output_prefix_ + ".popped.region-" + i
        }

        call SplitBubblesPanel { input:
            panel_bubble_vcf = panel_bubble_vcf_,
            panel_bubble_vcf_idx = panel_bubble_vcf_idx_,
            reference_fasta_fai = reference_fasta_fai,
            region = regions[i],
            output_prefix = output_prefix_ + ".bubble.split.region-" + i,
            leave_out_samples = []
        }
    }

    call ConcatVcfs.ConcatVcfs as ConcatPopped { input:
        vcfs = PopBubblesPanel.popped_vcf,
        vcf_idxs = PopBubblesPanel.popped_vcf_idx,
        output_prefix = output_prefix_ + ".popped"
    }

    call ConcatVcfs.ConcatVcfs as ConcatPoppedSitesOnly { input:
        vcfs = PopBubblesPanel.popped_sites_only_vcf,
        vcf_idxs = PopBubblesPanel.popped_sites_only_vcf_idx,
        output_prefix = output_prefix_ + ".popped.sites"
    }

    call ConcatVcfs.ConcatVcfs as ConcatBubbleSplit { input:
        vcfs = SplitBubblesPanel.split_bubbles_vcf,
        vcf_idxs = SplitBubblesPanel.split_bubbles_vcf_idx,
        output_prefix = output_prefix_ + ".bubble.split"
    }

    call ConcatVcfs.ConcatVcfs as ConcatBubbleSplitSitesOnly { input:
        vcfs = SplitBubblesPanel.split_bubbles_sites_only_vcf,
        vcf_idxs = SplitBubblesPanel.split_bubbles_sites_only_vcf_idx,
        output_prefix = output_prefix_ + ".bubble.split.sites"
    }

    if (defined(leave_out_samples)) {
        scatter (i in range(length(regions))) {
            call SplitBubblesPanel as SplitBubblesPanelLeaveOut { input:
                panel_bubble_vcf = panel_bubble_vcf_,
                panel_bubble_vcf_idx = panel_bubble_vcf_idx_,
                reference_fasta_fai = reference_fasta_fai,
                region = regions[i],
                output_prefix = output_prefix_ + ".bubble.split.leaveout.region-" + i,
                leave_out_samples = leave_out_samples
            }
        }

        call ConcatVcfs.ConcatVcfs as ConcatBubbleSplitLeaveOut { input:
            vcfs = SplitBubblesPanelLeaveOut.split_bubbles_vcf,
            vcf_idxs = SplitBubblesPanelLeaveOut.split_bubbles_vcf_idx,
            output_prefix = output_prefix_ + ".bubble.split.leaveout"
        }
        
        call ConcatVcfs.ConcatVcfs as ConcatBubbleSplitSitesOnlyLeaveOut { input:
            vcfs = SplitBubblesPanelLeaveOut.split_bubbles_sites_only_vcf,
            vcf_idxs = SplitBubblesPanelLeaveOut.split_bubbles_sites_only_vcf_idx,
            output_prefix = output_prefix_ + ".bubble.split.sites.leaveout"
        }
    }

    output {
        File panel_id_split_vcf_gz = ConvertToVcfGzIdSplit.vcf_gz
        File panel_id_split_vcf_gz_tbi = ConvertToVcfGzIdSplit.vcf_gz_tbi

        File? panel_id_split_sv_vcf_gz = ExtractSVIds.panel_id_split_sv_vcf_gz
        File? panel_id_split_sv_vcf_gz_tbi = ExtractSVIds.panel_id_split_sv_vcf_gz_tbi
        File? reduced_panel_bubble_vcf = ConcatReducedBubble.concatenated_vcf
        File? reduced_panel_bubble_vcf_idx = ConcatReducedBubble.concatenated_vcf_idx

        File panel_popped_vcf = ConcatPopped.concatenated_vcf
        File panel_popped_vcf_idx = ConcatPopped.concatenated_vcf_idx
        File panel_popped_sites_only_vcf = ConcatPoppedSitesOnly.concatenated_vcf
        File panel_popped_sites_only_vcf_idx = ConcatPoppedSitesOnly.concatenated_vcf_idx

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

task ConvertToVcfGz {
    input {
        File vcf
        File vcf_idx
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size(vcf, "GB"))

    command <<<
        set -euox pipefail

        bcftools view ~{vcf} -W=tbi -Oz -o ~{output_prefix}.vcf.gz
    >>>

    output {
        File vcf_gz = "~{output_prefix}.vcf.gz"
        File vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
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

task PopBubblesPanel {
    input {
        File panel_bubble_vcf
        File panel_bubble_vcf_idx
        File panel_id_split_vcf_gz
        File panel_id_split_vcf_gz_tbi
        File pop_python_script
        String region
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size([panel_bubble_vcf, panel_id_split_vcf_gz], "GB"))

    command <<<
        set -euox pipefail

        bcftools view -r ~{region} --regions-overlap 0 ~{panel_bubble_vcf} | \
            pypy ~{pop_python_script} ~{panel_id_split_vcf_gz} | \
            bcftools +fill-tags -Ou -- -t AC,AN,AF | \
            bcftools sort -W=csi -Ob -o ~{output_prefix}.bcf
        bcftools view ~{output_prefix}.bcf -G -W=tbi -Ob -o ~{output_prefix}.sites.bcf
    >>>

    output {
        File popped_vcf = "~{output_prefix}.bcf"
        File popped_vcf_idx = "~{output_prefix}.bcf.csi"
        File popped_sites_only_vcf = "~{output_prefix}.sites.bcf"
        File popped_sites_only_vcf_idx = "~{output_prefix}.sites.bcf.csi"
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

        bcftools view -r ~{region} --regions-overlap 0 ~{panel_bubble_vcf} -Ou | \
            bcftools norm -m-any -Ou -N | \
            bcftools +setGT -Ou -- -t . -n 0p | \
            bcftools +fill-AN-AC -Ou | \
            bcftools reheader -f ~{reference_fasta_fai} | \
            ~{if length(leave_out_samples_array) > 0 then "bcftools view -S ^" + leave_out_samples_list + " --force-samples -Ou |" else ""} \
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

task ExtractSVIds {
    input {
        File panel_id_split_vcf
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size(panel_id_split_vcf, "GB"))

    command <<<
        set -euox pipefail

        bcftools view -e 'ID ~ "^chr_*"' ~{panel_id_split_vcf} -W=tbi -Oz -o ~{output_prefix}.vcf.gz
    >>>

    output {
        File panel_id_split_sv_vcf_gz = "~{output_prefix}.vcf.gz"
        File panel_id_split_sv_vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
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

task ReduceBubblePanel {
    input {
        File panel_bubble_vcf
        File panel_id_split_sv_vcf_gz
        Int length_threshold = 20
        Float af_threshold = 0.001
        String region
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size([panel_bubble_vcf, panel_id_split_sv_vcf_gz], "GB"))

    command <<<
        set -euox pipefail

        # Create the inline Python filtering script
        cat << 'EOF' > reduce_panel.py
        import sys
        import gzip

        sv_vcf_path = sys.argv[1]
        length_thresh = int(sys.argv[2])
        af_thresh = float(sys.argv[3])

        # 1. Parse valid SV assigned IDs from the ID VCF
        sv_ids = set()
        opener = gzip.open if sv_vcf_path.endswith('.gz') else open
        with opener(sv_vcf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'): continue
                info_col = line.strip().split('\t')[7]
                for item in info_col.split(';'):
                    if item.startswith('ID='):
                        ids = item[3:].split(',')
                        sv_ids.update(ids)

        # 2. Filter the incoming VCF stream
        for line in sys.stdin:
            if line.startswith('#'):
                sys.stdout.write(line)
                continue
                
            fields = line.strip().split('\t')
            ref = fields[3]
            alt = fields[4]
            info_col = fields[7]
            
            alts = alt.split(',')
            is_biallelic = len(alts) == 1
            keep = False
            
            if is_biallelic:
                # Check Biallelic SNV/Indel Criteria
                is_snp = len(ref) == 1 and len(alt) == 1
                abs_ilen = abs(len(ref) - len(alts[0]))
                if is_snp or abs_ilen <= length_thresh:
                    af = 0.0
                    for item in info_col.split(';'):
                        if item.startswith('AF='):
                            try:
                                af = float(item[3:].split(',')[0])
                            except ValueError:
                                pass
                    if af >= af_thresh:
                        keep = True
                        
            if not keep:
                # Check SV Criteria
                # Retains isolated SVs AND multiallelic bubbles if *any* constituent ID is an SV.
                bubble_ids_str = ""
                for item in info_col.split(';'):
                    if item.startswith('ID='):
                        bubble_ids_str = item[3:]
                        break
                        
                if bubble_ids_str:
                    # Flatten all alleles and constituents into a single list by replacing ',' with ':'
                    all_constituent_ids = bubble_ids_str.replace(',', ':').split(':')
                    if any(cid in sv_ids for cid in all_constituent_ids):
                        keep = True
                        
            if keep:
                sys.stdout.write(line)
        EOF

        # Stream uncompressed VCF (-Ov) to the script, then sort and index as BCF
        bcftools +fill-tags -r ~{region} --regions-overlap 0 -Ov ~{panel_bubble_vcf} -- -t AC,AN,AF | \
            pypy reduce_panel.py ~{panel_id_split_sv_vcf_gz} ~{length_threshold} ~{af_threshold} | \
            bcftools view --threads $(nproc) -W=csi -Ob -o ~{output_prefix}.bcf
    >>>

    output {
        File reduced_vcf = "~{output_prefix}.bcf"
        File reduced_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
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
