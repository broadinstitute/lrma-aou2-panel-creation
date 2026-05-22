version 1.0

workflow HierarchicallyMergeVcfs {
    input {
        Array[File]? vcfs_array
        Array[File]? vcf_idxs_array
        File? vcfs_fofn
        File? vcf_idxs_fofn
        Array[String] regions   # bcftools regions, e.g. ["chr1,chr2,chr3", "chr4,chr5,chr6", ...]
        Array[Int] batch_sizes  # Parameterizable hierarchical levels, e.g., [100, 50]
        String output_prefix
        String extra_merge_args = "--threads $(nproc) --force-single --merge none --info-rules -"       # non-region args; note "--info-rules -" is needed to turn off DP summation, which can lead to MAX_INT overflows and bad VCF behavior
        String extra_concat_args = "--threads $(nproc) --naive"
    }

    Array[File] vcfs_in = if defined(vcfs_array) then select_first([vcfs_array]) else read_lines(select_first([vcfs_fofn]))
    Array[File] vcf_idxs_in = if defined(vcf_idxs_array) then select_first([vcf_idxs_array]) else read_lines(select_first([vcf_idxs_fofn]))

    # Scatter by region FIRST to isolate chunks and reduce combinatorial explosion
    scatter (j in range(length(regions))) {
        String region = regions[j]
        String region_prefix = output_prefix + ".region-" + j

        # ==========================================
        # LEVEL 0
        # ==========================================
        call CreateBatches as L0_Batches {
            input:
                vcfs = vcfs_in,
                vcf_idxs = vcf_idxs_in,
                batch_size = batch_sizes[0]
        }

        scatter (i in range(length(L0_Batches.vcf_batch_fofns))) {
            call MergeVcfs as L0_Merge {
                input:
                    vcfs = read_lines(L0_Batches.vcf_batch_fofns[i]),
                    vcf_idxs = read_lines(L0_Batches.vcf_idx_batch_fofns[i]),
                    output_prefix = region_prefix + ".L0-" + i,
                    extra_args = "--regions-overlap 0 -r " + region + " " + extra_merge_args
            }
        }

        # ==========================================
        # LEVEL 1
        # ==========================================
        if (length(batch_sizes) > 1) {
            call CreateBatches as L1_Batches {
                input:
                    vcfs = L0_Merge.merged_vcf,
                    vcf_idxs = L0_Merge.merged_vcf_idx,
                    batch_size = batch_sizes[1]
            }
            
            scatter (i in range(length(L1_Batches.vcf_batch_fofns))) {
                call MergeVcfs as L1_Merge {
                    input:
                        vcfs = read_lines(L1_Batches.vcf_batch_fofns[i]),
                        vcf_idxs = read_lines(L1_Batches.vcf_idx_batch_fofns[i]),
                        output_prefix = region_prefix + ".L1-" + i,
                        extra_args = extra_merge_args  # Region is already subset in L0
                }
            }
        }

        # ==========================================
        # LEVEL 2
        # ==========================================
        if (length(batch_sizes) > 2) {
            call CreateBatches as L2_Batches {
                input:
                    vcfs = select_first([L1_Merge.merged_vcf]),
                    vcf_idxs = select_first([L1_Merge.merged_vcf_idx]),
                    batch_size = batch_sizes[2]
            }
            
            scatter (i in range(length(L2_Batches.vcf_batch_fofns))) {
                call MergeVcfs as L2_Merge {
                    input:
                        vcfs = read_lines(L2_Batches.vcf_batch_fofns[i]),
                        vcf_idxs = read_lines(L2_Batches.vcf_idx_batch_fofns[i]),
                        output_prefix = region_prefix + ".L2-" + i,
                        extra_args = extra_merge_args
                }
            }
        }

        # ==========================================
        # FINAL REGION COLLAPSE
        # ==========================================
        # Capture the output from the deepest executed level
        Array[File] deepest_vcfs = select_first([L2_Merge.merged_vcf, L1_Merge.merged_vcf, L0_Merge.merged_vcf])
        Array[File] deepest_idxs = select_first([L2_Merge.merged_vcf_idx, L1_Merge.merged_vcf_idx, L0_Merge.merged_vcf_idx])

        # Safety Net: If the provided `batch_sizes` array wasn't steep enough to collapse 
        # the region down to 1 file, dynamically catch whatever is left and merge it.
        if (length(deepest_vcfs) > 1) {
            call MergeVcfs as FinalRegionMerge {
                input:
                    vcfs = deepest_vcfs,
                    vcf_idxs = deepest_idxs,
                    output_prefix = region_prefix + ".final",
                    extra_args = extra_merge_args
            }
        }

        # Select exactly 1 file for this region to pass to the final concat step
        File final_region_vcf = select_first([FinalRegionMerge.merged_vcf, deepest_vcfs[0]])
        File final_region_idx = select_first([FinalRegionMerge.merged_vcf_idx, deepest_idxs[0]])
    }

    # concatenate all regions together
    call ConcatVcfs {
        input:
            vcfs = final_region_vcf,
            vcf_idxs = final_region_idx,
            output_prefix = output_prefix,
            extra_args = extra_concat_args
    }

    output {
        File merged_vcf = ConcatVcfs.concatenated_vcf
        File merged_vcf_idx = ConcatVcfs.concatenated_vcf_idx
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

task CreateBatches {
    input {
        Array[String] vcfs
        Array[String] vcf_idxs
        Int batch_size

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euox pipefail

        cat ~{write_lines(vcfs)} | split -l ~{batch_size} - vcf_batch_
        cat ~{write_lines(vcf_idxs)} | split -l ~{batch_size} - vcf_idx_batch_
    >>>

    output {
        Array[File] vcf_batch_fofns = glob("vcf_batch_*")
        Array[File] vcf_idx_batch_fofns = glob("vcf_idx_batch_*")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            10,
        boot_disk_gb:       10,
        use_ssd:            false,
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

task MergeVcfs {
    input{
        Array[File] vcfs
        Array[File] vcf_idxs
        String output_prefix
        String? extra_args

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size(vcfs, "GiB"))

    command <<<
        set -euox pipefail

        bcftools merge \
            -l ~{write_lines(vcfs)} \
            ~{extra_args} \
            -W=csi -Ob -o ~{output_prefix}.bcf
    >>>

    output {
        File merged_vcf = "~{output_prefix}.bcf"
        File merged_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  3,
        max_retries:        0,
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

task ConcatVcfs {
    input{
        Array[File] vcfs
        Array[File] vcf_idxs
        String output_prefix
        String? extra_args

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size(vcfs, "GiB"))

    command <<<
        set -euox pipefail

        bcftools concat \
            -f ~{write_lines(vcfs)} \
            ~{extra_args} \
            -Ob -o ~{output_prefix}.bcf
        bcftools index ~{output_prefix}.bcf
    >>>

    output {
        File concatenated_vcf = "~{output_prefix}.bcf"
        File concatenated_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  3,
        max_retries:        0,
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
