version 1.0

workflow HierarchicallyMergeVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] regions   # bcftools regions, e.g. ["chr1,chr2,chr3", "chr4,chr5,chr6", ...]
        Int batch_size
        String output_prefix
        String extra_merge_args = "--threads $(nproc) --force-single --merge none"       # non-region args
        String extra_concat_args = "--threads $(nproc) --naive"
        Boolean use_ivcfmerge   # requires sample_names
        Array[String]? sample_names
    }

    call CreateBatches {
        input:
            vcfs = vcfs,
            vcf_idxs = vcf_idxs,
            batch_size = batch_size
    }

    if (use_ivcfmerge) {
        scatter (i in range(length(CreateBatches.vcf_batch_fofns))) {
            scatter (j in range(length(regions))) {
                call Ivcfmerge as IvcfmergeSingleBatchRegion {
                    input:
                        vcfs = read_lines(CreateBatches.vcf_batch_fofns[i]),
                        vcf_idxs = read_lines(CreateBatches.vcf_idx_batch_fofns[i]),
                        output_prefix = output_prefix + ".batch-" + i + ".region-" + i,
                        sample_names = select_first([sample_names]),
                        region_args = "-r " + regions[j]
                }
            }
        }
    }
    if (!use_ivcfmerge) {
        scatter (i in range(length(CreateBatches.vcf_batch_fofns))) {
            scatter (j in range(length(regions))) {
                call MergeVcfs as MergeVcfsSingleBatchRegion {
                    input:
                        vcfs = read_lines(CreateBatches.vcf_batch_fofns[i]),
                        vcf_idxs = read_lines(CreateBatches.vcf_idx_batch_fofns[i]),
                        output_prefix = output_prefix + ".batch-" + i + ".region-" + i,
                        extra_args = "-r " + regions[j] + " " + extra_merge_args
                }
            }
        }
    }

    Array[Array[File]] region_by_batch_vcfs = transpose(select_first([IvcfmergeSingleBatchRegion.merged_vcf, MergeVcfsSingleBatchRegion.merged_vcf]))
    Array[Array[File]] region_by_batch_vcf_idxs = transpose(select_first([IvcfmergeSingleBatchRegion.merged_vcf_idx, MergeVcfsSingleBatchRegion.merged_vcf_idx]))

    if (use_ivcfmerge) {
        # merge all samples in each region
        scatter (j in range(length(regions))) {
            call Ivcfmerge as IvcfmergeSingleRegion {
                input:
                    vcfs = region_by_batch_vcfs[j],
                    vcf_idxs = region_by_batch_vcf_idxs[j],
                    output_prefix = output_prefix + ".region-" + j,
                    sample_names = select_first([sample_names])
            }
        }
    }
    if (!use_ivcfmerge) {
        # merge all samples in each region
        scatter (j in range(length(regions))) {
            call MergeVcfs as MergeVcfsSingleRegion {
                input:
                    vcfs = region_by_batch_vcfs[j],
                    vcf_idxs = region_by_batch_vcf_idxs[j],
                    output_prefix = output_prefix + ".region-" + j,
                    extra_args = extra_merge_args
            }
        }
    }

    # concatenate all regions
    call ConcatVcfs {
        input:
            vcfs = select_first([IvcfmergeSingleRegion.merged_vcf, MergeVcfsSingleRegion.merged_vcf]),
            vcf_idxs = select_first([IvcfmergeSingleRegion.merged_vcf_idx, MergeVcfsSingleRegion.merged_vcf_idx]),
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

    Int disk_gb = 10 + 2 * ceil(size([vcfs], "GiB"))

    command <<<
        set -euox pipefail

        bcftools merge \
            -l ~{write_lines(vcfs)} \
            ~{extra_args} \
            -Oz -o ~{output_prefix}.vcf.gz
        bcftools index -t ~{output_prefix}.vcf.gz
    >>>

    output {
        File merged_vcf = "~{output_prefix}.vcf.gz"
        File merged_vcf_idx = "~{output_prefix}.vcf.gz.tbi"
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

# assumes all VCFs have identical variants
# TODO make this work for identical filenames
task Ivcfmerge {
    input{
        Array[File] vcfs
        Array[File] vcf_idxs
        Array[String] sample_names
        String output_prefix
        String? region_args

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size([vcfs], "GiB"))

    command <<<
        set -euox pipefail

        wget https://github.com/iqbal-lab-org/ivcfmerge/archive/refs/tags/v1.0.0.tar.gz
        tar -xvf v1.0.0.tar.gz

        mkdir compressed
        mv ~{sep=' ' vcfs} compressed
        mv ~{sep=' ' vcf_idxs} compressed

        if [ $(ls compressed/*.vcf.gz | wc -l) == 1 ]
        then
            cp $(ls compressed/*.vcf.gz) ~{output_prefix}.vcf.gz
            cp $(ls compressed/*.vcf.gz.tbi) ~{output_prefix}.vcf.gz.tbi
        else
            mkdir decompressed
            ls compressed/*.vcf.gz | xargs -I % sh -c 'bcftools annotate --no-version ~{region_args} -x INFO % --threads 2 -Ov -o decompressed/$(basename % .gz)'
            time python ivcfmerge-1.0.0/ivcfmerge.py <(ls decompressed/*.vcf) ~{output_prefix}.vcf
            bcftools annotate --no-version -S ~{write_lines(sample_names)} -x FORMAT/FT ~{output_prefix}.vcf --threads 2 -Oz -o ~{output_prefix}.vcf.gz
            bcftools index -t ~{output_prefix}.vcf.gz
        fi
    >>>

    output {
        File merged_vcf = "~{output_prefix}.vcf.gz"
        File merged_vcf_idx = "~{output_prefix}.vcf.gz.tbi"
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
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.11"
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

    Int disk_gb = 10 + 2 * ceil(size([vcfs], "GiB"))

    command <<<
        set -euox pipefail

        bcftools concat \
            -f ~{write_lines(vcfs)} \
            ~{extra_args} \
            -Oz -o ~{output_prefix}.vcf.gz
        bcftools index -t --threads $(nproc) ~{output_prefix}.vcf.gz
    >>>

    output {
        File concatenated_vcf = "~{output_prefix}.vcf.gz"
        File concatenated_vcf_idx = "~{output_prefix}.vcf.gz.tbi"
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
