version 1.0

# Gather collisionless shards, create a sites-only VCF, and create Shapeit4 chunks.
# TODO create Terra data table from chunks

workflow GatherSitesAndChunk {

    input {
        String region
        Array[File] collisionless_vcfs
        Array[File] collisionless_vcf_idxs
        String output_prefix

        String chunk_extra_args = "--thread $(nproc) --window-size 1000000 --buffer-size 200000 --window-count 50000 --buffer-count 500" # we want window counts to drive the constraints
    }

    call BcftoolsConcatNaive as CollisionlessPreShapeitConcat { input:
        vcfs = collisionless_vcfs,
        vcf_idxs = collisionless_vcf_idxs,
        output_prefix = output_prefix + ".collisionless"
    }

    call CreateSitesOnlyVCF as CollisionlessPreShapeitSitesOnly { input:
        vcf = CollisionlessPreShapeitConcat.concatenated_vcf,
        vcf_idx = CollisionlessPreShapeitConcat.concatenated_vcf_idx,
        output_prefix = output_prefix + ".collisionless.sites"
    }

    call CreateShapeitChunks { input:
        vcf = CollisionlessPreShapeitSitesOnly.sites_only_vcf,
        vcf_idx = CollisionlessPreShapeitSitesOnly.sites_only_vcf_idx,
        region = region,
        extra_args = chunk_extra_args
    }

    output {
        File collisionless_vcf = CollisionlessPreShapeitConcat.concatenated_vcf
        File collisionless_vcf_idx = CollisionlessPreShapeitConcat.concatenated_vcf_idx
        File sites_only_vcf = CollisionlessPreShapeitSitesOnly.sites_only_vcf
        File sites_only_vcf_idx = CollisionlessPreShapeitSitesOnly.sites_only_vcf_idx
        File chunks_tsv = CreateShapeitChunks.chunks_tsv
        File common_chunks = CreateShapeitChunks.common_chunks
        File rare_chunks = CreateShapeitChunks.rare_chunks
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

task BcftoolsConcatNaive {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 50 + 4 * ceil(size(vcfs, "GiB"))

    command <<<
        set -euxo pipefail

        bcftools concat --no-version ~{sep=" " vcfs} --naive -Ob -o ~{output_prefix}.bcf
        bcftools index ~{output_prefix}.bcf
    >>>

    output {
        File concatenated_vcf = "~{output_prefix}.bcf"
        File concatenated_vcf_idx = "~{output_prefix}.bcf.csi"
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

task CreateSitesOnlyVCF {
    input {
        File vcf
        File vcf_idx
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + ceil(size(vcf, "GiB"))

    command <<<
        set -euxo pipefail

        bcftools view --no-version --threads $(nproc) ~{vcf} -G -Ob -o ~{output_prefix}.bcf
        bcftools index ~{output_prefix}.bcf
    >>>

    output {
        File sites_only_vcf = "~{output_prefix}.bcf"
        File sites_only_vcf_idx = "~{output_prefix}.bcf.csi"
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

task CreateShapeitChunks {
    input {
        File vcf
        File vcf_idx
        String region
        String extra_args = "--thread $(nproc) --window-size 1000000 --buffer-size 200000 --window-count 50000 --buffer-count 500" # we want window counts to drive the constraints

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size([vcf, vcf_idx], "GiB")) + 1

    command <<<
        set -euxo pipefail

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_chunk_static
        chmod +x GLIMPSE_chunk_static

        ./GLIMPSE_chunk_static \
            -I ~{vcf} \
            --region ~{region} \
            ~{extra_args} \
            -O chunks.tsv

        # cut chunks + buffers
        cut -f 3 chunks.tsv > common.chunks.regions.txt
        cut -f 4 chunks.tsv > rare.chunks.regions.txt
    >>>

    output {
        File chunks_tsv = "chunks.tsv"
        File common_chunks = "common.chunks.regions.txt"
        File rare_chunks = "rare.chunks.regions.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
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
