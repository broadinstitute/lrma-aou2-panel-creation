version 1.0

# Gather collisionless shards, create a sites-only VCF, and create Shapeit4 chunks.

workflow GatherSitesAndChunk {

    input {
        String region
        Array[File] filter_and_concat_vcfs
        Array[File] filter_and_concat_vcf_idxs
        Array[File] collisionless_vcfs
        Array[File] collisionless_vcf_idxs
        String entity_name
        String output_prefix

        String chunk_extra_args = "--thread $(nproc) --window-size 1000000 --buffer-size 200000 --window-count 50000 --buffer-count 500" # we want window counts to drive the constraints
    }

    call BcftoolsConcatNaive as FilterAndConcatConcat { input:
        vcfs = filter_and_concat_vcfs,
        vcf_idxs = filter_and_concat_vcf_idxs,
        output_prefix = output_prefix + ".filterAndConcat"
    }

    call BcftoolsConcatNaive as CollisionlessPreShapeitConcat { input:
        vcfs = collisionless_vcfs,
        vcf_idxs = collisionless_vcf_idxs,
        output_prefix = output_prefix + ".collisionless"
    }

    call CreateSitesOnlyVCF as SitesOnly { input:
        vcf = FilterAndConcatConcat.concatenated_vcf,
        vcf_idx = FilterAndConcatConcat.concatenated_vcf_idx,
        output_prefix = output_prefix + ".collisionless.sites"
    }

    call CreateShapeitChunks { input:
        vcf = SitesOnly.sites_only_vcf,
        vcf_idx = SitesOnly.sites_only_vcf_idx,
        region = region,
        entity_name = entity_name,
        output_prefix = output_prefix,
        extra_args = chunk_extra_args
    }

    output {
        File filter_and_concat_vcf = FilterAndConcatConcat.concatenated_vcf
        File filter_and_concat_vcf_idx = FilterAndConcatConcat.concatenated_vcf_idx
        File collisionless_vcf = CollisionlessPreShapeitConcat.concatenated_vcf
        File collisionless_vcf_idx = CollisionlessPreShapeitConcat.concatenated_vcf_idx
        File sites_only_vcf = SitesOnly.sites_only_vcf
        File sites_only_vcf_idx = SitesOnly.sites_only_vcf_idx
        File common_chunks = CreateShapeitChunks.common_chunks
        File rare_chunks = CreateShapeitChunks.rare_chunks
        File chunks_terra_tsv = CreateShapeitChunks.chunks_terra_tsv
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    String? disk_type
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

        bcftools concat ~{sep=" " vcfs} --naive -Ob -o ~{output_prefix}.bcf
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
        disk_type:          "SSD",
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.23"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type])
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

        bcftools view --threads $(nproc) ~{vcf} -G -Ob -o ~{output_prefix}.bcf
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
        disk_type:          "SSD",
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.23"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type])
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
        String output_prefix
        String region
        String entity_name
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
        cut -f 3 chunks.tsv > ~{output_prefix}.chunks.regions.common.txt
        cut -f 4 chunks.tsv > ~{output_prefix}.chunks.regions.rare.txt
        
        # Generate Terra-compatible TSV
        python3 <<EOF
        import sys

        in_file = "chunks.tsv"
        out_file = "~{output_prefix}.chunks.tsv"

        with open(in_file, "r") as f_in, open(out_file, "w") as f_out:
            lines = [line.strip().split("\t") for line in f_in if line.strip()]

            # GLIMPSE chunk output has 6 columns
            header = [f"entity:~{entity_name}_id", "index", "chrom", "outer_region", "inner_region", "inner_region_length", "number_of_variants"]

            f_out.write("\t".join(header) + "\n")

            for i, row in enumerate(lines):
                # Construct an entity ID using the outer region
                outer_region = row[2].replace(':', '-')
                entity_id = f"~{output_prefix}.shard-{i:04d}.{outer_region}"
                f_out.write("\t".join([entity_id] + row) + "\n")
        EOF
    >>>

    output {
        File common_chunks = "~{output_prefix}.chunks.regions.common.txt"
        File rare_chunks = "~{output_prefix}.chunks.regions.rare.txt"
        File chunks_terra_tsv = "~{output_prefix}.chunks.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "SSD",
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.11"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type])
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
