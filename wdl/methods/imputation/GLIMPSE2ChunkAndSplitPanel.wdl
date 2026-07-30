version 1.0


workflow GLIMPSE2ChunkAndSplitPanel {
    input {
        Array[String] chromosomes
        File genetic_maps_tsv
        
        # per chromosome, in same order
        Array[File] panel_bubble_split_vcfs          # "split" here means "split to biallelic"; "split" just below means "chunked"
        Array[File] panel_bubble_split_vcf_idxs
        Array[File] panel_bubble_split_sites_only_vcfs
        Array[File] panel_bubble_split_sites_only_vcf_idxs

        String extra_chunk_args = "--thread $(nproc) --window-mb 5 --buffer-mb 0.5 --sequential"
        String extra_split_args = "--keep-monomorphic-ref-sites"
        String output_prefix

        String glimpse2_docker = "us.gcr.io/broad-gotc-prod/imputation-glimpse2:1.0.0-2cee597-1778869818"    # enables checkpointing, but note this contains bcftools/htslib 1.16!
    }

    Map[String, String] genetic_maps_dict = read_map(genetic_maps_tsv)

    scatter (i in range(length(chromosomes))) {
        String chromosome = chromosomes[i]

        call GLIMPSE2Chunk {
            input:
                vcf = panel_bubble_split_sites_only_vcfs[i],
                vcf_idx = panel_bubble_split_sites_only_vcf_idxs[i],
                region = chromosome,
                genetic_map = genetic_maps_dict[chromosome],
                output_prefix = output_prefix + "." + chromosome,
                extra_chunk_args = extra_chunk_args,
                docker = glimpse2_docker
        }

        Array[String] input_regions = read_lines(GLIMPSE2Chunk.input_regions)
        Array[String] output_regions = read_lines(GLIMPSE2Chunk.output_regions)

        scatter (k in range(length(output_regions))) {
            call GLIMPSE2SplitReference as ChunkedGLIMPSE2SplitReference {
                input:
                    panel_bubble_split_vcf = panel_bubble_split_vcfs[i],
                    panel_bubble_split_vcf_idx = panel_bubble_split_vcf_idxs[i],
                    input_region = input_regions[k],
                    output_region = output_regions[k],
                    genetic_map = genetic_maps_dict[chromosome],
                    output_prefix = output_prefix + "." + chromosome + ".shard-" + k + ".split",
                    extra_split_args = extra_split_args,
                    docker = glimpse2_docker
            }
        }

        ChunkedPanelChromosome chunked_panel_chromosome = object {
            input_regions: input_regions,
            output_regions: output_regions,
            panel_split_chunk_bins: ChunkedGLIMPSE2SplitReference.panel_split_chunk_bin
        }
        Pair[String, ChunkedPanelChromosome] chunked_panel_chromosome_pair = (chromosome, chunked_panel_chromosome)
    }

    call CoercePairsToMap {
        input:
            pair_array = chunked_panel_chromosome_pair,
            output_prefix = output_prefix
    }

    output {
        File chunked_panel_json = CoercePairsToMap.out_map_json
        Map[String, ChunkedPanelChromosome] chunked_panel = read_json(CoercePairsToMap.out_map_json)
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

struct ChunkedPanelChromosome {
    Array[String] input_regions
    Array[String] output_regions
    Array[String] panel_split_chunk_bins
}

task GLIMPSE2Chunk {
    input {
        File vcf
        File vcf_idx
        String region
        File genetic_map
        String output_prefix
        String? extra_chunk_args
        String docker

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 2 * ceil(size(vcf, "GB"))

    command <<<
        set -euxo pipefail

        /bin/GLIMPSE2_chunk \
            -I ~{vcf} \
            --region ~{region} \
            --map ~{genetic_map} \
            ~{extra_chunk_args} \
            -O ~{output_prefix}.chunks.tsv

        # cut chunks + buffers
        cut -f 3 ~{output_prefix}.chunks.tsv > ~{output_prefix}.input-regions.tsv
        cut -f 4 ~{output_prefix}.chunks.tsv > ~{output_prefix}.output-regions.tsv
    >>>

    output {
        File input_regions = "~{output_prefix}.input-regions.tsv"
        File output_regions = "~{output_prefix}.output-regions.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             docker
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

task GLIMPSE2SplitReference {
    input {
        File panel_bubble_split_vcf
        File panel_bubble_split_vcf_idx
        String input_region
        String output_region
        File genetic_map
        String output_prefix
        String? extra_split_args
        String docker

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 2 * ceil(size(panel_bubble_split_vcf, "GB"))

    command <<<
        set -euxo pipefail

        /bin/GLIMPSE2_split_reference \
            -R ~{panel_bubble_split_vcf} \
            --input-region ~{input_region} \
            --output-region ~{output_region} \
            --map ~{genetic_map} \
            --thread $(nproc) \
            ~{extra_split_args} \
            --output ~{output_prefix}
    >>>

    output {
        File panel_split_chunk_bin = glob("~{output_prefix}_*bin")[0]       # TODO parse input region and construct ~{output_prefix}_chr_start_end.bin filename
    }
    
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             docker
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

task CoercePairsToMap {
    input {
        Array[Pair[String, ChunkedPanelChromosome]] pair_array
        String output_prefix
        
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        python3 <<CODE
        import json

        # write_json is safe here because it executes inside the task container
        with open("~{write_json(pair_array)}", "r") as f:
            pairs = json.load(f)

        out_map = {item["left"]: item["right"] for item in pairs}

        with open("~{output_prefix}.json", "w") as f:
            json.dump(out_map, f, indent=2)
        CODE
    >>>

    output {
        File out_map_json = "~{output_prefix}.json"
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
