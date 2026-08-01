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

        call CountPanelVariantsPerShard {
            input:
                panel_bubble_split_sites_only_vcf = panel_bubble_split_sites_only_vcfs[i],
                panel_bubble_split_sites_only_vcf_idx = panel_bubble_split_sites_only_vcf_idxs[i],
                input_regions = input_regions,
                output_prefix = output_prefix + "." + chromosome,
                docker = glimpse2_docker
        }
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
            chunks_tsv: GLIMPSE2Chunk.chunks_tsv,
            input_regions: input_regions,
            output_regions: output_regions,
            panel_split_chunk_bins: ChunkedGLIMPSE2SplitReference.panel_split_chunk_bin,
            n_variants: CountPanelVariantsPerShard.n_variants
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
    String chunks_tsv
    Array[String] input_regions
    Array[String] output_regions
    Array[String] panel_split_chunk_bins
    # Panel variants in each shard's input (buffered) region -- GLIMPSE2's L, parallel to
    # input_regions. Computed once when the panel is chunked, so every batch reads it instead
    # of recounting the same numbers. NOT column 7 of chunks.tsv, which covers a different
    # region: 403,101 against an actual 965,039 for chr20 shard 4.
    Array[Int] n_variants
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
        File chunks_tsv = "~{output_prefix}.chunks.tsv"
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

# Counts panel variants per shard input (buffered) region -- exactly GLIMPSE2's L. Verified
# against its logged L on three shards: 657,784 / 965,039 / 1,346,888, all exact. NOT column 7
# of chunks.tsv, which covers a different region (403,101 vs an actual 965,039 for chr20 s4).
# One task per chromosome, run once when the panel is chunked. Every batch then reads the
# numbers out of chunked_panel.json instead of recounting them, which is what lets phase size
# each shard from its own L rather than from the densest one in the genome.
task CountPanelVariantsPerShard {
    input {
        File panel_bubble_split_sites_only_vcf
        File panel_bubble_split_sites_only_vcf_idx
        Array[String] input_regions
        String output_prefix

        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 2 * ceil(size(panel_bubble_split_sites_only_vcf, "GB"))


    #########################
    # 8 GiB, up from 4: the measured 2.49 GiB peak was a single reader and the loop runs
    # eff_cpu of them concurrently. Memory at that concurrency is unmeasured -- the readers
    # stream and should share page cache -- and the instrumentation below reports it.
    # MEASURED on chr20 (7 regions): 65 s wall, 2.49 GiB peak, 1 GB of disk.
    #
    # NOT preemptible: the whole panel build waits on these, and they cost cents.
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       0,
        use_ssd:            true,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             docker
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    # runtime_attr_override applies field by field, so overriding only mem_gb would keep a cpu
    # computed for the DEFAULT memory. Re-derive from whichever memory won; an explicit
    # cpu_cores still wins. Same pattern as the other tasks in this file.
    Float eff_mem_gb  = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
    Int eff_ratio_cpu = ceil(eff_mem_gb / 6.5)
    Int eff_unrounded = if eff_ratio_cpu > 4 then eff_ratio_cpu else 4
    Int eff_cpu       = select_first([runtime_attr.cpu_cores, eff_unrounded + (eff_unrounded % 2)])
    command <<<
        set -euxo pipefail

        # Peak RSS via EXIT trap on stderr. Must be a trap: under set -e a SIGKILLed child
        # aborts the script, so an end-of-script report would miss every OOM -- the only runs
        # worth measuring. Must be stderr: Cromwell delocalizes it for FAILED tasks, File
        # outputs do not exist. Guarded throughout; cannot fail the task it measures.
        # $1 may be the literal "max" (cgroup v2, no limit); emit NA rather than 0.00.
        _instr_gib() { for f in "$@"; do if [ -r "$f" ]; then awk '$1 ~ /^[0-9]+$/ {printf "%.2f", $1/1073741824; ok=1} END{if(!ok) printf "NA"}' "$f" 2>/dev/null && return 0; fi; done; printf NA; }
        _instr_peak()  { _instr_gib /sys/fs/cgroup/memory.peak /sys/fs/cgroup/memory/memory.max_usage_in_bytes; }
        _instr_limit() { _instr_gib /sys/fs/cgroup/memory.max  /sys/fs/cgroup/memory/memory.limit_in_bytes; }
        _instr_cur()   { _instr_gib /sys/fs/cgroup/memory.current /sys/fs/cgroup/memory/memory.usage_in_bytes; }
        _INSTR_T0=$SECONDS
        _INSTR_SAMPLER=""
        _instr_report() {
            _rc=$?
            set +e +x
            [ -n "$_INSTR_SAMPLER" ] && kill "$_INSTR_SAMPLER" 2>/dev/null
            # -Pk is POSIX; -BG is GNU-only and absent on busybox. 1K blocks -> GiB in awk.
            _du=$(df -Pk . 2>/dev/null | awk 'NR==2{printf "%.0f %.0f", $3/1048576, $2/1048576}')
            echo "[RESOURCE] task=CountPanelVariantsPerShard n_regions=~{length(input_regions)} requested_mem_gib=~{eff_mem_gb} requested_cpu=~{eff_cpu} rc=$_rc peak_rss_gib=$(_instr_peak) limit_gib=$(_instr_limit) disk_used_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f1) disk_total_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f2) wall_s=$((SECONDS-_INSTR_T0))" >&2
        }
        trap _instr_report EXIT
        # 10 s time series, so a spike can be located against GLIMPSE2's Cnk/Buf markers.
        ( set +x; while :; do
            echo "[RESOURCE-TS] t=$((SECONDS-_INSTR_T0)) rss_gib=$(_instr_cur)" >&2
            sleep 10
          done ) & _INSTR_SAMPLER=$!

        # Indexed read per region. The regions overlap, so this reads ~1.5x the file total,
        # and a miscount mis-sizes a
        # shard silently. Run in parallel because every phase shard for this chromosome blocks
        # on this task -- chr20's 7 regions took 65 s and chr2 has 47, so serial execution
        # would put ~7 minutes on the critical path of every chromosome. Parallel is also
        # cheaper as well as faster: cpu bills linearly but the memory reservation is charged
        # for the whole wall time, so finishing 4x sooner on 4x the cores costs slightly less.
        #
        # Counts go to per-index files and are assembled in order afterwards, because parallel
        # completion order is not region order and the array must line up with input_regions.
        # A per-job success sentinel distinguishes "bcftools failed" from "region genuinely has
        # no variants": `wait` with no arguments returns 0 regardless of what the background
        # jobs did, so without it a tool failure would surface later as a spurious empty-region
        # error pointing at the wrong cause.
        PAR=~{eff_cpu}
        REGIONS_FILE=~{write_lines(input_regions)}
        awk '{print NR-1 "\t" $0}' "$REGIONS_FILE" > indexed_regions.txt
        while IFS=$'\t' read -r IDX REGION; do
            (
                set -o pipefail
                # --regions-overlap 0 matches GLIMPSE2: a record counts if its POS is inside.
                # The default mode also counts records starting before the region, inflating it.
                bcftools view --no-version --threads 1 -H -r "$REGION" --regions-overlap 0 \
                    ~{panel_bubble_split_sites_only_vcf} \
                    | awk 'END{print NR+0}' > "count.$IDX" && touch "ok.$IDX"
            ) &
            # Exclude the instrumentation sampler, which is also a background job of this
            # shell. Counting it caps concurrency at PAR-1, and at PAR=1 -- reachable by a
            # cpu_cores override -- the sampler alone satisfies the condition and this loop
            # spins forever, hanging a non-preemptible task that every phase shard waits on.
            while [ "$(jobs -rp | grep -vx "${_INSTR_SAMPLER:-}" | awk 'END{print NR+0}')" -ge "$PAR" ]; do
                sleep 1
            done
        done < indexed_regions.txt
        wait

        # awk, because wc -l pads its output with spaces on some platforms, which would
        # break the numeric checks below.
        EXPECTED=$(awk 'END{print NR+0}' "$REGIONS_FILE")
        : > counts.txt
        for IDX in $(seq 0 $((EXPECTED-1))); do
            [ -f "ok.$IDX" ] || { echo "ERROR: bcftools failed on region index $IDX." >&2; exit 1; }
            grep -qE '^[0-9]+$' "count.$IDX" || { echo "ERROR: non-numeric count at index $IDX." >&2; exit 1; }
            cat "count.$IDX" >> counts.txt
        done

        if grep -qx '0' counts.txt; then
            echo "ERROR: at least one input region contains no panel variants." >&2
            paste -d' ' "$REGIONS_FILE" counts.txt >&2
            exit 1
        fi

        # read_json wants a JSON array; read_lines gives Array[String], which does not coerce.
        printf '[%s]\n' "$(paste -sd, counts.txt)" > ~{output_prefix}.n_variants.json
        cat ~{output_prefix}.n_variants.json
    >>>

    output {
        Array[Int] n_variants = read_json("~{output_prefix}.n_variants.json")
        File n_variants_json = "~{output_prefix}.n_variants.json"
    }
    runtime {
        cpu:                    eff_cpu
        memory:                 eff_mem_gb + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
