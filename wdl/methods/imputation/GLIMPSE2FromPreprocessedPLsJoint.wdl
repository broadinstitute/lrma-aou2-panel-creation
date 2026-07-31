version 1.0

import "../ConcatVcfs.wdl" as ConcatVcfs

workflow GLIMPSE2FromPreprocessedPLsJoint {
    input {
        File input_preprocessed_joint_vcf
        File input_preprocessed_joint_vcf_idx

        File? remap_sample_names_file    # TSV with old_name new_name mappings

        String chromosome
        File genetic_maps_tsv
        File chunked_panel_json

        String extra_phase_args = "--impute-reference-only-variants --keep-monomorphic-ref-sites --main 10 --burnin 5 --err-imp 1E-3"

        # Workflow-level, not call-qualified: the phase memory request is computed from both.
        Int phase_threads = 4
        Int phase_kpbwt = 1000

        String output_prefix

        # inputs for PopAndMarginalizeCollisions
        File pop_glimpse2_panel_resources_json
        File? pop_glimpse2_script               # heavily modified version of convert-to-biallelic.py
        File? pop_glimpse2_cargo_toml
        File? pop_glimpse2_binary

        String glimpse2_docker = "us.gcr.io/broad-gotc-prod/imputation-glimpse2:1.0.0-2cee597-1778869818"    # enables checkpointing, but note this contains bcftools/htslib 1.16!

        # Byte-identical to Cromwell's own ZonesDefaultValue: this RESTORES the stock four
        # zones, which the VWB backend overrode to pin us-central1-a (71% of 502 attempts were
        # preempted there). ASK VWB FIRST -- pinning may be a data-locality or compliance
        # decision, not an oversight. Set to "us-central1-a" to restore the pin.
        # ZonesValidation accepts String and Array[String] alike, so the form is style only.
        # Values must lie in the backend's Batch job region. Does not reach the imported
        # ConcatVcfs call, which takes no zones argument.
        # UNVERIFIED: whether a task-level value overrides the backend's allowedLocations.
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }

    Map[String, String] genetic_maps_dict = read_map(genetic_maps_tsv)
    String genetic_map = genetic_maps_dict[chromosome]

    Map[String, ChunkedPanelChromosome] chunked_panel = read_json(chunked_panel_json)
    Array[String] input_regions = chunked_panel[chromosome].input_regions
    Array[String] output_regions = chunked_panel[chromosome].output_regions
    Array[File] panel_split_chunk_bins = chunked_panel[chromosome].panel_split_chunk_bins

    Map[String, PopAndMarginalizePanelResourcesChromosome] pop_glimpse2_panel_resources = read_json(pop_glimpse2_panel_resources_json)
    File panel_bubble_split_sites_only_vcf = pop_glimpse2_panel_resources[chromosome].panel_bubble_split_sites_only_vcf
    File panel_bubble_split_sites_only_vcf_idx = pop_glimpse2_panel_resources[chromosome].panel_bubble_split_sites_only_vcf_idx
    File panel_id_split_vcf_gz = pop_glimpse2_panel_resources[chromosome].panel_id_split_vcf_gz
    File panel_id_split_vcf_gz_tbi = pop_glimpse2_panel_resources[chromosome].panel_id_split_vcf_gz_tbi
    Array[String] pop_regions = select_first([pop_glimpse2_panel_resources[chromosome].pop_regions, output_regions])
    

    # Per-shard L, so phase is sized per shard rather than for the worst one. Derived in-run
    # from this chromosome's panel over these regions in order, so it cannot be stale,
    # reordered or cross-chromosome -- which a JSON field of the right length could all be.
    call CountPanelVariantsPerShard {
        input:
            panel_bubble_split_sites_only_vcf = panel_bubble_split_sites_only_vcf,
            panel_bubble_split_sites_only_vcf_idx = panel_bubble_split_sites_only_vcf_idx,
            input_regions = input_regions,
            output_prefix = output_prefix,
            zones = zones,
            docker = glimpse2_docker
    }

    scatter (k in range(length(output_regions))) {
        call GLIMPSE2Phase as ChunkedGLIMPSE2Phase {
            input:
                input_vcf = input_preprocessed_joint_vcf,
                input_vcf_idx = input_preprocessed_joint_vcf_idx,
                panel_split_chunk_bin = panel_split_chunk_bins[k],
                input_region = input_regions[k],
                output_region = output_regions[k],
                genetic_map = genetic_map,
                output_prefix = output_prefix + ".shard-" + k + ".glimpse2.phased",
                extra_phase_args = extra_phase_args,
                phase_threads = phase_threads,
                phase_kpbwt = phase_kpbwt,
                n_variants = CountPanelVariantsPerShard.n_variants[k],
                zones = zones,
                docker = glimpse2_docker
        }
    }

    call GLIMPSE2Ligate {
        input:
            phased_vcfs = ChunkedGLIMPSE2Phase.phased_vcf,
            phased_vcf_idxs = ChunkedGLIMPSE2Phase.phased_vcf_idx,
            output_prefix = output_prefix + ".glimpse2.bubble",
            zones = zones,
            docker = glimpse2_docker
    }

    scatter (k in range(length(pop_regions))) {
        call PopAndMarginalizeCollisions { input:
            posteriors_vcf = GLIMPSE2Ligate.ligated_vcf,
            posteriors_vcf_idx = GLIMPSE2Ligate.ligated_vcf_idx,
            panel_bubble_split_sites_only_vcf = panel_bubble_split_sites_only_vcf,
            panel_bubble_split_sites_only_vcf_idx = panel_bubble_split_sites_only_vcf_idx,
            panel_id_split_vcf_gz = panel_id_split_vcf_gz,
            panel_id_split_vcf_gz_tbi = panel_id_split_vcf_gz_tbi,
            pop_glimpse2_script = pop_glimpse2_script,
            cargo_toml = pop_glimpse2_cargo_toml,
            pop_glimpse2_binary = pop_glimpse2_binary,
            region = pop_regions[k],
            zones = zones,
            output_prefix = output_prefix + ".glimpse2.popped"
        }
    }
    
    call ConcatVcfs.ConcatVcfs as ConcatPopAndMarginalizeCollisions { input:
        vcfs = PopAndMarginalizeCollisions.popped_vcf,
        vcf_idxs = PopAndMarginalizeCollisions.popped_vcf_idx,
        output_prefix = output_prefix + ".glimpse2.popped",
        do_bcf = true,
        do_sort = false,
        extra_args = "--threads $(nproc) --naive",
        regions = [],
        do_sort_shard = false,
        extra_args_shard = ""
    }

    # Conditionally trigger remapping tasks
    if (defined(remap_sample_names_file)) {
        call RemapSampleNames as RemapBubblePosteriors {
            input:
                vcf = GLIMPSE2Ligate.ligated_vcf,
                vcf_idx = GLIMPSE2Ligate.ligated_vcf_idx,
                remap_file = select_first([remap_sample_names_file]),
                zones = zones,
                output_prefix = output_prefix + ".glimpse2.bubble"
        }

        call RemapSampleNames as RemapPoppedPosteriors {
            input:
                vcf = ConcatPopAndMarginalizeCollisions.concatenated_vcf,
                vcf_idx = ConcatPopAndMarginalizeCollisions.concatenated_vcf_idx,
                remap_file = select_first([remap_sample_names_file]),
                zones = zones,
                output_prefix = output_prefix + ".glimpse2.popped"
        }
    }

    output {
        File glimpse2_bubble_posteriors_vcf = select_first([RemapBubblePosteriors.output_vcf, GLIMPSE2Ligate.ligated_vcf])
        File glimpse2_bubble_posteriors_vcf_idx = select_first([RemapBubblePosteriors.output_vcf_idx, GLIMPSE2Ligate.ligated_vcf_idx])
        File glimpse2_popped_posteriors_vcf = select_first([RemapPoppedPosteriors.output_vcf, ConcatPopAndMarginalizeCollisions.concatenated_vcf])
        File glimpse2_popped_posteriors_vcf_idx = select_first([RemapPoppedPosteriors.output_vcf_idx, ConcatPopAndMarginalizeCollisions.concatenated_vcf_idx])
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    # NOT consumed here, deliberately. Cromwell ADDS its 30 GB default to any bootDiskSizeGb
    # request, so the old value of 10 provisioned 40; omitting it yields the 30 GB minimum.
    # An explicit 0 would keep this override live, but it is unverified that
    # BootDiskSizeValidation accepts 0 -- and that would fail every task at once.
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

    # No per-shard variant count here, deliberately. These are unkeyed parallel arrays, so a
    # count that was stale but happened to have the right length would pair silently with the
    # wrong bins -- and underprovisioned memory is the failure the sizing exists to prevent.
    # It would also be dead: GLIMPSE2ChunkAndSplitPanel never emits such a field.
    # CountPanelVariantsPerShard derives the counts at runtime instead, from this chromosome's
    # panel over these regions in order -- correct by construction.
}

struct PopAndMarginalizePanelResourcesChromosome {
    String panel_bubble_split_sites_only_vcf
    String panel_bubble_split_sites_only_vcf_idx
    String panel_id_split_vcf_gz
    String panel_id_split_vcf_gz_tbi
    Array[String]? pop_regions              # non-overlapping, if not provided then GLIMPSE2 chunks will be used
}

# Counts panel variants per shard input (buffered) region -- exactly GLIMPSE2's L. Verified
# against its logged L on three shards: 657,784 / 965,039 / 1,346,888, all exact. NOT column 7
# of chunks.tsv, which covers a different region (403,101 vs an actual 965,039 for chr20 s4).
# One task per chromosome, cents, against ~$20/batch saved by not sizing for the worst shard.
task CountPanelVariantsPerShard {
    input {
        File panel_bubble_split_sites_only_vcf
        File panel_bubble_split_sites_only_vcf_idx
        Array[String] input_regions
        String output_prefix

        String docker
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 2 * ceil(size(panel_bubble_split_sites_only_vcf, "GB"))

    command <<<
        set -euxo pipefail

        # --regions-overlap 0 matches GLIMPSE2: a record counts if its POS is inside. The
        # default mode would also include records starting before the region and inflate.
        # Indexed read per region, not one bucketed stream: regions overlap so this reads
        # ~1.5x the file, and a miscount mis-sizes a shard silently. Correctness over speed.
        REGIONS_FILE=~{write_lines(input_regions)}
        : > counts.txt
        while read -r REGION; do
            bcftools view --no-version -H -r "$REGION" --regions-overlap 0 \
                ~{panel_bubble_split_sites_only_vcf} | wc -l >> counts.txt
        done < "$REGIONS_FILE"

        # Fail loudly rather than silently mis-sizing shards if a region yielded nothing.
        EXPECTED=$(wc -l < "$REGIONS_FILE")
        ACTUAL=$(wc -l < counts.txt)
        if [ "$ACTUAL" -ne "$EXPECTED" ]; then
            echo "ERROR: counted $ACTUAL regions, expected $EXPECTED." >&2
            exit 1
        fi
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

    #########################
    # NOT preemptible: every phase shard blocks on this, so it is a serialisation point and a
    # per-chromosome single point of failure. Minutes of runtime against a 7-minute median
    # preemption would often restart and delay all 523 shards. It costs cents.
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size_gb,
        use_ssd:            true,
        preemptible_tries:  0,
        max_retries:        1,
        docker:             docker
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        zones:                  zones
    }
}

# checkpoint implementation borrowed from https://github.com/broadinstitute/palantir-workflows/blob/main/GlimpseImputationPipeline/Glimpse2Imputation.wdl
task GLIMPSE2Phase {
    input {
        File input_vcf
        File input_vcf_idx
        File panel_split_chunk_bin
        String input_region
        String output_region
        File genetic_map
        String output_prefix
        String? extra_phase_args = "--impute-reference-only-variants --keep-monomorphic-ref-sites --main 10 --burnin 5 --err-imp 1E-3"

        String docker

        # Typed, not text in extra_phase_args, because the memory request is computed from
        # them: a caller replacing that string would otherwise drop --Kpbwt 1000, get this
        # build's default of 2000, and double the matrix while the request stayed put.
        # threads is pinned rather than $(nproc) because the matrix is per-thread, which under
        # $(nproc) makes memory depend on cpu while the N1 ratio makes cpu depend on memory --
        # unsolvable above L ~= 1.6M. cpu_cores may therefore exceed threads: ratio headroom,
        # not parallelism.
        Int phase_threads = 4
        Int phase_kpbwt = 1000

        # Panel variants in this shard's input region -- GLIMPSE2's L. Required, and supplied
        # by CountPanelVariantsPerShard in the same run rather than read from a resource file,
        # so it cannot be stale or paired with the wrong shard.
        Int n_variants

        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"


        RuntimeAttr? runtime_attr_override
    }

    # imputation_hmm.cpp allocates Alpha as polymorphic_sites * modK floats, modK = n_states
    # rounded to a multiple of 8 and bounded by Kpbwt, so the matrix costs at most
    # 4 bytes * L * Kpbwt per thread. UPPER BOUND, not a fit: n_states is often well below
    # Kpbwt and logged L counts all panel sites.
    #
    # Against the observed 16 GiB pass/fail boundary at 4 threads / Kpbwt 1000:
    #   L =   657,784 -> 19 GiB  (passed at 16, marginally)
    #   L =   965,039 -> 24 GiB  (OOMed at 16)
    #   L = 1,346,888 -> 30 GiB  (OOMed at 16)
    # Typical L ~= 400k -> 15 GiB / 4 cpu.
    #
    # Sizing per shard is what makes fixing the OOMs free: phase costs $18.4/batch against
    # $18.7 for the old OOM-prone 4/16 shape, and $40.6 for a flat worst-case 8/40.
    # Parenthesised to force Float promotion first: as integers, phase_kpbwt * n_variants is
    # 2.7e9 at Kpbwt 2000, past Int32. Do not reorder.
    Int final_mem_gb  = 8 + ceil((((4.0 * phase_threads) * phase_kpbwt) * n_variants) / 1000000000.0)
    # N1 allows <= 6.5 GB/cpu and rounds any cpu count but 1 up to even; doing both here keeps
    # the requested shape visible instead of letting Cromwell adjust it silently.
    Int ratio_min_cpu = ceil(final_mem_gb / 6.5)
    Int unrounded_cpu = if ratio_min_cpu > phase_threads then ratio_min_cpu else phase_threads
    Int final_cpu     = unrounded_cpu + (unrounded_cpu % 2)

    # Sized from the actual inputs rather than a flat 50 GB.
    #
    # Peak is not just the inputs: Cromwell's checkpoint sync keeps two copies (cp to -tmp) and
    # reheader writes the output twice, so chr7 s12 peaks ~16 GB against 32 provisioned.
    # The 30 GB floor is deliberate: localization was measured only at 50 GiB (161 MiB/s, ~7x
    # the documented pd-ssd per-GiB scaling), which says nothing about the curve at 16-20 GiB.
    # Sam's TODO stands -- only one shard of input_vcf is used, so pre-splitting it upstream
    # would cut ~1.4 TB of redundant localization per batch.
    Int computed_disk_gb = 10 + ceil(2.0 * (size(panel_split_chunk_bin, "GB") + size(input_vcf, "GB")))
    Int disk_size_gb = if computed_disk_gb > 30 then computed_disk_gb else 30

    command <<<
        set -euxo pipefail

        # Peak RSS via EXIT trap on stderr. Must be a trap: under set -e a SIGKILLed child
        # aborts the script, so an end-of-script report would miss every OOM -- the only runs
        # worth measuring. Must be stderr: Cromwell delocalizes it for FAILED tasks, File
        # outputs do not exist. Guarded throughout; cannot fail the task it measures.
        _instr_gib() { for f in "$@"; do if [ -r "$f" ]; then awk '{printf "%.2f", $1/1073741824}' "$f" 2>/dev/null && return 0; fi; done; printf NA; }
        _instr_peak()  { _instr_gib /sys/fs/cgroup/memory.peak /sys/fs/cgroup/memory/memory.max_usage_in_bytes; }
        _instr_limit() { _instr_gib /sys/fs/cgroup/memory.max  /sys/fs/cgroup/memory/memory.limit_in_bytes; }
        _instr_cur()   { _instr_gib /sys/fs/cgroup/memory.current /sys/fs/cgroup/memory/memory.usage_in_bytes; }
        _INSTR_T0=$SECONDS
        _INSTR_SAMPLER=""
        RESUMED=unknown
        _instr_report() {
            _rc=$?
            set +e +x
            [ -n "$_INSTR_SAMPLER" ] && kill "$_INSTR_SAMPLER" 2>/dev/null
            _du=$(df -P -BG . 2>/dev/null | awk 'NR==2{gsub(/G/,"",$3); gsub(/G/,"",$2); print $3" "$2}')
            echo "[RESOURCE] task=GLIMPSE2Phase region=~{input_region} n_variants=~{n_variants} requested_mem_gib=~{eff_mem_gb} requested_cpu=~{eff_cpu} threads=~{phase_threads} kpbwt=~{phase_kpbwt} resumed=$RESUMED rc=$_rc peak_rss_gib=$(_instr_peak) limit_gib=$(_instr_limit) disk_used_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f1) disk_total_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f2) wall_s=$((SECONDS-_INSTR_T0))" >&2
        }
        trap _instr_report EXIT
        # 10 s time series, so a spike can be located against GLIMPSE2's Cnk/Buf markers.
        ( set +x; while :; do
            echo "[RESOURCE-TS] t=$((SECONDS-_INSTR_T0)) rss_gib=$(_instr_cur)" >&2
            sleep 10
          done ) & _INSTR_SAMPLER=$!



        # Fail fast: GLIMPSE2 swallows program_options errors and exits 0, so a bad value
        # would surface much later as a missing output file.
        if [ "~{phase_threads}" -lt 1 ] || [ "~{phase_kpbwt}" -lt 1 ]; then
            echo "ERROR: phase_threads and phase_kpbwt must both be >= 1 (got ~{phase_threads}, ~{phase_kpbwt})." >&2
            exit 1
        fi

        # Strip both from extra_phase_args and inject the typed values. Passing either twice
        # is a duplicate option, which boost program_options rejects -- and GLIMPSE2 catches
        # parser errors and exits 0, so it surfaces later as a confusing missing-BCF failure.
        # Honouring the string's value instead would decouple it from the memory sized on it.
        # The placeholder is substituted as literal text, so bash expands "$(nproc)" here
        # before stripping sees it. Both spellings are handled and both are tested.
        EXTRA_PHASE_ARGS="~{extra_phase_args}"
        for OPT in thread Kpbwt; do
            if echo "$EXTRA_PHASE_ARGS" | grep -qE "(^|[[:space:]])--${OPT}([[:space:]]|=|$)"; then
                echo "WARNING: --${OPT} found in extra_phase_args; overriding with the typed input." >&2
                echo "WARNING: set phase_threads / phase_kpbwt instead -- memory is sized from them." >&2
                # Two passes. The first removes the option together with its value, where the
                # value is any token not itself starting with '-' -- it must match more than
                # [0-9]+, because the historical default string used "--thread $(nproc)" and
                # leaving that behind would reintroduce the duplicate this block exists to
                # prevent. The second removes a valueless leftover (e.g. a trailing "--thread"),
                # which would otherwise survive as a bare duplicate.
                # Signed numeric first (else "--thread -1" leaves a stray "-1"), then any
                # value not starting with '-' (must beat [0-9]+, the old default was
                # "--thread $(nproc)"), then a valueless leftover.
                EXTRA_PHASE_ARGS=$(echo "$EXTRA_PHASE_ARGS" \
                    | sed -E "s/(^|[[:space:]])--${OPT}([[:space:]]+|=)-?[0-9]+/ /g" \
                    | sed -E "s/(^|[[:space:]])--${OPT}([[:space:]]+|=)[^-[:space:]][^[:space:]]*/ /g" \
                    | sed -E "s/(^|[[:space:]])--${OPT}([[:space:]]|=|$)/ /g")
            fi
        done

        cmd="/bin/GLIMPSE2_phase \
                --input-gl ~{input_vcf} \
                -R ~{panel_split_chunk_bin} \
                --thread ~{phase_threads} \
                --Kpbwt ~{phase_kpbwt} \
                $EXTRA_PHASE_ARGS \
                --output ~{output_prefix}.raw.bcf \
                --checkpoint-file-out checkpoint.bin"

        RESUMED=no
        if [ -s "checkpoint.bin" ]; then
            cmd="$cmd --checkpoint-file-in checkpoint.bin"
            RESUMED=yes
        fi

        eval "$cmd"

        # take input VCF header and add GLIMPSE INFO and FORMAT lines (GLIMPSE header only contains a single chromosome and breaks bcftools concat --naive)
        bcftools view --no-version -h ~{input_vcf} | grep '^##' > input.header.txt
        bcftools view --no-version -h ~{output_prefix}.raw.bcf | grep -E '^##INFO|^##FORMAT|^##NMAIN|^##FPLOIDY' > glimpse2.header.txt
        bcftools view --no-version -h ~{input_vcf} | grep '^#CHROM' > input.columns.txt
        cat input.header.txt glimpse2.header.txt input.columns.txt > header.txt
        bcftools reheader -h header.txt ~{output_prefix}.raw.bcf -o ~{output_prefix}.bcf
        bcftools index ~{output_prefix}.bcf

    >>>

    output {
        File phased_vcf = "~{output_prefix}.bcf"
        File phased_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          final_cpu,
        mem_gb:             final_mem_gb,
        disk_gb:            disk_size_gb,
        use_ssd:            true,
        preemptible_tries:  10,
        max_retries:        1,
        docker:             docker
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    # runtime_attr_override applies field by field, so overriding only mem_gb would keep a cpu
    # computed for the DEFAULT memory (64 GiB with 4 cpu is 16 GB/cpu, silently widened by
    # Cromwell). Re-derive from whichever memory won; an explicit cpu_cores still wins.
    Float eff_mem_gb  = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
    Int eff_ratio_cpu = ceil(eff_mem_gb / 6.5)
    Int eff_unrounded = if eff_ratio_cpu > phase_threads then eff_ratio_cpu else phase_threads
    Int eff_cpu       = select_first([runtime_attr.cpu_cores, eff_unrounded + (eff_unrounded % 2)])
    runtime {
        cpu:                    eff_cpu
        memory:                 eff_mem_gb + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        zones:                  zones
        checkpointFile:         "checkpoint.bin"
    }
}

task GLIMPSE2Ligate {
    input {
        Array[File] phased_vcfs
        Array[File] phased_vcf_idxs
        String output_prefix

        String docker

        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"


        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 2 * ceil(size(phased_vcfs, "GB")) + 10

    command <<<
        set -euox pipefail

        # Peak RSS via EXIT trap on stderr. Must be a trap: under set -e a SIGKILLed child
        # aborts the script, so an end-of-script report would miss every OOM -- the only runs
        # worth measuring. Must be stderr: Cromwell delocalizes it for FAILED tasks, File
        # outputs do not exist. Guarded throughout; cannot fail the task it measures.
        _instr_gib() { for f in "$@"; do if [ -r "$f" ]; then awk '{printf "%.2f", $1/1073741824}' "$f" 2>/dev/null && return 0; fi; done; printf NA; }
        _instr_peak()  { _instr_gib /sys/fs/cgroup/memory.peak /sys/fs/cgroup/memory/memory.max_usage_in_bytes; }
        _instr_limit() { _instr_gib /sys/fs/cgroup/memory.max  /sys/fs/cgroup/memory/memory.limit_in_bytes; }
        _instr_cur()   { _instr_gib /sys/fs/cgroup/memory.current /sys/fs/cgroup/memory/memory.usage_in_bytes; }
        _INSTR_T0=$SECONDS
        _INSTR_SAMPLER=""
        _instr_report() {
            _rc=$?
            set +e +x
            [ -n "$_INSTR_SAMPLER" ] && kill "$_INSTR_SAMPLER" 2>/dev/null
            _du=$(df -P -BG . 2>/dev/null | awk 'NR==2{gsub(/G/,"",$3); gsub(/G/,"",$2); print $3" "$2}')
            echo "[RESOURCE] task=GLIMPSE2Ligate n_shards=~{length(phased_vcfs)} requested_mem_gib=~{eff_mem_gb} requested_cpu=~{eff_cpu} threads=2 rc=$_rc peak_rss_gib=$(_instr_peak) limit_gib=$(_instr_limit) disk_used_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f1) disk_total_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f2) wall_s=$((SECONDS-_INSTR_T0))" >&2
        }
        trap _instr_report EXIT
        # 10 s time series, so a spike can be located against GLIMPSE2's Cnk/Buf markers.
        ( set +x; while :; do
            echo "[RESOURCE-TS] t=$((SECONDS-_INSTR_T0)) rss_gib=$(_instr_cur)" >&2
            sleep 10
          done ) & _INSTR_SAMPLER=$!


        # Threads pinned rather than $(nproc): bcf_sr_set_threads creates a thread pool, so
        # letting the thread count follow cpu_cores would grow memory alongside the widening
        # below. Two is a conservative choice pending measurement, not a tuned value.
        /bin/GLIMPSE2_ligate --input ~{write_lines(phased_vcfs)} --output ~{output_prefix}.bcf --thread 2

        # the index generated by ligate appears to be corrupt for both bcf and vcf.gz output (possibly due to https://github.com/samtools/htslib/issues/1740), so we regenerate with bcftools
        bcftools index -f ~{output_prefix}.bcf
    >>>

    output {
        File ligated_vcf = "~{output_prefix}.bcf"
        File ligated_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    # 32 GiB is an empirical cap, NOT a validated model.
    #
    # Established: rc=137 at 12 GiB on chr2 and chr20 twice each, chr7 on three attempts,
    # chr22 passing -- correlating with each chromosome's densest seam, predicted 3/3 blind.
    # TWO mechanisms have been proposed and BOTH are dead; recorded so neither is re-derived:
    #  1. Linear in seam length (mem ~= 1.5 + 18*L_isec/1e6). HTSlib 1.16 _reader_fill_buffer
    #     only buffers records sharing a coordinate, so nothing accumulates across a seam.
    #  2. Max records per coordinate. chr20 peaked at 10,465 and FAILED, chr22 at 9,937 and
    #     PASSED -- 5% apart, opposite outcomes -- and 10,465 x ~9 KB x 2 readers is ~190 MB,
    #     ~100x short of the ~12 GB to explain.
    # The cause is unidentified. Do not tune this on a model; every model has been wrong.
    # The EXIT-trap peak RSS above is what will settle it.
    # cpu 6 only keeps 32/6 under the N1 6.5 GB/cpu limit; with --thread 2, four cores idle.
    RuntimeAttr default_attr = object {
        cpu_cores:          6,
        mem_gb:             32,
        disk_gb:            disk_size_gb,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             docker
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    # runtime_attr_override applies field by field, so overriding only mem_gb would keep a cpu
    # computed for the DEFAULT memory (64 GiB with 4 cpu is 16 GB/cpu, silently widened by
    # Cromwell). Re-derive from whichever memory won; an explicit cpu_cores still wins.
    Float eff_mem_gb  = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
    Int eff_ratio_cpu = ceil(eff_mem_gb / 6.5)
    Int eff_unrounded = if eff_ratio_cpu > 2 then eff_ratio_cpu else 2
    Int eff_cpu       = select_first([runtime_attr.cpu_cores, eff_unrounded + (eff_unrounded % 2)])
    runtime {
        cpu:                    eff_cpu
        memory:                 eff_mem_gb + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        zones:                  zones
    }
}


task PopAndMarginalizeCollisions {
    input {
        # all VCFs should be split to biallelic
        File posteriors_vcf
        File posteriors_vcf_idx
        File panel_bubble_split_sites_only_vcf          # for annotation of INFO fields
        File panel_bubble_split_sites_only_vcf_idx
        File panel_id_split_vcf_gz           # panel popping script currently requires vcf.gz, so we also use that here
        File panel_id_split_vcf_gz_tbi
        
        File? pop_glimpse2_script             # modified version of convert-to-biallelic.py translated to Rust
        File? cargo_toml
        File? pop_glimpse2_binary
        
        String region
        String output_prefix

        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"


        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 3 * ceil(size([posteriors_vcf, panel_bubble_split_sites_only_vcf, panel_id_split_vcf_gz], "GB"))

    command <<<
        set -euox pipefail

        # Peak RSS via EXIT trap on stderr. Must be a trap: under set -e a SIGKILLed child
        # aborts the script, so an end-of-script report would miss every OOM -- the only runs
        # worth measuring. Must be stderr: Cromwell delocalizes it for FAILED tasks, File
        # outputs do not exist. Guarded throughout; cannot fail the task it measures.
        _instr_gib() { for f in "$@"; do if [ -r "$f" ]; then awk '{printf "%.2f", $1/1073741824}' "$f" 2>/dev/null && return 0; fi; done; printf NA; }
        _instr_peak()  { _instr_gib /sys/fs/cgroup/memory.peak /sys/fs/cgroup/memory/memory.max_usage_in_bytes; }
        _instr_limit() { _instr_gib /sys/fs/cgroup/memory.max  /sys/fs/cgroup/memory/memory.limit_in_bytes; }
        _instr_cur()   { _instr_gib /sys/fs/cgroup/memory.current /sys/fs/cgroup/memory/memory.usage_in_bytes; }
        _INSTR_T0=$SECONDS
        _INSTR_SAMPLER=""
        _instr_report() {
            _rc=$?
            set +e +x
            [ -n "$_INSTR_SAMPLER" ] && kill "$_INSTR_SAMPLER" 2>/dev/null
            _du=$(df -P -BG . 2>/dev/null | awk 'NR==2{gsub(/G/,"",$3); gsub(/G/,"",$2); print $3" "$2}')
            echo "[RESOURCE] task=PopAndMarginalizeCollisions region=~{region} requested_mem_gib=~{eff_mem_gb} requested_cpu=~{eff_cpu} rc=$_rc peak_rss_gib=$(_instr_peak) limit_gib=$(_instr_limit) disk_used_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f1) disk_total_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f2) wall_s=$((SECONDS-_INSTR_T0))" >&2
        }
        trap _instr_report EXIT
        # 10 s time series, so a spike can be located against GLIMPSE2's Cnk/Buf markers.
        ( set +x; while :; do
            echo "[RESOURCE-TS] t=$((SECONDS-_INSTR_T0)) rss_gib=$(_instr_cur)" >&2
            sleep 10
          done ) & _INSTR_SAMPLER=$!


        if [ -n "~{pop_glimpse2_binary}" ]; then
            POP_BIN="~{pop_glimpse2_binary}"
            chmod +x $POP_BIN
        else
            mkdir -p pop-glimpse2/src/bin
            cp ~{pop_glimpse2_script} pop-glimpse2/src/bin/pop-glimpse2.rs
            cp ~{cargo_toml} pop-glimpse2
            cd pop-glimpse2
            cargo build --release
            cd ..
            POP_BIN="./pop-glimpse2/target/release/pop-glimpse2"
        fi

        # this now only works for pop-glimpse2-joint-opt.rs;
        # the sort may also be extraneous, but we keep it in to guard against getting out of sync with the popped panel
        bcftools view -r ~{region} --regions-overlap 0 ~{panel_bubble_split_sites_only_vcf} -Oz -o panel.bubble.split.sites.shard.vcf.gz
        bcftools view -r ~{region} --regions-overlap 0 ~{posteriors_vcf} | \
            $POP_BIN ~{panel_id_split_vcf_gz} panel.bubble.split.sites.shard.vcf.gz | \
            bcftools sort --max-mem=2G -W -Ob -o ~{output_prefix}.bcf

    >>>

    output {
        File popped_vcf = "~{output_prefix}.bcf"
        File popped_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    # DELIBERATELY UNCHANGED at 12 GiB / 2 cpu. chr2 and chr20 have never reached this stage
    # (ligate failed first), so it is untested at the sizes that matter -- but no pop task has
    # ever OOMed, and pop_regions defaults to output_regions, so this scatters ~523 times per
    # batch, not 22. Doubling it would cost ~$3.5/batch for no evidence. "Downstream of a
    # failure" is not evidence; it is true of every task after phase.
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             12,
        disk_gb:            disk_gb,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsde-methods/slee/lrma-aou2-panel-creation-rust:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    # runtime_attr_override applies field by field, so overriding only mem_gb would keep a cpu
    # computed for the DEFAULT memory (64 GiB with 4 cpu is 16 GB/cpu, silently widened by
    # Cromwell). Re-derive from whichever memory won; an explicit cpu_cores still wins.
    Float eff_mem_gb  = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
    Int eff_ratio_cpu = ceil(eff_mem_gb / 6.5)
    Int eff_unrounded = if eff_ratio_cpu > 2 then eff_ratio_cpu else 2
    Int eff_cpu       = select_first([runtime_attr.cpu_cores, eff_unrounded + (eff_unrounded % 2)])
    runtime {
        cpu:                    eff_cpu
        memory:                 eff_mem_gb + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        zones:                  zones
    }
}

task RemapSampleNames {
    input {
        File vcf
        File vcf_idx
        File remap_file
        String output_prefix

        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"


        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 2 * ceil(size(vcf, "GB"))

    command <<<
        set -euxo pipefail

        bcftools reheader --samples ~{remap_file} ~{vcf} -o ~{output_prefix}.bcf
        bcftools index ~{output_prefix}.bcf
    >>>

    output {
        File output_vcf = "~{output_prefix}.bcf"
        File output_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size_gb,
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
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        zones:                  zones
    }
}
