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
        Int ligate_threads = 4
        # NOTE: GLIMPSE2 is not deterministic under multithreading. chr22 run twice with
        # identical parameters gave identical AF at only 11.2% of sites (AF r = 0.999982,
        # max |dAF| = 0.0247, INFO r = 0.990) -- agreement in distribution, not bit-for-bit.
        # Pinning threads makes the thread COUNT reproducible, not the output.

        String output_prefix

        # inputs for PopAndMarginalizeCollisions
        File pop_glimpse2_panel_resources_json
        File? pop_glimpse2_script               # heavily modified version of convert-to-biallelic.py
        File? pop_glimpse2_cargo_toml
        File? pop_glimpse2_binary

        # Pinned by digest, not tag: tags in GCR are mutable, so a repush would silently
        # change what runs. This digest is tag 1.0.0-2cee597-1778869818 as of 2026-07-31.
        # Enables checkpointing; note it contains bcftools/htslib 1.16.
        #
        # NB the staged input CSVs override this with :tfenne-opt, which IS a mutable tag --
        # pinning it there matters more than pinning the default here, and that is a change to
        # the run config rather than to this file.
        String glimpse2_docker = "us.gcr.io/broad-gotc-prod/imputation-glimpse2@sha256:c3d64c5af3b8e789bcda94451931f09073d3a8750fc7375a8d44246d193e8a21"

        # Byte-identical to Cromwell's own ZonesDefaultValue: this RESTORES the stock four
        # zones, which the VWB backend overrode to pin us-central1-a (71% of 502 attempts were
        # preempted there). ASK VWB FIRST -- pinning may be a data-locality or compliance
        # decision, not an oversight. Set to "us-central1-a" to restore the pin.
        # ZonesValidation accepts String and Array[String] alike, so the form is style only.
        # Values must lie in the backend's Batch job region. Does not reach the imported
        # ConcatVcfs call, which takes no zones argument.
        # MEASURED 2026-07-31, and the answer is no: with these four zones set, a submitted
        # job still reported
        #     allowedLocations: ['regions/us-central1', 'zones/us-central1-a']
        # so on the VWB backend a task-level zones does NOT override the configured
        # allowedLocations. This attribute is therefore inert here. Kept because it is correct
        # on backends that do honour it and it documents the intent, but widening the spot pool
        # on VWB requires a backend config change, not a WDL change -- raise it with them
        # alongside checkpointingInterval.
        #
        # This matters more than the cost: measured preemption ran 41-50% (run1, 4cpu/16G) and
        # 63-72% (run2, 8cpu/40G) in the first hours of each run, decaying to ~17-21% later.
        # run2 had only 150 shards yet preempted harder than run1's 523, so machine shape and
        # clock both matter -- but 500+ shards of one shape against a single zone plausibly
        # congests that zone itself. Until the backend allows more zones, staggering
        # submissions is the only lever available from this side.
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
            ligate_threads = ligate_threads,
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
    # CONFIRMED 2026-07-31: a submitted job reported bootDiskMib 28611 (30 GB), against 38148
    # (40 GB) before -- 10 GB per shard, ~5.2 TB per batch. An explicit 0 would keep this
    # override live, but it is unverified that BootDiskSizeValidation accepts 0, and that
    # would fail every task at once.
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


    #########################
    # NOT preemptible: every phase shard blocks on this, so it is a serialisation point and a
    # per-chromosome single point of failure. MEASURED on chr20 (7 regions): 65 s wall,
    # 2.49 GiB peak against 4 requested, 1 GB of disk. Cheaper and faster than the "minutes"
    # estimated -- the larger chromosomes have ~47 regions so will take longer, but the
    # preemption exposure is smaller than assumed. Kept off spot anyway: it costs cents, and
    # one preemption here stalls every downstream shard.
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
    # runtime_attr_override applies field by field, so overriding only mem_gb would keep a cpu
    # computed for the DEFAULT memory. Re-derive from whichever memory won; an explicit
    # cpu_cores still wins. Same pattern as the other tasks in this file.
    Float eff_mem_gb  = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
    Int eff_ratio_cpu = ceil(eff_mem_gb / 6.5)
    Int eff_unrounded = if eff_ratio_cpu > 2 then eff_ratio_cpu else 2
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

        # --regions-overlap 0 matches GLIMPSE2: a record counts if its POS is inside. The
        # default mode would also include records starting before the region and inflate.
        # Indexed read per region, not one bucketed stream: regions overlap so this reads
        # ~1.5x the file, and a miscount mis-sizes a shard silently. Correctness over speed.
        REGIONS_FILE=~{write_lines(input_regions)}
        : > counts.txt
        while read -r REGION; do
            # --threads is ADDITIONAL worker threads, so 1 here is 2 total on 2 cpu.
            bcftools view --no-version --threads 1 -H -r "$REGION" --regions-overlap 0 \
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

    # FITTED FROM MEASURED PEAK RSS, superseding the earlier source-derived bound.
    #
    # That bound was 4 bytes per site per thread per state, from imputation_hmm.cpp allocating
    # Alpha as polymorphic_sites * modK floats. It was never actually a bound: it counted one
    # matrix, and the process allocates more than Alpha. Eight instrumented chr20 shards give
    #     peak_GiB = 0.37 + 2.17e-5 * L        (threads 4, Kpbwt 1000)
    # i.e. 5.83 bytes per site per thread per state -- the old slope was 46% short. It had not
    # bitten only because the +8 GiB intercept was generous while the true fixed cost is
    # ~0.37 GiB, and those two errors cancel at small L and stop cancelling as L grows:
    # measured utilisation ran 55% at L=331k to 90% at L=1.005M, monotonically. Extrapolated,
    # chr11 s9 would have sat at 95% of request and chr7 s12 at 99%, against a hard cgroup
    # limit.
    #
    # CAVEAT on the fit: cgroup memory.peak counts page cache as well as anonymous memory,
    # and these shards read 9-17 GB of input. The kernel reclaims cache before OOM-killing, so
    # part of the measured peak is reclaimable and the true allocation slope is below 5.82.
    # Sizing on what the cgroup accounts is the safe direction -- it is also what the OOM
    # killer sees -- but it means this is an upper bound on demand, not a measurement of it.
    # Anyone tightening these numbers should separate anon from cache first (memory.stat)
    # rather than trusting the slope.
    #
    # 8.0 with a 2 GiB constant gives ~64% utilisation flat across the range, never
    # under-provisions any measured shard, and costs about $18.6 per batch against $18.4 for
    # the under-sloped version -- so the correction is nearly free.
    #   L =   400,379 -> 15 GiB / 4 cpu   (measured peak  9.32, 62%)
    #   L = 1,005,104 -> 35 GiB / 6 cpu   (measured peak 22.44, 64%)
    #   L = 1,346,888 -> 46 GiB / 8 cpu   (projected     29.61, 64%)
    #
    # Sizing per shard is what makes fixing the OOMs free: phase costs $18.4/batch against
    # $18.7 for the old OOM-prone 4/16 shape, and $40.6 for a flat worst-case 8/40.
    # Parenthesised to force Float promotion first: as integers, phase_kpbwt * n_variants is
    # 2.7e9 at Kpbwt 2000, past Int32. Do not reorder.
    Int final_mem_gb  = 2 + ceil((((8.0 * phase_threads) * phase_kpbwt) * n_variants) / 1000000000.0)
    # N1 allows <= 6.5 GB/cpu and rounds any cpu count but 1 up to even; doing both here keeps
    # the requested shape visible instead of letting Cromwell adjust it silently.
    Int ratio_min_cpu = ceil(final_mem_gb / 6.5)
    Int unrounded_cpu = if ratio_min_cpu > phase_threads then ratio_min_cpu else phase_threads
    Int final_cpu     = unrounded_cpu + (unrounded_cpu % 2)

    # Sized from the actual inputs rather than a flat 50 GB.
    #
    # Peak is not just the inputs: Cromwell's checkpoint sync keeps two copies (cp to -tmp) and
    # reheader writes the output twice, so chr7 s12 peaks ~16 GB against 32 provisioned.
    # The 30 GB floor is now exercised rather than assumed. Every chr20 shard ran on it (df
    # reports 29 GB usable) and the worst used 17 GB -- 59%, on the L=1,005,104 shard -- with
    # wall times in the normal 842-3152 s range, so the floor is neither wasteful nor a
    # throughput cliff at this size. What remains untested is smaller: localization was only
    # ever timed at 50 GiB (161 MiB/s, ~7x the documented pd-ssd per-GiB scaling), so the curve
    # below 30 GB is still unknown and the floor is what keeps us out of it.
    # Sam's TODO stands -- only one shard of input_vcf is used, so pre-splitting it upstream
    # would cut ~1.4 TB of redundant localization per batch.
    Int computed_disk_gb = 10 + ceil(2.0 * (size(panel_split_chunk_bin, "GB") + size(input_vcf, "GB")))
    Int disk_size_gb = if computed_disk_gb > 30 then computed_disk_gb else 30


    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          final_cpu,
        mem_gb:             final_mem_gb,
        disk_gb:            disk_size_gb,
        use_ssd:            true,
        # 10 spot attempts, and the measurement says spot is worth it. Over the instrumented
        # chr20 run: 415 successful VM-minutes against 286 wasted across 32 preemptions, a 41%
        # waste fraction, i.e. 1.69x the VM-time of a clean run.
        #
        # That comes to ~0.58x the cost of running the same work non-preemptible. NOT 0.51x:
        # spot discounts vCPU and memory to 30% of list, but persistent disk bills the same
        # either way, and the wasted attempts pay for their disk too. On a 4 cpu / 15 GiB shard
        # with 30 GB of working disk and a 30 GB boot disk, compute is $0.1996/hr on demand
        # against $0.0599 spot, while disk is $0.0140/hr regardless:
        #     1.69 * (0.0599 + 0.0140) / (0.1996 + 0.0140) = 0.58
        # A compute-only calculation gives 0.507 and understates the true cost by about a
        # sixth. Avoiding checkpoint writes on demand recovers a further 0.8% (15 s of writes
        # on a ~1873 s shard), which is inside the rounding. Still a clear win for spot.
        #
        # Checkpointing is doing real work
        # here: 5 of 11 shards resumed, and resumed shards averaged 954 s against 1873 s fresh,
        # for 15 s of total checkpoint-write overhead. An earlier note in this branch claimed
        # resumes never fired; that was wrong.
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
        RESUMED=unknown
        _instr_report() {
            _rc=$?
            set +e +x
            [ -n "$_INSTR_SAMPLER" ] && kill "$_INSTR_SAMPLER" 2>/dev/null
            # -Pk is POSIX; -BG is GNU-only and absent on busybox. 1K blocks -> GiB in awk.
            _du=$(df -Pk . 2>/dev/null | awk 'NR==2{printf "%.0f %.0f", $3/1048576, $2/1048576}')
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

        # Typed, and used as the cpu floor below, so threads can never exceed the cores the
        # task is given. Hardcoding the floor at 2 while the command asked for 4 threads meant
        # an override of mem_gb to 12 produced eff_cpu 2 running --thread 4 -- oversubscribed.
        # phase avoids this by flooring on phase_threads; ligate now does the same.
        #
        # Left at 2, not raised to 4. The 12.39 GiB peak this task is sized against was
        # measured at 2 threads, and bcf_sr_set_threads allocates per-thread decompression
        # buffers, so 4 threads is a different configuration than the one that was measured.
        # An earlier revision did raise it, reasoning that ~1 GiB per extra thread would fit
        # in the headroom -- but that is a derived estimate, and every derived estimate in this
        # file has been wrong (phase by 46% on slope, ligate by 4.5x on coefficient). Raise it
        # once a 4-thread peak exists; the instrumentation will report one.
        Int ligate_threads = 2

        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"


        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 2 * ceil(size(phased_vcfs, "GB")) + 10


    #########################
    # 24 GiB / 4 cpu, from a measured peak.
    #
    # MEASURED on chr20: peak 12.39 GiB, rc=0, 618 s. The 10 s series puts the spike exactly
    # on Buf 4 [L_isec=744,316] -- the seam that took two rc=137 kills at 12 GiB. It needed
    # 12.39 and had 12.00. The margin was 0.39 GiB.
    #
    # Two anchors from that series (L_isec 117,072 -> 10.01 GiB, 744,316 -> 12.39) give
    #     peak_GiB ~= 9.57 + 3.79e-6 * L_isec        (~4.0 KB per intersecting site)
    # An earlier version of this comment declared BOTH proposed mechanisms dead. That was half
    # wrong: seam size is real, but at ~4 KB/site on a ~9.6 GiB baseline, not the 18 KB/site on
    # a 1.5 GiB baseline first claimed -- right variable, wrong intercept and coefficient, the
    # same failure as the phase formula. What is genuinely dead is max-records-per-coordinate:
    # chr20 peaked at 10,465 records and FAILED while chr22 peaked at 9,937 and PASSED.
    #
    # The spike is TRANSIENT -- roughly 20 s of the 618 s run, back to ~10.5 GiB immediately
    # after. That is why reading ligater_algorithm.cpp for a structure that grows across the
    # seam kept coming up empty: nothing accumulates, and the peak is a short-lived allocation
    # while the dense seam is processed. The static analysis was correct about the code and
    # the wrong conclusion was drawn from it. Note also that the 10 s sampler and the cgroup
    # exit trap agree to the digit (12.39), so the poller did not miss a taller, shorter spike.
    #
    # The model also explains why 12 GiB behaved like a coin flip. Projected peaks across the
    # genome's worst seams span 10.5-13.5 GiB, so chr2 (11.83), chr10 (11.91) and chr20 (12.39)
    # all sit within half a gigabyte of the old limit -- which is why chr10 passed, chr2 failed,
    # and chr7 needed three attempts.
    #
    # 24 GiB is 1.9x the ONE measured peak. That measurement is the basis for the number; the
    # model below is a secondary check and should not be read as genome-wide coverage.
    #
    # Applying it to every chromosome's worst seam (computed from chunks.tsv against the sites
    # BCF -- a method that reproduces chr2's 596,913, chr20's 744,316 and chr7's 1,023,060
    # exactly) puts the projected maximum at chr7, 13.45 GiB, 56% of this request. But that is
    # a two-point fit from a single chromosome extrapolated to twenty-one others, and one
    # circulating figure puts chr2 nearer 16 GiB. This model would under-predict that by ~39%.
    # 24 GiB still covers 16 at 68%, which is why the request is sized on the measurement plus
    # margin rather than on the fit. Treat any projection here as indicative until a second
    # chromosome reports a peak; chr2 and chr11 are the next to run and will settle it.
    #
    # chr1 is a separate gap: the position file used for that sweep held 370,261 entries
    # against the 10,112,850 records chr1 actually imputed, so it was truncated and chr1's
    # worst seam is unknown rather than zero.
    #
    # Down from 32, which measured 39% utilised. The reason to shrink is not
    # the $0.28 per batch: 24/4 is a materially smaller VM than 32/6, and wider shapes were
    # measured to preempt harder (63-72% on 8cpu/40G against 41-50% on 4cpu/16G), so the
    # narrower shape should place more easily on spot.
    #
    # Caveat as with phase: memory.peak includes page cache, so part of 12.39 is reclaimable.
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             24,
        disk_gb:            disk_size_gb,
        use_ssd:            true,
        # 2, unlike pop's 4, and for a reason rather than by inheritance. The spot arithmetic
        # that justified raising pop barely applies here: ligate is 22 tasks costing ~$0.26 per
        # batch in total, so the spread between 2 and 5 attempts is about $0.03. What does
        # matter is that ligate has no checkpoint and a failure costs the whole chromosome --
        # chr7 needed three attempts, which under tries=2 means two preemptions then a
        # guaranteed on-demand success. Falling back quickly is the reliable path, and at this
        # cost the reliability is effectively free.
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
    Int eff_unrounded = if eff_ratio_cpu > ligate_threads then eff_ratio_cpu else ligate_threads
    Int eff_cpu       = select_first([runtime_attr.cpu_cores, eff_unrounded + (eff_unrounded % 2)])
    command <<<
        set -euox pipefail

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
            echo "[RESOURCE] task=GLIMPSE2Ligate n_shards=~{length(phased_vcfs)} requested_mem_gib=~{eff_mem_gb} requested_cpu=~{eff_cpu} threads=~{ligate_threads} rc=$_rc peak_rss_gib=$(_instr_peak) limit_gib=$(_instr_limit) disk_used_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f1) disk_total_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f2) wall_s=$((SECONDS-_INSTR_T0))" >&2
        }
        trap _instr_report EXIT
        # 10 s time series, so a spike can be located against GLIMPSE2's Cnk/Buf markers.
        ( set +x; while :; do
            echo "[RESOURCE-TS] t=$((SECONDS-_INSTR_T0)) rss_gib=$(_instr_cur)" >&2
            sleep 10
          done ) & _INSTR_SAMPLER=$!


        # Pinned rather than $(nproc), because bcf_sr_set_threads creates a thread pool and
        # memory would otherwise track cpu_cores. See ligate_threads in the inputs for why it
        # is 2 and what would justify raising it.
        #
        # NOTE: eff_cpu is 4 here while only 2 threads run, and that is deliberate. The cpu
        # count is set by the N1 memory ratio -- 24 GiB needs ceil(24/6.5) = 4 cpu -- not by
        # any parallelism target. The two idle cores are the price of the memory, not an
        # oversight, and matching threads to them would change the configuration the 12.39 GiB
        # measurement was taken under.
        /bin/GLIMPSE2_ligate --input ~{write_lines(phased_vcfs)} --output ~{output_prefix}.bcf --thread ~{ligate_threads}

        # the index generated by ligate appears to be corrupt for both bcf and vcf.gz output (possibly due to https://github.com/samtools/htslib/issues/1740), so we regenerate with bcftools
        bcftools index -f ~{output_prefix}.bcf
    >>>

    output {
        File ligated_vcf = "~{output_prefix}.bcf"
        File ligated_vcf_idx = "~{output_prefix}.bcf.csi"
    }
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

    # Measured on chr20: 2 GB used of the 16 this produced. Trimmed 3x -> 2x, which still
    # leaves ~5x the observed high-water. Worth doing only because pop scatters ~523 times per
    # batch, so 2 GB per task is ~1 TB of the SSD quota that bounds how many batches run at once.
    Int disk_gb = 10 + 2 * ceil(size([posteriors_vcf, panel_bubble_split_sites_only_vcf, panel_id_split_vcf_gz], "GB"))


    #########################
    # 8 GiB / 2 cpu, from measurement. This task was previously left at 12 GiB on the grounds
    # that nothing was known about it either way; eleven instrumented chr20 shards now put the
    # worst peak at 4.68 GiB, 39% of that request.
    #
    # An earlier revision proposed RAISING this to 24 GiB because chr2 and chr20 had never
    # reached the stage. That was refused for lack of evidence, and the refusal was right for
    # the wrong reason -- the evidence, once collected, points the other way.
    #
    # 8 GiB is 1.7x the worst measured peak, the same margin used for ligate and phase. It
    # matters more here than anywhere else in the file because pop_regions defaults to
    # output_regions, so this scatters ~523 times per batch -- the same order as phase, not the
    # 22 that ligate runs. 12 -> 8 is ~$0.5 per batch, ~$200 across 200 batches.
    #
    # CAVEAT: all eleven measurements are chr20. The regions are chunked to a similar variant
    # count on every chromosome so the peak should not scale with chromosome size, but bubble
    # density does vary. The [RESOURCE] lines will say so if that is wrong.
    # bcftools sort --max-mem=2G in the command sets the floor.
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_gb,
        use_ssd:            true,
        # 4 rather than 2. Expected cost per success, spot at 30% of on-demand: at the
        # measured 40-60% preemption rate, 2 tries costs 0.58-0.84 of on-demand and 4-5 tries
        # costs 0.51-0.77, because falling through to on-demand is the expensive outcome. The
        # price of the extra attempts is wall-clock, and pop restarts cheaply -- 226-475 s
        # measured, no checkpoint to lose. Above ~70% preemption the ordering reverses and
        # fewer tries would win, so this is not worth pushing higher.
        preemptible_tries:  4,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsde-methods/slee/lrma-aou2-panel-creation-rust@sha256:0f25c4091c49d8eb0c3d8bcdb45e7680a093e0a691e74ccec2d3743758a7d22c"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    # runtime_attr_override applies field by field, so overriding only mem_gb would keep a cpu
    # computed for the DEFAULT memory (64 GiB with 4 cpu is 16 GB/cpu, silently widened by
    # Cromwell). Re-derive from whichever memory won; an explicit cpu_cores still wins.
    Float eff_mem_gb  = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
    Int eff_ratio_cpu = ceil(eff_mem_gb / 6.5)
    Int eff_unrounded = if eff_ratio_cpu > 2 then eff_ratio_cpu else 2
    Int eff_cpu       = select_first([runtime_attr.cpu_cores, eff_unrounded + (eff_unrounded % 2)])
    command <<<
        set -euox pipefail

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
        # --threads is ADDITIONAL worker threads, so 1 is 2 total on this task's 2 cpu.
        bcftools view --threads 1 -r ~{region} --regions-overlap 0 ~{panel_bubble_split_sites_only_vcf} -Oz -o panel.bubble.split.sites.shard.vcf.gz
        bcftools view --threads 1 -r ~{region} --regions-overlap 0 ~{posteriors_vcf} | \
            $POP_BIN ~{panel_id_split_vcf_gz} panel.bubble.split.sites.shard.vcf.gz | \
            bcftools sort --max-mem=2G -W -Ob -o ~{output_prefix}.bcf

    >>>

    output {
        File popped_vcf = "~{output_prefix}.bcf"
        File popped_vcf_idx = "~{output_prefix}.bcf.csi"
    }
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


    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size_gb,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools@sha256:f820f1708624242c9b35912be83446d502d36ddadfe0f3ccac6492591314454c"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    # Same partial-override fix as the other tasks: select_first is per field, so overriding
    # only mem_gb would keep a cpu computed for the default memory and breach the N1 ratio.
    Float eff_mem_gb  = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
    Int eff_ratio_cpu = ceil(eff_mem_gb / 6.5)
    Int eff_unrounded = if eff_ratio_cpu > 2 then eff_ratio_cpu else 2
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
            echo "[RESOURCE] task=RemapSampleNames requested_mem_gib=~{eff_mem_gb} requested_cpu=~{eff_cpu} rc=$_rc peak_rss_gib=$(_instr_peak) limit_gib=$(_instr_limit) disk_used_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f1) disk_total_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f2) wall_s=$((SECONDS-_INSTR_T0))" >&2
        }
        trap _instr_report EXIT
        # 10 s time series, so a spike can be located against GLIMPSE2's Cnk/Buf markers.
        ( set +x; while :; do
            echo "[RESOURCE-TS] t=$((SECONDS-_INSTR_T0)) rss_gib=$(_instr_cur)" >&2
            sleep 10
          done ) & _INSTR_SAMPLER=$!

        bcftools reheader --samples ~{remap_file} ~{vcf} -o ~{output_prefix}.bcf
        bcftools index ~{output_prefix}.bcf
    >>>

    output {
        File output_vcf = "~{output_prefix}.bcf"
        File output_vcf_idx = "~{output_prefix}.bcf.csi"
    }
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
