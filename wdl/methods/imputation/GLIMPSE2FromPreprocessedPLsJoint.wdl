version 1.0

# Resource sizing and cost rationale: docs/glimpse2-cost-and-resources.md
#
# Short version, so it is not lost: base compute is ~11 cents/sample, which matches the Terra
# figure and is effectively the floor. The observed 1.97 full-run-equivalents per success --
# preempted shards restarting from zero because the checkpoint sync interval (600 s) is longer
# than the median preemption (7 min) -- is worth more than every knob in this file combined,
# and it is backend config, not a WDL setting. Also noted there: tasks default to the N1
# machine family, and the levers deliberately NOT taken (thread count, batch size, Kpbwt).

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

        # First-class workflow inputs, not buried in extra_phase_args and not left as
        # call-qualified task inputs: GLIMPSE2Phase computes its memory request from both, so
        # setting them must not require reaching into
        # GLIMPSE2FromPreprocessedPLsJoint.ChunkedGLIMPSE2Phase.* -- the same fragile mechanism
        # this file avoids elsewhere.
        Int phase_threads = 4
        Int phase_kpbwt = 1000

        String output_prefix

        # inputs for PopAndMarginalizeCollisions
        File pop_glimpse2_panel_resources_json
        File? pop_glimpse2_script               # heavily modified version of convert-to-biallelic.py
        File? pop_glimpse2_cargo_toml
        File? pop_glimpse2_binary

        String glimpse2_docker = "us.gcr.io/broad-gotc-prod/imputation-glimpse2:1.0.0-2cee597-1778869818"    # enables checkpointing, but note this contains bcftools/htslib 1.16!

        # Widens the spot pool: the VWB backend pins allocation to us-central1-a, where 355 of
        # 502 observed attempts (71%) were preempted on the wider phase shape. Four zones cost
        # nothing here -- shards never communicate and in-region GCS traffic is free.
        #
        # HARD CONSTRAINT: these must lie inside the backend's configured Batch job region.
        # Google is enforcing that allowedLocations match the job region (from 2026-07-31 for
        # affected projects), so a backend configured for anything but us-central1 will reject
        # these values. Exposed as an input rather than hardcoded so such a backend can
        # override without editing the WDL.
        #
        # NOTE: this does not reach the imported ConcatVcfs.ConcatVcfs call, which takes no
        # zones argument, so that one task still runs wherever the backend places it.
        #
        # Space-separated String rather than Array[String]. To be clear about why, since an
        # earlier version of this comment gave the wrong reason: ZonesValidation accepts BOTH
        # forms, in the same PartialFunction --
        #     coercion = Set(WomStringType, WomArrayType(WomStringType))
        #     case WomString(s) => s.split("\\s+").toVector.validNel
        #     case WomArray(womType, value) if womType.memberType == WomStringType => ...
        # so neither is safer than the other and this is purely a style choice. Do not "fix"
        # it back to an array believing the string form is legacy, and do not assume the
        # reverse either.
        #
        # IMPORTANT -- this value is byte-identical to Cromwell's own ZonesDefaultValue:
        #     private val ZonesDefaultValue =
        #         WomString("us-central1-a us-central1-b us-central1-c us-central1-f")
        # so this is not introducing four-zone scheduling. It is RESTORING the stock default
        # that the VWB backend overrode via configDefaultWomValue, since the observed policy
        # was allowedLocations: ['regions/us-central1', 'zones/us-central1-a'].
        #
        # That someone deliberately pinned this workspace to a single zone is a reason to ask
        # before overriding it, not merely a performance detail: it could be a data-locality
        # or compliance decision. A task-level attribute does outrank a config default, so
        # this should take effect -- which is exactly why it should be raised with VWB first.
        # Set zones to "us-central1-a" to restore the pinned behaviour without editing this.
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
    

    # Per-shard variant counts (GLIMPSE2's L), so phase can be sized per shard rather than
    # every shard being sized for the worst one. Counting here rather than reading a field
    # from chunked_panel_json is what makes this safe: the counts are derived in this run,
    # from this chromosome's sites-only panel, over these input_regions in order, so they
    # cannot be stale, reordered, or copied from another chromosome. A JSON field could be
    # all three while still having the right length.
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
    # NOT consumed by the tasks in this file, deliberately. Cromwell adds its own 30 GB
    # default to whatever bootDiskSizeGb is requested, so the previous value of 10 provisioned
    # 40 GB; omitting the attribute entirely yields the 30 GB minimum and saves 10 GB across
    # every shard. Passing an explicit 0 would express the same intent and keep this override
    # live, but it is not established that BootDiskSizeValidation accepts 0 rather than
    # requiring a positive int -- and that failure mode would hit every task at once. A
    # documented no-op is the cheaper mistake. Retained in the struct so callers sharing this
    # RuntimeAttr shape do not break; re-wire it once 0 is known to validate.
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

    # NOTE: deliberately no per-shard variant count here.
    #
    # GLIMPSE2Phase is sized from each shard's L, but that count is produced at runtime by
    # CountPanelVariantsPerShard rather than carried in this struct, and that is a safety
    # decision rather than an oversight.
    #
    # These are unkeyed parallel arrays. A count array added here would have no binding to the
    # shard it describes, so one that was stale but happened to have the right length would
    # pair silently with the wrong bins -- and underprovisioned memory is the exact failure the
    # sizing exists to prevent. A length check catches a truncated array but not a reordered
    # one, a panel regenerated with the same shard count, or counts copied from another
    # chromosome. It would also be dead on arrival: GLIMPSE2ChunkAndSplitPanel builds this
    # struct from four fields and never emits a count, so the field could only ever be
    # populated by hand-editing a generated resource.
    #
    # Counting at runtime removes all of that. The counts come from this chromosome's
    # sites-only panel, over these input_regions, in order, in the same run -- correct by
    # construction, and costing one small task per chromosome.
}

struct PopAndMarginalizePanelResourcesChromosome {
    String panel_bubble_split_sites_only_vcf
    String panel_bubble_split_sites_only_vcf_idx
    String panel_id_split_vcf_gz
    String panel_id_split_vcf_gz_tbi
    Array[String]? pop_regions              # non-overlapping, if not provided then GLIMPSE2 chunks will be used
}

# Counts panel variants in each shard's *input* (buffered) region -- exactly GLIMPSE2's L.
#
# Verified against GLIMPSE2's own logged L on three shards spanning the range, matching
# exactly: chr22 shard 0 = 657,784; chr20 shard 4 = 965,039; chr7 shard 12 = 1,346,888.
#
# Note this is NOT column 7 of chunks.tsv, which counts a different region and disagrees badly
# (403,101 against an actual L of 965,039 for chr20 shard 4).
#
# One task per chromosome, reading a sites-only BCF -- a few minutes and a couple of cents,
# against roughly $20 per batch saved by not sizing every shard for the worst one.
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

        # --regions-overlap 0 matches how GLIMPSE2 selects sites for the input region: a
        # record counts if its POS falls inside, regardless of how far its REF span reaches.
        # Using the default overlap mode would include records starting before the region and
        # inflate every count.
        # One indexed read per region rather than a single stream bucketed in awk. The regions
        # are buffered and overlap, so the loop reads roughly 1.5x the file overall, not once
        # per region -- and the obviously-correct form is worth more here than the faster one,
        # because a miscount silently mis-sizes a shard rather than failing. If this ever
        # becomes a bottleneck on the larger chromosomes, a single `bcftools query -f '%POS\n'`
        # piped into an interval-bucketing awk is the drop-in replacement.
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

        # read_json wants a JSON array; read_lines would give Array[String], which does not
        # coerce to Array[Int] reliably across engines.
        printf '[%s]\n' "$(paste -sd, counts.txt)" > ~{output_prefix}.n_variants.json
        cat ~{output_prefix}.n_variants.json
    >>>

    output {
        Array[Int] n_variants = read_json("~{output_prefix}.n_variants.json")
        File n_variants_json = "~{output_prefix}.n_variants.json"
    }

    #########################
    # NOT preemptible, deliberately. Every phase shard for this chromosome blocks on this
    # task, so it is both a serialisation point and a single point of failure for the whole
    # chromosome. Streaming a ~1.2 GB sites BCF across ~47 regions on chr2 takes minutes,
    # against a median preemption at 7 minutes -- so on the big chromosomes a preemptible
    # instance would frequently restart from scratch and delay all 523 downstream shards.
    # This task costs cents; taking it off spot removes a new failure mode from the critical
    # path for less than the price of one preemption.
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

        # Both of these are typed inputs rather than text inside extra_phase_args because the
        # memory request below is computed from them. Anything that changes the size of the
        # forward matrix has to be visible to the sizing arithmetic; a caller who replaced
        # extra_phase_args to add one unrelated option would otherwise silently drop
        # --Kpbwt 1000 and get this build's default of 2000, doubling the matrix while the
        # memory request stayed put. The command strips both options from extra_phase_args
        # and injects these instead.
        #
        # phase_threads is pinned rather than $(nproc): the forward matrix is allocated per
        # thread, so under $(nproc) memory would be a function of cpu_cores while the N1
        # ratio limit makes cpu_cores a function of memory. Solving that coupling for a shard
        # of L variants gives cpu * (6.5 - 4*L/1e6) >= 5.5, which is unsolvable above
        # L ~= 1.6M. Pinning decouples them. It is also why cpu_cores below may exceed
        # phase_threads -- the extra cores buy ratio headroom, not parallelism.
        Int phase_threads = 4
        Int phase_kpbwt = 1000

        # Panel variants in this shard's input region -- GLIMPSE2's L. Required, and supplied
        # by CountPanelVariantsPerShard in the same run rather than read from a resource file,
        # so it cannot be stale or paired with the wrong shard.
        Int n_variants

        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"


        RuntimeAttr? runtime_attr_override
    }

    # Per-shard memory. From phase/src/models/imputation_hmm.cpp:
    #   Alpha.resize(polymorphic_sites.size() * modK)   // float
    # modK is n_states rounded up to a multiple of 8, and n_states is bounded above by Kpbwt,
    # so the matrix costs at most 4 bytes * L * Kpbwt per thread:
    #
    #   mem = 8 + 4 * threads * Kpbwt * L / 1e9    (GB)
    #
    # The 8 GB constant covers GLIMPSE2's fixed structures with margin -- the observed fixed
    # cost was nearer 5.5, and the fit came from a handful of points near a noisy threshold.
    #
    # Calibration against the observed pass/fail boundary at 16 GiB, four threads, Kpbwt 1000:
    #   chr22 s0   L =   657,784  -> 19 GiB   (passed at 16, marginally)
    #   chr20 s4   L =   965,039  -> 24 GiB   (OOMed at 16)
    #   chr7  s12  L = 1,346,888  -> 30 GiB   (OOMed at 16)
    # and a typical shard at L ~= 400k asks for 15 GiB / 4 cpu.
    #
    # This is an upper bound rather than an exact model: n_states is frequently well below
    # Kpbwt, and GLIMPSE2's logged L counts total panel sites. Bounding in the safe direction
    # is deliberate, but it is why this is not a tight fit.
    #
    # COST: sizing per shard is what keeps fixing the OOMs from costing anything. Most shards
    # land at 4 cpu / ~15 GiB; only the densest handful widen. At n1-custom us-central1 spot
    # rates and ~35 min per shard, across 523 shards:
    #
    #   old  4 cpu / 16 GiB flat (OOM-prone)   $18.7 per batch
    #   flat 8 cpu / 40 GiB (worst-case bound) $40.6
    #   per-shard, this change                 $18.4
    #
    # Whole batch including ligate, pop and the counting task: $22.4 -> $22.8, or ~$44 -> ~$45
    # effective at the measured 1.97 full-run-equivalents per success. Over 200 batches that
    # is about +$140 against the old shape -- and about $8.7k less than the flat bound would
    # have cost. The six OOMing shards get fixed essentially for free.
    # Parenthesised to force Float promotion before the multiplications, not for readability.
    # Evaluated as integers, phase_kpbwt * n_variants alone is 1000 * 1,346,888 = 1.35e9, and
    # at Kpbwt 2000 it is 2.7e9 -- past Int32. Promoting first keeps the whole product in
    # Float. Do not reorder this expression.
    Int final_mem_gb  = 8 + ceil((((4.0 * phase_threads) * phase_kpbwt) * n_variants) / 1000000000.0)
    # N1 allows at most 6.5 GB per cpu (GcpBatchCustomMachineType.scala) and rounds any cpu
    # count other than 1 up to an even number. Doing both here keeps the requested shape
    # visible instead of letting Cromwell silently adjust it.
    Int ratio_min_cpu = ceil(final_mem_gb / 6.5)
    Int unrounded_cpu = if ratio_min_cpu > phase_threads then ratio_min_cpu else phase_threads
    Int final_cpu     = unrounded_cpu + (unrounded_cpu % 2)

    # Sized from the actual inputs rather than a flat 50 GB.
    #
    # Peak is NOT just the localized inputs. Cromwell's checkpoint sync copies the checkpoint
    # before replacing it (cp checkpoint.bin checkpoint.bin-tmp), so two full copies coexist,
    # and the reheader step writes the output a second time. For chr7 s12 that is roughly
    #   10.65 (bin + PL VCF) + ~1.3 (output, written twice) + ~4 (1.39 GB checkpoint x2) ~= 16 GB
    # against 10 + ceil(2.0 * 10.61) = 32 GB provisioned, i.e. a 2.0x margin. (size() returns
    # decimal GB, so the provisioned figure is nearer 33 GiB in practice.)
    #
    # The 30 GB floor is deliberate and costs little. Localization was measured once, at 50 GiB:
    # 8.96 GiB in 57 s (161 MiB/s), about 7x what the documented pd-ssd scaling (0.48 MiB/s per
    # GiB) predicts. That shows per-GiB scaling is not binding *at 50 GiB*; it says nothing about
    # the curve at 16-20 GiB. If the scaling does bite there, a 16 GiB disk would localize at
    # ~8 MiB/s and turn a 57 s step into ~8 min, on every one of 523 shards. The floor keeps the
    # smallest disk within striking distance of the one size actually measured until someone
    # times a shard on a small disk; raising the risk to save 10 GB is a bad trade.
    #
    # Sam's original TODO also notes that only one shard of input_vcf is ever used; pre-splitting
    # the PL VCF per chromosome upstream would shrink this further and cut ~1.4 TB of redundant
    # localization per batch. Not done here -- it is a change to the preprocessing stage.
    Int computed_disk_gb = 10 + ceil(2.0 * (size(panel_split_chunk_bin, "GB") + size(input_vcf, "GB")))
    Int disk_size_gb = if computed_disk_gb > 30 then computed_disk_gb else 30

    command <<<
        set -euxo pipefail

        # ---- resource instrumentation (best effort; must survive an OOM kill) -------------
        # Two deliberate choices, both load-bearing:
        #
        #  * Emitted from an EXIT trap, not the end of the script. Under `set -e` a SIGKILLed
        #    child aborts the command, so an end-of-script report captures nothing from
        #    precisely the runs worth measuring. The shell outlives the killed child -- which
        #    is why "line 32: 16 Killed" appears in these logs at all -- so the trap fires and
        #    still reports the peak that caused the kill.
        #  * Written to stderr, not stdout. Cromwell delocalizes stderr for FAILED tasks;
        #    File outputs are not produced at all when a task fails.
        #
        # Every read is guarded and the trap disables errexit, so instrumentation can never
        # fail the task it measures.
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
            echo "[RESOURCE] task=GLIMPSE2Phase region=~{input_region} n_variants=~{n_variants} requested_mem_gib=~{final_mem_gb} requested_cpu=~{final_cpu} threads=~{phase_threads} kpbwt=~{phase_kpbwt} rc=$_rc peak_rss_gib=$(_instr_peak) limit_gib=$(_instr_limit) disk_used_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f1) disk_total_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f2) wall_s=$((SECONDS-_INSTR_T0))" >&2
        }
        trap _instr_report EXIT
        # Coarse time series so a spike can be located against the tool's own progress markers
        # (GLIMPSE2 prints Cnk/Buf lines as it advances). One stderr line per 10 s.
        ( set +x; while :; do
            echo "[RESOURCE-TS] t=$((SECONDS-_INSTR_T0)) rss_gib=$(_instr_cur)" >&2
            sleep 10
          done ) & _INSTR_SAMPLER=$!


        # ---- resource reporting (best effort; never fails the task) -----------------------
        # Every memory figure in this file is a bound or an empirical cap: no peak RSS has ever
        # been collected anywhere in this pipeline, three different models of the ligate OOM
        # have been proposed and all three were wrong, and phase is sized from an upper bound
        # that is known to overshoot. These lines make the next batch produce the missing data
        # at no extra cost and with no dedicated profiling run.
        #
        # Harvest with:  grep -h '\[RESOURCE\]' <cromwell logs> | scripts/harvest-resources.py
        #
        # Reads the container's own cgroup rather than /usr/bin/time, which is not guaranteed
        # to be in these images, and which would not capture children. Every read is guarded so
        # instrumentation can never fail the task it is measuring.
        RESOURCE_T0=$SECONDS
        report_resources() {
            set +e
            RR_PEAK=""; RR_SRC="unavailable"
            if [ -r /sys/fs/cgroup/memory.peak ]; then
                RR_PEAK=$(cat /sys/fs/cgroup/memory.peak 2>/dev/null); RR_SRC="cgroup-v2"
            elif [ -r /sys/fs/cgroup/memory/memory.max_usage_in_bytes ]; then
                RR_PEAK=$(cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes 2>/dev/null); RR_SRC="cgroup-v1"
            fi
            RR_DISK=$(df -P -BG . 2>/dev/null | awk 'NR==2{gsub(/G/,"",$3); gsub(/G/,"",$2); print $3" "$2}')
            echo "[RESOURCE] $1 peak_rss_bytes=${RR_PEAK:-NA} peak_rss_source=${RR_SRC}" \
                 "disk_used_gb=$(echo "${RR_DISK:-NA NA}" | cut -d' ' -f1)" \
                 "disk_total_gb=$(echo "${RR_DISK:-NA NA}" | cut -d' ' -f2)" \
                 "wall_s=$((SECONDS-RESOURCE_T0))"
            set -e
        }

        # Fail fast on values that would make the memory request meaningless. GLIMPSE2 would
        # otherwise either reject these itself -- and boost program_options errors are caught
        # and exit 0, so the failure would surface much later as a missing output file -- or
        # accept them and overrun the request.
        if [ "~{phase_threads}" -lt 1 ] || [ "~{phase_kpbwt}" -lt 1 ]; then
            echo "ERROR: phase_threads and phase_kpbwt must both be >= 1 (got ~{phase_threads}, ~{phase_kpbwt})." >&2
            exit 1
        fi

        # The memory request is computed from phase_threads and phase_kpbwt (see the sizing
        # block above), so both have to come from that one place. Existing input CSVs carry
        # --thread inside extra_phase_args, and the previous default string carried --Kpbwt.
        # Passing either twice hands GLIMPSE2 a duplicate option, which boost program_options
        # rejects -- and because GLIMPSE2 catches parser errors and exits 0, that surfaces later
        # as a confusing missing-BCF failure rather than a clean error. Honouring the string's
        # value instead is worse still: it decouples the parameter from the memory sized
        # against it. So strip both from extra_phase_args, warn, and inject the typed values.
        # NOTE: the placeholder is substituted into this script as literal text, so bash
        # expands anything expandable inside it at assignment time. The historical default
        # contained "--thread $(nproc)", which therefore becomes "--thread 8" here before the
        # stripping below ever sees it. Both forms are handled -- the value matcher accepts any
        # non-'-' token, so it removes "$(nproc)" and "8" alike -- but the distinction matters
        # when reasoning about this block, and the tests cover both spellings for that reason.
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
                # Three passes. The first takes a signed numeric value, so "--thread -1"
                # does not leave a stray "-1" behind for GLIMPSE2 to read as a positional.
                # The second takes any value not starting with '-' -- it must match more than
                # [0-9]+, because the historical default used "--thread $(nproc)". The third
                # removes a valueless leftover such as a trailing "--thread".
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

        # requested_* are what the sizing arithmetic asked for; peak_rss_bytes is what it
        # actually needed. The gap between them is the whole point of emitting this.
        report_resources "task=GLIMPSE2Phase region=~{input_region} n_variants=~{n_variants} \
requested_mem_gib=~{final_mem_gb} requested_cpu=~{final_cpu} threads=~{phase_threads} \
kpbwt=~{phase_kpbwt} resumed=$RESUMED"
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
    # runtime_attr_override is applied field by field, so a caller who overrides only mem_gb
    # would otherwise keep a cpu_cores computed for the DEFAULT memory -- e.g. mem_gb 64 with
    # the default 4 cpu is 16 GB/cpu, far past the N1 limit, and Cromwell would silently widen
    # the CPU count. That is the exact hidden adjustment the sizing arithmetic exists to
    # prevent, so the ratio is re-derived here from whatever memory actually won, and an
    # explicit cpu_cores override still takes precedence over the derivation.
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

        # ---- resource instrumentation (best effort; must survive an OOM kill) -------------
        # Two deliberate choices, both load-bearing:
        #
        #  * Emitted from an EXIT trap, not the end of the script. Under `set -e` a SIGKILLed
        #    child aborts the command, so an end-of-script report captures nothing from
        #    precisely the runs worth measuring. The shell outlives the killed child -- which
        #    is why "line 32: 16 Killed" appears in these logs at all -- so the trap fires and
        #    still reports the peak that caused the kill.
        #  * Written to stderr, not stdout. Cromwell delocalizes stderr for FAILED tasks;
        #    File outputs are not produced at all when a task fails.
        #
        # Every read is guarded and the trap disables errexit, so instrumentation can never
        # fail the task it measures.
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
            echo "[RESOURCE] task=GLIMPSE2Ligate n_shards=~{length(phased_vcfs)} requested_mem_gib=32 requested_cpu=6 threads=2 rc=$_rc peak_rss_gib=$(_instr_peak) limit_gib=$(_instr_limit) disk_used_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f1) disk_total_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f2) wall_s=$((SECONDS-_INSTR_T0))" >&2
        }
        trap _instr_report EXIT
        # Coarse time series so a spike can be located against the tool's own progress markers
        # (GLIMPSE2 prints Cnk/Buf lines as it advances). One stderr line per 10 s.
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
    # 32 GiB is an empirical cap, NOT a validated model. What is established: ligate was
    # SIGKILLed (rc=137) at 12 GiB on chr2 and chr20, twice each, and chr7 needed three
    # attempts; chr22 passed. The failures correlate with the number of variants in the
    # densest seam of each chromosome, and 12 GiB sits between chr22's worst seam (245k,
    # passed) and chr2's (597k, failed).
    #
    # An earlier version of this comment asserted a mechanism -- that the two synced readers
    # held every seam record, giving mem ~= 1.5 + 18*L_isec/1e6. That does not survive
    # reading HTSlib 1.16: _reader_fill_buffer only buffers records sharing a coordinate and
    # stops at the first record with a different one, reusing the buffer as it advances. So
    # the linear-in-seam-length model is wrong, and the correlation is probably confounded
    # (denser seams also carry more multiallelic sites, and the highly multiallelic bubbles
    # are a known sore spot in this panel).
    #
    # What survives is only that something seam-specific is real: the failing seam was
    # predicted correctly three times out of three (chr20 seam 4, chr2 seam 17, chr7 seam 11),
    # computed from chunks.tsv before the outcomes were known. TWO mechanisms have been
    # proposed to explain that and BOTH are dead. Recorded here so neither gets re-derived:
    #
    #  1. Linear in seam length, mem ~= 1.5 + 18*L_isec/1e6. Killed by HTSlib 1.16:
    #     _reader_fill_buffer only buffers records sharing a coordinate and stops at the first
    #     record with a different one, reusing the buffer as it advances. Nothing accumulates
    #     across the seam.
    #
    #  2. Max records per coordinate (large split multiallelic bubbles piling many ALT records
    #     at one position). Killed on both separation and magnitude:
    #       - chr20 peaked at 10,465 records at one position and FAILED; chr22 peaked at
    #         9,937 and PASSED. Five percent apart, opposite outcomes. Only chr2 (24,745)
    #         stands out at all.
    #       - magnitude is off by ~100x: 10,465 records x ~9 KB x 2 readers is ~190 MB, not
    #         the ~12 GB that would have to be explained.
    #
    # So the cause is genuinely unidentified. Remaining candidates -- buffer high-water marks,
    # index or output construction, allocator behaviour, something outside the ligater process
    # entirely -- are not distinguishable from rc=137. The only thing that would settle it is
    # peak RSS from /usr/bin/time -v around one known-failing seam.
    #
    # Until then this is a cap chosen to stop the bleeding, and it should not be tuned in
    # either direction on the strength of a model, because every model tried so far has been
    # wrong.
    #
    # cpu 6 exists only to keep 32/6 = 5.33 under the N1 6.5 GB/cpu limit; with --thread 2
    # four of those cores are deliberately idle. That is the price of memory on this shape.
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
    # runtime_attr_override is applied field by field, so a caller who overrides only mem_gb
    # would otherwise keep a cpu_cores computed for the DEFAULT memory -- e.g. mem_gb 64 with
    # the default 4 cpu is 16 GB/cpu, far past the N1 limit, and Cromwell would silently widen
    # the CPU count. That is the exact hidden adjustment the sizing arithmetic exists to
    # prevent, so the ratio is re-derived here from whatever memory actually won, and an
    # explicit cpu_cores override still takes precedence over the derivation.
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

        # ---- resource instrumentation (best effort; must survive an OOM kill) -------------
        # Two deliberate choices, both load-bearing:
        #
        #  * Emitted from an EXIT trap, not the end of the script. Under `set -e` a SIGKILLed
        #    child aborts the command, so an end-of-script report captures nothing from
        #    precisely the runs worth measuring. The shell outlives the killed child -- which
        #    is why "line 32: 16 Killed" appears in these logs at all -- so the trap fires and
        #    still reports the peak that caused the kill.
        #  * Written to stderr, not stdout. Cromwell delocalizes stderr for FAILED tasks;
        #    File outputs are not produced at all when a task fails.
        #
        # Every read is guarded and the trap disables errexit, so instrumentation can never
        # fail the task it measures.
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
            echo "[RESOURCE] task=PopAndMarginalizeCollisions region=~{region} requested_mem_gib=12 requested_cpu=2 rc=$_rc peak_rss_gib=$(_instr_peak) limit_gib=$(_instr_limit) disk_used_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f1) disk_total_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f2) wall_s=$((SECONDS-_INSTR_T0))" >&2
        }
        trap _instr_report EXIT
        # Coarse time series so a spike can be located against the tool's own progress markers
        # (GLIMPSE2 prints Cnk/Buf lines as it advances). One stderr line per 10 s.
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

        # Deliberately left at 12 GiB with no evidence either way. This is the evidence.
        report_resources "task=PopAndMarginalizeCollisions region=~{region} requested_mem_gib=12 requested_cpu=2"
    >>>

    output {
        File popped_vcf = "~{output_prefix}.bcf"
        File popped_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    # DELIBERATELY UNCHANGED at 12 GiB / 2 cpu (6.0 GB/cpu, within the N1 limit).
    #
    # An earlier revision raised this to 24 GiB on the reasoning that chr2 and chr20 -- the two
    # chromosomes with the densest seams -- have never reached this stage, ligate having failed
    # first, so it is untested at the sizes that matter. That observation stands. The
    # justification did not: it rested on pop being a couple of dozen tasks per batch, where
    # over-provisioning is free.
    #
    # It is not. pop_regions defaults to output_regions, so this scatters once per shard --
    # ~523 tasks per batch, the same order as phase, not the 22 that ligate runs. Doubling the
    # memory therefore costs roughly $3.5 per batch, ~$700 across the 50k, to guard against a
    # failure that has never been observed in any pop task on any chromosome.
    #
    # "Downstream of something that failed" is also not evidence of anything: it is equally
    # true of every task after phase. Raising this on that basis would mean raising all of
    # them. Left alone; if a pop task ever does OOM, that is the evidence to act on.
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
    # runtime_attr_override is applied field by field, so a caller who overrides only mem_gb
    # would otherwise keep a cpu_cores computed for the DEFAULT memory -- e.g. mem_gb 64 with
    # the default 4 cpu is 16 GB/cpu, far past the N1 limit, and Cromwell would silently widen
    # the CPU count. That is the exact hidden adjustment the sizing arithmetic exists to
    # prevent, so the ratio is re-derived here from whatever memory actually won, and an
    # explicit cpu_cores override still takes precedence over the derivation.
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
