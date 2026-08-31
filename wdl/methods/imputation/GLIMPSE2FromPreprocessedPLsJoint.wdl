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

        # Workflow-level, because the phase memory request is computed from both.
        Int phase_threads = 4
        Int phase_kpbwt = 1000
        # MUST match the GLIMPSE2Ligate task default. The workflow passes this through, so
        # this value is the one that runs and the task default is inert. 24 GiB is sized
        # against a peak measured at 2 threads.
        Int ligate_threads = 2

        # Phase now receives a pre-split regional PL BCF rather than localizing the full
        # chromosome BCF in every shard. The smaller input removes the reason for the old
        # 30-GB throughput floor. 20 GB retains space for the panel bin, phase output and
        # reheader copy while saving about 5.2 TB of provisioned work disk per genome batch.
        Int phase_disk_floor_gb = 20
        # GLIMPSE2 is not deterministic under multithreading; pinning threads fixes the
        # thread COUNT, not the output. The size of that non-determinism has now been measured
        # two ways, and the two disagree by three orders of magnitude.
        #
        # Controlled replicates, leaveout panel, --Kpbwt 1000 confirmed in the log, 8 threads:
        #     mean dDS      -1.18e-05     sd 1.26e-02     |mu|/sd 9.3e-04
        #     GT flips       0.014%
        #     bootstrap CI straddles zero
        # That is the real floor: replicate-to-replicate noise is centred on zero and tiny.
        #
        # The earlier chr22 pair gave a much larger number -- 4.671% non-ref discordance over
        # 1,896,739 sites x 250 samples, 99.822% hard-call concordance, 1.226% phase flips.
        # Do NOT read that as the determinism floor. Those two runs also differed in
        # preemptible_tries, so checkpoint-restart was confounded with thread scheduling, and
        # they did not pin Kpbwt, so they ran at the 2000-state default rather than 1000. It
        # bounds a combined effect, not thread non-determinism.
        #
        # The practical consequence: a Terra-vs-VWB comparison has far more resolution than
        # 4.7% implies. Differences at the percent level are signal, not noise.

        String output_prefix

        # inputs for PopAndMarginalizeCollisions
        File pop_glimpse2_panel_resources_json
        # Built once and reused. Compiling in each pop task would repeat the same Cargo build
        # thousands of times in production.
        File pop_glimpse2_binary

        # Pinned by digest. GCR tags are mutable, so a repush would silently
        # change what runs. This digest is tag 1.0.0-2cee597-1778869818 as of 2026-07-31.
        # Enables checkpointing; note it contains bcftools/htslib 1.16.
        #
        # NB the staged input CSVs override this with :tfenne-opt, which IS a mutable tag --
        # pinning it there matters more than pinning the default here, and that is a change to
        # the run config.
        String glimpse2_docker = "us.gcr.io/broad-gotc-prod/imputation-glimpse2@sha256:c3d64c5af3b8e789bcda94451931f09073d3a8750fc7375a8d44246d193e8a21"

        # VWB's Batch backend is restricted to us-central1-a. Task-level attempts to provide
        # a wider zone pool were measured to be inert, so state the location that actually
        # runs instead of advertising unavailable zones.
        String zones = "us-central1-a"
    }

    Map[String, String] genetic_maps_dict = read_map(genetic_maps_tsv)
    String genetic_map = genetic_maps_dict[chromosome]

    Map[String, ChunkedPanelChromosome] chunked_panel = read_json(chunked_panel_json)
    Array[String] input_regions = chunked_panel[chromosome].input_regions
    Array[String] output_regions = chunked_panel[chromosome].output_regions
    Array[File] panel_split_chunk_bins = chunked_panel[chromosome].panel_split_chunk_bins
    # Panel variants per shard, parallel to input_regions. Required: an absent field fails at
    # input evaluation and a short array fails on the index below, both before any VM starts.
    Array[Int] n_variants = chunked_panel[chromosome].n_variants

    Map[String, PopAndMarginalizePanelResourcesChromosome] pop_glimpse2_panel_resources = read_json(pop_glimpse2_panel_resources_json)
    File panel_bubble_split_sites_only_vcf = pop_glimpse2_panel_resources[chromosome].panel_bubble_split_sites_only_vcf
    File panel_bubble_split_sites_only_vcf_idx = pop_glimpse2_panel_resources[chromosome].panel_bubble_split_sites_only_vcf_idx
    File panel_id_split_vcf_gz = pop_glimpse2_panel_resources[chromosome].panel_id_split_vcf_gz
    File panel_id_split_vcf_gz_tbi = pop_glimpse2_panel_resources[chromosome].panel_id_split_vcf_gz_tbi
    Array[String] pop_regions = select_first([pop_glimpse2_panel_resources[chromosome].pop_regions, output_regions])

    call SplitPreprocessedPLsForPhase {
        input:
            input_vcf = input_preprocessed_joint_vcf,
            input_vcf_idx = input_preprocessed_joint_vcf_idx,
            input_regions = input_regions,
            zones = zones,
            docker = glimpse2_docker
    }

    scatter (k in range(length(output_regions))) {
        call GLIMPSE2Phase as ChunkedGLIMPSE2Phase {
            input:
                input_vcf = SplitPreprocessedPLsForPhase.sharded_vcfs[k],
                input_vcf_idx = SplitPreprocessedPLsForPhase.sharded_vcf_idxs[k],
                panel_split_chunk_bin = panel_split_chunk_bins[k],
                input_region = input_regions[k],
                output_region = output_regions[k],
                genetic_map = genetic_map,
                output_prefix = output_prefix + ".shard-" + k + ".glimpse2.phased",
                extra_phase_args = extra_phase_args,
                phase_threads = phase_threads,
                phase_kpbwt = phase_kpbwt,
                phase_disk_floor_gb = phase_disk_floor_gb,
                n_variants = n_variants[k],
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
    # EXTRA boot disk added to the backend default, so 0 is correct.
    # From BootDiskSizeValidation in Cromwell 91:
    #     override protected def staticDefaultOption: Option[WomInteger] = Option(WomInteger(0))
    #     case WomInteger(value) => (value + defaultBootDiskSize).validNel
    # so an absent attribute takes the static default of 0 and an explicit 0 takes the same
    # path -- both yield 0 + 30 = 30 GB, and there is no positive-int guard to trip. That is
    # why the old value of 10 provisioned 40. CONFIRMED on a live job: bootDiskMib 28611
    # (30 GB), against 38148 (40 GB) before -- 10 GB per shard, ~5.2 TB per batch.
    #
    # Passing 0 explicitly rather than omitting the attribute keeps this override usable by a
    # caller running a larger custom Docker image, at identical default behaviour.
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

struct PopAndMarginalizePanelResourcesChromosome {
    String panel_bubble_split_sites_only_vcf
    String panel_bubble_split_sites_only_vcf_idx
    String panel_id_split_vcf_gz
    String panel_id_split_vcf_gz_tbi
    Array[String]? pop_regions
}

task SplitPreprocessedPLsForPhase {
    input {
        File input_vcf
        File input_vcf_idx
        Array[String] input_regions
        String docker
        String zones = "us-central1-a"

        RuntimeAttr? runtime_attr_override
    }

    # The task localizes the chromosome BCF once and writes overlapping regional BCFs whose
    # total size is modestly larger than the input because GLIMPSE input regions include
    # buffers. Three input sizes plus 10 GB covers the localized input, all regional outputs,
    # their indexes, and compression overhead.
    Int disk_size_gb = 10 + 3 * ceil(size(input_vcf, "GB"))

    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       0,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             docker
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    Float eff_mem_gb  = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
    Int eff_ratio_cpu = ceil(eff_mem_gb / 6.5)
    Int eff_unrounded = if eff_ratio_cpu > 4 then eff_ratio_cpu else 4
    Int eff_cpu       = select_first([runtime_attr.cpu_cores, eff_unrounded + (eff_unrounded % 2)])

    command <<<
        set -euxo pipefail

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
            _du=$(df -Pk . 2>/dev/null | awk 'NR==2{printf "%.0f %.0f", $3/1048576, $2/1048576}')
            echo "[RESOURCE] task=SplitPreprocessedPLsForPhase requested_mem_gib=~{eff_mem_gb} requested_cpu=~{eff_cpu} rc=$_rc peak_rss_gib=$(_instr_peak) limit_gib=$(_instr_limit) disk_used_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f1) disk_total_gb=$(echo "${_du:-NA NA}" | cut -d' ' -f2) wall_s=$((SECONDS-_INSTR_T0))" >&2
        }
        trap _instr_report EXIT
        ( set +x; while :; do
            echo "[RESOURCE-TS] t=$((SECONDS-_INSTR_T0)) rss_gib=$(_instr_cur)" >&2
            sleep 10
          done ) & _INSTR_SAMPLER=$!

        mkdir phase-pl-shards
        : > phase-pl-shards/vcfs.list
        : > phase-pl-shards/indexes.list
        i=0
        while IFS= read -r REGION; do
            [ -n "$REGION" ] || continue
            PREFIX=$(printf "phase-pl-shards/shard-%04d" "$i")
            bcftools view \
                --threads 3 \
                --regions-overlap 0 \
                --regions "$REGION" \
                --output-type b \
                --output "$PREFIX.bcf" \
                ~{input_vcf}
            bcftools index --threads 3 --force "$PREFIX.bcf"
            printf '%s\n' "$PREFIX.bcf" >> phase-pl-shards/vcfs.list
            printf '%s\n' "$PREFIX.bcf.csi" >> phase-pl-shards/indexes.list
            i=$((i + 1))
        done < ~{write_lines(input_regions)}

        if [ "$i" -ne ~{length(input_regions)} ]; then
            echo "ERROR: wrote $i phase PL shards for ~{length(input_regions)} input regions" >&2
            exit 1
        fi
    >>>

    output {
        Array[File] sharded_vcfs = read_lines("phase-pl-shards/vcfs.list")
        Array[File] sharded_vcf_idxs = read_lines("phase-pl-shards/indexes.list")
    }

    runtime {
        cpu:                    eff_cpu
        memory:                 eff_mem_gb + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker, default_attr.docker])
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

        # Typed inputs, because the memory request is computed from
        # them. A caller replacing that string drops --Kpbwt 1000, gets this build's default
        # of 2000, and doubles the matrix while the request stays put. Threads are pinned
        # rather than $(nproc) because the matrix is per-thread: under $(nproc) memory depends
        # on cpu while the N1 ratio makes cpu depend on memory, unsolvable above L ~= 1.6M.
        # cpu_cores may exceed threads -- ratio headroom, not parallelism.
        #
        # 4 is the measured configuration. Modelling wall as
        # 415 + 0.00267*L s makes 8 threads clearly worse (+$2.11/batch) and 2 threads
        # marginally worse (+$1.23) -- but that assumes the parallel part scales perfectly,
        # which the memory axis does not. It flips to cheaper if 4-thread
        # efficiency is below ~92% of 2-thread. Worth one chromosome at phase_threads=2
        # before treating 4 as settled; it also gives a 2 cpu / 9 GiB shape, which places
        # better on spot.
        Int phase_threads = 4
        Int phase_kpbwt = 1000

        # The workflow pre-splits the chromosome PL BCF, so phase no longer needs the old
        # 30-GB floor to localize a multi-gigabyte chromosome input.
        Int phase_disk_floor_gb = 20

        # Panel variants in this shard's input region -- GLIMPSE2's L. Required, and supplied
        # by the panel, parallel to input_regions, so it is bound to the shard it describes.
        Int n_variants

        String zones = "us-central1-a"


        RuntimeAttr? runtime_attr_override
    }

    # Fitted from measured peak RSS. The chr20 sweep gave 11 points at 4 threads / Kpbwt 1000:
    #     peak_GiB = 0.42 + 2.17e-5 * L        (r2 ~ 0.99)
    # i.e. 5.82 bytes per site per thread per state. Requesting 8.0 on the strength of that
    # then OOMed every shard of the first chr2/chr11 run -- 3 of 3, all rc=137:
    #     L =  50,227  requested  4 GiB, killed at  3.37
    #     L = 351,339  requested 14 GiB, killed at 13.15
    #     L = 357,973  requested 14 GiB, killed at 13.15
    # Refitting on those: peak = 1.75 + 3.21e-5 * L, i.e. 8.62 bytes -- above the 8.0 they
    # were sized with. The chr20 fit understated the slope by a third, so a sweep of one
    # chromosome does not generalise to another.
    #
    # Two corrections, because the failures show two things at once. 8.0 -> 11.0 covers the
    # refit with margin, and 2 -> 4 GiB pays for what the task cannot use: every kill landed
    # BELOW its request (84%, 94%, 94%), so the guest OS, Docker, the monitoring script and
    # localization buffers take a cut of the VM that GLIMPSE2 never sees.
    #
    # Both peaks and the refit are lower bounds. peak_rss_gib on an rc=137 task is where the
    # kernel stopped it, not what it wanted, which biases the slope down -- the direction that
    # OOMs. Treat 11.0 as the floor of what is defensible until shards finish.
    #     L =   400,379 -> 22 GiB / 4 cpu
    #     L = 1,005,104 -> 49 GiB / 8 cpu
    #     L = 1,346,888 -> 64 GiB / 10 cpu
    # Per-shard sizing redistributes memory rather than saving it against the old flat 16. 490 shards get less, 33 get more, and the total is 3538 GiB-hours against 3552, a
    # 0.4% difference. Phase compute is $13.68/batch either way (+$0.09). What it does buy is
    # the freedom to size the dense shards correctly at all: a flat request safe for chr7 s12
    # would be 46 GiB / 8 cpu everywhere, $31.32/batch, so per-shard is $17.65 cheaper than
    # the only flat alternative that does not OOM.
    #
    # So the OOM fix is free. The batch does get ~$1.35 cheaper, and that comes from disk
    # (50 -> 30 GB working, 40 -> 30 boot).
    #
    # Two caveats. All 11 points are at threads=4 and Kpbwt=1000 and the formula multiplies
    # by both, so neither scaling is measured. The thread one is sublinear:
    # since linear predicts 58.8 GiB for chr7 s12 at 8 threads and it succeeded on 40. And
    # cgroup memory.peak counts page cache, so this is an upper bound on demand.
    #
    # Largest untaken dollar lever, not in this file: cpuPlatform is unset, so these tasks
    # land on N1CustomMachineType. N2D is ~14% cheaper (~$2.26/batch). It cannot be an
    # optional input -- Cromwell's validate keys off attribute presence, so a runtime key
    # wired to a String is always present and an empty string reaches Batch as
    # minCpuPlatform="". Testing it means editing the runtime block. Not default because
    # setMinCpuPlatform restricts placement, which at 41-72% preemption may cost more than
    # it saves.
    # Floored at 16 GiB. The regression cannot be trusted to set the low end: it was fit to
    # peak_rss values from runs killed at rc=137, which are censored observations -- a kill
    # point is where measurement stopped, not what the task wanted -- and those runs did not
    # pin Kpbwt, so they ran at the 2000-state default while the coefficient is applied here
    # at 1000. Both biases push the fit down. Against that, 517 of 523 shards completed at a
    # flat 16 GiB, which is a measurement rather than an extrapolation. Keep the regression
    # for the large shards it was built for and let the measured value hold the floor.
    Int computed_mem_gb = 4 + ceil((((11.0 * phase_threads) * phase_kpbwt) * n_variants) / 1000000000.0)
    Int final_mem_gb    = if computed_mem_gb > 16 then computed_mem_gb else 16
    # N1 allows <= 6.5 GB/cpu and rounds any cpu count but 1 up to even; doing both here keeps
    # the requested shape visible instead of letting Cromwell adjust it silently.
    Int ratio_min_cpu = ceil(final_mem_gb / 6.5)
    Int unrounded_cpu = if ratio_min_cpu > phase_threads then ratio_min_cpu else phase_threads
    Int final_cpu     = unrounded_cpu + (unrounded_cpu % 2)

    # Sized from the regional inputs rather than a flat 50 GB. The old 17-GB measured
    # high-water included localization of the full chromosome PL BCF into every phase shard.
    # SplitPreprocessedPLsForPhase removes that repeated input. The 20-GB floor retains room
    # for the panel bin, raw output, reheader copy, index, and checkpoint while reducing
    # provisioned SSD by about 5.2 TB per genome batch relative to 30 GB.
    Int computed_disk_gb = 10 + ceil(2.0 * (size(panel_split_chunk_bin, "GB") + size(input_vcf, "GB")))
    Int disk_size_gb = if computed_disk_gb > phase_disk_floor_gb then computed_disk_gb else phase_disk_floor_gb


    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          final_cpu,
        mem_gb:             final_mem_gb,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       0,
        use_ssd:            true,
        # Spot is worth it, measured: 415 successful VM-minutes against 286 wasted over 32
        # preemptions (1.69x VM-time) still costs ~0.58x on-demand. Not 0.51x -- spot
        # discounts vCPU and memory to 30% but disk bills the same either way and the wasted
        # attempts pay for their disk:
        #     1.69 * (0.0599 + 0.0140) / (0.1996 + 0.0140) = 0.58
        # Checkpointing works: 5 of 11 shards resumed, averaging 954 s against 1873 s fresh.
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
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
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
        # phase floors on phase_threads for the same reason.
        #
        # 2, matching the configuration the 12.39 GiB peak was measured under.
        # bcf_sr_set_threads allocates per-thread decompression buffers, so a different thread
        # count needs a fresh peak before the 24 GiB request applies. The instrumentation
        # reports one on every run.
        Int ligate_threads = 2

        String zones = "us-central1-a"


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
    # Seam size drives it, at ~4 KB per intersecting site on a ~9.6 GiB baseline.
    # Max-records-per-coordinate does not: chr20 peaked at 10,465 records and failed while
    # chr22 peaked at 9,937 and passed.
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
    # chromosome reports a peak.
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
    # memory.peak includes page cache, so part of 12.39 GiB is reclaimable.
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             24,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       0,
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
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
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
        
        File pop_glimpse2_binary
        
        String region
        String output_prefix

        String zones = "us-central1-a"


        RuntimeAttr? runtime_attr_override
    }

    # Measured on chr20: 2 GB used of the 16 this produced. Two input sizes plus 10 GB leaves
    # about 5x the observed high-water while keeping the short pop shards cheap to schedule.
    Int disk_gb = 10 + 2 * ceil(size([posteriors_vcf, panel_bubble_split_sites_only_vcf, panel_id_split_vcf_gz], "GB"))


    #########################
    # 8 GiB / 2 cpu, from eleven instrumented chr20 shards whose worst peak was 4.68 GiB.
    # Retain regional pop parallelism: collapsing to one chromosome task would serialize
    # roughly 20-30 short shards and turn a 5-8 minute Spot task into a long uncheckpointed
    # critical path. bcftools sort --max-mem=2G sets the other material memory floor.
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_gb,
        boot_disk_gb:       0,
        use_ssd:            true,
        # Short pop shards restart cheaply. Four Spot attempts was the measured cost optimum
        # at the observed 40-60% preemption rate before falling through to on-demand.
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


        POP_BIN="~{pop_glimpse2_binary}"
        chmod +x "$POP_BIN"

        # this now only works for pop-glimpse2-joint-opt.rs;
        # the sort may also be extraneous, but we keep it in to guard against getting out of sync with the popped panel
        # --threads is ADDITIONAL worker threads, so 1 is 2 total on this task's 2 cpu.
        bcftools view --threads 1 -r ~{region} --regions-overlap 0 ~{panel_bubble_split_sites_only_vcf} -Oz -o panel.bubble.split.sites.shard.vcf.gz
        bcftools view --threads 1 -r ~{region} --regions-overlap 0 ~{posteriors_vcf} | \
            "$POP_BIN" ~{panel_id_split_vcf_gz} panel.bubble.split.sites.shard.vcf.gz | \
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
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
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

        String zones = "us-central1-a"


        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 2 * ceil(size(vcf, "GB"))


    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       0,
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
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        zones:                  zones
    }
}
