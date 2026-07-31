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
        Array[String] zones = ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]
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
    # Extra boot disk *on top of* the backend default, not an absolute size. Cromwell adds its
    # 30 GB default to whatever is requested here, so 0 yields the 30 GB minimum and the old
    # value of 10 yielded 40. Defaulted to 0 in every task below rather than dropped, so that
    # an override still works for callers running a larger custom Docker image.
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

    # NOTE: there is deliberately no per-shard variant count here yet.
    #
    # Sizing GLIMPSE2Phase from each shard's real L (rather than from a worst-case bound)
    # would roughly halve per-batch vCPU, since only a handful of shards need the largest
    # shape. It was implemented and then removed, because doing it safely needs two things
    # this schema cannot express:
    #
    #  1. The producer must emit the counts. GLIMPSE2ChunkAndSplitPanel builds this struct
    #     with four fields and no count, so an optional count is dead on the normal
    #     producer-to-consumer path -- it can only be supplied by hand-editing a generated
    #     resource.
    #  2. The counts must be bound to the shards they describe. These are unkeyed parallel
    #     arrays, so a count array that is stale but happens to have the right length pairs
    #     silently with the wrong bins -- and underprovisioning memory is exactly the failure
    #     this was meant to prevent. A length check does not catch reordered shards, a
    #     regenerated panel with the same shard count, or counts copied from another
    #     chromosome.
    #
    # The right shape is one array of shard objects carrying region, bin and count together,
    # emitted atomically by the producer and required rather than optional. That is a change
    # to both files and to every staged chunked_panel.json, so it is deliberately not bundled
    # with the memory fixes.
}

struct PopAndMarginalizePanelResourcesChromosome {
    String panel_bubble_split_sites_only_vcf
    String panel_bubble_split_sites_only_vcf_idx
    String panel_id_split_vcf_gz
    String panel_id_split_vcf_gz_tbi
    Array[String]? pop_regions              # non-overlapping, if not provided then GLIMPSE2 chunks will be used
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

        Array[String] zones = ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]


        RuntimeAttr? runtime_attr_override
    }

    # 40 GiB / 8 cpu at the defaults, which is the shape that ran 217 shards without an OOM.
    # The expression exists so that raising phase_kpbwt or phase_threads scales the request
    # instead of quietly invalidating it: at phase_kpbwt=1000 and phase_threads=4 the second
    # term is exactly 32, giving 40.
    #
    # Upper bound, not an exact model. phase/src/models/imputation_hmm.cpp allocates
    #   Alpha.resize(polymorphic_sites.size() * modK)   // float
    # where modK is n_states rounded up to a multiple of 8 and n_states is the number of
    # states actually selected for the shard, bounded above by Kpbwt (see
    # containers/conditioning_set.cpp). GLIMPSE2's logged L is total panel sites, so
    # L * Kpbwt * 4 bytes * threads over-counts. That is deliberate -- it is the safe
    # direction -- but it means the coefficient is a bound, not an invariant.
    #
    # Per-shard sizing from the real L would cut this substantially (most shards need far
    # less than the worst one), but it requires the chunked-panel producer to emit a
    # per-shard variant count. It does not, so that is a separate change; see the note on
    # ChunkedPanelChromosome.
    Int final_mem_gb  = 8 + ceil(32.0 * phase_kpbwt * phase_threads / 4000.0)
    # N1 allows at most 6.5 GB per cpu (GcpBatchCustomMachineType.scala) and rounds odd cpu
    # counts up to even. Doing it here keeps the requested shape visible instead of letting
    # Cromwell silently adjust it.
    Int ratio_min_cpu = ceil(final_mem_gb / 6.5)
    Int final_cpu     = if ratio_min_cpu > phase_threads then ratio_min_cpu + (ratio_min_cpu % 2) else phase_threads

    # Sized from the actual inputs rather than a flat 50 GB.
    #
    # Peak is NOT just the localized inputs. Cromwell's checkpoint sync copies the checkpoint
    # before replacing it (cp checkpoint.bin checkpoint.bin-tmp), so two full copies coexist,
    # and the reheader step writes the output a second time. For chr7 s12 that is roughly
    #   10.65 (bin + PL VCF) + ~1.3 (output, written twice) + ~4 (1.39 GB checkpoint x2) ~= 16 GB
    # against 26 GB provisioned -- a 1.6x margin, not the 2x the input-only view suggests.
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
        EXTRA_PHASE_ARGS="~{extra_phase_args}"
        for OPT in thread Kpbwt; do
            if echo "$EXTRA_PHASE_ARGS" | grep -qE "(^|[[:space:]])--${OPT}([[:space:]]|=)"; then
                echo "WARNING: --${OPT} found in extra_phase_args; overriding with the typed input." >&2
                echo "WARNING: set phase_threads / phase_kpbwt instead -- memory is sized from them." >&2
                EXTRA_PHASE_ARGS=$(echo "$EXTRA_PHASE_ARGS" | sed -E "s/(^|[[:space:]])--${OPT}([[:space:]]+|=)[0-9]+/ /g")
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

        if [ -s "checkpoint.bin" ]; then
            cmd="$cmd --checkpoint-file-in checkpoint.bin" 
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
        boot_disk_gb:       0,
        use_ssd:            true,
        preemptible_tries:  10,
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

        Array[String] zones = ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]


        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 2 * ceil(size(phased_vcfs, "GB")) + 10

    command <<<
        set -euox pipefail

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
    # The real cause is unidentified: candidates include per-coordinate record pileups at
    # large bubbles, buffer high-water marks, index/output construction, or allocator
    # behaviour. rc=137 alone does not distinguish them. Treat this as "big enough to stop
    # the bleeding" and measure peak RSS on a known-failing seam before tuning further.
    #
    # cpu 6 exists only to keep 32/6 = 5.33 under the N1 6.5 GB/cpu limit; with --thread 2
    # four of those cores are deliberately idle. That is the price of memory on this shape.
    RuntimeAttr default_attr = object {
        cpu_cores:          6,
        mem_gb:             32,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       0,
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

        Array[String] zones = ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]


        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 3 * ceil(size([posteriors_vcf, panel_bubble_split_sites_only_vcf, panel_id_split_vcf_gz], "GB"))

    command <<<
        set -euox pipefail

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
    # NOTE: unlike phase and ligate, 24 is a round number, not a derivation. No pop task has
    # ever OOMed and there is no memory model for this step. It is raised only because chr2 and
    # chr20 -- the two chromosomes with the densest seams -- have never reached it, ligate
    # having failed first, so it is untested at exactly the sizes that matter.
    # cpu_cores is 4 rather than 2 to keep 24/4 = 6.0 under the N1 6.5 GB/cpu limit; at 2 cpu
    # Cromwell would compute ceil(24/6.5) = 4 and silently substitute the same shape anyway.
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             24,
        disk_gb:            disk_gb,
        boot_disk_gb:       0,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsde-methods/slee/lrma-aou2-panel-creation-rust:v1"
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
        zones:                  zones
    }
}

task RemapSampleNames {
    input {
        File vcf
        File vcf_idx
        File remap_file
        String output_prefix

        Array[String] zones = ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]


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
        boot_disk_gb:       0,
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
        zones:                  zones
    }
}
