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

        String extra_phase_args = "--impute-reference-only-variants --keep-monomorphic-ref-sites --Kpbwt 1000 --main 10 --burnin 5 --err-imp 1E-3"
        String output_prefix

        # inputs for PopAndMarginalizeCollisions
        File pop_glimpse2_panel_resources_json
        File? pop_glimpse2_script               # heavily modified version of convert-to-biallelic.py
        File? pop_glimpse2_cargo_toml
        File? pop_glimpse2_binary

        String glimpse2_docker = "us.gcr.io/broad-gotc-prod/imputation-glimpse2:1.0.0-2cee597-1778869818"    # enables checkpointing, but note this contains bcftools/htslib 1.16!
    }

    Map[String, String] genetic_maps_dict = read_map(genetic_maps_tsv)
    String genetic_map = genetic_maps_dict[chromosome]

    Map[String, ChunkedPanelChromosome] chunked_panel = read_json(chunked_panel_json)
    Array[String] input_regions = chunked_panel[chromosome].input_regions
    Array[String] output_regions = chunked_panel[chromosome].output_regions
    Array[File] panel_split_chunk_bins = chunked_panel[chromosome].panel_split_chunk_bins

    # Empty unless the panel JSON supplies n_variants. The length is checked against the shard
    # count rather than trusted, so a stale or partial array degrades to static sizing instead
    # of silently pairing shards with the wrong L.
    # select_all + flatten rather than `if defined(...) then select_first(...) else []`: the
    # ternary form relies on Cromwell not evaluating the untaken branch, and select_first on an
    # absent optional throws. select_all drops the None without ever unwrapping it, yielding a
    # 0- or 1-element Array[Array[Int]] that flattens to the values or to [].
    Array[Int] panel_n_variants = flatten(select_all([chunked_panel[chromosome].n_variants]))
    Boolean use_dynamic_phase_mem = length(panel_n_variants) == length(output_regions)

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
                n_variants = if use_dynamic_phase_mem then panel_n_variants[k] else 0,
                docker = glimpse2_docker
        }
    }

    call GLIMPSE2Ligate {
        input:
            phased_vcfs = ChunkedGLIMPSE2Phase.phased_vcf,
            phased_vcf_idxs = ChunkedGLIMPSE2Phase.phased_vcf_idx,
            output_prefix = output_prefix + ".glimpse2.bubble",
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
                output_prefix = output_prefix + ".glimpse2.bubble"
        }

        call RemapSampleNames as RemapPoppedPosteriors {
            input:
                vcf = ConcatPopAndMarginalizeCollisions.concatenated_vcf,
                vcf_idx = ConcatPopAndMarginalizeCollisions.concatenated_vcf_idx,
                remap_file = select_first([remap_sample_names_file]),
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
    Int? boot_disk_gb       # NOTE: no longer consumed by the tasks in this file; see the
                            # bootDiskSizeGb note in GLIMPSE2Phase. Retained so the struct
                            # stays compatible with the other WDLs that share it.
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

    # Optional, and absent from the panel JSONs shipped so far -- when it is missing every
    # shard falls back to static sizing, so existing inputs keep working unchanged.
    #
    # One entry per shard, parallel to input_regions: the number of panel variants in that
    # shard's *input* region (the buffered region, column 3 of chunks.tsv). This is exactly
    # GLIMPSE2's reported L. Note it is NOT column 7 of chunks.tsv, which counts a different
    # region and disagrees badly -- 403,101 vs an actual L of 965,039 for chr20 shard 4.
    #
    # Counting the sites-only panel BCF over each input region reproduced GLIMPSE2's L exactly
    # on all three shards checked (chr22 s0 = 657,784; chr20 s4 = 965,039; chr7 s12 =
    # 1,346,888), so this can be generated offline from resources already staged, with no VMs.
    Array[Int]? n_variants
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
        String? extra_phase_args = "--impute-reference-only-variants --keep-monomorphic-ref-sites --Kpbwt 1000 --main 10 --burnin 5 --err-imp 1E-3"

        String docker

        # Pinned rather than $(nproc). The forward matrix is allocated per thread, so with
        # --thread $(nproc) the memory requirement is a function of cpu_cores, while the N1
        # ratio limit makes cpu_cores a function of the memory requirement. Solving that
        # coupling for a shard of L variants (mem ~= 5.5 + 4*T*L/1e6, cpu >= mem/6.5) gives
        #   cpu * (6.5 - 4*L/1e6) >= 5.5
        # which has no solution above L ~= 1.6M and is already marginal at the genome-wide
        # max of L = 1.35M. Pinning T decouples the two and makes the sizing below well-posed.
        # It is also why cpu_cores below may exceed phase_threads: the extra cores buy ratio
        # headroom, not parallelism.
        Int phase_threads = 4

        # Number of variants (GLIMPSE2's L) in this shard's *input* region, if known.
        # 0 means unknown and selects the static fallback sizing. See the workflow-level
        # n_variants note on ChunkedPanelChromosome.
        Int n_variants = 0

        RuntimeAttr? runtime_attr_override
    }

    # Memory model, from phase/src/models/imputation_hmm.cpp:
    #   Alpha.resize(polymorphic_sites.size() * modK)   // float, modK = 1000 at --Kpbwt 1000
    # => 4 KB per site per thread, i.e. 4 GB per 1e6 sites per thread, plus ~5.5 GB fixed.
    # The +8 rather than +5.5 is deliberate headroom: the fit comes from a handful of points
    # near a noisy threshold and no peak RSS has actually been measured.
    Int auto_mem_gb   = 8 + ceil(4.0 * phase_threads * n_variants / 1000000.0)
    Int final_mem_gb  = if n_variants > 0 then auto_mem_gb else 40
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

        # The memory request is a function of phase_threads (see the sizing block above), so the
        # thread count has to come from that one place. Existing input CSVs carry --thread inside
        # extra_phase_args; passing both would hand GLIMPSE2 a duplicate option (boost
        # program_options rejects those), and honouring the CSV's value instead would silently
        # decouple threads from the memory it was sized against. So: strip any --thread from
        # extra_phase_args, warn loudly, and always inject phase_threads.
        EXTRA_PHASE_ARGS="~{extra_phase_args}"
        if echo "$EXTRA_PHASE_ARGS" | grep -qE '(^|[[:space:]])--thread([[:space:]]|=)'; then
            echo "WARNING: --thread found in extra_phase_args; overriding with phase_threads=~{phase_threads}." >&2
            echo "WARNING: set the phase_threads input instead -- memory is sized from it." >&2
            EXTRA_PHASE_ARGS=$(echo "$EXTRA_PHASE_ARGS" | sed -E 's/(^|[[:space:]])--thread([[:space:]]+|=)[0-9]+/ /g')
        fi

        cmd="/bin/GLIMPSE2_phase \
                --input-gl ~{input_vcf} \
                -R ~{panel_split_chunk_bin} \
                --thread ~{phase_threads} \
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
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        zones:                  ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]
        checkpointFile:         "checkpoint.bin"
    }
}

task GLIMPSE2Ligate {
    input {
        Array[File] phased_vcfs
        Array[File] phased_vcf_idxs
        String output_prefix

        String docker

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 2 * ceil(size(phased_vcfs, "GB")) + 10

    command <<<
        set -euox pipefail

        # Threads pinned rather than $(nproc): bcf_sr_set_threads allocates per-thread
        # decompression buffers, so scaling threads with cpu_cores would grow the very
        # memory this task is being widened to accommodate.
        /bin/GLIMPSE2_ligate --input ~{write_lines(phased_vcfs)} --output ~{output_prefix}.bcf --thread 2

        # the index generated by ligate appears to be corrupt for both bcf and vcf.gz output (possibly due to https://github.com/samtools/htslib/issues/1740), so we regenerate with bcftools
        bcftools index -f ~{output_prefix}.bcf
    >>>

    output {
        File ligated_vcf = "~{output_prefix}.bcf"
        File ligated_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
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
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        zones:                  ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]
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
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
        zones:                  ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]
    }
}

task RemapSampleNames {
    input {
        File vcf
        File vcf_idx
        File remap_file
        String output_prefix

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
        zones:                  ["us-central1-a", "us-central1-b", "us-central1-c", "us-central1-f"]
    }
}
