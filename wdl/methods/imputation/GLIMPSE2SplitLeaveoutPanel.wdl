version 1.0

import "GLIMPSE2ChunkAndSplitPanel.wdl" as ChunkAndSplit

# Build the .bin reference files for a leaveout panel, reusing the full panel's shard boundaries.
#
# The leaveout evaluation imputes held-out samples against a panel they are not in, then scores
# against the full panel, which still carries their assembly-based genotypes. The leaveout BCF
# exists; what does not is its .bin files. GLIMPSE2_phase consumes .bin, not BCF, and the
# chunked_panel that ships with the panel was split from the FULL panel -- using it would put
# the held-out samples' own haplotypes back into the reference and make the evaluation circular.
#
# WHY THIS IS NOT GLIMPSE2ChunkAndSplitPanel. That workflow runs GLIMPSE2_chunk first and splits
# against the boundaries it produces. Re-chunking the leaveout panel would produce DIFFERENT
# windows -- the panel has ~490 fewer samples, so some sites are no longer variant and the
# chunker's sequential window walk lands elsewhere. Different windows mean the memory numbers
# are not comparable to the full-panel run and the per-shard concordance no longer lines up
# shard-for-shard with Samuel's. So this reuses the full panel's chunks.tsv verbatim. Only the
# .bin contents and the variant counts change.
#
# It takes chunks.tsv rather than the full panel's chunked_panel JSON on purpose. Reading that
# JSON would mean coercing it to ChunkedPanelChromosome, whose n_variants field a panel chunked
# before that field existed does not have -- the exact gap scripts/add-panel-variant-counts.py
# was written to backfill. chunks.tsv has no such dependency, and columns 3 and 4 are where
# GLIMPSE2ChunkAndSplitPanel gets its own regions, so the boundaries are identical by
# construction.
#
# n_variants IS recounted, against the leaveout panel's own sites. Dropping samples drops sites,
# so the full panel's counts would over-size the phase memory request. The counts are emitted
# straight into the JSON, so no add-panel-variant-counts.py backfill is needed downstream.
#
# The output JSON is a drop-in for GLIMPSE2FromPreprocessedPLsJoint.chunked_panel_json.

workflow GLIMPSE2SplitLeaveoutPanel {
    input {
        Array[String] chromosomes
        File genetic_maps_tsv

        # The FULL panel's chunks.tsv, per chromosome, in the same order as `chromosomes`.
        # Reused verbatim; nothing here re-chunks.
        Array[File] full_panel_chunks_tsvs

        # per chromosome, in the same order as `chromosomes`
        Array[File] leaveout_panel_bubble_split_vcfs
        Array[File] leaveout_panel_bubble_split_vcf_idxs
        Array[File] leaveout_panel_bubble_split_sites_only_vcfs
        Array[File] leaveout_panel_bubble_split_sites_only_vcf_idxs

        # Sample count the leaveout panel is expected to have. Checked before anything expensive
        # runs, because pointing this at the full panel by mistake produces a complete, valid,
        # silently circular set of .bin files. Unset disables the check.
        Int? expected_n_samples

        # Must match what the full panel was split with, or the .bin files are not interchangeable.
        String extra_split_args = "--keep-monomorphic-ref-sites"
        String output_prefix

        String glimpse2_docker = "us.gcr.io/broad-gotc-prod/imputation-glimpse2:1.0.0-2cee597-1778869818"
    }

    Map[String, String] genetic_maps_dict = read_map(genetic_maps_tsv)

    scatter (i in range(length(chromosomes))) {
        String chromosome = chromosomes[i]

        call ParseChunksTsv {
            input:
                chunks_tsv = full_panel_chunks_tsvs[i],
                chromosome = chromosome,
                docker = glimpse2_docker
        }

        # Verbatim from the full panel's chunks.tsv. Not recomputed anywhere in this workflow.
        Array[String] input_regions = read_lines(ParseChunksTsv.input_regions)
        Array[String] output_regions = read_lines(ParseChunksTsv.output_regions)

        call AssertLeaveoutPanel {
            input:
                leaveout_panel_vcf = leaveout_panel_bubble_split_vcfs[i],
                leaveout_panel_vcf_idx = leaveout_panel_bubble_split_vcf_idxs[i],
                expected_n_samples = expected_n_samples,
                n_shards = length(output_regions),
                chromosome = chromosome,
                docker = glimpse2_docker
        }

        call ChunkAndSplit.CountPanelVariantsPerShard as CountLeaveoutPanelVariantsPerShard {
            input:
                panel_bubble_split_sites_only_vcf = leaveout_panel_bubble_split_sites_only_vcfs[i],
                panel_bubble_split_sites_only_vcf_idx = leaveout_panel_bubble_split_sites_only_vcf_idxs[i],
                input_regions = input_regions,
                output_prefix = output_prefix + "." + chromosome,
                docker = glimpse2_docker
        }

        scatter (k in range(length(output_regions))) {
            call ChunkAndSplit.GLIMPSE2SplitReference as ChunkedGLIMPSE2SplitReference {
                input:
                    panel_bubble_split_vcf = leaveout_panel_bubble_split_vcfs[i],
                    panel_bubble_split_vcf_idx = leaveout_panel_bubble_split_vcf_idxs[i],
                    input_region = input_regions[k],
                    output_region = output_regions[k],
                    genetic_map = genetic_maps_dict[chromosome],
                    output_prefix = output_prefix + "." + chromosome + ".shard-" + k + ".split",
                    extra_split_args = extra_split_args,
                    docker = glimpse2_docker
            }
        }

        # chunks_tsv is the full panel's file, carried through unchanged: same boundaries, so it
        # still describes these shards, and downstream reads it for provenance rather than for
        # control flow.
        ChunkedPanelChromosome leaveout_chunked_panel_chromosome = object {
            chunks_tsv: full_panel_chunks_tsvs[i],
            input_regions: input_regions,
            output_regions: output_regions,
            panel_split_chunk_bins: ChunkedGLIMPSE2SplitReference.panel_split_chunk_bin,
            n_variants: CountLeaveoutPanelVariantsPerShard.n_variants
        }
        Pair[String, ChunkedPanelChromosome] leaveout_chunked_panel_chromosome_pair = (chromosome, leaveout_chunked_panel_chromosome)
    }

    call ChunkAndSplit.CoercePairsToMap {
        input:
            pair_array = leaveout_chunked_panel_chromosome_pair,
            output_prefix = output_prefix + ".leaveout"
    }

    output {
        File leaveout_chunked_panel_json = CoercePairsToMap.out_map_json
        Map[String, ChunkedPanelChromosome] leaveout_chunked_panel = read_json(CoercePairsToMap.out_map_json)
        Array[File] sample_count_reports = AssertLeaveoutPanel.report
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

# Columns 3 and 4 of chunks.tsv, exactly as GLIMPSE2ChunkAndSplitPanel.GLIMPSE2Chunk cuts them.
# A task rather than transpose(read_tsv(...)) because transpose requires every row to have the
# same width and would abort on a ragged or trailing line, where cut does not.
task ParseChunksTsv {
    input {
        File chunks_tsv
        String chromosome

        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        # Drop blank lines so a trailing newline cannot become an empty region that
        # GLIMPSE2_split_reference would reject only after the scatter has fanned out.
        awk -F'\t' 'NF >= 4 && $3 != "" && $4 != ""' ~{chunks_tsv} > chunks.clean.tsv

        n=$(wc -l < chunks.clean.tsv)
        if [ "${n}" -eq 0 ]; then
            echo "ERROR: no usable rows in ~{chunks_tsv}; expected GLIMPSE2_chunk output with >=4 columns" >&2
            exit 1
        fi

        cut -f 3 chunks.clean.tsv > ~{chromosome}.input-regions.tsv
        cut -f 4 chunks.clean.tsv > ~{chromosome}.output-regions.tsv

        # Every region must name the chromosome being split. A mismatch means the per-chromosome
        # arrays are misaligned -- chr22's panel against chr21's chunks -- which otherwise
        # surfaces as an empty or wrong .bin much later.
        #
        # awk into a file, then test -s. NOT `grep -qv`: its exit status is inconsistent across
        # implementations (ugrep 7.5.0 returns 1 from `grep -vE` even while printing
        # non-matching lines), and a guard that silently never fires is worse than no guard.
        awk -v c="~{chromosome}" 'index($0, c ":") != 1' ~{chromosome}.input-regions.tsv > wrong-chrom.txt
        if [ -s wrong-chrom.txt ]; then
            echo "ERROR: ~{chunks_tsv} has regions outside ~{chromosome}:" >&2
            cat wrong-chrom.txt >&2
            exit 1
        fi

        echo "shards for ~{chromosome}: ${n}"
    >>>

    output {
        File input_regions = "~{chromosome}.input-regions.tsv"
        File output_regions = "~{chromosome}.output-regions.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            10,
        boot_disk_gb:       10,
        use_ssd:            false,
        preemptible_tries:  0,
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

# Reads the header only. Cheap, non-preemptible, and it runs before the split scatter so a
# mis-wired panel costs one task rather than a chromosome of .bin files that look fine.
task AssertLeaveoutPanel {
    input {
        File leaveout_panel_vcf
        File leaveout_panel_vcf_idx
        Int? expected_n_samples
        Int n_shards
        String chromosome

        String docker
        RuntimeAttr? runtime_attr_override
    }

    # 0 means "not supplied", so the check is skipped. Resolved here rather than inline in the
    # command, where `if defined(x) then x else ''` would mix Int and String.
    Int expected = select_first([expected_n_samples, 0])

    command <<<
        set -euxo pipefail

        bcftools query -l ~{leaveout_panel_vcf} > samples.txt
        n=$(wc -l < samples.txt)

        {
            echo "chromosome:       ~{chromosome}"
            echo "leaveout samples: ${n}"
            echo "shards:           ~{n_shards}"
        } | tee ~{chromosome}.leaveout-panel-check.txt

        if [ "~{expected}" -gt 0 ] && [ "${n}" -ne "~{expected}" ]; then
            echo "ERROR: expected ~{expected} samples in ~{chromosome}, found ${n}." >&2
            echo "       A count matching the FULL panel means this IS the full panel, and" >&2
            echo "       splitting it would make the leaveout evaluation circular." >&2
            exit 1
        fi
    >>>

    output {
        File report = "~{chromosome}.leaveout-panel-check.txt"
    }

    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            20,
        boot_disk_gb:       10,
        use_ssd:            false,
        preemptible_tries:  0,
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
