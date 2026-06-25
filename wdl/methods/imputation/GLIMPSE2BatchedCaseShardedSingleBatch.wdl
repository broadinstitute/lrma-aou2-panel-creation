version 1.0

import "../ConcatVcfs.wdl" as ConcatVcfs
import "../MultilevelHierarchicallyPasteVcfsStreaming.wdl" as MultilevelHierarchicallyPasteVcfsStreaming

workflow GLIMPSE2BatchedCaseShardedSingleBatch {
    input {
        # joint VCF, single-sample gVCF, and preprocessed joint VCF inputs are mutually exclusive
        File? input_joint_vcf
        File? input_joint_vcf_idx

        File? input_gvcfs_fofn
        File? input_gvcf_idxs_fofn
        File paste_vcfs_binary
        # TODO expose batch_sizes

        File? input_preprocessed_joint_vcf
        File? input_preprocessed_joint_vcf_idx

        File? sample_names_file          # in gVCF mode, order of sample names must match that of gVCFs
        File? remap_sample_names_file   # TSV with old_name new_name mappings

        String chromosome
        File genetic_maps_tsv
        File panel_bubble_split_vcf          # "split" here means "split to biallelic"; "split" just below means "chunked"
        File panel_bubble_split_vcf_idx
        File panel_bubble_split_sites_only_vcf
        File panel_bubble_split_sites_only_vcf_idx

        String extra_chunk_args = "--thread $(nproc) --window-mb 5 --buffer-mb 0.5 --sequential"
        String extra_split_args = "--keep-monomorphic-ref-sites"
        String extra_phase_args = "--impute-reference-only-variants --keep-monomorphic-ref-sites"
        String output_prefix

        # bypass for panel resources
        Array[String]? input_regions_bypass
        Array[String]? output_regions_bypass
        Array[File]? panel_split_chunk_bins_bypass

        # inputs for PreprocessPLs
        File? preprocess_panel_bubble_split_sites_only_vcf       # can be subset of panel, e.g., simple bubble alleles only
        File? preprocess_panel_bubble_split_sites_only_vcf_idx
        File? extract_bubble_likelihoods_script
        File? extract_bubble_likelihoods_cargo_toml
        File? extract_bubble_likelihoods_binary
        String? extract_bubble_likelihoods_extra_args

        # inputs for PopAndMarginalizeCollisions
        File panel_id_split_vcf_gz
        File panel_id_split_vcf_gz_tbi
        File? pop_glimpse2_script      # modified version of convert-to-biallelic.py
        File? pop_glimpse2_cargo_toml
        File? pop_glimpse2_binary

        String glimpse2_docker = "us.gcr.io/broad-gotc-prod/imputation-glimpse2:1.0.0-2cee597-1778869818"    # enables checkpointing, but note this contains bcftools/htslib 1.16!
    }

    Map[String, String] genetic_maps_dict = read_map(genetic_maps_tsv)

    if (!defined(input_regions_bypass) && !defined(output_regions_bypass) && !defined(panel_split_chunk_bins_bypass)) {
        call GLIMPSE2Chunk {
            input:
                vcf = panel_bubble_split_vcf,
                vcf_idx = panel_bubble_split_vcf_idx,
                region = chromosome,
                genetic_map = genetic_maps_dict[chromosome],
                output_prefix = output_prefix,
                extra_chunk_args = extra_chunk_args,
                docker = glimpse2_docker
        }

        Array[String] input_regions_generated = read_lines(GLIMPSE2Chunk.input_regions)
        Array[String] output_regions_generated = read_lines(GLIMPSE2Chunk.output_regions)

        scatter (k in range(length(output_regions_generated))) {
            call GLIMPSE2SplitReference as ChunkedGLIMPSE2SplitReference {
                input:
                    panel_bubble_split_vcf = panel_bubble_split_vcf,
                    panel_bubble_split_vcf_idx = panel_bubble_split_vcf_idx,
                    input_region = input_regions_generated[k],
                    output_region = output_regions_generated[k],
                    genetic_map = genetic_maps_dict[chromosome],
                    output_prefix = output_prefix + ".shard-" + k + ".split",
                    extra_split_args = extra_split_args,
                    docker = glimpse2_docker
            }
        }
    }

    Array[String] input_regions_ = select_first([input_regions_bypass, input_regions_generated])
    Array[String] output_regions_ = select_first([output_regions_bypass, output_regions_generated])
    Array[File] panel_split_chunk_bins_ = select_first([panel_split_chunk_bins_bypass, ChunkedGLIMPSE2SplitReference.panel_split_chunk_bin])

    # joint ###################################################################
    if (defined(input_joint_vcf) && defined(input_joint_vcf_idx) && !defined(input_gvcfs_fofn) && !defined(input_gvcf_idxs_fofn) && !defined(input_preprocessed_joint_vcf) && !defined(input_preprocessed_joint_vcf_idx) && 
        defined(preprocess_panel_bubble_split_sites_only_vcf) && defined(preprocess_panel_bubble_split_sites_only_vcf_idx) && defined(sample_names_file)) {
        Array[String] sample_names_joint = read_lines(select_first([sample_names_file]))
        
        scatter (k in range(length(output_regions_))) {
            call PreprocessPLs as ChunkedPreprocessPLsJoint {
                input:
                    input_vcf = select_first([input_joint_vcf]),
                    input_vcf_idx = select_first([input_joint_vcf_idx]),
                    mode = "joint",
                    panel_bubble_split_sites_only_vcf = select_first([preprocess_panel_bubble_split_sites_only_vcf]),
                    panel_bubble_split_sites_only_vcf_idx = select_first([preprocess_panel_bubble_split_sites_only_vcf_idx]),
                    output_region = output_regions_[k],
                    sample_names = sample_names_joint,
                    output_prefix = output_prefix + ".shard-" + k + ".preprocessedPLs",
                    extract_bubble_likelihoods_script = extract_bubble_likelihoods_script,
                    cargo_toml = extract_bubble_likelihoods_cargo_toml,
                    extract_bubble_likelihoods_binary = extract_bubble_likelihoods_binary,
                    extra_args = extract_bubble_likelihoods_extra_args
            }
        }
    }
    ###########################################################################

    # gVCF ####################################################################
    if (!defined(input_joint_vcf) && !defined(input_joint_vcf_idx) && defined(input_gvcfs_fofn) && defined(input_gvcf_idxs_fofn) && !defined(input_preprocessed_joint_vcf) && !defined(input_preprocessed_joint_vcf_idx) && 
        defined(preprocess_panel_bubble_split_sites_only_vcf) && defined(preprocess_panel_bubble_split_sites_only_vcf_idx) && defined(sample_names_file)) {
        Array[File] input_gvcfs = read_lines(select_first([input_gvcfs_fofn]))
        Array[File] input_gvcf_idxs = read_lines(select_first([input_gvcf_idxs_fofn]))
        Array[String] sample_names_gvcf = read_lines(select_first([sample_names_file]))
    
        # for each input, localize entire genome gVCF but genotype only requested chromosome (i.e., some redundant localization)
        scatter (j in range(length(select_first([input_gvcfs])))) {
            call PreprocessPLs as PreprocessPLsGVCF {
                input:
                    input_vcf = input_gvcfs[j],
                    input_vcf_idx = input_gvcf_idxs[j],
                    mode = "gvcf",
                    panel_bubble_split_sites_only_vcf = select_first([preprocess_panel_bubble_split_sites_only_vcf]),
                    panel_bubble_split_sites_only_vcf_idx = select_first([preprocess_panel_bubble_split_sites_only_vcf_idx]),
                    output_region = chromosome,
                    sample_names = [sample_names_gvcf[j]],
                    output_prefix = output_prefix + ".sample-" + j + "." + sample_names_gvcf[j] + ".preprocessedPLs",
                    extract_bubble_likelihoods_script = extract_bubble_likelihoods_script,
                    cargo_toml = extract_bubble_likelihoods_cargo_toml,
                    extract_bubble_likelihoods_binary = extract_bubble_likelihoods_binary,
                    extra_args = extract_bubble_likelihoods_extra_args
            }
        }

        # two-level localized hierarchical merge over entire chromosome
        call MultilevelHierarchicallyPasteVcfsStreaming.HierarchicallyMergeVcfs as PastePreprocessPLsGVCFs {
            input:
                vcfs_array = PreprocessPLsGVCF.preprocessed_pls_vcf,
                vcf_idxs_array = PreprocessPLsGVCF.preprocessed_pls_vcf_idx,
                regions = [chromosome],
                batch_sizes = [50, 50],
                do_localization = [true, true],
                timeouts_min = [0, 0],
                output_prefix = output_prefix + ".preprocessedPLs",
                paste_vcfs_binary = paste_vcfs_binary,
                extra_merge_args = "--threads $(nproc) --format GT,PL",
                extra_concat_args = "--threads $(nproc) --naive"
        }

        scatter (k in range(length(output_regions_))) {
            File preprocessed_pls_vcf = PastePreprocessPLsGVCFs.merged_vcf
            File preprocessed_pls_vcf_idx = PastePreprocessPLsGVCFs.merged_vcf_idx
        }

        Array[File] gvcf_preprocessed_pls_vcfs = preprocessed_pls_vcf
        Array[File] gvcf_preprocessed_pls_vcf_idxs = preprocessed_pls_vcf_idx
    }
    ###########################################################################

    # preprocessed joint#######################################################
    if (!defined(input_joint_vcf) && !defined(input_joint_vcf_idx) && !defined(input_gvcfs_fofn) && !defined(input_gvcf_idxs_fofn) && defined(input_preprocessed_joint_vcf) && defined(input_preprocessed_joint_vcf_idx) && 
        !defined(preprocess_panel_bubble_split_sites_only_vcf) && !defined(preprocess_panel_bubble_split_sites_only_vcf_idx) && !defined(sample_names_file)) {
        scatter (k in range(length(output_regions_))) {
            File joint_preprocessed_pls_vcf = select_first([input_preprocessed_joint_vcf])
            File joint_preprocessed_pls_vcf_idx = select_first([input_preprocessed_joint_vcf_idx])
        }

        Array[File] joint_preprocessed_pls_vcfs = joint_preprocessed_pls_vcf
        Array[File] joint_preprocessed_pls_vcf_idxs = joint_preprocessed_pls_vcf_idx
    }
    ###########################################################################

    # in joint mode the preprocessed VCFs are chunked, but in gVCF and joint preprocessed VCF modes we redundantly pass the same preprocessed VCF to all chunks
    Array[File] preprocessed_pls_vcfs = select_first([ChunkedPreprocessPLsJoint.preprocessed_pls_vcf, gvcf_preprocessed_pls_vcfs, joint_preprocessed_pls_vcfs])
    Array[File] preprocessed_pls_vcf_idxs = select_first([ChunkedPreprocessPLsJoint.preprocessed_pls_vcf_idx, gvcf_preprocessed_pls_vcf_idxs, joint_preprocessed_pls_vcf_idxs])

    scatter (k in range(length(output_regions_))) {
        call GLIMPSE2Phase as ChunkedGLIMPSE2Phase {
            input:
                input_vcf = preprocessed_pls_vcfs[k],
                input_vcf_idx = preprocessed_pls_vcf_idxs[k],
                panel_split_chunk_bin = panel_split_chunk_bins_[k],
                input_region = input_regions_[k],
                output_region = output_regions_[k],
                genetic_map = genetic_maps_dict[chromosome],
                output_prefix = output_prefix + ".shard-" + k + ".glimpse2.phased",
                extra_phase_args = extra_phase_args,
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

    scatter (k in range(length(output_regions_))) {
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
            region = output_regions_[k],
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
        Array[String] input_regions = input_regions_
        Array[String] output_regions = output_regions_
        Array[File] panel_split_chunk_bins = panel_split_chunk_bins_
        
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
    Int? boot_disk_gb
    Boolean? use_ssd
    Int? preemptible_tries
    Int? max_retries
    String? docker
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
            -O chunks.tsv

        # cut chunks + buffers
        cut -f 3 chunks.tsv > ~{output_prefix}.input-regions.tsv
        cut -f 4 chunks.tsv > ~{output_prefix}.output-regions.tsv
    >>>

    output {
        File input_regions = "~{output_prefix}.input-regions.tsv"
        File output_regions = "~{output_prefix}.output-regions.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             7,
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


task PreprocessPLs {
    input {
        File input_vcf
        File input_vcf_idx
        String mode     # joint or gvcf
        File panel_bubble_split_sites_only_vcf
        File panel_bubble_split_sites_only_vcf_idx
        String? output_region
        Array[String] sample_names
        String output_prefix

        File? extract_bubble_likelihoods_script
        File? cargo_toml
        File? extract_bubble_likelihoods_binary
        String? extra_args = "--window 15000 --cap-pl 30 --scale-pl 5.0 --threads $(nproc)"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 10 + 2 * ceil(size([input_vcf, panel_bubble_split_sites_only_vcf], "GB"))

    File sample_names_list = write_lines(sample_names)

    command <<<
        set -euxo pipefail

        if [ -n "~{extract_bubble_likelihoods_binary}" ]; then
            EXTRACT_BIN="~{extract_bubble_likelihoods_binary}"
            chmod +x $EXTRACT_BIN
        else
            mkdir -p extract-bubble-PLs/src
            cp ~{extract_bubble_likelihoods_script} extract-bubble-PLs/src/main.rs
            cp ~{cargo_toml} extract-bubble-PLs
            cd extract-bubble-PLs
            cargo build --release
            cd ..
            EXTRACT_BIN="./extract-bubble-PLs/target/release/extract-bubble-PLs"
        fi

        $EXTRACT_BIN ~{mode} \
            ~{panel_bubble_split_sites_only_vcf}##idx##~{panel_bubble_split_sites_only_vcf_idx} \
            ~{input_vcf}##idx##~{input_vcf_idx} \
            ~{output_prefix}.bcf \
            ~{"--region " + output_region} \
            --samples ~{sample_names_list} \
            ~{extra_args}
            
        bcftools index ~{output_prefix}.bcf

        echo "Number of bubble alleles extracted..."
        bcftools index -n ~{output_prefix}.bcf
    >>>

    output {
        File preprocessed_pls_vcf = "~{output_prefix}.bcf"
        File preprocessed_pls_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
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
        String? extra_phase_args

        String docker

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 50 + 3 * ceil(size([input_vcf, panel_split_chunk_bin], "GB"))

    command {
        set -euxo pipefail

        cmd="/bin/GLIMPSE2_phase \
                --input-gl ~{input_vcf} \
                -R ~{panel_split_chunk_bin} \
                --thread $(nproc) \
                ~{extra_phase_args} \
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
    }

    output {
        File phased_vcf = "~{output_prefix}.bcf"
        File phased_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  9,
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

        /bin/GLIMPSE2_ligate --input ~{write_lines(phased_vcfs)} --output ~{output_prefix}.bcf --thread $(nproc)

        # the index generated by ligate appears to be corrupt for both bcf and vcf.gz output (possibly due to https://github.com/samtools/htslib/issues/1740), so we regenerate with bcftools
        bcftools index -f ~{output_prefix}.bcf
    >>>

    output {
        File ligated_vcf = "~{output_prefix}.bcf"
        File ligated_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
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
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
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
        boot_disk_gb:       10,
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
    }
}
