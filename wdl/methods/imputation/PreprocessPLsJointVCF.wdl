version 1.0

import "GLIMPSE2BatchedCaseShardedSingleBatch.wdl" as GLIMPSE2BatchedCaseShardedSingleBatch
import "../ConcatVcfs.wdl" as ConcatVcfs

workflow PreprocessPLsJointVCF {
    input {
        File input_joint_vcf
        File input_joint_vcf_idx
        File sample_names_file

        String chromosome
        File? chunked_panel_json

        String output_prefix

        # inputs for PreprocessPLs
        File preprocess_panel_bubble_split_sites_only_vcf       # can be subset of panel, e.g., simple bubble alleles only
        File preprocess_panel_bubble_split_sites_only_vcf_idx
        File? extract_bubble_likelihoods_script
        File? extract_bubble_likelihoods_cargo_toml
        File? extract_bubble_likelihoods_binary
        String? extract_bubble_likelihoods_extra_args
    }

    if (defined(chunked_panel_json)) {
        Map[String, ChunkedPanelChromosome] chunked_panel = read_json(select_first([chunked_panel_json]))
        Array[String] output_regions = chunked_panel[chromosome].output_regions
    }

    Array[String] regions = select_first([output_regions, [chromosome]])

    scatter (k in range(length(regions))) {
        call GLIMPSE2BatchedCaseShardedSingleBatch.PreprocessPLs as ChunkedPreprocessPLsJoint {
            input:
                input_vcf = input_joint_vcf,
                input_vcf_idx = input_joint_vcf_idx,
                mode = "joint",
                panel_bubble_split_sites_only_vcf = preprocess_panel_bubble_split_sites_only_vcf,
                panel_bubble_split_sites_only_vcf_idx = preprocess_panel_bubble_split_sites_only_vcf_idx,
                output_region = regions[k],
                sample_names = read_lines(sample_names_file),
                output_prefix = output_prefix + ".shard-" + k + ".preprocessedPLs",
                extract_bubble_likelihoods_script = extract_bubble_likelihoods_script,
                cargo_toml = extract_bubble_likelihoods_cargo_toml,
                extract_bubble_likelihoods_binary = extract_bubble_likelihoods_binary,
                extra_args = extract_bubble_likelihoods_extra_args
        }
    }

    call ConcatVcfs.ConcatVcfs as ConcatPreprocessPLsJoint { input:
        vcfs = ChunkedPreprocessPLsJoint.preprocessed_pls_vcf,
        vcf_idxs = ChunkedPreprocessPLsJoint.preprocessed_pls_vcf_idx,
        output_prefix = output_prefix + ".preprocessedPLs",
        do_bcf = true,
        do_sort = false,
        extra_args = "--threads $(nproc) --naive",
        regions = [],
        do_sort_shard = false,
        extra_args_shard = ""
    }

    output {
        File preprocessed_pls_vcf = ConcatPreprocessPLsJoint.concatenated_vcf
        File preprocessed_pls_vcf_idx = ConcatPreprocessPLsJoint.concatenated_vcf_idx
    }
}

struct ChunkedPanelChromosome {
    String chunks_tsv
    Array[String] input_regions
    Array[String] output_regions
    Array[String] panel_split_chunk_bins
}
