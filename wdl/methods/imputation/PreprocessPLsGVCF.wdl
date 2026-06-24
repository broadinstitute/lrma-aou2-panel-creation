version 1.0

import "GLIMPSE2BatchedCaseShardedSingleBatch.wdl" as GLIMPSE2BatchedCaseShardedSingleBatch
import "../MultilevelHierarchicallyPasteVcfsStreaming.wdl" as MultilevelHierarchicallyPasteVcfsStreaming

workflow PreprocessPLsGVCF {
    input {
        File? input_gvcfs_fofn
        File? input_gvcf_idxs_fofn
        File? sample_names_file          # order of sample names must match that of gVCFs
        
        Array[File]? input_gvcfs
        Array[File]? input_gvcf_idxs
        Array[String]? sample_names

        String output_prefix

        # inputs for PreprocessPLs
        File preprocess_panel_bubble_split_sites_only_vcf       # can be subset of panel, e.g., simple bubble alleles only
        File preprocess_panel_bubble_split_sites_only_vcf_idx
        File? extract_bubble_likelihoods_script
        File? extract_bubble_likelihoods_cargo_toml
        File? extract_bubble_likelihoods_binary
        String? extract_bubble_likelihoods_extra_args

        File paste_vcfs_binary
        Array[String] paste_regions
    }

    Array[File] input_gvcfs_ = select_first([read_lines(select_first([input_gvcfs_fofn])), input_gvcfs])
    Array[File] input_gvcf_idxs_ = select_first([read_lines(select_first([input_gvcf_idxs_fofn])), input_gvcf_idxs])
    Array[String] sample_names_ = select_first([read_lines(select_first([sample_names_file])), sample_names])

    scatter (j in range(length(select_first([input_gvcfs_])))) {
        call GLIMPSE2BatchedCaseShardedSingleBatch.PreprocessPLs as PreprocessPLsGVCF {
            input:
                input_vcf = input_gvcfs_[j],
                input_vcf_idx = input_gvcf_idxs_[j],
                mode = "gvcf",
                panel_bubble_split_sites_only_vcf = preprocess_panel_bubble_split_sites_only_vcf,
                panel_bubble_split_sites_only_vcf_idx = preprocess_panel_bubble_split_sites_only_vcf_idx,
                sample_names = [sample_names_[j]],
                output_prefix = output_prefix + ".sample-" + j + "." + sample_names_[j] + ".preprocessedPLs",
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
            regions = paste_regions,
            batch_sizes = [50, 50],
            do_localization = [true, true],
            timeouts_min = [0, 0],
            output_prefix = output_prefix + ".preprocessedPLs",
            paste_vcfs_binary = paste_vcfs_binary,
            extra_merge_args = "--threads $(nproc) --format GT,PL",
            extra_concat_args = "--threads $(nproc) --naive"
    }

    output {
        File preprocessed_pls_vcf = PastePreprocessPLsGVCFs.merged_vcf
        File preprocessed_pls_vcf_idx = PastePreprocessPLsGVCFs.merged_vcf_idx
    }
}
