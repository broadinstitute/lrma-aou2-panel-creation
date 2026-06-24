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
        Array[String]? entity_ids
        File? sample_names_map_file           # TSV map of entity_id (research_id) to id2 for AoU DRAGEN gVCFs;
                                              # Terra struggles with id2 as they are parsed as mixed strings/numbers

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

    if (defined(input_gvcfs_fofn)) {
        Array[File] parsed_gvcfs = read_lines(select_first([input_gvcfs_fofn]))
    }
    Array[File] input_gvcfs_ = select_first([input_gvcfs, parsed_gvcfs])

    if (defined(input_gvcf_idxs_fofn)) {
        Array[File] parsed_gvcf_idxs = read_lines(select_first([input_gvcf_idxs_fofn]))
    }
    Array[File] input_gvcf_idxs_ = select_first([input_gvcf_idxs, parsed_gvcf_idxs])

    if (defined(sample_names_file)) {
        Array[String] parsed_sample_names = read_lines(select_first([sample_names_file]))
    }
    
    # Replaced map scatter with a bash task call
    if (defined(entity_ids) && defined(sample_names_map_file)) {
        call MapSampleNames {
            input:
                entity_ids = select_first([entity_ids]),
                sample_names_map_file = select_first([sample_names_map_file])
        }
    }
    
    Array[String] sample_names_ = select_first([MapSampleNames.mapped_sample_names, parsed_sample_names])

    scatter (j in range(length(input_gvcfs_))) {
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

task MapSampleNames {
    input {
        Array[String] entity_ids
        File sample_names_map_file
    }

    command <<<
        set -euo pipefail

        # Use awk to load the TSV map into memory, then translate the entity IDs array in order
        awk 'BEGIN {FS="\t"; OFS="\t"} 
        NR==FNR { 
            # First pass: read the map file into an array
            map[$1] = $2; 
            next 
        } 
        { 
            # Second pass: read the entity_ids file
            if ($1 in map) { 
                print map[$1] 
            } else { 
                print "Error: ID " $1 " not found in map file" > "/dev/stderr"
                exit 1 
            } 
        }' ~{sample_names_map_file} ~{write_lines(entity_ids)} > mapped_names.txt
    >>>

    output {
        Array[String] mapped_sample_names = read_lines("mapped_names.txt")
    }

    runtime {
        docker: "ubuntu:22.04"
        cpu: 1
        memory: "4 GB"
        disks: "local-disk 10 HDD"
        preemptible: 3
    }
}
