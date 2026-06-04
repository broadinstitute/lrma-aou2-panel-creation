version 1.0

import "../ConcatVcfs.wdl" as ConcatVcfs

workflow GLIMPSE2BatchedCaseShardedSingleBatch {
    input {
        File input_vcf
        File input_vcf_idx
        Array[String] sample_names = []     # empty list to select all

        String chromosome
        File genetic_maps_tsv
        File panel_bubble_split_vcf          # "split" here means "split to biallelic"; "split" below means "chunked"
        File panel_bubble_split_vcf_idx

        String extra_chunk_args = "--thread $(nproc) --window-mb 5 --buffer-mb 0.5 --sequential"
        String extra_split_args = "--keep-monomorphic-ref-sites"
        String extra_phase_args = "--impute-reference-only-variants --keep-monomorphic-ref-sites"
        String output_prefix

        # inputs for PreprocessPLs
        File remap_simple_bubble_likelihoods_python_script
        File swap_alleles_python_script
        File? preprocess_regions_bed
        String? preprocess_view_extra_args
        String? remap_simple_bubble_likelihoods_extra_args

        # inputs for PopAndMarkCollisions
        File panel_id_split_vcf
        File panel_id_split_vcf_idx
        File pop_python_script                      # modified version of convert-to-biallelic.py

        String glimpse2_docker = "us.gcr.io/broad-gotc-prod/imputation-glimpse2:1.0.0-2cee597-1778869818"    # enables checkpointing, but note this contains bcftools/htslib 1.16!
    }

    Map[String, String] genetic_maps_dict = read_map(genetic_maps_tsv)

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

    Array[String] input_regions_ = read_lines(GLIMPSE2Chunk.input_regions)
    Array[String] output_regions_ = read_lines(GLIMPSE2Chunk.output_regions)

    scatter (k in range(length(output_regions_))) {
        call GLIMPSE2SplitReference as ChunkedGLIMPSE2SplitReference {
            input:
                panel_bubble_split_vcf = panel_bubble_split_vcf,
                panel_bubble_split_vcf_idx = panel_bubble_split_vcf_idx,
                input_region = input_regions_[k],
                output_region = output_regions_[k],
                genetic_map = genetic_maps_dict[chromosome],
                output_prefix = output_prefix + ".shard-" + k + ".split",
                extra_split_args = extra_split_args,
                docker = glimpse2_docker
        }

        call PreprocessPLs as ChunkedPreprocessPLs {
            input:
                input_vcf = input_vcf,
                input_vcf_idx = input_vcf_idx,
                panel_bubble_split_vcf = panel_bubble_split_vcf,
                panel_bubble_split_vcf_idx = panel_bubble_split_vcf_idx,
                output_region = output_regions_[k],
                sample_names = sample_names,
                output_prefix = output_prefix + ".shard-" + k + ".preprocessedPLs",
                remap_simple_bubble_likelihoods_python_script = remap_simple_bubble_likelihoods_python_script,
                swap_alleles_python_script = swap_alleles_python_script,
                preprocess_regions_bed = preprocess_regions_bed,
                preprocess_view_extra_args = preprocess_view_extra_args,
                remap_simple_bubble_likelihoods_extra_args = remap_simple_bubble_likelihoods_extra_args
        }

        call GLIMPSE2Phase as ChunkedGLIMPSE2Phase {
            input:
                input_vcf = ChunkedPreprocessPLs.preprocessed_pls_bcf,
                input_vcf_idx = ChunkedPreprocessPLs.preprocessed_pls_bcf_csi,
                panel_split_chunk_bin = ChunkedGLIMPSE2SplitReference.panel_split_chunk_bin,
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
            phased_bcfs = ChunkedGLIMPSE2Phase.phased_bcf,
            phased_bcf_csis = ChunkedGLIMPSE2Phase.phased_bcf_csi,
            output_prefix = output_prefix + ".glimpse2.bubble",
            docker = glimpse2_docker
    }

    scatter (k in range(length(output_regions_))) {
        call PopAndMarkCollisions { input:
            posteriors_vcf = GLIMPSE2Ligate.ligated_vcf,
            posteriors_vcf_idx = GLIMPSE2Ligate.ligated_vcf_idx,
            panel_bubble_split_vcf = panel_bubble_split_vcf,
            panel_bubble_split_vcf_idx = panel_bubble_split_vcf_idx,
            panel_id_split_vcf = panel_id_split_vcf,
            panel_id_split_vcf_idx = panel_id_split_vcf_idx,
            pop_python_script = pop_python_script,
            region = output_regions_[k],
            output_prefix = output_prefix + ".glimpse2.popped"
        }
    }
    
    call ConcatVcfs.ConcatVcfs as ConcatPopAndMarkCollisions { input:
        vcfs = PopAndMarkCollisions.popped_vcf,
        vcf_idxs = PopAndMarkCollisions.popped_vcf_idx,
        output_prefix = output_prefix + ".popped",
        do_sort = false,
        extra_args = "--threads $(nproc) --naive"
    }

    output {
        Array[String] input_regions = input_regions_
        Array[String] output_regions = output_regions_
        Array[File] panel_split_chunk_bins = ChunkedGLIMPSE2SplitReference.panel_split_chunk_bin
        File glimpse2_bubble_posteriors_vcf = GLIMPSE2Ligate.ligated_vcf
        File glimpse2_bubble_posteriors_vcf_idx = GLIMPSE2Ligate.ligated_vcf_idx
        File glimpse2_popped_posteriors_vcf = ConcatPopAndMarkCollisions.concatenated_vcf
        File glimpse2_popped_posteriors_vcf_idx = ConcatPopAndMarkCollisions.concatenated_vcf_idx
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

    Int disk_size_gb = 2 * ceil(size(vcf, "GB")) + 10

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

    Int disk_size_gb = 2 * ceil(size(panel_bubble_split_vcf, "GB")) + 10

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
        File panel_bubble_split_vcf
        File panel_bubble_split_vcf_idx
        String output_region
        Array[String] sample_names
        String output_prefix

        File remap_simple_bubble_likelihoods_python_script
        File swap_alleles_python_script
        File? preprocess_regions_bed
        String? preprocess_view_extra_args
        String? remap_simple_bubble_likelihoods_extra_args

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 2 * ceil(size([input_vcf, panel_bubble_split_vcf], "GB")) + 10

    File sample_names_list = write_lines(sample_names)

    command {
        set -euxo pipefail

       # TODO add gcloud to Docker
#        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        # TODO stream; or, once shards are fixed, prepare beforehand?
        bcftools view --no-version -G ~{panel_bubble_split_vcf}##idx##~{panel_bubble_split_vcf_idx} \
            --regions-overlap pos -r ~{output_region} -Ou \
            ~{"-T " + preprocess_regions_bed} \
            ~{preprocess_view_extra_args} | \
        bcftools norm -m+any -N \
            --write-index=tbi -Oz -o panel.subset.sites.vcf.gz

        echo "Number of sites to preprocess..."
        bcftools index -n panel.subset.sites.vcf.gz

        pypy -m pip install --no-input tqdm

        # TODO stream
        bcftools view ~{input_vcf}##idx##~{input_vcf_idx} \
            -r ~{output_region} \
            --regions-overlap pos \
            ~{if length(sample_names) > 0 then "-S " + sample_names_list else ""} \
            --threads 2 | \
        pypy ~{remap_simple_bubble_likelihoods_python_script} \
            --bubble panel.subset.sites.vcf.gz \
            ~{remap_simple_bubble_likelihoods_extra_args} | \
        bcftools +tag2tag -Ou -- --LPL-to-PL | \
        bcftools norm -m-any -Ou | \
        bcftools filter -i 'INFO/BMAP != "."' | \
        pypy ~{swap_alleles_python_script} | \
        bcftools annotate -x INFO/BMAP,INFO/BUBBLE,INFO/BPOS,INFO/BREF,INFO/BALT | \
        bcftools sort --write-index=csi -Ob -o ~{output_prefix}.bcf
    }

    output {
        File preprocessed_pls_bcf = "~{output_prefix}.bcf"
        File preprocessed_pls_bcf_csi = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             6,
        disk_gb:            disk_size_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsde-methods/slee/lrma-aou2-panel-creation-pypy:v1"
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

    Int disk_size_gb = 10 + 2 * ceil(size([input_vcf, panel_split_chunk_bin], "GB"))

    command {
        set -euxo pipefail

        cmd="/bin/GLIMPSE2_phase \
                --input-gl ~{input_vcf} \
                -R ~{panel_split_chunk_bin} \
                --thread $(nproc) \
                ~{extra_phase_args} \
                --output ~{output_prefix}.raw.bcf \
                --checkpoint-file-out ~{output_prefix}.checkpoint.bin"

        if [ -s "~{output_prefix}.checkpoint.bin" ]; then
            cmd="$cmd --checkpoint-file-in ~{output_prefix}.checkpoint.bin" 
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
        File phased_bcf = "~{output_prefix}.bcf"
        File phased_bcf_csi = "~{output_prefix}.bcf.csi"
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
        docker:             docker,
        checkpointFile: "~{output_prefix}.checkpoint.bin"
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

task GLIMPSE2Ligate {
    input {
        Array[File] phased_bcfs
        Array[File] phased_bcf_csis
        String output_prefix

        String docker

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 2 * ceil(size(phased_bcfs, "GB")) + 10

    command <<<
        set -euox pipefail

        ./GLIMPSE2_ligate --input ~{write_lines(phased_bcfs)} --output ~{output_prefix}.bcf --thread $(nproc)

        # when generating BCF output, the index appears to be corrupt (possibly due to https://github.com/samtools/htslib/issues/1740), so we regenerate with bcftools
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


task PopAndMarkCollisions {
    input {
        # all VCFs should be split to biallelic
        File posteriors_vcf
        File posteriors_vcf_idx
        File panel_bubble_split_vcf
        File panel_bubble_split_vcf_idx
        File panel_id_split_vcf
        File panel_id_split_vcf_idx
        File pop_python_script                      # modified version of convert-to-biallelic.py
        String region
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 3 * ceil(size([posteriors_vcf, panel_bubble_split_vcf, panel_id_split_vcf], "GB"))

    command <<<
        set -euox pipefail

        pypy -m pip install tqdm

        # annotate bubble IDs
        bcftools annotate -r region --regions-overlap 0 -a ~{panel_bubble_split_vcf} ~{posteriors_vcf} \
            -c CHROM,POS,REF,ALT,ID:=INFO/ID,INFO/ID:=INFO/ID | \
        pypy ~{pop_python_script} \
            --panel_id_split_vcf ~{panel_id_split_vcf} \
            bcftools view -W -Ob -o ~{output_prefix}.bcf
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
        docker:             "us.gcr.io/broad-dsde-methods/slee/lrma-aou2-panel-creation-pypy:v1"
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
