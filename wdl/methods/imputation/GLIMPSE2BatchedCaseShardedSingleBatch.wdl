version 1.0

struct RuntimeAttributes {
    Int? cpu
    Int? command_mem_gb
    Int? additional_mem_gb
    Int? disk_size_gb
    Int? boot_disk_size_gb
    Boolean? use_ssd
    Int? preemptible
    Int? max_retries
}

workflow GLIMPSE2BatchedCaseShardedSingleBatch {
    input {
        File input_vcf
        File input_vcf_idx
        Array[String] sample_names

        # per chromosome
        Array[String]+ chromosomes
        Array[File]+ genetic_maps
        Array[File] panel_split_vcf          # "split" here means "split to biallelic"; "split" below means "chunked"
        Array[File] panel_split_vcf_idx

        String extra_chunk_args = "--thread $(nproc) --window-mb 5 --buffer-mb 0.5 --sequential"
        String extra_split_args = "--keep-monomorphic-ref-sites"
        String extra_phase_args = "--impute-reference-only-variants --keep-monomorphic-ref-sites"
        String output_prefix

        # inputs for PreprocessPLs
        File remap_simple_bubble_likelihoods_python_script
        File swap_alleles_python_script

        # inputs for FixVariantCollisions
        File annotations_vcf
        File annotations_vcf_idx
        File fix_variant_collisions_script
        Int operation
        String weight_tag
        Int is_weight_format_field

        String docker
    }

    scatter (j in range(length(chromosomes))) {
        String chromosome = chromosomes[j]

        call GLIMPSE2Chunk as ChromosomeGLIMPSE2Chunk {
            input:
                vcf = panel_split_vcf[j],
                vcf_idx = panel_split_vcf_idx[j],
                region = chromosome,
                genetic_map = genetic_maps[j],
                output_prefix = output_prefix + "." + chromosome,
                extra_chunk_args = extra_chunk_args,
                docker = docker
        }

        Array[String] input_regions = read_lines(ChromosomeGLIMPSE2Chunk.input_regions)
        Array[String] output_regions = read_lines(ChromosomeGLIMPSE2Chunk.output_regions)

        scatter (k in range(length(input_regions))) {
            call GLIMPSE2SplitReference as ChunkedGLIMPSE2SplitReference {
                input:
                    panel_split_vcf = panel_split_vcf[j],
                    panel_split_vcf_idx = panel_split_vcf_idx[j],
                    input_region = input_regions[k],
                    output_region = output_regions[k],
                    genetic_map = genetic_maps[j],
                    output_prefix = output_prefix + "." + chromosome + ".shard-" + k + ".split",
                    extra_split_args = extra_split_args,
                    docker = docker
            }

            call PreprocessPLs as ChunkedPreprocessPLs {
                input:
                    input_vcf = input_vcf,
                    input_vcf_idx = input_vcf_idx,
                    panel_split_vcf = panel_split_vcf[j],
                    panel_split_vcf_idx = panel_split_vcf_idx[j],
                    output_region = output_regions[k],
                    sample_names = sample_names,
                    output_prefix = output_prefix + "." + chromosome + ".shard-" + k + ".preprocessedPLs",
                    remap_simple_bubble_likelihoods_python_script = remap_simple_bubble_likelihoods_python_script,
                    swap_alleles_python_script = swap_alleles_python_script,
                    docker = docker
            }

            call GLIMPSE2Phase as ChunkedGLIMPSE2Phase {
                input:
                    input_vcf = ChunkedPreprocessPLs.preprocessed_pls_bcf,
                    input_vcf_idx = ChunkedPreprocessPLs.preprocessed_pls_bcf_csi,
                    panel_split_chunk_bin = ChunkedGLIMPSE2SplitReference.panel_split_chunk_bin,
                    input_region = input_regions[k],
                    output_region = output_regions[k],
                    sample_names = sample_names,
                    genetic_map = genetic_maps[j],
                    output_prefix = output_prefix + "." + chromosome + ".shard-" + k + ".phased",
                    extra_phase_args = extra_phase_args,
                    docker = docker
            }
        }

        call GLIMPSE2Ligate as ChromosomeGLIMPSE2Ligate {
            input:
                phased_bcfs = ChunkedGLIMPSE2Phase.phased_bcf,
                phased_bcf_csis = ChunkedGLIMPSE2Phase.phased_bcf_csi,
                output_prefix = output_prefix + "." + chromosome + ".ligated",
                docker = docker
        }

        call FixVariantCollisions as ChromosomeGLIMPSE2PosteriorsCollisionless { input:
            phased_vcf = ChromosomeGLIMPSE2Ligate.ligated_vcf_gz,
            phased_vcf_idx = ChromosomeGLIMPSE2Ligate.ligated_vcf_gz_tbi,
            annotations_vcf = annotations_vcf,
            annotations_vcf_idx = annotations_vcf_idx,
            fix_variant_collisions_script = fix_variant_collisions_script,
            operation = operation,
            weight_tag = weight_tag,
            is_weight_format_field = is_weight_format_field,
            output_prefix = output_prefix + "." + chromosome + ".glimpse2.collisionless"
        }
    }

    call ConcatVcfs as GLIMPSE2PosteriorsConcatVcfs {
        input:
            vcfs = ChromosomeGLIMPSE2Ligate.ligated_vcf_gz,
            vcf_idxs = ChromosomeGLIMPSE2Ligate.ligated_vcf_gz_tbi,
            output_prefix = output_prefix + ".glimpse2.posteriors",
            docker = docker
    }

    call ConcatVcfs as GLIMPSE2PosteriorsCollisionlessConcatVcfs {
        input:
            vcfs = ChromosomeGLIMPSE2PosteriorsCollisionless.phased_collisionless_vcf,
            vcf_idxs = ChromosomeGLIMPSE2PosteriorsCollisionless.phased_collisionless_vcf_idx,
            output_prefix = output_prefix + ".glimpse2.collisionless",
            docker = docker
    }

    output {
        Array[Array[File]] chromosome_panel_split_chunk_bins = ChunkedGLIMPSE2SplitReference.panel_split_chunk_bin
        File glimpse2_posteriors_vcf_gz = GLIMPSE2PosteriorsConcatVcfs.vcf_gz
        File glimpse2_posteriors_vcf_gz_tbi = GLIMPSE2PosteriorsConcatVcfs.vcf_gz_tbi
        File glimpse2_posteriors_collisionless_vcf_gz = GLIMPSE2PosteriorsCollisionlessConcatVcfs.vcf_gz
        File glimpse2_posteriors_collisionless_vcf_gz_tbi = GLIMPSE2PosteriorsCollisionlessConcatVcfs.vcf_gz_tbi
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

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_chunk_static
        chmod +x GLIMPSE2_chunk_static

        ./GLIMPSE2_chunk_static \
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
        File panel_split_vcf
        File panel_split_vcf_idx
        String input_region
        String output_region
        File genetic_map
        String output_prefix
        String? extra_split_args
        String docker

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 2 * ceil(size(panel_split_vcf, "GB")) + 10

    command <<<
        set -euxo pipefail

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_split_reference_static
        chmod +x GLIMPSE2_split_reference_static

        ./GLIMPSE2_split_reference_static \
            -R ~{panel_split_vcf} \
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


task PreprocessPLs {
    input {
        File input_vcf
        File input_vcf_idx
        File panel_split_vcf
        File panel_split_vcf_idx
        String output_region
        Array[String] sample_names
        String output_prefix

        File remap_simple_bubble_likelihoods_python_script
        File swap_alleles_python_script

        String docker

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 2 * ceil(size([input_vcf, panel_split_vcf], "GB")) + 10

    command {
        set -euxo pipefail

       # TODO add gcloud to Docker
#        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        # TODO stream; or, once shards are fixed, prepare beforehand?
        bcftools view --no-version -G ~{panel_split_vcf}##idx##~{panel_split_vcf_idx} \
            --regions-overlap pos -r ~{output_region} -Ou | \
        bcftools norm -m+any -N \
            --write-index=tbi -Oz -o panel.subset.sites.vcf.gz

        pypy -m pip install --no-input tqdm

        # TODO stream
        bcftools view ~{input_vcf}##idx##~{input_vcf_idx} \
            -r ~{output_region} \
            --regions-overlap pos \
            -S ~{write_lines(sample_names)} \
            --threads $(nproc) | \
        pypy ~{remap_simple_bubble_likelihoods_python_script} \
            --bubble panel.subset.sites.vcf.gz | \
        bcftools +tag2tag -Ou -- --LPL-to-PL | \
        bcftools norm -m-any -Ou | \
        bcftools filter -i 'INFO/BMAP != "."' | \
        pypy ~{swap_alleles_python_script} | \
        bcftools annotate -x INFO/BMAP,INFO/BUBBLE,INFO/BPOS,INFO/BREF,INFO/BALT \
            --write-index=csi -Ob -o ~{output_prefix}.bcf
    }

    output {
        File preprocessed_pls_bcf = "~{output_prefix}.bcf"
        File preprocessed_pls_bcf_csi = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
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

task GLIMPSE2Phase {
    input {
        File input_vcf
        File input_vcf_idx
        File panel_split_chunk_bin
        String input_region
        String output_region
        Array[String] sample_names
        File genetic_map
        String output_prefix
        String? extra_phase_args

        String docker

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 2 * ceil(size([input_vcf, panel_split_chunk_bin], "GB")) + 10

    command {
        set -euxo pipefail

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_phase_static
        chmod +x GLIMPSE2_phase_static

        ./GLIMPSE2_phase_static \
            --input-gl ~{input_vcf} \
            -R ~{panel_split_chunk_bin} \
            --thread $(nproc) \
            ~{extra_phase_args} \
            --output ~{output_prefix}.raw.bcf

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
        cpu_cores:          1,
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

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v2.0.1/GLIMPSE2_ligate_static
        chmod +x GLIMPSE2_ligate_static

        ./GLIMPSE2_ligate_static --input ~{write_lines(phased_bcfs)} --output ~{output_prefix}.vcf.gz --thread $(nproc)
    >>>

    output {
        File ligated_vcf_gz = "~{output_prefix}.vcf.gz"
        File ligated_vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
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

task FixVariantCollisions {
    input {
        File phased_vcf                     # biallelic
        File phased_vcf_idx
        File annotations_vcf
        File annotations_vcf_idx
        File fix_variant_collisions_script
        Int operation = 1                   # 0=can only remove an entire VCF record; 1=can remove single ones from a GT
        String weight_tag = "SCORE"         # ID of the weight field; weights are assumed to be non-negative; we set to SCORE to prefer kanpig records (and moreover, those with higher SCORE) over DeepVariant records (these should have no SCORE, and will be assigned the low default_weight below)
        Int is_weight_format_field = 0      # given a VCF record in a sample, assign it a weight encoded in the INFO field (0) or in the sample column (1)
        Float default_weight = 0.1          # default weight if the weight field is not found
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 50 + 4 * (ceil(size(phased_vcf, "GiB")))

    command <<<
        set -euxo pipefail

        rustc -O ~{fix_variant_collisions_script} -o FixVariantCollisions

        # after FixVariantCollisions, replace all missing alleles (correctly) emitted with reference alleles, since this is expected by PanGenie panel-creation script
        time bcftools annotate --no-version -c CHROM,POS,REF,ALT,ID,INFO/SCORE,INFO/SVLEN -a ~{annotations_vcf} ~{phased_vcf} --threads 2 | \
        ./FixVariantCollisions \
            ~{operation} \
            ~{weight_tag} \
            ~{is_weight_format_field} \
            ~{default_weight} \
            histogram.txt | \
        bcftools +setGT --no-version -Ou -- -t . -n 0p | \
            bcftools +fill-tags --no-version --threads 2 --write-index=csi -Ob -o ~{output_prefix}.phased.collisionless.bcf -- -t AF,AC,AN
    >>>

    output {
        File phased_collisionless_vcf = "~{output_prefix}.phased.collisionless.bcf"
        File phased_collisionless_vcf_idx = "~{output_prefix}.phased.collisionless.bcf.csi"
        File histogram = "histogram.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          6,
        mem_gb:             8,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsde-methods/slee/pangenie-panel-creation-rust:v1"
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

task ConcatVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String output_prefix
        String docker

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size_gb = 3 * ceil(size(vcfs, "GB")) + 10

    command {
        set -euox pipefail

        # TODO FIX LEXICOGRAPHICAL BUG!
        mkdir inputs
        mv ~{sep=' ' vcfs} inputs
        mv ~{sep=' ' vcf_idxs} inputs

        if [ $(ls inputs/*.vcf.gz | wc -l) == 1 ]
        then
            cp $(ls inputs/*.vcf.gz) ~{output_prefix}.vcf.gz
            cp $(ls inputs/*.vcf.gz.tbi) ~{output_prefix}.vcf.gz.tbi
        else
            bcftools concat $(ls inputs/*.vcf.gz) --naive -Oz -o ~{output_prefix}.vcf.gz
            bcftools index -t ~{output_prefix}.vcf.gz
        fi
    }

    output {
        File vcf_gz = "~{output_prefix}.vcf.gz"
        File vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
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
