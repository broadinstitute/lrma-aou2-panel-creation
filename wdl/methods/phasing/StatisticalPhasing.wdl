version 1.0


workflow StatisticalPhasing {

    input {

        File joint_short_vcf
        File joint_short_vcf_tbi
        File? joint_sv_vcf
        File? joint_sv_vcf_tbi
        File reference_fasta
        File reference_fasta_fai
        File genetic_mapping_tsv_for_shapeit
        String chromosome
        String region
        String output_prefix

        # inputs for FixVariantCollisions
        File fix_variant_collisions_java
        Int? operation
        String? weight_tag
        Int? is_weight_format_field

        Int bin_size = 1000000
        String extra_chunk_args = "--thread $(nproc) --window-size 2000000 --buffer-size 200000"

        String filter_and_concat_short_filter_args = "-i 'MAC>=2 && abs(strlen(ALT)-strlen(REF))<50'"
        String filter_and_concat_sv_filter_args = "-i 'MAC>=2 && abs(strlen(ALT)-strlen(REF))>=50'"

        Boolean shapeit5 = true
        Int shapeit4_cpu
        Int shapeit4_memory
        Int shapeit5_cpu
        Int shapeit5_memory
        String shapeit4_common_extra_args = "--thread $(nproc)"
        String shapeit5_rare_extra_args =  "--thread $(nproc)"
        String filter_common_args = "-i 'MAF>=0.01'"
        String filter_rare_args = "-i 'MAF<=0.01'"
        #String shapeit5_phase_rare_filter_args = "-e 'F_MISSING > 0.10 || ALT=\".\" || ALT=\"*\"'"
    }

    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv_for_shapeit)


    call SplitIntoShard as SubsetCreateChunks { input:
        locus = region,
        bin_size = bin_size,
        pad_size = 0,
        output_prefix = output_prefix + ".subset_create_chunks"
    }

    scatter (s_region in SubsetCreateChunks.locuslist) {
        call SubsetVCFStreaming as SubsetVcfShort { input:
            vcf_gz = joint_short_vcf,
            vcf_tbi = joint_short_vcf_tbi,
            locus = s_region
        }
        # if no SV VCF provided, use the short VCF as the input
        # concatenate the SV VCF if provided    
        if (defined(joint_sv_vcf) && defined(joint_sv_vcf_tbi)) {
            call SubsetVCF as SubsetVcfSV { input:
                vcf_gz = select_first([joint_sv_vcf ]),
                vcf_tbi = select_first([joint_sv_vcf_tbi]),
                locus = s_region
            }

            call FilterAndConcatVcfs { input:
                short_vcf = SubsetVcfShort.subset_vcf,
                short_vcf_tbi = SubsetVcfShort.subset_tbi,
                sv_vcf = SubsetVcfSV.subset_vcf,
                sv_vcf_tbi = SubsetVcfSV.subset_tbi,
                output_prefix = output_prefix + ".concat",
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                region = s_region,
                filter_and_concat_short_filter_args = filter_and_concat_short_filter_args,
                filter_and_concat_sv_filter_args = filter_and_concat_sv_filter_args
            }
        } 
    }

    call BcftoolsConcatNaive as ConcatSubsets { input:
        vcfs = select_all(select_first([FilterAndConcatVcfs.filter_and_concat_vcf,SubsetVcfShort.subset_vcf])),
        vcf_tbis = select_all(select_first([FilterAndConcatVcfs.filter_and_concat_vcf_tbi,SubsetVcfShort.subset_tbi])),
        output_prefix = output_prefix + ".subset.concat"
    }
    # added variant collision fix step
    call FixVariantCollisions { input:
        phased_bcf = ConcatSubsets.concatenated_vcf,
        fix_variant_collisions_java = fix_variant_collisions_java,
        operation = operation,
        weight_tag = weight_tag,
        is_weight_format_field = is_weight_format_field,
        output_prefix = output_prefix
    }

    call CreateChunks as CreateChunks { input:
        vcf = FixVariantCollisions.phased_collisionless_bcf,
        tbi = FixVariantCollisions.phased_collisionless_bcf_index,
        region = region,
        extra_chunk_args = extra_chunk_args
    }

    Array[String] region_list = read_lines(CreateChunks.chunks)

    scatter (i in range(length(region_list))) {
        # phase common using shapeit4
        if (!shapeit5) {
            call Shapeit4 as Shapeit4_all { input:
                vcf_input = FixVariantCollisions.phased_collisionless_bcf,
                vcf_index = FixVariantCollisions.phased_collisionless_bcf_index,
                mappingfile = genetic_mapping_dict[chromosome],
                region = region_list[i],
                output_prefix = output_prefix + ".filter_and_concat.phased",
                cpu = shapeit4_cpu,
                memory = shapeit4_memory,
                extra_args = shapeit4_common_extra_args
            }
        }
        if (shapeit5) {
            call FilterCommonandRareVariants { input:
                vcf_gz = FixVariantCollisions.phased_collisionless_bcf,
                vcf_gz_tbi = FixVariantCollisions.phased_collisionless_bcf_index,
                output_prefix = output_prefix + ".filter_common_and_rare",
                region = region_list[i],
                filter_common_args = filter_common_args,
                filter_rare_args = filter_rare_args
            }
            call Shapeit4 as Shapeit4_common { input:
                vcf_input = FilterCommonandRareVariants.filter_common_vcf,
                vcf_index = FilterCommonandRareVariants.filter_common_vcf_tbi,
                mappingfile = genetic_mapping_dict[chromosome],
                region = region_list[i],
                output_prefix = output_prefix + ".filter_and_concat.common",
                cpu = shapeit4_cpu,
                memory = shapeit4_memory,
                extra_args = shapeit4_common_extra_args
            }
        }
        
    }

    call LigateVcfs as LigateScaffold { input:
        vcfs = select_all(flatten([Shapeit4_common.phased_bcf,Shapeit4_all.phased_bcf])),
        vcf_idxs = select_all(flatten([Shapeit4_common.phased_bcf_index, Shapeit4_all.phased_bcf_index])),
        output_prefix = output_prefix + ".scaffold.ligated"
    }

    # phase rare
    if (shapeit5) {
        scatter (i in range(length(region_list))) {
            call Shapeit5PhaseRare as Shapeit5_phase_rare { input:
                vcf_input = FixVariantCollisions.phased_collisionless_bcf,
                vcf_index = FixVariantCollisions.phased_collisionless_bcf_index,
                scaffold_bcf = LigateScaffold.ligated_vcf_gz,
                scaffold_bcf_index = LigateScaffold.ligated_vcf_gz_tbi,
                mappingfile = genetic_mapping_dict[chromosome],
                chunk_region = region_list[i],
                scaffold_region = region,
                output_prefix = output_prefix + ".chunk.phase.rare.phased",
                chunknum = i,
                cpu = shapeit5_cpu,
                memory = shapeit5_memory,
                extra_args = shapeit5_rare_extra_args,
                #shapeit5_phase_rare_filter_args = shapeit5_phase_rare_filter_args
            }
        }

        call LigateVcfs as LigateRare { input:
            vcfs = select_all(flatten([Shapeit5_phase_rare.chunk_vcf])),
            vcf_idxs = select_all(flatten([Shapeit5_phase_rare.chunk_vcf_index])),
            output_prefix = output_prefix + ".phase.rare.concat"
        }

        # concat common and rare
        call BcftoolsConcatNaive as ConcatCommonRare { input:
            vcfs = [LigateScaffold.ligated_vcf_gz, LigateRare.ligated_vcf_gz],
            vcf_tbis = [LigateScaffold.ligated_vcf_gz_tbi, LigateRare.ligated_vcf_gz_tbi],
            output_prefix = output_prefix + ".phase.common.rare.concat"
        }

    }

    output {
        File phased_vcf = select_first([ConcatCommonRare.concatenated_vcf,LigateScaffold.ligated_vcf_gz])
        File phased_vcf_tbi = select_first([ConcatCommonRare.concatenated_vcf_tbi,LigateScaffold.ligated_vcf_gz_tbi])
    }
}


struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}


task SubsetVCF {

    meta {
        description: "Subset a VCF file to a given locus"
    }

    parameter_meta {
        vcf_gz: "VCF file to be subsetted"
        vcf_tbi: "Tabix index for the VCF file"
        locus: "Locus to be subsetted"
        output_prefix: "output_prefix for the output file"
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        File vcf_gz
        File? vcf_tbi
        String locus
        String output_prefix = "subset"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([vcf_gz, vcf_tbi], "GB")) + 100

    command <<<
        set -euxo pipefail
        if ! ~{defined(vcf_tbi)}; then
            bcftools index ~{vcf_gz}
        fi
        bcftools view ~{vcf_gz} --regions ~{locus} -O b -o ~{output_prefix}.bcf
        bcftools index ~{output_prefix}.bcf
    >>>

    output {
        File subset_vcf = "~{output_prefix}.bcf"
        File subset_tbi = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task SubsetVCFStreaming {

    meta {
        description: "Subset a VCF file to a given locus"
    }

    parameter_meta {
        vcf_gz: {
            description: "VCF file to be subsetted",
            localization_optional: true
        }
        vcf_tbi: {
            description: "Tabix index for the VCF file",
            localization_optional: true
        }
        locus: "Locus to be subsetted"
        output_prefix: "output_prefix for the output file"
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        File vcf_gz
        File vcf_tbi
        String locus
        String output_prefix = "subset"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([vcf_gz, vcf_tbi], "GB")) + 100

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        bcftools view --no-version ~{vcf_gz} --regions ~{locus} -Oz -o ~{output_prefix}.vcf.gz
        bcftools index -t ~{output_prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{output_prefix}.vcf.gz"
        File subset_tbi = "~{output_prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


task Shapeit4 {
    input{
        File vcf_input
        File vcf_index
        File mappingfile
        String region
        String output_prefix
        Int cpu
        Int memory
        String extra_args = "--thread $(nproc)"

        RuntimeAttr? runtime_attr_override
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
    command <<<
        set -euxo pipefail

        shapeit4.2 --input ~{vcf_input} \
                --map ~{mappingfile} \
                --region ~{region} \
                --sequencing \
                --output ~{output_prefix}.bcf \
                ~{extra_args}
        bcftools +fill-tags ~{output_prefix}.bcf -Ob -o ~{output_prefix}.tagged.bcf -- -t AN,AC
        bcftools index ~{output_prefix}.tagged.bcf
    >>>

    output{
        # File resouce_monitor_log = "resources.log"
        File phased_bcf = "~{output_prefix}.tagged.bcf"
        File phased_bcf_index = "~{output_prefix}.tagged.bcf.csi"
    }

    #Int disk_size = 100 + ceil(2 * size(vcf_input, "GiB"))

 #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpu,
        mem_gb:             memory,
        disk_gb:            100,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/shapeit4:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Shapeit5PhaseCommon{
    input{
        File vcf_input
        File vcf_index
        File mappingfile
        String region
        String output_prefix
        Int cpu
        Int memory
        String extra_args = "--thread $(nproc)"
        Float minimal_maf = 0.01

        RuntimeAttr? runtime_attr_override
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
    command <<<
        set -euxo pipefail
        # add AN AC tag
        bcftools +fill-tags ~{vcf_input} -Ob -o tmp.out.bcf -- -t AN,AC
        bcftools index tmp.out.bcf
        phase_common_static --input tmp.out.bcf \
                            --filter-maf ~{minimal_maf} \
                            --region ~{region} \
                            --map ~{mappingfile} \
                            --output scaffold.bcf \
                            ~{extra_args}
        bcftools +fill-tags scaffold.bcf -Ob -o ~{output_prefix}.scaffold.bcf -- -t AN,AC
        bcftools index ~{output_prefix}.scaffold.bcf
    >>>

    output{
        File scaffold_vcf = "~{output_prefix}.scaffold.bcf"
        File scaffold_vcf_index = "~{output_prefix}.scaffold.bcf.csi"
    }

    Int disk_size = 100 + ceil(2 * size(vcf_input, "GiB"))

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpu,
        mem_gb:             memory,
        disk_gb:            100,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "hangsuunc/shapeit5:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task CreateChunks {

    input {
        File vcf
        File tbi
        String region
        String? extra_chunk_args = "--thread $(nproc) --window-size 5000000 --buffer-size 500000"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size([vcf, tbi], "GB")) + 1

    command <<<
        set -euxo pipefail

        wget https://github.com/odelaneau/GLIMPSE/releases/download/v1.1.1/GLIMPSE_chunk_static
        chmod +x GLIMPSE_chunk_static

        ./GLIMPSE_chunk_static \
            -I ~{vcf} \
            --region ~{region} \
            ~{extra_chunk_args} \
            -O chunks.txt

        # cut chunks + buffers
        cut -f 3 chunks.txt > chunks.regions.txt
    >>>

    output {
        File chunks = "chunks.regions.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.2"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task LigateVcfs {

    input {
        Array[File] vcfs
        Array[File]? vcf_idxs
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(vcfs, "GB")) + 100

    command <<<
        set -euxo pipefail
        if ! ~{defined(vcf_idxs)}; then
            for ff in ~{sep=' ' vcfs}; do bcftools index $ff; done
        fi

        ligate_static --input ~{write_lines(vcfs)} --output ~{output_prefix}.vcf.gz
        bcftools index -t ~{output_prefix}.vcf.gz
    >>>

    output {
        File ligated_vcf_gz = "~{output_prefix}.vcf.gz"
        File ligated_vcf_gz_tbi = "~{output_prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:"hangsuunc/shapeit5:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


task Shapeit5PhaseRare{
    input{
        File vcf_input
        File vcf_index
        File scaffold_bcf
        File scaffold_bcf_index
        File mappingfile
        String chunk_region
        String scaffold_region
        String output_prefix
        Int chunknum
        Int cpu
        Int memory
        String extra_args = "--thread $(nproc)"
        #String shapeit5_phase_rare_filter_args = "-e 'F_MISSING > 0.10 || ALT=\".\" || ALT=\"*\"'"

        RuntimeAttr? runtime_attr_override
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
    command <<<
        set -euxo pipefail

        bcftools +fill-tags ~{scaffold_bcf} -Ob -o tmp.scaffold.out.bcf -- -t AN,AC
        bcftools index tmp.scaffold.out.bcf

        bcftools +fill-tags ~{vcf_input} -Ob -o tmp.out.bcf -- -t AN,AC
        bcftools index tmp.out.bcf
        
        # try to fix bugs in https://github.com/odelaneau/shapeit5/issues/33
        # replace filtering with setGT to set missing genotypes to 0|0
        # try to fix the issue of ID starts with numbers, will revisit later
        bcftools +setGT tmp.out.bcf -- -t . -n 0p | \
            bcftools annotate -x 'ID' \
                -Ob -o tmp.rare.out.bcf
        bcftools index tmp.rare.out.bcf
        
        phase_rare_static --input tmp.rare.out.bcf \
                    --scaffold tmp.scaffold.out.bcf \
                    --map ~{mappingfile} \
                    --input-region ~{chunk_region} \
                    --scaffold-region ~{scaffold_region} \
                    --output ~{output_prefix}.chunk.~{chunknum}.bcf \
                    ~{extra_args}

        bcftools +fill-tags ~{output_prefix}.chunk.~{chunknum}.bcf -Ob -o ~{output_prefix}.chunk.~{chunknum}.tagged.bcf -- -t AN,AC
        bcftools index ~{output_prefix}.chunk.~{chunknum}.tagged.bcf

    >>>

    output{
        File chunk_vcf = "~{output_prefix}.chunk.~{chunknum}.tagged.bcf"
        File chunk_vcf_index = "~{output_prefix}.chunk.~{chunknum}.tagged.bcf.csi"
    }

    Int disk_size = 100 + ceil(2 * size(vcf_input, "GiB"))

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpu,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "hangsuunc/shapeit5:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


task BcftoolsConcatBCFs {

    input {
        Array[File] vcfs
        Array[File]? vcf_idxs
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(vcfs, "GB")) + 1000

    command <<<
        set -euxo pipefail
        if ! ~{defined(vcf_idxs)}; then
            for ff in ~{sep=' ' vcfs}; do bcftools index $ff; done
        fi

        bcftools concat --allow-overlap --remove-duplicates -Ob -o ~{output_prefix}.bcf -f ~{write_lines(vcfs)} 
        bcftools sort ~{output_prefix}.bcf -Ob -o ~{output_prefix}.sorted.bcf
        bcftools index ~{output_prefix}.sorted.bcf
    >>>

    output {
        File concated_bcf = "~{output_prefix}.sorted.bcf"
        File concated_bcf_index = "~{output_prefix}.sorted.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:"hangsuunc/shapeit5:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FilterCommonandRareVariants {

    input {
        File vcf_gz        
        File vcf_gz_tbi
        String output_prefix
        String region
        String filter_common_args 
        String filter_rare_args

        RuntimeAttr? runtime_attr_override
    }

    command {
        set -euxo pipefail

        # filter common
        bcftools +fill-tags -r ~{region} ~{vcf_gz} -- -t AF,AC,AN | \
            bcftools view ~{filter_common_args} \
                -Oz -o ~{output_prefix}.common.vcf.gz
        bcftools index -t ~{output_prefix}.common.vcf.gz

        # filter rare
        bcftools +fill-tags -r ~{region} ~{vcf_gz} -- -t AF,AC,AN | \
            bcftools view ~{filter_rare_args} \
                -Oz -o ~{output_prefix}.rare.vcf.gz
        bcftools index -t ~{output_prefix}.rare.vcf.gz

    }

    output {
        File filter_common_vcf = "~{output_prefix}.common.vcf.gz"
        File filter_common_vcf_tbi = "~{output_prefix}.common.vcf.gz.tbi"
        File filter_rare_vcf = "~{output_prefix}.rare.vcf.gz"
        File filter_rare_vcf_tbi = "~{output_prefix}.rare.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            1000,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:"hangsuunc/shapeit5:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task FilterAndConcatVcfs {

    input {
        File short_vcf         # multiallelic
        File short_vcf_tbi
        File sv_vcf            # biallelic
        File sv_vcf_tbi
        String output_prefix
        String region
        File reference_fasta
        File reference_fasta_fai
        String? filter_and_concat_short_filter_args = "-i 'MAC>=2 && abs(strlen(ALT)-strlen(REF))<50'"
        String? filter_and_concat_sv_filter_args = "-i 'MAC>=2 && abs(strlen(ALT)-strlen(REF))>=50'"

        RuntimeAttr? runtime_attr_override
    }

    command {
        set -euxo pipefail

        # filter SV
        bcftools +fill-tags -r ~{region} ~{sv_vcf} -- -t AF,AC,AN | \
            bcftools view ~{filter_and_concat_sv_filter_args} \
                -Oz -o ~{output_prefix}.SV.vcf.gz
        bcftools index -t ~{output_prefix}.SV.vcf.gz

        # split to biallelic and filter short
        bcftools norm -r ~{region} -m-any -N -f ~{reference_fasta} ~{short_vcf} | \
            bcftools +fill-tags -- -t AF,AC,AN | \
            bcftools view ~{filter_and_concat_short_filter_args} | \
            bcftools sort -Oz -o ~{output_prefix}.short.vcf.gz
        bcftools index -t ~{output_prefix}.short.vcf.gz

        # concatenate with deduplication; providing SV VCF as first argument preferentially keeps those records
        bcftools concat \
            ~{output_prefix}.SV.vcf.gz \
            ~{output_prefix}.short.vcf.gz \
            --allow-overlaps --remove-duplicates | \
            bcftools sort -Oz -o ~{output_prefix}.vcf.gz
        bcftools index -t ~{output_prefix}.vcf.gz
    }

    output {
        File filter_and_concat_vcf = "~{output_prefix}.vcf.gz"
        File filter_and_concat_vcf_tbi = "~{output_prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            1000,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        1,
        docker:"hangsuunc/shapeit5:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task SplitIntoShard {
    input {
        String locus
        Int bin_size
        Int pad_size
        String output_prefix

        Int? preemptible_tries
    }


    command <<<
        set -eo pipefail

        python - --locus ~{locus} \
                 --bin_size ~{bin_size} \
                 --pad_size ~{pad_size} \
                 --output_file ~{output_prefix} \
                 <<-'EOF'
        import argparse

        def split_locus(locus):
            chromosome, span = locus.split(":")
            start, end = span.split("-")
            return(chromosome, int(start), int(end))

        def split_locus_to_intervals(locus, bin_size, pad_size):
            chromo, start, end = split_locus(locus)
            bin_num = (end - start)//bin_size
            intervals = [(chromo, start, start + bin_size + pad_size)]
            for i in range(1, bin_num):
                start_pos = start + i*bin_size - pad_size
                end_pos = start_pos + bin_size + pad_size
                intervals.append((chromo, start_pos, end_pos))
            if end > start + bin_num*bin_size:
                intervals.append((chromo, start + bin_num*bin_size - pad_size, end))
            return(intervals)

        def write_output_file(content, output_file):
            with open(output_file, "w") as f:
                for item in content:
                    l = "%s:%d-%d" % (item[0], item[1], item[2])
                    f.write(l+ "\n")

        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--locus',
                                type=str)

            parser.add_argument('--output_file',
                                type=str)

            parser.add_argument('--bin_size',
                    type=int)

            parser.add_argument('--pad_size',
                    type=int)

            args = parser.parse_args()

            intervals = split_locus_to_intervals(args.locus, args.bin_size, args.pad_size)
            write_output_file(intervals, args.output_file + ".txt")

        if __name__ == "__main__":
            main()
        EOF

    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/slee/kage-lite:pr_29"
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
    }

    output {
        Array[String] locuslist = read_lines("~{output_prefix}.txt")
    }
}

task BcftoolsConcatNaive {
    input {
        Array[File] vcfs
        Array[File]? vcf_tbis
        String output_prefix
    }

    command <<<
        set -euxo pipefail

        # Index all input VCF files if no precomputed tbis provided.
        if ! ~{defined(vcf_tbis)}; then
            for vcf in ~{sep=" " vcfs}; do
                bcftools index "$vcf"
            done
        fi

        bcftools concat \
            ~{sep=" " vcfs} \
            -n \
            --no-version | \
            bcftools sort -Oz -o ~{output_prefix}.vcf.gz
        bcftools index -t ~{output_prefix}.vcf.gz
    >>>

    output {
        File concatenated_vcf = "~{output_prefix}.vcf.gz"
        File concatenated_vcf_tbi = "~{output_prefix}.vcf.gz.tbi"
    }

    runtime {
        docker: "hangsuunc/shapeit5:v1"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 1000 SSD"
    }
}

task FixVariantCollisions {

    input {
        File phased_bcf                     # biallelic
        File fix_variant_collisions_java
        Int operation = 1                   # 0=can only remove an entire VCF record; 1=can remove single ones from a GT
        String weight_tag = "UNIT_WEIGHT"   # ID of the weight field; if this field is not found, all weights are set to one; weights are assumed to be non-negative
        Int is_weight_format_field = 0      # given a VCF record in a sample, assign it a weight encoded in the sample column (1) or in the INFO field (0)
        String output_prefix
    }

    command <<<
        set -euxo pipefail

        # convert bcf to vcf.gz
        bcftools view ~{phased_bcf} -Oz -o phased.vcf.gz
        bcftools index -t phased.vcf.gz

        java ~{fix_variant_collisions_java} \
            phased.vcf.gz \
            ~{operation} \
            ~{weight_tag} \
            ~{is_weight_format_field} \
            collisionless.vcf \
            windows.txt \
            histogram.txt \
            null                            # do not output figures

        # replace all missing alleles (correctly) emitted with reference alleles, since this is expected by PanGenie panel-creation script
        bcftools view collisionless.vcf | \
            sed -e 's/\.|0/0|0/g' | sed -e 's/0|\./0|0/g' | sed -e 's/\.|1/0|1/g' | sed -e 's/1|\./1|0/g' | sed -e 's/\.|\./0|0/g' | \
            bcftools view -Oz -o ~{output_prefix}.phased.collisionless.vcf.gz
        # index and convert via vcf.gz to avoid errors from missing header lines
        bcftools index -t ~{output_prefix}.phased.collisionless.vcf.gz
        bcftools view ~{output_prefix}.phased.collisionless.vcf.gz -Ob -o ~{output_prefix}.phased.collisionless.bcf
        bcftools index ~{output_prefix}.phased.collisionless.bcf
    >>>

    output {
        File phased_collisionless_bcf = "~{output_prefix}.phased.collisionless.bcf"
        File phased_collisionless_bcf_index = "~{output_prefix}.phased.collisionless.bcf.csi"
        File windows = "windows.txt"
        File histogram = "histogram.txt"
    }
    ###################
    runtime {
        cpu: 1
        memory:  "16 GiB"
        disks: "local-disk 1000 HDD"
        bootDiskSizeGb: 10
        preemptible_tries:     0
        max_retries:           0
        docker:"us.gcr.io/broad-gatk/gatk:4.6.0.0"
    }
}
