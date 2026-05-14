version 1.0

import "StatisticalPhasing_Stage0_CreateShards.wdl" as CreateShards
import "StatisticalPhasing_Stage1_FilterAndConcatVcfs.wdl" as FilterAndConcatVcfs
import "StatisticalPhasing_Stage2_FixVariantCollisions.wdl" as FixVariantCollisions
import "StatisticalPhasing_Stage3_GatherSitesAndChunk.wdl" as GatherSitesAndChunk

workflow StatisticalPhasing {

    input {
        File joint_short_vcf
        File joint_short_vcf_idx
        File joint_sv_vcf
        File joint_sv_vcf_idx
        File reference_fasta
        File reference_fasta_fai
        File genetic_maps_tsv
        String chromosome
        String region
        String output_prefix

        # CreateShards
        Int min_variants_per_shard = 450000
        Int max_variants_per_shard = 550000
        Int min_boundary_dist_bp = 15000
        Int min_sv_len = 50

        # FilterAndConcat
        String? filter_and_concat_short_filter_args
        String filter_and_concat_short_view_args = "-i 'MAC>=2 && F_MISSING<=0.4'"
        String filter_and_concat_sv_view_args = "-i 'MAC>=2'"

        # FixVariantCollisions (see documentation for arguments in task)
        File fix_variant_collisions_script
        Int operation = 1
        String weight_tag = "SCORE"
        Int is_weight_format_field = 0
        Float default_weight = 0.05
        String fix_variant_collisions_extra_args = "--use-gq --use-af --verbosity 1"

        # CreateShapeitChunks
        String chunk_extra_args = "--thread $(nproc) --window-size 1000000 --buffer-size 200000 --window-count 50000 --buffer-count 500" # we want window counts to drive the constraints

        Boolean do_shapeit5 = true
        String shapeit4_extra_args = "--thread $(nproc) --use-PS 0.0001"
        String shapeit5_extra_args =  "--thread $(nproc)"
        String filter_common_args = "-i 'MAF>=0.001'"
        
    }

    Map[String, String] genetic_maps_dict = read_map(genetic_maps_tsv)

    call CreateShards.CreateShards as CreateFilterAndConcatShards { input:
        region = region,
        output_prefix = output_prefix,
        entity_name = output_prefix,
        short_vcf = joint_short_vcf,
        short_vcf_idx = joint_short_vcf_idx,
        sv_vcf = joint_sv_vcf,
        sv_vcf_idx = joint_sv_vcf_idx,
        min_variants_per_shard = min_variants_per_shard,
        max_variants_per_shard = max_variants_per_shard,
        min_boundary_dist_bp = min_boundary_dist_bp,
        min_sv_len = min_sv_len
    }

    scatter (s in range(length(CreateFilterAndConcatShards.shard_regions))) {
        String shard_region = CreateFilterAndConcatShards.shard_regions[s]
        
        call FilterAndConcatVcfs.FilterAndConcatVcfs as FilterAndConcatVcfs { input:
            short_vcf = joint_short_vcf,
            short_vcf_idx = joint_short_vcf_idx,
            sv_vcf = joint_sv_vcf,
            sv_vcf_idx = joint_sv_vcf_idx,
            output_prefix = output_prefix + ".shard-" + s,
            reference_fasta = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
            region = shard_region,
            short_view_args = filter_and_concat_short_view_args,
            short_filter_args = filter_and_concat_short_filter_args,
            sv_view_args = filter_and_concat_sv_view_args
        }

        call FixVariantCollisions.FixVariantCollisions as FixVariantCollisions { input:
            phased_vcf = FilterAndConcatVcfs.filter_and_concat_vcf,
            fix_variant_collisions_script = fix_variant_collisions_script,
            output_prefix = output_prefix + ".shard-" + s,
            operation = operation,
            weight_tag = weight_tag,
            is_weight_format_field = is_weight_format_field,
            default_weight = default_weight,
            fix_variant_collisions_extra_args = fix_variant_collisions_extra_args
        }
    }

    call GatherSitesAndChunk.GatherSitesAndChunk as GatherSitesAndChunk { input:
        region = region,
        filter_and_concat_vcfs = FilterAndConcatVcfs.filter_and_concat_vcf,
        filter_and_concat_vcf_idxs = FilterAndConcatVcfs.filter_and_concat_vcf_idx,
        collisionless_vcfs = FixVariantCollisions.collisionless_vcf,
        collisionless_vcf_idxs = FixVariantCollisions.collisionless_vcf_idx,
        output_prefix = output_prefix,
        chunk_extra_args = chunk_extra_args
    }

    Array[String] common_regions = read_lines(GatherSitesAndChunk.common_chunks) # for shapeit4

    scatter (i in range(length(common_regions))) {
        if (!do_shapeit5) {
            # phase all using Shapeit4
            call Shapeit4 as Shapeit4All { input:
                vcf = GatherSitesAndChunk.collisionless_vcf,
                vcf_idx = GatherSitesAndChunk.collisionless_vcf_idx,
                genetic_map = genetic_maps_dict[chromosome],
                region = common_regions[i],
                output_prefix = output_prefix + ".phased",
                extra_args = shapeit4_extra_args
            }
        }
        if (do_shapeit5) {
            call FilterCommon { input:
                vcf = GatherSitesAndChunk.collisionless_vcf,
                vcf_idx = GatherSitesAndChunk.collisionless_vcf_idx,
                output_prefix = output_prefix + ".common",
                region = common_regions[i],
                filter_common_args = filter_common_args
            }
            # phase common scaffold using Shapeit4
            call Shapeit4 as Shapeit4Common { input:
                vcf = FilterCommon.common_vcf,
                vcf_idx = FilterCommon.common_vcf_idx,
                genetic_map = genetic_maps_dict[chromosome],
                region = common_regions[i],
                output_prefix = output_prefix + ".phased",
                extra_args = shapeit4_extra_args
            }
        }
    }

    call LigateVcfs as LigateScaffold { input:
        vcfs = select_all(flatten([Shapeit4All.phased_vcf, Shapeit4Common.phased_vcf])),
        vcf_idxs = select_all(flatten([Shapeit4All.phased_vcf_idx, Shapeit4Common.phased_vcf_idx])),
        output_prefix = output_prefix + ".phased.ligated"
    }

    if (do_shapeit5) {
        # phase rare using Shapeit5
        Array[String] rare_regions = read_lines(GatherSitesAndChunk.rare_chunks) # for shapeit5

        scatter (i in range(length(rare_regions))) {
            call Shapeit5Rare { input:
                vcf = GatherSitesAndChunk.collisionless_vcf,
                vcf_idx = GatherSitesAndChunk.collisionless_vcf_idx,
                scaffold_vcf = LigateScaffold.ligated_vcf,
                scaffold_vcf_idx = LigateScaffold.ligated_vcf_idx,
                genetic_map = genetic_maps_dict[chromosome],
                region = rare_regions[i],
                scaffold_region = common_regions[i],
                output_prefix = output_prefix + ".phased." + "chunk-" + i,
                extra_args = shapeit5_extra_args
            }
        }

        call GatherSitesAndChunk.BcftoolsConcatNaive as ConcatShapeit5 { input:
            vcfs = flatten([Shapeit5Rare.phased_vcf]),
            vcf_idxs = flatten([Shapeit5Rare.phased_vcf_idx]),
            output_prefix = output_prefix + ".phased.concat"
        }
    }

    output {
        File sites_only_vcf = GatherSitesAndChunk.sites_only_vcf
        File sites_only_vcf_idx = GatherSitesAndChunk.sites_only_vcf_idx
        File phased_vcf = select_first([ConcatShapeit5.concatenated_vcf, LigateScaffold.ligated_vcf])
        File phased_vcf_idx = select_first([ConcatShapeit5.concatenated_vcf_idx, LigateScaffold.ligated_vcf_idx])
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

task FilterCommon {
    input {
        File vcf
        File vcf_idx
        String output_prefix
        String region
        String filter_common_args 

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 50 + 4 * ceil(size(vcf, "GiB"))

    command <<<
        set -euxo pipefail

        # filter to common
        bcftools +fill-tags --no-version -r ~{region} ~{vcf} -Ou -- -t AF,AC,AN | \
            bcftools view ~{filter_common_args} \
                -Ob -o ~{output_prefix}.common.bcf
        bcftools index ~{output_prefix}.common.bcf
    >>>

    output {
        File common_vcf = "~{output_prefix}.common.bcf"
        File common_vcf_idx = "~{output_prefix}.common.bcf.csi"
    }
    
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_gb,
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

task Shapeit4 {
    input {
        File vcf
        File vcf_idx
        File genetic_map
        String region
        String output_prefix
        String extra_args = "--thread $(nproc) --use-PS 0.0001"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 4 * ceil(size(vcf, "GiB"))
    
    command <<<
        set -euxo pipefail

        shapeit4.2 --input ~{vcf} \
                --map ~{genetic_map} \
                --region ~{region} \
                --sequencing \
                --output ~{output_prefix}.bcf \
                ~{extra_args}
        bcftools index ~{output_prefix}.bcf
    >>>

    output{
        File phased_vcf = "~{output_prefix}.bcf"
        File phased_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             16,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/shapeit4:v1"
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

task LigateVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 4 * ceil(size(vcfs, "GiB"))

    command <<<
        set -euxo pipefail

        ligate_static --input ~{write_lines(vcfs)} --output ~{output_prefix}.ligate.bcf
        bcftools +fill-tags --no-version --threads $(nproc) ~{output_prefix}.ligate.bcf \
            -Ob -o ~{output_prefix}.bcf -- -t AF,AC,AN
        bcftools index ~{output_prefix}.bcf
    >>>

    output {
        File ligated_vcf = "~{output_prefix}.bcf"
        File ligated_vcf_idx = "~{output_prefix}.bcf.csi"
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
        docker:             "hangsuunc/shapeit5:v1"
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

task Shapeit5Rare {
    input {
        File vcf
        File vcf_idx
        File scaffold_vcf
        File scaffold_vcf_idx
        File genetic_map
        String region
        String scaffold_region
        String output_prefix
        String extra_args = "--thread $(nproc)"
        
        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 4 * ceil(size(vcf, "GiB")) + ceil(size(scaffold_vcf, "GiB"))

    command <<<
        set -euxo pipefail

        # we only need to fill rare in input (common in scaffold should have been imputed or filled previously);
        # this also fills common in input, but those records will be ignored by Shapeit5
        bcftools +setGT --no-version ~{vcf} -Ou -- -t . -n 0p | \
            bcftools +fill-tags --no-version -Ob -o input.bcf -- -t AF,AC,AN
        bcftools index input.bcf

        /shapeit5/phase_rare --input input.bcf \
            --scaffold ~{scaffold_vcf} \
            --map ~{genetic_map} \
            --input-region ~{region} \
            --scaffold-region ~{scaffold_region} \
            --output phased.bcf \
            ~{extra_args}

        bcftools +fill-tags phased.bcf --no-version -Ob -o ~{output_prefix}.bcf -- -t AF,AC,AN
        bcftools index ~{output_prefix}.bcf

    >>>

    output{
        File phased_vcf = "~{output_prefix}.bcf"
        File phased_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/shapeit5:develop"
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
