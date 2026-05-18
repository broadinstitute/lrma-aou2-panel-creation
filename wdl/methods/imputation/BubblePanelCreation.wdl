version 1.0

workflow BubblePanelCreation {
    input {
        File phased_vcf
        File phased_vcf_idx
        File annotations_vcf        # output of FixVariantCollisions step before Shapeit4
        File annotations_vcf_idx
        File reference_fasta
        File reference_fasta_fai
        String region
        String output_prefix

        # inputs for FixVariantCollisions (see documentation for arguments in task)
        File fix_variant_collisions_script
        Int operation = 1
        String weight_tag = "SCORE"
        Int is_weight_format_field = 0
        Float default_weight = 0.05

        File prepare_vcf_and_add_ids_script
        File merge_vcfs_script
        File cargo_toml
        Float frac_missing = 0.2
    }

    call FixVariantCollisions { input:
        phased_vcf = phased_vcf,
        phased_vcf_idx = phased_vcf_idx,
        annotations_vcf = annotations_vcf,
        annotations_vcf_idx = annotations_vcf_idx,
        fix_variant_collisions_script = fix_variant_collisions_script,
        operation = operation,
        weight_tag = weight_tag,
        is_weight_format_field = is_weight_format_field,
        default_weight = default_weight,
        region = region,
        output_prefix = output_prefix + ".ligated.collisionless"
    }

    call BubblePanelCreation { input:
        phased_vcf = FixVariantCollisions.collisionless_vcf,
        phased_vcf_idx = FixVariantCollisions.collisionless_vcf_idx,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        region = region,
        prepare_vcf_and_add_ids_script = prepare_vcf_and_add_ids_script,
        merge_vcfs_script = merge_vcfs_script,
        cargo_toml = cargo_toml,
        frac_missing = frac_missing,
        weight_tag = weight_tag,
        default_weight = default_weight,
        output_prefix = output_prefix + ".ligated.collisionless.bubble"
    }

    # make sure dict in header
    # TODO add preprocessing steps from KAGE Panel WDL

    output {
        File phased_collisionless_vcf = FixVariantCollisions.collisionless_vcf
        File phased_collisionless_vcf_idx = FixVariantCollisions.collisionless_vcf_idx
        File phased_collisionless_removed_counts_tsv = FixVariantCollisions.collisionless_removed_counts_tsv
        File phased_collisionless_histogram_tsv = FixVariantCollisions.collisionless_histogram_tsv
        File panel_vcf = BubblePanelCreation.panel_vcf
        File panel_vcf_idx = BubblePanelCreation.panel_vcf_idx
        File panel_id_split_vcf = BubblePanelCreation.panel_id_split_vcf
        File panel_id_split_vcf_idx = BubblePanelCreation.panel_id_split_vcf_idx
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

# TODO consolidate with Stage2; break out annotation step
task FixVariantCollisions {
    input {
        File phased_vcf                          # biallelic, can be locally or fully phased
        File phased_vcf_idx
        File annotations_vcf
        File annotations_vcf_idx
        File fix_variant_collisions_script
        Int operation = 1                        # 0=can only remove an entire VCF record; 1=can remove single ones from a GT
        String weight_tag = "SCORE"              # ID of the weight field; weights are assumed to be non-negative; we set to SCORE to prefer kanpig records (and moreover, those with higher SCORE) over DeepVariant records (these should have no SCORE, and will be assigned the low default_weight below)
        Int is_weight_format_field = 0           # given a VCF record in a sample, assign it a weight encoded in the INFO field (0) or in the sample column (1)
        Float default_weight = 0.05              # default weight if the weight field is not found
        String extra_args = "--use-gq --use-af --verbosity 1"  # use GQ then AF as tie breakers
        String region
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * (ceil(size(phased_vcf, "GiB")) + ceil(size(annotations_vcf, "GiB")))

    command <<<
        set -euxo pipefail

        rustc -O ~{fix_variant_collisions_script} -o FixVariantCollisions

        # after FixVariantCollisions, replace all missing alleles (correctly) emitted with reference alleles, since this is expected by PanGenie panel-creation script
        time bcftools annotate --no-version -c CHROM,POS,REF,ALT,ID,INFO/SCORE,INFO/SVLEN,INFO/AN,INFO/AC,INFO/AF -a ~{annotations_vcf} ~{phased_vcf} ~{region} --regions-overlap 0 --threads 2 | \
        ./FixVariantCollisions \
            ~{operation} \
            ~{weight_tag} \
            ~{is_weight_format_field} \
            ~{default_weight} \
            ~{extra_args} \
            --removed-counts-tsv ~{output_prefix}.removed.tsv \
            --histogram-tsv ~{output_prefix}.histogram.tsv | \
        bcftools +setGT --no-version -Ou -- -t . -n 0p | \
            bcftools +fill-tags --no-version --threads 2 --write-index=csi -Ob -o ~{output_prefix}.bcf -- -t AF,AC,AN
    >>>

    output {
        File collisionless_vcf = "~{output_prefix}.bcf"
        File collisionless_vcf_idx = "~{output_prefix}.bcf.csi"
        File collisionless_removed_counts_tsv = "~{output_prefix}.removed.tsv"
        File collisionless_histogram_tsv = "~{output_prefix}.histogram.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          6,
        mem_gb:             8,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  3,
        max_retries:        0,
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

task BubblePanelCreation {
    input {
        File phased_vcf
        File phased_vcf_idx
        File reference_fasta
        File reference_fasta_fai
        String region
        String output_prefix

        File prepare_vcf_and_add_ids_script
        File merge_vcfs_script
        File cargo_toml
        Float frac_missing
        String? weight_tag
        Float? default_weight

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 3 * ceil(size(phased_vcf, "GiB")) + ceil(size(reference_fasta, "GiB")) 

    command <<<
        set -euxo pipefail

        mkdir -p bubble-utils/src/bin
        cp ~{prepare_vcf_and_add_ids_script} bubble-utils/src/bin/prepare_vcf_and_add_ids.rs
        cp ~{merge_vcfs_script} bubble-utils/src/bin/merge_vcfs.rs
        cp ~{cargo_toml} bubble-utils
        cd bubble-utils
        cargo build --release
        cd ..

        bcftools stats -r ~{region} --regions-overlap 0 ~{phased_vcf} --threads 6 > ~{output_prefix}.stats.txt
        bcftools view --no-version -h ~{phased_vcf} > header.txt

        # validate variants against reference, run bubble prepare-vcf and add-ids scripts, split to biallelic, and run bubble merge script;
        # everything should be normalized or in the desired representation at this point
        time bcftools norm --no-version -r ~{region} --regions-overlap 0 --do-not-normalize --check-ref e --fasta-ref ~{reference_fasta} --threads 2 ~{phased_vcf} | \
            ./bubble-utils/target/release/prepare_vcf_and_add_ids --missing ~{frac_missing} | \
            bcftools norm --no-version -m-any --do-not-normalize | tee \
        >(  bcftools view --no-version --write-index=csi -Ob -o ~{output_prefix}.prepare.id.split.bcf ) | \
         (  ./bubble-utils/target/release/merge_vcfs merge \
                --header header.txt \
                -r ~{reference_fasta} \
                --ploidy 2 \
                ~{"--weight " + weight_tag} \
                ~{"--default-weight " + default_weight} | \
            bcftools view --no-version --threads 2 --write-index=csi -Ob -o ~{output_prefix}.prepare.id.split.mergehap.bcf )

        bcftools stats --threads 6 ~{output_prefix}.prepare.id.split.mergehap.bcf > ~{output_prefix}.prepare.id.split.mergehap.stats.txt
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          6,
        mem_gb:             6,
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

    output {
        File input_stats = "~{output_prefix}.stats.txt"
        File panel_stats = "~{output_prefix}.prepare.id.split.mergehap.stats.txt"
        File panel_vcf = "~{output_prefix}.prepare.id.split.mergehap.bcf"
        File panel_vcf_idx = "~{output_prefix}.prepare.id.split.mergehap.bcf.csi"
        File panel_id_split_vcf = "~{output_prefix}.prepare.id.split.bcf"
        File panel_id_split_vcf_idx = "~{output_prefix}.prepare.id.split.bcf.csi"
    }
}
