version 1.0

workflow PanGeniePanelCreation {
    input {
        File phased_vcf
        File phased_vcf_idx
        File reference_fasta
        String region
        String output_prefix

        # inputs for FixVariantCollisions (see documentation for arguments in task)
        File fix_variant_collisions_script
        Int operation = 1
        String weight_tag = "AF"
        Int is_weight_format_field = 0
        Float default_weight = 0.1

        File prepare_vcf_and_add_ids_script
        File merge_vcfs_script
        File cargo_toml
        Float frac_missing = 0.2
    }

    call FixVariantCollisions { input:
        phased_vcf = phased_vcf,
        fix_variant_collisions_script = fix_variant_collisions_script,
        operation = operation,
        weight_tag = weight_tag,
        is_weight_format_field = is_weight_format_field,
        default_weight = default_weight,
        output_prefix = output_prefix
    }

    call PanGeniePanelCreation { input:
        phased_vcf = FixVariantCollisions.phased_collisionless_vcf,
        phased_vcf_idx = FixVariantCollisions.phased_collisionless_vcf_idx,
        reference_fasta = reference_fasta,
        region = region,
        prepare_vcf_and_add_ids_script = prepare_vcf_and_add_ids_script,
        merge_vcfs_script = merge_vcfs_script,
        cargo_toml = cargo_toml,
        frac_missing = frac_missing,
        output_prefix = output_prefix
    }

    # make sure dict in header
    # TODO add preprocessing steps from KAGE Panel WDL

    output {
        File phased_collisionless_vcf = FixVariantCollisions.phased_collisionless_vcf
        File phased_collisionless_vcf_idx = FixVariantCollisions.phased_collisionless_vcf_idx
        File panel_vcf = PanGeniePanelCreation.panel_vcf
        File panel_vcf_idx = PanGeniePanelCreation.panel_vcf_idx
        File panel_id_split_vcf = PanGeniePanelCreation.panel_id_split_vcf
        File panel_id_split_vcf_idx = PanGeniePanelCreation.panel_id_split_vcf_idx
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

task FixVariantCollisions {
    input {
        File phased_vcf                     # biallelic
        File fix_variant_collisions_script
        Int operation = 1                   # 0=can only remove an entire VCF record; 1=can remove single ones from a GT
        String weight_tag = "AF"            # ID of the weight field; weights are assumed to be non-negative; TODO we set to AF for now, perhaps should annotate SVLEN or something else that would prefer SVs
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
        bcftools view ~{phased_vcf} --threads $(nproc) | \
        ./FixVariantCollisions \
            ~{operation} \
            ~{weight_tag} \
            ~{is_weight_format_field} \
            ~{default_weight} \
            histogram.txt | \
        bcftools +setGT --no-version -Ou -- -t . -n 0p | \
            bcftools +fill-tags --no-version --threads $(nproc) --write-index=csi -Ob -o ~{output_prefix}.phased.collisionless.bcf -- -t AF,AC,AN
    >>>

    output {
        File phased_collisionless_vcf = "~{output_prefix}.phased.collisionless.bcf"
        File phased_collisionless_vcf_idx = "~{output_prefix}.phased.collisionless.bcf.csi"
        File histogram = "histogram.txt"
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

task PanGeniePanelCreation {
    input {
        File phased_vcf
        File phased_vcf_idx
        File reference_fasta
        String region
        String output_prefix

        File prepare_vcf_and_add_ids_script
        File merge_vcfs_script
        File cargo_toml
        Float frac_missing

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 3 * ceil(size(phased_vcf, "GiB")) + ceil(size(reference_fasta, "GiB")) 

    command <<<
        set -euxo pipefail

        mkdir -p pangenie-utils/src/bin
        cp ~{prepare_vcf_and_add_ids_script} pangenie-utils/src/bin
        cp ~{merge_vcfs_script} pangenie-utils/src/bin
        cp ~{cargo_toml} pangenie-utils
        cd pangenie-utils
        cargo build --release
        cd ..

        bcftools stats -r ~{region} --regions-overlap 0 ~{phased_vcf} > ~{output_prefix}.stats.txt
        bcftools view --no-version -h ~{phased_vcf} > header.txt

        # validate variants against reference, run PanGenie prepare-vcf and add-ids scripts, split to biallelic, and run PanGenie merge script;
        # everything should be normalized or in the desired representation at this point
        bcftools norm --no-version -r ~{region} --regions-overlap 0 --do-not-normalize --check-ref e --fasta-ref ~{reference_fasta} ~{phased_vcf} | \
            ./pangenie-utils/target/release/prepare_vcf_and_add_ids --missing ~{frac_missing} | \
            bcftools norm --no-version -m-any --do-not-normalize | tee \
        >(  bcftools view --no-version --write-index=csi -Ob -o ~{output_prefix}.prepare.id.split.bcf ) | \
         (  ./pangenie-utils/target/release/merge_vcfs merge \
                --header header.txt \
                -r ~{reference_fasta} \
                --ploidy 2 | \
            bcftools view --no-version --write-index=csi -Ob -o ~{output_prefix}.prepare.id.split.mergehap.bcf )

        bcftools stats ~{output_prefix}.prepare.id.split.mergehap.bcf > ~{output_prefix}.prepare.id.split.mergehap.stats.txt
    >>>

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             3,
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

    output {
        File input_stats = "~{output_prefix}.stats.txt"
        File panel_stats = "~{output_prefix}.prepare.id.split.mergehap.stats.txt"
        File panel_vcf = "~{output_prefix}.prepare.id.split.mergehap.bcf"
        File panel_vcf_idx = "~{output_prefix}.prepare.id.split.mergehap.bcf.csi"
        File panel_id_split_vcf = "~{output_prefix}.prepare.id.split.bcf"
        File panel_id_split_vcf_idx = "~{output_prefix}.prepare.id.split.bcf.csi"
    }
}
