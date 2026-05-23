version 1.0

# Remove colliding variants on phased haplotypes based on specified weights.

workflow FixVariantCollisions {

    input {
        File phased_vcf                     # biallelic, can be locally or fully phased
        File fix_variant_collisions_script
        String output_prefix

        Int operation = 1
        String weight_tag = "SCORE"
        Int is_weight_format_field = 0
        Float default_weight = 0.05
        String fix_variant_collisions_extra_args = "--use-gq --use-af --verbosity 1"
    }

    call FixVariantCollisions { input:
        phased_vcf = phased_vcf,
        fix_variant_collisions_script = fix_variant_collisions_script,
        output_prefix = output_prefix + ".collisionless",
        operation = operation,
        weight_tag = weight_tag,
        is_weight_format_field = is_weight_format_field,
        default_weight = default_weight,
        extra_args = fix_variant_collisions_extra_args
    }

    output {
        File collisionless_vcf = FixVariantCollisions.collisionless_vcf
        File collisionless_vcf_idx = FixVariantCollisions.collisionless_vcf_idx
        File collisionless_removed_counts_tsv = FixVariantCollisions.collisionless_removed_counts_tsv
        File collisionless_histogram_tsv = FixVariantCollisions.collisionless_histogram_tsv
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    String? disk_type
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

task FixVariantCollisions {
    input {
        File phased_vcf                          # biallelic, can be locally or fully phased
        File fix_variant_collisions_script
        Int operation = 1                        # 0=can only remove an entire VCF record; 1=can remove single ones from a GT
        String weight_tag = "SCORE"              # ID of the weight field; weights are assumed to be non-negative; we set to SCORE to prefer kanpig records (and moreover, those with higher SCORE) over DeepVariant records (these should have no SCORE, and will be assigned the low default_weight below)
        Int is_weight_format_field = 0           # given a VCF record in a sample, assign it a weight encoded in the INFO field (0) or in the sample column (1)
        Float default_weight = 0.05              # default weight if the weight field is not found
        String extra_args = "--use-gq --use-af --verbosity 1"  # use GQ then AF as tie breakers
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * (ceil(size(phased_vcf, "GiB")))

    command <<<
        set -euxo pipefail

        rustc -O ~{fix_variant_collisions_script} -o FixVariantCollisions

        # after FixVariantCollisions, replace all missing alleles (correctly) emitted with reference alleles, since this is expected by PanGenie panel-creation script
        bcftools view ~{phased_vcf} | \
        ./FixVariantCollisions \
            ~{operation} \
            ~{weight_tag} \
            ~{is_weight_format_field} \
            ~{default_weight} \
            ~{extra_args} \
            --removed-counts-tsv ~{output_prefix}.removed.tsv \
            --histogram-tsv ~{output_prefix}.histogram.tsv | \
        bcftools +setGT --no-version -Ou -- -t . -n 0p | \
            bcftools +fill-tags --no-version --write-index=csi -Ob -o ~{output_prefix}.bcf -- -t AF,AC,AN
    >>>

    output {
        File collisionless_vcf = "~{output_prefix}.bcf"
        File collisionless_vcf_idx = "~{output_prefix}.bcf.csi"
        File collisionless_removed_counts_tsv = "~{output_prefix}.removed.tsv"
        File collisionless_histogram_tsv = "~{output_prefix}.histogram.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             6,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "SSD",
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsde-methods/slee/lrma-aou2-panel-creation-rust:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type])
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
