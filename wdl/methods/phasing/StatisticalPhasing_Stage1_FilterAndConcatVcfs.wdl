version 1.0

# Subset and filter short and SV VCFs separately, then concatenate.
# Run over a Terra data table of non-overlapping region shards containing roughly equal numbers of variants.

workflow FilterAndConcatVcfs {

    input {
        File short_vcf        # multiallelic
        File short_vcf_idx
        File sv_vcf           # biallelic
        File sv_vcf_idx
        File reference_fasta
        File reference_fasta_fai
        String region
        String output_prefix

        String? short_filter_args
        String short_view_args = "-i 'MAC>=2 && abs(strlen(ALT)-strlen(REF))<50'"
        String sv_view_args = "-i 'MAC>=2 && abs(strlen(ALT)-strlen(REF))>=50'"
    }

    call SubsetVCFStreaming as SubsetVcfShort { input:
        vcf = short_vcf,
        vcf_idx = short_vcf_idx,
        region = region,
        output_prefix = output_prefix + ".subsetShort"
    }

    call SubsetVCF as SubsetVcfSV { input:
        vcf = select_first([sv_vcf ]),
        vcf_idx = select_first([sv_vcf_idx]),
        region = region,
        output_prefix = output_prefix + ".subsetSV"
    }

    call FilterShortVcf { input:
        short_vcf = SubsetVcfShort.subset_vcf,
        short_vcf_idx = SubsetVcfShort.subset_idx,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        output_prefix = output_prefix + ".filterShort",
        short_filter_args = short_filter_args,
        short_view_args = short_view_args
    }

    call FilterSVVcf { input:
        sv_vcf = SubsetVcfSV.subset_vcf,
        sv_vcf_idx = SubsetVcfSV.subset_idx,
        output_prefix = output_prefix + ".filterSV",
        sv_view_args = sv_view_args
    }

    call ConcatShortAndSVVcfs { input:
        short_vcf = FilterShortVcf.filter_short_vcf,
        short_vcf_idx = FilterShortVcf.filter_short_vcf_idx,
        sv_vcf = FilterSVVcf.filter_sv_vcf,
        sv_vcf_idx = FilterSVVcf.filter_sv_vcf_idx,
        output_prefix = output_prefix + ".filterAndConcat"
    }

    output {
        File filter_and_concat_vcf = ConcatShortAndSVVcfs.concat_vcf
        File filter_and_concat_vcf_idx = ConcatShortAndSVVcfs.concat_vcf_idx
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

task SubsetVCF {
    input {
        File vcf
        File vcf_idx
        String region
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size([vcf, vcf_idx], "GiB"))

    command <<<
        set -euxo pipefail

        bcftools view --no-version ~{vcf}##idx##~{vcf_idx} \
            --regions ~{region} \
            --regions-overlap 0 \
            --write-index=csi -Ob -o ~{output_prefix}.bcf
    >>>

    output {
        File subset_vcf = "~{output_prefix}.bcf"
        File subset_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
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

task SubsetVCFStreaming {
    parameter_meta {
        vcf: {
            localization_optional: true
        }
        vcf_idx: {
            localization_optional: true
        }
    }

    input {
        File vcf
        File vcf_idx
        String region
        String output_prefix

        Int disk_gb = 10

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        bcftools view --no-version ~{vcf} \
            --regions ~{region} \
            --regions-overlap 0 \
            --write-index=csi -Ob -o ~{output_prefix}.bcf
    >>>

    output {
        File subset_vcf = "~{output_prefix}.bcf"
        File subset_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             3,
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

task FilterShortVcf {
    input {
        File short_vcf         # multiallelic
        File short_vcf_idx
        File reference_fasta
        File reference_fasta_fai
        String output_prefix

        String? short_filter_args
        String short_view_args = "-i 'MAC>=2 && abs(strlen(ALT)-strlen(REF))<50'"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size(short_vcf, "GiB"))

    command <<<
        set -euxo pipefail

        # normalize then split to biallelic, then filter short (re-fill tags when needed)
        # order of operations matters for final variant representations, only change with caution!
        bcftools norm --no-version -m-any -f ~{reference_fasta} ~{short_vcf} -Ou | \
            bcftools +fill-tags --no-version -Ou -- -t AF,AC,AN | \
            bcftools filter --no-version ~{short_filter_args} -Ou | \
            bcftools +fill-tags --no-version -Ou -- -t AF,AC,AN | \
            bcftools view --no-version ~{short_view_args} \
                --write-index=csi -Ob -o ~{output_prefix}.short.bcf
    >>>

    output {
        File filter_short_vcf = "~{output_prefix}.short.bcf"
        File filter_short_vcf_idx = "~{output_prefix}.short.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             6,
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

task FilterSVVcf {
    input {
        File sv_vcf            # biallelic
        File sv_vcf_idx
        String output_prefix

        String sv_view_args = "-i 'MAC>=2 && abs(strlen(ALT)-strlen(REF))>=50'"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size(sv_vcf, "GiB"))

    command <<<
        set -euxo pipefail

        # populate missing with hom-ref and filter SV
        bcftools +setGT ~{sv_vcf} --no-version -Ou -- -t . -n 0p | \
            bcftools +fill-tags --no-version -Ou -- -t AF,AC,AN | \
            bcftools view --no-version ~{sv_view_args} \
                --write-index=csi -Ob -o ~{output_prefix}.SV.bcf
    >>>

    output {
        File filter_sv_vcf = "~{output_prefix}.SV.bcf"
        File filter_sv_vcf_idx = "~{output_prefix}.SV.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             6,
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

task ConcatShortAndSVVcfs {
    input {
        File short_vcf
        File short_vcf_idx
        File sv_vcf
        File sv_vcf_idx
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size(short_vcf, "GiB") + size(sv_vcf, "GiB"))

    command <<<
        set -euxo pipefail

        # concatenate with deduplication; providing SV VCF as first argument preferentially keeps those records
        bcftools concat --no-version \
            ~{sv_vcf} \
            ~{short_vcf} \
            --allow-overlaps --remove-duplicates -Ou | \
            bcftools sort --write-index=csi -Ob -o ~{output_prefix}.bcf
    >>>

    output {
        File concat_vcf = "~{output_prefix}.bcf"
        File concat_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             6,
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
