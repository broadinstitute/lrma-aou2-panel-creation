version 1.0

workflow PhysicalPhasing {

    input {
        String sample_name
        String? region

        Array[String] short_vcf_paths_per_chromosome = [
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr1/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr2/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr3/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr4/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr5/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr6/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr7/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr8/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr9/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr10/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr11/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr12/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr13/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr14/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr15/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr16/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr17/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr18/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr19/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr20/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr21/shard-*/" + sample_name + ".short.shard-*.bcf",
            "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_short/test_chr22/shard-*/" + sample_name + ".short.shard-*.bcf"]

        String sv_vcf_shard_paths = "gs://fc-secure-e9018a40-98de-4d5d-8f40-16e60c8f8a0b/prod/split_sv/v3_main/shard-*/" + sample_name + ".sv.shard-*.bcf"

        File? trgt_vcf                  # untested
        File? trgt_vcf_idx

        String short_view_args = "-e 'QUAL<20 || abs(ILEN)>=20 || (FILTER!=\"PASS\" && FILTER!=\".\")'"
        String short_filter_args = "-S . -e 'GT=\"alt\" && ((TYPE=\"snp\" && GQ<15) || (TYPE!=\"snp\" && GQ<5))'"

        File bam
        File bam_idx
        File reference_fasta            # note that a preprocessing step that converted REF/ALT bases in all input VCFs from lower case to upper case (necessary for CHM13) was removed
        File reference_fasta_fai

        Boolean do_haplotagging = false
        String hiphase_extra_args = "--threads $(nproc) --global-realignment-cputime 300"
    }

    call PreProcessVCFs { input:
        sample_name = sample_name,
        region = region,
        short_vcf_paths_per_chromosome = short_vcf_paths_per_chromosome,
        sv_vcf_shard_paths = sv_vcf_shard_paths,
        trgt_vcf = trgt_vcf,
        trgt_vcf_idx = trgt_vcf_idx,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        short_view_args = short_view_args,
        short_filter_args = short_filter_args
    }

    call HiPhase { input:
        sample_name = sample_name,
        short_vcf = PreProcessVCFs.preprocessed_short_vcf,
        short_vcf_idx = PreProcessVCFs.preprocessed_short_vcf_idx,
        sv_vcf = PreProcessVCFs.preprocessed_sv_vcf,
        sv_vcf_idx = PreProcessVCFs.preprocessed_sv_vcf_idx,
        trgt_vcf = PreProcessVCFs.preprocessed_trgt_vcf,
        trgt_vcf_idx = PreProcessVCFs.preprocessed_trgt_vcf_idx,
        bam = bam,
        bam_idx = bam_idx,
        reference_fasta = reference_fasta,
        reference_fasta_fai = reference_fasta_fai,
        do_haplotagging = do_haplotagging,
        extra_args = hiphase_extra_args
    }

    output {
        File hiphase_short_vcf = HiPhase.hiphase_short_vcf
        File hiphase_short_vcf_idx = HiPhase.hiphase_short_vcf_idx
        File hiphase_sv_vcf = HiPhase.hiphase_sv_vcf
        File hiphase_sv_vcf_idx = HiPhase.hiphase_sv_vcf_idx
        File? hiphase_trgt_vcf = HiPhase.hiphase_trgt_vcf
        File? hiphase_trgt_vcf_idx = HiPhase.hiphase_trgt_vcf_idx
        Array[File] hiphase_metrics_files = HiPhase.hiphase_metrics_files
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

# separately naive concat short and SV VCFs, subset all VCFs, norm and filter short VCFs
task PreProcessVCFs {
    input {
        String sample_name
        String? region
        Array[String] short_vcf_paths_per_chromosome
        String sv_vcf_shard_paths
        File? trgt_vcf
        File? trgt_vcf_idx

        File reference_fasta
        File reference_fasta_fai

        String short_view_args = "-e 'QUAL<20 || abs(ILEN)>=20 || (FILTER!=\"PASS\" && FILTER!=\".\")'"
        String short_filter_args = "-S . -e 'GT=\"alt\" && ((TYPE=\"snp\" && GQ<15) || (TYPE!=\"snp\" && GQ<5))'"

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        # localize, concat, subset, norm, and filter short
        # -- localize and concat by chrom
        chrom_index=1
        for chromosome_short_vcf_paths in ~{sep=" " short_vcf_paths_per_chromosome}; do
            padded_chrom_index=$(printf "%02d" "$chrom_index")
            gsutil -m cp $chromosome_short_vcf_paths .
            bcftools concat ~{sample_name}.short.shard-*.bcf \
                --naive \
                -Ob -o ~{sample_name}.short.chr-$padded_chrom_index.bcf
            rm ~{sample_name}.short.shard-*.bcf
            chrom_index=$((chrom_index + 1))
        done
        # -- concat across chroms
        bcftools concat ~{sample_name}.short.chr-*.bcf\
            --naive \
            -Ob -o ~{sample_name}.short.bcf
        bcftools index ~{sample_name}.short.bcf
        rm ~{sample_name}.short.chr-*.bcf

        # note that norm splitting is critical for SNV/non-SNV filters to be properly applied to mixed sites
        bcftools norm ~{sample_name}.short.bcf \
                ~{"-r " + region} \
                -f ~{reference_fasta} \
                -m-any -Ou | \
            bcftools view ~{short_view_args} -Ou | \
            bcftools filter ~{short_filter_args} -Ou | \
            bcftools sort -W=csi -Ob -o ~{sample_name}.preprocessed.short.bcf

        # localize, concat, and subset sv
        gsutil -m cp ~{sv_vcf_shard_paths} .
        bcftools concat ~{sample_name}.sv.shard-*.bcf\
            --naive \
            -Ob -o ~{sample_name}.sv.bcf
        bcftools index ~{sample_name}.sv.bcf
        rm ~{sample_name}.sv.shard-*.bcf
        bcftools view ~{sample_name}.sv.bcf \
             ~{"-r " + region} \
             -W=csi -Ob -o ~{sample_name}.preprocessed.sv.bcf

        # subset trgt
        if [ -n "~{trgt_vcf}" ]; then
            bcftools view ~{"-r " + region} ~{trgt_vcf}##idx##~{trgt_vcf_idx} -W=csi -Ob -o ~{sample_name}.preprocessed.trgt.bcf
        fi
    >>>

    output {
        File preprocessed_short_vcf = "~{sample_name}.preprocessed.short.bcf"
        File preprocessed_short_vcf_idx = "~{sample_name}.preprocessed.short.bcf.csi"
        File preprocessed_sv_vcf = "~{sample_name}.preprocessed.sv.bcf"
        File preprocessed_sv_vcf_idx = "~{sample_name}.preprocessed.sv.bcf.csi"
        File? preprocessed_trgt_vcf = "~{sample_name}.preprocessed.trgt.bcf"
        File? preprocessed_trgt_vcf_idx = "~{sample_name}.preprocessed.trgt.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            30,
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

task HiPhase {
    input {
        String sample_name
        File short_vcf
        File short_vcf_idx
        File sv_vcf
        File sv_vcf_idx
        File? trgt_vcf
        File? trgt_vcf_idx
        File bam
        File bam_idx
        File reference_fasta            # note that a preprocessing step that converted REF/ALT bases in all input VCFs from lower case to upper case (necessary for CHM13) was removed
        File reference_fasta_fai

        Boolean do_haplotagging = false
        String extra_args = "--threads $(nproc) --global-realignment-cputime 300"

        RuntimeAttr? runtime_attr_override
    }

    # TODO adjust based on do_haplotagging?
    Int disk_gb = 10 + 2 * ceil(size(select_all([bam, short_vcf, sv_vcf, trgt_vcf, reference_fasta]), "GB"))

    String haplotagging_args = if do_haplotagging then "--output-bam ~{sample_name}.hiphase.haplotagged.bam --haplotag-file ~{sample_name}.hiphase.haplotagged.tsv" else ""

    command <<<
        set -euxo pipefail

        # touch indices to avoid "older than" error messages
        # TODO: update htslib to a more recent version that doesn't throw these?
        touch ~{short_vcf_idx}
        touch ~{sv_vcf_idx}
        touch ~{bam_idx}
        touch ~{reference_fasta_fai}
        ~{if defined(trgt_vcf_idx) then "touch " + trgt_vcf_idx else ""}

        hiphase \
        --bam ~{bam} \
        --reference ~{reference_fasta} \
        --vcf ~{short_vcf} \
        --output-vcf ~{sample_name}.hiphase.short.bcf \
        --vcf ~{sv_vcf} \
        --output-vcf ~{sample_name}.hiphase.sv.bcf \
        --csi-index \
        --stats-file ~{sample_name}.hiphase.stats.csv \
        --blocks-file ~{sample_name}.hiphase.blocks.tsv \
        --summary-file ~{sample_name}.hiphase.summary.tsv \
        ~{haplotagging_args} \
        ~{if defined(trgt_vcf) then "--vcf " + trgt_vcf + " --output-vcf " + sample_name + ".hiphase.trgt.bcf" else ""} \
        ~{extra_args}
    >>>

    output {
        File hiphase_short_vcf = "~{sample_name}.hiphase.short.bcf"
        File hiphase_short_vcf_idx = "~{sample_name}.hiphase.short.bcf.csi"
        File hiphase_sv_vcf = "~{sample_name}.hiphase.sv.bcf"
        File hiphase_sv_vcf_idx = "~{sample_name}.hiphase.sv.bcf.csi"
        File? hiphase_trgt_vcf = "~{sample_name}.hiphase.trgt.bcf"
        File? hiphase_trgt_vcf_idx = "~{sample_name}.hiphase.trgt.bcf.csi"
        File? haplotagged_bam = "~{sample_name}.hiphase.haplotagged.bam"
        File? haplotagged_bam_idx = "~{sample_name}.hiphase.haplotagged.bam.csi"
        Array[File] hiphase_metrics_files = glob("~{sample_name}.hiphase.*.*sv")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             12,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/hiphase:v1.5.0"
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
