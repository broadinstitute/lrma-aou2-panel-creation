version 1.0

workflow PopGLIMPSE2 {
    input {
        Array[File] posteriors_vcf_gzs          # whole-genome, per-batch
        Array[File] posteriors_vcf_gz_tbis
        Array[File] panel_split_vcf_gzs         # per-chromosome
        Array[File] panel_split_vcf_gz_tbis
        Array[File] panel_id_split_vcf_gzs      # per-chromosome
        Array[File] panel_id_split_vcf_gz_tbis
        Array[String]+ chromosomes

        Array[String] output_prefixes

        File pop_python_script
    }

    scatter (i in range(length(posteriors_vcf_gzs)))
    {
        scatter (j in range(length(chromosomes))) {
            call PopGLIMPSE2 as ChromosomePopGLIMPSE2 {
                input:
                    posteriors_vcf_gz = posteriors_vcf_gzs[i],
                    posteriors_vcf_gz_tbi = posteriors_vcf_gz_tbis[i],
                    panel_split_vcf_gz = panel_split_vcf_gzs[j],
                    panel_split_vcf_gz_tbi = panel_split_vcf_gz_tbis[j],
                    panel_id_split_vcf_gz = panel_id_split_vcf_gzs[j],
                    panel_id_split_vcf_gz_tbi = panel_id_split_vcf_gz_tbis[j],
                    pop_python_script = pop_python_script,
                    chromosome = chromosomes[j],
                    output_prefix = output_prefixes[i] + "." + chromosomes[j]
            }
        }

        # concat across chromosomes
        call ConcatVcfs {
            input:
                vcf_gzs = ChromosomePopGLIMPSE2.popped_vcf_gz,
                vcf_gz_tbis = ChromosomePopGLIMPSE2.popped_vcf_gz_tbi,
                output_prefix = output_prefixes[i]
        }
    }

    output {
        Array[File] popped_vcf_gzs = ConcatVcfs.vcf_gz
        Array[File] popped_vcf_gz_tbis = ConcatVcfs.vcf_gz_tbi
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

task PopGLIMPSE2 {
    input {
        # all VCFs should be split to biallelic
        File posteriors_vcf_gz
        File posteriors_vcf_gz_tbi
        File panel_split_vcf_gz
        File panel_split_vcf_gz_tbi
        File panel_id_split_vcf_gz
        File panel_id_split_vcf_gz_tbi
        File pop_python_script
        String chromosome
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 3 * ceil(size(posteriors_vcf_gz, "GB"))

    command <<<
        set -euox pipefail

        pypy -m pip install tqdm

        bcftools annotate -r ~{chromosome} -a ~{panel_split_vcf_gz} ~{posteriors_vcf_gz} \
            -c CHROM,POS,REF,ALT,ID:=INFO/ID,INFO/ID:=INFO/ID \
            --write-index=tbi \
            -Oz -o ~{output_prefix}.annotated.vcf.gz

        # modified version of convert-to-biallelic.py
        # DO NOT apply bcftools norm -m+ before using this, pass a biallelic VCF instead!
        pypy ~{pop_python_script} \
            --panel_id_split_vcf_gz ~{panel_id_split_vcf_gz} \
            --input_vcf_gz ~{output_prefix}.annotated.vcf.gz | \
            bcftools view -Oz -o ~{output_prefix}.popped.vcf.gz
        bcftools index -t ~{output_prefix}.popped.vcf.gz
    >>>

    output {
        File popped_vcf_gz = "~{output_prefix}.popped.vcf.gz"
        File popped_vcf_gz_tbi = "~{output_prefix}.popped.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             6,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsde-methods/slee/pangenie-panel-creation:v1"
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
        Array[File] vcf_gzs
        Array[File] vcf_gz_tbis
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 3 * ceil(size(vcf_gzs, "GB"))

    command {
        set -euox pipefail

        mkdir inputs
        mv ~{sep=' ' vcf_gzs} inputs
        mv ~{sep=' ' vcf_gz_tbis} inputs

        # TODO NOTE USE OF SORT, ENSURE THIS GIVES DESIRED CHROMOSOME ORDER (COULD TAKE IN ARRAY OF CHROMOSOMES INSTEAD)
        if [ $(ls inputs/*.vcf.gz | wc -l) == 1 ]
        then
            cp $(ls inputs/*.vcf.gz) ~{output_prefix}.vcf.gz
            cp $(ls inputs/*.vcf.gz.tbi) ~{output_prefix}.vcf.gz.tbi
        else
            bcftools concat $(ls inputs/*.vcf.gz | sort -V -d) --naive -Oz -o ~{output_prefix}.vcf.gz
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
        mem_gb:             6,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsde-methods/slee/pangenie-panel-creation:v1"
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
