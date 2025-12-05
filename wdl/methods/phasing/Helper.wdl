version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}


task HiPhase {

    meta {
        description: "Generates phased VCF. Note this runs fast so no need to parallize."
    }


    input {
        File bam
        File bai

        File snp_vcf_gz
        File snp__vcf_gz_tbi
        File sv__vcf_gz
        File sv__vcf_gz_tbi

        File ref_fasta
        File ref_fasta_fai
        String sample_name

        Int memory
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"

        String extra_args

        RuntimeAttr? runtime_attr_override
    }

    # Int bam_sz = ceil(size(bam, "GB"))
	Int disk_size = 30 # if bam_sz > 200 then 2*bam_sz else bam_sz + 200
    Int thread_num = memory/2

    command <<<
        set -euxo pipefail

        touch ~{bai}

        hiphase \
        --threads ~{thread_num} \
        --bam ~{bam} \
        --reference ~{ref_fasta} \
        --global-realignment-cputime 300 \
        --vcf ~{snp_vcf_gz} \
        --output-vcf ~{sample_name}_phased_snp.vcf.gz \
        --vcf ~{sv__vcf_gz} \
        --output-vcf ~{sample_name}_phased_sv.vcf.gz \
        --haplotag-file ~{sample_name}_phased_sv_haplotag.tsv \
        --stats-file ~{sample_name}.stats.csv \
        --blocks-file ~{sample_name}.blocks.tsv \
        --summary-file ~{sample_name}.summary.tsv \
        --verbose \
        ~{extra_args}

        bcftools sort ~{sample_name}_phased_snp.vcf.gz -O z -o ~{sample_name}_phased_snp.sorted.vcf.gz
        tabix -p vcf ~{sample_name}_phased_snp.sorted.vcf.gz

        bcftools sort ~{sample_name}_phased_sv.vcf.gz -O z -o ~{sample_name}_phased_sv.sorted.vcf.gz
        tabix -p vcf ~{sample_name}_phased_sv.sorted.vcf.gz
        
    >>>

    output {
        File phased_snp_vcf = "~{sample_name}_phased_snp.sorted.vcf.gz"
        File phased_snp_vcf_tbi = "~{sample_name}_phased_snp.sorted.vcf.gz.tbi"
        File phased_sv_vcf   = "~{sample_name}_phased_sv.sorted.vcf.gz"
        File phased_sv_vcf_tbi = "~{sample_name}_phased_sv.sorted.vcf.gz.tbi"
        File haplotag_file = "~{sample_name}_phased_sv_haplotag.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          thread_num,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "hangsuunc/hiphase:1.3.0"
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
        String extra_args

        RuntimeAttr? runtime_attr_override
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
    command <<<

        shapeit4.2 --input ~{vcf_input} \
                --map ~{mappingfile} \
                --region ~{region} \
                --sequencing \
                --output ~{output_prefix}.bcf \
                --thread $(nproc) \
                ~{extra_args}
        bcftools index ~{output_prefix}.bcf
    >>>

    output{
        # File resouce_monitor_log = "resources.log"
        File phased_bcf = "~{output_prefix}.bcf"
        File phased_bcf_index = "~{output_prefix}.bcf.csi"
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

task shapeit5_phase_common{
    input{
        File vcf_input
        File vcf_index
        File mappingfile
        String region
        String output_prefix
        Int cpu
        Int memory
        String extra_args
        Float minimal_maf = 0.01

        RuntimeAttr? runtime_attr_override
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
    command <<<
        # add AN AC tag
        bcftools +fill-tags ~{vcf_input} -Ob -o tmp.out.bcf -- -t AN,AC
        bcftools index tmp.out.bcf
        phase_common_static --input tmp.out.bcf \
                            --filter-maf ~{minimal_maf} \
                            --region ~{region} \
                            --map ~{mappingfile} \
                            --output scaffold.bcf \
                            --thread $(nproc) \
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

    Int disk_size = 2*ceil(size(vcfs, "GB")) + 1000

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


task shapeit5_phase_rare{
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
        String extra_args

        RuntimeAttr? runtime_attr_override
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
    }
    command <<<

        bcftools +fill-tags ~{vcf_input} -Ob -o tmp.out.bcf -- -t AN,AC
        bcftools index tmp.out.bcf

        phase_rare_static --input tmp.out.bcf \
                    --scaffold ~{scaffold_bcf} \
                    --map ~{mappingfile} \
                    --input-region ~{chunk_region} \
                    --scaffold-region ~{scaffold_region} \
                    --output ~{output_prefix}.chunk.~{chunknum}.bcf \
                    --thread $(nproc) \
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

# filter out singletons (i.e., keep MAC >= 2) and concatenate with deduplication
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

task split_into_shard {
    input {
        String locus
        Int bin_size
        Int pad_size
        String output_output_prefix

        Int? preemptible_tries
    }


    command <<<
        set -eo pipefail

        python - --locus ~{locus} \
                 --bin_size ~{bin_size} \
                 --pad_size ~{pad_size} \
                 --output_file ~{output_output_prefix} \
                 <<-'EOF'
        import gzip
        import argparse

        def split_locus(locus):
            chromosome, span = locus.split(":")
            start, end = span.split("-")
            return(chromosome, int(start), int(end))

        def split_locus_to_intervals(locus, bin_size=15000, pad_size=500):
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

        def write_bed_file(content, output_file):
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
            write_bed_file(intervals, args.output_file + ".txt")

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
        Array[String] locuslist = read_lines("~{output_output_prefix}.txt")
    }
}

task bcftools_concat_naive {
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
            --no-version \
            -Oz -o ~{output_prefix}.vcf.gz
        bcftools index -t ~{output_prefix}.vcf.gz
    >>>

    output {
        File concatenated_vcf = "~{output_prefix}.vcf.gz"
        File concatenated_vcf_tbi = "~{output_prefix}.vcf.gz.tbi"
    }

    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/lr-basic:0.1.1"
        memory: "32 GB"
        cpu: 8
        disks: "local-disk 1000 SSD"
    }
}
