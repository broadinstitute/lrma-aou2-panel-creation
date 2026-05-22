version 1.0

workflow ConcatVcfs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String output_prefix
        Boolean do_sort = false
        String extra_args = "--threads $(nproc) --naive"
    }

    call ConcatVcfs {
        input:
            vcfs = vcfs,
            vcf_idxs = vcf_idxs,
            output_prefix = output_prefix,
            do_sort = do_sort,
            extra_args = extra_args
    }

    output {
        File concatenated_vcf = ConcatVcfs.concatenated_vcf
        File concatenated_vcf_idx = ConcatVcfs.concatenated_vcf_idx
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

task ConcatVcfs {
    input{
        Array[File] vcfs
        Array[File] vcf_idxs
        String output_prefix
        Boolean do_sort = false
        String? extra_args

        RuntimeAttr? runtime_attr_override
    }

    # If sorting, provide extra disk space for the temporary sort shards
    Int disk_gb = if do_sort then 10 + 4 * ceil(size(vcfs, "GiB")) else 10 + 2 * ceil(size(vcfs, "GiB"))

    command <<<
        set -euox pipefail

        # Start zero-overhead background heartbeat monitor
        (
            echo "Starting concat monitoring..." >&2
            while true; do
                if [ -f "~{output_prefix}.bcf" ]; then
                    SIZE=$(ls -lh "~{output_prefix}.bcf" | awk '{print $5}')
                    echo "[Heartbeat] ~{output_prefix}.bcf is currently $SIZE..." >&2
                fi
                sleep 60
            done
        ) &
        HEARTBEAT_PID=$!

        if [ "~{do_sort}" == "true" ]; then
            echo "Concatenating and piping to bcftools sort..."
            # Output as uncompressed BCF (-Ou) to bypass intermediate compression overhead
            bcftools concat \
                -f ~{write_lines(vcfs)} \
                ~{extra_args} \
                -Ou | bcftools sort -m 2G -Ob -o ~{output_prefix}.bcf
        else
            echo "Concatenating directly to disk..."
            bcftools concat \
                -f ~{write_lines(vcfs)} \
                ~{extra_args} \
                -Ob -o ~{output_prefix}.bcf
        fi

        bcftools index ~{output_prefix}.bcf

        # Kill the background monitor the second the pipeline finishes
        kill $HEARTBEAT_PID || true
    >>>

    output {
        File concatenated_vcf = "~{output_prefix}.bcf"
        File concatenated_vcf_idx = "~{output_prefix}.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "SSD",
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.23"
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
