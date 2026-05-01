version 1.0

# Over a given region, create non-overlapping shards of a given genomic size for subsequent 
# 1) subsetting, filtering, short + SV concatenation, and collision removal
# 2) bubble creation

# TODO use numbers of variants instead
# TODO check no bp overlaps
# TODO use logic in https://app.terra.bio/#workspaces/allofus-drc-wgs-LR-prodPaper/AoU_DRC_WGS_LongReads_Imputation_Phase_2/analysis/launch/find-sv-windows-bedtools.ipynb
#   we want to create genomic shards that:
#   1) are above a minimum size, and
#   2) have boundaries that are well away from any SVs
#   this is so we can shard FixVariantCollisions/PanGenieBubbleCreation without worrying about splitting potential bubbles across shards
# TODO scatter over all chromosome regions and create a Terra data table with a row for each shard region

workflow CreateShards {

    input {
        String region
        String output_prefix

        Int shard_size = 2000000
    }

    call CreateShards { input:
        region = region,
        shard_size = shard_size,
        pad_size = 0,
        output_prefix = output_prefix + ".shards"
    }

    output {
        Array[String] shard_regions = CreateShards.shard_regions
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

task CreateShards {
    input {
        String region
        Int shard_size
        Int pad_size
        String output_prefix

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        python - --region ~{region} \
                 --shard_size ~{shard_size} \
                 --pad_size ~{pad_size} \
                 --output_file ~{output_prefix} \
                 <<-'EOF'
        import argparse

        def split_region(region):
            chromosome, span = region.split(":")
            start, end = span.split("-")
            return(chromosome, int(start), int(end))

        def split_region_to_intervals(region, shard_size, pad_size):
            chromo, start, end = split_region(region)
            shard_num = (end - start)//shard_size
            intervals = [(chromo, start, start + shard_size + pad_size)]
            for i in range(1, shard_num):
                start_pos = start + i*shard_size - pad_size
                end_pos = start_pos + shard_size + pad_size
                intervals.append((chromo, start_pos, end_pos))
            if end > start + shard_num*shard_size:
                intervals.append((chromo, start + shard_num*shard_size - pad_size, end))
            return(intervals)

        def write_output_file(content, output_file):
            with open(output_file, "w") as f:
                for item in content:
                    l = "%s:%d-%d" % (item[0], item[1], item[2])
                    f.write(l+ "\n")

        def main():
            parser = argparse.ArgumentParser()

            parser.add_argument('--region',
                                type=str)

            parser.add_argument('--output_file',
                                type=str)

            parser.add_argument('--shard_size',
                    type=int)

            parser.add_argument('--pad_size',
                    type=int)

            args = parser.parse_args()

            intervals = split_region_to_intervals(args.region, args.shard_size, args.pad_size)
            write_output_file(intervals, args.output_file + ".txt")

        if __name__ == "__main__":
            main()
        EOF
    >>>

    output {
        Array[String] shard_regions = read_lines("~{output_prefix}.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            10,
        boot_disk_gb:       10,
        use_ssd:            false,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.11"
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
