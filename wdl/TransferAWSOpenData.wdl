version 1.0

workflow TransferAWSOpenData {
    input {
        String s3_uri
    }

    call DownloadPublicS3 {
        input:
            s3_uri = s3_uri
    }

    output {
        File transferred_file = DownloadPublicS3.downloaded_file
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

task DownloadPublicS3 {
    input {
        String s3_uri

        RuntimeAttr? runtime_attr_override
    }

    String filename = basename(s3_uri)

    command <<<
        set -e
        
        echo "Downloading ~{s3_uri} from AWS Open Data..."
        
        # Use --no-sign-request to bypass the need for AWS credentials
        aws s3 cp --no-sign-request "~{s3_uri}" "./~{filename}"
        
        echo "Download complete."
    >>>

    output {
        File downloaded_file = "~{filename}"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            500,
        boot_disk_gb:       10,
        use_ssd:            true,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "amazon/aws-cli:latest"
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
