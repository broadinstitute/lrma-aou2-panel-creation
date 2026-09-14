version 1.0

# Read-only probes, change nothing:
#   Probe:     which bucket-level permissions the task VM's service account holds on the
#              workspace bucket (testIamPermissions), plus the bucket's soft-delete policy.
#   ProbeBoot: same task with bootDiskSizeGb set, so the resulting Batch job spec and VM disks
#              show whether the WDL controls the boot disk on this backend.

workflow ProbeBucketPerms {
    input {
        String bucket_name = "longreadsphase2imputation"
        String docker = "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.23"
    }
    call Probe { input: bucket_name = bucket_name, docker = docker }
    call ProbeBoot { input: bucket_name = bucket_name, docker = docker }
    output {
        File report = Probe.report
        File boot_report = ProbeBoot.report
    }
}

task Probe {
    input {
        String bucket_name
        String docker
    }
    command <<<
        set -euo pipefail
        { gcloud auth list 2>&1; df -h / /mnt/disks/cromwell_root 2>&1; lsblk 2>&1 || true; } | tee report.txt
        TOK=$(gcloud auth print-access-token)
        curl -s -H "Authorization: Bearer $TOK" \
          "https://storage.googleapis.com/storage/v1/b/~{bucket_name}/iam/testPermissions?permissions=storage.buckets.update&permissions=storage.buckets.get&permissions=storage.objects.delete" \
          | tee -a report.txt
        echo >> report.txt
        gcloud storage buckets describe "gs://~{bucket_name}" --format="yaml(soft_delete_policy,lifecycle)" 2>&1 | tee -a report.txt
    >>>
    output { File report = "report.txt" }
    runtime {
        docker: docker
        cpu: 1
        memory: "2 GiB"
        disks: "local-disk 10 HDD"
        preemptible: 1
        maxRetries: 1
    }
}

task ProbeBoot {
    input {
        String bucket_name
        String docker
    }
    command <<<
        set -euo pipefail
        { echo "bootDiskSizeGb 5 requested"; df -h / /mnt/disks/cromwell_root 2>&1; lsblk 2>&1 || true; } | tee report.txt
    >>>
    output { File report = "report.txt" }
    runtime {
        docker: docker
        cpu: 1
        memory: "2 GiB"
        disks: "local-disk 10 HDD"
        bootDiskSizeGb: 5
        preemptible: 1
        maxRetries: 1
    }
}
