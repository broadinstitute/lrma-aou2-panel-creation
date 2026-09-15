version 1.0

workflow Glimpse2SVImputationFixupInfo {
    input {
        Array[File] original_imputed_vcfs
        Array[File] original_imputed_vcf_idxs
        Array[File] panel_popped_sites_only_vcf_gzs
        Array[File] panel_popped_sites_only_vcf_gz_tbis
        Array[Int] batch_sizes
        Array[String] chromosomes
        Array[Array[String]] regions
        
        String output_basename
        String docker
    }

    scatter (contig_idx in range(length(chromosomes))) {
        Array[String] contig_regions = regions[contig_idx]

        scatter (shard_idx in range(length(contig_regions))) {
            call FixupInfo {
                input:
                    original_vcf = original_imputed_vcfs[contig_idx],
                    original_vcf_idx = original_imputed_vcf_idxs[contig_idx],
                    panel_popped_sites_only_vcf_gz = panel_popped_sites_only_vcf_gzs[contig_idx],
                    panel_popped_sites_only_vcf_gz_tbi = panel_popped_sites_only_vcf_gz_tbis[contig_idx],
                    region = contig_regions[shard_idx],
                    batch_sizes = batch_sizes,
                    output_prefix = output_basename + "." + chromosomes[contig_idx] + ".shard_" + shard_idx,
                    docker = docker
            }
        }

        call ConcatVCFs {
            input:
                vcfs = FixupInfo.fixed_vcf,
                vcf_idxs = FixupInfo.fixed_vcf_idx,
                output_name = output_basename + "." + chromosomes[contig_idx],
                docker = docker
        }
    }

    output {
        Array[File] fixed_vcfs = ConcatVCFs.concat_vcf
        Array[File] fixed_vcf_idxs = ConcatVCFs.concat_vcf_idx
        Array[Array[File]] batched_infos = FixupInfo.batched_info_tsv
    }
}

task FixupInfo {
    input {
        File original_vcf
        File original_vcf_idx
        File panel_popped_sites_only_vcf_gz
        File panel_popped_sites_only_vcf_gz_tbi
        String region
        Array[Int] batch_sizes
        String output_prefix
        File fixup_info_bin
        
        String docker
        Int mem_gb = 8
        Int cpu = 2
        Int disk_gb = ceil(3 * size(original_vcf, "GiB") + size(panel_popped_sites_only_vcf_gz, "GiB") + 20)
    }

    command <<<
        set -euox pipefail

        cp ~{fixup_info_bin} fixup-info
        chmod +x fixup-info

        bcftools view -r ~{region} --regions-overlap 0 ~{original_vcf} | \
            ./fixup-info ~{panel_popped_sites_only_vcf_gz} ~{sep="," batch_sizes} ~{output_prefix}.batched_info.tsv.gz | \
            bcftools view -W=tbi -Oz -o ~{output_prefix}.vcf.gz
    >>>

    runtime {
        docker: docker
        memory: mem_gb + " GiB"
        cpu: cpu
        disks: "local-disk " + disk_gb + " SSD"
        preemptible: 3
    }

    output {
        File fixed_vcf = "~{output_prefix}.vcf.gz"
        File fixed_vcf_idx = "~{output_prefix}.vcf.gz.tbi"
        File batched_info_tsv = "~{output_prefix}.batched_info.tsv.gz"
    }
}

task ConcatVCFs {
    input {
        Array[File] vcfs
        Array[File] vcf_idxs
        String output_name
        
        String docker
        Int mem_gb = 4
        Int cpu = 2
    }

    # Naive concatenation doesn't decompress, so footprint is small
    Int disk_gb = ceil(2 * size(vcfs, "GiB") + 20)

    command <<<
        set -euox pipefail
        
        # Write array to file to bypass maximum command-line length limits 
        # in case there are hundreds of shard inputs
        cat <<EOF > vcf_list.txt
        ~{sep='\n' vcfs}
        EOF

        bcftools concat \
            --naive \
            --file-list vcf_list.txt \
            -Oz \
            -o ~{output_name}.vcf.gz
        bcftools index -t ~{output_name}.vcf.gz
    >>>

    runtime {
        docker: docker
        memory: mem_gb + " GiB"
        cpu: cpu
        disks: "local-disk " + disk_gb + " SSD"
        preemptible: 3
    }

    output {
        File concat_vcf = "~{output_name}.vcf.gz"
        File concat_vcf_idx = "~{output_name}.vcf.gz.tbi"
    }
}
