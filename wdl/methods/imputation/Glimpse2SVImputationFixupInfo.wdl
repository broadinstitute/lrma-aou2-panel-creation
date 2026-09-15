version 1.0

workflow Glimpse2SVImputationFixupInfo {
    input {
        Array[File] original_imputed_vcfs
        Array[File] original_imputed_vcf_idxs
        File panel_popped_sites_vcf
        File panel_popped_sites_vcf_idx
        Array[Int] batch_sizes
        Array[String] chromosomes
        Array[Array[String]] regions
        
        String output_basename
        String rust_docker
    }

    scatter (contig_idx in range(length(chromosomes))) {
        File contig_vcf = original_imputed_vcfs[contig_idx]
        File contig_vcf_idx = original_imputed_vcf_idxs[contig_idx]
        String chromosome = chromosomes[contig_idx]
        Array[String] contig_regions = regions[contig_idx]

        scatter (shard_idx in range(length(contig_regions))) {
            String region = contig_regions[shard_idx]

            call FixupInfo {
                input:
                    original_vcf = contig_vcf,
                    original_vcf_idx = contig_vcf_idx,
                    panel_popped_sites_vcf = panel_popped_sites_vcf,
                    panel_popped_sites_vcf_idx = panel_popped_sites_vcf_idx,
                    region = region,
                    batch_sizes = batch_sizes,
                    output_prefix = output_basename + "." + chromosome + ".shard_" + shard_idx,
                    docker = rust_docker
            }
        }

        call ConcatVCFs {
            input:
                vcfs = FixupInfo.fixed_vcf,
                vcf_idxs = FixupInfo.fixed_vcf_idx,
                output_name = output_basename + "." + chromosome,
                docker = rust_docker
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
        File panel_popped_sites_vcf
        File panel_popped_sites_vcf_idx
        String region
        Array[Int] batch_sizes
        String output_prefix
        File fixup_info_bin
        
        String docker
        Int mem_gb = 8
        Int cpu = 2
        Int disk_gb = ceil(3 * size(original_vcf, "GiB") + size(panel_popped_sites_vcf, "GiB") + 20)
    }

    command <<<
        set -euox pipefail

        cp ~{fixup_info_bin} fixup-info
        chmod +x fixup-info

        bcftools view -r ~{region} --regions-overlap 0 ~{original_vcf} | \
            ./fixup-info ~{panel_popped_sites_vcf} ~{sep="," batch_sizes} ~{output_prefix}.batched_info.tsv.gz | \
            bcftools view -W=tbi -Oz -o ~{output_prefix}.vcf.gz
    >>>

    runtime {
        docker: docker
        memory: mem_gb + " GiB"
        cpu: cpu
        disks: "local-disk " + disk_gb + " HDD"
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
        disks: "local-disk " + disk_gb + " HDD"
        preemptible: 3
    }

    output {
        File concat_vcf = "~{output_name}.vcf.gz"
        File concat_vcf_idx = "~{output_name}.vcf.gz.tbi"
    }
}
