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
        Int mem_gb = 6
        Int cpu = 2
        Int disk_gb = ceil(2 * size(original_vcf, "GiB") + size(panel_popped_sites_only_vcf_gz, "GiB") + 20)
    }

    command <<<
        set -euox pipefail

        cp ~{fixup_info_bin} fixup-info
        chmod +x fixup-info

        # 1. Pre-shard the sites-only panel to drastically reduce Rust memory/startup time
        echo "Sharding panel VCF..."
        bcftools view -r ~{region} --regions-overlap 0 -Oz -o panel_shard.vcf.gz ~{panel_popped_sites_only_vcf_gz}
        bcftools index -t vcf panel_shard.vcf.gz

        # 2. Get expected site and sample counts
        IN_SITES=$(bcftools index -n panel_shard.vcf.gz)
        IN_SAMPLES=$(bcftools query -l ~{original_vcf} | wc -l)

        # 3. Execute processing pipeline (bgzip runs single-threaded to prevent CPU thrashing on 2 cores)
        echo "Running fixup-info..."
        bcftools view -r ~{region} --regions-overlap 0 ~{original_vcf} | \
            ./fixup-info panel_shard.vcf.gz ~{sep="," batch_sizes} ~{output_prefix}.batched_info.tsv.gz | \
            bgzip > ~{output_prefix}.vcf.gz
            
        # 4. Index the newly created bgzip output
        bcftools index -t ~{output_prefix}.vcf.gz
        
        # 5. Validation Check
        echo "Validating sites and samples..."
        OUT_SITES=$(bcftools index -n ~{output_prefix}.vcf.gz)
        OUT_SAMPLES=$(bcftools query -l ~{output_prefix}.vcf.gz | wc -l)
        
        if [ "$IN_SITES" -ne "$OUT_SITES" ]; then
            echo "FATAL ERROR: Site count mismatch. Panel region had $IN_SITES, but output has $OUT_SITES."
            exit 1
        fi
        
        if [ "$IN_SAMPLES" -ne "$OUT_SAMPLES" ]; then
            echo "FATAL ERROR: Sample count mismatch. Input had $IN_SAMPLES, but output has $OUT_SAMPLES."
            exit 1
        fi
        
        echo "Validation passed! $OUT_SITES sites and $OUT_SAMPLES samples successfully processed."
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
