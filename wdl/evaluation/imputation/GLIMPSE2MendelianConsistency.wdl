version 1.0

workflow MendelianConsistency {
    input {
        Array[File] panel_sites_only_vcfs
        Array[File] panel_sites_only_vcf_idxs
        Array[File] imputed_vcfs
        Array[File] imputed_vcf_idxs
        File trh_bed
        File trh_bed_idx
        File pedigree
        String output_prefix
    }

    scatter (idx in range(length(panel_sites_only_vcfs))) {
        call AnnotateVcf { input:
            panel_sites_only_vcf = panel_sites_only_vcfs[idx],
            panel_sites_only_vcf_idx = panel_sites_only_vcf_idxs[idx],
            imputed_vcf = imputed_vcfs[idx],
            imputed_vcf_idx = imputed_vcf_idxs[idx],
            trh_bed = trh_bed,
            trh_bed_idx = trh_bed_idx,
            output_prefix = output_prefix + "." + idx
        }

        call CalculateMendelianMetrics { input:
            annotated_vcf = AnnotateVcf.annotated_vcf,
            annotated_vcf_idx = AnnotateVcf.annotated_vcf_idx,
            pedigree = pedigree,
            output_prefix = output_prefix + "." + idx
        }

        # Generate Per-Chromosome Plots
        call PlotMendelianMetrics as PlotMetricsPerChrom { input:
            unfiltered_pkls = [CalculateMendelianMetrics.unfiltered_pkl],
            filtered_pkls = [CalculateMendelianMetrics.filtered_pkl],
            pedigree = pedigree,
            output_prefix = output_prefix + "." + idx
        }
    }

    # Generate Aggregate Plots
    call PlotMendelianMetrics as PlotMetricsAggregate { input:
        unfiltered_pkls = CalculateMendelianMetrics.unfiltered_pkl,
        filtered_pkls = CalculateMendelianMetrics.filtered_pkl,
        pedigree = pedigree,
        output_prefix = output_prefix + ".aggregate"
    }

    output {
        Array[File] mendelian_aggregate_plots_png = PlotMetricsAggregate.plots_png
        Array[File] mendelian_aggregate_plots_pdf = PlotMetricsAggregate.plots_pdf
        Array[File] mendelian_per_chrom_plots_png = flatten(PlotMetricsPerChrom.plots_png)
        Array[File] mendelian_per_chrom_plots_pdf = flatten(PlotMetricsPerChrom.plots_pdf)
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

task AnnotateVcf {
    input {
        File panel_sites_only_vcf
        File panel_sites_only_vcf_idx
        File imputed_vcf
        File imputed_vcf_idx
        File trh_bed
        File trh_bed_idx
        String output_prefix
        RuntimeAttr? runtime_attr_override
    }
    Int disk_gb = 20 + 3 * ceil(size(panel_sites_only_vcf, "GiB") + size(imputed_vcf, "GiB"))
    command <<<
        set -euxo pipefail
        bcftools annotate -a ~{trh_bed}##idx##~{trh_bed_idx} -c CHROM,FROM,TO -m +TRH --threads $(nproc) \
            ~{panel_sites_only_vcf}##idx##~{panel_sites_only_vcf_idx} -W -Ob -o panel.trh.bcf
        bcftools annotate -a panel.trh.bcf -c CHROM,POS,REF,ALT,INFO/AF,INFO/TRH --threads $(nproc) \
            ~{imputed_vcf}##idx##~{imputed_vcf_idx} -W -Ob -o ~{output_prefix}.annotated.bcf
    >>>
    output {
        File annotated_vcf = "~{output_prefix}.annotated.bcf"
        File annotated_vcf_idx = "~{output_prefix}.annotated.bcf.csi"
    }
    RuntimeAttr default_attr = object { cpu_cores: 2, mem_gb: 8, disk_gb: disk_gb, boot_disk_gb: 10, disk_type: "SSD", preemptible_tries: 1, max_retries: 0, docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.23" }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime { cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores]) memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB" disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type]) bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb]) preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries]) maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries]) docker: select_first([runtime_attr.docker, default_attr.docker]) }
}

task CalculateMendelianMetrics {
    input {
        File annotated_vcf
        File annotated_vcf_idx
        File pedigree
        String output_prefix
        Int chunk_size = 10000
        RuntimeAttr? runtime_attr_override
    }
    Int disk_gb = 20 + ceil(size(annotated_vcf, "GiB"))
    
    command <<<
        set -euxo pipefail
        conda install -y -c bioconda -c conda-forge bcftools scikit-allel pandas numpy pyarrow
        
        python - --input_path ~{annotated_vcf} --ped_path ~{pedigree} --output_prefix ~{output_prefix} --chunk_size ~{chunk_size} <<-'EOF'
        import argparse, numpy as np, pandas as pd, allel, sys, time, subprocess, warnings
        from datetime import timedelta

        def parse_pedigree(ped_path, vcf_samples):
            ped_df = pd.read_csv(ped_path, sep='\t', names=['familyID', 'sampleID', 'fatherID', 'motherID', 'sex', 'status'])
            samples_set = set(vcf_samples)
            valid_samples_set = set()
            valid_trios = []
            for _, row in ped_df.iterrows():
                if row['sampleID'] in samples_set and row['fatherID'] in samples_set and row['motherID'] in samples_set:
                    valid_trios.append({'sampleID': row['sampleID'], 'fatherID': row['fatherID'], 'motherID': row['motherID']})
                    valid_samples_set.update([row['sampleID'], row['fatherID'], row['motherID']])
            subset_samples = [s for s in vcf_samples if s in valid_samples_set]
            return valid_trios, subset_samples

        def process_vcf_in_chunks(input_path, valid_trios, subset_samples, chunk_size, num_trios):
            agg_results = {}
            fields_to_extract = ['calldata/GT', 'calldata/GP', 'variants/CHROM', 'variants/REF', 'variants/ALT', 'variants/AF', 'variants/TRH', 'variants/altlen', 'variants/is_snp']
            trio_samples_list = ",".join(subset_samples)
            bcftools_pipeline = f"bcftools view --threads 4 -Ov -s {trio_samples_list} {input_path}"
            proc = subprocess.Popen(bcftools_pipeline, shell=True, stdout=subprocess.PIPE)
            vcf_stream = proc.stdout
            fields, stream_samples, headers, it = allel.iter_vcf_chunks(vcf_stream, fields=fields_to_extract, chunk_length=chunk_size)
            
            stream_sample_to_idx = {s: i for i, s in enumerate(stream_samples)}
            proband_idx = [stream_sample_to_idx[t['sampleID']] for t in valid_trios]
            father_idx = [stream_sample_to_idx[t['fatherID']] for t in valid_trios]
            mother_idx = [stream_sample_to_idx[t['motherID']] for t in valid_trios]
            
            try:
                for chunk, actual_chunk_size, chrom_val, pos_val in it:
                    altlen_arr = chunk['variants/altlen']
                    length = altlen_arr[:, 0] if altlen_arr.ndim > 1 else altlen_arr
                    is_snp = chunk['variants/is_snp']
                    
                    len_bins = np.full(len(length), 'UNKNOWN', dtype=object)
                    len_bins[is_snp] = 'SNP'
                    len_bins[(50 <= length)] = '[50, inf)'
                    len_bins[(0 <= length) & (length < 50) & ~is_snp] = '[0, 50)'
                    len_bins[(-50 < length) & (length <= -1)] = '(-50, -1]'
                    len_bins[(length <= -50)] = '(-inf, -50]'
                    
                    afs = chunk['variants/AF']
                    afs = afs[:, 0] if afs.ndim > 1 else afs
                    af_bins = np.full(len(afs), 'UNKNOWN', dtype=object)
                    af_bins[(0 <= afs) & (afs < 0.01)] = '[0, 0.01)'
                    af_bins[(0.01 <= afs) & (afs < 0.1)] = '[0.01, 0.1)'
                    af_bins[(0.1 <= afs) & (afs <= 1.0)] = '[0.1, 1]'
                    
                    trhs = chunk['variants/TRH']
                    trhs = (trhs[:, 0] if trhs.ndim > 1 else trhs).astype(bool)

                    GT_base = chunk['calldata/GT']
                    if 'calldata/GP' in chunk:
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", category=RuntimeWarning)
                            max_gp = np.nanmax(chunk['calldata/GP'], axis=2) 
                    else:
                        max_gp = np.zeros(GT_base.shape[:2], dtype=np.float32) 

                    P = GT_base[:, proband_idx, :]
                    F = GT_base[:, father_idx, :]
                    M = GT_base[:, mother_idx, :]
                    P0, P1 = P[:, :, 0], P[:, :, 1]
                    F0, F1 = F[:, :, 0], F[:, :, 1]
                    M0, M1 = M[:, :, 0], M[:, :, 1]
                    
                    c1 = ((P0 == F0) | (P0 == F1)) & ((P1 == M0) | (P1 == M1))
                    c2 = ((P0 == M0) | (P0 == M1)) & ((P1 == F0) | (P1 == F1))
                    is_error_base = ~(c1 | c2)
                    
                    is_hom_ref_base = (P0 == 0) & (P1 == 0) & (F0 == 0) & (F1 == 0) & (M0 == 0) & (M1 == 0)
                    has_missing_base = (P0 == -1) | (P1 == -1) | (F0 == -1) | (F1 == -1) | (M0 == -1) | (M1 == -1)

                    for min_gp in [0.0, 0.9]:
                        if min_gp > 0:
                            p_mask = ~np.greater(max_gp[:, proband_idx], min_gp)
                            f_mask = ~np.greater(max_gp[:, father_idx], min_gp)
                            m_mask = ~np.greater(max_gp[:, mother_idx], min_gp)
                            has_missing = has_missing_base | p_mask | f_mask | m_mask
                        else:
                            has_missing = has_missing_base
                        
                        non_hom_ref = ~is_hom_ref_base & ~has_missing
                        errors = is_error_base & non_hom_ref 
                        
                        locus_errors = errors.sum(axis=1)
                        locus_nhr = non_hom_ref.sum(axis=1)
                        
                        chunk_df = pd.DataFrame({'AF_BIN': af_bins, 'LENGTH_BIN': len_bins, 'IN_TRH': trhs})
                        groups = chunk_df.groupby(['AF_BIN', 'LENGTH_BIN', 'IN_TRH'])
                        
                        for key, indices in groups.groups.items():
                            if 'UNKNOWN' in key: continue
                            full_key = (*key, min_gp)
                            if full_key not in agg_results:
                                agg_results[full_key] = {'errors_per_trio': np.zeros(num_trios, dtype=int), 'nhr_per_trio': np.zeros(num_trios, dtype=int), 'locus_rates': [], 'n_loci': 0}
                            agg_results[full_key]['errors_per_trio'] += errors[indices].sum(axis=0)
                            agg_results[full_key]['nhr_per_trio'] += non_hom_ref[indices].sum(axis=0)
                            
                            valid_loci = indices[locus_nhr[indices] > 0]
                            rates = locus_errors[valid_loci] / locus_nhr[valid_loci]
                            agg_results[full_key]['locus_rates'].extend(rates.tolist())
                            agg_results[full_key]['n_loci'] += len(indices)
            finally:
                if proc is not None:
                    proc.stdout.close()
                    proc.wait()
            return agg_results

        def main():
            parser = argparse.ArgumentParser()
            parser.add_argument("--input_path", required=True)
            parser.add_argument("--ped_path", required=True)
            parser.add_argument("--output_prefix", required=True)
            parser.add_argument("--chunk_size", type=int, default=10000)
            args = parser.parse_args()

            is_bcf = args.input_path.lower().endswith('.bcf')
            if is_bcf:
                samples_out = subprocess.run(['bcftools', 'query', '-l', args.input_path], stdout=subprocess.PIPE, text=True, check=True)
                vcf_samples = [s for s in samples_out.stdout.strip().split('\n') if s]
            else:
                headers = allel.read_vcf_headers(args.input_path)
                vcf_samples = headers.samples
            
            trios, subset_samples = parse_pedigree(args.ped_path, vcf_samples)
            num_trios = len(trios)
            
            agg_results = process_vcf_in_chunks(args.input_path, trios, subset_samples, args.chunk_size, num_trios)
            
            rows_0, rows_09 = [], []
            for key, res in agg_results.items():
                af_bin, len_bin, in_trh, min_gp = key
                row = {'AF_BIN': af_bin, 'LENGTH_BIN': len_bin, 'IN_TRH': in_trh, 'ERROR_VT': res['errors_per_trio'], 'NON_HOM_REF_VT': res['nhr_per_trio'], 'LOCUS_RATES': res['locus_rates'], 'NUM_LOCI': res['n_loci']}
                if min_gp == 0.0: rows_0.append(row)
                elif min_gp == 0.9: rows_09.append(row)
            
            cols = ['AF_BIN', 'LENGTH_BIN', 'IN_TRH', 'ERROR_VT', 'NON_HOM_REF_VT', 'LOCUS_RATES', 'NUM_LOCI']
            
            # Use explicit Dataframe instantiations to prevent syntax errors
            df_0 = pd.DataFrame(rows_0) if rows_0 else pd.DataFrame(columns=cols)
            df_09 = pd.DataFrame(rows_09) if rows_09 else pd.DataFrame(columns=cols)
            
            df_0.to_pickle(f"{args.output_prefix}-unfiltered.pkl")
            df_09.to_pickle(f"{args.output_prefix}-filtered-0.9.pkl")

        if __name__ == "__main__": main()
        EOF
    >>>

    output {
        File unfiltered_pkl = "~{output_prefix}-unfiltered.pkl"
        File filtered_pkl = "~{output_prefix}-filtered-0.9.pkl"
    }

    RuntimeAttr default_attr = object { cpu_cores: 4, mem_gb: 16, disk_gb: disk_gb, boot_disk_gb: 10, disk_type: "SSD", preemptible_tries: 1, max_retries: 0, docker: "continuumio/miniconda3:latest" }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    
    runtime { 
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores]) 
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB" 
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type]) 
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb]) 
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries]) 
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries]) 
        docker: select_first([runtime_attr.docker, default_attr.docker]) 
    }
}

task PlotMendelianMetrics {
    input {
        Array[File] unfiltered_pkls
        Array[File] filtered_pkls
        File pedigree
        String output_prefix
        RuntimeAttr? runtime_attr_override
    }
    Int disk_gb = 10 + ceil(size(unfiltered_pkls, "GiB") + size(filtered_pkls, "GiB"))
    command <<<
        set -euxo pipefail
        conda install -y -c bioconda -c conda-forge pandas numpy matplotlib seaborn pyarrow
        python - "~{sep=',' unfiltered_pkls}" "~{sep=',' filtered_pkls}" "~{pedigree}" "~{output_prefix}" <<-'EOF'
        import sys, pandas as pd, numpy as np, matplotlib.pyplot as plt, seaborn as sns, warnings

        unfiltered_files = [f for f in sys.argv[1].split(',') if f.strip()]
        filtered_files = [f for f in sys.argv[2].split(',') if f.strip()]
        pedigree_path = sys.argv[3]
        output_prefix = sys.argv[4]

        test_df = pd.read_pickle(unfiltered_files[0]) if unfiltered_files else pd.DataFrame()
        num_trios = len(test_df.iloc[0]['ERROR_VT']) if not test_df.empty else 0

        def aggregate_pkls(file_list, min_gp):
            combined_rows = []
            for f in file_list:
                df = pd.read_pickle(f)
                if df.empty: continue
                combined_rows.extend(df.to_dict('records'))
            
            if not combined_rows: return {}
            combined_df = pd.DataFrame(combined_rows)
            
            agg_results = {}
            for key, group in combined_df.groupby(['AF_BIN', 'LENGTH_BIN', 'IN_TRH']):
                agg_results[(*key, min_gp)] = {
                    'errors_per_trio': np.sum(group['ERROR_VT'].values, axis=0),
                    'nhr_per_trio': np.sum(group['NON_HOM_REF_VT'].values, axis=0),
                    'locus_rates': [rate for rates in group['LOCUS_RATES'] for rate in rates],
                    'n_loci': group['NUM_LOCI'].sum()
                }
            return agg_results

        agg_results = {**aggregate_pkls(unfiltered_files, 0.0), **aggregate_pkls(filtered_files, 0.9)}

        length_bin_labels = ['(-inf, -50]', '(-50, -1]', 'SNP', '[0, 50)', '[50, inf)']
        af_bin_labels = ['[0, 0.01)', '[0.01, 0.1)', '[0.1, 1]']
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            # 1. Per-Trio Boxplot
            for in_trh in [False, True]:
                fig, ax = plt.subplots(3, 1, figsize=(9, 10))
                tr_tag = 'non-TR/homopolymer' if not in_trh else 'TR/homopolymer'
                ax[0].set_title(f'{num_trios} trios, {tr_tag}')
                
                for i, af_bin_label in enumerate(af_bin_labels):
                    plt_df_values = []
                    for j, length_bin_label in enumerate(length_bin_labels):
                        for g, min_tar_gp in enumerate([0.0, 0.9]):
                            key = (af_bin_label, length_bin_label, in_trh, min_tar_gp)
                            if key in agg_results:
                                res = agg_results[key]
                                errs = res['errors_per_trio']
                                nhr = res['nhr_per_trio']
                                error_rates_per_trio = np.divide(errs, nhr, out=np.full(num_trios, np.nan), where=(nhr > 0))
                                mean_num_non_hom_ref = nhr.mean()
                                min_tar_gp_label = 'unfiltered' if min_tar_gp == 0 else f'GP > {min_tar_gp}'
                                len_text = g * '\n' + f'$\\langle N_{{l}} \\rangle={mean_num_non_hom_ref:.2f}$' + ('\n\n\n' + length_bin_label if g == 1 else '')
                                plt_df_values.extend([[min_tar_gp_label, len_text, error_rates_per_trio[t]] for t in range(num_trios) if not np.isnan(error_rates_per_trio[t])])
                                    
                    if plt_df_values:
                        plt_df = pd.DataFrame(plt_df_values, columns=['MIN_TAR_GP_TEXT', 'LENGTH_BIN_TEXT', 'ERROR_RATE'])
                        sns.boxplot(data=plt_df, x='LENGTH_BIN_TEXT', y='ERROR_RATE', hue='MIN_TAR_GP_TEXT', ax=ax[i], legend=i==2)
                        ax[i].set_xlabel('ALT length - REF length (bp)' if i == len(af_bin_labels) - 1 else None)
                        ax[i].set_ylabel(('panel allele frequency\n' if i == 1 else '\n\n') + f'{af_bin_label}\n\n' + ('Mendelian error rate per trio' if i == 1 else ''))
                        ax[i].set_yscale('symlog', linthresh=0.001)
                        ax[i].set_ylim([-1E-4, 1.001])
                        if i == 2:
                            handles, labels = ax[i].get_legend_handles_labels()
                            ax[i].legend(handles=handles, labels=labels, loc='upper center', fontsize=8)
                plt.tight_layout()
                plt.savefig(f"{output_prefix}.trio.{'inTRH' if in_trh else 'outTRH'}.png")
                plt.savefig(f"{output_prefix}.trio.{'inTRH' if in_trh else 'outTRH'}.pdf")
                plt.close()

            # 2. Per-Locus Boxplot
            for in_trh in [False, True]:
                fig, ax = plt.subplots(3, 1, figsize=(9, 12))
                tr_tag = 'non-TR/homopolymer' if not in_trh else 'TR/homopolymer'
                ax[0].set_title(f'{num_trios} trios, {tr_tag}')
                
                for i, af_bin_label in enumerate(af_bin_labels):
                    plt_df_values = []
                    for j, length_bin_label in enumerate(length_bin_labels):
                        for g, min_tar_gp in enumerate([0.0, 0.9]):
                            key = (af_bin_label, length_bin_label, in_trh, min_tar_gp)
                            if key in agg_results:
                                res = agg_results[key]
                                error_rates_per_locus = res['locus_rates']
                                num_loci = len(error_rates_per_locus)
                                mean_num_non_hom_ref_trios = res['nhr_per_trio'].sum() / max(num_loci, 1)
                                min_tar_gp_label = 'unfiltered' if min_tar_gp == 0 else f'GP > {min_tar_gp}'
                                len_text = g * '\n\n' + f'$\\langle N_{{t}} \\rangle$={mean_num_non_hom_ref_trios:.2f}' + f'\n$N_{{l}}={num_loci}$' + ('\n\n\n\n' + length_bin_label if g == 1 else '')
                                plt_df_values.extend([[min_tar_gp_label, len_text, rate] for rate in error_rates_per_locus])

                    if plt_df_values:
                        plt_df = pd.DataFrame(plt_df_values, columns=['MIN_TAR_GP_TEXT', 'LENGTH_BIN_TEXT', 'ERROR_RATE'])
                        sns.boxplot(data=plt_df, x='LENGTH_BIN_TEXT', y='ERROR_RATE', hue='MIN_TAR_GP_TEXT', ax=ax[i], legend=i==2)
                        ax[i].set_xlabel('ALT length - REF length (bp)' if i == len(af_bin_labels) - 1 else None)
                        ax[i].set_ylabel(('panel allele frequency\n' if i == 1 else '\n\n') + f'{af_bin_label}\n\n' + ('Mendelian error rate per locus' if i == 1 else ''))
                        ax[i].set_yscale('symlog', linthresh=0.001)
                        ax[i].set_ylim([-1E-5, 1])
                        if i == 2:
                            handles, labels = ax[i].get_legend_handles_labels()
                            ax[i].legend(handles=handles, labels=labels, loc='upper center', fontsize=8)
                plt.tight_layout()
                plt.savefig(f"{output_prefix}.locus.{'inTRH' if in_trh else 'outTRH'}.png")
                plt.savefig(f"{output_prefix}.locus.{'inTRH' if in_trh else 'outTRH'}.pdf")
                plt.close()
        EOF
    >>>
    output {
        Array[File] plots_png = glob("*.png")
        Array[File] plots_pdf = glob("*.pdf")
    }
    RuntimeAttr default_attr = object { cpu_cores: 4, mem_gb: 16, disk_gb: disk_gb, boot_disk_gb: 10, disk_type: "SSD", preemptible_tries: 1, max_retries: 0, docker: "continuumio/miniconda3:latest" }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime { cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores]) memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB" disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type]) bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb]) preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries]) maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries]) docker: select_first([runtime_attr.docker, default_attr.docker]) }
}
