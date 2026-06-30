version 1.0

workflow GLIMPSE2Summarize {
    input {
        Array[File] panel_vcfs
        Array[File] panel_vcf_idxs
        Array[File] imputed_vcfs
        Array[File] imputed_vcf_idxs
        File population_tsv
        String output_prefix
    }

    scatter (idx in range(length(panel_vcfs))) {
        call SummarizeMetrics { input:
            panel_vcf = panel_vcfs[idx],
            panel_vcf_idx = panel_vcf_idxs[idx],
            imputed_vcf = imputed_vcfs[idx],
            imputed_vcf_idx = imputed_vcf_idxs[idx],
            output_prefix = output_prefix + "." + idx
        }

        # Generate Per-Chromosome Plots
        call PlotSummaries as PlotSummariesPerChrom { input:
            summary_pkls = [SummarizeMetrics.summary_pkl],
            population_tsv = population_tsv,
            output_prefix = output_prefix + "." + idx
        }
    }

    # Generate Aggregate Plots
    call PlotSummaries as PlotSummariesAggregate { input:
        summary_pkls = SummarizeMetrics.summary_pkl,
        population_tsv = population_tsv,
        output_prefix = output_prefix + ".aggregate"
    }

    output {
        File summarize_pearson_tsv = PlotSummariesAggregate.summarize_pearson_tsv
        Array[File] summarize_aggregate_plots_png = PlotSummariesAggregate.plots_png
        Array[File] summarize_aggregate_plots_pdf = PlotSummariesAggregate.plots_pdf
        Array[File] summarize_per_chrom_plots_png = flatten(PlotSummariesPerChrom.plots_png)
        Array[File] summarize_per_chrom_plots_pdf = flatten(PlotSummariesPerChrom.plots_pdf)
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

task SummarizeMetrics {
    input {
        File panel_vcf
        File panel_vcf_idx
        File imputed_vcf
        File imputed_vcf_idx
        String output_prefix
        RuntimeAttr? runtime_attr_override
    }
    Int disk_gb = 20 + ceil(size(panel_vcf, "GiB") + size(imputed_vcf, "GiB"))
    
    command <<<
        set -euxo pipefail
        conda install -y -c bioconda -c conda-forge cyvcf2 numpy pandas
        
        python - ~{panel_vcf} ~{imputed_vcf} ~{output_prefix} <<-'EOF'
        import sys, numpy as np, pandas as pd, cyvcf2, pickle

        panel_vcf = cyvcf2.VCF(sys.argv[1])
        imputed_vcf = cyvcf2.VCF(sys.argv[2])

        num_p_samples, num_c_samples = len(panel_vcf.samples), len(imputed_vcf.samples)

        altlen_all, panel_af_all, target_af_all = [], [], []
        panel_hwe = {'hom_ref': [], 'het': [], 'hom_alt': []}
        target_hwe = {'hom_ref': [], 'het': [], 'hom_alt': []}
        panel_mean_alt_alleles_all, target_mean_alt_alleles_all = [], []

        p_het_all, p_hom_ref_all, p_hom_alt_all = np.zeros(num_p_samples, dtype=int), np.zeros(num_p_samples, dtype=int), np.zeros(num_p_samples, dtype=int)
        p_het_ins, p_het_del = np.zeros(num_p_samples, dtype=int), np.zeros(num_p_samples, dtype=int)
        c_het_all, c_hom_ref_all, c_hom_alt_all = np.zeros(num_c_samples, dtype=int), np.zeros(num_c_samples, dtype=int), np.zeros(num_c_samples, dtype=int)
        c_het_ins, c_het_del = np.zeros(num_c_samples, dtype=int), np.zeros(num_c_samples, dtype=int)

        p_is_het, p_is_hom_ref, p_is_hom_alt = np.empty(num_p_samples, dtype=bool), np.empty(num_p_samples, dtype=bool), np.empty(num_p_samples, dtype=bool)
        c_is_het, c_is_hom_ref, c_is_hom_alt = np.empty(num_c_samples, dtype=bool), np.empty(num_c_samples, dtype=bool), np.empty(num_c_samples, dtype=bool)

        for p_var, c_var in zip(panel_vcf, imputed_vcf):
            # Strict safety check to prevent mismatched chromosomes/positions
            if p_var.CHROM != c_var.CHROM or p_var.POS != c_var.POS:
                raise ValueError(
                    f"VCFs are out of sync! Check your input arrays.\n"
                    f"Panel : {p_var.CHROM}:{p_var.POS}\n"
                    f"Target: {c_var.CHROM}:{c_var.POS}"
                )

            alts = p_var.ALT
            if not alts: continue
            altlen = len(alts[0]) - len(p_var.REF)
            altlen_all.append(altlen)
            is_ins, is_del = altlen >= 50, altlen <= -50

            # --- Panel Extraction ---
            p_gt = p_var.gt_types
            np.equal(p_gt, 1, out=p_is_het); np.equal(p_gt, 0, out=p_is_hom_ref); np.equal(p_gt, 3, out=p_is_hom_alt)
            p_het_all += p_is_het; p_hom_ref_all += p_is_hom_ref; p_hom_alt_all += p_is_hom_alt
            if is_ins: p_het_ins += p_is_het
            elif is_del: p_het_del += p_is_het
            
            panel_hwe['hom_ref'].append(p_var.num_hom_ref); panel_hwe['het'].append(p_var.num_het); panel_hwe['hom_alt'].append(p_var.num_hom_alt)
            p_valid = (p_var.num_hom_ref + p_var.num_het + p_var.num_hom_alt) * 2
            p_alt_count = p_var.num_het + 2 * p_var.num_hom_alt
            panel_af_all.append(p_alt_count / p_valid if p_valid > 0 else 0.0)
            
            # Safe division to support sites-only VCFs
            panel_mean_alt_alleles_all.append(p_alt_count / num_p_samples if num_p_samples > 0 else 0.0)

            # --- Target Extraction ---
            c_gt = c_var.gt_types
            np.equal(c_gt, 1, out=c_is_het); np.equal(c_gt, 0, out=c_is_hom_ref); np.equal(c_gt, 3, out=c_is_hom_alt)
            c_het_all += c_is_het; c_hom_ref_all += c_is_hom_ref; c_hom_alt_all += c_is_hom_alt
            if is_ins: c_het_ins += c_is_het
            elif is_del: c_het_del += c_is_het
            
            target_hwe['hom_ref'].append(c_var.num_hom_ref); target_hwe['het'].append(c_var.num_het); target_hwe['hom_alt'].append(c_var.num_hom_alt)
            c_valid = (c_var.num_hom_ref + c_var.num_het + c_var.num_hom_alt) * 2
            c_alt_count = c_var.num_het + 2 * c_var.num_hom_alt
            target_af_all.append(c_alt_count / c_valid if c_valid > 0 else 0.0)
            
            # Safe division
            target_mean_alt_alleles_all.append(c_alt_count / num_c_samples if num_c_samples > 0 else 0.0)

        stats_dict = {
            'altlen_all': altlen_all, 'panel_af_all': panel_af_all, 'target_af_all': target_af_all,
            'panel_hwe': panel_hwe, 'target_hwe': target_hwe,
            'panel_mean_alt_alleles_all': panel_mean_alt_alleles_all, 'target_mean_alt_alleles_all': target_mean_alt_alleles_all,
            'sample_stats': {
                'panel': {'het_all': p_het_all, 'hom_ref_all': p_hom_ref_all, 'hom_alt_all': p_hom_alt_all, 'het_ins': p_het_ins, 'het_del': p_het_del},
                'target': {'het_all': c_het_all, 'hom_ref_all': c_hom_ref_all, 'hom_alt_all': c_hom_alt_all, 'het_ins': c_het_ins, 'het_del': c_het_del}
            },
            'panel_samples': list(panel_vcf.samples), 'target_samples': list(imputed_vcf.samples)
        }
        
        with open(f"{sys.argv[3]}.summary.pkl", "wb") as f:
            pickle.dump(stats_dict, f)
        EOF
    >>>

    output {
        File summary_pkl = "~{output_prefix}.summary.pkl"
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

task PlotSummaries {
    input {
        Array[File] summary_pkls
        File population_tsv
        String output_prefix
        RuntimeAttr? runtime_attr_override
    }
    Int disk_gb = 20 + ceil(size(summary_pkls, "GiB"))
    command <<<
        set -euxo pipefail
        conda install -y -c bioconda -c conda-forge numpy pandas matplotlib seaborn scipy
        python - "~{sep=',' summary_pkls}" "~{population_tsv}" "~{output_prefix}" <<-'EOF'
        import sys, pickle, numpy as np, pandas as pd, matplotlib.pyplot as plt, seaborn as sns, matplotlib.colors, warnings
        from scipy.stats import pearsonr

        warnings.filterwarnings('ignore', category=RuntimeWarning)
        pkl_files = [f for f in sys.argv[1].split(',') if f.strip()]
        population_tsv_path = sys.argv[2]
        output_prefix = sys.argv[3]

        # Combine all dictionaries
        agg = {'altlen_all': [], 'panel_af_all': [], 'target_af_all': [], 'panel_mean_alt_alleles_all': [], 'target_mean_alt_alleles_all': []}
        panel_hwe = {'hom_ref': [], 'het': [], 'hom_alt': []}
        target_hwe = {'hom_ref': [], 'het': [], 'hom_alt': []}
        
        sample_stats = None

        for f in pkl_files:
            with open(f, 'rb') as fp: data = pickle.load(fp)
            agg['altlen_all'].extend(data['altlen_all'])
            agg['panel_af_all'].extend(data['panel_af_all'])
            agg['target_af_all'].extend(data['target_af_all'])
            agg['panel_mean_alt_alleles_all'].extend(data['panel_mean_alt_alleles_all'])
            agg['target_mean_alt_alleles_all'].extend(data['target_mean_alt_alleles_all'])
            
            for k in ['hom_ref', 'het', 'hom_alt']:
                panel_hwe[k].extend(data['panel_hwe'][k])
                target_hwe[k].extend(data['target_hwe'][k])
            
            if sample_stats is None:
                sample_stats = data['sample_stats']
                panel_samples = data['panel_samples']
                target_samples = data['target_samples']
            else:
                for grp in ['panel', 'target']:
                    for k in sample_stats[grp]:
                        sample_stats[grp][k] += data['sample_stats'][grp][k]

        panel_af = np.array(agg['panel_af_all'])
        target_af = np.array(agg['target_af_all'])
        altlen = np.array(agg['altlen_all'])
        is_sv_ins = altlen >= 50
        is_sv_del = altlen <= -50

        # Pearson
        def safe_pearson(x, y):
            if len(x) > 1 and np.std(x) > 0 and np.std(y) > 0: return pearsonr(x, y)[0]
            return np.nan

        pearson_records = [{"VARIANT_TYPE": "ALL", "PEARSON_R": safe_pearson(panel_af, target_af)}]
        if np.sum(is_sv_ins) > 1: pearson_records.append({"VARIANT_TYPE": "SV_INS", "PEARSON_R": safe_pearson(panel_af[is_sv_ins], target_af[is_sv_ins])})
        if np.sum(is_sv_del) > 1: pearson_records.append({"VARIANT_TYPE": "SV_DEL", "PEARSON_R": safe_pearson(panel_af[is_sv_del], target_af[is_sv_del])})
        is_sv = is_sv_ins | is_sv_del
        if np.sum(is_sv) > 1: pearson_records.append({"VARIANT_TYPE": "SV", "PEARSON_R": safe_pearson(panel_af[is_sv], target_af[is_sv])})
        pd.DataFrame(pearson_records).to_csv(f'{output_prefix}.pearson.tsv', sep='\t', index=False)

        # Plots
        def plot_hist2d(p_af, c_af, title, outfile):
            if len(p_af) == 0: return
            plt.figure()
            plt.hist2d(p_af, c_af, bins=np.linspace(0, 1, 50), norm=matplotlib.colors.LogNorm())
            plt.title(title); plt.xlabel('AoU+HPRC2+HGSVC3 allele frequency'); plt.ylabel('Target allele frequency')
            plt.gca().set_aspect('equal')
            plt.colorbar().set_label('Number of variants', rotation=270, labelpad=10)
            plt.savefig(f'{outfile}.png', bbox_inches='tight')
            plt.savefig(f'{outfile}.pdf', bbox_inches='tight')
            plt.close()

        plot_hist2d(panel_af, target_af, 'All variants', f'{output_prefix}-AF-all')
        plot_hist2d(panel_af[is_sv_ins], target_af[is_sv_ins], 'SV-length insertions', f'{output_prefix}-AF-SV-ins')
        plot_hist2d(panel_af[is_sv_del], target_af[is_sv_del], 'SV-length deletions', f'{output_prefix}-AF-SV-del')

        population_df = pd.read_csv(population_tsv_path, sep='\t')
        pop_color_dict = {'ESN': '#ffcd00', 'GWD': '#ffb900', 'LWK': '#cc9933', 'MSL': '#e1b919', 'YRI': '#ffb933', 'ACB': '#ff9900', 'ASW': '#cc6600', 'CLM': '#cc3333', 'MXL': '#e10033', 'PEL': '#ff0000', 'PUR': '#cc3300', 'CDX': '#339900', 'CHB': '#adcd00', 'CHS': '#01ff00', 'JPT': '#008b00', 'KHV': '#00cc33', 'CEU': '#0000ff', 'GBR': '#00c5cd', 'FIN': '#00ebff', 'IBS': '#6495ed', 'TSI': '#00008b', 'BEB': '#8b008b', 'GIH': '#9400d3', 'ITU': '#b03060', 'PJL': '#e11289', 'STU': '#ff00ff'}

        def make_results_df(samples, stats):
            if stats is None: return pd.DataFrame()
            df = pd.DataFrame({'Sample': samples, 'Heterozygous variants per sample': stats['het_all'], 'Homozygous reference variants per sample': stats['hom_ref_all'], 'Homozygous alternate variants per sample': stats['hom_alt_all'], 'Heterozygous SV-length insertions per sample': stats['het_ins'], 'Heterozygous SV-length deletions per sample': stats['het_del']})
            return pd.merge(df, population_df[['Sample name', 'Population code']], left_on='Sample', right_on='Sample name').rename(columns={'Population code': 'Population'})

        panel_results_df = make_results_df(panel_samples, sample_stats['panel'] if sample_stats else None)
        target_results_df = make_results_df(target_samples, sample_stats['target'] if sample_stats else None)

        def plot_boxplot(df, x_col, title, outfile, xlim):
            if df.empty: return
            plt.figure(figsize=(4, 6))
            sns.boxplot(data=df, x=x_col, y='Population', hue='Population', order=pop_color_dict.keys(), hue_order=pop_color_dict.keys(), palette=pop_color_dict.values())
            plt.title(title); plt.xlim(xlim)
            plt.savefig(f'{outfile}.png', bbox_inches='tight')
            plt.savefig(f'{outfile}.pdf', bbox_inches='tight')
            plt.close()

        plot_boxplot(panel_results_df, 'Heterozygous variants per sample', 'HPRC2+HGSVC3 in AoU+HPRC2+HGSVC3', f'{output_prefix}-panel-het-all', [0, 6E4])
        plot_boxplot(target_results_df, 'Heterozygous variants per sample', 'Target', f'{output_prefix}-target-het-all', [0, 6E4])
        plot_boxplot(panel_results_df, 'Heterozygous SV-length insertions per sample', 'HPRC2+HGSVC3 in AoU+HPRC2+HGSVC3', f'{output_prefix}-panel-het-SV-ins', [0, 500])
        plot_boxplot(target_results_df, 'Heterozygous SV-length insertions per sample', 'Target', f'{output_prefix}-target-het-SV-ins', [0, 500])
        plot_boxplot(panel_results_df, 'Heterozygous SV-length deletions per sample', 'HPRC2+HGSVC3 in AoU+HPRC2+HGSVC3', f'{output_prefix}-panel-het-SV-del', [0, 500])
        plot_boxplot(target_results_df, 'Heterozygous SV-length deletions per sample', 'Target', f'{output_prefix}-target-het-SV-del', [0, 500])

        if len(altlen) > 0:
            bins = list(np.linspace(-10000, -100, 397)) + [-75, -50, -25, -1, -0.1, 0.1, 1, 25, 50, 75] + list(np.linspace(100, 10000, 397))
            plt.figure()
            plt.hist(altlen, bins=bins, label='AoU+HPRC2+HGSVC3', log=True, histtype='step', alpha=0.5, weights=agg['panel_mean_alt_alleles_all'])
            plt.hist(altlen, bins=bins, label='Target', log=True, histtype='step', alpha=0.5, weights=agg['target_mean_alt_alleles_all'])
            plt.xscale('symlog', linthresh=100)
            plt.xticks([-1E4, -1E3] + list(np.linspace(-100, 100, 9)) + [1E3, 1E4], labels=['$-10^4$', '$-10^3$', '$-10^2$'] + ['', '$-50$', '', '$0$', '', '$50$', ''] + ['$10^2$', '$10^3$', '$10^4$'])
            plt.ylabel('Number of ALT alleles per sample'); plt.xlabel('ALT length - REF length (bp)'); plt.legend()
            plt.savefig(f'{output_prefix}-alt-alleles-per-sample-hist.png', bbox_inches='tight')
            plt.savefig(f'{output_prefix}-alt-alleles-per-sample-hist.pdf', bbox_inches='tight')
            plt.close()

        ternary_to_cartesian = lambda a, b, c: (0.5 * (2 * b + c) / (a + b + c + 1E-10), 0.5 * np.sqrt(3) * c / (a + b + c + 1E-10))
        def calc_hwe_ternary(x, m=1, f=1): return np.array([1 - f, 0, 0]) + f * np.array([m * (1 - x)**2 + (1 - m) * (1 - x), m * x**2 + (1 - m) * x, 2 * m * x * (1 - x)])
        
        def plot_de_finetti(hom_ref_arr, hom_alt_arr, het_arr, title, outfile, gridsize=70):
            if len(hom_ref_arr) == 0: return
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.set_xlim([-0.1, 1.1]); ax.set_ylim([-0.1, np.sqrt(3) / 2 + 0.1]); ax.set_aspect(np.sqrt(3) / 2); ax.axis('off')
            ep = 0.02
            ax.plot([-2 * ep / np.sqrt(3), 1 + 2 * ep / np.sqrt(3)], [-ep, -ep], lw=3, c='k')
            ax.plot([-2 * ep / np.sqrt(3), 0.5], [-ep, 0.5 * np.sqrt(3) + ep], lw=3, c='k')
            ax.plot([1 + 2 * ep / np.sqrt(3), 0.5], [-ep, 0.5 * np.sqrt(3) + ep], lw=3, c='k')
            ax.text(-4 * ep / np.sqrt(3), -4 * ep, 'HOM\nREF ', fontsize=16, ha='right')
            ax.text(0.5, 0.5 * np.sqrt(3) + 3 * ep, 'HET', fontsize=16, ha='center')
            ax.text(1 + 4 * ep / np.sqrt(3), -4 * ep, 'HOM\n ALT', fontsize=16, ha='left')
            
            counts_arr = np.column_stack((hom_ref_arr, hom_alt_arr, het_arr))
            unique_counts, point_weights = np.unique(counts_arr, axis=0, return_counts=True)
            x_ternary_v, y_ternary_v = ternary_to_cartesian(unique_counts[:,0], unique_counts[:,1], unique_counts[:,2])
            hb = ax.hexbin(x_ternary_v, y_ternary_v, C=point_weights, reduce_C_function=np.sum, gridsize=gridsize, extent=[0, 1, 0, np.sqrt(3) / 2], norm=matplotlib.colors.LogNorm())
            plt.colorbar(hb, ax=ax, shrink=0.5).set_label('Number of variants', rotation=270, labelpad=10)
            
            x_values = np.linspace(0, 1, 50)
            cart_values = np.array([ternary_to_cartesian(*calc_hwe_ternary(x)) for x in x_values])
            ax.plot(cart_values[:, 0], cart_values[:, 1], c='C1', ls='solid', lw=3)
            ax.text(0.5, -0.3, title, fontsize=18, ha='center')
            plt.savefig(f'{outfile}.png', bbox_inches='tight')
            plt.savefig(f'{outfile}.pdf', bbox_inches='tight')
            plt.close()

        def slice_hwe(hwe_dict, mask): return {k: np.array(v)[mask] for k, v in hwe_dict.items()}

        p_hwe_ins, p_hwe_del = slice_hwe(panel_hwe, is_sv_ins), slice_hwe(panel_hwe, is_sv_del)
        c_hwe_ins, c_hwe_del = slice_hwe(target_hwe, is_sv_ins), slice_hwe(target_hwe, is_sv_del)

        plot_de_finetti(panel_hwe['hom_ref'], panel_hwe['hom_alt'], panel_hwe['het'], 'AoU+HPRC2+HGSVC3\nAll variants', f'{output_prefix}-panel-hwe-all')
        plot_de_finetti(p_hwe_ins['hom_ref'], p_hwe_ins['hom_alt'], p_hwe_ins['het'], 'AoU+HPRC2+HGSVC3\nSV-length insertions', f'{output_prefix}-panel-hwe-SV-ins')
        plot_de_finetti(p_hwe_del['hom_ref'], p_hwe_del['hom_alt'], p_hwe_del['het'], 'AoU+HPRC2+HGSVC3\nSV-length deletions', f'{output_prefix}-panel-hwe-SV-del')
        plot_de_finetti(target_hwe['hom_ref'], target_hwe['hom_alt'], target_hwe['het'], 'Target\nAll variants', f'{output_prefix}-target-hwe-all', gridsize=80)
        plot_de_finetti(c_hwe_ins['hom_ref'], c_hwe_ins['hom_alt'], c_hwe_ins['het'], 'Target\nSV-length insertions', f'{output_prefix}-target-hwe-SV-ins', gridsize=80)
        plot_de_finetti(c_hwe_del['hom_ref'], c_hwe_del['hom_alt'], c_hwe_del['het'], 'Target\nSV-length deletions', f'{output_prefix}-target-hwe-SV-del', gridsize=80)
        EOF
    >>>
    output {
        File summarize_pearson_tsv = "~{output_prefix}.pearson.tsv"
        Array[File] plots_png = glob("*.png")
        Array[File] plots_pdf = glob("*.pdf")
    }
    RuntimeAttr default_attr = object { cpu_cores: 4, mem_gb: 16, disk_gb: disk_gb, boot_disk_gb: 10, disk_type: "SSD", preemptible_tries: 1, max_retries: 0, docker: "continuumio/miniconda3:latest" }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime { cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores]) memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB" disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " " + select_first([runtime_attr.disk_type, default_attr.disk_type]) bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb]) preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries]) maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries]) docker: select_first([runtime_attr.docker, default_attr.docker]) }
}
