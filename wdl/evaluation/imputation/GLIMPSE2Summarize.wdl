version 1.0

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

workflow GLIMPSE2Summarize {
    input {
        File panel_vcf              # split to biallelic
        File panel_vcf_idx
        File imputed_vcf            # split to biallelic, variants in same order as in panel
        File imputed_vcf_idx
        File population_tsv         # should contain columns: Sample name, Population code (e.g., igsr_samples.tsv tables generated from https://www.internationalgenome.org/data/)
        String output_prefix
    }

    call SummarizeAndPlot { input:
        panel_vcf = panel_vcf,
        panel_vcf_idx = panel_vcf_idx,
        imputed_vcf = imputed_vcf,
        imputed_vcf_idx = imputed_vcf_idx,
        population_tsv = population_tsv,
        output_prefix = output_prefix
    }

    output {
        File summarize_pearson_tsv = SummarizeAndPlot.summarize_pearson_tsv
        Array[File] summarize_plots_png = SummarizeAndPlot.summarize_plots_png
        Array[File] summarize_plots_pdf = SummarizeAndPlot.summarize_plots_pdf
    }
}

task SummarizeAndPlot {
    input {
        File panel_vcf
        File panel_vcf_idx
        File imputed_vcf
        File imputed_vcf_idx
        File population_tsv
        String output_prefix
        
        Int chunk_size = 10000

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 20 + ceil(size(panel_vcf, "GiB") + size(imputed_vcf, "GiB"))

    command <<<
        set -euxo pipefail
        
        # Install dependencies for variant streaming, stats, and plotting
        conda install -y -c bioconda -c conda-forge cyvcf2 pandas numpy matplotlib seaborn scipy

        cat << 'EOF' > summarize.py
        import sys
        import time
        import numpy as np
        import pandas as pd
        import cyvcf2
        import matplotlib.pyplot as plt
        import seaborn as sns
        import matplotlib.colors
        from scipy.stats import pearsonr
        import warnings

        warnings.filterwarnings('ignore', category=RuntimeWarning)

        panel_vcf_path = sys.argv[1]
        imputed_vcf_path = sys.argv[2]
        population_tsv_path = sys.argv[3]
        output_prefix = sys.argv[4]
        chunk_size = int(sys.argv[5])

        # 1. Initialize VCF readers
        panel_vcf = cyvcf2.VCF(panel_vcf_path)
        imputed_vcf = cyvcf2.VCF(imputed_vcf_path)

        panel_samples = np.array(panel_vcf.samples)
        target_samples = np.array(imputed_vcf.samples)
        num_p_samples = len(panel_samples)
        num_c_samples = len(target_samples)

        # 2. Global Storage for Metrics
        altlen_all = []
        panel_af_all, target_af_all = [], []

        panel_hwe = {'hom_ref': [], 'het': [], 'hom_alt': []}
        target_hwe = {'hom_ref': [], 'het': [], 'hom_alt': []}

        panel_mean_alt_alleles_all, target_mean_alt_alleles_all = [], []

        # Per-sample counters
        sample_stats = {
            'panel': {
                'het_all': np.zeros(num_p_samples), 'hom_ref_all': np.zeros(num_p_samples), 'hom_alt_all': np.zeros(num_p_samples),
                'het_ins': np.zeros(num_p_samples), 'het_del': np.zeros(num_p_samples)
            },
            'target': {
                'het_all': np.zeros(num_c_samples), 'hom_ref_all': np.zeros(num_c_samples), 'hom_alt_all': np.zeros(num_c_samples),
                'het_ins': np.zeros(num_c_samples), 'het_del': np.zeros(num_c_samples)
            }
        }

        def process_chunk(p_gt_types, c_gt_types, altlens):
            global altlen_all, panel_af_all, target_af_all
            
            altlen_all.extend(altlens)
            l_arr = np.array(altlens)
            is_sv_ins = l_arr >= 50
            is_sv_del = l_arr <= -50

            for gt, prefix in [(p_gt_types, 'panel'), (c_gt_types, 'target')]:
                # In cyvcf2 gt_types: 0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT
                is_hom_ref = (gt == 0)
                is_het = (gt == 1)
                is_missing = (gt == 2)
                is_hom_alt = (gt == 3)

                # Calculate AF (total alt alleles / total valid alleles)
                valid_alleles = np.sum(~is_missing, axis=1) * 2
                alt_alleles = np.sum(is_het, axis=1) + 2 * np.sum(is_hom_alt, axis=1)
                af = np.divide(alt_alleles, valid_alleles, out=np.zeros_like(alt_alleles, dtype=float), where=valid_alleles!=0)
                
                if prefix == 'panel':
                    panel_af_all.extend(af)
                    panel_mean_alt_alleles_all.extend(np.mean(is_het * 1 + is_hom_alt * 2, axis=1))
                    hw_dict = panel_hwe
                else:
                    target_af_all.extend(af)
                    target_mean_alt_alleles_all.extend(np.mean(is_het * 1 + is_hom_alt * 2, axis=1))
                    hw_dict = target_hwe
                    
                # HWE Counts (per variant)
                hw_dict['hom_ref'].extend(np.sum(is_hom_ref, axis=1))
                hw_dict['het'].extend(np.sum(is_het, axis=1))
                hw_dict['hom_alt'].extend(np.sum(is_hom_alt, axis=1))

                # Update Sample Counters (per sample)
                sample_stats[prefix]['het_all'] += np.sum(is_het, axis=0)
                sample_stats[prefix]['hom_ref_all'] += np.sum(is_hom_ref, axis=0)
                sample_stats[prefix]['hom_alt_all'] += np.sum(is_hom_alt, axis=0)
                sample_stats[prefix]['het_ins'] += np.sum(is_het[is_sv_ins, :], axis=0)
                sample_stats[prefix]['het_del'] += np.sum(is_het[is_sv_del, :], axis=0)

        # 3. Stream Variants Out-of-Core
        print(f"Streaming variants with chunk size {chunk_size}...", flush=True)
        p_gt_chunk, c_gt_chunk, altlen_chunk = [], [], []
        variants_processed = 0
        start_time = time.time()

        for p_var, c_var in zip(panel_vcf, imputed_vcf):
            if p_var.POS != c_var.POS:
                raise ValueError(f"VCFs are out of sync at Panel POS: {p_var.POS} vs Target POS: {c_var.POS}")
            
            ref = p_var.REF
            alts = p_var.ALT
            if not alts: continue
            alt = alts[0]
            
            altlen = len(alt) - len(ref)
            altlen_chunk.append(altlen)
            
            # Use native C-array gt_types (drastically faster than extracting tuple genotypes)
            p_gt_chunk.append(p_var.gt_types)
            c_gt_chunk.append(c_var.gt_types)
            
            if len(p_gt_chunk) >= chunk_size:
                process_chunk(np.array(p_gt_chunk), np.array(c_gt_chunk), altlen_chunk)
                variants_processed += len(p_gt_chunk)
                elapsed = time.time() - start_time
                print(f"Processed {variants_processed} variants in {elapsed:.2f}s ({variants_processed/elapsed:.2f} var/s)...", flush=True)
                p_gt_chunk, c_gt_chunk, altlen_chunk = [], [], []

        if len(p_gt_chunk) > 0:
            process_chunk(np.array(p_gt_chunk), np.array(c_gt_chunk), altlen_chunk)
            variants_processed += len(p_gt_chunk)
            elapsed = time.time() - start_time
            print(f"Finished processing {variants_processed} total variants in {elapsed:.2f}s.", flush=True)

        # Convert globals to numpy arrays for fast indexing
        panel_af = np.array(panel_af_all)
        target_af = np.array(target_af_all)
        altlen = np.array(altlen_all)
        is_sv_ins = altlen >= 50
        is_sv_del = altlen <= -50

        # 4. Pearson Correlations (TSV Output)
        print("Calculating Pearson correlations...", flush=True)
        pearson_records = []
        
        # Helper to safely calculate pearsonr against edge-case 0-variance bins
        def safe_pearson(x, y):
            if len(x) > 1 and np.std(x) > 0 and np.std(y) > 0:
                return pearsonr(x, y)[0]
            return np.nan
        
        r_all = safe_pearson(panel_af, target_af)
        pearson_records.append({"VARIANT_TYPE": "ALL", "PEARSON_R": r_all})
            
        if np.sum(is_sv_ins) > 1:
            r_ins = safe_pearson(panel_af[is_sv_ins], target_af[is_sv_ins])
            pearson_records.append({"VARIANT_TYPE": "SV_INS", "PEARSON_R": r_ins})
                
        if np.sum(is_sv_del) > 1:
            r_del = safe_pearson(panel_af[is_sv_del], target_af[is_sv_del])
            pearson_records.append({"VARIANT_TYPE": "SV_DEL", "PEARSON_R": r_del})
                
        is_sv = is_sv_ins | is_sv_del
        if np.sum(is_sv) > 1:
            r_sv = safe_pearson(panel_af[is_sv], target_af[is_sv])
            pearson_records.append({"VARIANT_TYPE": "SV", "PEARSON_R": r_sv})

        pd.DataFrame(pearson_records).to_csv(f'{output_prefix}.pearson.tsv', sep='\t', index=False)

        # 5. AF Hist2D Plots
        print("Generating AF Hist2D Plots...", flush=True)
        def plot_hist2d(p_af, c_af, title, outfile):
            if len(p_af) == 0: return
            plt.figure()
            plt.hist2d(p_af, c_af, bins=50, norm=matplotlib.colors.LogNorm())
            plt.title(title)
            plt.xlabel('AoU+HPRC2+HGSVC3 allele frequency')
            plt.ylabel('Target allele frequency')
            plt.gca().set_aspect('equal')
            cbar = plt.colorbar()
            cbar.set_label('Number of variants', rotation=270, labelpad=10)
            plt.savefig(f'{outfile}.png', bbox_inches='tight')
            plt.savefig(f'{outfile}.pdf', bbox_inches='tight')
            plt.close()

        plot_hist2d(panel_af, target_af, 'All variants', f'{output_prefix}-AF-all')
        plot_hist2d(panel_af[is_sv_ins], target_af[is_sv_ins], 'SV-length insertion bubbles', f'{output_prefix}-AF-SV-ins')
        plot_hist2d(panel_af[is_sv_del], target_af[is_sv_del], 'SV-length deletion bubbles', f'{output_prefix}-AF-SV-del')

        # 6. Sample Metrics & Boxplots
        print("Generating population boxplots...", flush=True)
        population_df = pd.read_csv(population_tsv_path, sep='\t')

        pop_color_dict = {
            'ESN': '#ffcd00', 'GWD': '#ffb900', 'LWK': '#cc9933', 'MSL': '#e1b919', 'YRI': '#ffb933',
            'ACB': '#ff9900', 'ASW': '#cc6600', 'CLM': '#cc3333', 'MXL': '#e10033', 'PEL': '#ff0000',
            'PUR': '#cc3300', 'CDX': '#339900', 'CHB': '#adcd00', 'CHS': '#01ff00', 'JPT': '#008b00',
            'KHV': '#00cc33', 'CEU': '#0000ff', 'GBR': '#00c5cd', 'FIN': '#00ebff', 'IBS': '#6495ed',
            'TSI': '#00008b', 'BEB': '#8b008b', 'GIH': '#9400d3', 'ITU': '#b03060', 'PJL': '#e11289',
            'STU': '#ff00ff', 'UNK': 'grey'
        }

        def make_results_df(samples, stats):
            df = pd.DataFrame()
            df['Sample'] = samples
            df['Heterozygous variants per sample'] = stats['het_all']
            df['Homozygous reference variants per sample'] = stats['hom_ref_all']
            df['Homozygous alternate variants per sample'] = stats['hom_alt_all']
            df['Heterozygous SV-length insertion bubbles per sample'] = stats['het_ins']
            df['Heterozygous SV-length deletion bubbles per sample'] = stats['het_del']
            
            df = pd.merge(df, population_df[['Sample name', 'Population code']], left_on='Sample', right_on='Sample name')
            return df.rename(columns={'Population code': 'Population'})

        panel_results_df = make_results_df(panel_samples, sample_stats['panel'])
        target_results_df = make_results_df(target_samples, sample_stats['target'])

        def plot_boxplot(df, x_col, title, outfile, xlim):
            plt.figure(figsize=(4, 6))
            sns.boxplot(data=df, x=x_col, y='Population', hue='Population',
                        order=pop_color_dict.keys(), hue_order=pop_color_dict.keys(), 
                        palette=pop_color_dict.values())
            plt.title(title)
            plt.xlim(xlim)
            plt.savefig(f'{outfile}.png', bbox_inches='tight')
            plt.savefig(f'{outfile}.pdf', bbox_inches='tight')
            plt.close()

        plot_boxplot(panel_results_df, 'Heterozygous variants per sample', 'HPRC2+HGSVC3 in AoU+HPRC2+HGSVC3', f'{output_prefix}-panel-het-all', [0, 6E4])
        plot_boxplot(target_results_df, 'Heterozygous variants per sample', 'Target', f'{output_prefix}-target-het-all', [0, 6E4])
        plot_boxplot(panel_results_df, 'Heterozygous SV-length insertion bubbles per sample', 'HPRC2+HGSVC3 in AoU+HPRC2+HGSVC3', f'{output_prefix}-panel-het-SV-ins', [0, 500])
        plot_boxplot(target_results_df, 'Heterozygous SV-length insertion bubbles per sample', 'Target', f'{output_prefix}-target-het-SV-ins', [0, 500])
        plot_boxplot(panel_results_df, 'Heterozygous SV-length deletion bubbles per sample', 'HPRC2+HGSVC3 in AoU+HPRC2+HGSVC3', f'{output_prefix}-panel-het-SV-del', [0, 500])
        plot_boxplot(target_results_df, 'Heterozygous SV-length deletion bubbles per sample', 'Target', f'{output_prefix}-target-het-SV-del', [0, 500])

        # 7. ALT Length Weighted Histogram
        print("Generating ALT Length Histogram...", flush=True)
        bins = list(np.linspace(-10000, -100, 397)) + [-75, -50, -25, -1, -0.1, 0.1, 1, 25, 50, 75] + list(np.linspace(100, 10000, 397))
        plt.figure()
        plt.hist(altlen, bins=bins, label='AoU+HPRC2+HGSVC3', log=True, histtype='step', alpha=0.5, weights=panel_mean_alt_alleles_all)
        plt.hist(altlen, bins=bins, label='Target', log=True, histtype='step', alpha=0.5, weights=target_mean_alt_alleles_all)
        plt.xscale('symlog', linthresh=100)
        plt.xticks([-1E4, -1E3] + list(np.linspace(-100, 100, 9)) + [1E3, 1E4],
                   labels=['$-10^4$', '$-10^3$', '$-10^2$'] + ['', '$-50$', '', '$0$', '', '$50$', ''] + ['$10^2$', '$10^3$', '$10^4$'])
        plt.ylabel('Number of ALT alleles per sample')
        plt.xlabel('ALT length - REF length (bp)')
        plt.legend()
        plt.savefig(f'{output_prefix}-alt-alleles-per-sample-hist.png', bbox_inches='tight')
        plt.savefig(f'{output_prefix}-alt-alleles-per-sample-hist.pdf', bbox_inches='tight')
        plt.close()

        # 8. De Finetti Plots
        print("Generating De Finetti Plots...", flush=True)
        ternary_to_cartesian = lambda a, b, c: (0.5 * (2 * b + c) / (a + b + c + 1E-10), 0.5 * np.sqrt(3) * c / (a + b + c + 1E-10))

        def calc_hwe_ternary(x, m=1, f=1):
            return np.array([1 - f, 0, 0]) + f * np.array([m * (1 - x)**2 + (1 - m) * (1 - x), m * x**2 + (1 - m) * x, 2 * m * x * (1 - x)])

        def make_de_finetti_ax():
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.set_xlim([-0.1, 1.1])
            ax.set_ylim([-0.1, np.sqrt(3) / 2 + 0.1])
            ax.set_aspect(np.sqrt(3) / 2)
            ax.axis('off')
            ep = 0.02
            ax.plot([-2 * ep / np.sqrt(3), 1 + 2 * ep / np.sqrt(3)], [-ep, -ep], lw=3, c='k')
            ax.plot([-2 * ep / np.sqrt(3), 0.5], [-ep, 0.5 * np.sqrt(3) + ep], lw=3, c='k')
            ax.plot([1 + 2 * ep / np.sqrt(3), 0.5], [-ep, 0.5 * np.sqrt(3) + ep], lw=3, c='k')
            ax.text(-4 * ep / np.sqrt(3), -4 * ep, 'HOM\nREF ', fontsize=16, ha='right')
            ax.text(0.5, 0.5 * np.sqrt(3) + 3 * ep, 'HET', fontsize=16, ha='center')
            ax.text(1 + 4 * ep / np.sqrt(3), -4 * ep, 'HOM\n ALT', fontsize=16, ha='left')
            return fig, ax

        def plot_de_finetti(hom_ref_arr, hom_alt_arr, het_arr, title, outfile, gridsize=70):
            if len(hom_ref_arr) == 0: return
            fig, ax = make_de_finetti_ax()
            
            # 1. Combine arrays into a 2D matrix
            counts_arr = np.column_stack((hom_ref_arr, hom_alt_arr, het_arr))
            
            # 2. Pre-aggregate identical variant counts
            unique_counts, point_weights = np.unique(counts_arr, axis=0, return_counts=True)
            
            # 3. Calculate X/Y coordinates ONLY for unique combinations
            x_ternary_v, y_ternary_v = ternary_to_cartesian(unique_counts[:,0], unique_counts[:,1], unique_counts[:,2])
            
            # 4. Pass the point_weights directly to hexbin using 'C'
            hb = ax.hexbin(
                x_ternary_v, 
                y_ternary_v, 
                C=point_weights, 
                reduce_C_function=np.sum, 
                gridsize=gridsize, 
                extent=[0, 1, 0, np.sqrt(3) / 2], 
                norm=matplotlib.colors.LogNorm()
            )
            
            cbar = plt.colorbar(hb, ax=ax, shrink=0.5)
            cbar.set_label('Number of variants', rotation=270, labelpad=10)
            
            x_values = np.linspace(0, 1, 50)
            cart_values = np.array([ternary_to_cartesian(*calc_hwe_ternary(x)) for x in x_values])
            ax.plot(cart_values[:, 0], cart_values[:, 1], c='C1', ls='solid', lw=3)
            
            ax.text(0.5, -0.3, title, fontsize=18, ha='center')
            plt.savefig(f'{outfile}.png', bbox_inches='tight')
            plt.savefig(f'{outfile}.pdf', bbox_inches='tight')
            plt.close()

        # Slice HWE dictionaries for subsets
        def slice_hwe(hwe_dict, mask):
            mask = np.array(mask)
            return {k: np.array(v)[mask] for k, v in hwe_dict.items()}

        p_hwe_ins = slice_hwe(panel_hwe, is_sv_ins)
        p_hwe_del = slice_hwe(panel_hwe, is_sv_del)
        c_hwe_ins = slice_hwe(target_hwe, is_sv_ins)
        c_hwe_del = slice_hwe(target_hwe, is_sv_del)

        plot_de_finetti(panel_hwe['hom_ref'], panel_hwe['hom_alt'], panel_hwe['het'], 
                        'AoU+HPRC2+HGSVC3\nAll variants', f'{output_prefix}-panel-hwe-all')
        plot_de_finetti(p_hwe_ins['hom_ref'], p_hwe_ins['hom_alt'], p_hwe_ins['het'], 
                        'AoU+HPRC2+HGSVC3\nSV-length insertion bubbles', f'{output_prefix}-panel-hwe-SV-ins')
        plot_de_finetti(p_hwe_del['hom_ref'], p_hwe_del['hom_alt'], p_hwe_del['het'], 
                        'AoU+HPRC2+HGSVC3\nSV-length deletion bubbles', f'{output_prefix}-panel-hwe-SV-del')

        plot_de_finetti(target_hwe['hom_ref'], target_hwe['hom_alt'], target_hwe['het'], 
                        'Target\nAll variants', f'{output_prefix}-target-hwe-all', gridsize=80)
        plot_de_finetti(c_hwe_ins['hom_ref'], c_hwe_ins['hom_alt'], c_hwe_ins['het'], 
                        'Target\nSV-length insertion bubbles', f'{output_prefix}-target-hwe-SV-ins', gridsize=80)
        plot_de_finetti(c_hwe_del['hom_ref'], c_hwe_del['hom_alt'], c_hwe_del['het'], 
                        'Target\nSV-length deletion bubbles', f'{output_prefix}-target-hwe-SV-del', gridsize=80)

        print("All tasks completed successfully.", flush=True)
        EOF

        python3 summarize.py "~{panel_vcf}##idx##~{panel_vcf_idx}" "~{imputed_vcf}##idx##~{imputed_vcf_idx}" "~{population_tsv}" "~{output_prefix}" ~{chunk_size}
    >>>

    output {
        File summarize_pearson_tsv = "~{output_prefix}.pearson.tsv"
        Array[File] summarize_plots_png = glob("*.png")
        Array[File] summarize_plots_pdf = glob("*.pdf")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             16,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "SSD",
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "continuumio/miniconda3:latest"
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
