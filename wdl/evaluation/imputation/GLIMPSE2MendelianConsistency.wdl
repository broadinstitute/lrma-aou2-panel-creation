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

workflow MendelianConsistency {
    input {
        File panel_sites_only_vcf       # split to biallelic
        File panel_sites_only_vcf_idx
        File imputed_vcf                # split to biallelic, variants in same order as in panel
        File imputed_vcf_idx
        File trh_bed
        File trh_bed_idx
        File pedigree
        String output_prefix
    }

    call AnnotateVcf { input:
        panel_sites_only_vcf = panel_sites_only_vcf,
        panel_sites_only_vcf_idx = panel_sites_only_vcf_idx,
        imputed_vcf = imputed_vcf,
        imputed_vcf_idx = imputed_vcf_idx,
        trh_bed = trh_bed,
        trh_bed_idx = trh_bed_idx,
        output_prefix = output_prefix
    }

    call CalculateMendelianConsistency { input:
        annotated_vcf = AnnotateVcf.annotated_vcf,
        annotated_vcf_idx = AnnotateVcf.annotated_vcf_idx,
        pedigree = pedigree,
        output_prefix = output_prefix
    }

    output {
        Array[File] mendelian_pkls = [CalculateMendelianConsistency.unfiltered_pkl, CalculateMendelianConsistency.filtered_pkl]
        Array[File] mendelian_plots = [CalculateMendelianConsistency.trio_plot_inTRH, CalculateMendelianConsistency.trio_plot_outTRH, 
                                       CalculateMendelianConsistency.locus_plot_inTRH, CalculateMendelianConsistency.locus_plot_outTRH]
    }
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

    Int disk_gb = 20 + ceil(size(panel_sites_only_vcf, "GiB") + size(imputed_vcf, "GiB") * 3)

    command <<<
        set -euxo pipefail

        # 1. Annotate Panel VCF with TRH
        echo "Annotating Panel VCF with TRH regions..."
        bcftools annotate -a ~{trh_bed}##idx##~{trh_bed_idx} -c CHROM,FROM,TO -m +TRH --threads $(nproc) \
            ~{panel_sites_only_vcf}##idx##~{panel_sites_only_vcf_idx} -W -Ob -o panel.trh.bcf

        # 2. Transfer panel AF and TRH tags to Imputed VCF 
        echo "Transferring panel AF and TRH annotations to Imputed VCF..."
        bcftools annotate -a panel.trh.bcf -c INFO/AF,INFO/TRH --threads $(nproc) \
            ~{imputed_vcf}##idx##~{imputed_vcf_idx} -W -Ob -o ~{output_prefix}.annotated.bcf
    >>>

    output {
        File annotated_vcf = "~{output_prefix}.annotated.bcf"
        File annotated_vcf_idx = "~{output_prefix}.annotated.bcf.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        disk_type:          "SSD",
        preemptible_tries:  1,
        max_retries:        0,
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

task CalculateMendelianConsistency {
    input {
        File annotated_vcf
        File annotated_vcf_idx
        File pedigree
        String output_prefix
        
        Int chunk_size = 10000

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = ceil(size(annotated_vcf, "GiB")) + 20

    command <<<
        set -euxo pipefail

        # Install the dependencies for variant streaming and plotting (if using standard miniconda)
        conda install -y -c bioconda -c conda-forge cyvcf2 pandas numpy matplotlib seaborn

        cat << 'EOF' > mendelian.py
        import sys
        import time
        import numpy as np
        import pandas as pd
        import cyvcf2
        import matplotlib.pyplot as plt
        import seaborn as sns
        import warnings

        # Ignore division by zero warnings in visualization calculations
        warnings.filterwarnings('ignore', category=RuntimeWarning)

        vcf_path = sys.argv[1]
        ped_path = sys.argv[2]
        output_prefix = sys.argv[3]
        chunk_size = int(sys.argv[4])

        vcf = cyvcf2.VCF(vcf_path)
        ped_df = pd.read_csv(ped_path, sep='\t', names=['familyID', 'sampleID', 'fatherID', 'motherID', 'sex', 'status'])
        samples = np.array(vcf.samples)

        is_trio_complete = []
        valid_trios = []
        for _, row in ped_df.iterrows():
            if row['sampleID'] in samples and row['fatherID'] in samples and row['motherID'] in samples:
                is_trio_complete.append(True)
                valid_trios.append({
                    's': np.where(samples == row['sampleID'])[0][0],
                    'f': np.where(samples == row['fatherID'])[0][0],
                    'm': np.where(samples == row['motherID'])[0][0]
                })
            else:
                is_trio_complete.append(False)

        num_complete_trios = len(valid_trios)
        s_idx = [t['s'] for t in valid_trios]
        f_idx = [t['f'] for t in valid_trios]
        m_idx = [t['m'] for t in valid_trios]

        groups = {}

        def get_bins(chrom, ref, alt, af, trh):
            if len(ref) == 1 and len(alt) == 1:
                length_bin = 'SNP'
            else:
                l = len(alt) - len(ref)
                if l <= -50: length_bin = '(-inf, -50]'
                elif -50 < l <= -1: length_bin = '(-50, -1]'
                elif 1 <= l < 50: length_bin = '[1, 50)'
                else: length_bin = '[50, inf)'

            if af < 0.01: af_bin = '[0, 0.01)'
            elif af < 0.1: af_bin = '[0.01, 0.1)'
            else: af_bin = '[0.1, 1]'

            return chrom, af_bin, length_bin, bool(trh)

        def process_chunk(gt_arr, gp_arr, meta_list):
            p_gt = gt_arr[:, s_idx]
            f_gt = gt_arr[:, f_idx]
            m_gt = gt_arr[:, m_idx]

            gp_mask = (gp_arr >= 0.9)
            p_gp_mask = gp_mask[:, s_idx]
            f_gp_mask = gp_mask[:, f_idx]
            m_gp_mask = gp_mask[:, m_idx]

            for gp_threshold in [0, 0.9]:
                if gp_threshold > 0:
                    p_gt_filt = np.where(p_gp_mask, p_gt, 2)
                    f_gt_filt = np.where(f_gp_mask, f_gt, 2)
                    m_gt_filt = np.where(m_gp_mask, m_gt, 2)
                else:
                    p_gt_filt = p_gt
                    f_gt_filt = f_gt
                    m_gt_filt = m_gt

                miss = (p_gt_filt == 2) | (f_gt_filt == 2) | (m_gt_filt == 2)

                err_p0 = (p_gt_filt == 0) & ((f_gt_filt == 3) | (m_gt_filt == 3))
                err_p3 = (p_gt_filt == 3) & ((f_gt_filt == 0) | (m_gt_filt == 0))
                err_p1 = (p_gt_filt == 1) & (((f_gt_filt == 0) & (m_gt_filt == 0)) | ((f_gt_filt == 3) & (m_gt_filt == 3)))

                errors = (err_p0 | err_p3 | err_p1) & ~miss

                is_hom_ref = (p_gt_filt == 0) & (f_gt_filt == 0) & (m_gt_filt == 0)
                non_hom_ref = ~is_hom_ref & ~miss

                # Memory Efficient Aggregation: Calculate rates for plotting on-the-fly
                for i, meta in enumerate(meta_list):
                    key = meta
                    if key not in groups:
                        groups[key] = {
                            0: {'trio_err': np.zeros(num_complete_trios), 'trio_nhr': np.zeros(num_complete_trios), 'locus_rates': [], 'locus_nhr': []}, 
                            0.9: {'trio_err': np.zeros(num_complete_trios), 'trio_nhr': np.zeros(num_complete_trios), 'locus_rates': [], 'locus_nhr': []}, 
                            'count': 0
                        }

                    if gp_threshold == 0:
                        groups[key]['count'] += 1

                    err_v = errors[i] & non_hom_ref[i]
                    nhr_v = non_hom_ref[i]

                    # Accumulate for Per Trio Plot (axis=0 sums)
                    groups[key][gp_threshold]['trio_err'] += err_v
                    groups[key][gp_threshold]['trio_nhr'] += nhr_v

                    # Accumulate for Per Locus Plot (axis=1 sums)
                    locus_nhr_sum = nhr_v.sum()
                    if locus_nhr_sum > 0:
                        locus_err_sum = err_v.sum()
                        groups[key][gp_threshold]['locus_rates'].append(locus_err_sum / locus_nhr_sum)
                        groups[key][gp_threshold]['locus_nhr'].append(locus_nhr_sum)

        print(f"Streaming variants with chunk size {chunk_size}...", flush=True)
        gt_chunk, gp_chunk, meta_chunk = [], [], []
        variants_processed = 0
        start_time = time.time()

        for variant in vcf:
            chrom = variant.CHROM
            ref = variant.REF
            alts = variant.ALT
            if not alts: continue
            alt = alts[0]

            af = variant.INFO.get('AF')
            if af is None: af = 0.0
            trh = variant.INFO.get('TRH')
            if trh is None: trh = False

            meta = get_bins(chrom, ref, alt, af, trh)
            gt_chunk.append(variant.gt_types)

            gp = variant.format('GP')
            if gp is not None:
                gp_max = np.max(gp, axis=1)
                gp_chunk.append(gp_max)
            else:
                raise ValueError("FORMAT/GP is missing from the annotated VCF! Cannot build step plots.")

            meta_chunk.append(meta)

            if len(gt_chunk) >= chunk_size:
                process_chunk(np.array(gt_chunk), np.array(gp_chunk), meta_chunk)
                variants_processed += len(gt_chunk)
                elapsed = time.time() - start_time
                print(f"Processed {variants_processed} variants in {elapsed:.2f}s ({variants_processed/elapsed:.2f} var/s)...", flush=True)
                gt_chunk, gp_chunk, meta_chunk = [], [], []

        if len(gt_chunk) > 0:
            process_chunk(np.array(gt_chunk), np.array(gp_chunk), meta_chunk)
            variants_processed += len(gt_chunk)
            elapsed = time.time() - start_time
            print(f"Finished processing {variants_processed} total variants in {elapsed:.2f}s.", flush=True)

        print("Assembling tables...", flush=True)
        mv_df_columns = ['CHROM', 'AF_BIN', 'LENGTH_BIN', 'IN_TRH', 'TRIO_ERR', 'TRIO_NHR', 'LOCUS_RATES', 'LOCUS_NHR', 'NUM_LOCI']
        mv_df_0, mv_df_09 = [], []
        for key, data in groups.items():
            chrom, af_bin, length_bin, in_trh = key

            err_0 = data[0]['trio_err']
            nhr_0 = data[0]['trio_nhr']
            l_rates_0 = np.array(data[0]['locus_rates'])
            l_nhr_0 = np.array(data[0]['locus_nhr'])
            mv_df_0.append([chrom, af_bin, length_bin, in_trh, err_0, nhr_0, l_rates_0, l_nhr_0, data['count']])

            err_09 = data[0.9]['trio_err']
            nhr_09 = data[0.9]['trio_nhr']
            l_rates_09 = np.array(data[0.9]['locus_rates'])
            l_nhr_09 = np.array(data[0.9]['locus_nhr'])
            mv_df_09.append([chrom, af_bin, length_bin, in_trh, err_09, nhr_09, l_rates_09, l_nhr_09, data['count']])

        df_0 = pd.DataFrame(mv_df_0, columns=mv_df_columns)
        df_09 = pd.DataFrame(mv_df_09, columns=mv_df_columns)

        df_0.to_pickle(f'{output_prefix}-unfiltered.pkl')
        df_09.to_pickle(f'{output_prefix}-filter-0.9.pkl')

        print("Generating plots...", flush=True)
        length_bin_labels = ['(-inf, -50]', '(-50, -1]', 'SNP', '[1, 50)', '[50, inf)']
        af_bin_labels = ['[0, 0.01)', '[0.01, 0.1)', '[0.1, 1]']

        # PER TRIO PLOTTING
        for in_trh in [False, True]:
            fig, ax = plt.subplots(3, 1, figsize=(9, 10))
            trh_tag = 'all' if in_trh is None else ('non-TR/homopolymer' if not in_trh else 'TR/homopolymer')
            ax[0].set_title(f'{num_complete_trios} trios, {trh_tag}') 

            for i, af_bin_label in enumerate(af_bin_labels):
                plt_df_values = []
                for j, length_bin_label in enumerate(length_bin_labels):
                    for g, min_tar_gp in enumerate([0, 0.9]):
                        min_tar_gp_label = 'unfiltered' if min_tar_gp == 0 else f'GP > {min_tar_gp}'
                        trio_err_values = []
                        trio_nhr_values = []

                        df = df_0 if min_tar_gp == 0 else df_09
                        y_trh = (df['IN_TRH'] == in_trh) if in_trh is not None else True
                        bin_df = df[y_trh & (df['AF_BIN'] == af_bin_label) & (df['LENGTH_BIN'] == length_bin_label)]

                        if not bin_df.empty:
                            trio_err_values.append(bin_df['TRIO_ERR'].values[0])
                            trio_nhr_values.append(bin_df['TRIO_NHR'].values[0])

                        if not trio_err_values: continue

                        trio_err = np.sum(trio_err_values, axis=0)
                        trio_nhr = np.sum(trio_nhr_values, axis=0)

                        error_rates_per_trio = trio_err / trio_nhr
                        mean_num_non_hom_ref = trio_nhr.mean()

                        plt_df_values.extend(
                            [[min_tar_gp_label, 
                              g * '\n' + f'$\\langle N_{{l}} \\rangle={mean_num_non_hom_ref:.2f}$' + ('\n\n\n' + length_bin_label if g == 1 else ''), 
                              error_rates_per_trio[t]] 
                             for t in range(num_complete_trios)])

                if not plt_df_values: continue
                plt_df = pd.DataFrame(plt_df_values, columns=['MIN_TAR_GP_TEXT', 'LENGTH_BIN_TEXT', 'ERROR_RATE'])
                sns.boxplot(data=plt_df, x='LENGTH_BIN_TEXT', y='ERROR_RATE', hue='MIN_TAR_GP_TEXT', ax=ax[i], legend=(i==2))

                ax[i].set_xlabel('ALT length - REF length (bp)' if i == len(af_bin_labels) - 1 else None)
                ax[i].set_ylabel(('panel allele frequency\n' if i == 1 else '\n\n')
                                 + f'{af_bin_label}\n\n' + ('Mendelian error rate per trio' if i == 1 else ''))
                ax[i].set_yscale('symlog', linthresh=0.001)
                ax[i].set_ylim([-1E-4, 1.001])
                ax[i].set_yticks([k * 0.0001 for k in range(0, 10)] + 
                                 [k * 0.001 for k in range(0, 10)] + 
                                 [k * 0.01 for k in range(1, 10)] + 
                                 [k * 0.1 for k in range(1, 11)])
                if i == 2:
                    handles, labels = ax[i].get_legend_handles_labels()
                    ax[i].legend(handles=handles, labels=labels, loc='upper center', fontsize=8)
            plt.tight_layout()
            plt.savefig(f'{output_prefix}.trio.{"inTRH" if in_trh else "outTRH"}.png')

        # PER LOCUS PLOTTING
        for in_trh in [False, True]:
            fig, ax = plt.subplots(3, 1, figsize=(9, 12))
            trh_tag = 'all' if in_trh is None else ('non-TR/homopolymer' if not in_trh else 'TR/homopolymer')
            ax[0].set_title(f'{num_complete_trios} trios, {trh_tag}')
            for i, af_bin_label in enumerate(af_bin_labels):
                plt_df_values = []
                for j, length_bin_label in enumerate(length_bin_labels):
                    for g, min_tar_gp in enumerate([0, 0.9]):
                        min_tar_gp_label = 'unfiltered' if min_tar_gp == 0 else f'GP > {min_tar_gp}'
                        locus_rates_values = []
                        locus_nhr_values = []

                        df = df_0 if min_tar_gp == 0 else df_09
                        y_trh = (df['IN_TRH'] == in_trh) if in_trh is not None else True
                        bin_df = df[y_trh & (df['AF_BIN'] == af_bin_label) & (df['LENGTH_BIN'] == length_bin_label)]

                        if not bin_df.empty:
                            locus_rates_values.append(bin_df['LOCUS_RATES'].values[0])
                            locus_nhr_values.append(bin_df['LOCUS_NHR'].values[0])

                        if not locus_rates_values: continue

                        error_rates_per_locus = np.concatenate(locus_rates_values)
                        valid_loci_nhr = np.concatenate(locus_nhr_values)

                        if len(error_rates_per_locus) == 0: continue

                        mean_num_non_hom_ref_trios = valid_loci_nhr.mean()
                        N_l = len(error_rates_per_locus)

                        length_text = g * '\n\n' + f'$\\langle N_{{t}} \\rangle$={mean_num_non_hom_ref_trios:.2f}' + f'\n$N_{{l}}$={N_l}' + ('\n\n\n\n' + length_bin_label if g == 1 else '')

                        plt_df_values.extend([[min_tar_gp_label, length_text, error_rates_per_locus[v]] for v in range(len(error_rates_per_locus))])

                if not plt_df_values: continue
                plt_df = pd.DataFrame(plt_df_values, columns=['MIN_TAR_GP_TEXT', 'LENGTH_BIN_TEXT', 'ERROR_RATE'])
                sns.boxplot(data=plt_df, x='LENGTH_BIN_TEXT', y='ERROR_RATE', hue='MIN_TAR_GP_TEXT', ax=ax[i], legend=(i==2))

                ax[i].set_xlabel('ALT length - REF length (bp)' if i == len(af_bin_labels) - 1 else None)
                ax[i].set_ylabel(('panel allele frequency\n' if i == 1 else '\n\n') + f'{af_bin_label}\n\n' + ('Mendelian error rate per locus' if i == 1 else ''))
                ax[i].set_yscale('symlog', linthresh=0.001)
                ax[i].set_ylim([-1E-5, 1])
                ax[i].set_yticks([k * 0.0001 for k in range(0, 10)] + [k * 0.001 for k in range(0, 10)] + [k * 0.01 for k in range(1, 10)] + [k * 0.1 for k in range(1, 11)])
                if i == 2:
                    handles, labels = ax[i].get_legend_handles_labels()
                    ax[i].legend(handles=handles, labels=labels, loc='upper center', fontsize=8)
            plt.tight_layout()
            plt.savefig(f'{output_prefix}.locus.{"inTRH" if in_trh else "outTRH"}.png')
        EOF

        python3 mendelian.py "~{annotated_vcf}##idx##~{annotated_vcf_idx}" "~{pedigree}" "~{output_prefix}" "~{chunk_size}"
    >>>

    output {
        File unfiltered_pkl = "~{output_prefix}-unfiltered.pkl"
        File filtered_pkl = "~{output_prefix}-filtered-0.9.pkl"
        File trio_plot_inTRH = "~{output_prefix}.trio.inTRH.png"
        File trio_plot_outTRH = "~{output_prefix}.trio.outTRH.png"
        File locus_plot_inTRH = "~{output_prefix}.locus.inTRH.png"
        File locus_plot_outTRH = "~{output_prefix}.locus.outTRH.png"
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
