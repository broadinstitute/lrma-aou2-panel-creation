version 1.0

import "../MultilevelHierarchicallyPasteVcfsStreaming.wdl" as Paste

# CohortInfoFromBubbles
#
# One chromosome of the batched imputation: re-pop every batch's GLIMPSE2 bubble posteriors with
# the fixed pop-glimpse2, pool the per-batch AF/INFO into one cohort INFO, and paste the batches
# into one cohort BCF that carries it.
#
# The fixed pop-glimpse2 writes the same GT/DS/GP as the one that made the deliverables and
# computes AF/INFO from the bubble paths instead of copying the first path's. So the re-pop is kept
# only as two small tables per batch (AF/INFO, and the site list), and RepopBatch refuses unless
# CHROM, POS, ID, REF, ALT and every sample's GT:DS:GP match the existing popped BCF byte for byte.
# The merged genotypes are the existing popped BCFs pasted, never re-encoded.
#
# The paste copies INFO from its first input at every level, so only the first batch's records are
# rewritten (PoolCohortInfo) and that file is pasted in front of the other batches.
#
# Inputs never get written to. Every output goes to the run's own output path.
workflow CohortInfoFromBubbles {
    input {
        String chromosome
        # One row per batch, first batch first: batch name, bubble BCF, its index, popped BCF, its index.
        # The popped BCFs are the paste inputs and the genotype reference for the re-pop check.
        File batches_tsv
        File pop_glimpse2_panel_resources_json
        # pop-glimpse2 from broadinstitute/lrma-sv-imputation-utils main (PR #11)
        File pop_glimpse2_binary
        File cohort_info_script
        File paste_vcfs_binary
        Int expected_samples = 250
        # Paste regions for this chromosome; each is pasted separately and the regions concatenated.
        Array[String] paste_regions
        Array[Int] paste_batch_sizes = [20, 10]
        Array[Int] paste_timeouts_min = [720, 720]
        # Level 0 slices each popped BCF by region straight from GCS. Local test runs localize instead.
        Array[Boolean] paste_do_localization = [false, true, true]
        # CheckMerged streams the merged BCF from GCS, so a local test run checks it in its own script
        Boolean check_merged = true
        # Paste VM settings per level, e.g. {"preemptible_tries": 3}; unset fields keep MergeVcfs' defaults (including its
        # input-sized disk). Required: Cromwell 91 cannot evaluate these optional inputs of the paste WDL when a calling
        # workflow leaves them unset (its nested scatters become an implicit subworkflow).
        RuntimeAttr paste_l0_runtime_attr
        RuntimeAttr paste_l1_runtime_attr
        RuntimeAttr paste_l2_runtime_attr
        String output_prefix = "aou2_50k.cohort_info"
        String docker = "us.gcr.io/broad-dsde-methods/slee/lrma-aou2-panel-creation-rust@sha256:0f25c4091c49d8eb0c3d8bcdb45e7680a093e0a691e74ccec2d3743758a7d22c"
        String zones = "us-central1-a"
    }

    Array[Array[String]] batches = read_tsv(batches_tsv)
    Map[String, PanelResources] panel_by_chrom = read_json(pop_glimpse2_panel_resources_json)
    PanelResources panel = panel_by_chrom[chromosome]

    scatter (b in batches) {
        call RepopBatch { input:
            batch = b[0],
            bubble_bcf = b[1],
            bubble_bcf_idx = b[2],
            popped_bcf = b[3],
            popped_bcf_idx = b[4],
            panel_bubble_split_sites_only_vcf = panel.panel_bubble_split_sites_only_vcf,
            panel_bubble_split_sites_only_vcf_idx = panel.panel_bubble_split_sites_only_vcf_idx,
            panel_id_split_vcf_gz = panel.panel_id_split_vcf_gz,
            panel_id_split_vcf_gz_tbi = panel.panel_id_split_vcf_gz_tbi,
            pop_regions = select_first([panel.pop_regions]),
            pop_glimpse2_binary = pop_glimpse2_binary,
            expected_samples = expected_samples,
            docker = docker,
            zones = zones
        }
    }

    call PoolCohortInfo { input:
        base_bcf = batches[0][3],
        base_bcf_idx = batches[0][4],
        sites_tsv = RepopBatch.sites_tsv[0],
        sites_md5s = RepopBatch.sites_md5,
        info_tsvs = RepopBatch.info_tsv,
        n_samples = RepopBatch.n_samples,
        cohort_info_script = cohort_info_script,
        output_prefix = output_prefix + "." + chromosome + ".base",
        docker = docker,
        zones = zones
    }

    scatter (i in range(length(batches))) {
        String paste_vcf = if i == 0 then PoolCohortInfo.base_cohort_bcf else batches[i][3]
        String paste_idx = if i == 0 then PoolCohortInfo.base_cohort_bcf_idx else batches[i][4]
    }

    call Paste.HierarchicallyMergeVcfs as PasteBatches { input:
        vcfs_array = paste_vcf,
        vcf_idxs_array = paste_idx,
        regions = paste_regions,
        batch_sizes = paste_batch_sizes,
        do_localization = paste_do_localization,
        timeouts_min = paste_timeouts_min,
        output_prefix = output_prefix + "." + chromosome,
        paste_vcfs_binary = paste_vcfs_binary,
        l0_merge_runtime_attr_override = paste_l0_runtime_attr,
        l1_merge_runtime_attr_override = paste_l1_runtime_attr,
        l2_merge_runtime_attr_override = paste_l2_runtime_attr,
        extra_merge_args = "--threads $(nproc) --info ID,RAF,AF,INFO,N_PATHS,N_PATHS_TOTAL,DS_SUM,VAR_SUM,N_HAPS --format GT,DS,GP"
    }

    if (check_merged) {
        call CheckMerged { input:
            merged_bcf = PasteBatches.merged_vcf,
            merged_bcf_idx = PasteBatches.merged_vcf_idx,
            n_records = PoolCohortInfo.n_records,
            n_samples = RepopBatch.n_samples,
            sites_info_md5 = PoolCohortInfo.sites_info_md5,
            docker = docker,
            zones = zones
        }
    }

    output {
        Array[File] repop_reports = RepopBatch.report
        File cohort_report = PoolCohortInfo.report
        File base_cohort_bcf = PoolCohortInfo.base_cohort_bcf
        File merged_bcf = PasteBatches.merged_vcf
        File merged_bcf_idx = PasteBatches.merged_vcf_idx
        File? merged_report = CheckMerged.report
    }
}

struct PanelResources {
    String panel_bubble_split_sites_only_vcf
    String panel_bubble_split_sites_only_vcf_idx
    String panel_id_split_vcf_gz
    String panel_id_split_vcf_gz_tbi
    Array[String]? pop_regions
}

# Pops one batch-chromosome region by region, as GLIMPSE2FromPreprocessedPLsJoint did (same regions,
# --regions-overlap 0, per-region sort, naive concat), but with the regions run in parallel on one
# VM instead of one VM each.
task RepopBatch {
    input {
        String batch
        File bubble_bcf
        File bubble_bcf_idx
        File popped_bcf
        File popped_bcf_idx
        File panel_bubble_split_sites_only_vcf
        File panel_bubble_split_sites_only_vcf_idx
        File panel_id_split_vcf_gz
        File panel_id_split_vcf_gz_tbi
        Array[String] pop_regions
        File pop_glimpse2_binary
        Int expected_samples
        String docker
        String zones
        Int cpu = 8
        # 4.7 GiB measured peak per region plus bcftools sort --max-mem=2G, two cpu per region
        Int mem_gb = 32
    }

    Int jobs = cpu / 2
    Int disk_gb = 20 + ceil(size(bubble_bcf, "GB") + 3 * size(popped_bcf, "GB") + size([panel_bubble_split_sites_only_vcf, panel_id_split_vcf_gz], "GB"))

    command <<<
        set -euo pipefail
        cp ~{pop_glimpse2_binary} pop && chmod +x pop

        bcftools query -l ~{bubble_bcf} > samples.bubble
        bcftools query -l ~{popped_bcf} > samples.popped
        n=$(awk 'END{print NR}' samples.popped)
        [ "$n" = ~{expected_samples} ] || { echo "ERROR: ~{batch} popped has $n samples" >&2; exit 1; }
        cmp -s samples.bubble samples.popped || { echo "ERROR: ~{batch} bubble and popped samples differ" >&2; exit 1; }

        pop_region() {
            set -euo pipefail
            k=$1; r=$2
            bcftools view -r "$r" --regions-overlap 0 ~{panel_bubble_split_sites_only_vcf} -Oz -o sites.$k.vcf.gz
            bcftools view -r "$r" --regions-overlap 0 ~{bubble_bcf} | \
                ./pop ~{panel_id_split_vcf_gz} sites.$k.vcf.gz | \
                bcftools sort --max-mem=2G -T sort.$k -Ob -o shard.$k.bcf
            rm sites.$k.vcf.gz
        }
        export -f pop_region
        awk '{printf "%05d\t%s\n", NR, $0}' ~{write_lines(pop_regions)} > regions.tsv
        xargs -P ~{jobs} -L 1 bash -c 'pop_region "$@"' _ < regions.tsv
        cut -f1 regions.tsv | awk '{print "shard." $1 ".bcf"}' > shards.txt
        bcftools concat --naive -f shards.txt -Ob -o repop.bcf
        bcftools index repop.bcf
        xargs rm < shards.txt

        r_old=$(bcftools index -n ~{popped_bcf}); r_new=$(bcftools index -n repop.bcf)
        [ "$r_old" = "$r_new" ] || { echo "ERROR: ~{batch} records $r_old popped, $r_new re-popped" >&2; exit 1; }

        # Only INFO may differ from the deliverable.
        q='%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT:%DS:%GP]\n'
        bcftools query -f "$q" ~{popped_bcf} | md5sum | cut -d' ' -f1 > md5.old &
        bcftools query -f "$q" repop.bcf | md5sum | cut -d' ' -f1 > md5.new &
        bcftools query -f '%INFO/AF\t%INFO/INFO\n' repop.bcf > info.new &
        bcftools query -f '%INFO/INFO\n' ~{popped_bcf} > info.old &
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/ID\t%INFO/RAF\t%INFO/N_PATHS\t%INFO/N_PATHS_TOTAL\n' repop.bcf > sites.tsv &
        wait
        cmp -s md5.old md5.new || { echo "ERROR: ~{batch} sites or genotypes differ from the popped BCF" >&2; exit 1; }

        paste info.new info.old | gzip -1 > info.tsv.gz
        md5sum < sites.tsv | cut -d' ' -f1 > sites.md5
        gzip -1 sites.tsv
        echo "$n" > n_samples
        changed=$(paste info.new info.old | awk '$3 != "." {d=$2-$3; if (d<0) d=-d; if (d>0.1) c++} END{print c+0}')
        old_zero=$(awk '$1 == 0' info.old | awk 'END{print NR}')
        printf '%s samples=%s records=%s genotypes_md5_equal=yes info_changed_gt_0.1=%s old_info_zero=%s\n' \
            "~{batch}" "$n" "$r_new" "$changed" "$old_zero" | tee report.txt
    >>>

    output {
        File info_tsv = "info.tsv.gz"
        File sites_tsv = "sites.tsv.gz"
        String sites_md5 = read_string("sites.md5")
        Int n_samples = read_int("n_samples")
        File report = "report.txt"
    }
    runtime {
        cpu: cpu
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " SSD"
        bootDiskSizeGb: -20
        preemptible: 3
        maxRetries: 1
        docker: docker
        cpuPlatform: "AMD Rome"
        zones: zones
    }
}

task PoolCohortInfo {
    input {
        File base_bcf
        File base_bcf_idx
        File sites_tsv
        Array[String] sites_md5s
        Array[File] info_tsvs
        Array[Int] n_samples
        File cohort_info_script
        String output_prefix
        String docker
        String zones
    }

    Int disk_gb = 20 + ceil(3 * size(base_bcf, "GB") + size(info_tsvs, "GB"))

    command <<<
        set -euo pipefail
        n_lists=$(sort -u ~{write_lines(sites_md5s)} | awk 'END{print NR}')
        [ "$n_lists" = 1 ] || { echo "ERROR: batches have $n_lists different site lists" >&2; exit 1; }

        bcftools view ~{base_bcf} | \
            python3 ~{cohort_info_script} --n-samples ~{sep="," n_samples} --sites ~{sites_tsv} --report report.txt \
                ~{sep=" " info_tsvs} | \
            bcftools view -Ob -o ~{output_prefix}.bcf --write-index

        r_in=$(bcftools index -n ~{base_bcf}); r_out=$(bcftools index -n ~{output_prefix}.bcf)
        [ "$r_in" = "$r_out" ] || { echo "ERROR: records $r_in -> $r_out" >&2; exit 1; }
        q='%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT:%DS:%GP]\n'
        m_in=$(bcftools query -f "$q" ~{base_bcf} | md5sum | cut -d' ' -f1)
        m_out=$(bcftools query -f "$q" ~{output_prefix}.bcf | md5sum | cut -d' ' -f1)
        [ "$m_in" = "$m_out" ] || { echo "ERROR: sites or genotypes changed" >&2; exit 1; }
        echo "$r_out" > n_records
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/ID\t%INFO/RAF\t%INFO/AF\t%INFO/INFO\t%INFO/N_PATHS\t%INFO/N_PATHS_TOTAL\t%INFO/DS_SUM\t%INFO/VAR_SUM\t%INFO/N_HAPS\n' ~{output_prefix}.bcf | md5sum | cut -d' ' -f1 > sites_info.md5
        cat report.txt
    >>>

    output {
        File base_cohort_bcf = "~{output_prefix}.bcf"
        File base_cohort_bcf_idx = "~{output_prefix}.bcf.csi"
        Int n_records = read_int("n_records")
        String sites_info_md5 = read_string("sites_info.md5")
        File report = "report.txt"
    }
    runtime {
        cpu: 4
        memory: "8 GiB"
        disks: "local-disk " + disk_gb + " SSD"
        bootDiskSizeGb: -20
        preemptible: 2
        maxRetries: 1
        docker: docker
        zones: zones
    }
}

# Streams the merged BCF instead of localizing it.
task CheckMerged {
    input {
        String merged_bcf
        String merged_bcf_idx
        Int n_records
        Array[Int] n_samples
        String sites_info_md5
        String docker
        String zones
    }

    command <<<
        set -euo pipefail
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        f="~{merged_bcf}##idx##~{merged_bcf_idx}"
        want_s=$(echo ~{sep=" " n_samples} | tr ' ' '\n' | awk '{s+=$1} END{print s}')
        n_s=$(bcftools query -l "$f" | awk 'END{print NR}')
        n_u=$(bcftools query -l "$f" | sort -u | awk 'END{print NR}')
        [ "$n_s" = "$want_s" ] && [ "$n_u" = "$n_s" ] || { echo "ERROR: $n_s samples ($n_u unique), expected $want_s" >&2; exit 1; }
        n_r=$(bcftools index -n "$f")
        [ "$n_r" = ~{n_records} ] || { echo "ERROR: $n_r records, expected ~{n_records}" >&2; exit 1; }
        m=$(bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/ID\t%INFO/RAF\t%INFO/AF\t%INFO/INFO\t%INFO/N_PATHS\t%INFO/N_PATHS_TOTAL\t%INFO/DS_SUM\t%INFO/VAR_SUM\t%INFO/N_HAPS\n' "$f" | md5sum | cut -d' ' -f1)
        [ "$m" = ~{sites_info_md5} ] || { echo "ERROR: merged sites or INFO differ from the cohort base" >&2; exit 1; }
        printf 'samples=%s records=%s sites_and_info_equal_to_base=yes\n' "$n_s" "$n_r" | tee report.txt
    >>>

    output { File report = "report.txt" }
    runtime {
        cpu: 2
        memory: "4 GiB"
        disks: "local-disk 20 HDD"
        bootDiskSizeGb: -20
        preemptible: 2
        maxRetries: 1
        docker: docker
        zones: zones
    }
}
