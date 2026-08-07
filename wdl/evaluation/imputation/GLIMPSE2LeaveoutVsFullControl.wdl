version 1.0

# Per-window GLIMPSE2_split_reference + GLIMPSE2_phase against the leaveout panel and, as a
# control, the full panel; each arm ligated separately. Self-contained, no imports.
#
# The cohort samples are in the full panel, so the full arm imputes them partly from their own
# haplotypes. It is a ceiling, not a result.
#
# The cohort GL file was built by an inline PL extractor, not extract-bubble-PLs: exact REF/ALT
# matching plus ref-block hom-ref, so SV alleles carry no PL. SV concordance will look bad for
# that reason. (extract-bubble-PLs is in none of lrma-aou2-glimpse2:tfenne-opt,
# lrma-aou2-panel-creation-rust:v1, -python:v1.)
#
# shard-2's output region is chr22:23640007-28828225; the supplied windows cover a 2.4 Mb subset.

workflow GLIMPSE2LeaveoutVsFullControl {
    input {
        File leaveout_panel_bcf
        File leaveout_panel_bcf_idx
        File full_panel_bcf
        File full_panel_bcf_idx

        File cohort_gl_bcf
        File cohort_gl_bcf_idx

        # Parallel arrays: --input-region (buffered) and --output-region. output_regions must
        # tile contiguously or GLIMPSE2_ligate leaves a gap.
        Array[String] input_regions
        Array[String] output_regions

        # split_reference only; phase reads the map out of the .bin.
        File genetic_map

        String output_prefix = "chr22-shard2-leaveout-vs-full"

        Int phase_threads = 8

        # Conditioning states. --pbwt-depth does NOT set this: a run passing --pbwt-depth 1000
        # logged #states=2000.0, the default. Unset means 2000 states and ~2x the memory.
        Int phase_kpbwt = 1000

        # Verbatim from GLIMPSE2FromPreprocessedPLsJoint.wdl.
        String extra_phase_args = "--impute-reference-only-variants --keep-monomorphic-ref-sites --main 10 --burnin 5 --err-imp 1E-3"
        String extra_split_args = "--keep-monomorphic-ref-sites"

        Int split_mem_gb = 16
        Int phase_mem_gb = 32
        Int ligate_threads = 2
        Int ligate_mem_gb = 24

        String docker = "us.gcr.io/broad-dsde-methods/slee/lrma-aou2-glimpse2:tfenne-opt"
    }

    Array[File] panel_bcfs = [leaveout_panel_bcf, full_panel_bcf]
    Array[File] panel_bcf_idxs = [leaveout_panel_bcf_idx, full_panel_bcf_idx]
    Array[String] arm_labels = ["leaveout", "full"]

    # One scatter body for both arms, so only the panel differs between them.
    scatter (a in range(length(panel_bcfs))) {
        scatter (w in range(length(output_regions))) {
            call GLIMPSE2SplitReference {
                input:
                    panel_bcf = panel_bcfs[a],
                    panel_bcf_idx = panel_bcf_idxs[a],
                    input_region = input_regions[w],
                    output_region = output_regions[w],
                    genetic_map = genetic_map,
                    output_prefix = output_prefix + "." + arm_labels[a] + ".window-" + w + ".split",
                    extra_split_args = extra_split_args,
                    mem_gb = split_mem_gb,
                    docker = docker
            }

            call GLIMPSE2Phase {
                input:
                    cohort_gl_bcf = cohort_gl_bcf,
                    cohort_gl_bcf_idx = cohort_gl_bcf_idx,
                    panel_split_bin = GLIMPSE2SplitReference.panel_split_bin,
                    output_prefix = output_prefix + "." + arm_labels[a] + ".window-" + w + ".phased",
                    phase_threads = phase_threads,
                    phase_kpbwt = phase_kpbwt,
                    extra_phase_args = extra_phase_args,
                    mem_gb = phase_mem_gb,
                    docker = docker
            }
        }

        call GLIMPSE2Ligate {
            input:
                phased_bcfs = GLIMPSE2Phase.phased_bcf,
                phased_bcf_idxs = GLIMPSE2Phase.phased_bcf_idx,
                output_prefix = output_prefix + "." + arm_labels[a] + ".ligated",
                ligate_threads = ligate_threads,
                mem_gb = ligate_mem_gb,
                docker = docker
        }
    }

    output {
        # Indices follow panel_bcfs.
        File leaveout_ligated_bcf = GLIMPSE2Ligate.ligated_bcf[0]
        File leaveout_ligated_bcf_idx = GLIMPSE2Ligate.ligated_bcf_idx[0]
        File full_ligated_bcf = GLIMPSE2Ligate.ligated_bcf[1]
        File full_ligated_bcf_idx = GLIMPSE2Ligate.ligated_bcf_idx[1]

        Array[Array[File]] phased_bcfs_per_arm = GLIMPSE2Phase.phased_bcf
    }
}

task GLIMPSE2SplitReference {
    input {
        File panel_bcf
        File panel_bcf_idx
        String input_region
        String output_region
        File genetic_map
        String output_prefix
        String extra_split_args
        Int mem_gb
        String docker
    }

    # Measured on chr22:23800000-26200000 in one go: 60 s, 0.54 GB .bin leaveout / 0.56 GB full.
    Int disk_gb = 20 + 3 * ceil(size(panel_bcf, "GiB"))

    command <<<
        set -euxo pipefail

        /bin/GLIMPSE2_split_reference \
            -R ~{panel_bcf} \
            --input-region ~{input_region} \
            --output-region ~{output_region} \
            --map ~{genetic_map} \
            --thread $(nproc) \
            ~{extra_split_args} \
            --output ~{output_prefix}

        # Named <prefix>_<chr>_<start>_<end>.bin, so it must be globbed. Assert one match: a
        # stale second .bin would make the glob's choice arbitrary.
        ls -1 ~{output_prefix}_*bin > bins.txt
        n=$(wc -l < bins.txt)
        if [ "${n}" -ne 1 ]; then
            echo "ERROR: expected exactly 1 .bin, found ${n}:" >&2
            cat bins.txt >&2
            exit 1
        fi
    >>>

    output {
        File panel_split_bin = glob("~{output_prefix}_*bin")[0]
    }

    runtime {
        cpu: 4
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " SSD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: docker
    }
}

task GLIMPSE2Phase {
    input {
        File cohort_gl_bcf
        File cohort_gl_bcf_idx
        File panel_split_bin
        String output_prefix
        Int phase_threads
        Int phase_kpbwt
        String extra_phase_args
        Int mem_gb
        String docker
    }

    # The full 5.2 Mb shard (433,383 sites, 22 samples, 4 threads) was on iteration 5/15 at
    # 10 min; hence the scatter. Size memory off the windowed variant count: production's
    # 4 + ceil(11.0 * threads * Kpbwt * L / 1e9) gives ~14 GiB at threads=8, Kpbwt=1000,
    # L ~= 105k for an 800 kb window. 32 also covers Kpbwt unset at 2000 states.
    Int disk_gb = 20 + 3 * ceil(size(panel_split_bin, "GiB") + size(cohort_gl_bcf, "GiB"))

    command <<<
        set -euxo pipefail

        /bin/GLIMPSE2_phase \
            --input-gl ~{cohort_gl_bcf} \
            -R ~{panel_split_bin} \
            --thread ~{phase_threads} \
            --Kpbwt ~{phase_kpbwt} \
            ~{extra_phase_args} \
            --output ~{output_prefix}.bcf

        bcftools index -f ~{output_prefix}.bcf
    >>>

    output {
        File phased_bcf = "~{output_prefix}.bcf"
        File phased_bcf_idx = "~{output_prefix}.bcf.csi"
    }

    runtime {
        # cpu >= phase_threads, or GLIMPSE2 oversubscribes the container.
        cpu: phase_threads
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " SSD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: docker
    }
}

task GLIMPSE2Ligate {
    input {
        Array[File] phased_bcfs
        # Unreferenced on purpose: declaring it makes Cromwell localize each index next to its
        # BCF, which GLIMPSE2_ligate opens implicitly.
        Array[File] phased_bcf_idxs
        String output_prefix
        Int ligate_threads
        Int mem_gb
        String docker
    }

    Int disk_gb = 20 + 3 * ceil(size(phased_bcfs, "GiB"))

    command <<<
        set -euxo pipefail

        /bin/GLIMPSE2_ligate \
            --input ~{write_lines(phased_bcfs)} \
            --output ~{output_prefix}.bcf \
            --thread ~{ligate_threads}

        # GLIMPSE2_ligate's own index is corrupt (htslib#1740); regenerate. Same workaround as
        # GLIMPSE2FromPreprocessedPLsJoint.wdl.
        bcftools index -f ~{output_prefix}.bcf
    >>>

    output {
        File ligated_bcf = "~{output_prefix}.bcf"
        File ligated_bcf_idx = "~{output_prefix}.bcf.csi"
    }

    runtime {
        cpu: 4
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " SSD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: docker
    }
}
