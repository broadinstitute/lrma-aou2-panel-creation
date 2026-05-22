task CreateShards {
    input {
        String? region
        String output_prefix
        String entity_name
        
        File short_vcf
        File short_vcf_idx
        File sv_vcf
        File sv_vcf_idx

        Int min_variants_per_shard
        Int max_variants_per_shard
        Int min_boundary_dist_bp
        Int min_sv_len

        RuntimeAttr? runtime_attr_override
    }

    Int disk_gb = 10 + 2 * ceil(size([short_vcf, sv_vcf], "GiB"))

    command <<<
        set -euxo pipefail

        python3 -m pip install tqdm pysam

        python3 - \
            ~{"--region '" + region + "'"} \
            --vcf '~{short_vcf}##idx##~{short_vcf_idx}' \
            --sv_vcf '~{sv_vcf}##idx##~{sv_vcf_idx}' \
            --min_vars ~{min_variants_per_shard} \
            --max_vars ~{max_variants_per_shard} \
            --min_sv_dist ~{min_boundary_dist_bp} \
            --min_sv_len ~{min_sv_len} \
            --entity_name ~{entity_name} \
            --output_prefix ~{output_prefix} \
            <<-'EOF'
        import argparse
        import pysam
        import subprocess
        import bisect
        import sys
        from tqdm import tqdm

        def get_sv_intervals(sv_vcf_path, chrom, start, end, min_sv_len):
            svs = []
            try:
                vcf = pysam.VariantFile(sv_vcf_path, drop_samples=True)
                for rec in tqdm(vcf.fetch(chrom, start, end), desc=f"Reading SVs for {chrom}", leave=False, unit=" SVs"):
                    svlen = 0
                    if "SVLEN" in rec.info:
                        val = rec.info["SVLEN"]
                        svlen = abs(val[0]) if isinstance(val, (tuple, list)) else abs(val)
                    else:
                        svlen = abs(rec.stop - rec.pos)

                    if svlen >= min_sv_len:
                        svs.append((rec.pos, rec.stop))
                vcf.close()
            except Exception as e:
                print(f"Warning reading SVs for {chrom}: {e}", file=sys.stderr)
            return svs

        def get_all_short_variants(vcf_path, chrom, r_start, r_end):
            cmd = f"bcftools query -f '%POS\n' -r {chrom}:{r_start}-{r_end} '{vcf_path}'"

            # FIX 1: Send stderr directly to sys.stderr to prevent OS pipe buffer deadlocks
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=sys.stderr, text=True)

            pos_list = []
            for line in tqdm(proc.stdout, desc=f"Reading short vars for {chrom}", leave=False, unit=" vars"):
                if line.strip():
                    pos_list.append(int(line.strip()))

            proc.wait()

            # FIX 2: Gracefully handle missing sex chromosomes instead of hard crashing
            if proc.returncode != 0:
                print(f"Warning: bcftools exited with code {proc.returncode} for {chrom}. Proceeding with 0 short variants.", file=sys.stderr)
                return []

            return pos_list

        def calc_dist_to_svs(pos, sorted_svs, sv_starts, max_sv_len):
            if not sorted_svs:
                return float('inf')

            idx = bisect.bisect_left(sv_starts, pos)
            min_d = float('inf')

            for i in range(idx, len(sorted_svs)):
                s, e = sorted_svs[i]
                if s <= pos <= e:
                    return 0
                d = min(abs(pos - s), abs(pos - e))
                if d < min_d: min_d = d
                if s - pos >= min_d: 
                    break

            for i in range(idx - 1, -1, -1):
                s, e = sorted_svs[i]
                if s <= pos <= e:
                    return 0
                d = min(abs(pos - s), abs(pos - e))
                if d < min_d: min_d = d
                if (pos - s) >= min_d + max_sv_len:
                    break

            return min_d

        def process_region(chrom, r_start, r_end, args, f_tsv, f_txt, start_shard_idx):
            raw_svs = get_sv_intervals(args.sv_vcf, chrom, r_start, r_end, args.min_sv_len)
            short_vars = get_all_short_variants(args.vcf, chrom, r_start, r_end)

            sorted_svs = sorted(raw_svs, key=lambda x: x[0])
            sv_starts = [s[0] for s in sorted_svs]
            max_sv_len = max([e - s for s, e in sorted_svs]) if sorted_svs else 0

            all_vars = sorted(sv_starts + short_vars)

            shards = []
            current_start = r_start
            ideal_vars = (args.min_vars + args.max_vars) / 2
            shard_idx = start_shard_idx

            with tqdm(total=r_end - r_start + 1, desc=f"Sharding {chrom}", unit=" bp") as pbar:
                while current_start <= r_end:
                    remaining_bp = r_end - current_start + 1

                    idx = bisect.bisect_left(all_vars, current_start)
                    remaining_total_vars = len(all_vars) - idx

                    if remaining_total_vars <= args.max_vars:
                        d_start = calc_dist_to_svs(current_start, sorted_svs, sv_starts, max_sv_len)
                        d_end = calc_dist_to_svs(r_end, sorted_svs, sv_starts, max_sv_len)

                        shard_sv_count = bisect.bisect_right(sv_starts, r_end) - bisect.bisect_left(sv_starts, current_start)
                        shard_short_count = bisect.bisect_right(short_vars, r_end) - bisect.bisect_left(short_vars, current_start)

                        shards.append((current_start, r_end, remaining_bp, shard_sv_count, shard_short_count, d_start, d_end))
                        pbar.update(remaining_bp)
                        break

                    expected_shards_left = max(1, round(remaining_total_vars / ideal_vars))
                    target_total_vars = remaining_total_vars / expected_shards_left

                    k_min = idx + args.min_vars
                    k_max = min(idx + args.max_vars, len(all_vars))

                    if 0 < (len(all_vars) - k_max) < args.min_vars:
                        k_max = len(all_vars) - args.min_vars
                        if k_min > k_max: k_max = k_min

                    search_start = all_vars[k_min - 1]
                    search_end = all_vars[k_max] - 1 if k_max < len(all_vars) else r_end

                    search_start = max(current_start, search_start)
                    search_end = min(r_end, max(search_start, search_end))

                    def score_candidate(cand):
                        c_var_count = bisect.bisect_right(all_vars, cand) - idx
                        var_diff = abs(c_var_count - target_total_vars)
                        dist = calc_dist_to_svs(cand, sorted_svs, sv_starts, max_sv_len)
                        return (-dist, var_diff, -cand)

                    valid_candidates = []
                    expanded_search_start = search_start
                    expanded_search_end = search_end
                    expansion_step = 50000

                    while not valid_candidates:
                        candidates = set([expanded_search_start, expanded_search_end])

                        s_idx = bisect.bisect_left(sv_starts, expanded_search_start)
                        e_idx = bisect.bisect_right(sv_starts, expanded_search_end)
                        
                        for i in range(max(0, s_idx - 1), min(len(sorted_svs), e_idx + 1)):
                            s1, e1 = sorted_svs[i]
                            
                            safe_left = s1 - args.min_sv_dist
                            safe_right = e1 + args.min_sv_dist
                            
                            if expanded_search_start <= safe_left <= expanded_search_end:
                                candidates.add(safe_left)
                            if expanded_search_start <= safe_right <= expanded_search_end:
                                candidates.add(safe_right)
                                
                            if i < len(sorted_svs) - 1:
                                s2, e2 = sorted_svs[i+1]
                                mid = (e1 + s2) // 2
                                if expanded_search_start <= mid <= expanded_search_end:
                                    candidates.add(mid)

                        for p in range(expanded_search_start, expanded_search_end + 1, 1000):
                            candidates.add(p)

                        valid_candidates = [c for c in candidates if calc_dist_to_svs(c, sorted_svs, sv_starts, max_sv_len) >= args.min_sv_dist]

                        if valid_candidates:
                            break

                        prev_start, prev_end = expanded_search_start, expanded_search_end
                        expanded_search_start = max(current_start + 1, expanded_search_start - expansion_step)
                        expanded_search_end = min(r_end, expanded_search_end + expansion_step)

                        if expanded_search_start == prev_start and expanded_search_end == prev_end:
                            break

                    if not valid_candidates:
                        valid_candidates = list(candidates)

                    # FIX 3: Prevent the 1-bp crawling bug by strictly enforcing forward movement
                    valid_candidates = [c for c in valid_candidates if c > current_start]
                    
                    # Absolute failsafe if no candidates > current_start were generated
                    if not valid_candidates:
                        valid_candidates = [min(r_end, current_start + 500000)]

                    valid_candidates = sorted(valid_candidates, key=score_candidate)

                    best_boundary = valid_candidates[0]
                    best_dist = calc_dist_to_svs(best_boundary, sorted_svs, sv_starts, max_sv_len)

                    shard_sv_count = bisect.bisect_right(sv_starts, best_boundary) - bisect.bisect_left(sv_starts, current_start)
                    shard_short_count = bisect.bisect_right(short_vars, best_boundary) - bisect.bisect_left(short_vars, current_start)

                    d_start = calc_dist_to_svs(current_start, sorted_svs, sv_starts, max_sv_len)
                    size_bp = best_boundary - current_start + 1

                    shards.append((current_start, best_boundary, size_bp, shard_sv_count, shard_short_count, d_start, best_dist))

                    current_start = best_boundary + 1
                    pbar.update(size_bp)

            for s in shards:
                c_start, c_end, s_bp, n_svs, n_short, d_s, d_e = s
                reg_str = f"{chrom}:{c_start}-{c_end}"

                min_d = min(d_s, d_e)
                if min_d == float('inf'):
                    min_d = -1

                shard_id = f"{args.output_prefix}.shard-{shard_idx:04d}.{chrom}-{c_start}-{c_end}"
                shard_idx += 1

                f_tsv.write(f"{shard_id}\t{reg_str}\t{s_bp}\t{n_svs}\t{n_short}\t{min_d}\n")
                f_txt.write(reg_str + "\n")
            
            return shard_idx

        def main():
            parser = argparse.ArgumentParser()
            parser.add_argument('--region', type=str, required=False, default=None)
            parser.add_argument('--vcf', type=str, required=True)
            parser.add_argument('--sv_vcf', type=str, required=True)
            parser.add_argument('--min_vars', type=int, required=True)
            parser.add_argument('--max_vars', type=int, required=True)
            parser.add_argument('--min_sv_dist', type=int, required=True)
            parser.add_argument('--min_sv_len', type=int, required=True)
            parser.add_argument('--entity_name', type=str, required=True)
            parser.add_argument('--output_prefix', type=str, required=True)
            args = parser.parse_args()

            regions_to_process = []

            if args.region:
                if ":" in args.region and "-" in args.region:
                    chrom, span = args.region.split(":")
                    r_start, r_end = map(int, span.split("-"))
                    regions_to_process.append((chrom, r_start, r_end))
                else:
                    chrom = args.region
                    with pysam.VariantFile(args.sv_vcf, drop_samples=True) as vcf:
                        if chrom in vcf.header.contigs:
                            r_end = vcf.header.contigs[chrom].length
                            if r_end:
                                regions_to_process.append((chrom, 1, r_end))
                                print(f"Detected whole chromosome '{chrom}'. Extracted length: {r_end} bp")
                            else:
                                raise ValueError(f"Chromosome '{chrom}' missing length in VCF header.")
                        else:
                            raise ValueError(f"Chromosome '{chrom}' not found in VCF header.")
            else:
                print("No region provided. Extracting all contigs from SV VCF header...")
                with pysam.VariantFile(args.sv_vcf, drop_samples=True) as vcf:
                    for chrom, contig in vcf.header.contigs.items():
                        if contig.length:
                            regions_to_process.append((chrom, 1, contig.length))
                        else:
                            print(f"Skipping '{chrom}' as it has no defined length in the header.")

            tsv_file = args.output_prefix + ".shards.tsv"
            txt_file = args.output_prefix + ".shards.txt"

            start_shard_idx = 0
            with open(tsv_file, 'w') as f_tsv, open(txt_file, 'w') as f_txt:
                f_tsv.write(f"entity:{args.entity_name}_id\tregion\tsize_bp\tnumber_of_svs\tnumber_of_short_variants\tmin_sv_boundary_dist\n")

                for chrom, r_start, r_end in regions_to_process:
                    print(f"\nProcessing {chrom}:{r_start}-{r_end}")
                    start_shard_idx = process_region(chrom, r_start, r_end, args, f_tsv, f_txt, start_shard_idx)
                    
        if __name__ == "__main__":
            main()
        EOF
    >>>

    output {
        Array[String] shard_regions = read_lines("~{output_prefix}.shards.txt")
        File shard_regions_terra_tsv = "~{output_prefix}.shards.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_gb,
        boot_disk_gb:       10,
        use_ssd:            false,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsde-methods/slee/lrma-aou2-panel-creation-python:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + if select_first([runtime_attr.use_ssd, default_attr.use_ssd]) then " SSD" else " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
