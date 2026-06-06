import random

def generate_sparse_lpl(true_a, true_b):
    """Generates a DRAGEN-compliant sparse LA and LPL array."""
    alleles = sorted(list(set([0, true_a, true_b])))
    if len(alleles) == 1:
        alleles.append(1)
        
    la_str = ",".join(map(str, alleles))
    num_local = len(alleles)
    
    lpl_len = (num_local * (num_local + 1)) // 2
    lpls = [50] * lpl_len 
    
    loc_a = alleles.index(true_a)
    loc_b = alleles.index(true_b)
    
    idx = (max(loc_a, loc_b) * (max(loc_a, loc_b) + 1)) // 2 + min(loc_a, loc_b)
    lpls[idx] = 0  
    
    return la_str, ",".join(map(str, lpls))

def main():
    random.seed(42)
    
    num_samples = 10
    chrom = "chr1"
    sample_names = [f"Sample_{i}" for i in range(1, num_samples + 1)]
    
    # ---------------------------------------------------------
    # 1. Define the 6 Bubbles
    # ---------------------------------------------------------
    b1_snvs = [f"snv_{i}" for i in range(10)]
    b1_alleles = {0: set(), 1: set(b1_snvs)}
    
    b2_snvs = [f"snv_{i}" for i in range(10, 100)]
    b2_alleles = {0: set()}
    for i in range(1, 101):
        b2_alleles[i] = set(random.sample(b2_snvs, random.randint(5, 50)))
        
    b3_alleles = {0: set(), 1: {"snv_100"}}
    b4_alleles = {0: set(), 1: {"snv_101"}}
    
    # Bubble VCF Canonical ALTs
    b5_alts_bubble = ["A" + ("C" * i) for i in range(1, 11)] # AC, ACC, ACCC...
    b6_alts_bubble = ["C", "G", "T", "AA", "AC", "AG", "AT", "CC", "CG", "CT"]

    # Input VCF Shifted/Padded ALTs (Mathematically equivalent)
    b5_alts_input  = ["TA" + ("C" * i) for i in range(1, 11)] # TAC, TACC, TACCC...
    b6_alts_input  = [a + "T" for a in b6_alts_bubble]        # CT, GT, TT, AAT...

    # ---------------------------------------------------------
    # 2. Draw Haplotypes
    # ---------------------------------------------------------
    sample_snv_states = {snv: [0] * num_samples for snv in [f"snv_{i}" for i in range(102)]}
    sample_gts = []

    for s in range(num_samples):
        gts = {
            "b1": random.choice([(0,0), (0,1), (1,1)]),
            "b2": random.choice([(0, 0), (0, random.randint(1, 100)), (random.randint(1, 100), random.randint(1, 100))]),
            "b3": random.choice([(0,0), (0,1), (1,1)]),
            "b4": random.choice([(0,0), (0,1), (1,1)]),
            "b5": random.choice([(0, 0), (0, random.randint(1, 10)), (random.randint(1, 10), random.randint(1, 10))]),
            "b6": random.choice([(0, 0), (0, random.randint(1, 10)), (random.randint(1, 10), random.randint(1, 10))]),
        }
        for k in gts: gts[k] = tuple(sorted(gts[k]))
        sample_gts.append(gts)
        
        for allele in gts["b1"]:
            for snv in b1_alleles[allele]: sample_snv_states[snv][s] += 1
        for allele in gts["b2"]:
            for snv in b2_alleles[allele]: sample_snv_states[snv][s] += 1
        for allele in gts["b3"]:
            for snv in b3_alleles[allele]: sample_snv_states[snv][s] += 1
        for allele in gts["b4"]:
            for snv in b4_alleles[allele]: sample_snv_states[snv][s] += 1

    # ---------------------------------------------------------
    # 3. Generate Bubble VCF (Canonical Representations)
    # ---------------------------------------------------------
    with open("test_bubble.vcf", "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID={chrom},length=10000>\n")
        f.write('##INFO=<ID=ID,Number=A,Type=String,Description="Constituent SNV IDs">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
        # B1 & B2 (Complex, will correctly be skipped by the Remapper)
        b1_id_str = ":".join(sorted(list(b1_alleles[1]), key=lambda x: int(x.split('_')[1])))
        f.write(f"{chrom}\t100\tbub_1\tA\t<ALT1>\t50\tPASS\tID={b1_id_str}\n")
        
        b2_alts_str = ",".join(f"<ALT{i}>" for i in range(1, 101))
        b2_id_str = ",".join([":".join(sorted(list(b2_alleles[i]), key=lambda x: int(x.split('_')[1]))) for i in range(1, 101)])
        f.write(f"{chrom}\t200\tbub_2\tA\t{b2_alts_str}\t50\tPASS\tID={b2_id_str}\n")
        
        # B3 & B4 (Canonical)
        f.write(f"{chrom}\t1100\tbub_3\tA\tATGC\t50\tPASS\tID=snv_100\n")
        f.write(f"{chrom}\t1110\tbub_4\tATGC\tA\t50\tPASS\tID=snv_101\n")
        
        # B5 & B6 (Canonical)
        b5_alts_str = ",".join(b5_alts_bubble)
        b5_id_str = ",".join([f"snv_b5_{i}" for i in range(1, 11)])
        f.write(f"{chrom}\t4000\tbub_5\tA\t{b5_alts_str}\t50\tPASS\tID={b5_id_str}\n")
        
        b6_alts_str = ",".join(b6_alts_bubble)
        b6_id_str = ",".join([f"snv_b6_{i}" for i in range(1, 11)])
        f.write(f"{chrom}\t5000\tbub_6\tA\t{b6_alts_str}\t50\tPASS\tID={b6_id_str}\n")

    # ---------------------------------------------------------
    # 4. Generate Input VCF (Shifted/Padded Representations)
    # ---------------------------------------------------------
    with open("test_input.vcf", "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write(f"##contig=<ID={chrom},length=10000>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('##FORMAT=<ID=LA,Number=.,Type=Integer,Description="Local Alleles">\n')
        f.write('##FORMAT=<ID=LPL,Number=G,Type=Integer,Description="Sparse Likelihoods">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names) + "\n")
        
        for i in range(102):
            # Apply padding/shifting logic to B3 and B4 only
            if i == 100:   # B3 (Insertion prefix-padded)
                pos, ref, alt = 1099, "CA", "CATGC"
            elif i == 101: # B4 (Deletion prefix-padded)
                pos, ref, alt = 1109, "GATGC", "GA"
            else:          # Normal SNVs
                pos, ref, alt = 100 + (i * 10), "A", "T"
                
            snv_id = f"snv_{i}"
            row = [chrom, str(pos), ".", ref, alt, "50", "PASS", ".", "GT:LA:LPL"]
            
            for s in range(num_samples):
                true_state = sample_snv_states[snv_id][s]
                if true_state == 0: row.append("0/0:0,1:0,20,40")
                elif true_state == 1: row.append("0/1:0,1:20,0,20")
                else: row.append("1/1:0,1:40,20,0")
            f.write("\t".join(row) + "\n")
            
        # Apply padding/shifting logic to B5 and B6 multiallelic rows
        shifted_bubbles = [
            ("b5", 3999, "TA", b5_alts_input), 
            ("b6", 5000, "AT", b6_alts_input)
        ]
        
        for b_name, b_pos, b_ref, b_alts in shifted_bubbles:
            row = [chrom, str(b_pos), ".", b_ref, ",".join(b_alts), "50", "PASS", ".", "GT:LA:LPL"]
            for s in range(num_samples):
                a, b = sample_gts[s][b_name]
                la_str, lpl_str = generate_sparse_lpl(a, b)
                row.append(f"{a}/{b}:{la_str}:{lpl_str}")
            f.write("\t".join(row) + "\n")
            
    print("Test data generated successfully! The Input VCF contains shifted/padded representations.")
    print("-" * 75)
    print("GROUND TRUTH (Simulated Genotypes):")
    print(f"{'Sample':<10} | {'B1 (Cpx)':<8} | {'B2 (Cpx)':<8} | {'B3 (Bi)':<8} | {'B4 (Bi)':<8} | {'B5 (Multi)':<10} | {'B6 (Multi)':<10}")
    print("-" * 75)
    for s, name in enumerate(sample_names):
        g = sample_gts[s]
        print(f"{name:<10} | {g['b1'][0]}/{g['b1'][1]:<8} | {g['b2'][0]}/{g['b2'][1]:<8} | "
              f"{g['b3'][0]}/{g['b3'][1]:<8} | {g['b4'][0]}/{g['b4'][1]:<8} | "
              f"{g['b5'][0]}/{g['b5'][1]:<10} | {g['b6'][0]}/{g['b6'][1]:<10}")
    print("-" * 75)

if __name__ == "__main__":
    main()
