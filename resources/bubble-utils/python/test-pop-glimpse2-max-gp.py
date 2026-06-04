import random

def generate_test_vcf():
    with open("resource.vcf", "w") as res_out, open("input.vcf", "w") as in_out, open("expected.vcf", "w") as exp_out:
        
        # --- Write Headers ---
        res_out.write("##fileformat=VCFv4.2\n")
        res_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        
        in_out.write("##fileformat=VCFv4.2\n")
        in_out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        in_out.write('##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype Probabilities">\n')
        in_out.write('##FORMAT=<ID=CL_TRUTH,Number=2,Type=Integer,Description="Truth for Collision Boolean (Hap0,Hap1)">\n')
        in_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join([f"SAMP{i}" for i in range(1, 11)]) + "\n")
        
        exp_out.write("##fileformat=VCFv4.2\n")
        exp_out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        exp_out.write('##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype Probabilities">\n')
        exp_out.write('##FORMAT=<ID=CL_TRUTH,Number=2,Type=Integer,Description="Truth for Collision Boolean (Hap0,Hap1)">\n')
        exp_out.write('##FORMAT=<ID=CL,Number=2,Type=Integer,Description="Collision indicator array (Hap0,Hap1): 0=optimal, 1=lost by GP, 2=lost by parsimony, 3=lost by stable sort">\n')
        exp_out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join([f"SAMP{i}" for i in range(1, 11)]) + "\n")
        
        # --- Define bases for the atomic components ---
        comp_refs = ["A", "C", "G", "T"]
        comp_alts = ["T", "G", "A", "C"]
        
        # --- Generate 25 bubbles with 4 paths each ---
        for b in range(1, 26): 
            bubble_pos = b * 1000
            
            input_sample_data = [[] for _ in range(4)]
            expected_sample_data = [[] for _ in range(4)]
            
            comp_ids = [f"comp_{b}_{i}" for i in range(4)]
            
            bubble_ref = "".join(comp_refs)
            
            # --- Generate Resource VCF Records ---
            for i, comp_id in enumerate(comp_ids):
                orig_var_id = f"var_{b}_{i}"
                res_pos = bubble_pos + i 
                res_out.write(f"chr1\t{res_pos}\t{orig_var_id}\t{comp_refs[i]}\t{comp_alts[i]}\t.\t.\tID={comp_id}\n")

            # --- Define DAG Bubble Paths ---
            path_indices = [
                [0, 1],       
                [1, 2],       
                [2, 3],       
                [0, 2, 3]     
            ]
            
            paths = [[comp_ids[i] for i in indices] for indices in path_indices]
            path_info_ids = [":".join(p) for p in paths]
            
            path_indices_for_c = [
                [0, 3],       
                [0, 1],       
                [1, 2, 3],    
                [2, 3]        
            ]

            # --- Generate Samples, Collisions, and Expected Truth ---
            for s in range(10):
                hap0_calls = [0]*4
                hap1_calls = [0]*4
                
                hap0_calls[random.randint(0, 3)] = 1
                hap1_calls[random.randint(0, 3)] = 1
                
                if random.random() < 0.3: hap0_calls[random.randint(0, 3)] = 1
                if random.random() < 0.3: hap1_calls[random.randint(0, 3)] = 1
                    
                gps = []
                gp_strs = []
                for a in range(4):
                    # Calculate index for GP array: 0 for 0|0, 1 for 0|1 or 1|0, 2 for 1|1
                    gt_idx = hap0_calls[a] + hap1_calls[a]
                    
                    # Generate base probabilities favoring the correct index
                    probs = [random.uniform(0.01, 0.05) for _ in range(3)]
                    probs[gt_idx] = random.uniform(0.80, 0.98)
                    
                    # Normalize to 1.0
                    total = sum(probs)
                    probs = [p / total for p in probs]
                    
                    # Round and cleanly absorb floating point differences into the dominant trait
                    rounded = [round(p, 2) for p in probs]
                    diff = round(1.0 - sum(rounded), 2)
                    rounded[gt_idx] = round(rounded[gt_idx] + diff, 2)
                    
                    strs = [f"{x:.2f}" for x in rounded]
                    gp_strs.append(",".join(strs))
                    gps.append(rounded)
                
                hap0_ones = [(max(gps[a]), len(path_indices[a]), a) for a in range(4) if hap0_calls[a] == 1]
                hap1_ones = [(max(gps[a]), len(path_indices[a]), a) for a in range(4) if hap1_calls[a] == 1]
                
                hap0_ones.sort(key=lambda x: (x[0], -x[1]), reverse=True)
                hap1_ones.sort(key=lambda x: (x[0], -x[1]), reverse=True)
                
                winner_h0_tuple = hap0_ones[0] if hap0_ones else None
                winner_h1_tuple = hap1_ones[0] if hap1_ones else None
                
                # 1. Store input data (bubble level)
                for a in range(4):
                    cl0 = '0'
                    if hap0_calls[a] == 1 and winner_h0_tuple and a != winner_h0_tuple[2]:
                        if winner_h0_tuple[0] > max(gps[a]): cl0 = '1'
                        elif winner_h0_tuple[1] < len(path_indices[a]): cl0 = '2'
                        else: cl0 = '3'
                        
                    cl1 = '0'
                    if hap1_calls[a] == 1 and winner_h1_tuple and a != winner_h1_tuple[2]:
                        if winner_h1_tuple[0] > max(gps[a]): cl1 = '1'
                        elif winner_h1_tuple[1] < len(path_indices[a]): cl1 = '2'
                        else: cl1 = '3'
                        
                    cl_truth = f"{cl0},{cl1}"
                    
                    gt = f"{hap0_calls[a]}|{hap1_calls[a]}"
                    gp_str = gp_strs[a]
                    input_sample_data[a].append(f"{gt}:{gp_str}:{cl_truth}")
                
                # 2. Calculate expected output data (component level)
                for c in range(4):
                    paths_for_c = path_indices_for_c[c]
                    
                    h0_winner = winner_h0_tuple[2] in paths_for_c if winner_h0_tuple else False
                    h0_losers = []
                    for p in paths_for_c:
                        if hap0_calls[p] == 1 and (not winner_h0_tuple or p != winner_h0_tuple[2]):
                            if winner_h0_tuple[0] > max(gps[p]): h0_losers.append(1)
                            elif winner_h0_tuple[1] < len(path_indices[p]): h0_losers.append(2)
                            else: h0_losers.append(3)
                            
                    h0_gt = '1' if (h0_winner or h0_losers) else '0'
                    h0_cl = '0' if h0_winner or not h0_losers else str(max(h0_losers))
                    
                    h1_winner = winner_h1_tuple[2] in paths_for_c if winner_h1_tuple else False
                    h1_losers = []
                    for p in paths_for_c:
                        if hap1_calls[p] == 1 and (not winner_h1_tuple or p != winner_h1_tuple[2]):
                            if winner_h1_tuple[0] > max(gps[p]): h1_losers.append(1)
                            elif winner_h1_tuple[1] < len(path_indices[p]): h1_losers.append(2)
                            else: h1_losers.append(3)
                            
                    h1_gt = '1' if (h1_winner or h1_losers) else '0'
                    h1_cl = '0' if h1_winner or not h1_losers else str(max(h1_losers))
                    
                    expected_sample_data[c].append(f"{h0_gt}|{h1_gt}:{h0_cl},{h1_cl}")

            # --- Write Input VCF Records ---
            for a in range(4):
                info = f"ID={path_info_ids[a]}"
                fmt = "GT:GP:CL_TRUTH"
                
                path_alt = ""
                for i in range(4):
                    if i in path_indices[a]:
                        path_alt += comp_alts[i]
                    else:
                        path_alt += comp_refs[i]
                        
                row = ["chr1", str(bubble_pos), ".", bubble_ref, path_alt, ".", ".", info, fmt] + input_sample_data[a]
                in_out.write("\t".join(row) + "\n")
                
            # --- Write Expected VCF Records ---
            for c in range(4):
                res_pos = bubble_pos + c 
                orig_var_id = f"var_{b}_{c}"
                comp_id = f"comp_{b}_{c}"
                info = f"ID={comp_id}"
                fmt = "GT:CL"
                row = ["chr1", str(res_pos), orig_var_id, comp_refs[c], comp_alts[c], ".", ".", info, fmt] + expected_sample_data[c]
                exp_out.write("\t".join(row) + "\n")

if __name__ == "__main__":
    generate_test_vcf()
    print("Successfully generated 'resource.vcf', 'input.vcf', and 'expected.vcf'.")
    print("Test pipeline command:")
    print("cat input.vcf | python convert-to-biallelic.py resource.vcf > output.vcf")
    print("diff expected.vcf output.vcf")
