import random

def generate_test_vcf():
    print("##fileformat=VCFv4.2")
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print('##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype Probabilities">')
    print('##FORMAT=<ID=CID_TRUTH,Number=1,Type=String,Description="Truth for Collision ID">')
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join([f"SAMP{i}" for i in range(1, 11)]))
    
    # 25 bubbles with 4 alleles each = 100 records
    for b in range(1, 26): 
        pos = b * 100
        alleles = [f"node{b}_{a}" for a in range(1, 5)]
        
        sample_data_list = [[] for _ in range(4)]
        for s in range(10):
            hap0_calls = [0]*4
            hap1_calls = [0]*4
            
            # assign true calls
            true_h0 = random.randint(0, 3)
            true_h1 = random.randint(0, 3)
            hap0_calls[true_h0] = 1
            hap1_calls[true_h1] = 1
            
            # arbitrarily inject collisions
            if random.random() < 0.3:
                hap0_calls[random.randint(0, 3)] = 1
            if random.random() < 0.3:
                hap1_calls[random.randint(0, 3)] = 1
                
            # generate GPs
            gps = []
            for a in range(4):
                gp0 = random.random() * 0.1
                gp1 = random.random() * 0.9 + 0.1
                gp2 = random.random() * 0.1
                gps.append([gp0, gp1, gp2])
            
            # Determine truth CID
            hap0_ones = [(max(gps[a]), a) for a in range(4) if hap0_calls[a] == 1]
            hap1_ones = [(max(gps[a]), a) for a in range(4) if hap1_calls[a] == 1]
            
            cid0 = alleles[max(hap0_ones)[1]] if len(hap0_ones) > 1 else '.'
            cid1 = alleles[max(hap1_ones)[1]] if len(hap1_ones) > 1 else '.'
            cid_truth = f"{cid0}|{cid1}"
            
            for a in range(4):
                gt = f"{hap0_calls[a]}|{hap1_calls[a]}"
                gp_str = ",".join([f"{x:.2f}" for x in gps[a]])
                sample_data_list[a].append(f"{gt}:{gp_str}:{cid_truth}")
                
        for a in range(4):
            info = f"ID={alleles[a]}"
            fmt = "GT:GP:CID_TRUTH"
            row = ["chr1", str(pos), alleles[a], "A", "T", ".", ".", info, fmt] + sample_data_list[a]
            print("\t".join(row))

if __name__ == "__main__":
    generate_test_vcf()
