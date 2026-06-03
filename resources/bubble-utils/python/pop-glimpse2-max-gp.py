#!/usr/bin/env python

# a version of https://github.com/eblerjana/pangenie/blob/master/pipelines/run-from-callset/scripts/convert-to-biallelic.py with:
# 1) / -> | in GTs
# 2) resolution of bubble-level collisions by max GP before popping (the INFO/ID of the kept bubble allele is stored as FORMAT/CID)

import sys
import argparse
from collections import defaultdict
import gzip

parser = argparse.ArgumentParser(prog='convert-to-biallelic.py', description='cat <multiallelic VCF> | python convert-to-biallelic.py <biallelic VCF>')
parser.add_argument('vcf', metavar='VCF', help='original VCF containing REF/ALT of each Variant ID.')
args = parser.parse_args()

# chromosome ->  ID -> [start, REF, ALT] per chromosome
chrom_to_variants = defaultdict(lambda: defaultdict(list))

# read the biallelic VCF containing REF/ALT for all variant IDs and store them
for line in gzip.open(args.vcf, 'rt'):
    if line.startswith('#'):
        continue
    fields = line.split()
    info_field = { i.split('=')[0] : i.split('=')[1] for i in fields[7].split(';') if "=" in i}
    if 'ID' not in info_field: continue
    ids = info_field['ID'].split(',')
    chrom_to_variants[fields[0]][ids[0]] = [fields[1], fields[3], fields[4]]

def process_group(group_lines):
    if not group_lines: return
    
    parsed_lines = [line.strip().split('\t') for line in group_lines]
    if len(parsed_lines[0]) <= 9:
        for fields in parsed_lines: print('\t'.join(fields))
        return
        
    num_samples = len(parsed_lines[0]) - 9
    line_ids = []
    
    for fields in parsed_lines:
        info_field = { i.split('=')[0] : i.split('=')[1] for i in fields[7].split(';') if "=" in i}
        line_ids.append(info_field.get('ID', ''))
        
    # Inject CID to FORMAT if not present and pad sample fields
    for fields in parsed_lines:
        fmt = fields[8].split(':')
        if 'CID' not in fmt:
            fields[8] = fields[8] + ':CID'
            for s in range(num_samples):
                fields[9+s] += ':.'
                
    # Resolve collisions
    for s in range(num_samples):
        col = 9 + s
        hap0_ones = []
        hap1_ones = []
        
        for i, fields in enumerate(parsed_lines):
            fmt = fields[8].split(':')
            if 'GT' not in fmt: continue
            gt_idx = fmt.index('GT')
            sample_data = fields[col].split(':')
            gt = sample_data[gt_idx]
            
            if '|' not in gt: continue
            alleles = gt.split('|')
            
            gp_val = -1.0
            if 'GP' in fmt:
                gp_idx = fmt.index('GP')
                if gp_idx < len(sample_data) and sample_data[gp_idx] != '.':
                    gp_val = max([float(x) for x in sample_data[gp_idx].split(',')])
                    
            if alleles[0] == '1': hap0_ones.append((gp_val, i))
            if alleles[1] == '1': hap1_ones.append((gp_val, i))
            
        cid_hap = ['.', '.']
        
        # Resolve hap0 collision
        if len(hap0_ones) > 1:
            hap0_ones.sort(key=lambda x: x[0], reverse=True)
            kept_idx = hap0_ones[0][1]
            cid_hap[0] = line_ids[kept_idx]
            for _, lost_idx in hap0_ones[1:]:
                fields = parsed_lines[lost_idx]
                fmt = fields[8].split(':')
                gt_idx = fmt.index('GT')
                sample_data = fields[col].split(':')
                alleles = sample_data[gt_idx].split('|')
                alleles[0] = '0'
                sample_data[gt_idx] = '|'.join(alleles)
                parsed_lines[lost_idx][col] = ':'.join(sample_data)
                
        # Resolve hap1 collision
        if len(hap1_ones) > 1:
            hap1_ones.sort(key=lambda x: x[0], reverse=True)
            kept_idx = hap1_ones[0][1]
            cid_hap[1] = line_ids[kept_idx]
            for _, lost_idx in hap1_ones[1:]:
                fields = parsed_lines[lost_idx]
                fmt = fields[8].split(':')
                gt_idx = fmt.index('GT')
                sample_data = fields[col].split(':')
                alleles = sample_data[gt_idx].split('|')
                alleles[1] = '0'
                sample_data[gt_idx] = '|'.join(alleles)
                parsed_lines[lost_idx][col] = ':'.join(sample_data)
                
        # Update CID in sample format block
        sample_cid_str = '|'.join(cid_hap)
        for i, fields in enumerate(parsed_lines):
            fmt = fields[8].split(':')
            cid_idx = fmt.index('CID')
            sample_data = fields[col].split(':')
            sample_data[cid_idx] = sample_cid_str
            parsed_lines[i][col] = ':'.join(sample_data)

    # Pop bubbles to constituent variants
    for fields in parsed_lines:
        info_field = { i.split('=')[0] : i.split('=')[1] for i in fields[7].split(';') if "=" in i}
        if 'ID' not in info_field:
            print('\t'.join(fields))
            continue
            
        allele_to_ids = [''] + info_field['ID'].split(',')
        info_ids = info_field['ID'].split(',')
        if (len(info_ids) == 1) and any([x not in chrom_to_variants[fields[0]] for x in info_ids[0].split(':')]):
            print('\t'.join(fields))
            continue
            
        ids = set([])
        for i in info_field['ID'].split(','):
            for j in i.split(':'):
                if j in chrom_to_variants[fields[0]]:
                    ids.add((j, int(chrom_to_variants[fields[0]][j][0])))
        ids = list(ids)
        ids.sort(key=lambda x : x[1])
        
        for (var_id, coord) in ids:
            vcf_line = list(fields[:9])
            vcf_line[1] = str(coord)
            vcf_line[2] = var_id
            vcf_line[3] = chrom_to_variants[fields[0]][var_id][1]
            vcf_line[4] = chrom_to_variants[fields[0]][var_id][2]
            vcf_line[7] = 'ID=' + var_id
            
            for k,v in info_field.items():
                if k == 'ID': continue
                if k in ['MA', 'UK']:
                    vcf_line[7] = vcf_line[7] + ';' + k + '=' + v
                    
            format_field = fields[8].split(':')
            index_of_gt = format_field.index('GT')
            
            for sample_field in fields[9:]:
                genotype = sample_field.split(':')
                biallelic_genotype = []
                for allele in genotype[index_of_gt].split('|'):
                    if allele == '.':
                        biallelic_genotype.append('.')
                    else:
                        if var_id in allele_to_ids[int(allele)].split(':'):
                            biallelic_genotype.append('1')
                        else:
                            biallelic_genotype.append('0')
                genotype[index_of_gt] = '|'.join(biallelic_genotype)
                vcf_line.append(':'.join(genotype))
                
            print('\t'.join(vcf_line))

cid_header_added = False
current_pos = None
group = []

for line in sys.stdin:
    if line.startswith('#'):
        if any([i in line for i in ['INFO=<ID=AF', 'INFO=<ID=AK', 'FORMAT=<ID=GL', 'FORMAT=<ID=KC']]):
            continue
        if line.startswith('#CHROM') and not cid_header_added:
            print('##FORMAT=<ID=CID,Number=1,Type=String,Description="Collision ID: INFO/ID of the bubble allele kept after resolving collisions">')
            cid_header_added = True
        elif line.startswith('##FORMAT=') and not cid_header_added:
            print('##FORMAT=<ID=CID,Number=1,Type=String,Description="Collision ID: INFO/ID of the bubble allele kept after resolving collisions">')
            cid_header_added = True
        print(line[:-1])
        continue
        
    fields = line.split('\t')
    if len(fields) < 2: continue
    pos = fields[1]
    
    if current_pos is None:
        current_pos = pos
        
    if pos != current_pos:
        process_group(group)
        group = [line]
        current_pos = pos
    else:
        group.append(line)
        
if group:
    process_group(group)
