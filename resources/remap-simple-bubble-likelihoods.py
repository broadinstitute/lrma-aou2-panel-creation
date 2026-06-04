import sys
import gzip
import argparse
import math
from typing import Tuple, Optional, Dict, List
from tqdm import tqdm

def open_vcf(path: str):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')

def extract_info_tag(info_str: str, tag: str) -> Optional[str]:
    if info_str == '.': return None
    for part in info_str.split(';'):
        if part.startswith(tag + '='):
            return part.split('=', 1)[1]
        elif part == tag:
            return "true"
    return None

def get_minimal_representation(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    if alt.startswith('<'): return pos, ref, alt
    r, a = ref, alt
    p = pos
    while len(r) > 0 and len(a) > 0 and r[-1] == a[-1]:
        r, a = r[:-1], a[:-1]
    while len(r) > 0 and len(a) > 0 and r[0] == a[0]:
        r, a = r[1:], a[1:]
        p += 1
    return p, r, a

def get_rank(chrom: str, contig_rank: Dict[str, int]) -> int:
    if chrom not in contig_rank:
        contig_rank[chrom] = len(contig_rank)
    return contig_rank[chrom]

def peek_line(iterator):
    try: return next(iterator)
    except StopIteration: return None

def main():
    parser = argparse.ArgumentParser(description="Streaming Minimal Representation Remapper.")
    parser.add_argument('--bubble', required=True, help="Bubble VCF to match against.")
    parser.add_argument('--window', type=int, default=100, help="Coordinate search window for shifting indels.")
    parser.add_argument('--cap-pl', type=int, default=255, help="Maximum value to cap PL/LPL scores.")
    parser.add_argument('--scale-pl', type=float, default=1.0, help="Divide PL/LPL scores by this value before capping.")
    args = parser.parse_args()

    contig_rank = {}
    with open_vcf(args.bubble) as f:
        for line in f:
            if line.startswith('##contig=<ID='):
                cid = line.split('ID=')[1].split(',')[0].split('>')[0]
                if cid not in contig_rank: contig_rank[cid] = len(contig_rank)
            elif line.startswith('#CHROM'): break

    in_iter = iter(sys.stdin)
    bubble_f = open_vcf(args.bubble)
    bubble_iter = iter(bubble_f)

    for line in bubble_iter:
        if line.startswith('#CHROM'): break

    next_bubble_line = peek_line(bubble_iter)
    bubble_buffer: List[Tuple] = []
    out_f = sys.stdout
    
    with tqdm(desc="Remapping Input Variants", unit=" rows", file=sys.stderr, miniters=10000) as pbar:
        for line in in_iter:
            if line.startswith('##'):
                out_f.write(line)
                continue
            if line.startswith('#CHROM'):
                out_f.write('##INFO=<ID=BUBBLE,Number=.,Type=String,Description="Matched Bubble IDs">\n')
                out_f.write('##INFO=<ID=BMAP,Number=A,Type=String,Description="Mapped Bubble Constituent ID">\n')
                out_f.write('##INFO=<ID=BPOS,Number=A,Type=Integer,Description="Original Bubble POS">\n')
                out_f.write('##INFO=<ID=BREF,Number=A,Type=String,Description="Original Bubble REF">\n')
                out_f.write('##INFO=<ID=BALT,Number=A,Type=String,Description="Original Bubble ALT">\n')
                out_f.write(line)
                continue

            pbar.update(1)

            cols = line.rstrip('\n').split('\t')
            i_chrom, i_pos_str, _, i_ref, i_alts_str = cols[:5]
            i_pos = int(i_pos_str)
            i_alts = i_alts_str.split(',')

            while next_bubble_line:
                b_cols = next_bubble_line.split('\t', 2)
                b_chrom = b_cols[0]
                b_pos = int(b_cols[1])

                r_b = get_rank(b_chrom, contig_rank)
                r_i = get_rank(i_chrom, contig_rank)

                if r_b < r_i:
                    next_bubble_line = peek_line(bubble_iter) 
                    continue
                if r_b > r_i or b_pos > i_pos + args.window:
                    break 

                curr_bubble = next_bubble_line 
                next_bubble_line = peek_line(bubble_iter) 

                b_full_cols = curr_bubble.rstrip('\n').split('\t')
                b_pos_str, b_vid, b_ref, b_alts_str, _, _, b_info = b_full_cols[1:8]
                b_alts = b_alts_str.split(',')

                b_ids_raw = extract_info_tag(b_info, 'ID')
                if b_ids_raw and ':' in b_ids_raw:
                    continue 

                b_ids = b_ids_raw.split(',') if b_ids_raw else [b_vid] * len(b_alts)

                for idx, b_alt in enumerate(b_alts):
                    min_p, min_r, min_a = get_minimal_representation(b_pos, b_ref, b_alt)
                    b_cid = b_ids[idx] if idx < len(b_ids) else b_vid
                    bubble_buffer.append((b_chrom, b_pos, b_pos_str, b_vid, b_ref, b_alt, b_cid, min_p, min_r, min_a))

            bubble_buffer = [b for b in bubble_buffer if b[0] == i_chrom and b[1] >= i_pos - args.window]

            matched_bubbles = []
            bmap_array = ["."] * len(i_alts)
            bpos_array = ["."] * len(i_alts)
            bref_array = ["."] * len(i_alts)
            balt_array = ["."] * len(i_alts)

            for alt_idx, i_alt in enumerate(i_alts):
                min_p, min_r, min_a = get_minimal_representation(i_pos, i_ref, i_alt)

                for b_rec in bubble_buffer:
                    if b_rec[7] == min_p and b_rec[8] == min_r and b_rec[9] == min_a:
                        b_pos_str, b_vid, b_ref, b_alt, b_cid = b_rec[2], b_rec[3], b_rec[4], b_rec[5], b_rec[6]
                        
                        if b_vid not in matched_bubbles:
                            matched_bubbles.append(b_vid)
                        
                        bmap_array[alt_idx] = b_cid
                        bpos_array[alt_idx] = b_pos_str
                        bref_array[alt_idx] = b_ref
                        balt_array[alt_idx] = b_alt
                        break 

            if matched_bubbles:
                cols[2] = ",".join(matched_bubbles) 
                info = cols[7]
                new_info_tags = (
                    f"BUBBLE={','.join(matched_bubbles)};"
                    f"BMAP={','.join(bmap_array)};"
                    f"BPOS={','.join(bpos_array)};"
                    f"BREF={','.join(bref_array)};"
                    f"BALT={','.join(balt_array)}"
                )
                cols[7] = new_info_tags if info == '.' else f"{info};{new_info_tags}"

                # --- NEW PL/LPL CAPPING LOGIC ---
                if (args.cap_pl != 255 or not math.isclose(args.scale_pl, 1.0)) and len(cols) > 8:
                    fmt_keys = cols[8].split(':')
                    target_idx = -1
                    if 'LPL' in fmt_keys: target_idx = fmt_keys.index('LPL')
                    elif 'PL' in fmt_keys: target_idx = fmt_keys.index('PL')

                    if target_idx != -1:
                        for s_idx in range(9, len(cols)):
                            s_data = cols[s_idx].split(':')
                            if len(s_data) > target_idx and s_data[target_idx] not in ('.', ''):
                                vals = s_data[target_idx].split(',')
                                # Parse, scale, cap, and replace
                                capped_vals = [str(min(int(int(v) / args.scale_pl), args.cap_pl)) if v != '.' else '.' for v in vals]
                                s_data[target_idx] = ','.join(capped_vals)
                                cols[s_idx] = ':'.join(s_data)
                # --------------------------------

                out_f.write('\t'.join(cols) + '\n')

if __name__ == "__main__":
    main()
