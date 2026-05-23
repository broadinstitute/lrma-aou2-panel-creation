//! =====================================================================
//! VCF VARIANT COLLISION RESOLVER (GREEDY OCCUPANCY ALGORITHM)
//! =====================================================================
//!
//! SUMMARY:
//! This script resolves physical overlapping collisions between variants 
//! in a multi-sample VCF. It ensures that no two variants simultaneously 
//! occupy the same physical DNA space on the same haplotype strand.
//!
//! ALGORITHM:
//! 1. Windowing: The script reads the VCF and groups physically overlapping 
//!    variants into manageable memory "windows".
//! 2. Ranking: Inside each window, variants are ranked by their confidence 
//!    score (Weight -> GQ -> AF -> Variant Class).
//! 3. Anchoring: The highest-scoring phased variant establishes the "Anchor 
//!    Phase Set" (PS) for that window.
//! 4. Greedy Occupancy: Working from highest score to lowest, variants 
//!    attempt to claim physical space (start to end) on Haplotype 1 and/or 2.
//!    - "Aligned" variants (phased, matching the Anchor PS) claim their 
//!      specific strand.
//!    - "Unaligned" variants (unphased or mismatched PS) cannot trust their 
//!      strand, so they must find BOTH strands empty to safely place.
//! 5. Mutation: If a variant is blocked by a higher-scoring variant, its 
//!    ALT allele is degraded to the Reference allele ('0').
//! =====================================================================

use std::env;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::str::FromStr;
use std::cmp::Ordering;

// =====================================================================
// 1. DATA STRUCTURES & SCORING
// =====================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VariantType { Del, Inv, Dup, Ins, Snp, Replacement }

#[derive(Clone, Copy, Debug, PartialEq)]
struct Score {
    w: f64,
    gq: f64,
    af: f64,
    class: i32,
    input_index: usize,
}

impl Eq for Score {}

impl PartialOrd for Score {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Score {
    fn cmp(&self, other: &Self) -> Ordering {
        // Tie-breaker hierarchy: Weight -> GQ -> AF -> Class -> Input Index
        self.w.partial_cmp(&other.w).unwrap_or(Ordering::Equal)
            .then_with(|| self.gq.partial_cmp(&other.gq).unwrap_or(Ordering::Equal))
            .then_with(|| self.af.partial_cmp(&other.af).unwrap_or(Ordering::Equal))
            .then_with(|| self.class.cmp(&other.class))
            .then_with(|| self.input_index.cmp(&other.input_index))
    }
}

// =====================================================================
// 2. CORE INTERVAL LOGIC
// =====================================================================

struct Interval {
    line: String,
    tabs: [usize; 9],
    chr: i32,
    start: i32,
    end: i32,
    v_type: VariantType,
    v_len: i32,
    af: f64,
    class_score: i32,
    modified_gts: Vec<Option<String>>, 
}

impl Interval {
    fn new(line: String, sv_length: i32, n_samples: usize) -> Self {
        let mut tabs = [0; 9];
        let mut tab_idx = 0;
        for (i, &b) in line.as_bytes().iter().enumerate() {
            if b == b'\t' {
                tabs[tab_idx] = i;
                tab_idx += 1;
                if tab_idx == 9 { break; }
            }
        }
        
        let col = |i: usize| -> &str {
            let start = if i == 0 { 0 } else { tabs[i - 1] + 1 };
            &line[start..tabs[i]]
        };

        let chr_str = col(0).to_uppercase();
        let chr = match chr_str.trim_start_matches("CHR") {
            "X" => 23, "Y" => 24, "M" | "MT" => 25,
            s => i32::from_str(s).unwrap_or(-1),
        };

        let pos = i32::from_str(col(1)).unwrap_or(0);
        let info = col(7);
        let ref_allele = col(3);
        let alt_allele = col(4);
        
        let v_type = if let Some(svtype) = get_info_field(info, "SVTYPE") {
            match svtype.to_uppercase().as_str() {
                "DEL" | "DEL:ME" => VariantType::Del,
                "INV" => VariantType::Inv,
                "DUP" | "DUP:TANDEM" | "DUP:INT" | "CNV" => VariantType::Dup,
                "INS" | "INS:ME" | "INS:NOVEL" => VariantType::Ins,
                _ => VariantType::Replacement, 
            }
        } else if ref_allele.len() == 1 {
            if alt_allele.len() > 1 { VariantType::Ins } else { VariantType::Snp }
        } else if alt_allele.len() == 1 { VariantType::Del } else { VariantType::Replacement };
        
        let v_len = if let Some(svlen_str) = get_info_field(info, "SVLEN") {
            i32::from_str(svlen_str).unwrap_or(0).abs()
        } else if v_type == VariantType::Replacement {
            (ref_allele.len() - 1) as i32
        } else {
            (std::cmp::max(ref_allele.len(), alt_allele.len()) - 1) as i32
        };
        
        let (start, end) = match v_type {
            VariantType::Del | VariantType::Inv | VariantType::Dup => (pos + 1, pos + v_len),
            VariantType::Replacement => (pos, pos + v_len),
            VariantType::Ins => (pos, pos + 1),
            VariantType::Snp => (pos, pos),
        };

        let class_score = match v_type {
            VariantType::Del | VariantType::Inv => if v_len >= sv_length { 5 } else { 3 },
            VariantType::Ins | VariantType::Dup => if v_len >= sv_length { 4 } else { 2 },
            VariantType::Snp | VariantType::Replacement => 1,
        };

        let af = get_info_field(info, "AF")
            .and_then(|af_str| af_str.split(',').filter_map(|s| f64::from_str(s).ok()).max_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal)))
            .unwrap_or(0.0);
        
        Interval {
            line, tabs, chr, start, end,
            v_type, v_len, af, class_score,
            modified_gts: vec![None; n_samples],
        }
    }

    fn col(&self, i: usize) -> &str {
        let start = if i == 0 { 0 } else { self.tabs[i - 1] + 1 };
        &self.line[start..self.tabs[i]]
    }

    fn sample_str(&self, sample: usize) -> &str {
        let rest = &self.line[self.tabs[8] + 1..];
        let mut start = 0;
        let mut tabs_seen = 0;
        if sample > 0 {
            for (i, &b) in rest.as_bytes().iter().enumerate() {
                if b == b'\t' {
                    tabs_seen += 1;
                    if tabs_seen == sample {
                        start = i + 1;
                        break;
                    }
                }
            }
        }
        let end = rest[start..].find('\t').map(|i| start + i).unwrap_or(rest.len());
        &rest[start..end]
    }

    fn extract_format_field(&self, sample: usize, tag: &str) -> Option<&str> {
        let format_col = &self.line[self.tabs[7] + 1..self.tabs[8]];
        let target_idx = format_col.split(':').position(|t| t == tag)?;
        let gt = self.sample_str(sample);
        gt.split(':').nth(target_idx)
    }

    fn get_score(&self, sample: usize, weight_tag: &str, in_sample: bool, def_w: f64, idx: usize) -> Score {
        let w_str = if in_sample { self.extract_format_field(sample, weight_tag) } 
                    else { get_info_field(&self.line[self.tabs[6]+1..self.tabs[7]], weight_tag) };
        let w = w_str.and_then(|s| f64::from_str(s).ok()).unwrap_or(def_w);
        let gq = self.extract_format_field(sample, "GQ").and_then(|s| f64::from_str(s).ok()).unwrap_or(0.0);
        Score { w, gq, af: self.af, class: self.class_score, input_index: idx }
    }

    fn gt_bytes(&self, sample: usize) -> &[u8] {
        let s = self.sample_str(sample);
        let colon = s.find(':').unwrap_or(s.len());
        &s.as_bytes()[..colon]
    }

    fn is_missing_or_ref(&self, sample: usize) -> bool { !self.gt_bytes(sample).contains(&b'1') }
    fn is_phased(&self, sample: usize) -> bool { self.gt_bytes(sample).get(1) == Some(&b'|') }
    fn is_unphased(&self, sample: usize) -> bool { self.gt_bytes(sample).get(1) == Some(&b'/') }
    fn has_alt_on_hap1(&self, sample: usize) -> bool { self.gt_bytes(sample).get(0) == Some(&b'1') }
    fn has_alt_on_hap2(&self, sample: usize) -> bool { self.gt_bytes(sample).get(2) == Some(&b'1') }
    
    fn get_ps(&self, sample: usize) -> u32 {
        self.extract_format_field(sample, "PS").map(|s| hash_ps(s)).unwrap_or(0)
    }

    fn category(&self, sv_length: i32) -> usize {
        match self.v_type {
            VariantType::Del | VariantType::Inv => if self.v_len >= sv_length { 0 } else { 1 },
            VariantType::Snp | VariantType::Replacement => 2,
            VariantType::Ins | VariantType::Dup => if self.v_len >= sv_length { 4 } else { 3 },
        }
    }
}

fn hash_ps(s: &str) -> u32 {
    if s == "." || s.is_empty() { return 0; }
    s.as_bytes().iter().fold(5381u32, |hash, &b| hash.wrapping_mul(33).wrapping_add(b as u32)).max(1)
}

fn get_info_field<'a>(info: &'a str, field: &str) -> Option<&'a str> {
    for part in info.split(';') {
        if part.starts_with(field) && part.as_bytes().get(field.len()) == Some(&b'=') {
            return Some(&part[field.len() + 1..]);
        }
    }
    None
}

// =====================================================================
// 3. LOGGING & GRAPH SOLVER
// =====================================================================

fn log_overlap(kept_iv: &Interval, kept_score: &Score, skipped_iv: &Interval, skipped_score: &Score, sample_name: &str) {
    eprintln!("  Kept Allele    - CHROM: {} POS: {} REF: {} ALT: {} SCORE: ({:.3}, {:.1}, {:.3}, {}) ID: {}", 
        kept_iv.col(0), kept_iv.col(1), kept_iv.col(3), kept_iv.col(4), 
        kept_score.w, kept_score.gq, kept_score.af, kept_score.class, kept_iv.col(2));
    eprintln!("  Skipped Allele - CHROM: {} POS: {} REF: {} ALT: {} SCORE: ({:.3}, {:.1}, {:.3}, {}) ID: {}", 
        skipped_iv.col(0), skipped_iv.col(1), skipped_iv.col(3), skipped_iv.col(4), 
        skipped_score.w, skipped_score.gq, skipped_score.af, skipped_score.class, skipped_iv.col(2));
    eprintln!("  SAMPLES: ['{}']", sample_name);
}

struct Candidate {
    idx: usize,
    score: Score,
}

fn process_window<W: Write>(
    window: &mut Vec<Interval>, 
    weight_tag: &str, 
    weight_loc: bool, 
    default_w: f64, 
    output: &mut BufWriter<W>, 
    n_samples: usize,
    sample_names: &[String],
    removed_counts: &mut [[usize; 5]],
    verbosity: u8,
    sv_length: i32,
) -> io::Result<()> {
    if window.is_empty() { return Ok(()); }
    
    for sample in 0..n_samples {
        let mut candidates = Vec::new();
        
        for (idx, iv) in window.iter().enumerate() {
            if !iv.is_missing_or_ref(sample) {
                candidates.push(Candidate {
                    idx,
                    score: iv.get_score(sample, weight_tag, weight_loc, default_w, idx),
                });
            }
        }

        if candidates.is_empty() { continue; }

        candidates.sort_by(|a, b| b.score.cmp(&a.score));

        let mut anchor_ps = 0;
        for c in &candidates {
            let iv = &window[c.idx];
            let ps = iv.get_ps(sample);
            if ps != 0 && iv.is_phased(sample) {
                anchor_ps = ps;
                break;
            }
        }

        let mut hap1_occupied: Vec<(i32, i32, usize)> = Vec::new();
        let mut hap2_occupied: Vec<(i32, i32, usize)> = Vec::new();

        let overlaps = |occupied: &[(i32, i32, usize)], s: i32, e: i32| -> Option<usize> {
            occupied.iter().find(|&&(os, oe, _)| std::cmp::max(s, os) <= std::cmp::min(e, oe)).map(|&(_, _, c_idx)| c_idx)
        };

        for (c_idx, c) in candidates.iter().enumerate() {
            // Immutable borrow for checking conditions and logging
            let iv = &window[c.idx];
            let start = iv.start;
            let end = iv.end;
            
            let is_aligned = iv.is_phased(sample) && (iv.get_ps(sample) == anchor_ps || anchor_ps == 0);
            let has_alt_h1 = iv.has_alt_on_hap1(sample);
            let has_alt_h2 = iv.has_alt_on_hap2(sample);
            
            let hit1 = overlaps(&hap1_occupied, start, end);
            let hit2 = overlaps(&hap2_occupied, start, end);
            
            let mut keep_h1 = false;
            let mut keep_h2 = false;

            if is_aligned {
                keep_h1 = has_alt_h1 && hit1.is_none();
                keep_h2 = has_alt_h2 && hit2.is_none();
                
                if keep_h1 { hap1_occupied.push((start, end, c_idx)); }
                if keep_h2 { hap2_occupied.push((start, end, c_idx)); }
            } else {
                if has_alt_h1 && has_alt_h2 {
                    keep_h1 = hit1.is_none();
                    keep_h2 = hit2.is_none();
                    if keep_h1 { hap1_occupied.push((start, end, c_idx)); }
                    if keep_h2 { hap2_occupied.push((start, end, c_idx)); }
                } else if has_alt_h1 || has_alt_h2 {
                    if hit1.is_none() && hit2.is_none() {
                        keep_h1 = has_alt_h1;
                        keep_h2 = has_alt_h2;
                        hap1_occupied.push((start, end, c_idx));
                        hap2_occupied.push((start, end, c_idx));
                    }
                }
            }

            // Apply Mutations and Logging
            if (has_alt_h1 && !keep_h1) || (has_alt_h2 && !keep_h2) {
                
                removed_counts[sample][iv.category(sv_length)] += 1;
                
                if verbosity >= 2 {
                    if let Some(b_idx) = hit1.or(hit2) {
                        let blocker_c = &candidates[b_idx];
                        log_overlap(&window[blocker_c.idx], &blocker_c.score, iv, &c.score, &sample_names[sample]);
                    }
                }

                let gt = iv.gt_bytes(sample);
                let mut new_a = gt[0];
                let mut new_b = gt[2];

                if iv.is_unphased(sample) {
                    if !keep_h1 && !keep_h2 {
                        new_a = if gt[0] == b'.' { b'.' } else { b'0' };
                        new_b = if gt[2] == b'.' { b'.' } else { b'0' };
                    }
                } else {
                    if !keep_h1 { new_a = if gt[0] == b'.' { b'.' } else { b'0' }; }
                    if !keep_h2 { new_b = if gt[2] == b'.' { b'.' } else { b'0' }; }
                }

                if new_a != gt[0] || new_b != gt[2] {
                    let s_str = iv.sample_str(sample);
                    let new_gt = format!("{}{}{}{}", new_a as char, gt[1] as char, new_b as char, &s_str[gt.len()..]);
                    
                    // Single mutable borrow required exclusively to write the change
                    window[c.idx].modified_gts[sample] = Some(new_gt);
                }
            }
        }
    }

    // Output Generation
    for iv in window {
        output.write_all(&iv.line.as_bytes()[..iv.tabs[8]])?;
        let mut samples_iter = iv.line[iv.tabs[8] + 1..].split('\t');
        for sample in 0..n_samples {
            output.write_all(b"\t")?;
            if let Some(ref mod_gt) = iv.modified_gts[sample] {
                output.write_all(mod_gt.as_bytes())?;
                samples_iter.next(); 
            } else {
                output.write_all(samples_iter.next().unwrap().as_bytes())?;
            }
        }
        output.write_all(b"\n")?;
    }
    Ok(())
}

// =====================================================================
// 4. MAIN EXECUTION I/O
// =====================================================================

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    let mut sv_length = 50;
    let mut verbosity: u8 = 0;
    let mut clean_args = Vec::new();
    
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--sv-length" => { i += 1; sv_length = args[i].parse().unwrap_or(50); }
            "--verbosity" => { i += 1; verbosity = args[i].parse().unwrap_or(0); }
            _ => clean_args.push(args[i].clone()),
        }
        i += 1;
    }

    if clean_args.len() < 4 {
        eprintln!("Usage: {} [--sv-length <int>] [--verbosity <0|1|2>] <method> <weight_tag> <weight_loc> <default_w> < input.vcf > output.vcf", args[0]);
        return Ok(());
    }

    let weight_tag = &clean_args[1];
    let weight_loc = clean_args[2] == "1";
    let default_w = f64::from_str(&clean_args[3]).unwrap_or(0.0);

    let stdin = io::stdin();
    let mut reader = BufReader::new(stdin.lock());
    let stdout = io::stdout();
    let mut output = BufWriter::new(stdout.lock());

    let mut window: Vec<Interval> = Vec::new();
    let mut window_max_end = -1;
    let mut n_samples = 0;
    let mut sample_names: Vec<String> = Vec::new();
    let mut removed_per_sample: Vec<[usize; 5]> = Vec::new();
    let mut line = String::new();

    while reader.read_line(&mut line)? > 0 {
        if line.ends_with('\n') { line.pop(); }
        if line.ends_with('\r') { line.pop(); }

        if line.starts_with('#') {
            if line.starts_with("#CHROM") { 
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() > 9 {
                    sample_names = fields[9..].iter().map(|s| s.to_string()).collect();
                    n_samples = sample_names.len();
                    removed_per_sample = vec![[0; 5]; n_samples];
                }
            }
            writeln!(output, "{}", line)?;
        } else {
            let iv = Interval::new(line.clone(), sv_length, n_samples);
            
            if window.is_empty() || (iv.chr == window[0].chr && iv.start <= window_max_end) {
                window_max_end = window_max_end.max(iv.end);
                window.push(iv);
            } else {
                process_window(&mut window, weight_tag, weight_loc, default_w, &mut output, n_samples, &sample_names, &mut removed_per_sample, verbosity, sv_length)?;
                window.clear();
                window_max_end = iv.end;
                window.push(iv);
            }
        }
        line.clear();
    }
    
    process_window(&mut window, weight_tag, weight_loc, default_w, &mut output, n_samples, &sample_names, &mut removed_per_sample, verbosity, sv_length)?;
    output.flush()?;

    // Print summary table if Verbosity >= 1
    if verbosity >= 1 {
        eprintln!("\nSAMPLE\tSV_DEL\tDEL\tSNP\tINS\tSV_INS\tTOTAL");
        for (i, counts) in removed_per_sample.iter().enumerate() {
            let total: usize = counts.iter().sum();
            eprintln!("{}\t{}\t{}\t{}\t{}\t{}\t{}", sample_names[i], counts[0], counts[1], counts[2], counts[3], counts[4], total);
        }
    }

    Ok(())
}
