use std::collections::BTreeMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VariantType { Del, Inv, Dup, Ins, Snp, Replacement }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
enum Genotype {
    Phased00 = 0, Phased01 = 1, Phased10 = 2, Phased11 = 3,
    Unphased00 = 4, Unphased01 = 5, Unphased10 = 6, Unphased11 = 7,
    PhasedD0 = 8, Phased0D = 9, PhasedD1 = 10, Phased1D = 11, PhasedDD = 12,
    UnphasedD0 = 13, Unphased0D = 14, UnphasedD1 = 15, Unphased1D = 16, UnphasedDD = 17,
}

impl Genotype {
    fn from_str(gt: &str) -> Self {
        let bytes = gt.as_bytes();
        if bytes.len() < 3 { return Genotype::UnphasedDD; }
        let a = bytes[0];
        let b = bytes[1];
        let c = bytes[2];

        if b == b'/' {
            match (a, c) {
                (b'.', b'.') => Genotype::UnphasedDD,
                (b'.', b'0') => Genotype::UnphasedD0,
                (b'.', b'1') => Genotype::UnphasedD1,
                (b'0', b'.') => Genotype::Unphased0D,
                (b'0', b'0') => Genotype::Unphased00,
                (b'0', b'1') => Genotype::Unphased01,
                (b'1', b'.') => Genotype::Unphased1D,
                (b'1', b'0') => Genotype::Unphased10,
                (b'1', b'1') | _ => Genotype::Unphased11,
            }
        } else {
            match (a, c) {
                (b'.', b'.') => Genotype::PhasedDD,
                (b'.', b'0') => Genotype::PhasedD0,
                (b'.', b'1') => Genotype::PhasedD1,
                (b'0', b'.') => Genotype::Phased0D,
                (b'0', b'0') => Genotype::Phased00,
                (b'0', b'1') => Genotype::Phased01,
                (b'1', b'.') => Genotype::Phased1D,
                (b'1', b'0') => Genotype::Phased10,
                (b'1', b'1') | _ => Genotype::Phased11,
            }
        }
    }

    fn to_str(self) -> &'static str {
        match self {
            Genotype::Phased00 => "0|0", Genotype::Phased01 => "0|1",
            Genotype::Phased10 => "1|0", Genotype::Phased11 => "1|1",
            Genotype::Unphased00 => "0/0", Genotype::Unphased01 => "0/1",
            Genotype::Unphased10 => "1/0", Genotype::Unphased11 => "1/1",
            Genotype::PhasedD0 => ".|0", Genotype::Phased0D => "0|.",
            Genotype::PhasedD1 => ".|1", Genotype::Phased1D => "1|.", Genotype::PhasedDD => ".|.",
            Genotype::UnphasedD0 => "./0", Genotype::Unphased0D => "0/.",
            Genotype::UnphasedD1 => "./1", Genotype::Unphased1D => "1/.", Genotype::UnphasedDD => "./.",
        }
    }

    fn to_unphased(self) -> Self {
        match self {
            Genotype::Phased00 => Genotype::Unphased00,
            Genotype::Phased01 => Genotype::Unphased01,
            Genotype::Phased10 => Genotype::Unphased10,
            Genotype::Phased11 => Genotype::Unphased11,
            Genotype::PhasedD0 => Genotype::UnphasedD0,
            Genotype::Phased0D => Genotype::Unphased0D,
            Genotype::PhasedD1 => Genotype::UnphasedD1,
            Genotype::Phased1D => Genotype::Unphased1D,
            Genotype::PhasedDD => Genotype::UnphasedDD,
            _ => self,
        }
    }
}

const N_GT_COLLISIONS: [[u8; 18]; 18] = [
    [0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0],
    [0,1,0,1, 0,0,0,1, 0,0,1,0,0, 0,0,0,0,0],
    [0,0,1,1, 0,0,0,1, 0,0,0,1,0, 0,0,0,0,0],
    [0,1,1,2, 0,0,0,2, 0,0,1,1,0, 0,0,1,1,0],
    [0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0],
    [0,0,0,1, 0,0,0,1, 0,0,0,0,0, 0,0,0,0,0],
    [0,0,0,1, 0,0,0,1, 0,0,0,0,0, 0,0,0,0,0],
    [0,1,1,2, 0,0,0,2, 0,0,1,1,0, 0,0,1,1,0],
    [0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0],
    [0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0],
    [0,1,0,1, 0,0,0,1, 0,0,1,0,0, 0,0,0,0,0],
    [0,0,1,1, 0,0,0,1, 0,0,0,1,0, 0,0,0,0,0],
    [0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0],
    [0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0],
    [0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0],
    [0,0,0,1, 0,0,0,1, 0,0,0,0,0, 0,0,0,0,0],
    [0,0,0,1, 0,0,0,1, 0,0,0,0,0, 0,0,0,0,0],
    [0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0],
];

use Genotype::*;
const REMOVE_HAP1_0: [Genotype; 18] = [Phased00,Phased01,Phased00,Phased01,Unphased00,Unphased01,Unphased10,Unphased01,PhasedD0,Phased0D,PhasedD1,Phased0D,PhasedDD,UnphasedD0,Unphased0D,UnphasedD1,Unphased1D,UnphasedDD];
const REMOVE_HAP1_D: [Genotype; 18] = [Phased00,Phased01,PhasedD0,PhasedD1,Unphased00,Unphased01,Unphased10,UnphasedD1,PhasedD0,Phased0D,PhasedD1,PhasedDD,PhasedDD,UnphasedD0,Unphased0D,UnphasedD1,Unphased1D,UnphasedDD];
const REMOVE_HAP2_0: [Genotype; 18] = [Phased00,Phased00,Phased10,Phased10,Unphased00,Unphased01,Unphased10,Unphased01,PhasedD0,Phased0D,PhasedD0,Phased1D,PhasedDD,UnphasedD0,Unphased0D,UnphasedD1,Unphased1D,UnphasedDD];
const REMOVE_HAP2_D: [Genotype; 18] = [Phased00,Phased0D,Phased10,Phased1D,Unphased00,Unphased01,Unphased10,UnphasedD1,PhasedD0,Phased0D,PhasedDD,Phased1D,PhasedDD,UnphasedD0,Unphased0D,UnphasedD1,Unphased1D,UnphasedDD];

fn hash_ps(s: &str) -> u32 {
    if s == "." || s.is_empty() { return 0; }
    let mut hash = 5381u32;
    for &b in s.as_bytes() {
        hash = hash.wrapping_mul(33).wrapping_add(b as u32);
    }
    if hash == 0 { 1 } else { hash } 
}

#[derive(Clone)]
struct Interval {
    line: String,
    tabs: [usize; 9],
    chr: i32,
    first: i32,
    last: i32,
    input_index: usize,
    
    v_type: VariantType,
    v_len: i32,
    weight: f64,
    gq: f64,
    af: f64,
    genotypes: Vec<Genotype>,
    phase_sets: Vec<u32>,
    
    in_independent_set: bool,
    independent_set_weight: f64,
    independent_set_gq: f64,
    independent_set_af: f64,
    independent_set_previous: Option<usize>,
    overlaps_is_hap1: bool,
    overlaps_is_hap2: bool,
}

impl Interval {
    fn new(line: String) -> Self {
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
            let end = tabs[i];
            &line[start..end]
        };

        let chr = parse_chr(col(0));
        let pos = i32::from_str(col(1)).unwrap_or(0);
        let info = col(7);
        let ref_allele = col(3);
        let alt_allele = col(4);
        
        let variant_type = if let Some(svtype) = get_info_field(info, "SVTYPE") {
            parse_svtype(svtype)
        } else if ref_allele.len() == 1 {
            if alt_allele.len() > 1 { VariantType::Ins } else { VariantType::Snp }
        } else if alt_allele.len() == 1 {
            VariantType::Del
        } else {
            VariantType::Replacement
        };
        
        let length = if let Some(svlen_str) = get_info_field(info, "SVLEN") {
            i32::from_str(svlen_str).unwrap_or(0).abs()
        } else if variant_type == VariantType::Replacement {
            (ref_allele.len() - 1) as i32
        } else {
            (std::cmp::max(ref_allele.len(), alt_allele.len()) - 1) as i32
        };
        
        let (first, last) = match variant_type {
            VariantType::Del | VariantType::Inv | VariantType::Dup => (pos + 1, pos + length),
            VariantType::Replacement => (pos, pos + length),
            VariantType::Ins => (pos, pos + 1),
            VariantType::Snp => (pos, pos),
        };

        let af = if let Some(af_str) = get_info_field(info, "AF") {
            af_str.split(',')
                .filter_map(|s| f64::from_str(s).ok())
                .fold(0.0, f64::max)
        } else {
            0.0
        };
        
        let mut genotypes = Vec::new();
        let mut phase_sets = Vec::new();
        
        if tab_idx == 9 {
            let format_col = col(8);
            let mut ps_idx = None;
            for (i, t) in format_col.split(':').enumerate() {
                if t == "PS" { ps_idx = Some(i); break; }
            }

            for sample_str in line[tabs[8] + 1..].split('\t') {
                genotypes.push(Genotype::from_str(sample_str));
                let mut ps = 0;
                
                if let Some(target) = ps_idx {
                    let mut current_j = 0;
                    let mut start = 0;
                    let bytes = sample_str.as_bytes();
                    for i in 0..bytes.len() {
                        if bytes[i] == b':' {
                            if current_j == target {
                                ps = hash_ps(&sample_str[start..i]);
                                break;
                            }
                            current_j += 1;
                            start = i + 1;
                        }
                    }
                    if current_j == target && start <= sample_str.len() && ps == 0 {
                        ps = hash_ps(&sample_str[start..]);
                    }
                }
                phase_sets.push(ps);
            }
        }
        
        Interval {
            line, tabs,
            chr, first, last, input_index: 0, 
            v_type: variant_type, v_len: length,
            weight: 0.0, gq: 0.0, af, genotypes, phase_sets,
            in_independent_set: false, independent_set_weight: 0.0, independent_set_gq: 0.0, independent_set_af: 0.0,
            independent_set_previous: None, overlaps_is_hap1: false, overlaps_is_hap2: false,
        }
    }

    fn col(&self, i: usize) -> &str {
        if i == 0 {
            &self.line[..self.tabs[0]]
        } else if i < 9 {
            &self.line[self.tabs[i - 1] + 1..self.tabs[i]]
        } else {
            ""
        }
    }

    fn col_info(&self) -> &str {
        self.col(7)
    }

    fn col_format(&self) -> &str {
        self.col(8)
    }

    fn sample_str(&self, sample_idx: usize) -> &str {
        let rest = &self.line[self.tabs[8] + 1..];
        let mut start = 0;
        let mut tabs_seen = 0;
        if sample_idx > 0 {
            for (i, &b) in rest.as_bytes().iter().enumerate() {
                if b == b'\t' {
                    tabs_seen += 1;
                    if tabs_seen == sample_idx {
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
        let format_col = self.col_format();
        let mut tag_index = None;
        for (i, t) in format_col.split(':').enumerate() {
            if t == tag {
                tag_index = Some(i);
                break;
            }
        }
        let target_idx = tag_index?;
        
        let gt = self.sample_str(sample);
        let mut current_j = 0;
        let mut start = 0;
        let gt_bytes = gt.as_bytes();
        
        for i in 0..gt_bytes.len() {
            if gt_bytes[i] == b':' {
                if current_j == target_idx {
                    return Some(&gt[start..i]);
                }
                current_j += 1;
                start = i + 1;
            }
        }
        
        if current_j == target_idx && start <= gt.len() {
            return Some(&gt[start..]);
        }
        None
    }

    fn set_scores(&mut self, sample: usize, weight_tag: &str, in_sample: bool, default_weight: f64) {
        let val_opt = if in_sample {
            self.extract_format_field(sample, weight_tag)
        } else {
            get_info_field(self.col_info(), weight_tag)
        };
        
        self.weight = val_opt.and_then(|s| f64::from_str(s).ok()).unwrap_or(default_weight);
        self.gq = self.extract_format_field(sample, "GQ").and_then(|s| f64::from_str(s).ok()).unwrap_or(0.0);
    }

    fn is_present(&self, sample: usize) -> bool {
        let gt = self.genotypes[sample];
        !matches!(gt, Phased00 | Unphased00 | PhasedDD | UnphasedDD | PhasedD0 | UnphasedD0 | Phased0D | Unphased0D)
    }

    fn on_hap1(&self, sample: usize) -> bool {
        matches!(self.genotypes[sample], Phased10 | Phased1D | Phased11 | Unphased11)
    }

    fn on_hap2(&self, sample: usize) -> bool {
        matches!(self.genotypes[sample], Phased01 | PhasedD1 | Phased11 | Unphased11)
    }

    fn precedes(&self, next: &Interval, sample: usize) -> bool {
        let mut gt_self = self.genotypes[sample];
        let mut gt_next = next.genotypes[sample];
        if self.phase_sets[sample] != next.phase_sets[sample] {
            gt_self = gt_self.to_unphased();
            gt_next = gt_next.to_unphased();
        }
        N_GT_COLLISIONS[gt_self as usize][gt_next as usize] == 0 || self.last < next.first
    }

    fn clear_is_variables(&mut self) {
        self.in_independent_set = false;
        self.independent_set_weight = 0.0;
        self.independent_set_gq = 0.0;
        self.independent_set_af = 0.0;
        self.independent_set_previous = None;
        self.overlaps_is_hap1 = false;
        self.overlaps_is_hap2 = false;
    }

    fn write_vcf<W: Write>(&self, out: &mut W) -> io::Result<()> {
        let format_str = self.col_format();
        let mut gt_index = 0;
        if let Some(p) = format_str.find("GT") {
            for byte in format_str[..p].bytes() {
                if byte == b':' { gt_index += 1; }
            }
        }

        out.write_all(self.line[..self.tabs[8]].as_bytes())?;

        let mut samples_iter = self.line[self.tabs[8] + 1..].split('\t');
        for gt_enum in &self.genotypes {
            out.write_all(b"\t")?;
            let sample = samples_iter.next().unwrap();
            let sample_bytes = sample.as_bytes();
            let mut k = 0;
            let mut p2 = 0;
            
            for j in 0..sample_bytes.len() {
                if sample_bytes[j] == b':' {
                    k += 1;
                    if k == gt_index {
                        p2 = j + 1;
                        break;
                    }
                }
            }
            
            if p2 > 0 {
                out.write_all(&sample_bytes[0..p2])?;
            }
            out.write_all(gt_enum.to_str().as_bytes())?;
            
            if let Some(q) = sample[p2..].find(':') {
                out.write_all(sample[p2 + q..].as_bytes())?;
            }
        }
        out.write_all(b"\n")?;
        Ok(())
    }

    fn category(&self, sv_length: i32) -> usize {
        match self.v_type {
            VariantType::Snp | VariantType::Replacement => 2, // SNP
            VariantType::Del | VariantType::Inv => if self.v_len >= sv_length { 0 } else { 1 }, // SV_DEL or DEL
            VariantType::Ins | VariantType::Dup => if self.v_len >= sv_length { 4 } else { 3 }, // SV_INS or INS
        }
    }
}

fn is_better(w1: f64, gq1: f64, af1: f64, w2: f64, gq2: f64, af2: f64, use_gq: bool, use_af: bool) -> bool {
    if w1 > w2 + 1e-9 { return true; }
    if (w1 - w2).abs() <= 1e-9 {
        if use_gq {
            if gq1 > gq2 + 1e-9 { return true; }
            if (gq1 - gq2).abs() > 1e-9 { return false; }
        }
        if use_af {
            if af1 > af2 + 1e-9 { return true; }
        }
    }
    false
}

fn is_equal(w1: f64, gq1: f64, af1: f64, w2: f64, gq2: f64, af2: f64, use_gq: bool, use_af: bool) -> bool {
    if (w1 - w2).abs() > 1e-9 { return false; }
    if use_gq && (gq1 - gq2).abs() > 1e-9 { return false; }
    if use_af && (af1 - af2).abs() > 1e-9 { return false; }
    true
}

fn get_info_field<'a>(info: &'a str, field: &str) -> Option<&'a str> {
    let mut start = 0;
    let info_bytes = info.as_bytes();
    while start < info.len() {
        if info[start..].starts_with(field) {
            let after_field = start + field.len();
            if after_field < info.len() && info_bytes[after_field] == b'=' {
                let val_start = after_field + 1;
                let val_end = info[val_start..].find(';').map(|i| i + val_start).unwrap_or(info.len());
                return Some(&info[val_start..val_end]);
            }
        }
        start = info[start..].find(';').map(|i| start + i + 1).unwrap_or(info.len());
    }
    None
}

fn parse_chr(chr_str: &str) -> i32 {
    match chr_str.to_uppercase().trim_start_matches("CHR") {
        "X" => 23, "Y" => 24, "M" | "MT" => 25,
        s => i32::from_str(s).unwrap_or(-1),
    }
}

fn parse_svtype(t: &str) -> VariantType {
    match t.to_uppercase().as_str() {
        "DEL" | "DEL:ME" => VariantType::Del,
        "INV" => VariantType::Inv,
        "DUP" | "DUP:TANDEM" | "DUP:INT" | "CNV" => VariantType::Dup,
        "INS" | "INS:ME" | "INS:NOVEL" => VariantType::Ins,
        _ => VariantType::Replacement, 
    }
}

fn is_last_in_window(window: &[Interval], window_last_pos: i32) -> bool {
    if window.is_empty() || window_last_pos == -1 { return true; }
    let first = &window[0];
    let last = &window[window.len() - 1];
    last.chr == first.chr && last.first <= window_last_pos
}

fn count_collisions(window: &[Interval], sample: usize) -> usize {
    let mut out = 0;
    for i in (0..window.len()).rev() {
        if !window[i].is_present(sample) { continue; }
        for j in (0..i).rev() {
            if window[j].last < window[i].first { break; }
            let mut gt_i = window[i].genotypes[sample];
            let mut gt_j = window[j].genotypes[sample];
            if window[i].phase_sets[sample] != window[j].phase_sets[sample] {
                gt_i = gt_i.to_unphased();
                gt_j = gt_j.to_unphased();
            }
            out += N_GT_COLLISIONS[gt_j as usize][gt_i as usize] as usize;
        }
    }
    out
}

fn log_overlap(kept: &Interval, skipped: &Interval, sample_name: &str) {
    eprintln!("  Kept Allele    - CHROM: {} POS: {} REF: {} ALT: {} WEIGHT: ({}, {}, {}) ID: {}", kept.col(0), kept.col(1), kept.col(3), kept.col(4), kept.weight, kept.gq, kept.af, kept.col(2));
    eprintln!("  Skipped Allele - CHROM: {} POS: {} REF: {} ALT: {} WEIGHT: ({}, {}, {}) ID: {}", skipped.col(0), skipped.col(1), skipped.col(3), skipped.col(4), skipped.weight, skipped.gq, skipped.af, skipped.col(2));
    eprintln!("  SAMPLES: ['{}']", sample_name);
}

fn mark_is_overlaps(
    window: &mut [Interval], 
    sample_window: &[usize], 
    sample: usize, 
    hap1: bool, 
    hap2: bool,
    sample_name: &str,
    already_logged: &mut [bool],
    verbosity: u8,
    sv_length: i32
) -> [usize; 5] {
    already_logged.fill(false);
    let mut removed_counts = [0; 5];

    for &idx in sample_window {
        window[idx].overlaps_is_hap1 = false;
        window[idx].overlaps_is_hap2 = false;
    }

    for i in (0..sample_window.len()).rev() {
        let idx_i = sample_window[i];
        let first_i = window[idx_i].first;
        let on_hap1_i = window[idx_i].on_hap1(sample);
        let on_hap2_i = window[idx_i].on_hap2(sample);
        let in_is_i = window[idx_i].in_independent_set;

        for j in (0..i).rev() {
            let idx_j = sample_window[j];
            if window[idx_j].last < first_i { break; }
            let in_is_j = window[idx_j].in_independent_set;
            let same_ps = window[idx_i].phase_sets[sample] == window[idx_j].phase_sets[sample];

            let mut log_skip = |w: &[Interval], logged: &mut [bool], kept_idx: usize, skipped_idx: usize, counts: &mut [usize; 5]| {
                if !logged[skipped_idx] {
                    if verbosity >= 2 {
                        log_overlap(&w[kept_idx], &w[skipped_idx], sample_name);
                    }
                    logged[skipped_idx] = true;
                    counts[w[skipped_idx].category(sv_length)] += 1;
                }
            };

            if in_is_i && !in_is_j {
                if same_ps && hap1 && on_hap1_i && window[idx_j].on_hap1(sample) { 
                    window[idx_j].overlaps_is_hap1 = true; 
                    log_skip(window, already_logged, idx_i, idx_j, &mut removed_counts);
                }
                if same_ps && hap2 && on_hap2_i && window[idx_j].on_hap2(sample) { 
                    window[idx_j].overlaps_is_hap2 = true; 
                    log_skip(window, already_logged, idx_i, idx_j, &mut removed_counts);
                }
            } else if in_is_j && !in_is_i {
                if same_ps && hap1 && window[idx_j].on_hap1(sample) && on_hap1_i { 
                    window[idx_i].overlaps_is_hap1 = true; 
                    log_skip(window, already_logged, idx_j, idx_i, &mut removed_counts);
                }
                if same_ps && hap2 && window[idx_j].on_hap2(sample) && on_hap2_i { 
                    window[idx_i].overlaps_is_hap2 = true; 
                    log_skip(window, already_logged, idx_j, idx_i, &mut removed_counts);
                }
            }
        }
    }
    
    removed_counts
}

fn independent_set_1(
    window: &mut [Interval], 
    sample: usize, 
    weight_tag: &str, 
    weight_loc: bool, 
    default_w: f64, 
    sample_name: &str,
    already_logged: &mut [bool],
    use_gq: bool,
    use_af: bool,
    verbosity: u8,
    sv_length: i32
) -> [usize; 5] {
    let mut sample_window = Vec::new();
    for i in 0..window.len() {
        if window[i].is_present(sample) {
            window[i].clear_is_variables();
            window[i].set_scores(sample, weight_tag, weight_loc, default_w);
            sample_window.push(i);
        }
    }

    if sample_window.is_empty() { return [0; 5]; }

    let mut max_weight = f64::NEG_INFINITY;
    let mut max_gq = f64::NEG_INFINITY;
    let mut max_af = f64::NEG_INFINITY;
    let mut active = vec![true; sample_window.len()];
    let mut best_is = Vec::new();

    for i in (0..sample_window.len()).rev() {
        let mut ancestors = Vec::new();
        is1_dfs(
            i, &mut active, &mut ancestors, 
            0.0, 0.0, 0.0, 
            sample, &sample_window, window, 
            &mut max_weight, &mut max_gq, &mut max_af, &mut best_is,
            use_gq, use_af
        );
    }

    for &idx in &best_is {
        window[sample_window[idx]].in_independent_set = true;
    }

    let removed = mark_is_overlaps(window, &sample_window, sample, true, true, sample_name, already_logged, verbosity, sv_length);

    for &idx in &sample_window {
        if window[idx].in_independent_set { continue; }
        
        let mut gt = window[idx].genotypes[sample];
        gt = if window[idx].overlaps_is_hap1 { REMOVE_HAP1_D[gt as usize] } else { REMOVE_HAP1_0[gt as usize] };
        gt = if window[idx].overlaps_is_hap2 { REMOVE_HAP2_D[gt as usize] } else { REMOVE_HAP2_0[gt as usize] };
        window[idx].genotypes[sample] = gt;
    }
    
    removed
}

fn is1_dfs(
    id: usize, active: &mut [bool], ancestors: &mut Vec<usize>, 
    ancestors_weight: f64, ancestors_gq: f64, ancestors_af: f64,
    sample: usize, sample_window: &[usize], window: &[Interval],
    max_weight: &mut f64, max_gq: &mut f64, max_af: &mut f64, best_is: &mut Vec<usize>,
    use_gq: bool, use_af: bool
) {
    let curr_idx = sample_window[id];
    let weight_prime = ancestors_weight + window[curr_idx].weight;
    let gq_prime = ancestors_gq + window[curr_idx].gq;
    let af_prime = ancestors_af + window[curr_idx].af;
    let mut active_prime = active.to_vec();
    active_prime[id] = false;
    
    let mut has_child = false;
    let mut upper_bound = weight_prime;
    
    for i in (0..id).rev() {
        if active_prime[i] && !window[sample_window[i]].precedes(&window[curr_idx], sample) {
            active_prime[i] = false;
        }
        if active_prime[i] {
            has_child = true;
            upper_bound += window[sample_window[i]].weight;
        }
    }
    
    if upper_bound < *max_weight - 1e-9 { return; }
    
    if has_child {
        ancestors.push(id);
        for i in (0..id).rev() {
            if !active_prime[i] { continue; }
            if upper_bound < *max_weight - 1e-9 { break; }
            is1_dfs(
                i, &mut active_prime, ancestors, 
                weight_prime, gq_prime, af_prime, 
                sample, sample_window, window, 
                max_weight, max_gq, max_af, best_is,
                use_gq, use_af
            );
            upper_bound -= window[sample_window[i]].weight;
        }
        ancestors.pop();
    } else if is_better(weight_prime, gq_prime, af_prime, *max_weight, *max_gq, *max_af, use_gq, use_af) {
        *max_weight = weight_prime;
        *max_gq = gq_prime;
        *max_af = af_prime;
        best_is.clear();
        best_is.extend(ancestors.iter().copied());
        best_is.push(id);
    }
}

#[derive(Default)]
struct EndpointGroup {
    opens: Vec<usize>,
    closed: Vec<usize>,
}

fn independent_set_2(
    window: &mut [Interval], 
    sample: usize, 
    weight_tag: &str, 
    weight_loc: bool, 
    default_w: f64, 
    sample_name: &str,
    already_logged: &mut [bool],
    use_gq: bool,
    use_af: bool,
    verbosity: u8,
    sv_length: i32
) -> [usize; 5] {
    let mut total_removed = [0; 5];
    
    for hap in [1, 2] {
        let mut sample_window = Vec::new();
        let mut points_map: BTreeMap<i32, EndpointGroup> = BTreeMap::new();
        
        for i in 0..window.len() {
            let on_hap = if hap == 1 { window[i].on_hap1(sample) } else { window[i].on_hap2(sample) };
            if on_hap {
                window[i].clear_is_variables();
                window[i].set_scores(sample, weight_tag, weight_loc, default_w);
                sample_window.push(i);
                points_map.entry(window[i].first).or_default().opens.push(i);
                points_map.entry(window[i].last).or_default().closed.push(i);
            }
        }

        if sample_window.is_empty() { continue; }

        let mut available_weight = f64::NEG_INFINITY;
        let mut available_gq = f64::NEG_INFINITY;
        let mut available_af = f64::NEG_INFINITY;
        let mut previous: Option<usize> = None;

        for (_pos, group) in points_map {
            for &idx in &group.opens {
                let base_w = if available_weight == f64::NEG_INFINITY { 0.0 } else { available_weight };
                let base_gq = if available_gq == f64::NEG_INFINITY { 0.0 } else { available_gq };
                let base_af = if available_af == f64::NEG_INFINITY { 0.0 } else { available_af };
                
                window[idx].independent_set_weight = base_w + window[idx].weight;
                window[idx].independent_set_gq = base_gq + window[idx].gq;
                window[idx].independent_set_af = base_af + window[idx].af;
                window[idx].independent_set_previous = previous;
            }
            for &idx in &group.closed {
                if is_better(
                    window[idx].independent_set_weight, window[idx].independent_set_gq, window[idx].independent_set_af,
                    available_weight, available_gq, available_af, use_gq, use_af
                ) {
                    available_weight = window[idx].independent_set_weight;
                    available_gq = window[idx].independent_set_gq;
                    available_af = window[idx].independent_set_af;
                    previous = Some(idx);
                }
            }
        }

        let mut curr_trace = None;
        for i in (0..sample_window.len()).rev() {
            let idx = sample_window[i];
            if is_equal(
                window[idx].independent_set_weight, window[idx].independent_set_gq, window[idx].independent_set_af,
                available_weight, available_gq, available_af, use_gq, use_af
            ) {
                curr_trace = Some(idx);
                break;
            }
        }

        while let Some(idx) = curr_trace {
            window[idx].in_independent_set = true;
            curr_trace = window[idx].independent_set_previous;
        }

        if hap == 1 {
            let removed = mark_is_overlaps(window, &sample_window, sample, true, false, sample_name, already_logged, verbosity, sv_length);
            for i in 0..5 { total_removed[i] += removed[i]; }
            for &idx in &sample_window {
                if !window[idx].in_independent_set {
                    let gt = window[idx].genotypes[sample];
                    window[idx].genotypes[sample] = if window[idx].overlaps_is_hap1 { REMOVE_HAP1_D[gt as usize] } else { REMOVE_HAP1_0[gt as usize] };
                }
            }
        } else {
            let removed = mark_is_overlaps(window, &sample_window, sample, false, true, sample_name, already_logged, verbosity, sv_length);
            for i in 0..5 { total_removed[i] += removed[i]; }
            for &idx in &sample_window {
                if !window[idx].in_independent_set {
                    let gt = window[idx].genotypes[sample];
                    window[idx].genotypes[sample] = if window[idx].overlaps_is_hap2 { REMOVE_HAP2_D[gt as usize] } else { REMOVE_HAP2_0[gt as usize] };
                }
            }
        }
    }
    
    total_removed
}

fn process_window<W: Write>(
    window: &mut Vec<Interval>, 
    method_2: bool, 
    weight_tag: &str, 
    weight_loc: bool, 
    default_w: f64, 
    output: &mut BufWriter<W>, 
    histogram: &mut [usize],
    sample_names: &[String],
    already_logged: &mut Vec<bool>,
    removed_per_sample: &mut [[usize; 5]],
    use_gq: bool,
    use_af: bool,
    verbosity: u8,
    sv_length: i32
) {
    if window.is_empty() { return; }
    window.sort_by(|a, b| a.last.cmp(&b.last));
    let n_samples = window[0].genotypes.len();
    
    already_logged.clear();
    already_logged.resize(window.len(), false);

    for sample in 0..n_samples {
        let sample_name = sample_names.get(sample).map(|s| s.as_str()).unwrap_or("UNKNOWN");
        let cols = count_collisions(window, sample);
        let hist_idx = std::cmp::min(cols, histogram.len() - 1);
        histogram[hist_idx] += 1;
        
        if cols > 0 {
            let removed = if method_2 { 
                independent_set_2(window, sample, weight_tag, weight_loc, default_w, sample_name, already_logged, use_gq, use_af, verbosity, sv_length) 
            } else { 
                independent_set_1(window, sample, weight_tag, weight_loc, default_w, sample_name, already_logged, use_gq, use_af, verbosity, sv_length) 
            };
            for i in 0..5 {
                removed_per_sample[sample][i] += removed[i];
            }
        }
    }

    window.sort_by(|a, b| a.input_index.cmp(&b.input_index));
    for iv in window {
        let _ = iv.write_vcf(output);
    }
}

fn main() -> std::io::Result<()> {
    let all_args: Vec<String> = env::args().collect();
    let mut use_gq = false;
    let mut use_af = false;
    let mut verbosity: u8 = 2;
    let mut sv_length: i32 = 50;
    let mut removed_counts_tsv_path: Option<String> = None;
    let mut histogram_tsv_path: Option<String> = None;
    
    let mut args = vec![all_args[0].clone()];
    let mut i = 1;
    while i < all_args.len() {
        match all_args[i].as_str() {
            "--use-gq" => use_gq = true,
            "--use-af" => use_af = true,
            "--verbosity" => {
                i += 1;
                verbosity = all_args.get(i).and_then(|s| s.parse().ok()).unwrap_or(2);
            }
            "--sv-length" => {
                i += 1;
                sv_length = all_args.get(i).and_then(|s| s.parse().ok()).unwrap_or(50);
            }
            "--removed-counts-tsv" => {
                i += 1;
                if let Some(path) = all_args.get(i) {
                    removed_counts_tsv_path = Some(path.clone());
                }
            }
            "--histogram-tsv" => {
                i += 1;
                if let Some(path) = all_args.get(i) {
                    histogram_tsv_path = Some(path.clone());
                }
            }
            _ => args.push(all_args[i].clone()),
        }
        i += 1;
    }

    if args.len() < 5 {
        eprintln!("Usage: {} [--use-gq] [--use-af] [--verbosity <1|2>] [--sv-length <int>] [--removed-counts-tsv <file>] [--histogram-tsv <file>] <method> <weight_tag> <weight_loc> <default_w> < input.vcf > output.vcf", args[0]);
        return Ok(());
    }

    let method_2 = args[1] == "1";
    let weight_tag = &args[2];
    let weight_loc = args[3] == "1";
    let default_weight = f64::from_str(&args[4]).unwrap_or(0.0);
    
    let stdin = io::stdin();
    let mut reader = BufReader::new(stdin.lock());
    
    let stdout = io::stdout();
    let mut out_vcf = BufWriter::new(stdout.lock());
    
    let mut histogram = vec![0; 100];
    let mut already_logged: Vec<bool> = Vec::new();
    let mut removed_per_sample: Vec<[usize; 5]> = Vec::new();

    let mut window: Vec<Interval> = Vec::new();
    let mut window_last_pos = -1;
    let mut line = String::new();
    let mut n_records = 0;
    
    let mut sample_names: Vec<String> = Vec::new();

    while reader.read_line(&mut line)? > 0 {
        let mut trimmed_len = line.len();
        while trimmed_len > 0 && (line.as_bytes()[trimmed_len - 1] == b'\n' || line.as_bytes()[trimmed_len - 1] == b'\r') {
            trimmed_len -= 1;
        }
        line.truncate(trimmed_len);
        
        if line.starts_with('#') {
            if line.starts_with("#CHROM") {
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() > 9 {
                    sample_names = fields[9..].iter().map(|s| s.to_string()).collect();
                    removed_per_sample = vec![[0; 5]; sample_names.len()];
                }
            }
            let _ = writeln!(out_vcf, "{}", line);
        } else {
            let iv = Interval::new(line.clone());
            window.push(iv);
            
            if !is_last_in_window(&window, window_last_pos) {
                let next_iv = window.pop().unwrap();
                for (i, iv) in window.iter_mut().enumerate() { iv.input_index = i; }
                
                let prev_records = n_records;
                n_records += window.len();
                
                process_window(&mut window, method_2, weight_tag, weight_loc, default_weight, &mut out_vcf, &mut histogram, &sample_names, &mut already_logged, &mut removed_per_sample, use_gq, use_af, verbosity, sv_length);
                
                if verbosity >= 2 && n_records / 10_000 > prev_records / 10_000 {
                    eprintln!("Processed {} records", n_records);
                }
                
                window.clear();
                window.push(next_iv);
                window_last_pos = window[0].last;
            } else {
                let last_iv = &window[window.len() - 1];
                if last_iv.last > window_last_pos { window_last_pos = last_iv.last; }
            }
        }
        line.clear();
    }
    
    if !window.is_empty() { 
        for (i, iv) in window.iter_mut().enumerate() { iv.input_index = i; }
        process_window(&mut window, method_2, weight_tag, weight_loc, default_weight, &mut out_vcf, &mut histogram, &sample_names, &mut already_logged, &mut removed_per_sample, use_gq, use_af, verbosity, sv_length); 
    }
    out_vcf.flush()?; 

    if let Some(hist_path) = histogram_tsv_path {
        let mut hist_w = BufWriter::new(File::create(hist_path)?);
        let _ = writeln!(hist_w, "NUM_COLLISIONS\tNUM_HAPLOTYPES");
        for (i, &count) in histogram.iter().enumerate() {
            let _ = writeln!(hist_w, "{}\t{}", i, count);
        }
    }

    let mut tsv_writer: Box<dyn Write> = if let Some(path) = removed_counts_tsv_path {
        Box::new(BufWriter::new(File::create(path)?))
    } else {
        Box::new(io::stderr())
    };

    let _ = writeln!(tsv_writer, "SAMPLE\tSV_DEL\tDEL\tSNP\tINS\tSV_INS\tTOTAL");
    for (i, counts) in removed_per_sample.iter().enumerate() {
        let total: usize = counts.iter().sum();
        let _ = writeln!(tsv_writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}", sample_names[i], counts[0], counts[1], counts[2], counts[3], counts[4], total);
    }

    Ok(())
}
