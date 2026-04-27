use std::collections::BTreeMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VariantType { Del, Inv, Dup, Ins, Snp, Replacement }

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(usize)]
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

#[derive(Clone)]
struct Interval {
    chr: i32,
    first: i32,
    last: i32,
    input_index: usize,
    weight: f64,
    vcf_columns: Vec<String>,
    genotypes: Vec<Genotype>,
    
    in_independent_set: bool,
    independent_set_weight: f64,
    independent_set_previous: Option<usize>,
    overlaps_is_hap1: bool,
    overlaps_is_hap2: bool,
}

impl Interval {
    fn new(record: &str) -> Self {
        let columns: Vec<String> = record.split('\t').map(|s| s.to_string()).collect();
        let chr = parse_chr(&columns[0]);
        let pos = i32::from_str(&columns[1]).unwrap_or(0);
        let info = &columns[7];
        let ref_allele = &columns[3];
        let alt_allele = &columns[4];
        
        let variant_type = if let Some(svtype) = get_info_field(info, "SVTYPE") {
            parse_svtype(&svtype)
        } else if ref_allele.len() == 1 {
            if alt_allele.len() > 1 { VariantType::Ins } else { VariantType::Snp }
        } else if alt_allele.len() == 1 {
            VariantType::Del
        } else {
            VariantType::Replacement
        };
        
        let length = if let Some(svlen_str) = get_info_field(info, "SVLEN") {
            i32::from_str(&svlen_str).unwrap_or(0).abs()
        } else if variant_type == VariantType::Replacement {
            (ref_allele.len() - 1) as i32
        } else {
            (std::cmp::max(ref_allele.len(), alt_allele.len()) - 1) as i32
        };
        
        let (first, last) = match variant_type {
            VariantType::Del | VariantType::Inv | VariantType::Dup | VariantType::Replacement => (pos + 1, pos + length),
            VariantType::Ins => (pos, pos + 1),
            VariantType::Snp => (pos, pos),
        };
        
        let n_samples = if columns.len() > 9 { columns.len() - 9 } else { 0 };
        let mut genotypes = Vec::with_capacity(n_samples);
        for i in 0..n_samples { genotypes.push(Genotype::from_str(&columns[9 + i])); }
        
        Interval {
            chr, first, last, input_index: 0, weight: 0.0, vcf_columns: columns, genotypes,
            in_independent_set: false, independent_set_weight: 0.0,
            independent_set_previous: None, overlaps_is_hap1: false, overlaps_is_hap2: false,
        }
    }

    fn set_weight(&mut self, sample: usize, weight_tag: &str, in_sample: bool, default_weight: f64) {
        let mut value = None;
        if in_sample {
            let format_col = &self.vcf_columns[8];
            if let Some(p) = format_col.find(weight_tag) {
                let mut j = 0;
                for byte in format_col[..p].bytes() {
                    if byte == b':' { j += 1; }
                }
                
                let gt = &self.vcf_columns[9 + sample];
                let gt_bytes = gt.as_bytes();
                
                for i in 0..gt_bytes.len() {
                    if gt_bytes[i] != b':' { continue; }
                    j -= 1;
                    if j > 0 { continue; }
                    
                    if let Some(q) = gt[i + 1..].find(':') {
                        value = Some(gt[i + 1..i + 1 + q].to_string());
                    } else {
                        value = Some(gt[i + 1..].to_string());
                    }
                    break;
                }
            }
        } else {
            value = get_info_field(&self.vcf_columns[7], weight_tag);
        }
        self.weight = value.and_then(|s| f64::from_str(&s).ok()).unwrap_or(default_weight);
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
        N_GT_COLLISIONS[self.genotypes[sample] as usize][next.genotypes[sample] as usize] == 0 || self.last < next.first
    }

    fn clear_is_variables(&mut self) {
        self.in_independent_set = false;
        self.independent_set_weight = 0.0;
        self.independent_set_previous = None;
        self.overlaps_is_hap1 = false;
        self.overlaps_is_hap2 = false;
    }

    fn to_vcf_string(&self) -> String {
        let mut out = self.vcf_columns[0..9].join("\t");
        let format_str = &self.vcf_columns[8];
        
        let mut gt_index = 0;
        if let Some(p) = format_str.find("GT") {
            for byte in format_str[..p].bytes() {
                if byte == b':' { gt_index += 1; }
            }
        }

        for (i, gt_enum) in self.genotypes.iter().enumerate() {
            out.push('\t');
            let gt = &self.vcf_columns[9 + i];
            let gt_bytes = gt.as_bytes();
            let mut k = 0;
            let mut p2 = 0;
            
            for j in 0..gt_bytes.len() {
                if gt_bytes[j] == b':' {
                    k += 1;
                    if k == gt_index {
                        p2 = j + 1;
                        break;
                    }
                }
            }
            
            if p2 > 0 {
                out.push_str(&gt[0..p2]);
            }
            out.push_str(gt_enum.to_str());
            
            if let Some(q) = gt[p2..].find(':') {
                out.push_str(&gt[p2 + q..]);
            }
        }
        out
    }
}

fn get_info_field(info: &str, field: &str) -> Option<String> {
    let target = format!("{}=", field);
    info.find(&target).map(|start| {
        let val_start = start + target.len();
        let val_end = info[val_start..].find(';').map(|i| i + val_start).unwrap_or(info.len());
        info[val_start..val_end].to_string()
    })
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
            out += N_GT_COLLISIONS[window[j].genotypes[sample] as usize][window[i].genotypes[sample] as usize] as usize;
        }
    }
    out
}

fn mark_is_overlaps(window: &mut [Interval], sample_window: &[usize], sample: usize, hap1: bool, hap2: bool) {
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

            if in_is_i && !in_is_j {
                if hap1 && on_hap1_i && window[idx_j].on_hap1(sample) { window[idx_j].overlaps_is_hap1 = true; }
                if hap2 && on_hap2_i && window[idx_j].on_hap2(sample) { window[idx_j].overlaps_is_hap2 = true; }
            } else if in_is_j && !in_is_i {
                if hap1 && window[idx_j].on_hap1(sample) && on_hap1_i { window[idx_i].overlaps_is_hap1 = true; }
                if hap2 && window[idx_j].on_hap2(sample) && on_hap2_i { window[idx_i].overlaps_is_hap2 = true; }
            }
        }
    }
}

fn independent_set_1(window: &mut [Interval], sample: usize, weight_tag: &str, weight_loc: bool, default_w: f64) {
    let mut sample_window = Vec::new();
    for i in 0..window.len() {
        if window[i].is_present(sample) {
            window[i].clear_is_variables();
            window[i].set_weight(sample, weight_tag, weight_loc, default_w);
            sample_window.push(i);
        }
    }

    if sample_window.is_empty() { return; }

    let mut max_weight = 0.0;
    let mut active = vec![true; sample_window.len()];
    let mut best_is = Vec::new();

    for i in (0..sample_window.len()).rev() {
        let mut ancestors = Vec::new();
        is1_dfs(i, &mut active, &mut ancestors, 0.0, sample, &sample_window, window, &mut max_weight, &mut best_is);
    }

    for &idx in &best_is {
        window[sample_window[idx]].in_independent_set = true;
    }

    mark_is_overlaps(window, &sample_window, sample, true, true);

    for &idx in &sample_window {
        if window[idx].in_independent_set { continue; }
        
        let mut gt = window[idx].genotypes[sample];
        gt = if window[idx].overlaps_is_hap1 { REMOVE_HAP1_D[gt as usize] } else { REMOVE_HAP1_0[gt as usize] };
        gt = if window[idx].overlaps_is_hap2 { REMOVE_HAP2_D[gt as usize] } else { REMOVE_HAP2_0[gt as usize] };
        window[idx].genotypes[sample] = gt;
    }
}

fn is1_dfs(
    id: usize, active: &mut [bool], ancestors: &mut Vec<usize>, ancestors_weight: f64,
    sample: usize, sample_window: &[usize], window: &[Interval],
    max_weight: &mut f64, best_is: &mut Vec<usize>
) {
    let curr_idx = sample_window[id];
    let weight_prime = ancestors_weight + window[curr_idx].weight;
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
    
    if upper_bound < *max_weight { return; }
    
    if has_child {
        ancestors.push(id);
        for i in (0..id).rev() {
            if !active_prime[i] { continue; }
            if upper_bound < *max_weight { break; }
            is1_dfs(i, &mut active_prime, ancestors, weight_prime, sample, sample_window, window, max_weight, best_is);
            upper_bound -= window[sample_window[i]].weight;
        }
        ancestors.pop();
    } else if weight_prime > *max_weight {
        *max_weight = weight_prime;
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

fn independent_set_2(window: &mut [Interval], sample: usize, weight_tag: &str, weight_loc: bool, default_w: f64) {
    for hap in [1, 2] {
        let mut sample_window = Vec::new();
        let mut points_map: BTreeMap<i32, EndpointGroup> = BTreeMap::new();
        
        for i in 0..window.len() {
            let on_hap = if hap == 1 { window[i].on_hap1(sample) } else { window[i].on_hap2(sample) };
            if on_hap {
                window[i].clear_is_variables();
                window[i].set_weight(sample, weight_tag, weight_loc, default_w);
                sample_window.push(i);
                points_map.entry(window[i].first).or_default().opens.push(i);
                points_map.entry(window[i].last).or_default().closed.push(i);
            }
        }

        if sample_window.is_empty() { continue; }

        let mut available_weight = 0.0;
        let mut previous: Option<usize> = None;

        for (_pos, group) in points_map {
            for &idx in &group.opens {
                window[idx].independent_set_weight = available_weight + window[idx].weight;
                window[idx].independent_set_previous = previous;
            }
            for &idx in &group.closed {
                if window[idx].independent_set_weight > available_weight {
                    available_weight = window[idx].independent_set_weight;
                    previous = Some(idx);
                }
            }
        }

        let mut curr_trace = None;
        for i in (0..sample_window.len()).rev() {
            let idx = sample_window[i];
            if window[idx].independent_set_weight == available_weight {
                curr_trace = Some(idx);
                break;
            }
        }

        while let Some(idx) = curr_trace {
            window[idx].in_independent_set = true;
            curr_trace = window[idx].independent_set_previous;
        }

        if hap == 1 {
            mark_is_overlaps(window, &sample_window, sample, true, false);
            for &idx in &sample_window {
                if !window[idx].in_independent_set {
                    let gt = window[idx].genotypes[sample];
                    window[idx].genotypes[sample] = if window[idx].overlaps_is_hap1 { REMOVE_HAP1_D[gt as usize] } else { REMOVE_HAP1_0[gt as usize] };
                }
            }
        } else {
            mark_is_overlaps(window, &sample_window, sample, false, true);
            for &idx in &sample_window {
                if !window[idx].in_independent_set {
                    let gt = window[idx].genotypes[sample];
                    window[idx].genotypes[sample] = if window[idx].overlaps_is_hap2 { REMOVE_HAP2_D[gt as usize] } else { REMOVE_HAP2_0[gt as usize] };
                }
            }
        }
    }
}

fn process_window<W: Write>(window: &mut Vec<Interval>, method_2: bool, weight_tag: &str, weight_loc: bool, default_w: f64, output: &mut BufWriter<W>, histogram: &mut [usize]) {
    if window.is_empty() { return; }
    window.sort_by(|a, b| a.last.cmp(&b.last));
    let n_samples = window[0].genotypes.len();

    for sample in 0..n_samples {
        let cols = count_collisions(window, sample);
        let hist_idx = std::cmp::min(cols, histogram.len() - 1);
        histogram[hist_idx] += 1;
        
        if cols > 0 {
            if method_2 { independent_set_2(window, sample, weight_tag, weight_loc, default_w); }
            else { independent_set_1(window, sample, weight_tag, weight_loc, default_w); }
        }
    }

    window.sort_by(|a, b| a.input_index.cmp(&b.input_index));
    for iv in window {
        let _ = writeln!(output, "{}", iv.to_vcf_string());
    }
}

fn main() -> std::io::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() < 5 {
        eprintln!("Usage: ... <method> <weight_tag> <weight_loc> <default_w> [<out.hist>] < input.vcf > output.vcf");
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

    let mut window: Vec<Interval> = Vec::new();
    let mut window_last_pos = -1;
    let mut line = String::new();
    let mut n_records = 0;

    while reader.read_line(&mut line)? > 0 {
        let trimmed = line.trim_end();
        if trimmed.starts_with('#') {
            let _ = writeln!(out_vcf, "{}", trimmed);
        } else {
            let iv = Interval::new(trimmed);
            window.push(iv);
            
            if !is_last_in_window(&window, window_last_pos) {
                let next_iv = window.pop().unwrap();
                for (i, iv) in window.iter_mut().enumerate() { iv.input_index = i; }
                
                n_records += window.len();
                
                process_window(&mut window, method_2, weight_tag, weight_loc, default_weight, &mut out_vcf, &mut histogram);
                
                if n_records % 10_000 == 0 {
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
        process_window(&mut window, method_2, weight_tag, weight_loc, default_weight, &mut out_vcf, &mut histogram); 
    }
    out_vcf.flush()?; 

    if args.len() > 5 && args[5] != "null" {
        let mut hist_w = BufWriter::new(File::create(&args[5])?);
        let _ = writeln!(hist_w, "#nCollision \t nHaplotypes");
        for (i, &count) in histogram.iter().enumerate() {
            let _ = writeln!(hist_w, "{}\t{}", i, count);
        }
    }

    Ok(())
}
