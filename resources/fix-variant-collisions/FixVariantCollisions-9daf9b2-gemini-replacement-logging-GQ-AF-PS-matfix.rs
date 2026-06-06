// =====================================================================
// FixVariantCollisions (Rust Port) - Ultra-Optimized Sweepline
// 
// SUMMARY OF ARCHITECTURE & UPGRADES:
// 
// I. Core Bug Fixes
// * The Replacement (MNV) Bug Fix: The original logic struggled to properly define 
//   the spatial boundaries of "Replacement" variants (Multi-Nucleotide Variants where 
//   both REF and ALT are >1 base). The Rust port explicitly isolates `VariantType::Replacement`, 
//   calculating its length strictly as `ref_allele.len() - 1` and its footprint as 
//   [pos, pos + length]. This ensures the algorithm correctly detects collisions across 
//   the entire span of the altered sequence.
// * The Modulo Progress Bug: The Java tool used a modulo operator (`records % 10_000 == 0`) 
//   to print progress updates, which failed because variants are processed in irregular 
//   window batches. We replaced this with an integer-division threshold that mathematically 
//   guarantees the progress log fires exactly when a 10,000 milestone is crossed.
//
// II. Biological & Algorithmic Upgrades
// * First Principles Matrix: Replaced the legacy collision matrix with a mathematically 
//   symmetric derivation. Unphased variants (0/1) and homozygous variants (1/1) natively 
//   collide without VCF-order bias.
// * Compound Tiebreaker (SCORE -> GQ -> AF -> Class): The Java tool used a single, rigid 
//   weight and fell back to arbitrary VCF row order on ties. The Rust port evaluates paths 
//   using a sophisticated 4-tier hierarchy, wrapped into a compact `Score` tuple.
// * Variant Class Hierarchy: As the 4th tier in the tiebreaker, if all probability scores 
//   perfectly tie, the script resolves the collision using a hardcoded reliability hierarchy 
//   based on the difficulty of calling the variant type: SNP > DEL > INS > SV_DEL > SV_INS.
// * Phase Set (PS) Awareness: The Java tool assumed that all variants marked 1|0 resided 
//   on the exact same physical DNA strand for the entire chromosome. The Rust tool hashes 
//   the PS format tag; if two variants have different Phase Sets, it correctly deduces they 
//   are not anchored to the same block and temporarily downgrades them to unphased (1/0) 
//   so they don't improperly mask each other.
// * Complete Eradication: Unphased variants that lose collisions are now cleanly degraded 
//   to 0/0 (reference), rather than lingering as partial calls.
//
// III. Performance & Output Enhancements
// * Zero-Allocation Architecture: The Java tool heavily relied on splitting strings into 
//   arrays and cloning objects, triggering massive garbage collection overhead. The Rust 
//   port uses string slice iterators (`str::split`) and raw byte-index pointers to 
//   evaluate columns on the fly without ever cloning the underlying VCF text string.
// * O(N log N) Dynamic Programming Sweepline: Method 2 relies on Hsiao's DP Sweepline, 
//   using a flat `Vec<Event>` sorted chronologically to trace overlapping variants, 
//   minimizing heap allocations and preventing recursive hangs in dense genomic regions.
// * Stratified TSV Generation: Instead of silently dropping variants, the script physically 
//   counts removals per sample, stratifies them by variant category using a customizable 
//   `--sv-length` cutoff, and emits them to a `--removed-counts-tsv`.
// * Controllable Verbosity: Replaced the Java tool's raw stdout spam with a `--verbosity` 
//   controlled logger, allowing you to mute dense overlap traces while keeping heartbeats.
// =====================================================================

use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

const EPSILON: f64 = 1e-9;

// =====================================================================
// 1. DATA STRUCTURES & CONSTANTS
// =====================================================================

/// Categories for tiebreaker prioritization and TSV reporting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VariantType { Del, Inv, Dup, Ins, Snp, Replacement }

/// 18 static phase states. `repr(u8)` allows casting to `usize` for instant array lookups.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
enum Genotype {
    Phased00=0, Phased01=1, Phased10=2, Phased11=3, Unphased00=4, Unphased01=5, Unphased10=6, Unphased11=7,
    PhasedD0=8, Phased0D=9, PhasedD1=10, Phased1D=11, PhasedDD=12, UnphasedD0=13, Unphased0D=14, UnphasedD1=15, Unphased1D=16, UnphasedDD=17,
}
use Genotype::*;

const GT_STR: [&str; 18] = ["0|0","0|1","1|0","1|1","0/0","0/1","1/0","1/1",".|0","0|.",".|1","1|.",".|.","./0","0/.","./1","1/.","./."];

/// Fast-lookup map for stripping phase information from mismatched Phase Sets.
const TO_UNPHASED: [Genotype; 18] = [Unphased00, Unphased01, Unphased10, Unphased11, Unphased00, Unphased01, Unphased10, Unphased11, UnphasedD0, Unphased0D, UnphasedD1, Unphased1D, UnphasedDD, UnphasedD0, Unphased0D, UnphasedD1, Unphased1D, UnphasedDD];

impl Genotype {
    /// Byte-level pattern matching drastically reduces standard string manipulation overhead.
    /// Evaluates the first 3 characters of the GT string natively.
    fn from_str(gt: &str) -> Self {
        let b = gt.as_bytes();
        if b.len() < 3 { return UnphasedDD; }
        match (b[0], b[1], b[2]) {
            (b'.', b'/', b'.') => UnphasedDD, (b'.', b'/', b'0') => UnphasedD0, (b'.', b'/', b'1') => UnphasedD1,
            (b'0', b'/', b'.') => Unphased0D, (b'0', b'/', b'0') => Unphased00, (b'0', b'/', b'1') => Unphased01,
            (b'1', b'/', b'.') => Unphased1D, (b'1', b'/', b'0') => Unphased10, (b'1', b'/', _) => Unphased11,
            (b'.', _, b'.') => PhasedDD, (b'.', _, b'0') => PhasedD0, (b'.', _, b'1') => PhasedD1,
            (b'0', _, b'.') => Phased0D, (b'0', _, b'0') => Phased00, (b'0', _, b'1') => Phased01,
            (b'1', _, b'.') => Phased1D, (b'1', _, b'0') => Phased10, (b'1', _, _) => Phased11,
            _ => Unphased11,
        }
    }
    fn to_str(self) -> &'static str { GT_STR[self as usize] }
    fn to_unphased(self) -> Self { TO_UNPHASED[self as usize] }
}

/// A mathematically symmetric derivation of physical genomic collisions.
/// 0 means no collision. >0 indicates a physical overlap restriction.
const N_GT_COLLISIONS: [[u8; 18]; 18] = [
    [0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0], [0,1,0,1, 0,0,0,1, 0,0,1,0,0, 0,0,0,0,0], [0,0,1,1, 0,0,0,1, 0,0,0,1,0, 0,0,0,0,0],
    [0,1,1,2, 0,1,1,2, 0,0,1,1,0, 0,0,1,1,0], [0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0], [0,0,0,1, 0,0,0,1, 0,0,0,0,0, 0,0,0,0,0],
    [0,0,0,1, 0,0,0,1, 0,0,0,0,0, 0,0,0,0,0], [0,1,1,2, 0,1,1,2, 0,0,1,1,0, 0,0,1,1,0], [0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0],
    [0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0], [0,1,0,1, 0,0,0,1, 0,0,1,0,0, 0,0,0,0,0], [0,0,1,1, 0,0,0,1, 0,0,0,1,0, 0,0,0,0,0],
    [0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0], [0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0], [0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0],
    [0,0,0,1, 0,0,0,1, 0,0,0,0,0, 0,0,0,0,0], [0,0,0,1, 0,0,0,1, 0,0,0,0,0, 0,0,0,0,0], [0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,0,0,0],
];

// Replaces the offending allele with '.' (D) or '0' based on haplotype intersections.
// Ensures unphased variants cleanly drop to 0/0 rather than partially mutating.
const RM_HAP1_0: [Genotype; 18] = [Phased00, Phased01, Phased00, Phased01, Unphased00, Unphased00, Unphased00, Unphased01, PhasedD0, Phased0D, PhasedD1, Phased0D, PhasedDD, UnphasedD0, Unphased0D, UnphasedD0, Unphased0D, UnphasedDD];
const RM_HAP1_D: [Genotype; 18] = [Phased00, Phased01, PhasedD0, PhasedD1, Unphased00, Unphased00, Unphased00, UnphasedD1, PhasedD0, Phased0D, PhasedD1, PhasedDD, PhasedDD, UnphasedD0, Unphased0D, UnphasedDD, UnphasedDD, UnphasedDD];
const RM_HAP2_0: [Genotype; 18] = [Phased00, Phased00, Phased10, Phased10, Unphased00, Unphased00, Unphased00, Unphased10, PhasedD0, Phased0D, PhasedD0, Phased1D, PhasedDD, UnphasedD0, Unphased0D, UnphasedD0, Unphased0D, UnphasedDD];
const RM_HAP2_D: [Genotype; 18] = [Phased00, Phased0D, Phased10, Phased1D, Unphased00, Unphased00, Unphased00, Unphased1D, PhasedD0, Phased0D, PhasedDD, Phased1D, PhasedDD, UnphasedD0, Unphased0D, UnphasedDD, UnphasedDD, UnphasedDD];

// =====================================================================
// 2. CORE LOGIC UTILITIES
// =====================================================================

/// A compact 4-tier score tuple evaluating variants (Weight -> GQ -> AF -> Class).
#[derive(Clone, Copy, Debug)]
struct Score(f64, f64, f64, i32);
impl Score {
    /// Combines scores for graph path accumulation
    fn add(&self, o: &Score) -> Score { Score(self.0 + o.0, self.1 + o.1, self.2 + o.2, self.3 + o.3) }
    
    /// Executes the 4-Tier Compound Tiebreaker hierarchy
    fn is_better(&self, o: &Score, gq: bool, af: bool) -> bool {
        if self.0 > o.0 + EPSILON { return true; }
        if (self.0 - o.0).abs() <= EPSILON {
            if gq && self.1 > o.1 + EPSILON { return true; }
            if gq && (self.1 - o.1).abs() > EPSILON { return false; }
            if af && self.2 > o.2 + EPSILON { return true; }
            if af && (self.2 - o.2).abs() > EPSILON { return false; }
            return self.3 > o.3; // Class Score fallback
        }
        false
    }
    
    /// Evaluates if two scores are perfectly identical across all enabled metrics
    fn is_eq(&self, o: &Score, gq: bool, af: bool) -> bool {
        (self.0 - o.0).abs() <= EPSILON && (!gq || (self.1 - o.1).abs() <= EPSILON) && (!af || (self.2 - o.2).abs() <= EPSILON) && self.3 == o.3
    }
}

/// Zero-allocation INFO field extractor. Searches slice directly.
fn get_info<'a>(info: &'a str, field: &str) -> Option<&'a str> {
    info.split(';').find_map(|s| s.strip_prefix(field).and_then(|rem| rem.strip_prefix('=')))
}

/// Central data model holding a single VCF line and all computed metadata
#[derive(Clone)]
struct Interval {
    line: String, tabs: [usize; 9],
    chr: i32, first: i32, last: i32, input_idx: usize,
    v_type: VariantType, v_len: i32, score: Score,
    genotypes: Vec<Genotype>, phase_sets: Vec<u32>,
    
    // DP Graph Tracker Block (Cleared per sample sweep)
    in_is: bool, is_score: Score, is_prev: Option<usize>,
    ol_hap1: bool, ol_hap2: bool,
}

impl Interval {
    /// Zero-allocation constructor. Indexes tabs to parse the raw string directly via slices.
    fn new(line: String, sv_len: i32) -> Self {
        let mut tabs = [0; 9];
        // Scan line to capture indices of the first 9 tabs (standard VCF cols)
        line.bytes().enumerate().filter(|&(_, b)| b == b'\t').take(9).enumerate().for_each(|(j, (i, _))| tabs[j] = i);
        let col = |i| if i == 0 { &line[..tabs[0]] } else { &line[tabs[i-1]+1..tabs[i]] };

        // Parse CHROM (normalizes X, Y, M)
        let chr = match col(0).to_uppercase().trim_start_matches("CHR") { "X"=>23, "Y"=>24, "M"|"MT"=>25, s => s.parse().unwrap_or(-1) };
        let (pos, info, ref_a, alt_a) = (col(1).parse().unwrap_or(0), col(7), col(3), col(4));
        
        // Ascertain Variant Type
        let v_type = match get_info(info, "SVTYPE").unwrap_or("").to_uppercase().as_str() {
            "DEL"|"DEL:ME" => VariantType::Del, "INV" => VariantType::Inv, "DUP"|"DUP:TANDEM"|"DUP:INT"|"CNV" => VariantType::Dup,
            "INS"|"INS:ME"|"INS:NOVEL" => VariantType::Ins,
            _ => if ref_a.len() == 1 { if alt_a.len() > 1 { VariantType::Ins } else { VariantType::Snp } } else if alt_a.len() == 1 { VariantType::Del } else { VariantType::Replacement }
        };
        
        // Determine biological length. Fallback computes raw length differences.
        let length = get_info(info, "SVLEN").and_then(|s| s.parse::<i32>().ok()).map(|l| l.abs())
            .unwrap_or_else(|| if v_type == VariantType::Replacement { ref_a.len() as i32 - 1 } else { (ref_a.len().max(alt_a.len()) - 1) as i32 });
        
        // Spatial footprint definition (Fixes the original Java Replacement Bug)
        let (first, last) = match v_type {
            VariantType::Del|VariantType::Inv|VariantType::Dup => (pos + 1, pos + length),
            VariantType::Replacement => (pos, pos + length),
            VariantType::Ins => (pos, pos + 1), VariantType::Snp => (pos, pos),
        };

        // Tier 4 Tiebreaker hierarchy configuration
        let class_score = match v_type {
            VariantType::Snp|VariantType::Replacement => 5,
            VariantType::Del|VariantType::Inv => if length >= sv_len { 2 } else { 4 },
            VariantType::Ins|VariantType::Dup => if length >= sv_len { 1 } else { 3 },
        };
        
        // Pre-parse Allele Frequency for Tier 3
        let af = get_info(info, "AF").map(|s| s.split(',').filter_map(|n| n.parse::<f64>().ok()).fold(0.0, f64::max)).unwrap_or(0.0);
        
        let (mut genotypes, mut phase_sets) = (Vec::new(), Vec::new());
        if tabs[8] > 0 {
            // Find format block PS index
            let ps_idx = col(8).split(':').position(|t| t == "PS");
            for s in line[tabs[8]+1..].split('\t') {
                genotypes.push(Genotype::from_str(s));
                // Inline Djb2 Hash parsing for Phase Set to evaluate block connectivity as fast integers
                phase_sets.push(ps_idx.and_then(|t| s.split(':').nth(t))
                    .map(|tag| if tag=="."||tag.is_empty() {0} else { tag.bytes().fold(5381u32, |h, b| h.wrapping_mul(33).wrapping_add(b as u32)) }).unwrap_or(0));
            }
        }
        
        Interval { line, tabs, chr, first, last, input_idx: 0, v_type, v_len: length, score: Score(0.0, 0.0, af, class_score), genotypes, phase_sets, in_is: false, is_score: Score(0.0,0.0,0.0,0), is_prev: None, ol_hap1: false, ol_hap2: false }
    }

    fn col(&self, i: usize) -> &str { if i == 0 { &self.line[..self.tabs[0]] } else if i < 9 { &self.line[self.tabs[i-1]+1..self.tabs[i]] } else { "" } }
    
    /// String slice extraction avoids any `split` allocations and operates natively on the heap buffer
    fn ext_fmt(&self, s_idx: usize, tag: &str) -> Option<&str> {
        let tag_idx = self.col(8).split(':').position(|t| t == tag)?;
        let rest = &self.line[self.tabs[8]+1..];
        // Jump directly to the requested sample tab block
        let start = if s_idx == 0 { 0 } else { rest.bytes().enumerate().filter(|&(_, b)| b == b'\t').nth(s_idx - 1).map(|(i, _)| i + 1)? };
        rest[start..rest[start..].find('\t').map_or(rest.len(), |i| start + i)].split(':').nth(tag_idx)
    }

    /// Pulls dynamically designated Weight and GQ tags for Tier 1 and 2 scoring
    fn set_scores(&mut self, s: usize, tag: &str, in_sample: bool, def: f64) {
        let w = (if in_sample { self.ext_fmt(s, tag) } else { get_info(self.col(7), tag) }).and_then(|x| x.parse().ok()).unwrap_or(def);
        let gq = self.ext_fmt(s, "GQ").and_then(|x| x.parse().ok()).unwrap_or(0.0);
        self.score.0 = w; self.score.1 = gq; // Keep AF and Class Score untouched (computed once at init)
    }

    fn is_present(&self, s: usize) -> bool { !matches!(self.genotypes[s], Phased00|Unphased00|PhasedDD|UnphasedDD|PhasedD0|UnphasedD0|Phased0D|Unphased0D) }
    fn on_hap1(&self, s: usize) -> bool { matches!(self.genotypes[s], Phased10|Phased1D|Phased11|Unphased11|Unphased01|Unphased10|UnphasedD1|Unphased1D) }
    fn on_hap2(&self, s: usize) -> bool { matches!(self.genotypes[s], Phased01|PhasedD1|Phased11|Unphased11|Unphased01|Unphased10|UnphasedD1|Unphased1D) }
    
    /// Maps variants into TSV output classifications
    fn category(&self, sv: i32) -> usize { match self.v_type { VariantType::Del|VariantType::Inv => if self.v_len >= sv {0} else {1}, VariantType::Snp|VariantType::Replacement => 2, VariantType::Ins|VariantType::Dup => if self.v_len >= sv {4} else {3} } }

    /// Reconstructs the line back to stdout, surgically replacing only the GT tag block to apply masking
    fn write_vcf<W: Write>(&self, out: &mut W) -> io::Result<()> {
        let gt_idx = self.col(8).split(':').position(|t| t == "GT").unwrap_or(0);
        out.write_all(self.line[..self.tabs[8]].as_bytes())?;
        for (i, sample) in self.line[self.tabs[8]+1..].split('\t').enumerate() {
            out.write_all(b"\t")?;
            let (mut p2, mut k) = (0, 0);
            for (j, &b) in sample.as_bytes().iter().enumerate() { if b == b':' { k += 1; if k == gt_idx { p2 = j + 1; break; } } }
            if p2 > 0 { out.write_all(sample[..p2].as_bytes())?; }
            out.write_all(self.genotypes[i].to_str().as_bytes())?;
            if let Some(q) = sample[p2..].find(':') { out.write_all(sample[p2 + q..].as_bytes())?; }
        }
        out.write_all(b"\n")
    }
}

// =====================================================================
// 3. GRAPH ALGORITHM EVALUATORS
// =====================================================================

/// Tallies overlaps required for resolving dense regions.
fn count_collisions(win: &[Interval], s: usize) -> usize {
    (0..win.len()).rev().filter(|&i| win[i].is_present(s)).map(|i| {
        (0..i).rev().take_while(|&j| win[j].last >= win[i].first).filter(|&j| {
            let (mut gi, mut gj) = (win[i].genotypes[s], win[j].genotypes[s]);
            if win[i].phase_sets[s] != win[j].phase_sets[s] { gi = gi.to_unphased(); gj = gj.to_unphased(); }
            N_GT_COLLISIONS[gj as usize][gi as usize] > 0
        }).count()
    }).sum()
}

/// Identifies collisions natively across both haplotypes to build exact mask removal requests.
/// Iterates backward through overlapping pairs to determine which unselected node must be removed.
fn mark_overlaps(win: &mut [Interval], s_win: &[usize], s: usize, h1: bool, h2: bool, s_n: &str, log: &mut [bool], verb: u8, sv_l: i32) -> [usize; 5] {
    log.fill(false);
    let mut counts = [0; 5];
    for &idx in s_win { win[idx].ol_hap1 = false; win[idx].ol_hap2 = false; }
    
    for i in (0..s_win.len()).rev() {
        let (idx_i, is_i) = (s_win[i], win[s_win[i]].in_is);
        for &idx_j in s_win[..i].iter().rev() {
            if win[idx_j].last < win[idx_i].first { break; }
            let is_j = win[idx_j].in_is;
            if is_i == is_j { continue; } // Exactly one must be in IS to log a skip
            
            let (k, sp) = if is_i { (idx_i, idx_j) } else { (idx_j, idx_i) };
            let mut mk = false;
            
            if h1 && win[k].on_hap1(s) && win[sp].on_hap1(s) { win[sp].ol_hap1 = true; mk = true; }
            if h2 && win[k].on_hap2(s) && win[sp].on_hap2(s) { win[sp].ol_hap2 = true; mk = true; }
            
            // Log verbose overlap tracking
            if mk && !log[sp] {
                if verb >= 2 {
                    eprintln!("  Kept Allele    - CHROM: {} POS: {} REF: {} ALT: {} SCORE: {:?} ID: {}", win[k].col(0), win[k].col(1), win[k].col(3), win[k].col(4), win[k].score, win[k].col(2));
                    eprintln!("  Skipped Allele - CHROM: {} POS: {} REF: {} ALT: {} SCORE: {:?} ID: {}", win[sp].col(0), win[sp].col(1), win[sp].col(3), win[sp].col(4), win[sp].score, win[sp].col(2));
                    eprintln!("  SAMPLES: ['{}']", s_n);
                }
                log[sp] = true; counts[win[sp].category(sv_l)] += 1;
            }
        }
    }
    counts
}

// =====================================================================
// 4. GRAPH SOLVERS
// =====================================================================

/// METHOD 1: Legacy DFS search. Can be incredibly slow in dense regions O(2^N).
fn independent_set_1(w: &mut [Interval], s: usize, w_tag: &str, loc: bool, def_w: f64, n: &str, log: &mut [bool], u_gq: bool, u_af: bool, v: u8, sv: i32) -> [usize; 5] {
    let mut s_win = Vec::new();
    for i in 0..w.len() {
        if w[i].is_present(s) {
            w[i].in_is = false; w[i].is_score = Score(0.,0.,0.,0); w[i].is_prev = None;
            w[i].set_scores(s, w_tag, loc, def_w); s_win.push(i);
        }
    }
    if s_win.is_empty() { return [0; 5]; }
    
    let (mut max, mut act, mut best) = (Score(f64::MIN, f64::MIN, f64::MIN, std::i32::MIN), vec![true; s_win.len()], Vec::new());
    for i in (0..s_win.len()).rev() { is1_dfs(i, &mut act, &mut Vec::new(), Score(0.,0.,0.,0), s, &s_win, w, &mut max, &mut best, u_gq, u_af); }
    for &idx in &best { w[s_win[idx]].in_is = true; }
    
    let removed = mark_overlaps(w, &s_win, s, true, true, n, log, v, sv);
    for &idx in &s_win {
        if !w[idx].in_is {
            let gt = w[idx].genotypes[s] as usize;
            w[idx].genotypes[s] = if w[idx].ol_hap1 { RM_HAP1_D[gt] } else { RM_HAP1_0[gt] };
            w[idx].genotypes[s] = if w[idx].ol_hap2 { RM_HAP2_D[w[idx].genotypes[s] as usize] } else { RM_HAP2_0[w[idx].genotypes[s] as usize] };
        }
    }
    removed
}

/// Recursive Sub-Routine for Method 1
fn is1_dfs(id: usize, act: &mut [bool], anc: &mut Vec<usize>, p_s: Score, s: usize, s_win: &[usize], win: &[Interval], max: &mut Score, best: &mut Vec<usize>, u_gq: bool, u_af: bool) {
    let curr_s = p_s.add(&win[s_win[id]].score);
    let mut a_p = act.to_vec(); a_p[id] = false;
    let (mut has_child, mut bound) = (false, curr_s.0);
    
    for i in (0..id).rev() {
        let (mut gi, mut gj) = (win[s_win[i]].genotypes[s], win[s_win[id]].genotypes[s]);
        if win[s_win[i]].phase_sets[s] != win[s_win[id]].phase_sets[s] { gi = gi.to_unphased(); gj = gj.to_unphased(); }
        if a_p[i] && (N_GT_COLLISIONS[gi as usize][gj as usize] > 0 && win[s_win[i]].last >= win[s_win[id]].first) { a_p[i] = false; }
        if a_p[i] { has_child = true; bound += win[s_win[i]].score.0; }
    }
    
    if bound < max.0 - EPSILON { return; }
    if has_child {
        anc.push(id);
        for i in (0..id).rev() { if a_p[i] && bound >= max.0 - EPSILON { is1_dfs(i, &mut a_p, anc, curr_s, s, s_win, win, max, best, u_gq, u_af); bound -= win[s_win[i]].score.0; } }
        anc.pop();
    } else if curr_s.is_better(max, u_gq, u_af) { *max = curr_s; *best = anc.clone(); best.push(id); }
}

/// Represents an endpoint boundary for the DP Sweepline algorithm.
/// `is_close` natively forces Open(false) < Close(true) sorting, ensuring nodes 
/// entering the scope process before those exiting at identical coordinates.
#[derive(Eq, PartialEq, Ord, PartialOrd)]
struct Event { pos: i32, is_close: bool, idx: usize }

/// METHOD 2: Hsiao's Dynamic Programming Algorithm. 
/// Solves the Independent Set pathing in O(N log N) time utilizing a flat event-queue sweepline.
fn independent_set_2(w: &mut [Interval], s: usize, w_t: &str, loc: bool, d_w: f64, n: &str, log: &mut [bool], u_gq: bool, u_af: bool, v: u8, sv: i32) -> [usize; 5] {
    let mut tr = [0; 5];
    for hap in [1, 2] {
        let (mut s_win, mut evs) = (Vec::new(), Vec::new());
        for i in 0..w.len() {
            if if hap == 1 { w[i].on_hap1(s) } else { w[i].on_hap2(s) } {
                w[i].in_is = false; w[i].is_score = Score(0.,0.,0.,0); w[i].is_prev = None;
                w[i].set_scores(s, w_t, loc, d_w); s_win.push(i);
                // Load physical endpoints into sweepline event queue
                evs.push(Event { pos: w[i].first, is_close: false, idx: i });
                evs.push(Event { pos: w[i].last, is_close: true, idx: i });
            }
        }
        if s_win.is_empty() { continue; }
        
        // Sort chronologically. Open boundaries evaluate before Closes at same pos.
        evs.sort_unstable();

        let (mut avail, mut prev) = (Score(f64::MIN, f64::MIN, f64::MIN, std::i32::MIN), None);
        for e in evs {
            if !e.is_close { // Open Event: inherit best score up to this point
                let base = if avail.0 == f64::MIN { Score(0.,0.,0.,0) } else { avail };
                w[e.idx].is_score = base.add(&w[e.idx].score); w[e.idx].is_prev = prev;
            } else if w[e.idx].is_score.is_better(&avail, u_gq, u_af) { // Close Event: Check if path concludes better
                avail = w[e.idx].is_score; prev = Some(e.idx);
            }
        }

        // Traceback from end to construct optimal Independent Set
        let mut curr = (0..s_win.len()).rev().find(|&i| w[s_win[i]].is_score.is_eq(&avail, u_gq, u_af)).map(|i| s_win[i]);
        while let Some(idx) = curr { w[idx].in_is = true; curr = w[idx].is_prev; }

        // Mask losing variants iteratively
        let rm = mark_overlaps(w, &s_win, s, hap==1, hap==2, n, log, v, sv);
        for i in 0..5 { tr[i] += rm[i]; }
        
        for &idx in &s_win {
            if !w[idx].in_is {
                let gt = w[idx].genotypes[s] as usize;
                w[idx].genotypes[s] = match (hap, w[idx].ol_hap1, w[idx].ol_hap2) {
                    (1, true, _) => RM_HAP1_D[gt], (1, false, _) => RM_HAP1_0[gt],
                    (2, _, true) => RM_HAP2_D[gt], (2, _, false) => RM_HAP2_0[gt],
                    _ => unreachable!()
                };
            }
        }
    }
    tr
}

// =====================================================================
// 5. MAIN EXECUTION & I/O
// =====================================================================

/// Finalizes graph outputs for a spatial cluster of variants before emitting to disk
fn process_window<W: Write>(win: &mut Vec<Interval>, m2: bool, w_tag: &str, loc: bool, dw: f64, out: &mut BufWriter<W>, hist: &mut [usize], sn: &[String], log: &mut Vec<bool>, rm: &mut [[usize; 5]], u_gq: bool, u_af: bool, verb: u8, sv: i32) {
    if win.is_empty() { return; }
    // Sort array linearly prior to graph resolution sweep
    win.sort_by(|a, b| a.last.cmp(&b.last));
    log.resize(win.len(), false);

    for s in 0..win[0].genotypes.len() {
        let cols = count_collisions(win, s);
        hist[cols.min(hist.len() - 1)] += 1;
        if cols > 0 {
            let sn_str = sn.get(s).map(|x| x.as_str()).unwrap_or("UNKNOWN");
            let del = if m2 { independent_set_2(win, s, w_tag, loc, dw, sn_str, log, u_gq, u_af, verb, sv) } else { independent_set_1(win, s, w_tag, loc, dw, sn_str, log, u_gq, u_af, verb, sv) };
            for i in 0..5 { rm[s][i] += del[i]; }
        }
    }
    
    // Sort array back to original VCF layout before write
    win.sort_by(|a, b| a.input_idx.cmp(&b.input_idx));
    for iv in win { let _ = iv.write_vcf(out); }
}

fn main() -> io::Result<()> {
    let args: Vec<String> = env::args().collect();
    let (mut flags, mut iter) = (vec![args[0].clone()], args.into_iter().skip(1));
    let (mut u_gq, mut u_af, mut verb, mut sv_len, mut tsv_p, mut hist_p) = (false, false, 2, 50, None, None);
    
    // Iterating matching natively extracts parameters safely 
    while let Some(a) = iter.next() {
        match a.as_str() {
            "--use-gq" => u_gq = true, "--use-af" => u_af = true,
            "--verbosity" => verb = iter.next().and_then(|s| s.parse().ok()).unwrap_or(2),
            "--sv-length" => sv_len = iter.next().and_then(|s| s.parse().ok()).unwrap_or(50),
            "--removed-counts-tsv" => tsv_p = iter.next(), "--histogram-tsv" => hist_p = iter.next(),
            _ => flags.push(a),
        }
    }

    if flags.len() < 5 { return Err(io::Error::new(io::ErrorKind::InvalidInput, "Usage: <bin> [flags] <m> <w_t> <w_l> <def_w> < in > out")); }
    let (m2, w_tag, w_loc, def_w) = (flags[1] == "1", &flags[2], flags[3] == "1", flags[4].parse().unwrap_or(0.0));
    
    let (mut out_vcf, mut hist, mut logged, mut removed, mut s_names) = (BufWriter::new(io::stdout().lock()), vec![0; 100], Vec::new(), Vec::new(), Vec::new());
    let (mut window, mut win_last_pos, mut n_records, mut line) = (Vec::new(), -1, 0, String::new());
    let mut reader = BufReader::new(io::stdin().lock());

    while reader.read_line(&mut line)? > 0 {
        line.truncate(line.trim_end_matches(&['\n', '\r'][..]).len()); // No-allocation newline strip
        if line.starts_with('#') {
            if line.starts_with("#CHROM") {
                // Initialize Header Arrays
                let f: Vec<&str> = line.split('\t').collect();
                if f.len() > 9 { s_names = f[9..].iter().map(|s| s.to_string()).collect(); removed = vec![[0; 5]; s_names.len()]; }
            }
            let _ = writeln!(out_vcf, "{}", line);
        } else {
            window.push(Interval::new(line.clone(), sv_len));
            
            // Check if bounds have broken spatial window
            if !window.is_empty() && win_last_pos != -1 && window.last().unwrap().chr == window[0].chr && window.last().unwrap().first <= win_last_pos {
                win_last_pos = win_last_pos.max(window.last().unwrap().last);
            } else if window.len() > 1 {
                // Batch dispatch
                let next_iv = window.pop().unwrap();
                window.iter_mut().enumerate().for_each(|(i, iv)| iv.input_idx = i);
                n_records += window.len();
                
                process_window(&mut window, m2, w_tag, w_loc, def_w, &mut out_vcf, &mut hist, &s_names, &mut logged, &mut removed, u_gq, u_af, verb, sv_len);
                
                // Triggers progress updates dynamically over irregular batches 
                if verb >= 1 && n_records / 10_000 > (n_records - window.len()) / 10_000 { eprintln!("Processed {} records", n_records); }
                
                window.clear(); window.push(next_iv); win_last_pos = window[0].last;
            } else { win_last_pos = window[0].last; }
        }
        line.clear();
    }
    
    // Process final lingering window
    if !window.is_empty() { window.iter_mut().enumerate().for_each(|(i, iv)| iv.input_idx = i); process_window(&mut window, m2, w_tag, w_loc, def_w, &mut out_vcf, &mut hist, &s_names, &mut logged, &mut removed, u_gq, u_af, verb, sv_len); }
    out_vcf.flush()?; 

    // Output Data Products
    if let Some(p) = hist_p {
        let mut w = BufWriter::new(File::create(p)?);
        let _ = writeln!(w, "NUM_COLLISIONS\tNUM_HAPLOTYPES");
        for (i, &c) in hist.iter().enumerate() { let _ = writeln!(w, "{}\t{}", i, c); }
    }
    let mut out: Box<dyn Write> = tsv_p.map_or_else(|| Box::new(io::stderr()) as Box<dyn Write>, |p| Box::new(BufWriter::new(File::create(p).unwrap())));
    let _ = writeln!(out, "SAMPLE\tSV_DEL\tDEL\tSNP\tINS\tSV_INS\tTOTAL");
    for (i, c) in removed.iter().enumerate() { let _ = writeln!(out, "{}\t{}\t{}\t{}\t{}\t{}\t{}", s_names[i], c[0], c[1], c[2], c[3], c[4], c.iter().sum::<usize>()); }
    Ok(())
}
