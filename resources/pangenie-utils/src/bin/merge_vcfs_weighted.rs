use clap::{Parser, Subcommand};
use indexmap::IndexMap;
use bio::io::fasta::IndexedReader;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader};

#[derive(Parser)]
#[command(name = "merge_vcfs")]
#[command(about = "Merge variant calls into pangenome multi-sample VCF", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Merge {
        #[arg(short = 'r')]
        reference: String,
        #[arg(short = 'h', long = "header")]
        header: String,
        #[arg(short = 'p', long = "ploidy")]
        ploidy: usize,
        #[arg(short = 'c', long = "chromosomes", default_value = "")]
        chromosomes: String,
        #[arg(long = "weight", help = "INFO field to use for resolving conflicts (e.g., SCORE, AF)")]
        weight: String,
        #[arg(long = "default-weight", default_value_t = 0.0, help = "Fallback weight if INFO field is missing")]
        default_weight: f64,
    },
}

#[derive(Clone, Debug)]
struct Allele {
    seq: String,
    ref_seq: String,
    weight: f64,
    start: usize,
    end: usize,
    id: String,
}

#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
struct Haplotype {
    alleles: Vec<i32>,
}

impl Haplotype {
    fn is_undefined(&self) -> bool {
        self.alleles.iter().any(|&a| a == -1)
    }
    fn is_reference(&self) -> bool {
        self.alleles.iter().all(|&a| a == 0)
    }
}

struct HaplotypeTable {
    samples: Vec<String>,
    ploidy: usize,
    weight_field: String,
    default_weight: f64,
    chrom: Option<String>,
    haplotypes: BTreeMap<usize, Vec<i32>>,
    alleles: BTreeMap<usize, Vec<Allele>>,
    start: usize,
    end: usize,
}

impl HaplotypeTable {
    fn new(samples: Vec<String>, ploidy: usize, weight_field: String, default_weight: f64) -> Self {
        Self {
            samples,
            ploidy,
            weight_field,
            default_weight,
            chrom: None,
            haplotypes: BTreeMap::new(),
            alleles: BTreeMap::new(),
            start: usize::MAX,
            end: 0,
        }
    }

    fn is_empty(&self) -> bool {
        self.alleles.is_empty()
    }

    fn parse(&mut self, row_index: usize, line: &str) -> usize {
        let fields: Vec<&str> = line.split_whitespace().collect();
        let chrom = fields[0].to_string();
        let pos = fields[1].parse::<usize>().unwrap();
        let ref_seq = fields[3];
        let alt_fields: Vec<&str> = fields[4].split(',').collect();

        let info = fields[7];
        let mut ids = Vec::new();
        let mut weights = vec![self.default_weight; alt_fields.len()];

        for item in info.split(';') {
            if let Some((k, v)) = item.split_once('=') {
                if k == "ID" {
                    ids = v.split(',').map(|s| s.to_string()).collect();
                } else if k == self.weight_field {
                    let vals: Vec<&str> = v.split(',').collect();
                    if vals.len() == alt_fields.len() {
                        for (i, val) in vals.iter().enumerate() {
                            weights[i] = val.parse::<f64>().unwrap_or(self.default_weight);
                        }
                    } else if vals.len() == 1 {
                        let w = vals[0].parse::<f64>().unwrap_or(self.default_weight);
                        for i in 0..alt_fields.len() {
                            weights[i] = w;
                        }
                    }
                }
            }
        }

        if ids.is_empty() {
            ids = fields[2].split(',').map(|s| s.to_string()).collect();
        }

        if alt_fields.len() != ids.len() {
            panic!("An ID needs to be provided for each individual allele.");
        }

        for (i, alt_seq) in alt_fields.iter().enumerate() {
            let allele = Allele {
                seq: alt_seq.to_string(),
                ref_seq: ref_seq.to_string(),
                weight: weights[i],
                start: pos,
                end: pos + ref_seq.len(),
                id: ids[i].clone(),
            };
            self.alleles.entry(row_index).or_default().push(allele);
            self.start = self.start.min(pos);
            self.end = self.end.max(pos + ref_seq.len());
            self.chrom = Some(chrom.clone());
        }

        let format_fields: Vec<&str> = fields[8].split(':').collect();
        let gt_idx = format_fields.iter().position(|&r| r == "GT").expect("Genotype information missing.");

        let mut col_idx = 0;
        for gt_str in &fields[9..] {
            let gt_val = gt_str.split(':').nth(gt_idx).unwrap();
            if gt_val.contains('/') { panic!("Line: {} contains an unphased genotype.", line); }

            for a in gt_val.split('|') {
                let val = if a == "." { -1 } else { a.parse::<i32>().unwrap() };
                self.haplotypes.entry(col_idx).or_default().push(val);
                col_idx += 1;
            }
        }
        ids.len()
    }

    fn merge(&self, reader: &mut IndexedReader<File>) -> (Option<String>, usize) {
        if self.is_empty() { return (None, 0); }

        let chrom = self.chrom.as_ref().unwrap();
        reader.fetch(chrom, (self.start - 1) as u64, (self.end - 1) as u64).unwrap();
        let mut ref_allele = Vec::new();
        reader.read(&mut ref_allele).unwrap();
        ref_allele.make_ascii_uppercase();
        let ref_allele_str = std::str::from_utf8(&ref_allele).unwrap();

        let mut hap_to_columns: IndexMap<Haplotype, Vec<usize>> = IndexMap::new();
        for (col, alleles) in &self.haplotypes {
            hap_to_columns.entry(Haplotype { alleles: alleles.clone() }).or_default().push(*col);
        }

        let mut hap_to_sequence: HashMap<Haplotype, (String, Vec<Allele>, f64, String)> = HashMap::new();
        let mut written_ids = HashSet::new();

        for (hap, columns) in &hap_to_columns {
            if hap.is_undefined() { continue; }
            
            // 1. Gather all active alleles on this haplotype
            let mut active_alleles: Vec<&Allele> = Vec::new();
            for (row_idx, &allele_idx) in hap.alleles.iter().enumerate() {
                if allele_idx > 0 {
                    active_alleles.push(&self.alleles[&row_idx][(allele_idx - 1) as usize]);
                }
            }

            // 2. Sort by weight descending to prioritize heavier variants during overlaps
            active_alleles.sort_by(|a, b| b.weight.partial_cmp(&a.weight).unwrap_or(std::cmp::Ordering::Equal));

            // 3. Greedily build a non-overlapping set of kept alleles
            let mut kept_alleles: Vec<&Allele> = Vec::new();
            for allele in active_alleles {
                let mut overlapped_with = None;
                for kept in &kept_alleles {
                    // Check for overlap [start, end)
                    if std::cmp::max(allele.start, kept.start) < std::cmp::min(allele.end, kept.end) {
                        overlapped_with = Some(*kept);
                        break;
                    }
                }

                if let Some(kept) = overlapped_with {
                    eprintln!("Two overlapping variants at same haplotype at {}:{}, resolved by weight.", chrom, self.start);
                    eprintln!("  Kept Allele    - CHROM: {} POS: {} REF: {} ALT: {} WEIGHT: {} ID: {}", chrom, kept.start, kept.ref_seq, kept.seq, kept.weight, kept.id);
                    eprintln!("  Skipped Allele - CHROM: {} POS: {} REF: {} ALT: {} WEIGHT: {} ID: {}", chrom, allele.start, allele.ref_seq, allele.seq, allele.weight, allele.id);
                    let sample_repr = format!("[{}]", columns.iter().map(|&h| format!("'{}'", self.samples[h / self.ploidy])).collect::<Vec<_>>().join(", "));
                    eprintln!("  SAMPLES: {}", sample_repr);
                } else {
                    kept_alleles.push(allele);
                }
            }

            // 4. Re-sort the surviving alleles by position to build the sequence
            kept_alleles.sort_by_key(|a| a.start);

            let mut sequence = ref_allele_str.to_string();
            let mut id_seq = Vec::new();
            let mut components = Vec::new();
            let mut offset: isize = self.start as isize;
            let mut total_weight = 0.0;

            for allele_obj in kept_alleles {
                let s = (allele_obj.start as isize - offset) as usize;
                let e = (allele_obj.end as isize - offset) as usize;
                
                let mut new_seq = sequence[..s].to_string();
                new_seq.push_str(&allele_obj.seq);
                new_seq.push_str(&sequence[e..]);
                
                offset -= (new_seq.len() as isize) - (sequence.len() as isize);
                sequence = new_seq;
                id_seq.push(allele_obj.id.clone());
                components.push(allele_obj.clone());
                total_weight += allele_obj.weight;
            }

            hap_to_sequence.insert(hap.clone(), (sequence, components, total_weight, id_seq.join(":")));
        }

        let mut alt_alleles = Vec::new();
        let mut ids_list = Vec::new();
        let mut genotypes: Vec<i32> = vec![-1; self.haplotypes.len()];
        let mut seq_to_kept_components: HashMap<String, (Vec<Allele>, f64, usize)> = HashMap::new();

        for (hap, columns) in &hap_to_columns {
            if let Some((seq, components, weight, id_str)) = hap_to_sequence.get(hap) {
                let allele_idx: i32;
                
                if hap.is_reference() {
                    allele_idx = 0;
                } else if let Some((kept_comps, kept_weight, pos)) = seq_to_kept_components.get_mut(seq) {
                    allele_idx = (*pos + 1) as i32;
                    
                    if weight > kept_weight {
                        eprintln!("Different allele combinations lead to same sequence at {}:{}.", chrom, self.start);
                        eprintln!("  Replacing kept combination with heavier alternative.");
                        eprintln!("  New Kept combination components (Total Weight {}):", weight);
                        for c in components {
                            eprintln!("    CHROM: {} POS: {} REF: {} ALT: {} WEIGHT: {} ID: {}", chrom, c.start, c.ref_seq, c.seq, c.weight, c.id);
                        }
                        eprintln!("  Skipped combination components (Total Weight {}):", kept_weight);
                        for c in kept_comps.iter() {
                            eprintln!("    CHROM: {} POS: {} REF: {} ALT: {} WEIGHT: {} ID: {}", chrom, c.start, c.ref_seq, c.seq, c.weight, c.id);
                        }
                        
                        *kept_comps = components.clone();
                        *kept_weight = *weight;
                        ids_list[*pos] = id_str.clone(); 
                        
                    } else if ids_list[*pos] != *id_str {
                        eprintln!("Different allele combinations lead to same sequence at {}:{}.", chrom, self.start);
                        eprintln!("  Keeping existing combination. Skipping lighter alternative.");
                        eprintln!("  Kept combination components (Total Weight {}):", kept_weight);
                        for c in kept_comps.iter() {
                            eprintln!("    CHROM: {} POS: {} REF: {} ALT: {} WEIGHT: {} ID: {}", chrom, c.start, c.ref_seq, c.seq, c.weight, c.id);
                        }
                        eprintln!("  Skipped combination components (Total Weight {}):", weight);
                        for c in components {
                            eprintln!("    CHROM: {} POS: {} REF: {} ALT: {} WEIGHT: {} ID: {}", chrom, c.start, c.ref_seq, c.seq, c.weight, c.id);
                        }
                    }
                } else {
                    alt_alleles.push(seq.clone());
                    let pos = alt_alleles.len() - 1;
                    seq_to_kept_components.insert(seq.clone(), (components.clone(), *weight, pos));
                    
                    ids_list.push(id_str.clone());
                    for id in id_str.split(':') { written_ids.insert(id.to_string()); }
                    allele_idx = alt_alleles.len() as i32;
                }
                
                for &col in columns { genotypes[col] = allele_idx; }
            } else {
                for &col in columns { genotypes[col] = -1; }
            }
        }

        if alt_alleles.is_empty() { return (None, 0); }

        let vcf_alt = alt_alleles.join(",");
        if ref_allele_str.chars().any(|c| !"CAGTcagt,".contains(c)) || vcf_alt.chars().any(|c| !"CAGTcagt,".contains(c)) {
            return (None, 0);
        }

        let vcf_genotypes: Vec<String> = genotypes.chunks(self.ploidy)
            .map(|chunk| {
                chunk.iter().map(|&g| if g == -1 { ".".to_string() } else { g.to_string() })
                    .collect::<Vec<_>>().join("|")
            }).collect();

        let vcf_line = format!(
            "{}\t{}\t.\t{}\t{}\t.\tPASS\tID={}\tGT\t{}",
            chrom, self.start, ref_allele_str, vcf_alt, ids_list.join(","), vcf_genotypes.join("\t")
        );

        (Some(vcf_line), written_ids.len())
    }
}

fn print_header(header_path: &str) {
    let file = File::open(header_path).expect("Could not open header file");
    let reader = BufReader::new(file);
    for line in reader.lines() {
        let l = line.unwrap();
        if l.starts_with("##") {
            println!("{}", l);
        } else if l.starts_with('#') {
            println!("##INFO=<ID=ID,Number=A,Type=String,Description=\"Variant IDs.\">");
            println!("{}", l);
        }
    }
}

fn main() {
    let cli = Cli::parse();
    match cli.command {
        Commands::Merge { reference, header, ploidy, chromosomes, weight, default_weight } => {
            let chrom_filter: HashSet<String> = chromosomes.split(',')
                .filter(|s| !s.is_empty()).map(|s| s.to_string()).collect();
            
            let mut ref_reader = IndexedReader::from_file(&reference).expect("Could not open reference");
            
            let mut total_input = 0;
            let mut total_written = 0;
            let mut samples = Vec::new();
            let mut table: Option<HaplotypeTable> = None;
            
            let mut prev_chrom = String::new();
            let mut prev_end = 0;
            let mut row_index = 0;

            let stdin = io::stdin();
            for line in stdin.lock().lines() {
                let l = line.unwrap();
                if l.starts_with("##") { continue; }
                if l.starts_with('#') {
                    if !samples.is_empty() {
                        panic!("Input is not a valid VCF file.");
                    }
                    let fields: Vec<&str> = l.split_whitespace().collect();
                    if fields.len() < 10 {
                        panic!("Input does not contain any samples.");
                    }
                    samples = fields[9..].iter().map(|s| s.to_string()).collect();
                    print_header(&header);
                    table = Some(HaplotypeTable::new(samples.clone(), ploidy, weight.clone(), default_weight));
                    continue;
                }

                let fields: Vec<&str> = l.split_whitespace().collect();
                let chrom = fields[0];
                if !chrom_filter.is_empty() && !chrom_filter.contains(chrom) { continue; }

                let start = fields[1].parse::<usize>().unwrap();
                let end = start + fields[3].len();

                if let Some(ref mut t) = table {
                    if ((start >= prev_end) || (prev_chrom != chrom)) && !t.is_empty() {
                        let (vcf, written) = t.merge(&mut ref_reader);
                        if let Some(line) = vcf { println!("{}", line); }
                        total_written += written;
                        *t = HaplotypeTable::new(samples.clone(), ploidy, weight.clone(), default_weight);
                        row_index = 0;
                    }
                    total_input += t.parse(row_index, &l);
                    row_index += 1;
                }

                prev_chrom = chrom.to_string();
                prev_end = if prev_chrom == chrom { prev_end.max(end) } else { end };
            }

            if let Some(t) = table {
                let (vcf, written) = t.merge(&mut ref_reader);
                if let Some(line) = vcf { println!("{}", line); }
                total_written += written;
            }

            eprintln!("Total number of input alleles: {}", total_input);
            eprintln!("Total number of written alleles: {}", total_written);
        }
    }
}
