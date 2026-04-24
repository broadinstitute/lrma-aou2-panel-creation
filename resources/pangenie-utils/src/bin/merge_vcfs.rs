use clap::{Parser, Subcommand};
use indexmap::IndexMap;
use rust_htslib::faidx;
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
    },
}

#[derive(Clone, Debug)]
struct Allele {
    seq: String,
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
    chrom: Option<String>,
    // BTreeMap guarantees iteration in 0, 1, 2... column/row order (matching Python Dict insertion order)
    haplotypes: BTreeMap<usize, Vec<i32>>,
    alleles: BTreeMap<usize, Vec<Allele>>,
    start: usize,
    end: usize,
}

impl HaplotypeTable {
    fn new(samples: Vec<String>, ploidy: usize) -> Self {
        Self {
            samples,
            ploidy,
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
        for item in info.split(';') {
            if item.starts_with("ID=") {
                ids = item[3..].split(',').map(|s| s.to_string()).collect();
            }
        }

        if ids.is_empty() {
            ids = fields[2].split(',').map(|s| s.to_string()).collect();
        }

        if alt_fields.len() != ids.len() {
            panic!("An ID needs to be provided for each individual ID, either in ID colum or info column.");
        }

        for (i, alt_seq) in alt_fields.iter().enumerate() {
            let allele = Allele {
                seq: alt_seq.to_string(),
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
        let gt_idx = format_fields.iter().position(|&r| r == "GT").expect("Genotype information missing from VCF file.");

        let mut col_idx = 0;
        for gt_str in &fields[9..] {
            let gt_val = gt_str.split(':').nth(gt_idx).unwrap();
            
            if gt_val.contains('/') {
                panic!("Line: {} contains an unphased genotype.", line);
            }

            let sub_alleles: Vec<&str> = gt_val.split('|').collect();

            for a in sub_alleles {
                let val = if a == "." { -1 } else { a.parse::<i32>().unwrap() };
                self.haplotypes.entry(col_idx).or_default().push(val);
                col_idx += 1;
            }
        }
        ids.len()
    }

    fn merge(&self, reader: &faidx::Reader) -> (Option<String>, usize) {
        if self.is_empty() { return (None, 0); }

        let chrom = self.chrom.as_ref().unwrap();
        // FIX: end is inclusive in rust-htslib! To match Python's [start-1:end-1] exclusive slice, we must subtract 2.
        let ref_allele = reader.fetch_seq(chrom, self.start - 1, self.end - 2).unwrap().to_ascii_uppercase();
        let ref_allele_str = std::str::from_utf8(&ref_allele).unwrap();

        let mut hap_to_columns: IndexMap<Haplotype, Vec<usize>> = IndexMap::new();
        for (col, alleles) in &self.haplotypes {
            hap_to_columns.entry(Haplotype { alleles: alleles.clone() }).or_default().push(*col);
        }

        let mut hap_to_sequence = HashMap::new();
        let mut hap_to_id = HashMap::new();
        let mut written_ids = HashSet::new();

        for (hap, columns) in &hap_to_columns {
            if hap.is_undefined() { continue; }
            
            let mut sequence = ref_allele_str.to_string();
            let mut id_seq = Vec::new();
            let mut offset: isize = self.start as isize;
            let mut prev_end = 0;
            let mut skip_hap = false;

            for (row_idx, &allele_idx) in hap.alleles.iter().enumerate() {
                if allele_idx > 0 {
                    let allele_obj = &self.alleles[&row_idx][(allele_idx - 1) as usize];
                    if prev_end > allele_obj.start {
                        eprintln!("Two overlapping variants at same haplotype at {}:{}, set allele to missing.", chrom, self.start);
                        eprintln!("ALT ID: {}", allele_obj.id);
                        let sample_repr = format!("[{}]", columns.iter().map(|&h| format!("'{}'", self.samples[h / self.ploidy])).collect::<Vec<_>>().join(", "));
                        eprintln!("SAMPLES: {}", sample_repr);
                        
                        skip_hap = true;
                        break;
                    }
                    
                    let s = (allele_obj.start as isize - offset) as usize;
                    let e = (allele_obj.end as isize - offset) as usize;
                    let mut new_seq = sequence[..s].to_string();
                    new_seq.push_str(&allele_obj.seq);
                    new_seq.push_str(&sequence[e..]);
                    
                    offset -= (new_seq.len() as isize) - (sequence.len() as isize);
                    sequence = new_seq;
                    id_seq.push(allele_obj.id.clone());
                    prev_end = allele_obj.end;
                }
            }

            if !skip_hap {
                hap_to_sequence.insert(hap.clone(), sequence);
                hap_to_id.insert(hap.clone(), id_seq.join(":"));
            }
        }

        let mut alt_alleles = Vec::new();
        let mut ids_list = Vec::new();
        let mut genotypes: Vec<i32> = vec![-1; self.haplotypes.len()];

        for (hap, columns) in &hap_to_columns {
            if let Some(seq) = hap_to_sequence.get(hap) {
                let allele_idx: i32;
                if hap.is_reference() {
                    allele_idx = 0;
                } else if let Some(pos) = alt_alleles.iter().position(|s| s == seq) {
                    allele_idx = (pos + 1) as i32;
                    let current_id = hap_to_id.get(hap).unwrap();
                    if !ids_list.contains(current_id) {
                        eprintln!("Different allele combinations lead to same sequence at {}:{}.", chrom, self.start);
                    }
                } else {
                    alt_alleles.push(seq.clone());
                    let id_str = hap_to_id.get(hap).unwrap();
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
        Commands::Merge { reference, header, ploidy, chromosomes } => {
            let chrom_filter: HashSet<String> = chromosomes.split(',')
                .filter(|s| !s.is_empty()).map(|s| s.to_string()).collect();
            
            let ref_reader = faidx::Reader::from_path(&reference).expect("Could not open reference");
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
                    table = Some(HaplotypeTable::new(samples.clone(), ploidy));
                    continue;
                }

                let fields: Vec<&str> = l.split_whitespace().collect();
                let chrom = fields[0];
                if !chrom_filter.is_empty() && !chrom_filter.contains(chrom) { continue; }

                let start = fields[1].parse::<usize>().unwrap();
                let end = start + fields[3].len();

                if let Some(ref mut t) = table {
                    if ((start >= prev_end) || (prev_chrom != chrom)) && !t.is_empty() {
                        let (vcf, written) = t.merge(&ref_reader);
                        if let Some(line) = vcf { println!("{}", line); }
                        total_written += written;
                        *t = HaplotypeTable::new(samples.clone(), ploidy);
                        row_index = 0;
                    }
                    total_input += t.parse(row_index, &l);
                    row_index += 1;
                }

                prev_chrom = chrom.to_string();
                prev_end = if prev_chrom == chrom { prev_end.max(end) } else { end };
            }

            if let Some(t) = table {
                let (vcf, written) = t.merge(&ref_reader);
                if let Some(line) = vcf { println!("{}", line); }
                total_written += written;
            }

            eprintln!("Total number of input alleles: {}", total_input);
            eprintln!("Total number of written alleles: {}", total_written);
        }
    }
}
