use clap::Parser;
use std::collections::HashMap;
use std::io::{self, BufRead, Write};
use std::process;

#[derive(Parser, Debug)]
#[command(name = "prepare_vcf_and_add_ids")]
#[command(about = "Filter missing/N alleles, add IDs, and strip non-GT formats", long_about = None)]
struct Args {
    #[arg(long, default_value_t = 0.0, help = "Maximum allowed fraction of missing alleles per position.")]
    missing: f64,
}

fn main() {
    let args = Args::parse();
    
    let stdin = io::stdin();
    let mut reader = stdin.lock();
    let stdout = io::stdout();
    let mut stdout_lock = stdout.lock();
    let mut stderr = io::stderr();

    let mut total_records = 0;
    let mut total_alleles = 0;
    let mut missing_records = 0;
    let mut missing_alleles = 0;
    let mut ns_records = 0;
    let mut ns_alleles = 0;
    let mut written_records = 0;
    let mut written_alleles = 0;

    let mut allele_index = 0;
    let mut buffer = String::new();

    while reader.read_line(&mut buffer).unwrap() > 0 {
        let line = buffer.trim_end();
        if line.is_empty() {
            buffer.clear();
            continue;
        }

        if line.starts_with("##") {
            writeln!(stdout_lock, "{}", line).unwrap();
            buffer.clear();
            continue;
        }
        
        if line.starts_with('#') {
            writeln!(stdout_lock, "##INFO=<ID=ID,Number=A,Type=String,Description=\"Variant IDs.\">").unwrap();
            writeln!(stdout_lock, "{}", line).unwrap();
            buffer.clear();
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 {
            buffer.clear();
            continue;
        }

        total_records += 1;
        let chrom = fields[0];
        let position = fields[1];
        let alt_field = fields[4];
        let alts: Vec<&str> = alt_field.split(',').collect();
        let n_alt_alleles = alts.len();
        total_alleles += n_alt_alleles;

        let has_ns = alt_field.chars().any(|c| !"CAGTcagt,".contains(c));
        if has_ns {
            ns_records += 1;
            ns_alleles += n_alt_alleles;
            buffer.clear();
            continue;
        }

        // --- Calculate Missingness using EXACT Python mechanics ---
        // Python performs this logic blindly on the raw string before any formatting is stripped.
        let mut n_missing = 0.0;
        let mut n_total = 0.0;
        let mut unphased_prepare_error = false;

        for gt_str in fields.iter().skip(9) {
            n_total += 2.0;
            if *gt_str == "./." || *gt_str == "." {
                n_missing += 2.0;
            } else if gt_str.contains('|') {
                let alleles: Vec<&str> = gt_str.split('|').collect();
                n_missing += alleles.iter().filter(|&&a| a == ".").count() as f64;
            } else {
                unphased_prepare_error = true;
            }
        }

        if unphased_prepare_error {
            eprintln!("VCF contains unphased positions.");
            process::exit(1);
        }

        let frac_missing = n_missing / n_total;
        if frac_missing > args.missing {
            missing_records += 1;
            missing_alleles += n_alt_alleles;
            buffer.clear();
            continue;
        }

        // --- Format GT Fields (add-ids.py logic) ---
        let format_field = fields[8];
        let format_parts: Vec<&str> = format_field.split(':').collect();
        let gt_index = format_parts.iter().position(|&r| r == "GT").unwrap_or_else(|| {
            eprintln!("Input VCF does not contain GT field at position {}:{}.", chrom, position);
            process::exit(1);
        });

        let mut formatted_gts = Vec::with_capacity(fields.len() - 9);
        for gt_str in fields.iter().skip(9) {
            let gt_parts: Vec<&str> = gt_str.split(':').collect();
            let gt = gt_parts[gt_index];
            if !gt.contains('|') {
                eprintln!("Input VCF contains unphased genotype at position {}:{}.", chrom, position);
                process::exit(1);
            }
            formatted_gts.push(gt);
        }

        // --- Replicate Python Dict Update Behavior (preserves order) ---
        let info_string = fields[7];
        let mut info_keys = Vec::new();
        let mut info_map = HashMap::new();
        
        for pair in info_string.split(';') {
            if let Some((k, v)) = pair.split_once('=') {
                if !info_map.contains_key(k) {
                    info_keys.push(k);
                }
                info_map.insert(k, v);
            }
        }

        let mut ids = Vec::with_capacity(n_alt_alleles);
        for a in &alts {
            let var_len = a.len();
            let var_id = format!("{}-{}-allele{}-{}", chrom, position, allele_index, var_len);
            ids.push(var_id);
            allele_index += 1;
        }
        
        let ids_joined = ids.join(",");
        if !info_map.contains_key("ID") {
            info_keys.push("ID");
        }
        info_map.insert("ID", ids_joined.as_str());
        
        let updated_info = info_keys.iter()
            .map(|k| format!("{}={}", k, info_map[k]))
            .collect::<Vec<String>>()
            .join(";");

        write!(stdout_lock, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tGT", 
            chrom, position, fields[2], fields[3], alt_field, fields[5], fields[6], updated_info).unwrap();
        
        for gt in formatted_gts {
            write!(stdout_lock, "\t{}", gt).unwrap();
        }
        writeln!(stdout_lock).unwrap();

        written_records += 1;
        written_alleles += n_alt_alleles;
        buffer.clear();
    }

    write!(stderr, "skipped {} ({}) records (alleles) for which fraction of missing alleles exceeds threshold.\n", missing_records, missing_alleles).unwrap();
    write!(stderr, "skipped {} ({}) records (alleles) for which alternative alleles contained Ns\n.", ns_records, ns_alleles).unwrap();
    write!(stderr, "kept {} ({}) records (alleles) of {} ({}).\n", written_records, written_alleles, total_records, total_alleles).unwrap();
}
