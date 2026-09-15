use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::time::Instant;

fn smart_open(filename: &str) -> Box<dyn BufRead> {
    let file = File::open(filename).expect("Cannot open file");
    if filename.ends_with(".gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    }
}

// Emulates the exact f32 formatting of main_2.rs before Python sees it
fn format_float_3(val: f32) -> String {
    let s = format!("{:.3}", val);
    let trimmed = s.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
}

// Emulates Python numpy's vectorized_sig_fig_round on the final f64 values (Banker's Rounding)
fn round_to_3_sig_figs(x: f64) -> f64 {
    if x == 0.0 || !x.is_finite() {
        return x;
    }
    let power = 2.0 - x.abs().log10().floor();
    let factor = 10_f64.powf(power);
    
    let scaled = x * factor;
    let r = scaled.round();
    
    // Emulate numpy's round-half-to-even behavior
    let rounded = if (scaled - r).abs() == 0.5 {
        if r % 2.0 == 0.0 { r } else { r - 1.0 * scaled.signum() }
    } else {
        r
    };
    
    rounded / factor
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 4 {
        eprintln!("Usage: zcat <orig_vcf> | {} <panel_sites_vcf> <num_samples_csv> <out_batched_tsv_gz> [max_lines]", args[0]);
        std::process::exit(1);
    }

    let panel_sites_path = &args[1];
    let num_samples: Vec<usize> = args[2].split(',').map(|s| s.parse().unwrap()).collect();
    let total_samples: usize = num_samples.iter().sum();
    let batched_tsv_path = &args[3];
    
    let max_records: Option<usize> = if args.len() > 4 {
        Some(args[4].parse().expect("max_lines must be a valid integer"))
    } else {
        None
    };

    let mut tsv_out = GzEncoder::new(File::create(batched_tsv_path).unwrap(), Compression::default());
    
    let mut tsv_header = String::from("CHROM\tPOS\tREF\tALT");
    for i in 0..num_samples.len() {
        tsv_header.push_str(&format!("\tAF_{}\tINFO_{}", i, i));
    }
    writeln!(tsv_out, "{}", tsv_header).unwrap();

    eprintln!("Loading panel sites...");
    let mut panel_dict: HashMap<String, (String, String, String)> = HashMap::new();
    
    for line in smart_open(panel_sites_path).lines().filter_map(Result::ok) {
        if line.starts_with('#') { continue; }
        let fields: Vec<&str> = line.split('\t').collect();
        let key = format!("{}:{}:{}:{}", fields[0], fields[1], fields[3], fields[4]);
        
        let orig_id = fields[2].to_string();
        let mut atomic_id = String::from(".");
        let mut raf_val = String::new();
        let mut ac_val = String::new();
        let mut an_val = String::new();
        
        for item in fields[7].split(';') {
            if let Some(v) = item.strip_prefix("ID=") { atomic_id = v.to_string(); }
            else if let Some(v) = item.strip_prefix("AF=") { raf_val = v.to_string(); }
            else if let Some(v) = item.strip_prefix("AC=") { ac_val = v.to_string(); }
            else if let Some(v) = item.strip_prefix("AN=") { an_val = v.to_string(); }
        }
        
        if raf_val.is_empty() && !ac_val.is_empty() && !an_val.is_empty() {
            if let (Ok(ac), Ok(an)) = (ac_val.parse::<f64>(), an_val.parse::<f64>()) {
                if an > 0.0 { raf_val = format_float_3((ac / an) as f32); }
            }
        }
        panel_dict.insert(key, (atomic_id, orig_id, raf_val));
    }

    let stdout = io::stdout();
    let mut vcf_out = BufWriter::new(stdout.lock());
    let stdin = io::stdin();
    let mut stdin_lock = stdin.lock();
    
    let mut line = String::new();
    let mut records_processed: usize = 0;
    let start_time = Instant::now();

    eprintln!("Processing original VCF...");
    
    // Read directly into the reused 'line' buffer to avoid allocation overhead
    while stdin_lock.read_line(&mut line).unwrap() > 0 {
        // Strip the trailing newline safely
        if line.ends_with('\n') {
            line.pop();
            if line.ends_with('\r') {
                line.pop();
            }
        }

        if line.starts_with("#CHROM") {
            writeln!(vcf_out, "##INFO=<ID=ID,Number=1,Type=String,Description=\"Atomic variant ID\">").unwrap();
            writeln!(vcf_out, "##INFO=<ID=RAF,Number=A,Type=Float,Description=\"Panel reference allele frequency\">").unwrap();
            writeln!(vcf_out, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Recalculated allele frequency\">").unwrap();
            writeln!(vcf_out, "##INFO=<ID=INFO,Number=A,Type=Float,Description=\"Recalculated IMPUTE INFO score\">").unwrap();
            writeln!(vcf_out, "{}", line).unwrap();
            line.clear();
            continue;
        } else if line.starts_with('#') {
            if !line.starts_with("##INFO=<ID=RAF,") && 
               !line.starts_with("##INFO=<ID=AF,") && 
               !line.starts_with("##INFO=<ID=INFO,") &&
               !line.starts_with("##INFO=<ID=ID,") {
                writeln!(vcf_out, "{}", line).unwrap();
            }
            line.clear();
            continue;
        }

        let mut fields: Vec<&str> = line.split('\t').collect();
        let key = format!("{}:{}:{}:{}", fields[0], fields[1], fields[3], fields[4]);
        
        records_processed += 1;
        if records_processed % 10_000 == 0 {
            let elapsed = start_time.elapsed().as_secs();
            let hours = elapsed / 3600;
            let mins = (elapsed % 3600) / 60;
            let secs = elapsed % 60;
            eprintln!(
                "[{:02}:{:02}:{:02}] Processed {} input records... (Currently at {}:{})",
                hours, mins, secs, records_processed, fields[0], fields[1]
            );
        }
        
        let (atomic_id, orig_id, true_raf) = panel_dict.get(&key).cloned().unwrap_or((".".to_string(), ".".to_string(), "".to_string()));
        
        let orig_id_str = orig_id.clone();
        fields[2] = &orig_id_str;

        let fmt: Vec<&str> = fields[8].split(':').collect();
        let gp_idx = fmt.iter().position(|&x| x == "GP").expect("No GP format found");

        let mut batch_af_strs = Vec::with_capacity(num_samples.len());
        let mut batch_info_strs = Vec::with_capacity(num_samples.len());
        let mut sample_offset = 9;

        // 1. Calculate per-batch statistics exactly as f32 mirroring main_2.rs
        for &batch_size in &num_samples {
            let mut ds_sum = 0.0_f32;
            let mut ds2_sum = 0.0_f32;
            let mut ds4_sum = 0.0_f32;

            for s in 0..batch_size {
                let sample_data = fields[sample_offset + s];
                if sample_data != "." {
                    let parts: Vec<&str> = sample_data.split(':').collect();
                    if let Some(gp_str) = parts.get(gp_idx) {
                        let gps: Vec<f32> = gp_str.split(',').filter_map(|v| v.parse().ok()).collect();
                        if gps.len() == 3 {
                            let gp1_q = gps[1];
                            let gp2_q = gps[2];
                            let ds_q = gp1_q + 2.0 * gp2_q;
                            
                            ds_sum += ds_q;
                            ds2_sum += ds_q * ds_q;
                            ds4_sum += gp1_q + 4.0 * gp2_q;
                        }
                    }
                }
            }
            sample_offset += batch_size;

            let n_tar_haps = 2.0 * batch_size as f32;
            let safe_n = n_tar_haps.max(1e-9);
            let af = ds_sum / safe_n;
            let denom = n_tar_haps * af * (1.0 - af);
            
            let mut recalc_info = 1.0_f32;
            if af > 0.0 && af < 1.0 && denom > 0.0 {
                recalc_info = 1.0 - (ds4_sum - ds2_sum) / denom;
            }
            recalc_info = recalc_info.max(0.0);
            
            // Round INFO natively before truncating to mirror C++
            let recalc_info_rounded = (recalc_info * 1000.0).round() / 1000.0;
            
            batch_af_strs.push(format_float_3(af));
            batch_info_strs.push(format_float_3(recalc_info_rounded));
        }

        let mut tsv_line = format!("{}\t{}\t{}\t{}", fields[0], fields[1], fields[3], fields[4]);
        for i in 0..num_samples.len() {
            tsv_line.push_str(&format!("\t{}\t{}", batch_af_strs[i], batch_info_strs[i]));
        }
        writeln!(tsv_out, "{}", tsv_line).unwrap();

        let mut agg_af_sum = 0.0_f64;
        let mut numerator = 0.0_f64;

        // 2. Aggregate across batches using f64 parsed from the truncated strings, mirroring the Python RecomputeAndAnnotate
        for i in 0..num_samples.len() {
            let n = num_samples[i] as f64;
            let b_af = batch_af_strs[i].parse::<f64>().unwrap();
            let b_info = batch_info_strs[i].parse::<f64>().unwrap();

            agg_af_sum += b_af * n;
            numerator += (1.0 - b_info) * 2.0 * n * b_af * (1.0 - b_af);
        }

        let total_af = agg_af_sum / (total_samples as f64);
        let denom = 2.0 * (total_samples as f64) * total_af * (1.0 - total_af);
        
        let mut total_info = 1.0_f64;
        if total_af > 0.0 && total_af < 1.0 && denom > 0.0 {
            total_info = 1.0 - (numerator / denom);
        }

        let final_af = round_to_3_sig_figs(total_af);
        let final_info = round_to_3_sig_figs(total_info);

        let mut new_info = format!("ID={}", atomic_id);
        if !true_raf.is_empty() {
            new_info.push_str(&format!(";RAF={}", true_raf));
        }
        
        // 3. Format using standard {} serialization to preserve significant figures perfectly
        new_info.push_str(&format!(";AF={};INFO={}", final_af, final_info));
        
        let new_info_str = new_info.clone();
        fields[7] = &new_info_str;

        writeln!(vcf_out, "{}", fields.join("\t")).unwrap();

        // Check if we hit the requested limit
        if let Some(max) = max_records {
            if records_processed >= max {
                eprintln!("Hit the requested limit of {} records. Stopping early.", max);
                break;
            }
        }

        // Clear the buffer at the end of the loop so it can be reused for the next line
        line.clear();
    }
    
    let total_elapsed = start_time.elapsed().as_secs();
    let t_hours = total_elapsed / 3600;
    let t_mins = (total_elapsed % 3600) / 60;
    let t_secs = total_elapsed % 60;
    eprintln!(
        "Finished! Processed a total of {} input records in {:02}:{:02}:{:02}.",
        records_processed, t_hours, t_mins, t_secs
    );
}
