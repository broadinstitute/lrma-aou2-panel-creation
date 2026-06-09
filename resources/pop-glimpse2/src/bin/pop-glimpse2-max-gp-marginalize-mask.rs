use flate2::read::MultiGzDecoder; 
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write, BufWriter};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

struct Record {
    chrom: String,
    pos: String,
    id: String,
    ref_seq: String,
    alt: String,
    qual: String,
    filter: String,
    info_raw: String,
    info_map: HashMap<String, String>,
    fmt_raw: String,
    fmt_idx: HashMap<String, usize>,
    samples: Vec<String>,
    atomic_ids: HashSet<String>,
    num_const: usize,
    is_known: bool,
    gps: Vec<f32>,
}

fn smart_open(filename: &str) -> Box<dyn BufRead> {
    let file = File::open(filename).expect("Cannot open file");
    if filename.ends_with(".gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file))) 
    } else {
        Box::new(BufReader::new(file))
    }
}

fn parse_info(info_str: &str) -> HashMap<String, String> {
    let mut map = HashMap::new();
    for item in info_str.split(';') {
        if let Some((k, v)) = item.split_once('=') {
            map.insert(k.to_string(), v.to_string());
        }
    }
    map
}

// Formats floats strictly to 3 decimal places and trims trailing zeros to match htslib/GLIMPSE2 standard formatting.
fn format_float(val: f32) -> String {
    let s = format!("{:.3}", val);
    let trimmed = s.trim_end_matches('0').trim_end_matches('.');
    if trimmed.is_empty() {
        "0".to_string()
    } else {
        trimmed.to_string()
    }
}

fn process_group(
    group_lines: &[String],
    chrom_to_variants: &HashMap<String, HashMap<String, (u32, String, String, String)>>,
    printed_phys_vars: &mut HashMap<(String, u32, String, String, String), ()>,
    out_handle: &mut impl Write,
) {
    if group_lines.is_empty() {
        return;
    }

    let mut parsed_lines: Vec<Vec<&str>> = group_lines
        .iter()
        .map(|line| line.trim_end().split('\t').collect())
        .collect();

    if parsed_lines[0].len() <= 9 {
        for fields in parsed_lines {
            writeln!(out_handle, "{}", fields.join("\t")).unwrap();
        }
        return;
    }

    let num_samples = parsed_lines[0].len() - 9;
    let chrom = parsed_lines[0][0].to_string();

    let mut records: Vec<Record> = Vec::with_capacity(parsed_lines.len());

    for fields in parsed_lines.iter_mut() {
        let info_map = parse_info(fields[7]);
        
        let mut fmt: Vec<&str> = fields[8].split(':').collect();
        let added_cl = if !fmt.contains(&"CL") {
            fmt.push("CL");
            true
        } else {
            false
        };
        let fmt_str = fmt.join(":");
        let mut fmt_idx = HashMap::new();
        for (idx, k) in fmt.iter().enumerate() {
            fmt_idx.insert(k.to_string(), idx);
        }

        let mut is_known = true;
        let mut num_const = 0;
        let mut atomic_ids = HashSet::new();

        if let Some(id_str) = info_map.get("ID") {
            let replaced_id = id_str.replace(',', ":");
            let constituents: Vec<&str> = replaced_id.split(':').map(|s| s.trim()).collect();
            num_const = constituents.len();
            
            for j in &constituents {
                if let Some(chrom_map) = chrom_to_variants.get(&chrom) {
                    if !chrom_map.contains_key(*j) {
                        eprintln!("WARNING: Assigned ID '{}' not found in resource VCF. Path containing this variant at {}:{} will be passed through unpopped.", j, chrom, fields[1]);
                        is_known = false;
                        break;
                    }
                } else {
                    eprintln!("WARNING: Chromosome '{}' not found in resource VCF when looking up assigned ID '{}'. Path at {}:{} will be passed through unpopped.", chrom, j, chrom, fields[1]);
                    is_known = false;
                    break;
                }
            }
            
            if is_known {
                for j in constituents {
                    atomic_ids.insert(j.to_string());
                }
            }
        } else {
            is_known = false;
        }

        let mut samples = Vec::with_capacity(num_samples);
        let mut gps = vec![-1.0; num_samples];
        let gp_idx = fmt_idx.get("GP").copied();

        for s in 0..num_samples {
            let mut s_vals: Vec<&str> = fields[9 + s].split(':').collect();
            if added_cl && s_vals.len() < fmt.len() {
                s_vals.push("0,0");
            }
            if let Some(idx) = gp_idx {
                if idx < s_vals.len() && s_vals[idx] != "." {
                    let max_gp = s_vals[idx]
                        .split(',')
                        .filter_map(|x| x.parse::<f32>().ok())
                        .fold(-1.0, f32::max);
                    gps[s] = max_gp;
                }
            }
            samples.push(s_vals.join(":"));
        }

        records.push(Record {
            chrom: chrom.clone(),
            pos: fields[1].to_string(),
            id: fields[2].to_string(),
            ref_seq: fields[3].to_string(),
            alt: fields[4].to_string(),
            qual: fields[5].to_string(),
            filter: fields[6].to_string(),
            info_raw: fields[7].to_string(),
            info_map,
            fmt_raw: fmt_str,
            fmt_idx,
            samples,
            atomic_ids,
            num_const,
            is_known,
            gps,
        });
    }

    for s in 0..num_samples {
        let mut hap0_ones: Vec<(f32, usize, usize)> = Vec::new();
        let mut hap1_ones: Vec<(f32, usize, usize)> = Vec::new();

        for (i, rec) in records.iter().enumerate() {
            if let Some(&gt_idx) = rec.fmt_idx.get("GT") {
                let s_vals: Vec<&str> = rec.samples[s].split(':').collect();
                if gt_idx < s_vals.len() {
                    let gt = s_vals[gt_idx];
                    if gt.contains('|') {
                        let alleles: Vec<&str> = gt.split('|').collect();
                        if alleles[0] == "1" { hap0_ones.push((rec.gps[s], rec.num_const, i)); }
                        if alleles[1] == "1" { hap1_ones.push((rec.gps[s], rec.num_const, i)); }
                    }
                }
            }
        }

        let sort_logic = |a: &(f32, usize, usize), b: &(f32, usize, usize)| {
            let gp_cmp = b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal);
            if gp_cmp == std::cmp::Ordering::Equal {
                a.1.cmp(&b.1) 
            } else {
                gp_cmp
            }
        };

        if hap0_ones.len() > 1 {
            hap0_ones.sort_by(sort_logic);
            let (w_gp, w_const, _) = hap0_ones[0];
            for &(l_gp, l_const, i) in hap0_ones.iter().skip(1) {
                let reason = if w_gp > l_gp { "1" } else if w_const < l_const { "2" } else { "3" };
                let rec = &mut records[i];
                let cl_idx = rec.fmt_idx["CL"];
                let gt_idx = rec.fmt_idx["GT"];
                let mut s_vals: Vec<String> = rec.samples[s].split(':').map(|x| x.to_string()).collect();
                
                let mut cl_arr: Vec<String> = s_vals[cl_idx].split(',').map(|x| x.to_string()).collect();
                cl_arr[0] = reason.to_string();
                s_vals[cl_idx] = cl_arr.join(",");
                
                let mut gt_arr: Vec<String> = s_vals[gt_idx].split('|').map(|x| x.to_string()).collect();
                gt_arr[0] = "0".to_string();
                s_vals[gt_idx] = gt_arr.join("|");
                
                rec.samples[s] = s_vals.join(":");
            }
        }

        if hap1_ones.len() > 1 {
            hap1_ones.sort_by(sort_logic);
            let (w_gp, w_const, _) = hap1_ones[0];
            for &(l_gp, l_const, i) in hap1_ones.iter().skip(1) {
                let reason = if w_gp > l_gp { "1" } else if w_const < l_const { "2" } else { "3" };
                let rec = &mut records[i];
                let cl_idx = rec.fmt_idx["CL"];
                let gt_idx = rec.fmt_idx["GT"];
                let mut s_vals: Vec<String> = rec.samples[s].split(':').map(|x| x.to_string()).collect();
                
                let mut cl_arr: Vec<String> = s_vals[cl_idx].split(',').map(|x| x.to_string()).collect();
                cl_arr[1] = reason.to_string();
                s_vals[cl_idx] = cl_arr.join(",");
                
                let mut gt_arr: Vec<String> = s_vals[gt_idx].split('|').map(|x| x.to_string()).collect();
                gt_arr[1] = "0".to_string();
                s_vals[gt_idx] = gt_arr.join("|");
                
                rec.samples[s] = s_vals.join(":");
            }
        }
    }

    let mut all_atomic_vars = HashSet::new();
    for rec in &records {
        if !rec.is_known {
            let mut out = vec![
                rec.chrom.clone(), rec.pos.clone(), rec.id.clone(), 
                rec.ref_seq.clone(), rec.alt.clone(), rec.qual.clone(), 
                rec.filter.clone(), rec.info_raw.clone(), rec.fmt_raw.clone()
            ];
            out.extend(rec.samples.clone());
            writeln!(out_handle, "{}", out.join("\t")).unwrap();
        } else {
            for j in &rec.atomic_ids {
                if let Some(data) = chrom_to_variants.get(&rec.chrom).and_then(|m| m.get(j)) {
                    all_atomic_vars.insert((j.clone(), data.0));
                }
            }
        }
    }

    let mut sorted_atomic_vars: Vec<_> = all_atomic_vars.into_iter().collect();
    sorted_atomic_vars.sort_by_key(|x| x.1);

    for (assigned_id, coord) in sorted_atomic_vars {
        let template_idx = records.iter().position(|r| r.atomic_ids.contains(&assigned_id));
        if template_idx.is_none() { continue; }
        let t_idx = template_idx.unwrap();
        let t_rec = &records[t_idx];
        
        let var_data = &chrom_to_variants[&t_rec.chrom][&assigned_id];

        let phys_sig = (t_rec.chrom.clone(), coord, var_data.1.clone(), var_data.2.clone(), var_data.3.clone());
        if printed_phys_vars.contains_key(&phys_sig) {
            continue; 
        }
        printed_phys_vars.insert(phys_sig, ());

        let mut new_info = vec![format!("ID={}", assigned_id)];
        for k in ["MA", "UK", "RAF", "AF", "INFO"] {
            if let Some(v) = t_rec.info_map.get(k) {
                new_info.push(format!("{}={}", k, v));
            }
        }

        let mut fmt_out = vec!["GT"];
        if t_rec.fmt_idx.contains_key("DS") { fmt_out.push("DS"); }
        if t_rec.fmt_idx.contains_key("GP") { fmt_out.push("GP"); }
        if t_rec.fmt_idx.contains_key("GQ") { fmt_out.push("GQ"); }
        if t_rec.fmt_idx.contains_key("CL") { fmt_out.push("CL"); }

        let mut vcf_line = vec![
            t_rec.chrom.clone(),
            coord.to_string(),
            var_data.1.clone(), 
            var_data.2.clone(), 
            var_data.3.clone(), 
            ".".to_string(),
            ".".to_string(),
            new_info.join(";"),
            fmt_out.join(":"),
        ];

        let lines_with_var: Vec<usize> = records.iter()
            .enumerate()
            .filter(|(_, r)| r.atomic_ids.contains(&assigned_id))
            .map(|(i, _)| i)
            .collect();

        for s in 0..num_samples {
            let mut p0 = 0.0_f32;
            let mut p1 = 0.0_f32;
            let mut has_ds = false;
            
            let mut hap0_winner = false;
            let mut hap1_winner = false;
            let mut h0_losers = Vec::new();
            let mut h1_losers = Vec::new();

            for &i in &lines_with_var {
                let rec = &records[i];
                let s_vals: Vec<&str> = rec.samples[s].split(':').collect();

                let mut a_gt_0 = false;
                let mut a_gt_1 = false;
                let mut p0_frac = 1.0_f32;
                let mut p1_frac = 1.0_f32;

                if let Some(&gt_idx) = rec.fmt_idx.get("GT") {
                    let gt = s_vals[gt_idx];
                    if let Some(&cl_idx) = rec.fmt_idx.get("CL") {
                        if gt.contains('|') {
                            let a_gt: Vec<&str> = gt.split('|').collect();
                            let a_cl: Vec<&str> = s_vals[cl_idx].split(',').collect();

                            a_gt_0 = a_gt[0] == "1";
                            a_gt_1 = a_gt[1] == "1";

                            if a_gt_0 && a_cl[0] == "0" { hap0_winner = true; }
                            else if a_cl[0] != "0" { h0_losers.push(a_cl[0].parse::<u8>().unwrap_or(0)); }

                            if a_gt_1 && a_cl[1] == "0" { hap1_winner = true; }
                            else if a_cl[1] != "0" { h1_losers.push(a_cl[1].parse::<u8>().unwrap_or(0)); }

                            if (a_gt_0 || a_cl[0] != "0") && a_cl[0] != "0" { p0_frac = 0.0; }
                            if (a_gt_1 || a_cl[1] != "0") && a_cl[1] != "0" { p1_frac = 0.0; }
                        }
                    }
                }

                if let Some(&ds_idx) = rec.fmt_idx.get("DS") {
                    if let Ok(ds_val) = s_vals[ds_idx].parse::<f32>() {
                        has_ds = true;
                        if a_gt_0 && !a_gt_1 {
                            p0 += ds_val * p0_frac;
                        } else if !a_gt_0 && a_gt_1 {
                            p1 += ds_val * p1_frac;
                        } else {
                            p0 += (ds_val / 2.0) * p0_frac;
                            p1 += (ds_val / 2.0) * p1_frac;
                        }
                    }
                }
            }

            let hap0_val = if hap0_winner { "1" } else { "0" };
            let hap1_val = if hap1_winner { "1" } else { "0" };

            let hap0_cl = if hap0_winner || h0_losers.is_empty() { "0".to_string() } else { h0_losers.into_iter().max().unwrap().to_string() };
            let hap1_cl = if hap1_winner || h1_losers.is_empty() { "0".to_string() } else { h1_losers.into_iter().max().unwrap().to_string() };

            p0 = p0.clamp(0.0, 1.0);
            p1 = p1.clamp(0.0, 1.0);

            let mut s_out = vec![format!("{}|{}", hap0_val, hap1_val)];

            if fmt_out.contains(&"DS") {
                if has_ds {
                    let final_ds = p0 + p1;
                    s_out.push(format_float(final_ds));
                } else {
                    s_out.push(".".to_string());
                }
            }

            if fmt_out.contains(&"GP") || fmt_out.contains(&"GQ") {
                if has_ds {
                    let gp0_raw = (1.0 - p0) * (1.0 - p1);
                    let gp1_raw = p0 * (1.0 - p1) + (1.0 - p0) * p1;
                    let gp2_raw = p0 * p1;

                    if fmt_out.contains(&"GP") {
                        // Integer-based Largest Bucket algorithm to ensure the rounded percentages perfectly sum to 1.000
                        let mut v0 = (gp0_raw * 1000.0).round() as i32;
                        let mut v1 = (gp1_raw * 1000.0).round() as i32;
                        let mut v2 = (gp2_raw * 1000.0).round() as i32;
                        
                        let diff = 1000 - (v0 + v1 + v2);
                        if diff != 0 {
                            if v0 >= v1 && v0 >= v2 {
                                v0 += diff;
                            } else if v1 >= v0 && v1 >= v2 {
                                v1 += diff;
                            } else {
                                v2 += diff;
                            }
                        }
                        
                        s_out.push(format!("{},{},{}", 
                            format_float(v0 as f32 / 1000.0), 
                            format_float(v1 as f32 / 1000.0), 
                            format_float(v2 as f32 / 1000.0)
                        ));
                    }
                    if fmt_out.contains(&"GQ") {
                        let called_idx = if hap0_winner && hap1_winner { 2 } else if hap0_winner || hap1_winner { 1 } else { 0 };
                        let p_correct = if called_idx == 2 { gp2_raw } else if called_idx == 1 { gp1_raw } else { gp0_raw };
                        let p_error = (1.0 - p_correct).max(1e-10);
                        let gq = (-10.0 * p_error.log10()).clamp(0.0, 999.0);
                        s_out.push(format!("{:.0}", gq));
                    }
                } else {
                    if fmt_out.contains(&"GP") { s_out.push(".".to_string()); }
                    if fmt_out.contains(&"GQ") { s_out.push(".".to_string()); }
                }
            }

            if fmt_out.contains(&"CL") { s_out.push(format!("{},{}", hap0_cl, hap1_cl)); }

            vcf_line.push(s_out.join(":"));
        }
        writeln!(out_handle, "{}", vcf_line.join("\t")).unwrap();
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: cat <multiallelic VCF> | {} <biallelic VCF>", args[0]);
        std::process::exit(1);
    }
    
    let vcf_path = &args[1];
    
    eprintln!("Loading resource VCF into memory...");
    let mut chrom_to_variants: HashMap<String, HashMap<String, (u32, String, String, String)>> = HashMap::new();
    
    let reader = smart_open(vcf_path);
    for line_result in reader.lines() {
        let line = line_result.unwrap();
        if line.starts_with('#') { continue; }
        
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        if fields.len() < 8 { continue; }
        
        for item in fields[7].split(';') {
            if let Some(id_val) = item.strip_prefix("ID=") {
                let r_pos: u32 = fields[1].parse().unwrap_or(0);
                let r_chrom = fields[0].to_string();
                let original_id = fields[2].to_string(); 
                let val = (r_pos, original_id, fields[3].to_string(), fields[4].to_string());
                
                chrom_to_variants.entry(r_chrom)
                    .or_default()
                    .insert(id_val.trim().to_string(), val);
                break;
            }
        }
    }
    eprintln!("Finished loading resource VCF.");

    let stdout = io::stdout();
    let mut out_handle = BufWriter::new(stdout.lock());

    let mut printed_phys_vars: HashMap<(String, u32, String, String, String), ()> = HashMap::new();

    let stdin = io::stdin();
    let mut cl_header_added = false;
    let mut current_pos: Option<String> = None;
    let mut current_chrom: String = String::new();
    let mut group: Vec<String> = Vec::new(); 
    let mut records_processed: usize = 0;

    eprintln!("Starting projection from stdin...");
    for line_result in stdin.lock().lines() {
        let line = line_result.unwrap();
        if line.starts_with('#') {
            if line.contains("INFO=<ID=AK") || line.contains("FORMAT=<ID=GL") || line.contains("FORMAT=<ID=KC") {
                continue;
            }
            if line.starts_with("#CHROM") && !cl_header_added {
                writeln!(out_handle, "##FORMAT=<ID=CL,Number=2,Type=Integer,Description=\"Collision indicator array (Hap0,Hap1): 0=optimal, 1=lost by GP, 2=lost by parsimony, 3=lost by stable sort\">").unwrap();
                cl_header_added = true;
            } else if line.starts_with("##FORMAT=") && !cl_header_added {
                writeln!(out_handle, "##FORMAT=<ID=CL,Number=2,Type=Integer,Description=\"Collision indicator array (Hap0,Hap1): 0=optimal, 1=lost by GP, 2=lost by parsimony, 3=lost by stable sort\">").unwrap();
                cl_header_added = true;
            }
            writeln!(out_handle, "{}", line).unwrap();
            continue;
        }

        let fields: Vec<&str> = line.splitn(3, '\t').collect();
        if fields.len() < 2 { continue; }
        
        let chrom = fields[0].to_string();
        let pos = fields[1].to_string();

        records_processed += 1;
        if records_processed % 10_000 == 0 {
            eprintln!("Processed {} input records... (Currently at {}:{})", records_processed, chrom, pos);
        }

        if current_pos.is_none() {
            current_pos = Some(pos.clone());
            current_chrom = chrom.clone();
        }

        if pos != *current_pos.as_ref().unwrap() || chrom != current_chrom {
            process_group(
                &group, 
                &chrom_to_variants, 
                &mut printed_phys_vars, 
                &mut out_handle
            );
            
            group.clear();
            
            let current_pos_u32 = pos.parse::<u32>().unwrap_or(0);
            printed_phys_vars.retain(|(c, p, _, _, _), _| c == &chrom && *p >= current_pos_u32.saturating_sub(50_000));
            
            current_pos = Some(pos);
            current_chrom = chrom;
        }
        group.push(line);
    }

    if !group.is_empty() {
        process_group(
            &group, 
            &chrom_to_variants, 
            &mut printed_phys_vars, 
            &mut out_handle
        );
    }
    
    out_handle.flush().unwrap();
    eprintln!("Finished! Processed a total of {} input records.", records_processed);
}
