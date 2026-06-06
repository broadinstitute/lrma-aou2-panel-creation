use flate2::read::GzDecoder;
use std::collections::{HashMap, HashSet};
use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

// Struct to hold VCF record data for the current group
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
        Box::new(BufReader::new(GzDecoder::new(file)))
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

fn process_group(
    group_lines: &[String],
    chrom_to_variants: &HashMap<String, HashMap<String, (u32, String, String, String)>>,
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
            println!("{}", fields.join("\t"));
        }
        return;
    }

    let num_samples = parsed_lines[0].len() - 9;
    let chrom = parsed_lines[0][0].to_string();

    let mut records: Vec<Record> = Vec::with_capacity(parsed_lines.len());

    // 1. Pre-parse and initialize records
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

        let mut num_const = 0;
        let mut atomic_ids = HashSet::new();
        let mut is_known = false;

        if let Some(id_str) = info_map.get("ID") {
            // FIX: Bind the new String to a variable so it lives long enough for the splits to borrow it.
            let replaced_id = id_str.replace(',', ":");
            let constituents: Vec<&str> = replaced_id.split(':').map(|s| s.trim()).collect();
            num_const = constituents.len();
            for j in constituents {
                if let Some(chrom_map) = chrom_to_variants.get(&chrom) {
                    if chrom_map.contains_key(j) {
                        atomic_ids.insert(j.to_string());
                        is_known = true;
                    }
                }
            }
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

    // 2. Annotate and revert collisions
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
                a.1.cmp(&b.1) // Parsimony: lower constituent count first
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

    // 3. Output passthroughs and gather atomic targets
    let mut all_atomic_vars = HashSet::new();
    for rec in &records {
        if !rec.is_known {
            let mut out = vec![
                rec.chrom.clone(), rec.pos.clone(), rec.id.clone(), 
                rec.ref_seq.clone(), rec.alt.clone(), rec.qual.clone(), 
                rec.filter.clone(), rec.info_raw.clone(), rec.fmt_raw.clone()
            ];
            out.extend(rec.samples.clone());
            println!("{}", out.join("\t"));
        } else {
            for j in &rec.atomic_ids {
                if let Some(data) = chrom_to_variants.get(&rec.chrom).and_then(|m| m.get(j)) {
                    all_atomic_vars.insert((j.clone(), data.0));
                }
            }
        }
    }

    // 4. Pop Bubbles
    let mut sorted_atomic_vars: Vec<_> = all_atomic_vars.into_iter().collect();
    sorted_atomic_vars.sort_by_key(|x| x.1);

    for (var_id, coord) in sorted_atomic_vars {
        let template_idx = records.iter().position(|r| r.atomic_ids.contains(&var_id));
        if template_idx.is_none() { continue; }
        let t_idx = template_idx.unwrap();
        let t_rec = &records[t_idx];
        let var_data = &chrom_to_variants[&t_rec.chrom][&var_id];

        let mut new_info = vec![format!("ID={}", var_id)];
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
            .filter(|(_, r)| r.atomic_ids.contains(&var_id))
            .map(|(i, _)| i)
            .collect();

        for s in 0..num_samples {
            let mut max_gp_for_ds_gp = -1.0_f32;
            let mut best_ds = ".".to_string();
            let mut best_gp = ".".to_string();
            let mut best_gq = ".".to_string();
            
            let mut hap0_winner = false;
            let mut hap1_winner = false;
            let mut h0_losers = Vec::new();
            let mut h1_losers = Vec::new();

            for &i in &lines_with_var {
                let rec = &records[i];
                let s_vals: Vec<&str> = rec.samples[s].split(':').collect();

                if let Some(&gt_idx) = rec.fmt_idx.get("GT") {
                    let gt = s_vals[gt_idx];
                    if let Some(&cl_idx) = rec.fmt_idx.get("CL") {
                        if gt.contains('|') {
                            let a_gt: Vec<&str> = gt.split('|').collect();
                            let a_cl: Vec<&str> = s_vals[cl_idx].split(',').collect();

                            if a_gt[0] == "1" && a_cl[0] == "0" { hap0_winner = true; }
                            else if a_cl[0] != "0" { h0_losers.push(a_cl[0].parse::<u8>().unwrap_or(0)); }

                            if a_gt[1] == "1" && a_cl[1] == "0" { hap1_winner = true; }
                            else if a_cl[1] != "0" { h1_losers.push(a_cl[1].parse::<u8>().unwrap_or(0)); }
                        }
                    }
                }

                let curr_max_gp = rec.gps[s];
                if curr_max_gp > max_gp_for_ds_gp {
                    max_gp_for_ds_gp = curr_max_gp;
                    if let Some(&ds_idx) = rec.fmt_idx.get("DS") { best_ds = s_vals[ds_idx].to_string(); }
                    if let Some(&gp_idx) = rec.fmt_idx.get("GP") { best_gp = s_vals[gp_idx].to_string(); }
                }

                if let Some(&gq_idx) = rec.fmt_idx.get("GQ") {
                    let gq_str = s_vals[gq_idx];
                    if gq_str != "." {
                        let gq_val = gq_str.parse::<f32>().unwrap_or(-1.0);
                        let best_gq_val = best_gq.parse::<f32>().unwrap_or(-1.0);
                        if best_gq == "." || gq_val > best_gq_val {
                            best_gq = gq_str.to_string();
                        }
                    }
                }
            }

            let hap0_val = if hap0_winner { "1" } else { "0" };
            let hap1_val = if hap1_winner { "1" } else { "0" };

            let hap0_cl = if hap0_winner || h0_losers.is_empty() { "0".to_string() } else { h0_losers.into_iter().max().unwrap().to_string() };
            let hap1_cl = if hap1_winner || h1_losers.is_empty() { "0".to_string() } else { h1_losers.into_iter().max().unwrap().to_string() };

            let mut s_out = vec![format!("{}|{}", hap0_val, hap1_val)];
            if fmt_out.contains(&"DS") { s_out.push(best_ds); }
            if fmt_out.contains(&"GP") { s_out.push(best_gp); }
            if fmt_out.contains(&"GQ") { s_out.push(best_gq); }
            if fmt_out.contains(&"CL") { s_out.push(format!("{},{}", hap0_cl, hap1_cl)); }

            vcf_line.push(s_out.join(":"));
        }
        println!("{}", vcf_line.join("\t"));
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: cat <multiallelic VCF> | {} <biallelic VCF>", args[0]);
        std::process::exit(1);
    }
    
    let vcf_path = &args[1];
    
    // chrom -> ID -> (pos, orig_id, ref, alt)
    let mut chrom_to_variants: HashMap<String, HashMap<String, (u32, String, String, String)>> = HashMap::new();
    
    let reader = smart_open(vcf_path);
    for line_result in reader.lines() {
        let line = line_result.unwrap();
        if line.starts_with('#') { continue; }
        
        let fields: Vec<&str> = line.trim_end().split('\t').collect();
        if fields.len() < 8 { continue; }
        
        let mut id_val = "";
        for item in fields[7].split(';') {
            if item.starts_with("ID=") {
                id_val = &item[3..];
                break;
            }
        }
        
        if id_val.is_empty() { continue; }
        let id_keys: Vec<&str> = id_val.split(',').collect();
        let main_id = id_keys[0];
        
        let pos: u32 = fields[1].parse().unwrap();
        chrom_to_variants.entry(fields[0].to_string())
            .or_default()
            .insert(main_id.to_string(), (pos, fields[2].to_string(), fields[3].to_string(), fields[4].to_string()));
    }

    let stdin = io::stdin();
    let mut cl_header_added = false;
    let mut current_pos: Option<String> = None;
    let mut group = Vec::new();

    for line_result in stdin.lock().lines() {
        let line = line_result.unwrap();
        if line.starts_with('#') {
            if line.contains("INFO=<ID=AK") || line.contains("FORMAT=<ID=GL") || line.contains("FORMAT=<ID=KC") {
                continue;
            }
            if line.starts_with("#CHROM") && !cl_header_added {
                println!("##FORMAT=<ID=CL,Number=2,Type=Integer,Description=\"Collision indicator array (Hap0,Hap1): 0=optimal, 1=lost by GP, 2=lost by parsimony, 3=lost by stable sort\">");
                cl_header_added = true;
            } else if line.starts_with("##FORMAT=") && !cl_header_added {
                println!("##FORMAT=<ID=CL,Number=2,Type=Integer,Description=\"Collision indicator array (Hap0,Hap1): 0=optimal, 1=lost by GP, 2=lost by parsimony, 3=lost by stable sort\">");
                cl_header_added = true;
            }
            println!("{}", line);
            continue;
        }

        let fields: Vec<&str> = line.splitn(3, '\t').collect();
        if fields.len() < 2 { continue; }
        let pos = fields[1].to_string();

        if current_pos.is_none() {
            current_pos = Some(pos.clone());
        }

        // FIX: Dereferenced current_pos with *
        if pos != *current_pos.as_ref().unwrap() {
            process_group(&group, &chrom_to_variants);
            group.clear();
            current_pos = Some(pos);
        }
        group.push(line);
    }

    if !group.is_empty() {
        process_group(&group, &chrom_to_variants);
    }
}
