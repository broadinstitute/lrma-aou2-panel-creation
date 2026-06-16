use clap::Parser;
use rust_htslib::bcf::{self, header::{Header, TagType}, Read};
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(author, version, about = "Strictly pastes BCF formats across identical variants with strict POS subsetting")]
struct Args {
    /// Region to subset (e.g., chr1:10000-20000). Behaves like --regions-overlap pos
    #[arg(short, long)]
    region: Option<String>,

    /// Number of BGZF compression threads for the output writer
    #[arg(short, long, default_value_t = 1)]
    threads: usize,

    /// INFO tags to RETAIN and check for equality (all other INFO tags are permanently stripped)
    #[arg(long, value_delimiter = ',')]
    info: Vec<String>,

    /// FORMAT tags to dynamically extract and paste horizontally
    #[arg(long, value_delimiter = ',')]
    format: Vec<String>,

    /// Output BCF file path
    #[arg(short, long)]
    output: String,

    /// Input BCF files
    #[arg(required = true)]
    inputs: Vec<String>,
}

enum FormatData {
    Integer(Vec<i32>),
    Float(Vec<f32>),
    String(Vec<Vec<u8>>),
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // 1. Initialize readers
    let mut readers: Vec<bcf::IndexedReader> = args.inputs.iter()
        .map(|path| bcf::IndexedReader::from_path(path).expect("Failed to open BCF"))
        .collect();

    let num_files = readers.len();
    if num_files == 0 {
        panic!("No input files provided.");
    }

    // 2. Build the merged header natively
    let mut merged_header = Header::from_template(readers[0].header());
    for reader in readers.iter().skip(1) {
        for sample in reader.header().samples() {
            merged_header.push_sample(sample);
        }
    }

    // 3. Initialize Output Writer and apply threads
    let mut writer = bcf::Writer::from_path(&args.output, &merged_header, true, bcf::Format::Bcf)?;
    if args.threads > 1 {
        writer.set_threads(args.threads).expect("Failed to set writer threads");
    }

    // 4. Region Subsetting Setup & Coordinate Parsing
    let mut req_start = i64::MIN;
    let mut req_end = i64::MAX;

    if let Some(r) = &args.region {
        let (chrom, coords) = r.split_once(':').unwrap_or((r, ""));
        let rid = readers[0].header().name2rid(chrom.as_bytes())
            .unwrap_or_else(|_| panic!("Chromosome {} not found in header", chrom));

        let mut start: u64 = 0;
        let mut end: Option<u64> = None;

        if !coords.is_empty() {
            let parts: Vec<&str> = coords.split('-').collect();
            if !parts.is_empty() && !parts[0].is_empty() {
                let s = parts[0].replace(",", "").parse::<u64>().unwrap_or(1);
                start = s.saturating_sub(1);
                req_start = start as i64;
            }
            if parts.len() > 1 && !parts[1].is_empty() {
                let e = parts[1].replace(",", "").parse::<u64>().unwrap_or(u64::MAX);
                end = Some(e);
                req_end = e.saturating_sub(1) as i64;
            } else if parts.len() == 1 {
                req_end = start as i64; 
            }
        }

        for reader in readers.iter_mut() {
            reader.fetch(rid, start, end).unwrap_or_else(|_| panic!("Failed to fetch region {}", r));
        }
    }

    let mut base_record = writer.empty_record();
    let mut side_records: Vec<bcf::Record> = (1..num_files).map(|_| writer.empty_record()).collect();

    let mut record_count: usize = 0;
    let start_time = Instant::now();

    loop {
        // 5. Lockstep Read
        match readers[0].read(&mut base_record) {
            Some(Ok(_)) => {},
            Some(Err(e)) => panic!("Error reading base record: {}", e),
            None => break,
        }

        for (i, reader) in readers.iter_mut().skip(1).enumerate() {
            match reader.read(&mut side_records[i]) {
                Some(Ok(_)) => {},
                _ => panic!("FATAL: Side file {} ran out of records prematurely!", args.inputs[i + 1]),
            }
        }

        // 6. Strict POS filtering
        let pos = base_record.pos();
        if pos < req_start || pos > req_end {
            continue;
        }

        let rid = base_record.rid().unwrap_or(0);

        // --- PROGRESS METER ---
        record_count += 1;
        if record_count % 10_000 == 0 {
            let chrom_bytes = readers[0].header().rid2name(rid).unwrap_or(b"unknown");
            let chrom_str = String::from_utf8_lossy(chrom_bytes);
            eprintln!(
                "[Progress] Pasted {} records in {:.2?} | Current: {}:{}",
                record_count,
                start_time.elapsed(),
                chrom_str,
                pos + 1
            );
        }

        // 7. Core Equality Checks
        for side_record in &side_records {
            if base_record.rid() != side_record.rid() || pos != side_record.pos() {
                panic!("FATAL: Position mismatch! Base: {}:{}, Side: {}:{}", rid, pos, side_record.rid().unwrap_or(0), side_record.pos());
            }
            if base_record.id() != side_record.id() {
                panic!("FATAL: ID mismatch at {}:{}", rid, pos);
            }
            if base_record.alleles() != side_record.alleles() {
                panic!("FATAL: REF/ALT mismatch at {}:{}", rid, pos);
            }
        }

        // 8. THE PERFORMANCE HIT: Object Reconstruction
        // We spin up a completely blank record and manually copy the core data 
        let mut new_record = writer.empty_record();
        new_record.set_rid(base_record.rid());
        new_record.set_pos(pos);
        new_record.set_id(&base_record.id())?;
        new_record.set_alleles(&base_record.alleles())?;
        new_record.set_qual(base_record.qual());

        // 9. Strict INFO Checks & Selective Copying
        // Only the tags specified in `--info` are ported over to the new record
        for info_tag_str in &args.info {
            let tag = info_tag_str.as_bytes();
            let (tag_type, _) = readers[0].header().info_type(tag)
                .unwrap_or_else(|_| panic!("INFO tag {} not found in header", info_tag_str));

            match tag_type {
                TagType::Integer => {
                    let b_val = base_record.info(tag).integer().unwrap_or(None);
                    let b_slice = b_val.as_ref().map(|b| &**b);
                    for side_record in &mut side_records {
                        let s_val = side_record.info(tag).integer().unwrap_or(None);
                        let s_slice = s_val.as_ref().map(|b| &**b);
                        if b_slice != s_slice { panic!("INFO {} mismatch at {}:{}", info_tag_str, rid, pos); }
                    }
                    if let Some(val) = b_slice { new_record.push_info_integer(tag, val)?; }
                },
                TagType::Float => {
                    let b_val = base_record.info(tag).float().unwrap_or(None);
                    let b_slice = b_val.as_ref().map(|b| &**b);
                    for side_record in &mut side_records {
                        let s_val = side_record.info(tag).float().unwrap_or(None);
                        let s_slice = s_val.as_ref().map(|b| &**b);
                        let match_float = match (b_slice, s_slice) {
                            (Some(va), Some(vb)) => {
                                va.len() == vb.len() && va.iter().zip(vb.iter()).all(|(x, y)| (x.is_nan() && y.is_nan()) || x == y)
                            },
                            (None, None) => true,
                            _ => false,
                        };
                        if !match_float { panic!("INFO {} mismatch at {}:{}", info_tag_str, rid, pos); }
                    }
                    if let Some(val) = b_slice { new_record.push_info_float(tag, val)?; }
                },
                TagType::String => {
                    let b_val = base_record.info(tag).string().unwrap_or(None);
                    let b_slice = b_val.as_ref().map(|b| &**b);
                    for side_record in &mut side_records {
                        let s_val = side_record.info(tag).string().unwrap_or(None);
                        let s_slice = s_val.as_ref().map(|b| &**b);
                        if b_slice != s_slice { panic!("INFO {} mismatch at {}:{}", info_tag_str, rid, pos); }
                    }
                    if let Some(val) = b_slice { new_record.push_info_string(tag, val)?; }
                },
                TagType::Flag => {
                    let b_val = base_record.info(tag).flag().unwrap_or(false);
                    for side_record in &mut side_records {
                        let s_val = side_record.info(tag).flag().unwrap_or(false);
                        if b_val != s_val { panic!("INFO {} mismatch at {}:{}", info_tag_str, rid, pos); }
                    }
                    if b_val { new_record.push_info_flag(tag)?; }
                }
            }
        }

        // 10. Extract FORMAT arrays into memory
        let mut extracted_formats: Vec<(&[u8], FormatData)> = Vec::new();

        for tag_str in &args.format {
            let tag = tag_str.as_bytes();
            let tag_type = if tag == b"GT" {
                TagType::Integer
            } else {
                readers[0].header().format_type(tag).unwrap_or_else(|_| panic!("FORMAT tag {} not in header", tag_str)).0
            };

            match tag_type {
                TagType::Integer => {
                    let mut all_vals = Vec::new();
                    if let Ok(vals) = base_record.format(tag).integer() {
                        for v in vals.iter() { all_vals.extend_from_slice(v); }
                    }
                    for side_record in side_records.iter_mut() {
                        if let Ok(vals) = side_record.format(tag).integer() {
                            for v in vals.iter() { all_vals.extend_from_slice(v); }
                        }
                    }
                    extracted_formats.push((tag, FormatData::Integer(all_vals)));
                },
                TagType::Float => {
                    let mut all_vals = Vec::new();
                    if let Ok(vals) = base_record.format(tag).float() {
                        for v in vals.iter() { all_vals.extend_from_slice(v); }
                    }
                    for side_record in side_records.iter_mut() {
                        if let Ok(vals) = side_record.format(tag).float() {
                            for v in vals.iter() { all_vals.extend_from_slice(v); }
                        }
                    }
                    extracted_formats.push((tag, FormatData::Float(all_vals)));
                },
                TagType::String => {
                    let mut all_vals = Vec::new();
                    if let Ok(vals) = base_record.format(tag).string() {
                        for v in vals.iter() { all_vals.push(v.to_vec()); }
                    }
                    for side_record in side_records.iter_mut() {
                        if let Ok(vals) = side_record.format(tag).string() {
                            for v in vals.iter() { all_vals.push(v.to_vec()); }
                        }
                    }
                    extracted_formats.push((tag, FormatData::String(all_vals)));
                },
                _ => panic!("Unsupported FORMAT type for tag {}", tag_str),
            }
        }

        // 11. Push concatenated formats to the newly spawned record
        for (tag, data) in extracted_formats {
            match data {
                FormatData::Integer(vals) => { new_record.push_format_integer(tag, &vals)?; },
                FormatData::Float(vals) => { new_record.push_format_float(tag, &vals)?; },
                FormatData::String(vals) => {
                    let slices: Vec<&[u8]> = vals.iter().map(|v| v.as_slice()).collect();
                    new_record.push_format_string(tag, &slices)?;
                }
            }
        }

        // 12. Write the finalized clean record
        writer.write(&new_record)?;
    }

    Ok(())
}
