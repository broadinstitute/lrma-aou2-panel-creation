use clap::Parser;
use rust_htslib::bcf::{self, header::{Header, TagType}, Read};
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(author, version, about = "Strictly pastes BCF or VCF.GZ formats across identical variants")]
struct Args {
    #[arg(short, long)]
    region: Option<String>,

    #[arg(short, long, default_value_t = 1)]
    threads: usize,

    #[arg(long, value_delimiter = ',')]
    info: Vec<String>,

    #[arg(long, value_delimiter = ',')]
    format: Vec<String>,

    #[arg(short, long)]
    output: String,

    #[arg(required = true)]
    inputs: Vec<String>,
}

// Optimization 1: Cache the header dictionary lookups
struct TagSpec {
    bytes: Vec<u8>,
    ty: TagType,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    // Streaming Sequential Reader (No CSI index needed)
    let mut readers: Vec<bcf::Reader> = args.inputs.iter()
        .map(|path| bcf::Reader::from_path(path).expect("Failed to open input file"))
        .collect();

    let num_files = readers.len();
    if num_files == 0 {
        panic!("No input files provided.");
    }

    let mut merged_header = Header::from_template(readers[0].header());
    for reader in readers.iter().skip(1) {
        for sample in reader.header().samples() {
            merged_header.push_sample(sample);
        }
    }

    let mut writer = bcf::Writer::from_path(&args.output, &merged_header, false, bcf::Format::Bcf)?;
    if args.threads > 1 {
        writer.set_threads(args.threads).expect("Failed to set writer threads");
    }

    // Resolve Tag Types ONCE outside the loop
    let info_specs: Vec<TagSpec> = args.info.iter().map(|s| {
        let bytes = s.as_bytes().to_vec();
        let ty = readers[0].header().info_type(&bytes).unwrap().0;
        TagSpec { bytes, ty }
    }).collect();

    let format_specs: Vec<TagSpec> = args.format.iter().map(|s| {
        let bytes = s.as_bytes().to_vec();
        let ty = if bytes == b"GT" {
            TagType::Integer
        } else {
            readers[0].header().format_type(&bytes).unwrap().0
        };
        TagSpec { bytes, ty }
    }).collect();

    let mut req_rid: Option<u32> = None;
    let mut req_start = i64::MIN;
    let mut req_end = i64::MAX;

    if let Some(r) = &args.region {
        let (chrom, coords) = r.split_once(':').unwrap_or((r, ""));
        let rid = readers[0].header().name2rid(chrom.as_bytes())
            .unwrap_or_else(|_| panic!("Chromosome {} not found in header", chrom));
        
        req_rid = Some(rid);

        let mut start: u64 = 0;

        if !coords.is_empty() {
            let parts: Vec<&str> = coords.split('-').collect();
            if !parts.is_empty() && !parts[0].is_empty() {
                let s = parts[0].replace(",", "").parse::<u64>().unwrap_or(1);
                start = s.saturating_sub(1);
                req_start = start as i64;
            }
            if parts.len() > 1 && !parts[1].is_empty() {
                let e = parts[1].replace(",", "").parse::<u64>().unwrap_or(u64::MAX);
                req_end = e.saturating_sub(1) as i64;
            } else if parts.len() == 1 {
                req_end = start as i64; 
            }
        }
    }

    let mut base_record = writer.empty_record();
    let mut side_records: Vec<bcf::Record> = (1..num_files).map(|_| writer.empty_record()).collect();
    
    // Optimization 3: Hoist the record and scratch buffers out of the loop
    let mut new_record = writer.empty_record();
    let mut int_buffer: Vec<i32> = Vec::new();
    let mut float_buffer: Vec<f32> = Vec::new();
    let mut string_storage: Vec<Vec<u8>> = Vec::new();

    let mut record_count: usize = 0;
    let start_time = Instant::now();

    loop {
        match readers[0].read(&mut base_record) {
            Some(Ok(_)) => {},
            Some(Err(e)) => panic!("Error reading base record: {}", e),
            None => break,
        }

        for (i, reader) in readers.iter_mut().skip(1).enumerate() {
            match reader.read(&mut side_records[i]) {
                Some(Ok(_)) => {},
                _ => panic!("FATAL: Side file ran out of records prematurely!"),
            }
        }

        let rid = base_record.rid().unwrap_or(0);
        let pos = base_record.pos();

        // Streaming Region Check
        if let Some(target_rid) = req_rid {
            if rid < target_rid { continue; }
            if rid > target_rid || pos > req_end { break; }
            if pos < req_start { continue; }
        }

        record_count += 1;
        if record_count % 10_000 == 0 {
            let chrom_bytes = readers[0].header().rid2name(rid).unwrap_or(b"unknown".as_slice());
            let chrom_str = String::from_utf8_lossy(chrom_bytes);
            eprintln!("[Progress] Pasted {} records in {:.2?} | Current: {}:{}", record_count, start_time.elapsed(), chrom_str, pos + 1);
        }

        // Optimization 4: Cache the base ID and Alleles to avoid rebuilding strings
        let base_id = base_record.id();
        let base_alleles = base_record.alleles();
        
        for side_record in &side_records {
            if base_record.rid() != side_record.rid() || pos != side_record.pos() {
                panic!("FATAL: Position mismatch!");
            }
            if base_id != side_record.id() {
                panic!("FATAL: ID mismatch!");
            }
            if base_alleles != side_record.alleles() {
                panic!("FATAL: REF/ALT mismatch!");
            }
        }

        // Reuse the record structure
        new_record.clear();
        new_record.set_rid(Some(rid));
        new_record.set_pos(pos);
        new_record.set_id(&base_id)?;
        new_record.set_alleles(&base_alleles)?;
        new_record.set_qual(base_record.qual());

        // Process INFO tags using cached types
        for spec in &info_specs {
            let tag = spec.bytes.as_slice();
            match spec.ty {
                TagType::Integer => {
                    if let Some(val) = base_record.info(tag).integer().unwrap_or(None).as_ref().map(|b| &**b) {
                        new_record.push_info_integer(tag, val)?;
                    }
                },
                TagType::Float => {
                    if let Some(val) = base_record.info(tag).float().unwrap_or(None).as_ref().map(|b| &**b) {
                        new_record.push_info_float(tag, val)?;
                    }
                },
                TagType::String => {
                    if let Some(val) = base_record.info(tag).string().unwrap_or(None).as_ref().map(|b| &**b) {
                        new_record.push_info_string(tag, val)?;
                    }
                },
                TagType::Flag => {
                    if base_record.info(tag).flag().unwrap_or(false) {
                        new_record.push_info_flag(tag)?;
                    }
                }
            }
        }

        // Optimization 2: Push formats directly, using hoisted buffers
        for spec in &format_specs {
            let tag = spec.bytes.as_slice();
            match spec.ty {
                TagType::Integer => {
                    int_buffer.clear();
                    if let Ok(vals) = base_record.format(tag).integer() {
                        for v in vals.iter() { int_buffer.extend_from_slice(v); }
                    }
                    for side_record in side_records.iter_mut() {
                        if let Ok(vals) = side_record.format(tag).integer() {
                            for v in vals.iter() { int_buffer.extend_from_slice(v); }
                        }
                    }
                    new_record.push_format_integer(tag, &int_buffer)?;
                },
                TagType::Float => {
                    float_buffer.clear();
                    if let Ok(vals) = base_record.format(tag).float() {
                        for v in vals.iter() { float_buffer.extend_from_slice(v); }
                    }
                    for side_record in side_records.iter_mut() {
                        if let Ok(vals) = side_record.format(tag).float() {
                            for v in vals.iter() { float_buffer.extend_from_slice(v); }
                        }
                    }
                    new_record.push_format_float(tag, &float_buffer)?;
                },
                TagType::String => {
                    let mut n = 0usize;
                    if let Ok(vals) = base_record.format(tag).string() {
                        for v in vals.iter() {
                            if n == string_storage.len() { string_storage.push(Vec::new()); }
                            string_storage[n].clear();
                            string_storage[n].extend_from_slice(v);
                            n += 1;
                        }
                    }
                    for side_record in side_records.iter_mut() {
                        if let Ok(vals) = side_record.format(tag).string() {
                            for v in vals.iter() {
                                if n == string_storage.len() { string_storage.push(Vec::new()); }
                                string_storage[n].clear();
                                string_storage[n].extend_from_slice(v);
                                n += 1;
                            }
                        }
                    }
                    let slices: Vec<&[u8]> = string_storage[..n].iter().map(|v| v.as_slice()).collect();
                    new_record.push_format_string(tag, &slices)?;
                },
                _ => panic!("Unsupported FORMAT type"),
            }
        }

        writer.write(&new_record)?;
    }

    Ok(())
}
