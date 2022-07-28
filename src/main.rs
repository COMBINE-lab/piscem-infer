use ahash::AHashMap;
use anyhow::Context;
use clap::{Parser, Subcommand};
use libradicl::exit_codes;
use libradicl::rad_types;
use scroll::Pread;
use slog::{crit, info, o, Drain};
use std::fs::File;
use std::io::Write;
use std::io::{BufReader, Read};
use std::mem;
use std::path::PathBuf;

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum MappingType {
    Unmapped,
    SingleMapped,
    MappedFirstOrphan,
    MappedSecondOrphan,
    MappedPair,
}

impl MappingType {
    pub fn from_u8(t: u8) -> Self {
        match t {
            0 => MappingType::Unmapped,
            1 => MappingType::SingleMapped,
            2 => MappingType::MappedFirstOrphan,
            3 => MappingType::MappedSecondOrphan,
            4 => MappingType::MappedPair,
            _ => MappingType::Unmapped,
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub enum MappedFragmentOrientation {
    Reverse,
    Forward,
    ReverseReverse,
    ReverseForward,
    ForwardReverse,
    ForwardForward,
    Unknown,
}

impl MappedFragmentOrientation {
    pub fn from_u32_paired_status(n: u32, m: MappingType) -> Self {
        // if not paired, then we don't care about
        // the lowest order bit so shift it off
        if m == MappingType::SingleMapped {
            if (n >> 1) > 0 {
                MappedFragmentOrientation::Forward
            } else {
                MappedFragmentOrientation::Reverse
            }
        } else {
            match n {
                0 => MappedFragmentOrientation::ReverseReverse,
                1 => MappedFragmentOrientation::ReverseForward,
                2 => MappedFragmentOrientation::ForwardReverse,
                3 => MappedFragmentOrientation::ForwardForward,
                _ => MappedFragmentOrientation::Unknown,
            }
        }
    }
}

#[derive(Debug)]
pub struct MetaReadRecord {
    pub frag_type: u8,
    pub dirs: Vec<MappedFragmentOrientation>,
    pub refs: Vec<u32>,
    pub positions: Vec<u32>,
    pub frag_lengths: Vec<u16>,
}

#[derive(Debug)]
pub struct MetaChunk {
    pub nbytes: u32,
    pub nrec: u32,
    pub reads: Vec<MetaReadRecord>,
}

const MASK_LOWER_30_BITS: u32 = 0xC0000000;
const MASK_UPPER_2_BITS: u32 = 0x3FFFFFFF;

impl MetaChunk {
    pub fn read_header<T: Read>(reader: &mut T) -> (u32, u32) {
        let mut buf = [0u8; 8];

        reader.read_exact(&mut buf).unwrap();
        let nbytes = buf.pread::<u32>(0).unwrap();
        let nrec = buf.pread::<u32>(4).unwrap();
        (nbytes, nrec)
    }

    pub fn from_bytes<T: Read>(reader: &mut T, fmt: &rad_types::RadIntId) -> Self {
        let mut buf = [0u8; 8];

        reader.read_exact(&mut buf).unwrap();
        let nbytes = buf.pread::<u32>(0).unwrap();
        let nrec = buf.pread::<u32>(4).unwrap();
        let mut c = Self {
            nbytes,
            nrec,
            reads: Vec::with_capacity(nrec as usize),
        };

        for _ in 0..(nrec as usize) {
            c.reads.push(MetaReadRecord::from_bytes(reader, fmt));
        }

        c
    }

    /// peeks to the first record in the buffer `buf`, and returns
    /// the barcode and umi associated with this record.  It is assumed
    /// that there is at least one record present in the buffer.
    pub fn peek_record(
        buf: &[u8],
        bct: &rad_types::RadIntId,
        umit: &rad_types::RadIntId,
    ) -> (u64, u64) {
        let na_size = mem::size_of::<u32>();
        let bc_size = bct.bytes_for_type();

        let _na = buf.pread::<u32>(0).unwrap();

        let bc = match bct {
            rad_types::RadIntId::U8 => buf.pread::<u8>(na_size).unwrap() as u64,
            rad_types::RadIntId::U16 => buf.pread::<u16>(na_size).unwrap() as u64,
            rad_types::RadIntId::U32 => buf.pread::<u32>(na_size).unwrap() as u64,
            rad_types::RadIntId::U64 => buf.pread::<u64>(na_size).unwrap(),
        };
        let umi = match umit {
            rad_types::RadIntId::U8 => buf.pread::<u8>(na_size + bc_size).unwrap() as u64,
            rad_types::RadIntId::U16 => buf.pread::<u16>(na_size + bc_size).unwrap() as u64,
            rad_types::RadIntId::U32 => buf.pread::<u32>(na_size + bc_size).unwrap() as u64,
            rad_types::RadIntId::U64 => buf.pread::<u64>(na_size + bc_size).unwrap(),
        };
        (bc, umi)
    }
}

fn read_into_u64<T: Read>(reader: &mut T, rt: &rad_types::RadIntId) -> u64 {
    let mut rbuf = [0u8; 8];

    let v: u64 = match rt {
        rad_types::RadIntId::U8 => {
            reader.read_exact(&mut rbuf[0..1]).unwrap();
            rbuf.pread::<u8>(0).unwrap() as u64
        }
        rad_types::RadIntId::U16 => {
            reader.read_exact(&mut rbuf[0..2]).unwrap();
            rbuf.pread::<u16>(0).unwrap() as u64
        }
        rad_types::RadIntId::U32 => {
            reader.read_exact(&mut rbuf[0..4]).unwrap();
            rbuf.pread::<u32>(0).unwrap() as u64
        }
        rad_types::RadIntId::U64 => {
            reader.read_exact(&mut rbuf[0..8]).unwrap();
            rbuf.pread::<u64>(0).unwrap()
        }
    };
    v
}

impl MetaReadRecord {
    pub fn is_empty(&self) -> bool {
        self.refs.is_empty()
    }

    pub fn from_bytes<T: Read>(reader: &mut T, frag_map_t: &rad_types::RadIntId) -> Self {
        let mut rbuf = [0u8; 255];

        reader.read_exact(&mut rbuf[0..4]).unwrap();
        let na = rbuf.pread::<u32>(0).unwrap();
        let fmt = read_into_u64(reader, frag_map_t);
        let f = MappingType::from_u8(fmt as u8);

        let mut rec = Self {
            frag_type: fmt as u8,
            dirs: Vec::with_capacity(na as usize),
            refs: Vec::with_capacity(na as usize),
            positions: Vec::with_capacity(na as usize),
            frag_lengths: Vec::with_capacity(na as usize),
        };

        //println!("number of records : {:?}",na);

        for _ in 0..(na as usize) {
            reader.read_exact(&mut rbuf[0..4]).unwrap();
            let v = rbuf.pread::<u32>(0).unwrap();

            let dir_int = v & MASK_LOWER_30_BITS;
            let dir = MappedFragmentOrientation::from_u32_paired_status(dir_int, f);
            rec.dirs.push(dir);
            rec.refs.push(v & MASK_UPPER_2_BITS);
            // position
            reader.read_exact(&mut rbuf[0..4]).unwrap();
            let pos = rbuf.pread::<u32>(0).unwrap();
            rec.positions.push(pos);
            // length
            reader.read_exact(&mut rbuf[0..2]).unwrap();
            let flen = rbuf.pread::<u16>(0).unwrap();
            rec.frag_lengths.push(flen);
        }

        rec
    }
}

#[derive(Hash, PartialEq, Eq)]
struct EqLabel {
    pub targets: Vec<u32>,
}

fn conditional_means(freq: &[u32]) -> Vec<f64> {
    let mut cond_means = vec![0.0f64; freq.len()];
    let mut vals = vec![0.0f64; freq.len()];
    let mut multiplicities = vec![0.0f64; freq.len()];

    multiplicities[0] = freq[0] as f64;
    for i in 1..(freq.len()) {
        let v = freq[i] as f64;
        vals[i] = (v * i as f64) + vals[i - 1];
        multiplicities[i] = v + multiplicities[i - 1];
        if multiplicities[i] > 0.0f64 {
            cond_means[i] = vals[i] / multiplicities[i];
        }
    }

    cond_means
}

fn process<T: Read>(
    br: &mut BufReader<T>,
    nrec: usize,
    frag_map_t: &rad_types::RadIntId,
    tot_mappings: &mut usize,
    num_mapped_reads: &mut usize,
) -> (AHashMap<EqLabel, usize>, Vec<f64>) {
    let mut map: AHashMap<EqLabel, usize> = AHashMap::new();
    let mut frag_lengths = vec![0u32; 65536];
    let mut _unique_frags = 0u32;

    for _ in 0..nrec {
        let c = MetaChunk::from_bytes(br, frag_map_t);
        for mappings in &c.reads {
            let nm = mappings.positions.len();

            *tot_mappings += nm;
            *num_mapped_reads += 1;

            let eql = EqLabel {
                targets: mappings.refs.clone(),
            };
            map.entry(eql)
                .and_modify(|counter| *counter += 1)
                .or_insert(1);

            if nm == 1 {
                frag_lengths[mappings.frag_lengths[0] as usize] += 1;
                _unique_frags += 1;
            }
        }
    }

    let cond_means = conditional_means(&frag_lengths);
    (map, cond_means)
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// quantify from the rad file
    #[clap(arg_required_else_help = true)]
    Quant {
        /// input stem (i.e. without the .rad suffix)
        #[clap(short, long, value_parser)]
        input: PathBuf,
        /// where output should be written
        #[clap(short, long, value_parser)]
        output: PathBuf,
        /// max iterations to run the EM
        #[clap(short, long, default_value_t = 1500, value_parser)]
        max_iter: u32,
    },
}

/// quantifying from a metagenomic rad file
#[derive(Debug, Parser)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

fn read_ref_lengths<T: Read>(reader: &mut T) -> Result<Vec<u32>, std::io::Error> {
    let mut rbuf = [0u8; 8];

    // read type of length
    reader.read_exact(&mut rbuf[0..1])?;
    let lt = rbuf.pread::<u8>(0).unwrap() as u64;
    assert_eq!(lt, 3);

    // read length of the array
    reader.read_exact(&mut rbuf[0..4])?;
    let len = rbuf.pread::<u32>(0).unwrap() as usize;

    // read type of entry
    reader.read_exact(&mut rbuf[0..1])?;
    let et = rbuf.pread::<u8>(0).unwrap() as u64;
    assert_eq!(et, 3);

    let mut vec = Vec::<u32>::with_capacity(len);
    for _ in 0..len {
        reader.read_exact(&mut rbuf[0..4])?;
        let v = rbuf.pread::<u32>(0).unwrap();
        vec.push(v);
    }

    Ok(vec)
}

fn adjust_ref_lengths(ref_lens: &[u32], cond_means: &[f64]) -> Vec<f64> {
    let tmean = cond_means.last().unwrap();
    let el = ref_lens
        .iter()
        .map(|rli| {
            let rl = *rli as usize;
            let adj_len = if rl > cond_means.len() {
                (rl as f64) - tmean
            } else {
                (rl as f64) - cond_means[rl]
            };
            if adj_len >= 1.0 {
                adj_len
            } else {
                rl as f64
            }
        })
        .collect::<Vec<f64>>();

    el
}

#[inline]
fn m_step(
    eq_map: &AHashMap<EqLabel, usize>,
    prev_count: &[f64],
    eff_lens: &[f64],
    curr_counts: &mut [f64],
) {
    for (k, v) in eq_map {
        let mut denom = 0.0_f64;
        let count = *v as f64;
        for target_id in &k.targets {
            denom += prev_count[*target_id as usize] / eff_lens[*target_id as usize];
        }

        if denom > 1e-8 {
            for target_id in &k.targets {
                curr_counts[*target_id as usize] += count
                    * ((prev_count[*target_id as usize] / eff_lens[*target_id as usize]) / denom);
            }
        }
    }
}

fn em(eq_map: &AHashMap<EqLabel, usize>, eff_lens: &[f64], max_iter: u32) -> Vec<f64> {
    let total_weight: f64 = eq_map.values().sum::<usize>() as f64;

    // init
    let avg = total_weight / (eff_lens.len() as f64);
    let mut prev_counts = vec![avg; eff_lens.len()];
    let mut curr_counts = vec![0.0f64; eff_lens.len()];

    let mut rel_diff = 0.0_f64;
    let mut niter = 0_u32;

    while niter < max_iter {
        m_step(eq_map, &prev_counts, eff_lens, &mut curr_counts);

        //std::mem::swap(&)
        for i in 0..curr_counts.len() {
            if prev_counts[i] > 1e-8 {
                let rd = (curr_counts[i] - prev_counts[i]) / prev_counts[i];
                rel_diff = if rel_diff > rd { rel_diff } else { rd };
            }
        }

        std::mem::swap(&mut prev_counts, &mut curr_counts);
        curr_counts.fill(0.0_f64);

        if rel_diff < 1e-3 {
            break;
        }
        niter += 1;
        if niter % 100 == 0 {
            eprintln!("iteration {}; rel diff {}", niter, rel_diff);
        }
        rel_diff = 0.0_f64;
    }

    prev_counts.iter_mut().for_each(|x| if *x < 1e-8 { *x = 0.0} );
    m_step(eq_map, &prev_counts, eff_lens, &mut curr_counts);

    curr_counts
}

fn write_results(
    output: PathBuf,
    hdr: &rad_types::RadHeader,
    e_counts: &[f64],
    eff_lengths: &[f64],
) -> anyhow::Result<()> {
    let mut ofile = File::create(output)?;

    ofile.write_all("target_name\teelen\tecount\n".to_string().as_bytes())?;

    for (i, name) in hdr.ref_names.iter().enumerate() {
        let l = format!("{}\t{}\t{}\n", name, eff_lengths[i], e_counts[i]);
        ofile.write_all(l.as_bytes())?;
    }

    Ok(())
}

fn main() -> anyhow::Result<()> {
    let decorator = slog_term::TermDecorator::new().build();
    let drain = slog_term::CompactFormat::new(decorator)
        /*.use_custom_timestamp(|out: &mut dyn std::io::Write| {
            write!(out, "{}", chrono::Local::now().format("%Y-%m-%d %H:%M:%S")).unwrap();
            Ok(())
        })*/
        .build()
        .fuse();
    let drain = slog_async::Async::new(drain).build().fuse();

    let log = slog::Logger::root(drain, o!());

    let cli_args = Cli::parse();

    match cli_args.command {
        Commands::Quant { input, output, max_iter } => {
            info!(log, "path {:?}", input);
            let mut input_rad = input;
            input_rad.set_extension("rad");
            let i_file = File::open(input_rad).context("could not open input rad file")?;
            let mut br = BufReader::new(i_file);

            let hdr = rad_types::RadHeader::from_bytes(&mut br);
            info!(log, "read header!");
            // file-level
            let fl_tags = rad_types::TagSection::from_bytes(&mut br);
            info!(log, "read {:?} file-level tags", fl_tags.tags.len());

            // read-level
            let rl_tags = rad_types::TagSection::from_bytes(&mut br);
            info!(log, "read {:?} read-level tags", rl_tags.tags.len());

            const FRAG_TYPE_NAME: &str = "frag_map_type";
            let mut ftt: Option<u8> = None;
            // parse actual tags
            for rt in &rl_tags.tags {
                if rad_types::decode_int_type_tag(rt.typeid).is_none() {
                    crit!(
                        log,
                        "currently only RAD types 1--4 are supported for 'b' and 'u' tags."
                    );
                    std::process::exit(exit_codes::EXIT_UNSUPPORTED_TAG_TYPE);
                }
                info!(log, "\tread-level tag {}", rt.name);
                if rt.name == FRAG_TYPE_NAME {
                    ftt = Some(rt.typeid);
                }
            }

            // alignment-level
            let al_tags = rad_types::TagSection::from_bytes(&mut br);
            info!(log, "read {:?} alignemnt-level tags", al_tags.tags.len());

            const REF_ORI_NAME: &str = "compressed_ori_ref";
            const POS_NAME: &str = "pos";
            const FRAGLEN_NAME: &str = "frag_len";

            let mut ref_ori_t: Option<u8> = None;
            let mut pos_t: Option<u8> = None;
            let mut fraglen_t: Option<u8> = None;
            // parse actual tags
            for at in &al_tags.tags {
                if rad_types::decode_int_type_tag(at.typeid).is_none() {
                    crit!(
                        log,
                        "currently only RAD types 1--4 are supported for 'b' and 'u' tags."
                    );
                    std::process::exit(exit_codes::EXIT_UNSUPPORTED_TAG_TYPE);
                }
                info!(log, "\talignment-level tag {}", at.name);
                match at.name.as_str() {
                    REF_ORI_NAME => {
                        ref_ori_t = Some(at.typeid);
                    }
                    POS_NAME => {
                        pos_t = Some(at.typeid);
                    }
                    FRAGLEN_NAME => {
                        fraglen_t = Some(at.typeid);
                    }
                    _ => {
                        info!(log, "unknown alignment-level tag {}", at.name);
                    }
                }
            }
            // ensure the tags we expect are found
            assert!(ref_ori_t.is_some());
            assert!(pos_t.is_some());
            assert!(fraglen_t.is_some());

            // parse the file-level tags
            const REF_LENGTHS_NAME: &str = "ref_lengths";
            let mut ref_lengths: Option<Vec<u32>> = None;
            for ft in &fl_tags.tags {
                if ft.name == REF_LENGTHS_NAME {
                    if ft.typeid != 7 {
                        crit!(log, "expected array type (7) but found, {}", ft.typeid);
                    }
                    info!(log, "found frag_length file-level tag");

                    let v = read_ref_lengths(&mut br)?;
                    ref_lengths = Some(v);
                }
            }

            let ref_lengths =
                ref_lengths.expect("was not able to read reference lengths from file!");
            info!(log, "read {} reference lengths", ref_lengths.len());

            let frag_map_t = rad_types::decode_int_type_tag(
                ftt.expect("no fragment mapping type description present."),
            )
            .context("unknown fragment mapping type id.")?;
            let mut tot_mappings = 0usize;
            let mut num_mapped_reads = 0usize;
            let (eq_map, cond_means) = process(
                &mut br,
                hdr.num_chunks as usize,
                &frag_map_t,
                &mut tot_mappings,
                &mut num_mapped_reads,
            );

            let eff_lengths = adjust_ref_lengths(&ref_lengths, &cond_means);

            let em_res = em(&eq_map, &eff_lengths, max_iter);

            write_results(output, &hdr, &em_res, &eff_lengths)?;

            info!(log, "num mapped reads = {}", num_mapped_reads);
            info!(log, "total mappings = {}", tot_mappings);
            info!(log, "number of equivalence classes = {}", eq_map.len());
            let total_weight: usize = eq_map.values().sum();
            info!(log, "total equivalence map weight = {}", total_weight);
        }
    }
    Ok(())
}
