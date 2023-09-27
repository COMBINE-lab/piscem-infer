use anyhow::{bail, Context};
use clap::{Parser, Subcommand};
use libradicl::exit_codes;
use libradicl::rad_types;
use scroll::Pread;
use std::fs::File;
use std::io::Write;
use std::io::{BufReader, Read};
use std::path::PathBuf;
use tabled::{settings::Style, Table, Tabled};
use tracing::{error, info, warn, Level};

mod utils;
use utils::custom_rad_utils::MetaChunk;
use utils::em::{adjust_ref_lengths, conditional_means, em, EMInfo, EqLabel, EqMap};

use crate::utils::em::conditional_means_from_params;
use crate::utils::em::OrientationProperty;
use utils::map_record_types::{LibraryType, MappingType};

struct MappedFragStats {
    tot_mappings: usize,
    num_mapped_reads: usize,
    mapped_ori_count: [u32; 7],
    filtered_ori_count: [u32; 7],
}

impl MappedFragStats {
    pub fn new() -> Self {
        Self {
            tot_mappings: 0,
            num_mapped_reads: 0,
            mapped_ori_count: [0u32; 7],
            filtered_ori_count: [0u32; 7],
        }
    }
}

#[derive(Tabled)]
struct DirEntry {
    name: &'static str,
    count: u32,
}

fn build_ori_table(mapped_ori_count_global: &[u32]) -> Vec<DirEntry> {
    vec![
        DirEntry {
            name: "unknown",
            count: mapped_ori_count_global[0],
        },
        DirEntry {
            name: "f",
            count: mapped_ori_count_global[1],
        },
        DirEntry {
            name: "r",
            count: mapped_ori_count_global[2],
        },
        DirEntry {
            name: "fr",
            count: mapped_ori_count_global[3],
        },
        DirEntry {
            name: "rf",
            count: mapped_ori_count_global[4],
        },
        DirEntry {
            name: "ff",
            count: mapped_ori_count_global[5],
        },
        DirEntry {
            name: "rr",
            count: mapped_ori_count_global[6],
        },
    ]
}

fn process<T: Read>(
    br: &mut BufReader<T>,
    nrec: usize,
    frag_map_t: &rad_types::RadIntId,
    lib_type: LibraryType,
    mapped_stats: &mut MappedFragStats,
) -> (EqMap, Vec<u32>) {
    // true because it contains orientations
    let mut eqmap = EqMap::new(OrientationProperty::OrientationAware);

    let map = &mut eqmap.count_map;
    let mut frag_lengths = vec![0u32; 65_536];
    const TARGET_UNIQUE_FRAGS: u32 = 5_000;
    let mut unique_frags = 0u32;

    let mut mapped_ori_count = vec![0u32; 7];
    let mut filtered_ori_count = vec![0u32; 7];

    let mut label_ints = vec![];
    let mut dir_ints = vec![];
    //let mut dir_vec = vec![0u32, 64];
    for _ in 0..nrec {
        let c = MetaChunk::from_bytes(br, frag_map_t);
        for mappings in &c.reads {
            let ft = MappingType::from_u8(mappings.frag_type);
            let nm = mappings.positions.len();

            mapped_stats.tot_mappings += nm;
            mapped_stats.num_mapped_reads += 1;

            // reset the counter
            mapped_ori_count.fill(0);
            filtered_ori_count.fill(0);

            label_ints.clear();
            dir_ints.clear();

            for (r, o) in mappings.refs.iter().zip(mappings.dirs.iter()) {
                let y = u32::from(*o);
                if lib_type.is_compatible_with(*o) {
                    mapped_ori_count[y as usize] += 1;
                    label_ints.push(*r);
                    dir_ints.push(y);
                } else {
                    filtered_ori_count[y as usize] += 1;
                }
            }

            label_ints.append(&mut dir_ints);

            let eql = EqLabel {
                targets: label_ints.clone(),
            };

            map.entry(eql)
                .and_modify(|counter| *counter += 1)
                .or_insert(1);

            if nm == 1 && !ft.is_orphan() {
                if let Some(fl) = mappings.frag_lengths.first() {
                    frag_lengths[*fl as usize] += 1;
                }
                unique_frags += 1;
            }

            // update global orientations
            for i in 0..mapped_ori_count.len() {
                mapped_stats.mapped_ori_count[i] += if mapped_ori_count[i] > 0 { 1 } else { 0 };
                mapped_stats.filtered_ori_count[i] += if filtered_ori_count[i] > 0 { 1 } else { 0 };
            }
        }
    }

    let count_table_pass = build_ori_table(&mapped_stats.mapped_ori_count);
    info!(
        "mapping counts passing filtering\n{}\n",
        Table::new(count_table_pass)
            .with(Style::rounded())
            .to_string()
    );

    let count_table_filter = build_ori_table(&mapped_stats.filtered_ori_count);
    info!(
        "mapping counts failing filtering\n{}\n",
        Table::new(count_table_filter)
            .with(Style::rounded())
            .to_string()
    );

    if unique_frags < TARGET_UNIQUE_FRAGS {
        warn!("Only observed {} uniquely-mapped fragments (< threshold of {}), the fragment length distribution estimate may not be robust",
            unique_frags, TARGET_UNIQUE_FRAGS);
    }

    (eqmap, frag_lengths)
}

/// Read the lengths of the reference sequences from the RAD file header.
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

fn write_results(
    output: PathBuf,
    hdr: &rad_types::RadHeader,
    e_counts: &[f64],
    eff_lengths: &[f64],
) -> anyhow::Result<()> {
    let mut ofile = File::create(output)?;

    ofile.write_all("target_name\teelen\tecount\n".to_string().as_bytes())?;

    for (i, name) in hdr.ref_names.iter().enumerate() {
        let l = format!("{}\t{:.3}\t{:.3}\n", name, eff_lengths[i], e_counts[i]);
        ofile.write_all(l.as_bytes())?;
    }

    Ok(())
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// quantify from the rad file
    #[command(arg_required_else_help = true)]
    Quant {
        /// input stem (i.e. without the .rad suffix)
        #[arg(short, long)]
        input: PathBuf,
        /// the expected library type
        #[arg(short, long, value_parser = clap::value_parser!(LibraryType))]
        lib_type: LibraryType,
        /// where output should be written
        #[arg(short, long)]
        output: PathBuf,
        /// max iterations to run the EM
        #[arg(short, long, default_value_t = 1500)]
        max_iter: u32,
        /// convergence threshold for EM
        #[arg(long, default_value_t = 1e-3)]
        convergence_thresh: f64,
        /// mean of fragment length distribution mean
        /// (required, and used, only in the case of unpaired fragments).
        #[arg(long, requires = "fld_sd")]
        fld_mean: Option<f64>,
        /// mean of fragment length distribution standard deviation
        /// (required, and used, only in the case of unpaired fragments).
        #[arg(long, requires = "fld_mean")]
        fld_sd: Option<f64>,
    },
}

/// quantifying from a metagenomic rad file
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(propagate_version = true)]
struct Cli {
    #[arg(short, long)]
    quiet: bool,
    #[command(subcommand)]
    command: Commands,
}

fn main() -> anyhow::Result<()> {
    let cli_args = Cli::parse();
    if cli_args.quiet {
        tracing_subscriber::fmt().with_max_level(Level::WARN).init();
    } else {
        tracing_subscriber::fmt().with_max_level(Level::INFO).init();
    }

    match cli_args.command {
        Commands::Quant {
            input,
            lib_type,
            output,
            max_iter,
            convergence_thresh,
            fld_mean,
            fld_sd,
        } => {
            info!("path {:?}", input);
            let mut input_rad = input;
            input_rad.set_extension("rad");
            let i_file = File::open(&input_rad).context("could not open input rad file")?;
            let mut br = BufReader::new(i_file);

            let paired_end: bool;
            let mut fl_mean = 0_f64;
            let mut fl_sd = 0_f64;

            let hdr = rad_types::RadHeader::from_bytes(&mut br);
            info!("read header!");
            if hdr.is_paired > 0_u8 {
                info!("fragments paired in sequencing");
                paired_end = true;
                if let (Some(flm), Some(flsd)) = (fld_mean, fld_sd) {
                    warn!(concat!("provided fragment length distribution mean and sd ({}, {}), but ",
                                  "the RAD file contains paired-end fragments, so these will be ignored ",
                                  "and the fragment length distribution will be estimated."), flm, flsd);
                }
            } else {
                info!("fragments unpaired in sequencing");
                paired_end = false;
                if let (Some(flm), Some(flsd)) = (fld_mean, fld_sd) {
                    fl_mean = flm;
                    fl_sd = flsd;
                } else {
                    bail!(
                        concat!(
                            "The input RAD file {} was for unpaired reads, so ",
                            "a fragment length distribution mean and standard deviation ",
                            "must be provided."
                        ),
                        &input_rad.display()
                    );
                }
            }

            // file-level
            let fl_tags = rad_types::TagSection::from_bytes(&mut br);
            info!("read {:?} file-level tags", fl_tags.tags.len());

            // read-level
            let rl_tags = rad_types::TagSection::from_bytes(&mut br);
            info!("read {:?} read-level tags", rl_tags.tags.len());

            const FRAG_TYPE_NAME: &str = "frag_map_type";
            let mut ftt: Option<u8> = None;
            // parse actual tags
            for rt in &rl_tags.tags {
                if rad_types::decode_int_type_tag(rt.typeid).is_none() {
                    error!("currently only RAD types 1--4 are supported for 'b' and 'u' tags.");
                    std::process::exit(exit_codes::EXIT_UNSUPPORTED_TAG_TYPE);
                }
                info!("\tread-level tag {}", rt.name);
                if rt.name == FRAG_TYPE_NAME {
                    ftt = Some(rt.typeid);
                }
            }

            // alignment-level
            let al_tags = rad_types::TagSection::from_bytes(&mut br);
            info!("read {:?} alignemnt-level tags", al_tags.tags.len());

            const REF_ORI_NAME: &str = "compressed_ori_ref";
            const POS_NAME: &str = "pos";
            const FRAGLEN_NAME: &str = "frag_len";

            let mut ref_ori_t: Option<u8> = None;
            let mut pos_t: Option<u8> = None;
            let mut fraglen_t: Option<u8> = None;
            // parse actual tags
            for at in &al_tags.tags {
                if rad_types::decode_int_type_tag(at.typeid).is_none() {
                    error!("currently only RAD types 1--4 are supported for 'b' and 'u' tags.");
                    std::process::exit(exit_codes::EXIT_UNSUPPORTED_TAG_TYPE);
                }
                info!("\talignment-level tag {}", at.name);
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
                        info!("unknown alignment-level tag {}", at.name);
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
                        error!("expected array type (7) but found, {}", ft.typeid);
                    }
                    info!("found frag_length file-level tag");

                    let v = read_ref_lengths(&mut br)?;
                    ref_lengths = Some(v);
                }
            }

            let ref_lengths =
                ref_lengths.expect("was not able to read reference lengths from file!");
            info!("read {} reference lengths", ref_lengths.len());

            let frag_map_t = rad_types::decode_int_type_tag(
                ftt.expect("no fragment mapping type description present."),
            )
            .context("unknown fragment mapping type id.")?;

            let mut frag_stats = MappedFragStats::new();
            let (eq_map, frag_lengths) = process(
                &mut br,
                hdr.num_chunks as usize,
                &frag_map_t,
                lib_type,
                &mut frag_stats,
            );

            let cond_means = if paired_end {
                conditional_means(&frag_lengths)
            } else {
                conditional_means_from_params(fl_mean, fl_sd, 65_636_usize)
            };
            let eff_lengths = adjust_ref_lengths(&ref_lengths, &cond_means);

            let eminfo = EMInfo {
                eq_map: &eq_map,
                eff_lens: &eff_lengths,
                max_iter,
                convergence_thresh,
            };
            let em_res = em(&eminfo);

            write_results(output, &hdr, &em_res, &eff_lengths)?;

            info!("num mapped reads = {}", frag_stats.num_mapped_reads);
            info!("total mappings = {}", frag_stats.tot_mappings);
            info!("number of equivalence classes = {}", eq_map.len());
            let total_weight: usize = eq_map.count_map.values().sum();
            info!("total equivalence map weight = {}", total_weight);
        }
    }
    Ok(())
}
