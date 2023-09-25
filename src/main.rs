use anyhow::Context;
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

use crate::utils::em::OrientationProperty;

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
    tot_mappings: &mut usize,
    num_mapped_reads: &mut usize,
) -> (EqMap, Vec<f64>) {
    // true because it contains orientations
    let mut eqmap = EqMap::new(OrientationProperty::OrientationAware);

    let map = &mut eqmap.count_map;
    let mut frag_lengths = vec![0u32; 65_536];
    const TARGET_UNIQUE_FRAGS: u32 = 5_000;
    let mut unique_frags = 0u32;

    let mut mapped_ori_count_global = vec![0u32; 7];
    let mut mapped_ori_count = vec![0u32; 7];
    //let mut dir_vec = vec![0u32, 64];
    for _ in 0..nrec {
        let c = MetaChunk::from_bytes(br, frag_map_t);
        for mappings in &c.reads {
            let nm = mappings.positions.len();

            *tot_mappings += nm;
            *num_mapped_reads += 1;

            let mut label_ints = mappings.refs.clone();

            // reset the counter
            mapped_ori_count.fill(0);
            // extend with the info on the mapping
            // orientations
            label_ints.extend(mappings.dirs.iter().map(|x| {
                let y = u32::from(*x);
                mapped_ori_count[y as usize] += 1;
                y
            }));

            let eql = EqLabel {
                targets: label_ints,
            };

            map.entry(eql)
                .and_modify(|counter| *counter += 1)
                .or_insert(1);

            if nm == 1 {
                frag_lengths[mappings.frag_lengths[0] as usize] += 1;
                unique_frags += 1;
            }

            // update global orientations
            for i in 0..mapped_ori_count.len() {
                mapped_ori_count_global[i] += if mapped_ori_count[i] > 0 { 1 } else { 0 };
            }
        }
    }

    let count_table = build_ori_table(&mapped_ori_count_global);

    info!(
        "\n{}",
        Table::new(count_table).with(Style::rounded()).to_string()
    );

    if unique_frags < TARGET_UNIQUE_FRAGS {
        warn!("Only observed {} uniquely-mapped fragments (< threshold of {}), the fragment length distribution estimate may not be robust",
            unique_frags, TARGET_UNIQUE_FRAGS);
    }
    let cond_means = conditional_means(&frag_lengths);
    (eqmap, cond_means)
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
        let l = format!("{}\t{}\t{}\n", name, eff_lengths[i], e_counts[i]);
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
        /// where output should be written
        #[arg(short, long)]
        output: PathBuf,
        /// max iterations to run the EM
        #[arg(short, long, default_value_t = 1500)]
        max_iter: u32,
        /// convergence threshold for EM
        #[arg(long, default_value_t = 1e-3)]
        convergence_thresh: f64,
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
            output,
            max_iter,
            convergence_thresh,
        } => {
            info!("path {:?}", input);
            let mut input_rad = input;
            input_rad.set_extension("rad");
            let i_file = File::open(input_rad).context("could not open input rad file")?;
            let mut br = BufReader::new(i_file);

            let hdr = rad_types::RadHeader::from_bytes(&mut br);
            info!("read header!");
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

            let eminfo = EMInfo {
                eq_map: &eq_map,
                eff_lens: &eff_lengths,
                max_iter,
                convergence_thresh,
            };
            let em_res = em(&eminfo);

            write_results(output, &hdr, &em_res, &eff_lengths)?;

            info!("num mapped reads = {}", num_mapped_reads);
            info!("total mappings = {}", tot_mappings);
            info!("number of equivalence classes = {}", eq_map.len());
            let total_weight: usize = eq_map.count_map.values().sum();
            info!("total equivalence map weight = {}", total_weight);
        }
    }
    Ok(())
}
