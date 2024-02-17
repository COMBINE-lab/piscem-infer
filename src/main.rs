use anyhow::{bail, Context};
use arrow2::{
    array::{Float64Array, UInt32Array},
    chunk::Chunk,
    datatypes::Field,
};

use clap::Args;
use clap::{Parser, Subcommand};
use libradicl::exit_codes;
use libradicl::rad_types;
use path_tools::WithAdditionalExtension;
use scroll::Pread;
use serde::Serialize;
use serde_json::{json, Value};
use std::fs::{create_dir_all, File};
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};
use tabled::{settings::Style, Table, Tabled};
use tracing::{error, info, warn, Level};

mod utils;
use utils::custom_rad_utils::MetaChunk;
use utils::em::{adjust_ref_lengths, conditional_means, em, em_par, EMInfo, EqLabel, EqMap};

use crate::utils::em::conditional_means_from_params;
use crate::utils::em::do_bootstrap;
use crate::utils::em::OrientationProperty;
use crate::utils::em::PackedEqMap;
use crate::utils::io;
use utils::map_record_types::{LibraryType, MappingType};

#[derive(Serialize)]
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

    let mut mapped_ori_count = [0u32; 7];
    let mut filtered_ori_count = [0u32; 7];

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
fn read_ref_lengths<T: Read>(reader: &mut T) -> anyhow::Result<Vec<u32>> {
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

#[derive(Args, Serialize, Clone, Debug)]
pub struct QuantOpts {
    /// input stem (i.e. without the .rad suffix)
    #[arg(short, long)]
    pub input: PathBuf,
    /// the expected library type
    #[arg(short, long, value_parser = clap::value_parser!(LibraryType))]
    pub lib_type: LibraryType,
    /// output file prefix (multiple output files may be created, the main will have a `.quant` suffix)
    #[arg(short, long)]
    pub output: PathBuf,
    /// max iterations to run the EM
    #[arg(short, long, default_value_t = 1500)]
    pub max_iter: u32,
    /// convergence threshold for EM
    #[arg(long, default_value_t = 1e-3)]
    pub convergence_thresh: f64,
    /// mean of fragment length distribution mean
    /// (required, and used, only in the case of unpaired fragments).
    #[arg(long, requires = "fld_sd")]
    pub fld_mean: Option<f64>,
    /// mean of fragment length distribution standard deviation
    /// (required, and used, only in the case of unpaired fragments).
    #[arg(long, requires = "fld_mean")]
    pub fld_sd: Option<f64>,
    /// number of bootstrap replicates to perform.
    #[arg(long, default_value_t = 0)]
    pub num_bootstraps: usize,
    /// number of threads to use (used during the EM and for bootstrapping)
    #[arg(long, default_value_t = 16)]
    pub num_threads: usize,
}

#[derive(Debug, Subcommand)]
enum Commands {
    /// quantify from the rad file
    #[command(arg_required_else_help = true)]
    Quant(QuantOpts),
}

/// quantify target abundance from bulk-sequencing data
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(propagate_version = true)]
struct Cli {
    #[arg(short, long)]
    quiet: bool,
    #[command(subcommand)]
    command: Commands,
}

/// https://github.com/polymonster/hotline/blob/05afabb353355e22afea3884bdf3da92bd572da4/src/gfx.rs#L1384
/// Take any sized silce and convert to a slice of u8
pub fn as_u8_slice<T: Sized>(p: &[T]) -> &[u8] {
    let (prefix, body, suffix) = unsafe { p.align_to::<u8>() };
    assert!(prefix.is_empty());
    assert!(suffix.is_empty());
    body
}

fn main() -> anyhow::Result<()> {
    let cli_args = Cli::parse();
    if cli_args.quiet {
        tracing_subscriber::fmt().with_max_level(Level::WARN).init();
    } else {
        tracing_subscriber::fmt().with_max_level(Level::INFO).init();
    }

    match cli_args.command {
        Commands::Quant(quant_opts) => {
            let qo = quant_opts.clone();
            let input = qo.input;
            let lib_type = qo.lib_type;
            let output = qo.output;
            let max_iter = qo.max_iter;
            let convergence_thresh = qo.convergence_thresh;
            let fld_mean = qo.fld_mean;
            let fld_sd = qo.fld_sd;
            let num_bootstraps = qo.num_bootstraps;
            let num_threads = qo.num_threads;

            // if there is a parent directory
            if let Some(p) = output.parent() {
                // unless this was a relative path with one component,
                // which we should treat as the file prefix, then grab
                // the non-empty parent and create it.
                if p != Path::new("") {
                    create_dir_all(p)?;
                }
            }

            info!("path {:?}", input);
            let ref_sig_json;
            {
                let mut input_map_info = input.clone();
                input_map_info.set_extension("map_info.json");
                if !input_map_info.exists() {
                    ref_sig_json = None;
                    warn!(concat!("Expected the mapping info file {:?} to exist, but it doesn't. ",
                          "This is bad, and means that reference provenance signatures cannot be ",
                          "propagated to the output of piscem-infer. It is strongly recommended ",
                          "that you investigate why this file does not exist at the expected location."),
                        input_map_info);
                } else {
                    let map_info_str = std::fs::read_to_string(&input_map_info)
                        .unwrap_or_else(|_| panic!("Couldn't open {:?}.", &input_map_info));
                    let v: Value = serde_json::from_str(&map_info_str)?;
                    if let Some(sigs) = v.get("signatures") {
                        ref_sig_json = Some(sigs.clone());
                    } else {
                        warn!("The file {:?} exists, but has no \"signatures\" entry holding the reference provenance signatures",
                            input_map_info);
                        ref_sig_json = None;
                    }
                }
            }

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
                    warn!(
                        concat!(
                            "provided fragment length distribution mean and sd ({}, {}), but ",
                            "the RAD file contains paired-end fragments, so these will be ignored ",
                            "and the fragment length distribution will be estimated."
                        ),
                        flm, flsd
                    );
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

            let packed_eq_map = PackedEqMap::from_eq_map(&eq_map);

            let eminfo = EMInfo {
                eq_map: &packed_eq_map,
                eff_lens: &eff_lengths,
                max_iter,
                convergence_thresh,
            };

            let em_res = if num_threads > 1 {
                em_par(&eminfo, num_threads)
            } else {
                em(&eminfo)
            };

            let quant_output = output.with_additional_extension(".quant");
            io::write_results(&quant_output, &hdr, &em_res, &ref_lengths, &eff_lengths)?;

            info!("num mapped reads = {}", frag_stats.num_mapped_reads);
            info!("total mappings = {}", frag_stats.tot_mappings);
            info!("number of equivalence classes = {}", eq_map.len());
            let total_weight: usize = eq_map.count_map.values().sum();
            info!("total equivalence map weight = {}", total_weight);

            {
                let fld_array = UInt32Array::from_vec(frag_lengths);
                let field =
                    Field::new("fragment_length_dist", fld_array.data_type().clone(), false);

                let chunk = Chunk::new(vec![fld_array.boxed()]);
                let fields = vec![field];
                io::write_fld_file(&output, fields, chunk)?;
            }

            if num_bootstraps > 0 {
                info!("performing bootstraps");
                rayon::ThreadPoolBuilder::new()
                    .num_threads(num_threads)
                    .build_global()?;
                let bootstraps = do_bootstrap(&eminfo, num_bootstraps);

                let mut new_arrays = vec![];
                let mut bs_fields = vec![];
                for (i, b) in bootstraps.into_iter().enumerate() {
                    let bs_array = Float64Array::from_vec(b);
                    bs_fields.push(Field::new(
                        format!("bootstrap.{}", i),
                        bs_array.data_type().clone(),
                        false,
                    ));
                    new_arrays.push(bs_array.boxed());
                }
                let chunk = Chunk::new(new_arrays);
                io::write_infrep_file(&output, bs_fields, chunk)?;
            }

            let meta_info_output = output.with_additional_extension(".meta_info.json");
            let ofile = File::create(meta_info_output)?;
            let meta_info = json!({
                "quant_opts": quant_opts,
                "mapped_frag_stats": frag_stats,
                "num_bootstraps": num_bootstraps,
                "num_targets": eff_lengths.len(),
                "signatures": ref_sig_json
            });
            serde_json::to_writer_pretty(ofile, &meta_info)?;
        }
    }
    Ok(())
}
