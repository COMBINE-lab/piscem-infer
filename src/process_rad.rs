use anyhow::{Context, bail};
use arrow2::{
    array::{Float64Array, UInt32Array},
    chunk::Chunk,
    datatypes::Field,
};

use num_format::{Locale, ToFormattedString};
use path_tools::WithAdditionalExtension;
use serde::Serialize;
use serde_json::{Value, json};
use std::io::{BufReader, Read};
use std::path::Path;
use std::{
    fs::{File, create_dir_all},
    io::Seek,
};
use tabled::{Table, Tabled, settings::Style};
use tracing::{info, warn};

use crate::utils::eq_maps::{
    BasicEqLabel, BasicEqMap, EqLabel, EqMap, EqMapType, OrientationProperty, PackedEqMap,
    RangeFactorizedEqMap,
};
use crate::utils::io;
use crate::utils::map_record_types::LibraryType;
use crate::{fld::FldPDF, prog_opts::QuantOpts};
use crate::{
    fld::{EmpiricalFLD, FLD, ParametricFLD},
    utils::em::{
        EMInfo, adjust_ref_lengths, conditional_means, conditional_means_from_params, do_bootstrap,
        em, em_par,
    },
};

use libradicl::rad_types::{self};
use libradicl::{
    chunk,
    header::RadPrelude,
    record::{PiscemBulkReadRecord, PiscemBulkRecordContext},
};

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
struct DirectionalEntry {
    name: &'static str,
    count: u32,
}

fn build_ori_table(mapped_ori_count_global: &[u32]) -> Vec<DirectionalEntry> {
    vec![
        DirectionalEntry {
            name: "unknown",
            count: mapped_ori_count_global[0],
        },
        DirectionalEntry {
            name: "f",
            count: mapped_ori_count_global[1],
        },
        DirectionalEntry {
            name: "r",
            count: mapped_ori_count_global[2],
        },
        DirectionalEntry {
            name: "fr",
            count: mapped_ori_count_global[3],
        },
        DirectionalEntry {
            name: "rf",
            count: mapped_ori_count_global[4],
        },
        DirectionalEntry {
            name: "ff",
            count: mapped_ori_count_global[5],
        },
        DirectionalEntry {
            name: "rr",
            count: mapped_ori_count_global[6],
        },
    ]
}

fn compute_fld_from_sample<T: Read>(
    br: &mut BufReader<T>,
    nchunk: usize,
    record_context: &PiscemBulkRecordContext,
    lib_type: LibraryType,
    mut param_est_frags: isize,
) -> anyhow::Result<Vec<u32>> {
    let mut temp_frag_lengths = vec![0u32; 65_536];
    let mut sufficient_samples = false;
    let requested_samples = param_est_frags;

    'estimate_fld: for _ in 0..nchunk {
        let c = chunk::Chunk::<PiscemBulkReadRecord>::from_bytes(br, record_context);
        for mappings in &c.reads {
            let ft = rad_types::MappingType::from_u8(mappings.frag_type);
            let nm = mappings.positions.len();
            if nm == 1 && !ft.is_orphan() {
                let o = mappings.dirs.first().expect("at least one mapping");
                if lib_type.is_compatible_with(*o) {
                    if let Some(fl) = mappings.frag_lengths.first() {
                        temp_frag_lengths[*fl as usize] += 1;
                        param_est_frags -= 1;
                        if param_est_frags <= 0 {
                            sufficient_samples = true;
                            break 'estimate_fld;
                        }
                    }
                }
            }
        }
    }

    if sufficient_samples {
        let cmeans = conditional_means(&temp_frag_lengths);
        info!(
            "computed conditional means ... last is {}",
            cmeans.last().expect("present")
        );
    } else {
        let nseen = requested_samples - param_est_frags;
        warn!(
            "insufficient uniquely mapped reads from which to estimate the fragment length distribution. {requested_samples} requested but only {nseen} were observed!"
        );
    }
    Ok(temp_frag_lengths)
}

#[allow(dead_code)]
fn compute_fld_from_params(mu: f64, sigma: f64, weight: usize, upper: usize) -> Vec<u32> {
    let inv_sigma = 1.0 / sigma;
    let denom_b = distrs::Normal::cdf(upper as f64, mu, sigma);
    let denom_a = distrs::Normal::cdf(0.0_f64, mu, sigma);
    let denom = denom_b - denom_a;
    let inv_denom = 1.0_f64 / denom;

    let trunc_pdf = |i: usize| -> f64 {
        let x = i as f64;
        inv_sigma * (distrs::Normal::pdf(x, mu, sigma) * inv_denom)
    };

    (0..upper)
        .map(|i| (trunc_pdf(i) * weight as f64).round() as u32)
        .collect()
}

pub fn process_bulk(quant_opts: QuantOpts) -> anyhow::Result<()> {
    let qo = quant_opts.clone();
    let input = qo.input;
    let lib_type = qo.lib_type;
    let output = qo.output;
    let max_iter = qo.max_iter;
    let convergence_thresh = qo.convergence_thresh;
    let presence_thresh = qo.presence_thresh;
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
            warn!(
                "Expected the mapping info file {input_map_info:?} to exist, but it doesn't. \
                    This is bad, and means that reference provenance signatures cannot be \
                    propagated to the output of piscem-infer. It is strongly recommended \
                    that you investigate why this file does not exist at the expected location."
            );
        } else {
            let map_info_str = std::fs::read_to_string(&input_map_info)
                .unwrap_or_else(|_| panic!("Couldn't open {:?}.", &input_map_info));
            let v: Value = serde_json::from_str(&map_info_str)?;
            if let Some(sigs) = v.get("signatures") {
                ref_sig_json = Some(sigs.clone());
            } else {
                warn!(
                    "The file {input_map_info:?} exists, but has no \"signatures\" entry holding the reference provenance signatures"
                );
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

    // read the header and tag sections from the rad file
    let prelude = RadPrelude::from_bytes(&mut br)?;

    info!("read header!");
    if prelude.hdr.is_paired > 0_u8 {
        info!("fragments paired in sequencing");
        paired_end = true;
        if let (Some(flm), Some(flsd)) = (fld_mean, fld_sd) {
            warn!(
                "provided fragment length distribution mean and sd ({flm}, {flsd}), but \
                    the RAD file contains paired-end fragments, so these will be ignored \"
                    and the fragment length distribution will be estimated."
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
                "The input RAD file {} was for unpaired reads, so \
                    a fragment length distribution mean and standard deviation \
                    must be provided.",
                &input_rad.display()
            );
        }
    }

    // file-level
    info!("read {:?} file-level tags", prelude.file_tags.tags.len());
    // parse actual tags
    for ft in &prelude.file_tags.tags {
        info!("\tfile-level tag {}", ft.name);
    }

    // read-level
    info!("read {:?} read-level tags", prelude.read_tags.tags.len());

    // required read-level tag
    const FRAG_TYPE_NAME: &str = "frag_map_type";
    let mut had_frag_map_type = false;
    // parse actual tags
    for rt in &prelude.read_tags.tags {
        info!("\tread-level tag {}", rt.name);
        if rt.name == FRAG_TYPE_NAME {
            had_frag_map_type = true;
        }
    }
    if !had_frag_map_type {
        bail!(
            "read-level tag description missing required tag \"{FRAG_TYPE_NAME}\"; can't proceed."
        );
    }

    // alignment-level
    info!(
        "read {:?} alignemnt-level tags",
        prelude.aln_tags.tags.len()
    );

    // required alignment level tags
    const REF_ORI_NAME: &str = "compressed_ori_ref";
    const POS_NAME: &str = "pos";
    const FRAGLEN_NAME: &str = "frag_len";

    let mut found_ref_ori_t = false;
    let mut found_pos_t = false;
    let mut found_fraglen_t = false;
    // parse actual tags
    for at in &prelude.aln_tags.tags {
        info!("\talignment-level tag {}", at.name);
        match at.name.as_str() {
            REF_ORI_NAME => {
                found_ref_ori_t = true;
            }
            POS_NAME => {
                found_pos_t = true;
            }
            FRAGLEN_NAME => {
                found_fraglen_t = true;
            }
            _ => {
                info!("unknown alignment-level tag {}", at.name);
            }
        }
    }
    // ensure the tags we expect are found
    assert!(
        found_ref_ori_t,
        "required alignment-level tag \"{REF_ORI_NAME}\" is missing"
    );
    assert!(
        found_pos_t,
        "required alignment-level tag \"{POS_NAME}\" is missing"
    );
    assert!(
        found_fraglen_t,
        "required alignment-level tag \"{FRAGLEN_NAME}\" is missing"
    );

    // parse the actual file-level tags
    const REF_LENGTHS_NAME: &str = "ref_lengths";
    let file_tag_map = prelude.file_tags.parse_tags_from_bytes(&mut br)?;
    // get the reference lengths from the tag map
    let ref_lengths = match file_tag_map.get(REF_LENGTHS_NAME) {
        Some(rad_types::TagValue::ArrayU32(v)) => Some(v),
        _ => None,
    };

    let ref_lengths = ref_lengths.expect("was not able to read reference lengths from file!");
    info!(
        "read {} reference lengths",
        ref_lengths.len().to_formatted_string(&Locale::en)
    );

    // extract whatever context we'll need to read the records
    let tag_context = prelude.get_record_context::<PiscemBulkRecordContext>()?;

    let mut frag_stats = MappedFragStats::new();
    let est_frag_lengths: Option<Vec<u32>> = if paired_end {
        let file_offset = br.stream_position()?;
        let temp_frag_lens = compute_fld_from_sample(
            &mut br,
            prelude.hdr.num_chunks as usize,
            &tag_context,
            lib_type,
            qo.param_est_frags,
        )?;
        // reset the stream to the start of the chunks
        br.seek(std::io::SeekFrom::Start(file_offset))?;
        Some(temp_frag_lens)
    } else {
        None
    };

    let fld: FLD = if let Some(est_frag_lengths) = est_frag_lengths {
        FLD::Empirical(EmpiricalFLD::new(est_frag_lengths, f64::MIN_POSITIVE))
    } else {
        FLD::Parametric(ParametricFLD::new(fl_mean, fl_sd, 65_536_usize))
    };

    let (packed_eq_map, frag_lengths) = process(
        &mut br,
        prelude.hdr.num_chunks as usize,
        &tag_context,
        lib_type,
        &mut frag_stats,
        EqMapType::RangeFactorizedEqMap,
        fld,
    );

    let cond_means = if paired_end {
        conditional_means(&frag_lengths)
    } else {
        conditional_means_from_params(fl_mean, fl_sd, 65_536_usize)
    };
    let eff_lengths = adjust_ref_lengths(ref_lengths, &cond_means);

    let eminfo = EMInfo {
        eq_map: &packed_eq_map,
        eff_lens: &eff_lengths,
        max_iter,
        convergence_thresh,
        presence_thresh,
    };

    let em_res = if num_threads > 1 {
        em_par(&eminfo, num_threads)
    } else {
        em(&eminfo)
    };

    let quant_output = output.with_additional_extension(".quant");
    io::write_results(
        &quant_output,
        &prelude.hdr,
        &em_res,
        ref_lengths,
        &eff_lengths,
    )?;

    info!(
        "num mapped reads = {}",
        frag_stats.num_mapped_reads.to_formatted_string(&Locale::en)
    );
    info!(
        "total mappings = {}",
        frag_stats.tot_mappings.to_formatted_string(&Locale::en)
    );
    info!(
        "number of equivalence classes = {}",
        packed_eq_map.len().to_formatted_string(&Locale::en)
    );
    let total_weight: usize = packed_eq_map.total_weight();
    info!(
        "total equivalence map weight = {}",
        total_weight.to_formatted_string(&Locale::en)
    );

    {
        let fld_array = UInt32Array::from_vec(frag_lengths);
        let field = Field::new("fragment_length_dist", fld_array.data_type().clone(), false);

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
                format!("bootstrap.{i}"),
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
    Ok(())
}

fn process<T: Read>(
    br: &mut BufReader<T>,
    nrec: usize,
    record_context: &PiscemBulkRecordContext,
    lib_type: LibraryType,
    mapped_stats: &mut MappedFragStats,
    eq_map_t: EqMapType,
    fld_pdf: FLD,
) -> (PackedEqMap, Vec<u32>) {
    let eqmap_orientation_status = OrientationProperty::OrientationAware;
    match (fld_pdf, eq_map_t) {
        (FLD::Empirical(f), EqMapType::BasicEqMap) => process_dispatch(
            br,
            nrec,
            record_context,
            lib_type,
            mapped_stats,
            f,
            BasicEqMap::new(eqmap_orientation_status),
        ),
        (FLD::Empirical(f), EqMapType::RangeFactorizedEqMap) => process_dispatch(
            br,
            nrec,
            record_context,
            lib_type,
            mapped_stats,
            f,
            RangeFactorizedEqMap::new(eqmap_orientation_status),
        ),
        (FLD::Parametric(f), EqMapType::BasicEqMap) => process_dispatch(
            br,
            nrec,
            record_context,
            lib_type,
            mapped_stats,
            f,
            BasicEqMap::new(eqmap_orientation_status),
        ),
        (FLD::Parametric(f), EqMapType::RangeFactorizedEqMap) => process_dispatch(
            br,
            nrec,
            record_context,
            lib_type,
            mapped_stats,
            f,
            RangeFactorizedEqMap::new(eqmap_orientation_status),
        ),
    }
}

fn process_dispatch<T: Read, D: FldPDF, EqLabelT: EqLabel>(
    br: &mut BufReader<T>,
    nrec: usize,
    record_context: &PiscemBulkRecordContext,
    lib_type: LibraryType,
    mapped_stats: &mut MappedFragStats,
    fld_pdf: D,
    mut eqmap: EqMap<EqLabelT>,
) -> (PackedEqMap, Vec<u32>) {
    let eqmap_orientation_status = if eqmap.contains_ori {
        OrientationProperty::OrientationAware
    } else {
        OrientationProperty::OrientationAgnostic
    };

    let map = &mut eqmap.count_map;
    let mut frag_lengths = vec![0u32; 65_536];
    const TARGET_UNIQUE_FRAGS: u32 = 5_000;
    let mut unique_frags = 0u32;

    let mut mapped_ori_count = [0u32; 7];
    let mut filtered_ori_count = [0u32; 7];

    let mut label_ints = vec![];
    let mut dir_ints = vec![];
    let mut probs = vec![];
    //let mut dir_vec = vec![0u32, 64];
    for _ in 0..nrec {
        let c = chunk::Chunk::<PiscemBulkReadRecord>::from_bytes(br, record_context);
        for mappings in &c.reads {
            let ft = rad_types::MappingType::from_u8(mappings.frag_type);
            let nm = mappings.positions.len();

            mapped_stats.tot_mappings += nm;
            mapped_stats.num_mapped_reads += 1;

            // reset the counter
            mapped_ori_count.fill(0);
            filtered_ori_count.fill(0);

            label_ints.clear();
            dir_ints.clear();
            probs.clear();

            for ((r, o), l) in mappings
                .refs
                .iter()
                .zip(mappings.dirs.iter())
                .zip(mappings.frag_lengths.iter())
            {
                let y = u32::from(*o);
                if lib_type.is_compatible_with(*o) {
                    mapped_ori_count[y as usize] += 1;
                    label_ints.push(*r);
                    dir_ints.push(y);
                    probs.push(fld_pdf.pdf(*l as usize))
                } else {
                    filtered_ori_count[y as usize] += 1;
                }
            }

            label_ints.append(&mut dir_ints);
            let eql = EqLabelT::new(&label_ints, Some(&probs));

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
        warn!(
            "Only observed {} uniquely-mapped fragments (< threshold of {}), the fragment length distribution estimate may not be robust",
            unique_frags, TARGET_UNIQUE_FRAGS
        );
    }

    let packed_eq_map = PackedEqMap::from_eq_map(&eqmap);
    drop(eqmap);

    (packed_eq_map, frag_lengths)
}
