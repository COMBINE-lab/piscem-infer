use anyhow::{Context, bail};
use arrow2::{
    array::{Float64Array, UInt32Array},
    chunk::Chunk,
    datatypes::Field,
};

use indicatif::{HumanCount, ProgressBar, ProgressDrawTarget, ProgressStyle};
use num_format::{Locale, ToFormattedString};
use path_tools::WithAdditionalExtension;
use serde::Serialize;
use serde_json::{Value, json};
use std::io::{BufReader, Read};
use std::path::Path;
use std::time::Duration;
use std::{
    fs::{File, create_dir_all},
    io::Seek,
};
use tabled::{Table, Tabled, settings::Style};
use tracing::{info, warn};

use crate::utils::gibbs::do_gibbs;
use crate::utils::eq_maps::{
    BasicEqMap, EqLabel, EqMap, EqMapType, OrientationProperty, PackedEqMap, RangeFactorizedEqMap,
};
use crate::utils::io;
use crate::utils::map_record_types::{
    LibraryType, OrientationCounts, check_strand_warnings, detect_library_type,
};
use crate::{
    fld::FldPDF,
    prog_opts::{LibTypeArg, QuantOpts},
};
use crate::{
    fld::{EmpiricalFLD, Fld, ParametricFLD},
    utils::em::{
        EMInfo, adjust_ref_lengths, conditional_means, conditional_means_from_params, do_bootstrap,
        em, em_par,
    },
};

use libradicl::rad_types::{self, MappedFragmentOrientation};
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

fn detect_lib_type_from_sample<T: Read>(
    br: &mut BufReader<T>,
    nchunk: usize,
    record_context: &PiscemBulkRecordContext,
    paired_end: bool,
    max_samples: usize,
) -> anyhow::Result<LibraryType> {
    let mut counts = OrientationCounts::default();
    let mut sampled = 0usize;

    'sample: for _ in 0..nchunk {
        let c = chunk::Chunk::<PiscemBulkReadRecord>::from_bytes(br, record_context);
        for mappings in &c.reads {
            for o in &mappings.dirs {
                counts.add(*o);
            }
            sampled += 1;
            if sampled >= max_samples {
                break 'sample;
            }
        }
    }

    info!(
        "Auto-detection sampled {} reads: forward={}, reverse={}, FR={}, RF={}, FF={}, RR={}, unknown={}",
        sampled,
        counts.forward, counts.reverse,
        counts.forward_reverse, counts.reverse_forward,
        counts.forward_forward, counts.reverse_reverse,
        counts.unknown
    );

    if sampled == 0 {
        anyhow::bail!("No mapped reads found in sample for library type auto-detection");
    }

    let (detected, ratio) = detect_library_type(&counts, paired_end);
    check_strand_warnings(detected, ratio, paired_end);

    info!(
        "Auto-detected library type: {} (forward-strand ratio: {:.4})",
        detected, ratio
    );

    Ok(detected)
}

pub fn process_bulk(quant_opts: QuantOpts, eq_map_t: EqMapType) -> anyhow::Result<()> {
    let eqmap_orientation_status = OrientationProperty::OrientationAware;
    match eq_map_t {
        EqMapType::BasicEqMap => {
            process_bulk_dispatch(quant_opts, BasicEqMap::new(eqmap_orientation_status))
        }
        EqMapType::RangeFactorizedEqMap => process_bulk_dispatch(
            quant_opts,
            RangeFactorizedEqMap::new(eqmap_orientation_status),
        ),
    }
}

pub fn process_bulk_dispatch<EqLabelT: EqLabel>(
    quant_opts: QuantOpts,
    eqc_map: EqMap<EqLabelT>,
) -> anyhow::Result<()> {
    let qo = quant_opts.clone();
    let input = qo.input;
    let output = qo.output;
    let max_iter = qo.max_iter;
    let convergence_thresh = qo.convergence_thresh;
    let presence_thresh = qo.presence_thresh;
    let fld_mean = qo.fld_mean;
    let fld_sd = qo.fld_sd;
    let num_bootstraps = qo.num_bootstraps;
    let num_gibbs_samples = qo.num_gibbs_samples;
    let gibbs_thinning_factor = qo.gibbs_thinning_factor;
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

    // resolve library type (auto-detect if requested)
    let lib_type: LibraryType = match qo.lib_type {
        LibTypeArg::Explicit(lt) => {
            info!("Using user-specified library type: {}", lt);
            lt
        }
        LibTypeArg::Auto => {
            let file_offset = br.stream_position()?;
            let detected = detect_lib_type_from_sample(
                &mut br,
                prelude.hdr.num_chunks as usize,
                &tag_context,
                paired_end,
                qo.auto_detect_samples,
            )?;
            br.seek(std::io::SeekFrom::Start(file_offset))?;
            detected
        }
    };

    let mut frag_stats = MappedFragStats::new();
    let est_frag_lengths: Option<Vec<u32>> = if paired_end {
        // record the position in the file (right at the start)
        // of the first chunk
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

    let fld: Fld = if let Some(est_frag_lengths) = est_frag_lengths {
        Fld::Empirical(EmpiricalFLD::new(est_frag_lengths, f64::MIN_POSITIVE))
    } else {
        Fld::Parametric(ParametricFLD::new(fl_mean, fl_sd, 65_536_usize))
    };

    let (packed_eq_map, frag_lengths) = process(
        &mut br,
        prelude.hdr.num_chunks as usize,
        &tag_context,
        lib_type,
        &mut frag_stats,
        ref_lengths,
        eqc_map,
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

    if num_gibbs_samples > 0 {
        info!("performing Gibbs sampling ({num_gibbs_samples} samples, thinning factor {gibbs_thinning_factor})");
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()?;
        let gibbs_samples =
            do_gibbs(&eminfo, &em_res, num_gibbs_samples, gibbs_thinning_factor);

        let mut new_arrays = vec![];
        let mut gs_fields = vec![];
        for (i, g) in gibbs_samples.into_iter().enumerate() {
            let gs_array = Float64Array::from_vec(g);
            gs_fields.push(Field::new(
                format!("bootstrap.{i}"),
                gs_array.data_type().clone(),
                false,
            ));
            new_arrays.push(gs_array.boxed());
        }
        let chunk = Chunk::new(new_arrays);
        io::write_infrep_file(&output, gs_fields, chunk)?;
    }

    let infrep_method = if num_gibbs_samples > 0 {
        "gibbs"
    } else if num_bootstraps > 0 {
        "bootstrap"
    } else {
        "none"
    };

    let meta_info_output = output.with_additional_extension(".meta_info.json");
    let ofile = File::create(meta_info_output)?;
    let meta_info = json!({
        "quant_opts": quant_opts,
        "inferred_lib_type": lib_type.to_string(),
        "mapped_frag_stats": frag_stats,
        "num_bootstraps": num_bootstraps,
        "num_gibbs_samples": num_gibbs_samples,
        "infrep_method": infrep_method,
        "num_targets": eff_lengths.len(),
        "signatures": ref_sig_json
    });
    serde_json::to_writer_pretty(ofile, &meta_info)?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn process<T: Read, EqLabelT: EqLabel>(
    br: &mut BufReader<T>,
    nrec: usize,
    record_context: &PiscemBulkRecordContext,
    lib_type: LibraryType,
    mapped_stats: &mut MappedFragStats,
    ref_lengths: &[u32],
    eq_map: EqMap<EqLabelT>,
    fld_pdf: Fld,
) -> (PackedEqMap<EqLabelT>, Vec<u32>) {
    match fld_pdf {
        Fld::Empirical(f) => process_dispatch(
            br,
            nrec,
            record_context,
            lib_type,
            mapped_stats,
            ref_lengths,
            f,
            eq_map,
        ),
        Fld::Parametric(f) => process_dispatch(
            br,
            nrec,
            record_context,
            lib_type,
            mapped_stats,
            ref_lengths,
            f,
            eq_map,
        ),
    }
}

#[allow(clippy::too_many_arguments)]
fn process_dispatch<T: Read, D: FldPDF, EqLabelT: EqLabel>(
    br: &mut BufReader<T>,
    nrec: usize,
    record_context: &PiscemBulkRecordContext,
    lib_type: LibraryType,
    mapped_stats: &mut MappedFragStats,
    ref_lengths: &[u32],
    fld_pdf: D,
    mut eqmap: EqMap<EqLabelT>,
) -> (PackedEqMap<EqLabelT>, Vec<u32>) {
    /*
    let eqmap_orientation_status = if eqmap.contains_ori {
        OrientationProperty::OrientationAware
    } else {
        OrientationProperty::OrientationAgnostic
    };
    */

    //let map = &mut eqmap.count_map;
    let mut frag_lengths = vec![0u32; 65_536];
    const TARGET_UNIQUE_FRAGS: u32 = 5_000;
    let mut unique_frags = 0u32;

    let mut mapped_ori_count = [0u32; 7];
    let mut filtered_ori_count = [0u32; 7];

    let mut label_ints = vec![];
    let mut dir_ints = vec![];
    let mut probs = vec![];

    let pb = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr_with_hz(1));
    pb.set_style(
        ProgressStyle::with_template("{spinner:.green} Processed {human_pos} reads [{elapsed_precise}]")
            .unwrap()
            .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
    );
    pb.enable_steady_tick(Duration::from_secs(1));

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

            for (((r, pos), o), l) in mappings
                .refs
                .iter()
                .zip(mappings.positions.iter())
                .zip(mappings.dirs.iter())
                .zip(mappings.frag_lengths.iter())
            {
                let y = u32::from(*o);
                if lib_type.is_compatible_with(*o) {
                    mapped_ori_count[y as usize] += 1;
                    label_ints.push(*r);
                    dir_ints.push(y);
                    let frag_len_prob = match o {
                        MappedFragmentOrientation::Forward => {
                            let max_frag_len = (ref_lengths[*r as usize] - *pos) as usize;
                            fld_pdf.cdf(max_frag_len)
                        }
                        MappedFragmentOrientation::Reverse => {
                            // TODO: modify the bulk RAD creation upstream so that
                            // if a read is mapped as an orphan, then the length
                            // field encodes the read length (instead of the maximum
                            // fragment length).  This will allow us to compute the
                            // actual right-most end of the mapped fragment, rather than
                            // just assuming that the read is around length 100.
                            let max_frag_len = *pos as usize + 100;
                            fld_pdf.cdf(max_frag_len)
                        }
                        MappedFragmentOrientation::ForwardReverse
                        | MappedFragmentOrientation::ReverseForward => fld_pdf.pdf(*l as usize),
                        _ => 2.0 * f64::MIN_POSITIVE,
                    };
                    probs.push(frag_len_prob)
                } else {
                    filtered_ori_count[y as usize] += 1;
                }
            }

            label_ints.append(&mut dir_ints);
            let eql = EqLabelT::new(&label_ints, Some(&probs));
            eqmap.add(eql);

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

            pb.inc(1);
        }
    }

    pb.finish_with_message(format!(
        "Done — processed {} reads",
        HumanCount(mapped_stats.num_mapped_reads as u64)
    ));

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
