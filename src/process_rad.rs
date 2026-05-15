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
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::time::Duration;
use std::{
    fs::{File, create_dir_all},
    io::Seek,
};
use tabled::{Table, Tabled, settings::Style};
use tracing::{info, warn};

use crate::utils::collapsed_eq::build_collapsed;
use crate::utils::eq_maps::{
    BasicEqMap, EqLabel, EqMap, EqMapType, OrientationProperty, PackedEqMap,
    RangeFactorizedEqLabel, RangeFactorizedEqMap,
};
use crate::utils::gibbs::{do_gibbs, do_gibbs_with_pool};
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
        do_bootstrap_with_pool, em, em_par_with_pool, squarem_em, squarem_em_par_with_pool,
    },
};

/// Options for RAD file processing (shared between single-sample and multi-sample modes).
pub struct RadProcessingOpts {
    pub input: std::path::PathBuf,
    pub lib_type: LibTypeArg,
    pub param_est_frags: isize,
    pub fld_mean: Option<f64>,
    pub fld_sd: Option<f64>,
    pub auto_detect_samples: usize,
    /// Number of threads for parallel EQ map building (1 = sequential).
    pub num_threads: usize,
}

/// Bundle of results from building an EQ map from a RAD file.
pub struct EqMapBundle<EqLabelT: EqLabel> {
    pub packed_eq_map: PackedEqMap<EqLabelT>,
    pub eff_lengths: Vec<f64>,
    pub ref_lengths: Vec<u32>,
    pub ref_names: Vec<String>,
    pub frag_lengths: Vec<u32>,
    pub frag_stats: MappedFragStats,
    pub ref_sig_json: Option<Value>,
    pub lib_type: LibraryType,
}

use libradicl::rad_types::{self, MappedFragmentOrientation, TagMap};
use libradicl::{
    chunk,
    header::RadPrelude,
    readers::ParallelRadReader,
    record::{PiscemBulkReadRecord, PiscemBulkRecordContext},
};

#[derive(Serialize, Default)]
pub struct MappedFragStats {
    pub tot_mappings: usize,
    pub num_mapped_reads: usize,
    pub mapped_ori_count: [u32; 7],
    pub filtered_ori_count: [u32; 7],
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

fn read_transcript_mask(path: &Path) -> anyhow::Result<std::collections::HashSet<String>> {
    let file = File::open(path)
        .with_context(|| format!("failed to open transcript mask file: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut keep = std::collections::HashSet::new();
    for line in reader.lines() {
        let line = line?;
        let txp = line.trim();
        if txp.is_empty() || txp.starts_with('#') {
            continue;
        }
        keep.insert(txp.to_string());
    }
    Ok(keep)
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
                if lib_type.is_compatible_with(*o)
                    && let Some(fl) = mappings.frag_lengths.first()
                {
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
        counts.forward,
        counts.reverse,
        counts.forward_reverse,
        counts.reverse_forward,
        counts.forward_forward,
        counts.reverse_reverse,
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

fn detect_lib_type_and_fld_from_sample<T: Read>(
    br: &mut BufReader<T>,
    nchunk: usize,
    record_context: &PiscemBulkRecordContext,
    paired_end: bool,
    max_samples: usize,
    mut param_est_frags: isize,
) -> anyhow::Result<(LibraryType, Vec<u32>)> {
    let mut counts = OrientationCounts::default();
    let mut sampled = 0usize;
    let requested_samples = param_est_frags;
    let mut sampled_for_lib_type = false;
    let mut temp_frag_lengths_by_ori = vec![vec![0u32; 65_536]; 7];

    'sample: for _ in 0..nchunk {
        let c = chunk::Chunk::<PiscemBulkReadRecord>::from_bytes(br, record_context);
        for mappings in &c.reads {
            if sampled < max_samples {
                for o in &mappings.dirs {
                    counts.add(*o);
                }
                sampled += 1;
                if sampled >= max_samples {
                    sampled_for_lib_type = true;
                }
            }

            let ft = rad_types::MappingType::from_u8(mappings.frag_type);
            let nm = mappings.positions.len();
            if nm == 1
                && !ft.is_orphan()
                && let (Some(o), Some(fl)) = (mappings.dirs.first(), mappings.frag_lengths.first())
            {
                temp_frag_lengths_by_ori[u32::from(*o) as usize][*fl as usize] += 1;
                param_est_frags -= 1;
            }

            if sampled_for_lib_type && param_est_frags <= 0 {
                break 'sample;
            }
        }
    }

    info!(
        "Auto-detection sampled {} reads: forward={}, reverse={}, FR={}, RF={}, FF={}, RR={}, unknown={}",
        sampled,
        counts.forward,
        counts.reverse,
        counts.forward_reverse,
        counts.reverse_forward,
        counts.forward_forward,
        counts.reverse_reverse,
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

    let mut temp_frag_lengths = vec![0u32; 65_536];
    for (ori_idx, ori_counts) in temp_frag_lengths_by_ori.into_iter().enumerate() {
        let ori = MappedFragmentOrientation::from(ori_idx as u32);
        if detected.is_compatible_with(ori) {
            for (dst, src) in temp_frag_lengths.iter_mut().zip(ori_counts.into_iter()) {
                *dst += src;
            }
        }
    }

    if param_est_frags <= 0 {
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

    Ok((detected, temp_frag_lengths))
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

/// Build an equivalence class map from a RAD file, including FLD estimation
/// and effective length computation. This is the common "Phase A" work shared
/// between single-sample `quant` and multi-sample `multi-quant`.
pub fn build_eq_map_from_rad<EqLabelT: EqLabel + Send + 'static>(
    opts: &RadProcessingOpts,
    eqc_map: EqMap<EqLabelT>,
) -> anyhow::Result<EqMapBundle<EqLabelT>> {
    let input = &opts.input;
    let fld_mean = opts.fld_mean;
    let fld_sd = opts.fld_sd;

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

    let mut input_rad = input.clone();
    input_rad.set_extension("rad");
    let i_file = File::open(&input_rad).context("could not open input rad file")?;
    let mut br = BufReader::new(i_file);

    let paired_end: bool;
    let mut fl_mean = 0_f64;
    let mut fl_sd = 0_f64;

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
    for ft in &prelude.file_tags.tags {
        info!("\tfile-level tag {}", ft.name);
    }

    // read-level
    info!("read {:?} read-level tags", prelude.read_tags.tags.len());
    const FRAG_TYPE_NAME: &str = "frag_map_type";
    let mut had_frag_map_type = false;
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

    const REF_ORI_NAME: &str = "compressed_ori_ref";
    const POS_NAME: &str = "pos";
    const FRAGLEN_NAME: &str = "frag_len";

    let mut found_ref_ori_t = false;
    let mut found_pos_t = false;
    let mut found_fraglen_t = false;
    for at in &prelude.aln_tags.tags {
        info!("\talignment-level tag {}", at.name);
        match at.name.as_str() {
            REF_ORI_NAME => found_ref_ori_t = true,
            POS_NAME => found_pos_t = true,
            FRAGLEN_NAME => found_fraglen_t = true,
            _ => info!("unknown alignment-level tag {}", at.name),
        }
    }
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

    const REF_LENGTHS_NAME: &str = "ref_lengths";
    let file_tag_map = prelude.file_tags.parse_tags_from_bytes(&mut br)?;
    let ref_lengths = match file_tag_map.get(REF_LENGTHS_NAME) {
        Some(rad_types::TagValue::ArrayU32(v)) => v,
        _ => bail!("was not able to read reference lengths from file!"),
    };
    info!(
        "read {} reference lengths",
        ref_lengths.len().to_formatted_string(&Locale::en)
    );

    let tag_context = prelude.get_record_context::<PiscemBulkRecordContext>()?;

    let mut frag_stats = MappedFragStats::new();
    let (lib_type, est_frag_lengths): (LibraryType, Option<Vec<u32>>) = match &opts.lib_type {
        LibTypeArg::Explicit(lt) => {
            info!("Using user-specified library type: {}", lt);
            let est_frag_lengths = if paired_end {
                let file_offset = br.stream_position()?;
                let temp_frag_lens = compute_fld_from_sample(
                    &mut br,
                    prelude.hdr.num_chunks as usize,
                    &tag_context,
                    *lt,
                    opts.param_est_frags,
                )?;
                br.seek(std::io::SeekFrom::Start(file_offset))?;
                Some(temp_frag_lens)
            } else {
                None
            };
            (*lt, est_frag_lengths)
        }
        LibTypeArg::Auto => {
            let file_offset = br.stream_position()?;
            if paired_end {
                let (detected, temp_frag_lens) = detect_lib_type_and_fld_from_sample(
                    &mut br,
                    prelude.hdr.num_chunks as usize,
                    &tag_context,
                    paired_end,
                    opts.auto_detect_samples,
                    opts.param_est_frags,
                )?;
                br.seek(std::io::SeekFrom::Start(file_offset))?;
                (detected, Some(temp_frag_lens))
            } else {
                let detected = detect_lib_type_from_sample(
                    &mut br,
                    prelude.hdr.num_chunks as usize,
                    &tag_context,
                    paired_end,
                    opts.auto_detect_samples,
                )?;
                br.seek(std::io::SeekFrom::Start(file_offset))?;
                (detected, None)
            }
        }
    };

    let fld: Fld = if let Some(est_frag_lengths) = est_frag_lengths {
        Fld::Empirical(EmpiricalFLD::new(est_frag_lengths, f64::MIN_POSITIVE))
    } else {
        Fld::Parametric(ParametricFLD::new(fl_mean, fl_sd, 65_536_usize))
    };

    // Extract values needed after processing before we potentially move prelude/file_tag_map.
    let ref_names = prelude.hdr.ref_names.clone();
    let num_chunks = prelude.hdr.num_chunks as usize;
    let ref_lengths_owned = ref_lengths.to_vec();

    let n_threads = opts.num_threads.max(1);
    let (packed_eq_map, frag_lengths) = if n_threads > 1 {
        process_parallel(
            br,
            prelude,
            file_tag_map,
            lib_type,
            &mut frag_stats,
            &ref_lengths_owned,
            eqc_map,
            fld,
            n_threads,
        )
    } else {
        process(
            &mut br,
            num_chunks,
            &tag_context,
            lib_type,
            &mut frag_stats,
            &ref_lengths_owned,
            eqc_map,
            fld,
            1,
        )
    };

    let cond_means = if paired_end {
        conditional_means(&frag_lengths)
    } else {
        conditional_means_from_params(fl_mean, fl_sd, 65_536_usize)
    };
    let eff_lengths = adjust_ref_lengths(&ref_lengths_owned, &cond_means);

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
    info!(
        "total equivalence map weight = {}",
        packed_eq_map
            .total_weight()
            .to_formatted_string(&Locale::en)
    );

    Ok(EqMapBundle {
        packed_eq_map,
        eff_lengths,
        ref_lengths: ref_lengths_owned,
        ref_names,
        frag_lengths,
        frag_stats,
        ref_sig_json,
        lib_type,
    })
}

pub fn process_bulk_dispatch<EqLabelT: EqLabel + Send + 'static>(
    quant_opts: QuantOpts,
    eqc_map: EqMap<EqLabelT>,
) -> anyhow::Result<()> {
    let output = quant_opts.output.clone();
    let num_threads = quant_opts.num_threads;
    let em_pool = if num_threads > 1 {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build()?,
        )
    } else {
        None
    };

    // if there is a parent directory
    if let Some(p) = output.parent()
        && p != Path::new("")
    {
        create_dir_all(p)?;
    }

    // Build EQ map from RAD (Phase A work)
    let rad_opts = RadProcessingOpts {
        input: quant_opts.input.clone(),
        lib_type: quant_opts.lib_type.clone(),
        param_est_frags: quant_opts.param_est_frags,
        fld_mean: quant_opts.fld_mean,
        fld_sd: quant_opts.fld_sd,
        auto_detect_samples: quant_opts.auto_detect_samples,
        num_threads: quant_opts.num_threads,
    };
    let bundle = build_eq_map_from_rad(&rad_opts, eqc_map)?;

    let transcript_mask = if let Some(mask_path) = quant_opts.transcript_mask.as_deref() {
        info!("Applying transcript mask from {}", mask_path.display());
        let keep = read_transcript_mask(mask_path)?;
        let mask: Vec<bool> = bundle
            .ref_names
            .iter()
            .map(|txp| keep.contains(txp))
            .collect();
        let n_kept = mask.iter().filter(|&&b| b).count();
        let n_missing = keep
            .iter()
            .filter(|txp| !bundle.ref_names.iter().any(|name| name == *txp))
            .count();
        info!(
            "Transcript mask keeps {} / {} transcripts{}",
            n_kept,
            bundle.ref_names.len(),
            if n_missing > 0 {
                format!(", {} IDs not found in reference", n_missing)
            } else {
                String::new()
            }
        );
        Some(mask)
    } else {
        None
    };

    // Optional transcript variable selection runs on the positional/basic
    // packed map (before any EC collapsing).
    let selection_mask = if quant_opts.txp_selection {
        info!("Running transcript variable selection...");
        let stages = quant_opts.selection_stages.clone().unwrap_or_default();
        let result = crate::utils::txp_selection::run_selection_with_stages(
            &bundle.packed_eq_map,
            bundle.ref_names.len(),
            &stages,
        );
        Some(result.keep_mask)
    } else {
        None
    };

    // Decide whether to route EM/SQUAREM/Gibbs/bootstrap through a
    // collapsed (pos-bin-stripped) view of the equivalence class map.
    // The collapse is semantics-preserving for M-steps that read only
    // target_labels() + target_probs() (true for standard EM, bootstrap,
    // and the Gibbs multinomial reassignment). Coverage-smoothing EM
    // must see the positional map because `target_pos_bins()` is
    // required by its M-step — so we fall back to the positional map
    // whenever coverage smoothing is active.
    let cov_smoothing_active = quant_opts.coverage_smooth_rounds > 0 && quant_opts.pos_bins > 1;
    let is_range_factorized =
        std::any::TypeId::of::<EqLabelT>() == std::any::TypeId::of::<RangeFactorizedEqLabel>();
    let use_collapsed = is_range_factorized
        && !quant_opts.no_collapsed_ec_em
        && quant_opts.pos_bins > 1
        && !cov_smoothing_active;

    if use_collapsed {
        // SAFETY: is_range_factorized was verified via TypeId above.
        let pos_map: &PackedEqMap<RangeFactorizedEqLabel> = unsafe {
            &*(&bundle.packed_eq_map as *const PackedEqMap<EqLabelT>
                as *const PackedEqMap<RangeFactorizedEqLabel>)
        };
        info!(
            "building collapsed EC view over {} positional ECs",
            pos_map.len().to_formatted_string(&Locale::en)
        );
        let collapsed = build_collapsed(pos_map);
        info!(
            "collapsed view has {} ECs ({} → {})",
            collapsed.len().to_formatted_string(&Locale::en),
            pos_map.len().to_formatted_string(&Locale::en),
            collapsed.len().to_formatted_string(&Locale::en)
        );
        run_inference_and_output(
            &collapsed.packed,
            &bundle,
            transcript_mask.as_deref(),
            selection_mask.as_deref(),
            em_pool.as_ref(),
            &quant_opts,
            &output,
            /* allow_coverage_smoothing = */ false,
        )
    } else {
        run_inference_and_output(
            &bundle.packed_eq_map,
            &bundle,
            transcript_mask.as_deref(),
            selection_mask.as_deref(),
            em_pool.as_ref(),
            &quant_opts,
            &output,
            /* allow_coverage_smoothing = */ true,
        )
    }
}

/// Run the post-bundle inference (EM / bootstrap / Gibbs) and write all
/// outputs. `packed` is the EC map actually used for inference; it can
/// be either the positional/basic map from `bundle` or a collapsed view
/// of it (with a different `EqLabelT`). `bundle` is only read for its
/// per-transcript data (effective lengths, ref names, fragment stats).
///
/// When `allow_coverage_smoothing` is false, the caller has passed a
/// map without positional bins, so `em_with_coverage` is skipped even if
/// the flags would normally trigger it.
#[allow(clippy::too_many_arguments)]
fn run_inference_and_output<EqLabelT: EqLabel, BundleEqLabelT: EqLabel>(
    packed: &PackedEqMap<EqLabelT>,
    bundle: &EqMapBundle<BundleEqLabelT>,
    transcript_mask: Option<&[bool]>,
    selection_mask: Option<&[bool]>,
    em_pool: Option<&rayon::ThreadPool>,
    quant_opts: &QuantOpts,
    output: &Path,
    allow_coverage_smoothing: bool,
) -> anyhow::Result<()> {
    let max_iter = quant_opts.max_iter;
    let convergence_thresh = quant_opts.convergence_thresh;
    let presence_thresh = quant_opts.presence_thresh;
    let num_bootstraps = quant_opts.num_bootstraps;
    let num_gibbs_samples = quant_opts.num_gibbs_samples;
    let gibbs_thinning_factor = quant_opts.gibbs_thinning_factor;

    let mut eminfo = EMInfo::new(
        packed,
        bundle.eff_lengths.clone(),
        max_iter,
        convergence_thresh,
        presence_thresh,
    );
    if let Some(mask) = transcript_mask {
        eminfo.apply_mask(mask);
    }

    let em_res = if allow_coverage_smoothing
        && quant_opts.coverage_smooth_rounds > 0
        && quant_opts.pos_bins > 1
    {
        crate::utils::em::em_with_coverage(
            &eminfo,
            None,
            quant_opts.pos_bins as usize,
            quant_opts.coverage_smooth_rounds,
            quant_opts.coverage_epsilon,
        )
    } else if !quant_opts.no_squarem {
        if let Some(pool) = em_pool {
            squarem_em_par_with_pool(&eminfo, pool)
        } else {
            squarem_em(&eminfo)
        }
    } else if let Some(pool) = em_pool {
        em_par_with_pool(&eminfo, pool)
    } else {
        em(&eminfo)
    };

    // Apply selection mask: zero out structurally redundant transcripts
    let em_res = if let Some(mask) = selection_mask {
        em_res
            .iter()
            .enumerate()
            .map(|(t, &c)| if mask[t] { c } else { 0.0 })
            .collect()
    } else {
        em_res
    };

    let quant_output = output.to_path_buf().with_additional_extension(".quant");
    io::write_results(
        &quant_output,
        &bundle.ref_names,
        &em_res,
        &bundle.ref_lengths,
        &bundle.eff_lengths,
    )
    .context("failed to write quant output")?;

    {
        let fld_array = UInt32Array::from_vec(bundle.frag_lengths.clone());
        let field = Field::new("fragment_length_dist", fld_array.data_type().clone(), false);
        let chunk = Chunk::new(vec![fld_array.boxed()]);
        let fields = vec![field];
        io::write_fld_file(output, fields, chunk)?;
    }

    if num_bootstraps > 0 {
        info!("performing bootstraps");
        let bootstraps = if let Some(pool) = em_pool {
            do_bootstrap_with_pool(&eminfo, num_bootstraps, pool)
        } else {
            do_bootstrap(&eminfo, num_bootstraps)
        };

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
        io::write_infrep_file(output, bs_fields, chunk)?;
    }

    if num_gibbs_samples > 0 {
        info!(
            "performing Gibbs sampling ({num_gibbs_samples} samples, thinning factor {gibbs_thinning_factor})"
        );
        let gibbs_samples = if let Some(pool) = em_pool {
            do_gibbs_with_pool(
                &eminfo,
                &em_res,
                num_gibbs_samples,
                gibbs_thinning_factor,
                pool,
            )
        } else {
            do_gibbs(&eminfo, &em_res, num_gibbs_samples, gibbs_thinning_factor)
        };

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
        io::write_infrep_file(output, gs_fields, chunk)?;
    }

    let infrep_method = if num_gibbs_samples > 0 {
        "gibbs"
    } else if num_bootstraps > 0 {
        "bootstrap"
    } else {
        "none"
    };

    let meta_info_output = output
        .to_path_buf()
        .with_additional_extension(".meta_info.json");
    let ofile = File::create(meta_info_output)?;
    let meta_info = json!({
        "quant_opts": quant_opts,
        "inferred_lib_type": bundle.lib_type.to_string(),
        "mapped_frag_stats": bundle.frag_stats,
        "num_bootstraps": num_bootstraps,
        "num_gibbs_samples": num_gibbs_samples,
        "infrep_method": infrep_method,
        "num_targets": bundle.eff_lengths.len(),
        "signatures": bundle.ref_sig_json
    });
    serde_json::to_writer_pretty(ofile, &meta_info)?;
    Ok(())
}

/// Parallel EQ map building using libradicl's ParallelRadReader.
/// The reader thread fills a lock-free ArrayQueue with MetaChunks;
/// N worker threads pop MetaChunks, iterate their chunks, and build
/// thread-local EqMaps which are merged at the end.
#[allow(clippy::too_many_arguments)]
fn process_parallel<EqLabelT: EqLabel + Send + 'static>(
    reader: BufReader<File>,
    prelude: RadPrelude,
    file_tag_map: TagMap,
    lib_type: LibraryType,
    mapped_stats: &mut MappedFragStats,
    ref_lengths: &[u32],
    eqmap: EqMap<EqLabelT>,
    fld_pdf: Fld,
    num_threads: usize,
) -> (PackedEqMap<EqLabelT>, Vec<u32>) {
    let n_workers = num_threads;
    let contains_ori = eqmap.contains_ori;

    let mut rad_reader =
        ParallelRadReader::<PiscemBulkReadRecord, BufReader<File>>::from_prelude_and_file_tag_map(
            reader,
            prelude,
            file_tag_map,
            std::num::NonZeroUsize::new(n_workers).unwrap(),
        );

    let pb = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr_with_hz(1));
    pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} Processed {human_pos} reads [{elapsed_precise}]",
        )
        .unwrap()
        .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
    );
    pb.enable_steady_tick(Duration::from_secs(1));

    // Clone shared state for workers.
    let queue = rad_reader.get_queue();
    let done = rad_reader.is_done();

    // Spawn worker threads before starting chunk parsing (as required by ParallelRadReader).
    let ref_lengths_arc = std::sync::Arc::new(ref_lengths.to_vec());
    let fld_arc = std::sync::Arc::new(fld_pdf);
    let pb_clone = pb.clone();

    let handles: Vec<_> = (0..n_workers)
        .map(|worker_id| {
            let q = queue.clone();
            let rd = done.clone();
            let rl = ref_lengths_arc.clone();
            let fld = fld_arc.clone();
            let pbc = pb_clone.clone();
            std::thread::spawn(move || {
                let mut state = WorkerState {
                    eqmap: EqMap::new(if contains_ori {
                        OrientationProperty::OrientationAware
                    } else {
                        OrientationProperty::OrientationAgnostic
                    }),
                    frag_lengths: vec![0u32; 65_536],
                    stats: MappedFragStats::default(),
                    unique_frags: 0,
                };

                // Use the FLD enum to get the concrete FldPDF impl.
                loop {
                    while let Some(meta_chunk) = q.pop() {
                        for chunk in meta_chunk.iter() {
                            pbc.inc(chunk.nrec as u64);
                            match fld.as_ref() {
                                Fld::Empirical(f) => {
                                    process_chunk(&chunk, lib_type, &rl, f, &mut state)
                                }
                                Fld::Parametric(f) => {
                                    process_chunk(&chunk, lib_type, &rl, f, &mut state)
                                }
                            }
                        }
                    }
                    if rd.load(std::sync::atomic::Ordering::SeqCst) {
                        // Drain any remaining items.
                        while let Some(meta_chunk) = q.pop() {
                            for chunk in meta_chunk.iter() {
                                pbc.inc(chunk.nrec as u64);
                                match fld.as_ref() {
                                    Fld::Empirical(f) => {
                                        process_chunk(&chunk, lib_type, &rl, f, &mut state)
                                    }
                                    Fld::Parametric(f) => {
                                        process_chunk(&chunk, lib_type, &rl, f, &mut state)
                                    }
                                }
                            }
                        }
                        break;
                    }
                    std::hint::spin_loop();
                }

                tracing::debug!("EQ map worker {} finished", worker_id);
                state
            })
        })
        .collect();

    // Main thread: fill the work queue (blocks until all chunks are enqueued).
    let _ = rad_reader.start_chunk_parsing(libradicl::readers::EMPTY_METACHUNK_CALLBACK);

    // Collect and merge worker results.
    let mut merged = EqMap::new(if contains_ori {
        OrientationProperty::OrientationAware
    } else {
        OrientationProperty::OrientationAgnostic
    });
    let mut frag_lengths = vec![0u32; 65_536];
    let mut unique_frags = 0u32;

    for handle in handles {
        let state = handle.join().unwrap();
        mapped_stats.tot_mappings += state.stats.tot_mappings;
        mapped_stats.num_mapped_reads += state.stats.num_mapped_reads;
        for i in 0..7 {
            mapped_stats.mapped_ori_count[i] += state.stats.mapped_ori_count[i];
            mapped_stats.filtered_ori_count[i] += state.stats.filtered_ori_count[i];
        }
        for (dst, src) in frag_lengths.iter_mut().zip(state.frag_lengths.iter()) {
            *dst += src;
        }
        unique_frags += state.unique_frags;
        merged.merge(state.eqmap);
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

    const TARGET_UNIQUE_FRAGS: u32 = 5_000;
    if unique_frags < TARGET_UNIQUE_FRAGS {
        warn!(
            "Only observed {} uniquely-mapped fragments (< threshold of {}), the fragment length distribution estimate may not be robust",
            unique_frags, TARGET_UNIQUE_FRAGS
        );
    }

    let packed_eq_map = PackedEqMap::from_eq_map(&merged);
    (packed_eq_map, frag_lengths)
}

#[allow(clippy::too_many_arguments)]
fn process<T: Read + Send, EqLabelT: EqLabel + Send>(
    br: &mut BufReader<T>,
    nrec: usize,
    record_context: &PiscemBulkRecordContext,
    lib_type: LibraryType,
    mapped_stats: &mut MappedFragStats,
    ref_lengths: &[u32],
    eq_map: EqMap<EqLabelT>,
    fld_pdf: Fld,
    num_threads: usize,
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
            num_threads,
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
            num_threads,
        ),
    }
}

/// Per-worker accumulator for parallel EQ map building.
struct WorkerState<EqLabelT: EqLabel> {
    eqmap: EqMap<EqLabelT>,
    frag_lengths: Vec<u32>,
    stats: MappedFragStats,
    unique_frags: u32,
}

/// Process a single chunk's reads into a worker's local state.
fn process_chunk<D: FldPDF, EqLabelT: EqLabel>(
    chunk: &chunk::Chunk<PiscemBulkReadRecord>,
    lib_type: LibraryType,
    ref_lengths: &[u32],
    fld_pdf: &D,
    state: &mut WorkerState<EqLabelT>,
) {
    let num_pos_bins = crate::utils::eq_maps::NUM_POS_BINS
        .get()
        .copied()
        .unwrap_or(1.0) as u32;
    let mut mapped_ori_count = [0u32; 7];
    let mut filtered_ori_count = [0u32; 7];
    let mut label_ints = vec![];
    let mut dir_ints = vec![];
    let mut probs = vec![];
    let mut pos_bin_ints = vec![];

    for mappings in &chunk.reads {
        let ft = rad_types::MappingType::from_u8(mappings.frag_type);
        let nm = mappings.positions.len();

        state.stats.tot_mappings += nm;
        state.stats.num_mapped_reads += 1;

        mapped_ori_count.fill(0);
        filtered_ori_count.fill(0);
        label_ints.clear();
        dir_ints.clear();
        probs.clear();
        pos_bin_ints.clear();

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
                        let max_frag_len = *pos as usize + 100;
                        fld_pdf.cdf(max_frag_len)
                    }
                    MappedFragmentOrientation::ForwardReverse
                    | MappedFragmentOrientation::ReverseForward => fld_pdf.pdf(*l as usize),
                    _ => 2.0 * f64::MIN_POSITIVE,
                };
                probs.push(frag_len_prob);
                // Compute positional bin from relative position on transcript.
                if num_pos_bins > 1 {
                    let rel_pos = *pos as f64 / ref_lengths[*r as usize] as f64;
                    let bin = (rel_pos * num_pos_bins as f64) as u32;
                    pos_bin_ints.push(bin.min(num_pos_bins - 1));
                }
            } else {
                filtered_ori_count[y as usize] += 1;
            }
        }

        label_ints.append(&mut dir_ints);
        let pos_bins_arg = if num_pos_bins > 1 {
            Some(pos_bin_ints.as_slice())
        } else {
            None
        };
        let eql = EqLabelT::new(&label_ints, Some(&probs), pos_bins_arg);
        state.eqmap.add(eql);

        if nm == 1 && !ft.is_orphan() {
            if let Some(fl) = mappings.frag_lengths.first() {
                state.frag_lengths[*fl as usize] += 1;
            }
            state.unique_frags += 1;
        }

        for i in 0..mapped_ori_count.len() {
            state.stats.mapped_ori_count[i] += if mapped_ori_count[i] > 0 { 1 } else { 0 };
            state.stats.filtered_ori_count[i] += if filtered_ori_count[i] > 0 { 1 } else { 0 };
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn process_dispatch<T: Read + Send, D: FldPDF + Sync, EqLabelT: EqLabel + Send>(
    br: &mut BufReader<T>,
    nrec: usize,
    record_context: &PiscemBulkRecordContext,
    lib_type: LibraryType,
    mapped_stats: &mut MappedFragStats,
    ref_lengths: &[u32],
    fld_pdf: D,
    eqmap: EqMap<EqLabelT>,
    num_threads: usize,
) -> (PackedEqMap<EqLabelT>, Vec<u32>) {
    let pb = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr_with_hz(1));
    pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} Processed {human_pos} reads [{elapsed_precise}]",
        )
        .unwrap()
        .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
    );
    pb.enable_steady_tick(Duration::from_secs(1));

    let n_workers = num_threads.max(1);
    let contains_ori = eqmap.contains_ori;

    let (merged_eqmap, frag_lengths, unique_frags) = if n_workers <= 1 {
        // Sequential path: single worker, no threading overhead.
        let mut state = WorkerState {
            eqmap,
            frag_lengths: vec![0u32; 65_536],
            stats: MappedFragStats::default(),
            unique_frags: 0,
        };
        for _ in 0..nrec {
            let c = chunk::Chunk::<PiscemBulkReadRecord>::from_bytes(br, record_context);
            pb.inc(c.reads.len() as u64);
            process_chunk(&c, lib_type, ref_lengths, &fld_pdf, &mut state);
        }
        *mapped_stats = state.stats;
        (state.eqmap, state.frag_lengths, state.unique_frags)
    } else {
        // Parallel path: main thread reads chunks, dispatches to N workers
        // via round-robin channels. Each worker builds a local EqMap.
        let mut senders = Vec::with_capacity(n_workers);
        let mut receivers = Vec::with_capacity(n_workers);
        for _ in 0..n_workers {
            let (tx, rx) = std::sync::mpsc::sync_channel::<chunk::Chunk<PiscemBulkReadRecord>>(2);
            senders.push(tx);
            receivers.push(Some(rx));
        }

        std::thread::scope(|scope| {
            // Spawn worker threads — each takes ownership of its receiver.
            let worker_handles: Vec<_> = receivers
                .iter_mut()
                .enumerate()
                .map(|(worker_id, rx_slot)| {
                    let rx = rx_slot.take().unwrap();
                    let fld_ref = &fld_pdf;
                    scope.spawn(move || {
                        let mut state = WorkerState {
                            eqmap: EqMap::new(if contains_ori {
                                OrientationProperty::OrientationAware
                            } else {
                                OrientationProperty::OrientationAgnostic
                            }),
                            frag_lengths: vec![0u32; 65_536],
                            stats: MappedFragStats::default(),
                            unique_frags: 0,
                        };
                        while let Ok(chunk) = rx.recv() {
                            process_chunk(&chunk, lib_type, ref_lengths, fld_ref, &mut state);
                        }
                        tracing::debug!("EQ map worker {} finished", worker_id);
                        state
                    })
                })
                .collect();

            // Main thread: read chunks and dispatch round-robin.
            for chunk_idx in 0..nrec {
                let c = chunk::Chunk::<PiscemBulkReadRecord>::from_bytes(br, record_context);
                pb.inc(c.reads.len() as u64);
                let worker = chunk_idx % n_workers;
                let _ = senders[worker].send(c);
            }
            // Drop senders to signal workers to finish.
            drop(senders);

            // Collect results from worker threads.
            let mut merged = EqMap::new(if contains_ori {
                OrientationProperty::OrientationAware
            } else {
                OrientationProperty::OrientationAgnostic
            });
            let mut frag_lens = vec![0u32; 65_536];
            let mut total_unique = 0u32;

            for handle in worker_handles {
                let state = handle.join().unwrap();
                mapped_stats.tot_mappings += state.stats.tot_mappings;
                mapped_stats.num_mapped_reads += state.stats.num_mapped_reads;
                for i in 0..7 {
                    mapped_stats.mapped_ori_count[i] += state.stats.mapped_ori_count[i];
                    mapped_stats.filtered_ori_count[i] += state.stats.filtered_ori_count[i];
                }
                for (dst, src) in frag_lens.iter_mut().zip(state.frag_lengths.iter()) {
                    *dst += src;
                }
                total_unique += state.unique_frags;
                merged.merge(state.eqmap);
            }
            (merged, frag_lens, total_unique)
        })
    };

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

    const TARGET_UNIQUE_FRAGS: u32 = 5_000;
    if unique_frags < TARGET_UNIQUE_FRAGS {
        warn!(
            "Only observed {} uniquely-mapped fragments (< threshold of {}), the fragment length distribution estimate may not be robust",
            unique_frags, TARGET_UNIQUE_FRAGS
        );
    }

    let packed_eq_map = PackedEqMap::from_eq_map(&merged_eqmap);
    drop(merged_eqmap);

    (packed_eq_map, frag_lengths)
}
