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
use std::num::NonZeroUsize;
use std::path::Path;
use std::time::Duration;
use std::{
    fs::{File, create_dir_all},
    io::Seek,
};
use tabled::{Table, Tabled, settings::Style};
use tracing::{info, warn};

use crate::utils::bulk_record::{BulkAlignment, BulkRecord};
use crate::utils::eq_maps::{
    BasicEqMap, EqLabel, EqMap, EqMapType, OrientationProperty, PackedEqMap, RangeFactorizedEqMap,
};
use crate::utils::io;
use crate::utils::map_record_types::{
    LibraryType, OrientationCounts, check_strand_warnings, detect_library_type,
};
use crate::utils::rad_input::{self, ChunkReader};
use crate::{
    fld::FldPDF,
    prog_opts::{LibTypeArg, QuantOpts},
};
use crate::{
    fld::{EmpiricalFLD, Fld, ParametricFLD},
    utils::em::{
        adjust_ref_lengths, build_packed_eq_classes, conditional_means,
        conditional_means_from_params, em_options, run_em,
    },
};
use salmon_infer::{EffLens, EmResult, GibbsOptions, bootstrap, gibbs_sample};

use libradicl::rad_types::{self, MappedFragmentOrientation};
use libradicl::{
    header::RadPrelude,
    readers::{MetaChunk, ParallelChunkReader},
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

    /// add the counts of `other` to these
    fn merge(&mut self, other: &MappedFragStats) {
        self.tot_mappings += other.tot_mappings;
        self.num_mapped_reads += other.num_mapped_reads;
        for (total, n) in self
            .mapped_ori_count
            .iter_mut()
            .zip(&other.mapped_ori_count)
        {
            *total += n;
        }
        for (total, n) in self
            .filtered_ori_count
            .iter_mut()
            .zip(&other.filtered_ori_count)
        {
            *total += n;
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
    chunks: &mut ChunkReader,
    record_context: &PiscemBulkRecordContext,
    lib_type: LibraryType,
    mut param_est_frags: isize,
) -> anyhow::Result<Vec<u32>> {
    let mut temp_frag_lengths = vec![0u32; 65_536];
    let mut sufficient_samples = false;
    let requested_samples = param_est_frags;

    'estimate_fld: for i in 0..chunks.num_chunks() {
        for rec in chunks.read_records(br, record_context, i)? {
            if rec.num_alignments() == 1 && !rec.mapping_type().is_orphan() {
                let aln = rec.alignment(0);
                if lib_type.is_compatible_with(aln.ori) {
                    temp_frag_lengths[aln.frag_len as usize] += 1;
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
    chunks: &mut ChunkReader,
    record_context: &PiscemBulkRecordContext,
    paired_end: bool,
    max_samples: usize,
) -> anyhow::Result<LibraryType> {
    let mut counts = OrientationCounts::default();
    let mut sampled = 0usize;

    'sample: for i in 0..chunks.num_chunks() {
        for rec in chunks.read_records(br, record_context, i)? {
            for aln in rec.alignments() {
                counts.add(aln.ori);
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

pub fn process_bulk(quant_opts: QuantOpts, eq_map_t: EqMapType) -> anyhow::Result<()> {
    // The packed classes the EM consumes drop orientation (see
    // `extract_key_for_packed_map`), so keying the eq map by orientation as
    // well only splits each class into per-orientation entries that are
    // summed back together when the map is packed.
    let eqmap_orientation_status = OrientationProperty::OrientationAgnostic;
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
    let seed = qo.seed;

    // all parallel work (the EM's M-step, bootstrap replicates, Gibbs chains)
    // runs on the global rayon pool; size it once, up front.
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .context("could not initialize the rayon thread pool")?;

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
    let (input_rad, input_map_info) = rad_input::input_paths(&input);
    let ref_sig_json;
    {
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
                .unwrap_or_else(|_| panic!("Couldn't open {input_map_info:?}."));
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

    let i_file = File::open(&input_rad)
        .with_context(|| format!("could not open input RAD file {}", input_rad.display()))?;
    let mut br = BufReader::new(i_file);

    let paired_end: bool;
    let mut fl_mean = 0_f64;
    let mut fl_sd = 0_f64;

    // read the header and tag sections from the rad file
    let prelude = RadPrelude::from_bytes(&mut br)
        .with_context(|| format!("could not read the RAD header of {}", input_rad.display()))?;
    // before anything else, make sure this is a bulk RAD file we can read
    rad_input::check_bulk_layout(&prelude, &input_rad)?;

    info!("read header!");
    if prelude.hdr.is_paired > 0_u8 {
        info!("fragments paired in sequencing");
        paired_end = true;
        if let (Some(flm), Some(flsd)) = (fld_mean, fld_sd) {
            warn!(
                "provided fragment length distribution mean and sd ({flm}, {flsd}), but \
                    the RAD file contains paired-end fragments, so these will be ignored \
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
                input_rad.display()
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
    // parse actual tags (their layout was checked above)
    for rt in &prelude.read_tags.tags {
        info!("\tread-level tag {}", rt.name);
    }

    // alignment-level
    info!(
        "read {:?} alignment-level tags",
        prelude.aln_tags.tags.len()
    );
    // parse actual tags (their layout was checked above)
    for at in &prelude.aln_tags.tags {
        info!("\talignment-level tag {}", at.name);
    }

    // parse the actual file-level tags, and get the reference lengths from them
    let ref_lengths = rad_input::read_file_tags(&prelude, &mut br, &input_rad)?;
    let ref_lengths = &ref_lengths[..];
    rad_input::check_num_chunks(&prelude.hdr, &mut br, &input_rad)?;
    info!(
        "read {} reference lengths",
        ref_lengths.len().to_formatted_string(&Locale::en)
    );

    // extract whatever context we'll need to read the records
    let tag_context = prelude.get_record_context_prefer_roles::<PiscemBulkRecordContext>()?;
    let mut chunks = ChunkReader::new(&input_rad, prelude.hdr.num_chunks as usize);

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
                &mut chunks,
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
            &mut chunks,
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
        &prelude,
        &input_rad,
        num_threads,
        lib_type,
        &mut frag_stats,
        ref_lengths,
        eqc_map,
        fld,
    )?;

    let cond_means = if paired_end {
        conditional_means(&frag_lengths)
    } else {
        conditional_means_from_params(fl_mean, fl_sd, 65_536_usize)
    };
    let eff_lengths = adjust_ref_lengths(ref_lengths, &cond_means);

    // build the flat CSR equivalence classes consumed by salmon-infer; the raw
    // conditional weights are only needed by the Gibbs sampler.
    let packed = build_packed_eq_classes(&packed_eq_map, &eff_lengths, num_gibbs_samples > 0);
    // only the two summary values logged below are needed from the packed map
    // from here on; free it before the EM and the inferential replicates.
    let num_eq_classes = packed_eq_map.len();
    let total_weight: usize = packed_eq_map.total_weight();
    drop(packed_eq_map);
    let em_opts = em_options(
        max_iter,
        convergence_thresh,
        qo.alpha_check_cutoff,
        presence_thresh,
        qo.em_accel.into(),
    );

    let EmResult {
        alphas: em_res,
        iters: em_iters,
        converged: em_converged,
        ..
    } = run_em(&packed, &em_opts);

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
        num_eq_classes.to_formatted_string(&Locale::en)
    );
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
        info!("performing {num_bootstraps} bootstraps");
        let bootstraps = bootstrap(
            &packed,
            &em_opts,
            EffLens::new(&eff_lengths),
            num_bootstraps as u32,
            seed,
        );

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
        info!(
            "performing Gibbs sampling ({num_gibbs_samples} samples, thinning factor {gibbs_thinning_factor})"
        );
        let gibbs_opts = GibbsOptions {
            num_samples: num_gibbs_samples as u32,
            thinning: gibbs_thinning_factor as u32,
            ..GibbsOptions::default()
        };
        let gibbs_samples = gibbs_sample(&packed, &eff_lengths, &em_res, &gibbs_opts, seed);

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
        "em_iters": em_iters,
        "em_converged": em_converged,
        "num_targets": eff_lengths.len(),
        "signatures": ref_sig_json
    });
    serde_json::to_writer_pretty(ofile, &meta_info)?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn process<T: BufRead, EqLabelT: EqLabel>(
    br: T,
    prelude: &RadPrelude,
    input_rad: &Path,
    num_threads: usize,
    lib_type: LibraryType,
    mapped_stats: &mut MappedFragStats,
    ref_lengths: &[u32],
    eq_map: EqMap<EqLabelT>,
    fld_pdf: Fld,
) -> anyhow::Result<(PackedEqMap<EqLabelT>, Vec<u32>)> {
    match fld_pdf {
        Fld::Empirical(f) => process_dispatch(
            br,
            prelude,
            input_rad,
            num_threads,
            lib_type,
            mapped_stats,
            ref_lengths,
            f,
            eq_map,
        ),
        Fld::Parametric(f) => process_dispatch(
            br,
            prelude,
            input_rad,
            num_threads,
            lib_type,
            mapped_stats,
            ref_lengths,
            f,
            eq_map,
        ),
    }
}

/// What the main pass accumulates over the fragments it parses. Each parsing
/// thread fills its own `ParseState`, and they are merged once every chunk has
/// been parsed. All of it is integer counts (and a class -> count map), so the
/// merged result does not depend on which thread parsed which chunk.
struct ParseState<'a, EqLabelT, D> {
    eqmap: EqMap<EqLabelT>,
    frag_lengths: Vec<u32>,
    unique_frags: u32,
    stats: MappedFragStats,
    lib_type: LibraryType,
    ref_lengths: &'a [u32],
    fld_pdf: &'a D,
    // per-fragment scratch space, reused across fragments
    mapped_ori_count: [u32; 7],
    filtered_ori_count: [u32; 7],
    label_ints: Vec<u32>,
    dir_ints: Vec<u32>,
    probs: Vec<f64>,
    key: Vec<u32>,
    /// the first corrupt chunk this state met (its index and the error); the
    /// state parses nothing after it, but keeps taking meta-chunks so the
    /// reader is never left waiting on a full queue
    error: Option<(usize, anyhow::Error)>,
}

impl<'a, EqLabelT: EqLabel, D: FldPDF> ParseState<'a, EqLabelT, D> {
    fn new(
        eqmap: EqMap<EqLabelT>,
        lib_type: LibraryType,
        ref_lengths: &'a [u32],
        fld_pdf: &'a D,
    ) -> Self {
        Self {
            eqmap,
            frag_lengths: vec![0u32; 65_536],
            unique_frags: 0,
            stats: MappedFragStats::new(),
            lib_type,
            ref_lengths,
            fld_pdf,
            mapped_ori_count: [0u32; 7],
            filtered_ori_count: [0u32; 7],
            label_ints: vec![],
            dir_ints: vec![],
            probs: vec![],
            key: vec![],
            error: None,
        }
    }

    fn add_meta_chunk(
        &mut self,
        meta_chunk: &MetaChunk<PiscemBulkReadRecord>,
        ctx: &PiscemBulkRecordContext,
        path: &Path,
    ) {
        if self.error.is_some() {
            return;
        }
        let first = meta_chunk.first_chunk_index();
        for (k, c) in meta_chunk.raw_chunks().enumerate() {
            match rad_input::chunk_records(c.records, c.num_records, ctx, path, first + k) {
                Ok(records) => {
                    for rec in records {
                        self.add_fragment(rec);
                    }
                }
                Err(e) => {
                    self.error = Some((first + k, e));
                    return;
                }
            }
        }
    }

    fn add_fragment(&mut self, rec: BulkRecord<'_>) {
        self.add_alignments(rec.mapping_type(), rec.num_alignments(), rec.alignments());
    }

    /// Count one fragment, of mapping type `ft`, with its `nm` alignments.
    fn add_alignments(
        &mut self,
        ft: rad_types::MappingType,
        nm: usize,
        alignments: impl Iterator<Item = BulkAlignment>,
    ) {
        self.stats.tot_mappings += nm;
        self.stats.num_mapped_reads += 1;

        // reset the counter
        self.mapped_ori_count.fill(0);
        self.filtered_ori_count.fill(0);

        self.label_ints.clear();
        self.dir_ints.clear();
        self.probs.clear();

        let mut first_frag_len = None;
        for aln in alignments {
            let (r, pos, o, l) = (&aln.tid, &aln.pos, &aln.ori, &aln.frag_len);
            first_frag_len.get_or_insert(*l);
            let y = u32::from(*o);
            if self.lib_type.is_compatible_with(*o) {
                self.mapped_ori_count[y as usize] += 1;
                self.label_ints.push(*r);
                if self.eqmap.contains_ori {
                    self.dir_ints.push(y);
                }
                let frag_len_prob = match o {
                    MappedFragmentOrientation::Forward => {
                        let max_frag_len = (self.ref_lengths[*r as usize] - *pos) as usize;
                        self.fld_pdf.cdf(max_frag_len)
                    }
                    MappedFragmentOrientation::Reverse => {
                        // TODO: modify the bulk RAD creation upstream so that
                        // if a read is mapped as an orphan, then the length
                        // field encodes the read length (instead of the maximum
                        // fragment length).  This will allow us to compute the
                        // actual right-most end of the mapped fragment, rather than
                        // just assuming that the read is around length 100.
                        let max_frag_len = *pos as usize + 100;
                        self.fld_pdf.cdf(max_frag_len)
                    }
                    MappedFragmentOrientation::ForwardReverse
                    | MappedFragmentOrientation::ReverseForward => self.fld_pdf.pdf(*l as usize),
                    _ => 2.0 * f64::MIN_POSITIVE,
                };
                self.probs.push(frag_len_prob)
            } else {
                self.filtered_ori_count[y as usize] += 1;
            }
        }

        // A stranded library can filter out every mapping of a fragment;
        // an empty label would be a class with no targets to assign to.
        if !self.label_ints.is_empty() {
            if self.eqmap.contains_ori {
                self.label_ints.append(&mut self.dir_ints);
            }
            EqLabelT::write_key(&self.label_ints, Some(&self.probs), &mut self.key);
            self.eqmap.add_key(&self.key);
        }

        if nm == 1 && !ft.is_orphan() {
            if let Some(fl) = first_frag_len {
                self.frag_lengths[fl as usize] += 1;
            }
            self.unique_frags += 1;
        }

        // update global orientations
        for i in 0..self.mapped_ori_count.len() {
            self.stats.mapped_ori_count[i] += if self.mapped_ori_count[i] > 0 { 1 } else { 0 };
            self.stats.filtered_ori_count[i] += if self.filtered_ori_count[i] > 0 { 1 } else { 0 };
        }
    }
}

/// Parse the chunks of the RAD file `br` (positioned at its first chunk) on
/// `num_threads` threads, and return the state each of them accumulated.
///
/// The calling thread reads the chunks and queues them for `num_threads - 1`
/// worker threads. libradicl's reader waits for room in a full queue by
/// spinning, so instead of burning a core on that wait, the calling thread
/// parses a queued chunk itself whenever it finds the queue full, and helps
/// parse what is left queued once it has read the last chunk. With no worker
/// threads (`num_threads == 1`), it thus parses every chunk itself.
fn parse_chunks<'a, T: BufRead, D: FldPDF + Sync, EqLabelT: EqLabel>(
    br: T,
    prelude: &RadPrelude,
    path: &Path,
    num_threads: usize,
    mut new_state: impl FnMut() -> ParseState<'a, EqLabelT, D>,
    pb: &ProgressBar,
) -> anyhow::Result<Vec<ParseState<'a, EqLabelT, D>>> {
    let num_threads = NonZeroUsize::new(num_threads).unwrap_or(NonZeroUsize::MIN);
    let ctx = prelude.get_record_context_prefer_roles::<PiscemBulkRecordContext>()?;
    let ctx = &ctx;
    let mut reader = ParallelChunkReader::<PiscemBulkReadRecord>::new(prelude, num_threads);
    let queue = reader.get_queue();
    let mut own_state = new_state();

    std::thread::scope(|s| -> anyhow::Result<_> {
        let workers: Vec<_> = (1..num_threads.get())
            .map(|_| {
                let meta_chunks = reader.chunk_iter();
                let mut state = new_state();
                s.spawn(move || {
                    for meta_chunk in meta_chunks {
                        state.add_meta_chunk(&meta_chunk, ctx, path);
                        pb.inc(meta_chunk.num_records() as u64);
                    }
                    state
                })
            })
            .collect();

        // called by the reader each time it has queued a meta-chunk
        let help_when_full = |_nbytes: u64, _nchunks: u64| {
            while queue.is_full() {
                if let Some(meta_chunk) = queue.pop() {
                    own_state.add_meta_chunk(&meta_chunk, ctx, path);
                    pb.inc(meta_chunk.num_records() as u64);
                }
            }
        };
        let read_res = reader.start(br, Some(help_when_full));
        if read_res.is_ok() {
            // parse what is still queued alongside the workers
            while let Some(meta_chunk) = queue.pop() {
                own_state.add_meta_chunk(&meta_chunk, ctx, path);
                pb.inc(meta_chunk.num_records() as u64);
            }
        }

        let mut states: Vec<_> = workers
            .into_iter()
            .map(|w| w.join().expect("a RAD parsing thread panicked"))
            .collect();
        read_res?;
        states.push(own_state);
        // report the corrupt chunk that comes first in the file, whichever
        // thread met it
        if let Some(i) = (0..states.len())
            .filter(|&i| states[i].error.is_some())
            .min_by_key(|&i| states[i].error.as_ref().map(|e| e.0))
        {
            let (_, e) = states.swap_remove(i).error.expect("filtered on");
            return Err(e);
        }
        Ok(states)
    })
}

/// Merge the states of the parsing threads: sum their fragment length
/// histograms, unique-fragment counts and mapping stats (the latter into
/// `mapped_stats`), and collect their eq maps (which
/// `PackedEqMap::from_eq_maps` packs together).
fn merge_states<EqLabelT, D>(
    states: Vec<ParseState<'_, EqLabelT, D>>,
    mapped_stats: &mut MappedFragStats,
) -> (Vec<EqMap<EqLabelT>>, Vec<u32>, u32) {
    let mut frag_lengths = vec![0u32; 65_536];
    let mut unique_frags = 0u32;
    let mut eqmaps = Vec::with_capacity(states.len());
    for state in states {
        for (total, n) in frag_lengths.iter_mut().zip(&state.frag_lengths) {
            *total += n;
        }
        unique_frags += state.unique_frags;
        mapped_stats.merge(&state.stats);
        eqmaps.push(state.eqmap);
    }
    (eqmaps, frag_lengths, unique_frags)
}

#[allow(clippy::too_many_arguments)]
fn process_dispatch<T: BufRead, D: FldPDF + Sync, EqLabelT: EqLabel>(
    br: T,
    prelude: &RadPrelude,
    input_rad: &Path,
    num_threads: usize,
    lib_type: LibraryType,
    mapped_stats: &mut MappedFragStats,
    ref_lengths: &[u32],
    fld_pdf: D,
    eqmap: EqMap<EqLabelT>,
) -> anyhow::Result<(PackedEqMap<EqLabelT>, Vec<u32>)> {
    const TARGET_UNIQUE_FRAGS: u32 = 5_000;

    let pb = ProgressBar::with_draw_target(None, ProgressDrawTarget::stderr_with_hz(1));
    pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} Processed {human_pos} reads [{elapsed_precise}]",
        )
        .unwrap()
        .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"),
    );
    pb.enable_steady_tick(Duration::from_secs(1));

    // each parsing thread counts into its own eq map: the first one takes
    // over `eqmap`, the others start from empty maps set up like it
    let contains_ori = eqmap.contains_ori;
    let mut eqmap = Some(eqmap);
    let mut new_state = || {
        let eqmap = eqmap.take().unwrap_or_else(|| EqMap {
            count_map: Default::default(),
            contains_ori,
        });
        ParseState::new(eqmap, lib_type, ref_lengths, &fld_pdf)
    };

    // piscem leaves a chunk count of 0 in the header of a file without
    // chunks, and of one whose mapping run did not finish (it fills the count
    // in at the end). The FLD and library-type pre-passes read `num_chunks`
    // chunks, so read none here either, rather than letting the chunk reader
    // fall back to reading until EOF.
    let states = if prelude.hdr.num_chunks > 0 {
        parse_chunks(br, prelude, input_rad, num_threads, new_state, &pb)?
    } else {
        vec![new_state()]
    };

    let (eqmaps, frag_lengths, unique_frags) = merge_states(states, mapped_stats);

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

    let packed_eq_map = PackedEqMap::from_eq_maps(&eqmaps);
    drop(eqmaps);

    Ok((packed_eq_map, frag_lengths))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::eq_maps::{NUM_BINS, OrientationProperty, RangeFactorizedEqLabel};
    use libradicl::ChunkBuf;
    use libradicl::header::RadHeader;
    use libradicl::rad_types::{RadIntId, RadType, TagDesc, TagSection, TagSectionLabel};

    const NUM_REFS: u32 = 50;
    const FRAGS_PER_CHUNK: usize = 1_000;

    /// pseudo-random paired-end fragments (pairs and orphans) mapping to one
    /// or more targets, in several orientations
    fn fragments(n: usize) -> Vec<PiscemBulkReadRecord> {
        let mut state = 0x9E37_79B9u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as u32
        };
        (0..n)
            .map(|_| {
                let frag_type = if next() % 4 == 0 { 2 } else { 4 };
                let nm = 1 + (next() % 4) as usize;
                let mut rec = PiscemBulkReadRecord {
                    frag_type,
                    dirs: vec![],
                    refs: vec![],
                    positions: vec![],
                    frag_lengths: vec![],
                };
                for _ in 0..nm {
                    rec.dirs.push(match (frag_type, next() % 2) {
                        (4, 0) => MappedFragmentOrientation::ForwardReverse,
                        (4, _) => MappedFragmentOrientation::ReverseForward,
                        (_, 0) => MappedFragmentOrientation::Forward,
                        _ => MappedFragmentOrientation::Reverse,
                    });
                    rec.refs.push(next() % NUM_REFS);
                    rec.positions.push(next() % 500);
                    rec.frag_lengths.push((150 + next() % 300) as u16);
                }
                rec
            })
            .collect()
    }

    /// the prelude and chunk bytes of a bulk RAD file holding `frags`
    fn rad_file(frags: &[PiscemBulkReadRecord]) -> (RadPrelude, Vec<u8>) {
        let mut hdr = RadHeader::from_ref_names((0..NUM_REFS).map(|i| format!("t{i}")));
        hdr.is_paired = 1;
        hdr.num_chunks = frags.len().div_ceil(FRAGS_PER_CHUNK) as u64;
        let mut read_tags = TagSection::new_with_label(TagSectionLabel::ReadTags);
        read_tags.add_tag_desc(TagDesc::new("frag_map_type", RadType::Int(RadIntId::U8)));
        let prelude = RadPrelude::from_header_and_tag_sections(
            hdr,
            TagSection::new_with_label(TagSectionLabel::FileTags),
            read_tags,
            TagSection::new_with_label(TagSectionLabel::AlignmentTags),
        );
        let ctx = prelude
            .get_record_context_prefer_roles::<PiscemBulkRecordContext>()
            .unwrap();
        let mut bytes = vec![];
        for chunk_frags in frags.chunks(FRAGS_PER_CHUNK) {
            let mut c = ChunkBuf::new();
            for rec in chunk_frags {
                c.write_record(rec, &ctx).unwrap();
            }
            bytes.extend(c.into_bytes());
        }
        (prelude, bytes)
    }

    struct Parsed {
        eq_labels: Vec<u32>,
        eq_label_starts: Vec<usize>,
        counts: Vec<usize>,
        frag_lengths: Vec<u32>,
        unique_frags: u32,
        stats: Value,
    }

    fn summarize(states: Vec<ParseState<'_, RangeFactorizedEqLabel, ParametricFLD>>) -> Parsed {
        let mut stats = MappedFragStats::new();
        let (eqmaps, frag_lengths, unique_frags) = merge_states(states, &mut stats);
        let packed = PackedEqMap::from_eq_maps(&eqmaps);
        Parsed {
            eq_labels: packed.eq_labels,
            eq_label_starts: packed.eq_label_starts,
            counts: packed.counts,
            frag_lengths,
            unique_frags,
            stats: serde_json::to_value(&stats).unwrap(),
        }
    }

    fn new_eqmap(ori: bool) -> EqMap<RangeFactorizedEqLabel> {
        EqMap::new(if ori {
            OrientationProperty::OrientationAware
        } else {
            OrientationProperty::OrientationAgnostic
        })
    }

    /// The parallel main pass (any thread count, orientation-aware or not)
    /// must give exactly what parsing the chunks one by one into a single
    /// orientation-aware map gives.
    #[test]
    fn parallel_parse_matches_sequential_parse() {
        let _ = NUM_BINS.set(4.0);
        let frags = fragments(300_000);
        let (prelude, bytes) = rad_file(&frags);
        let ref_lengths = vec![1_000u32; NUM_REFS as usize];
        let fld = ParametricFLD::new(250.0, 25.0, 65_536);
        let ctx = prelude
            .get_record_context_prefer_roles::<PiscemBulkRecordContext>()
            .unwrap();

        for lib_type in [
            LibraryType::InwardUnstranded,
            LibraryType::InwardStrandedReverse,
        ] {
            let mut seq = ParseState::new(new_eqmap(true), lib_type, &ref_lengths, &fld);
            let mut reader = std::io::Cursor::new(&bytes[..]);
            for _ in 0..prelude.hdr.num_chunks {
                let c =
                    libradicl::chunk::Chunk::<PiscemBulkReadRecord>::from_bytes(&mut reader, &ctx);
                // the reference: libradicl's own decoding of each record
                for rec in &c.reads {
                    let ft = rad_types::MappingType::from_u8(rec.frag_type);
                    let alns = (0..rec.refs.len()).map(|i| BulkAlignment {
                        tid: rec.refs[i],
                        ori: rec.dirs[i],
                        pos: rec.positions[i],
                        frag_len: rec.frag_lengths[i],
                    });
                    seq.add_alignments(ft, rec.refs.len(), alns);
                }
            }
            let expected = summarize(vec![seq]);
            assert_eq!(expected.stats["num_mapped_reads"], frags.len());
            assert!(expected.counts.len() > 1_000);

            for ori in [true, false] {
                for num_threads in [1, 2, 3, 8] {
                    let states = parse_chunks(
                        &bytes[..],
                        &prelude,
                        Path::new("test.rad"),
                        num_threads,
                        || ParseState::new(new_eqmap(ori), lib_type, &ref_lengths, &fld),
                        &ProgressBar::hidden(),
                    )
                    .unwrap();
                    assert_eq!(states.len(), num_threads);
                    let got = summarize(states);
                    let what = format!("{lib_type}, ori {ori}, {num_threads} threads");
                    assert_eq!(got.eq_labels, expected.eq_labels, "{what}");
                    assert_eq!(got.eq_label_starts, expected.eq_label_starts, "{what}");
                    assert_eq!(got.counts, expected.counts, "{what}");
                    assert_eq!(got.frag_lengths, expected.frag_lengths, "{what}");
                    assert_eq!(got.unique_frags, expected.unique_frags, "{what}");
                    assert_eq!(got.stats, expected.stats, "{what}");
                }
            }
        }
    }

    /// A truncated file is an error (not a panic or a hang), whatever the
    /// number of threads.
    #[test]
    fn parallel_parse_reports_truncated_file() {
        let _ = NUM_BINS.set(4.0);
        let frags = fragments(50_000);
        let (prelude, bytes) = rad_file(&frags);
        let ref_lengths = vec![1_000u32; NUM_REFS as usize];
        let fld = ParametricFLD::new(250.0, 25.0, 65_536);
        for num_threads in [1, 4] {
            let res = parse_chunks(
                &bytes[..bytes.len() - 100],
                &prelude,
                Path::new("test.rad"),
                num_threads,
                || {
                    ParseState::new(
                        new_eqmap(false),
                        LibraryType::InwardUnstranded,
                        &ref_lengths,
                        &fld,
                    )
                },
                &ProgressBar::hidden(),
            );
            let err = res.err().expect("truncated input must be an error");
            assert!(format!("{err:#}").contains("truncated"), "{err:#}");
        }
    }
}
