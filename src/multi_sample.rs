use std::fs::{File, create_dir_all};
use std::io::BufReader;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use path_tools::WithAdditionalExtension;
use serde::{Deserialize, Serialize};
use serde_json::json;
use tracing::info;

use crate::process_rad::{RadProcessingOpts, build_eq_map_from_rad};
use crate::prog_opts::MultiQuantOpts;
use crate::utils::em::{self, EMInfo};
use crate::utils::eq_maps::{
    BasicEqLabel, BasicEqMap, EqMapType, OrientationProperty, PackedEqMap,
    RangeFactorizedEqMap,
};
use crate::utils::eq_serialize::{
    EqMapTypeTag, SampleMeta, deserialize_eq_map, serialize_eq_map,
};
use crate::utils::hierarchical;
use crate::utils::io;

/// A single entry from the sample manifest.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SampleEntry {
    pub sample_name: String,
    pub condition: String,
    pub rad_path: PathBuf,
    pub output_dir: PathBuf,
}

/// Wrapper for JSON/YAML manifest format.
#[derive(Debug, Deserialize)]
struct ManifestDoc {
    samples: Vec<SampleEntry>,
}

/// Parse the manifest file, auto-detecting format by extension.
pub fn parse_manifest(path: &Path) -> Result<Vec<SampleEntry>> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_lowercase();

    match ext.as_str() {
        "csv" => parse_csv_manifest(path),
        "json" => parse_json_manifest(path),
        "yaml" | "yml" => parse_yaml_manifest(path),
        _ => bail!(
            "Unrecognized manifest extension '.{}'. Expected .csv, .json, .yaml, or .yml",
            ext
        ),
    }
}

fn parse_csv_manifest(path: &Path) -> Result<Vec<SampleEntry>> {
    let mut reader = csv::Reader::from_path(path)
        .with_context(|| format!("Failed to open CSV manifest: {}", path.display()))?;

    let mut entries = Vec::new();
    for result in reader.deserialize() {
        let entry: SampleEntry = result
            .with_context(|| format!("Failed to parse CSV row in: {}", path.display()))?;
        entries.push(entry);
    }

    if entries.is_empty() {
        bail!("Manifest file is empty: {}", path.display());
    }
    Ok(entries)
}

fn parse_json_manifest(path: &Path) -> Result<Vec<SampleEntry>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open JSON manifest: {}", path.display()))?;
    let doc: ManifestDoc = serde_json::from_reader(BufReader::new(file))
        .with_context(|| format!("Failed to parse JSON manifest: {}", path.display()))?;

    if doc.samples.is_empty() {
        bail!("Manifest file has no samples: {}", path.display());
    }
    Ok(doc.samples)
}

fn parse_yaml_manifest(path: &Path) -> Result<Vec<SampleEntry>> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("Failed to read YAML manifest: {}", path.display()))?;
    use saphyr::LoadableYamlNode;
    let docs = saphyr::Yaml::load_from_str(&content)
        .map_err(|e| anyhow::anyhow!("Failed to parse YAML manifest {}: {}", path.display(), e))?;

    if docs.is_empty() {
        bail!("Empty YAML document: {}", path.display());
    }

    let doc = &docs[0];
    let samples = doc["samples"]
        .as_vec()
        .ok_or_else(|| anyhow::anyhow!("YAML manifest must have a 'samples' list"))?;

    let mut entries = Vec::new();
    for (i, sample) in samples.iter().enumerate() {
        let get_str = |key: &str| -> Result<String> {
            sample[key]
                .as_str()
                .map(|s| s.to_string())
                .ok_or_else(|| anyhow::anyhow!("Sample {} missing required field '{}'", i, key))
        };
        entries.push(SampleEntry {
            sample_name: get_str("sample_name")?,
            condition: get_str("condition")?,
            rad_path: PathBuf::from(get_str("rad_path")?),
            output_dir: PathBuf::from(get_str("output_dir")?),
        });
    }

    if entries.is_empty() {
        bail!("Manifest file has no samples: {}", path.display());
    }
    Ok(entries)
}

/// Build the sorted, deduplicated list of condition names from manifest entries.
fn condition_names(samples: &[SampleEntry]) -> Vec<String> {
    let mut conditions: Vec<String> = samples.iter().map(|s| s.condition.clone()).collect();
    conditions.sort();
    conditions.dedup();
    conditions
}

/// Map a condition string to its index in the sorted condition list.
fn condition_index(condition: &str, condition_names: &[String]) -> usize {
    condition_names
        .iter()
        .position(|c| c == condition)
        .expect("condition not found in condition_names list")
}

/// Main entry point for the multi-sample quantification workflow.
pub fn run(opts: &MultiQuantOpts) -> Result<()> {
    // Initialize rayon thread pool for parallel EM
    if opts.num_threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(opts.num_threads)
            .build_global()
            .unwrap_or_else(|_| {
                info!("Rayon global thread pool already initialized, using existing pool");
            });
    }

    let samples = parse_manifest(&opts.manifest)?;
    let conditions = condition_names(&samples);
    info!(
        "Loaded manifest with {} samples across {} conditions",
        samples.len(),
        conditions.len()
    );

    for s in &samples {
        info!(
            "  sample={}, condition={}, rad_path={}",
            s.sample_name,
            s.condition,
            s.rad_path.display()
        );
    }

    // Force BasicEqMap for multi-quant: the serialization/deserialization pipeline
    // in Phase B operates in BasicEqLabel space. RangeFactorized labels can't be
    // reinterpreted as Basic labels without losing correctness.
    if opts.factorized_eqc_bins > 1 {
        info!("Note: multi-quant uses basic equivalence classes (overriding factorized_eqc_bins={})", opts.factorized_eqc_bins);
    }
    let eq_map_type = EqMapType::BasicEqMap;

    if !opts.phase_b_only {
        info!("=== Phase A: Per-sample EQ class building ===");
        for (i, sample) in samples.iter().enumerate() {
            info!(
                "[{}/{}] Building EQ map for sample '{}'",
                i + 1,
                samples.len(),
                sample.sample_name
            );
            let rad_opts = RadProcessingOpts {
                input: sample.rad_path.clone(),
                lib_type: opts.lib_type.clone(),
                param_est_frags: opts.param_est_frags,
                fld_mean: opts.fld_mean,
                fld_sd: opts.fld_sd,
                auto_detect_samples: opts.auto_detect_samples,
            };

            let eqmap_ori = OrientationProperty::OrientationAware;

            match eq_map_type {
                EqMapType::BasicEqMap => {
                    let bundle =
                        build_eq_map_from_rad(&rad_opts, BasicEqMap::new(eqmap_ori))?;
                    let meta = SampleMeta {
                        sample_name: sample.sample_name.clone(),
                        condition: sample.condition.clone(),
                        num_targets: bundle.ref_names.len(),
                        num_eqcs: bundle.packed_eq_map.len(),
                        total_weight: bundle.packed_eq_map.total_weight(),
                        ref_names: bundle.ref_names,
                        ref_lengths: bundle.ref_lengths,
                        eff_lengths: bundle.eff_lengths,
                        eq_map_type: EqMapTypeTag::Basic,
                        num_bins: 1,
                        contains_ori: bundle.packed_eq_map.contains_ori,
                    };
                    serialize_eq_map(&bundle.packed_eq_map, &meta, &sample.output_dir)?;
                    info!(
                        "  Serialized {} EQCs ({} total weight) for '{}'",
                        meta.num_eqcs, meta.total_weight, sample.sample_name
                    );
                }
                EqMapType::RangeFactorizedEqMap => {
                    let bundle = build_eq_map_from_rad(
                        &rad_opts,
                        RangeFactorizedEqMap::new(eqmap_ori),
                    )?;
                    let meta = SampleMeta {
                        sample_name: sample.sample_name.clone(),
                        condition: sample.condition.clone(),
                        num_targets: bundle.ref_names.len(),
                        num_eqcs: bundle.packed_eq_map.len(),
                        total_weight: bundle.packed_eq_map.total_weight(),
                        ref_names: bundle.ref_names,
                        ref_lengths: bundle.ref_lengths,
                        eff_lengths: bundle.eff_lengths,
                        eq_map_type: EqMapTypeTag::RangeFactorized,
                        num_bins: opts.factorized_eqc_bins,
                        contains_ori: bundle.packed_eq_map.contains_ori,
                    };
                    serialize_eq_map(&bundle.packed_eq_map, &meta, &sample.output_dir)?;
                    info!(
                        "  Serialized {} EQCs ({} total weight) for '{}'",
                        meta.num_eqcs, meta.total_weight, sample.sample_name
                    );
                }
            }
        }
        info!("Phase A complete.");
    }

    if !opts.phase_a_only {
        info!("=== Phase B: Joint hierarchical inference ===");

        // Load all serialized EQ maps and metadata
        let mut all_meta: Vec<SampleMeta> = Vec::with_capacity(samples.len());
        let mut all_deser = Vec::with_capacity(samples.len());

        for sample in &samples {
            let pq_path = sample.output_dir.join(format!("{}.eqc.pq", sample.sample_name));
            let meta_path = sample
                .output_dir
                .join(format!("{}.eqmeta.json", sample.sample_name));
            let (deser, meta) = deserialize_eq_map(&pq_path, &meta_path)
                .with_context(|| format!("Failed to load EQ map for sample '{}'", sample.sample_name))?;
            all_deser.push(deser);
            all_meta.push(meta);
        }

        // Validate reference set compatibility across samples
        let ref_meta = &all_meta[0];
        for (i, meta) in all_meta.iter().enumerate().skip(1) {
            if meta.ref_names != ref_meta.ref_names || meta.ref_lengths != ref_meta.ref_lengths {
                bail!(
                    "Reference mismatch: sample '{}' has different references than sample '{}'. \
                     All samples must map against the same reference.",
                    samples[i].sample_name,
                    samples[0].sample_name
                );
            }
        }

        let num_targets = ref_meta.num_targets;
        info!(
            "All {} samples validated against same reference ({} targets)",
            samples.len(),
            num_targets
        );

        // Reconstruct PackedEqMaps (BasicEqLabel only for now)
        // Phase B operates in BasicEqLabel space regardless of the original EQ map type,
        // since we've already resolved conditional probabilities during Phase A.
        let packed_maps: Vec<PackedEqMap<BasicEqLabel>> = all_deser
            .into_iter()
            .map(|d| PackedEqMap::from_raw(d.eq_labels, d.eq_label_starts, d.counts, d.contains_ori))
            .collect();

        // Optional transcript variable selection on merged EC graph.
        let selection_mask: Option<Vec<bool>> = if opts.txp_selection {
            info!("Running multi-sample transcript variable selection...");
            use crate::utils::txp_selection;
            let stages = opts
                .selection_stages
                .clone()
                .unwrap_or_default();
            let indices: Vec<_> = packed_maps
                .iter()
                .map(|m| txp_selection::TranscriptEqIndex::from_packed_eq_map(m, num_targets))
                .collect();
            let eqc_counts: Vec<usize> = packed_maps.iter().map(|m| m.len()).collect();
            let total_eqcs: usize = eqc_counts.iter().sum();
            let merged_index =
                txp_selection::merge_transcript_indices(&indices, &eqc_counts, num_targets);
            let result = txp_selection::run_selection_from_index_with_stages(
                &merged_index,
                num_targets,
                total_eqcs,
                &stages,
            );
            Some(result.keep_mask)
        } else {
            None
        };

        // Initialize Dirichlet-Multinomial hyperparameters.
        // α₀ = prior_weight × average total reads across samples.
        let avg_total_reads: f64 = packed_maps.iter()
            .map(|m| m.total_weight() as f64)
            .sum::<f64>()
            / packed_maps.len() as f64;
        let alpha_0 = opts.prior_weight * avg_total_reads;
        let mut hyperparams =
            hierarchical::init_hyperparams(num_targets, conditions.clone(), alpha_0);

        let mode_label = if opts.spike_slab { "spike-and-slab" } else { "Dirichlet-Multinomial" };
        info!(
            "Starting hierarchical outer loop ({} iterations, α₀={:.1}, mode={})",
            opts.num_outer_iters, alpha_0, mode_label
        );

        // Track convergence metrics per iteration
        let mut convergence_history: Vec<serde_json::Value> = Vec::new();
        let mut converged_at: Option<u32> = None;

        // Per-sample cached counts from previous iteration
        let mut prev_counts: Vec<Vec<f64>> = Vec::new();

        // Per-sample presence masks (used in Dirichlet mode only).
        let mut presence_masks: Vec<Vec<bool>> = Vec::new();

        // Spike-and-slab inclusion probabilities (used in spike-slab mode only).
        let mut gamma: Vec<f64> = vec![1.0; num_targets];

        // Pre-compute condition indices for all samples
        let sample_cond_indices: Vec<usize> = samples
            .iter()
            .map(|s| condition_index(&s.condition, &conditions))
            .collect();

        // Total iterations = 1 (initialization with standard EM) + num_outer_iters
        let total_iters = 1 + opts.num_outer_iters;

        for outer_iter in 0..total_iters {
            let is_init_iter = outer_iter == 0;
            if is_init_iter {
                info!("--- Initialization iteration (standard EM, no prior) ---");
            } else {
                info!(
                    "--- Outer iteration {}/{} (penalized EM) ---",
                    outer_iter, opts.num_outer_iters
                );
            }

            // Snapshot current hyperparams for convergence check
            let old_pi = hyperparams.pi.clone();

            // Inner loop: per-sample EM (standard or penalized)
            let mut results: Vec<hierarchical::SampleResult> = Vec::with_capacity(samples.len());
            let mut new_counts: Vec<Vec<f64>> = Vec::with_capacity(samples.len());

            for (i, sample) in samples.iter().enumerate() {
                info!("  Sample '{}'", sample.sample_name);

                let em_info = EMInfo {
                    eq_map: &packed_maps[i],
                    eff_lens: &all_meta[i].eff_lengths,
                    max_iter: opts.max_em_iter,
                    convergence_thresh: opts.convergence_thresh,
                    presence_thresh: opts.presence_thresh,
                };

                let counts = if is_init_iter {
                    // Standard EM to establish presence mask and data-driven estimates
                    let c = if opts.num_threads > 1 {
                        em::em_par(&em_info, opts.num_threads)
                    } else {
                        em::em(&em_info)
                    };
                    let mask: Vec<bool> = c.iter().map(|&v| v > opts.presence_thresh).collect();
                    let n_present = mask.iter().filter(|&&b| b).count();
                    info!("    {} / {} transcripts present", n_present, num_targets);
                    presence_masks.push(mask);
                    c
                } else {
                    // Penalized EM with pseudo-counts from hierarchical prior
                    let cond_idx = sample_cond_indices[i];
                    let alpha = if opts.spike_slab {
                        hierarchical::compute_pseudo_counts_spike_slab(
                            &hyperparams,
                            cond_idx,
                            &gamma,
                        )
                    } else {
                        hierarchical::compute_pseudo_counts(
                            &hyperparams,
                            cond_idx,
                            &presence_masks[i],
                        )
                    };

                    if opts.num_threads > 1 {
                        em::em_penalized_par(&em_info, &alpha, opts.num_threads)
                    } else {
                        em::em_penalized(&em_info, &alpha)
                    }
                };

                let present: Vec<bool> = counts
                    .iter()
                    .map(|&c| c > opts.presence_thresh)
                    .collect();

                results.push(hierarchical::SampleResult {
                    counts: counts.clone(),
                    present,
                    sample_name: sample.sample_name.clone(),
                    condition_idx: sample_cond_indices[i],
                });

                new_counts.push(counts);
            }

            // M-step: update condition-level mean proportions
            hierarchical::update_condition_means(&results, &mut hyperparams);

            // Check convergence
            let pi_change = hierarchical::convergence_metric(&old_pi, &hyperparams);
            info!("  Convergence: max_pi_change={:.6e}", pi_change);

            convergence_history.push(json!({
                "iteration": if is_init_iter { "init".to_string() } else { outer_iter.to_string() },
                "max_pi_change": pi_change,
            }));

            prev_counts = new_counts;

            if opts.spike_slab && is_init_iter {
                // Empirical Bayes spike-and-slab: compute γ ONCE from init iteration
                // counts. Using penalized counts in later iterations creates a positive
                // feedback loop where pseudo-counts inflate presence → γ grows → repeat.
                let pi_0 = hierarchical::estimate_inclusion_rate(
                    &prev_counts,
                    num_targets,
                    opts.presence_thresh,
                );
                let n_expressed = (pi_0 * num_targets as f64).round() as usize;
                info!(
                    "  Estimated inclusion rate π₀ = {:.4} ({} / {} expressed in ≥1 sample)",
                    pi_0, n_expressed, num_targets
                );

                gamma = hierarchical::compute_inclusion_probabilities(
                    &prev_counts,
                    &sample_cond_indices,
                    conditions.len(),
                    num_targets,
                    opts.presence_thresh,
                    pi_0,
                );

                // Unique EQ class rule: transcripts with strong unique evidence
                // get γ = 1.0 regardless of replicate consistency.
                const UNIQUE_EQ_MIN_COUNT: usize = 10;
                let mut unique_eq_counts = vec![0usize; num_targets];
                for packed_map in packed_maps.iter() {
                    for eqc_idx in 0..packed_map.len() {
                        if packed_map.counts[eqc_idx] > 0
                            && packed_map.num_targets_in_eqc(eqc_idx) == 1
                        {
                            let s = packed_map.eq_label_starts[eqc_idx] as usize;
                            let target_id = packed_map.eq_labels[s] as usize;
                            unique_eq_counts[target_id] += packed_map.counts[eqc_idx];
                        }
                    }
                }
                let mut n_unique_added = 0usize;
                for (t, &count) in unique_eq_counts.iter().enumerate() {
                    if count >= UNIQUE_EQ_MIN_COUNT && gamma[t] < 0.5 {
                        gamma[t] = 1.0;
                        n_unique_added += 1;
                    }
                }
                if n_unique_added > 0 {
                    info!(
                        "  Added {} transcripts via unique EQ class rule (γ set to 1.0)",
                        n_unique_added
                    );
                }

                let n_included = gamma.iter().filter(|&&g| g >= 0.5).count();
                let n_excluded = gamma.iter().filter(|&&g| g < 0.5).count();
                info!(
                    "  Spike-and-slab support: {} included (γ≥0.5), {} excluded",
                    n_included, n_excluded
                );
            } else if !opts.spike_slab && is_init_iter {
                // Dirichlet mode: hard consensus filtering after init iteration
                info!("--- Consensus support filtering ---");

                let num_conditions = conditions.len();
                let mut cond_counts: Vec<Vec<u32>> =
                    vec![vec![0u32; num_targets]; num_conditions];
                let mut cond_n_reps: Vec<u32> = vec![0u32; num_conditions];

                for (i, _sample) in samples.iter().enumerate() {
                    let cond_idx = sample_cond_indices[i];
                    cond_n_reps[cond_idx] += 1;
                    for t in 0..num_targets {
                        if presence_masks[i][t] {
                            cond_counts[cond_idx][t] += 1;
                        }
                    }
                }

                let mut consensus_support = vec![false; num_targets];
                for c in 0..num_conditions {
                    let threshold = (opts.consensus_thresh * cond_n_reps[c] as f64).ceil() as u32;
                    let threshold = threshold.max(1);
                    for t in 0..num_targets {
                        if cond_counts[c][t] >= threshold {
                            consensus_support[t] = true;
                        }
                    }
                }

                // Unique EQ class rule
                const UNIQUE_EQ_MIN_COUNT: usize = 10;
                let n_before_unique = consensus_support.iter().filter(|&&b| b).count();
                let mut unique_eq_counts = vec![0usize; num_targets];
                for packed_map in packed_maps.iter() {
                    for eqc_idx in 0..packed_map.len() {
                        if packed_map.counts[eqc_idx] > 0
                            && packed_map.num_targets_in_eqc(eqc_idx) == 1
                        {
                            let s = packed_map.eq_label_starts[eqc_idx] as usize;
                            let target_id = packed_map.eq_labels[s] as usize;
                            unique_eq_counts[target_id] += packed_map.counts[eqc_idx];
                        }
                    }
                }
                for (t, &count) in unique_eq_counts.iter().enumerate() {
                    if count >= UNIQUE_EQ_MIN_COUNT {
                        consensus_support[t] = true;
                    }
                }
                let n_unique_added = consensus_support.iter().filter(|&&b| b).count()
                    - n_before_unique;
                if n_unique_added > 0 {
                    info!(
                        "  Added {} transcripts via unique EQ class rule",
                        n_unique_added
                    );
                }

                let n_consensus = consensus_support.iter().filter(|&&b| b).count();
                let n_before: usize = presence_masks
                    .iter()
                    .map(|m| m.iter().filter(|&&b| b).count())
                    .sum::<usize>()
                    / samples.len();
                let n_removed = n_before.saturating_sub(n_consensus);
                info!(
                    "  Consensus support: {} transcripts (was ~{}/sample, removed ~{})",
                    n_consensus, n_before, n_removed
                );

                for mask in presence_masks.iter_mut() {
                    *mask = consensus_support.clone();
                }
            }

            // Don't check convergence on the init iteration
            if !is_init_iter && pi_change < 1e-7 {
                converged_at = Some(outer_iter);
                info!("Hierarchical loop converged at iteration {}", outer_iter);
                break;
            }
        }

        // Write per-sample output using cached final counts.
        info!("Writing per-sample quantification results...");
        for (i, sample) in samples.iter().enumerate() {
            // Apply support filtering: hard mask (Dirichlet) or soft γ (spike-and-slab),
            // composed with structural variable selection if enabled.
            let e_counts: Vec<f64> = if opts.spike_slab {
                // MAP decision: include transcript if γ ≥ 0.5 (posterior
                // probability of expression exceeds 0.5). This is the
                // Bayes-optimal binary decision under 0-1 loss.
                prev_counts[i]
                    .iter()
                    .enumerate()
                    .map(|(t, &c)| {
                        let structural_ok = selection_mask.as_ref().is_none_or(|m| m[t]);
                        if gamma[t] >= 0.5 && structural_ok { c } else { 0.0 }
                    })
                    .collect()
            } else {
                prev_counts[i]
                    .iter()
                    .enumerate()
                    .map(|(t, &c)| {
                        let structural_ok = selection_mask.as_ref().is_none_or(|m| m[t]);
                        if presence_masks[i][t] && structural_ok { c } else { 0.0 }
                    })
                    .collect()
            };

            create_dir_all(&sample.output_dir)
                .with_context(|| format!("Failed to create output dir: {}", sample.output_dir.display()))?;

            let quant_path = sample
                .output_dir
                .join(&sample.sample_name)
                .with_additional_extension(".quant");
            io::write_results(
                &quant_path,
                &all_meta[i].ref_names,
                &e_counts,
                &all_meta[i].ref_lengths,
                &all_meta[i].eff_lengths,
            )
            .with_context(|| format!("Failed to write quant for '{}'", sample.sample_name))?;

            info!("  Wrote {}", quant_path.display());
        }

        // Write joint output files to the global output directory
        create_dir_all(&opts.output)?;

        // Hierarchical parameters
        let hp_path = opts.output.join("hierarchical_params.json");
        let hp_json = json!({
            "condition_names": hyperparams.condition_names,
            "pi": hyperparams.pi,
            "alpha_0": hyperparams.alpha_0,
        });
        let hf = File::create(&hp_path)?;
        serde_json::to_writer_pretty(hf, &hp_json)?;
        info!("Wrote hierarchical parameters to {}", hp_path.display());

        // Convergence history
        let conv_path = opts.output.join("convergence.json");
        let conv_json = json!({
            "converged": converged_at.is_some(),
            "converged_at_iteration": converged_at,
            "total_iterations": convergence_history.len(),
            "history": convergence_history,
        });
        let cf = File::create(&conv_path)?;
        serde_json::to_writer_pretty(cf, &conv_json)?;
        info!("Wrote convergence history to {}", conv_path.display());

        // Meta info
        let meta_path = opts.output.join("meta_info.json");
        let meta_json = json!({
            "mode": "hierarchical_multi_sample",
            "num_samples": samples.len(),
            "num_conditions": conditions.len(),
            "condition_names": conditions,
            "num_targets": num_targets,
            "num_outer_iters": opts.num_outer_iters,
            "max_em_iter": opts.max_em_iter,
            "prior_weight": opts.prior_weight,
            "convergence_thresh": opts.convergence_thresh,
            "converged": converged_at.is_some(),
            "converged_at_iteration": converged_at,
            "samples": samples.iter().map(|s| json!({
                "sample_name": s.sample_name,
                "condition": s.condition,
                "output_dir": s.output_dir.to_string_lossy(),
            })).collect::<Vec<_>>(),
        });
        let mf = File::create(&meta_path)?;
        serde_json::to_writer_pretty(mf, &meta_json)?;
        info!("Wrote meta info to {}", meta_path.display());

        info!("Phase B complete.");
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::eq_maps::EqLabel;
    use std::io::{BufRead, Write};
    use tempfile::TempDir;

    #[test]
    fn test_parse_csv_manifest() {
        let tmp = TempDir::new().unwrap();
        let path = tmp.path().join("manifest.csv");
        let mut f = File::create(&path).unwrap();
        writeln!(f, "sample_name,condition,rad_path,output_dir").unwrap();
        writeln!(f, "s1,control,/data/s1/out,/results/s1").unwrap();
        writeln!(f, "s2,treatment,/data/s2/out,/results/s2").unwrap();

        let entries = parse_manifest(&path).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].sample_name, "s1");
        assert_eq!(entries[0].condition, "control");
        assert_eq!(entries[1].sample_name, "s2");
        assert_eq!(entries[1].condition, "treatment");
    }

    #[test]
    fn test_parse_json_manifest() {
        let tmp = TempDir::new().unwrap();
        let path = tmp.path().join("manifest.json");
        let mut f = File::create(&path).unwrap();
        writeln!(
            f,
            r#"{{
                "samples": [
                    {{"sample_name": "s1", "condition": "ctrl", "rad_path": "/d/s1", "output_dir": "/r/s1"}},
                    {{"sample_name": "s2", "condition": "treat", "rad_path": "/d/s2", "output_dir": "/r/s2"}}
                ]
            }}"#
        )
        .unwrap();

        let entries = parse_manifest(&path).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].condition, "ctrl");
        assert_eq!(entries[1].condition, "treat");
    }

    #[test]
    fn test_parse_yaml_manifest() {
        let tmp = TempDir::new().unwrap();
        let path = tmp.path().join("manifest.yaml");
        let mut f = File::create(&path).unwrap();
        writeln!(
            f,
            "samples:\n  - sample_name: s1\n    condition: ctrl\n    rad_path: /d/s1\n    output_dir: /r/s1\n  - sample_name: s2\n    condition: treat\n    rad_path: /d/s2\n    output_dir: /r/s2"
        )
        .unwrap();

        let entries = parse_manifest(&path).unwrap();
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].sample_name, "s1");
        assert_eq!(entries[1].sample_name, "s2");
    }

    #[test]
    fn test_empty_manifest_errors() {
        let tmp = TempDir::new().unwrap();
        let path = tmp.path().join("empty.csv");
        let mut f = File::create(&path).unwrap();
        writeln!(f, "sample_name,condition,rad_path,output_dir").unwrap();

        let result = parse_manifest(&path);
        assert!(result.is_err());
    }

    #[test]
    fn test_unknown_extension_errors() {
        let tmp = TempDir::new().unwrap();
        let path = tmp.path().join("manifest.txt");
        File::create(&path).unwrap();

        let result = parse_manifest(&path);
        assert!(result.is_err());
    }

    /// Build a synthetic PackedEqMap with known ground truth abundances.
    ///
    /// The EQ class structure has:
    /// - Unique EQCs: {0}, {1}, ..., {num_targets-1}
    /// - Ambiguous EQCs: {0,1}, {2,3}, {4,5}, {6,7}
    ///
    /// Each transcript sends `unique_frac` of its reads to its unique EQC
    /// and `(1 - unique_frac)` to its paired ambiguous EQC.
    fn build_synthetic_eq_map(
        theta: &[f64],
        total_reads: usize,
        unique_frac: f64,
    ) -> PackedEqMap<BasicEqLabel> {
        use crate::utils::eq_maps::{BasicEqMap, OrientationProperty};

        let num_targets = theta.len();
        assert!(num_targets % 2 == 0, "need even number of targets for pairing");

        let mut eqm = BasicEqMap::new(OrientationProperty::OrientationAgnostic);

        for t in 0..num_targets {
            // Unique reads for this transcript
            let n_unique =
                (total_reads as f64 * theta[t] * unique_frac).round() as usize;
            for _ in 0..n_unique {
                eqm.add(BasicEqLabel::new(&[t as u32], None));
            }
        }

        // Ambiguous reads for paired transcripts
        for pair_idx in 0..(num_targets / 2) {
            let t1 = pair_idx * 2;
            let t2 = pair_idx * 2 + 1;
            let n_ambig_t1 =
                (total_reads as f64 * theta[t1] * (1.0 - unique_frac)).round() as usize;
            let n_ambig_t2 =
                (total_reads as f64 * theta[t2] * (1.0 - unique_frac)).round() as usize;
            for _ in 0..(n_ambig_t1 + n_ambig_t2) {
                eqm.add(BasicEqLabel::new(&[t1 as u32, t2 as u32], None));
            }
        }

        PackedEqMap::from_eq_map(&eqm)
    }

    /// Parse a .quant TSV file and return (target_name, ecount) pairs.
    fn read_quant_file(path: &Path) -> Vec<(String, f64)> {
        let file = File::open(path).unwrap();
        let reader = std::io::BufReader::new(file);
        let mut results = Vec::new();
        for (i, line) in reader.lines().enumerate() {
            let line = line.unwrap();
            if i == 0 {
                continue; // skip header
            }
            let fields: Vec<&str> = line.split('\t').collect();
            let name = fields[0].to_string();
            let ecount: f64 = fields[4].parse().unwrap();
            results.push((name, ecount));
        }
        results
    }

    /// Full integration test: build synthetic multi-sample data, serialize,
    /// run Phase B hierarchical inference, and verify output.
    ///
    /// Setup:
    /// - 8 transcripts, equal effective lengths
    /// - 2 conditions (control, treatment), 3 replicates each
    /// - Transcripts 0,1 are highly expressed in control (DE down in treatment)
    /// - Transcripts 2,3 are highly expressed in treatment (DE up in treatment)
    /// - Transcripts 4-7 are not differentially expressed
    /// - 30% of reads are ambiguous (shared between transcript pairs)
    #[test]
    fn test_hierarchical_phase_b_synthetic() {
        use crate::prog_opts::{LibTypeArg, MultiQuantOpts};
        use crate::utils::eq_serialize::{EqMapTypeTag, SampleMeta};
        use crate::utils::map_record_types::LibraryType;

        let tmp = TempDir::new().unwrap();
        let num_targets = 8;
        let eff_lengths = vec![200.0; num_targets];
        let ref_lengths = vec![250u32; num_targets];
        let ref_names: Vec<String> =
            (0..num_targets).map(|i| format!("tx{}", i)).collect();

        // Ground truth theta per condition (sums to 1.0)
        let control_theta = [0.30, 0.20, 0.10, 0.08, 0.10, 0.08, 0.07, 0.07];
        let treatment_theta = [0.10, 0.08, 0.30, 0.20, 0.08, 0.10, 0.07, 0.07];

        // Per-replicate log-perturbations to create biological variability
        let perturbations: [[f64; 8]; 3] = [
            [0.10, -0.05, 0.02, -0.03, 0.01, -0.02, 0.05, -0.01],
            [-0.05, 0.10, -0.03, 0.02, -0.01, 0.03, -0.02, 0.01],
            [0.02, -0.03, 0.05, 0.01, -0.02, 0.01, -0.03, 0.02],
        ];

        let total_reads = 500;
        let unique_frac = 0.7;

        let conditions = ["control", "treatment"];
        let base_thetas: [&[f64]; 2] = [&control_theta, &treatment_theta];

        // Build and serialize 6 samples (3 per condition)
        let mut manifest_lines = vec!["sample_name,condition,rad_path,output_dir".to_string()];

        for (cond_idx, &cond_name) in conditions.iter().enumerate() {
            for rep in 0..3 {
                let sample_name = format!("{}_{}", cond_name, rep);
                let sample_dir = tmp.path().join(&sample_name);

                // Compute per-replicate theta with biological variation
                let base = base_thetas[cond_idx];
                let mut theta_rep: Vec<f64> = (0..num_targets)
                    .map(|t| base[t] * perturbations[rep][t].exp())
                    .collect();
                let sum: f64 = theta_rep.iter().sum();
                for t in theta_rep.iter_mut() {
                    *t /= sum;
                }

                let packed = build_synthetic_eq_map(&theta_rep, total_reads, unique_frac);

                let meta = SampleMeta {
                    sample_name: sample_name.clone(),
                    condition: cond_name.to_string(),
                    num_targets,
                    num_eqcs: packed.len(),
                    total_weight: packed.total_weight(),
                    ref_names: ref_names.clone(),
                    ref_lengths: ref_lengths.clone(),
                    eff_lengths: eff_lengths.clone(),
                    eq_map_type: EqMapTypeTag::Basic,
                    num_bins: 1,
                    contains_ori: false,
                };

                serialize_eq_map(&packed, &meta, &sample_dir).unwrap();

                manifest_lines.push(format!(
                    "{},{},unused,{}",
                    sample_name,
                    cond_name,
                    sample_dir.display()
                ));
            }
        }

        // Write manifest CSV
        let manifest_path = tmp.path().join("manifest.csv");
        let mut f = File::create(&manifest_path).unwrap();
        for line in &manifest_lines {
            writeln!(f, "{}", line).unwrap();
        }

        // Run Phase B only (use small iteration counts for fast debug-mode testing)
        let output_dir = tmp.path().join("joint_output");
        let opts = MultiQuantOpts {
            manifest: manifest_path,
            lib_type: LibTypeArg::Explicit(LibraryType::InwardStrandedForward),
            output: output_dir.clone(),
            max_em_iter: 1500,
            convergence_thresh: 0.001,
            presence_thresh: 1e-8,
            num_outer_iters: 5,
            param_est_frags: 500_000,
            fld_mean: None,
            fld_sd: None,
            factorized_eqc_bins: 1,
            num_threads: 1,
            auto_detect_samples: 10_000,
            phase_a_only: false,
            phase_b_only: true,
            prior_weight: 0.25,
            consensus_thresh: 0.67,
            spike_slab: false,
            txp_selection: false,
            selection_stages: None,
        };

        run(&opts).expect("Phase B should succeed");

        // === Verify joint output files ===
        assert!(
            output_dir.join("hierarchical_params.json").exists(),
            "hierarchical_params.json should exist"
        );
        assert!(
            output_dir.join("convergence.json").exists(),
            "convergence.json should exist"
        );
        assert!(
            output_dir.join("meta_info.json").exists(),
            "meta_info.json should exist"
        );

        // === Verify hierarchical parameters show expected DE pattern ===
        let hp: serde_json::Value = serde_json::from_reader(std::io::BufReader::new(
            File::open(output_dir.join("hierarchical_params.json")).unwrap(),
        ))
        .unwrap();

        let pi = hp["pi"].as_array().unwrap();
        // Conditions are sorted alphabetically: control=0, treatment=1
        let pi_control: Vec<f64> = pi[0]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap())
            .collect();
        let pi_treatment: Vec<f64> = pi[1]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap())
            .collect();

        // pi should sum to 1 per condition
        assert!(
            (pi_control.iter().sum::<f64>() - 1.0).abs() < 1e-6,
            "pi_control should sum to 1, got {}",
            pi_control.iter().sum::<f64>()
        );

        // Transcript 0: higher in control than treatment (DE down)
        assert!(
            pi_control[0] > pi_treatment[0],
            "tx0 should have higher pi in control ({:.4}) than treatment ({:.4})",
            pi_control[0],
            pi_treatment[0]
        );

        // Transcript 2: higher in treatment than control (DE up)
        assert!(
            pi_treatment[2] > pi_control[2],
            "tx2 should have higher pi in treatment ({:.4}) than control ({:.4})",
            pi_treatment[2],
            pi_control[2]
        );

        // === Verify per-sample quant files exist and have plausible estimates ===
        use path_tools::WithAdditionalExtension;
        for (cond_idx, &cond_name) in conditions.iter().enumerate() {
            for rep in 0..3 {
                let sample_name = format!("{}_{}", cond_name, rep);
                let sample_dir = tmp.path().join(&sample_name);
                let quant_path = sample_dir
                    .join(&sample_name)
                    .with_additional_extension(".quant");
                assert!(
                    quant_path.exists(),
                    "Quant file should exist: {}",
                    quant_path.display()
                );

                let quant_data = read_quant_file(&quant_path);
                assert_eq!(quant_data.len(), num_targets);

                // Convert ecounts to proportions
                let total: f64 = quant_data.iter().map(|(_, c)| c).sum();
                assert!(total > 0.0, "Total estimated count should be positive");
                let est_theta: Vec<f64> =
                    quant_data.iter().map(|(_, c)| c / total).collect();

                // The top transcript should match ground truth direction:
                // In control, tx0 should be the most abundant
                // In treatment, tx2 should be the most abundant
                let gt = base_thetas[cond_idx];
                let gt_top = gt
                    .iter()
                    .enumerate()
                    .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                    .unwrap()
                    .0;
                let est_top = est_theta
                    .iter()
                    .enumerate()
                    .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                    .unwrap()
                    .0;

                assert_eq!(
                    est_top, gt_top,
                    "Sample '{}': top transcript should be tx{} (gt), got tx{} \
                     (est_theta={:.4?})",
                    sample_name, gt_top, est_top, est_theta
                );

                // Spearman-like check: rank correlation should be positive
                // (simplified: just check the top-2 match)
                let mut gt_ranked: Vec<(usize, f64)> =
                    gt.iter().copied().enumerate().collect();
                gt_ranked.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
                let mut est_ranked: Vec<(usize, f64)> =
                    est_theta.iter().copied().enumerate().collect();
                est_ranked.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

                // Top-2 ground truth transcripts should be in the top-4 estimated
                let top2_gt: Vec<usize> =
                    gt_ranked.iter().take(2).map(|(i, _)| *i).collect();
                let top4_est: Vec<usize> =
                    est_ranked.iter().take(4).map(|(i, _)| *i).collect();
                for &t in &top2_gt {
                    assert!(
                        top4_est.contains(&t),
                        "Sample '{}': ground truth top-2 tx{} should be in \
                         estimated top-4 {:?}",
                        sample_name,
                        t,
                        top4_est
                    );
                }
            }
        }

        // === Verify convergence.json ===
        let conv: serde_json::Value = serde_json::from_reader(std::io::BufReader::new(
            File::open(output_dir.join("convergence.json")).unwrap(),
        ))
        .unwrap();
        let history = conv["history"].as_array().unwrap();
        assert!(
            !history.is_empty(),
            "Convergence history should have at least one entry"
        );
        // Verify each iteration has valid metric fields
        for entry in history {
            assert!(entry["max_pi_change"].as_f64().unwrap() >= 0.0);
        }
    }
}
