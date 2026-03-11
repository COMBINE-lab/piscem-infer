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
use crate::utils::em::EMInfo;
use crate::utils::eq_maps::{
    BasicEqLabel, BasicEqMap, EqMapType, OrientationProperty, PackedEqMap,
    RangeFactorizedEqMap,
};
use crate::utils::eq_serialize::{
    EqMapTypeTag, SampleMeta, deserialize_eq_map, serialize_eq_map,
};
use crate::utils::gradient;
use crate::utils::hierarchical::{self, SamplePosterior as HierSamplePosterior};
use crate::utils::io;
use crate::utils::lbfgs::{PenalizedPrior, penalized_em};

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

        // Initialize hierarchical hyperparameters with vague prior
        let mut hyperparams = hierarchical::init_hyperparams(num_targets, conditions.clone());

        info!(
            "Starting hierarchical outer loop ({} iterations)",
            opts.num_outer_iters
        );

        // Track convergence metrics per iteration
        let mut convergence_history: Vec<serde_json::Value> = Vec::new();
        let mut final_posteriors: Vec<HierSamplePosterior> = Vec::new();
        let mut converged_at: Option<u32> = None;

        // Total iterations = 1 (initialization with flat prior) + num_outer_iters
        let total_iters = 1 + opts.num_outer_iters;

        for outer_iter in 0..total_iters {
            // Iteration 0 uses a flat prior to get unbiased initial estimates;
            // subsequent iterations use the learned hyperparameters.
            let is_init_iter = outer_iter == 0;
            if is_init_iter {
                info!("--- Initialization iteration (flat prior) ---");
            } else {
                info!("--- Outer iteration {}/{} ---", outer_iter, opts.num_outer_iters);
            }

            // Snapshot current hyperparams for convergence check
            let old_nu = hyperparams.nu.clone();
            let old_sigma_sq = hyperparams.sigma_sq.clone();

            // Inner loop: per-sample penalized MAP estimation
            let mut posteriors: Vec<HierSamplePosterior> = Vec::with_capacity(samples.len());

            for (i, sample) in samples.iter().enumerate() {
                info!("  Sample '{}'", sample.sample_name);

                let prior = if is_init_iter {
                    // Flat prior: let EM+L-BFGS find the MLE without shrinkage
                    PenalizedPrior {
                        nu: vec![0.0; num_targets],
                        sigma_sq: vec![1e30; num_targets],
                    }
                } else {
                    PenalizedPrior {
                        nu: hyperparams.nu[condition_index(&sample.condition, &conditions)].clone(),
                        sigma_sq: hyperparams.sigma_sq.clone(),
                    }
                };

                let em_info = EMInfo {
                    eq_map: &packed_maps[i],
                    eff_lens: &all_meta[i].eff_lengths,
                    max_iter: opts.em_warmstart_iters,
                    convergence_thresh: opts.convergence_thresh,
                    presence_thresh: opts.presence_thresh,
                };

                let lbfgs_result = penalized_em(
                    &em_info,
                    &prior,
                    opts.em_warmstart_iters,
                    opts.lbfgs_max_iters,
                    opts.lbfgs_history,
                    opts.num_threads,
                )?;

                posteriors.push(HierSamplePosterior {
                    phi_hat: lbfgs_result.phi_hat,
                    sigma_hat_sq: lbfgs_result.sigma_hat_sq,
                    sample_name: sample.sample_name.clone(),
                    condition_idx: condition_index(&sample.condition, &conditions),
                });
            }

            // M-step: update hyperparameters
            hierarchical::update_condition_means(&posteriors, &mut hyperparams);
            hierarchical::update_biological_variance(&posteriors, &mut hyperparams);

            // Check convergence
            let (nu_change, sigma_change) =
                hierarchical::convergence_metrics(&old_nu, &old_sigma_sq, &hyperparams);
            info!(
                "  Convergence: max_nu_change={:.6e}, max_sigma_sq_change={:.6e}",
                nu_change, sigma_change
            );

            convergence_history.push(json!({
                "iteration": if is_init_iter { "init".to_string() } else { outer_iter.to_string() },
                "max_nu_change": nu_change,
                "max_sigma_sq_change": sigma_change,
            }));

            final_posteriors = posteriors;

            // Don't check convergence on the init iteration
            if !is_init_iter && nu_change < 1e-5 && sigma_change < 1e-5 {
                converged_at = Some(outer_iter);
                info!("Hierarchical loop converged at iteration {}", outer_iter);
                break;
            }
        }

        // Write per-sample output using cached final posteriors
        info!("Writing per-sample quantification results...");
        for (i, sample) in samples.iter().enumerate() {
            let posterior = &final_posteriors[i];
            let theta = gradient::softmax(&posterior.phi_hat);
            let total_weight = packed_maps[i].total_weight() as f64;
            let e_counts: Vec<f64> = theta.iter().map(|&t| t * total_weight).collect();

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

            // Write posterior JSON
            let posterior_path = sample
                .output_dir
                .join(format!("{}.posterior.json", sample.sample_name));
            let posterior_json = json!({
                "phi_hat": posterior.phi_hat,
                "sigma_hat_sq": posterior.sigma_hat_sq,
            });
            let pf = File::create(&posterior_path)?;
            serde_json::to_writer_pretty(pf, &posterior_json)?;

            info!("  Wrote {}", quant_path.display());
        }

        // Write joint output files to the global output directory
        create_dir_all(&opts.output)?;

        // Hierarchical parameters
        let hp_path = opts.output.join("hierarchical_params.json");
        let hp_json = json!({
            "condition_names": hyperparams.condition_names,
            "nu": hyperparams.nu,
            "sigma_sq": hyperparams.sigma_sq,
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
            "em_warmstart_iters": opts.em_warmstart_iters,
            "lbfgs_max_iters": opts.lbfgs_max_iters,
            "lbfgs_history": opts.lbfgs_history,
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
            em_warmstart_iters: 10,
            convergence_thresh: 0.001,
            presence_thresh: 1e-8,
            lbfgs_max_iters: 100,
            lbfgs_history: 5,
            num_outer_iters: 5,
            param_est_frags: 500_000,
            fld_mean: None,
            fld_sd: None,
            factorized_eqc_bins: 1,
            num_threads: 1,
            auto_detect_samples: 10_000,
            phase_a_only: false,
            phase_b_only: true,
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

        let nu = hp["nu"].as_array().unwrap();
        // Conditions are sorted alphabetically: control=0, treatment=1
        let nu_control: Vec<f64> = nu[0]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap())
            .collect();
        let nu_treatment: Vec<f64> = nu[1]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_f64().unwrap())
            .collect();

        // Transcript 0: higher in control than treatment (DE down)
        assert!(
            nu_control[0] > nu_treatment[0],
            "tx0 should have higher nu in control ({:.4}) than treatment ({:.4})",
            nu_control[0],
            nu_treatment[0]
        );

        // Transcript 2: higher in treatment than control (DE up)
        assert!(
            nu_treatment[2] > nu_control[2],
            "tx2 should have higher nu in treatment ({:.4}) than control ({:.4})",
            nu_treatment[2],
            nu_control[2]
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

        // === Verify posterior files exist ===
        for cond_name in &conditions {
            for rep in 0..3 {
                let sample_name = format!("{}_{}", cond_name, rep);
                let sample_dir = tmp.path().join(&sample_name);
                let posterior_path =
                    sample_dir.join(format!("{}.posterior.json", sample_name));
                assert!(
                    posterior_path.exists(),
                    "Posterior file should exist: {}",
                    posterior_path.display()
                );

                let post: serde_json::Value = serde_json::from_reader(
                    std::io::BufReader::new(File::open(&posterior_path).unwrap()),
                )
                .unwrap();
                let phi_hat = post["phi_hat"].as_array().unwrap();
                let sigma_hat_sq = post["sigma_hat_sq"].as_array().unwrap();
                assert_eq!(phi_hat.len(), num_targets);
                assert_eq!(sigma_hat_sq.len(), num_targets);

                // All posterior variances should be positive
                for (t, v) in sigma_hat_sq.iter().enumerate() {
                    assert!(
                        v.as_f64().unwrap() > 0.0,
                        "Posterior variance for tx{} should be positive",
                        t
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
            assert!(entry["max_nu_change"].as_f64().unwrap() >= 0.0);
            assert!(entry["max_sigma_sq_change"].as_f64().unwrap() >= 0.0);
        }
        // Check that convergence metrics decrease over iterations
        if history.len() >= 2 {
            let first_nu = history[0]["max_nu_change"].as_f64().unwrap();
            let last_nu = history.last().unwrap()["max_nu_change"]
                .as_f64()
                .unwrap();
            assert!(
                last_nu <= first_nu,
                "Nu change should generally decrease: first={:.6e}, last={:.6e}",
                first_nu,
                last_nu
            );
        }
    }
}
