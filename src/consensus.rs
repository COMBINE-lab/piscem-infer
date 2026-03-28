//! Two-pass consensus-filtered multi-sample quantification.
//!
//! Phase 1: Build EQ maps and run per-sample EM to get initial abundance estimates.
//! Consensus filter: Keep transcripts expressed in ≥K of N samples.
//! Phase 2: Re-run per-sample EM with non-consensus transcripts masked out
//!          (effective length set to 0), redistributing reads to the consensus set.

use anyhow::{Context, Result};
use path_tools::WithAdditionalExtension;
use rayon::prelude::*;
use serde_json::json;
use std::fs::{File, create_dir_all};
use tracing::info;

use crate::multi_sample::{SampleEntry, parse_manifest};
use crate::process_rad::{EqMapBundle, RadProcessingOpts, build_eq_map_from_rad};
use crate::prog_opts::{ConsensusQuantOpts, FilterMode};
use crate::utils::em::{
    EMInfo, em, em_init, em_par, em_par_init, em_par_with_pool, em_par_with_pool_init,
    squarem_em, squarem_em_par, squarem_em_par_with_pool,
};
use crate::utils::eq_maps::{EqLabel, EqMap, OrientationProperty, PackedEqMap, TargetLabelsRef};
use crate::utils::io;

/// Per-transcript EC-based evidence metrics computed after EM convergence.
struct EvidenceMetrics {
    /// Unique Evidence Score: count-weighted average of posterior share per EC.
    /// High UES means the transcript dominates the ECs that contribute to it.
    ues: Vec<f64>,
    /// Effective support count: number of distinct ECs contributing ≥ min_count
    /// assigned fragments to this transcript.
    support: Vec<u32>,
    /// Count-weighted average EC size for this transcript's contributing ECs.
    /// High values indicate the transcript lives in highly ambiguous EC neighborhoods.
    mean_ec_size: Vec<f64>,
}

/// Compute per-transcript UES and effective support count from converged EM counts.
///
/// For each EC `e` containing transcript `t` with count `C_e`:
///   posterior_share[t,e] = (prob[t,e] * mu[t] / eff_len[t]) / denom_e
///   assigned[t,e] = C_e * posterior_share[t,e]
///
/// UES[t] = Σ_e (assigned[t,e] * share[t,e]) / Σ_e assigned[t,e]
/// support[t] = |{e : assigned[t,e] >= min_support_count}|
fn compute_evidence_metrics<EqLabelT: EqLabel>(
    packed_eq_map: &PackedEqMap<EqLabelT>,
    em_counts: &[f64],
    eff_lens: &[f64],
    min_support_count: f64,
) -> EvidenceMetrics {
    let n_targets = em_counts.len();
    let inv_eff_lens: Vec<f64> = eff_lens
        .iter()
        .map(|&l| {
            let inv = 1.0 / l;
            if inv.is_finite() { inv } else { 0.0 }
        })
        .collect();

    // Accumulators: UES numerator (sum of assigned * share) and denominator (sum of assigned)
    let mut ues_num = vec![0.0f64; n_targets];
    let mut ues_den = vec![0.0f64; n_targets];
    let mut support = vec![0u32; n_targets];
    // Accumulators for mean EC size: weighted sum of EC sizes and total weight
    let mut ec_size_sum = vec![0.0f64; n_targets];
    let mut ec_size_weight = vec![0.0f64; n_targets];

    let mut weights: Vec<f64> = Vec::with_capacity(64);

    for (label, &count) in packed_eq_map.iter_labels().zip(packed_eq_map.counts.iter()) {
        let ec_count = count as f64;
        if ec_count == 0.0 {
            continue;
        }
        let ec_size = label.target_labels().len() as f64;

        // Compute weights (same as EM M-step)
        let mut denom = 0.0f64;
        for (tid, cond_prob) in label.target_labels().iter().zip(label.target_probs()) {
            let w = cond_prob * em_counts[*tid as usize] * inv_eff_lens[*tid as usize];
            weights.push(w);
            denom += w;
        }

        if denom > 1e-8 {
            for (tid, &w) in label.target_labels().iter().zip(weights.iter()) {
                let t = *tid as usize;
                let share = w / denom;
                let assigned = ec_count * share;
                ues_num[t] += assigned * share;
                ues_den[t] += assigned;
                if assigned >= min_support_count {
                    support[t] += 1;
                }
                ec_size_sum[t] += assigned * ec_size;
                ec_size_weight[t] += assigned;
            }
        }
        weights.clear();
    }

    // Finalize UES and mean EC size
    let ues: Vec<f64> = ues_num
        .iter()
        .zip(ues_den.iter())
        .map(|(&num, &den)| if den > 0.0 { num / den } else { 0.0 })
        .collect();

    let mean_ec_size: Vec<f64> = ec_size_sum
        .iter()
        .zip(ec_size_weight.iter())
        .map(|(&s, &w)| if w > 0.0 { s / w } else { 0.0 })
        .collect();

    EvidenceMetrics { ues, support, mean_ec_size }
}

/// Compute TPM from estimated counts and effective lengths.
fn compute_tpm(e_counts: &[f64], eff_lengths: &[f64]) -> Vec<f64> {
    const ONE_MILLION: f64 = 1_000_000.0;
    let denom: f64 = e_counts
        .iter()
        .zip(eff_lengths.iter())
        .map(|(c, l)| if *l > 0.0 { c / l } else { 0.0 })
        .sum::<f64>();
    let inv_denom = if denom > 0.0 {
        ONE_MILLION / denom
    } else {
        0.0
    };
    e_counts
        .iter()
        .zip(eff_lengths.iter())
        .map(|(c, l)| if *l > 0.0 { inv_denom * (c / l) } else { 0.0 })
        .collect()
}

fn phase1_max_iter(opts: &ConsensusQuantOpts) -> u32 {
    opts.phase1_max_iter.unwrap_or(opts.max_iter)
}

fn phase1_convergence_thresh(opts: &ConsensusQuantOpts) -> f64 {
    opts.phase1_convergence_thresh
        .unwrap_or(opts.convergence_thresh)
}

fn phase2_max_iter(opts: &ConsensusQuantOpts) -> u32 {
    opts.phase2_max_iter.unwrap_or(opts.max_iter)
}

fn phase2_convergence_thresh(opts: &ConsensusQuantOpts) -> f64 {
    opts.phase2_convergence_thresh
        .unwrap_or(opts.convergence_thresh)
}

fn consensus_thread_split(opts: &ConsensusQuantOpts, n_samples: usize) -> (usize, usize) {
    let num_threads = opts.num_threads.max(1);
    let outer_threads = if opts.sample_parallelism == 0 {
        // Auto: use up to n_samples concurrent jobs, but keep ≥2 inner threads
        // for EM parallelism when possible.
        n_samples.min(num_threads / 2).max(1)
    } else {
        (opts.sample_parallelism as usize).max(1)
    };
    let outer_threads = outer_threads.min(num_threads).min(n_samples.max(1));
    let inner_threads = (num_threads / outer_threads).max(1);
    (outer_threads, inner_threads)
}

fn phase2_init_counts(init_counts: &[f64], active_mask: &[bool]) -> Vec<f64> {
    init_counts
        .iter()
        .zip(active_mask.iter())
        .map(|(&c, &keep)| if keep && c.is_finite() && c > 0.0 { c } else { 0.0 })
        .collect()
}

/// Core implementation generic over EQ label type.
fn run_dispatch<EqLabelT: EqLabel + Send + Sync + 'static>(
    opts: &ConsensusQuantOpts,
    samples: &[SampleEntry],
) -> Result<()> {
    let n_samples = samples.len();
    let (outer_threads, inner_threads) = consensus_thread_split(opts, n_samples);
    let serial_inner_pool = if outer_threads == 1 && inner_threads > 1 {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(inner_threads)
                .build()
                .unwrap(),
        )
    } else {
        None
    };
    // ====== Build EQ maps for all samples ======
    info!(
        "Building EQ maps for {} samples ({} job(s) in parallel, {} thread(s) per sample)",
        n_samples, outer_threads, inner_threads
    );

    let mut bundles: Vec<EqMapBundle<EqLabelT>> = if outer_threads > 1 {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(outer_threads)
            .build()
            .unwrap();
        let results: Vec<Result<(usize, EqMapBundle<EqLabelT>)>> = pool.install(|| {
            samples
                .par_iter()
                .enumerate()
                .map(|(i, sample)| -> Result<(usize, EqMapBundle<EqLabelT>)> {
                    info!("  [{}/{}] Building EQ map for {}", i + 1, n_samples, sample.sample_name);
                    let rad_opts = RadProcessingOpts {
                        input: sample.rad_path.clone(),
                        lib_type: opts.lib_type.clone(),
                        param_est_frags: opts.param_est_frags,
                        fld_mean: opts.fld_mean,
                        fld_sd: opts.fld_sd,
                        auto_detect_samples: opts.auto_detect_samples,
                        num_threads: 1, // EQ map building is I/O-bound; prefer sample parallelism
                    };
                    let bundle = build_eq_map_from_rad(
                        &rad_opts,
                        EqMap::<EqLabelT>::new(OrientationProperty::OrientationAware),
                    )?;
                    Ok((i, bundle))
                })
                .collect()
        });
        let mut results: Vec<(usize, EqMapBundle<EqLabelT>)> =
            results.into_iter().collect::<Result<Vec<_>>>()?;
        results.sort_by_key(|(i, _)| *i);
        results.into_iter().map(|(_, b)| b).collect()
    } else {
        let mut bs = Vec::with_capacity(n_samples);
        for (i, sample) in samples.iter().enumerate() {
            info!("  [{}/{}] Building EQ map for {}", i + 1, n_samples, sample.sample_name);
            let rad_opts = RadProcessingOpts {
                input: sample.rad_path.clone(),
                lib_type: opts.lib_type.clone(),
                param_est_frags: opts.param_est_frags,
                fld_mean: opts.fld_mean,
                fld_sd: opts.fld_sd,
                auto_detect_samples: opts.auto_detect_samples,
                num_threads: inner_threads,
            };
            let bundle = build_eq_map_from_rad(
                &rad_opts,
                EqMap::<EqLabelT>::new(OrientationProperty::OrientationAware),
            )?;
            bs.push(bundle);
        }
        bs
    };

    let n_targets = bundles[0].ref_names.len();

    // ====== Optional: structural transcript variable selection ======
    let selection_stats = if opts.txp_selection {
        use crate::utils::txp_selection;
        info!("Running structural transcript variable selection on merged EC graph...");
        let stages = opts.selection_stages.clone().unwrap_or_default();
        let indices: Vec<_> = bundles
            .iter()
            .map(|b| txp_selection::TranscriptEqIndex::from_packed_eq_map(&b.packed_eq_map, n_targets))
            .collect();
        let eqc_counts: Vec<usize> = bundles.iter().map(|b| b.packed_eq_map.len()).collect();
        let total_eqcs: usize = eqc_counts.iter().sum();
        let merged_index =
            txp_selection::merge_transcript_indices(&indices, &eqc_counts, n_targets);
        let result = txp_selection::run_selection_from_index_with_stages(
            &merged_index,
            n_targets,
            total_eqcs,
            &stages,
        );
        info!(
            "Structural selection: {} kept, {} removed (from {} transcripts)",
            result.num_kept, result.num_removed, n_targets
        );
        // Zero out effective lengths for removed transcripts in all bundles.
        for bundle in &mut bundles {
            for (t, &keep) in result.keep_mask.iter().enumerate() {
                if !keep {
                    bundle.eff_lengths[t] = 0.0;
                }
            }
        }
        Some((result.num_kept, result.num_removed))
    } else {
        None
    };

    // ====== Phase 1: Run initial per-sample EM ======
    info!("Phase 1: running initial EM for {} samples", n_samples);

    let phase1_counts: Vec<Vec<f64>> = if outer_threads > 1 {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(outer_threads)
            .build()
            .unwrap();
        let results: Vec<Result<(usize, Vec<f64>)>> = pool.install(|| {
            bundles
                .par_iter()
                .enumerate()
                .map(|(i, bundle)| -> Result<(usize, Vec<f64>)> {
                    info!("  [{}/{}] Phase 1 EM for {}", i + 1, n_samples, samples[i].sample_name);
                    let eminfo = EMInfo {
                        eq_map: &bundle.packed_eq_map,
                        eff_lens: &bundle.eff_lengths,
                        max_iter: phase1_max_iter(opts),
                        convergence_thresh: phase1_convergence_thresh(opts),
                        presence_thresh: opts.presence_thresh,
                    };
                    let counts = if !opts.no_phase1_squarem {
                        if inner_threads > 1 {
                            squarem_em_par(&eminfo, inner_threads)
                        } else {
                            squarem_em(&eminfo)
                        }
                    } else if inner_threads > 1 {
                        em_par(&eminfo, inner_threads)
                    } else {
                        em(&eminfo)
                    };
                    Ok((i, counts))
                })
                .collect()
        });
        let mut results: Vec<(usize, Vec<f64>)> =
            results.into_iter().collect::<Result<Vec<_>>>()?;
        results.sort_by_key(|(i, _)| *i);
        results.into_iter().map(|(_, c)| c).collect()
    } else {
        let mut counts_vec = Vec::with_capacity(n_samples);
        for (i, bundle) in bundles.iter().enumerate() {
            info!("  [{}/{}] Phase 1 EM for {}", i + 1, n_samples, samples[i].sample_name);
            let eminfo = EMInfo {
                eq_map: &bundle.packed_eq_map,
                eff_lens: &bundle.eff_lengths,
                max_iter: phase1_max_iter(opts),
                convergence_thresh: phase1_convergence_thresh(opts),
                presence_thresh: opts.presence_thresh,
            };
            let counts = if !opts.no_phase1_squarem {
                if let Some(pool) = serial_inner_pool.as_ref() {
                    squarem_em_par_with_pool(&eminfo, pool)
                } else if inner_threads > 1 {
                    squarem_em_par(&eminfo, inner_threads)
                } else {
                    squarem_em(&eminfo)
                }
            } else if let Some(pool) = serial_inner_pool.as_ref() {
                em_par_with_pool(&eminfo, pool)
            } else if inner_threads > 1 {
                em_par(&eminfo, inner_threads)
            } else {
                em(&eminfo)
            };
            counts_vec.push(counts);
        }
        counts_vec
    };

    // ====== Consensus filter ======
    let filter_mode = &opts.filter_mode;
    let base_min_ecs = opts.min_ec_support;

    // Compute per-transcript adaptive EC support threshold if requested.
    let adaptive_thresh: Option<Vec<u32>> = if opts.adaptive_ec_support
        && matches!(filter_mode, FilterMode::Support)
    {
        let mut max_ec_size = vec![0.0f64; n_targets];
        for (i, counts) in phase1_counts.iter().enumerate() {
            let metrics = compute_evidence_metrics(
                &bundles[i].packed_eq_map,
                counts,
                &bundles[i].eff_lengths,
                opts.min_support_count,
            );
            for (t, &sz) in metrics.mean_ec_size.iter().enumerate() {
                if sz > max_ec_size[t] {
                    max_ec_size[t] = sz;
                }
            }
        }
        let thresh: Vec<u32> = max_ec_size
            .iter()
            .map(|&sz| {
                if sz > 1.0 {
                    base_min_ecs.max(sz.log2().ceil() as u32)
                } else {
                    base_min_ecs
                }
            })
            .collect();
        let n_elevated = thresh.iter().filter(|&&t| t > base_min_ecs).count();
        info!(
            "Adaptive EC support: {} transcripts have elevated threshold (base={})",
            n_elevated, base_min_ecs
        );
        Some(thresh)
    } else {
        None
    };

    // Compute per-sample evidence and count how many samples pass for each transcript
    let mut express_count = vec![0u32; n_targets];

    let filter_label = match filter_mode {
        FilterMode::Tpm => {
            let tpm_threshold = opts.tpm_threshold;
            for (i, counts) in phase1_counts.iter().enumerate() {
                let tpms = compute_tpm(counts, &bundles[i].eff_lengths);
                for (t, &tpm) in tpms.iter().enumerate() {
                    if tpm > tpm_threshold {
                        express_count[t] += 1;
                    }
                }
            }
            format!("TPM > {}", tpm_threshold)
        }
        FilterMode::Ues => {
            let ues_threshold = opts.ues_threshold;
            for (i, counts) in phase1_counts.iter().enumerate() {
                let metrics = compute_evidence_metrics(
                    &bundles[i].packed_eq_map,
                    counts,
                    &bundles[i].eff_lengths,
                    opts.min_support_count,
                );
                for (t, &ues) in metrics.ues.iter().enumerate() {
                    if ues > ues_threshold {
                        express_count[t] += 1;
                    }
                }
            }
            format!("UES > {}", ues_threshold)
        }
        FilterMode::Support => {
            for (i, counts) in phase1_counts.iter().enumerate() {
                let metrics = compute_evidence_metrics(
                    &bundles[i].packed_eq_map,
                    counts,
                    &bundles[i].eff_lengths,
                    opts.min_support_count,
                );
                for (t, &sup) in metrics.support.iter().enumerate() {
                    let thresh = adaptive_thresh.as_ref().map_or(base_min_ecs, |v| v[t]);
                    if sup >= thresh {
                        express_count[t] += 1;
                    }
                }
            }
            format!(
                "EC support >= {}{}",
                base_min_ecs,
                if opts.adaptive_ec_support { " (adaptive)" } else { "" }
            )
        }
    };

    // Determine K threshold
    let min_fraction = opts.min_fraction.unwrap_or_else(|| {
        // Default: (N-1)/N
        (n_samples as f64 - 1.0) / n_samples as f64
    });
    let min_k = ((min_fraction * n_samples as f64).ceil() as u32).max(1);

    // Compute per-condition evidence counts (needed for condition_aware and condition_rescue).
    let has_conditions = samples.iter().any(|s| s.condition != samples[0].condition);
    let condition_data = if has_conditions && (opts.condition_aware_consensus || opts.condition_rescue) {
        let mut condition_names: Vec<String> = samples.iter().map(|s| s.condition.clone()).collect();
        condition_names.sort();
        condition_names.dedup();
        let condition_index = |condition: &str| -> usize {
            condition_names
                .iter()
                .position(|c| c == condition)
                .expect("condition should exist")
        };

        let mut cond_counts = vec![vec![0u32; n_targets]; condition_names.len()];
        let mut cond_reps = vec![0u32; condition_names.len()];
        for (sample_idx, sample) in samples.iter().enumerate() {
            let ci = condition_index(&sample.condition);
            cond_reps[ci] += 1;

            let mut sample_pass = vec![false; n_targets];
            match filter_mode {
                FilterMode::Tpm => {
                    let tpms = compute_tpm(&phase1_counts[sample_idx], &bundles[sample_idx].eff_lengths);
                    for (t, &tpm) in tpms.iter().enumerate() {
                        sample_pass[t] = tpm > opts.tpm_threshold;
                    }
                }
                FilterMode::Ues => {
                    let metrics = compute_evidence_metrics(
                        &bundles[sample_idx].packed_eq_map,
                        &phase1_counts[sample_idx],
                        &bundles[sample_idx].eff_lengths,
                        opts.min_support_count,
                    );
                    for (t, &ues) in metrics.ues.iter().enumerate() {
                        sample_pass[t] = ues > opts.ues_threshold;
                    }
                }
                FilterMode::Support => {
                    let metrics = compute_evidence_metrics(
                        &bundles[sample_idx].packed_eq_map,
                        &phase1_counts[sample_idx],
                        &bundles[sample_idx].eff_lengths,
                        opts.min_support_count,
                    );
                    for (t, &sup) in metrics.support.iter().enumerate() {
                        let thresh = adaptive_thresh.as_ref().map_or(base_min_ecs, |v| v[t]);
                        sample_pass[t] = sup >= thresh;
                    }
                }
            }

            for (t, pass) in sample_pass.into_iter().enumerate() {
                if pass {
                    cond_counts[ci][t] += 1;
                }
            }
        }

        let cond_k: Vec<u32> = cond_reps
            .iter()
            .map(|&nrep| ((min_fraction * nrep as f64).ceil() as u32).max(1))
            .collect();

        Some((condition_names, cond_counts, cond_k))
    } else {
        None
    };

    let consensus_mask: Vec<bool> = if opts.condition_aware_consensus {
        let (condition_names, cond_counts, cond_k) =
            condition_data.as_ref().expect("condition data required for condition-aware consensus");
        let mut mask = vec![false; n_targets];
        for t in 0..n_targets {
            mask[t] = (0..condition_names.len()).any(|ci| cond_counts[ci][t] >= cond_k[ci]);
        }
        info!(
            "Condition-aware consensus enabled across {} conditions",
            condition_names.len()
        );
        mask
    } else if opts.condition_rescue {
        // Strict global consensus + condition-specific rescue.
        let (condition_names, cond_counts, cond_k) =
            condition_data.as_ref().expect("condition data required for condition-rescue");

        let mut mask = vec![false; n_targets];
        let mut n_global = 0usize;
        let mut n_rescued = 0usize;
        for t in 0..n_targets {
            if express_count[t] >= min_k {
                // Passes strict global consensus.
                mask[t] = true;
                n_global += 1;
            } else if (0..condition_names.len()).any(|ci| cond_counts[ci][t] >= cond_k[ci]) {
                // Fails globally but passes within at least one condition.
                mask[t] = true;
                n_rescued += 1;
            }
        }
        info!(
            "Condition rescue: {} pass global, {} rescued from within-condition consensus ({} conditions)",
            n_global, n_rescued, condition_names.len()
        );
        mask
    } else {
        express_count
            .iter()
            .map(|&c| c >= min_k)
            .collect()
    };

    // Save the strict (pre-rescue) mask for Phase 2 EM — rescued transcripts
    // will use Phase 1 estimates instead of Phase 2 re-estimation.
    let strict_mask = consensus_mask.clone();
    let n_consensus = consensus_mask.iter().filter(|&&b| b).count();
    let n_filtered = n_targets - n_consensus;

    info!(
        "Consensus filter ({}{}): K={} (min_fraction={:.2}), {} transcripts pass, {} filtered out",
        filter_label,
        if opts.condition_aware_consensus { ", condition-aware" } else if opts.condition_rescue { ", condition-rescue" } else { "" },
        min_k,
        min_fraction,
        n_consensus,
        n_filtered
    );

    // ====== Gene-level rescue ======
    // For transcripts that fail consensus, check if their gene's total TPM
    // passes consensus. If so, rescue all transcripts of that gene.
    // Gene names are parsed from GENCODE-style pipe-delimited transcript names
    // (field 6, 0-indexed field 5).
    let consensus_mask = {
        let mut mask = consensus_mask;

        // Parse gene name from transcript name (field index 5 in pipe-delimited GENCODE IDs).
        let gene_names: Vec<Option<&str>> = bundles[0]
            .ref_names
            .iter()
            .map(|name| name.split('|').nth(5))
            .collect();

        // Build gene -> transcript indices mapping.
        let mut gene_to_txps: std::collections::HashMap<&str, Vec<usize>> =
            std::collections::HashMap::new();
        for (t, gene) in gene_names.iter().enumerate() {
            if let Some(g) = gene {
                gene_to_txps.entry(g).or_default().push(t);
            }
        }

        // Compute per-gene TPM in each sample from Phase 1 estimates.
        // A gene "passes" in a sample if its total TPM >= 10 (strong signal).
        const GENE_TPM_FLOOR: f64 = 3.0;
        let mut gene_express_count: std::collections::HashMap<&str, u32> =
            std::collections::HashMap::new();
        let mut gene_mean_tpm: std::collections::HashMap<&str, f64> =
            std::collections::HashMap::new();
        let n_samp = phase1_counts.len() as f64;
        for counts in &phase1_counts {
            let tpms = compute_tpm(counts, &bundles[0].eff_lengths);
            let mut gene_tpm: std::collections::HashMap<&str, f64> =
                std::collections::HashMap::new();
            for (t, &tpm) in tpms.iter().enumerate() {
                if let Some(g) = gene_names[t] {
                    *gene_tpm.entry(g).or_insert(0.0) += tpm;
                }
            }
            for (gene, &tpm) in &gene_tpm {
                *gene_mean_tpm.entry(gene).or_insert(0.0) += tpm / n_samp;
                if tpm >= GENE_TPM_FLOOR {
                    *gene_express_count.entry(gene).or_insert(0) += 1;
                }
            }
        }

        // Compute mean Phase 1 TPM per transcript (across all samples) for
        // ranking isoforms within a gene.
        let mean_phase1_tpm: Vec<f64> = {
            let n = phase1_counts.len() as f64;
            let mut mean_tpm = vec![0.0f64; n_targets];
            for counts in &phase1_counts {
                let tpms = compute_tpm(counts, &bundles[0].eff_lengths);
                for (t, &tpm) in tpms.iter().enumerate() {
                    mean_tpm[t] += tpm / n;
                }
            }
            mean_tpm
        };

        // Rescue transcripts of genes that pass gene-level consensus but have
        // no individual transcript in the consensus set. Only rescue the top
        // isoforms (by Phase 1 TPM) to avoid flooding the EM with weak isoforms.
        let mut n_gene_rescued = 0usize;
        let mut n_genes_rescued = 0usize;
        for (gene, txps) in &gene_to_txps {
            // Skip if any transcript already passes consensus.
            if txps.iter().any(|&t| mask[t]) {
                continue;
            }
            // Check gene-level consensus.
            let gene_count = gene_express_count.get(gene).copied().unwrap_or(0);
            if gene_count >= min_k {
                // Sort isoforms by Phase 1 TPM (descending).
                let mut ranked: Vec<(usize, f64)> = txps
                    .iter()
                    .filter(|&&t| bundles[0].eff_lengths[t] > 0.0)
                    .map(|&t| (t, mean_phase1_tpm[t]))
                    .collect();
                ranked.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

                // Rescue top isoforms that collectively capture >= 90% of the
                // gene's total Phase 1 TPM, or at least the top 1.
                let gene_total: f64 = ranked.iter().map(|(_, tpm)| tpm).sum();
                let mut cumulative = 0.0f64;
                let mut n_rescued_here = 0usize;
                for &(t, tpm) in &ranked {
                    mask[t] = true;
                    n_gene_rescued += 1;
                    n_rescued_here += 1;
                    cumulative += tpm;
                    if cumulative >= 0.9 * gene_total && n_rescued_here >= 1 {
                        break;
                    }
                }
                n_genes_rescued += 1;
            }
        }

        if n_gene_rescued > 0 {
            info!(
                "Gene-level rescue: {} genes ({} transcripts) rescued by gene-level TPM consensus",
                n_genes_rescued, n_gene_rescued
            );
        }

        let n_consensus_after = mask.iter().filter(|&&b| b).count();
        if n_consensus_after > n_consensus {
            info!(
                "After gene rescue: {} transcripts pass (was {})",
                n_consensus_after, n_consensus
            );
        }

        mask
    };

    // ====== Phase 2: Re-run EM with consensus-masked effective lengths ======
    info!("Phase 2: re-running EM with consensus-filtered transcript set");

    let phase2_results: Vec<(usize, Vec<f64>)> = if outer_threads > 1 {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(outer_threads)
            .build()
            .unwrap();
        let results: Vec<Result<(usize, Vec<f64>)>> = pool.install(|| {
            samples
                .par_iter()
                .enumerate()
                .map(|(i, sample)| -> Result<(usize, Vec<f64>)> {
                    info!(
                        "  [{}/{}] Re-estimating {}",
                        i + 1,
                        n_samples,
                        sample.sample_name
                    );
                    // Phase 2 EM uses strict consensus only (no rescued transcripts).
                    let masked_eff_lens: Vec<f64> = bundles[i]
                        .eff_lengths
                        .iter()
                        .enumerate()
                        .map(|(t, &el)| if strict_mask[t] { el } else { 0.0 })
                        .collect();
                    let init = if opts.no_phase2_warm_start {
                        None
                    } else {
                        Some(phase2_init_counts(&phase1_counts[i], &strict_mask))
                    };
                    let eminfo = EMInfo {
                        eq_map: &bundles[i].packed_eq_map,
                        eff_lens: &masked_eff_lens,
                        max_iter: phase2_max_iter(opts),
                        convergence_thresh: phase2_convergence_thresh(opts),
                        presence_thresh: opts.presence_thresh,
                    };
                    let em_res = if inner_threads > 1 {
                        em_par_init(&eminfo, init.as_deref(), inner_threads)
                    } else {
                        em_init(&eminfo, init.as_deref())
                    };
                    Ok((i, em_res))
                })
                .collect()
        });
        let mut results: Vec<(usize, Vec<f64>)> =
            results.into_iter().collect::<Result<Vec<_>>>()?;
        results.sort_by_key(|(i, _)| *i);
        results
    } else {
        let mut results = Vec::with_capacity(n_samples);
        for (i, sample) in samples.iter().enumerate() {
            info!(
                "  [{}/{}] Re-estimating {}",
                i + 1,
                n_samples,
                sample.sample_name
            );
            let masked_eff_lens: Vec<f64> = bundles[i]
                .eff_lengths
                .iter()
                .enumerate()
                .map(|(t, &el)| if strict_mask[t] { el } else { 0.0 })
                .collect();
            let init = if opts.no_phase2_warm_start {
                None
            } else {
                Some(phase2_init_counts(&phase1_counts[i], &strict_mask))
            };
            let eminfo = EMInfo {
                eq_map: &bundles[i].packed_eq_map,
                eff_lens: &masked_eff_lens,
                max_iter: phase2_max_iter(opts),
                convergence_thresh: phase2_convergence_thresh(opts),
                presence_thresh: opts.presence_thresh,
            };
            let em_res = if let Some(pool) = serial_inner_pool.as_ref() {
                em_par_with_pool_init(&eminfo, init.as_deref(), pool)
            } else if inner_threads > 1 {
                em_par_init(&eminfo, init.as_deref(), inner_threads)
            } else {
                em_init(&eminfo, init.as_deref())
            };
            results.push((i, em_res));
        }
        results
    };

    for (i, mut em_res) in phase2_results {
        // For rescued transcripts (in consensus_mask but not strict_mask),
        // use Phase 1 estimated counts instead of Phase 2. This avoids
        // FC distortion from read redistribution among rescued transcripts.
        for t in 0..n_targets {
            if consensus_mask[t] && !strict_mask[t] {
                em_res[t] = phase1_counts[i][t];
            }
        }

        let sample = &samples[i];
        create_dir_all(&sample.output_dir)?;
        let output_stem = sample.output_dir.join(&sample.sample_name);
        let quant_output = output_stem.with_additional_extension(".quant");

        io::write_results(
            &quant_output,
            &bundles[i].ref_names,
            &em_res,
            &bundles[i].ref_lengths,
            &bundles[i].eff_lengths, // original eff_lens for reporting
        )
        .with_context(|| {
            format!(
                "failed to write quant output for {}",
                sample.sample_name
            )
        })?;

        // Write per-sample meta_info
        let meta_info_output = output_stem.with_additional_extension(".meta_info.json");
        let ofile = File::create(&meta_info_output)?;
        let mut meta_info = json!({
            "inferred_lib_type": format!("{:?}", bundles[i].lib_type),
            "infrep_method": "none",
            "mapped_frag_stats": {
                "num_mapped_reads": bundles[i].frag_stats.num_mapped_reads,
                "tot_mappings": bundles[i].frag_stats.tot_mappings,
            },
            "num_targets": n_targets,
            "num_consensus_targets": n_consensus,
            "consensus_min_k": min_k,
            "consensus_min_fraction": min_fraction,
            "filter_mode": format!("{}", filter_label),
            "piscem_infer_version": env!("CARGO_PKG_VERSION"),
        });
        if let Some((kept, removed)) = selection_stats {
            meta_info["structural_selection"] = json!({
                "enabled": true,
                "num_kept": kept,
                "num_removed": removed,
            });
        }
        serde_json::to_writer_pretty(ofile, &meta_info)?;
    }

    info!("Done. Wrote output for {} samples.", n_samples);
    Ok(())
}

/// Entry point for consensus-quant subcommand.
pub fn run(opts: &ConsensusQuantOpts) -> Result<()> {
    let samples = parse_manifest(&opts.manifest)?;
    info!(
        "Parsed manifest with {} samples",
        samples.len()
    );

    if samples.is_empty() {
        anyhow::bail!("Consensus filtering requires at least 1 sample");
    }

    // Dispatch based on EQ class type
    if opts.factorized_eqc_bins > 1 {
        run_dispatch::<crate::utils::eq_maps::RangeFactorizedEqLabel>(opts, &samples)
    } else {
        run_dispatch::<crate::utils::eq_maps::BasicEqLabel>(opts, &samples)
    }
}
