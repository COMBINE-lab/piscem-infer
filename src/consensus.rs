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
use crate::utils::collapsed_eq::{CollapsedEqMap, build_collapsed};
use crate::utils::em::{
    EMInfo, em, em_init, em_par, em_par_init, em_par_with_pool, em_par_with_pool_init,
    em_penalized_init, em_penalized_par_init, em_penalized_par_with_pool_init, em_with_coverage,
    squarem_em, squarem_em_par, squarem_em_par_with_pool,
};
use crate::utils::eq_maps::{
    EqLabel, EqMap, OrientationProperty, PackedEqMap, RangeFactorizedEqLabel, TargetLabelsRef,
};
use crate::utils::hierarchical;
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

    EvidenceMetrics {
        ues,
        support,
        mean_ec_size,
    }
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

#[inline]
fn support_threshold_for_transcript(
    base_min_ecs: u32,
    adaptive_thresh: Option<&[u32]>,
    t: usize,
) -> u32 {
    adaptive_thresh.map_or(base_min_ecs, |v| v[t])
}

#[inline]
fn ambiguity_adjusted_ues_margin(ues: f64, mean_ec_size: f64) -> f64 {
    let ambiguity_baseline = if mean_ec_size > 1.0 {
        1.0 / mean_ec_size
    } else {
        0.0
    };
    (ues - ambiguity_baseline).max(0.0)
}

fn sample_expression_mask<EqLabelT: EqLabel>(
    packed_eq_map: &PackedEqMap<EqLabelT>,
    counts: &[f64],
    eff_lengths: &[f64],
    filter_mode: &FilterMode,
    opts: &ConsensusQuantOpts,
    base_min_ecs: u32,
    adaptive_thresh: Option<&[u32]>,
) -> Vec<bool> {
    match filter_mode {
        FilterMode::Tpm => compute_tpm(counts, eff_lengths)
            .into_iter()
            .map(|tpm| tpm > opts.tpm_threshold)
            .collect(),
        FilterMode::Ues => {
            let metrics = compute_evidence_metrics(
                packed_eq_map,
                counts,
                eff_lengths,
                opts.min_support_count,
            );
            metrics
                .ues
                .into_iter()
                .map(|ues| ues > opts.ues_threshold)
                .collect()
        }
        FilterMode::Support => {
            let metrics = compute_evidence_metrics(
                packed_eq_map,
                counts,
                eff_lengths,
                opts.min_support_count,
            );
            metrics
                .support
                .into_iter()
                .enumerate()
                .map(|(t, sup)| {
                    sup >= support_threshold_for_transcript(base_min_ecs, adaptive_thresh, t)
                })
                .collect()
        }
        FilterMode::Hybrid => {
            let metrics = compute_evidence_metrics(
                packed_eq_map,
                counts,
                eff_lengths,
                opts.min_support_count,
            );
            (0..counts.len())
                .map(|t| {
                    let support_thresh =
                        support_threshold_for_transcript(base_min_ecs, adaptive_thresh, t);
                    let relaxed_support = support_thresh.saturating_sub(1).max(1);
                    let dominance_margin =
                        ambiguity_adjusted_ues_margin(metrics.ues[t], metrics.mean_ec_size[t]);
                    metrics.support[t] >= support_thresh
                        || (metrics.support[t] >= relaxed_support
                            && dominance_margin > opts.ues_threshold)
                })
                .collect()
        }
    }
}

fn sample_tpm_rescue_mask(tpms: &[f64], tpm_threshold: f64) -> Vec<bool> {
    tpms.iter().map(|&tpm| tpm > tpm_threshold).collect()
}

fn sample_specific_phase2_masks(
    samples: &[SampleEntry],
    strict_global_mask: &[bool],
    pre_gene_consensus_mask: &[bool],
    condition_names: Option<&[String]>,
    condition_support_masks: Option<&[Vec<bool>]>,
    condition_aware_consensus: bool,
    use_condition_rescue: bool,
) -> Vec<Vec<bool>> {
    if !(condition_aware_consensus || use_condition_rescue) {
        return vec![pre_gene_consensus_mask.to_vec(); samples.len()];
    }

    let condition_names =
        condition_names.expect("condition names required for condition-aware phase-2 masks");
    let condition_support_masks = condition_support_masks
        .expect("condition support masks required for condition-aware phase-2 masks");

    samples
        .iter()
        .map(|sample| {
            let ci = condition_names
                .iter()
                .position(|c| c == &sample.condition)
                .expect("sample condition should exist in condition_names");
            let cond_mask = &condition_support_masks[ci];
            (0..pre_gene_consensus_mask.len())
                .map(|t| cond_mask[t] || (!condition_aware_consensus && strict_global_mask[t]))
                .collect()
        })
        .collect()
}

/// Build gene-like groups from EC graph structure (annotation-free).
///
/// Uses Jaccard similarity of EC signatures to identify transcripts that likely
/// belong to the same gene. Two transcripts are grouped if they share >= 10% of
/// their combined EC signatures (Jaccard >= 0.10). Connected components of this
/// graph form the groups.
///
/// Returns a mapping from transcript index to group ID. Transcripts with no
/// EC-sharing neighbors get their own singleton group.
fn build_ec_graph_groups(
    index: &crate::utils::txp_selection::TranscriptEqIndex,
    consensus_mask: &[bool],
    n_targets: usize,
) -> Vec<u32> {
    let jaccard_threshold = 0.10;

    // Build inverted index: EC → list of consensus transcripts
    let max_eqc = index
        .eqc_ids
        .iter()
        .copied()
        .max()
        .map(|m| m as usize + 1)
        .unwrap_or(0);
    let mut eqc_to_txps: Vec<Vec<u32>> = vec![Vec::new(); max_eqc];
    for t in 0..n_targets {
        if !consensus_mask[t] {
            continue;
        }
        for &eqc in index.signature(t) {
            eqc_to_txps[eqc as usize].push(t as u32);
        }
    }

    // For each transcript, find neighbors with Jaccard >= threshold.
    // Use union-find for efficient connected components.
    let mut parent: Vec<u32> = (0..n_targets as u32).collect();
    let mut rank = vec![0u32; n_targets];

    fn find(parent: &mut [u32], x: u32) -> u32 {
        let mut r = x;
        while parent[r as usize] != r {
            r = parent[r as usize];
        }
        // Path compression
        let mut c = x;
        while c != r {
            let next = parent[c as usize];
            parent[c as usize] = r;
            c = next;
        }
        r
    }
    fn union(parent: &mut [u32], rank: &mut [u32], a: u32, b: u32) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra == rb {
            return;
        }
        if rank[ra as usize] < rank[rb as usize] {
            parent[ra as usize] = rb;
        } else if rank[ra as usize] > rank[rb as usize] {
            parent[rb as usize] = ra;
        } else {
            parent[rb as usize] = ra;
            rank[ra as usize] += 1;
        }
    }

    // For each consensus transcript, compute pairwise Jaccard with EC-sharing neighbors
    for t in 0..n_targets {
        if !consensus_mask[t] {
            continue;
        }
        let sig_t = index.signature(t);
        if sig_t.is_empty() {
            continue;
        }
        let deg_t = sig_t.len() as u32;

        // Count shared ECs with each neighbor
        let mut neighbor_shared: ahash::AHashMap<u32, u32> = ahash::AHashMap::new();
        for &eqc in sig_t {
            for &u in &eqc_to_txps[eqc as usize] {
                if u as usize != t {
                    *neighbor_shared.entry(u).or_insert(0) += 1;
                }
            }
        }

        // Check Jaccard threshold and union
        for (&u, &shared) in &neighbor_shared {
            let deg_u = index.degree(u as usize) as u32;
            let union_size = deg_t + deg_u - shared;
            if union_size > 0 {
                let jaccard = shared as f64 / union_size as f64;
                if jaccard >= jaccard_threshold {
                    union(&mut parent, &mut rank, t as u32, u);
                }
            }
        }
    }

    // Flatten to group IDs
    for t in 0..n_targets {
        find(&mut parent, t as u32);
    }
    parent
}

/// Check if transcript t's position profile is near-identical to transcript dom's,
/// suggesting EM leakage. Returns true if correlation > 0.9.
fn profile_corr_is_leakage(
    pos_bin_profiles: &[Vec<u64>],
    t: usize,
    dom: usize,
    n_pos_bins: usize,
) -> bool {
    let prof_t = &pos_bin_profiles[t];
    let prof_dom = &pos_bin_profiles[dom];
    let total_t: f64 = prof_t.iter().map(|&c| c as f64).sum();
    let total_dom: f64 = prof_dom.iter().map(|&c| c as f64).sum();

    if total_t <= 10.0 || total_dom <= 10.0 {
        return false;
    }

    let dom_bin_t = prof_t
        .iter()
        .enumerate()
        .max_by_key(|&(_, &count)| count)
        .map(|(idx, _)| idx);
    let dom_bin_dom = prof_dom
        .iter()
        .enumerate()
        .max_by_key(|&(_, &count)| count)
        .map(|(idx, _)| idx);
    if dom_bin_t != dom_bin_dom {
        return false;
    }

    let mean_t = 1.0 / n_pos_bins as f64; // fracs sum to 1, mean = 1/n
    let mean_dom = mean_t;

    let mut cov = 0.0f64;
    let mut var_t = 0.0f64;
    let mut var_dom = 0.0f64;
    let mut l1 = 0.0f64;
    for b in 0..n_pos_bins {
        let ft = prof_t[b] as f64 / total_t - mean_t;
        let fd = prof_dom[b] as f64 / total_dom - mean_dom;
        cov += ft * fd;
        var_t += ft * ft;
        var_dom += fd * fd;
        l1 += ((prof_t[b] as f64 / total_t) - (prof_dom[b] as f64 / total_dom)).abs();
    }
    let denom = (var_t * var_dom).sqrt();
    let corr = if denom > 1e-12 { cov / denom } else { 0.0 };
    corr > 0.9 && l1 < 0.35
}

/// Compute per-transcript EC-neighborhood leakage scores.
///
/// For each transcript t, computes the total posterior-weighted count attributed
/// to all other transcripts from t's shared ECs. Returns:
/// - `nbr_frac[t]`: t's abundance as a fraction of (t + all competitors)
/// - `dominant[t]`: index of t's single most abundant competitor
///
/// This generalizes the gene-fraction filter to work without gene annotations:
/// the "gene" is replaced by the EC-neighborhood structure.
fn compute_neighborhood_leakage<EqLabelT: EqLabel>(
    packed_map: &PackedEqMap<EqLabelT>,
    em_counts: &[f64],
    eff_lens: &[f64],
) -> (Vec<f64>, Vec<usize>) {
    let n_targets = em_counts.len();
    let mut nbr_frac = vec![1.0f64; n_targets];
    let mut dominant = (0..n_targets).collect::<Vec<usize>>();

    let inv_eff_lens: Vec<f64> = eff_lens
        .iter()
        .map(|&l| {
            let inv = 1.0 / l;
            if inv.is_finite() { inv } else { 0.0 }
        })
        .collect();

    // For each transcript t, accumulate per-competitor posterior counts
    // across all shared ECs, then compute the neighborhood fraction as
    // t's count relative to the top-K competitors that collectively account
    // for most of the competition.
    let mut competitor_acc: Vec<Option<ahash::AHashMap<u32, f64>>> =
        (0..n_targets).map(|_| None).collect();

    let mut weights: Vec<f64> = Vec::with_capacity(64);

    for (label, &count) in packed_map.iter_labels().zip(packed_map.counts.iter()) {
        let ec_count = count as f64;
        if ec_count == 0.0 {
            continue;
        }

        weights.clear();
        let mut denom = 0.0f64;
        for (tid, cond_prob) in label.target_labels().iter().zip(label.target_probs()) {
            let w = cond_prob * em_counts[*tid as usize] * inv_eff_lens[*tid as usize];
            weights.push(w);
            denom += w;
        }
        if denom <= 1e-8 {
            continue;
        }

        let targets = label.target_labels();
        let n = targets.len();
        if n < 2 {
            continue;
        }

        for i in 0..n {
            let t = targets[i] as usize;
            if weights[i] < 1e-10 {
                continue;
            }

            let acc = competitor_acc[t].get_or_insert_with(ahash::AHashMap::new);
            for j in 0..n {
                if i == j {
                    continue;
                }
                let u = targets[j] as u32;
                let u_share = ec_count * weights[j] / denom;
                if u_share > 0.01 {
                    *acc.entry(u).or_insert(0.0) += u_share;
                }
            }
        }
    }

    // For each transcript, find the dominant competitor and compute
    // neighborhood fraction. Use the top competitors that collectively
    // account for >= 90% of the competition (like a gene's top isoforms).
    for t in 0..n_targets {
        if em_counts[t] <= 0.0 {
            continue;
        }

        let acc = match &competitor_acc[t] {
            Some(a) if !a.is_empty() => a,
            _ => continue,
        };

        // Sort competitors by accumulated count (descending)
        let mut sorted: Vec<(u32, f64)> = acc.iter().map(|(&u, &c)| (u, c)).collect();
        sorted.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

        // Top competitor is the dominant
        dominant[t] = sorted[0].0 as usize;

        // Sum top competitors until we cover 90% of total competition
        let total_competition: f64 = sorted.iter().map(|(_, c)| c).sum();
        if total_competition <= 0.0 {
            continue;
        }

        let mut cumulative = 0.0f64;
        let mut cluster_count = 0.0f64;
        for &(_, c) in &sorted {
            cluster_count += c;
            cumulative += c;
            if cumulative >= 0.9 * total_competition {
                break;
            }
        }

        let total = em_counts[t] + cluster_count;
        if total > 0.0 {
            nbr_frac[t] = em_counts[t] / total;
        }
    }

    (nbr_frac, dominant)
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
        .map(|(&c, &keep)| {
            if keep && c.is_finite() && c > 0.0 {
                c
            } else {
                0.0
            }
        })
        .collect()
}

fn sorted_condition_names(samples: &[SampleEntry]) -> Vec<String> {
    let mut condition_names: Vec<String> = samples.iter().map(|s| s.condition.clone()).collect();
    condition_names.sort();
    condition_names.dedup();
    condition_names
}

fn sample_condition_indices(samples: &[SampleEntry], condition_names: &[String]) -> Vec<usize> {
    samples
        .iter()
        .map(|sample| {
            condition_names
                .iter()
                .position(|c| c == &sample.condition)
                .expect("sample condition should exist in condition_names")
        })
        .collect()
}

/// Phase-1 EM dispatch over any `PackedEqMap<L>`. `allow_coverage_smoothing`
/// must be false when `packed` is a collapsed map (no positional bins).
fn phase1_em_step<L: EqLabel>(
    packed: &PackedEqMap<L>,
    eff_lengths: Vec<f64>,
    opts: &ConsensusQuantOpts,
    inner_threads: usize,
    serial_inner_pool: Option<&rayon::ThreadPool>,
    allow_coverage_smoothing: bool,
) -> Vec<f64> {
    let eminfo = EMInfo::new(
        packed,
        eff_lengths,
        phase1_max_iter(opts),
        phase1_convergence_thresh(opts),
        opts.presence_thresh,
    );
    if allow_coverage_smoothing && opts.coverage_smooth_rounds > 0 && opts.pos_bins > 1 {
        em_with_coverage(
            &eminfo,
            None,
            opts.pos_bins as usize,
            opts.coverage_smooth_rounds,
            opts.coverage_epsilon,
        )
    } else if !opts.no_phase1_squarem {
        if let Some(pool) = serial_inner_pool {
            squarem_em_par_with_pool(&eminfo, pool)
        } else if inner_threads > 1 {
            squarem_em_par(&eminfo, inner_threads)
        } else {
            squarem_em(&eminfo)
        }
    } else if let Some(pool) = serial_inner_pool {
        em_par_with_pool(&eminfo, pool)
    } else if inner_threads > 1 {
        em_par(&eminfo, inner_threads)
    } else {
        em(&eminfo)
    }
}

/// Phase-2 EM dispatch over any `PackedEqMap<L>`. Always uses standard or
/// penalized (non-SQUAREM) EM and requires a transcript-level mask to be
/// applied via `EMInfo::apply_mask` before the EM call.
#[allow(clippy::too_many_arguments)]
fn phase2_em_step<L: EqLabel>(
    packed: &PackedEqMap<L>,
    eff_lengths: Vec<f64>,
    phase2_mask: &[bool],
    init: Option<&[f64]>,
    alpha: Option<&[f64]>,
    opts: &ConsensusQuantOpts,
    inner_threads: usize,
    serial_inner_pool: Option<&rayon::ThreadPool>,
) -> Vec<f64> {
    let mut eminfo = EMInfo::new(
        packed,
        eff_lengths,
        phase2_max_iter(opts),
        phase2_convergence_thresh(opts),
        opts.presence_thresh,
    );
    eminfo.apply_mask(phase2_mask);
    if let Some(alpha) = alpha {
        if let Some(pool) = serial_inner_pool {
            em_penalized_par_with_pool_init(&eminfo, alpha, init, pool)
        } else if inner_threads > 1 {
            em_penalized_par_init(&eminfo, alpha, init, inner_threads)
        } else {
            em_penalized_init(&eminfo, alpha, init)
        }
    } else if let Some(pool) = serial_inner_pool {
        em_par_with_pool_init(&eminfo, init, pool)
    } else if inner_threads > 1 {
        em_par_init(&eminfo, init, inner_threads)
    } else {
        em_init(&eminfo, init)
    }
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
                    info!(
                        "  [{}/{}] Building EQ map for {}",
                        i + 1,
                        n_samples,
                        sample.sample_name
                    );
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
            info!(
                "  [{}/{}] Building EQ map for {}",
                i + 1,
                n_samples,
                sample.sample_name
            );
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

    // Build merged transcript-EC index (needed for structural selection and
    // gene-fraction filter EC uniqueness).
    use crate::utils::txp_selection;
    let indices: Vec<_> = bundles
        .iter()
        .map(|b| txp_selection::TranscriptEqIndex::from_packed_eq_map(&b.packed_eq_map, n_targets))
        .collect();
    let eqc_counts: Vec<usize> = bundles.iter().map(|b| b.packed_eq_map.len()).collect();
    let total_eqcs: usize = eqc_counts.iter().sum();
    let merged_index = txp_selection::merge_transcript_indices(&indices, &eqc_counts, n_targets);

    // ====== Optional: structural transcript variable selection ======
    let selection_stats = if opts.txp_selection {
        info!("Running structural transcript variable selection on merged EC graph...");
        let stages = opts.selection_stages.clone().unwrap_or_default();
        let packed_maps_ref: Vec<&_> = bundles.iter().map(|b| &b.packed_eq_map).collect();
        let result = txp_selection::run_selection_from_index_with_coverage(
            &merged_index,
            n_targets,
            total_eqcs,
            &stages,
            &packed_maps_ref,
            &eqc_counts,
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

    // ====== Optional: collapsed EC view for EM/SQUAREM ======
    // Collapse positional ECs that share the same (targets, prob_bins) key.
    // Semantics-preserving for the M-step (which reads only target_labels()
    // and target_probs()). We gate on RangeFactorizedEqLabel (the only
    // label type whose positional bins can be stripped), the
    // `--no-collapsed-ec-em` override, pos_bins > 1 (otherwise nothing to
    // collapse), and whether coverage-smoothing EM is active (it needs the
    // positional map and hence blocks collapse).
    let cov_smoothing_active = opts.coverage_smooth_rounds > 0 && opts.pos_bins > 1;
    let is_range_factorized = std::any::TypeId::of::<EqLabelT>()
        == std::any::TypeId::of::<RangeFactorizedEqLabel>();
    let use_collapsed = is_range_factorized
        && !opts.no_collapsed_ec_em
        && opts.pos_bins > 1
        && !cov_smoothing_active;
    let collapsed_maps: Vec<Option<CollapsedEqMap>> = if use_collapsed {
        info!(
            "building collapsed EC views for {} sample{} (EM will iterate over the collapsed map)",
            n_samples,
            if n_samples == 1 { "" } else { "s" }
        );
        bundles
            .iter()
            .map(|b| {
                // SAFETY: `is_range_factorized` above verified EqLabelT == RangeFactorizedEqLabel.
                let pos_map: &PackedEqMap<RangeFactorizedEqLabel> = unsafe {
                    &*(&b.packed_eq_map as *const PackedEqMap<EqLabelT>
                        as *const PackedEqMap<RangeFactorizedEqLabel>)
                };
                Some(build_collapsed(pos_map))
            })
            .collect()
    } else {
        (0..n_samples).map(|_| None).collect()
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
                    info!(
                        "  [{}/{}] Phase 1 EM for {}",
                        i + 1,
                        n_samples,
                        samples[i].sample_name
                    );
                    let counts = if let Some(cm) = collapsed_maps[i].as_ref() {
                        phase1_em_step(
                            &cm.packed,
                            bundle.eff_lengths.clone(),
                            opts,
                            inner_threads,
                            None,
                            false,
                        )
                    } else {
                        phase1_em_step(
                            &bundle.packed_eq_map,
                            bundle.eff_lengths.clone(),
                            opts,
                            inner_threads,
                            None,
                            true,
                        )
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
            info!(
                "  [{}/{}] Phase 1 EM for {}",
                i + 1,
                n_samples,
                samples[i].sample_name
            );
            let counts = if let Some(cm) = collapsed_maps[i].as_ref() {
                phase1_em_step(
                    &cm.packed,
                    bundle.eff_lengths.clone(),
                    opts,
                    inner_threads,
                    serial_inner_pool.as_ref(),
                    false,
                )
            } else {
                phase1_em_step(
                    &bundle.packed_eq_map,
                    bundle.eff_lengths.clone(),
                    opts,
                    inner_threads,
                    serial_inner_pool.as_ref(),
                    true,
                )
            };
            counts_vec.push(counts);
        }
        counts_vec
    };

    let phase1_tpms: Vec<Vec<f64>> = phase1_counts
        .iter()
        .enumerate()
        .map(|(i, counts)| compute_tpm(counts, &bundles[i].eff_lengths))
        .collect();

    // ====== Consensus filter ======
    let filter_mode = &opts.filter_mode;
    let base_min_ecs = opts.min_ec_support;

    // Compute per-transcript adaptive EC support threshold if requested.
    let adaptive_thresh: Option<Vec<u32>> = if opts.adaptive_ec_support
        && matches!(filter_mode, FilterMode::Support | FilterMode::Hybrid)
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

    // Compute per-sample evidence once and reuse it for global consensus,
    // condition rescue, and gene-level rescue masks.
    let sample_pass_masks: Vec<Vec<bool>> = phase1_counts
        .iter()
        .enumerate()
        .map(|(i, counts)| {
            sample_expression_mask(
                &bundles[i].packed_eq_map,
                counts,
                &bundles[i].eff_lengths,
                filter_mode,
                opts,
                base_min_ecs,
                adaptive_thresh.as_deref(),
            )
        })
        .collect();

    let rescue_sample_pass_masks: Vec<Vec<bool>> = phase1_tpms
        .iter()
        .map(|tpms| sample_tpm_rescue_mask(tpms, opts.tpm_threshold))
        .collect();

    let mut express_count = vec![0u32; n_targets];
    for sample_pass in &sample_pass_masks {
        for (t, &pass) in sample_pass.iter().enumerate() {
            if pass {
                express_count[t] += 1;
            }
        }
    }

    let filter_label = match filter_mode {
        FilterMode::Tpm => format!("TPM > {}", opts.tpm_threshold),
        FilterMode::Ues => format!("UES > {}", opts.ues_threshold),
        FilterMode::Support => format!(
            "EC support >= {}{}",
            base_min_ecs,
            if opts.adaptive_ec_support {
                " (adaptive)"
            } else {
                ""
            }
        ),
        FilterMode::Hybrid => format!(
            "hybrid support/UES (base_support={}, ues_margin>{}{})",
            base_min_ecs,
            opts.ues_threshold,
            if opts.adaptive_ec_support {
                ", adaptive"
            } else {
                ""
            }
        ),
    };

    // Determine K threshold
    let min_fraction = opts.min_fraction.unwrap_or_else(|| {
        // Default: (N-1)/N
        (n_samples as f64 - 1.0) / n_samples as f64
    });
    let min_k = ((min_fraction * n_samples as f64).ceil() as u32).max(1);

    // Auto-enable condition rescue when multiple conditions are present,
    // unless explicitly disabled or condition-aware mode is selected.
    let has_conditions = samples.iter().any(|s| s.condition != samples[0].condition);
    let use_condition_rescue = if opts.condition_aware_consensus || opts.no_condition_rescue {
        false
    } else if opts.condition_rescue {
        true
    } else {
        // Auto-enable when multiple conditions exist
        has_conditions
    };
    if use_condition_rescue && !opts.condition_rescue {
        info!(
            "Auto-enabling condition rescue (multiple conditions detected; use --no-condition-rescue to disable)"
        );
    }

    let condition_data =
        if has_conditions && (opts.condition_aware_consensus || use_condition_rescue) {
            let mut condition_names: Vec<String> =
                samples.iter().map(|s| s.condition.clone()).collect();
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

                let condition_pass_mask = if opts.condition_aware_consensus {
                    &sample_pass_masks[sample_idx]
                } else {
                    &rescue_sample_pass_masks[sample_idx]
                };

                for (t, &pass) in condition_pass_mask.iter().enumerate() {
                    if pass {
                        cond_counts[ci][t] += 1;
                    }
                }
            }

            let cond_k: Vec<u32> = cond_reps
                .iter()
                .map(|&nrep| ((min_fraction * nrep as f64).ceil() as u32).max(1))
                .collect();

            let cond_support_masks: Vec<Vec<bool>> = cond_counts
                .iter()
                .enumerate()
                .map(|(ci, counts)| counts.iter().map(|&c| c >= cond_k[ci]).collect())
                .collect();

            Some((condition_names, cond_counts, cond_k, cond_support_masks))
        } else {
            None
        };

    let strict_global_mask: Vec<bool> = express_count.iter().map(|&c| c >= min_k).collect();

    let pre_gene_consensus_mask: Vec<bool> = if opts.condition_aware_consensus {
        let (condition_names, cond_counts, cond_k, _) = condition_data
            .as_ref()
            .expect("condition data required for condition-aware consensus");
        let mut mask = vec![false; n_targets];
        for t in 0..n_targets {
            mask[t] = (0..condition_names.len()).any(|ci| cond_counts[ci][t] >= cond_k[ci]);
        }
        info!(
            "Condition-aware consensus enabled across {} conditions",
            condition_names.len()
        );
        mask
    } else if use_condition_rescue {
        // Strict global consensus + condition-specific rescue.
        let (condition_names, cond_counts, cond_k, _) = condition_data
            .as_ref()
            .expect("condition data required for condition-rescue");

        let mut mask = vec![false; n_targets];
        let mut n_global = 0usize;
        let mut n_rescued = 0usize;
        for t in 0..n_targets {
            if strict_global_mask[t] {
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
            "Condition rescue: {} pass global, {} rescued from within-condition TPM consensus ({} conditions, threshold={})",
            n_global,
            n_rescued,
            condition_names.len(),
            opts.tpm_threshold
        );
        mask
    } else {
        strict_global_mask.clone()
    };

    let n_consensus = pre_gene_consensus_mask.iter().filter(|&&b| b).count();
    let n_filtered = n_targets - n_consensus;

    info!(
        "Consensus filter ({}{}): K={} (min_fraction={:.2}), {} transcripts pass, {} filtered out",
        filter_label,
        if opts.condition_aware_consensus {
            ", condition-aware"
        } else if use_condition_rescue {
            ", condition-rescue"
        } else {
            ""
        },
        min_k,
        min_fraction,
        n_consensus,
        n_filtered
    );

    // Parse gene names from GENCODE-style pipe-delimited transcript names
    // (field 6, 0-indexed field 5). Only used when --use-gene-annotation is set.
    // Default: EC-graph-based grouping (annotation-free, generally more accurate).
    let gene_names: Vec<Option<&str>> = if opts.use_gene_annotation {
        bundles[0]
            .ref_names
            .iter()
            .map(|name| name.split('|').nth(5))
            .collect()
    } else {
        vec![None; n_targets]
    };

    let gene_to_txps: std::collections::HashMap<&str, Vec<usize>> = {
        let mut map: std::collections::HashMap<&str, Vec<usize>> = std::collections::HashMap::new();
        for (t, gene) in gene_names.iter().enumerate() {
            if let Some(g) = gene {
                map.entry(g).or_default().push(t);
            }
        }
        map
    };

    // ====== Gene-level rescue ======
    // For transcripts that fail consensus, check if their gene's total TPM
    // passes consensus. If so, rescue all transcripts of that gene.
    let consensus_mask = {
        let mut mask = pre_gene_consensus_mask.clone();

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

    let final_n_consensus = consensus_mask.iter().filter(|&&b| b).count();
    let phase2_sample_masks = sample_specific_phase2_masks(
        samples,
        &strict_global_mask,
        &pre_gene_consensus_mask,
        condition_data.as_ref().map(|(names, _, _, _)| names.as_slice()),
        condition_data
            .as_ref()
            .map(|(_, _, _, support_masks)| support_masks.as_slice()),
        opts.condition_aware_consensus,
        use_condition_rescue,
    );
    let phase2_prior_alpha: Option<Vec<Vec<f64>>> =
        if opts.condition_specific_prior_weight > 0.0 {
        let avg_total_reads: f64 = bundles
            .iter()
            .map(|b| b.packed_eq_map.total_weight() as f64)
            .sum::<f64>()
            / n_samples as f64;
        let alpha_0 = opts.condition_specific_prior_weight * avg_total_reads;
        let condition_names = sorted_condition_names(samples);
        let condition_indices = sample_condition_indices(samples, &condition_names);
        let mut hyperparams =
            hierarchical::init_hyperparams(n_targets, condition_names, alpha_0);
        let phase1_results: Vec<hierarchical::SampleResult> = samples
            .iter()
            .enumerate()
            .map(|(i, sample)| hierarchical::SampleResult {
                counts: phase1_counts[i].clone(),
                present: phase1_counts[i]
                    .iter()
                    .zip(phase2_sample_masks[i].iter())
                    .map(|(&count, &keep)| keep && count > opts.presence_thresh)
                    .collect(),
                sample_name: sample.sample_name.clone(),
                condition_idx: condition_indices[i],
            })
            .collect();
        hierarchical::update_condition_means(&phase1_results, &mut hyperparams);
        let alphas: Vec<Vec<f64>> = samples
            .iter()
            .enumerate()
            .map(|(i, _)| {
                hierarchical::compute_pseudo_counts(
                    &hyperparams,
                    condition_indices[i],
                    &phase2_sample_masks[i],
                )
            })
            .collect();
        info!(
            "Phase 2 condition-specific hierarchical prior enabled: prior_weight={:.3}, alpha_0={:.1}",
            opts.condition_specific_prior_weight, alpha_0
        );
        Some(alphas)
    } else {
        None
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
                    let phase2_mask = &phase2_sample_masks[i];
                    let init = if opts.no_phase2_warm_start {
                        None
                    } else {
                        Some(phase2_init_counts(&phase1_counts[i], phase2_mask))
                    };
                    let alpha = phase2_prior_alpha.as_ref().map(|a| a[i].as_slice());
                    let em_res = if let Some(cm) = collapsed_maps[i].as_ref() {
                        phase2_em_step(
                            &cm.packed,
                            bundles[i].eff_lengths.clone(),
                            phase2_mask,
                            init.as_deref(),
                            alpha,
                            opts,
                            inner_threads,
                            None,
                        )
                    } else {
                        phase2_em_step(
                            &bundles[i].packed_eq_map,
                            bundles[i].eff_lengths.clone(),
                            phase2_mask,
                            init.as_deref(),
                            alpha,
                            opts,
                            inner_threads,
                            None,
                        )
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
            let phase2_mask = &phase2_sample_masks[i];
            let init = if opts.no_phase2_warm_start {
                None
            } else {
                Some(phase2_init_counts(&phase1_counts[i], phase2_mask))
            };
            let alpha = phase2_prior_alpha.as_ref().map(|a| a[i].as_slice());
            let em_res = if let Some(cm) = collapsed_maps[i].as_ref() {
                phase2_em_step(
                    &cm.packed,
                    bundles[i].eff_lengths.clone(),
                    phase2_mask,
                    init.as_deref(),
                    alpha,
                    opts,
                    inner_threads,
                    serial_inner_pool.as_ref(),
                )
            } else {
                phase2_em_step(
                    &bundles[i].packed_eq_map,
                    bundles[i].eff_lengths.clone(),
                    phase2_mask,
                    init.as_deref(),
                    alpha,
                    opts,
                    inner_threads,
                    serial_inner_pool.as_ref(),
                )
            };
            results.push((i, em_res));
        }
        results
    };

    // ====== Gene-fraction filter ======
    // Precompute robust EC uniqueness per transcript within each gene.
    // A transcript's unique ECs are those not shared with any other
    // consensus transcript of the same gene. We weight uniqueness by
    // cross-sample consistency: unique ECs appearing in more samples
    // are stronger evidence of independent expression.
    //
    // The merged index's EC IDs encode sample origin via offset ranges.
    // ECs from sample s have IDs in [eqc_offsets[s], eqc_offsets[s+1]).
    let eqc_offset_boundaries: Vec<u32> = {
        let mut bounds = Vec::with_capacity(eqc_counts.len() + 1);
        let mut cum = 0u32;
        bounds.push(cum);
        for &c in &eqc_counts {
            cum += c as u32;
            bounds.push(cum);
        }
        bounds
    };
    let n_samples_total = eqc_counts.len();

    // Determine which sample an EC ID belongs to.
    let ec_to_sample = |eqc_id: u32| -> usize {
        eqc_offset_boundaries
            .partition_point(|&b| b <= eqc_id)
            .saturating_sub(1)
    };

    let ec_unique_frac: Vec<f64> = if !gene_to_txps.is_empty() {
        let mut frac = vec![0.0f64; n_targets];
        for txps in gene_to_txps.values() {
            let active: Vec<usize> = txps
                .iter()
                .filter(|&&t| consensus_mask[t])
                .copied()
                .collect();
            if active.len() <= 1 {
                for &t in &active {
                    frac[t] = 1.0;
                }
                continue;
            }
            for &t in &active {
                let sig_t = merged_index.signature(t);
                let n_total = sig_t.len();
                if n_total == 0 {
                    continue;
                }
                // Count unique ECs weighted by cross-sample robustness.
                // A unique EC counts as 1.0 only if it appears in multiple
                // samples. If it appears in only 1 sample, it counts as
                // 1/n_samples (fragile evidence).
                let mut unique_weight = 0.0f64;
                let mut total_weight = n_total as f64;
                let mut prev_eqc = u32::MAX;
                let mut prev_sample = usize::MAX;
                let mut unique_samples: Vec<usize> = Vec::new();

                for &eqc in sig_t {
                    let is_shared = active.iter().any(|&u| {
                        u != t && {
                            let sig_u = merged_index.signature(u);
                            sig_u.binary_search(&eqc).is_ok()
                        }
                    });
                    if !is_shared {
                        let sample = ec_to_sample(eqc);
                        unique_samples.push(sample);
                    }
                }

                if !unique_samples.is_empty() {
                    // Count distinct samples contributing unique ECs.
                    unique_samples.sort_unstable();
                    unique_samples.dedup();
                    let n_unique_samples = unique_samples.len();

                    // Robust uniqueness: require unique ECs in >= 2 samples
                    // (or >= 50% of samples for small N) to count as full
                    // evidence. Single-sample unique ECs get partial credit.
                    let min_robust = 2.min(n_samples_total);
                    if n_unique_samples >= min_robust {
                        frac[t] = 1.0; // robust unique evidence
                    } else {
                        // Fragile: unique ECs in only 1 sample
                        frac[t] = 0.5 / n_samples_total as f64;
                    }
                }
                // else: frac[t] remains 0.0 (no unique ECs)
            }
        }
        frac
    } else {
        vec![0.0f64; n_targets]
    };

    // Log EC uniqueness stats for expressed transcripts.
    {
        let active_fracs: Vec<f64> = (0..n_targets)
            .filter(|&t| consensus_mask[t])
            .map(|t| ec_unique_frac[t])
            .collect();
        if !active_fracs.is_empty() {
            let n_no_unique = active_fracs.iter().filter(|&&f| f == 0.0).count();
            let n_some_unique = active_fracs.iter().filter(|&&f| f > 0.0 && f < 1.0).count();
            let n_all_unique = active_fracs.iter().filter(|&&f| f == 1.0).count();
            info!(
                "EC uniqueness: {} consensus transcripts — {} with no unique ECs, {} with some, {} fully unique",
                active_fracs.len(),
                n_no_unique,
                n_some_unique,
                n_all_unique
            );
        }
    }

    // Diagnostic: dump pairwise EC graph metrics for expressed transcript pairs.
    if std::env::var("PISCEM_EC_GRAPH_DIAG").is_ok() {
        info!("Computing pairwise EC graph metrics (this may take a while)...");
        // Build inverted index: EC → list of consensus transcripts containing it
        let max_eqc = merged_index
            .eqc_ids
            .iter()
            .copied()
            .max()
            .map(|m| m as usize + 1)
            .unwrap_or(0);
        let mut eqc_to_txps: Vec<Vec<u32>> = vec![Vec::new(); max_eqc];
        for t in 0..n_targets {
            if !consensus_mask[t] {
                continue;
            }
            for &eqc in merged_index.signature(t) {
                eqc_to_txps[eqc as usize].push(t as u32);
            }
        }

        // For each pair of consensus transcripts sharing >= 1 EC, compute:
        // - Jaccard similarity of EC signatures
        // - Fraction of t's ECs shared with u (asymmetric containment)
        // - Mean EC size of shared ECs
        // - Whether they're in the same gene
        let mut pair_metrics: ahash::AHashMap<(u32, u32), (u32, u32, u32, f64)> =
            ahash::AHashMap::new();
        // (shared_count, sig_t_size, sig_u_size, sum_shared_ec_sizes)

        for t in 0..n_targets {
            if !consensus_mask[t] {
                continue;
            }
            let sig_t = merged_index.signature(t);
            if sig_t.is_empty() {
                continue;
            }

            // Find all transcripts sharing at least one EC with t
            let mut neighbor_shared: ahash::AHashMap<u32, (u32, f64)> = ahash::AHashMap::new();
            for &eqc in sig_t {
                let ec_size = eqc_to_txps[eqc as usize].len() as f64;
                for &u in &eqc_to_txps[eqc as usize] {
                    if u as usize == t {
                        continue;
                    }
                    let entry = neighbor_shared.entry(u).or_insert((0, 0.0));
                    entry.0 += 1;
                    entry.1 += ec_size;
                }
            }

            for (&u, &(shared, sum_ec_size)) in &neighbor_shared {
                let key = if (t as u32) < u {
                    (t as u32, u)
                } else {
                    (u, t as u32)
                };
                let sig_u_len = merged_index.degree(u as usize) as u32;
                pair_metrics.entry(key).or_insert((
                    shared,
                    sig_t.len() as u32,
                    sig_u_len,
                    sum_ec_size / shared as f64, // mean shared EC size
                ));
            }
        }

        // Dump to file
        if let Ok(mut f) = std::fs::File::create("ec_graph_pairs.tsv") {
            use std::io::Write;
            writeln!(f, "txp_t\ttxp_u\tgene_t\tgene_u\tsame_gene\tshared_ecs\tsig_t_size\tsig_u_size\tjaccard\tcontain_t\tcontain_u\tmean_shared_ec_size").ok();
            for (&(t, u), &(shared, sig_t, sig_u, mean_ec_sz)) in &pair_metrics {
                let gene_t = gene_names.get(t as usize).and_then(|g| *g).unwrap_or("NA");
                let gene_u = gene_names.get(u as usize).and_then(|g| *g).unwrap_or("NA");
                let same = if gene_t != "NA" && gene_t == gene_u {
                    1
                } else {
                    0
                };
                let union = sig_t + sig_u - shared;
                let jaccard = if union > 0 {
                    shared as f64 / union as f64
                } else {
                    0.0
                };
                let contain_t = if sig_t > 0 {
                    shared as f64 / sig_t as f64
                } else {
                    0.0
                };
                let contain_u = if sig_u > 0 {
                    shared as f64 / sig_u as f64
                } else {
                    0.0
                };
                writeln!(
                    f,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{:.1}",
                    t,
                    u,
                    gene_t,
                    gene_u,
                    same,
                    shared,
                    sig_t,
                    sig_u,
                    jaccard,
                    contain_t,
                    contain_u,
                    mean_ec_sz
                )
                .ok();
            }
            info!(
                "Wrote {} pairwise EC metrics to ec_graph_pairs.tsv",
                pair_metrics.len()
            );
        }
    }

    // Build EC-graph-based groups for annotation-free leakage filtering.
    // Groups transcripts with Jaccard(EC signatures) >= 0.10 into connected components.
    let ec_graph_groups = build_ec_graph_groups(&merged_index, &consensus_mask, n_targets);

    // Build group → members mapping
    let ec_graph_group_members: std::collections::HashMap<u32, Vec<usize>> = {
        let mut map: std::collections::HashMap<u32, Vec<usize>> = std::collections::HashMap::new();
        for t in 0..n_targets {
            if !consensus_mask[t] {
                continue;
            }
            map.entry(ec_graph_groups[t]).or_default().push(t);
        }
        map
    };

    {
        let n_groups = ec_graph_group_members.len();
        let n_multi = ec_graph_group_members
            .values()
            .filter(|v| v.len() > 1)
            .count();
        let max_size = ec_graph_group_members
            .values()
            .map(|v| v.len())
            .max()
            .unwrap_or(0);
        info!(
            "EC-graph groups: {} total ({} multi-member, max size {})",
            n_groups, n_multi, max_size
        );
    }

    // Compute per-transcript position bin profile from the merged index.
    // For each transcript, count fragments in each position bin (count-weighted).
    let n_pos_bins = crate::utils::eq_maps::NUM_POS_BINS
        .get()
        .copied()
        .unwrap_or(1.0) as usize;
    let pos_bin_profiles: Vec<Vec<u64>> =
        if n_pos_bins > 1 && merged_index.has_pos_bins() && merged_index.has_counts() {
            let mut profiles = vec![vec![0u64; n_pos_bins]; n_targets];
            for t in 0..n_targets {
                if !consensus_mask[t] {
                    continue;
                }
                let sig = merged_index.signature(t);
                let base = merged_index.offsets[t] as usize;
                for (j, _) in sig.iter().enumerate() {
                    let pb = merged_index.pos_bins[base + j] as usize;
                    let cnt = merged_index.ec_counts[base + j] as u64;
                    if pb < n_pos_bins {
                        profiles[t][pb] += cnt;
                    }
                }
            }
            profiles
        } else {
            Vec::new()
        };

    // Compute position bin CV for each consensus transcript.
    // High CV -> non-uniform coverage -> possible leakage signal.
    let pos_cv: Vec<f64> = if !pos_bin_profiles.is_empty() && n_pos_bins > 1 {
        let mut cv = vec![0.0f64; n_targets];
        for t in 0..n_targets {
            if !consensus_mask[t] {
                continue;
            }
            let total: f64 = pos_bin_profiles[t].iter().map(|&c| c as f64).sum();
            if total < 1.0 {
                continue;
            }
            let mean = total / n_pos_bins as f64;
            let var: f64 = pos_bin_profiles[t]
                .iter()
                .map(|&c| {
                    let d = c as f64 - mean;
                    d * d
                })
                .sum::<f64>()
                / n_pos_bins as f64;
            cv[t] = var.sqrt() / (mean + 1e-10);
        }
        cv
    } else {
        vec![0.0f64; n_targets]
    };

    // Log position bin coverage stats if available.
    if !pos_bin_profiles.is_empty() && n_pos_bins > 1 {
        let mut n_full_coverage = 0usize;
        let mut n_partial = 0usize;
        let mut n_single = 0usize;
        for t in 0..n_targets {
            if !consensus_mask[t] {
                continue;
            }
            let occupied = pos_bin_profiles[t].iter().filter(|&&c| c > 0).count();
            match occupied {
                0 => {}
                1 => n_single += 1,
                2..=4 => n_partial += 1,
                _ => n_full_coverage += 1,
            }
        }
        info!(
            "Position coverage: {} full (all {} bins), {} partial, {} single-bin",
            n_full_coverage, n_pos_bins, n_partial, n_single
        );

        // Dump per-transcript diagnostic if enabled via environment variable.
        if std::env::var("PISCEM_POS_DIAG").is_ok() {
            if let Ok(mut f) = std::fs::File::create("pos_diagnostic.tsv") {
                use std::io::Write;
                writeln!(
                    f,
                    "target_name\tec_unique_frac\tpos_cv\t{}\ttotal_count",
                    (0..n_pos_bins)
                        .map(|b| format!("bin_{}", b))
                        .collect::<Vec<_>>()
                        .join("\t")
                )
                .ok();
                for t in 0..n_targets {
                    if !consensus_mask[t] {
                        continue;
                    }
                    let bins_str = pos_bin_profiles[t]
                        .iter()
                        .map(|c| c.to_string())
                        .collect::<Vec<_>>()
                        .join("\t");
                    let total: u64 = pos_bin_profiles[t].iter().sum();
                    writeln!(
                        f,
                        "{}\t{:.4}\t{:.4}\t{}\t{}",
                        bundles[0].ref_names[t], ec_unique_frac[t], pos_cv[t], bins_str, total
                    )
                    .ok();
                }
                info!("Wrote position diagnostic to pos_diagnostic.tsv");
            }
        }
    }

    let gene_frac_threshold = opts.gene_fraction_filter;
    let has_gene_annot = !gene_to_txps.is_empty();
    let has_pos = !pos_bin_profiles.is_empty() && n_pos_bins > 1;
    let mut total_gene_frac_removed = 0usize;
    let mut total_nbr_removed = 0usize;
    let gene_rescue_only: Vec<bool> = consensus_mask
        .iter()
        .zip(pre_gene_consensus_mask.iter())
        .map(|(&final_keep, &pre_gene_keep)| final_keep && !pre_gene_keep)
        .collect();
    let condition_rescue_only: Vec<bool> = pre_gene_consensus_mask
        .iter()
        .zip(strict_global_mask.iter())
        .map(|(&pre_gene_keep, &strict_keep)| pre_gene_keep && !strict_keep)
        .collect();
    for (i, mut em_res) in phase2_results {
        // Rescue-only transcripts retain their phase-1 estimates after phase 2.
        for t in 0..n_targets {
            if gene_rescue_only[t] || condition_rescue_only[t] {
                em_res[t] = phase1_counts[i][t];
            }
        }

        // ---- Gene-annotated filtering (when gene names are available) ----
        if has_gene_annot && (gene_frac_threshold > 0.0 || has_pos) {
            for txps in gene_to_txps.values() {
                let gene_total: f64 = txps.iter().map(|&t| em_res[t]).sum();
                if gene_total <= 0.0 {
                    continue;
                }

                let dominant = txps
                    .iter()
                    .filter(|&&t| em_res[t] > 0.0)
                    .max_by(|&&a, &&b| em_res[a].partial_cmp(&em_res[b]).unwrap());
                let dom_t = match dominant {
                    Some(&t) => t,
                    None => continue,
                };

                for &t in txps {
                    if t == dom_t || em_res[t] <= 0.0 {
                        continue;
                    }
                    let gene_frac = em_res[t] / gene_total;

                    // Gene-fraction filter
                    if gene_frac_threshold > 0.0 && gene_frac < gene_frac_threshold {
                        let pos_uneven = if has_pos { pos_cv[t] > 0.5 } else { false };
                        if ec_unique_frac[t] < 1.0 || pos_uneven {
                            em_res[t] = 0.0;
                            if i == 0 {
                                total_gene_frac_removed += 1;
                            }
                            continue;
                        }
                    }

                    // Profile-correlation filter (gene-annotated version)
                    if has_pos && gene_frac < 0.05 {
                        if profile_corr_is_leakage(&pos_bin_profiles, t, dom_t, n_pos_bins) {
                            em_res[t] = 0.0;
                            if i == 0 {
                                total_gene_frac_removed += 1;
                            }
                        }
                    }
                }
            }
        }

        // ---- EC-graph-group filtering (annotation-free) ----
        // For transcripts without gene annotation, reconstruct gene-like groups
        // from EC graph structure (Jaccard similarity >= 0.10) and apply the
        // same fraction + profile correlation filters within each group.
        {
            for (&group_root, members) in &ec_graph_group_members {
                let active: Vec<usize> = members
                    .iter()
                    .filter(|&&t| em_res[t] > 0.0)
                    .copied()
                    .collect();
                if active.len() <= 1 {
                    continue;
                }

                // Skip groups where all members have gene annotation
                // (already handled by gene-annotated filtering)
                if has_gene_annot && active.iter().all(|&t| gene_names[t].is_some()) {
                    continue;
                }

                let group_total: f64 = active.iter().map(|&t| em_res[t]).sum();
                if group_total <= 0.0 {
                    continue;
                }

                let dom_t = *active
                    .iter()
                    .max_by(|&&a, &&b| em_res[a].partial_cmp(&em_res[b]).unwrap())
                    .unwrap();

                for &t in &active {
                    if t == dom_t {
                        continue;
                    }
                    // Skip if this transcript was handled by gene-annotated filter
                    if has_gene_annot && gene_names[t].is_some() {
                        continue;
                    }

                    let group_frac = em_res[t] / group_total;

                    // Group-fraction filter (same logic as gene-fraction)
                    if group_frac < gene_frac_threshold.max(0.01) {
                        let pos_uneven = if has_pos { pos_cv[t] > 0.5 } else { false };
                        if ec_unique_frac[t] < 1.0 || pos_uneven {
                            em_res[t] = 0.0;
                            if i == 0 {
                                total_nbr_removed += 1;
                            }
                            continue;
                        }
                    }

                    // Profile-correlation filter (group version)
                    if has_pos && group_frac < 0.05 {
                        if profile_corr_is_leakage(&pos_bin_profiles, t, dom_t, n_pos_bins) {
                            em_res[t] = 0.0;
                            if i == 0 {
                                total_nbr_removed += 1;
                            }
                        }
                    }
                }
            }
        }

        // Diagnostic: on sample 0, compare EC-graph groups vs gene annotation
        if i == 0 && has_gene_annot {
            // Run EC-graph group filter on ALL transcripts (pretend no annotation)
            let p1 = &phase1_counts[0];
            let mut ecg_would_remove = 0usize;
            let mut gene_did_remove = 0usize;
            let mut both = 0usize;

            for (&_group_root, members) in &ec_graph_group_members {
                let active: Vec<usize> = members
                    .iter()
                    .filter(|&&t| p1[t] > 0.0 && consensus_mask[t])
                    .copied()
                    .collect();
                if active.len() <= 1 {
                    continue;
                }
                let group_total: f64 = active.iter().map(|&t| p1[t]).sum();
                if group_total <= 0.0 {
                    continue;
                }
                let dom_t = *active
                    .iter()
                    .max_by(|&&a, &&b| p1[a].partial_cmp(&p1[b]).unwrap())
                    .unwrap();

                for &t in &active {
                    if t == dom_t {
                        continue;
                    }
                    let gf = p1[t] / group_total;
                    let gene_removed = em_res[t] == 0.0;
                    let mut ecg_remove = false;

                    if gf < gene_frac_threshold.max(0.01) {
                        let pos_uneven = if has_pos { pos_cv[t] > 0.5 } else { false };
                        if ec_unique_frac[t] < 1.0 || pos_uneven {
                            ecg_remove = true;
                        }
                    }
                    if !ecg_remove && has_pos && gf < 0.05 {
                        if profile_corr_is_leakage(&pos_bin_profiles, t, dom_t, n_pos_bins) {
                            ecg_remove = true;
                        }
                    }

                    if gene_removed {
                        gene_did_remove += 1;
                    }
                    if ecg_remove {
                        ecg_would_remove += 1;
                    }
                    if gene_removed && ecg_remove {
                        both += 1;
                    }
                }
            }
            info!(
                "EC-graph vs gene filter: gene={}, ecg={}, overlap={}, ecg_only={}, gene_only={}",
                gene_did_remove,
                ecg_would_remove,
                both,
                ecg_would_remove.saturating_sub(both),
                gene_did_remove.saturating_sub(both)
            );
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
        .with_context(|| format!("failed to write quant output for {}", sample.sample_name))?;

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
            "num_consensus_targets": final_n_consensus,
            "consensus_min_k": min_k,
            "consensus_min_fraction": min_fraction,
            "filter_mode": format!("{}", filter_label),
            "condition_specific_prior_weight": opts.condition_specific_prior_weight,
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

    if total_gene_frac_removed > 0 || total_nbr_removed > 0 {
        info!(
            "Leakage filter: zeroed {} (gene-annotated) + {} (EC-neighborhood) isoforms in sample 1",
            total_gene_frac_removed, total_nbr_removed
        );
    }

    info!("Done. Wrote output for {} samples.", n_samples);
    Ok(())
}

/// Entry point for consensus-quant subcommand.
pub fn run(opts: &ConsensusQuantOpts) -> Result<()> {
    let samples = parse_manifest(&opts.manifest)?;
    info!("Parsed manifest with {} samples", samples.len());

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

#[cfg(test)]
mod tests {
    use super::{profile_corr_is_leakage, sample_specific_phase2_masks};
    use crate::multi_sample::SampleEntry;
    use std::path::PathBuf;

    fn sample(sample_name: &str, condition: &str) -> SampleEntry {
        SampleEntry {
            sample_name: sample_name.to_string(),
            condition: condition.to_string(),
            rad_path: PathBuf::from("dummy.rad"),
            output_dir: PathBuf::from("dummy_out"),
        }
    }

    #[test]
    fn phase2_masks_use_condition_specific_rescue() {
        let samples = vec![sample("s1", "A"), sample("s2", "B")];
        let strict_global_mask = vec![true, false, false];
        let pre_gene_consensus_mask = vec![true, true, false];
        let condition_names = vec!["A".to_string(), "B".to_string()];
        let condition_support_masks = vec![vec![true, true, false], vec![true, false, false]];

        let masks = sample_specific_phase2_masks(
            &samples,
            &strict_global_mask,
            &pre_gene_consensus_mask,
            Some(&condition_names),
            Some(&condition_support_masks),
            false,
            true,
        );

        assert_eq!(masks[0], vec![true, true, false]);
        assert_eq!(masks[1], vec![true, false, false]);
    }

    #[test]
    fn phase2_masks_fall_back_to_pre_gene_consensus_without_condition_logic() {
        let samples = vec![sample("s1", "A"), sample("s2", "B")];
        let strict_global_mask = vec![true, false];
        let pre_gene_consensus_mask = vec![true, true];

        let masks = sample_specific_phase2_masks(
            &samples,
            &strict_global_mask,
            &pre_gene_consensus_mask,
            None,
            None,
            false,
            false,
        );

        assert_eq!(masks, vec![vec![true, true], vec![true, true]]);
    }

    #[test]
    fn profile_corr_requires_peak_alignment() {
        let profiles = vec![vec![2, 12, 6, 0, 0], vec![0, 2, 12, 6, 0]];
        assert!(!profile_corr_is_leakage(&profiles, 0, 1, 5));
    }

    #[test]
    fn profile_corr_accepts_nearly_identical_profiles() {
        let profiles = vec![vec![1, 4, 12, 5, 1], vec![1, 5, 11, 6, 1]];
        assert!(profile_corr_is_leakage(&profiles, 0, 1, 5));
    }
}
