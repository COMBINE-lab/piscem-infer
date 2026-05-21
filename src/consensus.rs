//! Two-pass consensus-filtered multi-sample quantification.
//!
//! Phase 1: Build EQ maps and run per-sample EM to get initial abundance estimates.
//! Consensus filter: Keep transcripts expressed in ≥K of N samples.
//! Phase 2: Re-run per-sample EM with non-consensus transcripts masked out
//!          (effective length set to 0), redistributing reads to the consensus set.

use anyhow::{Context, Result, bail};
use path_tools::WithAdditionalExtension;
use rayon::prelude::*;
use serde_json::json;
use std::collections::HashSet;
use std::fs::{File, create_dir_all, read_to_string};
use std::io::Write;
use tracing::info;

use crate::multi_sample::{SampleEntry, parse_manifest};
use crate::process_rad::{EqMapBundle, RadProcessingOpts, build_eq_map_from_rad};
use crate::prog_opts::{
    ConditionRescueLockMode, ConsensusQuantOpts, FilterMode, PostFilterRedistributeMode,
    StructuredRescueRankMode,
};
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
    if use_condition_rescue && !condition_aware_consensus {
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
                .map(|t| cond_mask[t] || strict_global_mask[t])
                .collect()
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn condition_rescue_free_variable(
    active: bool,
    rescue_only: bool,
    locked: bool,
    phase1_count: f64,
    free_min_count: f64,
    floor_barrier_mode: bool,
    fully_locked: bool,
    lock_mode: &ConditionRescueLockMode,
) -> bool {
    if !active {
        return false;
    }
    if floor_barrier_mode {
        return true;
    }
    if rescue_only
        && matches!(
            lock_mode,
            ConditionRescueLockMode::GuardedConfidence
                | ConditionRescueLockMode::Instability
                | ConditionRescueLockMode::TranscriptStability
                | ConditionRescueLockMode::EnrichmentCredibleStability
        )
        && free_min_count > 0.0
        && phase1_count < free_min_count
    {
        return false;
    }
    if matches!(
        lock_mode,
        ConditionRescueLockMode::GuardedConfidence
            | ConditionRescueLockMode::Instability
            | ConditionRescueLockMode::TranscriptStability
            | ConditionRescueLockMode::EnrichmentCredibleStability
    ) && locked
    {
        return free_min_count <= 0.0 || phase1_count >= free_min_count;
    }
    !(fully_locked && locked)
}

fn condition_rescue_residual_caps(
    condition_rescue_only: &[bool],
    phase1_counts: &[f64],
    locked_counts: &[f64],
    active_mask: &[bool],
    z: f64,
) -> Vec<f64> {
    let z = z.max(0.0);
    condition_rescue_only
        .iter()
        .zip(phase1_counts.iter())
        .zip(locked_counts.iter())
        .zip(active_mask.iter())
        .map(
            |(((&rescue_only, &phase1_count), &locked_count), &active)| {
                if rescue_only && active {
                    let phase1_count = phase1_count.max(0.0);
                    let upper = phase1_count + z * phase1_count.max(1.0).sqrt();
                    (upper - locked_count.max(0.0)).max(0.0)
                } else {
                    f64::INFINITY
                }
            },
        )
        .collect()
}

#[inline]
fn condition_rescue_residual_cap_enabled(opts: &ConsensusQuantOpts) -> bool {
    opts.condition_rescue_residual_cap && !opts.no_condition_rescue_residual_cap
}

#[inline]
fn condition_rescue_residual_cap_active_set_enabled(opts: &ConsensusQuantOpts) -> bool {
    opts.condition_rescue_residual_cap_active_set
        && !opts.no_condition_rescue_residual_cap_active_set
}

#[inline]
fn lock_condition_rescue_allocations_enabled(
    opts: &ConsensusQuantOpts,
    use_condition_rescue: bool,
) -> bool {
    use_condition_rescue
        && !opts.no_lock_condition_rescue_allocations
        && !opts.reestimate_condition_rescue
}

#[inline]
fn sorted_subset(a: &[u32], b: &[u32]) -> bool {
    if a.len() > b.len() {
        return false;
    }
    let mut bi = 0usize;
    for &val in a {
        while bi < b.len() && b[bi] < val {
            bi += 1;
        }
        if bi >= b.len() || b[bi] != val {
            return false;
        }
        bi += 1;
    }
    true
}

fn write_selection_audit(
    path: &std::path::Path,
    target_names: &[String],
    result: &crate::utils::txp_selection::SelectionResult,
    index: &crate::utils::txp_selection::TranscriptEqIndex,
) -> Result<()> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        create_dir_all(parent)?;
    }

    let n_targets = target_names.len();
    let mut target_to_group = vec![usize::MAX; n_targets];
    for gi in 0..result.groups.num_groups() {
        for &member in result.groups.group_members(gi) {
            target_to_group[member as usize] = gi;
        }
    }

    let max_eqc = index
        .eqc_ids
        .iter()
        .copied()
        .max()
        .map(|m| m as usize + 1)
        .unwrap_or(0);
    let mut eqc_to_groups: Vec<Vec<u32>> = vec![Vec::new(); max_eqc];
    for gi in 0..result.groups.num_groups() {
        let rep = result.groups.representatives[gi] as usize;
        for &eqc in index.signature(rep) {
            eqc_to_groups[eqc as usize].push(gi as u32);
        }
    }

    let mut group_kept = vec![false; result.groups.num_groups()];
    for (t, &keep) in result.keep_mask.iter().enumerate() {
        let gi = target_to_group[t];
        if keep && gi != usize::MAX {
            group_kept[gi] = true;
        }
    }

    let mut file = File::create(path)?;
    writeln!(
        file,
        "target_idx\ttarget_name\tstructurally_kept\tselection_reason\tgroup_idx\tgroup_size\tgroup_representative_idx\tgroup_representative_name\tgroup_required\tgroup_dominated\tdegree\ttotal_ec_count\tdominator_group_idx\tdominator_idx\tdominator_name\tdominator_degree\tdominator_extra_eqcs\tshared_eqcs\tshared_ec_count_sum\tdominator_total_ec_count\tdominator_extra_ec_count_sum\tdominator_shared_count_fraction\tdominator_shared_pos_bins"
    )?;

    for (t, target_name) in target_names.iter().enumerate() {
        let keep = result.keep_mask[t];
        let degree = index.degree(t);
        let total_ec_count = index.total_count(t);
        let gi = target_to_group[t];

        let (reason, group_idx, group_size, rep_idx, rep_name, required, dominated, dominator) =
            if gi == usize::MAX {
                (
                    "no_eqc",
                    String::from("NA"),
                    0usize,
                    String::from("NA"),
                    String::from("NA"),
                    false,
                    false,
                    None,
                )
            } else {
                let group_dominated = result.dominated[gi];
                let reason = if keep {
                    "kept"
                } else if group_dominated {
                    "subset_dominated"
                } else {
                    "removed_unknown"
                };
                let rep = result.groups.representatives[gi] as usize;
                let dominator = if keep {
                    None
                } else {
                    find_selection_dominator(gi, &group_kept, &eqc_to_groups, result, index)
                };
                (
                    reason,
                    gi.to_string(),
                    result.groups.group_members(gi).len(),
                    rep.to_string(),
                    target_names[rep].clone(),
                    result.required[gi],
                    group_dominated,
                    dominator,
                )
            };

        let (
            dom_group_idx,
            dom_idx,
            dom_name,
            dom_degree,
            dom_extra_eqcs,
            shared_eqcs,
            shared_count_sum,
            dom_total_count,
            dom_extra_count,
            dom_shared_fraction,
            dom_bins,
        ) = if let Some(dom_gi) = dominator {
            let dom = result.groups.representatives[dom_gi] as usize;
            let sig = index.signature(t);
            let dom_sig = index.signature(dom);
            let dom_degree = dom_sig.len();
            let dom_extra = dom_degree.saturating_sub(sig.len());
            let shared_count_sum = index
                .counts_for(t)
                .map(|counts| counts.iter().map(|&c| c as u64).sum::<u64>())
                .unwrap_or(0);
            let dom_total_count = index.total_count(dom);
            let dom_extra_count = dom_total_count.saturating_sub(shared_count_sum);
            let dom_shared_fraction = if dom_total_count > 0 {
                shared_count_sum as f64 / dom_total_count as f64
            } else {
                0.0
            };
            let dom_bins = index
                .pos_bins_for(dom)
                .map(|bins| {
                    let mut used = std::collections::BTreeSet::new();
                    for &eqc in sig {
                        if let Ok(pos) = dom_sig.binary_search(&eqc) {
                            used.insert(bins[pos]);
                        }
                    }
                    used.into_iter()
                        .map(|b| b.to_string())
                        .collect::<Vec<_>>()
                        .join(",")
                })
                .unwrap_or_else(|| String::from("NA"));
            (
                dom_gi.to_string(),
                dom.to_string(),
                target_names[dom].clone(),
                dom_degree.to_string(),
                dom_extra.to_string(),
                sig.len().to_string(),
                shared_count_sum.to_string(),
                dom_total_count.to_string(),
                dom_extra_count.to_string(),
                format!("{dom_shared_fraction:.6}"),
                dom_bins,
            )
        } else {
            (
                String::from("NA"),
                String::from("NA"),
                String::from("NA"),
                String::from("NA"),
                String::from("NA"),
                String::from("NA"),
                String::from("NA"),
                String::from("NA"),
                String::from("NA"),
                String::from("NA"),
                String::from("NA"),
            )
        };

        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            t,
            target_name,
            keep,
            reason,
            group_idx,
            group_size,
            rep_idx,
            rep_name,
            required,
            dominated,
            degree,
            total_ec_count,
            dom_group_idx,
            dom_idx,
            dom_name,
            dom_degree,
            dom_extra_eqcs,
            shared_eqcs,
            shared_count_sum,
            dom_total_count,
            dom_extra_count,
            dom_shared_fraction,
            dom_bins
        )?;
    }

    Ok(())
}

fn find_selection_dominator(
    gi: usize,
    group_kept: &[bool],
    eqc_to_groups: &[Vec<u32>],
    result: &crate::utils::txp_selection::SelectionResult,
    index: &crate::utils::txp_selection::TranscriptEqIndex,
) -> Option<usize> {
    let rep_i = result.groups.representatives[gi] as usize;
    let sig_i = index.signature(rep_i);
    if sig_i.is_empty() {
        return None;
    }
    let rarest_eqc = sig_i
        .iter()
        .min_by_key(|&&eqc| eqc_to_groups[eqc as usize].len())?;

    eqc_to_groups[*rarest_eqc as usize]
        .iter()
        .map(|&candidate| candidate as usize)
        .filter(|&gj| gj != gi && group_kept[gj])
        .filter(|&gj| {
            let rep_j = result.groups.representatives[gj] as usize;
            let sig_j = index.signature(rep_j);
            sig_j.len() > sig_i.len() && sorted_subset(sig_i, sig_j)
        })
        .min_by_key(|&gj| {
            let rep_j = result.groups.representatives[gj] as usize;
            let sig_j = index.signature(rep_j);
            (sig_j.len() - sig_i.len(), sig_j.len(), rep_j)
        })
}

fn position_effective_bins(counts: &[u64]) -> f64 {
    let total: u64 = counts.iter().sum();
    if total == 0 {
        return 0.0;
    }
    let sum_sq: u128 = counts
        .iter()
        .map(|&count| {
            let count = count as u128;
            count * count
        })
        .sum();
    if sum_sq == 0 {
        0.0
    } else {
        (total as f64 * total as f64) / sum_sq as f64
    }
}

fn max_position_bin_fraction(counts: &[u64]) -> f64 {
    let total: u64 = counts.iter().sum();
    if total == 0 {
        return 0.0;
    }
    counts.iter().copied().max().unwrap_or(0) as f64 / total as f64
}

fn transcript_position_profile(
    index: &crate::utils::txp_selection::TranscriptEqIndex,
    t: usize,
    n_pos_bins: usize,
) -> Vec<u64> {
    let mut profile = vec![0u64; n_pos_bins];
    if n_pos_bins == 0 || !index.has_pos_bins() || !index.has_counts() {
        return profile;
    }
    let sig = index.signature(t);
    let base = index.offsets[t] as usize;
    for (j, _) in sig.iter().enumerate() {
        let pb = index.pos_bins[base + j] as usize;
        if pb < n_pos_bins {
            profile[pb] += index.ec_counts[base + j] as u64;
        }
    }
    profile
}

#[derive(Clone)]
struct StructuralRepairCandidate {
    target: usize,
    dominator: usize,
    group_idx: usize,
    shared_count: u64,
    dominator_total_count: u64,
    dominator_extra_count: u64,
    shared_fraction: f64,
}

fn local_pair_em<EqLabelT: EqLabel>(
    packed_eq_map: &PackedEqMap<EqLabelT>,
    merged_index: &crate::utils::txp_selection::TranscriptEqIndex,
    eqc_offsets: &[u32],
    sample_idx: usize,
    target: usize,
    related: usize,
    eff_lengths: &[f64],
) -> Option<(f64, f64)> {
    let target_eff_len = eff_lengths[target];
    let related_eff_len = eff_lengths[related];
    if target_eff_len <= 0.0 || related_eff_len <= 0.0 {
        return None;
    }

    let start = eqc_offsets[sample_idx];
    let end = eqc_offsets[sample_idx + 1];
    let mut eqcs: Vec<u32> = merged_index
        .signature(target)
        .iter()
        .chain(merged_index.signature(related).iter())
        .copied()
        .filter(|&eqc| eqc >= start && eqc < end)
        .collect();
    eqcs.sort_unstable();
    eqcs.dedup();
    if eqcs.is_empty() {
        return None;
    }

    let mut local_eqcs = Vec::new();
    let mut total_count = 0.0f64;
    for eqc in eqcs {
        let local_eqc = (eqc - start) as usize;
        let label = packed_eq_map.refs_for_eqc(local_eqc);
        let mut target_prob = 0.0f64;
        let mut related_prob = 0.0f64;
        for (tid, prob) in label.target_labels().iter().zip(label.target_probs()) {
            let tid = *tid as usize;
            if tid == target {
                target_prob = prob;
            } else if tid == related {
                related_prob = prob;
            }
        }
        if target_prob > 0.0 || related_prob > 0.0 {
            let count = packed_eq_map.counts[local_eqc] as f64;
            if count > 0.0 {
                total_count += count;
                local_eqcs.push((count, target_prob, related_prob));
            }
        }
    }
    if total_count <= 0.0 {
        return None;
    }

    let mut target_count = total_count * 0.5;
    let mut related_count = total_count * 0.5;
    for _ in 0..100 {
        let mut next_target = 0.0f64;
        let mut next_related = 0.0f64;
        for &(count, target_prob, related_prob) in &local_eqcs {
            let target_weight = target_prob * target_count / target_eff_len;
            let related_weight = related_prob * related_count / related_eff_len;
            let denom = target_weight + related_weight;
            if denom <= 0.0 {
                continue;
            }
            next_target += count * target_weight / denom;
            next_related += count * related_weight / denom;
        }
        let diff = (next_target - target_count).abs() + (next_related - related_count).abs();
        target_count = next_target;
        related_count = next_related;
        if diff / total_count.max(1.0) < 1e-5 {
            break;
        }
    }

    Some((target_count, related_count))
}

#[allow(clippy::too_many_arguments)]
fn condition_rescue_leakage_exempt<EqLabelT: EqLabel>(
    packed_eq_map: &PackedEqMap<EqLabelT>,
    merged_index: &crate::utils::txp_selection::TranscriptEqIndex,
    eqc_offsets: &[u32],
    sample_idx: usize,
    target: usize,
    dominator: usize,
    eff_lengths: &[f64],
    pos_bin_profiles: &[Vec<u64>],
    n_pos_bins: usize,
    min_count: f64,
    min_fraction: f64,
) -> Option<(f64, f64)> {
    if !pos_bin_profiles.is_empty()
        && profile_corr_is_leakage(pos_bin_profiles, target, dominator, n_pos_bins)
    {
        return None;
    }
    let (target_count, dominator_count) = local_pair_em(
        packed_eq_map,
        merged_index,
        eqc_offsets,
        sample_idx,
        target,
        dominator,
        eff_lengths,
    )?;
    let total = target_count + dominator_count;
    if total <= 0.0 {
        return None;
    }
    let fraction = target_count / total;
    if target_count >= min_count && fraction >= min_fraction {
        Some((target_count, fraction))
    } else {
        None
    }
}

fn signatures_intersect(
    index: &crate::utils::txp_selection::TranscriptEqIndex,
    a: usize,
    b: usize,
) -> bool {
    let sig_a = index.signature(a);
    let sig_b = index.signature(b);
    let mut ia = 0usize;
    let mut ib = 0usize;
    while ia < sig_a.len() && ib < sig_b.len() {
        match sig_a[ia].cmp(&sig_b[ib]) {
            std::cmp::Ordering::Equal => return true,
            std::cmp::Ordering::Less => ia += 1,
            std::cmp::Ordering::Greater => ib += 1,
        }
    }
    false
}

struct SharedResponsibilityContext<'a, L: EqLabel> {
    packed: &'a PackedEqMap<L>,
    index: &'a crate::utils::txp_selection::TranscriptEqIndex,
    eff_lengths: &'a [f64],
}

fn redistribute_by_shared_responsibility<L: EqLabel>(
    ctx: &SharedResponsibilityContext<'_, L>,
    em_res: &mut [f64],
    target: usize,
    candidates: &[usize],
    removed: f64,
    fallback_recipient: usize,
) {
    if !removed.is_finite() || removed <= 0.0 {
        return;
    }
    if target >= em_res.len() || target >= ctx.eff_lengths.len() || ctx.eff_lengths[target] <= 0.0 {
        if target != fallback_recipient {
            em_res[fallback_recipient] += removed;
        }
        return;
    }

    let mut local_candidate = vec![false; em_res.len()];
    for &candidate in candidates {
        if candidate != target
            && candidate < em_res.len()
            && candidate < ctx.eff_lengths.len()
            && em_res[candidate].is_finite()
            && em_res[candidate] > 0.0
            && ctx.eff_lengths[candidate] > 0.0
        {
            local_candidate[candidate] = true;
        }
    }

    let mut contributions: Vec<(usize, f64)> = Vec::new();
    let mut total_reassigned = 0.0;
    for &eq_idx_u32 in ctx.index.signature(target) {
        let eq_idx = eq_idx_u32 as usize;
        if eq_idx >= ctx.packed.len() {
            continue;
        }
        let label = ctx.packed.refs_for_eqc(eq_idx);
        let targets = label.target_labels();
        let probs = label.target_probs().collect::<Vec<_>>();
        let count = ctx.packed.counts[eq_idx] as f64;
        if count <= 0.0 || targets.is_empty() || probs.len() != targets.len() {
            continue;
        }

        let mut denom_all = 0.0;
        let mut target_weight = 0.0;
        let mut recipient_weights: Vec<(usize, f64)> = Vec::new();
        let mut denom_recipients = 0.0;
        for (&tid_u32, &prob) in targets.iter().zip(probs.iter()) {
            let tid = tid_u32 as usize;
            if tid >= em_res.len() || tid >= ctx.eff_lengths.len() || ctx.eff_lengths[tid] <= 0.0 {
                continue;
            }
            let weight = prob * em_res[tid].max(0.0) / ctx.eff_lengths[tid];
            if !weight.is_finite() || weight <= 0.0 {
                continue;
            }
            denom_all += weight;
            if tid == target {
                target_weight += weight;
            } else if local_candidate[tid] {
                recipient_weights.push((tid, weight));
                denom_recipients += weight;
            }
        }
        if denom_all <= 0.0 || target_weight <= 0.0 || denom_recipients <= 0.0 {
            continue;
        }

        let target_assignment = count * target_weight / denom_all;
        if !target_assignment.is_finite() || target_assignment <= 0.0 {
            continue;
        }
        for (recipient, weight) in recipient_weights {
            let delta = target_assignment * weight / denom_recipients;
            if delta > 0.0 && delta.is_finite() {
                if let Some((_t, existing)) = contributions
                    .iter_mut()
                    .find(|(existing_recipient, _)| *existing_recipient == recipient)
                {
                    *existing += delta;
                } else {
                    contributions.push((recipient, delta));
                }
                total_reassigned += delta;
            }
        }
    }

    if total_reassigned > 0.0 {
        let scale = removed / total_reassigned;
        for (recipient, delta) in contributions {
            em_res[recipient] += delta * scale;
        }
    } else if target != fallback_recipient {
        em_res[fallback_recipient] += removed;
    }
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
    for (t, &active) in consensus_mask.iter().enumerate().take(n_targets) {
        if !active {
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
    for (t, &active) in consensus_mask.iter().enumerate().take(n_targets) {
        if !active {
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

struct StructuredRescueResult {
    mask: Vec<bool>,
    condition_support_masks: Vec<Vec<bool>>,
    raw_rescued: usize,
    selected_rescued: usize,
    admitted_groups: usize,
}

#[derive(Default)]
struct RescueAuditRecord {
    target: usize,
    group: u32,
    group_size: usize,
    has_strict_member: bool,
    admitted: bool,
    raw_rescue: bool,
    selected: bool,
    primary_supported: bool,
    score: f64,
    score_frac: f64,
    rank: Option<usize>,
    peak_rank: Option<usize>,
    breadth_rank: Option<usize>,
    support_conditions: usize,
    mean_score: f64,
    suppression_reason: String,
    raw_conditions: Vec<usize>,
    admitted_conditions: Vec<usize>,
    selected_conditions: Vec<usize>,
    condition_mean_tpms: Vec<f64>,
    group_condition_mean_tpms: Vec<f64>,
}

struct StructuredRescueInputs<'a> {
    samples: &'a [SampleEntry],
    condition_names: &'a [String],
    cond_k: &'a [u32],
    raw_condition_support_masks: &'a [Vec<bool>],
    sample_pass_masks: &'a [Vec<bool>],
    strict_global_mask: &'a [bool],
    phase1_tpms: &'a [Vec<f64>],
    index: &'a crate::utils::txp_selection::TranscriptEqIndex,
    ref_names: &'a [String],
    opts: &'a ConsensusQuantOpts,
}

struct RescueCandidate {
    target: usize,
    score: f64,
    primary_supported: bool,
    support_conditions: usize,
    mean_score: f64,
    raw_conditions: Vec<usize>,
    condition_mean_tpms: Vec<f64>,
    group_condition_mean_tpms: Vec<f64>,
    peak_rank: Option<usize>,
    breadth_rank: Option<usize>,
}

fn structured_condition_rescue_mask(input: StructuredRescueInputs<'_>) -> StructuredRescueResult {
    let StructuredRescueInputs {
        samples,
        condition_names,
        cond_k,
        raw_condition_support_masks,
        sample_pass_masks,
        strict_global_mask,
        phase1_tpms,
        index,
        ref_names,
        opts,
    } = input;
    let n_targets = strict_global_mask.len();
    let n_conditions = condition_names.len();
    let max_isoforms = opts.structured_rescue_max_isoforms as usize;
    let cumulative_frac = opts.structured_rescue_cumulative_frac.clamp(0.0, 1.0);

    let sample_condition_indices: Vec<usize> = samples
        .iter()
        .map(|sample| {
            condition_names
                .iter()
                .position(|c| c == &sample.condition)
                .expect("sample condition should exist in condition_names")
        })
        .collect();

    let mut condition_sample_counts = vec![0usize; n_conditions];
    let mut condition_mean_tpms = vec![vec![0.0f64; n_targets]; n_conditions];
    for (sample_idx, &ci) in sample_condition_indices.iter().enumerate() {
        condition_sample_counts[ci] += 1;
        for t in 0..n_targets {
            condition_mean_tpms[ci][t] += phase1_tpms[sample_idx][t];
        }
    }
    for ci in 0..n_conditions {
        if condition_sample_counts[ci] > 0 {
            let denom = condition_sample_counts[ci] as f64;
            for mean_tpm in &mut condition_mean_tpms[ci] {
                *mean_tpm /= denom;
            }
        }
    }

    let mut primary_cond_counts = vec![vec![0u32; n_targets]; n_conditions];
    for (sample_idx, sample_pass) in sample_pass_masks.iter().enumerate() {
        let ci = sample_condition_indices[sample_idx];
        for (t, &pass) in sample_pass.iter().enumerate() {
            if pass {
                primary_cond_counts[ci][t] += 1;
            }
        }
    }
    let primary_condition_support_masks: Vec<Vec<bool>> = primary_cond_counts
        .iter()
        .enumerate()
        .map(|(ci, counts)| counts.iter().map(|&c| c >= cond_k[ci]).collect())
        .collect();

    let mut raw_rescue_mask = vec![false; n_targets];
    for cond_mask in raw_condition_support_masks {
        for (t, &pass) in cond_mask.iter().enumerate() {
            raw_rescue_mask[t] |= pass;
        }
    }

    let raw_rescued = raw_rescue_mask
        .iter()
        .zip(strict_global_mask.iter())
        .filter(|&(&raw, &strict)| raw && !strict)
        .count();

    let candidate_mask: Vec<bool> = strict_global_mask
        .iter()
        .zip(raw_rescue_mask.iter())
        .map(|(&strict, &raw)| strict || raw)
        .collect();
    let groups = build_ec_graph_groups(index, &candidate_mask, n_targets);

    let mut group_members: std::collections::HashMap<u32, Vec<usize>> =
        std::collections::HashMap::new();
    for t in 0..n_targets {
        if candidate_mask[t] {
            group_members.entry(groups[t]).or_default().push(t);
        }
    }

    let mut mask = strict_global_mask.to_vec();
    let mut selected_condition_masks = vec![vec![false; n_targets]; n_conditions];
    let mut selected_rescued = 0usize;
    let mut admitted_groups = 0usize;
    let mut audit_records: Vec<RescueAuditRecord> = Vec::new();

    for (&group, members) in &group_members {
        let group_condition_mean_tpms: Vec<f64> = (0..n_conditions)
            .map(|ci| members.iter().map(|&t| condition_mean_tpms[ci][t]).sum())
            .collect();
        let has_strict_member = members.iter().any(|&t| strict_global_mask[t]);
        let strict_group_per_condition_escape =
            has_strict_member && opts.structured_rescue_strict_group_escape;
        let strict_group_balanced_escape =
            has_strict_member && opts.structured_rescue_strict_group_balanced_escape;
        let strict_group_escape = strict_group_per_condition_escape || strict_group_balanced_escape;
        if has_strict_member && !strict_group_escape {
            for &t in members {
                if raw_rescue_mask[t] && !strict_global_mask[t] {
                    let raw_conditions: Vec<usize> = (0..n_conditions)
                        .filter(|&ci| raw_condition_support_masks[ci][t])
                        .collect();
                    audit_records.push(RescueAuditRecord {
                        target: t,
                        group,
                        group_size: members.len(),
                        has_strict_member,
                        raw_rescue: true,
                        suppression_reason: "strict_group".to_string(),
                        raw_conditions,
                        condition_mean_tpms: condition_mean_tpms
                            .iter()
                            .map(|tpms| tpms[t])
                            .collect(),
                        group_condition_mean_tpms: group_condition_mean_tpms.clone(),
                        ..Default::default()
                    });
                }
            }
            continue;
        }

        let has_raw_rescue = members.iter().any(|&t| raw_rescue_mask[t]);
        if !has_raw_rescue {
            continue;
        }

        let mut admitted_conditions = Vec::new();
        for (ci, &condition_k) in cond_k.iter().enumerate().take(n_conditions) {
            let mut n_group_pass = 0u32;
            for (sample_idx, &sample_ci) in sample_condition_indices.iter().enumerate() {
                if sample_ci != ci {
                    continue;
                }
                let group_tpm: f64 = members.iter().map(|&t| phase1_tpms[sample_idx][t]).sum();
                if group_tpm >= opts.structured_rescue_group_tpm_floor {
                    n_group_pass += 1;
                }
            }
            if n_group_pass >= condition_k {
                admitted_conditions.push(ci);
            }
        }

        if admitted_conditions.is_empty() {
            for &t in members {
                if raw_rescue_mask[t] && !strict_global_mask[t] {
                    let raw_conditions: Vec<usize> = (0..n_conditions)
                        .filter(|&ci| raw_condition_support_masks[ci][t])
                        .collect();
                    audit_records.push(RescueAuditRecord {
                        target: t,
                        group,
                        group_size: members.len(),
                        has_strict_member,
                        raw_rescue: true,
                        suppression_reason: "group_not_admitted".to_string(),
                        raw_conditions,
                        condition_mean_tpms: condition_mean_tpms
                            .iter()
                            .map(|tpms| tpms[t])
                            .collect(),
                        group_condition_mean_tpms: group_condition_mean_tpms.clone(),
                        ..Default::default()
                    });
                }
            }
            continue;
        }

        let mut ranked = Vec::new();
        for &t in members {
            if strict_global_mask[t] {
                continue;
            }

            let raw_in_admitted_condition = admitted_conditions
                .iter()
                .any(|&ci| raw_condition_support_masks[ci][t]);
            if !raw_in_admitted_condition {
                continue;
            }

            let mut score = 0.0f64;
            let mut score_sum = 0.0f64;
            let mut support_conditions = 0usize;
            let mut primary_supported = false;
            let raw_conditions: Vec<usize> = (0..n_conditions)
                .filter(|&ci| raw_condition_support_masks[ci][t])
                .collect();
            for &ci in &admitted_conditions {
                let condition_score = condition_mean_tpms[ci][t];
                score = score.max(condition_score);
                score_sum += condition_score;
                if raw_condition_support_masks[ci][t] {
                    support_conditions += 1;
                }
                primary_supported |= primary_condition_support_masks[ci][t];
            }

            if score > 0.0 {
                ranked.push(RescueCandidate {
                    target: t,
                    score,
                    primary_supported,
                    support_conditions,
                    mean_score: score_sum / admitted_conditions.len() as f64,
                    raw_conditions,
                    condition_mean_tpms: condition_mean_tpms.iter().map(|tpms| tpms[t]).collect(),
                    group_condition_mean_tpms: group_condition_mean_tpms.clone(),
                    peak_rank: None,
                    breadth_rank: None,
                });
            }
        }

        if ranked.is_empty() {
            for &t in members {
                if raw_rescue_mask[t] && !strict_global_mask[t] {
                    let raw_conditions: Vec<usize> = (0..n_conditions)
                        .filter(|&ci| raw_condition_support_masks[ci][t])
                        .collect();
                    audit_records.push(RescueAuditRecord {
                        target: t,
                        group,
                        group_size: members.len(),
                        has_strict_member,
                        admitted: true,
                        raw_rescue: true,
                        suppression_reason: "not_raw_in_admitted_condition".to_string(),
                        raw_conditions,
                        admitted_conditions: admitted_conditions.clone(),
                        condition_mean_tpms: condition_mean_tpms
                            .iter()
                            .map(|tpms| tpms[t])
                            .collect(),
                        group_condition_mean_tpms: group_condition_mean_tpms.clone(),
                        ..Default::default()
                    });
                }
            }
            continue;
        }

        let mut peak_order: Vec<usize> = (0..ranked.len()).collect();
        peak_order.sort_by(|&a, &b| ranked[b].score.partial_cmp(&ranked[a].score).unwrap());
        for (rank, &idx) in peak_order.iter().enumerate() {
            ranked[idx].peak_rank = Some(rank + 1);
        }
        let mut breadth_order: Vec<usize> = (0..ranked.len()).collect();
        breadth_order.sort_by(|&a, &b| {
            ranked[b]
                .support_conditions
                .cmp(&ranked[a].support_conditions)
                .then_with(|| {
                    ranked[b]
                        .mean_score
                        .partial_cmp(&ranked[a].mean_score)
                        .unwrap()
                })
                .then_with(|| ranked[b].score.partial_cmp(&ranked[a].score).unwrap())
        });
        for (rank, &idx) in breadth_order.iter().enumerate() {
            ranked[idx].breadth_rank = Some(rank + 1);
        }

        match opts.structured_rescue_rank_mode {
            StructuredRescueRankMode::Peak => {
                ranked.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());
            }
            StructuredRescueRankMode::Breadth => {
                ranked.sort_by(|a, b| {
                    b.support_conditions
                        .cmp(&a.support_conditions)
                        .then_with(|| b.mean_score.partial_cmp(&a.mean_score).unwrap())
                        .then_with(|| b.score.partial_cmp(&a.score).unwrap())
                });
            }
        }
        let total_score: f64 = ranked.iter().map(|candidate| candidate.score).sum();
        if total_score <= 0.0 {
            for &t in members {
                if raw_rescue_mask[t] && !strict_global_mask[t] {
                    let raw_conditions: Vec<usize> = (0..n_conditions)
                        .filter(|&ci| raw_condition_support_masks[ci][t])
                        .collect();
                    audit_records.push(RescueAuditRecord {
                        target: t,
                        group,
                        group_size: members.len(),
                        has_strict_member,
                        admitted: true,
                        raw_rescue: true,
                        suppression_reason: "zero_score".to_string(),
                        raw_conditions,
                        admitted_conditions: admitted_conditions.clone(),
                        condition_mean_tpms: condition_mean_tpms
                            .iter()
                            .map(|tpms| tpms[t])
                            .collect(),
                        group_condition_mean_tpms: group_condition_mean_tpms.clone(),
                        ..Default::default()
                    });
                }
            }
            continue;
        }

        admitted_groups += 1;
        let target_score = cumulative_frac * total_score;
        let mut cumulative = 0.0f64;
        let mut kept = 0usize;
        let mut selected_targets = std::collections::HashSet::new();
        let mut per_condition_targets = Vec::new();
        if opts.structured_rescue_per_condition_representatives {
            for &ci in &admitted_conditions {
                if let Some(candidate) = ranked
                    .iter()
                    .filter(|candidate| raw_condition_support_masks[ci][candidate.target])
                    .max_by(|a, b| {
                        condition_mean_tpms[ci][a.target]
                            .partial_cmp(&condition_mean_tpms[ci][b.target])
                            .unwrap()
                    })
                {
                    per_condition_targets.push(candidate.target);
                }
            }
        }
        if strict_group_per_condition_escape {
            for &ci in &admitted_conditions {
                if let Some(candidate) = ranked
                    .iter()
                    .filter(|candidate| {
                        raw_condition_support_masks[ci][candidate.target]
                            && condition_mean_tpms[ci][candidate.target]
                                >= opts.structured_rescue_strict_group_candidate_tpm_floor
                    })
                    .max_by(|a, b| {
                        condition_mean_tpms[ci][a.target]
                            .partial_cmp(&condition_mean_tpms[ci][b.target])
                            .unwrap()
                    })
                {
                    per_condition_targets.push(candidate.target);
                }
            }
        }
        if strict_group_balanced_escape {
            let floor = opts.structured_rescue_strict_group_balanced_tpm_floor;
            let min_conditions =
                opts.structured_rescue_strict_group_balanced_min_conditions as usize;
            if let Some(candidate) = ranked
                .iter()
                .filter(|candidate| {
                    condition_mean_tpms
                        .iter()
                        .filter(|tpms| tpms[candidate.target] >= floor)
                        .count()
                        >= min_conditions
                })
                .max_by(|a, b| {
                    a.support_conditions
                        .cmp(&b.support_conditions)
                        .then_with(|| a.mean_score.partial_cmp(&b.mean_score).unwrap())
                        .then_with(|| a.score.partial_cmp(&b.score).unwrap())
                })
            {
                per_condition_targets.push(candidate.target);
            }
        }
        per_condition_targets.sort_unstable();
        per_condition_targets.dedup();
        let mut group_audit: Vec<RescueAuditRecord> = ranked
            .iter()
            .enumerate()
            .map(|(rank, candidate)| RescueAuditRecord {
                target: candidate.target,
                group,
                group_size: members.len(),
                has_strict_member,
                admitted: true,
                raw_rescue: true,
                primary_supported: candidate.primary_supported,
                score: candidate.score,
                score_frac: candidate.score / total_score,
                rank: Some(rank + 1),
                peak_rank: candidate.peak_rank,
                breadth_rank: candidate.breadth_rank,
                support_conditions: candidate.support_conditions,
                mean_score: candidate.mean_score,
                suppression_reason: "ranked_out".to_string(),
                raw_conditions: candidate.raw_conditions.clone(),
                admitted_conditions: admitted_conditions.clone(),
                condition_mean_tpms: candidate.condition_mean_tpms.clone(),
                group_condition_mean_tpms: candidate.group_condition_mean_tpms.clone(),
                ..Default::default()
            })
            .collect();

        for target in per_condition_targets {
            let candidate = ranked
                .iter()
                .find(|candidate| candidate.target == target)
                .expect("per-condition target should be present in ranked candidates");
            selected_targets.insert(candidate.target);
            mask[candidate.target] = true;
            if let Some(rec) = group_audit
                .iter_mut()
                .find(|rec| rec.target == candidate.target)
            {
                rec.selected = true;
                rec.suppression_reason = if strict_group_balanced_escape {
                    "selected_strict_group_balanced_escape".to_string()
                } else if strict_group_per_condition_escape {
                    "selected_strict_group_escape".to_string()
                } else {
                    "selected_condition_representative".to_string()
                };
                rec.selected_conditions = admitted_conditions
                    .iter()
                    .copied()
                    .filter(|&ci| raw_condition_support_masks[ci][candidate.target])
                    .collect();
            }
            selected_rescued += 1;
            kept += 1;
            cumulative += candidate.score;
            for &ci in &admitted_conditions {
                if raw_condition_support_masks[ci][candidate.target] {
                    selected_condition_masks[ci][candidate.target] = true;
                }
            }
        }

        if strict_group_escape {
            audit_records.append(&mut group_audit);
            continue;
        }

        for (rank, candidate) in ranked.iter().enumerate() {
            if selected_targets.contains(&candidate.target) {
                continue;
            }
            let score_frac = candidate.score / total_score;
            let keep = (rank == 0 && kept == 0)
                || (kept < max_isoforms
                    && cumulative < target_score
                    && (candidate.primary_supported || score_frac >= 0.05));
            if !keep {
                continue;
            }

            selected_targets.insert(candidate.target);
            mask[candidate.target] = true;
            if let Some(rec) = group_audit
                .iter_mut()
                .find(|rec| rec.target == candidate.target)
            {
                rec.selected = true;
                rec.suppression_reason = "selected".to_string();
                rec.selected_conditions = admitted_conditions
                    .iter()
                    .copied()
                    .filter(|&ci| raw_condition_support_masks[ci][candidate.target])
                    .collect();
            }
            selected_rescued += 1;
            kept += 1;
            cumulative += candidate.score;
            for &ci in &admitted_conditions {
                if raw_condition_support_masks[ci][candidate.target] {
                    selected_condition_masks[ci][candidate.target] = true;
                }
            }

            if kept >= max_isoforms || (kept >= 1 && cumulative >= target_score) {
                break;
            }
        }

        if opts.structured_rescue_phase1_condition_complements {
            let complement_floor = opts.structured_rescue_phase1_complement_tpm_floor;
            let mut complement_targets = Vec::new();
            for &ci in &admitted_conditions {
                let selected_condition_tpm: f64 = selected_targets
                    .iter()
                    .map(|&target| condition_mean_tpms[ci][target])
                    .sum();
                if selected_condition_tpm > 0.0 || group_condition_mean_tpms[ci] < complement_floor
                {
                    continue;
                }
                if let Some(candidate) = ranked
                    .iter()
                    .filter(|candidate| {
                        !selected_targets.contains(&candidate.target)
                            && candidate.primary_supported
                            && condition_mean_tpms[ci][candidate.target] >= complement_floor
                    })
                    .max_by(|a, b| {
                        condition_mean_tpms[ci][a.target]
                            .partial_cmp(&condition_mean_tpms[ci][b.target])
                            .unwrap()
                    })
                {
                    complement_targets.push((candidate.target, ci));
                }
            }
            complement_targets.sort_unstable();
            complement_targets.dedup();

            for (target, ci) in complement_targets {
                let Some(candidate) = ranked.iter().find(|candidate| candidate.target == target)
                else {
                    continue;
                };
                selected_targets.insert(candidate.target);
                mask[candidate.target] = true;
                selected_condition_masks[ci][candidate.target] = true;
                if let Some(rec) = group_audit
                    .iter_mut()
                    .find(|rec| rec.target == candidate.target)
                {
                    rec.selected = true;
                    rec.suppression_reason = "selected_phase1_condition_complement".to_string();
                    if !rec.selected_conditions.contains(&ci) {
                        rec.selected_conditions.push(ci);
                        rec.selected_conditions.sort_unstable();
                    }
                }
                selected_rescued += 1;
            }
        }

        audit_records.append(&mut group_audit);
    }

    if let Ok(path) = std::env::var("PISCEM_STRUCTURED_RESCUE_AUDIT")
        && let Ok(mut f) = std::fs::File::create(&path)
    {
        use std::io::Write;
        let condition_mean_headers = condition_names
            .iter()
            .map(|condition| format!("mean_tpm_{}", condition))
            .collect::<Vec<_>>()
            .join("\t");
        let group_mean_headers = condition_names
            .iter()
            .map(|condition| format!("group_mean_tpm_{}", condition))
            .collect::<Vec<_>>()
            .join("\t");
        writeln!(
            f,
            "target_idx\ttarget_name\tgroup\tgroup_size\thas_strict_member\tadmitted\traw_rescue\tselected\tprimary_supported\tscore\tscore_frac\trank\tpeak_rank\tbreadth_rank\tsupport_condition_count\tmean_score\tsuppression_reason\traw_conditions\tadmitted_conditions\tselected_conditions\t{}\t{}",
            condition_mean_headers, group_mean_headers
        )
        .ok();
        for rec in &audit_records {
            let format_conditions = |conditions: &[usize]| {
                conditions
                    .iter()
                    .map(|&ci| condition_names[ci].as_str())
                    .collect::<Vec<_>>()
                    .join(",")
            };
            let format_tpms = |tpms: &[f64]| {
                tpms.iter()
                    .map(|tpm| format!("{:.6}", tpm))
                    .collect::<Vec<_>>()
                    .join("\t")
            };
            let raw_conditions = format_conditions(&rec.raw_conditions);
            let admitted_conditions = rec
                .admitted_conditions
                .iter()
                .map(|&ci| condition_names[ci].as_str())
                .collect::<Vec<_>>()
                .join(",");
            let selected_conditions = format_conditions(&rec.selected_conditions);
            let rank = rec.rank.map(|r| r.to_string()).unwrap_or_default();
            let peak_rank = rec.peak_rank.map(|r| r.to_string()).unwrap_or_default();
            let breadth_rank = rec.breadth_rank.map(|r| r.to_string()).unwrap_or_default();
            let condition_mean_tpms = format_tpms(&rec.condition_mean_tpms);
            let group_condition_mean_tpms = format_tpms(&rec.group_condition_mean_tpms);
            writeln!(
                f,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}",
                rec.target,
                ref_names[rec.target],
                rec.group,
                rec.group_size,
                rec.has_strict_member,
                rec.admitted,
                rec.raw_rescue,
                rec.selected,
                rec.primary_supported,
                rec.score,
                rec.score_frac,
                rank,
                peak_rank,
                breadth_rank,
                rec.support_conditions,
                rec.mean_score,
                rec.suppression_reason,
                raw_conditions,
                admitted_conditions,
                selected_conditions,
                condition_mean_tpms,
                group_condition_mean_tpms
            )
            .ok();
        }
        info!("Wrote structured rescue audit to {}", path);
    }

    StructuredRescueResult {
        mask,
        condition_support_masks: selected_condition_masks,
        raw_rescued,
        selected_rescued,
        admitted_groups,
    }
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
#[allow(dead_code)]
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
                let u = targets[j];
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

fn project_counts_to_total(counts: &mut [f64], eff_lens: &[f64], total_weight: f64) {
    let mut sum = 0.0f64;
    for (count, &eff_len) in counts.iter_mut().zip(eff_lens.iter()) {
        if eff_len <= 0.0 || !count.is_finite() || *count < 0.0 {
            *count = 0.0;
        }
        sum += *count;
    }
    if sum > 0.0 && total_weight > 0.0 {
        let scale = total_weight / sum;
        for count in counts {
            *count *= scale;
        }
    }
}

fn project_counts_to_total_with_caps(
    counts: &mut [f64],
    eff_lens: &[f64],
    total_weight: f64,
    caps: Option<&[f64]>,
) {
    let Some(caps) = caps else {
        project_counts_to_total(counts, eff_lens, total_weight);
        return;
    };
    let mut fixed = vec![false; counts.len()];
    for (((count, &eff_len), &cap), is_fixed) in counts
        .iter_mut()
        .zip(eff_lens.iter())
        .zip(caps.iter())
        .zip(fixed.iter_mut())
    {
        if eff_len <= 0.0 || !count.is_finite() || *count < 0.0 {
            *count = 0.0;
            *is_fixed = true;
        } else if cap.is_finite() {
            let cap = cap.max(0.0);
            if *count > cap {
                *count = cap;
                *is_fixed = true;
            }
        }
    }
    loop {
        let sum: f64 = counts.iter().sum();
        let diff = total_weight - sum;
        if diff.abs() <= 1e-6 || total_weight <= 0.0 {
            break;
        }
        let free_sum: f64 = counts
            .iter()
            .zip(fixed.iter())
            .filter_map(|(&count, &is_fixed)| (!is_fixed && count > 0.0).then_some(count))
            .sum();
        if free_sum <= 0.0 {
            break;
        }
        let scale = (total_weight - (sum - free_sum)) / free_sum;
        if scale <= 0.0 || !scale.is_finite() {
            break;
        }
        let mut newly_fixed = false;
        for (((count, is_fixed), &eff_len), &cap) in counts
            .iter_mut()
            .zip(fixed.iter_mut())
            .zip(eff_lens.iter())
            .zip(caps.iter())
        {
            if *is_fixed || eff_len <= 0.0 {
                continue;
            }
            *count *= scale;
            if cap.is_finite() {
                let cap = cap.max(0.0);
                if *count > cap {
                    *count = cap;
                    *is_fixed = true;
                    newly_fixed = true;
                }
            }
        }
        if !newly_fixed {
            break;
        }
    }
}

fn compute_residual_max_rel_diff(prev: &[f64], curr: &[f64], presence_thresh: f64) -> f64 {
    let mut max_rel = 0.0f64;
    for (&p, &c) in prev.iter().zip(curr.iter()) {
        if p > presence_thresh || c > presence_thresh {
            let denom = p.abs().max(c.abs()).max(1.0);
            max_rel = max_rel.max((c - p).abs() / denom);
        }
    }
    max_rel
}

fn compute_residual_mean_rel_diff(prev: &[f64], curr: &[f64], presence_thresh: f64) -> f64 {
    let mut sum_abs_rel = 0.0f64;
    let mut n = 0u64;
    for (&p, &c) in prev.iter().zip(curr.iter()) {
        if p > presence_thresh {
            sum_abs_rel += ((c - p) / p.max(1.0)).abs();
            n += 1;
        }
    }
    if n > 0 { sum_abs_rel / n as f64 } else { 0.0 }
}

fn normal_cdf(z: f64) -> f64 {
    if !z.is_finite() {
        return if z.is_sign_positive() { 1.0 } else { 0.0 };
    }
    let x = z.abs() / 2.0f64.sqrt();
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;
    let erf_approx = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();
    let erf = if z >= 0.0 { erf_approx } else { -erf_approx };
    0.5 * (1.0 + erf)
}

fn rescued_target_audit_fields(
    target_labels: &[u32],
    condition_locked_mask: &[bool],
    ref_names: &[String],
) -> (String, String) {
    let mut indices = Vec::new();
    let mut names = Vec::new();
    for &tid in target_labels {
        let t = tid as usize;
        if condition_locked_mask[t] {
            indices.push(t.to_string());
            names.push(ref_names[t].as_str());
        }
    }
    (indices.join(","), names.join(","))
}

fn residual_squarem_alpha(x0: &[f64], x1: &[f64], x2: &[f64]) -> Option<f64> {
    let mut rr = 0.0f64;
    let mut vv = 0.0f64;
    for ((&a, &b), &c) in x0.iter().zip(x1.iter()).zip(x2.iter()) {
        let r = b - a;
        let v = c - (2.0 * b) + a;
        rr += r * r;
        vv += v * v;
    }
    if rr <= 0.0 || vv <= 0.0 {
        return None;
    }
    let alpha = -(rr / vv).sqrt();
    if alpha.is_finite() {
        Some(alpha.clamp(-10.0, -1.0))
    } else {
        None
    }
}

fn residual_em_step_f64_counts<L: EqLabel>(
    packed: &PackedEqMap<L>,
    eq_counts: &[f64],
    inv_eff_lens: &[f64],
    prev: &[f64],
    curr: &mut [f64],
) {
    curr.fill(0.0);
    let mut weights = Vec::with_capacity(64);
    for (label, &eq_count) in packed.iter_labels().zip(eq_counts.iter()) {
        if eq_count <= 0.0 {
            continue;
        }
        weights.clear();
        let mut denom = 0.0f64;
        for (tid, cond_prob) in label.target_labels().iter().zip(label.target_probs()) {
            let t = *tid as usize;
            let w = cond_prob * prev[t] * inv_eff_lens[t];
            weights.push(w);
            denom += w;
        }
        if denom <= 0.0 {
            continue;
        }
        let scale = eq_count / denom;
        for (tid, &w) in label.target_labels().iter().zip(weights.iter()) {
            curr[*tid as usize] += scale * w;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn subtract_capped_responsibilities<L: EqLabel>(
    packed: &PackedEqMap<L>,
    eq_counts: &mut [f64],
    inv_eff_lens: &[f64],
    estimate: &[f64],
    newly_capped: &[usize],
    caps: &[f64],
) -> Vec<f64> {
    let mut fixed_counts = vec![0.0f64; estimate.len()];
    if newly_capped.is_empty() {
        return fixed_counts;
    }
    let mut cap_scale = vec![0.0f64; estimate.len()];
    for &t in newly_capped {
        if estimate[t] > 0.0 && caps[t].is_finite() {
            cap_scale[t] = (caps[t].max(0.0) / estimate[t]).clamp(0.0, 1.0);
        }
    }
    let mut weights = Vec::with_capacity(64);
    for (label_idx, label) in packed.iter_labels().enumerate() {
        let eq_count = eq_counts[label_idx];
        if eq_count <= 0.0 {
            continue;
        }
        weights.clear();
        let mut denom = 0.0f64;
        for (tid, cond_prob) in label.target_labels().iter().zip(label.target_probs()) {
            let t = *tid as usize;
            let w = cond_prob * estimate[t] * inv_eff_lens[t];
            weights.push(w);
            denom += w;
        }
        if denom <= 0.0 {
            continue;
        }
        let mut fixed_sum = 0.0f64;
        for (tid, &w) in label.target_labels().iter().zip(weights.iter()) {
            let t = *tid as usize;
            if w > 0.0 && cap_scale[t] > 0.0 {
                let fixed = eq_count * w / denom * cap_scale[t];
                fixed_counts[t] += fixed;
                fixed_sum += fixed;
            }
        }
        eq_counts[label_idx] = (eq_count - fixed_sum).max(0.0);
    }
    fixed_counts
}

#[allow(clippy::too_many_arguments)]
fn phase2_em_step_f64_counts_active_set_caps<L: EqLabel>(
    packed: &PackedEqMap<L>,
    eq_counts: &[f64],
    eff_lens: &[f64],
    active_mask: &[bool],
    init_counts: Option<&[f64]>,
    opts: &ConsensusQuantOpts,
    residual_caps: &[f64],
    eq_index: Option<&crate::utils::txp_selection::TranscriptEqIndex>,
) -> Vec<f64> {
    let n_targets = eff_lens.len();
    let mut residual_eq_counts = eq_counts.to_vec();
    let mut active = active_mask.to_vec();
    let mut fixed_counts = vec![0.0f64; n_targets];
    let mut init = init_counts.map(|counts| counts.to_vec());
    let mut inline_opts = opts.clone();
    inline_opts.phase2_max_iter = Some(phase2_max_iter(opts).clamp(1, 24));
    let inv_eff_lens = eff_lens
        .iter()
        .map(|&eff_len| {
            let inv = 1.0 / eff_len;
            if inv.is_finite() { inv } else { 0.0 }
        })
        .collect::<Vec<_>>();
    let mut estimate;
    let mut rounds = 0usize;
    loop {
        rounds += 1;
        estimate = phase2_em_step_f64_counts(
            packed,
            &residual_eq_counts,
            eff_lens,
            &active,
            init.as_deref(),
            &inline_opts,
            None,
            None,
            false,
            eq_index,
        );
        let newly_capped = estimate
            .iter()
            .zip(active.iter())
            .zip(residual_caps.iter())
            .enumerate()
            .filter_map(|(t, ((&count, &is_active), &cap))| {
                (is_active && cap.is_finite() && count > cap.max(0.0)).then_some(t)
            })
            .collect::<Vec<_>>();
        if newly_capped.is_empty() {
            estimate = phase2_em_step_f64_counts(
                packed,
                &residual_eq_counts,
                eff_lens,
                &active,
                Some(&phase2_init_counts(&estimate, &active)),
                opts,
                None,
                None,
                false,
                eq_index,
            );
            break;
        }
        if rounds >= 32 {
            break;
        }
        let fixed_delta = subtract_capped_responsibilities(
            packed,
            &mut residual_eq_counts,
            &inv_eff_lens,
            &estimate,
            &newly_capped,
            residual_caps,
        );
        for &t in &newly_capped {
            active[t] = false;
            fixed_counts[t] += fixed_delta[t];
        }
        init = Some(phase2_init_counts(&estimate, &active));
    }
    for (count, fixed) in estimate.iter_mut().zip(fixed_counts.iter()) {
        *count += *fixed;
    }
    info!(
        "Residual active-set cap stats: rounds={} fixed_targets={} residual_weight={:.3}",
        rounds,
        fixed_counts.iter().filter(|&&x| x > 0.0).count(),
        residual_eq_counts.iter().sum::<f64>()
    );
    estimate
}

fn apply_one_sided_floor_barrier(
    counts: &mut [f64],
    floor_counts: &[f64],
    eff_lens: &[f64],
    total_weight: f64,
    weight: f64,
) {
    if weight <= 0.0 {
        return;
    }
    let shrink = weight / (1.0 + weight);
    let mut touched = false;
    for ((count, &floor), &eff_len) in counts
        .iter_mut()
        .zip(floor_counts.iter())
        .zip(eff_lens.iter())
    {
        if eff_len > 0.0 && floor > 0.0 && *count < floor {
            *count += shrink * (floor - *count);
            touched = true;
        }
    }
    if touched {
        project_counts_to_total(counts, eff_lens, total_weight);
    }
}

fn residual_tail_mask(
    prev: &[f64],
    curr: &[f64],
    presence_thresh: f64,
    threshold: f64,
) -> Vec<bool> {
    prev.iter()
        .zip(curr.iter())
        .map(|(&p, &c)| {
            if p > presence_thresh || c > presence_thresh {
                let denom = p.abs().max(c.abs()).max(1.0);
                ((c - p).abs() / denom) > threshold
            } else {
                false
            }
        })
        .collect()
}

struct TailExpansion {
    active: Vec<bool>,
    eq_indices: Vec<usize>,
    crossed_eqs: usize,
    crossed_weight: f64,
    boundary_eqs: usize,
    boundary_weight: f64,
}

fn expand_tail_component<L: EqLabel>(
    packed: &PackedEqMap<L>,
    eq_counts: &[f64],
    eff_lens: &[f64],
    initial_tail: &[bool],
    eq_index: Option<&crate::utils::txp_selection::TranscriptEqIndex>,
) -> TailExpansion {
    let n_targets = eff_lens.len();
    let mut active = vec![false; n_targets];
    let mut active_count = 0usize;
    for (t, (&tail, &eff_len)) in initial_tail.iter().zip(eff_lens.iter()).enumerate() {
        if tail && eff_len > 0.0 {
            active[t] = true;
            active_count += 1;
        }
    }

    if active_count == 0 {
        return TailExpansion {
            active,
            eq_indices: Vec::new(),
            crossed_eqs: 0,
            crossed_weight: 0.0,
            boundary_eqs: 0,
            boundary_weight: 0.0,
        };
    }

    if let Some(index) = eq_index {
        return expand_tail_component_indexed(
            packed,
            eq_counts,
            eff_lens,
            active,
            active_count,
            index,
        );
    }

    let mut frontier = active.clone();
    let mut hop = 0usize;
    let mut changed = true;
    let mut crossed_eqs_total = 0usize;
    let mut crossed_weight_total = 0.0f64;
    while changed {
        changed = false;
        let before_count = active_count;
        let mut bridge_eqs = 0usize;
        let mut bridge_weight = 0.0f64;
        let mut next_frontier = vec![false; n_targets];
        for (label, &eq_count) in packed.iter_labels().zip(eq_counts.iter()) {
            if eq_count <= 0.0 {
                continue;
            }
            let touches_frontier = label
                .target_labels()
                .iter()
                .any(|tid| frontier[*tid as usize]);
            if !touches_frontier {
                continue;
            }
            let mut added_from_eq = false;
            for tid in label.target_labels() {
                let t = *tid as usize;
                if eff_lens[t] > 0.0 && !active[t] {
                    active[t] = true;
                    next_frontier[t] = true;
                    active_count += 1;
                    changed = true;
                    added_from_eq = true;
                }
            }
            if added_from_eq {
                bridge_eqs += 1;
                bridge_weight += eq_count;
            }
        }
        let added = active_count - before_count;
        crossed_eqs_total += bridge_eqs;
        crossed_weight_total += bridge_weight;
        if added > 0 {
            info!(
                "Residual tail component hop {}: added_targets={} bridge_eqs={} bridge_weight={:.3}",
                hop + 1,
                added,
                bridge_eqs,
                bridge_weight
            );
        }
        frontier = next_frontier;
        hop += 1;
    }

    let mut eq_indices = Vec::new();
    let boundary_eqs = 0usize;
    let boundary_weight = 0.0f64;
    for (eq_idx, (label, &eq_count)) in packed.iter_labels().zip(eq_counts.iter()).enumerate() {
        if eq_count <= 0.0 {
            continue;
        }
        let touches_active = label
            .target_labels()
            .iter()
            .any(|tid| active[*tid as usize]);
        if touches_active {
            eq_indices.push(eq_idx);
        }
    }
    TailExpansion {
        active,
        eq_indices,
        crossed_eqs: crossed_eqs_total,
        crossed_weight: crossed_weight_total,
        boundary_eqs,
        boundary_weight,
    }
}

fn expand_tail_component_indexed<L: EqLabel>(
    packed: &PackedEqMap<L>,
    eq_counts: &[f64],
    eff_lens: &[f64],
    mut active: Vec<bool>,
    mut active_count: usize,
    index: &crate::utils::txp_selection::TranscriptEqIndex,
) -> TailExpansion {
    let n_targets = eff_lens.len();
    let mut frontier = active.clone();
    let mut visited_eq = vec![false; packed.len()];
    let mut eq_indices = Vec::new();
    let mut hop = 0usize;
    let mut crossed_eqs_total = 0usize;
    let mut crossed_weight_total = 0.0f64;

    loop {
        let before_count = active_count;
        let mut bridge_eqs = 0usize;
        let mut bridge_weight = 0.0f64;
        let mut next_frontier = vec![false; n_targets];

        for (t, &is_frontier) in frontier.iter().enumerate().take(n_targets) {
            if !is_frontier {
                continue;
            }
            for &eq_idx_u32 in index.signature(t) {
                let eq_idx = eq_idx_u32 as usize;
                if eq_idx >= packed.len() || visited_eq[eq_idx] || eq_counts[eq_idx] <= 0.0 {
                    continue;
                }
                visited_eq[eq_idx] = true;
                eq_indices.push(eq_idx);
                let label = packed.refs_for_eqc(eq_idx);
                let mut added_from_eq = false;
                for tid in label.target_labels() {
                    let target = *tid as usize;
                    if eff_lens[target] > 0.0 && !active[target] {
                        active[target] = true;
                        next_frontier[target] = true;
                        active_count += 1;
                        added_from_eq = true;
                    }
                }
                if added_from_eq {
                    bridge_eqs += 1;
                    bridge_weight += eq_counts[eq_idx];
                }
            }
        }

        let added = active_count - before_count;
        crossed_eqs_total += bridge_eqs;
        crossed_weight_total += bridge_weight;
        if added > 0 {
            info!(
                "Residual tail component hop {}: added_targets={} bridge_eqs={} bridge_weight={:.3}",
                hop + 1,
                added,
                bridge_eqs,
                bridge_weight
            );
        }
        if added == 0 {
            break;
        }
        frontier = next_frontier;
        hop += 1;
    }

    TailExpansion {
        active,
        eq_indices,
        crossed_eqs: crossed_eqs_total,
        crossed_weight: crossed_weight_total,
        boundary_eqs: 0,
        boundary_weight: 0.0,
    }
}

#[allow(clippy::too_many_arguments)]
fn residual_tail_em_step_f64_counts<L: EqLabel>(
    packed: &PackedEqMap<L>,
    eq_counts: &[f64],
    eq_indices: &[usize],
    inv_eff_lens: &[f64],
    active_tail: &[bool],
    fixed_counts: &[f64],
    prev: &[f64],
    curr: &mut [f64],
) {
    for (c, &active) in curr.iter_mut().zip(active_tail.iter()) {
        if active {
            *c = 0.0;
        }
    }
    let mut weights = Vec::with_capacity(64);
    for &eq_idx in eq_indices {
        let eq_count = eq_counts[eq_idx];
        if eq_count <= 0.0 {
            continue;
        }
        let label = packed.refs_for_eqc(eq_idx);
        weights.clear();
        let mut denom = 0.0f64;
        for (tid, cond_prob) in label.target_labels().iter().zip(label.target_probs()) {
            let t = *tid as usize;
            let count = if active_tail[t] {
                prev[t]
            } else {
                fixed_counts[t]
            };
            let w = cond_prob * count * inv_eff_lens[t];
            weights.push(w);
            denom += w;
        }
        if denom <= 0.0 {
            continue;
        }
        let scale = eq_count / denom;
        for (tid, &w) in label.target_labels().iter().zip(weights.iter()) {
            let t = *tid as usize;
            if active_tail[t] {
                curr[t] += scale * w;
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn refine_residual_tail_f64_counts<L: EqLabel>(
    packed: &PackedEqMap<L>,
    eq_counts: &[f64],
    eff_lens: &[f64],
    inv_eff_lens: &[f64],
    counts: &mut [f64],
    tail_seed: &[bool],
    eq_index: Option<&crate::utils::txp_selection::TranscriptEqIndex>,
    opts: &ConsensusQuantOpts,
    max_iter: u32,
) {
    let seed_count = tail_seed.iter().filter(|&&x| x).count();
    if seed_count == 0 || max_iter == 0 {
        return;
    }
    let expansion = expand_tail_component(packed, eq_counts, eff_lens, tail_seed, eq_index);
    let active_tail = expansion.active;
    let eq_indices = expansion.eq_indices;
    let active_count = active_tail.iter().filter(|&&x| x).count();
    if active_count == 0 || eq_indices.is_empty() {
        return;
    }

    let mut prev = counts.to_vec();
    let fixed_counts = counts.to_vec();
    let mut curr = counts.to_vec();
    let mut final_mean_rel = f64::INFINITY;
    let mut final_max_rel = f64::INFINITY;
    let tail_conv_thresh = phase2_convergence_thresh(opts);
    let mut steps = 0u32;
    while steps < max_iter {
        residual_tail_em_step_f64_counts(
            packed,
            eq_counts,
            &eq_indices,
            inv_eff_lens,
            &active_tail,
            &fixed_counts,
            &prev,
            &mut curr,
        );
        steps += 1;
        final_mean_rel = compute_residual_mean_rel_diff(&prev, &curr, opts.presence_thresh);
        final_max_rel = compute_residual_max_rel_diff(&prev, &curr, opts.presence_thresh);
        prev.clone_from_slice(&curr);
        if final_max_rel < tail_conv_thresh {
            break;
        }
    }

    counts.clone_from_slice(&prev);
    info!(
        "Residual tail refinement: seed_targets={} active_targets={} active_eqs={} crossed_eqs={} crossed_weight={:.3} boundary_eqs={} boundary_weight={:.3} steps={} final_mean_rel_diff={:.6} final_max_rel_diff={:.6}",
        seed_count,
        active_count,
        eq_indices.len(),
        expansion.crossed_eqs,
        expansion.crossed_weight,
        expansion.boundary_eqs,
        expansion.boundary_weight,
        steps,
        final_mean_rel,
        final_max_rel
    );
}

#[allow(clippy::too_many_arguments)]
fn phase2_em_step_f64_counts<L: EqLabel>(
    packed: &PackedEqMap<L>,
    eq_counts: &[f64],
    eff_lens: &[f64],
    active_mask: &[bool],
    init_counts: Option<&[f64]>,
    opts: &ConsensusQuantOpts,
    floor_barrier: Option<(&[f64], f64)>,
    residual_caps: Option<&[f64]>,
    active_set_caps: bool,
    eq_index: Option<&crate::utils::txp_selection::TranscriptEqIndex>,
) -> Vec<f64> {
    if active_set_caps
        && floor_barrier.is_none()
        && let Some(caps) = residual_caps
    {
        return phase2_em_step_f64_counts_active_set_caps(
            packed,
            eq_counts,
            eff_lens,
            active_mask,
            init_counts,
            opts,
            caps,
            eq_index,
        );
    }
    let mut eff_lens = eff_lens.to_vec();
    for (eff_len, &active) in eff_lens.iter_mut().zip(active_mask.iter()) {
        if !active {
            *eff_len = 0.0;
        }
    }
    let n_targets = eff_lens.len();
    let total_weight: f64 = eq_counts.iter().sum();
    if total_weight <= 0.0 {
        return vec![0.0; n_targets];
    }
    let mut prev = if let Some(init) = init_counts {
        init.iter()
            .zip(eff_lens.iter())
            .map(|(&count, &eff_len)| if eff_len > 0.0 { count.max(0.0) } else { 0.0 })
            .collect::<Vec<_>>()
    } else {
        let n_active = eff_lens.iter().filter(|&&eff_len| eff_len > 0.0).count();
        let avg = if n_active > 0 {
            total_weight / n_active as f64
        } else {
            0.0
        };
        eff_lens
            .iter()
            .map(|&eff_len| if eff_len > 0.0 { avg } else { 0.0 })
            .collect::<Vec<_>>()
    };
    project_counts_to_total_with_caps(&mut prev, &eff_lens, total_weight, residual_caps);

    let inv_eff_lens = eff_lens
        .iter()
        .map(|&eff_len| {
            let inv = 1.0 / eff_len;
            if inv.is_finite() { inv } else { 0.0 }
        })
        .collect::<Vec<_>>();
    let mut x0 = prev;
    let mut x1 = vec![0.0f64; n_targets];
    let mut x2 = vec![0.0f64; n_targets];
    let mut x_sq = vec![0.0f64; n_targets];
    let mut x_next = vec![0.0f64; n_targets];
    let max_iter = phase2_max_iter(opts);
    let conv_thresh = phase2_convergence_thresh(opts);
    let mut em_steps = 0u32;
    let mut final_mean_rel_diff = f64::INFINITY;
    let mut final_max_rel_diff = f64::INFINITY;
    let mut accel_attempts = 0u32;
    let mut accel_accepts = 0u32;
    let mut accel_alpha_none = 0u32;
    let mut accel_invalid_candidate = 0u32;
    let mut accel_not_improved = 0u32;
    let mut alpha_min = f64::INFINITY;
    let mut alpha_max = f64::NEG_INFINITY;
    let mut candidate_ratio_sum = 0.0f64;
    let mut candidate_ratio_n = 0u32;
    while em_steps < max_iter {
        residual_em_step_f64_counts(packed, eq_counts, &inv_eff_lens, &x0, &mut x1);
        if let Some((floor_counts, weight)) = floor_barrier {
            apply_one_sided_floor_barrier(&mut x1, floor_counts, &eff_lens, total_weight, weight);
        }
        project_counts_to_total_with_caps(&mut x1, &eff_lens, total_weight, residual_caps);
        em_steps += 1;
        let rel1 = compute_residual_mean_rel_diff(&x0, &x1, opts.presence_thresh);
        final_mean_rel_diff = rel1;
        final_max_rel_diff = compute_residual_max_rel_diff(&x0, &x1, opts.presence_thresh);
        if rel1 < conv_thresh || em_steps >= max_iter {
            x0.clone_from_slice(&x1);
            break;
        }

        // Keep the first few iterations plain EM, matching the shared SQUAREM
        // solver's conservative burn-in.
        if em_steps < 3 {
            x0.clone_from_slice(&x1);
            continue;
        }

        accel_attempts += 1;
        residual_em_step_f64_counts(packed, eq_counts, &inv_eff_lens, &x1, &mut x2);
        if let Some((floor_counts, weight)) = floor_barrier {
            apply_one_sided_floor_barrier(&mut x2, floor_counts, &eff_lens, total_weight, weight);
        }
        project_counts_to_total_with_caps(&mut x2, &eff_lens, total_weight, residual_caps);
        em_steps += 1;
        let ordinary_rel = compute_residual_mean_rel_diff(&x1, &x2, opts.presence_thresh);

        let use_candidate = if let Some(alpha) = residual_squarem_alpha(&x0, &x1, &x2) {
            alpha_min = alpha_min.min(alpha);
            alpha_max = alpha_max.max(alpha);
            for (((sq, &a), &b), &c) in x_sq.iter_mut().zip(x0.iter()).zip(x1.iter()).zip(x2.iter())
            {
                let r = b - a;
                let v = c - (2.0 * b) + a;
                *sq = a - (2.0 * alpha * r) + (alpha * alpha * v);
            }
            project_counts_to_total_with_caps(&mut x_sq, &eff_lens, total_weight, residual_caps);
            if em_steps < max_iter {
                residual_em_step_f64_counts(packed, eq_counts, &inv_eff_lens, &x_sq, &mut x_next);
                if let Some((floor_counts, weight)) = floor_barrier {
                    apply_one_sided_floor_barrier(
                        &mut x_next,
                        floor_counts,
                        &eff_lens,
                        total_weight,
                        weight,
                    );
                }
                project_counts_to_total_with_caps(
                    &mut x_next,
                    &eff_lens,
                    total_weight,
                    residual_caps,
                );
                em_steps += 1;
                let candidate_rel =
                    compute_residual_mean_rel_diff(&x_sq, &x_next, opts.presence_thresh);
                if ordinary_rel > 0.0 && candidate_rel.is_finite() {
                    candidate_ratio_sum += candidate_rel / ordinary_rel;
                    candidate_ratio_n += 1;
                }
                if !x_next.iter().all(|x| x.is_finite() && *x >= 0.0) {
                    accel_invalid_candidate += 1;
                    false
                } else if candidate_rel < ordinary_rel {
                    true
                } else {
                    accel_not_improved += 1;
                    false
                }
            } else {
                false
            }
        } else {
            accel_alpha_none += 1;
            false
        };

        if use_candidate {
            accel_accepts += 1;
            final_mean_rel_diff =
                compute_residual_mean_rel_diff(&x0, &x_next, opts.presence_thresh);
            final_max_rel_diff = compute_residual_max_rel_diff(&x0, &x_next, opts.presence_thresh);
            x0.clone_from_slice(&x_next);
        } else {
            final_mean_rel_diff = compute_residual_mean_rel_diff(&x0, &x2, opts.presence_thresh);
            final_max_rel_diff = compute_residual_max_rel_diff(&x0, &x2, opts.presence_thresh);
            x0.clone_from_slice(&x2);
        }
        if final_mean_rel_diff < conv_thresh {
            break;
        }
    }

    residual_em_step_f64_counts(packed, eq_counts, &inv_eff_lens, &x0, &mut x1);
    if let Some((floor_counts, weight)) = floor_barrier {
        apply_one_sided_floor_barrier(&mut x1, floor_counts, &eff_lens, total_weight, weight);
    }
    project_counts_to_total_with_caps(&mut x1, &eff_lens, total_weight, residual_caps);
    let tail_threshold = (20.0 * conv_thresh).max(0.01);
    let tail_seed = residual_tail_mask(&x0, &x1, opts.presence_thresh, tail_threshold);
    let tail_seed_count = tail_seed.iter().filter(|&&x| x).count();
    if tail_seed_count > 0 {
        let remaining_iter = max_iter.saturating_sub(em_steps);
        let tail_max_iter = remaining_iter.min(750);
        if tail_max_iter > 0 {
            refine_residual_tail_f64_counts(
                packed,
                eq_counts,
                &eff_lens,
                &inv_eff_lens,
                &mut x0,
                &tail_seed,
                eq_index,
                opts,
                tail_max_iter,
            );
            residual_em_step_f64_counts(packed, eq_counts, &inv_eff_lens, &x0, &mut x1);
            if let Some((floor_counts, weight)) = floor_barrier {
                apply_one_sided_floor_barrier(
                    &mut x1,
                    floor_counts,
                    &eff_lens,
                    total_weight,
                    weight,
                );
            }
            project_counts_to_total_with_caps(&mut x1, &eff_lens, total_weight, residual_caps);
            final_mean_rel_diff = compute_residual_mean_rel_diff(&x0, &x1, opts.presence_thresh);
            final_max_rel_diff = compute_residual_max_rel_diff(&x0, &x1, opts.presence_thresh);
        }
    }

    for count in &mut x0 {
        if *count < opts.presence_thresh {
            *count = 0.0;
        }
    }
    residual_em_step_f64_counts(packed, eq_counts, &inv_eff_lens, &x0, &mut x1);
    if let Some((floor_counts, weight)) = floor_barrier {
        apply_one_sided_floor_barrier(&mut x1, floor_counts, &eff_lens, total_weight, weight);
    }
    project_counts_to_total_with_caps(&mut x1, &eff_lens, total_weight, residual_caps);
    info!(
        "Residual locked-rescue SQUAREM stats: em_steps={} accel_attempts={} accel_accepts={} final_mean_rel_diff={:.6} final_max_rel_diff={:.6} alpha_none={} invalid_candidate={} not_improved={} alpha_range=[{:.3},{:.3}] mean_candidate_rel_ratio={:.3}",
        em_steps,
        accel_attempts,
        accel_accepts,
        final_mean_rel_diff,
        final_max_rel_diff,
        accel_alpha_none,
        accel_invalid_candidate,
        accel_not_improved,
        alpha_min,
        alpha_max,
        if candidate_ratio_n > 0 {
            candidate_ratio_sum / candidate_ratio_n as f64
        } else {
            f64::NAN
        }
    );
    x1
}

fn record_residual_eq_count(
    residual_eq_counts: &mut Vec<f64>,
    collapsed: Option<&CollapsedEqMap>,
    label_idx: usize,
    count: f64,
) {
    if let Some(cm) = collapsed {
        if count > 0.0 {
            residual_eq_counts[cm.pos_to_collapsed[label_idx] as usize] += count;
        }
    } else {
        residual_eq_counts.push(count);
    }
}

fn sample_cv_f64(values: &[f64]) -> f64 {
    if values.len() <= 1 {
        return 0.0;
    }
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    if mean <= 0.0 {
        return 0.0;
    }
    let var = values
        .iter()
        .map(|v| {
            let d = v - mean;
            d * d
        })
        .sum::<f64>()
        / (n - 1.0);
    var.sqrt() / (mean + 1e-6)
}

#[allow(clippy::too_many_arguments)]
fn condition_rescue_instability_masks(
    samples: &[SampleEntry],
    condition_names: &[String],
    condition_rescue_only: &[bool],
    condition_support_masks: &[Vec<bool>],
    rescue_sample_pass_masks: &[Vec<bool>],
    phase1_counts: &[Vec<f64>],
    cv_threshold: f64,
    min_pass_fraction: f64,
) -> Vec<Vec<bool>> {
    let n_conditions = condition_names.len();
    let n_targets = condition_rescue_only.len();
    let mut sample_indices_by_condition = vec![Vec::<usize>::new(); n_conditions];
    for (sample_idx, sample) in samples.iter().enumerate() {
        let ci = condition_names
            .iter()
            .position(|condition| condition == &sample.condition)
            .expect("sample condition should exist in condition_names");
        sample_indices_by_condition[ci].push(sample_idx);
    }

    let cv_threshold = cv_threshold.max(0.0);
    let min_pass_fraction = min_pass_fraction.clamp(0.0, 1.0);
    let mut unstable = vec![vec![false; n_targets]; n_conditions];
    for ci in 0..n_conditions {
        let sample_indices = &sample_indices_by_condition[ci];
        if sample_indices.is_empty() {
            continue;
        }
        let denom = sample_indices.len() as f64;
        for t in 0..n_targets {
            if !condition_rescue_only[t] || !condition_support_masks[ci][t] {
                continue;
            }
            let pass_count = sample_indices
                .iter()
                .filter(|&&sample_idx| rescue_sample_pass_masks[sample_idx][t])
                .count();
            let pass_fraction = pass_count as f64 / denom;
            let counts: Vec<f64> = sample_indices
                .iter()
                .map(|&sample_idx| phase1_counts[sample_idx][t])
                .collect();
            let count_cv = sample_cv_f64(&counts);
            unstable[ci][t] = pass_fraction < min_pass_fraction || count_cv >= cv_threshold;
        }
    }
    unstable
}

fn condition_rescue_guarded_confidence_masks(
    samples: &[SampleEntry],
    condition_names: &[String],
    condition_rescue_only: &[bool],
    condition_support_masks: &[Vec<bool>],
    phase1_counts: &[Vec<f64>],
    mean_count_threshold: f64,
) -> Vec<Vec<bool>> {
    let n_conditions = condition_names.len();
    let n_targets = condition_rescue_only.len();
    let mut sample_indices_by_condition = vec![Vec::<usize>::new(); n_conditions];
    for (sample_idx, sample) in samples.iter().enumerate() {
        let ci = condition_names
            .iter()
            .position(|condition| condition == &sample.condition)
            .expect("sample condition should exist in condition_names");
        sample_indices_by_condition[ci].push(sample_idx);
    }

    let mean_count_threshold = mean_count_threshold.max(0.0);
    let mut condition_mean_count = vec![vec![0.0f64; n_targets]; n_conditions];
    for ci in 0..n_conditions {
        let sample_indices = &sample_indices_by_condition[ci];
        if sample_indices.is_empty() {
            continue;
        }
        let denom = sample_indices.len() as f64;
        for &sample_idx in sample_indices {
            for t in 0..n_targets {
                condition_mean_count[ci][t] += phase1_counts[sample_idx][t] / denom;
            }
        }
    }

    let mut allow_relax = vec![vec![false; n_targets]; n_conditions];
    for t in 0..n_targets {
        if !condition_rescue_only[t] {
            continue;
        }
        let evidence_condition_count = (0..n_conditions)
            .filter(|&ci| condition_mean_count[ci][t] >= mean_count_threshold)
            .count();
        let condition_local = evidence_condition_count <= 1;
        if !condition_local {
            continue;
        }
        for ci in 0..n_conditions {
            allow_relax[ci][t] = condition_support_masks[ci][t];
        }
    }
    allow_relax
}

#[allow(clippy::too_many_arguments)]
fn condition_rescue_transcript_stability_masks(
    samples: &[SampleEntry],
    condition_names: &[String],
    condition_rescue_only: &[bool],
    condition_support_masks: &[Vec<bool>],
    rescue_sample_pass_masks: &[Vec<bool>],
    phase1_counts: &[Vec<f64>],
    mean_count_threshold: f64,
    cv_threshold: f64,
    min_pass_fraction: f64,
) -> Vec<Vec<bool>> {
    let n_conditions = condition_names.len();
    let n_targets = condition_rescue_only.len();
    let mut sample_indices_by_condition = vec![Vec::<usize>::new(); n_conditions];
    for (sample_idx, sample) in samples.iter().enumerate() {
        let ci = condition_names
            .iter()
            .position(|condition| condition == &sample.condition)
            .expect("sample condition should exist in condition_names");
        sample_indices_by_condition[ci].push(sample_idx);
    }

    let mean_count_threshold = mean_count_threshold.max(0.0);
    let cv_threshold = cv_threshold.max(0.0);
    let min_pass_fraction = min_pass_fraction.clamp(0.0, 1.0);
    let mut allow_relax = vec![vec![false; n_targets]; n_conditions];
    for ci in 0..n_conditions {
        let sample_indices = &sample_indices_by_condition[ci];
        if sample_indices.is_empty() {
            continue;
        }
        let denom = sample_indices.len() as f64;
        for t in 0..n_targets {
            if !condition_rescue_only[t] || !condition_support_masks[ci][t] {
                continue;
            }
            let pass_count = sample_indices
                .iter()
                .filter(|&&sample_idx| rescue_sample_pass_masks[sample_idx][t])
                .count();
            let pass_fraction = pass_count as f64 / denom;
            let counts: Vec<f64> = sample_indices
                .iter()
                .map(|&sample_idx| phase1_counts[sample_idx][t])
                .collect();
            let mean_count = counts.iter().sum::<f64>() / denom;
            let count_cv = sample_cv_f64(&counts);
            let stable = pass_fraction >= min_pass_fraction
                && mean_count >= mean_count_threshold
                && count_cv <= cv_threshold;
            allow_relax[ci][t] = !stable;
        }
    }
    allow_relax
}

#[allow(clippy::too_many_arguments)]
fn lock_condition_rescue_allocations_for_sample<L: EqLabel>(
    sample_name: &str,
    ref_names: &[String],
    packed: &PackedEqMap<L>,
    phase1_counts: &[f64],
    eff_lens: &[f64],
    condition_locked_mask: &[bool],
    condition_relax_mask: Option<&[bool]>,
    collapsed: Option<&CollapsedEqMap>,
    lock_fraction: f64,
    lock_mode: &ConditionRescueLockMode,
    full_lock_threshold: f64,
    min_lock_threshold: f64,
    credible_floor_z: f64,
    enrichment_posterior_threshold: f64,
) -> (Vec<f64>, Vec<f64>) {
    let mut residual_eq_counts = if let Some(cm) = collapsed {
        vec![0.0f64; cm.packed.len()]
    } else {
        Vec::with_capacity(packed.len())
    };
    let mut locked_counts = vec![0.0f64; phase1_counts.len()];
    let audit_targets: Option<HashSet<String>> = std::env::var("PISCEM_LOCK_AUDIT_TARGETS")
        .ok()
        .and_then(|path| read_to_string(path).ok())
        .map(|contents| {
            contents
                .lines()
                .map(str::trim)
                .filter(|line| !line.is_empty())
                .map(ToOwned::to_owned)
                .collect()
        });
    let mut audit = std::env::var("PISCEM_LOCK_AUDIT_DIR")
        .ok()
        .and_then(|dir| {
            create_dir_all(&dir).ok()?;
            let path = std::path::Path::new(&dir).join(format!("{sample_name}.lock_ec.tsv"));
            let mut file = File::create(path).ok()?;
            writeln!(
                file,
                "sample\teq_idx\teq_count\tresidual_count\tlocked_sum\trescue_fraction\tselected_lock_fraction\tec_size\trescued_targets\tlocked_positive_targets\trescued_target_indices\trescued_target_names\tlocked_target_indices"
            )
            .ok()?;
            Some(file)
        });
    let lock_fraction = lock_fraction.clamp(0.0, 1.0);
    let min_lock_threshold = min_lock_threshold.clamp(0.0, 1.0);
    let full_lock_threshold = full_lock_threshold.clamp(min_lock_threshold, 1.0);
    let credible_floor_z = credible_floor_z.max(0.0);
    let enrichment_posterior_threshold = enrichment_posterior_threshold.clamp(0.0, 1.0);
    let inv_eff_lens = eff_lens
        .iter()
        .map(|&eff_len| {
            let inv = 1.0 / eff_len;
            if inv.is_finite() { inv } else { 0.0 }
        })
        .collect::<Vec<_>>();
    let mut weights = Vec::with_capacity(64);

    let confidence_lock_fraction = |rescue_fraction: f64| -> f64 {
        if rescue_fraction >= full_lock_threshold {
            1.0
        } else if rescue_fraction >= min_lock_threshold {
            lock_fraction
        } else {
            0.0
        }
    };
    let floor_smooth_confidence_lock_fraction = |rescue_fraction: f64| -> f64 {
        if rescue_fraction >= full_lock_threshold {
            1.0
        } else if rescue_fraction < min_lock_threshold {
            0.0
        } else {
            let span = full_lock_threshold - min_lock_threshold;
            if span <= 0.0 {
                1.0
            } else {
                let x = ((rescue_fraction - min_lock_threshold) / span).clamp(0.0, 1.0);
                let smooth = x * x * x * (x * (x * 6.0 - 15.0) + 10.0);
                lock_fraction + (1.0 - lock_fraction) * smooth
            }
        }
    };
    let credible_floor_lock_fraction = |rescue_fraction: f64, eq_count: f64| -> f64 {
        if rescue_fraction <= 0.0 || eq_count <= 0.0 {
            return 0.0;
        }
        let se = (rescue_fraction * (1.0 - rescue_fraction) / eq_count).sqrt();
        let lower = (rescue_fraction - credible_floor_z * se).max(0.0);
        (lower / rescue_fraction).clamp(0.0, 1.0)
    };
    let enrichment_lock_fraction =
        |rescue_fraction: f64, rescued_targets: usize, ec_size: usize| {
            if rescue_fraction <= 0.0 || rescue_fraction >= 1.0 {
                return rescue_fraction.clamp(0.0, 1.0);
            }
            if rescued_targets == 0 || rescued_targets >= ec_size || ec_size == 0 {
                return rescue_fraction.clamp(0.0, 1.0);
            }
            let null_fraction = (rescued_targets as f64 / ec_size as f64).clamp(1e-12, 1.0 - 1e-12);
            if rescue_fraction <= null_fraction {
                return 0.0;
            }
            let posterior_odds = rescue_fraction / (1.0 - rescue_fraction);
            let null_odds = null_fraction / (1.0 - null_fraction);
            let enrichment = posterior_odds / null_odds;
            if enrichment.is_finite() && enrichment > 0.0 {
                enrichment / (1.0 + enrichment)
            } else {
                0.0
            }
        };
    let enrichment_credible_floor_lock_fraction =
        |rescue_fraction: f64, eq_count: f64, rescued_targets: usize, ec_size: usize| {
            if rescue_fraction <= 0.0 || eq_count <= 0.0 {
                return 0.0;
            }
            if rescued_targets == 0 || rescued_targets >= ec_size || ec_size == 0 {
                return floor_smooth_confidence_lock_fraction(rescue_fraction);
            }
            let null_fraction = (rescued_targets as f64 / ec_size as f64).clamp(1e-9, 1.0 - 1e-9);
            let variance = (rescue_fraction * (1.0 - rescue_fraction) / eq_count).max(1e-12);
            let enrichment_posterior =
                normal_cdf((rescue_fraction - null_fraction) / variance.sqrt());
            if enrichment_posterior >= enrichment_posterior_threshold {
                floor_smooth_confidence_lock_fraction(rescue_fraction)
            } else {
                0.0
            }
        };

    for label_idx in 0..packed.len() {
        let label = packed.refs_for_eqc(label_idx);
        let eq_count = packed.counts[label_idx] as f64;
        let audit_label_selected = audit_targets
            .as_ref()
            .map(|targets| {
                label.target_labels().iter().any(|tid| {
                    let t = *tid as usize;
                    targets.contains(&ref_names[t]) || targets.contains(&t.to_string())
                })
            })
            .unwrap_or(true);
        if eq_count <= 0.0 {
            record_residual_eq_count(&mut residual_eq_counts, collapsed, label_idx, 0.0);
            continue;
        }
        weights.clear();
        let mut denom = 0.0f64;
        for (tid, cond_prob) in label.target_labels().iter().zip(label.target_probs()) {
            let t = *tid as usize;
            let w = cond_prob * phase1_counts[t] * inv_eff_lens[t];
            weights.push(w);
            denom += w;
        }
        if denom <= 0.0 {
            if audit_label_selected && let Some(file) = audit.as_mut() {
                let _ = writeln!(
                    file,
                    "{sample_name}\t{label_idx}\t{eq_count:.6}\t{eq_count:.6}\t0.000000\t0.000000\t0.000000\t{}\t0\t0\t\t\t",
                    label.target_labels().len()
                );
            }
            record_residual_eq_count(&mut residual_eq_counts, collapsed, label_idx, eq_count);
            continue;
        }
        let mut locked_sum = 0.0f64;
        let rescue_weight_sum: f64 = label
            .target_labels()
            .iter()
            .zip(weights.iter())
            .filter_map(|(tid, &w)| {
                let t = *tid as usize;
                (condition_locked_mask[t] && w > 0.0).then_some(w)
            })
            .sum();
        let rescue_fraction = rescue_weight_sum / denom;
        let rescued_targets = label
            .target_labels()
            .iter()
            .filter(|tid| condition_locked_mask[**tid as usize])
            .count();
        let relaxed_ec_lock_fraction = match lock_mode {
            ConditionRescueLockMode::FloorSmoothConfidence => {
                floor_smooth_confidence_lock_fraction(rescue_fraction)
            }
            ConditionRescueLockMode::CredibleFloor => {
                credible_floor_lock_fraction(rescue_fraction, eq_count)
            }
            ConditionRescueLockMode::EnrichmentFloor => enrichment_lock_fraction(
                rescue_fraction,
                rescued_targets,
                label.target_labels().len(),
            ),
            ConditionRescueLockMode::EnrichmentCredibleFloor => {
                enrichment_credible_floor_lock_fraction(
                    rescue_fraction,
                    eq_count,
                    rescued_targets,
                    label.target_labels().len(),
                )
            }
            ConditionRescueLockMode::EnrichmentCredibleStability => {
                enrichment_credible_floor_lock_fraction(
                    rescue_fraction,
                    eq_count,
                    rescued_targets,
                    label.target_labels().len(),
                )
            }
            _ => confidence_lock_fraction(rescue_fraction),
        };
        let any_lockable = match lock_mode {
            ConditionRescueLockMode::Fixed | ConditionRescueLockMode::FloorBarrier => {
                lock_fraction > 0.0
            }
            ConditionRescueLockMode::Confidence => relaxed_ec_lock_fraction > 0.0,
            ConditionRescueLockMode::FloorSmoothConfidence => relaxed_ec_lock_fraction > 0.0,
            ConditionRescueLockMode::CredibleFloor => relaxed_ec_lock_fraction > 0.0,
            ConditionRescueLockMode::EnrichmentFloor => relaxed_ec_lock_fraction > 0.0,
            ConditionRescueLockMode::EnrichmentCredibleFloor => relaxed_ec_lock_fraction > 0.0,
            ConditionRescueLockMode::EnrichmentCredibleStability => {
                label.target_labels().iter().any(|tid| {
                    let t = *tid as usize;
                    condition_locked_mask[t]
                        && condition_relax_mask.map(|mask| mask[t]).unwrap_or(false)
                        && relaxed_ec_lock_fraction > 0.0
                }) || label.target_labels().iter().any(|tid| {
                    let t = *tid as usize;
                    condition_locked_mask[t]
                        && !condition_relax_mask.map(|mask| mask[t]).unwrap_or(false)
                        && floor_smooth_confidence_lock_fraction(rescue_fraction) > 0.0
                })
            }
            ConditionRescueLockMode::GuardedConfidence => {
                label.target_labels().iter().any(|tid| {
                    let t = *tid as usize;
                    condition_locked_mask[t]
                        && condition_relax_mask.map(|mask| mask[t]).unwrap_or(false)
                        && relaxed_ec_lock_fraction > 0.0
                }) || label.target_labels().iter().any(|tid| {
                    let t = *tid as usize;
                    condition_locked_mask[t]
                        && !condition_relax_mask.map(|mask| mask[t]).unwrap_or(false)
                })
            }
            ConditionRescueLockMode::Instability => {
                label.target_labels().iter().any(|tid| {
                    let t = *tid as usize;
                    condition_locked_mask[t]
                        && condition_relax_mask.map(|mask| mask[t]).unwrap_or(false)
                        && relaxed_ec_lock_fraction > 0.0
                }) || label.target_labels().iter().any(|tid| {
                    let t = *tid as usize;
                    condition_locked_mask[t]
                        && !condition_relax_mask.map(|mask| mask[t]).unwrap_or(false)
                })
            }
            ConditionRescueLockMode::TranscriptStability => {
                label.target_labels().iter().any(|tid| {
                    let t = *tid as usize;
                    condition_locked_mask[t]
                        && condition_relax_mask.map(|mask| mask[t]).unwrap_or(false)
                        && relaxed_ec_lock_fraction > 0.0
                }) || label.target_labels().iter().any(|tid| {
                    let t = *tid as usize;
                    condition_locked_mask[t]
                        && !condition_relax_mask.map(|mask| mask[t]).unwrap_or(false)
                })
            }
        };
        if !any_lockable {
            if audit_label_selected && let Some(file) = audit.as_mut() {
                let (rescued_target_indices, rescued_target_names) = rescued_target_audit_fields(
                    label.target_labels(),
                    condition_locked_mask,
                    ref_names,
                );
                let _ = writeln!(
                    file,
                    "{sample_name}\t{label_idx}\t{eq_count:.6}\t{eq_count:.6}\t0.000000\t{rescue_fraction:.6}\t0.000000\t{}\t{rescued_targets}\t0\t{rescued_target_indices}\t{rescued_target_names}\t",
                    label.target_labels().len(),
                );
            }
            record_residual_eq_count(&mut residual_eq_counts, collapsed, label_idx, eq_count);
            continue;
        }
        let mut locked_positive_targets = 0usize;
        let mut locked_target_indices = audit
            .as_ref()
            .filter(|_| audit_label_selected)
            .map(|_| Vec::new());
        for (tid, &w) in label.target_labels().iter().zip(weights.iter()) {
            let t = *tid as usize;
            if condition_locked_mask[t] && w > 0.0 {
                let target_lock_fraction = match lock_mode {
                    ConditionRescueLockMode::Fixed | ConditionRescueLockMode::FloorBarrier => {
                        lock_fraction
                    }
                    ConditionRescueLockMode::Confidence => relaxed_ec_lock_fraction,
                    ConditionRescueLockMode::FloorSmoothConfidence => relaxed_ec_lock_fraction,
                    ConditionRescueLockMode::CredibleFloor => relaxed_ec_lock_fraction,
                    ConditionRescueLockMode::EnrichmentFloor => relaxed_ec_lock_fraction,
                    ConditionRescueLockMode::EnrichmentCredibleFloor => relaxed_ec_lock_fraction,
                    ConditionRescueLockMode::EnrichmentCredibleStability => {
                        if condition_relax_mask.map(|mask| mask[t]).unwrap_or(false) {
                            relaxed_ec_lock_fraction
                        } else {
                            floor_smooth_confidence_lock_fraction(rescue_fraction)
                        }
                    }
                    ConditionRescueLockMode::GuardedConfidence => {
                        if condition_relax_mask.map(|mask| mask[t]).unwrap_or(false) {
                            relaxed_ec_lock_fraction
                        } else {
                            1.0
                        }
                    }
                    ConditionRescueLockMode::Instability => {
                        if condition_relax_mask.map(|mask| mask[t]).unwrap_or(false) {
                            relaxed_ec_lock_fraction
                        } else {
                            1.0
                        }
                    }
                    ConditionRescueLockMode::TranscriptStability => {
                        if condition_relax_mask.map(|mask| mask[t]).unwrap_or(false) {
                            relaxed_ec_lock_fraction
                        } else {
                            1.0
                        }
                    }
                };
                if target_lock_fraction <= 0.0 {
                    continue;
                }
                let locked = target_lock_fraction * eq_count * w / denom;
                if locked > 0.0 {
                    locked_positive_targets += 1;
                    if let Some(indices) = locked_target_indices.as_mut() {
                        indices.push(t.to_string());
                    }
                }
                locked_counts[t] += locked;
                locked_sum += locked;
            }
        }
        let residual_count = (eq_count - locked_sum).max(0.0);
        if audit_label_selected && let Some(file) = audit.as_mut() {
            let (rescued_target_indices, rescued_target_names) = rescued_target_audit_fields(
                label.target_labels(),
                condition_locked_mask,
                ref_names,
            );
            let _ = writeln!(
                file,
                "{sample_name}\t{label_idx}\t{eq_count:.6}\t{residual_count:.6}\t{locked_sum:.6}\t{rescue_fraction:.6}\t{relaxed_ec_lock_fraction:.6}\t{}\t{rescued_targets}\t{locked_positive_targets}\t{rescued_target_indices}\t{rescued_target_names}\t{}",
                label.target_labels().len(),
                locked_target_indices
                    .as_ref()
                    .map(|indices| indices.join(","))
                    .unwrap_or_default()
            );
        }
        record_residual_eq_count(
            &mut residual_eq_counts,
            collapsed,
            label_idx,
            residual_count,
        );
    }

    (residual_eq_counts, locked_counts)
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

#[allow(clippy::too_many_arguments)]
fn locked_rescue_phase2_for_sample<L: EqLabel>(
    sample_name: &str,
    ref_names: &[String],
    packed: &PackedEqMap<L>,
    collapsed: Option<&CollapsedEqMap>,
    collapsed_index: Option<&crate::utils::txp_selection::TranscriptEqIndex>,
    full_index: &crate::utils::txp_selection::TranscriptEqIndex,
    eff_lengths: &[f64],
    phase1_counts: &[f64],
    init_counts: Option<&[f64]>,
    active_mask: &[bool],
    strict_global_mask: &[bool],
    condition_rescue_only: &[bool],
    condition_locked_mask: &[bool],
    condition_relax_mask: Option<&[bool]>,
    opts: &ConsensusQuantOpts,
) -> Vec<f64> {
    let floor_barrier_mode = matches!(
        opts.condition_rescue_lock_mode,
        ConditionRescueLockMode::FloorBarrier
    );
    let fully_locked = matches!(
        opts.condition_rescue_lock_mode,
        ConditionRescueLockMode::Fixed
    ) && opts.condition_rescue_lock_fraction >= 1.0;
    let free_mask: Vec<bool> = active_mask
        .iter()
        .zip(strict_global_mask.iter())
        .zip(condition_locked_mask.iter())
        .zip(phase1_counts.iter())
        .map(|(((&active, &global), &locked), &phase1_count)| {
            condition_rescue_free_variable(
                active,
                active && !global,
                locked,
                phase1_count,
                opts.condition_rescue_free_min_count,
                floor_barrier_mode,
                fully_locked,
                &opts.condition_rescue_lock_mode,
            )
        })
        .collect();
    let (residual_eq_counts, locked_counts) = lock_condition_rescue_allocations_for_sample(
        sample_name,
        ref_names,
        packed,
        phase1_counts,
        eff_lengths,
        condition_locked_mask,
        condition_relax_mask,
        collapsed,
        if floor_barrier_mode {
            1.0
        } else {
            opts.condition_rescue_lock_fraction
        },
        if floor_barrier_mode {
            &ConditionRescueLockMode::Fixed
        } else {
            &opts.condition_rescue_lock_mode
        },
        opts.condition_rescue_full_lock_threshold,
        opts.condition_rescue_min_lock_threshold,
        opts.condition_rescue_credible_floor_z,
        opts.condition_rescue_enrichment_posterior_threshold,
    );
    let free_init = init_counts.map(|counts| phase2_init_counts(counts, &free_mask));
    let residual_caps = condition_rescue_residual_cap_enabled(opts).then(|| {
        condition_rescue_residual_caps(
            condition_rescue_only,
            phase1_counts,
            &locked_counts,
            active_mask,
            opts.condition_rescue_residual_cap_z,
        )
    });
    let mut em_res = if let Some(cm) = collapsed {
        let collapsed_counts = if floor_barrier_mode {
            cm.packed
                .counts
                .iter()
                .map(|&count| count as f64)
                .collect::<Vec<_>>()
        } else {
            residual_eq_counts
        };
        phase2_em_step_f64_counts(
            &cm.packed,
            &collapsed_counts,
            eff_lengths,
            &free_mask,
            free_init.as_deref(),
            opts,
            floor_barrier_mode.then_some((
                locked_counts.as_slice(),
                opts.condition_rescue_floor_barrier_weight,
            )),
            residual_caps.as_deref(),
            condition_rescue_residual_cap_active_set_enabled(opts),
            collapsed_index,
        )
    } else {
        let full_counts = if floor_barrier_mode {
            packed
                .counts
                .iter()
                .map(|&count| count as f64)
                .collect::<Vec<_>>()
        } else {
            residual_eq_counts
        };
        phase2_em_step_f64_counts(
            packed,
            &full_counts,
            eff_lengths,
            &free_mask,
            free_init.as_deref(),
            opts,
            floor_barrier_mode.then_some((
                locked_counts.as_slice(),
                opts.condition_rescue_floor_barrier_weight,
            )),
            residual_caps.as_deref(),
            condition_rescue_residual_cap_active_set_enabled(opts),
            Some(full_index),
        )
    };
    if !floor_barrier_mode {
        for (count, locked) in em_res.iter_mut().zip(locked_counts.iter()) {
            *count += *locked;
        }
    }
    em_res
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
    let mut selection_keep_mask: Option<Vec<bool>> = None;
    let mut selection_result: Option<txp_selection::SelectionResult> = None;
    let original_eff_lengths: Vec<Vec<f64>> = bundles
        .iter()
        .map(|bundle| bundle.eff_lengths.clone())
        .collect();
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
        if let Some(audit_path) = &opts.selection_audit_output {
            write_selection_audit(audit_path, &bundles[0].ref_names, &result, &merged_index)
                .with_context(|| {
                    format!(
                        "failed to write structural selection audit to {}",
                        audit_path.display()
                    )
                })?;
            info!(
                "Wrote structural selection audit to {}",
                audit_path.display()
            );
        }
        // Zero out effective lengths for removed transcripts in all bundles.
        for bundle in &mut bundles {
            for (t, &keep) in result.keep_mask.iter().enumerate() {
                if !keep {
                    bundle.eff_lengths[t] = 0.0;
                }
            }
        }
        let num_kept = result.num_kept;
        let num_removed = result.num_removed;
        selection_keep_mask = Some(result.keep_mask.clone());
        selection_result = Some(result);
        Some((num_kept, num_removed))
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
    let is_range_factorized =
        std::any::TypeId::of::<EqLabelT>() == std::any::TypeId::of::<RangeFactorizedEqLabel>();
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
    let collapsed_indices: Vec<Option<txp_selection::TranscriptEqIndex>> = collapsed_maps
        .iter()
        .map(|cm| {
            cm.as_ref().map(|cm| {
                txp_selection::TranscriptEqIndex::from_packed_eq_map(&cm.packed, n_targets)
            })
        })
        .collect();

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
    let lock_condition_rescue_allocations =
        lock_condition_rescue_allocations_enabled(opts, use_condition_rescue);
    if opts.lock_condition_rescue_allocations && !use_condition_rescue {
        bail!(
            "--lock-condition-rescue-allocations requires condition rescue; provide multiple manifest conditions or pass --condition-rescue"
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

    let mut structured_phase2_condition_masks: Option<Vec<Vec<bool>>> = None;
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
        let (condition_names, cond_counts, cond_k, cond_support_masks) = condition_data
            .as_ref()
            .expect("condition data required for condition-rescue");

        if opts.structured_condition_rescue {
            let structured = structured_condition_rescue_mask(StructuredRescueInputs {
                samples,
                condition_names,
                cond_k,
                raw_condition_support_masks: cond_support_masks,
                sample_pass_masks: &sample_pass_masks,
                strict_global_mask: &strict_global_mask,
                phase1_tpms: &phase1_tpms,
                index: &merged_index,
                ref_names: &bundles[0].ref_names,
                opts,
            });
            let n_global = strict_global_mask.iter().filter(|&&b| b).count();
            let suppressed = structured
                .raw_rescued
                .saturating_sub(structured.selected_rescued);
            info!(
                "Structured condition rescue: {} pass global, {} raw rescue candidates, {} selected across {} admitted EC-graph groups, {} suppressed",
                n_global,
                structured.raw_rescued,
                structured.selected_rescued,
                structured.admitted_groups,
                suppressed
            );
            structured_phase2_condition_masks = Some(structured.condition_support_masks);
            structured.mask
        } else {
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
        }
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

    if let Some(audit_path) = &opts.consensus_audit_output {
        if let Some(parent) = audit_path.parent()
            && !parent.as_os_str().is_empty()
        {
            create_dir_all(parent).with_context(|| {
                format!(
                    "failed to create consensus audit output directory {}",
                    parent.display()
                )
            })?;
        }

        let mut audit = File::create(audit_path).with_context(|| {
            format!(
                "failed to create consensus audit output {}",
                audit_path.display()
            )
        })?;
        writeln!(
            audit,
            "target_idx\ttarget_name\tstructurally_kept\tphase1_count_sum\tphase1_tpm_mean\tsample_pass_count\tmin_k\tstrict_global\tcondition_counts\tcondition_rescue\tpre_gene_consensus\tgene_rescue\tfinal_consensus\tdecision"
        )?;

        let structural_keep = selection_keep_mask.as_deref();
        for t in 0..n_targets {
            let structurally_kept = structural_keep.map(|m| m[t]).unwrap_or(true);
            let phase1_count_sum: f64 = phase1_counts.iter().map(|counts| counts[t]).sum();
            let phase1_tpm_mean: f64 =
                phase1_tpms.iter().map(|tpms| tpms[t]).sum::<f64>() / n_samples as f64;
            let condition_counts = condition_data
                .as_ref()
                .map(|(names, counts, cond_k, _)| {
                    names
                        .iter()
                        .enumerate()
                        .map(|(ci, name)| format!("{}:{}/{}", name, counts[ci][t], cond_k[ci]))
                        .collect::<Vec<_>>()
                        .join(";")
                })
                .unwrap_or_else(|| ".".to_string());
            let condition_rescue = pre_gene_consensus_mask[t] && !strict_global_mask[t];
            let gene_rescue = consensus_mask[t] && !pre_gene_consensus_mask[t];
            let decision = if !structurally_kept {
                "structural_selection_removed"
            } else if strict_global_mask[t] {
                "strict_global"
            } else if condition_rescue {
                "condition_rescue"
            } else if gene_rescue {
                "gene_rescue"
            } else if consensus_mask[t] {
                "final_consensus"
            } else if express_count[t] > 0 {
                "failed_consensus"
            } else {
                "no_sample_pass"
            };
            writeln!(
                audit,
                "{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                t,
                bundles[0].ref_names[t],
                structurally_kept,
                phase1_count_sum,
                phase1_tpm_mean,
                express_count[t],
                min_k,
                strict_global_mask[t],
                condition_counts,
                condition_rescue,
                pre_gene_consensus_mask[t],
                gene_rescue,
                consensus_mask[t],
                decision
            )?;
        }
        info!("Wrote consensus audit to {}", audit_path.display());
    }

    let phase2_condition_support_masks =
        structured_phase2_condition_masks.as_deref().or_else(|| {
            condition_data
                .as_ref()
                .map(|(_, _, _, support_masks)| support_masks.as_slice())
        });
    let phase2_sample_masks = sample_specific_phase2_masks(
        samples,
        &strict_global_mask,
        &pre_gene_consensus_mask,
        condition_data
            .as_ref()
            .map(|(names, _, _, _)| names.as_slice()),
        phase2_condition_support_masks,
        opts.condition_aware_consensus,
        use_condition_rescue,
    );
    let condition_rescue_only: Vec<bool> = pre_gene_consensus_mask
        .iter()
        .zip(strict_global_mask.iter())
        .map(|(&pre_gene_keep, &strict_keep)| pre_gene_keep && !strict_keep)
        .collect();
    let condition_rescue_lock_relax_masks: Option<Vec<Vec<bool>>> =
        if lock_condition_rescue_allocations
            && matches!(
                opts.condition_rescue_lock_mode,
                ConditionRescueLockMode::GuardedConfidence
                    | ConditionRescueLockMode::Instability
                    | ConditionRescueLockMode::TranscriptStability
                    | ConditionRescueLockMode::EnrichmentCredibleStability
            )
        {
            condition_data
            .as_ref()
            .map(|(condition_names, _, _, condition_support_masks)| {
                let masks = match opts.condition_rescue_lock_mode {
                    ConditionRescueLockMode::GuardedConfidence => {
                        condition_rescue_guarded_confidence_masks(
                            samples,
                            condition_names,
                            &condition_rescue_only,
                            condition_support_masks,
                            &phase1_counts,
                            opts.condition_rescue_guard_mean_count_threshold,
                        )
                    }
                    ConditionRescueLockMode::Instability => condition_rescue_instability_masks(
                        samples,
                        condition_names,
                        &condition_rescue_only,
                        condition_support_masks,
                        &rescue_sample_pass_masks,
                        &phase1_counts,
                        opts.condition_rescue_instability_cv_threshold,
                        opts.condition_rescue_instability_min_pass_fraction,
                    ),
                    ConditionRescueLockMode::TranscriptStability
                    | ConditionRescueLockMode::EnrichmentCredibleStability => {
                        condition_rescue_transcript_stability_masks(
                            samples,
                            condition_names,
                            &condition_rescue_only,
                            condition_support_masks,
                            &rescue_sample_pass_masks,
                            &phase1_counts,
                            opts.condition_rescue_stability_mean_count_threshold,
                            opts.condition_rescue_stability_cv_threshold,
                            opts.condition_rescue_stability_min_pass_fraction,
                        )
                    }
                    _ => unreachable!("relax masks are only used by guarded lock modes"),
                };
                let n_relaxable: usize = masks
                    .iter()
                    .map(|mask| mask.iter().filter(|&&relaxable| relaxable).count())
                    .sum();
                match opts.condition_rescue_lock_mode {
                    ConditionRescueLockMode::GuardedConfidence => info!(
                        "Condition-rescue guarded-confidence lock: {} condition/transcript pairs marked relaxable (mean_count_threshold={:.3})",
                        n_relaxable,
                        opts.condition_rescue_guard_mean_count_threshold
                    ),
                    ConditionRescueLockMode::Instability => info!(
                        "Condition-rescue instability lock: {} condition/transcript pairs marked relaxable (cv_threshold={:.3}, min_pass_fraction={:.3})",
                        n_relaxable,
                        opts.condition_rescue_instability_cv_threshold,
                        opts.condition_rescue_instability_min_pass_fraction
                    ),
                    ConditionRescueLockMode::TranscriptStability => info!(
                        "Condition-rescue transcript-stability lock: {} condition/transcript pairs marked relaxable (mean_count_threshold={:.3}, cv_threshold={:.3}, min_pass_fraction={:.3})",
                        n_relaxable,
                        opts.condition_rescue_stability_mean_count_threshold,
                        opts.condition_rescue_stability_cv_threshold,
                        opts.condition_rescue_stability_min_pass_fraction
                    ),
                    ConditionRescueLockMode::EnrichmentCredibleStability => info!(
                        "Condition-rescue enrichment-credible-stability lock: {} condition/transcript pairs require enrichment gate (mean_count_threshold={:.3}, cv_threshold={:.3}, min_pass_fraction={:.3})",
                        n_relaxable,
                        opts.condition_rescue_stability_mean_count_threshold,
                        opts.condition_rescue_stability_cv_threshold,
                        opts.condition_rescue_stability_min_pass_fraction
                    ),
                    _ => {}
                }
                masks
            })
        } else {
            None
        };
    let condition_rescue_lock_condition_indices: Option<Vec<usize>> =
        condition_rescue_lock_relax_masks.as_ref().map(|_| {
            let condition_names = condition_data
                .as_ref()
                .map(|(names, _, _, _)| names.as_slice())
                .expect("condition data required for condition rescue guarded lock");
            sample_condition_indices(samples, condition_names)
        });
    let phase2_prior_alpha: Option<Vec<Vec<f64>>> = if opts.condition_specific_prior_weight > 0.0 {
        let avg_total_reads: f64 = bundles
            .iter()
            .map(|b| b.packed_eq_map.total_weight() as f64)
            .sum::<f64>()
            / n_samples as f64;
        let alpha_0 = opts.condition_specific_prior_weight * avg_total_reads;
        let condition_names = sorted_condition_names(samples);
        let condition_indices = sample_condition_indices(samples, &condition_names);
        let mut hyperparams = hierarchical::init_hyperparams(n_targets, condition_names, alpha_0);
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
                    let em_res = if lock_condition_rescue_allocations {
                        let condition_relax_mask =
                            condition_rescue_lock_relax_masks.as_ref().map(|masks| {
                                let condition_idx =
                                    condition_rescue_lock_condition_indices.as_ref().unwrap()[i];
                                masks[condition_idx].as_slice()
                            });
                        locked_rescue_phase2_for_sample(
                            &sample.sample_name,
                            &bundles[i].ref_names,
                            &bundles[i].packed_eq_map,
                            collapsed_maps[i].as_ref(),
                            collapsed_indices[i].as_ref(),
                            &indices[i],
                            &bundles[i].eff_lengths,
                            &phase1_counts[i],
                            init.as_deref(),
                            phase2_mask,
                            &strict_global_mask,
                            &condition_rescue_only,
                            &condition_rescue_only,
                            condition_relax_mask,
                            opts,
                        )
                    } else if let Some(cm) = collapsed_maps[i].as_ref() {
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
            let em_res = if lock_condition_rescue_allocations {
                let condition_relax_mask =
                    condition_rescue_lock_relax_masks.as_ref().map(|masks| {
                        let condition_idx =
                            condition_rescue_lock_condition_indices.as_ref().unwrap()[i];
                        masks[condition_idx].as_slice()
                    });
                locked_rescue_phase2_for_sample(
                    &sample.sample_name,
                    &bundles[i].ref_names,
                    &bundles[i].packed_eq_map,
                    collapsed_maps[i].as_ref(),
                    collapsed_indices[i].as_ref(),
                    &indices[i],
                    &bundles[i].eff_lengths,
                    &phase1_counts[i],
                    init.as_deref(),
                    phase2_mask,
                    &strict_global_mask,
                    &condition_rescue_only,
                    &condition_rescue_only,
                    condition_relax_mask,
                    opts,
                )
            } else if let Some(cm) = collapsed_maps[i].as_ref() {
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
        for (t, &active) in consensus_mask.iter().enumerate().take(n_targets) {
            if !active {
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

        for (t, &active) in consensus_mask.iter().enumerate().take(n_targets) {
            if !active {
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
        if std::env::var("PISCEM_POS_DIAG").is_ok()
            && let Ok(mut f) = std::fs::File::create("pos_diagnostic.tsv")
        {
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
            for (t, &active) in consensus_mask.iter().enumerate().take(n_targets) {
                if !active {
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

    let mut ambiguity_repair_report = if let Some(report_path) = &opts.ambiguity_repair_report {
        if let Some(parent) = report_path.parent()
            && !parent.as_os_str().is_empty()
        {
            create_dir_all(parent)?;
        }
        let mut file = File::create(report_path).with_context(|| {
            format!(
                "failed to create ambiguity repair report at {}",
                report_path.display()
            )
        })?;
        writeln!(
            file,
            "sample_name\taction\tgroup_id\ttarget_idx\ttarget_name\trelated_idx\trelated_name\told_count\tgroup_total\tgroup_fraction\tphase1_count\tprivate_ec_fraction\tposition_effective_bins\tmax_position_bin_fraction\tshared_ec_count\trelated_total_ec_count\trelated_extra_ec_count\trelated_shared_fraction\tlocal_candidate_count\tlocal_related_count\tlocal_candidate_fraction\treason"
        )?;
        Some(file)
    } else {
        None
    };

    let mut structural_repair_candidates: Vec<StructuralRepairCandidate> = Vec::new();
    if let (Some(file), Some(selection), Some(structural_keep)) = (
        ambiguity_repair_report.as_mut(),
        selection_result.as_ref(),
        selection_keep_mask.as_deref(),
    ) {
        let mut target_to_group = vec![usize::MAX; n_targets];
        for gi in 0..selection.groups.num_groups() {
            for &member in selection.groups.group_members(gi) {
                target_to_group[member as usize] = gi;
            }
        }

        let max_eqc = merged_index
            .eqc_ids
            .iter()
            .copied()
            .max()
            .map(|m| m as usize + 1)
            .unwrap_or(0);
        let mut eqc_to_groups: Vec<Vec<u32>> = vec![Vec::new(); max_eqc];
        for gi in 0..selection.groups.num_groups() {
            let rep = selection.groups.representatives[gi] as usize;
            for &eqc in merged_index.signature(rep) {
                eqc_to_groups[eqc as usize].push(gi as u32);
            }
        }

        let mut group_kept = vec![false; selection.groups.num_groups()];
        for (t, &keep) in structural_keep.iter().enumerate() {
            let gi = target_to_group[t];
            if keep && gi != usize::MAX {
                group_kept[gi] = true;
            }
        }

        let report_pos_bins = n_pos_bins.max(1);
        let mut n_structural_candidates = 0usize;
        for t in 0..n_targets {
            if structural_keep[t] {
                continue;
            }
            let gi = target_to_group[t];
            if gi == usize::MAX || !selection.dominated[gi] {
                continue;
            }
            let Some(dom_gi) =
                find_selection_dominator(gi, &group_kept, &eqc_to_groups, selection, &merged_index)
            else {
                continue;
            };
            let dom_t = selection.groups.representatives[dom_gi] as usize;
            let shared_count = merged_index.total_count(t);
            let dom_total = merged_index.total_count(dom_t);
            if shared_count < 100 || dom_total == 0 {
                continue;
            }
            let dom_extra = dom_total.saturating_sub(shared_count);
            let shared_fraction = shared_count as f64 / dom_total as f64;
            if dom_extra > 100 && shared_fraction < 0.90 {
                continue;
            }
            let pos_profile = transcript_position_profile(&merged_index, t, report_pos_bins);
            let effective_bins = position_effective_bins(&pos_profile);
            let max_bin_frac = max_position_bin_fraction(&pos_profile);
            let private_ec_fraction = 0.0;
            writeln!(
                file,
                "NA\tstructural_injection_candidate\t{}\t{}\t{}\t{}\t{}\t0.000000\t0.000000\t0.000000\t0.000000\t{:.6}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{:.6}\tNA\tNA\tNA\tweak_dominator_private_support",
                gi,
                t,
                bundles[0].ref_names[t],
                dom_t,
                bundles[0].ref_names[dom_t],
                private_ec_fraction,
                effective_bins,
                max_bin_frac,
                shared_count,
                dom_total,
                dom_extra,
                shared_fraction
            )?;
            structural_repair_candidates.push(StructuralRepairCandidate {
                target: t,
                dominator: dom_t,
                group_idx: gi,
                shared_count,
                dominator_total_count: dom_total,
                dominator_extra_count: dom_extra,
                shared_fraction,
            });
            n_structural_candidates += 1;
        }
        info!(
            "Ambiguity repair report: wrote {} structural injection candidates",
            n_structural_candidates
        );
    }

    let gene_frac_threshold = opts.gene_fraction_filter;
    let has_gene_annot = !gene_to_txps.is_empty();
    let has_pos = !pos_bin_profiles.is_empty() && n_pos_bins > 1;
    let mut total_gene_frac_removed = 0usize;
    let mut total_nbr_removed = 0usize;
    let leakage_audit_path = std::env::var("PISCEM_LEAKAGE_AUDIT").ok();
    let mut leakage_audit_rows: Vec<String> = Vec::new();
    let gene_rescue_only: Vec<bool> = consensus_mask
        .iter()
        .zip(pre_gene_consensus_mask.iter())
        .map(|(&final_keep, &pre_gene_keep)| final_keep && !pre_gene_keep)
        .collect();
    let post_filter_redistribute_mass = !opts.no_post_filter_redistribute_mass;
    for (i, mut em_res) in phase2_results {
        let final_em_mode = post_filter_redistribute_mass
            && matches!(
                opts.post_filter_redistribute_mode,
                PostFilterRedistributeMode::FinalEm
            );
        let mut post_filter_mask = if final_em_mode {
            Some(phase2_sample_masks[i].clone())
        } else {
            None
        };
        let responsibility_index = if post_filter_redistribute_mass
            && matches!(
                opts.post_filter_redistribute_mode,
                PostFilterRedistributeMode::SharedResponsibility
            ) {
            Some(&indices[i])
        } else {
            None
        };
        let remove_leakage = |em_res: &mut [f64],
                              post_filter_mask: &mut Option<Vec<bool>>,
                              target: usize,
                              recipient: usize,
                              candidates: &[usize]| {
            if post_filter_redistribute_mass {
                let removed = em_res[target];
                if removed.is_finite() && removed > 0.0 {
                    match opts.post_filter_redistribute_mode {
                        PostFilterRedistributeMode::Dominator => {
                            if target != recipient {
                                em_res[recipient] += removed;
                            }
                        }
                        PostFilterRedistributeMode::SharedPosterior => {
                            let shared: Vec<usize> = candidates
                                .iter()
                                .copied()
                                .filter(|&u| {
                                    u != target
                                        && em_res[u].is_finite()
                                        && em_res[u] > 0.0
                                        && signatures_intersect(&merged_index, target, u)
                                })
                                .collect();
                            let denom: f64 = shared.iter().map(|&u| em_res[u]).sum();
                            if denom > 0.0 {
                                for u in shared {
                                    em_res[u] += removed * (em_res[u] / denom);
                                }
                            } else if target != recipient {
                                em_res[recipient] += removed;
                            }
                        }
                        PostFilterRedistributeMode::SharedResponsibility => {
                            if let Some(index) = responsibility_index.as_ref() {
                                let ctx = SharedResponsibilityContext {
                                    packed: &bundles[i].packed_eq_map,
                                    index,
                                    eff_lengths: &bundles[i].eff_lengths,
                                };
                                redistribute_by_shared_responsibility(
                                    &ctx, em_res, target, candidates, removed, recipient,
                                );
                            } else if target != recipient {
                                em_res[recipient] += removed;
                            }
                        }
                        PostFilterRedistributeMode::FinalEm => {
                            if let Some(mask) = post_filter_mask.as_mut() {
                                mask[target] = false;
                            }
                        }
                    }
                }
            }
            em_res[target] = 0.0;
        };

        // Rescue-only transcripts retain their phase-1 estimates after phase 2.
        for t in 0..n_targets {
            if gene_rescue_only[t]
                || (condition_rescue_only[t]
                    && !opts.reestimate_condition_rescue
                    && !lock_condition_rescue_allocations)
            {
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
                    if opts.preserve_condition_rescue_leakage && condition_rescue_only[t] {
                        continue;
                    }
                    let gene_frac = em_res[t] / gene_total;

                    // Gene-fraction filter
                    if gene_frac_threshold > 0.0 && gene_frac < gene_frac_threshold {
                        let pos_uneven = if has_pos { pos_cv[t] > 0.5 } else { false };
                        if ec_unique_frac[t] < 1.0 || pos_uneven {
                            if opts.condition_rescue_leakage_phase1_floor_count > 0.0
                                && condition_rescue_only[t]
                                && phase1_counts[i][t]
                                    >= opts.condition_rescue_leakage_phase1_floor_count
                            {
                                if leakage_audit_path.is_some() {
                                    leakage_audit_rows.push(format!(
                                        "{}\t{}\t{}\tgene_fraction_phase1_floor_exempt\t{:.6}\t{:.6}\t{}\t{}\t{:.6}\t0.000000",
                                        samples[i].sample_name,
                                        t,
                                        bundles[i].ref_names[t],
                                        em_res[t],
                                        gene_frac,
                                        dom_t,
                                        bundles[i].ref_names[dom_t],
                                        phase1_counts[i][t]
                                    ));
                                }
                                continue;
                            }
                            if opts.selective_condition_rescue_leakage
                                && condition_rescue_only[t]
                                && let Some((pair_count, pair_frac)) =
                                    condition_rescue_leakage_exempt(
                                        &bundles[i].packed_eq_map,
                                        &merged_index,
                                        &eqc_offset_boundaries,
                                        i,
                                        t,
                                        dom_t,
                                        &original_eff_lengths[i],
                                        &pos_bin_profiles,
                                        n_pos_bins,
                                        opts.selective_condition_rescue_leakage_min_count,
                                        opts.selective_condition_rescue_leakage_min_fraction,
                                    )
                            {
                                if leakage_audit_path.is_some() {
                                    leakage_audit_rows.push(format!(
                                        "{}\t{}\t{}\tgene_fraction_exempt\t{:.6}\t{:.6}\t{}\t{}\t{:.6}\t{:.6}",
                                        samples[i].sample_name,
                                        t,
                                        bundles[i].ref_names[t],
                                        em_res[t],
                                        gene_frac,
                                        dom_t,
                                        bundles[i].ref_names[dom_t],
                                        pair_count,
                                        pair_frac
                                    ));
                                }
                                continue;
                            }
                            if leakage_audit_path.is_some() {
                                leakage_audit_rows.push(format!(
                                    "{}\t{}\t{}\tgene_fraction\t{:.6}\t{:.6}\t{}\t{}",
                                    samples[i].sample_name,
                                    t,
                                    bundles[i].ref_names[t],
                                    em_res[t],
                                    gene_frac,
                                    dom_t,
                                    bundles[i].ref_names[dom_t]
                                ));
                            }
                            remove_leakage(&mut em_res, &mut post_filter_mask, t, dom_t, txps);
                            if i == 0 {
                                total_gene_frac_removed += 1;
                            }
                            continue;
                        }
                    }

                    // Profile-correlation filter (gene-annotated version)
                    if has_pos
                        && gene_frac < 0.05
                        && profile_corr_is_leakage(&pos_bin_profiles, t, dom_t, n_pos_bins)
                    {
                        if leakage_audit_path.is_some() {
                            leakage_audit_rows.push(format!(
                                "{}\t{}\t{}\tgene_profile\t{:.6}\t{:.6}\t{}\t{}",
                                samples[i].sample_name,
                                t,
                                bundles[i].ref_names[t],
                                em_res[t],
                                gene_frac,
                                dom_t,
                                bundles[i].ref_names[dom_t]
                            ));
                        }
                        remove_leakage(&mut em_res, &mut post_filter_mask, t, dom_t, txps);
                        if i == 0 {
                            total_gene_frac_removed += 1;
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
            for members in ec_graph_group_members.values() {
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
                    if opts.preserve_condition_rescue_leakage && condition_rescue_only[t] {
                        continue;
                    }

                    let group_frac = em_res[t] / group_total;

                    // Group-fraction filter (same logic as gene-fraction)
                    if group_frac < gene_frac_threshold.max(0.01) {
                        let pos_uneven = if has_pos { pos_cv[t] > 0.5 } else { false };
                        if ec_unique_frac[t] < 1.0 || pos_uneven {
                            if opts.condition_rescue_leakage_phase1_floor_count > 0.0
                                && condition_rescue_only[t]
                                && phase1_counts[i][t]
                                    >= opts.condition_rescue_leakage_phase1_floor_count
                            {
                                if leakage_audit_path.is_some() {
                                    leakage_audit_rows.push(format!(
                                        "{}\t{}\t{}\tec_group_fraction_phase1_floor_exempt\t{:.6}\t{:.6}\t{}\t{}\t{:.6}\t0.000000",
                                        samples[i].sample_name,
                                        t,
                                        bundles[i].ref_names[t],
                                        em_res[t],
                                        group_frac,
                                        dom_t,
                                        bundles[i].ref_names[dom_t],
                                        phase1_counts[i][t]
                                    ));
                                }
                                continue;
                            }
                            if opts.selective_condition_rescue_leakage
                                && condition_rescue_only[t]
                                && let Some((pair_count, pair_frac)) =
                                    condition_rescue_leakage_exempt(
                                        &bundles[i].packed_eq_map,
                                        &merged_index,
                                        &eqc_offset_boundaries,
                                        i,
                                        t,
                                        dom_t,
                                        &original_eff_lengths[i],
                                        &pos_bin_profiles,
                                        n_pos_bins,
                                        opts.selective_condition_rescue_leakage_min_count,
                                        opts.selective_condition_rescue_leakage_min_fraction,
                                    )
                            {
                                if leakage_audit_path.is_some() {
                                    leakage_audit_rows.push(format!(
                                        "{}\t{}\t{}\tec_group_fraction_exempt\t{:.6}\t{:.6}\t{}\t{}\t{:.6}\t{:.6}",
                                        samples[i].sample_name,
                                        t,
                                        bundles[i].ref_names[t],
                                        em_res[t],
                                        group_frac,
                                        dom_t,
                                        bundles[i].ref_names[dom_t],
                                        pair_count,
                                        pair_frac
                                    ));
                                }
                                continue;
                            }
                            if leakage_audit_path.is_some() {
                                leakage_audit_rows.push(format!(
                                    "{}\t{}\t{}\tec_group_fraction\t{:.6}\t{:.6}\t{}\t{}",
                                    samples[i].sample_name,
                                    t,
                                    bundles[i].ref_names[t],
                                    em_res[t],
                                    group_frac,
                                    dom_t,
                                    bundles[i].ref_names[dom_t]
                                ));
                            }
                            remove_leakage(&mut em_res, &mut post_filter_mask, t, dom_t, &active);
                            if i == 0 {
                                total_nbr_removed += 1;
                            }
                            continue;
                        }
                    }

                    // Profile-correlation filter (group version)
                    if has_pos
                        && group_frac < 0.05
                        && profile_corr_is_leakage(&pos_bin_profiles, t, dom_t, n_pos_bins)
                    {
                        if leakage_audit_path.is_some() {
                            leakage_audit_rows.push(format!(
                                "{}\t{}\t{}\tec_group_profile\t{:.6}\t{:.6}\t{}\t{}",
                                samples[i].sample_name,
                                t,
                                bundles[i].ref_names[t],
                                em_res[t],
                                group_frac,
                                dom_t,
                                bundles[i].ref_names[dom_t]
                            ));
                        }
                        remove_leakage(&mut em_res, &mut post_filter_mask, t, dom_t, &active);
                        if i == 0 {
                            total_nbr_removed += 1;
                        }
                    }
                }
            }
        }

        if let Some(file) = ambiguity_repair_report.as_mut() {
            let mut n_local_rows = 0usize;
            for candidate in &structural_repair_candidates {
                let Some((candidate_count, related_count)) = local_pair_em(
                    &bundles[i].packed_eq_map,
                    &merged_index,
                    &eqc_offset_boundaries,
                    i,
                    candidate.target,
                    candidate.dominator,
                    &original_eff_lengths[i],
                ) else {
                    continue;
                };
                let local_total = candidate_count + related_count;
                if local_total <= 0.0 {
                    continue;
                }
                let candidate_fraction = candidate_count / local_total;
                if candidate_count < 1.0 || candidate_fraction < 0.01 {
                    continue;
                }
                let pos_profile =
                    transcript_position_profile(&merged_index, candidate.target, n_pos_bins.max(1));
                writeln!(
                    file,
                    "{}\tlocal_pair_em_candidate\t{}\t{}\t{}\t{}\t{}\t0.000000\t{:.6}\t{:.6}\t{:.6}\t0.000000\t{:.6}\t{:.6}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\tlocal_pair_em_nonzero_candidate",
                    samples[i].sample_name,
                    candidate.group_idx,
                    candidate.target,
                    bundles[i].ref_names[candidate.target],
                    candidate.dominator,
                    bundles[i].ref_names[candidate.dominator],
                    local_total,
                    candidate_fraction,
                    phase1_counts[i][candidate.target],
                    position_effective_bins(&pos_profile),
                    max_position_bin_fraction(&pos_profile),
                    candidate.shared_count,
                    candidate.dominator_total_count,
                    candidate.dominator_extra_count,
                    candidate.shared_fraction,
                    candidate_count,
                    related_count,
                    candidate_fraction
                )?;
                n_local_rows += 1;
            }
            if i == 0 {
                info!(
                    "Ambiguity repair report: wrote local pair-EM candidates for {} sample-1 rows",
                    n_local_rows
                );
            }

            for (&group_id, members) in &ec_graph_group_members {
                let active: Vec<usize> = members
                    .iter()
                    .filter(|&&t| em_res[t] > 0.0)
                    .copied()
                    .collect();
                if active.len() <= 1 {
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
                    let group_frac = em_res[t] / group_total;
                    let weak_private = ec_unique_frac[t] < 1.0;
                    let pos_uneven = if has_pos { pos_cv[t] > 0.5 } else { false };
                    if group_frac >= 0.20 || !(weak_private || pos_uneven) {
                        continue;
                    }
                    let (effective_bins, max_bin_frac) = if !pos_bin_profiles.is_empty() {
                        (
                            position_effective_bins(&pos_bin_profiles[t]),
                            max_position_bin_fraction(&pos_bin_profiles[t]),
                        )
                    } else {
                        (0.0, 0.0)
                    };
                    let reason = match (weak_private, pos_uneven) {
                        (true, true) => "weak_private_and_uneven_position",
                        (true, false) => "weak_private_support",
                        (false, true) => "uneven_position",
                        (false, false) => "none",
                    };
                    writeln!(
                        file,
                        "{}\tprune_candidate\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t{}",
                        samples[i].sample_name,
                        group_id,
                        t,
                        bundles[i].ref_names[t],
                        dom_t,
                        bundles[i].ref_names[dom_t],
                        em_res[t],
                        group_total,
                        group_frac,
                        phase1_counts[i][t],
                        ec_unique_frac[t],
                        effective_bins,
                        max_bin_frac,
                        reason
                    )?;
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
                    if !ecg_remove
                        && has_pos
                        && gf < 0.05
                        && profile_corr_is_leakage(&pos_bin_profiles, t, dom_t, n_pos_bins)
                    {
                        ecg_remove = true;
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

        if let Some(final_mask) = post_filter_mask.as_ref() {
            info!(
                "  [{}/{}] Final post-filter EM for {}",
                i + 1,
                n_samples,
                samples[i].sample_name
            );
            em_res = if lock_condition_rescue_allocations {
                let condition_locked_mask: Vec<bool> = condition_rescue_only
                    .iter()
                    .zip(final_mask.iter())
                    .map(|(&locked, &active)| locked && active)
                    .collect();
                let condition_relax_mask =
                    condition_rescue_lock_relax_masks.as_ref().map(|masks| {
                        let condition_idx =
                            condition_rescue_lock_condition_indices.as_ref().unwrap()[i];
                        masks[condition_idx].as_slice()
                    });
                locked_rescue_phase2_for_sample(
                    &samples[i].sample_name,
                    &bundles[i].ref_names,
                    &bundles[i].packed_eq_map,
                    collapsed_maps[i].as_ref(),
                    collapsed_indices[i].as_ref(),
                    &indices[i],
                    &bundles[i].eff_lengths,
                    &phase1_counts[i],
                    Some(&em_res),
                    final_mask,
                    &strict_global_mask,
                    &condition_rescue_only,
                    &condition_locked_mask,
                    condition_relax_mask,
                    opts,
                )
            } else if let Some(cm) = collapsed_maps[i].as_ref() {
                phase2_em_step(
                    &cm.packed,
                    bundles[i].eff_lengths.clone(),
                    final_mask,
                    Some(&em_res),
                    None,
                    opts,
                    inner_threads,
                    serial_inner_pool.as_ref(),
                )
            } else {
                phase2_em_step(
                    &bundles[i].packed_eq_map,
                    bundles[i].eff_lengths.clone(),
                    final_mask,
                    Some(&em_res),
                    None,
                    opts,
                    inner_threads,
                    serial_inner_pool.as_ref(),
                )
            };
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
            "structured_condition_rescue": opts.structured_condition_rescue,
            "structured_rescue_group_tpm_floor": opts.structured_rescue_group_tpm_floor,
            "structured_rescue_cumulative_frac": opts.structured_rescue_cumulative_frac,
            "structured_rescue_max_isoforms": opts.structured_rescue_max_isoforms,
            "structured_rescue_rank_mode": opts.structured_rescue_rank_mode,
            "structured_rescue_per_condition_representatives": opts.structured_rescue_per_condition_representatives,
            "structured_rescue_strict_group_escape": opts.structured_rescue_strict_group_escape,
            "structured_rescue_strict_group_candidate_tpm_floor": opts.structured_rescue_strict_group_candidate_tpm_floor,
            "structured_rescue_strict_group_balanced_escape": opts.structured_rescue_strict_group_balanced_escape,
            "structured_rescue_strict_group_balanced_tpm_floor": opts.structured_rescue_strict_group_balanced_tpm_floor,
            "structured_rescue_strict_group_balanced_min_conditions": opts.structured_rescue_strict_group_balanced_min_conditions,
            "structured_rescue_phase1_condition_complements": opts.structured_rescue_phase1_condition_complements,
            "structured_rescue_phase1_complement_tpm_floor": opts.structured_rescue_phase1_complement_tpm_floor,
            "reestimate_condition_rescue": opts.reestimate_condition_rescue,
            "lock_condition_rescue_allocations": lock_condition_rescue_allocations,
            "condition_rescue_lock_fraction": opts.condition_rescue_lock_fraction,
            "condition_rescue_lock_mode": opts.condition_rescue_lock_mode,
            "condition_rescue_floor_barrier_weight": opts.condition_rescue_floor_barrier_weight,
            "condition_rescue_full_lock_threshold": opts.condition_rescue_full_lock_threshold,
            "condition_rescue_min_lock_threshold": opts.condition_rescue_min_lock_threshold,
            "condition_rescue_instability_cv_threshold": opts.condition_rescue_instability_cv_threshold,
            "condition_rescue_instability_min_pass_fraction": opts.condition_rescue_instability_min_pass_fraction,
            "condition_rescue_guard_mean_count_threshold": opts.condition_rescue_guard_mean_count_threshold,
            "condition_rescue_stability_mean_count_threshold": opts.condition_rescue_stability_mean_count_threshold,
            "condition_rescue_stability_cv_threshold": opts.condition_rescue_stability_cv_threshold,
            "condition_rescue_stability_min_pass_fraction": opts.condition_rescue_stability_min_pass_fraction,
            "condition_specific_prior_weight": opts.condition_specific_prior_weight,
            "preserve_condition_rescue_leakage": opts.preserve_condition_rescue_leakage,
            "selective_condition_rescue_leakage": opts.selective_condition_rescue_leakage,
            "selective_condition_rescue_leakage_min_count": opts.selective_condition_rescue_leakage_min_count,
            "selective_condition_rescue_leakage_min_fraction": opts.selective_condition_rescue_leakage_min_fraction,
            "condition_rescue_leakage_phase1_floor_count": opts.condition_rescue_leakage_phase1_floor_count,
            "piscem_infer_version": env!("CARGO_PKG_VERSION"),
        });
        meta_info["post_filter_redistribute_mass"] = json!(post_filter_redistribute_mass);
        meta_info["post_filter_redistribute_mode"] = json!(opts.post_filter_redistribute_mode);
        meta_info["no_lock_condition_rescue_allocations"] =
            json!(opts.no_lock_condition_rescue_allocations);
        meta_info["condition_rescue_free_min_count"] = json!(opts.condition_rescue_free_min_count);
        meta_info["condition_rescue_residual_cap"] =
            json!(condition_rescue_residual_cap_enabled(opts));
        meta_info["condition_rescue_residual_cap_z"] = json!(opts.condition_rescue_residual_cap_z);
        meta_info["condition_rescue_residual_cap_active_set"] =
            json!(condition_rescue_residual_cap_active_set_enabled(opts));
        meta_info["condition_rescue_credible_floor_z"] =
            json!(opts.condition_rescue_credible_floor_z);
        meta_info["condition_rescue_enrichment_posterior_threshold"] =
            json!(opts.condition_rescue_enrichment_posterior_threshold);
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
    if let Some(path) = leakage_audit_path
        && let Ok(mut f) = std::fs::File::create(&path)
    {
        use std::io::Write;
        writeln!(
            f,
            "sample_name\ttarget_idx\ttarget_name\treason\tcount_before\tfraction\tdominant_target_idx\tdominant_target_name\tlocal_pair_count\tlocal_pair_fraction"
        )
        .ok();
        for row in &leakage_audit_rows {
            writeln!(f, "{}", row).ok();
        }
        info!("Wrote leakage audit to {}", path);
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
    fn phase2_masks_use_experiment_wide_condition_rescue() {
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
        assert_eq!(masks[1], vec![true, true, false]);
    }

    #[test]
    fn phase2_masks_remain_condition_specific_for_condition_aware_consensus() {
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
            true,
            false,
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
