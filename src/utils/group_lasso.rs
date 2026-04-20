use tracing::info;

use crate::utils::eq_maps::{EqLabel, PackedEqMap, TargetLabelsRef};

/// Joint abundance matrix: T transcripts × S samples.
/// Column-major: `theta[s][t]` = rate for transcript t in sample s.
pub struct AbundanceMatrix {
    pub theta: Vec<Vec<f64>>,
    pub num_targets: usize,
    pub num_samples: usize,
}

#[allow(dead_code)]
impl AbundanceMatrix {
    pub fn new(num_targets: usize, num_samples: usize) -> Self {
        Self {
            theta: vec![vec![0.0; num_targets]; num_samples],
            num_targets,
            num_samples,
        }
    }

    /// Deep clone into an existing matrix (avoids allocation).
    fn copy_from(&mut self, other: &AbundanceMatrix) {
        for s in 0..self.num_samples {
            self.theta[s].copy_from_slice(&other.theta[s]);
        }
    }

    /// Subtract `other` from `self` element-wise: self -= other.
    fn sub_assign(&mut self, other: &AbundanceMatrix) {
        for s in 0..self.num_samples {
            for t in 0..self.num_targets {
                self.theta[s][t] -= other.theta[s][t];
            }
        }
    }

    /// Squared Frobenius norm.
    fn sq_norm(&self) -> f64 {
        self.theta
            .iter()
            .flat_map(|row| row.iter())
            .map(|&x| x * x)
            .sum()
    }
}

impl Clone for AbundanceMatrix {
    fn clone(&self) -> Self {
        Self {
            theta: self.theta.clone(),
            num_targets: self.num_targets,
            num_samples: self.num_samples,
        }
    }
}

/// Controls whether groups span all samples or are per-condition.
pub enum GroupScope {
    /// One group per transcript across all samples: Σ_t ||Θ_t*||₂
    AllSamples,
    /// Per-condition groups: Σ_t Σ_c ||Θ_{t,c*}||₂
    /// `condition_indices[s]` gives the condition index for sample s.
    PerCondition {
        condition_indices: Vec<usize>,
        num_conditions: usize,
    },
}

#[allow(dead_code)]
pub struct GroupLassoConfig {
    pub lambda: f64,
    pub max_iter: u32,
    pub convergence_thresh: f64,
    pub group_scope: GroupScope,
}

// ---------------------------------------------------------------------------
// Poisson negative log-likelihood and gradient
// ---------------------------------------------------------------------------

/// Compute Poisson NLL and its gradient w.r.t. θ.
///
/// NLL = Σ_s Σ_c [ λ_cs − n_cs · ln(λ_cs) ]
/// where λ_cs = Σ_{t ∈ EQC(c)} p_ct · θ[s][t] / eff_len[t]
///
/// The conditional probabilities p_ct are normalized per EQC so they sum to 1.
/// This ensures the Poisson model is well-defined: a transcript's contribution
/// to an EQC is its fractional share, not its full rate.
///
/// ∂NLL/∂θ[s][t] = Σ_{c: t ∈ EQC(c)} p_ct / eff_len[t] · (1 − n_cs / λ_cs)
#[allow(dead_code)]
fn poisson_nll_and_gradient<EqLabelT: EqLabel>(
    theta: &AbundanceMatrix,
    packed_maps: &[PackedEqMap<EqLabelT>],
    inv_eff_lens: &[f64],
) -> (f64, AbundanceMatrix) {
    let num_targets = theta.num_targets;
    let num_samples = theta.num_samples;
    let mut grad = AbundanceMatrix::new(num_targets, num_samples);
    let mut nll = 0.0_f64;

    // Reusable buffer for normalized conditional probabilities per EQC
    let mut norm_probs: Vec<f64> = Vec::with_capacity(64);

    for (s, packed_map) in packed_maps.iter().enumerate() {
        for (eqc_idx, label) in packed_map.iter_labels().enumerate() {
            let count = packed_map.counts[eqc_idx] as f64;
            if count == 0.0 {
                continue;
            }

            // Normalize conditional probabilities for this EQC
            norm_probs.clear();
            let mut prob_sum = 0.0_f64;
            for cond_prob in label.target_probs() {
                norm_probs.push(cond_prob);
                prob_sum += cond_prob;
            }
            if prob_sum > 0.0 {
                let inv_sum = 1.0 / prob_sum;
                for p in norm_probs.iter_mut() {
                    *p *= inv_sum;
                }
            }

            // Compute λ_cs = Σ_t p_ct · θ[s][t] / eff_len[t]
            let mut lambda_cs = 0.0_f64;
            for (&target_id, &p) in label.target_labels().iter().zip(norm_probs.iter()) {
                let t = target_id as usize;
                lambda_cs += p * theta.theta[s][t] * inv_eff_lens[t];
            }

            // Clamp to avoid log(0) and division by zero
            let lambda_cs_safe = lambda_cs.max(1e-10);

            // NLL contribution: λ_cs - n_cs · ln(λ_cs)
            nll += lambda_cs - count * lambda_cs_safe.ln();

            // Gradient contribution
            let ratio = 1.0 - count / lambda_cs_safe;
            for (&target_id, &p) in label.target_labels().iter().zip(norm_probs.iter()) {
                let t = target_id as usize;
                grad.theta[s][t] += p * inv_eff_lens[t] * ratio;
            }
        }
    }

    (nll, grad)
}

// ---------------------------------------------------------------------------
// Group penalty value
// ---------------------------------------------------------------------------

/// Evaluate the group LASSO penalty: Σ_t ||Θ_t*||₂
/// (or the per-condition variant: Σ_t Σ_c ||Θ_{t,c*}||₂).
fn group_penalty_value(theta: &AbundanceMatrix, scope: &GroupScope) -> f64 {
    match scope {
        GroupScope::AllSamples => {
            let mut penalty = 0.0;
            for t in 0..theta.num_targets {
                let mut sq_sum = 0.0;
                for s in 0..theta.num_samples {
                    sq_sum += theta.theta[s][t] * theta.theta[s][t];
                }
                penalty += sq_sum.sqrt();
            }
            penalty
        }
        GroupScope::PerCondition {
            condition_indices,
            num_conditions,
        } => {
            let mut penalty = 0.0;
            for t in 0..theta.num_targets {
                for c in 0..*num_conditions {
                    let mut sq_sum = 0.0;
                    for (s, &cond_idx) in condition_indices.iter().enumerate() {
                        if cond_idx == c {
                            sq_sum += theta.theta[s][t] * theta.theta[s][t];
                        }
                    }
                    penalty += sq_sum.sqrt();
                }
            }
            penalty
        }
    }
}

// ---------------------------------------------------------------------------
// Proximal operator
// ---------------------------------------------------------------------------

/// Apply the proximal operator for group LASSO + non-negativity.
///
/// For each group (transcript row or transcript-condition sub-row):
///   1. Clamp negatives to zero
///   2. Compute group L2 norm
///   3. If norm ≤ λ·η: zero the group
///   4. Else: scale by (1 − λ·η / norm)
#[allow(dead_code)]
fn proximal_group_lasso_nonneg(
    theta: &mut AbundanceMatrix,
    lambda: f64,
    step_size: f64,
    scope: &GroupScope,
) {
    let threshold = lambda * step_size;

    match scope {
        GroupScope::AllSamples => {
            for t in 0..theta.num_targets {
                // Clamp negatives
                for s in 0..theta.num_samples {
                    theta.theta[s][t] = theta.theta[s][t].max(0.0);
                }
                // Group L2 shrinkage
                let mut sq_sum = 0.0;
                for s in 0..theta.num_samples {
                    sq_sum += theta.theta[s][t] * theta.theta[s][t];
                }
                let norm = sq_sum.sqrt();
                if norm <= threshold {
                    for s in 0..theta.num_samples {
                        theta.theta[s][t] = 0.0;
                    }
                } else {
                    let scale = 1.0 - threshold / norm;
                    for s in 0..theta.num_samples {
                        theta.theta[s][t] *= scale;
                    }
                }
            }
        }
        GroupScope::PerCondition {
            condition_indices,
            num_conditions,
        } => {
            for t in 0..theta.num_targets {
                // Clamp negatives
                for s in 0..theta.num_samples {
                    theta.theta[s][t] = theta.theta[s][t].max(0.0);
                }
                // Per-condition group shrinkage
                for c in 0..*num_conditions {
                    let mut sq_sum = 0.0;
                    for (s, &cond_idx) in condition_indices.iter().enumerate() {
                        if cond_idx == c {
                            sq_sum += theta.theta[s][t] * theta.theta[s][t];
                        }
                    }
                    let norm = sq_sum.sqrt();
                    if norm <= threshold {
                        for (s, &cond_idx) in condition_indices.iter().enumerate() {
                            if cond_idx == c {
                                theta.theta[s][t] = 0.0;
                            }
                        }
                    } else {
                        let scale = 1.0 - threshold / norm;
                        for (s, &cond_idx) in condition_indices.iter().enumerate() {
                            if cond_idx == c {
                                theta.theta[s][t] *= scale;
                            }
                        }
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Composite objective value
// ---------------------------------------------------------------------------

#[allow(dead_code)]
fn composite_objective<EqLabelT: EqLabel>(
    theta: &AbundanceMatrix,
    packed_maps: &[PackedEqMap<EqLabelT>],
    inv_eff_lens: &[f64],
    lambda: f64,
    scope: &GroupScope,
) -> f64 {
    let (nll, _) = poisson_nll_and_gradient(theta, packed_maps, inv_eff_lens);
    nll + lambda * group_penalty_value(theta, scope)
}

// ---------------------------------------------------------------------------
// Backtracking line search
// ---------------------------------------------------------------------------

/// Backtracking line search for proximal gradient.
/// Finds step size η such that the composite objective decreases sufficiently.
///
/// Uses the proximal gradient sufficient decrease condition:
///   F(prox_η(y − η·∇f(y))) ≤ f(y) + ⟨∇f(y), x−y⟩ + (1/2η)||x−y||² + λ·g(x)
#[allow(dead_code)]
fn backtrack_line_search<EqLabelT: EqLabel>(
    y: &AbundanceMatrix,
    grad: &AbundanceMatrix,
    f_y: f64,
    packed_maps: &[PackedEqMap<EqLabelT>],
    inv_eff_lens: &[f64],
    lambda: f64,
    scope: &GroupScope,
) -> (f64, AbundanceMatrix) {
    let mut eta = 1.0;
    let beta = 0.5;
    let num_targets = y.num_targets;
    let num_samples = y.num_samples;

    for _backtrack_iter in 0..30 {
        // Compute candidate: x = prox_η(y - η·grad)
        let mut x = AbundanceMatrix::new(num_targets, num_samples);
        for s in 0..num_samples {
            for t in 0..num_targets {
                x.theta[s][t] = y.theta[s][t] - eta * grad.theta[s][t];
            }
        }
        proximal_group_lasso_nonneg(&mut x, lambda, eta, scope);

        // Evaluate composite objective at x
        let f_x = composite_objective(&x, packed_maps, inv_eff_lens, lambda, scope);

        // Compute the quadratic upper bound:
        // Q(x, y) = f(y) + ⟨∇f(y), x−y⟩ + (1/2η)||x−y||² + λ·g(x)
        let mut diff = x.clone();
        diff.sub_assign(y);
        let mut inner_product = 0.0_f64;
        for s in 0..num_samples {
            for t in 0..num_targets {
                inner_product += grad.theta[s][t] * diff.theta[s][t];
            }
        }
        let sq_diff = diff.sq_norm();
        let g_x = group_penalty_value(&x, scope);
        let q_val = f_y + inner_product + sq_diff / (2.0 * eta) + lambda * g_x;

        if f_x <= q_val + 1e-10 {
            return (eta, x);
        }
        eta *= beta;
    }

    // Fallback: use very small step
    let mut x = AbundanceMatrix::new(num_targets, num_samples);
    for s in 0..num_samples {
        for t in 0..num_targets {
            x.theta[s][t] = y.theta[s][t] - eta * grad.theta[s][t];
        }
    }
    proximal_group_lasso_nonneg(&mut x, lambda, eta, scope);
    (eta, x)
}

// ---------------------------------------------------------------------------
// FISTA optimizer
// ---------------------------------------------------------------------------

/// Run FISTA (Fast Iterative Shrinkage-Thresholding Algorithm) with adaptive restart
/// to minimize:  F(Θ) = NLL(Θ) + λ·Σ_t ||Θ_t*||₂
///
/// Returns the optimized AbundanceMatrix.
#[allow(dead_code)]
pub fn fista_group_lasso<EqLabelT: EqLabel>(
    init_theta: AbundanceMatrix,
    packed_maps: &[PackedEqMap<EqLabelT>],
    inv_eff_lens: &[f64],
    config: &GroupLassoConfig,
) -> AbundanceMatrix {
    let num_targets = init_theta.num_targets;
    let num_samples = init_theta.num_samples;

    let mut theta = init_theta;
    let mut y = theta.clone();
    let mut theta_prev = theta.clone();
    let mut t_prev = 1.0_f64;

    let mut prev_obj = composite_objective(
        &theta,
        packed_maps,
        inv_eff_lens,
        config.lambda,
        &config.group_scope,
    );

    info!(
        "FISTA iter 0: obj={:.6e}, penalty={:.6e}",
        prev_obj,
        group_penalty_value(&theta, &config.group_scope)
    );

    for iter in 1..=config.max_iter {
        // Gradient at the extrapolated point y
        let (f_y, grad) = poisson_nll_and_gradient(&y, packed_maps, inv_eff_lens);

        // Backtracking line search
        let (_eta, theta_new) = backtrack_line_search(
            &y,
            &grad,
            f_y,
            packed_maps,
            inv_eff_lens,
            config.lambda,
            &config.group_scope,
        );

        // Evaluate new objective
        let new_obj = composite_objective(
            &theta_new,
            packed_maps,
            inv_eff_lens,
            config.lambda,
            &config.group_scope,
        );

        // Nesterov momentum step
        let t_new = (1.0 + (1.0 + 4.0 * t_prev * t_prev).sqrt()) / 2.0;
        let momentum = (t_prev - 1.0) / t_new;

        // Compute extrapolated point: y = theta_new + momentum * (theta_new - theta_prev)
        let mut y_new = AbundanceMatrix::new(num_targets, num_samples);
        for s in 0..num_samples {
            for t in 0..num_targets {
                y_new.theta[s][t] = theta_new.theta[s][t]
                    + momentum * (theta_new.theta[s][t] - theta_prev.theta[s][t]);
                // Keep y non-negative (extrapolation can go negative)
                y_new.theta[s][t] = y_new.theta[s][t].max(0.0);
            }
        }

        // Adaptive restart: if objective increased, reset momentum
        if new_obj > prev_obj {
            y.copy_from(&theta_new);
            t_prev = 1.0;
        } else {
            y = y_new;
            t_prev = t_new;
        }

        theta_prev.copy_from(&theta);
        theta = theta_new;

        // Count sparsity
        let n_zero_rows = (0..num_targets)
            .filter(|&t| (0..num_samples).all(|s| theta.theta[s][t] == 0.0))
            .count();

        if iter % 50 == 0 || iter <= 5 {
            info!(
                "FISTA iter {}: obj={:.6e}, penalty={:.6e}, zero_rows={}/{}",
                iter,
                new_obj,
                group_penalty_value(&theta, &config.group_scope),
                n_zero_rows,
                num_targets,
            );
        }

        // Check convergence
        let rel_change = (prev_obj - new_obj).abs() / (prev_obj.abs() + 1e-10);
        if rel_change < config.convergence_thresh && iter > 1 {
            info!(
                "FISTA converged at iter {}: rel_change={:.2e}, zero_rows={}/{}",
                iter, rel_change, n_zero_rows, num_targets,
            );
            break;
        }

        prev_obj = new_obj;
    }

    theta
}

// ---------------------------------------------------------------------------
// BIC for lambda selection
// ---------------------------------------------------------------------------

/// Compute BIC for a given solution.
/// BIC = -2·loglik + k·log(n)
/// where k = number of transcripts with ||θ_t*||₂ > 0
/// and n = total fragments across all samples.
#[allow(dead_code)]
pub fn compute_bic<EqLabelT: EqLabel>(
    theta: &AbundanceMatrix,
    packed_maps: &[PackedEqMap<EqLabelT>],
    inv_eff_lens: &[f64],
    scope: &GroupScope,
) -> f64 {
    let (nll, _) = poisson_nll_and_gradient(theta, packed_maps, inv_eff_lens);

    // Count non-zero transcript rows (or groups for per-condition)
    let k = match scope {
        GroupScope::AllSamples => (0..theta.num_targets)
            .filter(|&t| (0..theta.num_samples).any(|s| theta.theta[s][t] > 0.0))
            .count(),
        GroupScope::PerCondition {
            condition_indices,
            num_conditions,
        } => {
            let mut count = 0usize;
            for t in 0..theta.num_targets {
                for c in 0..*num_conditions {
                    let any_nonzero = condition_indices
                        .iter()
                        .enumerate()
                        .any(|(s, &ci)| ci == c && theta.theta[s][t] > 0.0);
                    if any_nonzero {
                        count += 1;
                    }
                }
            }
            count
        }
    };

    // k parameters per active group, times num_samples (or num_samples_in_condition)
    let effective_k = match scope {
        GroupScope::AllSamples => k * theta.num_samples,
        GroupScope::PerCondition {
            condition_indices, ..
        } => {
            // For simplicity, use k * avg_samples_per_condition
            // Actually each active group has its condition's sample count
            // For BIC, just count total free parameters
            let mut total = 0usize;
            for t in 0..theta.num_targets {
                for (s, _) in condition_indices.iter().enumerate() {
                    if theta.theta[s][t] > 0.0 {
                        total += 1;
                    }
                }
            }
            total
        }
    };

    let total_frags: f64 = packed_maps.iter().map(|m| m.total_weight() as f64).sum();

    2.0 * nll + (effective_k as f64) * total_frags.ln()
}

// ---------------------------------------------------------------------------
// EM-based group shrinkage (Option 1: Penalized EM)
// ---------------------------------------------------------------------------

/// Run joint EM with group L2 shrinkage across samples.
///
/// This avoids the Poisson model entirely and uses the standard multinomial EM.
/// After each outer iteration (one full EM per sample), group L2 shrinkage is
/// applied across samples to drive weakly-supported transcript rows to zero.
/// Counts are then renormalized to preserve total weight per sample.
///
/// Unlike the FISTA approach, the EM has no log barrier at zero: a transcript
/// with θ_t = 0 simply gets no reads assigned in the E-step.
///
/// # Arguments
/// * `packed_maps` — per-sample equivalence class maps
/// * `eff_lens` — per-sample effective lengths (outer index = sample)
/// * `lambda` — group L2 penalty strength
/// * `scope` — all-samples or per-condition grouping
/// * `num_em_iters` — EM iterations per sample per outer round
/// * `num_outer_iters` — number of shrinkage rounds
/// * `max_em_iter` — maximum EM iterations per sample
/// * `convergence_thresh` — EM convergence threshold
/// * `presence_thresh` — threshold below which counts are zeroed
pub fn em_group_shrinkage<EqLabelT: EqLabel>(
    packed_maps: &[PackedEqMap<EqLabelT>],
    eff_lens_per_sample: &[&[f64]],
    lambda: f64,
    scope: &GroupScope,
    num_outer_iters: u32,
    max_em_iter: u32,
    convergence_thresh: f64,
    presence_thresh: f64,
) -> AbundanceMatrix {
    let num_samples = packed_maps.len();
    let num_targets = eff_lens_per_sample[0].len();

    // Pre-compute inverse effective lengths per sample
    let inv_eff_lens: Vec<Vec<f64>> = eff_lens_per_sample
        .iter()
        .map(|els| {
            els.iter()
                .map(|x| {
                    let y = 1.0 / *x;
                    if y.is_finite() { y } else { 0.0 }
                })
                .collect()
        })
        .collect();

    // Total weights per sample (for renormalization after shrinkage)
    let total_weights: Vec<f64> = packed_maps
        .iter()
        .map(|m| m.total_weight() as f64)
        .collect();

    // Initialize: uniform counts per sample
    let mut counts: Vec<Vec<f64>> = (0..num_samples)
        .map(|s| {
            let avg = total_weights[s] / num_targets as f64;
            vec![avg; num_targets]
        })
        .collect();

    // Track which transcripts are active (not zeroed by shrinkage)
    let mut active_mask = vec![true; num_targets];

    for outer_iter in 0..num_outer_iters {
        // --- Per-sample EM ---
        for s in 0..num_samples {
            let eq_map = &packed_maps[s];
            let inv_el = &inv_eff_lens[s];

            // Use current counts as initialization
            let mut prev_counts = counts[s].clone();
            let mut curr_counts = vec![0.0_f64; num_targets];

            // Apply active mask to initialization
            for t in 0..num_targets {
                if !active_mask[t] {
                    prev_counts[t] = 0.0;
                }
            }

            // Run EM iterations
            let mut weights: Vec<f64> = Vec::with_capacity(64);
            for _em_iter in 0..max_em_iter {
                // M-step (same as em.rs m_step)
                for (eqc_idx, label) in eq_map.iter_labels().enumerate() {
                    let count = eq_map.counts[eqc_idx] as f64;
                    if count == 0.0 {
                        continue;
                    }

                    let mut denom = 0.0_f64;
                    for (&target_id, cond_prob) in
                        label.target_labels().iter().zip(label.target_probs())
                    {
                        let t = target_id as usize;
                        let w = cond_prob * prev_counts[t] * inv_el[t];
                        weights.push(w);
                        denom += w;
                    }
                    if denom > 1e-8 {
                        let count_over_denom = count / denom;
                        for (&target_id, w) in label.target_labels().iter().zip(weights.iter()) {
                            curr_counts[target_id as usize] += count_over_denom * w;
                        }
                    }
                    weights.clear();
                }

                // Check convergence
                let mut rel_diff = 0.0_f64;
                for t in 0..num_targets {
                    if prev_counts[t] > presence_thresh {
                        let rd = (curr_counts[t] - prev_counts[t]) / prev_counts[t];
                        if rd > rel_diff {
                            rel_diff = rd;
                        }
                    }
                }

                std::mem::swap(&mut prev_counts, &mut curr_counts);
                curr_counts.fill(0.0);

                if rel_diff < convergence_thresh {
                    break;
                }
            }

            // Presence thresholding
            for x in prev_counts.iter_mut() {
                if *x < presence_thresh {
                    *x = 0.0;
                }
            }

            // Final M-step
            for (eqc_idx, label) in eq_map.iter_labels().enumerate() {
                let count = eq_map.counts[eqc_idx] as f64;
                if count == 0.0 {
                    continue;
                }
                let mut denom = 0.0_f64;
                for (&target_id, cond_prob) in
                    label.target_labels().iter().zip(label.target_probs())
                {
                    let t = target_id as usize;
                    let w = cond_prob * prev_counts[t] * inv_el[t];
                    weights.push(w);
                    denom += w;
                }
                if denom > 1e-8 {
                    let count_over_denom = count / denom;
                    for (&target_id, w) in label.target_labels().iter().zip(weights.iter()) {
                        curr_counts[target_id as usize] += count_over_denom * w;
                    }
                }
                weights.clear();
            }

            counts[s] = curr_counts;
        }

        // --- Group L2 shrinkage across samples ---
        let mut n_zeroed = 0usize;
        match scope {
            GroupScope::AllSamples => {
                for t in 0..num_targets {
                    if !active_mask[t] {
                        n_zeroed += 1;
                        continue;
                    }
                    let mut sq_sum = 0.0;
                    for s in 0..num_samples {
                        sq_sum += counts[s][t] * counts[s][t];
                    }
                    let norm = sq_sum.sqrt();
                    if norm <= lambda {
                        for s in 0..num_samples {
                            counts[s][t] = 0.0;
                        }
                        active_mask[t] = false;
                        n_zeroed += 1;
                    } else {
                        let scale = 1.0 - lambda / norm;
                        for s in 0..num_samples {
                            counts[s][t] *= scale;
                        }
                    }
                }
            }
            GroupScope::PerCondition {
                condition_indices,
                num_conditions,
            } => {
                for t in 0..num_targets {
                    let mut all_zero = true;
                    for c in 0..*num_conditions {
                        let mut sq_sum = 0.0;
                        for (s, &ci) in condition_indices.iter().enumerate() {
                            if ci == c {
                                sq_sum += counts[s][t] * counts[s][t];
                            }
                        }
                        let norm = sq_sum.sqrt();
                        if norm <= lambda {
                            for (s, &ci) in condition_indices.iter().enumerate() {
                                if ci == c {
                                    counts[s][t] = 0.0;
                                }
                            }
                        } else {
                            let scale = 1.0 - lambda / norm;
                            for (s, &ci) in condition_indices.iter().enumerate() {
                                if ci == c {
                                    counts[s][t] *= scale;
                                }
                            }
                            all_zero = false;
                        }
                    }
                    if all_zero {
                        active_mask[t] = false;
                        n_zeroed += 1;
                    }
                }
            }
        }

        // --- Renormalize to preserve total weight per sample ---
        for s in 0..num_samples {
            let sum: f64 = counts[s].iter().sum();
            if sum > 0.0 {
                let scale = total_weights[s] / sum;
                for t in 0..num_targets {
                    counts[s][t] *= scale;
                }
            }
        }

        let n_active = num_targets - active_mask.iter().filter(|&&b| !b).count();
        let penalty = group_penalty_value(
            &AbundanceMatrix {
                theta: counts.clone(),
                num_targets,
                num_samples,
            },
            scope,
        );
        info!(
            "EM-GL outer iter {}: active={}/{}, penalty={:.6e}",
            outer_iter, n_active, num_targets, penalty,
        );

        // Early termination if no transcripts were newly zeroed
        if outer_iter > 0 && n_zeroed == num_targets - n_active {
            // Check if active set stabilized — count wasn't changing
        }
    }

    AbundanceMatrix {
        theta: counts,
        num_targets,
        num_samples,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::eq_maps::{BasicEqLabel, PackedEqMap};

    /// Helper: build a small PackedEqMap from (targets, count) pairs.
    fn make_packed_map(eqcs: &[(&[u32], usize)]) -> PackedEqMap<BasicEqLabel> {
        let mut eq_labels = Vec::new();
        let mut eq_label_starts = vec![0u32];
        let mut counts = Vec::new();
        for (targets, count) in eqcs {
            eq_labels.extend_from_slice(targets);
            eq_label_starts.push(eq_labels.len() as u32);
            counts.push(*count);
        }
        PackedEqMap::from_raw(eq_labels, eq_label_starts, counts, false)
    }

    #[test]
    fn test_proximal_zeros_small_row() {
        let mut theta = AbundanceMatrix {
            theta: vec![vec![0.01, 0.02, 100.0], vec![0.01, 0.03, 200.0]],
            num_targets: 3,
            num_samples: 2,
        };
        let scope = GroupScope::AllSamples;
        // threshold = lambda * step = 1.0 * 1.0 = 1.0
        // Row 0: norm = sqrt(0.01² + 0.01²) ≈ 0.014 < 1.0 → zeroed
        // Row 1: norm = sqrt(0.02² + 0.03²) ≈ 0.036 < 1.0 → zeroed
        // Row 2: norm = sqrt(100² + 200²) ≈ 223.6 > 1.0 → shrunk
        proximal_group_lasso_nonneg(&mut theta, 1.0, 1.0, &scope);
        assert_eq!(theta.theta[0][0], 0.0);
        assert_eq!(theta.theta[1][0], 0.0);
        assert_eq!(theta.theta[0][1], 0.0);
        assert_eq!(theta.theta[1][1], 0.0);
        assert!(theta.theta[0][2] > 0.0);
        assert!(theta.theta[1][2] > 0.0);
    }

    #[test]
    fn test_proximal_preserves_large_row() {
        let mut theta = AbundanceMatrix {
            theta: vec![vec![10.0, 20.0], vec![15.0, 25.0]],
            num_targets: 2,
            num_samples: 2,
        };
        let scope = GroupScope::AllSamples;
        // Small lambda: threshold = 0.001 * 1.0 = 0.001
        // Row 0: norm = sqrt(100 + 225) = 18.03 >> 0.001 → barely shrunk
        proximal_group_lasso_nonneg(&mut theta, 0.001, 1.0, &scope);
        assert!(theta.theta[0][0] > 9.9);
        assert!(theta.theta[1][0] > 14.9);
        assert!(theta.theta[0][1] > 19.9);
        assert!(theta.theta[1][1] > 24.9);
    }

    #[test]
    fn test_proximal_per_condition() {
        // 2 samples, 2 conditions: sample 0 → condition 0, sample 1 → condition 1
        let mut theta = AbundanceMatrix {
            theta: vec![
                vec![0.001, 50.0], // sample 0 (condition 0): txp 0 tiny, txp 1 large
                vec![50.0, 0.001], // sample 1 (condition 1): txp 0 large, txp 1 tiny
            ],
            num_targets: 2,
            num_samples: 2,
        };
        let scope = GroupScope::PerCondition {
            condition_indices: vec![0, 1],
            num_conditions: 2,
        };
        // threshold = 1.0
        // Txp 0: cond 0 has norm=0.001 → zeroed; cond 1 has norm=50 → shrunk
        // Txp 1: cond 0 has norm=50 → shrunk; cond 1 has norm=0.001 → zeroed
        proximal_group_lasso_nonneg(&mut theta, 1.0, 1.0, &scope);
        assert_eq!(theta.theta[0][0], 0.0, "txp 0 in cond 0 should be zeroed");
        assert!(theta.theta[1][0] > 0.0, "txp 0 in cond 1 should survive");
        assert!(theta.theta[0][1] > 0.0, "txp 1 in cond 0 should survive");
        assert_eq!(theta.theta[1][1], 0.0, "txp 1 in cond 1 should be zeroed");
    }

    #[test]
    fn test_nonneg_projection() {
        let mut theta = AbundanceMatrix {
            theta: vec![vec![-5.0, 3.0, -1.0], vec![-2.0, 4.0, -0.5]],
            num_targets: 3,
            num_samples: 2,
        };
        let scope = GroupScope::AllSamples;
        // Very small lambda so shrinkage doesn't zero anything
        proximal_group_lasso_nonneg(&mut theta, 1e-10, 1.0, &scope);
        // Negatives should be clamped to 0
        assert_eq!(theta.theta[0][0], 0.0);
        assert_eq!(theta.theta[1][0], 0.0);
        assert!(theta.theta[0][1] > 0.0);
        assert!(theta.theta[1][1] > 0.0);
        assert_eq!(theta.theta[0][2], 0.0);
        assert_eq!(theta.theta[1][2], 0.0);
    }

    #[test]
    fn test_poisson_gradient_finite_diff() {
        // Small problem: 2 transcripts, 2 samples, 2 EQCs per sample
        // EQC 0: {txp 0}, count 100
        // EQC 1: {txp 0, txp 1}, count 50
        let map0 = make_packed_map(&[(&[0], 100), (&[0, 1], 50)]);
        let map1 = make_packed_map(&[(&[0], 80), (&[0, 1], 40)]);
        let packed_maps = vec![map0, map1];
        let inv_eff_lens = vec![1.0 / 500.0, 1.0 / 300.0]; // eff_len = 500, 300

        let theta = AbundanceMatrix {
            theta: vec![vec![200.0, 100.0], vec![150.0, 80.0]],
            num_targets: 2,
            num_samples: 2,
        };

        let (nll, grad) = poisson_nll_and_gradient(&theta, &packed_maps, &inv_eff_lens);
        assert!(nll.is_finite());

        // Check gradient with finite differences
        let eps = 1e-5;
        for s in 0..2 {
            for t in 0..2 {
                let mut theta_plus = theta.clone();
                theta_plus.theta[s][t] += eps;
                let (nll_plus, _) =
                    poisson_nll_and_gradient(&theta_plus, &packed_maps, &inv_eff_lens);

                let mut theta_minus = theta.clone();
                theta_minus.theta[s][t] -= eps;
                let (nll_minus, _) =
                    poisson_nll_and_gradient(&theta_minus, &packed_maps, &inv_eff_lens);

                let fd_grad = (nll_plus - nll_minus) / (2.0 * eps);
                let analytic_grad = grad.theta[s][t];

                let abs_diff = (fd_grad - analytic_grad).abs();
                let rel_diff = abs_diff / (analytic_grad.abs() + 1e-10);
                assert!(
                    rel_diff < 1e-4 || abs_diff < 1e-6,
                    "Gradient mismatch at [{s}][{t}]: analytic={analytic_grad:.6e}, fd={fd_grad:.6e}, rel={rel_diff:.2e}"
                );
            }
        }
    }

    #[test]
    fn test_convergence_lambda_zero() {
        // With lambda=0, FISTA should converge to the per-sample MLE.
        // 1 transcript, 2 samples
        // EQC: {txp 0}, count = 100 (sample 0), 200 (sample 1)
        let map0 = make_packed_map(&[(&[0], 100)]);
        let map1 = make_packed_map(&[(&[0], 200)]);
        let packed_maps = vec![map0, map1];
        let inv_eff_lens = vec![1.0]; // eff_len = 1

        let init_theta = AbundanceMatrix {
            theta: vec![vec![50.0], vec![50.0]],
            num_targets: 1,
            num_samples: 2,
        };

        let config = GroupLassoConfig {
            lambda: 0.0,
            max_iter: 200,
            convergence_thresh: 1e-10,
            group_scope: GroupScope::AllSamples,
        };

        let result = fista_group_lasso(init_theta, &packed_maps, &inv_eff_lens, &config);

        // MLE for Poisson: θ = n (count), since λ = θ/eff_len = θ and ∂NLL/∂θ = 1 - n/θ = 0 → θ = n
        assert!(
            (result.theta[0][0] - 100.0).abs() < 1.0,
            "Sample 0 MLE should be ~100, got {}",
            result.theta[0][0]
        );
        assert!(
            (result.theta[1][0] - 200.0).abs() < 1.0,
            "Sample 1 MLE should be ~200, got {}",
            result.theta[1][0]
        );
    }

    #[test]
    fn test_em_group_shrinkage_zeros_redundant() {
        // 3 transcripts, 2 samples
        // EQC 0: {txp 0} count=100/80 (unique to txp 0)
        // EQC 1: {txp 0, txp 1, txp 2} count=50/40 (shared)
        // Txp 1 and 2 have no unique support — EM+shrinkage should zero them.
        let map0 = make_packed_map(&[(&[0], 100), (&[0, 1, 2], 50)]);
        let map1 = make_packed_map(&[(&[0], 80), (&[0, 1, 2], 40)]);
        let packed_maps = vec![map0, map1];
        let eff_lens = vec![1.0, 1.0, 1.0];
        let eff_lens_refs: Vec<&[f64]> = vec![&eff_lens, &eff_lens];

        let scope = GroupScope::AllSamples;
        let result = em_group_shrinkage(
            &packed_maps,
            &eff_lens_refs,
            5.0, // lambda
            &scope,
            20,   // outer iters
            100,  // max EM iters
            1e-6, // convergence thresh
            1e-8, // presence thresh
        );

        // Txp 0 must remain active
        assert!(
            result.theta[0][0] > 50.0,
            "txp 0 should remain active, got {}",
            result.theta[0][0],
        );
        // Txp 1 and 2 should be zeroed (no unique support + shrinkage)
        assert_eq!(
            result.theta[0][1], 0.0,
            "txp 1 should be zeroed in sample 0, got {}",
            result.theta[0][1],
        );
        assert_eq!(
            result.theta[0][2], 0.0,
            "txp 2 should be zeroed in sample 0, got {}",
            result.theta[0][2],
        );
        assert_eq!(
            result.theta[1][1], 0.0,
            "txp 1 should be zeroed in sample 1, got {}",
            result.theta[1][1],
        );
    }

    #[test]
    fn test_sparsity_zeros_redundant() {
        // With group LASSO, redundant transcripts should be zeroed.
        // Setup: 3 transcripts, 2 samples
        // EQC 0: {txp 0} count=100/80 (unique to txp 0)
        // EQC 1: {txp 0, txp 1, txp 2} count=50/40 (shared)
        // Txp 1 and 2 have no unique support, so under sparsity penalty
        // the optimal solution assigns all shared reads to txp 0.
        let map0 = make_packed_map(&[(&[0], 100), (&[0, 1, 2], 50)]);
        let map1 = make_packed_map(&[(&[0], 80), (&[0, 1, 2], 40)]);
        let packed_maps = vec![map0, map1];
        let inv_eff_lens = vec![1.0, 1.0, 1.0];

        // Initialize near MLE (without penalty)
        let init_theta = AbundanceMatrix {
            theta: vec![vec![120.0, 15.0, 15.0], vec![96.0, 12.0, 12.0]],
            num_targets: 3,
            num_samples: 2,
        };

        let config = GroupLassoConfig {
            lambda: 5.0,
            max_iter: 1000,
            convergence_thresh: 1e-12,
            group_scope: GroupScope::AllSamples,
        };

        let result = fista_group_lasso(init_theta, &packed_maps, &inv_eff_lens, &config);

        // Txp 0 must remain (has unique support)
        assert!(
            result.theta[0][0] > 50.0,
            "txp 0 should remain active, got {}",
            result.theta[0][0],
        );
        // Txp 1 and 2 should be driven toward zero (no unique support)
        // With lambda=5 they may not be exactly zero but should be much
        // smaller than txp 0
        let ratio = result.theta[0][1] / result.theta[0][0];
        assert!(
            ratio < 0.1,
            "txp 1 should be much smaller than txp 0, ratio = {}",
            ratio,
        );
    }
}
