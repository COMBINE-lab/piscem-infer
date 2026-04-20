/// Hierarchical Dirichlet-Multinomial model for multi-sample quantification.
///
/// Model:
///   π_c           = condition-level mean proportions (simplex)
///   α₀            = concentration parameter (controls shrinkage strength)
///   θ_s ~ Dir(α₀ · π_c)   for sample s in condition c
///   x_s | θ_s ~ Multinomial(N_s, θ_s)   (via EQ class likelihood)
///
/// MAP-EM pseudo-counts: α_t = α₀ · π_c,t for present transcripts.

/// Per-sample result from (penalized) EM.
#[allow(dead_code)]
pub struct SampleResult {
    /// Estimated counts per transcript from EM
    pub counts: Vec<f64>,
    /// Per-transcript presence mask (true if EM count > presence threshold)
    pub present: Vec<bool>,
    /// Sample name
    pub sample_name: String,
    /// Index into condition list
    pub condition_idx: usize,
}

/// Hyperparameters for the Dirichlet-Multinomial hierarchical model.
pub struct DirichletHyperparams {
    /// Condition-level mean proportions: `pi[condition_idx][transcript_idx]`
    /// Each pi[c] sums to 1 (or 0 if no transcripts are present).
    pub pi: Vec<Vec<f64>>,
    /// Concentration parameter: total pseudo-count strength.
    /// Pseudo-counts for sample s are α_t = alpha_0 * pi[c][t].
    pub alpha_0: f64,
    /// Condition names (for indexing into `pi`)
    pub condition_names: Vec<String>,
}

/// Initialize hyperparameters with a uniform prior.
///
/// - `pi` initialized to uniform over all transcripts (1/T)
/// - `alpha_0` set to `prior_weight * avg_total_reads` (will be updated after init iteration)
pub fn init_hyperparams(
    num_targets: usize,
    condition_names: Vec<String>,
    alpha_0: f64,
) -> DirichletHyperparams {
    let num_conditions = condition_names.len();
    let uniform = 1.0 / num_targets as f64;
    DirichletHyperparams {
        pi: vec![vec![uniform; num_targets]; num_conditions],
        alpha_0,
        condition_names,
    }
}

/// Update condition-level mean proportions from per-sample EM results.
///
/// For each condition c:
///   π̂_c,t = mean(θ̂_s,t for s in condition c)
///   then normalize so Σ_t π̂_c,t = 1
///
/// Only present transcripts contribute to the mean.
pub fn update_condition_means(results: &[SampleResult], hyperparams: &mut DirichletHyperparams) {
    let num_targets = hyperparams.pi[0].len();
    let num_conditions = hyperparams.condition_names.len();

    for c in 0..num_conditions {
        let mut mean_counts = vec![0.0_f64; num_targets];
        let mut n_samples = 0u32;

        for res in results.iter().filter(|r| r.condition_idx == c) {
            n_samples += 1;
            for t in 0..num_targets {
                if res.present[t] {
                    mean_counts[t] += res.counts[t];
                }
            }
        }

        if n_samples == 0 {
            continue;
        }

        // Normalize to proportions
        let total: f64 = mean_counts.iter().sum();
        if total > 0.0 {
            for t in 0..num_targets {
                hyperparams.pi[c][t] = mean_counts[t] / total;
            }
        }
    }
}

/// Compute pseudo-count vector for a sample given the current hyperparameters.
///
/// α_t = α₀ · π_c,t   for present transcripts
/// α_t = 0             for absent transcripts
pub fn compute_pseudo_counts(
    hyperparams: &DirichletHyperparams,
    condition_idx: usize,
    present: &[bool],
) -> Vec<f64> {
    let pi = &hyperparams.pi[condition_idx];
    let alpha_0 = hyperparams.alpha_0;

    // Only distribute pseudo-counts among present transcripts,
    // renormalizing pi over the support set.
    let support_sum: f64 = pi
        .iter()
        .zip(present.iter())
        .filter(|&(_, &p)| p)
        .map(|(&v, _)| v)
        .sum();

    if support_sum <= 0.0 {
        return vec![0.0; pi.len()];
    }

    pi.iter()
        .zip(present.iter())
        .map(|(&p, &is_present)| {
            if is_present {
                alpha_0 * p / support_sum
            } else {
                0.0
            }
        })
        .collect()
}

/// Estimate the fraction of expressed transcripts from initial EM counts.
///
/// A transcript is considered expressed if it exceeds `presence_thresh`
/// in at least one sample. Returns the fraction of such transcripts.
pub fn estimate_inclusion_rate(
    counts: &[Vec<f64>],
    num_targets: usize,
    presence_thresh: f64,
) -> f64 {
    let mut any_present = vec![false; num_targets];
    for sample_counts in counts {
        for (t, &c) in sample_counts.iter().enumerate() {
            if c > presence_thresh {
                any_present[t] = true;
            }
        }
    }
    let n_present = any_present.iter().filter(|&&b| b).count();
    // Clamp to [0.001, 0.999] to avoid degenerate Beta priors
    (n_present as f64 / num_targets as f64).clamp(0.001, 0.999)
}

/// Compute spike-and-slab inclusion probabilities per transcript.
///
/// For each transcript t, condition c:
///   γ_t^(c) = Beta posterior mean = (n_present + a) / (n_reps + a + b)
/// where n_present = number of replicates in c with count > threshold.
///
/// Global inclusion: γ_t = max_c(γ_t^(c))
///
/// Uses an empirical Bayes prior: Beta(π₀, 1-π₀) where π₀ is the
/// estimated fraction of expressed transcripts. This encodes sparsity:
/// with π₀ = 0.02 (sparse transcriptome), the prior strongly favors
/// exclusion and transcripts need replicate consistency to pass.
pub fn compute_inclusion_probabilities(
    counts: &[Vec<f64>],
    condition_indices: &[usize],
    num_conditions: usize,
    num_targets: usize,
    presence_thresh: f64,
    prior_inclusion_rate: f64,
) -> Vec<f64> {
    let beta_a = prior_inclusion_rate;
    let beta_b = 1.0 - prior_inclusion_rate;

    // Count replicates per condition
    let mut cond_n_reps = vec![0u32; num_conditions];
    for &c in condition_indices {
        cond_n_reps[c] += 1;
    }

    // Count present replicates per transcript per condition
    let mut cond_present: Vec<Vec<u32>> = vec![vec![0u32; num_targets]; num_conditions];
    for (i, &cond_idx) in condition_indices.iter().enumerate() {
        for t in 0..num_targets {
            if counts[i][t] > presence_thresh {
                cond_present[cond_idx][t] += 1;
            }
        }
    }

    // Compute gamma as max over conditions of Beta posterior mean
    let mut gamma = vec![0.0f64; num_targets];
    for c in 0..num_conditions {
        let n_reps = cond_n_reps[c] as f64;
        if n_reps == 0.0 {
            continue;
        }
        for t in 0..num_targets {
            let n_present = cond_present[c][t] as f64;
            let gamma_ct = (n_present + beta_a) / (n_reps + beta_a + beta_b);
            gamma[t] = gamma[t].max(gamma_ct);
        }
    }

    gamma
}

/// Compute pseudo-counts weighted by spike-and-slab inclusion probability.
///
/// α_t = γ_t × α₀ × π_c,t / Σ_{t': γ_t' > 0} (γ_t' × π_c,t')
pub fn compute_pseudo_counts_spike_slab(
    hyperparams: &DirichletHyperparams,
    condition_idx: usize,
    gamma: &[f64],
) -> Vec<f64> {
    let pi = &hyperparams.pi[condition_idx];
    let alpha_0 = hyperparams.alpha_0;

    let weighted_sum: f64 = pi.iter().zip(gamma.iter()).map(|(&p, &g)| g * p).sum();

    if weighted_sum <= 0.0 {
        return vec![0.0; pi.len()];
    }

    pi.iter()
        .zip(gamma.iter())
        .map(|(&p, &g)| alpha_0 * g * p / weighted_sum)
        .collect()
}

/// Compute convergence metric: max absolute change in pi across all conditions.
pub fn convergence_metric(old_pi: &[Vec<f64>], new: &DirichletHyperparams) -> f64 {
    let mut max_change = 0.0_f64;
    for (c, old_pi_c) in old_pi.iter().enumerate() {
        for (t, &old_val) in old_pi_c.iter().enumerate() {
            let new_val = new.pi[c][t];
            let change = (new_val - old_val).abs();
            max_change = max_change.max(change);
        }
    }
    max_change
}

// ---------------------------------------------------------------------------
// Special functions: digamma, trigamma, inverse trigamma
// ---------------------------------------------------------------------------

/// Digamma function ψ(x) = d/dx ln Γ(x).
///
/// Uses recurrence ψ(x) = ψ(x+1) - 1/x for x < 6, then asymptotic expansion.
fn digamma(mut x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    let mut result = 0.0;
    // Recurrence: shift x up to ≥ 6
    while x < 6.0 {
        result -= 1.0 / x;
        x += 1.0;
    }
    // Asymptotic (Stirling) series for x ≥ 6
    let inv_x = 1.0 / x;
    let inv_x2 = inv_x * inv_x;
    result +=
        x.ln() - 0.5 * inv_x - inv_x2 * (1.0 / 12.0 - inv_x2 * (1.0 / 120.0 - inv_x2 / 252.0));
    result
}

/// Trigamma function ψ₁(x) = d²/dx² ln Γ(x).
///
/// Uses recurrence ψ₁(x) = ψ₁(x+1) + 1/x² for x < 6, then asymptotic expansion.
fn trigamma(mut x: f64) -> f64 {
    if x <= 0.0 {
        return f64::NAN;
    }
    let mut result = 0.0;
    while x < 6.0 {
        result += 1.0 / (x * x);
        x += 1.0;
    }
    // Asymptotic series
    let inv_x = 1.0 / x;
    let inv_x2 = inv_x * inv_x;
    result +=
        inv_x + inv_x2 * 0.5 + inv_x2 * inv_x * (1.0 / 6.0 - inv_x2 * (1.0 / 30.0 - inv_x2 / 42.0));
    result
}

/// Inverse trigamma: find x > 0 such that trigamma(x) = y.
///
/// Uses Newton's method. For the derivative, d/dx trigamma(x) = -tetragamma(x),
/// approximated via finite differences.
fn inv_trigamma(y: f64) -> f64 {
    if y <= 0.0 || !y.is_finite() {
        return f64::INFINITY;
    }
    // Initial guess
    let mut x = if y > 1.0 {
        1.0 / y.sqrt()
    } else {
        1.0 / (y - 1e-10)
    };
    x = x.max(0.5);

    for _ in 0..25 {
        let tx = trigamma(x);
        let residual = tx - y;
        if residual.abs() < 1e-12 {
            break;
        }
        // Numerical derivative: d(trigamma)/dx ≈ (trigamma(x+h) - trigamma(x)) / h
        let h = x * 1e-6;
        let dtx = (trigamma(x + h) - tx) / h;
        if dtx.abs() < 1e-30 {
            break;
        }
        x -= residual / dtx;
        x = x.max(1e-10);
    }
    x
}

// ---------------------------------------------------------------------------
// Moderated variance estimation (limma squeezeVar)
// ---------------------------------------------------------------------------

/// Per-transcript moderated variance estimates and adaptive concentration parameters.
#[allow(dead_code)]
pub struct ModeratedVariances {
    /// Per-transcript sample variance of log-counts
    pub sample_var: Vec<f64>,
    /// Per-transcript residual degrees of freedom (n_present - 1)
    pub df: Vec<f64>,
    /// Prior degrees of freedom (empirical Bayes)
    pub d0: f64,
    /// Prior variance (empirical Bayes)
    pub s0_sq: f64,
    /// Per-transcript moderated (posterior) variance
    pub moderated_var: Vec<f64>,
    /// Per-transcript concentration parameter
    pub alpha_0_t: Vec<f64>,
}

/// Compute per-transcript sample variance of log(count + pseudo) across samples.
///
/// Returns `(sample_var, df, mean_log_expr)`. Transcripts present in fewer than
/// 2 samples get `sample_var = NaN` and `df = 0`.
pub fn compute_log_count_variances(
    counts: &[Vec<f64>],
    presence: &[Vec<bool>],
    num_targets: usize,
    log_pseudo: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let num_samples = counts.len();
    let mut sample_var = vec![f64::NAN; num_targets];
    let mut df = vec![0.0; num_targets];
    let mut mean_log = vec![0.0; num_targets];

    for t in 0..num_targets {
        // Collect log-counts for present samples
        let mut vals = Vec::new();
        for s in 0..num_samples {
            if presence[s][t] {
                vals.push((counts[s][t] + log_pseudo).ln());
            }
        }
        let n = vals.len();
        if n < 2 {
            continue;
        }
        let mean: f64 = vals.iter().sum::<f64>() / n as f64;
        let var: f64 = vals.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / (n - 1) as f64;
        sample_var[t] = var;
        df[t] = (n - 1) as f64;
        mean_log[t] = mean;
    }

    (sample_var, df, mean_log)
}

/// Fit scaled inverse chi-squared prior (d₀, s₀²) from per-transcript sample variances.
///
/// Uses method of moments on log(s²_t), following limma's `squeezeVar`:
/// - `Var(log s²) = trigamma(d/2) + trigamma(d₀/2)` → solve for d₀
/// - `E(log s²) = log(s₀²) + ψ(d/2) - log(d/2) + ψ(d₀/2) - log(d₀/2)` → solve for s₀²
///
/// Falls back to (∞, median(s²)) if fitting fails.
pub fn fit_variance_prior(sample_var: &[f64], df: &[f64]) -> (f64, f64) {
    // Collect valid (finite, positive) log-variances
    let valid: Vec<(f64, f64)> = sample_var
        .iter()
        .zip(df.iter())
        .filter(|(s2, d)| **d >= 1.0 && s2.is_finite() && **s2 > 0.0)
        .map(|(s2, d)| (s2.ln(), *d))
        .collect();

    if valid.len() < 3 {
        // Not enough data to fit prior; use median variance as fallback
        let mut finite_vars: Vec<f64> = sample_var
            .iter()
            .filter(|v| v.is_finite() && **v > 0.0)
            .copied()
            .collect();
        finite_vars.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median = if finite_vars.is_empty() {
            1.0
        } else {
            finite_vars[finite_vars.len() / 2]
        };
        return (f64::INFINITY, median);
    }

    let n = valid.len() as f64;
    let mean_log_s2: f64 = valid.iter().map(|(z, _)| z).sum::<f64>() / n;
    let var_log_s2: f64 = valid
        .iter()
        .map(|(z, _)| (z - mean_log_s2).powi(2))
        .sum::<f64>()
        / (n - 1.0);

    // Most transcripts share the same df; use the median df
    let mut dfs: Vec<f64> = valid.iter().map(|(_, d)| *d).collect();
    dfs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median_df = dfs[dfs.len() / 2];
    let half_d = median_df / 2.0;

    // trigamma(d₀/2) = Var(log s²) - trigamma(d/2)
    let tg_residual = var_log_s2 - trigamma(half_d);
    if tg_residual <= 0.0 {
        // Prior is uninformative (d₀ → ∞)
        let mut finite_vars: Vec<f64> = sample_var
            .iter()
            .filter(|v| v.is_finite() && **v > 0.0)
            .copied()
            .collect();
        finite_vars.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median = finite_vars[finite_vars.len() / 2];
        return (f64::INFINITY, median);
    }

    let half_d0 = inv_trigamma(tg_residual);
    let d0 = 2.0 * half_d0;

    // log(s₀²) = E(log s²) - ψ(d/2) + log(d/2) - ψ(d₀/2) + log(d₀/2)
    let log_s0_sq = mean_log_s2 - digamma(half_d) + half_d.ln() - digamma(half_d0) + half_d0.ln();
    let s0_sq = log_s0_sq.exp();

    (d0, s0_sq)
}

/// Compute moderated variances and per-transcript concentration parameters.
///
/// σ̃²_t = (d₀·s₀² + d_t·s²_t) / (d₀ + d_t)
/// α₀_t = base_α₀ · s₀² / σ̃²_t, clamped to [min_alpha, max_mult × base_α₀]
pub fn compute_moderated_variances(
    sample_var: &[f64],
    df: &[f64],
    d0: f64,
    s0_sq: f64,
    base_alpha_0: f64,
    min_alpha: f64,
    max_mult: f64,
) -> ModeratedVariances {
    let num_targets = sample_var.len();
    let mut moderated_var = vec![s0_sq; num_targets];
    let mut alpha_0_t = vec![base_alpha_0; num_targets];

    for t in 0..num_targets {
        if df[t] < 1.0 || !sample_var[t].is_finite() {
            // Insufficient data: use prior as-is
            continue;
        }

        let s2 = sample_var[t];
        let dt = df[t];

        let sigma2_tilde = if d0.is_infinite() {
            s0_sq
        } else {
            (d0 * s0_sq + dt * s2) / (d0 + dt)
        };

        moderated_var[t] = sigma2_tilde;

        // Scale concentration inversely with moderated variance.
        // One-sided: only reduce shrinkage for high-variance transcripts,
        // never increase it for stable ones (min(1.0, ...)).
        // Squared ratio for stronger modulation of high-variance transcripts.
        let ratio = if sigma2_tilde > 1e-30 {
            (s0_sq / sigma2_tilde).min(1.0)
        } else {
            1.0
        };
        let ratio_sq = ratio * ratio;
        alpha_0_t[t] = (base_alpha_0 * ratio_sq).clamp(min_alpha, max_mult * base_alpha_0);
    }

    ModeratedVariances {
        sample_var: sample_var.to_vec(),
        df: df.to_vec(),
        d0,
        s0_sq,
        moderated_var,
        alpha_0_t,
    }
}

/// Compute pseudo-counts with per-transcript adaptive concentration.
///
/// α_t = α₀_t · π_c,t / support_sum
pub fn compute_pseudo_counts_adaptive(
    hyperparams: &DirichletHyperparams,
    condition_idx: usize,
    present: &[bool],
    alpha_0_t: &[f64],
) -> Vec<f64> {
    let pi = &hyperparams.pi[condition_idx];

    let support_sum: f64 = pi
        .iter()
        .zip(present.iter())
        .filter(|&(_, &p)| p)
        .map(|(&v, _)| v)
        .sum();

    if support_sum <= 0.0 {
        return vec![0.0; pi.len()];
    }

    pi.iter()
        .zip(present.iter())
        .zip(alpha_0_t.iter())
        .map(|((&p, &is_present), &a0)| {
            if is_present {
                a0 * p / support_sum
            } else {
                0.0
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_results(counts: Vec<Vec<f64>>, conditions: Vec<usize>) -> Vec<SampleResult> {
        counts
            .into_iter()
            .zip(conditions)
            .enumerate()
            .map(|(i, (c, cond))| {
                let present = c.iter().map(|&v| v > 0.0).collect();
                SampleResult {
                    counts: c,
                    present,
                    sample_name: format!("sample_{}", i),
                    condition_idx: cond,
                }
            })
            .collect()
    }

    #[test]
    fn test_init_hyperparams() {
        let hp = init_hyperparams(100, vec!["ctrl".into(), "treat".into()], 50.0);
        assert_eq!(hp.pi.len(), 2);
        assert_eq!(hp.pi[0].len(), 100);
        assert!((hp.pi[0].iter().sum::<f64>() - 1.0).abs() < 1e-10);
        assert_eq!(hp.alpha_0, 50.0);
    }

    #[test]
    fn test_update_condition_means_single_condition() {
        let num_targets = 3;
        let mut hp = init_hyperparams(num_targets, vec!["ctrl".into()], 10.0);

        let results = make_results(
            vec![vec![10.0, 20.0, 30.0], vec![30.0, 40.0, 50.0]],
            vec![0, 0],
        );

        update_condition_means(&results, &mut hp);

        // Mean counts: [20, 30, 40] → proportions: [2/9, 3/9, 4/9]
        let total = 20.0 + 30.0 + 40.0;
        assert!((hp.pi[0][0] - 20.0 / total).abs() < 1e-10);
        assert!((hp.pi[0][1] - 30.0 / total).abs() < 1e-10);
        assert!((hp.pi[0][2] - 40.0 / total).abs() < 1e-10);
    }

    #[test]
    fn test_update_condition_means_two_conditions() {
        let num_targets = 2;
        let mut hp = init_hyperparams(num_targets, vec!["ctrl".into(), "treat".into()], 10.0);

        let results = make_results(
            vec![
                vec![80.0, 20.0], // ctrl sample 1
                vec![60.0, 40.0], // ctrl sample 2
                vec![20.0, 80.0], // treat sample 1
            ],
            vec![0, 0, 1],
        );

        update_condition_means(&results, &mut hp);

        // ctrl: [140, 60] → [0.7, 0.3]
        assert!((hp.pi[0][0] - 0.7).abs() < 1e-10);
        assert!((hp.pi[0][1] - 0.3).abs() < 1e-10);
        // treat: [20, 80] → [0.2, 0.8]
        assert!((hp.pi[1][0] - 0.2).abs() < 1e-10);
        assert!((hp.pi[1][1] - 0.8).abs() < 1e-10);
    }

    #[test]
    fn test_compute_pseudo_counts() {
        let hp = DirichletHyperparams {
            pi: vec![vec![0.6, 0.3, 0.1]],
            alpha_0: 100.0,
            condition_names: vec!["ctrl".into()],
        };

        let present = vec![true, true, true];
        let alpha = compute_pseudo_counts(&hp, 0, &present);
        assert!((alpha[0] - 60.0).abs() < 1e-10);
        assert!((alpha[1] - 30.0).abs() < 1e-10);
        assert!((alpha[2] - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_pseudo_counts_with_absent() {
        let hp = DirichletHyperparams {
            pi: vec![vec![0.6, 0.3, 0.1]],
            alpha_0: 100.0,
            condition_names: vec!["ctrl".into()],
        };

        // Transcript 2 is absent — pseudo-counts should be renormalized over present
        let present = vec![true, true, false];
        let alpha = compute_pseudo_counts(&hp, 0, &present);
        assert_eq!(alpha[2], 0.0);
        // 0.6/(0.6+0.3) = 2/3, 0.3/(0.6+0.3) = 1/3
        assert!((alpha[0] - 100.0 * 2.0 / 3.0).abs() < 1e-10);
        assert!((alpha[1] - 100.0 * 1.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_estimate_inclusion_rate() {
        // 5 targets, 2 samples: targets 0,1,2 present in at least one sample
        let counts = vec![
            vec![100.0, 50.0, 0.0, 0.0, 0.0],
            vec![0.0, 30.0, 20.0, 0.0, 0.0],
        ];
        let rate = estimate_inclusion_rate(&counts, 5, 1e-8);
        assert!((rate - 0.6).abs() < 1e-10); // 3/5

        // All present
        let counts_all = vec![vec![1.0, 1.0, 1.0]];
        let rate_all = estimate_inclusion_rate(&counts_all, 3, 1e-8);
        assert!((rate_all - 0.999).abs() < 1e-10); // clamped to 0.999

        // None present
        let counts_none = vec![vec![0.0, 0.0, 0.0]];
        let rate_none = estimate_inclusion_rate(&counts_none, 3, 1e-8);
        assert!((rate_none - 0.001).abs() < 1e-10); // clamped to 0.001
    }

    #[test]
    fn test_inclusion_probabilities_sparse_prior() {
        // 4 targets, 2 conditions (2 reps each), sparse transcriptome (π₀ = 0.02)
        let counts = vec![
            // condition 0, rep 0: target 0,1 present
            vec![100.0, 50.0, 0.0, 0.0],
            // condition 0, rep 1: target 0 present
            vec![80.0, 0.0, 0.0, 0.0],
            // condition 1, rep 0: target 2,3 present
            vec![0.0, 0.0, 60.0, 40.0],
            // condition 1, rep 1: target 2 present
            vec![0.0, 0.0, 70.0, 0.0],
        ];
        let cond_indices = vec![0, 0, 1, 1];
        let gamma = compute_inclusion_probabilities(&counts, &cond_indices, 2, 4, 1e-8, 0.02);

        // Target 0: 2/2 in cond 0 → γ = (2 + 0.02) / (2 + 1) = 0.673
        assert!(
            gamma[0] > 0.5,
            "target 0 (2/2 present) should have γ > 0.5, got {}",
            gamma[0]
        );
        // Target 1: 1/2 in cond 0 → γ = (1 + 0.02) / (2 + 1) = 0.34
        assert!(
            gamma[1] < 0.5,
            "target 1 (1/2 present) should have γ < 0.5, got {}",
            gamma[1]
        );
        // Target 2: 2/2 in cond 1 → γ > 0.5
        assert!(gamma[2] > 0.5);
        // Target 3: 1/2 in cond 1 → γ < 0.5
        assert!(gamma[3] < 0.5);
    }

    #[test]
    fn test_inclusion_probabilities_dense_prior() {
        // Same data but with dense prior (π₀ = 0.5, like Jeffreys)
        let counts = vec![
            vec![100.0, 50.0, 0.0, 0.0],
            vec![80.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 60.0, 40.0],
            vec![0.0, 0.0, 70.0, 0.0],
        ];
        let cond_indices = vec![0, 0, 1, 1];
        let gamma = compute_inclusion_probabilities(&counts, &cond_indices, 2, 4, 1e-8, 0.5);

        // With symmetric prior, 1/2 present → γ = (1 + 0.5) / (2 + 1) = 0.5
        // So targets with 1/2 present are borderline (exactly 0.5)
        assert!(
            gamma[1] >= 0.49,
            "with dense prior, 1/2 present should be ~0.5, got {}",
            gamma[1]
        );
    }

    #[test]
    fn test_convergence_metric() {
        let old_pi = vec![vec![0.5, 0.3, 0.2]];
        let new = DirichletHyperparams {
            pi: vec![vec![0.6, 0.25, 0.15]],
            alpha_0: 10.0,
            condition_names: vec!["ctrl".into()],
        };
        let change = convergence_metric(&old_pi, &new);
        assert!((change - 0.1).abs() < 1e-10);
    }

    // --- Moderated variance tests ---

    #[test]
    fn test_digamma_trigamma_known_values() {
        // ψ(1) = -γ (Euler-Mascheroni constant)
        let euler = 0.5772156649015329;
        assert!(
            (digamma(1.0) + euler).abs() < 1e-8,
            "ψ(1) = {}, expected {}",
            digamma(1.0),
            -euler
        );

        // ψ'(1) = π²/6
        let pi_sq_6 = std::f64::consts::PI.powi(2) / 6.0;
        assert!(
            (trigamma(1.0) - pi_sq_6).abs() < 1e-8,
            "ψ'(1) = {}, expected {}",
            trigamma(1.0),
            pi_sq_6
        );

        // ψ(5) = 1 + 1/2 + 1/3 + 1/4 - γ = 25/12 - γ
        let psi5 = 1.0 + 0.5 + 1.0 / 3.0 + 0.25 - euler;
        assert!(
            (digamma(5.0) - psi5).abs() < 1e-8,
            "ψ(5) = {}, expected {}",
            digamma(5.0),
            psi5
        );
    }

    #[test]
    fn test_inv_trigamma_roundtrip() {
        for &x in &[0.5, 1.0, 2.0, 5.0, 10.0, 50.0] {
            let y = trigamma(x);
            let x_recovered = inv_trigamma(y);
            assert!(
                (x_recovered - x).abs() < 1e-6,
                "inv_trigamma(trigamma({})) = {}, expected {}",
                x,
                x_recovered,
                x
            );
        }
    }

    #[test]
    fn test_log_count_variances_basic() {
        // 3 transcripts, 4 samples
        let counts = vec![
            vec![100.0, 10.0, 50.0],
            vec![110.0, 12.0, 50.0],
            vec![90.0, 8.0, 50.0],
            vec![105.0, 11.0, 50.0],
        ];
        let presence = vec![
            vec![true, true, true],
            vec![true, true, true],
            vec![true, true, true],
            vec![true, true, true],
        ];

        let (s2, df, mean_log) = compute_log_count_variances(&counts, &presence, 3, 1.0);

        // All present in 4 samples → df = 3
        assert_eq!(df[0], 3.0);
        assert_eq!(df[1], 3.0);
        assert_eq!(df[2], 3.0);

        // Transcript 2 (constant at 50) should have very low variance
        assert!(
            s2[2] < 1e-10,
            "constant transcript should have ~0 variance, got {}",
            s2[2]
        );
        // Transcript 0 and 1 should have positive variance
        assert!(s2[0] > 0.0);
        assert!(s2[1] > 0.0);
        // Mean log should be reasonable
        assert!(mean_log[0] > 4.0); // ln(101) ≈ 4.62
    }

    #[test]
    fn test_log_count_variances_single_sample() {
        let counts = vec![vec![100.0, 50.0]];
        let presence = vec![vec![true, true]];

        let (s2, df, _) = compute_log_count_variances(&counts, &presence, 2, 1.0);

        // Only 1 sample → df = 0, s2 = NaN
        assert_eq!(df[0], 0.0);
        assert_eq!(df[1], 0.0);
        assert!(s2[0].is_nan());
        assert!(s2[1].is_nan());
    }

    #[test]
    fn test_fit_variance_prior_fallback() {
        // Only 2 valid transcripts — not enough for fitting
        let s2 = vec![0.1, 0.2, f64::NAN];
        let df = vec![3.0, 3.0, 0.0];

        let (d0, s0_sq) = fit_variance_prior(&s2, &df);
        // Should fallback: d0 = INF, s0_sq = median of valid s2
        assert!(d0.is_infinite());
        assert!(s0_sq > 0.0);
    }

    #[test]
    fn test_moderated_variance_shrinkage_direction() {
        let s0_sq = 1.0;
        let d0 = 4.0;
        let base_alpha = 100.0;

        // Transcript with high variance (s² = 10) and low variance (s² = 0.1)
        let sample_var = vec![10.0, 0.1, 1.0]; // high, low, typical
        let df = vec![7.0, 7.0, 7.0];

        let mv = compute_moderated_variances(&sample_var, &df, d0, s0_sq, base_alpha, 1.0, 1.0);

        // High-variance transcript: σ̃² between s0_sq and s²
        assert!(
            mv.moderated_var[0] > s0_sq,
            "high-var σ̃² should exceed s0_sq"
        );
        assert!(
            mv.moderated_var[0] < 10.0,
            "high-var σ̃² should be less than s²"
        );

        // Low-variance transcript: σ̃² between s² and s0_sq
        assert!(
            mv.moderated_var[1] < s0_sq,
            "low-var σ̃² should be less than s0_sq"
        );
        assert!(mv.moderated_var[1] > 0.1, "low-var σ̃² should exceed s²");

        // Typical transcript: σ̃² ≈ s0_sq
        assert!((mv.moderated_var[2] - s0_sq).abs() < 0.5);

        // Concentration: high-var → less shrinkage (lower α₀_t, squared ratio)
        assert!(mv.alpha_0_t[0] < base_alpha);
        // Low-var → NO increase (one-sided: capped at base)
        assert_eq!(
            mv.alpha_0_t[1], base_alpha,
            "stable transcripts should not get increased shrinkage"
        );
        // Typical → unchanged
        assert_eq!(mv.alpha_0_t[2], base_alpha);
    }

    #[test]
    fn test_adaptive_pseudo_counts_de_less_shrinkage() {
        let hp = DirichletHyperparams {
            pi: vec![vec![0.5, 0.5]],
            alpha_0: 100.0,
            condition_names: vec!["all".into()],
        };
        let present = vec![true, true];

        // Transcript 0: DE (high variance → low α₀_t)
        // Transcript 1: stable (low variance → high α₀_t)
        let alpha_0_t = vec![20.0, 200.0];

        let alpha = compute_pseudo_counts_adaptive(&hp, 0, &present, &alpha_0_t);

        // Both have equal π (0.5), so ratio is purely from α₀_t
        assert!(
            alpha[0] < alpha[1],
            "DE transcript should get less shrinkage"
        );
        assert_eq!(alpha[0], 20.0 * 0.5 / 1.0); // 10
        assert_eq!(alpha[1], 200.0 * 0.5 / 1.0); // 100
    }

    #[test]
    fn test_adaptive_matches_global_when_uniform() {
        let hp = DirichletHyperparams {
            pi: vec![vec![0.6, 0.3, 0.1]],
            alpha_0: 100.0,
            condition_names: vec!["ctrl".into()],
        };
        let present = vec![true, true, true];

        // All α₀_t = base_α₀ → should match compute_pseudo_counts exactly
        let alpha_0_t = vec![100.0, 100.0, 100.0];
        let adaptive = compute_pseudo_counts_adaptive(&hp, 0, &present, &alpha_0_t);
        let global = compute_pseudo_counts(&hp, 0, &present);

        for (a, g) in adaptive.iter().zip(global.iter()) {
            assert!((a - g).abs() < 1e-10, "adaptive={} vs global={}", a, g);
        }
    }
}
