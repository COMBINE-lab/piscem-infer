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
pub fn update_condition_means(
    results: &[SampleResult],
    hyperparams: &mut DirichletHyperparams,
) {
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

    let weighted_sum: f64 = pi
        .iter()
        .zip(gamma.iter())
        .map(|(&p, &g)| g * p)
        .sum();

    if weighted_sum <= 0.0 {
        return vec![0.0; pi.len()];
    }

    pi.iter()
        .zip(gamma.iter())
        .map(|(&p, &g)| alpha_0 * g * p / weighted_sum)
        .collect()
}

/// Compute convergence metric: max absolute change in pi across all conditions.
pub fn convergence_metric(
    old_pi: &[Vec<f64>],
    new: &DirichletHyperparams,
) -> f64 {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn make_results(
        counts: Vec<Vec<f64>>,
        conditions: Vec<usize>,
    ) -> Vec<SampleResult> {
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
                vec![80.0, 20.0],  // ctrl sample 1
                vec![60.0, 40.0],  // ctrl sample 2
                vec![20.0, 80.0],  // treat sample 1
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
        let gamma = compute_inclusion_probabilities(
            &counts, &cond_indices, 2, 4, 1e-8, 0.02,
        );

        // Target 0: 2/2 in cond 0 → γ = (2 + 0.02) / (2 + 1) = 0.673
        assert!(gamma[0] > 0.5, "target 0 (2/2 present) should have γ > 0.5, got {}", gamma[0]);
        // Target 1: 1/2 in cond 0 → γ = (1 + 0.02) / (2 + 1) = 0.34
        assert!(gamma[1] < 0.5, "target 1 (1/2 present) should have γ < 0.5, got {}", gamma[1]);
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
        let gamma = compute_inclusion_probabilities(
            &counts, &cond_indices, 2, 4, 1e-8, 0.5,
        );

        // With symmetric prior, 1/2 present → γ = (1 + 0.5) / (2 + 1) = 0.5
        // So targets with 1/2 present are borderline (exactly 0.5)
        assert!(gamma[1] >= 0.49, "with dense prior, 1/2 present should be ~0.5, got {}", gamma[1]);
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
}
