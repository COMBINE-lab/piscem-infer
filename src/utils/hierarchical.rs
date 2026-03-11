/// Per-sample posterior result from penalized EM + L-BFGS.
/// This is the output of the inner loop (Phase 1) for each sample.
pub struct SamplePosterior {
    /// MAP estimate in log-space (softmax parameterization)
    pub phi_hat: Vec<f64>,
    /// Diagonal Laplace posterior variance per transcript
    pub sigma_hat_sq: Vec<f64>,
    /// Sample name (for identification)
    pub sample_name: String,
    /// Index into `HierarchicalHyperparams::condition_names`
    pub condition_idx: usize,
}

/// Hyperparameters shared across samples in the hierarchical model.
pub struct HierarchicalHyperparams {
    /// Condition-level means per transcript: `nu[condition_idx][transcript_idx]`
    pub nu: Vec<Vec<f64>>,
    /// Biological variance per transcript (shared across conditions)
    pub sigma_sq: Vec<f64>,
    /// Condition names (for indexing into `nu`)
    pub condition_names: Vec<String>,
}

/// Initialize hyperparameters with a vague prior.
///
/// - `nu` initialized to 0 for all conditions and transcripts
/// - `sigma_sq` initialized to 1.0 (vague prior)
pub fn init_hyperparams(
    num_targets: usize,
    condition_names: Vec<String>,
) -> HierarchicalHyperparams {
    let num_conditions = condition_names.len();
    HierarchicalHyperparams {
        nu: vec![vec![0.0; num_targets]; num_conditions],
        sigma_sq: vec![1.0; num_targets],
        condition_names,
    }
}

/// Update condition-level means (precision-weighted).
///
/// For each condition c and transcript t:
///   ν̂_t^(c) = [ Σ_{s: c(s)=c} φ̂_t^(s) / (σ_t² + Σ̂_t^(s)) ]
///            / [ Σ_{s: c(s)=c} 1 / (σ_t² + Σ̂_t^(s)) ]
pub fn update_condition_means(
    posteriors: &[SamplePosterior],
    hyperparams: &mut HierarchicalHyperparams,
) {
    let num_targets = hyperparams.sigma_sq.len();
    let num_conditions = hyperparams.condition_names.len();

    for c in 0..num_conditions {
        let mut numer = vec![0.0_f64; num_targets];
        let mut denom = vec![0.0_f64; num_targets];

        for post in posteriors.iter().filter(|p| p.condition_idx == c) {
            for t in 0..num_targets {
                let total_var = hyperparams.sigma_sq[t] + post.sigma_hat_sq[t];
                if total_var > 0.0 {
                    let precision = 1.0 / total_var;
                    numer[t] += post.phi_hat[t] * precision;
                    denom[t] += precision;
                }
            }
        }

        for t in 0..num_targets {
            hyperparams.nu[c][t] = if denom[t] > 0.0 {
                numer[t] / denom[t]
            } else {
                0.0
            };
        }
    }
}

/// Update biological variance (method-of-moments with max(0, ...) clamp).
///
///   σ̂_t² = max(0, (1/S) Σ_s (φ̂_t^(s) - ν̂_t^(c(s)))² - (1/S) Σ_s Σ̂_t^(s))
///
/// The max(0, ·) clamp ensures that when sample-to-sample variation is
/// dominated by quantification noise, the biological variance goes to zero
/// and the prior loses strength.
pub fn update_biological_variance(
    posteriors: &[SamplePosterior],
    hyperparams: &mut HierarchicalHyperparams,
) {
    let num_targets = hyperparams.sigma_sq.len();
    let num_samples = posteriors.len() as f64;

    if num_samples < 1.0 {
        return;
    }

    for t in 0..num_targets {
        // Compute (1/S) Σ_s (φ̂_t^(s) - ν̂_t^(c(s)))²
        let mut sum_sq_dev = 0.0_f64;
        let mut sum_sigma_hat = 0.0_f64;

        for post in posteriors {
            let nu_t = hyperparams.nu[post.condition_idx][t];
            let dev = post.phi_hat[t] - nu_t;
            sum_sq_dev += dev * dev;
            sum_sigma_hat += post.sigma_hat_sq[t];
        }

        let mean_sq_dev = sum_sq_dev / num_samples;
        let mean_sigma_hat = sum_sigma_hat / num_samples;

        // max(floor, ...) clamp — floor prevents division by zero in the prior term
        // when measurement noise dominates (common with few samples)
        const SIGMA_SQ_FLOOR: f64 = 1e-4;
        hyperparams.sigma_sq[t] = (mean_sq_dev - mean_sigma_hat).max(SIGMA_SQ_FLOOR);
    }
}

/// Compute convergence metrics for the outer loop.
///
/// Returns (max_nu_change, max_sigma_sq_change) as relative changes.
pub fn convergence_metrics(
    old_nu: &[Vec<f64>],
    old_sigma_sq: &[f64],
    new: &HierarchicalHyperparams,
) -> (f64, f64) {
    let mut max_nu_change = 0.0_f64;
    for (c, old_nu_c) in old_nu.iter().enumerate() {
        for (t, &old_val) in old_nu_c.iter().enumerate() {
            let new_val = new.nu[c][t];
            let rel_change = (new_val - old_val).abs() / (1.0 + old_val.abs());
            max_nu_change = max_nu_change.max(rel_change);
        }
    }

    let mut max_sigma_change = 0.0_f64;
    for (t, &old_val) in old_sigma_sq.iter().enumerate() {
        let new_val = new.sigma_sq[t];
        let rel_change = (new_val - old_val).abs() / (1.0 + old_val);
        max_sigma_change = max_sigma_change.max(rel_change);
    }

    (max_nu_change, max_sigma_change)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_posteriors(
        phi_hats: Vec<Vec<f64>>,
        sigma_hats: Vec<Vec<f64>>,
        conditions: Vec<usize>,
    ) -> Vec<SamplePosterior> {
        phi_hats
            .into_iter()
            .zip(sigma_hats)
            .zip(conditions)
            .enumerate()
            .map(|(i, ((phi, sigma), cond))| SamplePosterior {
                phi_hat: phi,
                sigma_hat_sq: sigma,
                sample_name: format!("sample_{}", i),
                condition_idx: cond,
            })
            .collect()
    }

    #[test]
    fn test_init_hyperparams() {
        let hp = init_hyperparams(100, vec!["ctrl".into(), "treat".into()]);
        assert_eq!(hp.nu.len(), 2);
        assert_eq!(hp.nu[0].len(), 100);
        assert_eq!(hp.sigma_sq.len(), 100);
        assert!(hp.nu[0].iter().all(|&v| v == 0.0));
        assert!(hp.sigma_sq.iter().all(|&v| v == 1.0));
    }

    #[test]
    fn test_update_condition_means_single_condition() {
        let num_targets = 3;
        let mut hp = init_hyperparams(num_targets, vec!["ctrl".into()]);

        // Two samples in the same condition
        let posteriors = make_posteriors(
            vec![vec![1.0, 2.0, 3.0], vec![3.0, 4.0, 5.0]],
            vec![vec![0.1, 0.1, 0.1], vec![0.1, 0.1, 0.1]],
            vec![0, 0],
        );

        update_condition_means(&posteriors, &mut hp);

        // With equal measurement variance, nu should be the mean of phi_hats
        for t in 0..num_targets {
            let expected = (posteriors[0].phi_hat[t] + posteriors[1].phi_hat[t]) / 2.0;
            assert!(
                (hp.nu[0][t] - expected).abs() < 1e-6,
                "Target {}: expected nu={:.4}, got {:.4}",
                t, expected, hp.nu[0][t]
            );
        }
    }

    #[test]
    fn test_update_condition_means_precision_weighting() {
        let num_targets = 1;
        let mut hp = init_hyperparams(num_targets, vec!["ctrl".into()]);
        hp.sigma_sq[0] = 0.0; // No biological variance → pure precision weighting

        // Sample 0 has small variance (precise), sample 1 has large variance (imprecise)
        let posteriors = make_posteriors(
            vec![vec![1.0], vec![5.0]],
            vec![vec![0.1], vec![10.0]],
            vec![0, 0],
        );

        update_condition_means(&posteriors, &mut hp);

        // Nu should be closer to the precise sample (1.0) than to the imprecise one (5.0)
        assert!(
            hp.nu[0][0] < 3.0,
            "Precision weighting should pull nu toward the precise sample, got {:.4}",
            hp.nu[0][0]
        );
    }

    #[test]
    fn test_update_biological_variance_max_zero_clamp() {
        let num_targets = 1;
        let mut hp = init_hyperparams(num_targets, vec!["ctrl".into()]);

        // Samples very close together with high measurement uncertainty
        // → biological variance should clamp to 0
        let posteriors = make_posteriors(
            vec![vec![1.0], vec![1.01]],
            vec![vec![10.0], vec![10.0]],
            vec![0, 0],
        );

        // First update nu
        update_condition_means(&posteriors, &mut hp);
        // Then update sigma_sq
        update_biological_variance(&posteriors, &mut hp);

        assert!(
            hp.sigma_sq[0] <= 1e-4,
            "When noise dominates, sigma_sq should clamp to floor, got {}",
            hp.sigma_sq[0]
        );
    }

    #[test]
    fn test_update_biological_variance_positive() {
        let num_targets = 1;
        let mut hp = init_hyperparams(num_targets, vec!["ctrl".into()]);

        // Samples far apart with low measurement uncertainty
        // → biological variance should be positive
        let posteriors = make_posteriors(
            vec![vec![0.0], vec![10.0]],
            vec![vec![0.01], vec![0.01]],
            vec![0, 0],
        );

        update_condition_means(&posteriors, &mut hp);
        update_biological_variance(&posteriors, &mut hp);

        assert!(
            hp.sigma_sq[0] > 0.0,
            "With real variation, sigma_sq should be positive, got {:.6}",
            hp.sigma_sq[0]
        );
    }

    #[test]
    fn test_convergence_metrics() {
        let old_nu = vec![vec![1.0, 2.0, 3.0]];
        let old_sigma_sq = vec![1.0, 1.0, 1.0];

        let new = HierarchicalHyperparams {
            nu: vec![vec![1.1, 2.0, 3.0]],
            sigma_sq: vec![1.0, 1.5, 1.0],
            condition_names: vec!["ctrl".into()],
        };

        let (nu_change, sigma_change) = convergence_metrics(&old_nu, &old_sigma_sq, &new);

        assert!(nu_change > 0.0);
        assert!(sigma_change > 0.0);
    }
}
