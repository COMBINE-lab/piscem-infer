use crate::utils::eq_maps::{EqLabel, PackedEqMap, TargetLabelsRef};

/// Minimum threshold for theta values to avoid log(0) and division by zero.
const THETA_MIN: f64 = 1e-300;

/// Minimum threshold for mu_e (per-EQC denominator) to skip degenerate EQCs.
const MU_MIN: f64 = 1e-300;

/// Numerically stable softmax: theta_t = exp(phi_t) / sum_t' exp(phi_t').
/// Subtracts max(phi) before exponentiating for numerical stability.
pub fn softmax(phi: &[f64]) -> Vec<f64> {
    let max_phi = phi.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut theta: Vec<f64> = phi.iter().map(|&p| (p - max_phi).exp()).collect();
    let sum: f64 = theta.iter().sum();
    let inv_sum = 1.0 / sum;
    for t in theta.iter_mut() {
        *t *= inv_sum;
    }
    theta
}

/// Inverse of softmax: phi_t = log(theta_t), with a floor to avoid log(0).
/// Note: the resulting phi is only unique up to a constant shift.
pub fn log_transform(theta: &[f64], floor: f64) -> Vec<f64> {
    theta
        .iter()
        .map(|&t| if t > floor { t.ln() } else { floor.ln() })
        .collect()
}

/// Compute the log-likelihood of the RFEQ model:
///   ℓ(θ) = Σ_e c_e · log(μ_e)
/// where μ_e = Σ_{t ∈ T_e} θ_t · w[t,e]
///
/// This is used as part of the penalized objective F(φ).
pub fn log_likelihood<EqLabelT: EqLabel>(
    eq_map: &PackedEqMap<EqLabelT>,
    theta: &[f64],
    inv_eff_lens: &[f64],
) -> f64 {
    let mut ll = 0.0_f64;
    for (k, &count) in eq_map.iter_labels().zip(eq_map.counts.iter()) {
        let c_e = count as f64;
        let mut mu_e = 0.0_f64;
        for (&target, cond_prob) in k.target_labels().iter().zip(k.target_probs()) {
            mu_e += cond_prob * theta[target as usize] * inv_eff_lens[target as usize];
        }
        if mu_e > MU_MIN {
            ll += c_e * mu_e.ln();
        }
    }
    ll
}

/// Compute the penalized objective value:
///   F(φ) = ℓ(softmax(φ)) - (1/2) Σ_t (φ_t - ν_t)² / σ_t²
pub fn penalized_objective<EqLabelT: EqLabel>(
    eq_map: &PackedEqMap<EqLabelT>,
    phi: &[f64],
    nu: &[f64],
    sigma_sq: &[f64],
    inv_eff_lens: &[f64],
) -> f64 {
    let theta = softmax(phi);
    let ll = log_likelihood(eq_map, &theta, inv_eff_lens);
    let penalty: f64 = phi
        .iter()
        .zip(nu.iter())
        .zip(sigma_sq.iter())
        .map(|((&p, &n), &s)| {
            let diff = p - n;
            0.5 * diff * diff / s
        })
        .sum();
    ll - penalty
}

/// Compute the gradient of the penalized objective F(φ) in φ-space.
///
/// The log-likelihood gradient through the softmax parameterization is:
///   ∂ℓ/∂φ_t = θ_t · (score_t - N)
///
/// where score_t = Σ_{e: t∈T_e} c_e · w_te / μ_e, N = Σ_e c_e (total fragments),
/// w_te = cond_prob[t,e] · inv_eff_len[t], and μ_e = Σ_{t ∈ T_e} θ_t · w_te.
///
/// This follows from the softmax Jacobian: ∂θ_{t'}/∂φ_t = θ_{t'}(δ_{tt'} - θ_t),
/// which gives ∂ℓ/∂φ_t = θ_t · (∂ℓ/∂θ_t - Σ_{t'} θ_{t'} · ∂ℓ/∂θ_{t'}).
/// The weighted sum Σ_{t'} θ_{t'} · ∂ℓ/∂θ_{t'} telescopes to N because
/// Σ_{t'} θ_{t'} · w_{t'e} = μ_e.
///
/// The full penalized gradient adds the prior term:
///   ∂F/∂φ_t = θ_t · (score_t - N) - (φ_t - ν_t) / σ_t²
pub fn compute_gradient_phi<EqLabelT: EqLabel>(
    eq_map: &PackedEqMap<EqLabelT>,
    theta: &[f64],
    phi: &[f64],
    nu: &[f64],
    sigma_sq: &[f64],
    inv_eff_lens: &[f64],
) -> Vec<f64> {
    let num_targets = theta.len();
    // First compute score_t = Σ_{e: t∈T_e} c_e · w_te / μ_e
    let mut score = vec![0.0_f64; num_targets];
    let total_n: f64 = eq_map.counts.iter().sum::<usize>() as f64;

    for (k, &count) in eq_map.iter_labels().zip(eq_map.counts.iter()) {
        let c_e = count as f64;

        // Compute μ_e = Σ_t θ_t · w_te
        let mut mu_e = 0.0_f64;
        for (&target, cond_prob) in k.target_labels().iter().zip(k.target_probs()) {
            let t = target as usize;
            mu_e += cond_prob * theta[t] * inv_eff_lens[t];
        }

        if mu_e < MU_MIN {
            continue;
        }

        let c_over_mu = c_e / mu_e;

        for (&target, cond_prob) in k.target_labels().iter().zip(k.target_probs()) {
            let t = target as usize;
            let w_te = cond_prob * inv_eff_lens[t];
            score[t] += c_over_mu * w_te;
        }
    }

    // Gradient: θ_t · (score_t - N) - prior term
    let mut grad = vec![0.0_f64; num_targets];
    for t in 0..num_targets {
        grad[t] = theta[t] * (score[t] - total_n) - (phi[t] - nu[t]) / sigma_sq[t];
    }

    grad
}

/// Compute the gradient of the log-likelihood only (no prior term).
///   ∂ℓ/∂φ_t = θ_t · (score_t - N)
/// Used for verifying convergence of EM to a stationary point.
#[cfg(test)]
pub fn compute_gradient_ll<EqLabelT: EqLabel>(
    eq_map: &PackedEqMap<EqLabelT>,
    theta: &[f64],
    inv_eff_lens: &[f64],
) -> Vec<f64> {
    let num_targets = theta.len();
    let mut score = vec![0.0_f64; num_targets];
    let total_n: f64 = eq_map.counts.iter().sum::<usize>() as f64;

    for (k, &count) in eq_map.iter_labels().zip(eq_map.counts.iter()) {
        let c_e = count as f64;

        let mut mu_e = 0.0_f64;
        for (&target, cond_prob) in k.target_labels().iter().zip(k.target_probs()) {
            let t = target as usize;
            mu_e += cond_prob * theta[t] * inv_eff_lens[t];
        }

        if mu_e < MU_MIN {
            continue;
        }

        let c_over_mu = c_e / mu_e;

        for (&target, cond_prob) in k.target_labels().iter().zip(k.target_probs()) {
            let t = target as usize;
            let w_te = cond_prob * inv_eff_lens[t];
            score[t] += c_over_mu * w_te;
        }
    }

    let mut grad = vec![0.0_f64; num_targets];
    for t in 0..num_targets {
        grad[t] = theta[t] * (score[t] - total_n);
    }
    grad
}

/// Compute the diagonal Fisher information in φ-space.
///
/// In θ-space: I_tt^(θ) = Σ_{e: t∈T_e} c_e · (w[t,e] · inv_eff_len[t])² / μ_e²
/// Transform to φ-space: I_tt^(φ) = θ_t² · (1 - θ_t)² · I_tt^(θ)
///
/// Uses exact Jacobian factor θ_t(1-θ_t) rather than the θ_t approximation.
pub fn compute_fisher_diag<EqLabelT: EqLabel>(
    eq_map: &PackedEqMap<EqLabelT>,
    theta: &[f64],
    inv_eff_lens: &[f64],
) -> Vec<f64> {
    let num_targets = theta.len();
    let mut fisher_theta = vec![0.0_f64; num_targets];

    // Compute Fisher in θ-space
    for (k, &count) in eq_map.iter_labels().zip(eq_map.counts.iter()) {
        let c_e = count as f64;

        // Compute μ_e
        let mut mu_e = 0.0_f64;
        for (&target, cond_prob) in k.target_labels().iter().zip(k.target_probs()) {
            let t = target as usize;
            mu_e += cond_prob * theta[t] * inv_eff_lens[t];
        }

        if mu_e < MU_MIN {
            continue;
        }

        let inv_mu_sq = 1.0 / (mu_e * mu_e);

        for (&target, cond_prob) in k.target_labels().iter().zip(k.target_probs()) {
            let t = target as usize;
            let w_t_e = cond_prob * inv_eff_lens[t];
            fisher_theta[t] += c_e * w_t_e * w_t_e * inv_mu_sq;
        }
    }

    // Transform to φ-space: I_tt^(φ) = θ_t² · (1 - θ_t)² · I_tt^(θ)
    let mut fisher_phi = vec![0.0_f64; num_targets];
    for t in 0..num_targets {
        let th = theta[t].max(THETA_MIN);
        let jacobian_sq = th * th * (1.0 - th) * (1.0 - th);
        fisher_phi[t] = jacobian_sq * fisher_theta[t];
    }

    fisher_phi
}

/// Compute the Laplace posterior variance per transcript:
///   Σ̂_t = 1 / (1/σ_t² + I_tt^(φ))
///
/// When I_tt^(φ) = 0 (unidentifiable transcript), Σ̂_t = σ_t² (prior dominates).
/// When σ_t² is very large (vague prior), Σ̂_t ≈ 1/I_tt^(φ) (data dominates).
pub fn laplace_posterior_variance(fisher_diag: &[f64], sigma_sq: &[f64]) -> Vec<f64> {
    fisher_diag
        .iter()
        .zip(sigma_sq.iter())
        .map(|(&fish, &sig)| {
            let precision = 1.0 / sig + fish;
            if precision > 0.0 {
                1.0 / precision
            } else {
                sig
            }
        })
        .collect()
}

/// Compute inverse effective lengths, mapping non-finite values to 0.
pub fn compute_inv_eff_lens(eff_lens: &[f64]) -> Vec<f64> {
    eff_lens
        .iter()
        .map(|&x| {
            let y = 1.0 / x;
            if y.is_finite() { y } else { 0.0 }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::eq_maps::{BasicEqLabel, BasicEqMap, OrientationProperty, PackedEqMap};

    /// Build a small synthetic EqMap for testing.
    /// 5 targets, several EQCs with known structure.
    fn build_test_eq_map() -> PackedEqMap<BasicEqLabel> {
        let mut eqm = BasicEqMap::new(OrientationProperty::OrientationAgnostic);

        // EQC 0: targets {0, 1}, count = 10
        for _ in 0..10 {
            eqm.add(BasicEqLabel::new(&[0, 1], None));
        }
        // EQC 1: targets {1, 2}, count = 15
        for _ in 0..15 {
            eqm.add(BasicEqLabel::new(&[1, 2], None));
        }
        // EQC 2: targets {0}, count = 20
        for _ in 0..20 {
            eqm.add(BasicEqLabel::new(&[0], None));
        }
        // EQC 3: targets {2, 3, 4}, count = 5
        for _ in 0..5 {
            eqm.add(BasicEqLabel::new(&[2, 3, 4], None));
        }
        // EQC 4: targets {3}, count = 8
        for _ in 0..8 {
            eqm.add(BasicEqLabel::new(&[3], None));
        }

        PackedEqMap::from_eq_map(&eqm)
    }

    #[test]
    fn test_softmax_basic() {
        let phi = vec![1.0, 2.0, 3.0];
        let theta = softmax(&phi);

        // Check sums to 1
        let sum: f64 = theta.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);

        // Check ordering preserved
        assert!(theta[2] > theta[1]);
        assert!(theta[1] > theta[0]);
    }

    #[test]
    fn test_softmax_numerical_stability() {
        // Very large values should not overflow
        let phi = vec![1000.0, 1001.0, 1002.0];
        let theta = softmax(&phi);
        let sum: f64 = theta.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);

        // Very negative values should not underflow to all zeros
        let phi2 = vec![-1000.0, -1001.0, -1002.0];
        let theta2 = softmax(&phi2);
        let sum2: f64 = theta2.iter().sum();
        assert!((sum2 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_log_transform_roundtrip() {
        let phi_orig = vec![1.0, 2.0, 3.0, 0.5, -1.0];
        let theta = softmax(&phi_orig);
        let phi_recovered = log_transform(&theta, THETA_MIN);
        let theta_recovered = softmax(&phi_recovered);

        // softmax(log(theta)) should recover theta
        for i in 0..theta.len() {
            assert!(
                (theta[i] - theta_recovered[i]).abs() < 1e-12,
                "Mismatch at index {}: {} vs {}",
                i,
                theta[i],
                theta_recovered[i]
            );
        }
    }

    #[test]
    fn test_finite_difference_gradient() {
        let eq_map = build_test_eq_map();
        let num_targets = 5;
        let eff_lens = vec![100.0, 200.0, 150.0, 300.0, 250.0];
        let inv_eff_lens = compute_inv_eff_lens(&eff_lens);

        // Set up a test point
        let phi = vec![0.5, -0.3, 0.1, 0.8, -0.5];
        let theta = softmax(&phi);

        // Prior parameters (moderate prior)
        let nu = vec![0.0; num_targets];
        let sigma_sq = vec![1.0; num_targets];

        // Compute analytic gradient
        let grad_analytic =
            compute_gradient_phi(&eq_map, &theta, &phi, &nu, &sigma_sq, &inv_eff_lens);

        // Compute finite-difference gradient
        let h = 1e-7;
        let mut grad_fd = vec![0.0_f64; num_targets];
        for t in 0..num_targets {
            let mut phi_plus = phi.clone();
            let mut phi_minus = phi.clone();
            phi_plus[t] += h;
            phi_minus[t] -= h;

            let f_plus =
                penalized_objective(&eq_map, &phi_plus, &nu, &sigma_sq, &inv_eff_lens);
            let f_minus =
                penalized_objective(&eq_map, &phi_minus, &nu, &sigma_sq, &inv_eff_lens);
            grad_fd[t] = (f_plus - f_minus) / (2.0 * h);
        }

        // Check agreement
        for t in 0..num_targets {
            let abs_diff = (grad_analytic[t] - grad_fd[t]).abs();
            let scale = grad_fd[t].abs().max(1.0);
            let rel_diff = abs_diff / scale;
            assert!(
                rel_diff < 1e-5,
                "Gradient mismatch at target {}: analytic={:.10}, fd={:.10}, rel_diff={:.2e}",
                t,
                grad_analytic[t],
                grad_fd[t],
                rel_diff
            );
        }
    }

    #[test]
    fn test_fisher_non_negative() {
        let eq_map = build_test_eq_map();
        let eff_lens = vec![100.0, 200.0, 150.0, 300.0, 250.0];
        let inv_eff_lens = compute_inv_eff_lens(&eff_lens);
        let phi = vec![0.5, -0.3, 0.1, 0.8, -0.5];
        let theta = softmax(&phi);

        let fisher = compute_fisher_diag(&eq_map, &theta, &inv_eff_lens);

        for (t, &f) in fisher.iter().enumerate() {
            assert!(
                f >= 0.0,
                "Fisher information negative at target {}: {}",
                t,
                f
            );
        }
    }

    #[test]
    fn test_laplace_posterior_variance() {
        // With zero Fisher (no data), posterior variance = prior variance
        let fisher = vec![0.0, 0.0, 0.0];
        let sigma_sq = vec![1.0, 2.0, 0.5];
        let post_var = laplace_posterior_variance(&fisher, &sigma_sq);
        for i in 0..3 {
            assert!(
                (post_var[i] - sigma_sq[i]).abs() < 1e-12,
                "With zero Fisher, posterior variance should equal prior"
            );
        }

        // With large Fisher (lots of data), posterior variance << prior
        let fisher_large = vec![1e6, 1e6, 1e6];
        let post_var_large = laplace_posterior_variance(&fisher_large, &sigma_sq);
        for i in 0..3 {
            assert!(
                post_var_large[i] < sigma_sq[i] * 0.01,
                "With large Fisher, posterior variance should be much smaller than prior"
            );
        }
    }

    #[test]
    fn test_gradient_near_zero_at_em_optimum() {
        // Run a simple manual EM to find the MLE, then verify that the
        // log-likelihood gradient in phi-space is approximately zero.
        let eq_map = build_test_eq_map();
        let num_targets = 5;
        let eff_lens = vec![100.0, 200.0, 150.0, 300.0, 250.0];
        let inv_eff_lens = compute_inv_eff_lens(&eff_lens);

        // Run EM manually (no prior)
        let total_weight: f64 = eq_map.counts.iter().sum::<usize>() as f64;
        let avg = total_weight / (num_targets as f64);
        let mut counts = vec![avg; num_targets];
        let mut new_counts = vec![0.0_f64; num_targets];

        for _ in 0..500 {
            new_counts.fill(0.0);
            for (k, &count) in eq_map.iter_labels().zip(eq_map.counts.iter()) {
                let c_e = count as f64;
                let mut denom = 0.0_f64;
                let mut weights = Vec::new();
                for (&target, cond_prob) in k.target_labels().iter().zip(k.target_probs()) {
                    let w = cond_prob * counts[target as usize] * inv_eff_lens[target as usize];
                    weights.push((target, w));
                    denom += w;
                }
                if denom > 1e-8 {
                    for &(target, w) in &weights {
                        new_counts[target as usize] += c_e * w / denom;
                    }
                }
            }
            std::mem::swap(&mut counts, &mut new_counts);
        }

        // Convert EM counts to theta (probability simplex)
        let total: f64 = counts.iter().sum();
        let theta: Vec<f64> = counts.iter().map(|&c| c / total).collect();

        // Compute gradient of log-likelihood only (no prior)
        let grad = compute_gradient_ll(&eq_map, &theta, &inv_eff_lens);

        let grad_norm: f64 = grad.iter().map(|g| g * g).sum::<f64>().sqrt();
        assert!(
            grad_norm < 1e-4,
            "Gradient norm at EM optimum should be near zero, got {:.6e}",
            grad_norm
        );
    }
}
