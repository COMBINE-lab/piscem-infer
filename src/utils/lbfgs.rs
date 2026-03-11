use anyhow::Result;
use argmin::core::{CostFunction, Executor, Gradient, State, TerminationStatus};
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::quasinewton::LBFGS;
use tracing::info;

use crate::utils::em::{em, em_par, EMInfo};
use crate::utils::eq_maps::{EqLabel, PackedEqMap};
use crate::utils::gradient;

/// Prior parameters for the penalized MAP objective.
pub struct PenalizedPrior {
    /// Condition-level mean per transcript (log-space)
    pub nu: Vec<f64>,
    /// Biological variance per transcript
    pub sigma_sq: Vec<f64>,
}

/// Per-sample posterior result from penalized EM + L-BFGS.
pub struct SamplePosterior {
    /// MAP estimate in log-space (softmax parameterization)
    pub phi_hat: Vec<f64>,
    /// Diagonal Laplace posterior variance per transcript
    pub sigma_hat_sq: Vec<f64>,
}

/// The negated penalized objective for use with argmin (which minimizes).
///
/// F(φ) = ℓ(softmax(φ)) - (1/2) Σ_t (φ_t - ν_t)² / σ_t²
///
/// We negate because argmin minimizes, and we want to maximize F.
struct PenalizedProblem<'a, EqLabelT: EqLabel> {
    eq_map: &'a PackedEqMap<EqLabelT>,
    inv_eff_lens: Vec<f64>,
    prior: &'a PenalizedPrior,
}

impl<EqLabelT: EqLabel> CostFunction for PenalizedProblem<'_, EqLabelT> {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, phi: &Self::Param) -> std::result::Result<Self::Output, argmin::core::Error> {
        let val = gradient::penalized_objective(
            self.eq_map,
            phi,
            &self.prior.nu,
            &self.prior.sigma_sq,
            &self.inv_eff_lens,
        );
        // Negate: argmin minimizes, we want to maximize F
        Ok(-val)
    }
}

impl<EqLabelT: EqLabel> Gradient for PenalizedProblem<'_, EqLabelT> {
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;

    fn gradient(
        &self,
        phi: &Self::Param,
    ) -> std::result::Result<Self::Gradient, argmin::core::Error> {
        let theta = gradient::softmax(phi);
        let grad = gradient::compute_gradient_phi(
            self.eq_map,
            &theta,
            phi,
            &self.prior.nu,
            &self.prior.sigma_sq,
            &self.inv_eff_lens,
        );
        // Negate for minimization
        Ok(grad.into_iter().map(|g| -g).collect())
    }
}

/// Run the three-phase per-sample inference:
///
/// Phase 1a: EM warm-start for `em_warmstart_iters` iterations.
/// Phase 1b: L-BFGS on the penalized objective F(φ) from the EM solution.
/// Phase 1c: Compute diagonal Fisher information and Laplace posterior variance.
///
/// Returns a `SamplePosterior` with the MAP estimate and uncertainty.
pub fn penalized_em<EqLabelT: EqLabel>(
    em_info: &EMInfo<'_, '_, EqLabelT>,
    prior: &PenalizedPrior,
    em_warmstart_iters: u32,
    lbfgs_max_iters: u64,
    lbfgs_history: usize,
    num_threads: usize,
) -> Result<SamplePosterior> {
    let num_targets = em_info.eff_lens.len();
    let inv_eff_lens = gradient::compute_inv_eff_lens(em_info.eff_lens);

    // Phase 1a: EM warm-start
    info!("Phase 1a: Running EM warm-start for {} iterations", em_warmstart_iters);
    let warmstart_info = EMInfo {
        eq_map: em_info.eq_map,
        eff_lens: em_info.eff_lens,
        max_iter: em_warmstart_iters,
        convergence_thresh: em_info.convergence_thresh,
        presence_thresh: em_info.presence_thresh,
    };

    let em_counts = if num_threads > 1 {
        em_par(&warmstart_info, num_threads)
    } else {
        em(&warmstart_info)
    };

    // Convert EM counts to theta (probability simplex), then to phi (log-space).
    // Use a moderate floor (1e-8) to keep phi in a reasonable range for L-BFGS,
    // then center phi to mean zero (softmax is shift-invariant).
    let total_counts: f64 = em_counts.iter().sum();
    let theta_floor = 1e-8;
    let theta: Vec<f64> = em_counts.iter().map(|&c| {
        let t = c / total_counts;
        if t > theta_floor { t } else { theta_floor }
    }).collect();
    let mut init_phi = gradient::log_transform(&theta, theta_floor);
    // Center phi so mean = 0 (stabilizes L-BFGS and prior interaction)
    let phi_mean = init_phi.iter().sum::<f64>() / init_phi.len() as f64;
    for p in init_phi.iter_mut() {
        *p -= phi_mean;
    }

    // Phase 1b: L-BFGS on penalized objective
    info!("Phase 1b: Running L-BFGS optimization (max {} iters, history {})",
          lbfgs_max_iters, lbfgs_history);

    let problem = PenalizedProblem {
        eq_map: em_info.eq_map,
        inv_eff_lens: inv_eff_lens.clone(),
        prior,
    };

    let linesearch = MoreThuenteLineSearch::new();
    let solver = LBFGS::new(linesearch, lbfgs_history)
        .with_tolerance_grad(1e-4 * (num_targets as f64).sqrt())
        .map_err(|e| anyhow::anyhow!("L-BFGS config error: {}", e))?
        .with_tolerance_cost(1e-10)
        .map_err(|e| anyhow::anyhow!("L-BFGS config error: {}", e))?;

    let executor = Executor::new(problem, solver)
        .configure(|state| state.param(init_phi).max_iters(lbfgs_max_iters));

    let result = executor.run()?;

    let final_phi = result
        .state()
        .get_best_param()
        .ok_or_else(|| anyhow::anyhow!("L-BFGS did not produce a result"))?
        .clone();

    let term_status = result.state().get_termination_status();
    let iters = result.state().get_iter();
    match term_status {
        TerminationStatus::Terminated(reason) => {
            info!("L-BFGS converged after {} iterations: {:?}", iters, reason);
        }
        TerminationStatus::NotTerminated => {
            info!("L-BFGS reached max iterations ({})", iters);
        }
    }

    // Phase 1c: Compute Fisher diagonal and Laplace posterior variance
    info!("Phase 1c: Computing Fisher information and posterior variance");
    let final_theta = gradient::softmax(&final_phi);
    let fisher_diag = gradient::compute_fisher_diag(em_info.eq_map, &final_theta, &inv_eff_lens);
    let sigma_hat_sq = gradient::laplace_posterior_variance(&fisher_diag, &prior.sigma_sq);

    Ok(SamplePosterior {
        phi_hat: final_phi,
        sigma_hat_sq,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::eq_maps::{BasicEqLabel, BasicEqMap, OrientationProperty, PackedEqMap};
    use crate::utils::gradient;

    /// Build the same synthetic EqMap used in gradient tests.
    fn build_test_eq_map() -> PackedEqMap<BasicEqLabel> {
        let mut eqm = BasicEqMap::new(OrientationProperty::OrientationAgnostic);
        for _ in 0..10 { eqm.add(BasicEqLabel::new(&[0, 1], None)); }
        for _ in 0..15 { eqm.add(BasicEqLabel::new(&[1, 2], None)); }
        for _ in 0..20 { eqm.add(BasicEqLabel::new(&[0], None)); }
        for _ in 0..5 { eqm.add(BasicEqLabel::new(&[2, 3, 4], None)); }
        for _ in 0..8 { eqm.add(BasicEqLabel::new(&[3], None)); }
        PackedEqMap::from_eq_map(&eqm)
    }

    #[test]
    fn test_flat_prior_matches_em() {
        let eq_map = build_test_eq_map();
        let num_targets = 5;
        let eff_lens = vec![100.0, 200.0, 150.0, 300.0, 250.0];

        // Run full EM (many iterations, to convergence)
        let em_info = EMInfo {
            eq_map: &eq_map,
            eff_lens: &eff_lens,
            max_iter: 1500,
            convergence_thresh: 1e-10,
            presence_thresh: 1e-8,
        };
        let em_counts = em(&em_info);
        let em_total: f64 = em_counts.iter().sum();
        let em_theta: Vec<f64> = em_counts.iter().map(|&c| c / em_total).collect();

        // Run penalized EM with very flat prior (sigma_sq = 1e30)
        let prior = PenalizedPrior {
            nu: vec![0.0; num_targets],
            sigma_sq: vec![1e30; num_targets],
        };

        let posterior = penalized_em(&em_info, &prior, 20, 200, 7, 1).unwrap();
        let lbfgs_theta = gradient::softmax(&posterior.phi_hat);

        // Compare theta values
        for t in 0..num_targets {
            let rel_diff = if em_theta[t] > 1e-10 {
                (em_theta[t] - lbfgs_theta[t]).abs() / em_theta[t]
            } else {
                (em_theta[t] - lbfgs_theta[t]).abs()
            };
            assert!(
                rel_diff < 1e-3,
                "Target {}: EM theta={:.8}, L-BFGS theta={:.8}, rel_diff={:.2e}",
                t, em_theta[t], lbfgs_theta[t], rel_diff
            );
        }
    }

    #[test]
    fn test_tight_prior_pulls_toward_nu() {
        let eq_map = build_test_eq_map();
        let num_targets = 5;
        let eff_lens = vec![100.0, 200.0, 150.0, 300.0, 250.0];

        let em_info = EMInfo {
            eq_map: &eq_map,
            eff_lens: &eff_lens,
            max_iter: 1500,
            convergence_thresh: 1e-10,
            presence_thresh: 1e-8,
        };

        // First run with flat prior to get the unconstrained solution
        let prior_flat = PenalizedPrior {
            nu: vec![0.0; num_targets],
            sigma_sq: vec![1e30; num_targets],
        };
        let posterior_flat = penalized_em(&em_info, &prior_flat, 20, 200, 7, 1).unwrap();

        // Now run with tight prior pulling toward a specific target
        let nu_target = vec![1.0, -1.0, 0.5, -0.5, 0.0];
        let prior_tight = PenalizedPrior {
            nu: nu_target.clone(),
            sigma_sq: vec![1e-3; num_targets],
        };
        let posterior_tight = penalized_em(&em_info, &prior_tight, 20, 200, 7, 1).unwrap();

        // With tight prior, phi should be much closer to nu than the flat solution
        let mut dist_to_nu_flat = 0.0_f64;
        let mut dist_to_nu_tight = 0.0_f64;
        for t in 0..num_targets {
            dist_to_nu_flat += (posterior_flat.phi_hat[t] - nu_target[t]).powi(2);
            dist_to_nu_tight += (posterior_tight.phi_hat[t] - nu_target[t]).powi(2);
        }

        assert!(
            dist_to_nu_tight < dist_to_nu_flat,
            "Tight prior should pull phi closer to nu: tight_dist={:.4}, flat_dist={:.4}",
            dist_to_nu_tight.sqrt(),
            dist_to_nu_flat.sqrt()
        );
    }

    #[test]
    fn test_posterior_variance_decreases_with_data() {
        let num_targets = 5;
        let eff_lens = vec![100.0, 200.0, 150.0, 300.0, 250.0];

        let prior = PenalizedPrior {
            nu: vec![0.0; num_targets],
            sigma_sq: vec![1.0; num_targets],
        };

        // Build a small dataset
        let mut eqm_small = BasicEqMap::new(OrientationProperty::OrientationAgnostic);
        for _ in 0..5 { eqm_small.add(BasicEqLabel::new(&[0, 1], None)); }
        for _ in 0..3 { eqm_small.add(BasicEqLabel::new(&[2], None)); }
        let eq_map_small = PackedEqMap::from_eq_map(&eqm_small);

        let em_info_small = EMInfo {
            eq_map: &eq_map_small,
            eff_lens: &eff_lens,
            max_iter: 1500,
            convergence_thresh: 1e-10,
            presence_thresh: 1e-8,
        };
        let post_small = penalized_em(&em_info_small, &prior, 20, 200, 7, 1).unwrap();

        // Build a larger dataset (10x counts)
        let mut eqm_large = BasicEqMap::new(OrientationProperty::OrientationAgnostic);
        for _ in 0..50 { eqm_large.add(BasicEqLabel::new(&[0, 1], None)); }
        for _ in 0..30 { eqm_large.add(BasicEqLabel::new(&[2], None)); }
        let eq_map_large = PackedEqMap::from_eq_map(&eqm_large);

        let em_info_large = EMInfo {
            eq_map: &eq_map_large,
            eff_lens: &eff_lens,
            max_iter: 1500,
            convergence_thresh: 1e-10,
            presence_thresh: 1e-8,
        };
        let post_large = penalized_em(&em_info_large, &prior, 20, 200, 7, 1).unwrap();

        // For target 2 (has its own unambiguous EQC), posterior variance
        // should be smaller with more data. Targets 0 and 1 share an
        // ambiguous EQC so the relationship is less clean.
        assert!(
            post_large.sigma_hat_sq[2] < post_small.sigma_hat_sq[2],
            "Target 2: more data should give smaller variance: large={:.6}, small={:.6}",
            post_large.sigma_hat_sq[2], post_small.sigma_hat_sq[2]
        );

        // All targets with data should have posterior variance < prior variance (1.0)
        for t in 0..3 {
            assert!(
                post_large.sigma_hat_sq[t] < 1.0,
                "Target {}: posterior variance ({:.6}) should be < prior (1.0) with data",
                t, post_large.sigma_hat_sq[t]
            );
        }
    }
}
