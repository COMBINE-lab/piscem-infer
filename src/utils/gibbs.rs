use rand::prelude::*;
use rand::rng;
use rand_distr::weighted::WeightedIndex;
use rand_distr::Gamma;
use rayon::prelude::*;
use tracing::info;

use crate::utils::em::EMInfo;
use crate::utils::eq_maps::{EqLabel, PackedEqMap, TargetLabelsRef};

pub fn do_gibbs<EqLabelT: EqLabel>(
    em_info: &EMInfo<EqLabelT>,
    em_result: &[f64],
    num_samples: usize,
    thinning_factor: usize,
) -> Vec<Vec<f64>> {
    let eq_map = em_info.eq_map;
    let eff_lens = em_info.eff_lens;
    let num_targets = eff_lens.len();

    // Per-nucleotide Dirichlet-like prior: alpha_prior / effLen[i]
    const ALPHA_PRIOR: f64 = 1e-3;
    let prior: Vec<f64> = eff_lens
        .iter()
        .map(|el| if *el > 0.0 { ALPHA_PRIOR / *el } else { 0.0 })
        .collect();

    // Determine number of chains (matching salmon's adaptive scheme)
    let num_chains = if num_samples >= 200 {
        8
    } else if num_samples >= 100 {
        4
    } else if num_samples >= 50 {
        2
    } else {
        1
    };
    let samples_per_chain = num_samples / num_chains;
    let remainder = num_samples % num_chains;

    // Run chains in parallel, each chain produces samples sequentially
    let all_samples: Vec<Vec<Vec<f64>>> = (0..num_chains)
        .into_par_iter()
        .map(|chain_id| {
            let my_samples = samples_per_chain + if chain_id < remainder { 1 } else { 0 };
            if my_samples == 0 {
                return vec![];
            }

            info!("Gibbs chain {} collecting {} samples", chain_id, my_samples);
            let mut rng = rng();
            let mut counts: Vec<f64> = em_result.to_vec();
            let mut mu = vec![0.0_f64; num_targets];
            let mut weights: Vec<f64> = Vec::with_capacity(64);
            let mut chain_samples: Vec<Vec<f64>> = Vec::with_capacity(my_samples);

            for sample_idx in 0..my_samples {
                // Run thinning_factor internal Gibbs iterations per collected sample
                for _ in 0..thinning_factor {
                    gibbs_iteration(
                        eq_map,
                        eff_lens,
                        &prior,
                        &mut counts,
                        &mut mu,
                        &mut weights,
                        &mut rng,
                    );
                }

                chain_samples.push(counts.clone());

                if (sample_idx + 1) % 50 == 0 {
                    info!(
                        "Gibbs chain {}: collected {}/{} samples",
                        chain_id,
                        sample_idx + 1,
                        my_samples
                    );
                }
            }

            chain_samples
        })
        .collect();

    // Flatten: collect all chain samples into a single vector
    let total_samples: usize = all_samples.iter().map(|c| c.len()).sum();
    let mut result: Vec<Vec<f64>> = Vec::with_capacity(total_samples);
    for chain_samples in all_samples {
        result.extend(chain_samples);
    }

    result
}

/// Single Gibbs iteration: Gamma step + multinomial reassignment
fn gibbs_iteration<EqLabelT: EqLabel>(
    eq_map: &PackedEqMap<EqLabelT>,
    eff_lens: &[f64],
    prior: &[f64],
    counts: &mut [f64],
    mu: &mut [f64],
    weights: &mut Vec<f64>,
    rng: &mut impl Rng,
) {
    const BETA: f64 = 0.1;

    // Gamma step: draw new transcript fractions
    let mut mu_sum = 0.0_f64;
    for i in 0..counts.len() {
        let shape = prior[i] + counts[i];
        if shape > 0.0 {
            let scale = 1.0 / (BETA + eff_lens[i]);
            mu[i] = Gamma::new(shape, scale).unwrap().sample(rng);
            mu_sum += mu[i];
        } else {
            mu[i] = 0.0;
        }
    }

    // Normalize mu for numerical stability
    if mu_sum > 0.0 {
        let inv_sum = 1.0 / mu_sum;
        mu.iter_mut().for_each(|x| *x *= inv_sum);
    }

    // Multinomial reassignment step
    counts.fill(0.0);
    for (label, eq_count) in eq_map.iter_labels().zip(eq_map.counts.iter()) {
        let targets = label.target_labels();
        let count = *eq_count;

        if targets.len() == 1 {
            // Unambiguous: assign all reads directly
            counts[targets[0] as usize] += count as f64;
            continue;
        }

        // Compute weights for ambiguous equivalence class
        weights.clear();
        let mut denom = 0.0_f64;
        for (tid, cp) in targets.iter().zip(label.target_probs()) {
            let w = mu[*tid as usize] * cp;
            weights.push(w);
            denom += w;
        }

        if denom <= f64::MIN_POSITIVE {
            // All weights essentially zero — skip this eq class
            continue;
        }

        // Sample from categorical distribution for each read
        if let Ok(dist) = WeightedIndex::new(weights.as_slice()) {
            for _ in 0..count {
                let chosen = dist.sample(rng);
                counts[targets[chosen] as usize] += 1.0;
            }
        }

        weights.clear();
    }
}
