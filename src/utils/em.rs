use atomic_float::AtomicF64;
use rand::prelude::*;
use rand::rng;
use rand_distr::weighted::{WeightedAliasIndex, WeightedIndex};
use rand_distr::Gamma;
use rayon::prelude::*;
use std::sync::atomic::Ordering;
use tracing::info;

use crate::utils::eq_maps::{EqLabel, PackedEqMap, TargetLabelsRef};

/// Takes as input `freq`, the frequency of occurrence of
/// fragments of each observed length, and returns
/// the conditional mean of the fragment length distribution
/// at every value 1 <= i <= freq.len()
pub fn conditional_means(freq: &[u32]) -> Vec<f64> {
    let mut cond_means = vec![0.0f64; freq.len()];
    let mut vals = vec![0.0f64; freq.len()];
    let mut multiplicities = vec![0.0f64; freq.len()];

    multiplicities[0] = freq[0] as f64;
    for i in 1..(freq.len()) {
        let v = freq[i] as f64;
        vals[i] = (v * i as f64) + vals[i - 1];
        multiplicities[i] = v + multiplicities[i - 1];
        if multiplicities[i] > 0.0f64 {
            cond_means[i] = vals[i] / multiplicities[i];
        }
    }

    cond_means
}

/// Takes as input the parameters of a truncated Normal distribution
/// with lower bound 0, upper bound `upper`, mean `mu` and standard
/// deviation `sigma`, and returns the conditional mean of the fragment
/// length distribution at every value 1 <= i <= `upper`
pub fn conditional_means_from_params(mu: f64, sigma: f64, upper: usize) -> Vec<f64> {
    let mut cond_means = vec![0.0f64; upper];
    let mut vals = vec![0.0f64; upper];
    let mut multiplicities = vec![0.0f64; upper];

    let inv_sigma = 1.0 / sigma;
    let denom_b = distrs::Normal::cdf(upper as f64, mu, sigma);
    let denom_a = distrs::Normal::cdf(0.0_f64, mu, sigma);
    let denom = denom_b - denom_a;
    let inv_denom = 1.0_f64 / denom;

    let trunc_pdf = |i: usize| -> f64 {
        let x = i as f64;
        inv_sigma * (distrs::Normal::pdf(x, mu, sigma) * inv_denom)
    };

    multiplicities[0] = trunc_pdf(0);
    for i in 1..upper {
        let v = trunc_pdf(i);
        vals[i] = (v * i as f64) + vals[i - 1];
        multiplicities[i] = v + multiplicities[i - 1];
        if multiplicities[i] > 0.0f64 {
            cond_means[i] = vals[i] / multiplicities[i];
        }
    }

    cond_means
}

/// Go through the set of references (`ref_lens`), and adjust their lengths according to
/// the computed conditional means `cond_means` of the fragment length distribution.
pub fn adjust_ref_lengths(ref_lens: &[u32], cond_means: &[f64]) -> Vec<f64> {
    let tmean = cond_means.last().unwrap();
    ref_lens
        .iter()
        .map(|rli| {
            let rl = *rli as usize;
            let adj_len = if rl > cond_means.len() {
                (rl as f64) - tmean
            } else {
                (rl as f64) - cond_means[rl]
            };
            if adj_len >= 1.0 { adj_len } else { rl as f64 }
        })
        .collect::<Vec<f64>>()
}

#[inline]
fn m_step_par<EqLabelT: EqLabel>(
    eq_iterates: &[(EqLabelT::LabelRefT<'_>, &usize)],
    prev_count: &mut [AtomicF64],
    inv_eff_lens: &[f64],
    curr_counts: &mut [AtomicF64],
) {
    // TODO: is there a better way to set the capacity on
    // this Vec?
    eq_iterates.par_iter().for_each_with(
        (&curr_counts, Vec::with_capacity(64)),
        |(curr_counts, weights), (k, v)| {
            let count = **v as f64;

            let mut denom = 0.0_f64;
            for (e, cond_prob) in k.target_labels().iter().zip(k.target_probs()) {
                let w = cond_prob
                    * prev_count[*e as usize].load(Ordering::Relaxed)
                    * inv_eff_lens[*e as usize];
                weights.push(w);
                denom += w;
            }
            if denom > 1e-8 {
                let count_over_denom = count / denom;
                for (target_id, w) in k.target_labels().iter().zip(weights.iter()) {
                    let inc = count_over_denom * w;
                    curr_counts[*target_id as usize].fetch_add(inc, Ordering::AcqRel);
                }
            }
            weights.clear();
        },
    );
}

#[inline]
fn m_step<EqLabelT: EqLabel>(
    eq_map: &PackedEqMap<EqLabelT>,
    eq_counts: &[usize],
    prev_count: &[f64],
    inv_eff_lens: &[f64],
    curr_counts: &mut [f64],
) {
    // TODO: is there a better way to set the capacity on
    // this Vec?
    let mut weights: Vec<f64> = Vec::with_capacity(64);

    for (k, v) in eq_map.iter_labels().zip(eq_counts.iter()) {
        let count = *v as f64;

        let mut denom = 0.0_f64;
        for (e, cond_prob) in k.target_labels().iter().zip(k.target_probs()) {
            let w = cond_prob * prev_count[*e as usize] * inv_eff_lens[*e as usize];
            weights.push(w);
            denom += w;
        }
        if denom > 1e-8 {
            let count_over_denom = count / denom;
            for (target_id, w) in k.target_labels().iter().zip(weights.iter()) {
                curr_counts[*target_id as usize] += count_over_denom * w;
            }
        }
        weights.clear();
    }
}

/// Holds the info relevant for running the EM algorithm
pub struct EMInfo<'eqm, 'el, EqLabelT> {
    pub eq_map: &'eqm PackedEqMap<EqLabelT>,
    pub eff_lens: &'el [f64],
    pub max_iter: u32,
    pub convergence_thresh: f64,
    pub presence_thresh: f64,
}

pub fn do_bootstrap<EqLabelT: EqLabel>(
    em_info: &EMInfo<EqLabelT>,
    num_boot: usize,
) -> Vec<Vec<f64>> {
    let converge_thresh = em_info.convergence_thresh;
    let presence_thresh = em_info.presence_thresh;
    let eq_map = em_info.eq_map;
    let max_iter = em_info.max_iter;
    let eff_lens = em_info.eff_lens;
    let inv_eff_lens = eff_lens
        .iter()
        .map(|x| {
            let y = 1.0_f64 / *x;
            if y.is_finite() { y } else { 0_f64 }
        })
        .collect::<Vec<f64>>();
    let total_weight = em_info.eq_map.counts.iter().sum::<usize>();
    // init
    let avg = (total_weight as f64) / (eff_lens.len() as f64);
    let dist = WeightedAliasIndex::new(em_info.eq_map.counts.clone()).unwrap();

    (0..num_boot)
        .into_par_iter()
        .map(|i| {
            info!("evaluating bootstrap replicate {}", i);
            let mut prev_counts = vec![avg; eff_lens.len()];
            let mut curr_counts = vec![0.0f64; eff_lens.len()];

            let mut rel_diff = 0.0_f64;
            let mut niter = 0_u32;
            let mut rng = rng();
            let mut base_counts = vec![0_usize; eq_map.counts.len()];
            for _s in 0..total_weight {
                base_counts[dist.sample(&mut rng)] += 1;
            }

            while niter < max_iter {
                m_step(
                    eq_map,
                    &base_counts,
                    &prev_counts,
                    &inv_eff_lens,
                    &mut curr_counts,
                );

                //std::mem::swap(&)
                for i in 0..curr_counts.len() {
                    if prev_counts[i] > presence_thresh {
                        let rd = (curr_counts[i] - prev_counts[i]) / prev_counts[i];
                        rel_diff = if rel_diff > rd { rel_diff } else { rd };
                    }
                }

                std::mem::swap(&mut prev_counts, &mut curr_counts);
                curr_counts.fill(0.0_f64);

                if rel_diff < converge_thresh {
                    break;
                }
                niter += 1;
                rel_diff = 0.0_f64;
            }

            prev_counts.iter_mut().for_each(|x| {
                if *x < presence_thresh {
                    *x = 0.0
                }
            });
            m_step(
                eq_map,
                &base_counts,
                &prev_counts,
                &inv_eff_lens,
                &mut curr_counts,
            );

            curr_counts
        })
        .collect()
}

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

    // Flatten: interleave samples from chains for better mixing in output
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

pub fn em_par<EqLabelT: EqLabel>(em_info: &EMInfo<EqLabelT>, nthreads: usize) -> Vec<f64> {
    let converge_thresh = em_info.convergence_thresh;
    let presence_thresh = em_info.presence_thresh;
    let eq_map = em_info.eq_map;
    let eff_lens = em_info.eff_lens;
    let inv_eff_lens = eff_lens
        .iter()
        .map(|x| {
            let y = 1.0_f64 / *x;
            if y.is_finite() { y } else { 0_f64 }
        })
        .collect::<Vec<f64>>();
    let max_iter = em_info.max_iter;
    let total_weight: f64 = eq_map.counts.iter().sum::<usize>() as f64;

    // init
    let avg = total_weight / (eff_lens.len() as f64);
    let mut prev_counts: Vec<AtomicF64> = vec![avg; eff_lens.len()]
        .iter()
        .map(|x| AtomicF64::new(*x))
        .collect();
    let mut curr_counts: Vec<AtomicF64> = vec![0.0f64; eff_lens.len()]
        .iter()
        .map(|x| AtomicF64::new(*x))
        .collect();

    let eq_iterates: Vec<(EqLabelT::LabelRefT<'_>, &usize)> =
        eq_map.iter_labels().zip(&eq_map.counts).collect();

    let mut rel_diff = 0.0_f64;
    let mut niter = 0_u32;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();

    pool.install(|| {
        while niter < max_iter {
            m_step_par::<EqLabelT>(
                &eq_iterates,
                &mut prev_counts,
                &inv_eff_lens,
                &mut curr_counts,
            );

            //std::mem::swap(&)
            for i in 0..curr_counts.len() {
                let pci = prev_counts[i].load(Ordering::Relaxed);
                if pci > presence_thresh {
                    let cci = curr_counts[i].load(Ordering::Relaxed);
                    let rd = (cci - pci) / pci;
                    rel_diff = if rel_diff > rd { rel_diff } else { rd };
                }
            }

            std::mem::swap(&mut prev_counts, &mut curr_counts);
            // zero out the vector
            curr_counts
                .par_iter()
                .for_each(|x| x.store(0.0f64, Ordering::Relaxed));

            if rel_diff < converge_thresh {
                break;
            }
            niter += 1;
            if niter % 100 == 0 {
                info!("iteration {}; rel diff {:.3}", niter, rel_diff);
            }
            rel_diff = 0.0_f64;
        }

        prev_counts.iter_mut().for_each(|x| {
            if x.load(Ordering::Relaxed) < presence_thresh {
                x.store(0.0, Ordering::Relaxed);
            }
        });
        m_step_par::<EqLabelT>(
            &eq_iterates,
            &mut prev_counts,
            &inv_eff_lens,
            &mut curr_counts,
        );
    });

    curr_counts
        .iter()
        .map(|x| x.load(Ordering::Relaxed))
        .collect::<Vec<f64>>()
}

pub fn em<EqLabelT: EqLabel>(em_info: &EMInfo<EqLabelT>) -> Vec<f64> {
    let converge_thresh = em_info.convergence_thresh;
    let presence_thresh = em_info.presence_thresh;
    let eq_map = em_info.eq_map;
    let eff_lens = em_info.eff_lens;
    let inv_eff_lens = eff_lens
        .iter()
        .map(|x| {
            let y = 1.0_f64 / *x;
            if y.is_finite() { y } else { 0_f64 }
        })
        .collect::<Vec<f64>>();
    let max_iter = em_info.max_iter;
    let total_weight: f64 = eq_map.counts.iter().sum::<usize>() as f64;

    // init
    let avg = total_weight / (eff_lens.len() as f64);
    let mut prev_counts = vec![avg; eff_lens.len()];
    let mut curr_counts = vec![0.0f64; eff_lens.len()];

    let mut rel_diff = 0.0_f64;
    let mut niter = 0_u32;

    while niter < max_iter {
        m_step(
            eq_map,
            &eq_map.counts,
            &prev_counts,
            &inv_eff_lens,
            &mut curr_counts,
        );

        //std::mem::swap(&)
        for i in 0..curr_counts.len() {
            if prev_counts[i] > presence_thresh {
                let rd = (curr_counts[i] - prev_counts[i]) / prev_counts[i];
                rel_diff = if rel_diff > rd { rel_diff } else { rd };
            }
        }

        std::mem::swap(&mut prev_counts, &mut curr_counts);
        curr_counts.fill(0.0_f64);

        if rel_diff < converge_thresh {
            break;
        }
        niter += 1;
        if niter % 100 == 0 {
            info!("iteration {}; rel diff {:.3}", niter, rel_diff);
        }
        rel_diff = 0.0_f64;
    }

    prev_counts.iter_mut().for_each(|x| {
        if *x < presence_thresh {
            *x = 0.0
        }
    });
    m_step(
        eq_map,
        &eq_map.counts,
        &prev_counts,
        &inv_eff_lens,
        &mut curr_counts,
    );

    curr_counts
}
