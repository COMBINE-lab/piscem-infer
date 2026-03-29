use atomic_float::AtomicF64;
use rand::prelude::*;
use rand::rng;
use rand_distr::weighted::WeightedAliasIndex;
use rayon::prelude::*;
use rayon::ThreadPool;
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
fn install_in_pool<R: Send>(pool: Option<&ThreadPool>, f: impl FnOnce() -> R + Send) -> R {
    if let Some(pool) = pool {
        pool.install(f)
    } else {
        f()
    }
}

#[inline]
fn m_step_par<EqLabelT: EqLabel>(
    eq_iterates: &[(EqLabelT::LabelRefT<'_>, &usize)],
    prev_count: &mut [AtomicF64],
    inv_eff_lens: &[f64],
    curr_counts: &mut [AtomicF64],
) {
    eq_iterates.par_iter().for_each_with(
        (&curr_counts, Vec::with_capacity(64)),
        |(curr_counts, weights), (k, v)| {
            let count = **v as f64;
            weights.clear();
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

#[derive(Clone, Copy, Debug)]
pub struct SquaremOptions {
    pub min_step: f64,
    pub max_step: f64,
    pub burn_in_steps: u32,
    pub acceleration_interval: u32,
}

impl Default for SquaremOptions {
    fn default() -> Self {
        Self {
            min_step: 1e-4,
            max_step: 10.0,
            burn_in_steps: 25,
            acceleration_interval: 4,
        }
    }
}

fn initial_counts(eff_lens: &[f64], total_weight: f64, init_counts: Option<&[f64]>) -> Vec<f64> {
    if let Some(init) = init_counts {
        let counts = init
            .iter()
            .zip(eff_lens.iter())
            .map(|(&c, &el)| if el > 0.0 && c.is_finite() && c > 0.0 { c } else { 0.0 })
            .collect::<Vec<f64>>();
        if counts.iter().any(|&x| x > 0.0) {
            return counts;
        }
    }

    // Uniform over active transcripts (eff_len > 0); zero for masked ones.
    let n_active = eff_lens.iter().filter(|&&el| el > 0.0).count() as f64;
    if n_active > 0.0 {
        let avg = total_weight / n_active;
        eff_lens
            .iter()
            .map(|&el| if el > 0.0 { avg } else { 0.0 })
            .collect()
    } else {
        vec![0.0; eff_lens.len()]
    }
}

#[inline]
fn compute_inv_eff_lens(eff_lens: &[f64]) -> Vec<f64> {
    eff_lens
        .iter()
        .map(|x| {
            let y = 1.0_f64 / *x;
            if y.is_finite() { y } else { 0_f64 }
        })
        .collect::<Vec<f64>>()
}

#[inline]
fn compute_rel_diff(prev_counts: &[f64], curr_counts: &[f64], presence_thresh: f64) -> f64 {
    let mut sum_abs_rel = 0.0_f64;
    let mut n = 0_u64;
    for (&prev, &curr) in prev_counts.iter().zip(curr_counts.iter()) {
        if prev > presence_thresh {
            sum_abs_rel += ((curr - prev) / prev).abs();
            n += 1;
        }
    }
    if n > 0 { sum_abs_rel / n as f64 } else { 0.0 }
}

#[inline]
fn compute_rel_diff_atomic(
    prev_counts: &[AtomicF64],
    curr_counts: &[AtomicF64],
    presence_thresh: f64,
) -> f64 {
    let mut sum_abs_rel = 0.0_f64;
    let mut n = 0_u64;
    for (prev, curr) in prev_counts.iter().zip(curr_counts.iter()) {
        let p = prev.load(Ordering::Relaxed);
        let c = curr.load(Ordering::Relaxed);
        if p > presence_thresh {
            sum_abs_rel += ((c - p) / p).abs();
            n += 1;
        }
    }
    if n > 0 { sum_abs_rel / n as f64 } else { 0.0 }
}

#[inline]
fn project_counts(counts: &mut [f64], eff_lens: &[f64], total_weight: f64) {
    let mut sum = 0.0f64;
    for (count, &el) in counts.iter_mut().zip(eff_lens.iter()) {
        if !count.is_finite() || *count < 0.0 || el <= 0.0 {
            *count = 0.0;
        }
        sum += *count;
    }
    if sum > 0.0 {
        let scale = total_weight / sum;
        for count in counts.iter_mut() {
            *count *= scale;
        }
    }
}

#[inline]
fn squarem_alpha(x0: &[f64], x1: &[f64], x2: &[f64], opts: SquaremOptions) -> Option<f64> {
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
    let mut alpha = -(rr / vv).sqrt();
    if !alpha.is_finite() {
        return None;
    }
    alpha = alpha.clamp(-opts.max_step, -opts.min_step);
    Some(alpha)
}

#[inline]
fn em_step_plain<EqLabelT: EqLabel>(
    eq_map: &PackedEqMap<EqLabelT>,
    inv_eff_lens: &[f64],
    prev_counts: &[f64],
    out: &mut [f64],
) {
    out.fill(0.0);
    m_step(eq_map, &eq_map.counts, prev_counts, inv_eff_lens, out);
}

#[inline]
fn m_step_par_from_slice<EqLabelT: EqLabel>(
    eq_iterates: &[(EqLabelT::LabelRefT<'_>, &usize)],
    prev_count: &[f64],
    inv_eff_lens: &[f64],
    curr_counts: &mut [AtomicF64],
) {
    eq_iterates.par_iter().for_each_with(
        (&curr_counts, Vec::with_capacity(64)),
        |(curr_counts, weights), (k, v)| {
            let count = **v as f64;
            weights.clear();
            let mut denom = 0.0_f64;
            for (e, cond_prob) in k.target_labels().iter().zip(k.target_probs()) {
                let w = cond_prob * prev_count[*e as usize] * inv_eff_lens[*e as usize];
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
fn em_step_plain_par_in_pool<EqLabelT: EqLabel>(
    eq_iterates: &[(EqLabelT::LabelRefT<'_>, &usize)],
    inv_eff_lens: &[f64],
    prev_counts: &[f64],
    curr_counts: &mut [AtomicF64],
    pool: &ThreadPool,
) -> Vec<f64> {
    install_in_pool(Some(pool), || {
        curr_counts
            .par_iter()
            .for_each(|x| x.store(0.0f64, Ordering::Relaxed));
        m_step_par_from_slice::<EqLabelT>(eq_iterates, prev_counts, inv_eff_lens, curr_counts);
    });
    curr_counts
        .iter()
        .map(|x| x.load(Ordering::Relaxed))
        .collect::<Vec<f64>>()
}

#[inline]
fn should_try_squarem(step_count: u32, rel_diff: f64, opts: SquaremOptions) -> bool {
    step_count >= opts.burn_in_steps
        && (step_count - opts.burn_in_steps).is_multiple_of(opts.acceleration_interval)
        && rel_diff > 10.0 * f64::EPSILON
}


pub fn do_bootstrap<EqLabelT: EqLabel>(
    em_info: &EMInfo<EqLabelT>,
    num_boot: usize,
) -> Vec<Vec<f64>> {
    do_bootstrap_in_pool(em_info, num_boot, None)
}

pub fn do_bootstrap_with_pool<EqLabelT: EqLabel>(
    em_info: &EMInfo<EqLabelT>,
    num_boot: usize,
    pool: &ThreadPool,
) -> Vec<Vec<f64>> {
    do_bootstrap_in_pool(em_info, num_boot, Some(pool))
}

fn do_bootstrap_in_pool<EqLabelT: EqLabel>(
    em_info: &EMInfo<EqLabelT>,
    num_boot: usize,
    pool: Option<&ThreadPool>,
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

    install_in_pool(pool, || {
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

                    rel_diff = compute_rel_diff(&prev_counts, &curr_counts, presence_thresh);

                    std::mem::swap(&mut prev_counts, &mut curr_counts);
                    curr_counts.fill(0.0_f64);

                    if rel_diff < converge_thresh {
                        break;
                    }
                    niter += 1;
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
    })
}

pub fn em_par<EqLabelT: EqLabel>(em_info: &EMInfo<EqLabelT>, nthreads: usize) -> Vec<f64> {
    em_par_init(em_info, None, nthreads)
}

pub fn em_par_with_pool<EqLabelT: EqLabel>(
    em_info: &EMInfo<EqLabelT>,
    pool: &ThreadPool,
) -> Vec<f64> {
    em_par_with_pool_init(em_info, None, pool)
}

pub fn em_par_init<EqLabelT: EqLabel>(
    em_info: &EMInfo<EqLabelT>,
    init_counts: Option<&[f64]>,
    nthreads: usize,
) -> Vec<f64> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();
    em_par_with_pool_init(em_info, init_counts, &pool)
}

pub fn em_par_with_pool_init<EqLabelT: EqLabel>(
    em_info: &EMInfo<EqLabelT>,
    init_counts: Option<&[f64]>,
    pool: &ThreadPool,
) -> Vec<f64> {
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

    let init = initial_counts(eff_lens, total_weight, init_counts);
    let mut prev_counts: Vec<AtomicF64> = init.iter().map(|x| AtomicF64::new(*x)).collect();
    let mut curr_counts: Vec<AtomicF64> = vec![0.0f64; eff_lens.len()]
        .iter()
        .map(|x| AtomicF64::new(*x))
        .collect();
    let eq_iterates: Vec<(EqLabelT::LabelRefT<'_>, &usize)> =
        eq_map.iter_labels().zip(&eq_map.counts).collect();

    let mut rel_diff = 0.0_f64;
    let mut last_rel_diff = f64::INFINITY;
    let mut niter = 0_u32;

    install_in_pool(Some(pool), || {
        while niter < max_iter {
            m_step_par::<EqLabelT>(
                &eq_iterates,
                &mut prev_counts,
                &inv_eff_lens,
                &mut curr_counts,
            );

            rel_diff = compute_rel_diff_atomic(&prev_counts, &curr_counts, presence_thresh);
            last_rel_diff = rel_diff;

            std::mem::swap(&mut prev_counts, &mut curr_counts);
            curr_counts
                .par_iter()
                .for_each(|x| x.store(0.0f64, Ordering::Relaxed));

            if rel_diff < converge_thresh {
                break;
            }
            niter += 1;
            if niter.is_multiple_of(100) {
                info!("iteration {}; rel diff {:.3}", niter, rel_diff);
            }
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

    let final_counts = curr_counts
        .iter()
        .map(|x| x.load(Ordering::Relaxed))
        .collect::<Vec<f64>>();
    info!(
        "EM stats: em_steps={} final_rel_diff={:.6}",
        niter, last_rel_diff
    );
    final_counts
}

pub fn squarem_em<EqLabelT: EqLabel>(em_info: &EMInfo<EqLabelT>) -> Vec<f64> {
    squarem_em_init(em_info, None)
}

pub fn squarem_em_init<EqLabelT: EqLabel>(
    em_info: &EMInfo<EqLabelT>,
    init_counts: Option<&[f64]>,
) -> Vec<f64> {
    let opts = SquaremOptions::default();
    let presence_thresh = em_info.presence_thresh;
    let eq_map = em_info.eq_map;
    let eff_lens = em_info.eff_lens;
    let inv_eff_lens = compute_inv_eff_lens(eff_lens);
    let max_iter = em_info.max_iter;
    let total_weight: f64 = eq_map.counts.iter().sum::<usize>() as f64;

    let mut x0 = initial_counts(eff_lens, total_weight, init_counts);
    project_counts(&mut x0, eff_lens, total_weight);

    let mut x1 = vec![0.0f64; eff_lens.len()];
    let mut x2 = vec![0.0f64; eff_lens.len()];
    let mut x_sq = vec![0.0f64; eff_lens.len()];
    let mut x_next = vec![0.0f64; eff_lens.len()];
    let mut em_steps = 0_u32;
    let mut last_rel_diff = f64::INFINITY;
    let mut accel_attempts = 0_u32;
    let mut accel_accepts = 0_u32;

    while em_steps < max_iter {
        em_step_plain(eq_map, &inv_eff_lens, &x0, &mut x1);
        em_steps += 1;
        let rel1 = compute_rel_diff(&x0, &x1, presence_thresh);
        last_rel_diff = rel1;
        if rel1 < em_info.convergence_thresh || em_steps >= max_iter {
            x0 = x1.clone();
            break;
        }

        if !should_try_squarem(em_steps, last_rel_diff, opts) || em_steps >= max_iter {
            x0.clone_from_slice(&x1);
            continue;
        }

        accel_attempts += 1;
        em_step_plain(eq_map, &inv_eff_lens, &x1, &mut x2);
        em_steps += 1;

        let ordinary_rel = compute_rel_diff(&x1, &x2, presence_thresh);
        let candidate = if let Some(alpha) = squarem_alpha(&x0, &x1, &x2, opts) {
            for (((sq, &a), &b), &c) in x_sq.iter_mut().zip(x0.iter()).zip(x1.iter()).zip(x2.iter())
            {
                let r = b - a;
                let v = c - (2.0 * b) + a;
                *sq = a - (2.0 * alpha * r) + (alpha * alpha * v);
            }
            project_counts(&mut x_sq, eff_lens, total_weight);
            if em_steps < max_iter {
                em_step_plain(eq_map, &inv_eff_lens, &x_sq, &mut x_next);
                em_steps += 1;
                let candidate_rel = compute_rel_diff(&x_sq, &x_next, presence_thresh);
                if x_next.iter().all(|x| x.is_finite() && *x >= 0.0)
                    && candidate_rel < ordinary_rel
                {
                    accel_accepts += 1;
                    &x_next
                } else {
                    &x2
                }
            } else {
                &x2
            }
        } else {
            &x2
        };

        let rel_diff = compute_rel_diff(&x0, candidate, presence_thresh);
        x0.clone_from_slice(candidate);
        if rel_diff < em_info.convergence_thresh {
            break;
        }
    }

    x0.iter_mut().for_each(|x| {
        if *x < presence_thresh {
            *x = 0.0;
        }
    });
    em_step_plain(eq_map, &inv_eff_lens, &x0, &mut x1);
    info!(
        "SQUAREM stats: em_steps={} accel_attempts={} accel_accepts={} final_rel_diff={:.6}",
        em_steps, accel_attempts, accel_accepts, last_rel_diff
    );
    x1
}

pub fn squarem_em_par<EqLabelT: EqLabel>(em_info: &EMInfo<EqLabelT>, nthreads: usize) -> Vec<f64> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();
    squarem_em_par_with_pool_init(em_info, None, &pool)
}

pub fn squarem_em_par_with_pool<EqLabelT: EqLabel>(
    em_info: &EMInfo<EqLabelT>,
    pool: &ThreadPool,
) -> Vec<f64> {
    squarem_em_par_with_pool_init(em_info, None, pool)
}

pub fn squarem_em_par_with_pool_init<EqLabelT: EqLabel>(
    em_info: &EMInfo<EqLabelT>,
    init_counts: Option<&[f64]>,
    pool: &ThreadPool,
) -> Vec<f64> {
    let opts = SquaremOptions::default();
    let presence_thresh = em_info.presence_thresh;
    let eq_map = em_info.eq_map;
    let eff_lens = em_info.eff_lens;
    let inv_eff_lens = compute_inv_eff_lens(eff_lens);
    let max_iter = em_info.max_iter;
    let total_weight: f64 = eq_map.counts.iter().sum::<usize>() as f64;
    let eq_iterates: Vec<(EqLabelT::LabelRefT<'_>, &usize)> =
        eq_map.iter_labels().zip(&eq_map.counts).collect();
    let mut curr_counts: Vec<AtomicF64> = vec![0.0f64; eff_lens.len()]
        .iter()
        .map(|x| AtomicF64::new(*x))
        .collect();

    let mut x0 = initial_counts(eff_lens, total_weight, init_counts);
    project_counts(&mut x0, eff_lens, total_weight);
    let mut em_steps = 0_u32;
    let mut last_rel_diff = f64::INFINITY;
    let mut accel_attempts = 0_u32;
    let mut accel_accepts = 0_u32;

    while em_steps < max_iter {
        let x1 = em_step_plain_par_in_pool::<EqLabelT>(
            &eq_iterates,
            &inv_eff_lens,
            &x0,
            &mut curr_counts,
            pool,
        );
        em_steps += 1;
        let rel1 = compute_rel_diff(&x0, &x1, presence_thresh);
        last_rel_diff = rel1;
        if rel1 < em_info.convergence_thresh || em_steps >= max_iter {
            x0 = x1;
            break;
        }

        if !should_try_squarem(em_steps, last_rel_diff, opts) || em_steps >= max_iter {
            x0 = x1;
            continue;
        }

        accel_attempts += 1;
        let x2 = em_step_plain_par_in_pool::<EqLabelT>(
            &eq_iterates,
            &inv_eff_lens,
            &x1,
            &mut curr_counts,
            pool,
        );
        em_steps += 1;

        let ordinary_rel = compute_rel_diff(&x1, &x2, presence_thresh);
        let candidate = if let Some(alpha) = squarem_alpha(&x0, &x1, &x2, opts) {
            let mut x_sq = vec![0.0f64; eff_lens.len()];
            for (((sq, &a), &b), &c) in x_sq.iter_mut().zip(x0.iter()).zip(x1.iter()).zip(x2.iter())
            {
                let r = b - a;
                let v = c - (2.0 * b) + a;
                *sq = a - (2.0 * alpha * r) + (alpha * alpha * v);
            }
            project_counts(&mut x_sq, eff_lens, total_weight);
            if em_steps < max_iter {
                let x_next = em_step_plain_par_in_pool::<EqLabelT>(
                    &eq_iterates,
                    &inv_eff_lens,
                    &x_sq,
                    &mut curr_counts,
                    pool,
                );
                em_steps += 1;
                let candidate_rel = compute_rel_diff(&x_sq, &x_next, presence_thresh);
                if x_next.iter().all(|x| x.is_finite() && *x >= 0.0)
                    && candidate_rel < ordinary_rel
                {
                    accel_accepts += 1;
                    x_next
                } else {
                    x2
                }
            } else {
                x2
            }
        } else {
            x2
        };

        let rel_diff = compute_rel_diff(&x0, &candidate, presence_thresh);
        x0 = candidate;
        if rel_diff < em_info.convergence_thresh {
            break;
        }
    }

    x0.iter_mut().for_each(|x| {
        if *x < presence_thresh {
            *x = 0.0;
        }
    });
    let final_counts =
        em_step_plain_par_in_pool::<EqLabelT>(&eq_iterates, &inv_eff_lens, &x0, &mut curr_counts, pool);
    info!(
        "SQUAREM stats: em_steps={} accel_attempts={} accel_accepts={} final_rel_diff={:.6}",
        em_steps, accel_attempts, accel_accepts, last_rel_diff
    );
    final_counts
}

/// Run EM with pseudo-count regularization from a hierarchical prior.
///
/// After each M-step, adds `alpha[t]` pseudo-counts to `curr_counts[t]`.
/// This causes the prior to influence read assignments in the next E-step,
/// unlike post-hoc L-BFGS which can only reweight the point estimate.
///
/// Returns the final estimated counts (not normalized).
pub fn em_penalized<EqLabelT: EqLabel>(
    em_info: &EMInfo<EqLabelT>,
    alpha: &[f64],
) -> Vec<f64> {
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
    let mut last_rel_diff = f64::INFINITY;
    let mut niter = 0_u32;

    while niter < max_iter {
        m_step(
            eq_map,
            &eq_map.counts,
            &prev_counts,
            &inv_eff_lens,
            &mut curr_counts,
        );

        // Add pseudo-counts from hierarchical prior
        for (c, &a) in curr_counts.iter_mut().zip(alpha.iter()) {
            *c += a;
        }

        rel_diff = compute_rel_diff(&prev_counts, &curr_counts, presence_thresh);
        last_rel_diff = rel_diff;

        std::mem::swap(&mut prev_counts, &mut curr_counts);
        curr_counts.fill(0.0_f64);

        if rel_diff < converge_thresh {
            break;
        }
        niter += 1;
        if niter.is_multiple_of(100) {
            info!("iteration {}; rel diff {:.3}", niter, rel_diff);
        }
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

    // Add pseudo-counts to the final step too
    for (c, &a) in curr_counts.iter_mut().zip(alpha.iter()) {
        *c += a;
    }

    curr_counts
}

/// Parallel version of penalized EM with pseudo-count regularization.
pub fn em_penalized_par<EqLabelT: EqLabel>(
    em_info: &EMInfo<EqLabelT>,
    alpha: &[f64],
    nthreads: usize,
) -> Vec<f64> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();
    em_penalized_par_with_pool(em_info, alpha, &pool)
}

pub fn em_penalized_par_with_pool<EqLabelT: EqLabel>(
    em_info: &EMInfo<EqLabelT>,
    alpha: &[f64],
    pool: &ThreadPool,
) -> Vec<f64> {
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
    let mut last_rel_diff = f64::INFINITY;
    let mut niter = 0_u32;

    install_in_pool(Some(pool), || {
        while niter < max_iter {
            m_step_par::<EqLabelT>(
                &eq_iterates,
                &mut prev_counts,
                &inv_eff_lens,
                &mut curr_counts,
            );

            // Add pseudo-counts from hierarchical prior
            for (c, &a) in curr_counts.iter().zip(alpha.iter()) {
                c.fetch_add(a, Ordering::AcqRel);
            }

            rel_diff = compute_rel_diff_atomic(&prev_counts, &curr_counts, presence_thresh);

            std::mem::swap(&mut prev_counts, &mut curr_counts);
            curr_counts
                .par_iter()
                .for_each(|x| x.store(0.0f64, Ordering::Relaxed));

            if rel_diff < converge_thresh {
                break;
            }
            niter += 1;
            if niter.is_multiple_of(100) {
                info!("iteration {}; rel diff {:.3}", niter, rel_diff);
            }
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

        // Add pseudo-counts to the final step too
        for (c, &a) in curr_counts.iter().zip(alpha.iter()) {
            c.fetch_add(a, Ordering::AcqRel);
        }
    });

    curr_counts
        .iter()
        .map(|x| x.load(Ordering::Relaxed))
        .collect::<Vec<f64>>()
}

pub fn em<EqLabelT: EqLabel>(em_info: &EMInfo<EqLabelT>) -> Vec<f64> {
    em_init(em_info, None)
}

pub fn em_init<EqLabelT: EqLabel>(
    em_info: &EMInfo<EqLabelT>,
    init_counts: Option<&[f64]>,
) -> Vec<f64> {
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

    let mut prev_counts = initial_counts(eff_lens, total_weight, init_counts);
    let mut curr_counts = vec![0.0f64; eff_lens.len()];

    let mut rel_diff = 0.0_f64;
    let mut last_rel_diff = f64::INFINITY;
    let mut niter = 0_u32;

    while niter < max_iter {
        m_step(
            eq_map,
            &eq_map.counts,
            &prev_counts,
            &inv_eff_lens,
            &mut curr_counts,
        );

        rel_diff = compute_rel_diff(&prev_counts, &curr_counts, presence_thresh);

        last_rel_diff = rel_diff;
        std::mem::swap(&mut prev_counts, &mut curr_counts);
        curr_counts.fill(0.0_f64);

        if rel_diff < converge_thresh {
            break;
        }
        niter += 1;
        if niter.is_multiple_of(100) {
            info!("iteration {}; rel diff {:.3}", niter, rel_diff);
        }
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

    info!(
        "EM stats: em_steps={} final_rel_diff={:.6}",
        niter, last_rel_diff
    );
    curr_counts
}

/// Compute per-transcript coverage profile from converged EM counts.
/// Returns a flattened `n_targets × n_pos_bins` array where
/// `profile[t * n_pos_bins + b]` = assigned reads to transcript t from pos_bin b.
pub fn compute_coverage_profile<EqLabelT: EqLabel>(
    eq_map: &PackedEqMap<EqLabelT>,
    em_counts: &[f64],
    eff_lens: &[f64],
    n_targets: usize,
    n_pos_bins: usize,
) -> Vec<f64> {
    let inv_eff_lens: Vec<f64> = eff_lens
        .iter()
        .map(|x| { let y = 1.0 / *x; if y.is_finite() { y } else { 0.0 } })
        .collect();

    let mut profile = vec![0.0f64; n_targets * n_pos_bins];

    for (label, &count) in eq_map.iter_labels().zip(eq_map.counts.iter()) {
        let ec_count = count as f64;
        if ec_count == 0.0 { continue; }

        let pos_bins = label.target_pos_bins();

        // Compute posterior shares (same as M-step)
        let mut denom = 0.0f64;
        let mut weights: Vec<f64> = Vec::new();
        for (tid, cond_prob) in label.target_labels().iter().zip(label.target_probs()) {
            let w = cond_prob * em_counts[*tid as usize] * inv_eff_lens[*tid as usize];
            weights.push(w);
            denom += w;
        }

        if denom > 1e-8 {
            for (i, (tid, &w)) in label.target_labels().iter().zip(weights.iter()).enumerate() {
                let t = *tid as usize;
                let assigned = ec_count * w / denom;
                let b = pos_bins.map_or(0, |pb| pb[i] as usize).min(n_pos_bins - 1);
                profile[t * n_pos_bins + b] += assigned;
            }
        }
    }

    profile
}

/// Compute coverage-consistency weights from a coverage profile.
/// For each transcript, bins with below-average coverage get higher weight,
/// bins with above-average coverage get lower weight.
/// Returns flattened `n_targets × n_pos_bins` weights, normalized per transcript.
pub fn coverage_weights_from_profile(
    profile: &[f64],
    n_targets: usize,
    n_pos_bins: usize,
    epsilon: f64,
) -> Vec<f64> {
    let mut weights = vec![1.0f64; n_targets * n_pos_bins];

    for t in 0..n_targets {
        let base = t * n_pos_bins;
        let total: f64 = profile[base..base + n_pos_bins].iter().sum();
        if total < 1.0 {
            // Transcript has negligible coverage — uniform weights
            continue;
        }
        let expected_per_bin = total / n_pos_bins as f64;

        for b in 0..n_pos_bins {
            // Inverse weighting: under-covered bins get higher weight
            weights[base + b] = expected_per_bin / (profile[base + b] + epsilon);
        }

        // Normalize weights to sum to n_pos_bins (preserves total mass)
        let w_sum: f64 = weights[base..base + n_pos_bins].iter().sum();
        if w_sum > 0.0 {
            let scale = n_pos_bins as f64 / w_sum;
            for b in 0..n_pos_bins {
                weights[base + b] *= scale;
            }
        }
    }

    weights
}

/// Run a coverage-aware EM: standard EM followed by coverage profile
/// estimation and a second EM pass with coverage weights.
pub fn em_with_coverage<EqLabelT: EqLabel>(
    em_info: &EMInfo<EqLabelT>,
    init_counts: Option<&[f64]>,
    n_pos_bins: usize,
    n_coverage_rounds: usize,
    coverage_epsilon: f64,
) -> Vec<f64> {
    let n_targets = em_info.eff_lens.len();

    // Initial EM (SQUAREM-accelerated)
    let mut counts = squarem_em_init(em_info, init_counts);

    if n_pos_bins <= 1 || n_coverage_rounds == 0 {
        return counts;
    }

    for round in 0..n_coverage_rounds {
        // Compute coverage profile from current estimates
        let profile = compute_coverage_profile(
            em_info.eq_map, &counts, em_info.eff_lens, n_targets, n_pos_bins,
        );

        // Compute coverage weights
        let cov_weights = coverage_weights_from_profile(&profile, n_targets, n_pos_bins, coverage_epsilon);

        // Run coverage-weighted EM
        counts = em_coverage_weighted(em_info, Some(&counts), &cov_weights, n_pos_bins);

        info!("Coverage EM round {}/{} complete", round + 1, n_coverage_rounds);
    }

    counts
}

/// One coverage-weighted M-step: redistribute reads using coverage weights.
#[inline]
fn m_step_coverage_weighted<EqLabelT: EqLabel>(
    eq_map: &PackedEqMap<EqLabelT>,
    prev_counts: &[f64],
    inv_eff_lens: &[f64],
    cov_weights: &[f64],
    n_pos_bins: usize,
    out: &mut [f64],
) {
    out.fill(0.0);
    let mut weights: Vec<f64> = Vec::with_capacity(64);
    for (label, &count) in eq_map.iter_labels().zip(eq_map.counts.iter()) {
        let ec_count = count as f64;
        let pos_bins = label.target_pos_bins();
        let mut denom = 0.0f64;
        for (i, (tid, cond_prob)) in label.target_labels().iter()
            .zip(label.target_probs())
            .enumerate()
        {
            let t = *tid as usize;
            let b = pos_bins.map_or(0, |pb| pb[i] as usize).min(n_pos_bins - 1);
            let cw = cov_weights[t * n_pos_bins + b];
            let w = cond_prob * cw * prev_counts[t] * inv_eff_lens[t];
            weights.push(w);
            denom += w;
        }
        if denom > 1e-8 {
            let count_over_denom = ec_count / denom;
            for (tid, &w) in label.target_labels().iter().zip(weights.iter()) {
                out[*tid as usize] += count_over_denom * w;
            }
        }
        weights.clear();
    }
}

/// SQUAREM-accelerated EM with coverage weights.
fn em_coverage_weighted<EqLabelT: EqLabel>(
    em_info: &EMInfo<EqLabelT>,
    init_counts: Option<&[f64]>,
    cov_weights: &[f64],
    n_pos_bins: usize,
) -> Vec<f64> {
    let opts = SquaremOptions::default();
    let presence_thresh = em_info.presence_thresh;
    let eq_map = em_info.eq_map;
    let eff_lens = em_info.eff_lens;
    let inv_eff_lens = compute_inv_eff_lens(eff_lens);
    let max_iter = em_info.max_iter;
    let total_weight: f64 = eq_map.counts.iter().sum::<usize>() as f64;

    let mut x0 = initial_counts(eff_lens, total_weight, init_counts);
    project_counts(&mut x0, eff_lens, total_weight);

    let mut x1 = vec![0.0f64; eff_lens.len()];
    let mut x2 = vec![0.0f64; eff_lens.len()];
    let mut x_sq = vec![0.0f64; eff_lens.len()];
    let mut x_next = vec![0.0f64; eff_lens.len()];
    let mut em_steps = 0_u32;
    let mut last_rel_diff = f64::INFINITY;
    let mut accel_attempts = 0_u32;
    let mut accel_accepts = 0_u32;

    while em_steps < max_iter {
        m_step_coverage_weighted(eq_map, &x0, &inv_eff_lens, cov_weights, n_pos_bins, &mut x1);
        em_steps += 1;
        let rel1 = compute_rel_diff(&x0, &x1, presence_thresh);
        last_rel_diff = rel1;
        if rel1 < em_info.convergence_thresh || em_steps >= max_iter {
            x0 = x1.clone();
            break;
        }

        if !should_try_squarem(em_steps, last_rel_diff, opts) || em_steps >= max_iter {
            x0.clone_from_slice(&x1);
            continue;
        }

        accel_attempts += 1;
        m_step_coverage_weighted(eq_map, &x1, &inv_eff_lens, cov_weights, n_pos_bins, &mut x2);
        em_steps += 1;

        let ordinary_rel = compute_rel_diff(&x1, &x2, presence_thresh);
        let candidate = if let Some(alpha) = squarem_alpha(&x0, &x1, &x2, opts) {
            for (((sq, &a), &b), &c) in x_sq.iter_mut().zip(x0.iter()).zip(x1.iter()).zip(x2.iter()) {
                let r = b - a;
                let v = c - (2.0 * b) + a;
                *sq = a - (2.0 * alpha * r) + (alpha * alpha * v);
            }
            project_counts(&mut x_sq, eff_lens, total_weight);
            if em_steps < max_iter {
                m_step_coverage_weighted(eq_map, &x_sq, &inv_eff_lens, cov_weights, n_pos_bins, &mut x_next);
                em_steps += 1;
                let candidate_rel = compute_rel_diff(&x_sq, &x_next, presence_thresh);
                if x_next.iter().all(|x| x.is_finite() && *x >= 0.0) && candidate_rel < ordinary_rel {
                    accel_accepts += 1;
                    &x_next
                } else {
                    &x2
                }
            } else {
                &x2
            }
        } else {
            &x2
        };

        let rel_diff = compute_rel_diff(&x0, candidate, presence_thresh);
        x0.clone_from_slice(candidate);
        if rel_diff < em_info.convergence_thresh {
            break;
        }
    }

    x0.iter_mut().for_each(|x| { if *x < presence_thresh { *x = 0.0 } });
    m_step_coverage_weighted(eq_map, &x0, &inv_eff_lens, cov_weights, n_pos_bins, &mut x1);

    info!(
        "Coverage SQUAREM: em_steps={} accel_attempts={} accel_accepts={} final_rel_diff={:.6}",
        em_steps, accel_attempts, accel_accepts, last_rel_diff
    );
    x1
}

