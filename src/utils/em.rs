//! Effective-length computation and the bridge from piscem-infer's
//! [`PackedEqMap`] to `salmon-infer`'s [`PackedEqClasses`].
//!
//! The numeric core (EM / VBEM, SQUAREM / DAAREM acceleration, bootstrap and
//! Gibbs sampling) lives in the `salmon-infer` crate. Its parallel M-step
//! partitions the work over a *data-derived* shard plan, so the point estimate
//! is bit-identical regardless of the thread count (the previous hand-rolled
//! shared-atomic M-step was not). This module only builds the flat CSR input it
//! consumes and maps piscem-infer's options onto its [`EmOptions`].

use salmon_infer::{EmAccel, EmOptions, EmResult, PackedEqClasses};
use tracing::{info, warn};

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

/// Build `salmon-infer`'s flat CSR [`PackedEqClasses`] from piscem-infer's
/// [`PackedEqMap`] and the per-target effective lengths.
///
/// For every equivalence class `i` (in `eq_map` order) and every target `j`
/// in that class:
/// * `labels[..]`   = the target id (orientation is already stripped by
///   [`PackedEqMap`]),
/// * `weights[..]`  = the conditional probability of the fragment given the
///   target (`1.0` for basic equivalence classes; the bin-center probability
///   for range-factorized ones) — read by the Gibbs sampler,
/// * `combined[..]` = `weights / eff_len[target]` — read by the EM and the
///   bootstrap. This is exactly the `cond_prob * inv_eff_len` factor of the
///   previous piscem-infer M-step, so the model is unchanged.
///
/// `keep_weights` controls whether the raw `weights` array is materialized; only
/// Gibbs sampling reads it.
pub fn build_packed_eq_classes<EqLabelT: EqLabel>(
    eq_map: &PackedEqMap<EqLabelT>,
    eff_lens: &[f64],
    keep_weights: bool,
) -> PackedEqClasses {
    let num_classes = eq_map.len();
    // `eq_labels` may also carry probability-bin ids (range factorized), so
    // count the incidences (class, target) to size the arrays exactly.
    let cap: usize = eq_map
        .iter_labels()
        .map(|label| label.target_labels().len())
        .sum();
    let mut labels = Vec::<u32>::with_capacity(cap);
    let mut combined = Vec::<f64>::with_capacity(cap);
    let mut weights = Vec::<f64>::with_capacity(if keep_weights { cap } else { 0 });
    let mut starts = Vec::<u64>::with_capacity(num_classes + 1);
    let mut counts = Vec::<u64>::with_capacity(num_classes);

    let inv_eff_lens: Vec<f64> = eff_lens
        .iter()
        .map(|x| {
            let y = 1.0_f64 / *x;
            if y.is_finite() { y } else { 0_f64 }
        })
        .collect();

    starts.push(0u64);
    let mut total_count = 0u64;
    for (label, count) in eq_map.iter_labels().zip(eq_map.counts.iter()) {
        // A class with no targets carries no assignable evidence, and the EM
        // kernels index a class's first target unconditionally.
        if label.target_labels().is_empty() {
            continue;
        }
        for (tid, cond_prob) in label.target_labels().iter().zip(label.target_probs()) {
            labels.push(*tid);
            combined.push(cond_prob * inv_eff_lens[*tid as usize]);
            if keep_weights {
                weights.push(cond_prob);
            }
        }
        starts.push(labels.len() as u64);
        counts.push(*count as u64);
        total_count += *count as u64;
    }

    PackedEqClasses {
        labels,
        starts,
        combined,
        weights,
        counts,
        num_txps: eff_lens.len(),
        total_count,
    }
}

/// Map piscem-infer's EM options onto `salmon-infer`'s [`EmOptions`].
///
/// Convergence follows salmon: `alpha_check_cutoff` is the abundance at or
/// below which a target is excluded from the relative-difference convergence
/// check, and `presence_thresh` (salmon's `min_alpha`) is the threshold below
/// which a target is truncated to zero on output.
pub fn em_options(
    max_iter: u32,
    convergence_thresh: f64,
    alpha_check_cutoff: f64,
    presence_thresh: f64,
    accel: EmAccel,
) -> EmOptions {
    EmOptions {
        max_iter,
        rel_diff_tol: convergence_thresh,
        alpha_check_cutoff,
        min_alpha: presence_thresh,
        accel,
        ..EmOptions::default()
    }
}

/// Run the EM to convergence over `packed`, returning the per-target
/// fractional counts (`alphas`) along with the number of M-steps run and
/// whether the convergence criterion was met before `--max-iter`.
///
/// The sharded M-step is always used (over the current rayon pool). Its work
/// partition and summation order are derived from the data alone, so the
/// result is bit-identical for any thread count — including a single
/// thread, where the shards simply run one after another. salmon-infer's
/// purely sequential M-step accumulates in a different order and is
/// therefore *not* used even when only one thread is available.
pub fn run_em(packed: &PackedEqClasses, opts: &EmOptions) -> EmResult {
    let res = salmon_infer::optimize_packed(packed, opts, true);
    info!(
        "EM finished after {} M-steps (converged = {}, accel = {:?})",
        res.iters, res.converged, opts.accel
    );
    if !res.converged {
        warn!(
            "EM stopped at --max-iter ({}) without converging: some abundances were still \
            changing by more than --convergence-thresh ({}) between M-steps. The estimates \
            are reported as they stand; consider raising --max-iter.",
            opts.max_iter, opts.rel_diff_tol
        );
    }
    if res.dropped_mass > 0.0 {
        info!(
            "{:.3} fragments belonged to equivalence classes whose every target was truncated below the presence threshold",
            res.dropped_mass
        );
    }
    res
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::eq_maps::{BasicEqLabel, EqMap, OrientationProperty, RangeFactorizedEqLabel};

    /// Fragments whose every mapping failed the library-type filter must not
    /// become target-less classes: the EM kernels index a class's first target.
    #[test]
    fn empty_labels_are_not_packed() {
        let mut eqm = EqMap::<BasicEqLabel>::new(OrientationProperty::OrientationAware);
        eqm.add(BasicEqLabel::new(&[], Some(&[])));
        eqm.add(BasicEqLabel::new(&[0, 1, 1, 1], Some(&[0.5, 0.5])));
        let packed_map = PackedEqMap::from_eq_map(&eqm);
        let p = build_packed_eq_classes(&packed_map, &[100.0, 200.0], false);
        assert_eq!(p.counts, vec![1]);
        assert_eq!(p.starts, vec![0, 2]);
        assert_eq!(p.total_count, 1);
        let opts = em_options(100, 1e-2, 1e-2, 1e-8, EmAccel::None);
        let alphas = run_em(&p, &opts).alphas;
        assert!((alphas.iter().sum::<f64>() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn packed_from_basic_eq_map_matches_expected_layout() {
        let mut eqm = EqMap::<BasicEqLabel>::new(OrientationProperty::OrientationAware);
        // labels followed by orientations (contains_ori = true)
        eqm.add(BasicEqLabel::new(&[0, 1, 1, 1], Some(&[0.5, 0.5])));
        eqm.add(BasicEqLabel::new(&[0, 1, 1, 1], Some(&[0.5, 0.5])));
        eqm.add(BasicEqLabel::new(&[2, 1], Some(&[1.0])));
        let packed_map = PackedEqMap::from_eq_map(&eqm);
        let eff_lens = [100.0, 200.0, 50.0];
        let p = build_packed_eq_classes(&packed_map, &eff_lens, true);

        assert_eq!(p.num_txps, 3);
        assert_eq!(p.total_count, 3);
        assert_eq!(p.counts.len(), 2);
        assert_eq!(p.starts.len(), 3);
        assert_eq!(*p.starts.last().unwrap() as usize, p.labels.len());
        assert_eq!(p.labels.len(), 3);
        assert_eq!(p.weights.len(), 3);
        assert_eq!(p.combined.len(), 3);
        for (i, (&tid, &w)) in p.labels.iter().zip(&p.weights).enumerate() {
            assert_eq!(w, 1.0);
            assert_eq!(p.combined[i], 1.0 / eff_lens[tid as usize]);
        }
    }

    #[test]
    fn packed_from_range_factorized_eq_map_uses_bin_probs() {
        let _ = crate::utils::eq_maps::NUM_BINS.set(4.0);
        let mut eqm = EqMap::<RangeFactorizedEqLabel>::new(OrientationProperty::OrientationAware);
        eqm.add(RangeFactorizedEqLabel::new(
            &[0, 1, 1, 1],
            Some(&[0.75, 0.25]),
        ));
        let packed_map = PackedEqMap::from_eq_map(&eqm);
        let eff_lens = [10.0, 20.0];
        let p = build_packed_eq_classes(&packed_map, &eff_lens, true);
        assert_eq!(p.labels, vec![0, 1]);
        // bins: 0.75 -> bin 3 (center 0.875), 0.25 -> bin 1 (center 0.375)
        assert!((p.weights[0] - 0.875).abs() < 1e-12);
        assert!((p.weights[1] - 0.375).abs() < 1e-12);
        assert!((p.combined[0] - 0.875 / 10.0).abs() < 1e-12);
        assert!((p.combined[1] - 0.375 / 20.0).abs() < 1e-12);
        assert_eq!(p.counts, vec![1]);
    }

    /// Entries that differ only in orientation coalesce into one packed class,
    /// and the (`usize`) offsets delimit the full range-factorized keys.
    #[test]
    fn packed_range_factorized_coalesces_orientations() {
        let _ = crate::utils::eq_maps::NUM_BINS.set(4.0);
        let mut eqm = EqMap::<RangeFactorizedEqLabel>::new(OrientationProperty::OrientationAware);
        eqm.add(RangeFactorizedEqLabel::new(
            &[0, 1, 1, 1],
            Some(&[0.75, 0.25]),
        ));
        eqm.add(RangeFactorizedEqLabel::new(
            &[0, 1, 2, 2],
            Some(&[0.75, 0.25]),
        ));
        eqm.add(RangeFactorizedEqLabel::new(&[2, 1], Some(&[1.0])));
        let packed_map = PackedEqMap::from_eq_map(&eqm);
        // keys (targets, then bin ids): [0, 1, 3, 1] and [2, 3]
        assert_eq!(packed_map.eq_label_starts, vec![0, 4, 6]);
        assert_eq!(packed_map.counts, vec![2, 1]);
        let p = build_packed_eq_classes(&packed_map, &[10.0, 20.0, 30.0], true);
        assert_eq!(p.labels, vec![0, 1, 2]);
        assert_eq!(p.starts, vec![0, 2, 3]);
        assert_eq!(p.counts, vec![2, 1]);
        assert_eq!(p.weights.len(), 3);
        assert_eq!(p.combined.len(), 3);
    }

    #[test]
    fn em_conserves_total_count() {
        let mut eqm = EqMap::<BasicEqLabel>::new(OrientationProperty::OrientationAware);
        for _ in 0..30 {
            eqm.add(BasicEqLabel::new(&[0, 1], Some(&[1.0])));
        }
        for _ in 0..70 {
            eqm.add(BasicEqLabel::new(&[1, 1], Some(&[1.0])));
        }
        for _ in 0..100 {
            eqm.add(BasicEqLabel::new(&[0, 1, 1, 1], Some(&[1.0, 1.0])));
        }
        let packed_map = PackedEqMap::from_eq_map(&eqm);
        let eff_lens = [1.0, 1.0];
        let p = build_packed_eq_classes(&packed_map, &eff_lens, false);
        let opts = em_options(1500, 1e-3, 1e-8, 1e-8, EmAccel::Daarem);
        let res = run_em(&p, &opts);
        assert!(res.converged);
        let total: f64 = res.alphas.iter().sum();
        assert!((total - 200.0).abs() < 1e-6, "total = {total}");
        assert!(res.alphas[1] > res.alphas[0]);
    }

    /// Hitting `--max-iter` is reported, not hidden: `converged` is false and
    /// `iters` is the cap (recorded in meta_info.json as em_converged /
    /// em_iters).
    #[test]
    fn em_reports_non_convergence() {
        let mut eqm = EqMap::<BasicEqLabel>::new(OrientationProperty::OrientationAware);
        for _ in 0..30 {
            eqm.add(BasicEqLabel::new(&[0, 1], Some(&[1.0])));
        }
        for _ in 0..100 {
            eqm.add(BasicEqLabel::new(&[0, 1, 1, 1], Some(&[1.0, 1.0])));
        }
        let packed_map = PackedEqMap::from_eq_map(&eqm);
        let p = build_packed_eq_classes(&packed_map, &[1.0, 1.0], false);
        // fewer M-steps than salmon-infer's `min_iter`, so the convergence
        // check never runs
        let opts = em_options(5, 1e-2, 1e-2, 1e-8, EmAccel::None);
        let res = run_em(&p, &opts);
        assert!(!res.converged);
        assert_eq!(res.iters, 5);
        assert!((res.alphas.iter().sum::<f64>() - 130.0).abs() < 1e-6);
    }
}

#[cfg(test)]
mod reproducibility {
    use super::*;
    use crate::utils::eq_maps::{BasicEqLabel, EqMap, OrientationProperty};

    /// A large, heavily ambiguous problem: many multi-target classes with
    /// non-trivial (non-integer) combined weights, so the EM's floating-point
    /// summation order actually matters.
    fn ambiguous_problem() -> (PackedEqClasses, Vec<f64>) {
        const NUM_TXPS: u32 = 5_000;
        let mut eqm = EqMap::<BasicEqLabel>::new(OrientationProperty::OrientationAware);
        // simple LCG so the construction is fixed but irregular
        let mut state = 0x9E37_79B9u64;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as u32
        };
        let mut buf = Vec::new();
        for _ in 0..60_000 {
            let card = 1 + (next() % 6) as usize;
            let anchor = next() % NUM_TXPS;
            buf.clear();
            for j in 0..card {
                // cluster the members so classes overlap heavily
                buf.push((anchor + (j as u32) * (1 + next() % 3)) % NUM_TXPS);
            }
            buf.sort_unstable();
            buf.dedup();
            let n = buf.len();
            buf.extend(std::iter::repeat_n(1u32, n)); // orientations
            let probs = vec![1.0; n];
            let reps = 1 + next() % 20;
            for _ in 0..reps {
                eqm.add(BasicEqLabel::new(&buf, Some(&probs)));
            }
        }
        let packed_map = PackedEqMap::from_eq_map(&eqm);
        let eff_lens: Vec<f64> = (0..NUM_TXPS)
            .map(|i| 200.0 + ((i * 7919) % 2_000) as f64 + 0.37)
            .collect();
        let p = build_packed_eq_classes(&packed_map, &eff_lens, false);
        (p, eff_lens)
    }

    fn in_pool<T: Send>(threads: usize, f: impl FnOnce() -> T + Send) -> T {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(f)
    }

    fn run_in_pool(threads: usize, p: &PackedEqClasses, opts: &EmOptions) -> Vec<f64> {
        in_pool(threads, || run_em(p, opts).alphas)
    }

    /// The headline property of the port: the point estimate is bit-identical
    /// across thread counts on data where summation order matters.
    #[test]
    fn em_point_estimate_is_bit_identical_across_thread_counts() {
        let (p, _) = ambiguous_problem();
        let ambiguous_classes = (0..p.num_classes())
            .filter(|&i| p.starts[i + 1] - p.starts[i] > 1)
            .count();
        assert!(
            ambiguous_classes > 10_000,
            "problem is not ambiguous enough"
        );
        for accel in [EmAccel::None, EmAccel::Squarem, EmAccel::Daarem] {
            let opts = em_options(1500, 1e-3, 1e-8, 1e-8, accel);
            let reference = run_in_pool(1, &p, &opts);
            let total: f64 = reference.iter().sum();
            assert!((total - p.total_count as f64).abs() < 1e-6 * total);
            for threads in [2, 4, 16, 64] {
                let other = run_in_pool(threads, &p, &opts);
                let mismatches = reference
                    .iter()
                    .zip(&other)
                    .filter(|(a, b)| a.to_bits() != b.to_bits())
                    .count();
                assert_eq!(
                    mismatches, 0,
                    "{accel:?}: {mismatches} targets differ between 1 and {threads} threads"
                );
            }
        }
    }

    /// A small ambiguous problem for the inferential-replicate tests.
    fn replicate_problem(keep_weights: bool) -> (PackedEqClasses, Vec<f64>) {
        let mut eqm = EqMap::<BasicEqLabel>::new(OrientationProperty::OrientationAware);
        for i in 0..2000u32 {
            let a = i % 37;
            let b = (i * 7) % 41;
            if i % 3 == 0 {
                eqm.add(BasicEqLabel::new(&[a, 1], Some(&[1.0])));
            } else {
                eqm.add(BasicEqLabel::new(&[a, b + 40, 1, 1], Some(&[1.0, 1.0])));
            }
        }
        let packed_map = PackedEqMap::from_eq_map(&eqm);
        let eff_lens: Vec<f64> = (0..81).map(|i| 100.0 + i as f64).collect();
        let p = build_packed_eq_classes(&packed_map, &eff_lens, keep_weights);
        (p, eff_lens)
    }

    /// Pins the property piscem-infer relies on for `--seed`: the same packed
    /// classes and seed give bit-identical bootstrap replicates whatever the
    /// thread count, and a different seed gives different replicates.
    #[test]
    fn bootstrap_is_reproducible_for_fixed_seed() {
        let (p, eff_lens) = replicate_problem(false);
        let opts = em_options(1500, 1e-3, 1e-8, 1e-8, EmAccel::None);
        let run = |threads: usize, seed: u64| {
            in_pool(threads, || {
                salmon_infer::bootstrap(&p, &opts, salmon_infer::EffLens::new(&eff_lens), 8, seed)
            })
        };
        let a = run(1, 42);
        assert_eq!(a.len(), 8);
        let b = run(8, 42);
        assert_eq!(a, b, "bootstrap replicates differ between 1 and 8 threads");
        let c = run(8, 43);
        assert_ne!(a, c, "bootstrap replicates do not depend on the seed");
    }

    /// The same for Gibbs sampling: the chain layout follows the sample count
    /// (not the thread count) and each chain is seeded from (seed, chain).
    #[test]
    fn gibbs_is_reproducible_for_fixed_seed() {
        let (p, eff_lens) = replicate_problem(true);
        let opts = em_options(1500, 1e-3, 1e-8, 1e-8, EmAccel::None);
        let alphas = run_em(&p, &opts).alphas;
        // 100 samples run as 4 chains, so 8 threads do run them concurrently
        let gibbs_opts = salmon_infer::GibbsOptions {
            num_samples: 100,
            thinning: 2,
            ..salmon_infer::GibbsOptions::default()
        };
        let run = |threads: usize, seed: u64| {
            in_pool(threads, || {
                salmon_infer::gibbs_sample(&p, &eff_lens, &alphas, &gibbs_opts, seed)
            })
        };
        let a = run(1, 42);
        assert_eq!(a.len(), 100);
        let b = run(8, 42);
        assert_eq!(a, b, "Gibbs samples differ between 1 and 8 threads");
        let c = run(8, 43);
        assert_ne!(a, c, "Gibbs samples do not depend on the seed");
    }
}
