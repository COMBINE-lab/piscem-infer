use ahash::AHashMap;
use atomic_float::AtomicF64;
use rand::prelude::*;
use rand::rng;
use rand_distr::weighted::WeightedAliasIndex;
use rayon::prelude::*;
use std::sync::atomic::Ordering;
use tracing::info;

pub enum OrientationProperty {
    OrientationAware,
    #[allow(dead_code)]
    OrientationAgnostic,
}

#[derive(Hash, PartialEq, Eq)]
pub struct EqLabel {
    pub targets: Vec<u32>,
}

impl EqLabel {
    // return the slice of identifiers (u32s) that correspond
    // to the target ids. If the EqLabel was built without orientations
    // this is the whole vector, otherwise it's the first half.
    #[inline]
    pub fn target_labels(&self, with_ori: bool) -> &[u32] {
        // number of targets is total length / 2
        let nt = if with_ori {
            self.targets.len() >> 1
        } else {
            self.targets.len()
        };
        &self.targets[0..nt]
    }
}

/// An equivalence class map that maps target equivalence
/// classes to their counts.
pub struct EqMap {
    pub count_map: AHashMap<EqLabel, usize>,
    pub contains_ori: bool,
}

pub struct PackedEqMap {
    /// the packed list of all equivalence class labels
    pub eq_labels: Vec<u32>,
    /// vector that deliniates where each equivalence class label
    /// begins and ends.  The label for equivalence class i begins
    /// at offset eq_label_starts[i], and it ends at
    /// eq_label_starts[i+1].  The length of this vector is 1 greater
    /// than the number of equivalence classes.
    pub eq_label_starts: Vec<u32>,
    /// the vector of counts for each equivalence class
    pub counts: Vec<usize>,
    /// whether or not the underlying equivalence map was built
    /// with orientation information encoded or not.
    pub contains_ori: bool,
}

impl PackedEqMap {
    pub fn from_eq_map(eqm: &EqMap) -> Self {
        let mut eq_labels = Vec::<u32>::with_capacity(eqm.len() * 5);
        let mut counts = Vec::<usize>::with_capacity(eqm.len());
        let mut eq_label_starts = Vec::<u32>::with_capacity(eqm.len() + 1);

        eq_label_starts.push(0);
        for (eq_lab, count) in eqm.iter() {
            eq_labels.extend_from_slice(eq_lab);
            eq_label_starts.push(eq_labels.len() as u32);
            counts.push(*count);
        }

        Self {
            eq_labels,
            eq_label_starts,
            counts,
            contains_ori: eqm.contains_ori,
        }
    }

    pub fn refs_for_eqc(&self, idx: usize) -> &[u32] {
        let s: usize = self.eq_label_starts[idx] as usize;
        let e: usize = self.eq_label_starts[idx + 1] as usize;
        // if we encode orientation, then it's the first half of the
        // label, otherwise it's the whole label.
        let l = if self.contains_ori {
            e - s
        } else {
            (e - s) >> 1
        };
        &self.eq_labels[s..(s + l)]
    }

    pub fn len(&self) -> usize {
        self.counts.len()
    }

    #[allow(dead_code)]
    pub fn iter(&self) -> PackedEqEntryIter {
        PackedEqEntryIter {
            counter: 0,
            underlying_packed_map: self,
        }
    }

    pub fn iter_labels(&self) -> PackedEqLabelIter {
        PackedEqLabelIter {
            counter: 0,
            underlying_packed_map: self,
        }
    }
}

pub struct PackedEqLabelIter<'a> {
    counter: u32,
    underlying_packed_map: &'a PackedEqMap,
}

impl<'a> Iterator for PackedEqLabelIter<'a> {
    type Item = &'a [u32];

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let c = self.counter as usize;
        if c < self.underlying_packed_map.len() {
            self.counter += 1;
            Some(self.underlying_packed_map.refs_for_eqc(c))
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let rem = self.underlying_packed_map.len() - self.counter as usize;
        (rem, Some(rem))
    }
}

pub struct PackedEqEntryIter<'a> {
    counter: u32,
    underlying_packed_map: &'a PackedEqMap,
}

impl<'a> Iterator for PackedEqEntryIter<'a> {
    type Item = (&'a [u32], &'a usize);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        let c = self.counter as usize;
        if c < self.underlying_packed_map.len() {
            self.counter += 1;
            Some((
                self.underlying_packed_map.refs_for_eqc(c),
                &self.underlying_packed_map.counts[c],
            ))
        } else {
            None
        }
    }
}

impl<'a> ExactSizeIterator for PackedEqEntryIter<'a> {
    fn len(&self) -> usize {
        self.underlying_packed_map.len() - self.counter as usize
    }
}

impl EqMap {
    /// Create a new equivalence class map, if
    /// `cotntains_ori` is true, the equivalence class
    /// definitions will include the orientation flags,
    /// if false, they will not.
    pub fn new(ori_prop: OrientationProperty) -> Self {
        Self {
            count_map: AHashMap::<EqLabel, usize>::new(),
            contains_ori: match ori_prop {
                OrientationProperty::OrientationAware => true,
                OrientationProperty::OrientationAgnostic => false,
            },
        }
    }

    /// The number of equivalence classes
    pub fn len(&self) -> usize {
        self.count_map.len()
    }

    /// Return an iterator over the equivalence class
    /// map iterator.
    pub fn iter(&self) -> EqEntryIter {
        EqEntryIter {
            underlying_iter: self.count_map.iter(),
            contains_ori: self.contains_ori,
        }
    }
}

pub struct EqEntryIter<'a> {
    underlying_iter: std::collections::hash_map::Iter<'a, EqLabel, usize>,
    contains_ori: bool,
}

impl<'a> Iterator for EqEntryIter<'a> {
    type Item = (&'a [u32], &'a usize);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        match self.underlying_iter.next() {
            Some((k, v)) => Some((k.target_labels(self.contains_ori), v)),
            None => None,
        }
    }
}

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

const ABSENCE_THRESH: f64 = 1e-8;
const RELDIFF_THRESH: f64 = 1e-3;

#[inline]
fn m_step_par(
    eq_iterates: &[(&[u32], &usize)],
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
            for e in k.iter() {
                let w = prev_count[*e as usize].load(Ordering::Relaxed) * inv_eff_lens[*e as usize];
                weights.push(w);
                denom += w;
            }
            if denom > 1e-8 {
                let count_over_denom = count / denom;
                for (target_id, w) in k.iter().zip(weights.iter()) {
                    let inc = count_over_denom * w;
                    curr_counts[*target_id as usize].fetch_add(inc, Ordering::AcqRel);
                }
            }
            weights.clear();
        },
    );
}

#[inline]
fn m_step(
    eq_map: &PackedEqMap,
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
        for e in k.iter() {
            let w = prev_count[*e as usize] * inv_eff_lens[*e as usize];
            weights.push(w);
            denom += w;
        }
        if denom > 1e-8 {
            let count_over_denom = count / denom;
            for (target_id, w) in k.iter().zip(weights.iter()) {
                curr_counts[*target_id as usize] += count_over_denom * w;
            }
        }
        weights.clear();
    }
}

/// Holds the info relevant for running the EM algorithm
pub struct EMInfo<'eqm, 'el> {
    pub eq_map: &'eqm PackedEqMap,
    pub eff_lens: &'el [f64],
    pub max_iter: u32,
    pub convergence_thresh: f64,
}

pub fn do_bootstrap(em_info: &EMInfo, num_boot: usize) -> Vec<Vec<f64>> {
    let converge_thresh = em_info.convergence_thresh;
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
                    if prev_counts[i] > ABSENCE_THRESH {
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
                if *x < ABSENCE_THRESH {
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

pub fn em_par(em_info: &EMInfo, nthreads: usize) -> Vec<f64> {
    let converge_thresh = em_info.convergence_thresh;
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

    let eq_iterates: Vec<(&[u32], &usize)> = eq_map.iter_labels().zip(&eq_map.counts).collect();

    let mut rel_diff = 0.0_f64;
    let mut niter = 0_u32;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthreads)
        .build()
        .unwrap();

    pool.install(|| {
        while niter < max_iter {
            m_step_par(
                &eq_iterates,
                &mut prev_counts,
                &inv_eff_lens,
                &mut curr_counts,
            );

            //std::mem::swap(&)
            for i in 0..curr_counts.len() {
                let pci = prev_counts[i].load(Ordering::Relaxed);
                if pci > ABSENCE_THRESH {
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
            if x.load(Ordering::Relaxed) < ABSENCE_THRESH {
                x.store(0.0, Ordering::Relaxed);
            }
        });
        m_step_par(
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

pub fn em(em_info: &EMInfo) -> Vec<f64> {
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
            if prev_counts[i] > ABSENCE_THRESH {
                let rd = (curr_counts[i] - prev_counts[i]) / prev_counts[i];
                rel_diff = if rel_diff > rd { rel_diff } else { rd };
            }
        }

        std::mem::swap(&mut prev_counts, &mut curr_counts);
        curr_counts.fill(0.0_f64);

        if rel_diff < RELDIFF_THRESH {
            break;
        }
        niter += 1;
        if niter % 100 == 0 {
            info!("iteration {}; rel diff {:.3}", niter, rel_diff);
        }
        rel_diff = 0.0_f64;
    }

    prev_counts.iter_mut().for_each(|x| {
        if *x < ABSENCE_THRESH {
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
