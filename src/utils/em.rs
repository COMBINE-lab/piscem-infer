use ahash::AHashMap;

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

/// Go through the set of references (`ref_lens`), and adjust their lengths according to
/// the computed conditional means `cond_means` of the fragment length distribution.
pub fn adjust_ref_lengths(ref_lens: &[u32], cond_means: &[f64]) -> Vec<f64> {
    let tmean = cond_means.last().unwrap();
    let el = ref_lens
        .iter()
        .map(|rli| {
            let rl = *rli as usize;
            let adj_len = if rl > cond_means.len() {
                (rl as f64) - tmean
            } else {
                (rl as f64) - cond_means[rl]
            };
            if adj_len >= 1.0 {
                adj_len
            } else {
                rl as f64
            }
        })
        .collect::<Vec<f64>>();

    el
}

#[inline]
fn m_step(eq_map: &EqMap, prev_count: &[f64], eff_lens: &[f64], curr_counts: &mut [f64]) {
    for (k, v) in eq_map.iter() {
        let mut denom = 0.0_f64;
        let count = *v as f64;
        for target_id in k {
            denom += prev_count[*target_id as usize] / eff_lens[*target_id as usize];
        }

        if denom > 1e-8 {
            for target_id in k {
                curr_counts[*target_id as usize] += count
                    * ((prev_count[*target_id as usize] / eff_lens[*target_id as usize]) / denom);
            }
        }
    }
}

/// Holds the info relevant for running the EM algorithm
pub struct EMInfo<'eqm, 'el> {
    pub eq_map: &'eqm EqMap,
    pub eff_lens: &'el [f64],
    pub max_iter: u32,
    pub convergence_thresh: f64,
}

pub fn em(em_info: &EMInfo) -> Vec<f64> {
    let eq_map = em_info.eq_map;
    let eff_lens = em_info.eff_lens;
    let max_iter = em_info.max_iter;
    let total_weight: f64 = eq_map.count_map.values().sum::<usize>() as f64;

    // init
    let avg = total_weight / (eff_lens.len() as f64);
    let mut prev_counts = vec![avg; eff_lens.len()];
    let mut curr_counts = vec![0.0f64; eff_lens.len()];

    let mut rel_diff = 0.0_f64;
    let mut niter = 0_u32;

    while niter < max_iter {
        m_step(eq_map, &prev_counts, eff_lens, &mut curr_counts);

        //std::mem::swap(&)
        for i in 0..curr_counts.len() {
            if prev_counts[i] > 1e-8 {
                let rd = (curr_counts[i] - prev_counts[i]) / prev_counts[i];
                rel_diff = if rel_diff > rd { rel_diff } else { rd };
            }
        }

        std::mem::swap(&mut prev_counts, &mut curr_counts);
        curr_counts.fill(0.0_f64);

        if rel_diff < 1e-3 {
            break;
        }
        niter += 1;
        if niter % 100 == 0 {
            eprintln!("iteration {}; rel diff {}", niter, rel_diff);
        }
        rel_diff = 0.0_f64;
    }

    prev_counts.iter_mut().for_each(|x| {
        if *x < 1e-8 {
            *x = 0.0
        }
    });
    m_step(eq_map, &prev_counts, eff_lens, &mut curr_counts);

    curr_counts
}
