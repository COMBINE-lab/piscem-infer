use ahash::AHashMap;

pub type BasicEqMap = EqMap<BasicEqLabel>;
pub type RangeFactorizedEqMap = EqMap<RangeFactorizedEqLabel>;

const NUM_BINS: f64 = 32_f64;

pub enum OrientationProperty {
    OrientationAware,
    #[allow(dead_code)]
    OrientationAgnostic,
}

pub enum EqMapType {
    BasicEqMap,
    RangeFactorizedEqMap,
}

pub trait EqLabel: TargetLabels + std::hash::Hash + PartialEq + Eq {
    type LabelRefT<'a>;
    fn new(labels: &[u32], probs: Option<&[f64]>) -> Self;
    fn new_ref(labels: &[u32], has_ori: bool) -> Self::LabelRefT<'_>;
}

#[derive(Hash, PartialEq, Eq)]
pub struct RangeFactorizedEqLabelRef<'a> {
    pub targets_and_bins: &'a [u32],
    pub contains_ori: bool,
}

impl<'a> TargetLabels for RangeFactorizedEqLabelRef<'a> {
    #[inline]
    fn target_labels(&self, with_ori: bool) -> &[u32] {
        let l = self.targets_and_bins.len() / if with_ori { 3 } else { 2 };
        &self.targets_and_bins[..l]
    }
    #[inline]
    fn target_probs(&self, with_ori: bool) -> impl Iterator<Item = f64> {
        let l = self.targets_and_bins.len() / if with_ori { 3 } else { 2 };
        let num_bins = NUM_BINS;
        let half_bin_width = 0.5 / num_bins;
        RangeFactorizedBinIterator {
            bin_iterator: self.targets_and_bins[l..2 * l].iter(),
            num_bins,
            half_bin_width,
        }
    }
}

#[derive(Hash, PartialEq, Eq)]
pub struct RangeFactorizedEqLabel {
    pub targets_and_bins: Vec<u32>,
}

impl EqLabel for RangeFactorizedEqLabel {
    type LabelRefT<'a> = RangeFactorizedEqLabelRef<'a>;

    fn new_ref(labels: &[u32], has_ori: bool) -> RangeFactorizedEqLabelRef {
        RangeFactorizedEqLabelRef {
            targets_and_bins: labels,
            contains_ori: has_ori
        }
    }

    fn new(labels: &[u32], probs: Option<&[f64]>) -> Self {
        let probs = probs.expect("probs *must* be present for range factorized equivalence class");
        // first, ensure probs are normalized
        let tot_prob: f64 = probs.iter().sum();
        let num_labels = probs.len();
        let num_bins = NUM_BINS as usize;//4_usize + ((num_labels as f64).sqrt()).ceil() as usize;

        let (just_labels, oris) = labels.split_at(num_labels);
        let mut targets_and_bins: Vec<u32> = just_labels.into();
        targets_and_bins.extend(probs.iter().map(|&prob| {
            let p: f64 = prob / tot_prob;
            // Handle edge case where p = 1.0 (should go to last bin)
            if p >= 1.0 {
                (num_bins - 1) as u32
            } else {
                // Calculate bin index: floor(p * num_bins)
                (p * num_bins as f64) as u32
            }
        }));
        targets_and_bins.extend_from_slice(oris);
        Self { targets_and_bins }
    }
}

struct RangeFactorizedBinIterator<'a> {
    bin_iterator: std::slice::Iter<'a, u32>,
    num_bins: f64,
    half_bin_width: f64,
}

impl<'a> Iterator for RangeFactorizedBinIterator<'a> {
    type Item = f64;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(x) = self.bin_iterator.next() {
            Some((*x as f64) / self.num_bins + self.half_bin_width)
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.bin_iterator.size_hint()
    }
}

impl<'a> ExactSizeIterator for RangeFactorizedBinIterator<'a> {}

impl TargetLabels for RangeFactorizedEqLabel {
    /// NOTE: The 2 * l below here is a huge hack. We want to include the bins in the
    /// contents of the packed map, and this is how we do it, but we certainly should
    /// find a more elegant and less hacky way.
    #[inline]
    fn target_labels(&self, with_ori: bool) -> &[u32] {
        // if this is an orientation-aware equivalence class factorization, then
        // the targets_and_bins vector contains (labels, oris, probs), otherwise
        // it contains just (labels, probs)
        let l = self.targets_and_bins.len() / if with_ori { 3 } else { 2 };
        &self.targets_and_bins[..2 * l]
    }

    #[inline]
    fn target_probs(&self, with_ori: bool) -> impl Iterator<Item = f64> {
        // if this is an orientation-aware equivalence class factorization, then
        // the targets_and_bins vector contains (labels, oris, probs), otherwise
        // it contains just (labels, probs)
        let l = self.targets_and_bins.len() / if with_ori { 3 } else { 2 };
        let num_bins = NUM_BINS;//4_f64 + ((l as f64).sqrt()).ceil();
        let half_bin_width = 0.5 / num_bins;
        RangeFactorizedBinIterator {
            bin_iterator: self.targets_and_bins[l..2 * l].iter(),
            num_bins,
            half_bin_width,
        }
    }
}

#[derive(Hash, PartialEq, Eq)]
pub struct BasicEqLabelRef<'a> {
    pub targets: &'a [u32],
    pub contains_ori: bool,
}

#[derive(Hash, PartialEq, Eq)]
pub struct BasicEqLabel {
    pub targets: Vec<u32>,
}

impl EqLabel for BasicEqLabel {
    type LabelRefT<'a> = BasicEqLabelRef<'a>;

    fn new_ref(labels: &[u32], has_ori: bool) -> BasicEqLabelRef {
        BasicEqLabelRef {
            targets: labels,
            contains_ori: has_ori
        }
    }

    fn new(targets: &[u32], probs: Option<&[f64]>) -> Self {
        Self {
            targets: targets.into(),
        }
    }
}

pub trait TargetLabels {
    fn target_labels(&self, with_ori: bool) -> &[u32];
    fn target_probs(&self, with_ori: bool) -> impl Iterator<Item = f64>;
}

impl TargetLabels for BasicEqLabel {
    // return the slice of identifiers (u32s) that correspond
    // to the target ids. If the EqLabel was built without orientations
    // this is the whole vector, otherwise it's the first half.
    #[inline]
    fn target_labels(&self, with_ori: bool) -> &[u32] {
        // number of targets is total length / 2
        let nt = if with_ori {
            self.targets.len() >> 1
        } else {
            self.targets.len()
        };
        &self.targets[0..nt]
    }

    #[inline]
    fn target_probs(&self, with_ori: bool) -> impl Iterator<Item = f64> {
        // number of targets is total length / 2
        let nt = if with_ori {
            self.targets.len() >> 1
        } else {
            self.targets.len()
        };
        std::iter::repeat_n(1.0_f64, nt)
    }
}

impl<'a> TargetLabels for BasicEqLabelRef<'a> {
    // return the slice of identifiers (u32s) that correspond
    // to the target ids. If the EqLabel was built without orientations
    // this is the whole vector, otherwise it's the first half.
    #[inline]
    fn target_labels(&self, with_ori: bool) -> &[u32] {
        // number of targets is total length / 2
        let nt = if with_ori {
            self.targets.len() >> 1
        } else {
            self.targets.len()
        };
        &self.targets[0..nt]
    }

    #[inline]
    fn target_probs(&self, with_ori: bool) -> impl Iterator<Item = f64> {
        // number of targets is total length / 2
        let nt = if with_ori {
            self.targets.len() >> 1
        } else {
            self.targets.len()
        };
        std::iter::repeat_n(1.0_f64, nt)
    }
}

/// An equivalence class map that maps target equivalence
/// classes to their counts.
pub struct EqMap<EqLabelT: TargetLabels> {
    pub count_map: AHashMap<EqLabelT, usize>,
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
    // phantom: std::marker::PhantomData<EqLabelT>
}

impl PackedEqMap {

    pub fn from_eq_map<EqLabelT: EqLabel>(eqm: &EqMap<EqLabelT>) -> Self {
        let mut eq_labels = Vec::<u32>::with_capacity(eqm.len() * 5);
        let mut counts = Vec::<usize>::with_capacity(eqm.len());
        let mut eq_label_starts = Vec::<u32>::with_capacity(eqm.count_map.len() + 1);

        eq_label_starts.push(0);
        for (eq_lab, count) in eqm.iter() {
            //eprintln!("eq_lab = {eq_lab:#?}");
            eq_labels.extend_from_slice(eq_lab);
            eq_label_starts.push(eq_labels.len() as u32);
            counts.push(*count);
        }

        Self {
            eq_labels,
            eq_label_starts,
            counts,
            contains_ori: eqm.contains_ori,
            //phantom: std::marker::PhantomData
        }
    }

    pub fn refs_for_eqc(&self, idx: usize) -> RangeFactorizedEqLabelRef {
        //<EqLabelT as EqLabel>::LabelRefT<'_> {
        let s: usize = self.eq_label_starts[idx] as usize;
        let e: usize = self.eq_label_starts[idx + 1] as usize;
        // if we encode orientation, then it's the first half of the
        // label, otherwise it's the whole label.
        let l = e - s;
        // right now contains_ori is alwayws false becuase
        // we have stripped the orientations from the
        // label vector when building the
        // PackedEqMap.
        RangeFactorizedEqLabel::new_ref(&self.eq_labels[s..(s+l)], false)
    }

    pub fn len(&self) -> usize {
        self.counts.len()
    }

    /*
    #[allow(dead_code)]
    pub fn iter(&self) -> PackedEqEntryIter {
        PackedEqEntryIter {
            counter: 0,
            underlying_packed_map: self,
        }
    }
    */

    pub fn iter_labels(&self) -> PackedEqLabelIter {
        PackedEqLabelIter {
            counter: 0,
            underlying_packed_map: self,
        }
    }

    pub fn total_weight(&self) -> usize {
        self.counts.iter().sum()
    }
}

/// An iterator over the labels of the
/// `PackedEqMap`.
pub struct PackedEqLabelIter<'a> {
    counter: u32,
    underlying_packed_map: &'a PackedEqMap,
}

impl<'a> Iterator for PackedEqLabelIter<'a> {
    type Item = RangeFactorizedEqLabelRef<'a>;

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

impl<'a> ExactSizeIterator for PackedEqLabelIter<'a> {}

/*
/// An iterator over the equivalence classes of the
/// `PackedEqMap`.
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

    fn size_hint(&self) -> (usize, Option<usize>) {
        let rem = self.underlying_packed_map.len() - self.counter as usize;
        (rem, Some(rem))
    }
}

impl<'a> ExactSizeIterator for PackedEqEntryIter<'a> {}
*/

impl<EqLabelT> EqMap<EqLabelT>
where
    EqLabelT: EqLabel,
{
    /// Create a new equivalence class map, if
    /// `cotntains_ori` is true, the equivalence class
    /// definitions will include the orientation flags,
    /// if false, they will not.
    pub fn new(ori_prop: OrientationProperty) -> Self {
        Self {
            count_map: AHashMap::<EqLabelT, usize>::new(),
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
    pub fn iter(&self) -> EqEntryIter<EqLabelT> {
        EqEntryIter {
            underlying_iter: self.count_map.iter(),
            contains_ori: self.contains_ori,
        }
    }
}

pub struct EqEntryIter<'a, EqLabelT: EqLabel> {
    underlying_iter: std::collections::hash_map::Iter<'a, EqLabelT, usize>,
    contains_ori: bool,
}

impl<'a, EqLabelT: EqLabel> Iterator for EqEntryIter<'a, EqLabelT> {
    type Item = (&'a [u32], &'a usize);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        match self.underlying_iter.next() {
            Some((k, v)) => Some((k.target_labels(self.contains_ori), v)),
            None => None,
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.underlying_iter.size_hint()
    }
}

impl<'a, EqLabelT: EqLabel> ExactSizeIterator for EqEntryIter<'a, EqLabelT> {}
