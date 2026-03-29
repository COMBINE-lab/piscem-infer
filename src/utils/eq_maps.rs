use ahash::AHashMap;

/// aliases for referring to the different basic types of eq class maps
pub type BasicEqMap = EqMap<BasicEqLabel>;
pub type RangeFactorizedEqMap = EqMap<RangeFactorizedEqLabel>;

/// the default number of bins to use for range-factorized equivalence classes
pub static NUM_BINS: std::sync::OnceLock<f64> = std::sync::OnceLock::new();

/// number of positional bins for positional equivalence classes (default 1 = disabled)
pub static NUM_POS_BINS: std::sync::OnceLock<f64> = std::sync::OnceLock::new();

/// Helper: number of segments in a RangeFactorized packed label.
fn rf_num_segments(has_ori: bool) -> usize {
    let has_pos = NUM_POS_BINS.get().is_some_and(|&n| n > 1.0);
    match (has_pos, has_ori) {
        (false, false) => 2, // targets, prob_bins
        (false, true) => 3,  // targets, prob_bins, oris
        (true, false) => 3,  // targets, prob_bins, pos_bins
        (true, true) => 4,   // targets, prob_bins, pos_bins, oris
    }
}

/// whether or not the equivalence classes differentiate between fragments
/// mapping in different orientations
pub enum OrientationProperty {
    OrientationAware,
    #[allow(dead_code)]
    OrientationAgnostic,
}

/// discern betweeen the types of equivalence class maps
pub enum EqMapType {
    BasicEqMap,
    RangeFactorizedEqMap,
}

/// The `EqLabel` trait gives us the ability to create things that are equivalence class labels.
/// These things can be created with or without associated probabilities.
/// Also, any type implementing this trait will have an associated `LabelRefT` type that will allow
/// the creation of borrowed references that can enumerate labels and probabilities.
/// These types must also be hashable and comparable to be used as keys in HashMaps.
pub trait EqLabel: TargetLabels + std::hash::Hash + PartialEq + Eq + Sync {
    type LabelRefT<'a>: TargetLabelsRef;
    fn new(labels: &[u32], probs: Option<&[f64]>, pos_bins: Option<&[u32]>) -> Self;
    fn new_ref(labels: &[u32], has_ori: bool) -> Self::LabelRefT<'_>;
}

/// This trait ensures that we can get a list of the labels of this equivalence class
pub trait TargetLabels {
    fn target_labels(&self, with_ori: bool) -> &[u32];
    fn extract_key_for_packed_map(&self, has_ori: bool) -> &[u32];
}

/// This trait ensures we can get a list of the lables (and the probabilities of this
/// equivalence class label refererence)
pub trait TargetLabelsRef: Sync {
    fn target_labels(&self) -> &[u32];
    fn target_probs(&self) -> impl Iterator<Item = f64>;
}

// === basic equivalence classes

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

    /// we create a new reference by passing along the label refs we
    /// are given
    fn new_ref(labels: &[u32], has_ori: bool) -> BasicEqLabelRef<'_> {
        BasicEqLabelRef {
            targets: labels,
            contains_ori: has_ori,
        }
    }

    /// we ignore probabilities and position bins in the basic equivalence class,
    /// and just pass along the targets we are given as the labels
    fn new(targets: &[u32], _probs: Option<&[f64]>, _pos_bins: Option<&[u32]>) -> Self {
        Self {
            targets: targets.into(),
        }
    }
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

    // return the slice of identifiers (u32s) that corresponds
    // to just the labels
    #[inline]
    fn extract_key_for_packed_map(&self, has_ori: bool) -> &[u32] {
        // number of targets is total length / 2
        let nt = if has_ori {
            self.targets.len() >> 1
        } else {
            self.targets.len()
        };
        &self.targets[0..nt]
    }
}

impl<'a> TargetLabelsRef for BasicEqLabelRef<'a> {
    /// return the slice of identifiers (u32s) that correspond
    /// to the target ids. If the EqLabel was built without orientations
    /// this is the whole vector, otherwise it's the first half.
    #[inline]
    fn target_labels(&self) -> &[u32] {
        let with_ori = self.contains_ori;
        // number of targets is total length / 2
        let nt = if with_ori {
            self.targets.len() >> 1
        } else {
            self.targets.len()
        };
        &self.targets[0..nt]
    }

    /// since there were no target probabilities provided when
    /// creating this equivalence class, we will instead  
    /// simply produce a stream of 1s
    #[inline]
    fn target_probs(&self) -> impl Iterator<Item = f64> {
        let with_ori = self.contains_ori;
        let nt = if with_ori {
            self.targets.len() >> 1
        } else {
            self.targets.len()
        };
        std::iter::repeat_n(1.0_f64, nt)
    }
}

// ===== Range factorized equivalence classes

#[derive(Hash, PartialEq, Eq)]
pub struct RangeFactorizedEqLabel {
    pub targets_and_bins: Vec<u32>,
}

impl EqLabel for RangeFactorizedEqLabel {
    type LabelRefT<'a> = RangeFactorizedEqLabelRef<'a>;

    fn new_ref(labels: &[u32], has_ori: bool) -> RangeFactorizedEqLabelRef<'_> {
        RangeFactorizedEqLabelRef {
            targets_and_bins: labels,
            contains_ori: has_ori,
        }
    }

    /// create the range factorized equivalence class from labels,
    /// associated probabilities, and optional position bins.
    fn new(labels: &[u32], probs: Option<&[f64]>, pos_bins: Option<&[u32]>) -> Self {
        let probs = probs.expect("probs *must* be present for range factorized equivalence class");

        let tot_prob: f64 = probs.iter().sum();
        let num_labels = probs.len();
        let num_bins = *NUM_BINS.get().unwrap() as usize;
        let num_pos_bins = NUM_POS_BINS.get().copied().unwrap_or(1.0) as usize;

        // Layout: [targets..., prob_bins..., pos_bins..., oris...]
        // where pos_bins are only present when NUM_POS_BINS > 1
        let (just_labels, oris) = labels.split_at(num_labels);
        let mut targets_and_bins: Vec<u32> = just_labels.into();

        // Probability bins
        targets_and_bins.extend(probs.iter().map(|&prob| {
            let p: f64 = prob / tot_prob;
            if p >= 1.0 {
                (num_bins - 1) as u32
            } else {
                (p * num_bins as f64) as u32
            }
        }));

        // Position bins (only when enabled)
        if num_pos_bins > 1 {
            if let Some(pb) = pos_bins {
                targets_and_bins.extend_from_slice(pb);
            } else {
                // Default to bin 0 if not provided
                targets_and_bins.extend(std::iter::repeat(0u32).take(num_labels));
            }
        }

        // Orientations
        targets_and_bins.extend_from_slice(oris);
        Self { targets_and_bins }
    }
}

/// allows us to iterate over the bins of the range factorized equivalence
/// classes and, for each bin, return the associated conditional probability
struct RangeFactorizedBinIterator<'a> {
    bin_iterator: std::slice::Iter<'a, u32>,
    num_bins: f64,
    half_bin_width: f64,
}

impl<'a> Iterator for RangeFactorizedBinIterator<'a> {
    type Item = f64;

    /// for the next bin, returns the probability associated with
    /// the center of the bin
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
    #[inline]
    fn target_labels(&self, with_ori: bool) -> &[u32] {
        let nseg = rf_num_segments(with_ori);
        let l = self.targets_and_bins.len() / nseg;
        &self.targets_and_bins[..l]
    }

    /// Extract the key for the packed map: targets + prob_bins + pos_bins
    /// (everything except orientations).
    fn extract_key_for_packed_map(&self, has_ori: bool) -> &[u32] {
        let nseg = rf_num_segments(has_ori);
        let l = self.targets_and_bins.len() / nseg;
        // Everything except the last segment (oris) if has_ori, else everything
        let key_segments = if has_ori { nseg - 1 } else { nseg };
        &self.targets_and_bins[..key_segments * l]
    }
}

#[derive(Hash, PartialEq, Eq)]
pub struct RangeFactorizedEqLabelRef<'a> {
    pub targets_and_bins: &'a [u32],
    pub contains_ori: bool,
}

/// treat a slice of u32s as labels and or probabilities
impl<'a> TargetLabelsRef for RangeFactorizedEqLabelRef<'a> {
    /// gets the labels associated with this reference
    #[inline]
    fn target_labels(&self) -> &[u32] {
        let nseg = rf_num_segments(self.contains_ori);
        let l = self.targets_and_bins.len() / nseg;
        &self.targets_and_bins[..l]
    }

    /// returns an iterator over the conditional probabilities from the
    /// probability bins of this reference (position weight = 1.0 under
    /// uniform coverage, folded in for future bias models)
    #[inline]
    fn target_probs(&self) -> impl Iterator<Item = f64> {
        let nseg = rf_num_segments(self.contains_ori);
        let l = self.targets_and_bins.len() / nseg;
        let num_bins = *NUM_BINS.get().unwrap();
        let half_bin_width = 0.5 / num_bins;
        RangeFactorizedBinIterator {
            bin_iterator: self.targets_and_bins[l..2 * l].iter(),
            num_bins,
            half_bin_width,
        }
    }
}

// equivalence class maps

/// An equivalence class map that maps target equivalence
/// classes to their counts.
pub struct EqMap<EqLabelT> {
    pub count_map: AHashMap<EqLabelT, usize>,
    pub contains_ori: bool,
}

impl<EqLabelT: EqLabel> EqMap<EqLabelT> {
    /// add the equivalence class label to this `EQMap`'s
    /// count map with a count of 1 if it hasen't yet been seen,
    /// or increment the count if it already exists.
    pub fn add(&mut self, lab: EqLabelT) -> &'_ mut usize {
        self.count_map
            .entry(lab)
            .and_modify(|counter| *counter += 1)
            .or_insert(1)
    }

    /// Merge another EqMap into this one by summing counts for shared labels.
    pub fn merge(&mut self, other: EqMap<EqLabelT>) {
        for (lab, count) in other.count_map {
            *self.count_map.entry(lab).or_insert(0) += count;
        }
    }
}

pub struct PackedEqMap<EqLabelT> {
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
    #[allow(dead_code)]
    pub contains_ori: bool,
    /// we need to be able to hold the type of the label of this
    /// equivalence class, so that we can generated the associated
    /// `LabelRefT`
    phantom: std::marker::PhantomData<EqLabelT>,
}

impl<EqLabelT: EqLabel> PackedEqMap<EqLabelT> {
    /// Reconstruct a `PackedEqMap` from raw vectors (e.g. after deserialization).
    pub fn from_raw(
        eq_labels: Vec<u32>,
        eq_label_starts: Vec<u32>,
        counts: Vec<usize>,
        contains_ori: bool,
    ) -> Self {
        Self {
            eq_labels,
            eq_label_starts,
            counts,
            contains_ori,
            phantom: std::marker::PhantomData,
        }
    }

    pub fn from_eq_map(eqm: &EqMap<EqLabelT>) -> Self {
        let mut eq_labels = Vec::<u32>::with_capacity(eqm.len() * 5);
        let mut counts = Vec::<usize>::with_capacity(eqm.len());
        let mut eq_label_starts = Vec::<u32>::with_capacity(eqm.count_map.len() + 1);

        eq_label_starts.push(0);
        for (eq_lab, count) in eqm.full_key_iter() {
            eq_labels.extend_from_slice(eq_lab);
            eq_label_starts.push(eq_labels.len() as u32);
            counts.push(*count);
        }

        Self {
            eq_labels,
            eq_label_starts,
            counts,
            contains_ori: eqm.contains_ori,
            phantom: std::marker::PhantomData,
        }
    }

    pub fn refs_for_eqc(&self, idx: usize) -> <EqLabelT as EqLabel>::LabelRefT<'_> {
        let s: usize = self.eq_label_starts[idx] as usize;
        let e: usize = self.eq_label_starts[idx + 1] as usize;
        // if we encode orientation, then it's the first half of the
        // label, otherwise it's the whole label.
        let l = e - s;
        // right now contains_ori is alwayws false becuase
        // we have stripped the orientations from the
        // label vector when building the
        // PackedEqMap.
        EqLabelT::new_ref(&self.eq_labels[s..(s + l)], false)
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

    pub fn iter_labels(&self) -> PackedEqLabelIter<'_, EqLabelT> {
        PackedEqLabelIter {
            counter: 0,
            underlying_packed_map: self,
        }
    }

    /// Returns the number of distinct targets in equivalence class `idx`.
    /// Works for both `BasicEqLabel` (label = targets) and
    /// `RangeFactorizedEqLabel` (label = targets + bins).
    pub fn num_targets_in_eqc(&self, idx: usize) -> usize {
        self.refs_for_eqc(idx).target_labels().len()
    }

    pub fn total_weight(&self) -> usize {
        self.counts.iter().sum()
    }
}

/// An iterator over the labels of the
/// `PackedEqMap`.
pub struct PackedEqLabelIter<'a, EqLabelT> {
    counter: u32,
    underlying_packed_map: &'a PackedEqMap<EqLabelT>,
}

impl<'a, EqLabelT: EqLabel> Iterator for PackedEqLabelIter<'a, EqLabelT> {
    type Item = <EqLabelT as EqLabel>::LabelRefT<'a>;

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

impl<'a, EqLabelT: EqLabel> ExactSizeIterator for PackedEqLabelIter<'a, EqLabelT> {}

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
    #[allow(dead_code)]
    pub fn iter(&self) -> EqEntryIter<'_, EqLabelT> {
        EqEntryIter {
            underlying_iter: self.count_map.iter(),
            contains_ori: self.contains_ori,
        }
    }

    /// Return an iterator over the equivalence class
    /// map iterator. Unlike the `iter` method, this iterator
    /// yields the "full" key, and primarily designed for filling
    /// in the `PackedEqMap`. For the `BasicEqLabel` the "full"
    /// key is just the target ids, while for the `RangeFactorizedEqLabel`
    /// it is the target ids and the conditional probability bin ids.
    pub fn full_key_iter(&self) -> EqEntryKeyIter<'_, EqLabelT> {
        EqEntryKeyIter {
            underlying_iter: self.count_map.iter(),
            contains_ori: self.contains_ori,
        }
    }
}

pub struct EqEntryKeyIter<'a, EqLabelT: EqLabel> {
    underlying_iter: std::collections::hash_map::Iter<'a, EqLabelT, usize>,
    contains_ori: bool,
}

impl<'a, EqLabelT: EqLabel> Iterator for EqEntryKeyIter<'a, EqLabelT> {
    type Item = (&'a [u32], &'a usize);

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        match self.underlying_iter.next() {
            Some((k, v)) => Some((k.extract_key_for_packed_map(self.contains_ori), v)),
            None => None,
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.underlying_iter.size_hint()
    }
}

impl<'a, EqLabelT: EqLabel> ExactSizeIterator for EqEntryKeyIter<'a, EqLabelT> {}

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
