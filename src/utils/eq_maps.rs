use ahash::AHashMap;

/// aliases for referring to the different basic types of eq class maps
pub type BasicEqMap = EqMap<BasicEqLabel>;
pub type RangeFactorizedEqMap = EqMap<RangeFactorizedEqLabel>;

/// the default number of bins to use for range-factorized equivalence classes
pub static NUM_BINS: std::sync::OnceLock<f64> = std::sync::OnceLock::new();

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
    fn new(labels: &[u32], probs: Option<&[f64]>) -> Self;
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
    fn new_ref(labels: &[u32], has_ori: bool) -> BasicEqLabelRef {
        BasicEqLabelRef {
            targets: labels,
            contains_ori: has_ori,
        }
    }

    /// we ignore probabiltiies in the basic equivalence class, and just
    /// pass along the targets we are given as the labels
    fn new(targets: &[u32], _probs: Option<&[f64]>) -> Self {
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

    fn new_ref(labels: &[u32], has_ori: bool) -> RangeFactorizedEqLabelRef {
        RangeFactorizedEqLabelRef {
            targets_and_bins: labels,
            contains_ori: has_ori,
        }
    }

    /// create the range factorized equivalence class the labels and
    /// the associated probabilities.
    fn new(labels: &[u32], probs: Option<&[f64]>) -> Self {
        // the probability vector must not be `None` when creating a
        // range-factorized equivalence class
        let probs = probs.expect("probs *must* be present for range factorized equivalence class");

        // compute the total probability so we can normalize later
        let tot_prob: f64 = probs.iter().sum();
        // get the number of labels (because the labels slice could contain either
        // just labels or labels and encoded orientations)
        let num_labels = probs.len();
        // get the number of bins we are using for the range factorized equivalence class
        // representation
        let num_bins = *NUM_BINS.get().unwrap() as usize;

        // we will structure the contents of the range factorized equivalence class, that
        // contains k labels as
        //
        // [label_1],[label_2],...,[label_k],[bin_id_1],[bin_id_2],...,[bin_id_k],[ori_1],[ori_2],...,[ori_k]
        //
        // where the `ori` components may be optional
        let (just_labels, oris) = labels.split_at(num_labels);
        let mut targets_and_bins: Vec<u32> = just_labels.into();
        // map the probabilities to the appropriate bin IDs
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
        // add back the encoding of the orientations if we have them
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
        // if this is an orientation-aware equivalence class factorization, then
        // the targets_and_bins vector contains (labels, oris, probs), otherwise
        // it contains just (labels, probs)
        let l = self.targets_and_bins.len() / if with_ori { 3 } else { 2 };
        &self.targets_and_bins[..l]
    }

    /// NOTE: The 2 * l below here is a huge hack. We want to include the bins in the
    /// contents of the packed map, and this is how we do it, but we certainly should
    /// find a more elegant and less hacky way.

    fn extract_key_for_packed_map(&self, has_ori: bool) -> &[u32] {
        let l = self.targets_and_bins.len() / if has_ori { 3 } else { 2 };
        &self.targets_and_bins[..2 * l]
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
        let with_ori = self.contains_ori;
        let l = self.targets_and_bins.len() / if with_ori { 3 } else { 2 };
        &self.targets_and_bins[..l]
    }

    /// returns an iterator over the probabilities assocaited with the
    /// bins of this reference
    #[inline]
    fn target_probs(&self) -> impl Iterator<Item = f64> {
        let with_ori = self.contains_ori;
        let l = self.targets_and_bins.len() / if with_ori { 3 } else { 2 };
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

    pub fn iter_labels(&self) -> PackedEqLabelIter<EqLabelT> {
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
    pub fn iter(&self) -> EqEntryIter<EqLabelT> {
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
    pub fn full_key_iter(&self) -> EqEntryKeyIter<EqLabelT> {
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
