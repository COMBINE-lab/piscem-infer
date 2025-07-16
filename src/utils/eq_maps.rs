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

/// An iterator over the labels of the
/// `PackedEqMap`.
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

impl<'a> ExactSizeIterator for PackedEqLabelIter<'a> {}

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

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.underlying_iter.size_hint()
    }
}

impl<'a> ExactSizeIterator for EqEntryIter<'a> {}
