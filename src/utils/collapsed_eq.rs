//! Collapsed equivalence class view for EM/SQUAREM/Gibbs/bootstrap.
//!
//! Positional bins (when enabled) subdivide equivalence classes by the
//! positional bin of each fragment. Downstream filters (position-CV
//! leakage, coverage smoothing) consume this subdivision, but the
//! standard M-step does not: `m_step` / `m_step_par` read only
//! `target_labels()` and `target_probs()`, and `target_probs()` depends
//! on the range-factorization prob bin only, not on the positional bin.
//!
//! This module provides a "collapsed" variant of the range-factorized
//! equivalence class label whose byte layout is `[targets, prob_bins]`
//! (no position bins, no orientations). A `CollapsedEqMap` wrapping a
//! `PackedEqMap<CollapsedRangeFactorizedEqLabel>` is built by scanning
//! the positional packed map, stripping the positional-bin segment, and
//! summing counts of ECs that share the same `(target_set,
//! prob_bin_vector)` key.
//!
//! The collapsed label type is intentionally distinct from
//! `RangeFactorizedEqLabel` so that both maps can coexist: the
//! collapsed label hardcodes a 2-segment layout and never consults the
//! `NUM_POS_BINS` global during parsing, while the positional label
//! retains its existing layout-via-global behavior. EM/SQUAREM/Gibbs
//! and bootstrap paths see the collapsed type through the existing
//! `EqLabel` generic; no trait changes are required.

use ahash::AHashMap;

use crate::utils::eq_maps::{
    EqLabel, NUM_BINS, NUM_POS_BINS, PackedEqMap, RangeFactorizedEqLabel, TargetLabels,
    TargetLabelsRef,
};

/// A collapsed range-factorized equivalence-class label. Layout is
/// `[targets, prob_bins]` — always 2 segments, regardless of whether
/// positional binning is enabled globally.
#[derive(Hash, PartialEq, Eq)]
pub struct CollapsedRangeFactorizedEqLabel {
    pub targets_and_bins: Vec<u32>,
}

#[derive(Hash, PartialEq, Eq)]
pub struct CollapsedRangeFactorizedEqLabelRef<'a> {
    pub targets_and_bins: &'a [u32],
}

impl EqLabel for CollapsedRangeFactorizedEqLabel {
    type LabelRefT<'a> = CollapsedRangeFactorizedEqLabelRef<'a>;

    fn new_ref(labels: &[u32], _has_ori: bool) -> CollapsedRangeFactorizedEqLabelRef<'_> {
        CollapsedRangeFactorizedEqLabelRef {
            targets_and_bins: labels,
        }
    }

    fn new(labels: &[u32], probs: Option<&[f64]>, _pos_bins: Option<&[u32]>) -> Self {
        let probs = probs.expect("probs *must* be present for collapsed range factorized EC");
        let tot_prob: f64 = probs.iter().sum();
        let num_labels = probs.len();
        let num_bins = *NUM_BINS.get().unwrap() as usize;

        // Input `labels` may carry trailing orientations; take only the
        // first `num_labels` target ids and drop the rest.
        let mut targets_and_bins: Vec<u32> = labels[..num_labels].into();
        targets_and_bins.extend(probs.iter().map(|&prob| {
            let p: f64 = prob / tot_prob;
            if p >= 1.0 {
                (num_bins - 1) as u32
            } else {
                (p * num_bins as f64) as u32
            }
        }));
        Self { targets_and_bins }
    }
}

impl TargetLabels for CollapsedRangeFactorizedEqLabel {
    #[inline]
    fn target_labels(&self, _with_ori: bool) -> &[u32] {
        let l = self.targets_and_bins.len() / 2;
        &self.targets_and_bins[..l]
    }

    #[inline]
    fn extract_key_for_packed_map(&self, _has_ori: bool) -> &[u32] {
        &self.targets_and_bins
    }
}

impl<'a> TargetLabelsRef for CollapsedRangeFactorizedEqLabelRef<'a> {
    #[inline]
    fn target_labels(&self) -> &[u32] {
        let l = self.targets_and_bins.len() / 2;
        &self.targets_and_bins[..l]
    }

    #[inline]
    fn target_probs(&self) -> impl Iterator<Item = f64> {
        let l = self.targets_and_bins.len() / 2;
        let num_bins = *NUM_BINS.get().unwrap();
        let half_bin_width = 0.5 / num_bins;
        self.targets_and_bins[l..2 * l]
            .iter()
            .map(move |&x| (x as f64) / num_bins + half_bin_width)
    }

    #[inline]
    fn target_pos_bins(&self) -> Option<&[u32]> {
        None
    }
}

/// A collapsed equivalence-class map plus the index from the source
/// positional map to the collapsed map. `pos_to_collapsed[i]` is the
/// index in `packed` of the collapsed EC that absorbed source EC `i`.
pub struct CollapsedEqMap {
    pub packed: PackedEqMap<CollapsedRangeFactorizedEqLabel>,
    /// Bijection from positional-EC index → collapsed-EC index. Kept so
    /// that bootstrap/Gibbs samples (or future diagnostics) can lift
    /// per-EC statistics back to the positional granularity.
    #[allow(dead_code)]
    pub pos_to_collapsed: Vec<u32>,
}

impl CollapsedEqMap {
    #[inline]
    pub fn len(&self) -> usize {
        self.packed.len()
    }

    #[inline]
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.packed.len() == 0
    }
}

/// Build a collapsed equivalence-class view by grouping all positional
/// ECs that share the same `(targets, prob_bins)` key and summing their
/// counts. When positional binning is disabled globally, the positional
/// map is already collapsed and this is a cheap identity rebuild.
pub fn build_collapsed(pos_map: &PackedEqMap<RangeFactorizedEqLabel>) -> CollapsedEqMap {
    build_collapsed_with_target_mask(pos_map, None)
}

/// Build a collapsed equivalence-class view, optionally dropping target labels
/// whose mask entry is false. EC counts are preserved even when all targets in
/// an EC are removed, so EM initialization and projection see the same total
/// input weight as the unfiltered map.
pub fn build_collapsed_with_target_mask(
    pos_map: &PackedEqMap<RangeFactorizedEqLabel>,
    target_mask: Option<&[bool]>,
) -> CollapsedEqMap {
    let n_pos = pos_map.len();
    let has_pos = NUM_POS_BINS.get().is_some_and(|&n| n > 1.0);

    if !has_pos && target_mask.is_none() {
        // Positional map is already [targets, prob_bins]; reinterpret
        // its bytes as a collapsed-label packed map. Identity mapping.
        let packed = PackedEqMap::<CollapsedRangeFactorizedEqLabel>::from_raw(
            pos_map.eq_labels.clone(),
            pos_map.eq_label_starts.clone(),
            pos_map.counts.clone(),
            false,
        );
        let pos_to_collapsed: Vec<u32> = (0..n_pos as u32).collect();
        return CollapsedEqMap {
            packed,
            pos_to_collapsed,
        };
    }

    // General case: strip the positional-bin segment and merge duplicates.
    // With a target mask, also drop structurally removed target/bin pairs
    // without renormalizing the retained probability bins.
    let mut pos_to_collapsed: Vec<u32> = Vec::with_capacity(n_pos);
    let mut key_to_idx: AHashMap<Vec<u32>, u32> = AHashMap::with_capacity(n_pos);
    let mut collapsed_labels: Vec<u32> = Vec::with_capacity(pos_map.eq_labels.len() * 2 / 3);
    let mut collapsed_starts: Vec<u32> = Vec::with_capacity(n_pos + 1);
    let mut collapsed_counts: Vec<usize> = Vec::with_capacity(n_pos);
    let mut filtered_key: Vec<u32> = Vec::new();
    collapsed_starts.push(0);

    for (i, count) in pos_map.counts.iter().enumerate() {
        let s = pos_map.eq_label_starts[i] as usize;
        let e = pos_map.eq_label_starts[i + 1] as usize;
        let full = &pos_map.eq_labels[s..e];
        let key = if let Some(mask) = target_mask {
            filtered_key.clear();
            let n_segments = if has_pos { 3 } else { 2 };
            let n = full.len() / n_segments;
            filtered_key.reserve(2 * n);
            for &target in &full[..n] {
                if mask.get(target as usize).copied().unwrap_or(false) {
                    filtered_key.push(target);
                }
            }
            for (&target, &bin) in full[..n].iter().zip(full[n..2 * n].iter()) {
                if mask.get(target as usize).copied().unwrap_or(false) {
                    filtered_key.push(bin);
                }
            }
            filtered_key.as_slice()
        } else {
            // 3-segment layout [targets, prob_bins, pos_bins] -> key length 2/3.
            let n = full.len() / 3;
            &full[..2 * n]
        };

        if let Some(&idx) = key_to_idx.get(key) {
            pos_to_collapsed.push(idx);
            collapsed_counts[idx as usize] += count;
        } else {
            let idx = collapsed_counts.len() as u32;
            key_to_idx.insert(key.to_vec(), idx);
            pos_to_collapsed.push(idx);
            collapsed_labels.extend_from_slice(key);
            collapsed_starts.push(collapsed_labels.len() as u32);
            collapsed_counts.push(*count);
        }
    }

    let packed = PackedEqMap::<CollapsedRangeFactorizedEqLabel>::from_raw(
        collapsed_labels,
        collapsed_starts,
        collapsed_counts,
        false,
    );
    CollapsedEqMap {
        packed,
        pos_to_collapsed,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::eq_maps::{EqMap, OrientationProperty, PackedEqMap};

    // Tests in this module share a single process-wide OnceLock, so we
    // fix NUM_BINS and NUM_POS_BINS to constants usable by every test in
    // the crate. `.ok()` swallows the "already set" error when another
    // test already installed the same values.
    const TEST_NUM_BINS: f64 = 10.0;
    const TEST_NUM_POS_BINS: f64 = 5.0;

    fn init_globals() {
        NUM_BINS.set(TEST_NUM_BINS).ok();
        NUM_POS_BINS.set(TEST_NUM_POS_BINS).ok();
    }

    /// Build a small positional PackedEqMap with three positional ECs
    /// that share targets and prob_bins but differ in pos_bins. The
    /// collapsed view should merge them into one EC with summed counts.
    #[test]
    fn collapse_merges_ecs_differing_only_by_pos_bins() {
        init_globals();

        let mut eqm =
            EqMap::<RangeFactorizedEqLabel>::new(OrientationProperty::OrientationAgnostic);

        // Three ECs all with targets = [3, 7] and prob_bins that will
        // hash to the same pair; they differ only in pos_bins.
        // Probabilities chosen so prob_bin = floor(p * 10).
        // (0.3, 0.7) → prob_bins [3, 7]; use pos_bins [0,0], [1,2], [4,4].
        let labels = [3u32, 7u32];
        let probs = [0.3_f64, 0.7_f64];

        let lab_a = RangeFactorizedEqLabel::new(&labels, Some(&probs), Some(&[0u32, 0u32]));
        let lab_b = RangeFactorizedEqLabel::new(&labels, Some(&probs), Some(&[1u32, 2u32]));
        let lab_c = RangeFactorizedEqLabel::new(&labels, Some(&probs), Some(&[4u32, 4u32]));

        *eqm.add(lab_a) = 5;
        *eqm.add(lab_b) = 3;
        *eqm.add(lab_c) = 2;

        let pos_packed = PackedEqMap::<RangeFactorizedEqLabel>::from_eq_map(&eqm);
        assert_eq!(pos_packed.len(), 3);

        let collapsed = build_collapsed(&pos_packed);
        assert_eq!(collapsed.len(), 1, "all 3 positional ECs share key");
        assert_eq!(
            collapsed.packed.counts[0], 10,
            "counts should sum: 5 + 3 + 2"
        );
        assert_eq!(
            collapsed.pos_to_collapsed,
            vec![0, 0, 0],
            "all source ECs map to collapsed idx 0"
        );
    }

    /// ECs that differ in target set or prob bins must NOT collapse.
    #[test]
    fn collapse_preserves_distinct_keys() {
        init_globals();

        let mut eqm =
            EqMap::<RangeFactorizedEqLabel>::new(OrientationProperty::OrientationAgnostic);
        let lab_a = RangeFactorizedEqLabel::new(&[1u32, 2u32], Some(&[0.5, 0.5]), Some(&[0, 0]));
        let lab_b = RangeFactorizedEqLabel::new(&[1u32, 3u32], Some(&[0.5, 0.5]), Some(&[0, 0]));
        let lab_c = RangeFactorizedEqLabel::new(&[1u32, 2u32], Some(&[0.2, 0.8]), Some(&[0, 0]));

        *eqm.add(lab_a) = 4;
        *eqm.add(lab_b) = 2;
        *eqm.add(lab_c) = 1;

        let pos_packed = PackedEqMap::<RangeFactorizedEqLabel>::from_eq_map(&eqm);
        let collapsed = build_collapsed(&pos_packed);

        assert_eq!(collapsed.len(), 3, "distinct keys stay distinct");
        let total: usize = collapsed.packed.counts.iter().sum();
        assert_eq!(total, 7, "total count preserved");
    }

    /// Exercise the serial M-step on both the positional and collapsed
    /// maps starting from the same transcript abundance vector and
    /// verify the resulting counts match exactly. This is the core
    /// correctness contract: EM on collapsed = EM on positional.
    #[test]
    fn m_step_bit_identical_serial() {
        init_globals();

        let mut eqm =
            EqMap::<RangeFactorizedEqLabel>::new(OrientationProperty::OrientationAgnostic);

        // A richer set: transcripts 0..5; four distinct target/prob
        // keys, each split across 2-3 positional bins.
        let t1 = [0u32, 1u32, 2u32];
        let p1 = [0.2_f64, 0.3, 0.5];
        *eqm.add(RangeFactorizedEqLabel::new(
            &t1,
            Some(&p1),
            Some(&[0, 0, 0]),
        )) = 6;
        *eqm.add(RangeFactorizedEqLabel::new(
            &t1,
            Some(&p1),
            Some(&[1, 2, 3]),
        )) = 4;
        *eqm.add(RangeFactorizedEqLabel::new(
            &t1,
            Some(&p1),
            Some(&[4, 4, 4]),
        )) = 2;

        let t2 = [1u32, 3u32];
        let p2 = [0.4_f64, 0.6];
        *eqm.add(RangeFactorizedEqLabel::new(&t2, Some(&p2), Some(&[0, 0]))) = 5;
        *eqm.add(RangeFactorizedEqLabel::new(&t2, Some(&p2), Some(&[2, 2]))) = 7;

        let t3 = [2u32, 4u32];
        let p3 = [0.7_f64, 0.3];
        *eqm.add(RangeFactorizedEqLabel::new(&t3, Some(&p3), Some(&[0, 0]))) = 3;

        let t4 = [4u32];
        let p4 = [1.0_f64];
        *eqm.add(RangeFactorizedEqLabel::new(&t4, Some(&p4), Some(&[0]))) = 8;

        let pos_packed = PackedEqMap::<RangeFactorizedEqLabel>::from_eq_map(&eqm);
        let collapsed = build_collapsed(&pos_packed);

        // Non-uniform prior to avoid trivial uniform propagation.
        let prev: Vec<f64> = vec![10.0, 5.0, 2.0, 8.0, 1.0];
        let inv_eff_lens: Vec<f64> = vec![0.01, 0.02, 0.015, 0.025, 0.03];
        let ntx = prev.len();

        let mut out_pos = vec![0.0_f64; ntx];
        serial_m_step(&pos_packed, &prev, &inv_eff_lens, &mut out_pos);

        let mut out_col = vec![0.0_f64; ntx];
        serial_m_step(&collapsed.packed, &prev, &inv_eff_lens, &mut out_col);

        for (i, (a, b)) in out_pos.iter().zip(out_col.iter()).enumerate() {
            assert!(
                (a - b).abs() < 1e-12,
                "transcript {i}: positional={a} collapsed={b}"
            );
        }

        // Total mapped reads must also match across the two maps.
        let total_pos: f64 = out_pos.iter().sum();
        let total_col: f64 = out_col.iter().sum();
        assert!((total_pos - total_col).abs() < 1e-12);
    }

    /// Dropping masked targets from the collapsed labels should be equivalent
    /// to keeping the full EC map but giving those targets zero effective
    /// length: they contribute zero weight to every denominator.
    #[test]
    fn masked_collapse_matches_zero_effective_length() {
        init_globals();

        let mut eqm =
            EqMap::<RangeFactorizedEqLabel>::new(OrientationProperty::OrientationAgnostic);

        *eqm.add(RangeFactorizedEqLabel::new(
            &[0u32, 1u32, 2u32],
            Some(&[0.2, 0.3, 0.5]),
            Some(&[0, 1, 2]),
        )) = 11;
        *eqm.add(RangeFactorizedEqLabel::new(
            &[1u32, 3u32],
            Some(&[0.4, 0.6]),
            Some(&[2, 3]),
        )) = 7;
        *eqm.add(RangeFactorizedEqLabel::new(
            &[2u32],
            Some(&[1.0]),
            Some(&[4]),
        )) = 5;

        let pos_packed = PackedEqMap::<RangeFactorizedEqLabel>::from_eq_map(&eqm);
        let keep = vec![true, false, true, false];
        let collapsed = build_collapsed_with_target_mask(&pos_packed, Some(&keep));

        let prev = vec![10.0, 5.0, 2.0, 8.0];
        let inv_eff_lens = vec![0.01, 0.0, 0.015, 0.0];

        let mut out_pos = vec![0.0_f64; prev.len()];
        serial_m_step(&pos_packed, &prev, &inv_eff_lens, &mut out_pos);

        let mut out_col = vec![0.0_f64; prev.len()];
        serial_m_step(&collapsed.packed, &prev, &inv_eff_lens, &mut out_col);

        for (i, (a, b)) in out_pos.iter().zip(out_col.iter()).enumerate() {
            assert!(
                (a - b).abs() < 1e-12,
                "transcript {i}: positional_zero_eff={a} masked_collapsed={b}"
            );
        }
        assert_eq!(collapsed.packed.counts.iter().sum::<usize>(), 23);
    }

    /// Local replica of the serial M-step in `em.rs:m_step`. We inline
    /// it here (rather than call the crate one) so this test stays a
    /// pure unit test of the collapse invariant without dragging in
    /// EMInfo construction.
    fn serial_m_step<EqLabelT>(
        packed: &PackedEqMap<EqLabelT>,
        prev: &[f64],
        inv_eff_lens: &[f64],
        curr: &mut [f64],
    ) where
        EqLabelT: EqLabel,
    {
        let mut weights: Vec<f64> = Vec::with_capacity(64);
        for (k, v) in packed.iter_labels().zip(packed.counts.iter()) {
            let count = *v as f64;
            let mut denom = 0.0_f64;
            for (e, cp) in k.target_labels().iter().zip(k.target_probs()) {
                let w = cp * prev[*e as usize] * inv_eff_lens[*e as usize];
                weights.push(w);
                denom += w;
            }
            if denom > 1e-8 {
                let c_over_d = count / denom;
                for (tid, w) in k.target_labels().iter().zip(weights.iter()) {
                    curr[*tid as usize] += c_over_d * w;
                }
            }
            weights.clear();
        }
    }
}
