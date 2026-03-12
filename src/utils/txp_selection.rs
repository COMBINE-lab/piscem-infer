//! Transcript variable selection via EC graph analysis.
//!
//! Exploits the geometry of the equivalence class compatibility matrix
//! to identify structurally redundant transcripts before EM.
//!
//! Pipeline:
//!   Stage 1: Signature collapsing — group transcripts with identical EC neighborhoods
//!   Stage 2: Unique EC peeling — iteratively identify required transcripts
//!   Stage 3: Subset dominance — remove transcripts whose EC support is a subset of another's

#![allow(dead_code)]

use ahash::AHashMap;
use serde::Serialize;
use tracing::info;

use crate::utils::eq_maps::{EqLabel, PackedEqMap, TargetLabelsRef};

/// Configuration for which selection stages to run.
#[derive(Clone, Debug, Serialize)]
pub struct SelectionStages {
    /// Stage 1: Collapse transcripts with identical EC signatures into groups.
    /// Always runs (required for the pipeline to function).
    pub collapse: bool,
    /// Stage 2: Iteratively peel degree-1 EQ classes to find required groups.
    pub peeling: bool,
    /// Stage 3: Remove groups whose EC signature is a subset of another's.
    pub dominance: bool,
}

impl Default for SelectionStages {
    fn default() -> Self {
        Self {
            collapse: true,
            peeling: true,
            dominance: true,
        }
    }
}

impl SelectionStages {
    /// All stages enabled (default).
    pub fn all() -> Self {
        Self::default()
    }

    /// Only collapse + peeling (no subset dominance).
    pub fn no_dominance() -> Self {
        Self {
            collapse: true,
            peeling: true,
            dominance: false,
        }
    }

    /// Parse from a comma-separated string like "collapse,peeling,dominance".
    pub fn from_str_list(s: &str) -> Result<Self, String> {
        let mut stages = Self {
            collapse: false,
            peeling: false,
            dominance: false,
        };
        for part in s.split(',') {
            match part.trim().to_lowercase().as_str() {
                "collapse" => stages.collapse = true,
                "peeling" => stages.peeling = true,
                "dominance" => stages.dominance = true,
                other => return Err(format!("unknown selection stage: '{other}'. Valid stages: collapse, peeling, dominance")),
            }
        }
        // Collapse is always required
        stages.collapse = true;
        Ok(stages)
    }
}

/// Packed reverse index mapping each transcript to its EQ class neighborhood.
///
/// Uses flat arrays with offsets (same layout as PackedEqMap) for cache locality.
/// `eqc_ids[offsets[t]..offsets[t+1]]` gives the sorted EQ class indices for transcript t.
#[allow(dead_code)]
pub struct TranscriptEqIndex {
    /// Flat packed list of EQ class indices, sorted per transcript.
    pub eqc_ids: Vec<u32>,
    /// Offset array: transcript t's EQ classes span eqc_ids[offsets[t]..offsets[t+1]].
    /// Length = num_targets + 1.
    pub offsets: Vec<u32>,
}

impl TranscriptEqIndex {
    /// Build reverse index from a PackedEqMap.
    ///
    /// Time: O(E) where E = total EC-transcript edges.
    pub fn from_packed_eq_map<EqLabelT: EqLabel>(
        packed_map: &PackedEqMap<EqLabelT>,
        num_targets: usize,
    ) -> Self {
        // First pass: count how many EQ classes each transcript belongs to.
        let mut counts = vec![0u32; num_targets];
        for eqc_idx in 0..packed_map.len() {
            let label = packed_map.refs_for_eqc(eqc_idx);
            for &tid in label.target_labels() {
                counts[tid as usize] += 1;
            }
        }

        // Build offset array from counts.
        let mut offsets = Vec::with_capacity(num_targets + 1);
        offsets.push(0u32);
        for &c in &counts {
            offsets.push(offsets.last().unwrap() + c);
        }
        let total_edges = *offsets.last().unwrap() as usize;

        // Second pass: fill in EQ class IDs using write cursors.
        let mut eqc_ids = vec![0u32; total_edges];
        let mut cursors = vec![0u32; num_targets]; // current write position per transcript
        for eqc_idx in 0..packed_map.len() {
            let label = packed_map.refs_for_eqc(eqc_idx);
            for &tid in label.target_labels() {
                let t = tid as usize;
                let pos = offsets[t] + cursors[t];
                eqc_ids[pos as usize] = eqc_idx as u32;
                cursors[t] += 1;
            }
        }

        // Sort each transcript's EQ class list (needed for signature comparison).
        for t in 0..num_targets {
            let s = offsets[t] as usize;
            let e = offsets[t + 1] as usize;
            eqc_ids[s..e].sort_unstable();
        }

        Self { eqc_ids, offsets }
    }

    /// Returns the number of transcripts in this index.
    #[inline]
    pub fn num_targets(&self) -> usize {
        self.offsets.len() - 1
    }

    /// Returns the sorted EQ class signature for transcript `t`.
    #[inline]
    pub fn signature(&self, t: usize) -> &[u32] {
        let s = self.offsets[t] as usize;
        let e = self.offsets[t + 1] as usize;
        &self.eqc_ids[s..e]
    }

    /// Returns the number of EQ classes containing transcript `t`.
    #[inline]
    pub fn degree(&self, t: usize) -> usize {
        (self.offsets[t + 1] - self.offsets[t]) as usize
    }
}

/// Result of signature collapsing (Stage 1).
///
/// Packed representation: each group has a representative and members stored
/// in flat arrays with offsets. Groups are sorted by representative ID.
#[allow(dead_code)]
pub struct SignatureGroups {
    /// For each group, the representative transcript ID.
    pub representatives: Vec<u32>,
    /// Flat packed member transcript IDs (including the representative).
    pub members: Vec<u32>,
    /// Offset array: group i's members span members[offsets[i]..offsets[i+1]].
    /// Length = num_groups + 1.
    pub offsets: Vec<u32>,
}

impl SignatureGroups {
    /// Number of groups.
    #[inline]
    pub fn num_groups(&self) -> usize {
        self.representatives.len()
    }

    /// Returns the member transcript IDs for group `i`.
    #[inline]
    pub fn group_members(&self, i: usize) -> &[u32] {
        let s = self.offsets[i] as usize;
        let e = self.offsets[i + 1] as usize;
        &self.members[s..e]
    }
}

/// Stage 1: Group transcripts with identical EC signatures.
///
/// Two transcripts are indistinguishable if they appear in exactly the same
/// set of equivalence classes. We hash each transcript's sorted EQ class list
/// and group by hash, then verify equality.
///
/// Returns the collapsed groups. The representative of each group is the
/// member with the smallest transcript ID.
pub fn signature_collapse(index: &TranscriptEqIndex) -> SignatureGroups {
    let num_targets = index.num_targets();

    // Hash each transcript's signature → group by hash.
    // Use the signature slice directly as the map key (via a u64 hash for speed).
    let mut groups: AHashMap<&[u32], Vec<u32>> = AHashMap::new();

    for t in 0..num_targets {
        let sig = index.signature(t);
        // Skip transcripts with no EQ classes (they can't be quantified).
        if sig.is_empty() {
            continue;
        }
        groups.entry(sig).or_default().push(t as u32);
    }

    // Build packed representation, sorted by representative (smallest ID in each group).
    let mut group_list: Vec<Vec<u32>> = groups.into_values().collect();
    for group in &mut group_list {
        group.sort_unstable();
    }
    group_list.sort_unstable_by_key(|g| g[0]);

    let mut representatives = Vec::with_capacity(group_list.len());
    let mut members = Vec::new();
    let mut offsets = Vec::with_capacity(group_list.len() + 1);
    offsets.push(0u32);

    for group in &group_list {
        representatives.push(group[0]);
        members.extend_from_slice(group);
        offsets.push(members.len() as u32);
    }

    SignatureGroups {
        representatives,
        members,
        offsets,
    }
}

/// Stage 2: Unique EC peeling — iteratively identify required transcripts.
///
/// A transcript (or signature group) is "required" if it is the sole member
/// of some equivalence class. After marking it required and removing its edges,
/// other EQ classes may become degree-1, cascading.
///
/// Operates on signature groups (treats each group as one node).
///
/// Returns a boolean mask over groups: `true` = required (must keep).
pub fn unique_ec_peeling<EqLabelT: EqLabel>(
    packed_map: &PackedEqMap<EqLabelT>,
    groups: &SignatureGroups,
    num_targets: usize,
) -> Vec<bool> {
    let num_eqcs = packed_map.len();
    let num_groups = groups.num_groups();

    // Map transcript ID → group index.
    let mut tx_to_group = vec![u32::MAX; num_targets];
    for (gi, &rep) in groups.representatives.iter().enumerate() {
        for &member in groups.group_members(gi) {
            tx_to_group[member as usize] = gi as u32;
        }
        // Ensure representative is also mapped (it's in members too,
        // but let's be explicit).
        tx_to_group[rep as usize] = gi as u32;
    }

    // Compute degree of each EQ class in group-space.
    // degree[k] = number of distinct groups in EQ class k.
    let mut eqc_degree = vec![0u32; num_eqcs];
    // Also build EQ class → group adjacency (packed).
    let mut eqc_groups: Vec<Vec<u32>> = vec![Vec::new(); num_eqcs];

    for eqc_idx in 0..num_eqcs {
        let label = packed_map.refs_for_eqc(eqc_idx);
        let mut seen_groups: Vec<u32> = label
            .target_labels()
            .iter()
            .map(|&tid| tx_to_group[tid as usize])
            .filter(|&g| g != u32::MAX)
            .collect();
        seen_groups.sort_unstable();
        seen_groups.dedup();
        eqc_degree[eqc_idx] = seen_groups.len() as u32;
        eqc_groups[eqc_idx] = seen_groups;
    }

    // Build group → EQ class reverse adjacency.
    // (We need this to know which EQ classes to update when removing a group.)
    let mut group_to_eqcs: Vec<Vec<u32>> = vec![Vec::new(); num_groups];
    for (eqc_idx, grps) in eqc_groups.iter().enumerate() {
        for &g in grps {
            group_to_eqcs[g as usize].push(eqc_idx as u32);
        }
    }

    let mut required = vec![false; num_groups];
    let mut removed = vec![false; num_groups];

    // Initialize worklist with groups that have a degree-1 EQ class.
    let mut worklist: Vec<u32> = Vec::new();
    for (eqc_idx, &deg) in eqc_degree.iter().enumerate() {
        if deg == 1 && !eqc_groups[eqc_idx].is_empty() {
            let g = eqc_groups[eqc_idx][0];
            if !required[g as usize] {
                required[g as usize] = true;
                worklist.push(g);
            }
        }
    }

    // Peel: remove required groups and propagate.
    while let Some(g) = worklist.pop() {
        removed[g as usize] = true;

        // Decrement degree of all EQ classes this group belongs to.
        for &eqc_idx in &group_to_eqcs[g as usize] {
            let eqc = eqc_idx as usize;
            if eqc_degree[eqc] > 0 {
                eqc_degree[eqc] -= 1;
            }
            // If this EQ class just became degree-1, mark the remaining group.
            if eqc_degree[eqc] == 1 {
                for &remaining_g in &eqc_groups[eqc] {
                    if !removed[remaining_g as usize] && !required[remaining_g as usize] {
                        required[remaining_g as usize] = true;
                        worklist.push(remaining_g);
                    }
                }
            }
        }
    }

    required
}

/// Stage 3: Subset dominance test.
///
/// Group i is redundant if its EQ signature is a subset of group j's signature
/// (for some j ≠ i). In that case, group j can always explain group i's reads.
///
/// Operates on signature groups via the TranscriptEqIndex (using representatives'
/// signatures, which are identical for all group members).
///
/// Returns a boolean mask over groups: `true` = dominated (can be removed).
pub fn subset_dominance(
    index: &TranscriptEqIndex,
    groups: &SignatureGroups,
) -> Vec<bool> {
    let num_groups = groups.num_groups();
    let mut dominated = vec![false; num_groups];

    // For each group, get the representative's signature.
    // Sort groups by signature length (ascending) — smaller signatures are
    // more likely to be subsets.
    let mut group_order: Vec<usize> = (0..num_groups).collect();
    group_order.sort_unstable_by_key(|&gi| {
        index.degree(groups.representatives[gi] as usize)
    });

    // Build an inverted index: EQ class → list of groups containing it.
    // This allows us to quickly find candidate supersets.
    let max_eqc = index
        .eqc_ids
        .iter()
        .copied()
        .max()
        .map(|m| m as usize + 1)
        .unwrap_or(0);
    let mut eqc_to_groups: Vec<Vec<u32>> = vec![Vec::new(); max_eqc];
    for gi in 0..num_groups {
        let rep = groups.representatives[gi] as usize;
        for &eqc in index.signature(rep) {
            eqc_to_groups[eqc as usize].push(gi as u32);
        }
    }

    // For each group (smallest signature first), check if any other group's
    // signature is a superset.
    for &gi in &group_order {
        if dominated[gi] {
            continue;
        }
        let rep_i = groups.representatives[gi] as usize;
        let sig_i = index.signature(rep_i);
        if sig_i.is_empty() {
            continue;
        }
        let deg_i = sig_i.len();

        // Candidate supersets: groups that share the first EQ class with gi
        // and have degree >= deg_i.
        // We check all groups sharing any EQ class, then verify the full subset relation.
        // Use the rarest EQ class (fewest groups) for initial candidate generation.
        let rarest_eqc = sig_i
            .iter()
            .min_by_key(|&&eqc| eqc_to_groups[eqc as usize].len())
            .unwrap();

        for &candidate_gj in &eqc_to_groups[*rarest_eqc as usize] {
            let gj = candidate_gj as usize;
            if gj == gi || dominated[gj] {
                continue;
            }
            let rep_j = groups.representatives[gj] as usize;
            let sig_j = index.signature(rep_j);
            if sig_j.len() <= deg_i {
                // Can't be a strict superset if same size or smaller.
                // (Same size = same signature = same group, already handled by Stage 1.)
                continue;
            }

            // Check if sig_i ⊆ sig_j (both sorted).
            if is_sorted_subset(sig_i, sig_j) {
                dominated[gi] = true;
                break;
            }
        }
    }

    dominated
}

/// Check if sorted slice `a` is a subset of sorted slice `b`.
/// Both must be sorted in ascending order.
fn is_sorted_subset(a: &[u32], b: &[u32]) -> bool {
    if a.len() > b.len() {
        return false;
    }
    let mut bi = 0;
    for &val in a {
        // Advance b until we find val or pass it.
        while bi < b.len() && b[bi] < val {
            bi += 1;
        }
        if bi >= b.len() || b[bi] != val {
            return false;
        }
        bi += 1;
    }
    true
}

/// Result of the full variable selection pipeline.
#[allow(dead_code)]
pub struct SelectionResult {
    /// Per-transcript mask: `true` = keep, `false` = removed.
    pub keep_mask: Vec<bool>,
    /// The signature groups (for reporting/output).
    pub groups: SignatureGroups,
    /// Per-group: `true` = required (from peeling).
    pub required: Vec<bool>,
    /// Per-group: `true` = dominated (from subset test).
    pub dominated: Vec<bool>,
    /// Number of transcripts kept.
    pub num_kept: usize,
    /// Number of transcripts removed.
    pub num_removed: usize,
}

/// Run the full variable selection pipeline (Stages 1-3).
///
/// Returns a per-transcript keep/remove mask and diagnostic information.
pub fn run_selection<EqLabelT: EqLabel>(
    packed_map: &PackedEqMap<EqLabelT>,
    num_targets: usize,
) -> SelectionResult {
    run_selection_with_stages(packed_map, num_targets, &SelectionStages::all())
}

/// Run variable selection with configurable stages.
pub fn run_selection_with_stages<EqLabelT: EqLabel>(
    packed_map: &PackedEqMap<EqLabelT>,
    num_targets: usize,
    stages: &SelectionStages,
) -> SelectionResult {
    // Stage 1: Build reverse index and collapse signatures.
    let index = TranscriptEqIndex::from_packed_eq_map(packed_map, num_targets);
    let groups = signature_collapse(&index);
    let n_no_eqc = (0..num_targets)
        .filter(|&t| index.degree(t) == 0)
        .count();
    info!(
        "  Stage 1 (signature collapse): {} transcripts → {} groups ({} with no EQ classes)",
        num_targets,
        groups.num_groups(),
        n_no_eqc
    );

    // Stage 2: Unique EC peeling.
    let required = if stages.peeling {
        let req = unique_ec_peeling(packed_map, &groups, num_targets);
        let n_required = req.iter().filter(|&&r| r).count();
        info!(
            "  Stage 2 (unique EC peeling): {} / {} groups are structurally required",
            n_required,
            groups.num_groups()
        );
        req
    } else {
        info!("  Stage 2 (unique EC peeling): skipped");
        vec![false; groups.num_groups()]
    };

    // Stage 3: Subset dominance.
    let dominated = if stages.dominance {
        let dom = subset_dominance(&index, &groups);
        let n_dominated = dom.iter().filter(|&&d| d).count();
        info!(
            "  Stage 3 (subset dominance): {} / {} groups are dominated (removable)",
            n_dominated,
            groups.num_groups()
        );
        dom
    } else {
        info!("  Stage 3 (subset dominance): skipped");
        vec![false; groups.num_groups()]
    };

    // Build per-transcript keep mask.
    // When dominance is enabled: keep groups that are NOT dominated.
    // When dominance is skipped: all groups with EQ classes are kept
    //   (the mask removes only transcripts with no EQ classes).
    let mut keep_mask = vec![false; num_targets];
    for (gi, &is_dominated) in dominated.iter().enumerate() {
        if !is_dominated {
            for &member in groups.group_members(gi) {
                keep_mask[member as usize] = true;
            }
        }
    }

    let num_kept = keep_mask.iter().filter(|&&k| k).count();
    let num_removed = num_targets - num_kept;
    info!(
        "  Variable selection result: {} kept, {} removed ({:.1}% reduction)",
        num_kept,
        num_removed,
        100.0 * num_removed as f64 / num_targets as f64
    );

    SelectionResult {
        keep_mask,
        groups,
        required,
        dominated,
        num_kept,
        num_removed,
    }
}

/// Merge multiple per-sample TranscriptEqIndexes into a union index.
///
/// The union signature for a transcript is the union of its per-sample
/// EC neighborhoods. EQ class IDs are remapped to a global namespace
/// (sample 0's EQCs: 0..n0, sample 1's: n0..n0+n1, etc.).
///
/// This gives maximum distinguishing power: two transcripts that look
/// identical in one sample may be distinguishable in another.
pub fn merge_transcript_indices(
    indices: &[TranscriptEqIndex],
    eqc_counts: &[usize], // number of EQ classes per sample
    num_targets: usize,
) -> TranscriptEqIndex {
    // Compute EQ class ID offsets per sample.
    let mut eqc_offsets = Vec::with_capacity(indices.len());
    let mut cumulative = 0u32;
    for &count in eqc_counts {
        eqc_offsets.push(cumulative);
        cumulative += count as u32;
    }

    // First pass: compute total degree per transcript.
    let mut total_degree = vec![0u32; num_targets];
    for index in indices {
        for (t, deg) in total_degree.iter_mut().enumerate() {
            *deg += index.degree(t) as u32;
        }
    }

    // Build offset array.
    let mut offsets = Vec::with_capacity(num_targets + 1);
    offsets.push(0u32);
    for &d in &total_degree {
        offsets.push(offsets.last().unwrap() + d);
    }
    let total_edges = *offsets.last().unwrap() as usize;

    // Second pass: fill in remapped EQ class IDs.
    let mut eqc_ids = vec![0u32; total_edges];
    let mut cursors = vec![0u32; num_targets];
    for (sample_idx, index) in indices.iter().enumerate() {
        let offset = eqc_offsets[sample_idx];
        for t in 0..num_targets {
            for &eqc in index.signature(t) {
                let pos = (offsets[t] + cursors[t]) as usize;
                eqc_ids[pos] = eqc + offset;
                cursors[t] += 1;
            }
        }
    }

    // Sort each transcript's merged signature.
    for t in 0..num_targets {
        let s = offsets[t] as usize;
        let e = offsets[t + 1] as usize;
        eqc_ids[s..e].sort_unstable();
    }

    TranscriptEqIndex { eqc_ids, offsets }
}

/// Peeling algorithm that works from TranscriptEqIndex + SignatureGroups,
/// without needing a PackedEqMap. Derives EQ class → group adjacency from
/// the reverse index.
///
/// Returns a boolean mask over groups: `true` = required.
pub fn unique_ec_peeling_from_index(
    index: &TranscriptEqIndex,
    groups: &SignatureGroups,
    num_eqcs: usize,
) -> Vec<bool> {
    let num_groups = groups.num_groups();

    // Build group → EQ class mapping from representatives' signatures.
    // Also build the reverse: EQ class → set of groups.
    let mut eqc_groups: Vec<Vec<u32>> = vec![Vec::new(); num_eqcs];
    for gi in 0..num_groups {
        let rep = groups.representatives[gi] as usize;
        for &eqc in index.signature(rep) {
            eqc_groups[eqc as usize].push(gi as u32);
        }
    }

    // Deduplicate (a group may appear multiple times if multiple members map to same EQ)
    for grps in &mut eqc_groups {
        grps.sort_unstable();
        grps.dedup();
    }

    let mut eqc_degree: Vec<u32> = eqc_groups.iter().map(|g| g.len() as u32).collect();

    // Build group → EQ class reverse adjacency.
    let mut group_to_eqcs: Vec<Vec<u32>> = vec![Vec::new(); num_groups];
    for (eqc_idx, grps) in eqc_groups.iter().enumerate() {
        for &g in grps {
            group_to_eqcs[g as usize].push(eqc_idx as u32);
        }
    }

    let mut required = vec![false; num_groups];
    let mut removed = vec![false; num_groups];
    let mut worklist: Vec<u32> = Vec::new();

    // Initialize with degree-1 EQ classes.
    for (eqc_idx, &deg) in eqc_degree.iter().enumerate() {
        if deg == 1 && !eqc_groups[eqc_idx].is_empty() {
            let g = eqc_groups[eqc_idx][0];
            if !required[g as usize] {
                required[g as usize] = true;
                worklist.push(g);
            }
        }
    }

    // Peel.
    while let Some(g) = worklist.pop() {
        removed[g as usize] = true;
        for &eqc_idx in &group_to_eqcs[g as usize] {
            let eqc = eqc_idx as usize;
            if eqc_degree[eqc] > 0 {
                eqc_degree[eqc] -= 1;
            }
            if eqc_degree[eqc] == 1 {
                for &remaining_g in &eqc_groups[eqc] {
                    if !removed[remaining_g as usize] && !required[remaining_g as usize] {
                        required[remaining_g as usize] = true;
                        worklist.push(remaining_g);
                    }
                }
            }
        }
    }

    required
}

/// Run variable selection from a pre-built TranscriptEqIndex (e.g., a merged index).
///
/// This variant doesn't need a PackedEqMap — it derives all structure from the index.
/// `num_eqcs` is the total number of distinct EQ classes in the index namespace.
pub fn run_selection_from_index(
    index: &TranscriptEqIndex,
    num_targets: usize,
    num_eqcs: usize,
) -> SelectionResult {
    run_selection_from_index_with_stages(index, num_targets, num_eqcs, &SelectionStages::all())
}

/// Run variable selection from a pre-built index with configurable stages.
pub fn run_selection_from_index_with_stages(
    index: &TranscriptEqIndex,
    num_targets: usize,
    num_eqcs: usize,
    stages: &SelectionStages,
) -> SelectionResult {
    let groups = signature_collapse(index);
    let n_no_eqc = (0..num_targets)
        .filter(|&t| index.degree(t) == 0)
        .count();
    info!(
        "  Stage 1 (signature collapse): {} transcripts → {} groups ({} with no EQ classes)",
        num_targets,
        groups.num_groups(),
        n_no_eqc
    );

    let required = if stages.peeling {
        let req = unique_ec_peeling_from_index(index, &groups, num_eqcs);
        let n_required = req.iter().filter(|&&r| r).count();
        info!(
            "  Stage 2 (unique EC peeling): {} / {} groups are structurally required",
            n_required,
            groups.num_groups()
        );
        req
    } else {
        info!("  Stage 2 (unique EC peeling): skipped");
        vec![false; groups.num_groups()]
    };

    let dominated = if stages.dominance {
        let dom = subset_dominance(index, &groups);
        let n_dominated = dom.iter().filter(|&&d| d).count();
        info!(
            "  Stage 3 (subset dominance): {} / {} groups are dominated (removable)",
            n_dominated,
            groups.num_groups()
        );
        dom
    } else {
        info!("  Stage 3 (subset dominance): skipped");
        vec![false; groups.num_groups()]
    };

    let mut keep_mask = vec![false; num_targets];
    for (gi, &is_dominated) in dominated.iter().enumerate() {
        if !is_dominated {
            for &member in groups.group_members(gi) {
                keep_mask[member as usize] = true;
            }
        }
    }

    let num_kept = keep_mask.iter().filter(|&&k| k).count();
    let num_removed = num_targets - num_kept;
    info!(
        "  Variable selection result: {} kept, {} removed ({:.1}% reduction)",
        num_kept,
        num_removed,
        100.0 * num_removed as f64 / num_targets as f64
    );

    SelectionResult {
        keep_mask,
        groups,
        required,
        dominated,
        num_kept,
        num_removed,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::eq_maps::{BasicEqLabel, BasicEqMap, OrientationProperty};

    /// Helper: build a PackedEqMap from a list of (target_ids, count) pairs.
    fn build_packed_map(eqcs: &[(&[u32], usize)]) -> PackedEqMap<BasicEqLabel> {
        let mut eqm = BasicEqMap::new(OrientationProperty::OrientationAgnostic);
        for &(targets, count) in eqcs {
            for _ in 0..count {
                eqm.add(BasicEqLabel::new(targets, None));
            }
        }
        PackedEqMap::from_eq_map(&eqm)
    }

    #[test]
    fn test_reverse_index() {
        // 4 targets, 3 EQ classes:
        //   EQC 0: {0, 1}   count=10
        //   EQC 1: {1, 2}   count=5
        //   EQC 2: {3}      count=8
        let packed = build_packed_map(&[(&[0, 1], 10), (&[1, 2], 5), (&[3], 8)]);
        let index = TranscriptEqIndex::from_packed_eq_map(&packed, 4);

        assert_eq!(index.num_targets(), 4);

        // Target 0: appears in EQC with {0,1}
        assert_eq!(index.degree(0), 1);
        // Target 1: appears in both EQCs with {0,1} and {1,2}
        assert_eq!(index.degree(1), 2);
        // Target 2: appears in EQC with {1,2}
        assert_eq!(index.degree(2), 1);
        // Target 3: appears in EQC with {3}
        assert_eq!(index.degree(3), 1);
    }

    #[test]
    fn test_signature_collapse() {
        // Targets 0 and 2 have identical signatures (same EQ classes).
        //   EQC 0: {0, 1, 2}  count=10
        //   EQC 1: {0, 2}     count=5
        //   EQC 2: {1}        count=3
        // Signatures: 0→{0,1}, 1→{0,2}, 2→{0,1}
        // So targets 0 and 2 collapse into one group.
        let packed = build_packed_map(&[(&[0, 1, 2], 10), (&[0, 2], 5), (&[1], 3)]);
        let index = TranscriptEqIndex::from_packed_eq_map(&packed, 3);
        let groups = signature_collapse(&index);

        // Should have 2 groups: {0,2} and {1}
        assert_eq!(groups.num_groups(), 2);

        // Find the group containing target 0.
        let g0 = (0..groups.num_groups())
            .find(|&i| groups.group_members(i).contains(&0))
            .unwrap();
        assert!(groups.group_members(g0).contains(&2));
        assert_eq!(groups.group_members(g0).len(), 2);
    }

    #[test]
    fn test_unique_ec_peeling() {
        // Target 0: unique EQC {0} (degree 1) → required
        // Target 1: shared EQCs only → not required initially
        // Target 2: unique EQC {2} (degree 1) → required
        //   EQC 0: {0}      count=5
        //   EQC 1: {0, 1}   count=10
        //   EQC 2: {1, 2}   count=8
        //   EQC 3: {2}      count=3
        let packed = build_packed_map(&[
            (&[0], 5),
            (&[0, 1], 10),
            (&[1, 2], 8),
            (&[2], 3),
        ]);
        let index = TranscriptEqIndex::from_packed_eq_map(&packed, 3);
        let groups = signature_collapse(&index);

        // Each target has a unique signature, so 3 groups.
        assert_eq!(groups.num_groups(), 3);

        let required = unique_ec_peeling(&packed, &groups, 3);

        // Target 0 (group 0): required via EQC {0}
        // Target 2 (group 2): required via EQC {2}
        // After removing 0 and 2, EQC {0,1} becomes degree-1 for target 1,
        // and EQC {1,2} becomes degree-1 for target 1.
        // So target 1 also becomes required via cascading.
        assert!(required.iter().all(|&r| r), "all groups should be required after cascading");
    }

    #[test]
    fn test_subset_dominance() {
        // Target 0: EQCs {0, 1, 2}
        // Target 1: EQCs {0, 1}      ← subset of target 0's signature
        // Target 2: EQCs {2, 3}      ← not a subset of anyone
        //   EQC 0: {0, 1}   count=10
        //   EQC 1: {0, 1}   count=5   (duplicate EQ label, but different instance)
        //   EQC 2: {0, 2}   count=8
        //   EQC 3: {2}      count=3
        // After building PackedEqMap, {0,1} with count=10 and {0,1} with count=5
        // will be merged into one EQC with count=15.
        //
        // Let's set it up more carefully:
        //   EQC A: {0, 1}  → targets 0 and 1 share this
        //   EQC B: {0}     → only target 0
        //   EQC C: {2}     → only target 2
        // Signatures: 0→{A,B}, 1→{A}, 2→{C}
        // Target 1's signature {A} ⊂ target 0's signature {A,B}
        // So target 1 is dominated.
        let packed = build_packed_map(&[(&[0, 1], 10), (&[0], 5), (&[2], 3)]);
        let index = TranscriptEqIndex::from_packed_eq_map(&packed, 3);
        let groups = signature_collapse(&index);

        // 3 groups (all unique signatures)
        assert_eq!(groups.num_groups(), 3);

        let dominated = subset_dominance(&index, &groups);

        // Find group for target 1 — it should be dominated.
        let g1 = (0..groups.num_groups())
            .find(|&i| groups.group_members(i).contains(&1))
            .unwrap();
        assert!(dominated[g1], "target 1 should be dominated by target 0");

        // Target 0 and target 2 should NOT be dominated.
        let g0 = (0..groups.num_groups())
            .find(|&i| groups.group_members(i).contains(&0))
            .unwrap();
        let g2 = (0..groups.num_groups())
            .find(|&i| groups.group_members(i).contains(&2))
            .unwrap();
        assert!(!dominated[g0], "target 0 should not be dominated");
        assert!(!dominated[g2], "target 2 should not be dominated");
    }

    #[test]
    fn test_is_sorted_subset() {
        assert!(is_sorted_subset(&[1, 3], &[1, 2, 3, 4]));
        assert!(is_sorted_subset(&[1, 2, 3], &[1, 2, 3]));
        assert!(!is_sorted_subset(&[1, 5], &[1, 2, 3, 4]));
        assert!(!is_sorted_subset(&[1, 2, 3], &[1, 2]));
        assert!(is_sorted_subset(&[], &[1, 2, 3]));
    }

    #[test]
    fn test_full_selection_pipeline() {
        // 5 targets:
        //   Target 0: unique evidence (EQC {0})
        //   Target 1: shared with 0 (EQC {0,1}), no unique evidence
        //   Target 2: unique evidence (EQC {2})
        //   Target 3: identical signature to target 4 (will collapse)
        //   Target 4: identical signature to target 3 (will collapse)
        //
        //   EQC 0: {0}      count=20
        //   EQC 1: {0, 1}   count=10
        //   EQC 2: {2}      count=15
        //   EQC 3: {3, 4}   count=8
        let packed = build_packed_map(&[
            (&[0], 20),
            (&[0, 1], 10),
            (&[2], 15),
            (&[3, 4], 8),
        ]);
        let result = run_selection(&packed, 5);

        // Targets 0 and 2: required (unique EQCs), kept
        assert!(result.keep_mask[0]);
        assert!(result.keep_mask[2]);

        // Target 1: signature {EQC1} ⊂ target 0's signature {EQC0, EQC1}
        // So target 1 is dominated and removed.
        assert!(!result.keep_mask[1]);

        // Targets 3 and 4: same signature group, have a shared EQC but
        // no unique evidence. After peeling removes 0 and 2, EQC {0,1}
        // has degree 1 → target 1 is required → but target 1 is dominated.
        // The group {3,4} has EQC {3,4} which has degree 1 (only that group),
        // so the group is required. Both members are kept.
        assert!(result.keep_mask[3]);
        assert!(result.keep_mask[4]);
    }

    #[test]
    fn test_merge_indices() {
        // Sample 0: 3 targets, 2 EQCs
        //   EQC 0: {0, 1}  EQC 1: {2}
        // Sample 1: 3 targets, 2 EQCs
        //   EQC 0: {0}     EQC 1: {1, 2}
        //
        // After merge with remapped IDs:
        // Target 0: sample0_EQC0(→0), sample1_EQC0(→2) = {0, 2}
        // Target 1: sample0_EQC0(→0), sample1_EQC1(→3) = {0, 3}
        // Target 2: sample0_EQC1(→1), sample1_EQC1(→3) = {1, 3}
        //
        // In sample 0 alone, targets 0 and 1 are indistinguishable.
        // After merge, they're distinguishable (different global signatures).

        let packed0 = build_packed_map(&[(&[0, 1], 10), (&[2], 5)]);
        let packed1 = build_packed_map(&[(&[0], 8), (&[1, 2], 3)]);

        let idx0 = TranscriptEqIndex::from_packed_eq_map(&packed0, 3);
        let idx1 = TranscriptEqIndex::from_packed_eq_map(&packed1, 3);

        let merged = merge_transcript_indices(&[idx0, idx1], &[2, 2], 3);

        // All three targets should now have distinct signatures.
        let sig0 = merged.signature(0);
        let sig1 = merged.signature(1);
        let sig2 = merged.signature(2);
        assert_ne!(sig0, sig1, "targets 0 and 1 should be distinguishable after merge");
        assert_ne!(sig1, sig2);
        assert_ne!(sig0, sig2);
    }
}
