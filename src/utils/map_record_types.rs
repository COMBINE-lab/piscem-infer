use anyhow::bail;
use libradicl::rad_types::MappedFragmentOrientation;
use serde::Serialize;
use std::fmt;
use std::str::FromStr;
use tracing::warn;

// const MASK_LOWER_30_BITS: u32 = 0xC0000000;
// const MASK_UPPER_2_BITS: u32 = 0x3FFFFFFF;

#[derive(Debug, Copy, Clone, Eq, PartialEq, Serialize)]
pub enum LibraryType {
    StrandedForward,
    InwardStrandedForward,
    StrandedReverse,
    InwardStrandedReverse,
    Unstranded,
    InwardUnstranded,
    Any,
}

impl FromStr for LibraryType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match (s.to_uppercase()).as_str() {
            "SF" => Ok(Self::StrandedForward),
            "ISF" => Ok(Self::InwardStrandedForward),
            "SR" => Ok(Self::StrandedReverse),
            "ISR" => Ok(Self::InwardStrandedReverse),
            "U" => Ok(Self::Unstranded),
            "IU" => Ok(Self::InwardUnstranded),
            "ANY" => Ok(Self::Any),
            _ => bail!("{} is an unknown library type!", s),
        }
    }
}

impl fmt::Display for LibraryType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::StrandedForward => write!(f, "SF"),
            Self::InwardStrandedForward => write!(f, "ISF"),
            Self::StrandedReverse => write!(f, "SR"),
            Self::InwardStrandedReverse => write!(f, "ISR"),
            Self::Unstranded => write!(f, "U"),
            Self::InwardUnstranded => write!(f, "IU"),
            Self::Any => write!(f, "Any"),
        }
    }
}

impl LibraryType {
    pub fn is_compatible_with(&self, fo: MappedFragmentOrientation) -> bool {
        match (&self, fo) {
            (Self::StrandedForward, MappedFragmentOrientation::Forward) => true,
            (Self::InwardStrandedForward, MappedFragmentOrientation::Forward) => true,
            (Self::InwardStrandedForward, MappedFragmentOrientation::ForwardReverse) => true,

            (Self::StrandedReverse, MappedFragmentOrientation::Reverse) => true,
            (Self::InwardStrandedReverse, MappedFragmentOrientation::Reverse) => true,
            (Self::InwardStrandedReverse, MappedFragmentOrientation::ReverseForward) => true,

            (Self::Unstranded, MappedFragmentOrientation::Forward) => true,
            (Self::Unstranded, MappedFragmentOrientation::Reverse) => true,

            (Self::InwardUnstranded, MappedFragmentOrientation::Forward) => true,
            (Self::InwardUnstranded, MappedFragmentOrientation::Reverse) => true,
            (Self::InwardUnstranded, MappedFragmentOrientation::ForwardReverse) => true,
            (Self::InwardUnstranded, MappedFragmentOrientation::ReverseForward) => true,

            (Self::Any, _) => true,

            (_, _) => false,
        }
    }
}

#[derive(Debug, Default)]
pub struct OrientationCounts {
    pub forward: u64,
    pub reverse: u64,
    pub forward_reverse: u64,
    pub reverse_forward: u64,
    pub forward_forward: u64,
    pub reverse_reverse: u64,
    pub unknown: u64,
}

impl OrientationCounts {
    pub fn add(&mut self, o: MappedFragmentOrientation) {
        match o {
            MappedFragmentOrientation::Unknown => self.unknown += 1,
            MappedFragmentOrientation::Forward => self.forward += 1,
            MappedFragmentOrientation::Reverse => self.reverse += 1,
            MappedFragmentOrientation::ForwardReverse => self.forward_reverse += 1,
            MappedFragmentOrientation::ReverseForward => self.reverse_forward += 1,
            MappedFragmentOrientation::ForwardForward => self.forward_forward += 1,
            MappedFragmentOrientation::ReverseReverse => self.reverse_reverse += 1,
        }
    }
}

/// Detect library type from orientation counts using salmon-style heuristics.
/// Returns the detected `LibraryType` and the forward-strand ratio.
pub fn detect_library_type(counts: &OrientationCounts, paired_end: bool) -> (LibraryType, f64) {
    let (fw_like, rv_like) = if paired_end {
        // Combine orphan + paired signals, matching salmon's approach
        let fr = counts.forward as f64 + counts.forward_reverse as f64;
        let rf = counts.reverse as f64 + counts.reverse_forward as f64;
        (fr, rf)
    } else {
        (counts.forward as f64, counts.reverse as f64)
    };

    let total = fw_like + rv_like;
    if total == 0.0 {
        let lt = if paired_end {
            LibraryType::InwardUnstranded
        } else {
            LibraryType::Unstranded
        };
        return (lt, 0.5);
    }

    let ratio = fw_like / total;

    let lib_type = if paired_end {
        if ratio < 0.3 {
            LibraryType::InwardStrandedReverse
        } else if ratio > 0.7 {
            LibraryType::InwardStrandedForward
        } else {
            LibraryType::InwardUnstranded
        }
    } else if ratio < 0.3 {
        LibraryType::StrandedReverse
    } else if ratio > 0.7 {
        LibraryType::StrandedForward
    } else {
        LibraryType::Unstranded
    };

    (lib_type, ratio)
}

/// Emit warnings about strand bias based on detected library type.
pub fn check_strand_warnings(
    detected: LibraryType,
    ratio: f64,
    paired_end: bool,
) {
    let (fw_like, rv_like_frac) = (ratio, 1.0 - ratio);
    let _ = paired_end; // used only for context in the message

    match detected {
        LibraryType::StrandedForward | LibraryType::InwardStrandedForward => {
            if rv_like_frac >= 0.05 {
                warn!(
                    "Detected forward-stranded library ({}), but {:.1}% of sampled mappings \
                     have reverse-strand orientation. This may indicate library preparation issues.",
                    detected,
                    rv_like_frac * 100.0
                );
            }
        }
        LibraryType::StrandedReverse | LibraryType::InwardStrandedReverse => {
            if fw_like >= 0.05 {
                warn!(
                    "Detected reverse-stranded library ({}), but {:.1}% of sampled mappings \
                     have forward-strand orientation. This may indicate library preparation issues.",
                    detected,
                    fw_like * 100.0
                );
            }
        }
        LibraryType::Unstranded | LibraryType::InwardUnstranded => {
            let bias = (ratio - 0.5).abs();
            if bias > 0.05 {
                warn!(
                    "Detected unstranded library ({}), but forward/reverse ratio is {:.3} \
                     (expected ~0.5). Strand bias of {:.1}% detected.",
                    detected,
                    ratio,
                    bias * 100.0
                );
            }
        }
        _ => {}
    }
}

#[cfg(test)]
mod tests {
    use super::{LibraryType, OrientationCounts, detect_library_type};

    #[test]
    fn lib_types_parse() {
        let lss = ["SF", "ISF", "SR", "ISR", "U", "IU", "ANY"];
        let lts = [
            LibraryType::StrandedForward,
            LibraryType::InwardStrandedForward,
            LibraryType::StrandedReverse,
            LibraryType::InwardStrandedReverse,
            LibraryType::Unstranded,
            LibraryType::InwardUnstranded,
            LibraryType::Any,
        ];

        for (ls, lt) in lss.iter().zip(lts.iter()) {
            if let Ok(pt) = ls.parse::<LibraryType>() {
                assert_eq!(*lt, pt);
            }
        }
    }

    #[test]
    fn bad_lib_fails() {
        match "ABC".parse::<LibraryType>() {
            Err(x) => {
                assert_eq!("ABC is an unknown library type!", format!("{x}"));
            }
            _ => {
                panic!("ABC should not parse!");
            }
        }
    }

    #[test]
    fn auto_detect_paired_isf() {
        let mut counts = OrientationCounts::default();
        counts.forward_reverse = 8000;
        counts.reverse_forward = 200;
        let (lt, ratio) = detect_library_type(&counts, true);
        assert_eq!(lt, LibraryType::InwardStrandedForward);
        assert!(ratio > 0.7);
    }

    #[test]
    fn auto_detect_paired_isr() {
        let mut counts = OrientationCounts::default();
        counts.forward_reverse = 200;
        counts.reverse_forward = 8000;
        let (lt, ratio) = detect_library_type(&counts, true);
        assert_eq!(lt, LibraryType::InwardStrandedReverse);
        assert!(ratio < 0.3);
    }

    #[test]
    fn auto_detect_paired_iu() {
        let mut counts = OrientationCounts::default();
        counts.forward_reverse = 5000;
        counts.reverse_forward = 5000;
        let (lt, ratio) = detect_library_type(&counts, true);
        assert_eq!(lt, LibraryType::InwardUnstranded);
        assert!((0.3..=0.7).contains(&ratio));
    }

    #[test]
    fn auto_detect_single_sf() {
        let mut counts = OrientationCounts::default();
        counts.forward = 9000;
        counts.reverse = 1000;
        let (lt, _) = detect_library_type(&counts, false);
        assert_eq!(lt, LibraryType::StrandedForward);
    }

    #[test]
    fn auto_detect_single_sr() {
        let mut counts = OrientationCounts::default();
        counts.forward = 1000;
        counts.reverse = 9000;
        let (lt, _) = detect_library_type(&counts, false);
        assert_eq!(lt, LibraryType::StrandedReverse);
    }

    #[test]
    fn auto_detect_single_u() {
        let mut counts = OrientationCounts::default();
        counts.forward = 5000;
        counts.reverse = 5000;
        let (lt, ratio) = detect_library_type(&counts, false);
        assert_eq!(lt, LibraryType::Unstranded);
        assert!((0.3..=0.7).contains(&ratio));
    }

    #[test]
    fn auto_detect_zero_counts_paired() {
        let counts = OrientationCounts::default();
        let (lt, ratio) = detect_library_type(&counts, true);
        assert_eq!(lt, LibraryType::InwardUnstranded);
        assert!((ratio - 0.5).abs() < f64::EPSILON);
    }

    #[test]
    fn auto_detect_zero_counts_single() {
        let counts = OrientationCounts::default();
        let (lt, ratio) = detect_library_type(&counts, false);
        assert_eq!(lt, LibraryType::Unstranded);
        assert!((ratio - 0.5).abs() < f64::EPSILON);
    }

    #[test]
    fn auto_detect_boundary_030() {
        // Exactly at 0.3 boundary — should be unstranded (not < 0.3)
        let mut counts = OrientationCounts::default();
        counts.forward = 3000;
        counts.reverse = 7000;
        let (lt, _) = detect_library_type(&counts, false);
        assert_eq!(lt, LibraryType::Unstranded);
    }

    #[test]
    fn auto_detect_boundary_070() {
        // Exactly at 0.7 boundary — should be unstranded (not > 0.7)
        let mut counts = OrientationCounts::default();
        counts.forward = 7000;
        counts.reverse = 3000;
        let (lt, _) = detect_library_type(&counts, false);
        assert_eq!(lt, LibraryType::Unstranded);
    }

    #[test]
    fn lib_type_display() {
        assert_eq!(format!("{}", LibraryType::StrandedForward), "SF");
        assert_eq!(format!("{}", LibraryType::InwardStrandedForward), "ISF");
        assert_eq!(format!("{}", LibraryType::InwardUnstranded), "IU");
    }
}
