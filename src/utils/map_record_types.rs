use anyhow::bail;
use libradicl::rad_types::MappedFragmentOrientation;
use serde::Serialize;
use std::str::FromStr;

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

#[cfg(test)]
mod tests {
    use super::LibraryType;

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
}
