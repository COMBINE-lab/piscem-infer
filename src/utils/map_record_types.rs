use libradicl::rad_types;
use scroll::Pread;
use std::io::Read;

use crate::utils::custom_rad_utils::*;

const MASK_LOWER_30_BITS: u32 = 0xC0000000;
const MASK_UPPER_2_BITS: u32 = 0x3FFFFFFF;

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum MappingType {
    Unmapped,
    SingleMapped,
    MappedFirstOrphan,
    MappedSecondOrphan,
    MappedPair,
}

impl MappingType {
    pub fn from_u8(t: u8) -> Self {
        match t {
            0 => MappingType::Unmapped,
            1 => MappingType::SingleMapped,
            2 => MappingType::MappedFirstOrphan,
            3 => MappingType::MappedSecondOrphan,
            4 => MappingType::MappedPair,
            _ => MappingType::Unmapped,
        }
    }

    #[inline]
    pub fn get_mask(&self) -> u32 {
        match &self {
            // if unmapped, we ignore flags all together.
            MappingType::Unmapped => 0b00,
            // if paired we care about both read and mate.
            MappingType::MappedPair => 0b11,
            // if orphan or single, we care only about read
            // mate flag should be ignored.
            _ => 0b10,
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub enum MappedFragmentOrientation {
    Reverse,
    Forward,
    ReverseReverse,
    ReverseForward,
    ForwardReverse,
    ForwardForward,
    Unknown,
}

impl MappedFragmentOrientation {
    pub fn from_u32_paired_status(n: u32, m: MappingType) -> Self {
        // if not paired, then we don't care about
        // the lowest order bit so shift it off
        if m == MappingType::SingleMapped {
            if (n & 0b10) == 2 {
                MappedFragmentOrientation::Forward
            } else {
                MappedFragmentOrientation::Reverse
            }
        } else {
            match n {
                0 => MappedFragmentOrientation::ReverseReverse,
                1 => MappedFragmentOrientation::ReverseForward,
                2 => MappedFragmentOrientation::ForwardReverse,
                3 => MappedFragmentOrientation::ForwardForward,
                _ => MappedFragmentOrientation::Unknown,
            }
        }
    }
}

impl From<MappedFragmentOrientation> for u32 {
    fn from(item: MappedFragmentOrientation) -> Self {
        match item {
            MappedFragmentOrientation::ForwardReverse => 0b011,
            MappedFragmentOrientation::ForwardForward => 0b101,
            MappedFragmentOrientation::ReverseReverse => 0b110,
            MappedFragmentOrientation::ReverseForward => 0b100,
            MappedFragmentOrientation::Forward => 0b1,
            MappedFragmentOrientation::Reverse => 0b10,
            MappedFragmentOrientation::Unknown => 0b0,
        }
    }
}

impl From<u32> for MappedFragmentOrientation {
    fn from(item: u32) -> Self {
        match item {
            0b011 => MappedFragmentOrientation::ForwardReverse,
            0b101 => MappedFragmentOrientation::ForwardForward,
            0b110 => MappedFragmentOrientation::ReverseReverse,
            0b100 => MappedFragmentOrientation::ReverseForward,
            0b1 => MappedFragmentOrientation::Forward,
            0b10 => MappedFragmentOrientation::Reverse,
            _ => MappedFragmentOrientation::Unknown,
        }
    }
}

#[derive(Debug)]
pub struct MetaReadRecord {
    pub frag_type: u8,
    pub dirs: Vec<MappedFragmentOrientation>,
    pub refs: Vec<u32>,
    pub positions: Vec<u32>,
    pub frag_lengths: Vec<u16>,
}

impl MetaReadRecord {
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.refs.is_empty()
    }

    pub fn from_bytes<T: Read>(reader: &mut T, frag_map_t: &rad_types::RadIntId) -> Self {
        let mut rbuf = [0u8; 255];

        reader.read_exact(&mut rbuf[0..4]).unwrap();
        let na = rbuf.pread::<u32>(0).unwrap();
        let fmt = read_into_u64(reader, frag_map_t);
        let f = MappingType::from_u8(fmt as u8);

        let mut rec = Self {
            frag_type: fmt as u8,
            dirs: Vec::with_capacity(na as usize),
            refs: Vec::with_capacity(na as usize),
            positions: Vec::with_capacity(na as usize),
            frag_lengths: Vec::with_capacity(na as usize),
        };

        //println!("number of records : {:?}",na);

        for _ in 0..(na as usize) {
            reader.read_exact(&mut rbuf[0..4]).unwrap();
            let v = rbuf.pread::<u32>(0).unwrap();

            let dir_int = (v & MASK_LOWER_30_BITS) >> 30;
            let dir = MappedFragmentOrientation::from_u32_paired_status(dir_int, f);
            rec.dirs.push(dir);
            rec.refs.push(v & MASK_UPPER_2_BITS);
            // position
            reader.read_exact(&mut rbuf[0..4]).unwrap();
            let pos = rbuf.pread::<u32>(0).unwrap();
            rec.positions.push(pos);
            // length
            reader.read_exact(&mut rbuf[0..2]).unwrap();
            let flen = rbuf.pread::<u16>(0).unwrap();
            rec.frag_lengths.push(flen);
        }

        rec
    }
}
