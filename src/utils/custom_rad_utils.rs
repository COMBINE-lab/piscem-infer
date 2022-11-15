use std::io::Read;
use scroll::Pread;

use crate::map_record_types::MetaReadRecord;

#[derive(Debug)]
pub struct MetaChunk {
    pub nbytes: u32,
    pub nrec: u32,
    pub reads: Vec<MetaReadRecord>,
}


impl MetaChunk {
    pub fn read_header<T: Read>(reader: &mut T) -> (u32, u32) {
        let mut buf = [0u8; 8];

        reader.read_exact(&mut buf).unwrap();
        let nbytes = buf.pread::<u32>(0).unwrap();
        let nrec = buf.pread::<u32>(4).unwrap();
        (nbytes, nrec)
    }

    pub fn from_bytes<T: Read>(reader: &mut T, fmt: &rad_types::RadIntId) -> Self {
        let mut buf = [0u8; 8];

        reader.read_exact(&mut buf).unwrap();
        let nbytes = buf.pread::<u32>(0).unwrap();
        let nrec = buf.pread::<u32>(4).unwrap();
        let mut c = Self {
            nbytes,
            nrec,
            reads: Vec::with_capacity(nrec as usize),
        };

        for _ in 0..(nrec as usize) {
            c.reads.push(MetaReadRecord::from_bytes(reader, fmt));
        }

        c
    }

    /// peeks to the first record in the buffer `buf`, and returns
    /// the barcode and umi associated with this record.  It is assumed
    /// that there is at least one record present in the buffer.
    pub fn peek_record(
        buf: &[u8],
        bct: &rad_types::RadIntId,
        umit: &rad_types::RadIntId,
    ) -> (u64, u64) {
        let na_size = mem::size_of::<u32>();
        let bc_size = bct.bytes_for_type();

        let _na = buf.pread::<u32>(0).unwrap();

        let bc = match bct {
            rad_types::RadIntId::U8 => buf.pread::<u8>(na_size).unwrap() as u64,
            rad_types::RadIntId::U16 => buf.pread::<u16>(na_size).unwrap() as u64,
            rad_types::RadIntId::U32 => buf.pread::<u32>(na_size).unwrap() as u64,
            rad_types::RadIntId::U64 => buf.pread::<u64>(na_size).unwrap(),
        };
        let umi = match umit {
            rad_types::RadIntId::U8 => buf.pread::<u8>(na_size + bc_size).unwrap() as u64,
            rad_types::RadIntId::U16 => buf.pread::<u16>(na_size + bc_size).unwrap() as u64,
            rad_types::RadIntId::U32 => buf.pread::<u32>(na_size + bc_size).unwrap() as u64,
            rad_types::RadIntId::U64 => buf.pread::<u64>(na_size + bc_size).unwrap(),
        };
        (bc, umi)
    }
}

pub fn read_into_u64<T: Read>(reader: &mut T, rt: &rad_types::RadIntId) -> u64 {
    let mut rbuf = [0u8; 8];

    let v: u64 = match rt {
        rad_types::RadIntId::U8 => {
            reader.read_exact(&mut rbuf[0..1]).unwrap();
            rbuf.pread::<u8>(0).unwrap() as u64
        }
        rad_types::RadIntId::U16 => {
            reader.read_exact(&mut rbuf[0..2]).unwrap();
            rbuf.pread::<u16>(0).unwrap() as u64
        }
        rad_types::RadIntId::U32 => {
            reader.read_exact(&mut rbuf[0..4]).unwrap();
            rbuf.pread::<u32>(0).unwrap() as u64
        }
        rad_types::RadIntId::U64 => {
            reader.read_exact(&mut rbuf[0..8]).unwrap();
            rbuf.pread::<u64>(0).unwrap()
        }
    };
    v
}
