//! Bulk RAD records read in place.
//!
//! libradicl's `PiscemBulkReadRecord` decodes every record into four freshly
//! allocated vectors, which is most of the cost of reading a bulk RAD file. The
//! views here borrow the chunk's bytes instead and decode each field as it is
//! asked for, exactly as `PiscemBulkReadRecord::from_bytes_with_context` does.
//! The records must have been checked to fill their chunk (see
//! [`crate::utils::rad_input`]) before they are walked.

use libradicl::rad_types::{MappedFragmentOrientation, MappingType};

/// Bytes per alignment in the bulk layout: `compressed_ori_ref: u32`,
/// `pos: u32`, `frag_len: u16`.
pub(crate) const ALN_BYTES: usize = 10;

#[inline]
fn u32_at(b: &[u8], off: usize) -> u32 {
    u32::from_le_bytes(b[off..off + 4].try_into().expect("slice of length 4"))
}

/// One alignment of a bulk record.
#[derive(Clone, Copy, Debug)]
pub(crate) struct BulkAlignment {
    pub tid: u32,
    pub ori: MappedFragmentOrientation,
    pub pos: u32,
    pub frag_len: u16,
}

/// One bulk record, borrowing its alignments' bytes.
#[derive(Clone, Copy, Debug)]
pub(crate) struct BulkRecord<'a> {
    pub frag_type: u8,
    alns: &'a [u8],
}

impl<'a> BulkRecord<'a> {
    pub(crate) fn num_alignments(&self) -> usize {
        self.alns.len() / ALN_BYTES
    }

    pub(crate) fn mapping_type(&self) -> MappingType {
        MappingType::from_u8(self.frag_type)
    }

    /// The `i`-th alignment (in record order).
    #[inline]
    pub(crate) fn alignment(&self, i: usize) -> BulkAlignment {
        let a = &self.alns[i * ALN_BYTES..(i + 1) * ALN_BYTES];
        let v = u32_at(a, 0);
        BulkAlignment {
            tid: v & 0x3FFF_FFFF,
            ori: MappedFragmentOrientation::from_u32_paired_status(
                (v & 0xC000_0000) >> 30,
                self.mapping_type(),
            ),
            pos: u32_at(a, 4),
            frag_len: u16::from_le_bytes([a[8], a[9]]),
        }
    }

    pub(crate) fn alignments(&self) -> impl ExactSizeIterator<Item = BulkAlignment> + 'a {
        let rec = *self;
        (0..rec.num_alignments()).map(move |i| rec.alignment(i))
    }
}

/// The records of one chunk, in order.
pub(crate) struct BulkRecords<'a> {
    data: &'a [u8],
    remaining: u32,
    frag_type_bytes: usize,
}

impl<'a> BulkRecords<'a> {
    /// The `nrec` records in `payload` (a chunk without its header), whose
    /// mapping type takes `frag_type_bytes` bytes. `payload` must already have
    /// been checked to hold exactly those records.
    pub(crate) fn new(payload: &'a [u8], nrec: u32, frag_type_bytes: usize) -> Self {
        Self {
            data: payload,
            remaining: nrec,
            frag_type_bytes,
        }
    }
}

impl<'a> Iterator for BulkRecords<'a> {
    type Item = BulkRecord<'a>;

    #[inline]
    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            return None;
        }
        self.remaining -= 1;
        let na = u32_at(self.data, 0) as usize;
        // The mapping type is a little-endian integer of `frag_type_bytes`
        // bytes, of which libradicl keeps the low byte.
        let frag_type = self.data[4];
        let start = 4 + self.frag_type_bytes;
        let end = start + na * ALN_BYTES;
        let alns = &self.data[start..end];
        self.data = &self.data[end..];
        Some(BulkRecord { frag_type, alns })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.remaining as usize, Some(self.remaining as usize))
    }
}

impl ExactSizeIterator for BulkRecords<'_> {}
