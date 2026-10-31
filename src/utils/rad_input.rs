//! Locating and validating a bulk RAD input, for both the legacy (v1) prelude
//! and the versioned (v2, `RAD_FILE` magic + tag roles) prelude.
//!
//! libradicl's bulk record parser reads records positionally and unwraps every
//! read, and the release profile uses `panic = "abort"`. So anything about the
//! input that would make it misread the records, or run off the end of the
//! file, is checked here first and reported as an error naming the file.

use anyhow::{Context, anyhow, bail, ensure};
use libradicl::codec::{ChunkCodec, chunk_codec_from_tag_map};
use libradicl::header::{RadHeader, RadPrelude};
use libradicl::rad_types::{RadIntId, RadType, TagMap, TagRole, TagSection, TagValue};
use libradicl::record::PiscemBulkRecordContext;
use libradicl::utils::{has_data_left, read_chunk_header};

use crate::utils::bulk_record::{ALN_BYTES, BulkRecords};
use std::io::{BufRead, Read};
use std::path::{Path, PathBuf};

/// Resolve the `-i` argument to the RAD file and its `map_info.json`.
///
/// `piscem map-bulk -o <stem>` writes `<stem>.rad` and `<stem>.map_info.json`
/// by appending to the stem's file name, so do the same here; replacing the
/// extension instead would turn the stem `sample.rep1` into `sample.rad`. The
/// path of the `.rad` file itself is accepted as well.
pub(crate) fn input_paths(input: &Path) -> (PathBuf, PathBuf) {
    let stem = if input.extension().is_some_and(|e| e == "rad") {
        input.with_extension("")
    } else {
        input.to_path_buf()
    };
    (
        append_to_file_name(&stem, ".rad"),
        append_to_file_name(&stem, ".map_info.json"),
    )
}

/// Append `suffix` to the last component of `stem` (as `PathBuf::add_extension`
/// does), so that e.g. `out/` becomes `out.rad` rather than `out/.rad`.
fn append_to_file_name(stem: &Path, suffix: &str) -> PathBuf {
    let mut name = stem.file_name().unwrap_or_default().to_os_string();
    name.push(suffix);
    stem.with_file_name(name)
}

/// The read-level tag of a bulk record: its name, or (in a versioned file) its
/// role.
const BULK_READ_TAG: (&str, TagRole) = ("frag_map_type", TagRole::MappingType);

/// The alignment-level tags of a bulk record, in on-disk order: the name (or,
/// in a versioned file, the role) and the width the record parser reads.
const BULK_ALN_TAGS: [(&str, RadIntId, TagRole); 3] = [
    ("compressed_ori_ref", RadIntId::U32, TagRole::Orientation),
    ("pos", RadIntId::U32, TagRole::MappingPosition),
    ("frag_len", RadIntId::U16, TagRole::FragmentLength),
];

/// Check that `prelude` describes a bulk RAD file, as written by `piscem
/// map-bulk`, in the record layout libradicl's `PiscemBulkReadRecord` parser
/// assumes. That parser ignores the declared tag types (except for the width of
/// the mapping type), so a file with other tags would be silently misread.
pub(crate) fn check_bulk_layout(prelude: &RadPrelude, path: &Path) -> anyhow::Result<()> {
    let (read_name, read_role) = BULK_READ_TAG;
    let read_ok = match prelude.read_tags.tags.as_slice() {
        [t] => {
            (t.name == read_name || t.role == read_role)
                && matches!(
                    t.typeid,
                    RadType::Int(RadIntId::U8 | RadIntId::U16 | RadIntId::U32 | RadIntId::U64)
                )
        }
        _ => false,
    };
    if !read_ok {
        bail!(
            "{} is not a bulk RAD file that piscem-infer can read: its read-level tags are [{}], \
             but bulk RAD written by `piscem map-bulk` has the single read-level tag {read_name}. \
             (Single-cell RAD, from `piscem map-sc` or `piscem map-sc-atac`, is processed with \
             alevin-fry, not piscem-infer.)",
            path.display(),
            describe_tags(&prelude.read_tags)
        );
    }

    let aln = &prelude.aln_tags.tags;
    let aln_ok = aln.len() == BULK_ALN_TAGS.len()
        && aln
            .iter()
            .zip(BULK_ALN_TAGS)
            .all(|(t, (name, width, role))| {
                (t.name == name || t.role == role) && t.typeid == RadType::Int(width)
            });
    if !aln_ok {
        let expected = BULK_ALN_TAGS
            .iter()
            .map(|(name, width, _)| format!("{name} ({})", type_name(RadType::Int(*width))))
            .collect::<Vec<_>>()
            .join(", ");
        bail!(
            "the alignment-level tags of {} are [{}], but piscem-infer reads bulk RAD as \
             written by `piscem map-bulk`, whose alignment-level tags are exactly [{expected}]",
            path.display(),
            describe_tags(&prelude.aln_tags)
        );
    }
    Ok(())
}

fn describe_tags(section: &TagSection) -> String {
    section
        .tags
        .iter()
        .map(|t| format!("{} ({})", t.name, type_name(t.typeid)))
        .collect::<Vec<_>>()
        .join(", ")
}

fn type_name(t: RadType) -> String {
    match t {
        RadType::Int(i) => format!("{i:?}").to_lowercase(),
        other => format!("{other:?}"),
    }
}

/// Read the file-level tag values, which directly follow the prelude, and
/// return the reference lengths from them. On success `reader` is left at the
/// first chunk.
pub(crate) fn read_file_tags<T: BufRead>(
    prelude: &RadPrelude,
    reader: &mut T,
    path: &Path,
) -> anyhow::Result<Vec<u32>> {
    let file_tag_map = prelude
        .file_tags
        .parse_tags_from_bytes(reader)
        .with_context(|| format!("could not read the file-level tags of {}", path.display()))?;

    // the records are parsed straight from the file, so the chunks must not be
    // compressed
    let codec = chunk_codec_from_tag_map(&file_tag_map)
        .with_context(|| format!("invalid file-level tags in {}", path.display()))?;
    ensure!(
        codec == ChunkCodec::None,
        "the chunks of {} are {}-compressed, but piscem-infer reads uncompressed bulk RAD, \
         as written by `piscem map-bulk`",
        path.display(),
        codec.as_str()
    );

    ref_lengths(&file_tag_map, &prelude.hdr.ref_names)
        .with_context(|| format!("invalid file-level tags in {}", path.display()))
}

/// The `ref_lengths` file-level tag, one length per reference. piscem writes
/// it as an array of u32, but the element width is self-described, so accept
/// any unsigned width as long as every length fits in a u32.
fn ref_lengths(file_tags: &TagMap, ref_names: &[String]) -> anyhow::Result<Vec<u32>> {
    const NAME: &str = "ref_lengths";
    let lens = match file_tags.get(NAME) {
        Some(TagValue::ArrayU8(v)) => lengths_to_u32(v, ref_names)?,
        Some(TagValue::ArrayU16(v)) => lengths_to_u32(v, ref_names)?,
        Some(TagValue::ArrayU32(v)) => v.clone(),
        Some(TagValue::ArrayU64(v)) => lengths_to_u32(v, ref_names)?,
        Some(TagValue::ArrayU128(v)) => lengths_to_u32(v, ref_names)?,
        Some(other) => bail!(
            "the file-level tag \"{NAME}\" has type {}, not an array of unsigned integers",
            type_name(other.rad_type())
        ),
        None => bail!("the required file-level tag \"{NAME}\" is missing"),
    };
    ensure!(
        lens.len() == ref_names.len(),
        "the file-level tag \"{NAME}\" holds {} lengths, but the header names {} references",
        lens.len(),
        ref_names.len()
    );
    Ok(lens)
}

fn lengths_to_u32<V>(lens: &[V], ref_names: &[String]) -> anyhow::Result<Vec<u32>>
where
    V: Copy + TryInto<u32> + std::fmt::Display,
{
    lens.iter()
        .enumerate()
        .map(|(i, &l)| {
            l.try_into().map_err(|_| {
                anyhow!(
                    "reference {} has length {l}, which does not fit in a u32",
                    ref_names.get(i).map_or("?", String::as_str)
                )
            })
        })
        .collect()
}

/// piscem writes the chunk count into the header only once mapping has
/// finished, so a header that records no chunks but is followed by chunk data
/// was left by a `piscem map-bulk` run that was killed or failed; quantifying
/// it would yield all-zero estimates. (A run that mapped nothing writes no
/// chunks, and is still accepted.) `reader` must be at the first chunk.
pub(crate) fn check_num_chunks<T: BufRead>(
    hdr: &RadHeader,
    reader: &mut T,
    path: &Path,
) -> anyhow::Result<()> {
    if hdr.num_chunks().is_none()
        && has_data_left(reader).with_context(|| format!("could not read {}", path.display()))?
    {
        bail!(
            "the header of {} records 0 chunks, but chunk data follows it. piscem writes the \
             chunk count only when mapping finishes, so the `piscem map-bulk` run that wrote \
             this file did not complete (it was killed or failed); please re-run it.",
            path.display()
        );
    }
    Ok(())
}

/// Bytes in a chunk header: `nbytes: u32` (which counts the header too) and
/// `nrec: u32`.
const CHUNK_HEADER_BYTES: usize = 8;

/// Reads the chunks of a bulk RAD file. Each chunk's bytes are read first, and
/// its records are checked to fill exactly those bytes before they are parsed,
/// so a truncated or corrupt file is reported as an error naming the file and
/// chunk rather than aborting inside libradicl's `Chunk::from_bytes`.
pub(crate) struct ChunkReader {
    path: PathBuf,
    num_chunks: usize,
    buf: Vec<u8>,
}

impl ChunkReader {
    pub(crate) fn new(path: &Path, num_chunks: usize) -> Self {
        Self {
            path: path.to_path_buf(),
            num_chunks,
            buf: Vec::new(),
        }
    }

    /// The number of chunks recorded in the header.
    pub(crate) fn num_chunks(&self) -> usize {
        self.num_chunks
    }

    /// Read the next chunk from `reader`; `idx` is its 0-based index in the
    /// file, used in error messages. (The tests check [`Self::read_records`]
    /// against the records libradicl decodes from the same chunk.)
    #[cfg(test)]
    pub(crate) fn read<T: BufRead>(
        &mut self,
        reader: &mut T,
        ctx: &PiscemBulkRecordContext,
        idx: usize,
    ) -> anyhow::Result<libradicl::chunk::Chunk<libradicl::record::PiscemBulkReadRecord>> {
        self.read_checked(reader, ctx, idx)?;
        Ok(libradicl::chunk::Chunk::from_bytes(
            &mut std::io::Cursor::new(&self.buf[..]),
            ctx,
        ))
    }

    /// As [`Self::read`], but walks the chunk's records in place rather than
    /// decoding each into an owned record.
    pub(crate) fn read_records<T: BufRead>(
        &mut self,
        reader: &mut T,
        ctx: &PiscemBulkRecordContext,
        idx: usize,
    ) -> anyhow::Result<BulkRecords<'_>> {
        let nrec = self.read_checked(reader, ctx, idx)?;
        Ok(BulkRecords::new(
            &self.buf[CHUNK_HEADER_BYTES..],
            nrec,
            ctx.frag_map_t.bytes_for_type(),
        ))
    }

    /// Read the next chunk into `self.buf` and check its records fill it;
    /// returns its record count.
    fn read_checked<T: BufRead>(
        &mut self,
        reader: &mut T,
        ctx: &PiscemBulkRecordContext,
        idx: usize,
    ) -> anyhow::Result<u32> {
        let (path, n) = (self.path.display(), self.num_chunks);
        let (nbytes, nrec) = read_chunk_header(reader).with_context(|| {
            format!("{path} is truncated: its header records {n} chunks, but it ends after {idx}")
        })?;
        ensure!(
            nbytes as usize >= CHUNK_HEADER_BYTES,
            "chunk {} of {n} in {path} is corrupt: its header records a size of {nbytes} bytes",
            idx + 1
        );

        self.buf.clear();
        self.buf.extend_from_slice(&nbytes.to_le_bytes());
        self.buf.extend_from_slice(&nrec.to_le_bytes());
        // read (rather than pre-allocate) the payload, so that a corrupt size
        // in a truncated file doesn't allocate up to 4 GiB before failing
        let payload_bytes = (nbytes as usize - CHUNK_HEADER_BYTES) as u64;
        let got = reader
            .take(payload_bytes)
            .read_to_end(&mut self.buf)
            .with_context(|| format!("could not read chunk {} of {n} from {path}", idx + 1))?;
        ensure!(
            got as u64 == payload_bytes,
            "{path} is truncated: it ends inside chunk {} of {n}, which should hold {nbytes} bytes",
            idx + 1
        );
        check_record_sizes(&self.buf[CHUNK_HEADER_BYTES..], nrec, ctx).with_context(|| {
            format!(
                "chunk {} of {n} in {path} is corrupt or not in the bulk RAD layout",
                idx + 1
            )
        })?;
        Ok(nrec)
    }
}

/// The records of a chunk handed over by libradicl's parallel reader (its
/// records, without the chunk header), checked to fill the chunk; `chunk` is
/// the chunk's 0-based index in the file, for the error message.
pub(crate) fn chunk_records<'a>(
    records: &'a [u8],
    nrec: u32,
    ctx: &PiscemBulkRecordContext,
    path: &Path,
    chunk: usize,
) -> anyhow::Result<BulkRecords<'a>> {
    check_record_sizes(records, nrec, ctx).with_context(|| {
        format!(
            "chunk {} in {} is corrupt or not in the bulk RAD layout",
            chunk + 1,
            path.display()
        )
    })?;
    Ok(BulkRecords::new(
        records,
        nrec,
        ctx.frag_map_t.bytes_for_type(),
    ))
}

/// Check that the `nrec` records of a chunk payload occupy exactly its bytes
/// under the bulk layout (`na: u32`, the mapping type, then `na` alignments of
/// [`BULK_ALN_TAGS`]), so that parsing them can't run past the chunk.
fn check_record_sizes(
    payload: &[u8],
    nrec: u32,
    ctx: &PiscemBulkRecordContext,
) -> anyhow::Result<()> {
    let rec_header = std::mem::size_of::<u32>() + ctx.frag_map_t.bytes_for_type();
    let aln_bytes: usize = BULK_ALN_TAGS
        .iter()
        .map(|(_, width, _)| width.bytes_for_type())
        .sum();
    debug_assert_eq!(aln_bytes, ALN_BYTES);
    let mut off = 0usize;
    for r in 0..nrec {
        ensure!(
            off + rec_header <= payload.len(),
            "record {r} of {nrec} runs past the end of the chunk"
        );
        let na = u32::from_le_bytes(payload[off..off + 4].try_into().expect("slice of length 4"))
            as usize;
        off += rec_header + na * aln_bytes;
        ensure!(
            off <= payload.len(),
            "record {r} of {nrec} ({na} alignments) runs past the end of the chunk"
        );
    }
    ensure!(
        off == payload.len(),
        "the chunk holds {} bytes after its {nrec} records",
        payload.len() - off
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use libradicl::chunk::Chunk;
    use libradicl::header::SpecVersion;
    use libradicl::rad_types::{MappedFragmentOrientation, TagDesc, TagSectionLabel};
    use libradicl::record::{MappedRecord, PiscemBulkReadRecord};

    #[test]
    fn input_paths_append_to_the_stem() {
        let p = |s: &str| PathBuf::from(s);
        let (rad, info) = input_paths(Path::new("out/SRR1.rep1"));
        assert_eq!(rad, p("out/SRR1.rep1.rad"));
        assert_eq!(info, p("out/SRR1.rep1.map_info.json"));
        // the .rad file itself resolves the same way (and never to e.g. a
        // sibling out/SRR1.rad)
        assert_eq!(input_paths(Path::new("out/SRR1.rep1.rad")), (rad, info));
        assert_eq!(
            input_paths(Path::new("sample")),
            (p("sample.rad"), p("sample.map_info.json"))
        );
        assert_eq!(
            input_paths(Path::new("sample.rad")),
            (p("sample.rad"), p("sample.map_info.json"))
        );
        assert_eq!(
            input_paths(Path::new("bulk_0.23.0")),
            (p("bulk_0.23.0.rad"), p("bulk_0.23.0.map_info.json"))
        );
        // dots in a parent directory don't matter
        assert_eq!(
            input_paths(Path::new("run.1/sample")),
            (p("run.1/sample.rad"), p("run.1/sample.map_info.json"))
        );
        // like `piscem map-bulk -o out/`, which writes out.rad
        assert_eq!(
            input_paths(Path::new("out/")),
            (p("out.rad"), p("out.map_info.json"))
        );
    }

    fn bulk_aln_tags(with_roles: bool) -> TagSection {
        let mut aln_tags = TagSection::new_with_label(TagSectionLabel::AlignmentTags);
        for (name, width, role) in BULK_ALN_TAGS {
            let d = TagDesc::new(name, RadType::Int(width));
            aln_tags.add_tag_desc(if with_roles { d.with_role(role) } else { d });
        }
        aln_tags
    }

    fn bulk_read_tags(with_roles: bool) -> TagSection {
        let mut read_tags = TagSection::new_with_label(TagSectionLabel::ReadTags);
        let d = TagDesc::new(BULK_READ_TAG.0, RadType::Int(RadIntId::U8));
        read_tags.add_tag_desc(if with_roles {
            d.with_role(BULK_READ_TAG.1)
        } else {
            d
        });
        read_tags
    }

    fn records() -> Vec<PiscemBulkReadRecord> {
        vec![
            PiscemBulkReadRecord {
                frag_type: 4, // paired
                dirs: vec![MappedFragmentOrientation::ForwardReverse],
                refs: vec![1],
                positions: vec![10],
                frag_lengths: vec![300],
            },
            PiscemBulkReadRecord {
                frag_type: 4,
                dirs: vec![
                    MappedFragmentOrientation::ReverseForward,
                    MappedFragmentOrientation::ForwardReverse,
                ],
                refs: vec![0, 2],
                positions: vec![5, 7],
                frag_lengths: vec![250, 251],
            },
            PiscemBulkReadRecord {
                frag_type: 1, // single
                dirs: vec![MappedFragmentOrientation::Forward],
                refs: vec![2],
                positions: vec![42],
                frag_lengths: vec![0],
            },
        ]
    }

    /// A bulk RAD file in the layout `piscem map-bulk` writes (C++ piscem and
    /// piscem-rs 0.10 for `Legacy`, piscem-rs 0.11 for v2), with the records
    /// split into chunks of `per_chunk`. Roles are only written to v2 files.
    struct BulkRad {
        version: SpecVersion,
        ref_lengths: TagValue,
        extra_file_tags: Vec<(&'static str, TagValue)>,
        read_tags: TagSection,
        aln_tags: TagSection,
        num_chunks: Option<u64>,
        per_chunk: usize,
    }

    impl BulkRad {
        fn new(version: SpecVersion) -> Self {
            let roles = version.is_versioned();
            Self {
                version,
                ref_lengths: TagValue::ArrayU32(vec![1000, 2000, 3000]),
                extra_file_tags: vec![],
                read_tags: bulk_read_tags(roles),
                aln_tags: bulk_aln_tags(roles),
                num_chunks: None,
                per_chunk: 2,
            }
        }

        /// The file's bytes, and the offset at which its first chunk starts.
        fn bytes(self) -> (Vec<u8>, usize) {
            let recs = records();
            let chunks: Vec<_> = recs.chunks(self.per_chunk).collect();
            let mut hdr = RadHeader::from_ref_names(["t0", "t1", "t2"]);
            hdr.version = self.version;
            hdr.is_paired = 1;
            hdr.num_chunks = self.num_chunks.unwrap_or(chunks.len() as u64);
            let mut file_tag_values = vec![
                ("known_rad_type", TagValue::String("bulk_with_pos".into())),
                ("ref_lengths", self.ref_lengths),
            ];
            file_tag_values.extend(self.extra_file_tags);
            let (file_tags, file_vals) =
                TagSection::from_tag_values(TagSectionLabel::FileTags, &file_tag_values);
            let prelude = RadPrelude::from_header_and_tag_sections(
                hdr,
                file_tags,
                self.read_tags,
                self.aln_tags,
            );
            let mut out = Vec::new();
            prelude.write(&mut out).unwrap();
            file_vals.write_values(&mut out).unwrap();
            let first_chunk = out.len();

            let ctx = PiscemBulkRecordContext {
                frag_map_t: RadIntId::U8,
            };
            for c in chunks {
                let mut body = Vec::new();
                for r in c {
                    r.write(&mut body, &ctx).unwrap();
                }
                out.extend_from_slice(&((body.len() + CHUNK_HEADER_BYTES) as u32).to_le_bytes());
                out.extend_from_slice(&(c.len() as u32).to_le_bytes());
                out.extend_from_slice(&body);
            }
            (out, first_chunk)
        }
    }

    /// What process_bulk_dispatch reads from a file: the reference lengths
    /// and, for each record, its references, positions and fragment lengths.
    type Contents = (Vec<u32>, Vec<(Vec<u32>, Vec<u32>, Vec<u16>)>);

    /// Read `bytes` the way process_bulk_dispatch does.
    fn read_bulk(mut bytes: &[u8]) -> anyhow::Result<Contents> {
        let path = Path::new("x.rad");
        let r = &mut bytes;
        let prelude = RadPrelude::from_bytes(r)?;
        check_bulk_layout(&prelude, path)?;
        let lens = read_file_tags(&prelude, r, path)?;
        check_num_chunks(&prelude.hdr, r, path)?;
        let ctx = prelude.get_record_context_prefer_roles::<PiscemBulkRecordContext>()?;
        let mut chunks = ChunkReader::new(path, prelude.hdr.num_chunks as usize);
        let mut recs = vec![];
        for i in 0..chunks.num_chunks() {
            for rec in chunks.read(r, &ctx, i)?.reads {
                recs.push((rec.refs, rec.positions, rec.frag_lengths));
            }
        }
        Ok((lens, recs))
    }

    fn expected() -> Contents {
        let recs = records()
            .into_iter()
            .map(|r| (r.refs, r.positions, r.frag_lengths))
            .collect();
        (vec![1000, 2000, 3000], recs)
    }

    fn err_of(bytes: &[u8]) -> String {
        format!("{:#}", read_bulk(bytes).unwrap_err())
    }

    #[test]
    fn reads_v1_and_v2_bulk_rad() {
        for v in [SpecVersion::Legacy, SpecVersion::current()] {
            let (bytes, _) = BulkRad::new(v).bytes();
            assert_eq!(read_bulk(&bytes).unwrap(), expected(), "{v:?}");
        }
    }

    #[test]
    fn libradicl_reads_the_test_files_the_same_way() {
        // the buffered reader must parse exactly what Chunk::from_bytes does
        let (bytes, first) = BulkRad::new(SpecVersion::current()).bytes();
        let ctx = PiscemBulkRecordContext {
            frag_map_t: RadIntId::U8,
        };
        let mut direct = &bytes[first..];
        let mut buffered = &bytes[first..];
        let mut chunks = ChunkReader::new(Path::new("x.rad"), 2);
        for i in 0..2 {
            // (records only derive Debug)
            let want = Chunk::<PiscemBulkReadRecord>::from_bytes(&mut direct, &ctx);
            let got = chunks.read(&mut buffered, &ctx, i).unwrap();
            assert_eq!(format!("{got:?}"), format!("{want:?}"));
        }
        assert!(direct.is_empty() && buffered.is_empty());
    }

    #[test]
    fn records_read_in_place_match_libradicl() {
        // Hand-encoded records covering every orientation code under each
        // mapping type, at every width of the mapping type; the in-place view
        // must decode exactly what libradicl's own record parser does.
        for frag_map_t in [RadIntId::U8, RadIntId::U16, RadIntId::U32, RadIntId::U64] {
            let ctx = PiscemBulkRecordContext { frag_map_t };
            let mut payload = Vec::new();
            let mut nrec = 0u32;
            for frag_type in 0u64..=5 {
                for na in 0u32..=4 {
                    payload.extend_from_slice(&na.to_le_bytes());
                    payload
                        .extend_from_slice(&frag_type.to_le_bytes()[..frag_map_t.bytes_for_type()]);
                    for i in 0..na {
                        let tid = 1_000_000 * frag_type as u32 + 7 * i;
                        payload.extend_from_slice(&((i << 30) | tid).to_le_bytes());
                        payload.extend_from_slice(&(10 * i + 3).to_le_bytes());
                        payload.extend_from_slice(&(250 + i as u16).to_le_bytes());
                    }
                    nrec += 1;
                }
            }
            let mut framed = Vec::new();
            framed.extend_from_slice(&((payload.len() + 8) as u32).to_le_bytes());
            framed.extend_from_slice(&nrec.to_le_bytes());
            framed.extend_from_slice(&payload);
            let want = Chunk::<PiscemBulkReadRecord>::from_bytes(&mut &framed[..], &ctx);

            let got = chunk_records(&payload, nrec, &ctx, Path::new("x.rad"), 0).unwrap();
            assert_eq!(got.len(), want.reads.len());
            for (g, w) in got.zip(&want.reads) {
                assert_eq!(g.frag_type, w.frag_type);
                assert_eq!(g.num_alignments(), w.refs.len());
                for (i, a) in g.alignments().enumerate() {
                    assert_eq!(a.tid, w.refs[i]);
                    assert_eq!(u32::from(a.ori), u32::from(w.dirs[i]));
                    assert_eq!(a.pos, w.positions[i]);
                    assert_eq!(a.frag_len, w.frag_lengths[i]);
                }
            }

            // a record that runs past its chunk is an error, not a panic
            let err = chunk_records(
                &payload[..payload.len() - 1],
                nrec,
                &ctx,
                Path::new("x.rad"),
                6,
            )
            .err()
            .expect("short chunk");
            assert!(format!("{err:#}").contains("chunk 7 in x.rad"), "{err:#}");
        }
    }

    #[test]
    fn accepts_any_unsigned_ref_lengths_width() {
        let lens = [
            TagValue::ArrayU16(vec![1000, 2000, 3000]),
            TagValue::ArrayU64(vec![1000, 2000, 3000]),
            TagValue::ArrayU128(vec![1000, 2000, 3000]),
        ];
        for l in lens {
            let mut f = BulkRad::new(SpecVersion::current());
            f.ref_lengths = l.clone();
            assert_eq!(read_bulk(&f.bytes().0).unwrap(), expected(), "{l:?}");
        }
        let mut f = BulkRad::new(SpecVersion::Legacy);
        f.ref_lengths = TagValue::ArrayU8(vec![10, 20, 30]);
        assert_eq!(read_bulk(&f.bytes().0).unwrap().0, vec![10, 20, 30]);
    }

    #[test]
    fn rejects_bad_ref_lengths() {
        let cases = [
            (
                TagValue::ArrayU64(vec![1000, 1 << 32, 3000]),
                "reference t1 has length 4294967296",
            ),
            (TagValue::ArrayU32(vec![1000, 2000]), "holds 2 lengths"),
            (TagValue::ArrayI32(vec![1000, 2000, 3000]), "has type Array"),
            (TagValue::U32(1000), "has type u32"),
        ];
        for (lens, msg) in cases {
            let mut f = BulkRad::new(SpecVersion::current());
            f.ref_lengths = lens;
            let e = err_of(&f.bytes().0);
            assert!(e.contains("x.rad") && e.contains(msg), "{e}");
        }
    }

    #[test]
    fn rejects_compressed_chunks() {
        let mut f = BulkRad::new(SpecVersion::current());
        f.extra_file_tags = vec![("chunk_codec", TagValue::U8(ChunkCodec::Lz4.as_u8()))];
        let e = err_of(&f.bytes().0);
        assert!(e.contains("lz4-compressed"), "{e}");
    }

    #[test]
    fn layout_accepts_v2_roles_under_other_names() {
        let mut f = BulkRad::new(SpecVersion::current());
        f.read_tags.tags[0].name = "frag_map_typX".into();
        for t in &mut f.aln_tags.tags {
            t.name.push('X');
        }
        assert_eq!(read_bulk(&f.bytes().0).unwrap(), expected());
    }

    #[test]
    fn layout_rejects_renamed_tags_without_roles() {
        let mut f = BulkRad::new(SpecVersion::Legacy);
        f.aln_tags.tags[2].name = "frag_lex".into();
        let e = err_of(&f.bytes().0);
        assert!(e.contains("frag_lex (u16)"), "{e}");

        let mut f = BulkRad::new(SpecVersion::Legacy);
        f.read_tags.tags[0].name = "frag_map_typX".into();
        let e = err_of(&f.bytes().0);
        assert!(e.contains("frag_map_typX (u8)"), "{e}");
    }

    #[test]
    fn layout_rejects_other_widths_order_and_extra_tags() {
        // frag_len declared as u32 (the parser would still read a u16)
        let mut f = BulkRad::new(SpecVersion::current());
        f.aln_tags.tags[2].typeid = RadType::Int(RadIntId::U32);
        assert!(err_of(&f.bytes().0).contains("frag_len (u32)"));

        // pos and frag_len swapped
        let mut f = BulkRad::new(SpecVersion::Legacy);
        f.aln_tags.tags.swap(1, 2);
        assert!(err_of(&f.bytes().0).contains("[compressed_ori_ref (u32), frag_len (u16)"));

        // an extra (e.g. score) alignment tag
        let mut f = BulkRad::new(SpecVersion::Legacy);
        f.aln_tags
            .add_tag_desc(TagDesc::new("alignment_score", RadType::Int(RadIntId::I32)));
        assert!(err_of(&f.bytes().0).contains("alignment_score (i32)"));

        // a second read-level tag
        let mut f = BulkRad::new(SpecVersion::Legacy);
        f.read_tags
            .add_tag_desc(TagDesc::new("extra", RadType::Int(RadIntId::U8)));
        assert!(err_of(&f.bytes().0).contains("[frag_map_type (u8), extra (u8)]"));
    }

    #[test]
    fn single_cell_rad_is_named_as_such() {
        // the tags `piscem map-sc` writes
        let mut prelude = RadPrelude::from_header_and_tag_sections(
            RadHeader::from_ref_names(["t0"]),
            TagSection::new_with_label(TagSectionLabel::FileTags),
            TagSection::new_with_label(TagSectionLabel::ReadTags),
            TagSection::new_with_label(TagSectionLabel::AlignmentTags),
        );
        prelude.hdr.version = SpecVersion::current();
        for name in ["b", "u"] {
            prelude
                .read_tags
                .add_tag_desc(TagDesc::new(name, RadType::Int(RadIntId::U32)));
        }
        prelude.aln_tags.add_tag_desc(
            TagDesc::new("compressed_ori_refid", RadType::Int(RadIntId::U32))
                .with_role(TagRole::Orientation),
        );
        let e = check_bulk_layout(&prelude, Path::new("map/map.rad"))
            .unwrap_err()
            .to_string();
        assert!(
            e.contains("map/map.rad is not a bulk RAD file")
                && e.contains("[b (u32), u (u32)]")
                && e.contains("alevin-fry"),
            "{e}"
        );
    }

    #[test]
    fn truncation_anywhere_is_an_error() {
        for v in [SpecVersion::Legacy, SpecVersion::current()] {
            let (bytes, first) = BulkRad::new(v).bytes();
            for cut in 0..bytes.len() {
                let e = err_of(&bytes[..cut]);
                if cut >= first {
                    assert!(e.contains("x.rad is truncated"), "cut at {cut}: {e}");
                }
            }
        }
        let (bytes, first) = BulkRad::new(SpecVersion::current()).bytes();
        // at a chunk boundary
        let second =
            first + u32::from_le_bytes(bytes[first..first + 4].try_into().unwrap()) as usize;
        let e = err_of(&bytes[..second]);
        assert!(e.contains("records 2 chunks, but it ends after 1"), "{e}");
        // inside a chunk
        let e = err_of(&bytes[..bytes.len() - 1]);
        assert!(e.contains("it ends inside chunk 2 of 2"), "{e}");
    }

    #[test]
    fn corrupt_chunk_sizes_are_errors() {
        let (bytes, first) = BulkRad::new(SpecVersion::current()).bytes();
        // an nbytes smaller than the header itself
        let mut b = bytes.clone();
        b[first..first + 4].copy_from_slice(&4u32.to_le_bytes());
        assert!(err_of(&b).contains("records a size of 4 bytes"));
        // an nbytes far beyond the end of the file
        let mut b = bytes.clone();
        b[first..first + 4].copy_from_slice(&u32::MAX.to_le_bytes());
        let e = err_of(&b);
        assert!(e.contains("ends inside chunk 1 of 2"), "{e}");
        // a record claiming more alignments than fit in its chunk
        let mut b = bytes.clone();
        b[first + 8..first + 12].copy_from_slice(&1000u32.to_le_bytes());
        assert!(err_of(&b).contains("runs past the end of the chunk"));
        // a chunk claiming fewer records than it holds, leaving bytes over
        let mut b = bytes;
        b[first + 4..first + 8].copy_from_slice(&1u32.to_le_bytes());
        assert!(err_of(&b).contains("bytes after its 1 records"));
    }

    #[test]
    fn zero_chunk_header_is_an_error_only_if_chunks_follow() {
        for v in [SpecVersion::Legacy, SpecVersion::current()] {
            let mut f = BulkRad::new(v);
            f.num_chunks = Some(0);
            let (bytes, first) = f.bytes();
            // an interrupted `map-bulk` run: chunks written, count never patched
            let e = err_of(&bytes);
            assert!(
                e.contains("records 0 chunks, but chunk data follows"),
                "{e}"
            );
            // a run that mapped nothing: no chunks at all
            let (lens, recs) = read_bulk(&bytes[..first]).unwrap();
            assert_eq!(lens, vec![1000, 2000, 3000]);
            assert!(recs.is_empty());
        }
    }
}
