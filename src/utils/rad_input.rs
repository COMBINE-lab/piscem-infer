//! Locating and validating a bulk RAD input, for both the legacy (v1) prelude
//! and the versioned (v2, `RAD_FILE` magic + tag roles) prelude.
//!
//! libradicl's bulk record parser reads records positionally, so anything in
//! the header that would make it misread the records is checked here first
//! and reported as an error naming the file.

use anyhow::{Context, anyhow, bail, ensure};
use libradicl::codec::{ChunkCodec, chunk_codec_from_tag_map};
use libradicl::header::RadPrelude;
use libradicl::rad_types::{RadIntId, RadType, TagMap, TagRole, TagSection, TagValue};
use std::io::BufRead;
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

#[cfg(test)]
mod tests {
    use super::*;
    use libradicl::chunk::Chunk;
    use libradicl::header::{RadHeader, SpecVersion};
    use libradicl::rad_types::{MappedFragmentOrientation, TagDesc, TagSectionLabel};
    use libradicl::record::{MappedRecord, PiscemBulkReadRecord, PiscemBulkRecordContext};

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
                out.extend_from_slice(&((body.len() + 8) as u32).to_le_bytes());
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
        let ctx = prelude.get_record_context_prefer_roles::<PiscemBulkRecordContext>()?;
        let mut recs = vec![];
        for _ in 0..prelude.hdr.num_chunks {
            for rec in Chunk::<PiscemBulkReadRecord>::from_bytes(r, &ctx).reads {
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
}
