use std::fs::File;
use std::io::{BufReader, Write as IoWrite};
use std::path::Path;

use anyhow::{Context, Result};
use arrow2::{
    array::{ListArray, PrimitiveArray, UInt32Array, UInt64Array},
    chunk::Chunk,
    datatypes::{DataType, Field, Schema},
    io::parquet::read::{self as pq_read, FileReader},
    offset::OffsetsBuffer,
};
use serde::{Deserialize, Serialize};

use crate::utils::eq_maps::{EqLabel, PackedEqMap};
use crate::utils::parquet_utils;

/// Metadata about a sample's equivalence class map, serialized as JSON.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct SampleMeta {
    pub sample_name: String,
    pub condition: String,
    pub num_targets: usize,
    pub num_eqcs: usize,
    pub total_weight: usize,
    pub ref_names: Vec<String>,
    pub ref_lengths: Vec<u32>,
    pub eff_lengths: Vec<f64>,
    pub eq_map_type: EqMapTypeTag,
    pub num_bins: u32,
    pub contains_ori: bool,
    /// Version of piscem-infer that produced this EQ map.
    /// Used to detect stale Phase A output when code changes affect EQ class structure.
    #[serde(default)]
    pub piscem_infer_version: Option<String>,
}

/// Tag to distinguish equivalence class map types during deserialization.
#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Eq)]
pub enum EqMapTypeTag {
    Basic,
    RangeFactorized,
}

/// Deserialized equivalence class map data (type-erased).
/// Holds the raw packed data without the generic EqLabel type parameter.
pub struct DeserializedEqMap {
    pub eq_labels: Vec<u32>,
    pub eq_label_starts: Vec<u32>,
    pub counts: Vec<usize>,
    pub contains_ori: bool,
}

/// Serialize a `PackedEqMap` and associated metadata to Parquet + JSON files.
///
/// Writes two files:
/// - `{sample_name}.eqc.pq` — Parquet with columns: eq_labels (List<UInt32>), counts (UInt64)
/// - `{sample_name}.eqmeta.json` — JSON metadata
pub fn serialize_eq_map<EqLabelT: EqLabel>(
    eq_map: &PackedEqMap<EqLabelT>,
    meta: &SampleMeta,
    output_dir: &Path,
) -> Result<()> {
    std::fs::create_dir_all(output_dir)
        .with_context(|| format!("Failed to create output directory: {}", output_dir.display()))?;

    // Write Parquet file
    let pq_path = output_dir.join(format!("{}.eqc.pq", meta.sample_name));
    write_eq_parquet(eq_map, pq_path.to_str().unwrap())?;

    // Write JSON metadata
    let meta_path = output_dir.join(format!("{}.eqmeta.json", meta.sample_name));
    let json = serde_json::to_string_pretty(meta)?;
    let mut file = File::create(&meta_path)
        .with_context(|| format!("Failed to create metadata file: {}", meta_path.display()))?;
    file.write_all(json.as_bytes())?;

    Ok(())
}

/// Deserialize a `PackedEqMap` and metadata from Parquet + JSON files.
pub fn deserialize_eq_map(
    eqc_path: &Path,
    meta_path: &Path,
) -> Result<(DeserializedEqMap, SampleMeta)> {
    // Read JSON metadata
    let meta_file = File::open(meta_path)
        .with_context(|| format!("Failed to open metadata file: {}", meta_path.display()))?;
    let meta: SampleMeta = serde_json::from_reader(BufReader::new(meta_file))
        .with_context(|| format!("Failed to parse metadata: {}", meta_path.display()))?;

    // Read Parquet file
    let eq_map = read_eq_parquet(eqc_path, meta.contains_ori)?;

    Ok((eq_map, meta))
}

/// Write PackedEqMap data as Parquet with a List<UInt32> column for labels
/// and a UInt64 column for counts.
fn write_eq_parquet<EqLabelT: EqLabel>(
    eq_map: &PackedEqMap<EqLabelT>,
    path: &str,
) -> Result<()> {
    let num_eqcs = eq_map.len();

    // Build the offsets for the list array from eq_label_starts
    let offsets: Vec<i32> = eq_map.eq_label_starts.iter().map(|&s| s as i32).collect();
    let offsets_buf = OffsetsBuffer::try_from(offsets)
        .map_err(|e| anyhow::anyhow!("Invalid offsets: {}", e))?;

    // Build the flat values array from eq_labels
    let values = UInt32Array::from_vec(eq_map.eq_labels.clone());

    // Build the List<UInt32> array
    let list_field = Field::new("item", DataType::UInt32, false);
    let list_data_type = DataType::List(Box::new(list_field));
    let labels_array = ListArray::new(list_data_type.clone(), offsets_buf, values.boxed(), None);

    // Build the counts array
    let counts: Vec<u64> = eq_map.counts.iter().map(|&c| c as u64).collect();
    let counts_array = UInt64Array::from_vec(counts);

    assert_eq!(labels_array.len(), num_eqcs);
    assert_eq!(counts_array.len(), num_eqcs);

    // Create schema and chunk
    let fields = vec![
        Field::new("eq_labels", list_data_type, false),
        Field::new("counts", DataType::UInt64, false),
    ];
    let schema = Schema::from(fields);
    let chunk = Chunk::new(vec![
        labels_array.boxed(),
        counts_array.boxed(),
    ]);

    parquet_utils::write_chunk_to_file(path, schema, chunk)
}

/// Read a PackedEqMap from a Parquet file.
fn read_eq_parquet(path: &Path, contains_ori: bool) -> Result<DeserializedEqMap> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open EQ class file: {}", path.display()))?;
    let reader = BufReader::new(file);

    let metadata = pq_read::read_metadata(&mut BufReader::new(File::open(path)?))?;
    let schema = pq_read::infer_schema(&metadata)?;

    let mut file_reader = FileReader::new(
        reader,
        metadata.row_groups,
        schema,
        None, // read all columns
        None, // no row limit
        None, // no page limit
    );

    // Read the single row group
    let chunk = file_reader
        .next()
        .ok_or_else(|| anyhow::anyhow!("Empty parquet file: {}", path.display()))?
        .with_context(|| format!("Failed to read row group from: {}", path.display()))?;

    let columns = chunk.columns();
    if columns.len() != 2 {
        anyhow::bail!(
            "Expected 2 columns (eq_labels, counts), got {}",
            columns.len()
        );
    }

    // Extract eq_labels (List<UInt32>)
    let labels_array = columns[0]
        .as_any()
        .downcast_ref::<ListArray<i32>>()
        .ok_or_else(|| anyhow::anyhow!("Column 0 is not a List<i32> array"))?;

    let values = labels_array
        .values()
        .as_any()
        .downcast_ref::<PrimitiveArray<u32>>()
        .ok_or_else(|| anyhow::anyhow!("List values are not UInt32"))?;

    let eq_labels: Vec<u32> = values.values_iter().copied().collect();
    let eq_label_starts: Vec<u32> = labels_array
        .offsets()
        .iter()
        .map(|&o| o as u32)
        .collect();

    // Extract counts (UInt64)
    let counts_array = columns[1]
        .as_any()
        .downcast_ref::<PrimitiveArray<u64>>()
        .ok_or_else(|| anyhow::anyhow!("Column 1 is not a UInt64 array"))?;

    let counts: Vec<usize> = counts_array.values_iter().map(|&c| c as usize).collect();

    Ok(DeserializedEqMap {
        eq_labels,
        eq_label_starts,
        counts,
        contains_ori,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::eq_maps::{BasicEqLabel, BasicEqMap, OrientationProperty, PackedEqMap};
    use tempfile::TempDir;

    fn build_test_eq_map() -> PackedEqMap<BasicEqLabel> {
        let mut eqm = BasicEqMap::new(OrientationProperty::OrientationAgnostic);
        for _ in 0..10 { eqm.add(BasicEqLabel::new(&[0, 1], None)); }
        for _ in 0..15 { eqm.add(BasicEqLabel::new(&[1, 2], None)); }
        for _ in 0..20 { eqm.add(BasicEqLabel::new(&[0], None)); }
        for _ in 0..5 { eqm.add(BasicEqLabel::new(&[2, 3, 4], None)); }
        for _ in 0..8 { eqm.add(BasicEqLabel::new(&[3], None)); }
        PackedEqMap::from_eq_map(&eqm)
    }

    #[test]
    fn test_roundtrip_basic() {
        let eq_map = build_test_eq_map();
        let meta = SampleMeta {
            sample_name: "test_sample".to_string(),
            condition: "control".to_string(),
            num_targets: 5,
            num_eqcs: eq_map.len(),
            total_weight: eq_map.total_weight(),
            ref_names: vec![
                "tx0".into(), "tx1".into(), "tx2".into(), "tx3".into(), "tx4".into(),
            ],
            ref_lengths: vec![100, 200, 150, 300, 250],
            eff_lengths: vec![90.5, 190.2, 140.8, 290.1, 240.3],
            eq_map_type: EqMapTypeTag::Basic,
            num_bins: 1,
            contains_ori: false,
            piscem_infer_version: Some(env!("CARGO_PKG_VERSION").to_string()),
        };

        let tmp_dir = TempDir::new().unwrap();
        let output_dir = tmp_dir.path();

        // Serialize
        serialize_eq_map(&eq_map, &meta, output_dir).unwrap();

        // Check files exist
        let pq_path = output_dir.join("test_sample.eqc.pq");
        let meta_path = output_dir.join("test_sample.eqmeta.json");
        assert!(pq_path.exists());
        assert!(meta_path.exists());

        // Deserialize
        let (deser_map, deser_meta) = deserialize_eq_map(&pq_path, &meta_path).unwrap();

        // Verify metadata
        assert_eq!(deser_meta.sample_name, "test_sample");
        assert_eq!(deser_meta.condition, "control");
        assert_eq!(deser_meta.num_targets, 5);
        assert_eq!(deser_meta.num_eqcs, eq_map.len());
        assert_eq!(deser_meta.total_weight, eq_map.total_weight());
        assert_eq!(deser_meta.ref_names.len(), 5);
        assert_eq!(deser_meta.ref_lengths, vec![100, 200, 150, 300, 250]);
        assert_eq!(deser_meta.eq_map_type, EqMapTypeTag::Basic);

        // Verify EQ map data
        assert_eq!(deser_map.eq_labels, eq_map.eq_labels);
        assert_eq!(deser_map.eq_label_starts, eq_map.eq_label_starts);
        assert_eq!(deser_map.counts, eq_map.counts);
        assert_eq!(deser_map.contains_ori, eq_map.contains_ori);
    }

    #[test]
    fn test_roundtrip_preserves_order() {
        // Verify that the label-to-count mapping is preserved exactly
        let eq_map = build_test_eq_map();
        let meta = SampleMeta {
            sample_name: "order_test".to_string(),
            condition: "treatment".to_string(),
            num_targets: 5,
            num_eqcs: eq_map.len(),
            total_weight: eq_map.total_weight(),
            ref_names: (0..5).map(|i| format!("tx{}", i)).collect(),
            ref_lengths: vec![100, 200, 150, 300, 250],
            eff_lengths: vec![90.0; 5],
            eq_map_type: EqMapTypeTag::Basic,
            num_bins: 1,
            contains_ori: false,
            piscem_infer_version: Some(env!("CARGO_PKG_VERSION").to_string()),
        };

        let tmp_dir = TempDir::new().unwrap();
        serialize_eq_map(&eq_map, &meta, tmp_dir.path()).unwrap();

        let pq_path = tmp_dir.path().join("order_test.eqc.pq");
        let meta_path = tmp_dir.path().join("order_test.eqmeta.json");
        let (deser, _) = deserialize_eq_map(&pq_path, &meta_path).unwrap();

        // For each EQC, verify the label slice matches
        for i in 0..eq_map.len() {
            let orig_start = eq_map.eq_label_starts[i] as usize;
            let orig_end = eq_map.eq_label_starts[i + 1] as usize;
            let deser_start = deser.eq_label_starts[i] as usize;
            let deser_end = deser.eq_label_starts[i + 1] as usize;

            assert_eq!(
                &eq_map.eq_labels[orig_start..orig_end],
                &deser.eq_labels[deser_start..deser_end],
                "EQC {} labels mismatch", i
            );
            assert_eq!(
                eq_map.counts[i], deser.counts[i],
                "EQC {} count mismatch", i
            );
        }
    }
}
