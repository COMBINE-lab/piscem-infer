use arrow_array::{ArrayRef, RecordBatch};
use arrow_schema::{Field, Schema};
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, Encoding, ZstdLevel};
use parquet::file::properties::{EnabledStatistics, WriterProperties, WriterVersion};
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

/// Write `columns` (name, values) as one table to the Parquet file `path`:
/// non-nullable columns, plain-encoded and zstd-compressed, with statistics.
pub(crate) fn write_columns(path: &Path, columns: Vec<(String, ArrayRef)>) -> anyhow::Result<()> {
    let fields: Vec<Field> = columns
        .iter()
        .map(|(name, a)| Field::new(name, a.data_type().clone(), false))
        .collect();
    let schema = Arc::new(Schema::new(fields));
    let batch = RecordBatch::try_new(
        schema.clone(),
        columns.into_iter().map(|(_, a)| a).collect(),
    )?;

    let props = WriterProperties::builder()
        .set_writer_version(WriterVersion::PARQUET_2_0)
        .set_compression(Compression::ZSTD(ZstdLevel::default()))
        .set_dictionary_enabled(false)
        .set_encoding(Encoding::PLAIN)
        .set_statistics_enabled(EnabledStatistics::Page)
        .build();

    let file = File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(props))?;
    writer.write(&batch)?;
    writer.close()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow_array::{Array, Float64Array, UInt32Array};
    use arrow_schema::DataType;
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

    #[test]
    fn columns_read_back_unchanged() {
        let dir = std::env::temp_dir().join(format!("pi-pq-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("t.pq");
        let a: Vec<f64> = (0..100_000).map(|i| i as f64 * 0.25).collect();
        let b: Vec<f64> = (0..100_000).map(|i| (i % 7) as f64).collect();
        let c: Vec<u32> = (0..100_000).map(|i| i * 3).collect();
        write_columns(
            &path,
            vec![
                (
                    "bootstrap.0".into(),
                    Arc::new(Float64Array::from(a.clone())),
                ),
                (
                    "bootstrap.1".into(),
                    Arc::new(Float64Array::from(b.clone())),
                ),
                (
                    "fragment_length_dist".into(),
                    Arc::new(UInt32Array::from(c.clone())),
                ),
            ],
        )
        .unwrap();

        let reader = ParquetRecordBatchReaderBuilder::try_new(File::open(&path).unwrap())
            .unwrap()
            .build()
            .unwrap();
        let (mut got_a, mut got_b, mut got_c) = (vec![], vec![], vec![]);
        for batch in reader {
            let batch = batch.unwrap();
            let s = batch.schema();
            assert_eq!(s.field(0).name(), "bootstrap.0");
            assert_eq!(s.field(1).name(), "bootstrap.1");
            assert_eq!(s.field(2).name(), "fragment_length_dist");
            assert_eq!(s.field(0).data_type(), &DataType::Float64);
            assert_eq!(s.field(2).data_type(), &DataType::UInt32);
            assert!(s.fields().iter().all(|f| !f.is_nullable()));
            let col = |i: usize| batch.column(i).clone();
            got_a.extend(
                col(0)
                    .as_any()
                    .downcast_ref::<Float64Array>()
                    .unwrap()
                    .values()
                    .iter(),
            );
            got_b.extend(
                col(1)
                    .as_any()
                    .downcast_ref::<Float64Array>()
                    .unwrap()
                    .values()
                    .iter(),
            );
            got_c.extend(
                col(2)
                    .as_any()
                    .downcast_ref::<UInt32Array>()
                    .unwrap()
                    .values()
                    .iter(),
            );
        }
        assert_eq!((got_a, got_b, got_c), (a, b, c));
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
