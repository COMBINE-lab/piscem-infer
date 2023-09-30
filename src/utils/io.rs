use crate::utils::parquet_utils;
use anyhow::Result;
use arrow2::{
    array::Array,
    chunk::Chunk,
    datatypes::{Field, Schema},
};
use libradicl::rad_types;
use path_tools::WithAdditionalExtension;
use std::fs::File;
use std::io::Write;
use std::path::Path;

pub(crate) fn write_results(
    output: &Path,
    hdr: &rad_types::RadHeader,
    e_counts: &[f64],
    eff_lengths: &[f64],
) -> anyhow::Result<()> {
    let mut ofile = File::create(output)?;

    ofile.write_all("target_name\teelen\tecount\n".to_string().as_bytes())?;

    for (i, name) in hdr.ref_names.iter().enumerate() {
        let l = format!("{}\t{:.3}\t{:.3}\n", name, eff_lengths[i], e_counts[i]);
        ofile.write_all(l.as_bytes())?;
    }

    Ok(())
}
pub(crate) fn write_fld_file(
    output_path: &Path,
    fields: Vec<Field>,
    chunk: Chunk<Box<dyn Array>>,
) -> Result<()> {
    let output_path = output_path
        .clone()
        .to_path_buf()
        .with_additional_extension(".fld.pq");
    let schema = Schema::from(fields);
    parquet_utils::write_chunk_to_file(output_path.to_str().unwrap(), schema, chunk)
}

pub(crate) fn write_infrep_file(
    output_path: &Path,
    fields: Vec<Field>,
    chunk: Chunk<Box<dyn Array>>,
) -> Result<()> {
    let output_path = output_path
        .clone()
        .to_path_buf()
        .with_additional_extension(".infreps.pq");
    let schema = Schema::from(fields);
    parquet_utils::write_chunk_to_file(output_path.to_str().unwrap(), schema, chunk)
}
