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

use tracing::warn;

pub(crate) fn write_results(
    output: &Path,
    hdr: &rad_types::RadHeader,
    e_counts: &[f64],
    lengths: &[u32],
    eff_lengths: &[f64],
) -> anyhow::Result<()> {
    let mut ofile = File::create(output)?;

    ofile.write_all(
        "target_name\tlen\teelen\ttpm\tecount\n"
            .to_string()
            .as_bytes(),
    )?;

    const ONE_MILLION: f64 = 1_000_000.0;
    let denom: f64 = e_counts
        .iter()
        .zip(eff_lengths.iter())
        .map(|(c, l)| if *l > 0.0 { c / l } else { 0.0_f64 })
        .sum::<f64>();
    let inv_denom: f64 = if denom > 0.0 {
        ONE_MILLION / denom
    } else {
        warn!("The sum of ecount / eeln for all transcripts was 0. It seems likely that no fragments were quantified. Please check the input sample!");
        0.0
    };
    let tpms: Vec<f64> = e_counts
        .iter()
        .zip(eff_lengths.iter())
        .map(|(c, l)| {
            if *l > 0.0 {
                inv_denom * (c / l)
            } else {
                0.0_f64
            }
        })
        .collect();

    for (i, name) in hdr.ref_names.iter().enumerate() {
        let l = format!(
            "{}\t{}\t{:.3}\t{:.3}\t{:.3}\n",
            name, lengths[i], eff_lengths[i], tpms[i], e_counts[i]
        );
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
        .to_path_buf()
        .with_additional_extension(".infreps.pq");
    let schema = Schema::from(fields);
    parquet_utils::write_chunk_to_file(output_path.to_str().unwrap(), schema, chunk)
}
