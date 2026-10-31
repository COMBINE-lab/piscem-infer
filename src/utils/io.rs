use crate::utils::parquet_utils;
use anyhow::Result;
use arrow_array::{ArrayRef, Float64Array, UInt32Array};
use libradicl::header;
use path_tools::WithAdditionalExtension;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::Arc;

use tracing::warn;

pub(crate) fn write_results(
    output: &Path,
    hdr: &header::RadHeader,
    e_counts: &[f64],
    lengths: &[u32],
    eff_lengths: &[f64],
) -> anyhow::Result<()> {
    let mut ofile = BufWriter::new(File::create(output)?);

    ofile.write_all(b"target_name\tlen\teelen\ttpm\tecount\n")?;

    const ONE_MILLION: f64 = 1_000_000.0;
    let denom: f64 = e_counts
        .iter()
        .zip(eff_lengths.iter())
        .map(|(c, l)| if *l > 0.0 { c / l } else { 0.0_f64 })
        .sum::<f64>();
    let inv_denom: f64 = if denom > 0.0 {
        ONE_MILLION / denom
    } else {
        warn!(
            "The sum of ecount / eelen for all transcripts was 0. It seems likely that no fragments were quantified. Please check the input sample!"
        );
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
        writeln!(
            ofile,
            "{}\t{}\t{:.3}\t{:.3}\t{:.3}",
            name, lengths[i], eff_lengths[i], tpms[i], e_counts[i]
        )?;
    }
    // flush explicitly, since dropping a BufWriter silently discards write errors
    ofile.flush()?;

    Ok(())
}
/// Write the fragment length distribution (`<output>.fld.pq`): one column,
/// `fragment_length_dist`, holding the count at each length.
pub(crate) fn write_fld_file(output_path: &Path, frag_lengths: Vec<u32>) -> Result<()> {
    let output_path = output_path
        .to_path_buf()
        .with_additional_extension(".fld.pq");
    let col: ArrayRef = Arc::new(UInt32Array::from(frag_lengths));
    parquet_utils::write_columns(
        &output_path,
        vec![("fragment_length_dist".to_string(), col)],
    )
}

/// Write the inferential replicates (`<output>.infreps.pq`): one column per
/// replicate, `bootstrap.<i>`, whether drawn by bootstrap or Gibbs sampling.
pub(crate) fn write_infrep_file(output_path: &Path, replicates: Vec<Vec<f64>>) -> Result<()> {
    let output_path = output_path
        .to_path_buf()
        .with_additional_extension(".infreps.pq");
    let columns = replicates
        .into_iter()
        .enumerate()
        .map(|(i, r)| {
            (
                format!("bootstrap.{i}"),
                Arc::new(Float64Array::from(r)) as ArrayRef,
            )
        })
        .collect();
    parquet_utils::write_columns(&output_path, columns)
}
