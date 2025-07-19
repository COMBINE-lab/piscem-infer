use clap::Args;
use clap::{Parser, Subcommand};
use clap_num::number_range;

use serde::Serialize;
use std::path::PathBuf;

use crate::utils::map_record_types::LibraryType;

const PRESENCE_THRESH: f64 = 1e-8;
const RELDIFF_THRESH: f64 = 1e-3;
const MAX_EM_ITER: u32 = 1500;

fn greater_than_0(s: &str) -> std::result::Result<u32, String> {
    number_range(s, 1, u32::MAX)
}

#[derive(Args, Serialize, Clone, Debug)]
pub struct QuantOpts {
    /// input stem (i.e. without the .rad suffix)
    #[arg(short, long)]
    pub input: PathBuf,
    /// the expected library type
    #[arg(short, long, value_parser = clap::value_parser!(LibraryType))]
    pub lib_type: LibraryType,
    /// output file prefix (multiple output files may be created, the main will have a `.quant` suffix)
    #[arg(short, long)]
    pub output: PathBuf,
    /// max iterations to run the EM
    #[arg(short, long, default_value_t = MAX_EM_ITER)]
    pub max_iter: u32,
    /// convergence threshold for EM
    #[arg(long, default_value_t = RELDIFF_THRESH)]
    pub convergence_thresh: f64,
    /// presence threshold for EM
    #[arg(long, default_value_t = PRESENCE_THRESH)]
    pub presence_thresh: f64,
    /// number of (unique) mappings to use to perform initial coarse-grained
    /// estimation of the fragment length distribution. These fragments will have
    /// to be read from the file and interrogated twice.
    #[arg(long, default_value_t = 500_000_isize)]
    pub param_est_frags: isize,
    /// number of probability bins to use in RangeFactorized equivalence classes.
    /// If this value is set to 1, then basic equivalence classes are used.
    #[arg(long, default_value_t = 64_u32, value_parser=greater_than_0)]
    pub factorized_eqc_bins: u32,
    /// mean of fragment length distribution mean
    /// (required, and used, only in the case of unpaired fragments).
    #[arg(long, requires = "fld_sd")]
    pub fld_mean: Option<f64>,
    /// mean of fragment length distribution standard deviation
    /// (required, and used, only in the case of unpaired fragments).
    #[arg(long, requires = "fld_mean")]
    pub fld_sd: Option<f64>,
    /// number of bootstrap replicates to perform.
    #[arg(long, default_value_t = 0)]
    pub num_bootstraps: usize,
    /// number of threads to use (used during the EM and for bootstrapping)
    #[arg(long, default_value_t = 16)]
    pub num_threads: usize,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// quantify from the rad file
    #[command(arg_required_else_help = true)]
    Quant(QuantOpts),
}

/// quantify target abundance from bulk-sequencing data
#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(propagate_version = true)]
pub struct Cli {
    /// be quiet while processing and only report errors or
    /// critical log messages
    #[arg(short, long)]
    pub quiet: bool,
    #[command(subcommand)]
    pub command: Commands,
}
