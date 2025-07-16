use clap::Args;
use clap::{Parser, Subcommand};

use serde::Serialize;
use std::path::PathBuf;

use crate::utils::map_record_types::LibraryType;

const PRESENCE_THRESH: f64 = 1e-8;
const RELDIFF_THRESH: f64 = 1e-3;
const MAX_EM_ITER: u32 = 1500;

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
