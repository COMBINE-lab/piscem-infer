use clap::Args;
use clap::{Parser, Subcommand};

use serde::Serialize;
use std::path::PathBuf;

use crate::utils::map_record_types::LibraryType;

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
    #[arg(short, long, default_value_t = 1500)]
    pub max_iter: u32,
    /// convergence threshold for EM
    #[arg(long, default_value_t = 1e-3)]
    pub convergence_thresh: f64,
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
    #[arg(short, long)]
    pub quiet: bool,
    #[command(subcommand)]
    pub command: Commands,
}
