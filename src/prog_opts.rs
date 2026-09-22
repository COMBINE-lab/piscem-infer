use anyhow::bail;
use clap::Args;
use clap::{Parser, Subcommand};
use clap_num::number_range;

use serde::{Serialize, Serializer};
use std::path::PathBuf;
use std::str::FromStr;

use crate::utils::map_record_types::LibraryType;

// EM convergence defaults follow salmon (`salmon_infer::EmOptions::default()`).
const PRESENCE_THRESH: f64 = 1e-8;
const RELDIFF_THRESH: f64 = 1e-2;
const ALPHA_CHECK_CUTOFF: f64 = 1e-2;
const MAX_EM_ITER: u32 = 10_000;
const DEFAULT_SEED: u64 = 0x5A15_0EED;

/// Convergence acceleration applied on top of the EM fixed-point iteration.
/// All three reach the same fixpoint; they differ in how many M-steps it takes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, clap::ValueEnum)]
#[serde(rename_all = "lowercase")]
pub enum EmAccelArg {
    /// plain fixed-point iteration (default; matches prior releases and salmon)
    None,
    /// SQUAREM (SqS3) extrapolation
    Squarem,
    /// damped Anderson acceleration with restarts
    Daarem,
}

impl From<EmAccelArg> for salmon_infer::EmAccel {
    fn from(a: EmAccelArg) -> Self {
        match a {
            EmAccelArg::None => Self::None,
            EmAccelArg::Squarem => Self::Squarem,
            EmAccelArg::Daarem => Self::Daarem,
        }
    }
}

fn greater_than_0(s: &str) -> std::result::Result<u32, String> {
    number_range(s, 1, u32::MAX)
}

#[derive(Debug, Clone)]
pub enum LibTypeArg {
    Explicit(LibraryType),
    Auto,
}

impl FromStr for LibTypeArg {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_uppercase().as_str() {
            "AUTO" => Ok(Self::Auto),
            other => match other.parse::<LibraryType>() {
                Ok(lt) => Ok(Self::Explicit(lt)),
                Err(e) => bail!("{e}"),
            },
        }
    }
}

impl Serialize for LibTypeArg {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self {
            Self::Explicit(lt) => lt.serialize(serializer),
            Self::Auto => serializer.serialize_str("Auto"),
        }
    }
}

#[derive(Args, Serialize, Clone, Debug)]
pub struct QuantOpts {
    // --- Input / Output ---
    /// input stem (i.e. without the .rad suffix)
    #[arg(short, long, help_heading = "Input / Output")]
    pub input: PathBuf,
    /// the expected library type (or 'auto' for automatic detection)
    #[arg(short, long, value_parser = clap::value_parser!(LibTypeArg), help_heading = "Input / Output")]
    pub lib_type: LibTypeArg,
    /// output file prefix (multiple output files may be created, the main will have a `.quant` suffix)
    #[arg(short, long, help_heading = "Input / Output")]
    pub output: PathBuf,

    // --- EM Algorithm ---
    /// max iterations to run the EM
    #[arg(short, long, default_value_t = MAX_EM_ITER, help_heading = "EM Algorithm")]
    pub max_iter: u32,
    /// convergence threshold for EM: stop once every target with abundance above
    /// `--alpha-check-cutoff` changes by less than this relative amount between
    /// iterations
    #[arg(long, default_value_t = RELDIFF_THRESH, help_heading = "EM Algorithm")]
    pub convergence_thresh: f64,
    /// targets whose abundance is at or below this value are ignored when
    /// checking convergence (salmon's `alphaCheckCutoff`)
    #[arg(long, default_value_t = ALPHA_CHECK_CUTOFF, help_heading = "EM Algorithm")]
    pub alpha_check_cutoff: f64,
    /// presence threshold for EM: abundances below this are truncated to zero
    /// (with mass-preserving redistribution) in the final estimate
    #[arg(long, default_value_t = PRESENCE_THRESH, help_heading = "EM Algorithm")]
    pub presence_thresh: f64,
    /// convergence acceleration scheme for the EM. `squarem`/`daarem` reach the
    /// same fixpoint in fewer M-steps but are not byte-identical to `none`.
    #[arg(long, value_enum, default_value_t = EmAccelArg::None, help_heading = "EM Algorithm")]
    pub em_accel: EmAccelArg,

    // --- Fragment Length Distribution ---
    /// number of (unique) mappings to use to perform initial coarse-grained
    /// estimation of the fragment length distribution. These fragments will have
    /// to be read from the file and interrogated twice.
    #[arg(
        long,
        default_value_t = 500_000_isize,
        help_heading = "Fragment Length Distribution"
    )]
    pub param_est_frags: isize,
    /// mean of fragment length distribution mean
    /// (required, and used, only in the case of unpaired fragments).
    #[arg(
        long,
        requires = "fld_sd",
        help_heading = "Fragment Length Distribution"
    )]
    pub fld_mean: Option<f64>,
    /// mean of fragment length distribution standard deviation
    /// (required, and used, only in the case of unpaired fragments).
    #[arg(
        long,
        requires = "fld_mean",
        help_heading = "Fragment Length Distribution"
    )]
    pub fld_sd: Option<f64>,

    // --- Inferential Replicates ---
    /// number of bootstrap replicates to perform.
    /// Mutually exclusive with --num-gibbs-samples.
    #[arg(
        long,
        default_value_t = 0,
        conflicts_with = "num_gibbs_samples",
        help_heading = "Inferential Replicates"
    )]
    pub num_bootstraps: usize,
    /// number of Gibbs samples to draw for posterior uncertainty estimation.
    /// Mutually exclusive with --num-bootstraps.
    #[arg(
        long,
        default_value_t = 0,
        conflicts_with = "num_bootstraps",
        help_heading = "Inferential Replicates"
    )]
    pub num_gibbs_samples: usize,
    /// number of internal Gibbs iterations between collected samples (thinning).
    /// Only used when --num-gibbs-samples > 0.
    #[arg(long, default_value_t = 5, help_heading = "Inferential Replicates")]
    pub gibbs_thinning_factor: usize,
    /// random seed used for bootstrap resampling and Gibbs sampling
    /// (replicates are reproducible for a fixed seed, independent of thread count)
    #[arg(long, default_value_t = DEFAULT_SEED, help_heading = "Inferential Replicates")]
    pub seed: u64,

    // --- Advanced ---
    /// number of probability bins to use in RangeFactorized equivalence classes.
    /// If this value is set to 1, then basic equivalence classes are used.
    #[arg(long, default_value_t = 64_u32, value_parser=greater_than_0, help_heading = "Advanced")]
    pub factorized_eqc_bins: u32,
    /// number of threads to use (used during the EM and for bootstrapping)
    #[arg(long, default_value_t = 16, help_heading = "Advanced")]
    pub num_threads: usize,
    /// number of mapped reads to sample for automatic library type detection
    /// (only used when --lib-type is set to 'auto')
    #[arg(long, default_value_t = 10_000, help_heading = "Advanced")]
    pub auto_detect_samples: usize,
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
