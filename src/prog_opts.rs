use anyhow::bail;
use clap::Args;
use clap::{Parser, Subcommand};
use clap_num::number_range;

use serde::{Serialize, Serializer};
use std::path::PathBuf;
use std::str::FromStr;

use crate::utils::map_record_types::LibraryType;
use crate::utils::txp_selection::SelectionStages;

/// Filter mode for consensus-quant: how to decide if a transcript is "expressed" in a sample.
#[derive(Debug, Clone, Default)]
pub enum FilterMode {
    /// Simple TPM threshold (default)
    #[default]
    Tpm,
    /// Unique Evidence Score: count-weighted average posterior share per EC
    Ues,
    /// Effective EC support count: number of ECs contributing non-trivially
    Support,
}

impl FromStr for FilterMode {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "tpm" => Ok(Self::Tpm),
            "ues" => Ok(Self::Ues),
            "support" => Ok(Self::Support),
            other => bail!("unknown filter mode '{}'; expected tpm, ues, or support", other),
        }
    }
}

impl Serialize for FilterMode {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self {
            Self::Tpm => serializer.serialize_str("tpm"),
            Self::Ues => serializer.serialize_str("ues"),
            Self::Support => serializer.serialize_str("support"),
        }
    }
}

const PRESENCE_THRESH: f64 = 1e-8;
const RELDIFF_THRESH: f64 = 5e-4;
const MAX_EM_ITER: u32 = 1500;

fn greater_than_0(s: &str) -> std::result::Result<u32, String> {
    number_range(s, 1, u32::MAX)
}

fn parse_selection_stages(s: &str) -> std::result::Result<SelectionStages, String> {
    SelectionStages::from_str_list(s)
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
    /// convergence threshold for EM
    #[arg(long, default_value_t = RELDIFF_THRESH, help_heading = "EM Algorithm")]
    pub convergence_thresh: f64,
    /// presence threshold for EM
    #[arg(long, default_value_t = PRESENCE_THRESH, help_heading = "EM Algorithm")]
    pub presence_thresh: f64,
    /// disable SQUAREM acceleration for the EM solver
    #[arg(long, help_heading = "EM Algorithm")]
    pub no_squarem: bool,
    /// number of coverage-smoothing EM rounds (0 = disabled).
    /// After the initial EM, re-estimates with position-dependent weights
    /// that encourage uniform transcript coverage.
    #[arg(long, default_value_t = 0, help_heading = "EM Algorithm")]
    pub coverage_smooth_rounds: usize,
    /// smoothing parameter for coverage weights (higher = gentler smoothing).
    /// Controls how aggressively non-uniform coverage is penalized.
    #[arg(long, default_value_t = 1.0, help_heading = "EM Algorithm")]
    pub coverage_epsilon: f64,
    // --- Fragment Length Distribution ---
    /// number of (unique) mappings to use to perform initial coarse-grained
    /// estimation of the fragment length distribution. These fragments will have
    /// to be read from the file and interrogated twice.
    #[arg(long, default_value_t = 500_000_isize, help_heading = "Fragment Length Distribution")]
    pub param_est_frags: isize,
    /// mean of fragment length distribution mean
    /// (required, and used, only in the case of unpaired fragments).
    #[arg(long, requires = "fld_sd", help_heading = "Fragment Length Distribution")]
    pub fld_mean: Option<f64>,
    /// mean of fragment length distribution standard deviation
    /// (required, and used, only in the case of unpaired fragments).
    #[arg(long, requires = "fld_mean", help_heading = "Fragment Length Distribution")]
    pub fld_sd: Option<f64>,

    // --- Inferential Replicates ---
    /// number of bootstrap replicates to perform.
    /// Mutually exclusive with --num-gibbs-samples.
    #[arg(long, default_value_t = 0, conflicts_with = "num_gibbs_samples", help_heading = "Inferential Replicates")]
    pub num_bootstraps: usize,
    /// number of Gibbs samples to draw for posterior uncertainty estimation.
    /// Mutually exclusive with --num-bootstraps.
    #[arg(long, default_value_t = 0, conflicts_with = "num_bootstraps", help_heading = "Inferential Replicates")]
    pub num_gibbs_samples: usize,
    /// number of internal Gibbs iterations between collected samples (thinning).
    /// Only used when --num-gibbs-samples > 0.
    #[arg(long, default_value_t = 5, help_heading = "Inferential Replicates")]
    pub gibbs_thinning_factor: usize,

    // --- Advanced ---
    /// number of probability bins to use in RangeFactorized equivalence classes.
    /// If this value is set to 1, then basic equivalence classes are used.
    #[arg(long, default_value_t = 64_u32, value_parser=greater_than_0, help_heading = "Advanced")]
    pub factorized_eqc_bins: u32,
    /// number of positional bins for positional equivalence classes.
    /// Fragments at different relative positions on the same transcripts
    /// become different ECs, improving isoform disambiguation. (1 = disabled)
    #[arg(long, default_value_t = 5_u32, value_parser=greater_than_0, help_heading = "Advanced")]
    pub pos_bins: u32,
    /// number of threads to use (used during the EM and for bootstrapping)
    #[arg(long, default_value_t = 16, help_heading = "Advanced")]
    pub num_threads: usize,
    /// number of mapped reads to sample for automatic library type detection
    /// (only used when --lib-type is set to 'auto')
    #[arg(long, default_value_t = 10_000, help_heading = "Advanced")]
    pub auto_detect_samples: usize,
    /// optional file listing transcript IDs to keep during inference
    /// (one transcript ID per line). Transcripts not listed are masked out
    /// during EM by setting their effective lengths to 0.
    #[arg(long, help_heading = "Advanced")]
    pub transcript_mask: Option<PathBuf>,
    /// enable transcript variable selection (removes structurally redundant
    /// transcripts before EM using EC graph analysis)
    #[arg(long, help_heading = "Advanced")]
    pub txp_selection: bool,
    /// which selection stages to run (comma-separated: collapse,peeling,dominance).
    /// Only used when --txp-selection is enabled. Default: all stages.
    #[arg(long, requires = "txp_selection", value_parser = parse_selection_stages, help_heading = "Advanced")]
    pub selection_stages: Option<SelectionStages>,
}

#[derive(Args, Serialize, Clone, Debug)]
pub struct MultiQuantOpts {
    // --- Input / Output ---
    /// path to a manifest file listing samples (CSV, JSON, or YAML).
    /// Format detected by extension (.csv, .json, .yaml/.yml).
    /// CSV columns: sample_name, condition, rad_path, output_dir
    #[arg(short, long, help_heading = "Input / Output")]
    pub manifest: PathBuf,
    /// the expected library type (or 'auto' for automatic detection)
    #[arg(short, long, value_parser = clap::value_parser!(LibTypeArg), help_heading = "Input / Output")]
    pub lib_type: LibTypeArg,
    /// global output directory for joint results and intermediate files
    #[arg(short, long, help_heading = "Input / Output")]
    pub output: PathBuf,

    // --- EM Algorithm ---
    /// max EM iterations per sample (determines presence mask and L-BFGS init)
    #[arg(long, default_value_t = MAX_EM_ITER, help_heading = "EM Algorithm")]
    pub max_em_iter: u32,
    /// convergence threshold for EM warm-start
    #[arg(long, default_value_t = RELDIFF_THRESH, help_heading = "EM Algorithm")]
    pub convergence_thresh: f64,
    /// presence threshold for EM
    #[arg(long, default_value_t = PRESENCE_THRESH, help_heading = "EM Algorithm")]
    pub presence_thresh: f64,

    // --- Hierarchical ---
    /// number of outer hierarchical EM iterations
    #[arg(long, default_value_t = 7, help_heading = "Hierarchical")]
    pub num_outer_iters: u32,
    /// prior weight (κ): fraction of total sample reads added as pseudo-counts
    /// from the hierarchical prior. E.g. 0.1 = 10% of reads as pseudo-counts.
    /// Higher values give the prior more influence.
    #[arg(long, default_value_t = 0.25, help_heading = "Hierarchical")]
    pub prior_weight: f64,
    /// minimum fraction of replicates (within any condition) in which a transcript
    /// must be present to enter the consensus support set. Transcripts below this
    /// threshold are zeroed out. E.g. 1.0 = all replicates, 0.67 = 2/3 majority.
    #[arg(long, default_value_t = 0.67, help_heading = "Hierarchical")]
    pub consensus_thresh: f64,
    /// use spike-and-slab prior instead of hard consensus filtering.
    /// Computes soft inclusion probabilities per transcript, updated each iteration.
    #[arg(long, help_heading = "Hierarchical")]
    pub spike_slab: bool,
    /// use per-transcript variance-adaptive shrinkage (moderated variances).
    /// Estimates cross-sample variance for each transcript and reduces
    /// shrinkage for high-variance (potentially DE) transcripts.
    #[arg(long, help_heading = "Hierarchical")]
    pub adaptive_variance: bool,

    // --- Fragment Length Distribution ---
    /// number of (unique) mappings to use for fragment length distribution estimation
    #[arg(long, default_value_t = 500_000_isize, help_heading = "Fragment Length Distribution")]
    pub param_est_frags: isize,
    /// mean of fragment length distribution
    /// (required, and used, only for unpaired fragments)
    #[arg(long, requires = "fld_sd", help_heading = "Fragment Length Distribution")]
    pub fld_mean: Option<f64>,
    /// standard deviation of fragment length distribution
    /// (required, and used, only for unpaired fragments)
    #[arg(long, requires = "fld_mean", help_heading = "Fragment Length Distribution")]
    pub fld_sd: Option<f64>,

    // --- Advanced ---
    /// number of probability bins for RangeFactorized equivalence classes (1 = basic)
    #[arg(long, default_value_t = 64_u32, value_parser = greater_than_0, help_heading = "Advanced")]
    pub factorized_eqc_bins: u32,
    /// number of positional bins for positional equivalence classes.
    /// Fragments at different relative positions on the same transcripts
    /// become different ECs, improving isoform disambiguation. (1 = disabled)
    #[arg(long, default_value_t = 5_u32, value_parser = greater_than_0, help_heading = "Advanced")]
    pub pos_bins: u32,
    /// number of threads to use
    #[arg(long, default_value_t = 16, help_heading = "Advanced")]
    pub num_threads: usize,
    /// number of mapped reads to sample for automatic library type detection
    #[arg(long, default_value_t = 10_000, help_heading = "Advanced")]
    pub auto_detect_samples: usize,
    /// enable transcript variable selection (removes structurally redundant
    /// transcripts before EM using EC graph analysis)
    #[arg(long, help_heading = "Advanced")]
    pub txp_selection: bool,
    /// which selection stages to run (comma-separated: collapse,peeling,dominance).
    /// Only used when --txp-selection is enabled. Default: all stages.
    #[arg(long, requires = "txp_selection", value_parser = parse_selection_stages, help_heading = "Advanced")]
    pub selection_stages: Option<SelectionStages>,

    // --- Group LASSO ---
    /// use group LASSO sparse inference instead of hierarchical EM.
    /// Encourages entire transcript rows to go to zero across all samples
    /// via a convex L2-norm penalty on the abundance matrix.
    #[arg(long, conflicts_with_all = ["spike_slab", "prior_weight"], help_heading = "Group LASSO")]
    pub group_lasso: bool,
    /// regularization parameter λ for group LASSO penalty.
    /// If not set, selected automatically via BIC.
    #[arg(long, requires = "group_lasso", help_heading = "Group LASSO")]
    pub gl_lambda: Option<f64>,
    /// max iterations for FISTA proximal gradient descent
    #[arg(long, default_value_t = 500, requires = "group_lasso", help_heading = "Group LASSO")]
    pub gl_max_iter: u32,
    /// convergence threshold for FISTA (relative objective change)
    #[arg(long, default_value_t = 1e-6, requires = "group_lasso", help_heading = "Group LASSO")]
    pub gl_convergence_thresh: f64,
    /// apply group sparsity per-condition instead of across all samples.
    /// Allows condition-specific transcript expression patterns.
    #[arg(long, requires = "group_lasso", help_heading = "Group LASSO")]
    pub gl_per_condition: bool,

    // --- Workflow ---
    /// only run Phase A (per-sample EQ class building + serialization)
    #[arg(long, conflicts_with = "phase_b_only", help_heading = "Workflow")]
    pub phase_a_only: bool,
    /// only run Phase B (joint hierarchical inference from serialized EQ classes)
    #[arg(long, conflicts_with = "phase_a_only", help_heading = "Workflow")]
    pub phase_b_only: bool,
}

#[derive(Args, Serialize, Clone, Debug)]
pub struct ConsensusQuantOpts {
    // --- Input / Output ---
    /// path to a manifest file listing samples (CSV, JSON, or YAML).
    /// Format detected by extension (.csv, .json, .yaml/.yml).
    /// CSV columns: sample_name, condition, rad_path, output_dir
    #[arg(short, long, help_heading = "Input / Output")]
    pub manifest: PathBuf,
    /// the expected library type (or 'auto' for automatic detection)
    #[arg(short, long, value_parser = clap::value_parser!(LibTypeArg), help_heading = "Input / Output")]
    pub lib_type: LibTypeArg,

    // --- EM Algorithm ---
    /// max iterations to run the EM
    #[arg(long, default_value_t = MAX_EM_ITER, help_heading = "EM Algorithm")]
    pub max_iter: u32,
    /// convergence threshold for EM
    #[arg(long, default_value_t = RELDIFF_THRESH, help_heading = "EM Algorithm")]
    pub convergence_thresh: f64,
    /// presence threshold for EM
    #[arg(long, default_value_t = PRESENCE_THRESH, help_heading = "EM Algorithm")]
    pub presence_thresh: f64,
    /// phase-1 override for the EM iteration cap. If unset, uses --max-iter.
    #[arg(long, help_heading = "EM Algorithm")]
    pub phase1_max_iter: Option<u32>,
    /// phase-1 override for the EM convergence threshold. If unset, uses
    /// --convergence-thresh.
    #[arg(long, help_heading = "EM Algorithm")]
    pub phase1_convergence_thresh: Option<f64>,
    /// disable SQUAREM acceleration for the phase-1 EM solver
    #[arg(long, help_heading = "EM Algorithm")]
    pub no_phase1_squarem: bool,
    /// number of coverage-smoothing EM rounds in Phase 1 (0 = disabled).
    #[arg(long, default_value_t = 0, help_heading = "EM Algorithm")]
    pub coverage_smooth_rounds: usize,
    /// smoothing parameter for coverage weights (higher = gentler).
    #[arg(long, default_value_t = 5.0, help_heading = "EM Algorithm")]
    pub coverage_epsilon: f64,
    /// phase-2 override for the EM iteration cap. If unset, uses --max-iter.
    #[arg(long, help_heading = "EM Algorithm")]
    pub phase2_max_iter: Option<u32>,
    /// phase-2 override for the EM convergence threshold. If unset, uses
    /// --convergence-thresh.
    #[arg(long, help_heading = "EM Algorithm")]
    pub phase2_convergence_thresh: Option<f64>,

    // --- Consensus Filter ---
    /// filter mode: how to decide if a transcript is "expressed" in a sample.
    /// tpm = simple TPM threshold; ues = unique evidence score;
    /// support = effective EC support count
    #[arg(long, default_value = "tpm", value_parser = clap::value_parser!(FilterMode), help_heading = "Consensus Filter")]
    pub filter_mode: FilterMode,
    /// minimum fraction of samples in which a transcript must be expressed
    /// to pass the consensus filter. Default: (N-1)/N.
    #[arg(long, help_heading = "Consensus Filter")]
    pub min_fraction: Option<f64>,
    /// require consensus within any one condition instead of globally across
    /// all samples. A transcript passes if it is supported in enough
    /// replicates within at least one condition.
    #[arg(long, help_heading = "Consensus Filter")]
    pub condition_aware_consensus: bool,
    /// use strict global consensus, then rescue transcripts that fail globally
    /// but pass within at least one condition. Combines high precision of
    /// global filtering with preservation of condition-specific expression.
    #[arg(long, conflicts_with = "condition_aware_consensus", help_heading = "Consensus Filter")]
    pub condition_rescue: bool,
    /// TPM threshold above which a transcript is considered expressed
    /// in a given sample. Only used with --filter-mode tpm. (default: 0.0)
    #[arg(long, default_value_t = 0.0, help_heading = "Consensus Filter")]
    pub tpm_threshold: f64,
    /// UES threshold above which a transcript is considered expressed.
    /// Only used with --filter-mode ues. (default: 0.01)
    #[arg(long, default_value_t = 0.01, help_heading = "Consensus Filter")]
    pub ues_threshold: f64,
    /// minimum number of ECs with non-trivial contribution for a transcript
    /// to be considered expressed. Only used with --filter-mode support. (default: 2)
    #[arg(long, default_value_t = 2, help_heading = "Consensus Filter")]
    pub min_ec_support: u32,
    /// scale the EC support threshold per transcript based on its ambiguity
    /// (average EC size). Transcripts in highly shared EC neighborhoods require
    /// more supporting ECs. Threshold = max(min_ec_support, ceil(log2(avg_ec_size))).
    #[arg(long, help_heading = "Consensus Filter")]
    pub adaptive_ec_support: bool,
    /// minimum assigned fragment count from an EC for it to count toward
    /// a transcript's support. Used by ues and support modes. (default: 0.5)
    #[arg(long, default_value_t = 0.5, help_heading = "Consensus Filter")]
    pub min_support_count: f64,
    /// disable phase-2 warm starts from the phase-1 abundance estimates.
    #[arg(long, help_heading = "Consensus Filter")]
    pub no_phase2_warm_start: bool,
    /// within-gene isoform fraction filter: after Phase 2 EM, zero out isoforms
    /// contributing less than this fraction of their gene's total estimated count.
    /// Removes EM leakage into sibling isoforms. Gene names are parsed from
    /// pipe-delimited transcript IDs (GENCODE format, field 6). Set to 0 to disable.
    #[arg(long, default_value_t = 0.01, help_heading = "Consensus Filter")]
    pub gene_fraction_filter: f64,
    /// disable gene name parsing from transcript IDs, forcing the annotation-free
    /// EC-graph-based leakage filter for all transcripts. Useful when transcript
    /// names don't follow GENCODE pipe-delimited format.
    #[arg(long, help_heading = "Consensus Filter")]
    pub no_gene_annotation: bool,

    // --- Fragment Length Distribution ---
    /// number of (unique) mappings to use for fragment length distribution estimation
    #[arg(long, default_value_t = 500_000_isize, help_heading = "Fragment Length Distribution")]
    pub param_est_frags: isize,
    /// mean of fragment length distribution
    /// (required, and used, only for unpaired fragments)
    #[arg(long, requires = "fld_sd", help_heading = "Fragment Length Distribution")]
    pub fld_mean: Option<f64>,
    /// standard deviation of fragment length distribution
    /// (required, and used, only for unpaired fragments)
    #[arg(long, requires = "fld_mean", help_heading = "Fragment Length Distribution")]
    pub fld_sd: Option<f64>,

    // --- Advanced ---
    /// number of probability bins for RangeFactorized equivalence classes (1 = basic)
    #[arg(long, default_value_t = 64_u32, value_parser = greater_than_0, help_heading = "Advanced")]
    pub factorized_eqc_bins: u32,
    /// number of positional bins for positional equivalence classes.
    /// Fragments at different relative positions on the same transcripts
    /// become different ECs, improving isoform disambiguation. (1 = disabled)
    #[arg(long, default_value_t = 5_u32, value_parser = greater_than_0, help_heading = "Advanced")]
    pub pos_bins: u32,
    /// number of threads to use
    #[arg(long, default_value_t = 16, help_heading = "Advanced")]
    pub num_threads: usize,
    /// number of samples to process concurrently in consensus-quant.
    /// The total thread budget from --num-threads is split across these jobs.
    /// Default (0) auto-selects based on sample count and thread budget.
    #[arg(long, default_value_t = 0, help_heading = "Advanced")]
    pub sample_parallelism: u32,
    /// number of mapped reads to sample for automatic library type detection
    #[arg(long, default_value_t = 10_000, help_heading = "Advanced")]
    pub auto_detect_samples: usize,
    /// enable transcript variable selection on the merged EC graph across all
    /// samples. Removes structurally redundant transcripts before Phase 1 EM
    /// using signature collapse, unique-EC peeling, and subset dominance.
    #[arg(long, help_heading = "Advanced")]
    pub txp_selection: bool,
    /// which selection stages to run (comma-separated: collapse,peeling,dominance).
    /// Only used when --txp-selection is enabled. Default: all stages.
    #[arg(long, requires = "txp_selection", value_parser = parse_selection_stages, help_heading = "Advanced")]
    pub selection_stages: Option<SelectionStages>,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// quantify from the rad file (single sample)
    #[command(arg_required_else_help = true)]
    Quant(QuantOpts),
    /// hierarchical multi-sample quantification
    #[command(arg_required_else_help = true)]
    MultiQuant(MultiQuantOpts),
    /// consensus-filtered multi-sample quantification (two-pass EM)
    #[command(arg_required_else_help = true)]
    ConsensusQuant(ConsensusQuantOpts),
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
