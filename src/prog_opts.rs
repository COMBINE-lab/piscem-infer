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
    /// Hybrid evidence rule: EC support or ambiguity-adjusted UES dominance
    Hybrid,
}

impl FromStr for FilterMode {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "tpm" => Ok(Self::Tpm),
            "ues" => Ok(Self::Ues),
            "support" => Ok(Self::Support),
            "hybrid" => Ok(Self::Hybrid),
            other => bail!(
                "unknown filter mode '{}'; expected tpm, ues, support, or hybrid",
                other
            ),
        }
    }
}

impl Serialize for FilterMode {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self {
            Self::Tpm => serializer.serialize_str("tpm"),
            Self::Ues => serializer.serialize_str("ues"),
            Self::Support => serializer.serialize_str("support"),
            Self::Hybrid => serializer.serialize_str("hybrid"),
        }
    }
}

/// How structured condition rescue ranks candidate representatives in each EC group.
#[derive(Debug, Clone, Default)]
pub enum StructuredRescueRankMode {
    /// Rank by the largest mean TPM in any admitted condition.
    #[default]
    Peak,
    /// Prefer transcripts supported across more admitted conditions, then rank by mean TPM.
    Breadth,
}

impl FromStr for StructuredRescueRankMode {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "peak" => Ok(Self::Peak),
            "breadth" => Ok(Self::Breadth),
            other => bail!(
                "unknown structured rescue rank mode '{}'; expected peak or breadth",
                other
            ),
        }
    }
}

impl Serialize for StructuredRescueRankMode {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self {
            Self::Peak => serializer.serialize_str("peak"),
            Self::Breadth => serializer.serialize_str("breadth"),
        }
    }
}

/// How locked condition rescue chooses the fraction of Phase-1 allocation to lock.
#[derive(Debug, Clone, Default)]
pub enum ConditionRescueLockMode {
    /// Lock the same fraction for every rescued posterior allocation.
    #[default]
    Fixed,
    /// Choose full, partial, or no lock from the per-EC rescued posterior fraction.
    Confidence,
    /// Use the partial lock fraction immediately above the minimum confidence
    /// threshold, then smoothly increase to full preservation.
    FloorSmoothConfidence,
    /// Lock the lower confidence bound of the Phase-1 rescued posterior mass
    /// in each EC, using the EC count as the effective sample size.
    CredibleFloor,
    /// Lock rescued mass in proportion to its posterior enrichment over the
    /// uniform rescued-target share of the EC.
    EnrichmentFloor,
    /// Penalize phase-2 estimates that fall below the phase-1 rescued
    /// allocation floor without subtracting locked mass from EC counts.
    FloorBarrier,
    /// Apply confidence-style relaxation only to condition-local rescued
    /// transcripts; rescued transcripts with evidence in multiple conditions
    /// stay fully locked.
    GuardedConfidence,
    /// Apply confidence-style relaxation only to rescued transcripts with
    /// replicate-level instability in the current condition.
    Instability,
    /// Apply confidence-style relaxation unless a rescued transcript has
    /// stable aggregate Phase-1 count support in the current condition.
    TranscriptStability,
}

impl FromStr for ConditionRescueLockMode {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "fixed" => Ok(Self::Fixed),
            "confidence" => Ok(Self::Confidence),
            "floor-smooth-confidence" | "floor_smooth_confidence" => {
                Ok(Self::FloorSmoothConfidence)
            }
            "credible-floor" | "credible_floor" => Ok(Self::CredibleFloor),
            "enrichment-floor" | "enrichment_floor" => Ok(Self::EnrichmentFloor),
            "floor-barrier" | "floor_barrier" => Ok(Self::FloorBarrier),
            "guarded-confidence" | "guarded_confidence" => Ok(Self::GuardedConfidence),
            "instability" => Ok(Self::Instability),
            "transcript-stability" | "transcript_stability" => Ok(Self::TranscriptStability),
            other => bail!(
                "unknown condition rescue lock mode '{}'; expected fixed, confidence, floor-smooth-confidence, credible-floor, enrichment-floor, floor-barrier, guarded-confidence, instability, or transcript-stability",
                other
            ),
        }
    }
}

impl Serialize for ConditionRescueLockMode {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        match self {
            Self::Fixed => serializer.serialize_str("fixed"),
            Self::Confidence => serializer.serialize_str("confidence"),
            Self::FloorSmoothConfidence => serializer.serialize_str("floor-smooth-confidence"),
            Self::CredibleFloor => serializer.serialize_str("credible-floor"),
            Self::EnrichmentFloor => serializer.serialize_str("enrichment-floor"),
            Self::FloorBarrier => serializer.serialize_str("floor-barrier"),
            Self::GuardedConfidence => serializer.serialize_str("guarded-confidence"),
            Self::Instability => serializer.serialize_str("instability"),
            Self::TranscriptStability => serializer.serialize_str("transcript-stability"),
        }
    }
}

const PRESENCE_THRESH: f64 = 1e-8;
const RELDIFF_THRESH: f64 = 5e-4;
const MAX_EM_ITER: u32 = 1500;

fn greater_than_0(s: &str) -> std::result::Result<u32, String> {
    number_range(s, 1, u32::MAX)
}

fn fraction_0_to_1(s: &str) -> std::result::Result<f64, String> {
    let value = s
        .parse::<f64>()
        .map_err(|e| format!("failed to parse fraction '{}': {}", s, e))?;
    if value.is_finite() && (0.0..=1.0).contains(&value) {
        Ok(value)
    } else {
        Err(format!("fraction must be finite and in [0, 1], got {}", s))
    }
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
    /// disable the collapsed-EC optimization for EM/SQUAREM/Gibbs/bootstrap.
    /// When set, EM iterates over the full positional equivalence-class map.
    /// Provided for A/B timing and output comparison.
    #[arg(long, help_heading = "Advanced")]
    pub no_collapsed_ec_em: bool,
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
    #[arg(
        long,
        default_value_t = 500_000_isize,
        help_heading = "Fragment Length Distribution"
    )]
    pub param_est_frags: isize,
    /// mean of fragment length distribution
    /// (required, and used, only for unpaired fragments)
    #[arg(
        long,
        requires = "fld_sd",
        help_heading = "Fragment Length Distribution"
    )]
    pub fld_mean: Option<f64>,
    /// standard deviation of fragment length distribution
    /// (required, and used, only for unpaired fragments)
    #[arg(
        long,
        requires = "fld_mean",
        help_heading = "Fragment Length Distribution"
    )]
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
    /// disable the collapsed-EC optimization for EM/SQUAREM/Gibbs/bootstrap.
    /// When set, EM iterates over the full positional equivalence-class map.
    /// Provided for A/B timing and output comparison.
    #[arg(long, help_heading = "Advanced")]
    pub no_collapsed_ec_em: bool,

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
    #[arg(
        long,
        default_value_t = 500,
        requires = "group_lasso",
        help_heading = "Group LASSO"
    )]
    pub gl_max_iter: u32,
    /// convergence threshold for FISTA (relative objective change)
    #[arg(
        long,
        default_value_t = 1e-6,
        requires = "group_lasso",
        help_heading = "Group LASSO"
    )]
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
    /// optional TSV audit path recording per-transcript consensus decisions.
    /// Useful for diagnosing why transcripts passed or failed structural
    /// selection, global consensus, condition rescue, and gene rescue.
    #[arg(long, help_heading = "Input / Output")]
    pub consensus_audit_output: Option<PathBuf>,
    /// optional TSV audit path recording structural selection decisions.
    /// Includes the removal reason and a retained superset competitor when one
    /// is found.
    #[arg(long, requires = "txp_selection", help_heading = "Input / Output")]
    pub selection_audit_output: Option<PathBuf>,
    /// optional TSV report path for non-mutating ambiguity-group repair
    /// candidates. This computes local repair/pruning features but does not
    /// change quantification output.
    #[arg(long, help_heading = "Input / Output")]
    pub ambiguity_repair_report: Option<PathBuf>,
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
    /// support = effective EC support count;
    /// hybrid = support or ambiguity-adjusted UES dominance
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
    /// but have reproducible phase-1 TPM evidence within at least one
    /// condition. By default, rescue-only transcripts retain their phase-1
    /// estimates after phase 2. Enabled automatically when the manifest has
    /// multiple conditions. Use --no-condition-rescue to disable.
    #[arg(
        long,
        conflicts_with = "condition_aware_consensus",
        help_heading = "Consensus Filter"
    )]
    pub condition_rescue: bool,
    /// disable automatic condition rescue when multiple conditions are present.
    #[arg(long, conflicts_with_all = ["condition_rescue", "condition_aware_consensus"], help_heading = "Consensus Filter")]
    pub no_condition_rescue: bool,
    /// re-estimate condition-rescued transcripts during phase 2 instead of
    /// restoring their phase-1 estimates. The phase-2 active set remains
    /// sample-specific: strict global consensus plus transcripts rescued for
    /// the sample's condition.
    #[arg(long, requires = "condition_rescue", help_heading = "Consensus Filter")]
    pub reestimate_condition_rescue: bool,
    /// lock condition-rescued transcripts to their phase-1 per-EC posterior
    /// allocations, subtract that locked mass from each EC, then run phase-2 EM
    /// on the residual EC counts. By default this uses the recommended
    /// targeted confidence rule: fully lock ECs where rescued transcripts have
    /// at least half the posterior mass, lock 75% for moderate rescue support,
    /// and leave very weak rescue support unlocked. With a multi-condition
    /// manifest, condition rescue is enabled automatically, so this flag is the
    /// only extra flag needed to request the locked-rescue path.
    #[arg(
        long,
        conflicts_with_all = [
            "reestimate_condition_rescue",
            "no_condition_rescue",
            "condition_aware_consensus"
        ],
        help_heading = "Consensus Filter"
    )]
    pub lock_condition_rescue_allocations: bool,
    /// fraction of each condition-rescued transcript's phase-1 per-EC posterior
    /// allocation to lock. A value below 1 leaves the un-locked residual mass
    /// available to phase-2 EM over the full active set.
    #[arg(long, default_value_t = 0.75, requires = "lock_condition_rescue_allocations", value_parser = fraction_0_to_1, help_heading = "Consensus Filter")]
    pub condition_rescue_lock_fraction: f64,
    /// rule for choosing how much condition-rescued posterior allocation to
    /// lock. `fixed` uses --condition-rescue-lock-fraction everywhere;
    /// `confidence` locks fully when rescued posterior mass dominates an EC,
    /// partially when it is moderate, and not at all when it is tiny;
    /// `floor-smooth-confidence` uses the partial lock fraction just above the
    /// minimum threshold, then smoothly interpolates to full preservation;
    /// `credible-floor` locks a lower confidence bound on the rescued posterior
    /// EC mass, replacing fixed posterior-share thresholds with an EC-count
    /// uncertainty adjustment;
    /// `enrichment-floor` locks rescued mass according to its posterior
    /// enrichment over the rescued targets' uniform EC share;
    /// `floor-barrier` keeps all condition-rescued transcripts in phase 2 and
    /// applies a one-sided penalty when estimates fall below their Phase-1
    /// rescued allocation floor;
    /// `guarded-confidence` uses that rule only for condition-local rescued
    /// transcripts and fully locks rescued transcripts with evidence in multiple
    /// conditions;
    /// `instability` applies the confidence rule only to rescued transcripts
    /// with replicate-level instability in the current condition, keeping
    /// stable rescued transcripts fully locked;
    /// `transcript-stability` applies the confidence rule unless the rescued
    /// transcript has stable aggregate Phase-1 count support in the current
    /// condition.
    #[arg(long, default_value = "floor-smooth-confidence", requires = "lock_condition_rescue_allocations", value_parser = clap::value_parser!(ConditionRescueLockMode), help_heading = "Consensus Filter")]
    pub condition_rescue_lock_mode: ConditionRescueLockMode,
    /// z-score used by --condition-rescue-lock-mode credible-floor. Larger
    /// values protect only rescue mass with stronger per-EC posterior support.
    #[arg(
        long,
        default_value_t = 1.96,
        requires = "lock_condition_rescue_allocations",
        help_heading = "Consensus Filter"
    )]
    pub condition_rescue_credible_floor_z: f64,
    /// strength of the one-sided rescued-allocation floor penalty used by
    /// --condition-rescue-lock-mode floor-barrier. Larger values approximate a
    /// hard floor; 0 disables the penalty.
    #[arg(
        long,
        default_value_t = 10.0,
        requires = "lock_condition_rescue_allocations",
        help_heading = "Consensus Filter"
    )]
    pub condition_rescue_floor_barrier_weight: f64,
    /// in confidence lock mode, fully lock rescued allocation when rescued
    /// posterior mass is at least this fraction of the EC.
    #[arg(long, default_value_t = 0.5, requires = "lock_condition_rescue_allocations", value_parser = fraction_0_to_1, help_heading = "Consensus Filter")]
    pub condition_rescue_full_lock_threshold: f64,
    /// in confidence lock mode, do not lock rescued allocation when rescued
    /// posterior mass is below this fraction of the EC.
    #[arg(long, default_value_t = 0.1, requires = "lock_condition_rescue_allocations", value_parser = fraction_0_to_1, help_heading = "Consensus Filter")]
    pub condition_rescue_min_lock_threshold: f64,
    /// in instability lock mode, relax rescued transcripts whose within-condition
    /// Phase-1 count CV is at least this value.
    #[arg(
        long,
        default_value_t = 1.0,
        requires = "lock_condition_rescue_allocations",
        help_heading = "Consensus Filter"
    )]
    pub condition_rescue_instability_cv_threshold: f64,
    /// in instability lock mode, relax rescued transcripts whose within-condition
    /// sample pass fraction is below this value. The default relaxes rescued
    /// transcripts with any replicate dropout in their rescued condition.
    #[arg(long, default_value_t = 1.0, requires = "lock_condition_rescue_allocations", value_parser = fraction_0_to_1, help_heading = "Consensus Filter")]
    pub condition_rescue_instability_min_pass_fraction: f64,
    /// in guarded-confidence lock mode, count a transcript as having condition
    /// evidence when its mean Phase-1 count in that condition is at least this
    /// value. Rescued transcripts with evidence in more than one condition are
    /// fully locked.
    #[arg(
        long,
        default_value_t = 1.0,
        requires = "lock_condition_rescue_allocations",
        help_heading = "Consensus Filter"
    )]
    pub condition_rescue_guard_mean_count_threshold: f64,
    /// in transcript-stability lock mode, fully lock rescued transcripts only
    /// when their mean Phase-1 count in the current condition is at least this
    /// value.
    #[arg(
        long,
        default_value_t = 5.0,
        requires = "lock_condition_rescue_allocations",
        help_heading = "Consensus Filter"
    )]
    pub condition_rescue_stability_mean_count_threshold: f64,
    /// in transcript-stability lock mode, fully lock rescued transcripts only
    /// when their within-condition Phase-1 count CV is at most this value.
    #[arg(
        long,
        default_value_t = 0.5,
        requires = "lock_condition_rescue_allocations",
        help_heading = "Consensus Filter"
    )]
    pub condition_rescue_stability_cv_threshold: f64,
    /// in transcript-stability lock mode, fully lock rescued transcripts only
    /// when their sample pass fraction in the current condition is at least
    /// this value.
    #[arg(long, default_value_t = 1.0, requires = "lock_condition_rescue_allocations", value_parser = fraction_0_to_1, help_heading = "Consensus Filter")]
    pub condition_rescue_stability_min_pass_fraction: f64,
    /// experimental rescue mode: admit rescue at an EC-graph-group level, then
    /// emit only the dominant rescued isoforms within each admitted group.
    #[arg(long, requires = "condition_rescue", help_heading = "Consensus Filter")]
    pub structured_condition_rescue: bool,
    /// minimum summed Phase-1 TPM for an EC-graph group to be rescued within a
    /// condition. Only used with --structured-condition-rescue.
    #[arg(
        long,
        default_value_t = 3.0,
        requires = "structured_condition_rescue",
        help_heading = "Consensus Filter"
    )]
    pub structured_rescue_group_tpm_floor: f64,
    /// fraction of an admitted group's rescue TPM captured by selected
    /// isoforms. Only used with --structured-condition-rescue.
    #[arg(
        long,
        default_value_t = 0.9,
        requires = "structured_condition_rescue",
        help_heading = "Consensus Filter"
    )]
    pub structured_rescue_cumulative_frac: f64,
    /// maximum number of rescue-only isoforms emitted per admitted EC-graph
    /// group. Strict-global transcripts are not counted against this cap.
    #[arg(long, default_value_t = 3_u32, value_parser = greater_than_0, requires = "structured_condition_rescue", help_heading = "Consensus Filter")]
    pub structured_rescue_max_isoforms: u32,
    /// ranking rule for rescue candidates within each admitted EC-graph group.
    #[arg(long, default_value = "peak", value_parser = clap::value_parser!(StructuredRescueRankMode), requires = "structured_condition_rescue", help_heading = "Consensus Filter")]
    pub structured_rescue_rank_mode: StructuredRescueRankMode,
    /// include the top rescue representative for each admitted condition in
    /// each EC-graph group before applying the general group-level ranking.
    #[arg(
        long,
        requires = "structured_condition_rescue",
        help_heading = "Consensus Filter"
    )]
    pub structured_rescue_per_condition_representatives: bool,
    /// allow rescue from EC-graph groups that also contain strict-global
    /// transcripts, but only as one top raw representative per admitted
    /// condition.
    #[arg(
        long,
        requires = "structured_condition_rescue",
        help_heading = "Consensus Filter"
    )]
    pub structured_rescue_strict_group_escape: bool,
    /// minimum condition mean Phase-1 TPM for a non-strict candidate rescued
    /// from a strict-containing EC-graph group.
    #[arg(
        long,
        default_value_t = 1.0,
        requires = "structured_rescue_strict_group_escape",
        help_heading = "Consensus Filter"
    )]
    pub structured_rescue_strict_group_candidate_tpm_floor: f64,
    /// allow at most one balanced non-strict rescue candidate from an EC-graph
    /// group that also contains strict-global transcripts.
    #[arg(
        long,
        requires = "structured_condition_rescue",
        help_heading = "Consensus Filter"
    )]
    pub structured_rescue_strict_group_balanced_escape: bool,
    /// minimum condition mean Phase-1 TPM for balanced strict-group escape.
    #[arg(
        long,
        default_value_t = 0.1,
        requires = "structured_rescue_strict_group_balanced_escape",
        help_heading = "Consensus Filter"
    )]
    pub structured_rescue_strict_group_balanced_tpm_floor: f64,
    /// minimum number of conditions with signal for balanced strict-group
    /// escape.
    #[arg(
        long,
        default_value_t = 2_u32,
        value_parser = greater_than_0,
        requires = "structured_rescue_strict_group_balanced_escape",
        help_heading = "Consensus Filter"
    )]
    pub structured_rescue_strict_group_balanced_min_conditions: u32,
    /// add one Phase-1 complement representative for admitted conditions where
    /// selected rescue representatives have no signal but the EC-graph group
    /// does.
    #[arg(
        long,
        requires = "structured_condition_rescue",
        help_heading = "Consensus Filter"
    )]
    pub structured_rescue_phase1_condition_complements: bool,
    /// minimum condition mean Phase-1 TPM for a complement representative.
    #[arg(
        long,
        default_value_t = 0.1,
        requires = "structured_rescue_phase1_condition_complements",
        help_heading = "Consensus Filter"
    )]
    pub structured_rescue_phase1_complement_tpm_floor: f64,
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
    /// add a condition-specific hierarchical Dirichlet prior during phase 2,
    /// with pseudo-count strength equal to this fraction of the average sample
    /// read count. 0 disables the prior and preserves the current phase-2
    /// behavior.
    #[arg(long, default_value_t = 0.0, help_heading = "Consensus Filter")]
    pub condition_specific_prior_weight: f64,
    /// within-gene isoform fraction filter: after Phase 2 EM, zero out isoforms
    /// contributing less than this fraction of their gene's total estimated count.
    /// Removes EM leakage into sibling isoforms. Applied within EC-graph-derived
    /// groups (default) or within annotated genes (with --use-gene-annotation).
    /// Set to 0 to disable.
    #[arg(long, default_value_t = 0.01, help_heading = "Consensus Filter")]
    pub gene_fraction_filter: f64,
    /// use gene names parsed from pipe-delimited transcript IDs (GENCODE format,
    /// field 6) for within-gene leakage filtering instead of the default
    /// EC-graph-based grouping. The EC-graph filter is generally more accurate
    /// as it also captures cross-gene leakage.
    #[arg(long, help_heading = "Consensus Filter")]
    pub use_gene_annotation: bool,
    /// diagnostic mode: do not zero condition-rescue-only transcripts in the
    /// post-EM leakage filter. This tests whether rescue dropout is being
    /// introduced downstream of locked/partial allocation.
    #[arg(long, help_heading = "Consensus Filter")]
    pub preserve_condition_rescue_leakage: bool,
    /// diagnostic mode: exempt condition-rescue-only transcripts from the
    /// post-EM fraction leakage filter only when local pairwise EM against the
    /// dominant competitor gives direct sample-local support. Position-profile
    /// leakage is still removed.
    #[arg(long, help_heading = "Consensus Filter")]
    pub selective_condition_rescue_leakage: bool,
    /// minimum local pairwise EM count required for
    /// --selective-condition-rescue-leakage.
    #[arg(
        long,
        default_value_t = 1.0,
        requires = "selective_condition_rescue_leakage",
        help_heading = "Consensus Filter"
    )]
    pub selective_condition_rescue_leakage_min_count: f64,
    /// minimum local pairwise EM fraction required for
    /// --selective-condition-rescue-leakage.
    #[arg(long, default_value_t = 0.05, requires = "selective_condition_rescue_leakage", value_parser = fraction_0_to_1, help_heading = "Consensus Filter")]
    pub selective_condition_rescue_leakage_min_fraction: f64,
    /// diagnostic mode: exempt condition-rescue-only transcripts from the
    /// post-EM fraction leakage filter when their sample-local Phase-1 count is
    /// at least this value. Position-profile leakage is still removed. Set to
    /// 0 to disable.
    #[arg(long, default_value_t = 0.0, help_heading = "Consensus Filter")]
    pub condition_rescue_leakage_phase1_floor_count: f64,

    // --- Fragment Length Distribution ---
    /// number of (unique) mappings to use for fragment length distribution estimation
    #[arg(
        long,
        default_value_t = 500_000_isize,
        help_heading = "Fragment Length Distribution"
    )]
    pub param_est_frags: isize,
    /// mean of fragment length distribution
    /// (required, and used, only for unpaired fragments)
    #[arg(
        long,
        requires = "fld_sd",
        help_heading = "Fragment Length Distribution"
    )]
    pub fld_mean: Option<f64>,
    /// standard deviation of fragment length distribution
    /// (required, and used, only for unpaired fragments)
    #[arg(
        long,
        requires = "fld_mean",
        help_heading = "Fragment Length Distribution"
    )]
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
    /// disable the collapsed-EC optimization for EM/SQUAREM/Gibbs/bootstrap.
    /// When set, EM iterates over the full positional equivalence-class map.
    /// Provided for A/B timing and output comparison.
    #[arg(long, help_heading = "Advanced")]
    pub no_collapsed_ec_em: bool,
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
