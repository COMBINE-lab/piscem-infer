# Feature Tracking

## Automatic Library Type Detection

**Status:** Implemented (not yet released)

**Summary:** Users can specify `--lib-type auto` to automatically detect the library type from a sample of mapped reads, instead of specifying it explicitly (SF, ISF, SR, ISR, U, IU).

**How it works:**
- Samples up to 10,000 mapped reads (configurable via `--auto-detect-samples`)
- Counts orientation frequencies (Forward, Reverse, FR, RF, etc.)
- Uses salmon-style heuristics (30%/70% thresholds) to classify as stranded or unstranded
- Paired-end: ISF / ISR / IU; Single-end: SF / SR / U
- Warns if ≥5% wrong-strand mappings (stranded) or >5% strand bias (unstranded)
- Detected type recorded in `.meta_info.json` as `inferred_lib_type`

**Files modified:**
- `src/utils/map_record_types.rs` — `OrientationCounts`, `detect_library_type()`, `check_strand_warnings()`, `Display` for `LibraryType`
- `src/prog_opts.rs` — `LibTypeArg` enum, `auto_detect_samples` CLI arg
- `src/process_rad.rs` — `detect_lib_type_from_sample()`, resolution in `process_bulk_dispatch()`

## Gibbs Sampling for Posterior Uncertainty

**Status:** Implemented (not yet released)

**Summary:** Optional Gibbs sampling as an alternative to bootstrapping for estimating inferential uncertainty in transcript abundance. Mutually exclusive with `--num-bootstraps`.

**How it works:**
- `--num-gibbs-samples N` activates Gibbs sampling (default 0, conflicts with `--num-bootstraps`)
- `--gibbs-thinning-factor T` controls internal iterations between samples (default 5)
- Initializes from EM point estimate
- Each iteration: Gamma step (draw transcript fractions) + multinomial reassignment of reads across equivalence classes
- Per-nucleotide Dirichlet prior (α=1e-3/effLen), Gamma rate = 1/(β+effLen) with β=0.1 (matching salmon)
- Adaptive multi-chain: 1 chain (<50 samples), 2 (≥50), 4 (≥100), 8 (≥200); chains run in parallel
- Output uses `bootstrap.N` column names in `.infreps.pq` for downstream compatibility
- `meta_info.json` includes `infrep_method` field ("gibbs", "bootstrap", or "none")

**Dependencies upgraded:**
- `rand` 0.9.1 → 0.10, `rand_distr` 0.5.1 → 0.6

**Files modified:**
- `src/utils/em.rs` — `do_gibbs()`, `gibbs_iteration()` functions
- `src/prog_opts.rs` — `num_gibbs_samples`, `gibbs_thinning_factor` CLI args
- `src/process_rad.rs` — Gibbs invocation block, `infrep_method` in meta_info.json
- `Cargo.toml` — Bumped rand/rand_distr versions

## Hierarchical Multi-Sample Quantification

**Status:** Implemented (not yet released, pending end-to-end validation with simulated data)

**Summary:** Multi-sample RNA-seq quantification that shares information across samples via a hierarchical empirical Bayes model. Improves abundance estimates for ambiguous and low-coverage transcripts by leveraging biological replicates within conditions.

**Algorithm:**
- **Phase A (per-sample):** Parse RAD files, build equivalence class maps, serialize to Parquet + JSON metadata
- **Phase B (joint inference):** Iterates between:
  - Per-sample penalized MAP estimation: warm-start EM (configurable iterations) followed by L-BFGS on the penalized log-likelihood F(phi) = l(softmax(phi)) - (1/2) sum_t (phi_t - nu_t)^2 / sigma_t^2
  - Hierarchical M-step: precision-weighted condition means + method-of-moments biological variance estimation with max(0,...) clamp
- Laplace approximation for per-sample posterior uncertainty (diagonal Fisher information)
- Two-phase workflow (`--phase-a-only` / `--phase-b-only`) enables distributed processing

**CLI:** `piscem-infer multi-quant --manifest <file> --output <dir> [options]`

**Manifest formats:** CSV, JSON, YAML (auto-detected by extension). Each entry specifies: sample_name, condition, rad_path, output_dir.

**Output:**
- Per-sample: `.quant` TSV (standard format), `.posterior.json` (phi_hat, sigma_hat_sq)
- Joint: `hierarchical_params.json` (nu, sigma_sq per condition/transcript), `convergence.json` (per-iteration metrics), `meta_info.json`

**Dependencies added:**
- `argmin` 0.11.0 (L-BFGS solver), `argmin-math` 0.5.1
- `csv` 1.4.0 (manifest parsing), `saphyr` 0.0.6 (YAML parsing)
- `tempfile` 3.27.0 (dev-dependency for tests)

**New files:**
- `src/utils/gradient.rs` — Softmax, gradient of penalized objective, Fisher information, Laplace variance
- `src/utils/lbfgs.rs` — L-BFGS wrapper (argmin), penalized_em() three-phase per-sample inference
- `src/utils/eq_serialize.rs` — EQ class serialization to/from Parquet + JSON metadata
- `src/utils/hierarchical.rs` — Hierarchical M-step: condition means, biological variance, convergence metrics
- `src/multi_sample.rs` — Manifest parsing (CSV/JSON/YAML), Phase A/B orchestration

**Files modified:**
- `src/process_rad.rs` — Extracted `build_eq_map_from_rad()` from `process_bulk_dispatch()`, added `RadProcessingOpts`, `EqMapBundle`
- `src/utils/eq_maps.rs` — Added `PackedEqMap::from_raw()` constructor for deserialization
- `src/utils/io.rs` — Changed `write_results()` to accept `&[String]` instead of `&RadHeader`
- `src/prog_opts.rs` — Added `MultiQuantOpts`, `Commands::MultiQuant`
- `src/main.rs` — Added `mod multi_sample` and dispatch
- `src/utils.rs` — Added module declarations for gradient, lbfgs, eq_serialize, hierarchical
