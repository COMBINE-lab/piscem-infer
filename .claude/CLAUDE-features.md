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
