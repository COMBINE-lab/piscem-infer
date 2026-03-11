# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

piscem-infer is a Rust CLI tool that consumes bulk-RAD files (produced by piscem/piscem-rs) and produces abundance estimates of mapped targets. It implements EM-based statistical inference for quantification of bulk sequencing data (RNA-seq transcript quantification, metagenomic abundance estimation, etc.).

## Build & Test Commands

```bash
cargo build                    # Debug build
cargo build --release          # Release build (thin LTO, panic=abort)
cargo test                     # Run all tests
cargo test -- --nocapture      # Run tests with stdout visible
cargo clippy                   # Lint
```

The binary is named `piscem-infer`. Rust edition 2024 is required.

## Architecture

**Data flow:** RAD file + map_info.json → parse header/metadata → (optional auto library type detection) → estimate fragment length distribution → build equivalence class maps → run EM algorithm → output abundance estimates (TSV .quant file) → (optional inferential replicates via bootstrap or Gibbs sampling)

**Key modules in `src/`:**

- `main.rs` — Entry point, dispatches to `quant` and `multi-quant` subcommands
- `prog_opts.rs` — CLI argument definitions using clap derive (`Cli`, `Commands::Quant`, `Commands::MultiQuant`, `QuantOpts`, `MultiQuantOpts`)
- `process_rad.rs` — Core quantification workflow: `build_eq_map_from_rad()` (shared EQ map building), `process_bulk()` / `process_bulk_dispatch()` (single-sample EM + output)
- `multi_sample.rs` — Multi-sample hierarchical quantification: manifest parsing (CSV/JSON/YAML), Phase A (per-sample EQ map building + serialization), Phase B (joint hierarchical inference)
- `fld.rs` — Fragment length distribution models (empirical and parametric), `FldPDF` trait
- `utils/em.rs` — EM algorithm implementation with Rayon parallelization, bootstrap replicates, conditional means
- `utils/gibbs.rs` — Gibbs sampler (`do_gibbs`/`gibbs_iteration`): Gamma step + multinomial reassignment over equivalence classes
- `utils/eq_maps.rs` — Equivalence class representations: `EqLabel` trait with `BasicEqMap` and `RangeFactorizedEqMap` (probability binning) implementations
- `utils/gradient.rs` — Softmax, gradient of penalized log-likelihood, diagonal Fisher information, Laplace posterior variance
- `utils/lbfgs.rs` — L-BFGS wrapper (argmin crate) for penalized MAP estimation; `penalized_em()` runs EM warm-start → L-BFGS → Fisher/Laplace
- `utils/eq_serialize.rs` — EQ class serialization to/from Parquet + JSON metadata; `serialize_eq_map()`, `deserialize_eq_map()`
- `utils/hierarchical.rs` — Hierarchical empirical Bayes M-step: precision-weighted condition means, method-of-moments biological variance
- `utils/io.rs` — Output writing (TSV results, Parquet FLD/bootstrap files via arrow2)
- `utils/map_record_types.rs` — Library type enums (SF, ISF, SR, ISR, U, IU) and compatibility checking; contains unit tests

**Inferential uncertainty:** Two mutually exclusive methods (`--num-bootstraps` vs `--num-gibbs-samples`):
- **Bootstrap**: Resamples equivalence class counts via `WeightedAliasIndex`, runs EM on each replicate. Parallelized across replicates.
- **Gibbs sampling**: Gamma step (transcript fractions from Gamma(prior+count, 1/(β+effLen))) + multinomial reassignment over equivalence classes. Per-nucleotide Dirichlet prior (α=1e-3/effLen), β=0.1. Adaptive multi-chain (1/2/4/8 chains parallelized via rayon), configurable thinning factor. Output uses salmon-style scaled Gamma fractions (`output[i] = μ[i] * effLen[i] * totalMapped / Σ(μ[j] * effLen[j])`) rather than raw count assignments, with values below 1e-8 truncated to zero (Bray et al. 2016).

Both output to `.infreps.pq` with `bootstrap.N` columns; `meta_info.json` records `infrep_method`.

**Key design patterns:**
- Generic trait-based equivalence classes (`EqLabel` trait) allow swappable strategies
- `FldPDF` trait abstracts over empirical vs parametric fragment length distributions
- Parallel EM using Rayon for multi-threaded abundance estimation
- `OnceLock<f64>` for lazy initialization of `NUM_BINS`
- `anyhow::Result` throughout for error propagation

**Core dependency:** `libradicl` (0.10.0) handles RAD file reading/parsing.
  - **NOTE**: This `libradicl` dependency is another tool written by us, expect the version to update frequently and reflect that here when it does.

## CI/CD

- **release.yml** — cargo-dist cross-platform binary builds on version tags
- **release-please.yml** — Automated versioning and crates.io publishing
- **sanitize-cargo.yml** — Runs cargo-sanitize when commit contains `[do_tag]`
- Targets: macOS (x86_64, aarch64), Linux (x86_64, aarch64)

## Documentation

Sphinx-based docs in `docs/` published to [ReadTheDocs](https://piscem-infer.readthedocs.io/). Build with `make html` from `docs/`.

## Feature Tracking

See [CLAUDE-features.md](CLAUDE-features.md) for the status of features in development.
