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

**Data flow:** RAD file + map_info.json → parse header/metadata → estimate fragment length distribution → build equivalence class maps → run EM algorithm → output abundance estimates (TSV .quant file)

**Key modules in `src/`:**

- `main.rs` — Entry point, dispatches to the `quant` subcommand
- `prog_opts.rs` — CLI argument definitions using clap derive (`Cli`, `Commands::Quant`, `QuantOpts`)
- `process_rad.rs` — Core quantification workflow: RAD parsing, FLD estimation, EM coordination via `process_bulk()` / `process_bulk_dispatch()`
- `fld.rs` — Fragment length distribution models (empirical and parametric), `FldPDF` trait
- `utils/em.rs` — EM algorithm implementation with Rayon parallelization, bootstrap replicates, conditional means
- `utils/eq_maps.rs` — Equivalence class representations: `EqLabel` trait with `BasicEqMap` and `RangeFactorizedEqMap` (probability binning) implementations
- `utils/io.rs` — Output writing (TSV results, Parquet FLD/bootstrap files via arrow2)
- `utils/map_record_types.rs` — Library type enums (SF, ISF, SR, ISR, U, IU) and compatibility checking; contains unit tests

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
