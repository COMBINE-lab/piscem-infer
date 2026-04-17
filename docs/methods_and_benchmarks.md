# piscem-infer: Methods and Benchmark Summary

This document summarizes the inference methods implemented in piscem-infer's multi-sample consensus quantification pipeline (`consensus-quant`) and their performance on standard RNA-seq benchmarks.

## Overview

piscem-infer implements a multi-sample consensus quantification pipeline that combines positional equivalence classes, structural variable selection, EC-graph-based leakage detection, and SQUAREM-accelerated EM inference to produce transcript-level abundance estimates with substantially improved accuracy and reproducibility compared to single-sample methods.

The pipeline operates in two phases:
1. **Phase 1** — Independent per-sample EM on positional equivalence classes
2. **Consensus filtering** — Multi-sample evidence integration with structural selection and leakage detection
3. **Phase 2** — Re-estimation on the filtered transcript set with warm-started EM

---

## Methods

### 1. Positional Equivalence Classes

Standard equivalence classes group fragments by their compatible transcript set. Positional ECs further partition by the fragment's relative position within each transcript, discretized into bins.

- **Default bins:** 5 (`--pos-bins N`)
- For transcript of length L with N bins, bin b spans positions `[b*L/N, (b+1)*L/N)`
- Fragments with the same transcript membership but different position bins become separate ECs
- Provides finer-grained evidence for distinguishing transcripts that share sequence

### 2. Positional Effective Lengths

When positional ECs are enabled, a single effective length per transcript is insufficient. Edge bins have fewer valid fragment start positions because fragments near transcript ends may extend beyond the boundary.

- For each (transcript, bin) pair: integrate over the fragment length distribution (FLD) to compute the number of valid start positions
- **Formula:** `pos_eff_len[t][b] = Σ_f FLD(f) * max(0, min(bin_end, L-f+1) - bin_start)`
- The M-step uses `inv_pos_eff_lens[t * n_bins + b]` instead of a single `inv_eff_lens[t]`
- Sum across bins approximates the standard effective length

### 3. Structural Variable Selection (3 Stages)

Before Phase 2, the pipeline removes transcripts that cannot be distinguished from the observed data.

**Stage 1 — Signature Collapse:** Group transcripts with identical EC signatures (same set of equivalence classes). Within each group, a single representative suffices.

**Stage 2 — Unique EC Peeling:** Iteratively identify transcripts that are the sole member of at least one EC. These are "required" — they have unique evidence. Remove their edges and cascade: newly-unique transcripts become required. Transcripts never identified as required may be removable.

**Stage 3 — Subset Dominance with Coverage Plausibility:** Transcript i is dominated if its EC signature is a strict subset of transcript j's. When positional ECs are enabled, a coverage plausibility check prevents false dominance: transcript j must use ≥2 distinct position bins across i's ECs. Single-bin coverage suggests leakage artifacts rather than true dominance.

### 4. SQUAREM EM Acceleration

The EM algorithm is accelerated using the SQUAREM method (Varadhan & Roland, 2008).

| Parameter | Default | Description |
|-----------|---------|-------------|
| Burn-in | 25 iterations | Plain EM before acceleration |
| Acceleration interval | Every 4 iterations | After burn-in |
| Step bounds | [1e-4, 10.0] | Clamp range for alpha |
| Max iterations | 1500 | Hard stop (`--max-iter`) |
| Convergence threshold | 5e-4 | Relative change (`--convergence-thresh`) |

The algorithm computes a quadratic extrapolation from three consecutive EM iterates, accepting the accelerated step only if it converges faster than ordinary EM. Parallel EM (Rayon) is used for samples with >1000 equivalence classes.

### 5. EC-Graph Leakage Detection

Annotation-free grouping of transcripts into gene-like clusters, followed by within-group leakage detection.

**Jaccard Grouping:** Two transcripts are connected if their EC signature overlap (Jaccard index) ≥ 0.10. Connected components form groups analogous to genes, without requiring annotation.

**Fraction Filter:** For each transcript, compute its abundance as a fraction of its group total. Low-fraction transcripts are candidates for leakage.

**Position CV Filter:** Compare the position-bin profile of a low-fraction transcript against the dominant group member. If Pearson correlation > 0.9 and the transcript contributes < 5% of group abundance, the transcript is flagged as leakage — its expression profile is indistinguishable from the dominant member.

### 6. Condition Rescue

When multiple experimental conditions are present (detected automatically from the sample manifest), condition rescue augments the global consensus filter:

- **Global consensus:** Transcript passes the primary cross-sample filter globally, typically support-based with adaptive EC support
- **Condition rescue:** Also keep transcripts that fail the global filter but are repeatedly present by phase-1 TPM within any single condition

This recovers condition-specific transcripts that may not meet the global threshold but have consistent within-condition abundance evidence. Rescue-only transcripts keep their phase-1 estimates after phase 2, which avoids re-suppressing rescued signal during re-estimation. Auto-enabled when multiple conditions are detected; disable with `--no-condition-rescue`.

### 7. Fragment Length Distribution Estimation

Two modes depending on data type:

- **Empirical** (default for paired-end): Estimate from the first 500K mapped fragments, with pseudocount smoothing
- **Parametric** (for single-end or manual override): Truncated normal distribution with user-specified mean and SD (`--fld-mean`, `--fld-sd`)

The FLD is used for effective length computation and positional effective length integration.

### 8. Inferential Replicates

Two mutually exclusive methods for quantifying estimation uncertainty:

**Bootstrap** (`--num-bootstraps N`): Resample EC counts via weighted alias sampling, re-run EM on each replicate. Parallelized across replicates.

**Gibbs Sampling** (`--num-gibbs-samples N`): Alternating Gamma step (draw transcript fractions from Gamma(α + count, 1/(β + effLen))) and multinomial reassignment of reads to transcripts. Per-nucleotide Dirichlet prior (α = 1e-3), β = 0.1, adaptive multi-chain (1–8 chains based on sample count), configurable thinning factor (default 5).

---

## Benchmark Results

All results generated from the committed codebase on the `pos-eq-class` branch. "pims" refers to piscem-infer multi-sample (`consensus-quant`).

### GENCODE v49 Simulation

Simulated paired-end RNA-seq from GENCODE v49 annotation using polyester, with 5000 expressed transcripts and known fold changes across 3 conditions (2 replicates each).

| Method | TP | FP | FN | Precision | Recall | F1 |
|--------|---:|---:|---:|----------:|-------:|---:|
| piscem-infer (single) | 4870 | 4436 | 130 | 0.523 | 0.974 | 0.681 |
| Salmon EM | 4778 | 2161 | 222 | 0.689 | 0.978 | 0.808 |
| Salmon VBEM | 4709 | 506 | 291 | 0.903 | 0.964 | 0.932 |
| Kallisto | 4897 | 5827 | 103 | 0.457 | 0.980 | 0.623 |
| **pims (strict)** | **4800** | **461** | **200** | **0.912** | **0.960** | **0.936** |
| pims (+rescue) | 4810 | 493 | 190 | 0.907 | 0.962 | 0.934 |

Key: pims achieves F1 = 0.936, matching Salmon VBEM (0.932) through structural filtering rather than a sparsity-inducing prior, with higher precision (0.912 vs 0.903).

### Airway Dataset (Himes et al. 2014)

8 samples of human airway smooth muscle cells (4 untreated, 4 dexamethasone-treated). Replicate concordance measured by coefficient of variation across untreated replicates.

| Method | Median CV | Transcripts with CV > 1 | N expressed (TPM ≥ 1) |
|--------|----------:|------------------------:|----------------------:|
| piscem-infer (single) | 0.391 | 5,546 | 45,682 |
| Salmon VBEM | 0.374 | 6,656 | 41,680 |
| Salmon EM | 0.376 | 5,025 | 43,810 |
| Kallisto | 0.385 | 5,402 | 45,461 |
| **pims (strict)** | **0.295** | **3,040** | **36,572** |
| pims (+rescue) | 0.309 | 3,543 | 37,552 |

Key: 21% improvement in median CV over Salmon VBEM (0.295 vs 0.374), with substantially fewer high-variance transcripts.

### SEQC/MAQC-III Titration

16 samples from the BGI site: 4 RNA mixtures (A=UHRR, B=HBRR, C=0.75A+0.25B, D=0.25A+0.75B) × 4 replicates. Tests whether predicted mixture proportions match known titration ratios.

| Method | Pearson(C) | Pearson(D) | FC slope | Mean CV |
|--------|----------:|----------:|---------:|--------:|
| piscem-infer (single) | 0.769 | 0.767 | 0.555 | 0.687 |
| Salmon VBEM | 0.758 | 0.758 | 0.572 | 0.720 |
| Salmon EM | 0.774 | 0.772 | 0.550 | 0.689 |
| Kallisto | 0.775 | 0.772 | 0.565 | 0.683 |
| **pims (strict)** | **0.995** | **0.995** | **0.991** | **0.192** |
| pims (+rescue) | 0.970 | 0.970 | 1.031 | 0.403 |

Key: Titration Pearson jumps from ~0.77 (all single-sample methods) to 0.995, with near-unity fold-change slope (0.991). Multi-sample consensus eliminates the inter-replicate noise that dominates single-sample estimates.

### TaqMan qRT-PCR Validation

Gene-level correlation against 900 TaqMan-validated genes from SEQC samples A and B.

| Method | Pearson(A) | Spearman(A) | FC Pearson | Genes detected |
|--------|----------:|----------:|----------:|---------------:|
| piscem-infer (single) | 0.531 | 0.688 | 0.925 | 711 |
| Salmon VBEM | 0.542 | 0.690 | 0.923 | 707 |
| Kallisto | 0.530 | 0.688 | 0.925 | 711 |
| **pims (+rescue)** | **0.563** | **0.693** | **0.910** | **689** |

Key: Modest improvement at gene level (Pearson 0.563 vs 0.542 for Salmon). Gene-level aggregation obscures most transcript-level differences; the strict consensus detects fewer genes (550) but rescue recovers to 689.

### Cross-Tool EC Analysis

Feeding kallisto's equivalence classes into piscem-infer's EM (to isolate mapping vs inference differences):

- **EM inference engines are equivalent:** Pearson correlation 0.999971 between piscem-infer and kallisto on identical ECs
- **SEQC titration gap** (piscem 0.6745 vs kallisto 0.6784): ~60% convergence criterion, ~40% mapping differences, 0% EM inference
- **Kallisto's extra mappings are false positives:** On simulation ground truth, kallisto has 5,827 FPs vs piscem's 4,582 (+27%). The marginal SEQC correlation advantage reflects more expressed transcripts (including FPs) adding correlated data points, not better accuracy.

---

## Default Configuration

The default `consensus-quant` configuration (as of April 2026):

- Positional ECs with 5 bins
- EC-graph leakage detection (Jaccard ≥ 0.10, position CV > 0.5, profile correlation > 0.9)
- Full structural selection (collapse + peeling + dominance with coverage plausibility)
- Condition rescue auto-enabled when multiple conditions detected
- SQUAREM acceleration (burn-in 25, threshold 5e-4)
- `--use-gene-annotation` available for legacy gene-name-based filtering (opt-in)
