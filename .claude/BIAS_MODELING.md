# Bias Modeling Implementation Plan

## Context

NG-CVB inference has been implemented and evaluated. At transcriptome scale (245K transcripts, 5K expressed), plain EM remains the best balanced option for point estimation accuracy. NG-CVB provides calibrated posterior uncertainty (Dirichlet posterior) and better global ranking metrics (Spearman), but at the cost of increased false negatives for low-abundance transcripts due to the Dirichlet prior's "rich get richer" dynamic.

The main conclusion: **improvements in quantification accuracy will come from modeling improvements (bias correction), not from changing the inference algorithm.** The current model assumes uniform fragment coverage and no library-prep biases — these are the biggest sources of systematic error on real data.

## NG-CVB Status (Implemented)

- Core algorithm in `src/utils/ngcvb.rs` — VBEM-style collapsed VB with ELBO monitoring
- CLI flags: `--use-ngcvb`, `--ngcvb-alpha` (default 1e-5)
- Convergence threshold automatically tightened to min(user_thresh, 1e-6) for ELBO-based criterion
- 5 unit tests passing
- Evaluated on sim_hard (200 txp) and sim_data_gencode (245K txp)
- **sim_hard**: NG-CVB competitive with EM (better Spearman +0.017, MARD -0.15, slightly worse FN +1%)
- **sim_data_gencode**: NG-CVB better global Spearman (+0.05), worse FN rate (5.9% vs 2.9%), worse expressed-Pearson (0.872 vs 0.925)
- Selection + NG-CVB compounds FN; α sweep showed no single α balances FP/FN at transcriptome scale
- Keeping NG-CVB as opt-in alternative; EM remains default

## Planned Implementation Order

1. Positional bias (`--pos-bias`)
2. GC bias (`--gc-bias --ref-seqs <fasta>`)
3. Sequence-context bias via VLMM (`--seq-bias --ref-seqs <fasta>`)

All bias models are optional CLI flags that modify the likelihood before inference (EM or NG-CVB).

## 1. Positional Bias

**Goal:** Model non-uniform fragment coverage along transcripts.

**What it corrects:** Library prep protocols (e.g. polyA selection, fragmentation) create systematic positional biases — fragments are more likely to originate from certain relative positions within a transcript (e.g. 3' bias in polyA-selected RNA-seq, or depletion near transcript ends).

**Approach:**
- Learn an observed/expected positional distribution from uniquely-mapped fragments
- Bin relative positions (e.g. 5' and 3' positions as fraction of transcript length) into ~20 bins
- Compute obs/exp ratio per bin from uniquely-mapped reads
- Apply as multiplicative weights to conditional probabilities p_t|k before EQ class construction
- No external data needed — uses only RAD file information

**Key design decisions:**
- Use relative position (0.0–1.0) to handle varying transcript lengths
- Separate 5' and 3' models for paired-end reads
- Smooth the obs/exp curve (kernel smoothing or spline) to avoid overfitting
- Short transcripts (< ~200bp) should be excluded from learning (positional bins too coarse)
- The bias model modifies effective lengths and conditional probabilities, not the inference algorithm

**Complexity:** Low — no external data, straightforward obs/exp ratio learning.

## 2. GC Bias

**Goal:** Correct for GC-content-dependent fragment representation.

**What it corrects:** PCR amplification and other library prep steps preferentially amplify or deplete fragments based on their GC content. This creates systematic over/under-counting that is transcript-dependent (since transcripts have different GC content profiles).

**Approach:**
- Requires reference FASTA (`--ref-seqs <fasta>`)
- For each possible fragment position on each transcript, compute the fragment GC content
- Learn obs/exp GC distribution from uniquely-mapped fragments:
  - Observed: GC content distribution of mapped fragments
  - Expected: GC content distribution of all possible fragments (weighted by transcript abundance from initial EM round)
- Compute per-bin weight w(gc) = obs(gc) / exp(gc)
- Apply weights to conditional probabilities: p_t|k → p_t|k · w(gc_fragment)
- Re-run EM with corrected probabilities (iterative: estimate abundance → update expected GC → re-estimate weights → re-run EM)

**Key design decisions:**
- GC bins: ~100 bins covering 0%–100% GC, or continuous via polynomial/spline fit
- Fragment GC computed from reference sequence at fragment start/end positions (from RAD alignment coordinates)
- Iterative refinement: 2–3 rounds of bias estimation + EM typically suffices
- Handle edge cases: very short fragments, fragments near transcript boundaries
- The salmon approach (Love et al. / Patro et al.) uses a conditional GC model that accounts for transcript-level GC and fragment-level GC jointly

**Complexity:** Medium — requires reference FASTA parsing, GC content computation per fragment position, iterative re-estimation.

## 3. Sequence-Context Bias via VLMM

**Goal:** Model sequence-specific biases at fragment start and end positions.

**What it corrects:** Enzymatic steps in library preparation (e.g. random hexamer priming, transposase insertion) have sequence preferences that create biased fragment start/end positions depending on the local sequence context.

**Approach (VLMM, NOT k-mer counting):**
- Requires reference FASTA (`--ref-seqs <fasta>`)
- Variable-Length Markov Model (VLMM) captures context-dependent nucleotide preferences at fragment boundaries
- Learn from uniquely-mapped fragments:
  - Extract reference sequence context windows (e.g. ±10bp) around fragment 5' start and 3' end positions
  - Build foreground VLMM from observed fragment start/end contexts
  - Build background VLMM from all possible start/end positions (weighted by transcript abundance)
  - Bias weight = foreground probability / background probability
- Apply weights to conditional probabilities before EM
- Context is looked up from the **reference sequence**, not from read sequences

**Why VLMM, not k-mer counting:**
- Simple k-mer counting is known to perform poorly (too rigid, doesn't capture variable-length dependencies)
- VLMM adapts context length based on data: uses longer contexts where the data supports it, shorter where it doesn't
- Captures the hierarchical nature of sequence preferences (e.g. strong preference for certain dinucleotides, weaker but real trinucleotide preferences)
- This is the approach used in salmon's sequence bias model

**Key design decisions:**
- Separate models for 5' (fragment start) and 3' (fragment end) biases
- Context window size: typically 10–20bp on each side
- VLMM depth: typically 5–7 levels
- Pruning: remove contexts with insufficient observations (use parent context instead)
- Like GC bias, this is iterative: estimate abundance → update background model → re-estimate bias → re-run EM

**Complexity:** High — VLMM implementation, reference sequence lookups, iterative re-estimation. Most complex of the three bias models.

## Integration Architecture

All bias models modify the **conditional probabilities** p_t|k and/or **effective lengths** before or during inference:

```
RAD file → parse → build initial EQ map → [optional: initial EM for abundance estimates]
  → [pos bias: learn obs/exp from unique mappers, apply weights]
  → [GC bias: compute fragment GC, learn obs/exp, apply weights]
  → [seq bias: build VLMM from unique mappers, apply weights]
  → rebuild EQ map with corrected probabilities
  → run final inference (EM or NG-CVB)
  → [optional: iterate bias estimation 2–3 times]
```

The bias weights are multiplicative on p_t|k:
```
p_t|k_corrected = p_t|k · w_pos(rel_position) · w_gc(gc_content) · w_seq(context)
```

This means the inference algorithm (EM or NG-CVB) doesn't change — only its inputs are corrected.

## Bias-Aware Simulations

Current simulations (polyester) generate reads without any biases, so they're the best case for plain EM. To validate bias models, we need:

1. **Positional bias simulation:** Generate reads with 3' bias (exponential decay from 3' end)
2. **GC bias simulation:** Generate reads with GC-dependent sampling probabilities (e.g. bell curve centered at 50% GC)
3. **Sequence-context bias simulation:** Generate reads with context-dependent start/end preferences

These should be built after the bias models are implemented, as the simulation parameters should match the model assumptions for validation. The simulation pipeline in `simulation/` can be extended for this.

## Novel Integration Ideas (C1–C4)

These are more speculative ideas that may be worth exploring once the core bias models are in place:

### C1: Bias-Aware RF Equivalence Classes
Currently, RangeFactorized EQ classes bin conditional probabilities. With bias correction, the probabilities change across iterations. RF bins could be made bias-aware by recomputing bins after bias correction, or by designing bins that are stable under expected bias perturbations.

### C2: Structural + Bayesian Sparsity
Combine transcript selection (structural filtering based on EQ graph analysis) with NG-CVB's Dirichlet posterior (statistical filtering based on evidence). Could use selection to set informative per-transcript priors (strong prior for structurally supported transcripts, weak for others) rather than post-hoc masking.

### C3: Iterative Bias Re-estimation
Rather than fixed-point bias correction followed by EM, interleave bias re-estimation with EM iterations. Each EM iteration updates abundance estimates, which update the expected distributions for bias models, which update the conditional probabilities. Converge jointly.

### C4: Grouped Priors from Signature Groups
Transcript selection identifies signature groups (transcripts with identical EQ class signatures). For NG-CVB, use group-level priors: transcripts in the same signature group share a prior, with the group-level α learned empirically. This provides better regularization than a uniform per-transcript α.

## CLI Design

```
# Positional bias only
piscem-infer quant -i sample -l ISF -o out --pos-bias

# GC bias (requires reference)
piscem-infer quant -i sample -l ISF -o out --gc-bias --ref-seqs transcripts.fa

# All bias models
piscem-infer quant -i sample -l ISF -o out --pos-bias --gc-bias --seq-bias --ref-seqs transcripts.fa

# Bias models + NG-CVB
piscem-infer quant -i sample -l ISF -o out --pos-bias --gc-bias --seq-bias --ref-seqs transcripts.fa --use-ngcvb
```

The `--ref-seqs` flag is required by `--gc-bias` and `--seq-bias` but not by `--pos-bias`.
