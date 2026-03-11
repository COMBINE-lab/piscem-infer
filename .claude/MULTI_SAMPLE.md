# piscem-infer: Hierarchical Quantification via Penalized EM + L-BFGS

## Project Context

This document captures the design decisions and mathematical foundations for adding
**hierarchical, multi-sample transcript quantification** to `piscem-infer`. It is
intended as persistent context for implementation sessions. Read it fully before
writing or modifying any code related to this feature.

---

## Goal

Extend piscem-infer's existing per-sample RFEQ-based EM to share information
across samples in an RNA-seq experiment. The aim is to:

1. Improve quantification accuracy for ambiguous/low-coverage transcripts by
   regularizing estimates toward experiment- and condition-level means.
2. Produce calibrated per-sample posterior uncertainty estimates suitable for
   downstream differential expression.
3. Ensure downstream DE test statistics remain valid (no double-shrinkage,
   correct variance propagation).

---

## Background: The Existing EM

piscem-infer builds **range-factorized equivalence classes (RFEQs)**. Each RFEQ `e`
has:
- A set of compatible transcripts `T_e`
- Data-driven conditional weights `w[t,e]` per transcript (range-factorized)
- A fragment count `c_e`

The per-sample log-likelihood is:

```
ℓ(θ) = Σ_e  c_e · log( Σ_{t ∈ T_e}  θ_t · w[t,e] )
```

where `θ` lives on the probability simplex (Σ θ_t = 1). The EM finds the MLE by
iterating:

```
E-step:  q[t,e] = θ_t · w[t,e] / μ_e          (μ_e = Σ_t θ_t · w[t,e])
M-step:  θ_t    ∝ Σ_e  c_e · q[t,e]
```

---

## Chosen Approach: Warm-Start L-BFGS (Option C)

### Rationale

Adding a hierarchical Gaussian prior in log-space breaks the EM M-step's closed
form. Three options were evaluated:

- **Option A** (EM + gradient correction): ad hoc, requires step-size tuning.
- **Option B** (Generalized EM / MM): cleaner but more invasive to implement.
- **Option C** (EM warm-start → L-BFGS): **selected**. Exact MAP estimation,
  integrates naturally with analytic Fisher information, handles the softmax
  Jacobian correctly, and requires minimal changes to the existing EM loop.

### The Penalized Objective

Work in **softmax / log-space** parameterization:

```
θ_t = softmax(φ)_t = exp(φ_t) / Σ_{t'} exp(φ_{t'})
```

The per-sample penalized objective is:

```
F(φ) = ℓ(softmax(φ))  -  (1/2) Σ_t  (φ_t - ν_t)² / σ_t²
```

where:
- `ν_t` = condition-level mean for transcript `t` (from outer hierarchical M-step)
- `σ_t²` = biological variance for transcript `t` (shared across samples, estimated
  empirically from data — see Outer Loop below)

### Gradient of F(φ)

```
∂F/∂φ_t = Σ_{e: t∈T_e}  c_e · θ_t · (w[t,e] - μ_e) / μ_e
           - (φ_t - ν_t) / σ_t²
```

The first term is the score of the RFEQ log-likelihood in φ-space (one pass over
EQCs, same cost as one EM iteration). The second term is the Gaussian prior score
(O(T), negligible).

**Implementation note**: compute `μ_e = Σ_t θ_t · w[t,e]` once per EQC and reuse
for both the gradient and Fisher diagonal.

### Two-Phase Per-Sample Inference

```
Phase 1a: Run standard EM for ~20 iterations (warm-start).
           No prior needed here — EM gets into the right basin cheaply.

Phase 1b: Switch to L-BFGS on F(φ) using EM solution as initial point.
           Convergence criterion: ||∇F|| < ε  (e.g. ε = 1e-6 · T)
           Typical iterations to convergence from warm-start: 10–30.

Phase 1c: At L-BFGS convergence, compute diagonal Fisher information
           (one additional pass over EQCs — see below).
```

---

## Analytic Fisher Information (Diagonal Approximation)

After convergence to `φ̂` (equivalently `θ̂ = softmax(φ̂)`), compute:

### In θ-space (direct from RFEQ structure):

```
I_tt = Σ_{e: t∈T_e}  c_e · w[t,e]² / μ_e²
```

### Transform to φ-space via delta method:

```
I_tt^(φ) = θ̂_t² · (1 - θ̂_t)²  · I_tt^(θ)
           ≈ θ̂_t²  · I_tt^(θ)      for θ̂_t << 1  (most transcripts)
```

Use the exact form `θ_t(1-θ_t)` for the Jacobian factor; the approximation only
matters for very highly expressed transcripts.

### Laplace Posterior Variance:

```
Σ̂_t = 1 / ( 1/σ_t²  +  I_tt^(φ) )
```

This is the per-sample quantification uncertainty that gets propagated to the
outer hierarchical M-step.

### Edge Cases to Handle:

- `θ̂_t = 0` (unobserved transcript): `I_tt = 0`, `Σ̂_t = σ_t²` (prior dominates).
  Use a small floor `θ_min = 1e-300` before computing logs.
- `μ_e ≈ 0`: skip EQC or add numerical guard (`if mu_e < 1e-300 { continue }`).
- `I_tt^(φ) = 0` (unidentifiable): `Σ̂_t = σ_t²`. Correct — prior is all we know.

---

## Outer Hierarchical Loop (Empirical Bayes M-step)

This is the cross-sample information sharing. After Phase 1 produces
`{(φ̂_t^(s), Σ̂_t^(s))}` for all samples `s`:

### Update condition-level means:

```
ν̂_t^(c) = [ Σ_{s: c(s)=c}  φ̂_t^(s) / (σ_t² + Σ̂_t^(s)) ]
           / [ Σ_{s: c(s)=c}  1 / (σ_t² + Σ̂_t^(s)) ]
```

### Update biological variance (method-of-moments, closed form):

```
σ̂_t² = max(0,  (1/S) Σ_s (φ̂_t^(s) - ν̂_t^(c(s)))²  -  (1/S) Σ_s Σ̂_t^(s) )
```

The `max(0, ·)` clamp is critical: if sample-to-sample variation is dominated by
quantification noise, the biological variance goes to zero, and the prior loses
strength. This is the self-regulating property — regularization strength is
data-driven, not fixed.

### Shrink σ̂_t² toward a global trend (recommended):

Fit a smooth function `σ²(μ)` of mean expression level (analogous to DESeq2
dispersion fitting). Use a moderated estimator:

```
σ̂_t²_moderated = w_t · σ̂_t²_local  +  (1 - w_t) · σ²(μ̂_t)
```

where `w_t` is a weight based on the reliability of the local estimate (e.g.,
effective sample size for transcript `t`). This prevents noisy estimates for
low-coverage transcripts from destabilizing the prior.

### Convergence:

Alternate Phase 1 and the M-step. In practice 5–10 outer iterations suffice.
Warm-start inner L-BFGS from the previous solution — later outer iterations
converge in very few L-BFGS steps.

---

## Data Structures to Add / Modify

### New struct: `SamplePosterior`

```rust
pub struct SamplePosterior {
    /// MAP estimate in log-space (softmax parameterization)
    pub phi_hat: Vec<f64>,
    /// Diagonal Laplace posterior variance per transcript
    pub sigma_hat_sq: Vec<f64>,
}
```

### New struct: `HierarchicalHyperparams`

```rust
pub struct HierarchicalHyperparams {
    /// Experiment-level mean per transcript (log-space)
    pub mu: Vec<f64>,
    /// Condition-level means per transcript per condition
    pub nu: Vec<Vec<f64>>,   // nu[condition_idx][transcript_idx]
    /// Biological variance per transcript (shared across conditions)
    pub sigma_sq: Vec<f64>,
}
```

### Modify existing EM loop:

Add a `prior: Option<&HierarchicalHyperparams>` parameter. When `None`, run
standard EM as today. When `Some(h)`, run EM warm-start followed by L-BFGS on
`F(φ)` using `h.nu[condition]` and `h.sigma_sq`.

---

## Downstream DE Validity

**The double-shrinkage problem**: if you hand the posterior means `φ̂_t^(s)` to
DESeq2/edgeR as if they were observed counts, the hierarchical shrinkage has
already reduced apparent inter-condition variation. A naive DE test will be
anti-conservative.

**The solution**: the total variance seen by the DE model must be:

```
Var_total = σ_t²  +  Σ̂_t^(s)
```

not just `σ_t²`. Concretely:

1. **Preferred**: Compute DE within the hierarchical model itself. The posterior
   over `ν_t^(c1) - ν_t^(c2)` has variance `Var(ν^(c1)) + Var(ν^(c2))` where
   each condition variance accounts for both `σ_t²` and `Σ̂_t^(s)`.

2. **Compatibility mode (Swish)**: Generate bootstrap resamples of the RFEQ
   counts, re-run the penalized EM on each resample, pass the resulting
   distribution of `φ̂^(b)` to Swish. This correctly propagates both biological
   and quantification uncertainty into the nonparametric test.

3. **Compatibility mode (Sleuth-style)**: Pass `φ̂_t^(s)` as the "observed"
   estimate and `Σ̂_t^(s)` as the measurement error variance to a regression model.
   The model must include `Σ̂_t^(s)` in its variance structure.

---

## L-BFGS Implementation Notes

- **Do not implement L-BFGS from scratch.** Use the `argmin` crate
  (https://argmin-rs.org/) which has a production-quality L-BFGS implementation
  and a clean Rust API. Alternatively, `lbfgs` crate if you prefer a smaller
  dependency.

- **History size**: m=5 to m=10 is standard and sufficient. Memory cost is
  `O(m · T)` per sample, which is negligible.

- **Convergence criterion**: gradient norm `||∇F||₂ < 1e-6 · sqrt(T)` is a
  reasonable relative criterion. Also cap at 200 L-BFGS iterations.

- **Line search**: L-BFGS requires a Wolfe-condition line search. Both `argmin`
  and `lbfgs` provide this. Do not implement your own.

- **Numerical gradient check**: Before integrating, validate the analytic gradient
  against finite differences on a small synthetic EQC set (T=10, |E|=20).
  This is essential — the softmax Jacobian is easy to get wrong.

---

## Files Expected to Change

Based on the repo structure (`src/` directory, Rust codebase):

| File | Change |
|------|--------|
| `src/main.rs` or CLI entry | Add `--hierarchical` flag, `--num-outer-iter N` |
| `src/em.rs` (or equivalent) | Add `compute_fisher_diag()`, modify EM to accept prior params |
| `src/lib.rs` or new `src/hierarchical.rs` | `SamplePosterior`, `HierarchicalHyperparams`, outer M-step logic |
| `src/optimizer.rs` (new) | L-BFGS wrapper, penalized objective `F(φ)` and gradient |
| `Cargo.toml` | Add `argmin` (or `lbfgs`) dependency |

---

## Key Mathematical References

The three papers in the project knowledge base are directly relevant:

- **rfeq.pdf** (Zakeri et al. 2017, Bioinformatics): The RFEQ likelihood model
  that piscem-infer is based on. Section 2 defines the likelihood and factorization
  — the Fisher information derivation follows directly from Equation 1 there.

- **isolator.pdf** (Jones et al. 2016): The hierarchical model design (three-level:
  experiment / condition / sample) and the Gibbs sampling approach. We adopt the
  same hierarchical structure but replace MCMC with the empirical Bayes EM +
  L-BFGS approach for scalability.

- **polee.pdf** (Jones & Ruzzo 2021): The Pólya tree transformation for likelihood
  approximation. Directly relevant to the question of whether the diagonal Laplace
  approximation is sufficient (it may not be for highly ambiguous loci — Polee's
  approach is the more principled alternative if Laplace proves inadequate).

---

## Open Questions / Deferred Decisions

1. **Off-diagonal curvature**: The diagonal Fisher approximation ignores transcript
   co-correlations within EQCs. For isoforms of the same gene that are highly
   co-ambiguous, the off-diagonal terms `I_tt'` may be substantial. A block-diagonal
   approximation (one block per gene) may be warranted. Defer until accuracy
   evaluation on real data.

2. **ILR transform**: Working in ILR coordinates that respect the gene-isoform
   tree structure might give better-conditioned priors. Deferred.

3. **Pólya tree vs. Laplace**: If diagonal Laplace posteriors prove inadequate
   (check via comparison with Polee outputs on the same data), consider
   implementing the Pólya tree approximation from polee.pdf as the per-sample
   posterior representation feeding into the outer M-step.

4. **Scalability**: With S > 500 samples, the outer M-step matrix
   `φ̂[T × S]` may require care in memory layout. Rust's ndarray crate or a
   simple transposed Vec<Vec<f64>> should be evaluated.

---

## Implementation Order (Suggested)

1. **Gradient + Fisher** (no prior yet): Add `compute_gradient_phi()` and
   `compute_fisher_diag()` functions. Validate gradient numerically. This is
   self-contained and testable without the outer loop.

2. **L-BFGS wrapper**: Wire up `argmin` with the gradient function. Test on a
   single sample with a flat prior (`σ_t² = ∞`) — should converge to the same
   answer as the current EM.

3. **Single-sample penalized MAP**: Add prior parameters, test that increasing
   prior strength (`σ_t²` → 0) pulls estimates toward `ν_t`.

4. **Outer M-step**: Implement variance estimation across a small synthetic
   multi-sample dataset (3 conditions × 3 replicates, simulated from known
   abundances).

5. **Full integration**: CLI flags, file I/O for `SamplePosterior` outputs,
   convergence diagnostics.

6. **Evaluation**: Compare against per-sample EM on simulated data with known
   ground truth. Check calibration of `Σ̂_t` against bootstrap variance.
