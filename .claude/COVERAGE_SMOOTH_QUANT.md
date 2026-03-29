# Coverage-Smooth Transcript Abundance Estimation

## Variance-Normalized Smoothness with Equivalence Classes

This document describes a proposed extension to Salmon-style transcript
abundance inference that incorporates **coverage smoothness
constraints** directly into the optimization objective. The method is
designed to discourage solutions that produce unrealistic coverage
spikes along transcripts while remaining statistically consistent with
RNA-seq sampling noise.

The formulation integrates:

-   equivalence class likelihood (as used in Salmon)
-   bin-level expected coverage modeling
-   **variance-normalized smoothness penalties**
-   optional **sparse transcript selection**

The approach is designed to remain **computationally scalable** and
compatible with existing equivalence-class pipelines.

------------------------------------------------------------------------

# 1. Overview

Most RNA-seq transcript quantification models assume that reads are
generated independently from transcripts. As a result, the likelihood
function does not distinguish between solutions that produce:

-   smooth coverage profiles
-   highly spiky coverage patterns

However, RNA fragmentation and sequencing typically produce **locally
smooth coverage distributions** after correcting for known biases.
Introducing a smoothness constraint allows the model to prefer
biologically plausible solutions.

To do this properly, the smoothness penalty must account for the fact
that **coverage variance scales with expected coverage**. A naive
penalty on absolute differences would incorrectly penalize
high-expression transcripts more strongly.

Instead, the smoothness penalty is **normalized by the expected sampling
variance**, yielding a statistically principled objective.

------------------------------------------------------------------------

# 2. Data Representation

## 2.1 Transcripts

t ∈ {{1,...,T}}

indexes transcripts.

------------------------------------------------------------------------

## 2.2 Transcript Binning

Each transcript is partitioned into bins of fixed width.

Typical choice:

bin_size ≈ 50--100 bp

Let

i ∈ {{1,...,B_t}}

index bins for transcript t.

Typical transcripts therefore have roughly:

10--40 bins

This binning allows spatial coverage information to be modeled without
storing per-base coverage.

------------------------------------------------------------------------

## 2.3 Equivalence Classes

Reads are grouped into equivalence classes using lightweight mapping.

For each equivalence class c we store:

n_c\
= number of fragments in the class

Compatibility weights:

w_ct

which represent the contribution of transcript t to equivalence class c.

Expected class counts:

λ_c = Σ_t w_ct θ_t

------------------------------------------------------------------------

# 3. Transcript Abundance Parameters

Each transcript has a non-negative abundance parameter:

θ_t ≥ 0

These parameters determine the expected contribution of transcripts to
fragment generation.

------------------------------------------------------------------------

# 4. Bias-Corrected Bin Coverage

For each transcript bin we compute a bias-corrected effective weight:

b_ti

This represents the expected probability that a fragment originates from
bin i of transcript t.

These weights incorporate standard Salmon bias corrections, including:

-   fragment length distribution
-   positional bias
-   sequence bias
-   GC bias

Expected bin coverage is therefore:

λ_ti = θ_t b_ti

------------------------------------------------------------------------

# 5. Likelihood (Equivalence Class Formulation)

Assume Poisson sampling for equivalence class counts:

n_c \~ Poisson(λ_c)

Negative log likelihood:

L(θ) = Σ_c ( λ_c − n_c log λ_c )

where

λ_c = Σ_t w_ct θ_t

Constants independent of θ are omitted.

------------------------------------------------------------------------

# 6. Coverage Smoothness Prior

Coverage changes between neighboring bins should typically be gradual.

However, observed differences increase naturally with expression level
due to sampling variance.

Under a Poisson model:

Var(y_i) = λ_i

For neighboring bins:

Var(y\_{{i+1}} − y_i) = λ\_{{i+1}} + λ_i

Therefore the correct normalization for differences is:

(y\_{{i+1}} − y_i) / sqrt(λ\_{{i+1}} + λ_i)

This ensures that large absolute changes are tolerated when coverage is
high but penalized when coverage is low.

------------------------------------------------------------------------

# 7. Variance-Normalized Smoothness Penalty

Define the smoothness penalty:

S(θ) = Σ_t Σ\_{{i=1}}\^{{B_t−1}} (λ\_{{t,i+1}} − λ\_{{t,i}})² /
(λ\_{{t,i+1}} + λ\_{{t,i}} + ε)

where ε is a small constant for numerical stability.

Substituting expected coverage:

λ_ti = θ_t b_ti

This penalty discourages abrupt changes in expected coverage that are
inconsistent with the variance of the sampling process.

------------------------------------------------------------------------

# 8. Optimization Objectives

## 8.1 Coverage-Smooth Abundance Estimation

The basic optimization problem is:

minimize over θ ≥ 0

L(θ) + α S(θ)

Expanded:

minimize over θ ≥ 0

Σ_c (λ_c − n_c log λ_c)

-   α Σ_t Σ_i (λ\_{{t,i+1}} − λ\_{{t,i}})² / (λ\_{{t,i+1}} +
    λ\_{{t,i}} + ε)

The parameter α controls the strength of the smoothness constraint.

------------------------------------------------------------------------

## 8.2 Coverage-Smooth Abundance with Sparse Transcript Selection

Sparse transcript usage can be encouraged using an L1 penalty.

Add:

β Σ_t θ_t

The full objective becomes:

minimize over θ ≥ 0

L(θ) + α S(θ) + β Σ_t θ_t

This encourages a smaller number of transcripts to explain the observed
reads.

------------------------------------------------------------------------

# 9. Interpretation of the Objective

The three components of the objective serve complementary purposes.

Likelihood term\
Explains observed fragment counts.

Smoothness term\
Discourages unrealistic coverage spikes.

L1 penalty\
Encourages sparse transcript usage.

Together these help resolve isoform ambiguity when multiple transcripts
share sequence.

------------------------------------------------------------------------

# 10. Efficient Implementation Strategy

## Step 1 --- Preprocessing

Partition transcripts into bins.

Compute bias-corrected bin weights:

b_ti

------------------------------------------------------------------------

## Step 2 --- Lightweight Mapping

Construct equivalence classes and compatibility weights:

n_c\
w_ct

using the existing mapping pipeline.

------------------------------------------------------------------------

## Step 3 --- Bin Data Structures

For each transcript store the vector:

b_t = \[b_t1, b_t2, ..., b_tB\]

These allow fast computation of expected bin coverage.

------------------------------------------------------------------------

## Step 4 --- Optimization

Optimize θ using:

-   projected gradient descent
-   proximal gradient methods
-   coordinate descent
-   L-BFGS with projection

Each iteration requires:

-   computing equivalence class expectations λ_c
-   computing expected bin coverage λ_ti
-   evaluating the smoothness penalty and gradient

------------------------------------------------------------------------

# 11. Computational Complexity

Per iteration cost is approximately proportional to:

O(number_of_equivalence_classes + number_of_transcript_bins)

Because the number of bins per transcript is small, this remains close
to Salmon's existing complexity.

------------------------------------------------------------------------

# 12. Practical Parameter Choices

Suggested starting values:

bin_size = 75 bp\
α ≈ 0.1 -- 10\
β ≈ 0 -- 0.01\
ε ≈ 1e-6

Hyperparameters can be tuned using simulation or held-out likelihood.

------------------------------------------------------------------------

# 13. Advantages of the Approach

Compared with standard transcript quantification models, this approach:

-   incorporates **spatial coverage structure**
-   penalizes **coverage spikes caused by multi-mapping**
-   remains compatible with **equivalence class compression**
-   preserves **convex optimization properties**
-   adds minimal additional computational overhead

------------------------------------------------------------------------

# 14. Expected Benefits

This approach should improve:

-   isoform disambiguation
-   robustness in highly overlapping transcript models
-   resistance to pathological EM solutions
-   stability of abundance estimates

while maintaining near-Salmon computational efficiency.

------------------------------------------------------------------------

# 15. Potential Extensions

Possible future directions include:

-   multi-sample joint inference
-   splice-graph Laplacian smoothness constraints
-   spike-and-slab transcript activity models
-   adaptive binning based on coverage density

------------------------------------------------------------------------

# 16. Implementation Goal

The goal of the first implementation prototype is to:

1.  load equivalence classes and transcript models\
2.  construct transcript bins\
3.  compute bias-corrected bin weights\
4.  optimize the objective described above\
5.  output transcript abundance estimates

This document is intended to serve as the **specification for
implementing the algorithm in code**.
