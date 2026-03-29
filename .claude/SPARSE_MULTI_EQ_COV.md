# Sparse Multi-Sample Transcript Inference over Equivalence Classes

## Overview

This document describes a proposed optimization framework for **joint
transcript abundance inference across multiple RNA‑seq samples** using
**equivalence classes** and **group sparsity**. The goal is to improve
transcript identifiability and reduce false transcript usage by
exploiting the fact that the set of expressed transcripts is typically
**sparse and largely shared across related samples**.

The formulation is designed to integrate naturally with Salmon-style
pipelines that already compute equivalence classes, while remaining
computationally tractable for tens to hundreds of samples.

------------------------------------------------------------------------

# Motivation

Transcript quantification from RNA‑seq data is fundamentally an
underdetermined inverse problem. Reads often map compatibly to multiple
transcripts, leading to ambiguous assignment.

Current quantifiers generally solve the problem **independently for each
sample**. However, biological datasets typically exhibit strong
structure:

-   Only a small fraction of transcripts are expressed.
-   Replicates and related samples share most expressed transcripts.
-   Cross-sample patterns provide information that can help resolve
    ambiguous isoforms.

Ignoring this shared structure leads to solutions where abundance mass
spreads across many transcripts that are not truly expressed.

The proposed model introduces **shared sparsity across samples**,
allowing the algorithm to jointly determine:

1.  Which transcripts are expressed.
2.  Their abundances in each sample.

------------------------------------------------------------------------

# Data Representation

## Samples

Let:

S = number of samples

## Transcripts

Let:

T = number of transcripts

## Equivalence Classes

Reads are grouped into equivalence classes based on transcript
compatibility.

For each equivalence class c we observe:

n_cs

= fragment count for class c in sample s.

Compatibility weights:

w_ct

represent the probability that transcript t generated reads in class c.

------------------------------------------------------------------------

# Abundance Matrix

Define the abundance matrix:

Θ_ts

where:

Θ_ts ≥ 0

represents the abundance of transcript t in sample s.

Expected equivalence class counts:

λ_cs = Σ_t w_ct Θ_ts

------------------------------------------------------------------------

# Likelihood

Assume Poisson sampling:

n_cs \~ Poisson(λ_cs)

Negative log-likelihood:

L(Θ) = Σ_c Σ_s ( λ_cs − n_cs log λ_cs )

Constants independent of Θ are omitted.

------------------------------------------------------------------------

# Shared Transcript Activity

The core assumption is that the set of expressed transcripts is sparse
and largely shared across samples.

Introduce conceptual activity variables:

z_t ∈ {0,1}

where:

z_t = 1 → transcript potentially expressed\
z_t = 0 → transcript inactive across all samples

Transcript abundances factorize as:

Θ_ts = z_t η_ts

However, discrete indicators make optimization difficult.

------------------------------------------------------------------------

# Convex Relaxation via Group Sparsity

Instead of binary indicators, impose a **group sparsity penalty** on
rows of Θ.

Define the row vector:

Θ_t\* = (Θ_t1, Θ_t2, ..., Θ_tS)

Add the penalty:

Σ_t \|\|Θ_t\*\|\|₂

This is the **group LASSO** penalty.

Interpretation:

-   If transcript t is inactive, the entire row Θ_t\* shrinks to zero.
-   If active, abundances across samples remain nonzero.

------------------------------------------------------------------------

# Optimization Objective

The final convex optimization problem is:

minimize over Θ ≥ 0

L(Θ) + λ Σ_t \|\|Θ_t\*\|\|₂

Expanded:

minimize

Σ_c Σ_s ( λ_cs − n_cs log λ_cs ) + λ Σ_t \|\|Θ_t\*\|\|₂

subject to:

Θ_ts ≥ 0

where:

λ_cs = Σ_t w_ct Θ_ts

λ controls the strength of sparsity.

------------------------------------------------------------------------

# Interpretation

The model simultaneously:

-   explains fragment counts via equivalence classes
-   identifies a **shared sparse set of expressed transcripts**
-   estimates transcript abundances across samples

Transcripts that only weakly explain ambiguous reads across samples are
eliminated automatically.

------------------------------------------------------------------------

# Relationship to Single-Sample Quantification

Traditional quantifiers solve independent optimization problems:

minimize L_s(θ_s)

for each sample s.

The proposed model instead solves a **single coupled problem**:

minimize L(Θ)

with shared sparsity across samples.

This effectively converts many small inverse problems into **one
structured sparse inverse problem**.

------------------------------------------------------------------------

# Benefits

This approach can improve:

-   isoform disambiguation
-   detection of lowly expressed transcripts
-   robustness to multi-mapping ambiguity
-   stability across replicates

It also reduces the number of transcripts used to explain reads.

------------------------------------------------------------------------

# Computational Considerations

The formulation is compatible with equivalence class compression.

Per-iteration cost is approximately:

O(#equivalence_classes × #samples)

Optimization can be performed using:

-   proximal gradient descent
-   coordinate descent
-   FISTA-style accelerated methods

The proximal operator for the group L2 penalty is well known and
efficient.

------------------------------------------------------------------------

# Implementation Outline

## Step 1 --- Input

Load:

-   equivalence classes
-   compatibility weights w_ct
-   fragment counts n_cs for each sample

## Step 2 --- Initialize Abundances

Initialize Θ_ts using independent single-sample quantification.

## Step 3 --- Iterative Optimization

Repeat until convergence:

1.  Compute expected class counts λ_cs
2.  Evaluate gradient of likelihood
3.  Apply group-sparsity proximal update
4.  Enforce non-negativity constraint

## Step 4 --- Output

Return transcript abundance estimates Θ_ts.

------------------------------------------------------------------------

# Future Extensions

The framework can naturally incorporate additional constraints,
including:

-   coverage smoothness penalties
-   splice graph regularization
-   multi-condition hierarchical priors
-   single-cell group inference

These can be added as additional regularization terms.

------------------------------------------------------------------------

# Summary

Sparse multi-sample inference over equivalence classes reframes
transcript quantification as a **joint sparse optimization problem
across samples**.

The method leverages:

-   equivalence class compression
-   shared transcript activity
-   convex group sparsity

to improve transcript identification and abundance estimation while
remaining scalable for large RNA‑seq datasets.
