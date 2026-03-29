# TXPSELECTION: Sparse Transcript Variable Selection for RNA‑seq Quantification

## Overview

This document summarizes a proposed research and engineering direction for improving RNA‑seq transcript quantification through **principled variable selection**. The focus is identifying which transcripts are truly expressed before or during abundance estimation.

The central observation is that the transcriptome is extremely sparse relative to the full annotation (typically a few thousand expressed transcripts among hundreds of thousands). Additionally, RNA‑seq likelihood models are highly non‑identifiable because many transcripts share reads.

Therefore the goal is to develop algorithms that:

- Identify the **minimal identifiable transcript set**
- Remove structurally redundant transcripts
- Improve abundance estimation stability
- Enable stronger **multi‑sample inference**

The ideas described here leverage the **equivalence‑class (EC) representation** already used in modern quantifiers (alevin / salmon / piscem).

---

# Core Model

Let:

- T = number of transcripts
- K = number of equivalence classes
- y_k = count of reads/fragments in EC k
- A_{kt} = compatibility probability between EC k and transcript t

Expected counts:

μ = Aθ

where

θ ∈ R⁺^T

is the transcript abundance vector.

Likelihood (Poisson form):

ℓ(θ) = Σ_k y_k log( Σ_t A_{kt} θ_t ) − Σ_t θ_t

The variable selection problem is identifying

S = { t | θ_t > 0 }

with |S| ≪ T.

---

# Key Observations

1. **Likelihood depends only on Aθ**.
2. Many transcripts have identical EC compatibility patterns.
3. The matrix A is extremely rank‑deficient.
4. Many transcripts are structurally redundant.

Thus inference should operate on a **reduced basis of transcripts** rather than the full annotation.

---

# EC Bipartite Graph

Equivalence classes define a bipartite graph:

EC nodes ↔ transcript nodes

Edge (k, t) exists if transcript t is compatible with EC k.

The EC neighborhood of a transcript defines its **signature**:

S_t = { k | transcript t compatible with EC k }

These signatures fully determine the columns of A.

---

# Single‑Sample Variable Selection

## Stage 1 — Signature Collapsing

Transcripts with identical EC signatures are indistinguishable.

Algorithm:

1. For each transcript compute sorted EC list
2. Hash the list
3. Group transcripts with identical signatures
4. Keep one representative

Typical reduction:

200k transcripts → 5k–15k groups

---

## Stage 2 — Unique EC Peeling

Define

unique_ec(t) = { k ∈ S_t | deg(k) = 1 }

If a transcript contains an EC that no other transcript contains, it must be present.

Algorithm (peeling):

repeat

    find EC nodes with degree 1

    mark their transcript as required

    remove transcript and its edges

until no EC degree = 1

This identifies **structurally required transcripts**.

---

## Stage 3 — Subset Dominance Test

Transcript t is redundant if

S_t ⊆ S_u

for some transcript u.

This means transcript u can always explain the EC support of t.

Algorithm:

- sort transcripts by signature size
- test subset relations using bitsets or inverted index

---

## Stage 4 — Compatibility Cone Decomposition

Even after subset pruning, some transcripts remain redundant.

We interpret transcript columns as vectors generating a **nonnegative cone**.

Transcript t is redundant if:

A_t = Σ_u w_u A_u , w_u ≥ 0

where u ≠ t.

Only transcripts that generate **extreme rays** of this cone are identifiable.

### Local NNLS Test

Only transcripts sharing ECs with t can represent it.

Define candidate set:

N(t) = ⋃_{k ∈ S_t} neighbors(k)

Solve NNLS:

A_t = A_{N(t)\t} w

If feasible → transcript redundant.

These problems are very small (typically < 20 variables).

---

## Result

Remaining transcripts form the **minimal identifiable transcript basis**.

Expected reduction:

200k transcripts
↓
~10k signature groups
↓
~6k after structural pruning
↓
~3k–4k extreme transcripts

---

# Sparse Likelihood Refinement

After basis computation, perform statistical inference.

Possible approaches:

### Fisher‑information pruning

Compute importance score

I_t = Σ_k y_k A_{kt}² / ( Σ_j A_{kj} θ_j )²

Remove transcripts with low information.

---

### Hierarchical shrinkage

Use automatic relevance determination:

θ_t ~ Gamma(α, λ_t)

λ_t ~ Gamma(a, b)

Large λ_t shrinks transcripts toward zero.

---

### Final EM estimation

Run EM or variational inference on reduced transcript set.

---

# Multi‑Sample Variable Selection

Joint inference across samples can greatly improve support recovery.

Let:

θ_{t,s} = abundance of transcript t in sample s

Counts:

y_{k,s} ~ Poisson( Σ_t A_{kt} θ_{t,s} )

---

## Shared Support Model

Introduce transcript inclusion variable:

z_t ∈ {0,1}

z_t ~ Bernoulli(π)

θ_{t,s} | z_t =

0 if z_t = 0

Gamma(α, β) if z_t = 1

This encourages transcripts to be either globally present or absent.

---

## Multi‑Sample Information Aggregation

Compute Fisher information per sample:

I_{t,s}

Aggregate:

I_t = Σ_s I_{t,s}

Prune transcripts with low global information.

---

## Multi‑Sample EC Graph

Construct a joint EC graph across samples.

Edge weights correspond to total counts across samples.

Cone decomposition is run once globally.

All samples then estimate abundances over the same transcript basis.

---

# Final Pipeline

Single‑sample pipeline:

1. Build equivalence classes
2. Signature collapsing
3. Unique EC peeling
4. Subset dominance pruning
5. Cone decomposition (NNLS tests)
6. Sparse inference / EM


Multi‑sample pipeline:

1. Build ECs for all samples
2. Merge EC graph
3. Compute global transcript basis
4. Run multi‑sample inference
5. Apply hierarchical sparsity prior

---

# Computational Complexity

Let

E = number of EC‑transcript edges

Then:

Signature collapsing: O(E)

Peeling: O(E)

Subset pruning: O(T log T)

NNLS tests: O(T × small constant)

Overall complexity is essentially **linear in EC graph size**.

---

# Expected Benefits

1. Major reduction in optimization dimension
2. Removal of structurally unidentifiable transcripts
3. Improved abundance estimation stability
4. Better isoform detection accuracy
5. Natural support for multi‑sample inference

---

# Potential Extensions

### Integration with molecule resolution

Combine transcript selection with molecule resolution frameworks such as the **monochromatic arborescence cover** used for UMI resolution.

### Compressed sensing formulations

Treat transcript inference as sparse recovery problem over compatibility matrix.

### Identifiability diagnostics

Expose transcript basis structure to users for interpretability.

---

# Implementation Goals for Coding Agent

1. Implement EC signature hashing
2. Implement peeling algorithm
3. Implement subset dominance tests
4. Implement local NNLS redundancy tests
5. Produce transcript basis
6. Integrate with existing quantification pipeline
7. Add multi‑sample support

---

# Summary

The central idea is to exploit the **geometry of the compatibility matrix** implied by equivalence classes to compute a **minimal identifiable transcript basis** before statistical inference.

This approach connects RNA‑seq quantification with concepts from:

- convex geometry
- sparse inference
- compressed sensing
- mixture identifiability

and has the potential to significantly improve both **single‑sample and multi‑sample transcript quantification accuracy**.

