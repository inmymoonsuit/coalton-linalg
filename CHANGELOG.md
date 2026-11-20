# Changelog

Changes a caller can meet. Each release's evidence is in the batch record it links.

## 0.7.2 — 2026-11-20

Changed: the fixed-pattern Gram analysis is faster and returns exactly the pattern it
returned before; refills and refusals are unchanged. On thirteen public systems
`analyze-ata` and `analyze-aat` take between 1.7 and 13 times less.
[Record](docs/batch8.md#addition-the-fixed-pattern-gram-analysis).

## 0.7.1 — 2026-11-20

Changed: sparse analysis is faster and returns exactly what it returned before. On
thirteen public patterns `analyze` takes between 1.7 and 25 times less; the permutation,
the factor structure and every refusal are unchanged for every input, and no call changes
meaning. [Record](docs/batch8.md).

## 0.7.0 — 2026-11-20

Added:

- `factor-signed!` in `coalton-linalg/sparse/spd-ldlt`: L D Lᵀ with pivots of either sign
  under the order fixed at analysis, for quasi-definite matrices, with `factor-growth`,
  an exact measure of how much accuracy the factorization kept, and
  `positive-pivot-count` and `negative-pivot-count`. **Read the growth before using the
  factor**: no fixed order is accurate on every matrix.
- `analyze-quasidefinite` and `analyze-quasidefinite-with-limit`: minimum-degree analysis
  that holds a trailing block of indices back, each released by an eliminated leading
  index.
- `coalton-linalg/sparse/augmented`: the matrix `[[top I, Aᵀ], [A, −bottom I]]` of a sparse
  A on a fixed pattern, analyzed once and refilled without allocating. With the two
  entries above it gives a regularized minimum-norm solve whose cost follows A's sparsity
  where A Aᵀ would be dense.
- A compiled example, `examples/augmented-minimum-norm.ct`.

[Contract](docs/quasidefinite.md) · [record](docs/batch7.md). No existing call changes
meaning. `factor!` and `factor-cholesky!` now also set the pivot counts and clear the
growth figure.

## 0.6.3 — 2026-11-20

Added: `factor-cholesky!` and `cholesky-ready?` in `coalton-linalg/sparse/spd-ldlt`, a
Cholesky-form numeric entry over the existing analysis, for proofs of positive
definiteness that are stated for Cholesky's arithmetic. It overwrites the Symbolic's one
factor and no solve accepts its result. [Contract](docs/cholesky-form.md) ·
[record](docs/batch6a.md). No existing call changes meaning.

## 0.6.2 — 2026-11-20

Fixed: **`eigen!` could return wrong eigenvalues as `JacobiSuccess`** for a matrix whose two
diagonal entries differ by more than the largest double, for example 10³⁰⁸ and −10³⁰⁸ with
an off-diagonal entry of 10³⁰⁷. The difference overflowed and the off-diagonal entry was
dropped without a rotation. The routine now scales the matrix by an exact power of two
first. Only a matrix with an entry above half the largest double, about 9 × 10³⁰⁷ in
magnitude, could be affected.
[Record](docs/batch6.md#after-review).

## 0.6.1 — 2026-11-20

Added: an artifact driver for the Jacobi routines, `scripts/check.ps1 -Task Jacobi`.
[Protocol](docs/jacobi.md#artifact-driver). No library call changed.

## 0.6.0 — 2026-11-19

Added:

- `coalton-linalg/dense/jacobi`: all eigenvalues and eigenvectors of a dense symmetric
  matrix by cyclic Jacobi rotations, and all singular values and vectors of a dense matrix
  of any shape by one-sided Jacobi rotations, with a complete right factor whose trailing
  columns are a nullspace basis, a rank count against the caller's threshold, and a sweep
  limit whose exhaustion is a named outcome. [Contract](docs/jacobi.md) ·
  [record](docs/batch6.md).
- A compiled example reading a rank from singular values, `examples/singular-rank.ct`.

Corrected: the orthogonal decomposition's contract now names a refusal it did not: a
column whose norm exceeds half the largest double. No code changed.

## 0.5.1 — 2026-11-19

Added: an artifact driver for the orthogonal decompositions, `scripts/check.ps1 -Task
Orthogonal`. [Protocol](docs/orthogonal.md#artifact-driver). No library call changed.

## 0.5.0 — 2026-11-19

Added:

- `coalton-linalg/dense/orthogonal`: column-pivoted Householder QR of a dense matrix of
  any shape, a rank chosen by the caller from the diagonal of R, the complete orthogonal
  decomposition for that rank, the minimum-norm least-squares solution with its residual
  norm, an orthonormal nullspace basis, and Q applied without being formed.
  [Contract](docs/orthogonal.md) · [record](docs/batch5.md).
- A compiled example with a repeated equation, `examples/dependent-constraints.ct`.

Corrected: the dense indefinite documents said a singular matrix completes with a nonzero
zero count. That holds when a computed pivot is exactly zero; an exactly singular input
can show a pivot at rounding level instead. No code changed.

## 0.4.1 — 2026-11-19

Added: `make-factor-with-limit` and `default-entry-limit` in
`coalton-linalg/dense/indefinite-ldlt`, and an artifact driver, `scripts/check.ps1 -Task
Indefinite`. [Record](docs/batch4.md#after-review).

Changed, and able to break a caller of 0.4.0: **`make-factor` refuses n above 4096**, where
n² exceeds the default budget of 16,777,216 entries, with `AnalysisStorageExceeded`. Pass a
larger budget to `make-factor-with-limit`.

## 0.4.0 — 2026-11-19

Added:

- `coalton-linalg/dense/indefinite-ldlt`: dense symmetric indefinite LDLᵀ with 1-by-1 and
  2-by-2 pivots, solves, and the inertia of the computed block diagonal, counting both
  eigenvalues of each 2-by-2 block. A singular matrix completes with a zero count and its
  solve is refused by name. [Contract](docs/indefinite-ldlt.md) · [record](docs/batch4.md).
- A compiled minimum-norm example through the augmented system, `examples/minimum-norm.ct`.

Changed: `AnalysisLimitExceeded` is now also raised by `make-factor` for a dense dimension
above `maximum-dimension`. No existing call changes meaning.

## 0.3.0 — 2026-11-19

Added:

- A condition estimate in `coalton-linalg/sparse/spd-ldlt`: `estimate-inverse-norm1!` gives a
  lower estimate of the 1-norm of the inverse of the factored, shifted matrix in at most 11
  solves without allocating; `shifted-norm1` gives the matching matrix norm and
  `reciprocal-condition1` combines them. Outcomes are a new `ConditionStatus`.
  [Contract](docs/condition.md) · [record](docs/batch3d.md).

Unchanged: every earlier result and artifact. A `factor!` call now also clears the stored
estimate; no existing call changes meaning.

## 0.2.0 — 2026-11-19

Added:

- `coalton-linalg/sparse/gram`: fixed-pattern AᵀA and AAᵀ. Analyze a checked CSC matrix once,
  then `refill!` caller-owned values without allocating. The output is the lower triangle with
  every diagonal, the pattern SPD analysis takes. [Contract](docs/gram.md) ·
  [record](docs/batch3c.md).
- `matvec!` and `transpose-matvec!` in `coalton-linalg/sparse/csc`, writing into caller storage.
  [Record](docs/batch3b.md).
- `analyze-with-limit` and `analyze-csc-with-limit` in `coalton-linalg/sparse/spd-ldlt`, taking
  an explicit budget of total factor entries. [Record](docs/batch3a.md).
- A compiled normal-equations example, `examples/normal-equations.ct`.

Changed, and able to break a caller:

- **SPD analysis no longer refuses dimensions above 4096.** It refuses instead when the factor
  would exceed an entry budget, 1,048,576 total factor entries by default. A matrix below 4096
  whose factor fills beyond that budget, which 0.1.0 analyzed, now needs an explicit larger budget.
- **`LinalgError` gained two constructors**, `InvalidAliasing` and `AnalysisStorageExceeded`.
  An exhaustive match on the exception family must add them.
- `maximum-analysis-size` is now the representability bound for dimensions and budgets, not 4096.

Unchanged: the elimination order, every numeric factor and solve result on the published test
rosters (byte-identical output files), and the positive-pivot contract.

## 0.1.0 — 2026-11-19

The first published snapshot: validated CSC and COO construction, transpose, permutation,
allocating matrix-vector products, stable norms, and sparse SPD LDLᵀ with reusable symbolic
analysis, typed factor and solve statuses, and typed construction errors.
