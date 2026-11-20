# Changelog

Changes a caller can meet. Each release's evidence is in the batch record it links.

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
