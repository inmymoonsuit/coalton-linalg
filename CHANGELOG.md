# Changelog

Changes a caller can meet. Each release's evidence is in the batch record it links.

## 0.19.0 — 2026-11-23

Changed, in `coalton-linalg/sparse/rank-deficient`:

- `refresh!` first looks for a vector's new multipliers on the rows the vector already
  holds, and solves on all the kept rows only where those do not serve. So it replaces
  such vectors where 0.18.0 returned `RefreshNoShift` and replaced none: **that status,
  and `RefreshFactorRefused`, now mean that a vector needed the kept rows' solve and it
  could not be had, the other vectors having been replaced.** After a refresh the gap
  bound is under the target as before, and can be nearer to it than 0.18.0 left it.
- A refresh no longer factors the system of the rows left out again where no vector
  changed its rows. The step solves it by conjugate gradients preconditioned with the
  factor it has, and makes the factor again only when a solve does not converge in
  eight iterations. A step taken with a factor of the vectors as they stand is
  unchanged, bit for bit.

Added: `refresh-to!`, a refresh with the fraction of its share a replaced vector is
corrected to, and `default-refresh-fraction`; `step-norm-term` and `step-outside-term`,
the two terms `step-bound` adds to `step-solve-bound`; `system-iterations`,
`system-factor-count` and `system-iteration-limit`.
[Contract](docs/rank-deficient.md#refreshing-the-kept-vectors) · [record](docs/batch18d.md).

## 0.18.0 — 2026-11-23

Added to `coalton-linalg/sparse/rank-deficient`: `refresh!(kept, matrix, gap-target)`,
which proves the kept vectors' bounds for new values as `reprove!` does and then replaces
each vector whose bound is above its share of the target by the row's multipliers for
those values; and `refreshed-count`. A redundant row whose values change at every step
no longer needs a search at every step. Six named outcomes.
[Contract](docs/rank-deficient.md#refreshing-the-kept-vectors) · [record](docs/batch18c.md).

No result of 0.17.0 changes.

## 0.17.0 — 2026-11-23

Added to `coalton-linalg/sparse/rank-deficient`: `step!(kept, r, out)`, the step on the
matrix a certificate stands on. It writes `−B⁺ r`, B being the matrix with each row left
out replaced by its kept combination of the kept rows, and leaves two proved bounds:
`step-solve-bound`, on the distance from the exact step for B, and `step-bound`, on the
distance from the step of the matrix truncated at the kept count. Seven named outcomes;
a refusal writes nothing. Also `step-taken?`, `step-corrections`, `correction-limit`, and
three diagnostics, `solve-growth`, `solve-backward-bound` and `predicted-rate`. A steady
loop of reproof, certificate and step allocates nothing.
[Contract](docs/rank-deficient.md#the-step) · [record](docs/batch18b.md).

Changed: `find-kept-rows!` also analyzes and factors the system of the rows left out, and
so can raise the typed errors of those analyses where 0.16.0 did not. No result of
0.16.0 changes.

## 0.16.0 — 2026-11-23

Added: `coalton-linalg/sparse/rank-deficient`, the kept rows of a rank-deficient matrix
for a loop that takes many steps on one pattern. `find-kept-rows!` searches once and
keeps the vector of each row left out, entries below a drop tolerance set to zero and
the bound proved for the vector as kept; `reprove!` proves every bound again for other
values with no search and no allocation; `certify!` proves, or refuses by name, that
the gap bound is below one floor and every singular value of the kept rows above
another; `largest-singular-bounds!` bounds the largest singular value of any stored
matrix from both sides. Nothing is solved by this package yet.
[Contract](docs/rank-deficient.md) · [record](docs/batch18a.md).

Added: `dependent-rows:check-pattern!`, the pattern check of a search by itself. No
result of an earlier package changes.

## 0.15.0 — 2026-11-22

Added: `coalton-linalg/sparse/damped`, damped least squares on the normal equations. One
analysis for a pattern, then `set-matrix!`, a damping by `set-damping!` or
`set-column-damping!`, and for each μ a `factor!` and a `step!` giving
`−(Jᵀ J + μ D)⁻¹ Jᵀ r`. Refusals are by name; a steady loop allocates nothing; the step at
any scale the format holds is the step at ordinary scale.
[Contract](docs/damped.md) · [record](docs/batch17.md).

## 0.14.0 — 2026-11-22

Changed, and a caller that matches `RowStatus` exhaustively must add a case:

- `RowStatus` gains `RowsUnderflow`: under a threshold of zero, a row holding a nonzero
  entry that is too small beside the matrix's largest entry for its square to be held.
  Version 0.13.4 reported this as `RowsNonFinite`, which again means only that an
  intermediate result of the factorization was not finite.
- `find-dependent-rows!` works on the matrix divided by the power of two at its largest
  entry, so it no longer reaches only half the exponent range. Entries of any size the
  format holds give the rows left out, the vectors and the values of the same matrix at
  ordinary size, bit for bit. A matrix of large entries no longer raises the Gram
  product's `NonFiniteResult`. `RowsBadThreshold` for a small positive threshold is now
  relative: the threshold is too small beside the largest entry. A threshold whose own
  square is not finite is no longer refused; it is a distance beyond every row, and
  leaves every row out.
- `transpose-product-bound` is formed the same way, and is finite wherever the norm it
  bounds can be held. Two products too large to hold that cancel now give a finite bound.
- `cg!` and `pcg!` carry the iteration on the residual divided by a power of two. A
  right-hand side from 2^510 no longer ends `CgNonFinite`, and an operator with entries
  far from 1 no longer ends `CgNonFinite` or, wrongly, `CgNotPositive`. Iterates on data
  that worked before are unchanged bit for bit. The operator and the preconditioner are
  now called with the divided vectors.

Added: `gram:check-csc!`, the source pattern check of `refill-csc!` by itself.
[Record](docs/batch16.md).

## 0.13.4 — 2026-11-22

Fixed: `find-dependent-rows!` on a matrix of entries near 10⁻¹⁸⁰ and below returned
`RowsFound` with every row left out, the squares it compares having underflowed to zero.
A positive threshold whose square is too small to be a normal number is now
`RowsBadThreshold`, and under a threshold of zero a row whose square is too small is
`RowsNonFinite`. A caller that passed a positive threshold below about 1.5 × 10⁻¹⁵⁴ and
took the answer now gets the refusal. [Contract](docs/dependent-rows.md#interface).

Added: two rosters in the test suite that every exported call must answer, one for two
arguments that are one array and one for the same problem at every scale, and the
measured range of each routine in its contract. No other result changes.
[The rosters](docs/rosters.md) · [record](docs/batch15.md).

## 0.13.3 — 2026-11-22

Changed: the pinned Coalton compiler moves from `b6e284b` to `deb300be`, which fixes the
miscompilation of a cell-holding function inlined from another package. No source of
the library changes and no result changes; `toolchain.lock.json` and
`toolchain.sources.json` describe the new compiler, and a build against the old one no
longer matches the lock. [Record](docs/batch14.md#the-pin-afterwards).

## 0.13.2 — 2026-11-22

Fixed, by an outside audit contributed as pull request 1 (Robert Smith), each with a
test that fails on the earlier sources. [Record](docs/batch14.md).

- `cg!`, `pcg!`, and `refine!` preserve the original right-hand side when it shares
  the iterate's array. CG no longer mistakes an underflowed square for convergence;
  PCG reports the final residual without asking for a step beyond the limit.
- Dense indefinite block solves no longer divide by an overflowed determinant scale.
  Jacobi SVD retains rank across the former square-underflow boundary, and QR forms
  orthogonal reflectors even for columns of subnormal entries.
- Sparse LU refuses nonfinite multipliers. Augmented analysis refuses an insufficient
  entry budget before copying the input pattern.

`make-cg` and `make-refiner` each allocate one more vector than before. On data away
from the ends of the floating-point range every result is the same, bit for bit. The
build now fails on a compiler warning in the library.

## 0.13.1 — 2026-11-20

Fixed, in `coalton-linalg/sparse/scaling`, after a held-back review:

- `scale!` and `scale-symmetric!` are now exact for powers of two whenever the result is
  zero or a normal number. Before, a value multiplied by one factor could leave the
  range on the way to a normal result and come back rounded: `[[1, 1e-310]]` scaled its
  second entry to a value wrong in its last digits.
- `equilibrate!` gives a column whose entries fall below the normal range under the row
  factors the largest factor, 2¹⁰²³. Before, such a column could keep the factor 1, so
  that `scale!` wrote zero for an entry it can keep exactly, or a factor one step short.
- The contract names the rows and columns that no factor brings to 1/2.

For a matrix whose entries span less than about 2¹⁰²² within every row, no factor and no
scaled value changes. With factors that are not powers of two, which the contract makes
no claim about, a scaled value may differ in its last bit from before.
[Contract](docs/scaling.md) · [record](docs/batch9a.md).

## 0.13.0 — 2026-11-20

Added: `coalton-linalg/conjugate-gradient`. `cg!` and `pcg!` solve a symmetric positive
definite system through products with its operator, plain or with a preconditioner the
caller supplies as a function, to an absolute tolerance on the 2-norm of the residual,
which is formed from the operator before convergence is reported.
[Contract](docs/conjugate-gradient.md) · [record](docs/batch13.md).

No existing call changes.

## 0.12.0 — 2026-11-20

Added: `coalton-linalg/sparse/lu`. `factor!` computes `P A Q = L U` for a square sparse
matrix with partial pivoting, the pivots chosen from the values and the diagonal
preferred within a threshold; `solve!` and `solve-transpose!` solve with the matrix and
with its transpose; `lu-growth` and `lu-backward-bound` report on the factors that
stand; `make-lu`, `make-lu-with-order` and `column-order` make the handle, whose storage
grows up to an entry budget. [Contract](docs/lu.md) · [record](docs/batch12.md).

No existing call changes.

## 0.11.0 — 2026-11-20

Added:

- `coalton-linalg/sparse/dependent-rows`: `find-dependent-rows!` gives the rows of a
  sparse matrix to keep and those that depend on them, `null-vector!` a left null vector
  for each row left out, `residual-bound` a bound on what that row leaves that holds for
  the exact product, and `gap-bound` a bound on the singular values beyond the kept
  count. `kept-values!` supplies what the existing rank certificate needs to prove the
  kept rows independent, and `transpose-product-bound` is the bound by itself.
  [Contract](docs/dependent-rows.md) · [record](docs/batch11.md).
- In `coalton-linalg/sparse/spd-ldlt`: `factor-skipping!` leaves out each index whose
  pivot is at or below a floor and factors the rest; `solve-skipping!` solves with the
  kept principal submatrix; `skipped-null-vector!` gives an index left out its
  combination of the kept ones; `skipping-ready?` and `skipped-count` report the state.
  [Contract](docs/skipping.md).

No existing call changes meaning. `factor!` and `factor-signed!` are unchanged,
operation for operation, and a factor from `factor-skipping!` is refused by `solve!`,
`solve-cholesky!` and `estimate-inverse-norm1!`.

## 0.10.0 — 2026-11-20

Added:

- `coalton-linalg/refinement`: `refine!` improves a solution through any approximate solve
  and an operator the caller applies, stopping on an absolute tolerance in the infinity
  norm, a caller's stall factor or a step limit, and returning the best iterate seen with
  its residual history. [Contract](docs/refinement.md) · [record](docs/batch10.md).
- `solve-cholesky!` in `coalton-linalg/sparse/spd-ldlt`: a solve with the Cholesky-form
  factor, of the shifted matrix that factor is of. [Contract](docs/cholesky-form.md).

No existing call changes meaning. `factor-cholesky!` is unchanged, operation for
operation.

## 0.9.1 — 2026-11-20

Changed: every exported symbol now carries a docstring, including each status value,
error constructor and exported constant, and the sources follow the Coalton developers'
published style rules. No call changes meaning and no result changes.

## 0.9.0 — 2026-11-20

Added: `coalton-linalg/sparse/scaling`. `equilibrate!` gives power-of-two row and column
factors that bring every row's and column's largest magnitude into [1/2, 1);
`equilibrate-symmetric!` scales a symmetric positive definite matrix by its diagonal;
`scale!`, `scale-symmetric!` and `scale-vector!` apply the factors, exactly; `row-maxima!`
and `column-maxima!` report the result. A compiled example, `examples/scaled-solve.ct`.
[Contract](docs/scaling.md) · [record](docs/batch9.md). No existing call changes meaning.

## 0.8.0 — 2026-11-20

Added: `factor-backward-bound` in `coalton-linalg/sparse/spd-ldlt`, an upper bound on the
norm of the signed factorization's backward error, with what underflow can add included,
and with it a documented certificate that every singular value of a sparse matrix exceeds
a threshold, from one signed factorization. [Contract](docs/quasidefinite.md#a-rank-certificate-from-the-inertia)
· [record](docs/batch7a.md).

Corrected: the signed factorization's entrywise error bound was stated without its
condition. It holds barring underflow; the contract now says so and gives the term.
No call changes meaning.

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
