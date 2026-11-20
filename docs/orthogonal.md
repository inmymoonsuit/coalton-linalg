# Orthogonal decompositions and minimum-norm least squares

Package: `coalton-linalg/dense/orthogonal`. Source: [orthogonal.ct](../src/dense/orthogonal.ct).
State: completed as Batch 5; verified at safety 1 and safety 0. The [batch record](batch5.md)
holds the measurements and the controls.

## What it computes

For a dense m-by-n matrix A of any shape, in three steps the caller drives:

1. `factor!` computes a column-pivoted Householder QR, `A P = Q R`.
2. The caller reads the diagonal of R, decides a rank r, and calls `set-rank!`. That keeps
   the leading r rows of R and reduces them from the right to a triangle:
   `A P ≈ Q [[T, 0], [0, 0]] Zᵀ`, a complete orthogonal decomposition of the matrix with
   the trailing block of R left out.
3. `solve!` returns the **minimum-norm least-squares solution of that rank-r problem**, and
   `nullspace-vector!` returns an orthonormal basis of its nullspace.

With r = n the solution is the ordinary least-squares solution of a full-column-rank
matrix; with r = m < n it is the minimum-norm solution of a full-row-rank system; with
r below both it is the minimum-norm solution among the least-squares minimizers of the
truncated matrix.

## Interface

A matrix is an F64 array of m n elements in column-major order: entry (i,j) is at
`i + j m`. Every entry is read.

- `make-decomposition(m, n) -> Decomposition` allocates reusable storage for one shape:
  about 2 m n doubles when m ≤ n, m n + n² when m > n, plus O(m + n). It refuses by name
  before allocating: `AnalysisStorageExceeded(m n, budget)` above `default-entry-limit`
  (16,777,216), `AnalysisLimitExceeded` for a dimension above `maximum-dimension`
  (32,768). `make-decomposition-with-limit(m, n, maximum-entries)` takes the budget from
  the caller. Either dimension may be zero.
- `factor!(decomposition, a) -> FactorStatus`. The input is copied, not retained.
- `pivot-magnitude(decomposition, k) -> F64` is |R[k,k]| for k below `step-count`, the
  smaller of m and n; it is 0 for any other k and when no factor stands.
- `rank-above(decomposition, threshold) -> UFix` is the number of leading diagonal
  entries of R **strictly larger in magnitude than threshold**. An entry equal to the
  threshold is not counted. The threshold is absolute; a caller wanting a relative one
  multiplies by `pivot-magnitude(decomposition, 0)`, the largest.
- `set-rank!(decomposition, rank) -> RankStatus` selects any rank from 0 to `step-count`.
  It does not need a new factorization, and a different rank may be selected afterward.
- `solve!(decomposition, b, x) -> SolveStatus`, with b of m entries and x of n. On any
  refusal x is unchanged. `solve-residual-norm` then gives the Euclidean norm of the
  residual of the latest successful solve, or -1 when the latest solve did not succeed.
- `nullspace-vector!(decomposition, index, out) -> SolveStatus` writes basis vector
  `index`, from 0 to n - rank - 1, of the selected nullspace into out of n entries.
- `apply-q!(decomposition, v)` and `apply-qt!(decomposition, v) -> SolveStatus` replace a
  vector of m entries by Q v or Qᵀ v, without forming Q. Q is m-by-m.
- `row-count`, `column-count`, `step-count`, `factor-ready?`, `rank-set?`,
  `selected-rank`, and the inspection views below.

| FactorStatus | Meaning |
| --- | --- |
| FactorSuccess | A factorization stands; no rank is selected yet |
| FactorNonFinite | A non-finite input, column norm, or computed entry of R |
| FactorBadSize | The array does not have m n elements |

| RankStatus | Meaning |
| --- | --- |
| RankSuccess | The decomposition is completed for this rank |
| RankInvalidFactor | No factorization stands |
| RankTooLarge | The rank exceeds `step-count` |
| RankSingular | A diagonal entry of R within the rank is exactly zero |
| RankNonFinite | Non-finite arithmetic while completing |

| SolveStatus | Meaning |
| --- | --- |
| SolveSuccess | The output holds the result |
| SolveInvalidFactor | No factorization stands |
| SolveNoRank | No rank is selected |
| SolveBadSize | A vector has the wrong length |
| SolveBadIndex | The nullspace index is not below n - rank |
| SolveNonFiniteRhs | Non-finite input vector |
| SolveNonFiniteResult | Non-finite computed result; an in-place application leaves its vector unspecified |

Every `factor!` call withdraws the selected rank, and a refused one withdraws the
factorization too. Every `set-rank!` call withdraws the earlier rank, whatever it returns.
One Decomposition owns its storage and scratch, so calls sharing it are not reentrant or
concurrent. Successful calls allocate nothing after warm-up on the measured runtime.

Inspection views are borrowed and read-only: `factor-permutation`, where `permutation[k]`
is the original column at pivoted position k; `factor-work`, m-by-n column-major, holding
R on and above the diagonal and each reflector's tail below it; `factor-tau`, one scalar a
reflector.

## What the rank means

**The caller decides the rank; the library reports what it was given to decide with.**
The diagonal of R is non-increasing in magnitude, and its entry k is the largest
Euclidean norm among the columns still remaining at step k. `rank-above` counts against
one explicit threshold with one explicit boundary rule, and nothing else is built in.

**Diagonal entries of R are not singular values.** The first is within a factor √n of
the largest singular value. The last can exceed the smallest singular value by a large
factor: on the Kahan matrix of order 30 the tests measure a last diagonal entry more than
three thousand times a proven upper bound on the smallest singular value. Column pivoting reveals rank reliably on most matrices and is
not guaranteed to on all. A caller whose decision turns on the smallest singular value
itself needs the singular value decomposition.

**A selected rank below the exact one solves a nearby problem.** The trailing block of R
is dropped, so the solution is the minimum-norm least-squares solution of a matrix that
differs from A P by that block, whose columns have norm at most the first dropped
diagonal entry. `solve-residual-norm` is the residual against that truncated matrix.

**A selected rank above the exact one is refused only when it meets an exact zero.**
Rounding usually leaves tiny nonzero diagonal entries where exact arithmetic would have
zeros; selecting them divides by them. Choosing the threshold is how a caller avoids that.

## Numerical assumptions and limits

Arithmetic is IEEE double with traps disabled. A sum of squares is formed unscaled while
its largest entry lies between 2⁻⁵⁰⁰ and 2⁵⁰⁰ and otherwise after an exact scaling by
2^±600, so a matrix scaled by a power of two gives the same pivoting and the same
reflectors bit for bit as long as its entries stay well inside the representable range. A column norm that overflows is a
refusal. No scaling or equilibration of rows or columns is applied, and the result
depends on column scaling, as pivoting by column norm does.

Work is about 2 m n² − 2n³/3 multiply-adds for the factorization when m ≥ n. The norms
of the remaining columns are measured in the same pass that updates each column, so
measuring adds no pass over the matrix. Completing for rank r costs about r²(n − r), and a
solve O(m n + n r). The code is unblocked and sequential.

## Artifact driver

The optional tools system writes reviewable outputs for a folder of systems:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/check.ps1 -Task Orthogonal `
  -InputFolder <folder> -OutputFolder build/<name> -AllocationRepeats 20
```

Each `<name>.M.txt` holds m and n, then m n values in column-major order. `<name>.b.txt`
holds m values. `<name>.meta.txt` is optional; each of its lines names one run, in order:
`threshold: <value>` selects the rank counted strictly above that value times the largest
diagonal magnitude of R, and `rank: <r>` selects exactly r. With no such line there is one
run at a threshold of 1e-10.

Once a system, `<name>.qr.out.txt` holds m, n, the permutation, the upper trapezoid of R as
k n values column-major with k the smaller dimension and its zeros written out, then the
first k columns of Q as m k values column-major, read back through `apply-q!`. Every value
is on its own line, floating values with 17 significant digits.

For run index i from 0: `<name>.<i>.x.out.txt` holds the n values of the solution;
`<name>.<i>.null.out.txt` the n - rank basis vectors of the selected nullspace, n values
each, one after another, and nothing at full column rank; `<name>.<i>.run.out.txt` the two
lines `rank <r>` and `residual <value>`. A run whose rank or solve is refused writes only
`<name>.<i>.refused.out.txt` holding the status name. A system refused at the factorization
writes only `<name>.refused.out.txt`.

`orthogonal-table.txt` and standard output carry one line a run:
`name run m n rank factor_ms solve_ms factor_bytes solve_bytes`, or
`name run m n refused <Status>`, with `factor` in place of the run index for a refused
factorization. Reruns replace every artifact of a system.

## Derivation

**A reflector.** For a vector with head α and tail of norm σ, the matrix
`H = I − τ u uᵀ` with u = (1, tail/(α − β)) maps it to (β, 0), where
`β = −sign(α) √(α² + σ²)` and `τ = (β − α)/β`. To check: `uᵀx = α + σ²/(α − β)`, and
with α² + σ² = β² the product τ uᵀx is α − β, so the head becomes β and each tail entry
`t − (α − β) t/(α − β) = 0`. H is symmetric and orthogonal because
`τ (uᵀu) = 2`. The sign of β is opposite to α's so that α − β adds magnitudes. When σ is 0
the vector is already of the wanted form and τ is 0. √(α² + σ²) is formed from the two
values divided by the larger.

**Pivoted QR.** At step k, the remaining column of largest norm, over rows k and below,
is exchanged into position k; ties take the first. A reflector built from that column's
entries k and below is applied to the columns after it. After min(m, n) steps the
reflectors' product is Q, with `Qᵀ = H_last ⋯ H₀`, and what remains on and above the
diagonal is R. Each remaining column's norm is measured from its current entries as they
are written, at every step. No formula updates an old norm, so nothing can drift.

**Completion.** With the leading r rows of R as `[R₁₁ R₁₂]`, work upward from row r − 1.
For row i, a reflector acting on column i and columns r to n − 1 folds that row's entries
beyond the triangle into its diagonal entry. Rows below i have zeros in all those
columns and are untouched; rows above are updated. After row 0 the result is `[T 0]`
with T upper triangular, and `[R₁₁ R₁₂] = [T 0] Zᵀ` with `Z = G₍r−1₎ ⋯ G₀` as the
reflections were applied.

**Solve.** Minimizing `‖A x − b‖` over the truncated matrix with `x = P Z y`: the
residual is `‖[[T, 0], [0, 0]] y − Qᵀb‖`. Its first r entries are made zero by
`y₁ = T⁻¹ (Qᵀb)₁`; the rest of Qᵀb is the residual whatever y is. The remaining n − r
entries of y do not affect the residual, and since `‖x‖ = ‖y‖`, zero is the choice of
least norm. So: form Qᵀb, back-substitute with T, pad with zeros, apply G₀ then G₁ up to
G₍r−1₎, and undo the column exchange. The output is written only after every value is
finite.

**Nullspace.** Column r + j of P Z is mapped by the truncated matrix to
`Q [[T, 0], [0, 0]] e₍r+j₎ = 0`. These n − r columns are orthonormal because P and Z are
orthogonal, and they are orthogonal to every solution `P Z [y₁; 0]`.

## Scope decisions

The planned stages named plain QR, full-rank least squares, column pivoting and the
complete decomposition. They are delivered as one factorization with pivoting always on:
a selected rank of n is the full-column-rank case, a selected rank of m < n is the
full-row-rank case that a transpose or LQ treatment would give, and neither needs a
separate entry point. An unpivoted variant was not added.

## Representation choices

`Decomposition` is a private nominal handle over a private record, as in the other
factorizations. The three status enums are distinct types; applications of Q and
nullspace requests share `SolveStatus` with solves. External compile probes cover
construction, destructuring, field access and confusion of two status types.
