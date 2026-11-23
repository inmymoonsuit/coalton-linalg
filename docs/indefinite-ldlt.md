# Dense symmetric indefinite LDLᵀ and inertia

Package: `coalton-linalg/dense/indefinite-ldlt`. Source:
[indefinite-ldlt.ct](../src/dense/indefinite-ldlt.ct). State: completed as Batch 4;
verified at safety 1 and safety 0. The [batch record](batch4.md) holds the measurements
and the controls.

## Interface

A symmetric n-by-n matrix is passed as an F64 array of n² elements in column-major
order: entry (i,j) is at index `i + j n`. **Only the lower triangle, i ≥ j, is read.**
The strict upper triangle may hold anything, including non-finite values.

- `make-factor(n) -> Factor` allocates reusable storage for n-by-n problems: n² + O(n)
  doubles. It refuses by name before allocating:
  `AnalysisStorageExceeded(n², budget)` when n² exceeds `default-entry-limit`
  (16,777,216 entries, n = 4096, 128 MB), and
  `AnalysisLimitExceeded(n, maximum-dimension)` when n exceeds `maximum-dimension`
  (32,768), which keeps n² representable.
- `make-factor-with-limit(n, maximum-entries) -> Factor` takes the entry budget from the
  caller. The budget counts the n² entries of L, not bytes; a host can still fail to
  allocate inside any budget.
- `factor!(factor, a, shift) -> FactorStatus` computes
  `P (A - shift I) Pᵀ = L D Lᵀ`, with L unit lower triangular and D block diagonal
  with 1-by-1 and 2-by-2 blocks. The input is copied; `a` is not retained or changed.
- `solve!(factor, b, out) -> SolveStatus` solves `(A - shift I) x = b` in the caller's
  ordering. b and out may be the same array. On any refusal out is unchanged.
- `positive-count`, `negative-count`, `zero-count` give the inertia of the computed D:
  the numbers of its positive, negative and exactly zero eigenvalues. They sum to n
  after a success and are all 0 when no factor stands.
- `dimension`, `factor-ready?`, and the inspection views below.

| FactorStatus | Meaning |
| --- | --- |
| FactorSuccess | A factor stands, **including one with exactly zero pivots**; read `zero-count` |
| FactorNonFinite | A non-finite lower-triangle input, shifted diagonal, pivot or multiplier |
| FactorBadSize | The array does not have n² elements |
| FactorBadShift | Non-finite shift |

| SolveStatus | Meaning |
| --- | --- |
| SolveSuccess | Output contains the solution |
| SolveInvalidFactor | No successful factor stands |
| SolveBadSize | Input or output vector length is not n |
| SolveSingular | The factor has at least one zero pivot; no solve is attempted |
| SolveNonFiniteRhs | Non-finite right-hand side |
| SolveNonFiniteResult | Non-finite computed solution |

Checks are made in the order of that table: an invalid factor is reported before a
size mismatch, a singular factor before the right-hand side is read.

Every `factor!` call overwrites the factor; a refusal leaves none standing and reports
no inertia. One Factor owns its storage and scratch, so calls sharing it are not
reentrant or concurrent. Successful `factor!` and `solve!` calls allocate nothing after
warm-up on the measured runtime.

## The factor

Inspection functions return borrowed, read-only views that the next `factor!` rewrites:

- `factor-permutation`: `permutation[k]` is the original index at pivoted position k.
- `factor-lower`: L as a full n-by-n column-major array, with its unit diagonal and
  zeros above it stored explicitly.
- `factor-diagonal`: the n diagonal entries of D.
- `factor-subdiagonal`: entry k is `D[k+1,k]` where a 2-by-2 block starts at k, and 0
  elsewhere.
- `factor-block-sizes`: 1 at a 1-by-1 pivot, 2 at the first row of a 2-by-2 block and 0
  at its second row.

Inside a 2-by-2 block `L[k+1,k]` is 0. Every multiplier in L is at most 1/(1 - α) ≈ 2.78
in magnitude, with α = `pivot-threshold` = (1 + √17)/8; a multiplier of a 1-by-1 pivot
is at most 1/α ≈ 1.56.

## Inertia, and what it does and does not say

By Sylvester's law of inertia, a congruence `Sᵀ M S` with nonsingular S has as many
positive, negative and zero eigenvalues as M. `A - shift I` is congruent to D, so in
exact arithmetic the inertia of D is the inertia of `A - shift I`: its negative count is
the number of eigenvalues of A below the shift, its zero count the number equal to it.

**Each 2-by-2 block contributes one positive and one negative eigenvalue**, whatever
the signs on its diagonal. A block is taken only when both diagonal entries are below
α times its off-diagonal entry in magnitude, so its determinant is negative.
`[[1,2],[2,1]]` is factored as one such block: both diagonals are positive, and the
matrix has eigenvalues 3 and -1. Counting blocks, or signs of diagonal entries, would be
wrong.

**The counts describe the computed D.** The computed factor reproduces a matrix within
the reconstruction error the batch record measures, and a count is reliable only for
eigenvalues further from the shift than that error. An eigenvalue at the shift in exact
arithmetic usually becomes a tiny pivot of either sign, not a zero. `zero-count` counts
pivots that are exactly zero in the computed arithmetic, which happens for a zero
matrix, for exact small-integer cancellations, and rarely otherwise. A caller deciding
rank near a threshold chooses the shift and reads the counts on both sides of it; this
package makes no rank decision.

**An exactly zero pivot completes.** When the remaining matrix has an entirely zero pivot
column, that pivot is recorded as zero with zero multipliers and elimination continues.
The result is `FactorSuccess` with `zero-count` above 0, and `solve!` then returns
`SolveSingular`. **An exactly singular input does not always get there**: rounding in
earlier steps usually leaves a pivot of rounding size instead, the zero count is 0, and a
solve succeeds with a very large answer. Reading the counts either side of a threshold of
the caller's choosing is the reliable test, and
[the orthogonal decomposition](orthogonal.md) is the tool when rank itself is the question. The factorization does not solve singular systems: whether a singular
system is consistent, and which solution is wanted, are different problems.

For an m-by-n matrix J of rank r, the augmented matrix `[[I, Jᵀ], [J, 0]]` is congruent
to `diag(I, -J Jᵀ)` and has inertia (n, r, m - r). With full row rank it is nonsingular
and the first n unknowns of `[[I, Jᵀ], [J, 0]] [x; y] = [0; b]` are the minimum-norm
solution of `J x = b`; [a compiled example](../examples/minimum-norm.ct) does this. The
inertia of `J Jᵀ` concerns m eigenvalues and is not the nullity n - r of J.

## Numerical assumptions and limits

The diagonal is formed as the rounded difference `a[j,j] - shift`, one rounding per
entry before factoring. No scaling, equilibration, iterative refinement, tolerance or
perturbation of small pivots is applied. Arithmetic is IEEE double with traps disabled;
overflow is detected by finiteness checks on every pivot and multiplier and is a
refusal. Underflow follows F64 and is not diagnosed.

A 2-by-2 block `[[a, b], [b, c]]` is inverted in the form divided by b: its scaled
determinant `(a/b)(c/b) - 1` lies strictly between -(1 + α²) and -(1 - α²), since the two
diagonals may differ in sign. At the ends of the exponent range, a common power of two
scales b and the trailing entries or right-hand side before the divisions. This avoids
both an overflowed denominator and a prematurely overflowed right-hand side divided
by b, even on perfectly conditioned blocks with representable solutions.

Storage is dense: n² doubles for L plus O(n). Work is about n³/3 multiply-adds for the
factorization, plus the pivot search described below, and n² for a solve. The code is a
plain column-oriented loop, with no blocking and no parallelism. Dense delivery says
nothing about large sparse augmented systems.

## Artifact driver

The optional tools system writes reviewable outputs for a folder of systems:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/check.ps1 -Task Indefinite `
  -InputFolder <folder> -OutputFolder build/<name> -AllocationRepeats 20
```

Each `<name>.S.txt` holds n, then n² values in column-major order, of which only the lower
triangle is read. `<name>.b.txt` holds n values. `<name>.meta.txt` is optional; each of its
lines `shift: <value>` names one run, in order, and with none there is one run at shift 0.

For run index i from 0, a successful factorization writes `<name>.<i>.factor.out.txt`: n,
the permutation, the block sizes, the diagonal, the subdiagonal, L column-major, each value
on its own line with 17 significant digits, then `counts <positive> <negative> <zero>`. The
solve writes `<name>.<i>.x.out.txt`, or `<name>.<i>.solve.out.txt` holding the SolveStatus
name when it is refused. A refused factorization writes only `<name>.<i>.refused.out.txt`
holding the FactorStatus name. `indefinite-table.txt` and standard output carry one line a
run: `name shift n positive negative zero factor_ms solve_ms factor_bytes`, with `NA` for the
time of a refused solve. Reruns replace every indexed artifact of a system.

## Derivation

Take the trailing symmetric matrix at some step and write its entries a.

**A 1-by-1 pivot.** Choosing diagonal entry d = a[i,i] and exchanging i to the front
gives multipliers `a[t,i] / d` and the Schur complement
`a[t,u] - a[t,i] a[u,i] / d`. If ω is the largest off-diagonal magnitude in column i,
every multiplier is at most ω/|d|. Requiring |d| ≥ α ω bounds them by 1/α.

**A 2-by-2 pivot.** Choosing rows i and r with the block `E = [[a, b], [b, c]]` gives
the multiplier rows `[a[t,i], a[t,r]] E⁻¹` and the Schur complement
`a[t,u] - [a[t,i], a[t,r]] E⁻¹ [a[u,i], a[u,r]]ᵀ`. Suppose |b| = ω is the largest
off-diagonal magnitude of column i and also of column r, and |a| < α ω, |c| < α ω. Then
`det E = ac - b²` is at most `-(1 - α²) ω²`, and each multiplier is at most
`(α ω · ω + ω · ω) / ((1 - α²) ω²) = 1/(1 - α)`.

**The search.** Start at the first remaining column i with largest off-diagonal ω at
row r. If ω is 0 or |a[i,i]| ≥ α ω, take i as a 1-by-1 pivot. Otherwise look at column
r, whose largest off-diagonal ω' is at least ω because it includes a[r,i]. If
|a[r,r]| ≥ α ω', take r as a 1-by-1 pivot. If ω' = ω, the entry a[r,i] is the largest
off-diagonal of both columns and both diagonals have failed their tests, which are
exactly the conditions above: take the pair. Otherwise ω' > ω; move to column r and
repeat. The off-diagonal maximum strictly increases at each move, so no column is
visited twice and the search ends within the remaining dimension. A zero pivot can only
be taken when ω is 0, that is, from a column that is entirely zero.

**The threshold.** With μ the largest magnitude in the matrix, a 1-by-1 step under the
test |d| ≥ α ω leaves entries of at most μ(1 + 1/α), so two such steps at most
μ(1 + 1/α)². A 2-by-2 step whose off-diagonal is the largest entry of the whole matrix
leaves at most μ(1 + 2/(1 - α)). Equating the two bounds gives 4α² - α - 1 = 0, whose
positive root is (1 + √17)/8. This is the classical choice, derived for a search over
the whole matrix; here it is used as the threshold of the column search above, for
which the multiplier bounds hold as stated and no bound on entry growth is claimed.

**Exchanges on a lower triangle.** Exchanging indices p < q of a symmetric matrix held
in its lower triangle swaps rows p and q in the columns before p, the two diagonal
entries, the column segment below q, and the entries (m,p) with (q,m) for p < m < q;
entry (q,p) stays. The already computed columns of L take the same row exchange.

**Solve.** Gather Pb, solve Ly = Pb, solve Dz = y block by block, solve Lᵀw = z, then
scatter x = Pᵀw. The output is written only after every computed value is finite.

## Representation choices

`Factor` is a private nominal handle over a private record, `(repr :transparent)` as in
the sparse factorization; external compile probes cover construction, destructuring,
field access and confusion of this package's `FactorStatus` with the sparse one. The
shift argument and the status enums follow [the sparse SPD contract](spd-ldlt.md) so
the two factorizations read alike; the types are distinct.
