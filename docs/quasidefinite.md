# Signed pivots and augmented systems

Packages: `coalton-linalg/sparse/spd-ldlt` (the signed entry points) and
`coalton-linalg/sparse/augmented` (the plan). Sources:
[spd-ldlt.ct](../src/sparse/spd-ldlt.ct), [augmented.ct](../src/sparse/augmented.ct),
[ordering.ct](../src/sparse/ordering.ct). State: completed as Batch 7; verified at safety 1
and safety 0. The [batch record](batch7.md) holds the measurements and the controls.

## What it is for

A wide matrix A with one dense column has a dense A Aᵀ, so the minimum-norm step
`x = Aᵀ (A Aᵀ)⁻¹ b` through a positive definite factorization costs dense work on a sparse
problem. The augmented matrix

    K = [[top I,  Aᵀ       ],
         [A,      −bottom I]]

is as sparse as A, and `K [x; y] = [0; b]` gives `x = Aᵀ (A Aᵀ + top·bottom I)⁻¹ b`. K is
symmetric and indefinite: with top and bottom positive it has n positive and m negative
eigenvalues for an m-by-n A. Factoring it needs pivots of both signs.

**This route is not accurate on every matrix, and it says when it is not.** No elimination
order fixed in advance keeps the factorization of K accurate for all values. Every signed
factorization therefore measures its own growth, exactly, and a caller reads that number
before trusting the factor. See *Accuracy, growth and the elimination order*.

## Interface

In `coalton-linalg/sparse/spd-ldlt`, beside the positive definite entry points and over the
same `Symbolic` and `NumericFactor`:

- `factor-signed!(symbolic, values, shift) -> SignedFactorStatus` computes
  `P (A − shift I) Pᵀ = L D Lᵀ` with D diagonal, **accepting pivots of either sign and
  refusing only a pivot that is exactly zero**. There is no numerical pivoting: the order
  is the one fixed at analysis.
- `factor-growth(factor) -> F64`: the largest entry of `|L| |D| |Lᵀ|` divided by the largest
  magnitude of the matrix that was factored, as measured by the latest successful
  `factor-signed!`; −1 when none stands, and after `factor!`, which does not measure it.
- `positive-pivot-count(factor)`, `negative-pivot-count(factor) -> UFix`: the signs of D in
  the factor as it now stands, 0 and 0 when none stands. After `factor!` they are n and 0.
  By Sylvester's law of inertia they are the numbers of positive and negative eigenvalues
  of the matrix the factor reproduces.
- `analyze-quasidefinite(pattern, leading) -> Symbolic`, and
  `analyze-quasidefinite-with-limit(pattern, leading, maximum-entries)`: minimum-degree
  analysis in which **the indices from `leading` on are held back, and each eliminated
  index below `leading` releases one of them**. `leading` above the size raises
  `InvalidDimensions`. With `leading` 0 or equal to the size the order is plain minimum
  degree.
- `solve!`, `refused-column`, `estimate-inverse-norm1!` and the inspection views work on a
  signed factor exactly as on a positive definite one.

| SignedFactorStatus | Meaning |
| --- | --- |
| SignedSuccess | A factor stands; read the growth and the pivot counts |
| SignedZeroPivot | A computed pivot was exactly zero; `refused-column` names it |
| SignedNonFinite | Non-finite input or arithmetic |
| SignedBadValues | Wrong value count |
| SignedBadShift | Non-finite shift |

`SignedFactorStatus` and `FactorStatus` are distinct types. Every factorization call on a
Symbolic, of any of the three kinds, overwrites its one factor, and a refusal withdraws the
factor, the pivot counts, the growth and the stored condition estimate. Successful calls
allocate nothing after warm-up on the measured runtime, and `solve!` may be called any
number of times on fresh right-hand sides against one factor.

In `coalton-linalg/sparse/augmented`:

- `analyze(matrix) -> AugmentedPlan` for a checked CSC matrix A, m-by-n, and
  `analyze-with-limit(matrix, maximum-entries)`. The plan owns a copy of A's pattern and
  the lower CSC pattern of K, of order n + m with n + nnz(A) + m stored entries: every
  diagonal, and A's entries below the leading block. Stored zeros of A stay structural.
  The budget counts those entries; `AnalysisStorageExceeded` and `AnalysisLimitExceeded`
  are raised as elsewhere.
- `dimension(plan)` is n + m; `leading-dimension(plan)` is n, the argument
  `analyze-quasidefinite` wants; `column-pointers`, `row-indices` and `entry-count` give
  the pattern for `ldlt:Pattern`.
- `refill!(plan, values, top, bottom, out) -> Void` writes K's values in the plan's order
  from A's values in its analyzed CSC order and the two scalars. It checks both lengths,
  that `out` is not the input array, and that the scalars and every value are finite,
  before the first write; no arithmetic is done, so a call either writes everything or
  nothing.
- `refill-csc!(plan, matrix, top, bottom, out)` first compares the matrix's dimensions and
  whole pattern with the analyzed one.

The unknowns of the leading block come first in every vector: for `K [x; y] = [0; b]` the
right-hand side holds n zeros and then b, and x is the first n entries of the solution.

A caller's loop: `analyze` the plan and `analyze-quasidefinite` once for a pattern; then
for each set of values `refill!`, `factor-signed!`, read the growth and the counts, and
`solve!`. [A compiled example](../examples/augmented-minimum-norm.ct) does this.

## What the factor is a factor of

**The regularized matrix.** The solution is `Aᵀ (A Aᵀ + top·bottom I)⁻¹ b`, which differs
from the minimum-norm solution `Aᵀ (A Aᵀ)⁻¹ b` by a relative amount of order
`top·bottom / σ_min²`, with σ_min the smallest singular value of A. A caller that needs
the unregularized solution corrects it outside this library, for example by iterative
refinement with residuals from `matvec!` and `transpose-matvec!` and corrections from
repeated `solve!` calls. With bottom zero K is not quasi-definite, a dependent or empty
row gives a zero pivot, and `SignedZeroPivot` is returned.

**The inertia is a count, not a rank decision.** In exact arithmetic the counts are n and m
for positive top and bottom and any A, dependent rows included. Computed, they are n and m
whenever the factorization behaved as the theory requires. With dependent rows and a
bottom at rounding level the last such row's pivot is what cancellation leaves of bottom:
the call may then refuse a zero pivot, or complete with other counts, and counts other than
n and m are the warning. Neither outcome says anything about the rank of A.

## Accuracy, growth and the elimination order

**The bound.** Write K̃ for the matrix factored, with its diagonal as rounded when the shift
was subtracted, u = 2⁻⁵³, `γₘ = m u / (1 − m u)`, and w for the longest row or column of L,
counting entries off the diagonal. The computed factor and every computed solution x̂
satisfy, entry by entry,

    | P K̃ Pᵀ − L D Lᵀ |  ≤  γ₍p+2₎ |L| |D| |Lᵀ|,     p ≤ w the products that entry received,
    | K̃ x̂ − b |          ≤  γ₍3w+5₎ |L| |D| |Lᵀ| |x̂|   (rows and columns as permuted).

Since no entry of `|L| |D| |Lᵀ|` exceeds `factor-growth · max|K̃|`,

    max | K̃ x̂ − b |  ≤  γ₍3w+5₎ · N · factor-growth · max|K̃| · max|x̂|,

with N the order of K. Both entrywise bounds are checked exactly, in rational arithmetic,
on every entry and every solve of the test rosters. So **growth is the number that says
how much accuracy the factorization kept**: near 1, the solve is as accurate as the
problem allows; at 10⁸, eight digits of the residual are gone.

**Growth is measured, not estimated.** `|L| |D| |Lᵀ|` is the Gram matrix of the rows of
`|L| |D|^½`, so its largest entry is on its diagonal, and its diagonal is
`Σₖ l[i,k]² |d[k]|`: one multiply-add for each entry of L as it is formed. The reported
figure is that maximum over the largest magnitude of K̃, in floating point.

**Why the order matters.** A matrix `[[H, Aᵀ], [A, −G]]` with H and G positive definite is
called quasi-definite, and every symmetric permutation of it has an L D Lᵀ factorization
with D diagonal. That guarantees nonzero pivots, not large enough ones. Under plain
minimum degree a row of A with few entries is eliminated before any of its columns; its
pivot is exactly −bottom and its multipliers are A's entries over bottom, and growth is
about `2 / bottom` relative to A's scale.

**What `analyze-quasidefinite` does about it.** Rows are held back, and each eliminated
column releases one row among those it has an entry in or has come to be joined to, the
one of least degree. So no row is eliminated before a column of its own, and no two rows
are released by the same column: a row's pivot is then `−(bottom + a²/p + …)` with a an
entry of A and p a column's pivot, and no longer shrinks with bottom. Among released rows
and columns the order is still minimum degree with smallest-index ties, so a dense column
is still left until late and its rows are not joined into a dense block. If no free vertex
is left, held ones are taken, least degree first.

**What it does not do.** The rule reads the pattern, never the values. Rows can be
independent and still nearly dependent when restricted to the columns eliminated before
them, and geometric constraint matrices are often exactly that. On the record's grid
system growth rises as one over bottom under this rule too; on its hub, chain and frame
systems growth is between 6 and 244 at every bottom down to rounding level. **A caller
decides per factorization, from `factor-growth`, whether to use the factor.** The
positive definite route through A Aᵀ and `factor!` is the accurate alternative, and is
costly only where A has a dense column, which is where the held order works.

`factor-signed!` over a plain `analyze` remains available and is appropriate when the two
diagonal blocks are comparable in size, or for a symmetric matrix that is safely
factorizable in any order, such as one with a strictly dominant diagonal of mixed signs.

## Artifact driver

The optional tools system writes reviewable outputs for a folder of systems:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/check.ps1 -Task Augmented `
  -InputFolder <folder> -OutputFolder build/<name>
```

Every `<name>.J.txt` is a matrix A: `rows columns nnz`, then zero-based `row column value`
triplets in any order. Each line `augmented: <top> <bottom> <held|plain>` of
`<name>.meta.txt` is one run, in order from index 0: `held` analyzes with
`analyze-quasidefinite` at the plan's leading dimension, `plain` with `analyze`. A system
with a run needs `<name>.b.txt`, one value a row of A; the right-hand side is n zeros and
then those values. A system whose meta file names no run is not run, and the table says
`<name> skipped`, so a folder can be shared with the Gram driver.

For run index i the driver writes:

- `<name>.<i>.K.out.txt`: the lower triangle of K as refilled, `size nnz` and then
  `row column value` by column.
- `<name>.<i>.signed.out.txt`: six lines, `status`, `refused-column`, `positive`,
  `negative`, `growth` and `solve`, each followed by its value; `solve none` when no factor
  stood.
- After `SignedSuccess`: `<name>.<i>.perm.out.txt`, `.D.out.txt` and `.L.out.txt`, the
  permutation, the pivots in elimination order and the strict lower triangle of L as
  `size nnz` and triplets in permuted numbering; and after a successful solve
  `<name>.<i>.x.out.txt`, the n + m solution values in original order.

Floating values carry 17 significant digits. A system the library refuses with a typed
error, at construction or analysis, writes only `<name>.refused.out.txt` holding the
constructor's name. `augmented-table.txt` carries one line a run:
`name run m n order top bottom nnzK nnzL status positive negative growth analyse_ms
refill_ms factor_ms solve_ms bytes`, with nnzL counting the unit diagonal. Times are the
best of the requested repeats; bytes are allocated by at least 1,000 repetitions of
refill, factorization and solve after warm-up, and a nonzero count fails the run. Reruns
replace every artifact of a system.

## Derivation

**Nonzero pivots.** Let K be quasi-definite and `K z = 0` with z = (x, y). Then
`H x + Aᵀ y = 0` and `A x − G y = 0`; multiplying the first by xᵀ and the second by yᵀ and
subtracting gives `xᵀ H x + yᵀ G y = 0`, so x and y are zero and K is nonsingular. A
symmetric permutation of K restricted to any leading set of indices is again of that form,
with sub-blocks of H and G, so every leading principal submatrix is nonsingular, which is
the condition for L D Lᵀ with diagonal D to exist.

**Inertia.** K is congruent to `diag(H, −(G + A H⁻¹ Aᵀ))`, which has n positive and m
negative eigenvalues; a congruence preserves the counts, and so does the factorization.

**A pivot depends only on what was eliminated before it.** The pivot of an index is the
ratio of two leading principal minors, so it is fixed by the set S of indices eliminated
earlier, not by their order. For a row r of A, with S made of columns C and rows R and with
top = 1, it is `−(bottom + ‖a‖²)` up to terms that vanish with bottom, where a is the part
of row r, restricted to the columns C, that is orthogonal to the rows R restricted to C. So
a row's pivot is small exactly when the row, cut down to the columns already eliminated, is
nearly a combination of the rows already eliminated cut down the same way. Holding a row
until it has a column of its own removes the structural cases of that; the numerical ones
remain, and growth reports them.

**The rounding bounds.** In the kernel a strict lower entry is
`l[i,k] = (a[i,k] − Σⱼ (l[k,j] d[j]) l[i,j]) / d[k]`: each of its p terms is rounded twice
as a product of three numbers, then subtracted, then one division follows. Counting the
roundings each term passes through gives at most p + 2, hence γ₍p+2₎; the diagonal is the
same without the division. In a solve, forward substitution with a unit triangle perturbs
row i of L by at most γ of its length, the division by D by u, and back substitution
column k of L by at most γ of its length; adding the three to the factorization's bound
and absorbing second-order terms gives γ₍3w+5₎.

**The augmented solve.** From `top x + Aᵀ y = 0` and `A x − bottom y = b`:
`x = −Aᵀ y / top` and `−(A Aᵀ / top + bottom I) y = b`, so
`x = Aᵀ (A Aᵀ + top·bottom I)⁻¹ b`.

## Representation choices

The signed entry shares one numeric kernel with `factor!`; the two differ in one test a
column and in whether growth is accumulated, and cannot drift apart. They live in the
package named for the positive definite case because they share its `Symbolic`, its
`NumericFactor` and its solve; a package name that fits all three numeric entries is a
naming decision left open. The plan follows [the Gram plan](gram.md): a private handle,
borrowed pattern views, caller-owned values.
