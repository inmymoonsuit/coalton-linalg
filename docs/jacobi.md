# Symmetric eigenvalues and singular values by Jacobi rotations

Package: `coalton-linalg/dense/jacobi`. Source: [jacobi.ct](../src/dense/jacobi.ct).
State: completed as Batch 6; verified at safety 1 and safety 0. The [batch record](batch6.md)
holds the measurements and the controls.

Two computations share one package and one outcome type:

- `eigen!`: all eigenvalues and eigenvectors of a dense symmetric matrix, `A = V Λ Vᵀ`.
- `svd!`: all singular values and vectors of a dense matrix of any shape, `A V = U Σ`,
  computed from the columns of A directly. AᵀA is never formed.

Both are for small and moderate sizes. Their work is bounded by a sweep limit, and
reaching it is a named outcome.

| JacobiStatus | Meaning |
| --- | --- |
| JacobiSuccess | Results stand |
| JacobiNotConverged | The sweep limit was reached while rotations were still being made; no result stands |
| JacobiNonFinite | A non-finite input or computed value |
| JacobiBadSize | The array has the wrong number of elements |

Every call overwrites what stood before, and every outcome but success leaves nothing
standing. One handle owns its storage, so calls sharing it are not reentrant or concurrent.
Successful calls allocate nothing after warm-up on the measured runtime. Storage is refused
by name before allocation, as in the other dense packages: `AnalysisStorageExceeded` above
`default-entry-limit` (16,777,216 entries in any one array), `AnalysisLimitExceeded` above
`maximum-dimension`; the `-with-limit` constructors take the budget from the caller.

## Symmetric eigenvalues

A symmetric n-by-n matrix is an F64 array of n² elements in column-major order, of which
**only the lower triangle is read**.

- `make-symmetric-eigen(n) -> SymmetricEigen`: storage of 2n² + n doubles.
- `eigen!(eigen, a) -> JacobiStatus`, and `eigen-with-sweep-limit!(eigen, a, sweeps)` with an
  explicit bound in place of `default-sweep-limit`, 60.
- `eigenvalues(eigen)`: n values in **ascending** order.
- `eigenvectors(eigen)`: n-by-n column-major; column k is a unit eigenvector for value k.
- `eigen-sweeps(eigen)`: sweeps made by the latest call, counting the final sweep that
  rotated nothing. `eigen-ready?`, `eigen-dimension`.

The result arrays are borrowed, read-only views that the next call rewrites.

A sweep visits every pair p < q once, row by row. A pair is rotated when
`|a[p,q]| > 2⁻⁵² √|a[p,p]| √|a[q,q]|`, and the computation has converged when a whole sweep
rotates nothing. An entry that is exactly zero is never rotated.

**What is promised and tested.** The computed V is orthogonal and A V = V Λ to within small
multiples of rounding, and each eigenvalue is within the same of an exact eigenvalue of
the stored matrix, in the absolute sense, relative to the size of A. Eigenvectors of
repeated or tightly clustered eigenvalues are not unique; the tests compare residuals and
orthogonality, never vectors.

**What is observed and not promised.** The convergence test is relative to the diagonal.
For a positive definite matrix that is the test under which small eigenvalues are known to
come out to high relative accuracy, and the tests hold one graded matrix whose five
eigenvalues span 80 binary orders to a relative 2⁻⁴⁰ each. No claim is made for an
arbitrary symmetric matrix: an indefinite one can have tiny eigenvalues that no algorithm
working from its entries determines to any relative accuracy.

## Singular values

An m-by-n matrix is an F64 array of m n elements in column-major order. Every entry is read.

- `make-svd(m, n) -> Svd`: storage of m n + n² + n doubles. Either dimension may be zero.
- `svd!(svd, a) -> JacobiStatus`, and `svd-with-sweep-limit!(svd, a, sweeps)`.
- `singular-values(svd)`: **n** values in **descending** order, all non-negative.
- `right-vectors(svd)`: V, n-by-n column-major, a complete orthogonal matrix.
- `left-vectors(svd)`: U, m-by-n column-major; column j is `A vⱼ / σⱼ`, and is a zero column
  where σⱼ is exactly zero.
- `rank-above(svd, threshold) -> UFix`: the number of singular values **strictly larger than
  threshold**. The threshold is absolute; multiply by the first value for a relative one.
- `svd-sweeps`, `svd-ready?`, `svd-rows`, `svd-columns`.

**There are n singular values, not min(m, n).** They are the square roots of the n
eigenvalues of AᵀA. For a wide matrix at least n − m of them are zero in exact arithmetic;
they are computed at or below rounding level and counted out by any sensible threshold.
This is what makes V complete for every shape:

**The nullspace.** With r = `rank-above` at the caller's threshold, columns r to n − 1 of V
are an orthonormal basis of the nullspace of the rank-r matrix, of the full dimension
n − r, for tall, square and wide matrices alike. The first r columns of U are an orthonormal
basis of its range. Columns of U beyond r are whatever `A vⱼ / σⱼ` rounds to and mean
nothing.

**The rank is the caller's decision**, as in [the orthogonal decomposition](orthogonal.md),
with the same strict boundary rule. The difference is what is being compared with the
threshold: here it is the singular values themselves, the quantity a rank decision is
properly about, where there it is the diagonal of R, which can be off by a large factor.

A sweep visits every pair of columns once. With α and β the squared norms of the two
columns and γ their inner product, the pair is rotated when
`|γ| > 2⁻⁵² √α √β`, unless either column's norm is below 2⁻¹⁰⁴ times the largest entry of
A, in which case it is left alone. The computation has converged when a whole sweep rotates
nothing.

**What is promised and tested.** V is orthogonal, A V = U Σ, and the kept columns of U are
orthonormal, to within small multiples of rounding; each singular value is within the same
of an exact one, in the absolute sense relative to the largest column norm.

**What is observed and not promised.** The test is relative to the two columns. The tests
hold one matrix with graded columns whose three singular values span 40 binary orders to a
relative 2⁻³⁰ each. No general relative claim is made, and none at all for values below
2⁻¹⁰⁴ of the largest entry, where the floor stops the rotations.

**Cost.** Each pair costs three inner products of length m, and two column updates of
length m and two of length n when rotated; a sweep has n(n − 1)/2 pairs. The work therefore
grows with n², so **a wide matrix is far more expensive than its transpose**: a caller that
needs only the singular values of a wide matrix, or its left vectors, should pass the
transpose, whose singular values are the same and whose left and right vectors change
places. The complete nullspace basis of a wide matrix is what the direct orientation buys.

## Numerical assumptions and limits

Arithmetic is IEEE double with traps disabled. Both routines first scale the matrix by an
exact power of two when its largest entry is above 2⁵⁰⁰ or below 2⁻⁵⁰⁰, and scale the values
back. For singular values that keeps sums of squares representable; entries more than
about 500 binary orders below the largest are lost to underflow. For eigenvalues it keeps
the difference of two diagonal entries and every rotated entry representable: a symmetric
matrix with any finite entries is accepted, and the only overflow left is an eigenvalue
that is itself too large to represent, which is `JacobiNonFinite`. Scaling a matrix by a
power of two (an even power for the eigenvalue routine) scales the values exactly and
leaves the vectors unchanged bit for bit. A value below the normal range is returned on
the subnormal grid, so its relative accuracy is that grid's.

Rank-deficient and wide matrices take more sweeps than well-conditioned ones, because
columns at rounding level are also orthogonalized among themselves; the record has the
counts. The largest count on any tested matrix was 30, against the default limit of 60. The code is
unblocked and sequential, and it is slower than a tridiagonalizing or bidiagonalizing
method by a large factor that grows with the dimension.

## Artifact driver

The optional tools system writes reviewable outputs for a folder of systems:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/check.ps1 -Task Jacobi `
  -InputFolder <folder> -OutputFolder build/<name> -Repeats 1 -AllocationRepeats 1
```

Every `<name>.S.txt` is a symmetric eigenvalue problem: n, then n² values in column-major
order, of which the lower triangle is read. Every `<name>.M.txt` is a singular value
problem: m and n, then m n values in column-major order. One name may not be both.
`<name>.meta.txt` is optional. Each of its lines `sweep limit: <k>` is one run, in order
from index 0; with none there is one run at the default limit. For a singular value
problem each line `threshold: <value>` is not a run: every successful run reports the
count of singular values strictly above that value times the largest.

For run index i, a successful eigenvalue run writes `<name>.<i>.eigen.out.txt`: n, the n
values ascending, the n² vector entries column-major, then `sweeps <k>`. A successful
singular value run writes `<name>.<i>.svd.out.txt`: m, n, the n values descending, the
left vectors as m n entries and the right vectors as n² entries, column-major, then
`sweeps <k>` and one line `rank <count>` a threshold, in meta order. Every value is on its
own line, floating values with 17 significant digits. A refused run writes only
`<name>.<i>.refused.out.txt` holding the JacobiStatus name.

`eigen-table.txt` carries `name run n sweeps eigen_ms eigen_bytes` and `svd-table.txt`
carries `name run m n sweeps svd_ms svd_bytes`, one line a run, with `refused <Status>` in
place of the last three fields for a refusal. Times are the best of the requested repeats
after one warm-up call, and bytes are over the requested allocation repeats after one
warm-up call; both honour a count of 1, for matrices that take seconds a call. Reruns
replace every artifact of a system.

## Derivation

**The rotation.** For a symmetric 2-by-2 matrix `[[a, b], [b, d]]` and
`J = [[c, s], [−s, c]]` with c = cos θ, s = sin θ, the off-diagonal entry of `JᵀMJ` is
`(c² − s²) b + c s (a − d)`. It vanishes when `cot 2θ = (d − a)/(2b) = ζ`, that is, when
t = tan θ satisfies `t² + 2ζt − 1 = 0`. The root of smaller magnitude,
`t = sign(ζ) / (|ζ| + √(1 + ζ²))`, gives |θ| ≤ π/4, and with it the new diagonal entries
are `a − t b` and `d + t b`. Then `c = 1/√(1 + t²)` and `s = t c`. Where ζ is so large that
ζ² would overflow, t is `1/(2ζ)` to every representable digit.

**Eigenvalues.** Applying J to rows and columns p and q of the whole matrix, with the same
J accumulated into V, keeps `A = V W Vᵀ` and removes `2b²` from the sum of squares of the
off-diagonal entries: the two changed rows and columns keep their combined sums of squares,
and the entry b itself goes to zero. Repeating over all pairs drives the off-diagonal part
to zero; the diagonal of W is then Λ. The rotated entry is set to exactly zero and the two
new diagonal entries by the formulas above.

**Singular values.** Apply the same rotation to columns p and q of G, which starts as A,
with a, d and b taken as α, β and γ: the two new columns are orthogonal. Accumulating the
rotations in V keeps `G = A V`. When every pair of columns is orthogonal, the column norms
are the singular values, the normalized columns are left singular vectors, and
`A V = U Σ`. The two-sided computation on AᵀA would square the condition number; this one
only ever touches A's own columns.

## Representation choices

`SymmetricEigen` and `Svd` are private nominal handles over private records, as in the
other factorizations, and share one status enum because their outcomes are the same four.
External compile probes cover construction, destructuring, field access and confusion of
the two handles.
