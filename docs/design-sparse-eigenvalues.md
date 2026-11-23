# Design: a few eigenvalues of a sparse symmetric matrix

State: **design; nothing is built.** It was read once by an independent reviewer, on
2026-11-23, who re-derived the theorem of *Several pairs* and the count and found them
sound, and found one decision wrong and eight smaller things; all are worked in below
and marked where the first version said otherwise. It answers the roadmap item "a few eigenvalues of a
sparse symmetric matrix: Lanczos, shift-invert through the signed factorization, and the
inertia as a count that none was missed". No consumer asked for it; it is for the
library's own users, and it changes no existing call. Where it says "in exact arithmetic"
it means a property of the method that rounding weakens and that the tests would measure,
not one a call proves. Where it says "proved" it says what from.

## The question

A real symmetric matrix A of order n is too large for
[the dense routine](jacobi.md), which computes all n eigenvalues at a cost that grows
with n³. A caller wants a few of them, with their vectors: the largest, the smallest, or
those nearest a number σ it names. And it wants to know that the ones it was given are
all there are in some interval, which no iteration can tell it.

Three things are wanted of the answer, in this library's usual order:

1. **What a call reports of a pair is formed from the operator**: for a unit vector y
   and a number θ, the residual `A y − θ y`. An eigenvalue of A lies within its norm of θ,
   and that is the whole of what an iteration can say.
2. **That none was missed is a separate statement with a separate proof**, from the
   signs of the pivots of a factorization, and it is either made or refused by name.
3. **Every outcome has a name, a refusal is made before anything is written, and the
   size of the problem does not limit it**: the same problem at another power of two is
   solved bit for bit the same.

## The pieces, and the order they would be built in

Four pieces, each usable without the ones after it, each its own batch.

1. **The eigenvalues of a symmetric tridiagonal matrix**, with as much of the vectors as
   is asked for. New package `coalton-linalg/dense/tridiagonal`. The iteration needs it
   at every step; it is also the second half of a dense symmetric routine faster than
   Jacobi rotations, whose first half is not in this design.
2. **The iteration on an operator**, package `coalton-linalg/lanczos`: a few pairs at
   either end of the spectrum of a symmetric operator, each with its residual formed.
3. **What a set of pairs proves together**, in the same package: one number δ, from the
   residuals and from how far the vectors are from orthonormal, such that as many
   eigenvalues as there are pairs lie within δ of the values, counted with multiplicity.
4. **The count**, package `coalton-linalg/sparse/inertia-count`: how many eigenvalues of
   a stored matrix lie below a number, from [the signed factorization](quasidefinite.md),
   with a margin that is proved; and the comparison of a count with the pairs found.

Shift and invert is not a piece. It is the iteration of piece 2 given the solve of an
existing factorization as its operator, and the contract page would show it.

## The iteration

The method Lanczos published in 1950, with every vector kept and made orthogonal to all
the earlier ones. It is set out here from the mathematics, because the code is to be
written from this page.

### The recurrence

B is the symmetric operator iterated on: A itself, or a function of it. From a unit
vector q₁, unit vectors q₂, q₃, … and numbers α_j, β_j ≥ 0 are defined by

    β_j q_{j+1} = B q_j − α_j q_j − β_{j−1} q_{j−1},        α_j = q_jᵀ B q_j,

β_j being the norm that makes q_{j+1} a unit vector, and β₀ q₀ = 0. In exact arithmetic the
q are orthonormal, and with `Q_j = [q₁ … q_j]` and T_j the symmetric tridiagonal matrix
with α₁ … α_j on its diagonal and β₁ … β_{j−1} beside it,

    B Q_j = Q_j T_j + β_j q_{j+1} e_jᵀ.

Each step costs one product with B.

### Pairs from the small matrix

Let θ be an eigenvalue of T_j and s its unit eigenvector, and `y = Q_j s`. Multiplying
the relation above by s,

    B y − θ y = β_j s_j q_{j+1},        so        ‖B y − θ y‖ = β_j |s_j|,

s_j being the last entry of s. So the eigenvalues of T_j approximate eigenvalues of B,
and how well is read from the last entries of its eigenvectors, without forming y. Those
at the two ends of the spectrum are approximated first and best. This estimate is exact
in exact arithmetic; here it only decides when to form the true quantity.

### Keeping the vectors orthogonal

In floating point the three-term recurrence loses the orthogonality of the q as soon as
a pair has converged, and then finds the same eigenvalue again, any number of times
(Paige, 1971). A count of eigenvalues cannot be compared with pairs that may be copies.
So the basis is kept, and what is left of `B q_j` after the two terms of the recurrence
is made orthogonal to every earlier q, by a rule in three branches:

1. One pass. If what is left is at least half as long as what the pass was given, it is
   taken: the pass removed little, and its own rounding is small beside what is left.
2. If not, a second pass, on what the first left. If what that leaves is at least half
   as long as what it was given, it is taken.
3. If not, **the vector is in the span of the basis as far as the arithmetic can tell,
   and the sequence has ended.** Each pass has taken away most of what it was given; what
   is left is the rounding of rounding, and has no direction of its own.

This is Kahan's rule as Parlett gives it, as far as memory serves: the book was not at
hand, and the reviewer's recollection of it is the same. What the passes subtract is
discarded, not added to T: with an orthonormal basis it is of the size of rounding, since
`q_iᵀ B q_j = (B q_i)ᵀ q_j` and `B q_i` lies along q_{i−1}, q_i and q_{i+1}. The third
branch is what keeps the basis orthonormal, and so keeps that sentence true.

The cost is storage of the basis, n numbers a step, and work of 4 j n a step. That is
the price of a count that means something, and it is why the step limit is the caller's
and is meant to be tens or a few hundreds, not thousands.

### Where the sequence ends

If `β_j = 0` the span of q₁ … q_j is invariant: every eigenvalue of T_j is an eigenvalue
of B exactly, and there is nothing to continue with. The iteration ends there, and it
ends by the third branch above as well as by an exact zero.

**The first version of this page ended the sequence on an exact zero alone, and that was
wrong.** Its reason was that a remainder at rounding level, made a unit vector, is
orthogonal to the basis like any other. That holds only while rounding can still leave
the span. Once the basis spans all the start can reach, which is the whole space, a
block of a reducible operator, or the complement of the pairs that stand, the remainder
is inside the span and so is the vector made of it. The reviewer's case: the matrix with
rows (2, 1, 0, 0), (1, 3, 1, 0), (0, 1, 5, 1), (0, 0, 1, 7) from the start (1, 2, 3, 4).
After the fourth step the norm is 4.7 × 10⁻¹⁶ after the recurrence, 3.7 × 10⁻³² after one
pass and 9.5 × 10⁻⁴⁸ after two. Made a unit vector, the fifth has inner products of order
1 with the first four, and T₅ has an eigenvalue 3.1561 that the matrix does not have,
with an estimate of 1.4 × 10⁻¹⁶ and a formed residual of 1.7. A call for the three
smallest at a tolerance under rounding would have returned that value twice and left out
two pairs it had held to rounding. No status would have been false, the residuals being
formed; the pairs would have been inventions.

The test is of the passes and not of the norm's size: no threshold beside the operator
decides that a small norm is zero. A norm at rounding level that a pass leaves nearly
whole is a direction, and is taken.

### One copy of a multiple eigenvalue

Every q is a polynomial in B applied to q₁, so the vectors span no more of an eigenspace
of B than q₁ has in it: one direction. **A call finds one copy of a multiple eigenvalue**,
in exact arithmetic, and in floating point a second only late and by rounding, which is
not promised. A second copy is found by a second run from another starting vector, kept
orthogonal to the vectors already found; the interface has that call. Whether there is a
second copy is what the count says.

**In that second run the estimate has a second term.** Keeping every vector orthogonal to
the vectors Y of the pairs that stand is the iteration on `(I − P) B (I − P)`, P the
projector onto their span, and what is discarded along Y is not rounding: with
`R = B Y − Y Θ` their residuals, and y orthogonal to Y,

    B y − θ y = β_j s_j q_{j+1} + P B y,        P B y = Y (Yᵀ Y)⁻¹ Rᵀ y,

so `‖B y − θ y‖² = (β_j s_j)² + ‖P B y‖²`, and the second term does not fall as j grows.
It is at most `ρ_Y / √(1 − γ_Y)`, ρ_Y the root of the sum of the squared residual norms of
the pairs that stand and γ_Y their Gram matrix's distance from the identity. The estimate
of the second run carries that bound as its second term. So **a tolerance below the
residuals of the pairs kept is not met by the estimate**, the run goes to its limit, and
the residuals formed there say what was reached: pairs to be kept for a second run are
first dropped down to those at least as good as the tolerance to be asked. The reviewer
measured a standing residual of 1.3 × 10⁻³ leaving the second run's best pair at
3.5 × 10⁻⁴ under an estimate, without the term, of 10⁻¹⁵.

### The small problem

All eigenvalues of T_j are wanted at every step, with the last entries of their
vectors; and at a confirmation the vectors whole. Both are one computation: shifted QR
steps on the tridiagonal matrix, applied implicitly. A step with shift μ is a sequence
of plane rotations applied to T as similarities. The first is the rotation that would
make the second entry of the first column of `T − μ I` zero; applied on both sides of T
it removes nothing and creates one entry outside the three diagonals. Each rotation
after it is the one that removes that entry, and moves it one place down, until it
leaves the matrix. The shift is the eigenvalue of the trailing two-by-two block nearer
its last entry (Wilkinson's): with a and b the block's diagonal entries, b the last, e
its off-diagonal entry and `d = (a − b) / 2`, it is `b − e² / (d + sign(d) √(d² + e²))`,
**and where the two diagonal entries are equal, so that both eigenvalues are equally
near, it is the smaller, `b − |e|`.** A rule that fell back to the last entry itself there
would leave the block `[[a, e], [e, a]]` where it stood for ever. The rotations are also applied, from the right, to a
matrix C of j columns that the caller supplies, which then holds C Z, Z being the matrix
of eigenvectors. With C the one row e_jᵀ that is the last entries, at a cost of j² for
the whole decomposition; with C the identity it is Z, at j³.

An off-diagonal entry is taken to be zero when its magnitude is at most 2⁻⁵³ times the
sum of the magnitudes of its two diagonal neighbours, and the matrix then splits there.
Between two diagonal entries that are both zero that test never passes for an entry that
is not zero, and convergence there rests on the steps bringing it to zero exactly or on
the bound.
The number of steps is bounded, thirty for each eigenvalue, and reaching the bound is a
named outcome. The method is that of Francis (1961) in the form used for symmetric
tridiagonal matrices; it converges for every such matrix with this shift, in exact
arithmetic, and the bound is there because a proof about exact arithmetic is not a proof
about this arithmetic.

**Why not the dense routine that exists.** Jacobi rotations on T_j as a dense matrix
would give the same numbers at j³ times a number of sweeps, where the last entries
alone cost j² here. Over a run of m steps that is m⁴ against m³; neither is measured.
The alternative, to solve the small problem only every few steps, makes every run take
steps it did not need, each a product or a solve.

## Shift and invert

For the eigenvalues of A nearest σ, iterate on `B = (A − σ I)⁻¹`: an eigenvalue λ of A is
the eigenvalue `1 / (λ − σ)` of B with the same vector, so those nearest σ are the
largest in magnitude, at the two ends of B's spectrum, where the iteration is best, and
well separated there **from the rest of the spectrum**. They are not separated from one
another by it: the gap between two of them, against their own size in B, is
`|λ₁ − λ₂| / |λ₂ − σ|`, small whenever the two are closer to each other than to the shift,
and the iteration resolves such a pair slowly. (The first version said "well separated
however close together they are in A", which is false of the pair.)

B is applied by `solve!` on the factor that `factor-signed!(symbolic, values, σ)` leaves:
one factorization, then one solve a step. Nothing new is built for it. What has to be
said of it:

- **The pairs wanted are A's.** After the iteration the vectors are measured again on A:
  each value becomes the Rayleigh quotient `yᵀ A y`, each residual `A y − θ y`, by one
  product with A a vector. The tolerance of the iteration was on B's residual, and
  `A y − λ y = −(λ − σ) (A − σ I) (B y − θ y)` for `λ = σ + 1/θ`, so it bounds A's only
  through the size of `A − σ I`, which the call does not know. The residuals on A are
  what is reported and what the bound below is formed from.
- **The factorization does not pivot.** Its [growth](quasidefinite.md#accuracy-growth-and-the-elimination-order)
  says how accurate the solve is, and a σ at or beside an eigenvalue of a leading block
  is refused as a zero pivot or factors with large growth. Either is answered by moving
  σ. An inaccurate solve is an operator that is not quite symmetric and not quite the
  inverse; the iteration is then slower or stalls, and the residuals on A, which do not
  use the solve, say what was reached.

## What a set of pairs proves

**One pair.** For any symmetric A, any vector y ≠ 0 and any number θ, an eigenvalue of A
lies within `‖A y − θ y‖ / ‖y‖` of θ: if θ is not an eigenvalue,
`‖y‖ ≤ ‖(A − θ I)⁻¹‖ ‖A y − θ y‖`, and the norm of that inverse is one over the distance
from θ to the nearest eigenvalue.

Two pairs with close values may be telling of one eigenvalue. That they are not takes
the vectors together.

**Several pairs.** Let Y be n by k with columns y₁ … y_k, Θ the diagonal matrix of any k
real numbers θ_i, and

    R = A Y − Y Θ,   ρ = ‖R‖_F,        G = Yᵀ Y − I,   γ = ‖G‖_F,        w = (max θ_i − min θ_i) / 2.

If `γ ≤ 1/4`, then **A has k eigenvalues, counted with multiplicity, of which the i-th
smallest is within `δ = 3 (ρ + γ w)` of the i-th smallest θ.**

The derivation, to be carried in full on the contract page:

1. Subtracting a multiple c of the identity from A and from Θ changes neither R nor the
   statement. Take c the midpoint of the θ, so that `‖Θ‖₂ = w`.
2. `S = (I + G)^{1/2}` exists, and `Q = Y S⁻¹` has orthonormal columns. Let `H = Qᵀ A Q`,
   with eigenvalues μ₁ ≤ … ≤ μ_k.
3. In an orthonormal basis that begins with the columns of Q, A is `[[H, Cᵀ], [C, W]]`.
   Setting C to zero moves no eigenvalue by more than `‖C‖₂` (Weyl), and leaves a matrix
   whose eigenvalues are those of H and those of W. So k eigenvalues of A, of distinct
   indices and in order, are each within `‖C‖₂` of a μ_i. This is Kahan's theorem (1967).
4. `‖C‖₂ = ‖A Q − Q H‖₂`, and H is the matrix X that makes `‖A Q − Q X‖₂` least, since
   `A Q − Q X` is `Q (H − X)` plus a part orthogonal to Q. With `X = S Θ S⁻¹`,
   `A Q − Q X = R S⁻¹`. So `‖C‖₂ ≤ ρ / √(1 − γ)`.
5. `M = Yᵀ A Y = (I + G) Θ + Yᵀ R`, and M is symmetric, so it is also the mean of that
   and its transpose: `M = Θ + (G Θ + Θ G)/2 + (Yᵀ R + Rᵀ Y)/2`. The i-th eigenvalue of M
   is therefore within `e = γ w + √(1 + γ) ρ` of the i-th smallest θ (Weyl again).
6. `H = S⁻¹ M S⁻¹`, and the i-th eigenvalue of such a product is the i-th of M times a
   number between the least and the greatest eigenvalue of `S⁻²` (Ostrowski, 1959), which
   lie between `1/(1 + γ)` and `1/(1 − γ)`. So `|μ_i − λ_i(M)| ≤ (w + e) γ / (1 − γ)`.
7. Adding the three distances with `γ ≤ 1/4` gives at most `2.66 ρ + 2.67 γ w`.

The constant is not sharp. With exactly orthonormal vectors the argument as written
gives `2 ρ`; one more step gives `√2 ρ`, since then `R = Q (H − Θ) + Q⊥ C` with the two
parts orthogonal, so `‖H − Θ‖_F² + ‖C‖_F² = ρ²` and the two distances add to at most `√2 ρ`.
With the vectors of this iteration γ is a small multiple of rounding, the third branch
of the orthogonalization seeing to it, and the term γ w is far below any residual.

**What δ is formed from.** ρ and γ are computed in floating point. On an operator their
rounding cannot be bounded, since nothing is known of the operator's own; there δ is a
theorem about the numbers formed and is offered as that. On a matrix stored as its lower
triangle the product's rounding is bounded entry by entry by the product of the
magnitudes, at the cost of a second product a vector, and the Gram matrix's by its own
count of operations; with both allowed for and rounded upward δ is an upper bound, and
that is the δ the count below may be compared with.

## What the inertia counts

`factor-signed!(symbolic, values, s)` leaves `P K̃ Pᵀ ≈ L D Lᵀ`, K̃ being `A − s I` with
its diagonal as rounded, and `factor-backward-bound` is a proved bound β on the 2-norm
of what the stored factor is out by. Let

    e = β + u · max_i |fl(a_ii − s)|,        u = 2⁻⁵³,

rounded upward, the second term being what the subtraction on the diagonal can have
lost. Then `L D Lᵀ = P (A − s I + E) Pᵀ` exactly, for a symmetric E with `‖E‖₂ ≤ e`. L is
unit triangular, so by Sylvester's law `A − s I + E` has exactly ν negative eigenvalues,
ν the count of negative pivots, and none that is zero. No eigenvalue of `A − s I + E` is
further than e from the eigenvalue of `A − s I` of the same index. So

    #{ eigenvalues of A below s − e }   ≤   ν   ≤   #{ eigenvalues of A below s + e },

both counted with multiplicity. **The count is of a band, not of a point**: of an
eigenvalue within e of s it says nothing, and e is a number the call returns.

For two numbers a < b with counts ν_a, ν_b and margins e_a, e_b, where
`a + e_a ≤ b − e_b`, subtracting the two statements:

    #{ eigenvalues of A in [a + e_a, b − e_b) }   ≤   ν_b − ν_a   ≤   #{ eigenvalues of A in [a − e_a, b + e_b) }.

The second inequality, which the first version of this page did not draw, is from the
other halves of the same two statements. The counts bound the inner interval from above
and the outer one from below.

## None was missed

Put the two together. Suppose k pairs stand, measured on the stored matrix, with a bound
δ; every interval `[θ_i − δ, θ_i + δ]` lies inside `[a + e_a, b − e_b)`; and
`k = ν_b − ν_a`. The pairs prove at least k eigenvalues of A in that interval and the
counts prove at most k. So **A has exactly k eigenvalues in `[a + e_a, b − e_b)`,
counted with multiplicity, and the i-th is within δ of the i-th θ.** Without a lower
end, `k = ν_b` gives the same of everything below `b − e_b`; without an upper end,
`k = n − ν_a` of everything from `a + e_a` up.

When k is less than `ν_b − ν_a`, something is proved, of the outer interval: if no
other pair's interval meets `[a − e_a, b + e_b)`, then **at least `ν_b − ν_a − k`
eigenvalues of A in that outer interval have no pair.** They are proved unfound; they are
not proved to lie in the inner interval, which the counts bound only from above. The
outcome is named for that, with the number, and the remedy is a second run from another
starting vector, kept orthogonal to the pairs that stand. When a pair's interval crosses
an end of the inner interval, or meets the outer one without lying in the inner, nothing
is said, and the end is moved.

## Interface

### `coalton-linalg/dense/tridiagonal`

- `eigen!(order, diagonal, off-diagonal, rows, c) -> TridiagonalStatus`. `diagonal` holds
  at least `order` numbers and `off-diagonal` at least `order − 1`. On success the first
  `order` of `diagonal` are the eigenvalues, ascending, and `off-diagonal` is destroyed.
  `c` is `rows` by `order`, column-major, and is overwritten by C Z with its columns in
  the eigenvalues' order; `rows` may be 0. No workspace: the call is in place and
  allocates nothing.

| TridiagonalStatus | Meaning |
| --- | --- |
| TridiagonalSuccess | The eigenvalues and C Z stand |
| TridiagonalNotConverged | The step bound was reached; the arrays are as the iteration left them and are not a result |
| TridiagonalNonFinite | An entry of the matrix or of C is not finite, nothing written; or an eigenvalue is too large for the format |
| TridiagonalBadSize | An array is too short for the order and the rows; nothing written |

### `coalton-linalg/lanczos`

- `make-lanczos(dimension, step-limit, pair-limit) -> LanczosWorkspace` allocates once:
  the basis, `dimension` by `step-limit + 1`; the vectors of the pairs, `dimension` by
  `pair-limit`; and the small problem, of the order of `step-limit²`. An entry budget
  refuses a workspace too large, by the error the other packages raise.
- `lanczos!(workspace, apply, start, wanted, which, tolerance) -> LanczosStatus`.
  - `apply` is a function of two arrays `(v, out)` that writes `B v` into out.
  - `start` is any vector that is not zero; it is read and not written. The library
    makes no random numbers. A start with nothing of a wanted eigenvector in it does
    not find that eigenvector.
  - `wanted` is how many pairs; `which` is `Largest`, `Smallest` or `LargestMagnitude`.
  - `tolerance` is absolute, on `‖B y − θ y‖₂` for each wanted pair with y a unit
    vector, formed from the operator. It may be zero.
- `lanczos-more!(workspace, apply, start, wanted, which, tolerance)`: the same, with the
  pairs that stand kept, the start and every vector of the run made orthogonal to their
  vectors, and the pairs it finds added to them. It is how a second copy of a multiple
  eigenvalue is found. It keeps every pair that stands, whatever its residual:
  `drop-above!(workspace, bound)` first removes those whose residual norm, as last
  formed, is above a bound the caller gives, and the bound to give is the tolerance about
  to be asked, since the run's estimate cannot fall below the residuals of the pairs
  kept. The values it adds are the operator's it ran on; where the pairs kept were
  measured on another, `measure!` afterwards makes them all one operator's. **The δ that
  stood is withdrawn when a pair is added**: it was of a subset. After `drop-above!` it
  stands, being a bound for any subset of the pairs it was formed for.
- `measure!(workspace, apply) -> MeasureStatus`: for the vectors that stand and this
  operator, which may be another than the one iterated on, replaces each value by the
  Rayleigh quotient and each residual norm by the one formed, one product a vector, and
  forms γ and δ. `measure-symmetric!(workspace, pattern, values)`: the same on a matrix
  stored as its lower triangle, with the rounding of its own arithmetic allowed for, so
  that the residual norms and δ are upper bounds.
- `pair-count`, `pair-values`, `pair-vectors`, `residual-norms`: the pairs that stand,
  ascending by value, the vectors as columns.
- `orthogonality-defect`, `error-bound`: γ and δ as last formed; `error-bound` is −1
  where γ is above 1/4, where nothing has been measured, and where a pair has been added
  since it was formed.
- `step-count`, `confirmation-count`, `estimates`: of the latest run; the estimates are
  `β_j |s_j|` for the wanted pairs at its last step, as they were.
- `asymmetry-estimate`: over the latest run, the largest of `|q_{j−1}ᵀ B q_j − β_{j−1}|`
  and of `|q_iᵀ B q_j|` for every i below j − 1, all of which the first pass computes
  anyway. For a symmetric operator and an orthonormal basis each is rounding. Far above
  rounding it says the operator is not symmetric on the vectors met; small, it proves
  nothing.
- `certify-between(workspace, a, count-a, margin-a, b, count-b, margin-b)`,
  `certify-below(workspace, b, count-b, margin-b)` and
  `certify-above(workspace, a, count-a, margin-a) -> CompletenessStatus`: the comparison
  under *None was missed*, on the pairs and the δ that stand, with every sum and
  difference rounded against the statement. Arithmetic only. `unfound-count` is
  `ν_b − ν_a − k` where that is positive: a proved lower bound on the eigenvalues of the
  outer interval that have no pair.

| LanczosStatus | Meaning |
| --- | --- |
| LanczosConverged | Every wanted pair has a residual, formed from the operator, within the tolerance |
| LanczosEnded | The sequence ended, by a remainder of zero or by the third branch of the orthogonalization, so that no more can come from this start, with fewer pairs than wanted or a residual above the tolerance |
| LanczosSmallProblem | The small problem did not converge within its bound, or an eigenvalue of it is too large for the format: the pairs that stood before the call stand and none is added |
| LanczosLimit | The step limit was reached with a wanted residual above the tolerance |
| LanczosNonFinite | The start, an image or an update was not finite |
| LanczosOutOfRange | An image is too large, in the units of the first, to be squared |
| LanczosBadStart | The start is zero, or is in the span of the pairs that stand as far as the arithmetic can tell, by the same three branches |
| LanczosBadCount | `wanted` is zero, or with the pairs that stand exceeds the pair limit or the dimension, or exceeds the step limit |
| LanczosBadSize | The start does not have the workspace's dimension |
| LanczosBadTolerance | The tolerance is negative or not a number |

After the first three the pairs found stand, each with its formed residual, and the
status says only whether all are within the tolerance. The four `Bad` statuses are
returned with nothing written. After `LanczosSmallProblem`, `LanczosNonFinite` and
`LanczosOutOfRange` the pairs that stood before the call stand and none is added.

| CompletenessStatus | Meaning |
| --- | --- |
| Complete | The interval holds exactly the pairs inside it, each eigenvalue within δ of its value |
| Unfound | Fewer pairs are inside than the counts give: at least `unfound-count` eigenvalues of the outer interval have no pair. Not proved to lie in the inner interval |
| Straddles | The interval of a pair crosses an end of the inner interval, or meets the outer one without lying in the inner; nothing is said |
| NoBound | No δ stands |
| Inconsistent | More pairs are inside than the counts allow, which the two theorems exclude: the matrix measured is not the one counted, or the δ that stands is not an upper bound |
| BadInput | A number is not finite, a margin is negative, the ends are out of order, or the count at the upper end is below the count at the lower |

### `coalton-linalg/sparse/inertia-count`

- `make-count() -> Count`; `count-below!(count, symbolic, values, shift) -> CountStatus`
  factors `A − shift I` by `factor-signed!`, which **replaces the factor that stood on
  that symbolic analysis**, and forms the margin. `below`, `margin`: ν and e. The
  statement is the one under *What the inertia counts*, of the matrix as stored.

| CountStatus | Meaning |
| --- | --- |
| CountSuccess | `below` and `margin` stand |
| CountZeroPivot | A pivot was exactly zero; move the shift |
| CountUnbounded | The factorization's bound is not finite, so no margin is |
| CountNonFinite, CountBadValues, CountBadShift | As the factorization's |

The factor left by `count-below!` at a is also the operator for shift and invert at
σ = a, so an interval certified at both ends costs two factorizations.

## Decisions

**D1. The basis is kept and every vector made orthogonal to it.** The alternatives keep
less and pay in copies of eigenvalues that have to be recognized afterwards. A count
compared with pairs needs the pairs to be distinct by construction, and the bound δ
needs γ small. Selective and partial schemes exist and are left out.

**D2. No restart inside a call.** A run that reaches its step limit returns what it has.
With shift and invert a run is expected to take tens of steps where the eigenvalues
wanted are not closer to one another than to the shift, which is not measured here; one
that needs more wants its shift nearer, or has such a pair. A caller continues with
`drop-above!` and `lanczos-more!`, from a start of its own.

**D3. Confirmation instead of trust**, as in [least squares on an operator](lsqr.md).
When every wanted estimate is within the tolerance the vectors are formed and their
residuals taken from the operator, a product each. If one is above the tolerance the
iteration goes on from where it was, and the trigger falls by the factor the worst
estimate was out by; a trigger that would be zero is switched off. At the limit and at
the end of the sequence the residuals are formed whatever the estimates say.

**D4. Absolute tolerances, the caller's**, on the operator iterated on. No relative test
against a norm estimated on the way.

**D5. The starting vector is the caller's.** A library that made one would have to make
random numbers or choose a fixed vector that some matrix defeats. The contract would
say what a good one is.

**D6. The sequence ends on a remainder of zero, or on one that two passes of
orthogonalization have each more than halved.** No threshold beside the operator's size
decides that a small norm is zero; the test is of what the passes removed. The first
version had an exact zero alone, and *Where the sequence ends* says what that would have
returned.

**D7. Everything iterated is near 1, where the first image has the operator's size.**
Where it has not, an image far above the first is `LanczosOutOfRange`. The start is divided by the power of two at its
largest entry before it is normalized, and every image by the power of two at the
largest entry of the first. T is then in those units and the values are multiplied back
once. The same operator multiplied by a power of two gives the same vectors, bit for
bit, and values scaled exactly, wherever its own products stayed in range. The small
problem scales its own matrix the same way. A first image of zero ends the sequence at
its first step, with the eigenvalue zero.

**D8. No squared norm is formed** where its root would do: the rotations' lengths as in
least squares on an operator, a norm of a vector that is all tiny taken again with its
entries lifted.

**D9. The count is its own call and its own package.** It needs a stored matrix and a
factorization; the iteration needs neither. A caller with an operator alone has pairs
and residuals and no count, and is told so.

**D10. Proved for a stored matrix, measured for an operator.** See *What δ is formed
from*. `certify-between` and its two companions take whatever δ stands; the contract
would say that only the δ of `measure-symmetric!` makes their statement a proof.

## What would be claimed, and what not

- **Each residual reported is of the vector returned**, formed from the operator, to the
  rounding of one product; from `measure-symmetric!` it is an upper bound.
- **`Complete` is a theorem** when δ is from `measure-symmetric!` and the counts from
  `count-below!` on the same stored values: of the matrix as stored, in exact
  arithmetic, with every rounding on the way allowed for.
- **`Unfound` proves pairs missing from the outer interval, and does not prove them in
  the inner one.** No number of runs proves that none is missing. Only `Complete` says
  none is.
- **No rate and no step count is promised.** Eigenvalues inside the spectrum, or close
  together against its width, converge slowly without shift and invert.
- **One copy of a multiple eigenvalue a run**, in exact arithmetic. A second copy within
  one run is rounding's doing and not promised.
- **A vector is not determined** where eigenvalues are equal or closer than the
  residual. What is claimed of vectors is their residuals and their Gram matrix, never
  the vectors.
- **An operator that is not symmetric is not detected**, though `asymmetry-estimate` may
  show it. The theorems are of symmetric matrices.
- **The small problem's eigenvalues are accurate against the size of its matrix**, not
  each against its own size: a tridiagonal matrix with entries of very different sizes
  can have small eigenvalues this method does not resolve relatively.

## Verification planned

For the small problem:

1. **Every eigenvalue bracketed by an exact count**: the number of eigenvalues of the
   stored tridiagonal matrix below a rational number, by the signs of an exact
   elimination, on either side of each computed value at a distance of a small multiple
   of rounding times the matrix's size.
2. **`T Z − Z Λ` and `Zᵀ Z − I`** exactly, in rational arithmetic, for C the identity; C
   a single row and C a general matrix against that Z.
3. **Matrices made for it**: one that splits, by a zero and by an entry at the threshold
   on either side; equal diagonals, the two-by-two `[[a, e], [e, a]]` among them, which a
   shift without its tie rule leaves standing; a zero diagonal with an off-diagonal entry
   that is not zero; Wilkinson's W₂₁⁺, whose top pairs agree to many
   digits; graded entries; order 1 and 2; a step bound of one, for the named outcome.

For the iteration:

4. **Against exact arithmetic.** Small integer matrices with known spectra, among them
   multiple eigenvalues, a zero eigenvalue and an indefinite spectrum: each value
   reported is the exact Rayleigh quotient of its vector and each residual norm the
   exact one, to rounding; each interval `θ ± δ` holds an exact eigenvalue as claimed,
   the claim of *Several pairs* checked against an exact count on both sides.
5. **An independent recomputation** in the test language, with the library's product,
   reproduces T, the estimates, the steps, the confirmations and the pairs bit for bit.
6. **The orthogonality kept**, and what one pass alone would have left on a case made to
   show it. **A run past what the start can reach**: the reviewer's four by four matrix
   with a step limit of seven, a block-diagonal operator with the start inside one block,
   and a second run whose start is in the span of the pairs that stand. Each must end by
   name, with no value the operator does not have. On small matrices a step limit at or
   above the reachable dimension is the ordinary case, and every exact test meets it.
7. **A multiple eigenvalue**: one copy from a run, the second from `lanczos-more!`, and
   the count naming the one unfound in between. **The second term of the second run's
   estimate**: pairs kept with residuals above the tolerance asked, where the estimate
   must not trigger, and the same after `drop-above!`, where it must.
8. **Every status and every completeness outcome**, each by a case made for it, with
   nothing written where that is said.
9. **A confirmation that fails**, by a product wrong in one step, as in least squares on
   an operator.

For the count and the statement:

10. **The band is true and is needed.** Exact counts of exact eigenvalues, by rational
    arithmetic, on both sides of every margin; and a matrix with an eigenvalue at the
    shift to rounding, on which the pivots' signs alone would miscount.
11. **`Complete` is never issued for a false statement**: on every test matrix the exact
    number of eigenvalues in the certified interval equals the pairs', including cases
    with one withheld, which must come out `PossiblyMissing`.
12. **Shift and invert end to end** on a sparse matrix of known spectrum, in the
    examples: the eigenvalues nearest a shift, measured on the matrix, certified between
    two counts.

For all of it: scale, with entries in [the scale-sweep roster](rosters.md); aliasing, in
[the aliasing roster](rosters.md); zero bytes a call after warm-up through compiled
callers; planted source defects, each applied alone; and one timed problem that is not
small, with the time of each part.

## Not in this design

- The generalized problem `A x = λ M x`, which is the form vibration and buckling take.
  It is this iteration in the inner product of M, and is the natural next step; it
  doubles what has to be said about orthogonality and is left out to say this first.
- Restarting with the basis compressed; a block of starting vectors; orthogonalizing
  only against the pairs that have converged.
- The reduction of a dense symmetric matrix to tridiagonal form.
- Eigenvalues of the small problem to high relative accuracy, and counts on it by
  bisection.
- Singular values of a sparse matrix, which are this iteration on a product or the
  bidiagonalization of [least squares on an operator](lsqr.md) with its basis kept.
- Matrices that are not symmetric.

## Sources

The iteration is that of Lanczos (1950); the account of how it loses orthogonality is
Paige's (1971); the rule that two passes of orthogonalization suffice is Kahan's, as
given in Parlett's *The Symmetric Eigenvalue Problem*, which also gives Kahan's theorem
on a set of orthonormal vectors and their residuals. Shift and invert for this iteration
is that of Ericsson and Ruhe (1980). The shifted QR step is that of Francis (1961), with
Wilkinson's shift and his proof of convergence (1968); Trefethen and Bau, *Numerical
Linear Algebra*, lectures 28 to 30 and 36, describe the iteration and the shifted QR
algorithm in its explicit form; the implicit form used here, and the conditions under
which an entry is taken to be zero, are not in those lectures and are set out above. The
theorems of Weyl, Sylvester and Ostrowski are used as stated. **Cited from memory, with
no text at hand that carries them**, and so to be checked before they are relied on for
more than a name: Ericsson and Ruhe (1980), Kahan's theorem (1967) and its place in
Parlett's book, Kahan's rule for the passes as given there, Wilkinson's proof of
convergence, and Ostrowski (1959). The reviewer could confirm Lanczos, Paige, Francis,
Wilkinson's shift and its convergence, Parlett's book and Sylvester's theorem from the
texts available, and recalls the others the same way. The extension of the bound to vectors that
are only nearly orthonormal, the band about a count, the comparison of the two and the
rules of confirmation, scaling and ending are worked out on this page, and the code
would be written from it. No implementation and no listing would be opened.
