# A factorization that leaves small pivots out

Package: `coalton-linalg/sparse/spd-ldlt`. Source: [spd-ldlt.ct](../src/sparse/spd-ldlt.ct).
State: completed as Batch 11; verified at safety 1 and safety 0. The [batch
record](batch11.md) holds the controls.

## What it is for

A symmetric positive semidefinite matrix has no L D Lᵀ factorization that `factor!` will
give: the first dependent index has a zero pivot and the call refuses it. Often the
dependent indices are the answer wanted. In `A Aᵀ` they are the dependent rows of A; in a
stiffness matrix they are the directions a structure is free to move in; in a Gram matrix
of functions they are the redundant ones.

`factor-skipping!` factors in the analyzed order and **leaves out every index whose pivot
is at or below a floor**. What stands afterwards is the ordinary factor of the principal
submatrix of the indices kept, in the storage and the structure already analyzed. Two
calls read it: `solve-skipping!` solves with the kept submatrix, and
`skipped-null-vector!` gives, for an index left out, its combination of the kept indices
before it.

[Dependent rows](dependent-rows.md) is this, applied to `A Aᵀ`, with the bounds that turn
what it finds into something proved.

## Interface

- `factor-skipping!(symbolic, values, shift, floor) -> SkipFactorStatus`. The values are
  those of `factor!`, the lower triangle in the analyzed order, and the shift is
  subtracted from the diagonal first. A computed pivot `p` with `p <= floor` leaves its
  index out; any other is kept and is positive, since the floor is not negative.
- `skipping-ready?(factor)`: whether this kind of factor stands.
- `skipped-count(factor)`: how many indices it left out. `positive-pivot-count` is the
  number kept and `negative-pivot-count` is zero.
- `factor-diagonal(factor)` holds **exactly zero at each index left out**, in elimination
  order, and the pivot at each kept one; with `factor-permutation` that names them.
- `solve-skipping!(factor, b, out) -> SolveStatus`: with S the kept indices and M the
  matrix factored, solves `M[S,S] x[S] = b[S]` and writes zero at the indices left out.
  The entries of b at those indices do not enter the result, though every entry of b must
  be finite. `out` may be `b`.
- `skipped-null-vector!(symbolic, index, out) -> NullVectorStatus`: for an index left out,
  in original numbering, writes the vector v with

  - `v[index] = 1`,
  - `v = 0` at every other index left out and at every kept index eliminated after it,
  - `(M v)[j] = 0`, to rounding, at every kept index j eliminated before it.

  Then `vᵀ M v` is the pivot that was judged, the Schur complement of the index on the
  kept indices before it.

| SkipFactorStatus | Meaning |
| --- | --- |
| SkipSuccess | A factor of the kept principal submatrix stands |
| SkipNonFinite | A value or an intermediate result was not finite; `refused-column` names the column |
| SkipBadValues | The value count differs from the analyzed pattern's |
| SkipBadShift | The shift was not finite |
| SkipBadFloor | The floor was negative, infinite or not a number |

| NullVectorStatus | Meaning |
| --- | --- |
| NullVectorSuccess | The output holds the vector |
| NullVectorNoFactor | No factor of this kind stands |
| NullVectorNotSkipped | The index is out of range, or was kept |
| NullVectorBadSize | The output's length differs from the order |
| NullVectorNonFinite | An entry of the vector was not finite; the output is unchanged |

`solve-skipping!` returns the outcomes of `solve!`, with `SolveInvalidFactor` when no
factor of this kind stands. Every refusal leaves its output untouched. A refused
`factor-skipping!` leaves no factor of any kind standing. None of these calls allocates.

**The kinds of factor do not mix.** A factor from `factor-skipping!` is refused by
`solve!`, by `solve-cholesky!` and by `estimate-inverse-norm1!`; a factor from `factor!`
or `factor-signed!` is refused by `solve-skipping!` and `skipped-null-vector!`. Each
factorization call on an analysis replaces the one factor it owns.

## What is and is not claimed

- **The kept factor is an ordinary one.** For kept i and j the computed factor satisfies
  the entrywise bound of [the signed factorization](quasidefinite.md#accuracy-growth-and-the-elimination-order),
  `|M[i,j] − (L D Lᵀ)[i,j]| ≤ γ₍p+2₎ (|L| |D| |Lᵀ|)[i,j]`, with p the products that entry
  received from kept columns, and a solve satisfies the solve bound on the kept rows. Both
  are checked exactly in the tests. With nothing at or below the floor the factor is the
  one `factor!` computes, bit for bit.
- **The floor is compared with a computed pivot.** A pivot is a difference: the diagonal
  entry less what the kept columns before it account for. Its rounding error is of the
  order of the unit roundoff times the diagonal entry, however small the exact pivot is.
  **A floor below that level decides on noise.** For a matrix with diagonal entries near
  1, floors below about 10⁻¹⁵ are not meaningful, and a caller separating "zero" from "not
  zero" should set the floor well above the noise and well below the pivots it means to
  keep. Where the matrix is a Gram matrix `A Aᵀ` the pivot is a *squared* distance, which
  is what limits the method there; [dependent rows](dependent-rows.md#the-floor-of-the-method)
  states the consequence.
- **Which indices are left out depends on the order.** The factorization meets the indices
  in the fill-reducing order and keeps each one that is independent of those kept before
  it. Another order keeps another, equally valid, set. The number kept is the same
  wherever the decisions are not near the floor.
- **No pivoting for rank.** A kept set can be badly conditioned though every kept pivot
  cleared the floor: the pivots measure each index against those before it, not the set
  as a whole. A caller that needs the kept submatrix to be well conditioned checks it, by
  [`factor-cholesky!`](cholesky-form.md) at a shift on a matrix of the kept indices or by
  the route [dependent rows](dependent-rows.md#proving-the-kept-rows-independent) gives.
- **A negative pivot is left out, not refused.** It is at or below every floor. For a
  matrix that is not positive semidefinite the result is still the factor of the kept
  principal submatrix, which is then positive definite, and it says nothing more about
  the matrix.

## The null vector

Let k be an index left out and S the kept indices eliminated before it. The factor of
`M[S,S]` is `L[S,S] D[S] L[S,S]ᵀ`, and row k of L, computed before the pivot of k was
judged and left in place, satisfies `M[S,k] = L[S,S] D[S] L[k,S]ᵀ`. So the y with
`M[S,S] y = M[S,k]` is the solution of

    L[S,S]ᵀ y = L[k,S]ᵀ,

one substitution from the last kept index up, and the vector is `e_k − y`. A column left
out is empty, so its entry in row k is zero and it drops out of the substitution by
itself.

Its error: writing B for `|L| |D| |Lᵀ|` and w for the longest row or column of L,

    | (M v)[j] |  ≤  γ₍2w+4₎ · Σ_l B[j,l] |v[l]|        for kept j before k,

which follows from the factorization bound on `M[S,S]`, the same bound on the entries
`M[S,k]` that row k accounts for, and the substitution's own error. The tests check this
inequality exactly on every vector of the rosters; the largest observed was 0.061 of it.

A vector computed this way inherits the conditioning of `M[S,S]`. Where M is `A Aᵀ` that
is the conditioning of the kept rows *squared*, and [dependent rows](dependent-rows.md)
refines the vector against A itself.

## Representation choices

The kernel is the one `factor!` and `factor-signed!` use, with one more branch at each
pivot; a kept index is treated exactly as `factor!` treats it. An index left out is
represented by a zero pivot and an empty column and not by a separate list, so every
later column's update from it is a multiplication by zero that the existing loop already
performs, and the solve needs only to clear its entry. The row of L it leaves behind is
found through the row incidence the analysis already holds, which is why
`skipped-null-vector!` takes the analysis and not the factor.

The index is located by a pass over the permutation, linear in the order. The vectors
asked for are as many as the indices left out, which is few where the method applies.
