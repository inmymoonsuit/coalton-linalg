# Cholesky-form factorization, for a proof

Package: `coalton-linalg/sparse/spd-ldlt`. Source: [spd-ldlt.ct](../src/sparse/spd-ldlt.ct).
State: completed as Batch 6a, with its own solve added in Batch 10; verified at safety 1
and safety 0. The batch records, [6a](batch6a.md) and [10](batch10.md), hold the tests and
the controls.

## What it is for

A floating-point Cholesky factorization that runs to completion on a suitably shifted
matrix proves the unshifted matrix positive definite: that is the method of verified
positive definiteness (Rump, *Verification of positive definiteness*, BIT 46, 2006). The
theorems are stated for Cholesky's own arithmetic, in which each entry of the factor is one
value minus a sum of products of two numbers, divided once or square-rooted once.

`factor!` computes L D Lᵀ, and each of its products multiplies three numbers. This entry
performs Cholesky's arithmetic, over the same analysis and in the same elimination order,
so that such a theorem applies to what the code does. **The library states the arithmetic.
Choosing the shift and checking a theorem's hypotheses is the caller's.**

## Interface

- `factor-cholesky!(symbolic, values, shift) -> FactorStatus` computes the lower triangular
  G with `P (A − shift I) Pᵀ = G Gᵀ` and returns one of the outcomes of
  [`factor!`](spd-ldlt.md): `FactorSuccess`; `FactorNonPositive` when a pivot is not
  strictly positive, with `refused-column` naming it; `FactorNonFinite`; `FactorBadValues`;
  `FactorBadShift`.
- `cholesky-ready?(factor) -> Boolean` says whether the storage holds the G of the latest
  successful call.
- `solve-cholesky!(factor, b, out) -> SolveStatus` solves **`(A − shift I) x = b`** with
  that G, the shift being the one the factorization was given. It has the outcomes and
  the rules of [`solve!`](spd-ldlt.md): `SolveInvalidFactor` when no Cholesky-form factor
  stands, `SolveBadSize`, `SolveNonFiniteRhs`, `SolveNonFiniteResult`; b and out may be the
  same array; on any refusal out is unchanged.

**It is a factor of the shifted matrix, and it has its own solve.** A Symbolic owns one
numeric factor, and this call overwrites it. Afterwards `factor-ready?` is false, `solve!`
returns `SolveInvalidFactor`, and `estimate-inverse-norm1!` returns
`ConditionInvalidFactor`, until `factor!` runs again; `solve-cholesky!` in turn refuses a
solving factor. A refused call leaves nothing standing. A successful factorization or
solve allocates nothing after warm-up on the measured runtime.

**The solve does not change the proof.** `factor-cholesky!` is what it was, operation for
operation. The solve reads G and writes only its own scratch and the caller's output. A
caller who has made the proof can therefore solve with the same factor, and where the
unshifted solution is wanted, [refine](refinement.md) against the unshifted matrix:
one factorization serves both.

G is read through the existing views, in elimination order: `factor-diagonal` holds the
diagonal `g[k,k]`, `factor-values` the entries below it, column by column as
`factor-columns` and `factor-rows` lay them out, and `factor-permutation[k]` is the
original index at position k.

## The arithmetic, step by step

All values are IEEE doubles under round-to-nearest, each operation rounded once: no fused
multiply-add, no extended precision, no reordering. Write fl for one rounded operation and
work in elimination order. For column k from first to last:

1. w takes the stored entries of column k of the permuted matrix, at rows k and below; the
   positions of the factor's structure that the matrix does not store take 0.
2. `w[k] ← fl(w[k] − shift)`.
3. For each earlier column j that has an entry in row k, **in increasing j**, and for each
   row i of column j from row k downward, in increasing i:
   `w[i] ← fl(w[i] − fl(g[k,j] · g[i,j]))`. The term for i = k is `g[k,j]²`.
4. If `w[k]` is not finite the call returns `FactorNonFinite`; if it is not greater than 0
   it returns `FactorNonPositive`.
5. `g[k,k] ← fl(√w[k])`.
6. For each row i of column k's structure below k: `g[i,k] ← fl(w[i] / g[k,k])`; a
   non-finite quotient returns `FactorNonFinite`.

So every entry of G is one stored value, or zero, less a chain of rounded products
subtracted one at a time in a fixed order, then divided once or square-rooted once. A
position that is in the structure but numerically zero takes part as an exact zero.
Underflow follows the arithmetic and is not diagnosed; overflow ends in a refusal.

## The counts a proof uses

For entry (i, j) with i ≥ j, the number of products subtracted is the number of earlier
columns k < j in which both row i and row j have a structural entry.

**t′ⱼ, the count in Rump's bound, is the number of entries of G in row j to the left of the
diagonal.** It is a row count, not the length of column j: `factor-columns[j+1] −
factor-columns[j]` counts the entries below the diagonal in column j, which is a different
number. No entry receives more products than t′ⱼ of its column index j. A caller tallies
the row counts from the row indices of the stored strict lower triangle:

```lisp
(let counts = (arr:make n (the UFix 0)))
(let rows = (ldlt:factor-rows factor))
(for ((q 0 (+ q 1))) :while (< q (arr:length rows))
  (let j = (arr:aref rows q))
  (arr:set! counts j (+ (arr:aref counts j) 1)))
;; counts[j] is t'_j for elimination position j.
```

The indices are elimination positions; the structure and the counts are fixed by the
analysis and are the same for `factor!`.

## What is tested, and what is not claimed

For each stored entry, with p its number of products, the standard rounding analysis of
this arithmetic gives

    | ã[i,j] − Σₖ g[i,k] g[j,k] |  ≤  γ₍p+1₎ Σₖ |g[i,k]| |g[j,k]|   off the diagonal,
    | ã[j,j] − Σₖ g[j,k]² |        ≤  γ₍p+2₎ Σₖ g[j,k]²             on it,

with `γₘ = m u / (1 − m u)`, u = 2⁻⁵³, and ã the permuted matrix with its diagonal as
rounded in step 2. The tests compute the left side exactly, in rational arithmetic on the
stored values, for every stored entry of every test matrix, and assert it against
`(p + 2) u` times the right side's sum. That is a check that the code performs the
arithmetic the analysis assumes. It is not a proof of the analysis, which is textbook
material, and it is not a proof that any particular matrix is positive definite.

Step 2 rounds once. Where the caller's theorem is about `A − shift I` itself and not about
the rounded diagonal, that rounding is the caller's to account for, or the caller passes a
shift of 0 with the diagonal already formed.

## The solve

`P (A − shift I) Pᵀ = G Gᵀ` to rounding, so the solve is two triangular substitutions
between a gather and a scatter by P: `G z = P b` a column at a time, each diagonal
dividing before its column is taken off the rest, then `Gᵀ y = z` from the last row up,
and `x = Pᵀ y`. With w the longest row or column of G off the diagonal, the computed x
satisfies, row by row and barring underflow,

    | (A − shift I) x − b |  ≤  γ₍3w+5₎ |G| |Gᵀ| |x|     (rows and columns as permuted),

the factorization's bound and one triangular solve's on each side. The tests check that
row by row in rational arithmetic.

## Representation choices

The entry reuses the factor's storage and views instead of adding a second factor to every
analysis: one state value distinguishes the three conditions of the storage, nothing, a
solving factor and a proof factor, and each solve accepts only its own.
