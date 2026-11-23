# Dependent rows of a sparse matrix

Package: `coalton-linalg/sparse/dependent-rows`. Source:
[dependent-rows.ct](../src/sparse/dependent-rows.ct).
State: completed as Batch 11; verified at safety 1 and safety 0. The [batch
record](batch11.md) holds the controls, and [the design](design-dependent-rows.md) the
reasoning that came before the code.

## What it is for

The rows of A are constraints, measurements or equations, and some of them say nothing
the others do not. This module answers three questions about them:

1. **Which rows to keep.** A set S of rows, each independent of those kept before it, and
   the rows left out.
2. **What each row left out depends on.** A vector v for it, with 1 at that row, the
   negated multipliers of its combination at the kept rows, and zeros at the other rows
   left out, so that `vᵀ A` is the row less its combination: a left null vector of A, to
   the accuracy the next item states.
3. **What of that is proved.** An upper bound on `‖Aᵀ v‖₂` that holds for the exact product
   of the stored values, and from the bounds together an upper bound on the singular
   values of A beyond the number of rows kept.

**The finding and the proving are separate.** The rows are found by a factorization of
`A Aᵀ` that leaves small pivots out, which is fast and limited in what it can resolve. The
bounds are formed afterwards from A itself and hold whatever the search did, so a search
that chose badly is seen to have.

## Interface

- `analyze(matrix) -> RowAnalysis`, `analyze-with-limit(matrix, maximum-entries)`: once for
  a pattern. It plans `A Aᵀ`, analyzes its factorization and allocates every workspace.
  The errors are those of [the Gram analysis](gram.md) and [the factor
  analysis](spd-ldlt.md); the budget applies to each.
- `find-dependent-rows!(analysis, matrix, threshold) -> RowStatus`. The matrix must have
  the analyzed pattern; another pattern or shape raises the typed error of
  `refill-csc!`, and a value that is not finite that of `refill!`. The threshold is a distance, in the
  units of the entries: a row is left out when its squared distance from the span of the
  rows kept before it, as computed, is at most the threshold squared.
- `rows-found?(analysis)`, `kept-count(analysis)`, `left-out-count(analysis)`.
- `left-out-rows(analysis)`: a borrowed array as long as the row count, whose first
  `left-out-count` entries are the rows left out. **A position in it** is what the next
  two calls take.
- `null-vector!(analysis, position, out) -> NullStatus`: writes v for the row at that
  position and records its bound.
- `residual-bound(analysis, position)`: the bound that call recorded on `‖Aᵀ v‖₂`; −1 when
  there is none since the latest search.
- `gap-bound(analysis)`: with r rows kept, an upper bound on singular value r + 1 of A.
  Zero when no row is left out; −1 until every position has a bound.
- `kept-values!(analysis, out) -> Boolean`: the matrix's values as the latest search read
  them, in its stored order, with zeros for the entries of the rows left out. It is what
  [the proof below](#proving-the-kept-rows-independent) is given.
- `transpose-product-bound(matrix, v) -> F64`: the bound by itself, formed on the matrix and
  the vector each divided by the power of two at its largest entry, for any stored matrix
  and any vector as long as it has rows. A caller that zeroes small multipliers and wants
  to know what the reduced combination leaves asks this.
- `row-count`, `column-count`, `refinement-step-limit`, `default-analysis-entry-limit`.

| RowStatus | Meaning |
| --- | --- |
| RowsFound | A kept set stands |
| RowsNonFinite | An intermediate result of the factorization was not finite; no set stands |
| RowsUnderflow | Under a threshold of zero, a row holding a nonzero entry is so small beside the matrix's largest entry that its square cannot be held; no set stands |
| RowsBadThreshold | The threshold was negative or not finite, or positive and so small beside the matrix's largest entry that its square cannot be held; no set stands |

| NullStatus | Meaning |
| --- | --- |
| NullFound | The output holds the vector, and its bound is recorded |
| NullNoRows | No kept set stands |
| NullBadPosition | The position is not below `left-out-count` |
| NullBadSize | The output's length differs from the row count |
| NullNonFinite | The vector or its bound was not finite; the output is unchanged |

Every search withdraws the earlier set and its bounds before it does anything else, so a
refused search leaves none. After the analysis nothing allocates: a compiled caller that
searches, takes every vector and reads the gap bound is measured at zero bytes.
`transpose-product-bound` is the exception, a convenience that returns a boxed number.
One analysis serves one call at a time.

## What is proved

Write `A_S` for the kept rows, r for their number, and `v_i` for the vector of the i-th
row left out, with bound `ε_i`.

**Each row left out is within ε_i of the span of the kept rows.** `‖Aᵀ v_i‖₂ ≤ ε_i`, and
`Aᵀ v_i` is that row less a combination of kept rows, since `v_i` is zero at the other
rows left out. The bound is on the exact product: for a column with c stored entries the
computed sum `d̂` of the rounded products and the computed sum `ŝ` of their magnitudes give

    | (Aᵀ v)[j] |  ≤  |d̂| + γ₍2c+2₎ ŝ,

because the exact entry is within `γ_c` of the exact sum of those magnitudes, which is
within `γ₍c−1₎` of the computed one. The squares are summed, one smallest positive number
is added for each column to cover what underflow can remove from a product or a square,
and the root is multiplied by `1 + 2u(n + 10)` for the rounding of that arithmetic
itself. A sum that is not finite gives infinity, never a number that is not one. The
tests compare the value returned with the exact norm in rational arithmetic on every
vector, with cases in which the plain computed norm is zero and the exact one is not, in
which a product underflows, and in which every addition of the sum rounds down.

**A has at most r singular values above the gap bound.** Replace each row left out by its
combination of kept rows. The result has rank at most r, and differs from A only in those
rows, by `Aᵀ v_i` in each: a Frobenius distance of `√(Σ ε_i²)`. So

    σ₍r+1₎(A)  ≤  √(Σ ε_i²)  ≤  gap-bound.

This holds whether or not the kept rows are independent.

**That the kept rows are independent is proved separately**, as the next section sets
out, and with it

    σ_r(A)  ≥  σ_min(A_S)  >  τ,

since adding rows cannot lower the singular values a matrix already has. Where the gap
bound is below τ the matrix has **numerical rank exactly r, with a proved gap**: r
singular values above τ and the rest at or below the gap bound.

**What is not proved**: that the kept set is the best conditioned one, that a row's
combination uses the fewest rows possible, or anything about a row whose distance lies
below what the search can resolve, which is the subject of [its own
section](#the-floor-of-the-method).

## Proving the kept rows independent

The search keeps a row when it is far enough from the rows *kept before it*. That does
not make the kept set well conditioned: every row of the set can clear the threshold
while the smallest singular value of the set is far below it. The tests hold such a
matrix, ten rows that are the columns of a Kahan matrix, each at distance at least 0.13
from the rows before it and with a smallest singular value near 0.003; at a threshold of
0.05 every row is kept.

The proof is [the rank certificate](quasidefinite.md#a-rank-certificate-from-the-inertia)
of the signed factorization, and **it needs no second analysis**:

1. Analyze the augmented matrix of A once, with `analyze-quasidefinite`, as for any use
   of that certificate.
2. After a search, `kept-values!` gives A's values with zeros for the rows left out.
   Refill the augmented plan with those values and with the certificate's two diagonal
   scalars.
3. Run the certificate's steps as written there, with one change: **the pivot counts to
   require are n + (m − r) positive and r negative**, m being the rows of A and n its
   columns.

Why that is the same proof. With those rows zero the matrix `[[t I, Aᵀ], [A, s I]]` is,
after a symmetric permutation, the matrix `[[t I, A_Sᵀ], [A_S, s I]]` of the kept rows
beside a block `s I` for the rows left out. That block contributes m − r positive
eigenvalues whatever A is, so the whole has r negative eigenvalues exactly when the kept
rows' matrix has, which is when every singular value of `A_S` exceeds `√(s t)`. The
argument from the margins and the backward bound is unchanged.

A certificate withheld is not a disproof. Its outcomes are the certificate's: on the
tests' roster it was issued for every true statement, 253 times, 177 of them with rows
left out, and never for a false one.

## The floor of the method

In exact arithmetic the pivot of a row in the factorization of `A Aᵀ` is the *square* of
its distance from the rows kept before it. In floating point that square is a difference
of numbers the size of the row's squared norm, each rounded, so it is lost below about
the unit roundoff times that squared norm:

**a distance below about 1.5 × 10⁻⁸ of the row's size cannot be told from zero.**

A row whose distance lies between the caller's threshold and that floor is left out
though it is independent at the caller's scale, and the search cannot do otherwise:
keeping it would make the kept rows' product singular to working precision. **The bound
shows it.** `residual-bound` is formed from A, not from `A Aᵀ`, and resolves the distance
to rounding, so the row is reported with a residual above the threshold it was judged
against. The tests hold this case: a distance of 2⁻³⁰ and a threshold of 2⁻³², a row left
out, and a bound that reports 2⁻³⁰.

So:

- A threshold of 10⁻⁶ of the rows' size or coarser is resolved, with a factor of about 60
  to spare, and this module is the right tool.
- A threshold finer than about 10⁻⁷ of the rows' size is not. The method for that is a
  triangular factor of Aᵀ by rotations, which works with distances and not their squares;
  [the design](design-dependent-rows.md) describes it and it is not built.
- After any search, compare each residual bound with the threshold. A bound above it
  marks a row the search could not resolve.

## The vectors, and their refinement

The factorization gives each row left out a vector directly, by one triangular solve
([how](skipping.md#the-null-vector)). Its multipliers solve the normal equations of a
least squares problem, the row against the kept rows, and carry the conditioning of the
kept rows *squared*: with kept rows conditioned at 10⁵ they are good to five or six
digits.

`null-vector!` then refines them **against A, not against `A Aᵀ`**. With y the multipliers
and k the row, it forms the residual `r = a_k − A_Sᵀ y` in the column space, where it is
small and computed to rounding, and corrects y by the solution δ of
`(A_S A_Sᵀ) δ = A_S r` from the standing factor. Each correction multiplies the error by
about the squared condition number times the unit roundoff, down to a floor of about the
condition number times the unit roundoff, which is what the rounding of r allows.

**The corrections are judged by their own size.** They stop when one is more than half the
one before, in the infinity norm, or after `refinement-step-limit` of them, and the
iterate kept is the one whose next correction was smallest. `‖A_S r‖` was tried as the
measure and is the wrong one: its rounding is dominated by the directions in which the
kept rows are large, so it reaches its floor after one correction while the error along
the direction in which they are nearly dependent, the one that matters, is still falling.
A planted control found this. This is [`refine!`](refinement.md) run on the correction:
its operator returns −δ against a right-hand side of zeros and its solve hands δ back.

In the tests, with kept rows conditioned near 10⁵ the factor's own multipliers are off by
1.9 × 10⁻⁶ and the refined ones by 6 × 10⁻¹²; with kept rows conditioned near 3 × 10⁶ and
the row left out at a real distance, 2 × 10⁻³ and 1 × 10⁻¹⁰, which takes three corrections.

The vector is unique for the kept set: independent rows have one nearest combination.
Another kept set, from another order or threshold, gives other vectors, equally valid.
Vectors for all the rows left out are linearly independent, each having a 1 where the
others have 0, so **where the gap is proved they are a basis of the left null space** in
the numerical sense: m − r vectors, each with `‖Aᵀ v‖₂` under its bound.

## Using it

    analysis = analyze(A)                              once for the pattern
    find-dependent-rows!(analysis, A, threshold)       each time the values change
    for position below left-out-count(analysis):
        null-vector!(analysis, position, v)            v is that row's vector
        residual-bound(analysis, position)             compare with the threshold
    gap-bound(analysis)                                σ₍r+1₎(A) is at most this
    kept-values!(analysis, values)                     then the certificate, for σ_r

[The example](../examples/redundant-rows.ct) does this for three constraints of which one
is the sum of the others.

## Representation choices

The analysis owns a copy of the matrix's pattern and, after each search, of its values,
so that the vectors and bounds of a search refer to the values that search read even if
the caller's matrix has changed since. The operator and the solve the refinement needs
are made once, at analysis, over that storage. The kept set is recorded as the rows left
out, which are few, and a flag for each row.

**The size of the matrix does not limit the search.** It compares squares, those of the
rows in `A Aᵀ` with that of the threshold, and a square leaves the range at half the
exponent of what is squared. So the search is made on the matrix divided by the power of
two at its largest entry, with the threshold divided by the same power. The division is
exact. The multipliers of the factorization are those of the undivided matrix, so the rows
left out, the vectors and the values handed back are the same at every scale, bit for
bit; a bound is formed on the divided matrix and multiplied back when it is read, raised
by the smallest positive number where the product falls below the normal range and could
have been rounded down. [The scale-sweep roster](rosters.md#what-is-declared) holds this
from 2^-1000 to 2^1020.

What a square cannot hold is now a matter of sizes inside one matrix, and is refused by
name:

- A positive threshold more than about 2^511 times smaller than the matrix's largest
  entry is `RowsBadThreshold`: squared, it is zero or has lost its digits, and would be
  compared as another threshold. A threshold beyond every row, whose own square may not be
  finite, is taken, and leaves every row out.
- Under a threshold that can be held, a row too small to square beside the largest entry
  is shorter than the threshold and is left out, as the contract says.
- Under a threshold of zero there is no such statement to fall back on, and such a row,
  if it holds a nonzero entry, is `RowsUnderflow`. The search would otherwise leave it out
  as dependent.

The threshold is absolute. A caller with rows of very different sizes scales them first,
by [`equilibrate!`](scaling.md) or its own normalization, and should know that scaling
rows changes which combinations are small: the statement proved is about the matrix
given.
