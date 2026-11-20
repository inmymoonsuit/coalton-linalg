# Batch 11: dependent rows, and a factorization that leaves small pivots out

State: completed; verified at safety 1 and safety 0 on 2026-11-20, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contracts are in
[dependent-rows.md](dependent-rows.md) and [skipping.md](skipping.md), and the design that
preceded the code in [design-dependent-rows.md](design-dependent-rows.md).

## What was delivered, and why

- In `coalton-linalg/sparse/spd-ldlt`: `factor-skipping!`, which leaves out each index
  whose pivot is at or below a floor and factors the principal submatrix of the rest;
  `solve-skipping!`, a solve with that submatrix; `skipped-null-vector!`, the combination
  of the kept indices that an index left out amounts to; `skipping-ready?` and
  `skipped-count`.
- A new package, `coalton-linalg/sparse/dependent-rows`: `find-dependent-rows!` for the
  rows of a sparse matrix to keep and to leave out; `null-vector!` for a left null vector
  of each row left out, refined against the matrix; `residual-bound` and `gap-bound`, which
  are proved; `kept-values!`, which feeds the existing rank certificate; and
  `transpose-product-bound`, the bound by itself.

A consumer of this library locates redundant and conflicting constraints from a basis of
the left null space of its constraint matrix, which it gets from a dense singular value
decomposition. It ranked a sparse route to the same thing first among the things it would
use, on one condition: that what the route says be a proof, or be handed back.

**A written design came before any code**, and the consumer read it. The design separates
finding from proving: three statements that hold whatever found the rows, and two ways of
finding them. The consumer's answers chose between the two. Its tolerance is coarse enough
for the cheaper method to resolve, by a factor of about 60; it needs the rank as a proof and the
vectors accurate enough to reproduce the sets it shows now; and it needs no smallest
set, since it reduces any basis to one canonical form. The design page records the
answers and where the result departs from the interface it proposed.

## How it was reached

**The refinement's stopping rule was wrong in the first version, and a planted control
found it.** The multipliers of a row left out come from the factor of `A Aᵀ` and are then
refined against A. The first version ran the library's refinement on the least squares
gradient `A_S r` and stopped when that stopped halving. A control that removed a masking
step from that computation was not caught, and building the case that would catch it
showed why: with kept rows conditioned near 3 × 10⁶ the gradient reaches its rounding
floor after one correction, because its rounding comes from the directions in which the
kept rows are large, while the error along the direction in which they are nearly
dependent is still falling by three digits a step. The refinement stopped at an error of
3 × 10⁻⁶ where 10⁻¹⁰ was available. It is now run on the correction itself, which is the
error's own estimate: the masking step is gone, and on the exact roster the worst
multiplier went from 8 × 10⁻¹⁴ to 1 × 10⁻¹⁵ of its size.

**A per-product underflow term was derived and then left out of the bound.** A product
can be off by half the smallest positive number. One smallest positive number for each
column, added under the root, exceeds what that amounts to by more than a hundred orders
of magnitude and is needed anyway for the squares, so the separate term could not be
made to matter in a test. It is not in the code. The term that is there has its test.

**Two tests were missing and two controls were not controls.** A negative floor was
refused in the tests only on a matrix whose zero pivot produced the same outcome by
another path; the tests now refuse one on a definite matrix. Two planted edits left a
variable unused and did not compile, and one was equivalent to the original, since the
entry it added to a sum is always zero. All three were rewritten.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`, with every earlier roster unchanged. For this batch:

- **The skipping factorization, exactly.** 91 symmetric matrices: named cases (a pivot
  equal to the floor and one just above it, a shift that makes a zero pivot, a negative
  pivot, a zero matrix) and Gram matrices of integer factors of orders 2 to 7 and every
  rank, on their own patterns and on full ones. Each is first factored as a definite
  matrix, so that every pivot and entry of the storage holds a number the skipping call
  must replace. Then: the indices left out are those an exact elimination in the same
  order leaves out; each has a zero pivot and an empty column; every entry of the kept
  principal submatrix is within `(p + 2) u |L| |D| |Lᵀ|` of the exact one, the largest at
  0.411 of that; every kept row of a solve's exact residual is within
  `(3w + 5) u |L| |D| |Lᵀ| |x|`, the largest at 0.112; the solution is zero at each index
  left out, is unchanged bit for bit by the right-hand side's entries there, and agrees
  with itself in place; and each null vector has its ones and zeros where the contract
  puts them and `M v` within `(2w + 4) u |L| |D| |Lᵀ| |v|` at the kept indices before it,
  the largest at 0.061. Indices were left out in 74 of the cases, six at most. With
  nothing at or below the floor the factor equals `factor!`'s bit for bit.
- **Refusals of the skipping calls.** A floor that is negative, negative and subnormal,
  infinite or not a number; a shift that is not finite; a wrong value count; a value that
  is not a number: each leaves no factor of any kind. The solve and the null vector
  refuse the other kinds of factor, wrong sizes, a right-hand side that is not finite, an
  index that was kept and one out of range, with the output untouched. A four-by-four
  matrix with pivots 10⁻³²⁰, 10²⁸⁰ and 10³⁰⁶ gives a null vector that overflows; it is
  reported and not written. `solve!`, `solve-cholesky!` and the condition estimate refuse
  a skipping factor.
- **Dependent rows, exactly.** 118 integer matrices, most with up to three rows that are
  small integer combinations of the others, and named cases: a repeated row, a row with
  no stored entry, a row of stored zeros, no columns, no rows, more rows than columns, a
  hub column. In each: the rows kept number the exact rank and are independent, by exact
  elimination; each vector has 1 at its row and zeros at the other rows left out; its
  multipliers agree with the exact combination to 1.3 × 10⁻¹⁵ of their size at worst; its
  residual bound is at least the exact norm of the exact product; the gap bound is at
  least the root of the sum of the exact squares; the exact inertia of `A Aᵀ − g² I`
  confirms that no more singular values than rows kept exceed the gap bound g; and the
  kept values are the source's with zeros for the rows left out. 218 rows were left out
  in 95 of the cases.
- **Near the threshold.** A row at distance 2⁻⁴ down to 2⁻²⁴ from another, with the
  threshold at half and at twice the distance: kept in the first, left out in the second,
  and the residual bound within a part in a million above the exact distance. Two rows
  near a third, both left out, with a gap bound that carries both. **Below the method's
  floor**: a distance of 2⁻³⁰ and a threshold of 2⁻³²; the row is left out though it is
  independent at that threshold, and its residual bound reports the distance, above the
  threshold.
- **Refinement.** With kept rows conditioned near 10⁵ the factor's own multipliers are off
  by 1.9 × 10⁻⁶ and the refined ones by 6.1 × 10⁻¹². With kept rows conditioned near
  3 × 10⁶ and the row left out 6 × 10⁻⁸ off the span: 2.1 × 10⁻³ and 1.0 × 10⁻¹⁰. The
  factor's own vector is measured on a separate analysis, and each case asserts that it
  is as poor as the case needs it to be.
- **The kept rows proved independent, on the whole matrix's plan.** 364 statements about
  the kept rows' smallest singular value, at thresholds from 2⁻⁴⁰ to 40, through the
  certificate with zeros written for the rows left out and the pivot counts the contract
  gives. 253 were certified, 177 of them with rows left out; none of those is false by
  exact inertia; 111 were withheld and all 111 are false. **A matrix that fools the
  search**: ten rows that are the columns of a Kahan matrix; every row is kept at a
  threshold of 0.05, the smallest singular value is below it by exact inertia, and the
  certificate at the threshold is withheld.
- **The product bound.** 77 matrices and vectors, the value returned against the exact
  norm of the exact product: a column in which the plain computed sum is zero and the
  exact one is 1; a product that underflows to zero; 61 squares of which every addition
  after the first rounds down; a case with no rounding at all, bounded within 10⁻¹³ of
  the exact norm; and seeded cases with magnitudes over eight orders, in 31 of which the
  plain computed norm is below the exact one. An overflow, and an overflow of both signs
  in one column, give infinity.
- **Lifetime.** Before any search every reader says there is nothing. A threshold that is
  negative, not a number, infinite, or too large to square leaves no set, whatever stood
  before, and so does a matrix of another pattern or shape, which is a typed error. The
  gap bound is withheld until every vector has been taken and withdrawn by the next
  search. The same analysis on other values and back gives the same vector bit for bit.
  A product that is finite with a factor that is not is `RowsNonFinite`.
- **Encapsulation.** A compiled external caller reaches the analysis through its public
  operations; six probes that reach into it, or pass one status type for another, are
  rejected by the compiler.
- **Allocation.** A compiled caller that searches, takes every null vector, reads the gap
  bound and takes the kept values, and one that factors with skipping, solves and takes a
  null vector: zero bytes over 10,000 rounds after 10,000 warm-up rounds, at both
  policies.

**Planted source controls.** Each of these edits was applied alone to the final source
and the focused tests run. Every planted source compiled, every one failed the tests, and
the restored source passed.

| Source | Planted defect | Caught |
| --- | --- | --- |
| spd-ldlt | A pivot equal to the floor is kept | yes |
| spd-ldlt | A skipped column is not emptied | yes |
| spd-ldlt | A skipped pivot is not set to zero | yes |
| spd-ldlt | A skipping factor is marked as a solving factor | yes |
| spd-ldlt | A negative floor is accepted | yes |
| spd-ldlt | An infinite floor is accepted | yes |
| spd-ldlt | The skipping solve does not clear a skipped entry | yes |
| spd-ldlt | The skipping solve accepts a solving factor | yes |
| spd-ldlt | The skipping solve divides by a skipped pivot | yes |
| spd-ldlt | The skipping solve accepts a right-hand side that is not finite | yes |
| spd-ldlt | The null vector does not read the skipped row | yes |
| spd-ldlt | The null vector's substitution runs past the skipped index | yes |
| spd-ldlt | A kept index is given a null vector | yes |
| spd-ldlt | A null vector that is not finite is written | yes |
| spd-ldlt | The null vector's output is not cleared first | yes |
| spd-ldlt | The null vector has the combination's sign | yes |
| spd-ldlt | The null vector's size is not checked | yes |
| dependent-rows | The threshold is not squared | yes |
| dependent-rows | A negative threshold is accepted | yes |
| dependent-rows | A refused search leaves the earlier set standing | yes |
| dependent-rows | A row left out is recorded by its place in the elimination | yes |
| dependent-rows | A row left out is flagged as kept | yes |
| dependent-rows | A new search keeps the earlier bounds | yes |
| dependent-rows | The factor's vector is not refined | yes |
| dependent-rows | The refinement stops after one correction | yes |
| dependent-rows | The refinement's correction has the wrong sign | yes |
| dependent-rows | The refinement's correction is not applied | yes |
| dependent-rows | The refinement's residual is formed without the kept rows' combination | yes |
| dependent-rows | The refinement leaves the row itself out of its residual | yes |
| dependent-rows | The gap bound ignores a missing residual bound | yes |
| dependent-rows | The gap bound takes the last residual bound alone | yes |
| dependent-rows | The product bound leaves out the rounding of the sum | yes |
| dependent-rows | The product bound leaves out the underflow term | yes |
| dependent-rows | The product bound leaves out the margin for its own rounding | yes |
| dependent-rows | The product bound returns what is not a number | yes |
| dependent-rows | The kept values keep the rows left out | yes |
| dependent-rows | A position past the rows left out is accepted | yes |
| dependent-rows | The null vector's size is not checked | yes |
| dependent-rows | The residual bound is recorded for another position | yes |

## Limits

- **The search does not resolve a distance below about 1.5 × 10⁻⁸ of a row's size.** It is
  a factorization of `A Aᵀ`, in which distances are squared. A row in that range is left
  out and its residual bound says so. The method that resolves it is designed and not
  built.
- **The kept set depends on the elimination order**, and is not chosen for its
  conditioning. Its independence at a tolerance is a separate proof that the caller runs,
  and that can be withheld.
- **Two outcomes have no test that reaches them through this module.** `NullNonFinite`
  needs a vector or a bound that overflows; the kernel outcome it passes on is tested
  directly, and no matrix was found whose product is finite and whose null vector is
  not. A correction that cannot be computed inside the refinement ends it with the best
  iterate, by the path that a non-finite value takes; that path is `refine!`'s, tested
  there, and is not reached here.
- **No timing is reported.** The cost is one factorization of `A Aᵀ`, with a triangular
  solve, at most `refinement-step-limit + 1` pairs of products with A and as many solves
  for each row left out. What that saves a consumer against a dense decomposition is the
  consumer's to measure.
- The threshold is absolute, and a dense column of A makes `A Aᵀ` dense, as it does for
  every use of the Gram route in this library.
- One host and finite rosters establish neither universal behaviour nor portability beyond
  the supported runtime.
