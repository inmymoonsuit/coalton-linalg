# Batch 12: sparse LU with partial pivoting

State: completed; verified at safety 1 and safety 0 on 2026-11-20, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract is in [lu.md](lu.md).

## What was delivered, and why

A new package, `coalton-linalg/sparse/lu`:

- `factor!`: `P A Q = L U` for a square sparse matrix, each pivot the largest magnitude in
  its column, or the matrix's diagonal entry where that is within a threshold of the
  largest.
- `solve!` and `solve-transpose!`: solves with the matrix and with its transpose.
- `make-lu`, `make-lu-with-order`, `column-order`: the handle, for matrices of one order,
  with a column order from the pattern or from the caller and a budget of factor entries.
- `lu-growth` and `lu-backward-bound`; the row order chosen, the column of a singular
  refusal, the entry counts, and the factors as borrowed views.

Every sparse factorization this library had fixes its order from the pattern and never
exchanges rows. That is the right design for a positive definite matrix and a checked one
for a quasi-definite matrix, and it leaves a square matrix of no special structure with
nothing. This batch is built ahead of any consumer's measurement: it is the second item
of a consumer's ranking of what it would use, and the general tool a library of this kind
is expected to have.

## How it was reached

**The first version passed every test it was given, and the planted controls then
removed code.** Three things came out of the first control run:

- The transposed solve checked each value for finiteness after its first pass and again
  after its second. A control that removed the first check was not caught, and cannot
  be: a value that is not finite after the first pass is not finite after the second,
  which looks at every entry. The first check is gone.
- A refused call cleaned its work column, and every call also zeroes that column when it
  starts. One of the two was enough. The cleaning is gone, since the zeroing costs the
  order of the matrix in a call that costs at least that.
- The underflow term of the backward bound had a test that the relative term passed by
  itself, because the bound sums `|L| |U|` and any entry of ordinary size dominates it.
  The test is now a matrix whose every entry is below 2⁻¹⁰³⁹, where the relative term is
  2¹⁶ times too small and the underflow term carries the bound.

A column order longer than the matrix's order was accepted and silently cut short until
a control showed that the test offered only a shorter one. Six of the first forty
controls did not compile, each for a variable left unused, and were rewritten.

**A claim was withdrawn from a docstring.** The column order was described as bounding
the fill of U under any row exchanges. That is a known result and nothing here tests it,
so the text now says only what the code relies on: the pattern of `Aᵀ A` does not depend
on the order of the rows.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`, with every earlier roster unchanged. For this batch:

- **The factorization, exactly.** 180 matrices: named cases (a cyclic permutation, where
  every pivot is an exchange; triangular matrices; stored zeros; a dense matrix of order
  8, whose factors outgrow their first allocation) and seeded integer matrices of orders
  2 to 8 at three densities, with and without a full diagonal, with a dense row, each in
  the order `column-order` gives, in the natural order, and in the reversed order at a
  threshold of 0.1. For each that factors: the factors have the structure the contract
  states; every entry of `P A Q − L U`, computed exactly, is within
  `γ₍r_j+1₎ (|L| |U|)[i,j]`, the largest at 0.800 of that; no multiplier exceeds one over
  the threshold; at a threshold of 1 every multiplier of magnitude exactly 1 is a tie
  that went to the diagonal or to the smaller row; the growth reported is the exact one;
  the backward bound is at least the exact Frobenius norm and is the contract's formula
  on the exported factors, rounded up; each row of a solve's exact residual is within
  `γ₍3n₎ |L| |U| |x|`, the largest at 0.167, and likewise each column for the transposed
  solve; and both solves give the same bits in place. Rows were exchanged in 88 of the
  cases. 27 were refused as singular, and each of those is singular by exact elimination.
- **The pivot rule, by case.** A column holding 1 on the diagonal and 2 below it: the
  diagonal at a threshold of one half, the larger entry at six tenths and at one. A tie
  between two rows off the diagonal goes to the smaller; a tie with the diagonal goes to
  the diagonal though its row is the larger.
- **Outcomes.** Singular: a column with no entry, a column of stored zeros, a repeated
  column in either order, a row with no entry, each naming the column where it showed.
  Not finite: a candidate for the pivot that overflows, and an entry of U that overflows
  while every candidate stays finite. The budget: a dense matrix of order 4 has twelve
  entries off the diagonals; budgets of 0 to 11 are refused and 12 is not, and a handle
  that ran out factors a smaller matrix afterwards. A threshold of 0, a negative one, one
  above 1, infinity and not-a-number are refused, and the smallest positive number is
  accepted. A matrix of another order or shape is refused. Each refusal leaves no
  factors, whatever stood before. After a refused call the same handle factors a sound
  matrix as a fresh handle does, bit for bit.
- **Solves.** Wrong sizes and right-hand sides that are not finite are refused with the
  output untouched, by both solves; a solution that overflows, from a pivot of 10⁻³²⁰, is
  reported and not written.
- **Construction.** A matrix that is not square, and an order that is too short, too long,
  out of range or repeating, are typed errors; a budget too small for the pattern of
  `Aᵀ A` is the analysis's.
- **Underflow.** The matrix described above: its exact squared error is 2⁻²¹⁵², above the
  square of the relative term, and at most the square of the bound reported.
- **Lifetime.** One handle on a matrix, then on a matrix of another pattern, then on the
  first again: rows, pivots and both factors agree bit for bit with a fresh handle each
  time.
- **Encapsulation.** A compiled external caller reaches the handle through its public
  operations; four probes that reach into it or pass one status type for another are
  rejected by the compiler.
- **Allocation.** A compiled caller that factors a dense matrix of order 8, reads the
  bound, solves and solves with the transpose in place: zero bytes over 10,000 rounds
  after 10,000 warm-up rounds, at both policies.

**Planted source controls.** Each of these edits was applied alone to the final source
and the focused tests run. Every planted source compiled, every one failed the tests, and
the restored source passed.

| Planted defect | Caught |
| --- | --- |
| A tie between rows goes to the larger row | yes |
| The diagonal is never preferred | yes |
| The diagonal is preferred whatever the threshold | yes |
| A diagonal row already pivoted is preferred again | yes |
| The column's update is added | yes |
| A row met in the reach is not marked | yes |
| A position is written to the order when first met, not when finished | yes |
| A candidate that is not finite is accepted | yes |
| An entry of U that is not finite is accepted | yes |
| The singular column is not recorded | yes |
| The budget is not checked | yes |
| The budget refuses an exact fit | yes |
| Grown value storage loses what it held | yes |
| Grown index storage loses what it held | yes |
| The work column keeps the entries of U | yes |
| The work column keeps the candidates | yes |
| A refused call leaves the earlier factors standing | yes |
| A refused call leaves the earlier counts standing | yes |
| A matrix of another order is accepted | yes |
| A threshold above one is accepted | yes |
| A threshold of zero is accepted | yes |
| A column order that repeats an index is accepted | yes |
| A column order of the wrong length is accepted | yes |
| Growth leaves the pivots out | yes |
| The backward bound leaves out the relative term | yes |
| The backward bound leaves out the underflow term | yes |
| The backward bound leaves out the unit diagonal of L | yes |
| The backward bound leaves the pivots out of the row sums | yes |
| The backward bound counts one product too few | yes |
| The solve does not follow the row exchanges | yes |
| The solve does not undo the column order | yes |
| The solve does not divide by the pivot | yes |
| The solve accepts a right-hand side that is not finite | yes |
| The solve writes a solution that is not finite | yes |
| The transposed solve does not follow the column order | yes |
| The transposed solve reads L by row of the matrix, not by step | yes |
| The transposed solve does not undo the row exchanges | yes |
| The transposed solve accepts a right-hand side that is not finite | yes |
| The transposed solve writes a solution that is not finite | yes |
| A solve is given without factors | yes |

## Limits

- **The matrices of the roster are small**, orders 1 to 8, because every entry is checked
  in exact arithmetic. No factorization of a large matrix was run, and no timing is
  reported: how this compares with the symmetric factorizations on matrices both can
  take, and how its fill compares with other column orders, are not measured.
- **The column order forms the pattern of `Aᵀ A`**, which is dense where A has a dense row.
  A caller with such a matrix supplies an order.
- **Each call searches for its pivots.** There is no refactorization on a fixed sequence,
  no condition estimate, no scaling inside and no dense kernel.
- **Singular is exact.** A matrix that is singular in exact arithmetic is usually factored,
  with a pivot at rounding level; the growth and the backward bound do not flag that, and
  a rank is not this package's to decide.
- The one-number bound sums `|L| |U|` in place of taking a norm of it, and so can exceed
  the Frobenius norm of that product by up to the root of the number of its entries. On
  the roster the exact error squared was at most 0.011 of the bound squared.
- One host and finite rosters establish neither universal behaviour nor portability beyond
  the supported runtime.
