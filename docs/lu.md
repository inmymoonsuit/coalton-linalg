# Sparse LU with partial pivoting

Package: `coalton-linalg/sparse/lu`. Source: [lu.ct](../src/sparse/lu.ct).
State: completed as Batch 12; verified at safety 1 and safety 0. The [batch
record](batch12.md) holds the controls.

## What it is for

A square system `A x = b` whose matrix is not symmetric, or is symmetric and not safely
factorizable in a fixed order: the Jacobian of a square nonlinear system, a matrix with a
border, a transport or circuit matrix, a saddle-point matrix with nothing on part of its
diagonal. The other sparse factorizations of this library fix their elimination order
from the pattern and never exchange rows; that is right for a positive definite matrix
and is checked by a growth figure for a quasi-definite one, and it is wrong in general.
This one **chooses each pivot from the values**, so it needs no property of the matrix
beyond being nonsingular.

The price is that the fill depends on the values and is not known in advance. The
factors therefore live in storage that grows, up to a budget the caller sets.

## Interface

- `make-lu(matrix) -> Lu`: a handle for square matrices of this matrix's order, with a
  column order chosen for its pattern and the default budget.
- `make-lu-with-order(n, order, maximum-entries) -> Lu`: the order given, as the column
  eliminated at each step, and the budget given.
- `column-order(matrix, maximum-entries) -> order`: the order `make-lu` uses.
- `factor!(lu, matrix, threshold) -> LuStatus`: factors any square matrix of the handle's
  order. The pattern need not be the one the order was chosen for; another pattern fills
  differently and is otherwise the same call.
- `solve!(lu, b, out)`, `solve-transpose!(lu, b, out) -> LuSolveStatus`: `A x = b` and
  `Aᵀ x = b` with the standing factors. `out` may be `b`.
- `lu-ready?`, `dimension`, `entry-limit`, `singular-column`, `lower-entry-count`,
  `upper-entry-count`, `lu-growth`, `lu-backward-bound`.
- The factors, as borrowed views: `row-permutation`, `column-permutation`,
  `lower-pointers`, `lower-rows`, `lower-values`, `upper-pointers`, `upper-rows`,
  `upper-values`, `upper-diagonal`.

| LuStatus | Meaning |
| --- | --- |
| LuSuccess | Factors stand |
| LuSingular | A column had no nonzero entry left to pivot on; `singular-column` names it, in the matrix's numbering |
| LuNonFinite | An entry of a factor was not finite |
| LuStorageExceeded | The factors need more entries than the handle's budget |
| LuBadMatrix | The matrix is not square of the handle's order |
| LuBadThreshold | The threshold is not in (0, 1] |

| LuSolveStatus | Meaning |
| --- | --- |
| LuSolveSuccess | The output holds the solution |
| LuSolveNoFactor | No factors stand |
| LuSolveBadSize | A vector's length differs from the order |
| LuSolveNonFiniteRhs | The right-hand side holds an entry that is not finite |
| LuSolveNonFiniteResult | The computed solution was not finite; the output is unchanged |

`make-lu` raises `InvalidDimensions` for a matrix that is not square, and
`make-lu-with-order` raises `InvalidPermutation` for an order that is not one.
`column-order` raises what the Gram analysis and the ordering raise under their budget.

Every refused `factor!` leaves no factors, whatever stood before, and every refused solve
leaves its output untouched. One handle serves one call at a time.

## The factorization

    P A Q = L U

Q is the handle's column order, fixed. L is unit lower triangular, U upper triangular,
and P is the row order the pivoting produced: `row-permutation[k]` is the row of A taken
as pivot at step k and `column-permutation[k]` the column eliminated then.

**The pivot rule.** At step k the column is first brought up to date with the steps
before it. Among its entries in rows not yet taken:

1. the largest magnitude is found, the smallest row among equals;
2. if the entry on the diagonal of A in this column, the one in the row of the column's
   own index, is among them, is not zero, and is at least `threshold` times that largest
   magnitude, it is the pivot instead.

With a threshold of 1 this is partial pivoting, with ties going to the diagonal: no
multiplier exceeds 1 in magnitude. A smaller threshold keeps the diagonal where it is
nearly as good, which preserves sparsity for a matrix whose order was chosen with its
diagonal in mind, and allows multipliers up to `1 / threshold`.

**Singular means exactly singular in the arithmetic performed.** `LuSingular` is returned
when every remaining entry of a column is exactly zero, a column with no entries
included. A matrix that is singular in exact arithmetic usually leaves a pivot at
rounding level and is factored; nothing here decides a rank. The growth and the backward
bound say how far the factors can be trusted, and a rank is the business of
[the orthogonal decompositions](orthogonal.md) or of [dependent rows](dependent-rows.md).

**The column order.** `column-order` is minimum degree on the pattern of `Aᵀ A`, the
library's own exact ordering applied to the Gram pattern. That pattern does not depend on
the order of the rows, which is why it can be used before the pivoting is known. It is
costly where A has a dense row, since `Aᵀ A` is then dense; a caller with a better order
passes it.

## Accuracy

Write u = 2⁻⁵³, `γₘ = m u / (1 − m u)`, and r_j for the number of entries of U above the
diagonal in column j.

**Entry by entry, barring underflow:**

    | (P A Q − L U)[i,j] |  ≤  γ₍r_j+1₎ (|L| |U|)[i,j].

An entry of column j is the matrix's entry less at most r_j products, and an entry of L
is divided once more; that is the whole computation, so the bound follows from the
rounding of one sum. The tests check it exactly, in rational arithmetic, on every entry
of every factorization of the roster.

**As one number, with underflow.** `lu-backward-bound` is

    ‖ P A Q − L U ‖_F  ≤  γ₍r+1₎ · Σ(|L| |U|)  +  η · E · (r + 2 + max|U|),

rounded upward, with r the largest r_j, `Σ(|L| |U|)` the sum of every entry of that
product, η = 2⁻¹⁰⁷⁴ and E the number of stored entries, diagonal included. The sum is
computed as the sum over the steps of a column sum of |L| times a row sum of |U|, which
costs one addition for each stored entry. The Frobenius norm is at most the sum of the
magnitudes, so this is coarser than a norm of `|L| |U|` would give and needs no pass over
the product. The underflow term covers a product or a quotient that rounds to a multiple
of η. For a matrix of ordinary size the relative term is larger by hundreds of orders of
magnitude and the underflow term changes nothing; it is what makes the statement true
for a matrix whose entries are themselves near the bottom of the range, and the tests
hold one, two by two with every entry below 2⁻¹⁰³⁹, whose error is about 2¹⁶ times what
the relative term allows.

**A solve** satisfies `(A + ΔA) x̂ = b` with `|P ΔA Q| ≤ γ₍3n₎ |L| |U|`, the classical bound
for two triangular solves after the factorization, and the tests check each row of the
exact residual against it; `solve-transpose!` likewise by columns. For a tighter answer
the solve is a correction step for [`refine!`](refinement.md), with the matrix's own
product as the operator.

**Growth.** `lu-growth` is the largest magnitude in U over the largest in A, measured.
With partial pivoting it is at most 2ⁿ⁻¹ and almost always small; it is the figure to
read when a threshold below 1 is used, since the multipliers are then not bounded by 1.

## Storage

L and U are each held in a pair of arrays that **grow by doubling** when a column does
not fit, up to the handle's budget. The budget counts the entries of both factors off
their diagonals; a column that would pass it ends the call with `LuStorageExceeded`.

Growth allocates. Once the storage has reached what a matrix needs, factoring that
matrix again, or any matrix that fills no more, allocates nothing, and a solve never
does; a compiled caller that factors and solves twice is measured at zero bytes after
its warm-up. `lower-rows`, `lower-values`, `upper-rows` and `upper-values` are views of
the storage as it is when asked for, and are replaced when it grows.

## How it works

The factorization goes column by column. Column j of U and of L is the solution of a
triangular system with the columns of L already computed and column j of A as the
right-hand side. Solving it densely would cost the order of the matrix for every column.
Two facts make it cost what it changes:

- **Which steps the column reaches is a graph question.** An entry of the column in the
  row pivoted at step k will spread, through column k of L, to every row that column
  holds; those of them that are themselves pivot rows spread further. The steps reached
  are the nodes reachable from the column's own pivot rows in the graph whose edges are
  the entries of L, found by a depth-first search over that structure alone.
- **The order in which to apply them comes with the search.** A step is written down when
  everything it spreads to has been written down, filling a list from its end, so the
  list holds each step before those it feeds. Applying the columns of L in that order
  brings every entry up to date before it is used.

The rows met that are not yet pivoted are the candidates for the pivot. The search is
written without recursion, with a stack of steps and a cursor for each.

The transposed solve reads the same storage the other way: column j of U is row j of
`Uᵀ`, and column k of L, whose rows were all pivoted after step k, is what row k of `Lᵀ`
needs.

## What is not here

- **No refactorization on a fixed pivot sequence.** Each `factor!` searches for its
  pivots. Where the values change little from one matrix to the next, reusing the
  sequence and the structure would be faster; it is not built.
- **No scaling inside.** Partial pivoting compares magnitudes within a column, which
  means what the caller's units make it mean. [Equilibration](scaling.md) first is the
  caller's choice.
- **No condition estimate.** `solve!` and `solve-transpose!` are what one is made of.
- **No supernodes and no dense kernels.** The inner loop is one multiply-add for each
  entry.

## Representation choices

The factors are stored by columns because the factorization produces them by columns. L
keeps the row of the matrix for each entry and not its pivot step, since a row's step is
not known when its entry is stored; the solves translate. U keeps the pivot step. Within
a column the entries of either factor are in the order the search met them, which no
solve depends on. Storage grows in place of being sized from a bound, because the
useful bounds on the fill of pivoted factors are loose and a caller should not have to
supply one.
