# Scaling by powers of two

Package: `coalton-linalg/sparse/scaling`. Source: [scaling.ct](../src/sparse/scaling.ct).
State: completed as Batch 9 and corrected at the edges of the floating-point range as
Batch 9a; verified at safety 1 and safety 0. The [batch record](batch9.md) and [the
correction's](batch9a.md) hold the controls.

## What it is for

A matrix whose rows or columns differ by many orders of magnitude is hard on every method
that follows. A geometric constraint Jacobian is the standing example: a row that
measures a length in millimetres beside one that measures an angle in radians, a column
for a coordinate beside one for a rotation. Forming AᵀA or A Aᵀ squares the spread, a
threshold chosen for one row means something else for the next, and a pivot or a
regularization sized for the large entries swamps the small ones.

Scaling the rows and columns first removes the part of that which is only a choice of
units. **Here every scale factor is a power of two, so scaling changes no digit of any
entry**: the scaled matrix is exactly `R A C`, and a solution carried back through the
factors is exactly the solution of the system that was posed.

## Interface

A matrix is a checked [`csc:Matrix`](primitives.md), which supplies the pattern, together
with a values array in that matrix's stored order. The array may be the matrix's own
`entries`, or a caller's array of fresh values for the same pattern.

- `equilibrate!(matrix, values, row-scales, column-scales) -> Void` writes powers of two r
  and c such that **every row and every column that has a nonzero has its largest scaled
  magnitude `|r[i] a[i,j] c[j]|` in [1/2, 1)**. A row or column with no nonzero gets the
  factor 1. **The exception is a row or a column that no factor reaches**: the largest
  factor there is, 2¹⁰²³, is not enough to bring it to 1/2. It gets that factor, and its
  largest scaled magnitude is then below 1/2, possibly far below. That is
  - a row whose largest magnitude is below 2⁻¹⁰²⁴, and
  - a column whose largest magnitude *after the row factors* is below 2⁻¹⁰²⁴, which takes
    an entry more than 2¹⁰²³ times smaller than the largest in its own row and no larger
    entry elsewhere in the column.
- `scale!(matrix, values, row-scales, column-scales, out) -> Void` writes
  `r[i] a[i,j] c[j]` for every stored entry, in the matrix's order. The factors need not
  come from `equilibrate!`; what is claimed of the result is claimed for any factors that
  are powers of two.
- `scale-vector!(scales, x, out) -> Void` writes `scales[i] x[i]`. The output may be x
  itself.
- `row-maxima!(matrix, values, out)` and `column-maxima!(matrix, values, out)` write the
  largest magnitude of each row or column, 0 where there is no entry.
- `equilibrate-symmetric!(matrix, values, scales) -> DiagonalStatus`, for the lower
  triangle of a symmetric matrix with a positive diagonal, writes powers of two d such
  that **every scaled diagonal entry `d[i]² a[i,i]` is in [1/2, 2)**. If the matrix is
  positive semidefinite, every scaled entry is then below 2 in magnitude. A diagonal entry
  that is zero or negative returns `DiagonalNotPositive` and writes nothing. The exception
  is a subnormal diagonal entry, whose factor is held to 2⁵¹¹.
- `scale-symmetric!(matrix, values, scales, out) -> Void` writes `d[i] a[i,j] d[j]` for
  every stored entry of that lower triangle.

| DiagonalStatus | Meaning |
| --- | --- |
| DiagonalScaled | The scales were written |
| DiagonalNotPositive | Some diagonal entry is zero or negative; nothing was written |

Misuse is an exception, as elsewhere: a length that does not match the matrix is
`InvalidDimensions`; a symmetric call on a matrix that is not square the same; a column of
a symmetric matrix that does not begin with its diagonal is `InvalidStructure`; a
non-finite value or scale factor is `NonFiniteInput`; an output that is one of the inputs,
other than `scale-vector!`'s x, is `InvalidAliasing`. Each is raised before anything is
written. A result that is not finite, which caller-chosen factors can produce, is
`NonFiniteResult` and leaves the output unspecified. Nothing is allocated by any of these
calls.

## Using it

**A linear system.** To solve `A x = b`: equilibrate, scale A to `B = R A C`, scale b to
`R b`, solve `B y = R b`, and return `x = C y`. Both vector steps are `scale-vector!`.
Nothing is approximated by the scaling, so x solves the original system as well as y
solves the scaled one.

**A symmetric positive definite system.** Scale to `D A D`, solve `(D A D) y = D b`, and
return `x = D y`. [A compiled example](../examples/scaled-solve.ct) does this for a
matrix whose diagonal spans twelve orders of magnitude.

**What scaling leaves unchanged, and what it changes.** R and C are invertible, so
`A x = b` and `R A C y = R b` with `x = C y` have the same solutions. What can change is
which x is chosen when there is not exactly one:

| The system | Row scaling | Column scaling |
| --- | --- | --- |
| Has exactly one solution | Unchanged | Unchanged |
| Consistent with many solutions: the minimum-norm one | **Unchanged**: the set of solutions is the same and so is its shortest member | **Changed**: the shortest y minimizes `‖C⁻¹ x‖`, not `‖x‖` |
| Inconsistent: the least-squares one | **Changed**: it minimizes `‖R (A x − b)‖`, not `‖A x − b‖` | Unchanged if it is unique; otherwise changed as above |

So a caller who needs the minimum-norm solution of a consistent system in the original
variables may scale rows freely and must leave the columns alone, passing a column-scale
vector of ones to `scale!`; the row factors of `equilibrate!` are computed before the
column factors and do not depend on them. A caller for whom the unscaled norm adds
millimetres to radians may want the column scaling for its own sake. They are different
problems, and this library does not choose between them. A proof about the singular
values of A is a statement about A: scaling changes them, so such a proof is made on the
matrix it is about.

**Thresholds.** After `equilibrate!` every nonzero entry is below 1 in magnitude and
every row and column reaches 1/2, so a threshold or a regularization can be stated
against 1 for every row alike. **The rows and columns of the exception above do not
reach 1/2**, and a threshold against 1 treats them as negligible. They are found by
their factor: a row or column whose factor is 2¹⁰²³ should be looked at with
`row-maxima!` or `column-maxima!` on the scaled values before a threshold is trusted.
A matrix has such a column only if the entries of some row span more than 2¹⁰²³, three
hundred orders of magnitude, and such a row only if every entry of the row is below
2⁻¹⁰²⁴.

## What is and is not claimed

- **Exactness.** Each scaled value is the exact product `r a c` whenever that product is
  zero or a normal number, **whatever the sizes of the three**: a value or a factor may be
  subnormal, and the product of the value with either factor alone may underflow or
  overflow, and the result is still exact. This holds for `scale!` and `scale-symmetric!`
  with any factors that are powers of two, the caller's own included. A product below
  the normal range is rounded, possibly to zero, as any floating-point product would be;
  after `equilibrate!` that happens only to an entry more than about 2¹⁰²² times smaller
  than the largest in its row and also far below the largest in its column. The pattern
  is not changed; a stored entry that became zero is still stored.
- **One pass each way.** `equilibrate!` scales rows and then columns once. That is enough
  for the property stated, by the argument below. It does not make all rows alike beyond
  that: a row can have one entry near 1 and the rest tiny.
- **Not optimal.** These scalings do not minimize a condition number. For a positive
  definite matrix, scaling the diagonal to 1 is known to be within a factor of the order
  of the best diagonal scaling, a result of van der Sluis; a power of two differs from
  the exact factor by at most √2 on each side. No such statement is made for the
  rectangular case.
- **Rows first.** Scaling columns first gives a different pair of factors with the same
  property. A caller who wants that order passes the transpose.

## Derivation

**The factor.** A positive number x is `f 2ᵉ` with f in [1/2, 1) and e an integer, in
exactly one way. The row factor is 2⁻ᵉ, which brings x to f. It is found without
decoding the number: divide by 2⁵¹², 2²⁵⁶, …, 2 in turn whenever the number is at least
that large, which leaves it below 2; halve it once if it is at least 1; then multiply by
the same powers in turn whenever the product stays below 1. Every step multiplies or
divides by a power of two and is exact, and the ten powers together reach 2¹⁰²³, which is
where the bound on the factor comes from.

**The product.** `r a c` is three numbers and floating point multiplies two. Either order
can fail though the result is an ordinary number: for a = 2.5, r = 2⁻¹⁰⁷⁴ and c = 2¹⁰²³
the result is 2.5 × 2⁻⁵¹, while `r a` is not representable and rounds to 2 × 2⁻¹⁰⁷⁴, and
`a c` overflows. The two factors are therefore multiplied with each other first. Their
product is a power of two and is exact unless it leaves the range, and the value is then
multiplied once, with no intermediate to lose anything in. Where the factors' product
does leave the range, both factors lie on the same side of 1: both below it if the
product underflowed to zero, both above it if it overflowed. The value is then
multiplied by them in turn, and each step moves it the same way, so every intermediate
lies between the value and the result: it is normal if they are, and exact.

**Rows, then columns.** After the row pass every row with a nonzero has its largest
magnitude in [1/2, 1), so every entry is below 1. Take a column with a nonzero and let x
be its largest magnitude after the row pass. Its factor is 2⁻ᵉ for x = f 2ᵉ, and x < 1
gives e ≤ 0, so the factor is at least 1: no entry shrinks, and every row keeps a largest
magnitude of at least 1/2. Every entry of the column is at most x in magnitude and x is
brought to f < 1, so every entry stays below 1, in every row. So after the column pass
each column's largest magnitude is in [1/2, 1) and each row's still is.

**A column's largest row-scaled magnitude** can be far below the normal range, where it
could not be formed: it would be rounded, across a power of two in the worst case, or
lost. Each row-scaled magnitude is therefore formed 2¹⁰²³ times larger, by the product
above, which is exact down to 2⁻²⁰⁴⁵ of the scale. A largest lifted magnitude `f 2ᵉ` with
e from 0 to 1023 asks for the factor 2^(1023 − e), which is at least 1 and at most 2¹⁰²³;
one below 1/2 asks for more than there is and is given 2¹⁰²³, and that is the column of
the exception. A column is empty, with the factor 1, only if every stored value in it is
zero, and not because its values vanished under the row factors.

**The symmetric factor.** For a diagonal entry `a = f 2ᵉ` the factor is `2⁻ᵏ` with
k = ⌊e/2⌋, the only power of two for which `2⁻²ᵏ a = f 2^(e − 2k)` is in [1/2, 2), since
e − 2k is 0 or 1. It is found as above with the even powers 2⁵¹², …, 2², halving the
exponent of the factor at each. Scaling both sides by the same D keeps the matrix
symmetric and, because `D A D` is a congruence, positive definite if it was. For a
positive semidefinite matrix `|a[i,j]|² ≤ a[i,i] a[j,j]`, so every scaled entry is at most
the geometric mean of two scaled diagonal entries, each below 2.

## Representation choices

There is no plan and no handle: a `csc:Matrix` already holds the row of every entry, the
factors are the caller's arrays, and the row pass uses the row-scale array as its own
workspace. The symmetric functions take a lower-triangular `csc:Matrix`, the form
`analyze-csc` accepts, and not a second representation of a symmetric pattern.
