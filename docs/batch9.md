# Batch 9: scaling by powers of two

State: completed; verified at safety 1 and safety 0 on 2026-11-20, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract and the derivation are in
[scaling.md](scaling.md).

## What was delivered, and why

A new package, `coalton-linalg/sparse/scaling`:

- `equilibrate!`: power-of-two row and column factors for a sparse matrix that bring the
  largest magnitude of every row and every column into [1/2, 1).
- `scale!`: the scaled values `r[i] a[i,j] c[j]`, each the exact product.
- `scale-vector!`: a vector through a set of factors, for right-hand sides and solutions.
- `row-maxima!`, `column-maxima!`: the largest magnitude by row and by column.
- `equilibrate-symmetric!`, `scale-symmetric!`: for the lower triangle of a symmetric
  matrix with a positive diagonal, one set of factors applied on both sides, bringing the
  diagonal into [1/2, 2).

This is the first batch built ahead of a consumer's measurement. A constraint Jacobian
mixes units, lengths beside angles, and everything downstream of it in this library,
the Gram products, the regularized augmented matrix, a threshold on a pivot, is stated
against a scale that such a matrix does not have. Scaling by powers of two removes the
part of the difficulty that is only a choice of units and perturbs nothing.

**The consumer that will use it set a limit on what it is for**, and the contract carries
that limit. Scaling changes singular values, so a proof about the singular values of a
matrix is made on that matrix and not on a scaled one. And scaling can change which
solution is chosen where there is not exactly one: row scaling leaves the minimum-norm
solution of a consistent system unchanged and changes a least-squares solution; column
scaling does the reverse. The contract page has this as a table.

## How it was reached

The first version found a number's binary exponent by asking the host to decode the
number. It passed every exact test and allocated 576 bytes for each round of the
compiled caller, 32 for each factor, where the library's rule is that a repeated call
allocates nothing. The factor is now found by dividing or multiplying by 2⁵¹², 2²⁵⁶, …, 2
in turn, ten exact steps down and ten up, in arithmetic the compiler keeps unboxed. The
bound on a factor, 2¹⁰²³, is the product of those ten powers.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`, with every earlier roster unchanged. For this batch:

- **Rectangular, exactly.** 173 matrices: named cases (a single entry, an identity, a
  matrix already in range, mixed units, a column that is small in every row, stored
  zeros, an empty row and column, no entries, one row, one column, the largest and
  smallest normal numbers together) and 162 seeded matrices up to 8 by 8 at three
  densities with magnitudes spread over 2⁶, 2⁸⁰ and 2⁶⁰⁰. In each: every factor is a power
  of two; every scaled value equals the exact rational product; every row and column with
  a nonzero has its largest scaled magnitude in [1/2, 1) and every other has the factor
  1; the reported maxima are the exact ones; the values passed in are unchanged; and
  equilibrating the scaled matrix again gives factors of exactly 1.
- **Symmetric, exactly.** 73 positive definite matrices: Gram matrices of integer columns
  scaled by powers of two up to 2¹⁰⁰ either way, orders 1 to 8, and one named case. In
  each: every factor is a power of two, every scaled value is the exact product, every
  scaled diagonal entry is in [1/2, 2), and every scaled off-diagonal entry is below 2 in
  magnitude.
- **The limits of the range.** A row whose largest magnitude is 2⁻¹⁰⁷⁰ gets the factor
  2¹⁰²³ and its column the rest. A diagonal entry of 2⁻¹⁰⁴⁰ gets 2⁵¹¹ and is not brought
  into range. An entry 2²⁰⁰⁰ times smaller than the largest in its row is rounded to zero
  when scaled, without a refusal, and stays in the pattern.
- **Refusals.** Each wrong length, each forbidden aliasing, each non-finite value or
  factor and each overflow from caller-chosen factors by its own typed error,
  with the outputs untouched where the contract says so. A zero and a negative diagonal
  return `DiagonalNotPositive` and write nothing; a column that does not begin with its
  diagonal and a matrix that is not square are refused by name.
- **Vectors.** The output of `scale-vector!` may be its input and may not be its scales.
- **Allocation.** A compiled Coalton caller making all eight calls: zero bytes over
  10,000 rounds after 10,000 warm-up rounds, at both policies.
- **The example.** A positive definite system whose diagonal spans 2⁻³⁹ to 2⁴¹, solved
  through its scaled form, returns its exact solution to thirteen digits.

**Planted source controls.** Each of these edits to `scaling.ct` was applied alone and the
focused tests run. Every planted source compiled, every one failed the tests, and the
restored source passed.

| Planted defect | Caught |
| --- | --- |
| A largest magnitude in [1, 2) is left there | yes |
| A largest magnitude may be scaled to exactly 1 | yes |
| A row or column with no nonzero is scaled by 0 | yes |
| The column pass reads the unscaled entries | yes |
| The row pass takes the first entry's magnitude, not the largest | yes |
| `scale!` leaves out the column factor | yes |
| `scale!` takes the first row's factor for every row | yes |
| `scale!` does not check its result | yes |
| Non-finite values are accepted | yes |
| The two scale vectors may be one array | yes |
| `row-maxima!` reports a signed entry | yes |
| `column-maxima!` reports a signed entry | yes |
| `scale-vector!` may write over its scales | yes |
| A scaled diagonal in [2, 4) is left there | yes |
| The symmetric factor falls by the whole power, not its root | yes |
| A zero diagonal is scaled | yes |
| A negative diagonal is scaled | yes |
| `scale-symmetric!` leaves out the column factor | yes |
| A column need not begin with its diagonal | yes |
| A refused diagonal still writes the scales before it | yes |

## Limits

- One pass each way. The property is the one stated; the scaling does not minimize a
  condition number, and for the rectangular case no near-optimality is claimed.
- A symmetric matrix that is not positive definite has no routine here. Scaling by the
  diagonal needs a positive diagonal; an indefinite matrix needs the rows-and-columns
  form applied symmetrically over several passes, which was not built.
- Scaling is not timed in this record. It is two passes over the stored entries.
- One host and finite rosters establish neither universal behaviour nor portability beyond
  the supported runtime.
