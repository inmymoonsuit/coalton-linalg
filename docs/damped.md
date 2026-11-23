# Damped least squares

State: implemented in Batch 17 and reworked in Batch 21, after a review, to divide each
column by its own power of two. The package is `coalton-linalg/sparse/damped`;
[the first record](batch17.md) and [the second](batch21.md) have the evidence.

## What it is for

A nonlinear least squares loop linearizes its residual, `r + J p`, and asks for a step
`p` that reduces it without leaving the region where the linearization can be trusted.
The damped step is the solution of

    (Jᵀ J + μ D) p = −Jᵀ r,        D a positive diagonal, μ ≥ 0.

With μ zero it is the least squares step. As μ grows the step shortens and turns toward
the scaled steepest descent direction, `−D⁻¹ Jᵀ r / μ`. A loop tries a μ, judges the step,
and tries another: on one matrix, many values of μ.

This package is the machinery for that: one analysis for the pattern of J, the values of
J for an iterate, a damping, and then any number of factorizations and steps, each for
its own μ, with nothing analyzed again and nothing allocated. **It chooses no μ and
judges no step.** The policy of the loop, when to raise μ and when to accept, is the
caller's.

## Interface

- `analyze(matrix) -> DampedPlan`, `analyze-with-limit(matrix, maximum-entries)`: once for
  a pattern. It plans `Jᵀ J`, analyzes its factorization and allocates every workspace;
  no value is read. The errors are those of [the Gram analysis](gram.md) and [the factor
  analysis](spd-ldlt.md); the budget applies to each.
- `set-matrix!(plan, matrix)`: the values of J for this iterate. The matrix must have the
  analyzed pattern; another pattern or shape raises the typed error of `refill-csc!` and
  leaves no matrix set. A standing factor does not survive.
- `set-damping!(plan, d) -> DampingStatus`: a positive diagonal, one entry a column, in
  the units of the matrix squared. It stays through later `set-matrix!` calls.
- `set-column-damping!(plan, floor) -> DampingStatus`: `D = diag(Jᵀ J) + floor`, each
  column's squared length plus a floor in the same units, formed again from every later
  matrix. A matrix must be set.
- `factor!(plan, mu) -> DampedStatus`: the factorization of `Jᵀ J + μ D` for the standing
  matrix and damping. Each call overwrites the factor.
- `step!(plan, r, out) -> StepStatus`: `out = −(Jᵀ J + μ D)⁻¹ Jᵀ r` with the standing
  factor. `r` has a row's length and `out` a column's; for a square matrix they may be
  one array.
- `matrix-set?`, `damping-set?`, `factor-ready?`, `refused-column`, `row-count`,
  `column-count`, `default-analysis-entry-limit`.

| DampingStatus | Meaning |
| --- | --- |
| DampingSet | The diagonal stands; a standing factor does not |
| DampingBadSize | The diagonal's length differs from the column count; nothing changed |
| DampingNotPositive | An entry was not positive or not finite; or the floor was negative or not finite; or a column of zeros met a floor of zero. Nothing changed |
| DampingNoMatrix | The diagonal is formed from the matrix and none is set; nothing changed |

| DampedStatus | Meaning |
| --- | --- |
| DampedSuccess | The factor stands |
| DampedNotPositive | A computed pivot was not positive, and `refused-column` names its column. With μ zero it is what a matrix without full column rank gives, and what one of full rank gives when the normal equations lose its smallest singular value to rounding. No factor stands |
| DampedZeroDamping | The damping is the columns' own with a floor of zero, μ is positive, and the standing matrix has a column of zeros: its damping entry is zero, and no μ makes its pivot positive. `refused-column` names it. No factor stands |
| DampedNonFinite | A damped diagonal entry or an intermediate result was not finite, [a damping a column's scale took out of range](#what-is-and-is-not-claimed) among them; no factor stands |
| DampedBadDamping | μ was negative or not finite; no factor stands |
| DampedNotReady | No matrix is set, or no damping; no factor stands |

| StepStatus | Meaning |
| --- | --- |
| StepSuccess | The output holds the step |
| StepNoFactor | No factor stands; the output is unchanged |
| StepBadSize | A length is wrong; the output is unchanged |
| StepNonFiniteRhs | An entry of the residual was not finite; the output is unchanged |
| StepNonFiniteResult | The step cannot be held; the output is unchanged |

Setting the matrix or the damping withdraws the factor: `factor-ready?` is true only
while the factor is that of the standing matrix and damping. One plan serves one call at
a time.

## What is and is not claimed

- **The step is that of the stated matrix**, computed by the
  [positive-pivot factorization](spd-ldlt.md) and its solve, with their backward error.
  The tests compare it with the step in exact rational arithmetic on 337 cases, with
  both kinds of damping.
- **The normal equations square the condition of J.** With μ zero and a J of condition κ
  the step is accurate to about κ² rounding errors, as for any solve through `Jᵀ J`. A
  positive μ bounds the condition of what is factored by `(‖J‖² + μ d_max) / (μ d_min)`.
- **The sizes in the problem do not limit it, within two stated conditions.** Each column
  of J is divided by the power of two at its own largest entry and r by the one at its
  own, so `Jᵀ J` and `Jᵀ r` are formed from numbers near 1 whatever the columns' sizes, and
  each entry of the step is multiplied back by its column's power and the residual's.
  Columns 2^900 apart are each held in full; a column far below the largest entry is not
  a column of zeros. The divisions are by powers of two, and so the step for a matrix, a
  column or a residual multiplied by a power of two is the step at ordinary scale, scaled,
  bit for bit, wherever nothing leaves the range.
  [The scale-sweep roster](rosters.md#what-is-declared) holds that from 2^-1020 to 2^1020
  for the columns' own damping, and from 2^-510 to 2^510 for a damping given, which has to
  be scaled as the matrix squared. The two conditions:
  - **An entry more than 2^1022 below the largest of its own column is held as rounded**,
    by at most half of 2^-1074 in the column's divided units.
  - **An entry of the residual more than 2^1022 below its largest is not seen.** One power
    of two cannot carry both, and the part of the step that entry alone would give is lost:
    for the identity and `r = (2^600, 2^-500)` the step is `(−2^600, 0)`.
- **μ times a damping entry is formed once, at the factor**, as one product of μ, the
  entry and the square of its column's power, in an order whose intermediates the format
  holds. A damping far below the matrix squared under a μ large enough to make it count is
  held in full: `J = [2^520]`, `D = [1]`, `μ = 2^1023` gives the step `−1/(1 + 2^-17)`.
- **One stated deviation: the pivot of a column of zeros is 1.** Where a column's stored
  entries are all zero and μ times its damping entry is positive, the matrix factored
  differs from `Jᵀ J + μ D` at that one diagonal entry and nowhere else, and the step is
  not affected. A column of zeros makes row and column j of `Jᵀ J` zero, so the system
  decouples there: entry j's equation is `pivot · p_j = −(Jᵀ r)_j = 0`, and no other entry
  reads that pivot. Any positive pivot gives the same step, bit for bit, with entry j
  exactly zero, and 1 is the one that cannot leave the range: no floor and no μ can be
  refused, or lost, on such a column's account. Where its damping entry is
  zero, which is the columns' own damping with a floor of zero, `factor!` answers
  `DampedZeroDamping` for every positive μ; and with μ zero its pivot is zero, which is
  `DampedNotPositive`.
- **The three ways a step can fail stay apart by name**: a pivot that is not positive, a
  result that is not finite, a step that cannot be held. A loop that treats all three as
  "raise μ" can; the names are for diagnosis.
- **No predicted reduction and no step norm are computed.** A loop that accepts a step on
  the residual's norm alone needs neither, and that is the loop this was built for;
  [`norm2`](primitives.md) gives a norm that does not overflow early.
- **A damping given in the caller's units must be one the format holds.** For a matrix of
  entries near 2^600 the columns' squared lengths are not; `set-column-damping!` forms
  them in the divided units and has no such limit.
- **A damping that a column's scale takes out of range is refused at the factor, not
  where it was given.** μ times a damping entry is wanted in its column's divided units:
  times the square of the power of two that column was divided by. Under a small enough
  column that is not a number the format holds: a floor of 10⁻¹² under a column of entries
  near 2^-700, or a given entry of 10³⁰⁰ under one near 2^-600. Such a damping is more
  than 2^1023 times the column's own squared length, and the entry of the step it would
  give is below anything the format resolves beside the others.
  `set-damping!` and `set-column-damping!` still answer `DampingSet`, because the damping
  stays through later matrices and may be in range for the next, and `set-matrix!`
  refuses nothing. `factor!` with μ positive then answers `DampedNonFinite`, and `step!`
  `StepNoFactor`. With μ zero the damping is not used and the factor stands. The tests
  hold both cases.
- **`DampedNotPositive` is a finding about a computed pivot.** With μ zero on a matrix
  that is rank deficient only to rounding, the factorization can succeed with a pivot at
  rounding level and the step is then large and meaningless; the package does not judge
  rank. [Dependent rows](dependent-rows.md) does.

## The method

`set-matrix!` copies column j divided by `2^(e_j)`, e_j chosen so that the column's
largest magnitude falls in [1, 2), and refills the fixed-pattern product `Ĝ = Ĵᵀ Ĵ`, which
is `C G C` with C the diagonal matrix of the `2^(−e_j)`. `factor!` copies `Ĝ`, adds to
diagonal entry j the damping in the same units, `μ D_j / 2^(2 e_j)`, and factors. `step!`
divides r by the power of two at its largest entry, forms `ĝ = Ĵᵀ r̂` down each column,
solves `(Ĝ + μ C D C) p̂ = −ĝ`, and returns `p_j = p̂_j · 2^f / 2^(e_j)`, f the residual's
exponent: substituting `p = 2^f C p̂`, `(G + μ D) p = −Jᵀ r` exactly as stated. The two
powers are multiplied together first, which is exact, unless that product leaves the
range while the entry need not; then they are applied one after the other.

Dividing a symmetric matrix on both sides by powers of two, and its right-hand side with
it, changes no significand anywhere in a factorization whose order is fixed by the
pattern: every pivot, multiplier and solved entry comes out the same number times a power
of two. So where one power for the whole matrix was enough, as it was before Batch 21,
the step is bit for bit what it was. [The record](batch21.md) holds 9,000 steps to that.

## Using it

```
plan = analyze(J)                       once for the pattern
each iterate:
    set-matrix!(plan, J)
    set-column-damping!(plan, floor)    once; it follows the matrix after that
    for mu in the loop's own sequence:
        factor!(plan, mu)               DampedNotPositive: raise mu; DampedZeroDamping: none helps
        step!(plan, r, p)
        judge p; accept or try another mu
```

[The example](../examples/damped-step.ct) takes the least squares step and a damped one
on one analysis.

## Representation choices

A handle, like [dependent rows](dependent-rows.md), and not free functions over a
caller's Gram values: the division by a power of two has to be applied to the matrix,
the damping, the residual and the step consistently, and a caller composing the Gram
refill and the factorization itself would have to carry it. The plan keeps the damping
as given beside the divided copy, so that a diagonal in the caller's units survives a
change of matrix.

`J Jᵀ` and the augmented matrix are not offered in damped form; no caller has asked.
