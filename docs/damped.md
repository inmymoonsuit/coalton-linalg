# Damped least squares

State: implemented in Batch 17. The package is `coalton-linalg/sparse/damped`;
[the batch record](batch17.md) has the evidence.

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
| DampedNotPositive | A computed pivot was not positive, and `refused-column` names its column. With μ zero this is a matrix without full column rank. No factor stands |
| DampedNonFinite | A damped diagonal entry or an intermediate result was not finite, [a damping the matrix's scale took out of range](#what-is-and-is-not-claimed) among them; no factor stands |
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
- **The size of the problem does not limit it.** J is divided by the power of two at its
  largest entry and r by the one at its own, so `Jᵀ J` and `Jᵀ r` are formed from numbers
  near 1, and the step is multiplied back. The divisions are exact: the step for a matrix
  and a residual at any scale the format holds is the step at ordinary scale times the
  residual's power over the matrix's, bit for bit.
  [The scale-sweep roster](rosters.md#what-is-declared) holds that from 2^-1020 to 2^1020.
- **One stated deviation.** A damping entry that the division takes below the least
  normal number is raised to it. That changes the matrix factored by at most
  `μ · 2^-1022` on the diagonal, in units of the largest entry squared, far below its
  rounding. It is what keeps the pivot of a column of zeros positive when its damping is
  only a floor: such a column's entry of the step is then exactly zero. In the caller's
  own units it moves that column's damping from the floor given to `2^-1022` times the
  square of the power of two the matrix was divided by. It matters only for a column that
  is exactly zero, and there both give a positive pivot and a zero entry of the step.
- **The three ways a step can fail stay apart by name**: a pivot that is not positive, a
  result that is not finite, a step that cannot be held. A loop that treats all three as
  "raise μ" can; the names are for diagnosis.
- **No predicted reduction and no step norm are computed.** A loop that accepts a step on
  the residual's norm alone needs neither, and that is the loop this was built for;
  [`norm2`](primitives.md) gives a norm that does not overflow early.
- **A damping given in the caller's units must be one the format holds.** For a matrix of
  entries near 2^600 the columns' squared lengths are not; `set-column-damping!` forms
  them in the divided units and has no such limit.
- **A damping that the matrix's scale takes out of range is refused at the factor, not
  where it was given.** The damping is carried in the divided matrix's units: a given
  entry, or the floor, times the square of the power of two the matrix was divided by.
  Under a small enough matrix that is not a number the format holds: a floor of 10⁻¹²
  under entries near 2^-700, or a given entry of 10³⁰⁰ under entries near 2^-600.
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

`set-matrix!` copies the values divided by `2^e`, e chosen so that the largest magnitude
falls in [1, 2), and refills the fixed-pattern product `Ĝ = Ĵᵀ Ĵ`. The damping is held in
the same units, `D̂ = D / 2^(2e)`. `factor!` copies `Ĝ`, adds `μ D̂` to its diagonal, and
factors. `step!` divides r by the power of two at its largest entry, forms `ĝ = Ĵᵀ r̂` down
each column, solves `(Ĝ + μ D̂) p̂ = −ĝ`, and returns `p = p̂ · 2^f / 2^e`, f the
residual's exponent: substituting, `(G + μ D) p = −Jᵀ r` exactly as stated. The two
powers are multiplied together first, which is exact, unless that product leaves the
range while the step need not; then they are applied one after the other.

## Using it

```
plan = analyze(J)                       once for the pattern
each iterate:
    set-matrix!(plan, J)
    set-column-damping!(plan, floor)    once; it follows the matrix after that
    for mu in the loop's own sequence:
        factor!(plan, mu)               DampedNotPositive: raise mu
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
