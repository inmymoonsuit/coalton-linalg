# Batch 17: damped least squares

State: completed; verified at safety 1 and safety 0 on 2026-11-23, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract is [damped.md](damped.md).

## What this is

The first item of a list a consumer sent of what a nonlinear solver built on this library
will need, with the order and the interface agreed with it before any code. Its loop
takes the damped least squares step

    (Jᵀ J + μ D) p = −Jᵀ r,        D = diag(Jᵀ J) + a floor,

accepts it if the residual's norm fell, and otherwise raises μ and tries again on the
same matrix. So it wants one analysis for a pattern, and then a factorization and a step
for each μ with nothing analyzed again. It asked for nothing else: it reads two norms of
the residual and decides for itself.

## What was built

`coalton-linalg/sparse/damped`, a handle in the manner of
[dependent rows](dependent-rows.md):

- `analyze` plans `Jᵀ J`, analyzes its factorization and allocates every workspace.
- `set-matrix!` takes the values of J for an iterate.
- `set-damping!` takes a positive diagonal; `set-column-damping!` takes the columns' own
  squared lengths plus a floor, and forms it again from each later matrix.
- `factor!(plan, mu)` and `step!(plan, r, out)`, any number of times in turn.
- Three status types, so that a damping, a factorization and a step that cannot be made
  are each refused by its own name, with the output unwritten.

Two things were proposed and withdrawn before they were built: a predicted reduction and
a step norm. The consumer's loop uses neither, and a quantity nobody reads has no
contract to be held to.

**The arithmetic is that of Batch 16 applied from the start.** The normal equations
square J, and a square leaves the range at half the exponent of what is squared. So J is
divided by the power of two at its largest entry and r by the one at its own; the
damping is held in the divided units; and the step is multiplied back by the two powers.
A step at any scale the format holds is the step at ordinary scale, bit for bit.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`.

- **Against exact arithmetic.** 337 steps on eight shapes from 2 by 5 to 7 by 7, seeded
  integer matrices, four values of μ and both kinds of damping, each compared with
  `−(Jᵀ J + μ D)⁻¹ Jᵀ r` in rational arithmetic from the stored doubles. The largest
  difference, relative to the step, is 4.9e-13; the bound asserted is 2^-30.
- **Across the range.** Eighteen cases with the matrix or the residual or both at powers
  from 2^-1020 to 2^1020 give the step at scale 1 times the residual's power over the
  matrix's, bit for bit. [The scale-sweep roster](rosters.md) has three entries for the
  package and holds the step over the whole sweep wherever its inputs are numbers.
- **The two-stage paths.** A step multiplied back by two powers of 2^520 whose product
  cannot be held, to −2^988 exactly; and a damping of the smallest numbers the format
  has, divided by a power whose square cannot be held, against the exact step.
- **Refusals.** Each value of each status is produced, with the plan's state and the
  output checked after it: before anything is set; dampings and floors that cannot be
  used; values of μ that cannot; a matrix of another pattern or shape; a matrix without
  full column rank at μ zero, with the column named; a column of zeros under a floor of
  zero, and under the smallest normal number at scale 1 and at 2^600; a step that cannot
  be held.
- **Allocation.** A compiled caller that sets the matrix and then factors and steps for
  four values of μ: zero bytes over 10,000 calls after warm-up.
- **Aliasing.** `step!` joins [the aliasing roster](rosters.md): with a square matrix the
  residual and the output may be one array, and the step is the same.
- **Planted source defects.** Twenty-nine edits to the new source, each applied alone,
  with the package's suite and both rosters run. Every planted source compiled, every one
  failed a test, and the restored source passed. This was the first run of them; none was
  missed.

| Planted defect | Number | Caught |
| --- | --- | --- |
| The matrix, a given damping, the floor of the columns' damping, or the residual is not divided | 4 | yes |
| A damping divided below the least normal number is not raised; a damping is always multiplied by the square of the power; the floor is not added to the columns' squared lengths | 3 | yes |
| Setting the matrix leaves a factor standing; a refused matrix leaves the earlier one set; a new matrix leaves the damping as it was divided for the old one | 3 | yes |
| A damping entry that is not positive, or a damping of the wrong length, is taken; setting a damping leaves a factor standing | 3 | yes |
| The columns' damping is taken with no matrix set, with a floor that is negative or not finite, or with a column of zeros under a floor of zero | 3 | yes |
| A factorization is tried with nothing set; a mu that cannot be used is taken; the damping is not added; a mu of zero is multiplied into the damping; a refused factorization leaves the earlier factor standing | 5 | yes |
| A step is taken with no factor, with a residual or an output of the wrong length, or with a residual that is not finite | 4 | yes |
| The step is toward the gradient; it is not multiplied back; it is always multiplied by the two powers' product; a step that cannot be held is written | 4 | yes |

## Limits

- **The step is through the normal equations**, so with μ zero its accuracy is that of a
  condition squared. The package judges no rank; a pivot at rounding level is factored.
- **A damping in the caller's units must be one the format holds.** The columns' own is
  formed in the divided units and has no such limit.
- **A damping entry below the least normal number, once divided, is raised to it.** The
  contract states what that changes and for which columns.
- **`J Jᵀ` and the augmented matrix have no damped form here.** No caller has asked.
- **The problems of the tests are small**, up to 8 by 5 and 7 by 7. Nothing is timed.
- **No predicted reduction, trust-region radius or step norm is computed.** A loop that
  needs them forms them from the step.
