# Batch 18b: the step on a rank-deficient matrix

State: completed; verified at safety 1 and safety 0 on 2026-11-23, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract is
[rank-deficient.md](rank-deficient.md#the-step) and the design
[design-rank-deficient-step.md](design-rank-deficient-step.md).

## What this is

The second of the two parts of a step on a rank-deficient matrix.
[Batch 18a](batch18a.md) kept the rows, proved their bounds and certified the rank. This
batch takes the step on what was certified: the minimum-norm least squares step of the
matrix with its redundant rows put back as their kept combinations, and a proved bound on
how far the step written is from the one a singular value decomposition cut at the
certified floor would give.

The solve was measured on the consumer's seven patterns before any code, and the
measurement changed one agreed term, with the consumer's agreement:

**The bound is read from the residuals, and no rate of convergence is proved or used.** A
refusal had been agreed for a solve whose predicted contraction, with the factor's
backward bound in it, was not below 1. On the largest pattern that prediction is 0.9
where the contraction measured is 2 × 10⁻⁴: the backward bound is a bound through a trace,
grows with the order, and would have refused solves that reach rounding in two
corrections. For any pair the error of the kept rows' solve is an exact expression in its
two residuals, so the bound needs the residuals, proved, and the bound the certificate
proved at the same call, and nothing else. The prediction is kept as a number a caller
can read.

## What was built

In `coalton-linalg/sparse/rank-deficient`:

- `step!(kept, r, out)` writes `−B⁺ r` for the matrix the standing certificate is of, and
  refuses by one of six names, writing nothing.
- `step-solve-bound` and `step-bound`: what the two solves left, and the distance to the
  truncated matrix's step.
- `step-corrections`, `correction-limit`, `step-taken?`; and three diagnostics,
  `solve-growth`, `solve-backward-bound`, `predicted-rate`.
- `find-kept-rows!` now also analyzes and factors the system of the rows left out,
  `I + W Wᵀ`, since W changes only when the vectors do.

Two things were settled while building.

**The limit on corrections is eight, not the six first written.** A correction contracts
by 2⁻¹⁰ or better while it is above rounding, so six exhaust the precision; on small
matrices whose smallest singular value is near the bound proved, five were needed. Two
are left for corrections that halve the bound by chance at the bottom, so that a
refinement that has finished is not refused as unfinished.

**The system of the rows left out is refined by the same rule**, not by one correction.
With a multiplier of 2²⁰ the inversion lemma loses digits to cancellation, and one
correction left more of it than several do.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`.

- **Against exact arithmetic.** Eighty steps on forty seeded matrices of small integers on
  eight shapes, from 3 by 4 to 13 by 10 and one with no redundant row, each against a
  residual anywhere and a residual in the range. The step written is compared with the
  exact minimum-norm least squares solution, in rational arithmetic, twice: for B, formed
  from the kept vectors' stored doubles, against `step-solve-bound`; and for the matrix
  itself against `step-bound`. Both hold exactly on every case, and the whole bound is at
  most 1.4 × 10⁻¹² of the sizes of the step and the residual together.
- **Kept vectors out of date on purpose.** On the thirty-five of those matrices that leave
  a row out, every row left out is then multiplied by 1 + 2⁻¹², which keeps the rank and
  the row space and makes each kept combination wrong by that much; the bounds are proved
  again, the certificate is asked at twice the gap, and seventy steps are taken. Both
  bounds hold exactly on all seventy, and on sixty-eight the step is further from the
  exact step than its solve bound allows: the distance to the truncated matrix's step is
  doing the work it is there for.
- **By hand.** Two equal rows against `[1, 3]` give −2, minus the mean, where the kept row
  alone would give −1. Rows `[1; 1.01]` with a vector found when they were equal, against
  `[1, −1]`: the step for B is zero, the exact step is not, and the bound covers it by the
  residual outside the range alone.
- **Across the range.** Nine pairs of scales for the matrix and the residual, from 2^-1000
  to 2^900 each: the step is the step at scale 1 times the residual's power over the
  matrix's, bit for bit, and the bounds are those times the same, compared from both
  sides where they fall below the normal range.
  [The scale-sweep roster](rosters.md#what-is-declared) has two entries for the step.
- **Refusals and edges.** Each status an input reaches is produced, with the output
  checked untouched and the handle's state after it: no certificate, before a search,
  after one, after a refused certificate and after a reproof; four sizes; three residuals
  that are not finite; `4γ/τ ≥ 1`, with a gap bound of 0.3 certified under a proved 0.35;
  a step of 2²⁰⁰⁰. A residual of zeros; no row left out, on more columns than rows; a row
  with no entry and a row of stored zeros; one entry, with eight quotients the format
  cannot hold; a matrix of zeros; two rows a part in a million apart; multipliers of
  2¹⁰ to 2²¹. A second step on one certificate is the first
  bit for bit; a step does not survive a certificate, a reproof, a search or a refused
  step, and neither does the factor it solved with.
- **Allocation.** A compiled caller that proves the bounds again, certifies and takes the
  step twice: zero bytes over 10,000 calls after warm-up.
- **Aliasing.** `step!` joins [the aliasing roster](rosters.md): on a square matrix with a
  row left out, the residual and the output may be one array and the step is the same.
- **Planted source defects.** Forty edits to the step's source, each applied alone,
  with the package's suite and both rosters run. Every planted source compiled.
  Thirty-seven failed a test on the first run. Of the three that did not, two showed a
  gap and are caught now: no test gave a residual that was too long, only ones too
  short; and zeroing the rounding allowance of one of the kept rows' two residuals is
  covered by the other's, so that edit became both at once, and eight steps on one
  entry whose quotient the format cannot hold were added, where a residual computes to
  exactly zero and the error is not zero. **The third no test fails**: the term that
  carries the error of y into the solve's bound. In every case tried, with multipliers
  up to 2²¹, y is accurate far inside what the kept rows' own residuals already allow.
  The term stays because the proof needs it. The restored source passes.

| Planted defect | Number | Caught |
| --- | --- | --- |
| The system of the rows left out: the identity is not added to `W Wᵀ`; there is taken to be nothing to solve whenever a row is left out | 2 | yes |
| What withdraws a step: a search, a reproof or a certificate leaves the step and the solve's factor standing; the withdrawal leaves the factor marked, or the step | 5 | yes |
| Refusals: a step with no certificate; an output too long; a residual that is not finite; a refused step leaves the earlier one standing | 4 | yes |
| A residual that is too long is taken | 1 | on the second run, after the test was added |
| The solve's factor: its bottom block is its top block; its top block is 1 and not the bound proved | 2 | yes |
| `y = C⁺ r`: the residual is not divided; the residual on the rows left out does not reach it; W or `Wᵀ` is applied with the other sign; the residual of y leaves out the rows left out's own, or takes `W y` with the other sign | 6 | yes |
| The error of y does not reach the solve's bound | 1 | **no**: nothing observable |
| `x = J_S⁺ y`: either product in the residual has the other sign; the bound is the residual and not the residual over the bound proved; the iterate with the smallest bound is not the one kept; a correction is never added | 5 | yes |
| The rounding of the kept rows' two residuals is not allowed for | 1 | on the second run, as one edit of both, after the tests were added |
| The truncated matrix: its distance is left out of the whole bound; a gap bound a quarter of the bound proved is taken; the residual outside the range, or the step's own norm, is left out; the residual over every row is taken to be zero | 5 | yes |
| Units: the step has the sign of `B⁺ r`; it is not multiplied back by the residual's power, or the matrix's; a step that cannot be held is written; either bound is not in the caller's units | 6 | yes |
| A bound, or the count of corrections, is read when no step stands | 2 | yes |

## What it costs, and what the diagnostics read

Timed once on the consumer's seven patterns, read in place and not part of the repository,
on a machine otherwise idle, safety 1, one process. The search is at 10⁻⁶ of the upper
bound on the largest singular value, the drop at 10⁻¹², the floors at 2 × 10⁻⁸ of the
lower bound and 10⁻⁶ of the upper. Two residuals, one with every entry seeded in
[−1, 1] and one in the range of the matrix. A later step, on the factor a first step
left, is the mean of twenty; a whole iteration is a reproof, a certificate and a step,
also the mean of twenty.

| Pattern | Left out | Search, vectors and the left-out rows' system | A first step | A later step | A whole iteration | Corrections |
| --- | --- | --- | --- | --- | --- | --- |
| wheel, well constrained | 2 | 0.24 s | 0.8 ms | 0.7 to 0.8 ms | 1.1 ms | 3 |
| wheel, redundant | 402 | 0.88 s | 1.1 ms | 1.0 ms | 1.3 ms | 3 |
| wheel, redundant and loose | 399 | 0.84 s | 1.1 to 1.5 ms | 1.0 to 1.2 ms | 1.2 ms | 3 |
| lattice 15, redundant | 814 | 0.59 s | 2.3 ms | 2.0 ms | 3.1 ms | 3 |
| lattice 15, redundant and loose | 812 | 0.46 s | 2.1 to 2.2 ms | 1.8 ms | 2.9 ms | 3 |
| lattice 21, redundant | 1,642 | 2.64 s | 4.8 to 5.1 ms | 3.9 to 4.0 ms | 6.8 ms | 3 |
| lattice 21, redundant and loose | 1,640 | 2.02 s | 4.3 to 5.0 ms | 3.5 to 4.2 ms | 6.6 ms | 3 and 4 |

Every step was taken. **The solve's own factorization is the difference between a first
step and a later one: a tenth to four tenths of a millisecond on the wheels and the
lattice of 15, and about one millisecond on the lattice of 21.** The rest of a step is
the system of the rows left out, the first solve, three corrections and the proved
residuals. Putting that system into the search added a few hundredths of a second to it.

The bounds, over the norm of the step written, and the diagnostics of the factor; ρ₀ is
2⁻¹⁰ on every pattern, and η and the predicted rate are the terms the first design would
have refused on:

| Pattern | Solve bound, residual anywhere; in the range | Whole bound, anywhere; in the range | Growth | Backward bound | η | Predicted rate |
| --- | --- | --- | --- | --- | --- | --- |
| wheel, well constrained | 5e-8; 4e-7 | 6e-8; 4e-7 | 3e4 | 5e-6 | 0.09 | 0.10 |
| wheel, redundant | 6e-8; 4e-7 | 6e-7; 4e-7 | 3e4 | 5e-6 | 0.09 | 0.10 |
| wheel, redundant and loose | 1e-10; 2e-11 | 6e-9; 1e-10 | 2e3 | 2e-7 | 2e-4 | 1e-3 |
| lattice 15, redundant | 8e-9; 8e-8 | 1e-6; 8e-8 | 2e7 | 7e-6 | 0.11 | 0.12 |
| lattice 15, redundant and loose | 7e-12; 4e-11 | 2e-8; 8e-10 | 4e6 | 2e-6 | 0.006 | 0.007 |
| lattice 21, redundant | 9e-9; 6e-8 | 2e-6; 6e-8 | 2e7 | 3e-5 | 0.48 | 0.91 |
| lattice 21, redundant and loose | 7e-12; 4e-11 | 3e-8; 1e-9 | 4e6 | 5e-6 | 0.02 | 0.02 |

**The predicted rate is 0.91 on the lattice of 21, and its refinement stalled at rounding
after three corrections**, as on every other pattern: a correction contracts by decades,
not by a tenth. The prediction is a bound, and a loose one.

The bounds are wider than the steps are wrong. They are dominated by the rounding allowed
for in the residuals, over the bound proved, which the ladder leaves up to four times
below the smallest singular value; and for a residual outside the range, whose step is
thousands of times longer here, the whole bound also carries that residual over the
square of the bound proved.

## Limits

- **The bound on the rows left out's system widens with their multipliers.** It is the
  norm of a residual, and the residual's own rounding grows with the square of the
  largest multiplier. Over the step's norm the solve bound was 9 × 10⁻⁹ at a multiplier
  of 2¹⁰, 4 × 10⁻⁵ at 2¹⁶ and 10⁻² at 2²⁰, while the step itself stayed far more
  accurate than that. On the consumer's patterns the multipliers are of order 1.
- **`StepFactorRefused` and a refinement still halving at its eighth correction are
  produced by no test.** A certificate bounds the multipliers by the matrix's size over
  the bound proved, far below where their squares cannot be held; and no matrix tried
  contracts slower than a half a correction.
- **The rounding allowances are not each tested alone.** One planted defect zeroes the
  kept rows' two together and is caught; one carries the error of y and is not; the
  others are covered only by the bounds holding exactly on every case.
- **A kept vector that has gone out of date is not refreshed.** The certificate refuses
  and the caller searches again. For a matrix whose redundant rows change at every call
  that is a search at every call; refreshing one vector from current values is designed
  and not built.
- **The kept rows alone, as an option, are not built**; the consumer will not call it.
- **The exact tests are small**, up to 13 by 10. The consumer's patterns are timed and
  their bounds reported, not checked against exact arithmetic.
