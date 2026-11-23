# Batch 18d: a refresh on a vector's own rows, and a factor kept

State: completed; verified at safety 1 and safety 0 on 2026-11-23, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract is
[rank-deficient.md](rank-deficient.md#refreshing-the-kept-vectors) and the design
[design-rank-deficient-step.md](design-rank-deficient-step.md).

## What this is

[Batch 18c](batch18c.md) made a drag through redundant rows that change possible, and
measured two costs it left: ten milliseconds at every refresh for the system of the rows
left out, factored again whole, and a solve on all the kept rows for every vector
replaced. This batch removes each where it was not needed. It was planned as one change
and became three, because two measurements made before any code said so.

**Plain conjugate gradients were measured and not built.** The agreed plan was to stop
factoring the system of the rows left out and iterate on it. Replayed on the consumer's
patterns, the iteration reaches a residual of 10⁻¹⁴ in 3 steps on the loose wheel and 12
on the fully defined one, and in 56 to 58 on the lattice of 21, where a multiplier of
100 puts the largest eigenvalue at 3 × 10⁴. There the sparse factor costs some tens of
microseconds a solve and the iteration would cost about a millisecond, on rows that
never change. So the factor is kept and used as a preconditioner, which costs nothing
where nothing changed.

**A refreshed vector was two hundred entries of noise.** Counted while measuring: after
a refresh the fully defined wheel's vectors held 81,858 entries where the search's held
1,602. A solve on all the kept rows, stopped at a fraction of the vector's share, leaves
small multipliers at every kept row, well above the relative drop, and they are part of
what makes the residual small, so they cannot be dropped. A vector that is still a
combination of the rows it holds does not need that solve at all.

## What was built

In `coalton-linalg/sparse/rank-deficient`:

- **A refresh on a vector's own rows.** For a vector of at most 16 rows the multipliers
  are found against those rows alone, by their normal equations and one correction, and
  taken when the residual is under the fraction of the share asked for, or at rounding.
  The solve on all the kept rows is made only for a vector this does not serve, and its
  factorization only at the first such vector.
- **The factor of the rows left out's system kept.** A refresh that changes no vector's
  rows leaves it, and the step's solves of that system become conjugate gradients
  preconditioned with it, at most `system-iteration-limit` iterations, 8, to a residual of
  2⁻⁴⁰ of the right-hand side. A solve that has not converged by then makes the factor
  again and is one solve. With a factor of the vectors as they stand, the code that runs
  is the code of 18b.
- **`refresh-to!`**, with the fraction of its share a replaced vector is corrected to,
  which the consumer asked for so that it can tighten only when its own comparison would
  otherwise refuse; `refresh!` is the same at a sixteenth.
- **`step-norm-term` and `step-outside-term`**, the two terms `step-bound` adds to
  `step-solve-bound`; and `system-iterations`, `system-factor-count`.

`RefreshNoShift` and `RefreshFactorRefused` change their meaning with the first: they now
say that a vector needed the solve on all the kept rows and it could not be had, the
vectors that did not need it having been replaced.

**Not built, and why.** The consumer asked that a search whose system of the rows left
out exceeds the entry budget refuse by a name of its own. On the patterns where that
system is dense, a hub, the search's own product of the rows is dense too and exceeds
the same budget first, at the analysis, as the typed error that has been there since
Batch 11; no input was found that reaches the other. It is left as the typed error
`find-kept-rows!` is documented to raise.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`.

- **The first route, against exact arithmetic.** The twenty-eight matrices of 18c, now
  refreshed with no certificate ever issued, which the first route does not need: every
  vector replaced, every bound exact, the gap bound under 10⁻¹⁰. Then as in 18c: the rows
  left out moved, one row moved, every row moved; 93, 24 and 86 vectors.
- **The second route.** Twenty matrices with no zero entry whose rows left out are
  multiples of one kept row and become sums of two or three: the rows a vector holds no
  longer serve. With no certificate issued, 19 of the 20 are `RefreshNoShift`, with the
  vectors that needed the solve left and the certificate refused on the gap; with one,
  every vector is brought under its share, 19 gain rows, and the bounds, the certificate
  and the step are exact.
- **The factor as a preconditioner.** On the twenty-eight, the steps after a refresh made
  137 preconditioned iterations and no factorization. On four matrices with twelve rows
  left out, each multiplied by a different power of two, the iteration does not converge
  in eight steps, the factor is made again in all four, the step is exact, and the next
  step makes no iteration and no factorization.
- **A factor of the vectors as they stand.** Every test of 18b and 18c that takes a step
  without a refresh between runs the code of 18b; a second step on one certificate is the
  first bit for bit, as there.
- **The stop target and the terms.** Four fractions that cannot be used change nothing;
  0, a sixteenth and 1 are taken. The whole bound is the sum of its three terms to
  rounding (corrected 2026-11-23: it can exceed that sum by a few units of rounding, the
  first term being returned without the margin the whole carries; the test allows 10⁻¹²)
  and at least 0.99 of it; the terms read −1 when no step stands.
- **The rosters.** [The scale sweep](rosters.md#what-is-declared) answers for `refresh-to!`
  and the two terms, at the ranges already declared.
- **Allocation.** The loops of 18b and 18c, refreshing, certifying and stepping: zero
  bytes over 10,000 calls.
- **Planted source defects.** Nineteen edits to the new source, each applied alone,
  with the package's suite and both rosters run. Every planted source compiled and every
  one failed a test, on the first run; the restored source passes. Three tests were
  added while the edits were being chosen, because an edit had nothing to fail: that a
  factor a part in four thousand out of date is not made again, that a fraction of zero
  still takes a vector on its own rows at rounding, and four matrices on which the
  factor must be made again.

| Planted defect | Number | Caught |
| --- | --- | --- |
| The few-rows fit: it is taken whatever it leaves; it is not taken at rounding when no fraction is given; its multipliers are written at their order and not at their rows; the other kept rows are not cleared first; its normal matrix is filled in the other triangle; it is found and not used | 6 | yes |
| `refresh-to!`: any fraction is taken; the kept rows' solve is stopped far above the fraction; a vector that needed a shift and had none is not reported | 3 | yes |
| The system of the rows left out: a factor of earlier vectors is taken for one of the vectors as they stand; an iteration that did not converge is taken as converged; an iteration is never seen to converge; its operator is the identity; the factor is made again from the values it had; a factor that did not serve is not made again; a factorization made again is not counted; a step's iterations are added to the last step's | 8 | yes |
| The step's terms: one is read where the other is; one is not in the caller's units | 2 | yes |

## What it costs

Timed once along the two drags of [Batch 18c](batch18c.md#what-it-costs), with the same
floors, search and target, on a machine otherwise idle, safety 1, one process. Each wheel
was run twice: with a replaced vector corrected to a sixteenth of its share, the default,
and corrected as far as it goes, a fraction of zero.

**The loose wheel**, where one vector of 399 changes a step and keeps its rows:

| | Batch 18c | This batch |
| --- | --- | --- |
| The refresh | 10.6 to 13.7 ms | 0.07 to 0.09 ms |
| The certificate after it | issued, 0.15 ms | issued, 0.15 to 0.33 ms |
| The step after it | taken, 1.0 to 1.2 ms | taken, 1.5 to 1.8 ms |
| A whole iteration | 11.8 to 15.0 ms | 1.7 to 2.1 ms |
| Preconditioned iterations a step | | 6, over the step's three solves |
| Factorizations of the rows left out's system | one at every refresh | the search's, and none after |
| Entries the vectors hold | | 1,193, as the search left them |
| The step's bound, over the step | 3 × 10⁻¹⁰ | 8 × 10⁻¹¹ |

The two fractions give the same figures: on its own rows the vector comes out at
rounding either way. **An iteration there now costs what it costs where no row changes**,
1.2 ms on this wheel in [Batch 18b](batch18b.md#what-it-costs-and-what-the-diagnostics-read),
plus half a millisecond for the six iterations.

**The fully defined wheel**, where all 402 vectors change at every step:

| | Corrected to a sixteenth of the share | Corrected as far as it goes |
| --- | --- | --- |
| The refresh | 66 to 70 ms | 81 to 121 ms |
| for each vector replaced | 0.16 to 0.18 ms | 0.20 to 0.30 ms |
| Gap bound after it, over the floor | 1 × 10⁻³ to 2 × 10⁻³ | 1 × 10⁻⁶ to 2 × 10⁻⁵ |
| The certificate, the step | issued, 0.2 ms; taken, 1.9 to 2.0 ms | issued, 0.2 ms; taken, 1.9 to 5.3 ms |
| A whole iteration | 68 to 73 ms | 87 to 123 ms |
| Entries the vectors hold | 80,880 to 81,780 | 81,858 to 81,984 |
| Factorizations of the rows left out's system | one a step: a vector changed its rows | at three steps of five |

The step's bound there, in its terms, each over the norm of the step:

| | A sixteenth of the share | As far as it goes |
| --- | --- | --- |
| The solve | 9 × 10⁻⁹ | 9 × 10⁻⁹ |
| The term in the step's norm | 1 × 10⁻⁵ to 3 × 10⁻⁵ | 2 × 10⁻⁸ to 2 × 10⁻⁷ |
| The term in the residual outside the range | 7 × 10⁻³ to 9 × 10⁻³ | 9 × 10⁻⁶ to 8 × 10⁻⁵ |
| The whole bound | 7 × 10⁻³ to 9 × 10⁻³ | 9 × 10⁻⁶ to 8 × 10⁻⁵ |

**What this wheel says.**

- **Its rows are not combinations of a few rows.** The kept rows span every column, so
  each row left out is a combination of all of them; of the three its vector held at the
  search it is one only to 4 × 10⁻¹⁰, a hundred times its share, at the values the drag
  passes through. So the first route serves about twenty-five vectors of 402 and the
  second the rest, which take two hundred rows each. That is the matrix, not the method,
  and the cost of a refresh stays what 18c measured.
- **The bound is the term in the residual outside the range**, the gap bound times that
  residual over the square of the bound proved on the kept rows. Correcting further
  lowers the gap bound a hundredfold and the bound with it, for half as much time again;
  it does not go lower: 10⁻⁵ of the floor is where the proved bound of vectors of two
  hundred entries stands, their rounding allowed for. **What is left to tighten is the
  bound proved on the kept
  rows**, which the ladder leaves up to four times below their smallest singular value
  and which enters squared.
- **The default fraction leaves the gap bound a hundred times looser than 18c did**,
  still five hundred times and more under the floor: the vectors the first route takes sit just
  under their sixteenth, where the second route's overshoot by decades. A caller whose
  own comparison needs the tighter bound asks for it with `refresh-to!`.

## Limits

- **A row that is not a combination of a few kept rows stays expensive.** Where the kept
  rows span every column, as on the fully defined wheel, each row left out is a
  combination of them all, and of its own three only to about 4 × 10⁻¹⁰ at the values
  met along the drag, a hundred times a vector's share. Its vector then takes two hundred
  rows, by the second route, at every step.
- **The first route is for vectors of at most 16 rows**, and uses the normal equations of
  those rows; rows that are nearly dependent among themselves send it to the second.
- **`RefreshFactorRefused` and `RefreshNonFinite` are produced by no test**, and neither
  is a search refused for the budget of the rows left out's system.
- **How far to correct is not chosen for the caller.** At a sixteenth of its share a
  vector taken by the first route can sit just under it; the gap bound is then a
  hundredth of the floor's share and not a millionth. `refresh-to!` is the answer.
