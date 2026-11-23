# Batch 18c: kept vectors refreshed

State: completed; verified at safety 1 and safety 0 on 2026-11-23, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract is
[rank-deficient.md](rank-deficient.md#refreshing-the-kept-vectors) and the design
[design-rank-deficient-step.md](design-rank-deficient-step.md).

## What this is

[Batch 18a](batch18a.md) kept, for each redundant row, the combination of kept rows that
gives it, and [18b](batch18b.md) took the step on them. A kept combination stays true only
while every row in it keeps its values. The consumer that asked for the step answered a
question put to it at 18b's landing: among its rows, those that say two things coincide
or are level are constant, and a distance is not, nor will the rows of tangency or of a
point on a curve be. A problem with a redundant distance on every spoke of a wheel would
have been refused by the certificate at every step of a drag, and searched again at 0.9 s
a step. So it asked for this before the null-space bases.

## What was built

`refresh!(kept, matrix, gap-target)` and `refreshed-count`, in
`coalton-linalg/sparse/rank-deficient`.

- Every bound is proved for the new values, as `reprove!` does.
- Each vector whose bound is above its share of the target is replaced by the row's
  multipliers against the kept rows for those values, found by a solve on the kept rows'
  augmented matrix and corrected until its residual is a sixteenth of the share.
- The new vector is stored with the search's drop tolerance and its bound proved as
  stored. The system of the rows left out is filled and factored again, and analyzed
  again only where a vector changed the rows it holds.

**Nothing is proved by the solves.** Each new vector's bound is proved by the product
bound for the vector as stored, so a solve that went badly gives a wide bound and a
refused certificate, never a false one. The solve's shift is therefore taken from a
bound proved at an earlier call, as a shift only.

**A vector is refreshed below the floor**, as the consumer asked: at its share of a
target the caller gives under the floor, `gap-target / (2 (d + 1)^½)`. The measurements
below give half the lower floor as the target, so with about 400 rows left out a vector
is replaced when its bound passes an eightieth of the floor.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`.

- **Against exact arithmetic.** Twenty-eight seeded matrices of small integers on seven
  shapes, each with a certificate and a step standing. Every row left out is multiplied
  by 1 + 2⁻¹², which keeps the rank and the span and makes each kept combination wrong by
  that much; the refresh replaces every vector, 93 in all; each new bound holds exactly
  for the new values and the vector as stored; the gap bound is back under 10⁻¹⁰; each new
  vector is the old one with its multipliers scaled, to nine digits; the certificate that
  would have been refused is issued; and the step is exact again to a bound of rounding
  size. Then one row alone is moved, on 24 of them: one vector is replaced and the others
  are unchanged bit for bit. Then every row is multiplied by a factor of its own, from
  1.125 to 1.5, so that the kept rows change too: 86 vectors replaced, every bound exact,
  the certificate issued and the step exact.
- **A vector that changes its rows.** A row that is one kept row, then the sum of two, then
  one again: its vector gains an entry and loses it, the system of the rows left out is
  analyzed again each time, and the bounds and the step are exact each time. A row that
  is one kept row and then another, beside a second row left out that stays the first:
  as many entries, at another row, and a system of another pattern.
- **A row that leaves the span.** Its vector becomes the best combination there is, its
  bound is its distance from the span, 0.5 to six digits, and the certificate refuses on
  the gap.
- **Refusals and state.** Before a search and after a refused one; four targets that
  cannot be used, which change nothing; with no certificate ever issued, where the bounds
  are proved and no vector replaced; another pattern or shape. With nothing out of date
  nothing changes, bit for bit. A refresh withdraws the certificate, the bound it proved,
  and the step.
- **Across the range.** The same refresh at five scales from 2^-900 to 2^900 gives the same
  vectors bit for bit and a gap bound times the scale;
  [the scale-sweep roster](rosters.md#what-is-declared) has an entry for it.
- **Allocation.** A compiled caller that alternates two matrices, refreshing two vectors,
  certifying and stepping for each: zero bytes over 10,000 calls after warm-up.
- **Planted source defects.** Twenty-four edits to the refresh's source, each
  applied alone, with the package's suite and both rosters run. Twenty-two failed a test
  on the first run. One was missed: with a vector that changed its row not seen to have,
  the test had a single row left out, whose system is one number whichever row its
  vector names. The test now has two rows left out that share a kept row and then do
  not, and the edit is caught. One did not compile on the first run, because it left a
  variable unused, and so showed nothing; written to compile, it is caught. The restored
  source passes.

| Planted defect | Number | Caught |
| --- | --- | --- |
| A refresh is made where no kept set stands; a target that is not positive and finite is taken | 2 | yes |
| A refresh leaves the certificate standing, or the step, or the bound the certificate proved | 3 | yes |
| A refresh does not check the pattern | 1 | yes |
| Which vectors: none is ever taken to be out of date; a vector's share is not in the divided matrix's units; every vector is replaced, out of date or not; a refresh is made with no shift to solve with | 4 | yes |
| The solve's factor is of the kept rows' earlier values | 1 | yes |
| The multipliers: the row itself is left out of the right-hand side; a vector is taken as good after its first solve; the iterate with the smallest residual is not the one kept; the entries are the multipliers and not minus them; a multiplier is stored at its number among the kept rows and not at its row; nothing is dropped | 6 | yes |
| The vectors afterwards: those not out of date are lost; the count replaced is not kept | 2 | yes |
| The new vectors are written and the old ones left in place | 1 | on the second run: the first edit did not compile |
| The system of the rows left out is not filled again; the bounds are not proved for the new vectors; the count replaced is read where no kept set stands | 3 | yes |
| A vector that changed its rows is not seen to have | 1 | on the second run, after the test was strengthened |

## What it costs

Timed once along two drags the consumer recorded from its own solver: a wheel of 400
spokes with a redundant distance on each, fully defined and loose, at six consecutive
steps of one drag, each step's matrix as the solver's step met it. Read in place and not
part of the repository; a machine otherwise idle, safety 1, one process. The search is at
step 0, at 10⁻⁶ of the upper bound on the largest singular value; the floors are 2 × 10⁻⁸
of the lower bound and 10⁻⁶ of the upper; **the refresh's target is half the lower floor,
so a vector is replaced when its bound passes 1.25 × 10⁻² of the floor.**

| | Loose wheel | Fully defined wheel |
| --- | --- | --- |
| Rows by columns, rows left out | 1,803 by 1,602, 399 | 2,004 by 1,602, 402 |
| The search, at step 0 | 0.88 s | 0.89 s |
| Without a refresh, at each later step | gap 12,000 times the floor; certificate refused | gap 3,900 times the floor; certificate refused |
| Vectors above their share, a step | 1 | 402 |
| Vectors that moved at all, a step | 1 | 402 |
| The refresh | 10.6 to 13.7 ms | 54 to 71 ms |
| for each vector replaced | | 0.13 to 0.18 ms |
| Gap bound after it, over the floor | 3 × 10⁻⁷ | 2 × 10⁻⁶ to 2 × 10⁻⁵ |
| The certificate after it | issued, 0.15 ms | issued, 0.19 to 0.26 ms |
| The step after it | taken, 1.0 to 1.2 ms | taken, 1.9 to 2.2 ms |
| The step's bound, over the step | 3 × 10⁻¹⁰ | 1 × 10⁻⁵ to 8 × 10⁻⁵ |
| A whole iteration | 11.8 to 15.0 ms | 56 to 73 ms |
| A new search at the last step | 0.88 s | 1.28 s |

**The refresh turns a search at every step into an iteration sixty times cheaper on the
loose wheel and about twenty on the fully defined one.** It is still ten and fifty times
what an iteration costs where no row changes, 1.2 to 1.3 ms on these wheels
([Batch 18b](batch18b.md#what-it-costs-and-what-the-diagnostics-read)).

Where the time goes:

- **About ten milliseconds of every refresh on these wheels is the system of the rows
  left out.** One dependency runs through the hub and touches every spoke, so `I + W Wᵀ`
  is dense at order 400, and it is factored again whole when one vector changes. That is
  the whole cost on the loose wheel, where one vector changes a step.
- **The rest is the solves**, one for each vector replaced with about three corrections at
  the solve's contraction of 2⁻¹⁰. Stopping a vector once its residual is a sixteenth of
  its share, and not at rounding, took the fully defined wheel's refresh from 105 to 122 ms,
  in a development run before this timing, to the times above.

Both were reported to the consumer before anything further was built. Agreed with it: the
system of the rows left out will be solved by conjugate gradients on `I + Wᵀ W`, whose
eigenvalues are at least 1 and whose products are sparse, as a batch of its own; and
nothing is built to make the solves faster until a measurement on real problems asks for
it.

## Limits

- **The kept set never changes.** A row that has left the span is the caller's signal to
  search again, given by the certificate's refusal on the gap after a refresh.
- **A refresh needs a certificate to have been issued once** for the kept set, for the
  shift of its solve. `RefreshNoShift` otherwise.
- **The system of the rows left out is factored again whole** whenever a vector was
  replaced. See the times above.
- **`RefreshFactorRefused` and `RefreshNonFinite` are produced by no test.**
- **A vector is refined only as far as its share needs**, so after a refresh the gap bound
  is well under the target and not at rounding.
- **The exact tests are small**, up to 13 by 10. The consumer's drags are timed and their
  certificates and steps reported, not checked against exact arithmetic.
