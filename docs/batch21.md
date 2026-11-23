# Batch 21: an independent review, and what it found

State: **in progress.** Parts 1 and 2, the kept rows of a rank-deficient matrix and
dependent rows, are fixed and verified at safety 1 and safety 0 on 2026-11-23,
Windows/x86-64, SBCL 2.6.10 and the pinned Coalton compiler. The other parts are listed under [what is still open](#what-is-still-open)
with what was reproduced of each: **those defects are in this release.**

## What this is

Batches 15 to 20 were read by five independent reviewers, three of the code and two of
the mathematics, each given the sources, the contracts and the records as claims to test
and none of the author's conclusions. They read and traced by hand, and ran nothing of
the library. A consumer's own review of Batches 15 to 19, 6,316 checks against exact
rational arithmetic, had found no defect earlier the same day; its inputs were well scaled, and
most of what follows is at the edges of the floating-point range or in what a status says.

Every finding acted on here was first reproduced through the public interface, and every
fix has a test that was watched to fail on the sources as they were.

## Part 1: `coalton-linalg/sparse/rank-deficient`

What the central results are worth was also asked. A reviewer re-derived the residual
bounds, the gap bound, the certificate, both bounds of a step and the bound on a
vector's multipliers, and tried them on 46,000 exact or emulated cases, among them 6,000
at the certificate's limit. They held at ordinary scale. What did not hold:

| Finding | Before | Now |
| --- | --- | --- |
| The step's two bounds below the normal range | `2^1000 · I` of order 5 with a residual of `(1 + 2⁻¹⁵) · 2⁻⁶⁰`: both bounds 2^-1074, the step 1.118 of that from the exact one. The rounding of the step's own entries was in neither | With k entries written below the normal range, `(k + 1) / 2` of 2^-1074 is added to both |
| The lower bound on the largest singular value | For the 3 by 1 matrix of the smallest positive number: 2 of it, above the 1.732 that is true | A lower bound has 2^-1074 taken from it below the normal range, as an upper one has it added; `proved-lower-bound` too |
| A floor that is a singular value exactly | `[[5]]` asked at 5, and the identity at 1: `CertificateUnresolved`, from a pivot of zero | `CertificateKeptNotAboveFloor`: the trial is made again with a margin written |
| A floor whose product with the scale overflows | `[[2^-600]]` asked at 2^500: `CertificateUnresolved` | `CertificateKeptNotAboveFloor` |
| No row left out, and a lower floor of the smallest positive number | `CertificateGapAboveFloor`, with a gap bound of zero | Certified: zero is below every floor |
| A vector that holds every kept row, with a bound above the cut | `SupportsUndecided`, naming a kept row the vector does not hold. There is none. The same where no row is kept | Classed. The refusal is kept for a kept row that exists |
| A refresh that raises after its vectors were exchanged | Sixteen rows kept, sixteen copies changed, an entry budget of 120: `AnalysisStorageExceeded`, and then a gap bound of 1.0 and residual bounds of 0.25, proved for the old vectors, standing beside the new | Nothing proved stands from the exchange until the bounds are proved again |
| A typed error for another pattern | `reprove!` and `refresh!` withdrew every bound and the certificate, then raised | Raised first; what was proved of the values held stands |
| A step that cannot reach its bound | Refused after the factorization, the solves and the corrections | Refused on the certificate's two numbers, before anything is factored |

**A sentence that was false, with no defect in the code.** The contract, the design and
Batch 19's record said the dependent sets do not depend on which rows the search kept.
That is a theorem of supports that are exact. With a cut, a multiplier that is not zero
and is below the cut is classed out, and its size belongs to the kept rows: the rows
(1, 0), (0, 1), (1, 2⁻⁸) at a cut of 2⁻⁴ are a set of two and a row in no dependency in
that order, and one set of three as (1, 0), (1, 2⁻⁸), (0, 1). Each classing is proved as
stated. [The contract](rank-deficient.md#dependent-sets) now says what the sets are; the
design and Batch 19's record carry dated corrections; the tests hold both orders.

**Smaller corrections.** `predicted-rate` is `(√2 · 2⁻¹⁰ + η) / (1 − η)`: 2⁻¹⁰ bounded the
iteration matrix's two blocks and not its norm. `step-bound` is the sum of its three
terms to rounding and can exceed it by a few units of rounding. The drop of a search is
relative to the largest entry at a kept row, the vector's own 1 not counted. A refused
`class-supports!` withdraws a standing classing. A division by the matrix's power of two
rounds nothing for an entry within 2^1022 of the largest, and the contract now says so.
Six docstrings, a sentence of the interface index and one of an example disagreed with
the code.

### Verification of part 1

Both complete serial checks pass on the tree of this part, `scripts/check.ps1 -Safety 1`
then `-Safety 0`.

- **Each fix by a test seen to fail first.** The nine rows of the table: the inputs
  named, at orders 1, 2, 5 and 9 for the step and four shapes for the singular value's
  bound, each asserted in exact rational arithmetic.
- **Two tests the review asked for, of what was not broken.** The whole bound of a step
  where something is truncated: eight steps on matrices of rank two with one row kept,
  against the step of the matrix truncated to rank one, found from the eigenvalues of
  its two-by-two product with a root taken to eighty digits. The step is within 0.056 of
  its bound of that step. And the bound on a vector's multipliers where it decides
  something, the vectors a hundredth out of date: the exact distance to the least squares
  multipliers is at most 0.2 of the bound, on five vectors.
- **A test that asserted a refusal for a reason that was absent** is corrected: two kept
  rows a part in a million apart, both held by the vector, now class at a cut of 10⁻¹².
  A counter that the assertion above it never let move is gone.
- **Planted source defects.** Fifteen edits to the source, each taking one fix back out or
  breaking what a new test holds, each applied alone with the package's suite run. Every
  planted source compiled and every one failed a test: the nine fixes of the table, the
  bound proved rounded as it was, a vector holding no kept row taken to hold every one,
  the rate's constant, the whole bound without what truncation adds, and the multipliers'
  bound without the bound proved under it. The restored source passes.
- The whole suite of the package, its allocation checks among it: zero bytes a call after
  warm-up for the search's loop, the step, the refresh and the classing, as before.

### Left as it was, in part 1

- **`FindVectorRefused` and `StepFactorRefused` are still reached by no test.**
- **`undecided-vector` and `undecided-entry` are not reset by a later search**, against
  their docstrings' "the row count otherwise".
- **A ladder walked to a tiny floor** costs a factorization a rung, about 500 for a floor
  of 10⁻³⁰⁰ of the matrix's size on rows it cannot resolve.
- **`CertificateKeptNotAboveFloor` can in principle be returned with the counts given at
  every round**, where the margin never clears the factorization's error. Twenty thousand
  emulated trials did not produce it.
- **Entries more than 2^1022 below a matrix's largest are held as rounded.** A reviewer
  followed where that rounding can reach and found no false statement, and no proof that
  there is none.
- The style guide's rule for a private constructor's name, and three central
  abbreviations, are not changed here: they are the same in every package and are decided
  once, for all of them.

## Part 2: `coalton-linalg/sparse/dependent-rows`

- **A matrix whose stored values are all zero** was `RowsBadThreshold` under a positive
  threshold below about 10⁻¹⁵⁴, "so small beside the matrix's largest entry that its square
  cannot be held"; it has no largest entry, and `[10⁻³⁰⁰]` under the same threshold was
  `RowsFound`. It is now `RowsFound` with every row left out at every valid threshold. The
  test, seen to fail first, holds three shapes of zeros, one of them with no stored entry,
  at five thresholds.
- **"The division is exact" has its condition.** A reviewer could not break the bounds
  and noted that the written proof takes the divided matrix to be the matrix over a power
  of two, which holds only for entries within 2^1022 of the largest. The contract now says
  what the bounds are of: the matrix as held, entries further below than that rounded by
  at most half of 2^-1074 in the divided units. No covering argument for the stored
  matrix itself is claimed.

## What is still open

Reproduced through the public interface and **not fixed in this release**:

- **Conjugate gradients** reports `CgConverged` with a residual norm of 0 for a right-hand
  side whose entries are more than 2^1022 apart: `A = I`, `b = (2^100, 2^-980)`, tolerance
  0, returns `x = (2^100, 0)`. The same for `(10³⁰⁸, 10⁻¹⁷)`.
- **Least squares on an operator**: the same false `LsqrResidualMet` on that input;
  `LsqrNonFinite` for solvable problems whose iterates or vectors exceed 2^512 in the
  iteration's units, in one case with the answer written into x; a gradient the format
  holds reported as zero and called met, where one scale factor is small and the other
  large; a right-hand side orthogonal to the range refused under a damping above 2^500
  in the caller's units. Its operator norm estimate exceeds the Frobenius norm after four
  steps on a matrix with one dominant singular value, where the contract says it stays
  below until past the rank.
- **Damped least squares**: a damping entry raised to the least normal number changes
  the step under a large μ (−1/3 where −0.99999 is exact); a column 2^540 below the
  largest entry is reported as a column of zeros, or as rank deficiency; the columns' own
  damping with a floor of zero, followed into a later matrix's column of zeros, is refused
  at every μ; `refused-column` is not −1 after a refusal that was not of a pivot.
- **The two rosters** miss a call that returns nothing and holds its numbers in a handle,
  count any condition as a refusal, and hold three outputs to expectations that cannot
  fail.
