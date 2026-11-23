# Batch 21: an independent review, and what it found

State: **completed.** Parts 1 to 6, the kept rows of a rank-deficient matrix, dependent
rows, damped least squares, conjugate gradients, least squares on an operator and the two
rosters, are fixed and verified at safety 1 and safety 0 on 2026-11-23,
Windows/x86-64, SBCL 2.6.10 and the pinned Coalton compiler. No defect the review found in the numerical packages is known to remain
open; what is left is [the style rules](#what-is-still-open).

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

## Part 3: `coalton-linalg/sparse/damped`

Three of the findings had one cause: the matrix was divided by one power of two, the one
at its largest entry, and then squared. A column far below that entry had a square that
could not be held, and a damping far below the matrix squared was moved to the least
normal number so that it would not vanish.

**Each column is now divided by the power of two at its own largest entry**, and μ times
a damping entry is formed once, at the factor, in its column's units.

| Finding | Before | Now |
| --- | --- | --- |
| A damping far below the matrix squared, under a large μ | `J = [2^520]`, `D = [1]`, `μ = 2^1023`: the step −1/3, where −1/(1 + 2⁻¹⁷) is exact. The entry was raised to 2^-1022 of the matrix squared, which μ multiplied | Exact to rounding. Nothing is raised |
| A column far below the largest entry | `diag(1, 2^-540)`: `DampingNotPositive`, "a column of zeros", for a floor of zero, and `DampedNotPositive` with μ zero. Neither is true | Both succeed, and the step is exact; columns 2^900 apart also |
| The columns' own damping with a floor of zero, followed into a later matrix with a column of zeros | `DampedNotPositive` at every μ, which a loop answers by raising μ | `DampedZeroDamping`, a new status: no μ helps, and `refused-column` names the column |
| `refused-column` after a refusal that was not of a pivot | The column of an earlier refusal, or 0 | −1 |
| A column of zeros under a floor, or a damping, near the top or bottom of the range | `DampedNonFinite`, or a pivot lost | Its pivot is 1: the one stated deviation, since any positive pivot gives a zero entry of the step there |

**The gate: nothing that was right moved.** Dividing a symmetric matrix on both sides by
powers of two changes no significand in a factorization whose order the pattern fixes.
That was held to, not assumed: 9,000 solves, on forty seeded matrices of ten shapes, eight
of them with a column of zeros, at five scales of the matrix from 2^-500 to 2^500 and
three of the residual, with a damping given, the columns' own, and the columns' own with
a floor, at five values of μ from 0 to 10⁸. Every outcome and every entry of every step
was written out bit for bit before the change and after. The two files are identical:
SHA-256 `167b3cfe2eb1627913745e8fbfdeeaff6d219542b01b733c1fc7369135ae4bc7` of 9,000
lines, both times. The package's exact tests report the same worst difference as before,
4.872 × 10⁻¹³ on 337 steps, and the scale sweep the same 925 runs.

**Stated, not changed.** An entry more than 2^1022 below the largest of its own column
is held as rounded. (Part 3 as first landed also stated, and did not refuse, a residual
entry that far below the residual's largest; [it is refused now](#part-3-continued-an-entry-lost-inside-the-solve).) `DampedNotPositive` with μ zero is described
as a finding about a computed pivot, which a matrix of full rank can also give. The
roster's page said the step with a given damping was held everywhere; the roster holds
it from 2^-510 to 2^510, where the damping can be scaled as the matrix squared.
Batch 17's "relative to the step" was relative to the larger of 1 and the step.

### Verification of part 3

Both complete serial checks pass on the tree of this part.

- **Each fix by a test seen to fail first**: fourteen diagonal problems whose steps are known
  exactly, with a damping given and with the columns' own, the numbers powers of two and
  not, each entry held to a relative 2⁻⁴⁸; `refused-column` after each kind of refusal; a
  column of zeros under each damping at μ and floors from 10⁻³⁰⁰ to 10³⁰⁰; and the new
  status at four values of μ, with the factor and the step that follow it.
- **The gate above.**
- **Planted source defects.** Eleven edits to the source, each applied alone with the package's suite run. Every planted source compiled. Ten failed a test at once: the step or the damping divided by another column's power, the damping entry divided before μ meets it whatever that leaves, the zero column's pivot not made positive, its zero damping not named, `refused-column` kept from an earlier refusal or not recorded, the floor not divided as its column is, the columns' own damping without the column's square, and a μ of zero adding its damping. The eleventh did not fail, and could not: it changed a branch that the value it tested for never reaches. Replaced by an edit that puts the old raise back where the division takes an entry below the normal range, it fails. The restored source passes.
- The aliasing roster and the scale sweep pass with their counts unchanged, and the
  allocation check of the package's loop reads zero bytes a call.

### Part 3, continued: an entry lost inside the solve

The consumer reviewed part 3 and found one more thing of the same family, which the
first change had left standing. The residual is divided by one power of two, and under a
large damping the solve's unknown is smaller than its right-hand side by the damping:
for the identity with `D = (1, 2^700)`, `μ = 1` and `r = (2^200, 2^-200)` the second
unknown was 2^-1100 in those units and came back as zero, under `StepSuccess`, where the
step's second entry is −1.18 × 10⁻²⁷¹.

**The damped matrix is now balanced at the factor**: multiplied on both sides by a
diagonal of powers of two, the power at one over the root of each diagonal entry, so
that what is factored has a diagonal between 1 and 4 whatever the damping and an unknown
is of the size of its right-hand side. Each entry of the step is multiplied back by its
three powers of two in the order that keeps the partial products in range.

**What that still cannot carry is refused by name, with a new status, `StepOutOfRange`,
and the output as it was**: a residual, or a gradient in the units of the solve, whose
entries lie more than 2^1022 apart, seen by an entry not coming back from its division as
it went in; and an entry of the solve below the normal range where the step's own entry
is not. The first of those had been a stated condition in part 3, under which the
identity and `r = (2^600, 2^-500)` returned `(−2^600, 0)` with `StepSuccess`; it is the
refusal now. A number that was given below the normal range, and a step entry that is
itself below it, are not refused.

The gate again: the same 9,000 solves, written out bit for bit, give the same SHA-256 as
before either change, `167b3cfe…4bc7`. The exact tests' worst difference, the scale
sweep's 925 runs and the zero bytes a call are as they were.

Tests seen to fail first: the consumer's probe; three more diagonal problems with numbers
that are not powers of two; and a matrix of two columns on separate rows against the step
in rational arithmetic. With the fix: three steps refused by the new status, with the
output untouched; a residual entry given below the normal range, solved exactly; and a
residual spanning 2^1000, solved. Planted source defects: nine edits, each applied alone with the package's suite run. Eight compiled and failed a test at once: the matrix not balanced, balanced on one side only, the gradient not brought to the solve's units, the step not multiplied back by the balance, each of the three refusals not made, and a number given below the normal range refused. The ninth, the three powers applied in the order they come, did not compile as first written, an edit that left a binding unread; repaired, it fails a test. The restored source passes.

## Part 4: `coalton-linalg/conjugate-gradient`

- **A false `CgConverged`.** The iteration is carried on the residual divided by a power
  of two. That power was chosen once, from the first residual, and every residual formed
  from the operator to confirm a tolerance was divided by it too. An entry more than
  2^1022 below the first residual's largest was lost to that division both times: for the
  identity and `b = (2^100, 2^-980)` at a tolerance of zero the call returned
  `(2^100, 0)`, `CgConverged`, and a residual norm of 0, where the residual is 2^-980. So
  for `(10³⁰⁸, 10⁻¹⁷)`. It was a regression of Batch 16, which brought the division in;
  [that record](batch16.md) now says so. **The power is now chosen again at every
  residual formed**, which is where the method starts again anyway. The residual formed
  is then measured at its own size, the tolerance is tested on that, and the next steps
  solve what the first power could not carry: five right-hand sides with entries up to
  2^1900 apart are solved to a residual of exactly zero.
- **A residual whose norm is above the largest number is not a refusal.** With every
  entry finite the norm is infinity, meets no finite tolerance, and the method goes on, as
  it already did for a running norm: `b = (1.5 × 10³⁰⁸, 1.5 × 10³⁰⁸)` under the identity
  was `CgNonFinite` and is solved.
- **Nothing else moved**: the ninety-six solves that the tests recompute are bit for bit
  what they were, the choice of a power of two changing no significand, and the scale
  sweep's three entries for the two calls pass as declared.
- The contract's "the size of the system does not limit it" has its condition, and no
  longer says the preconditioned iterates are identical at 2^1020, where the roster
  declares them within the tolerance.

Two planted defects, each putting one of the two behaviours back, both caught; the
restored source passes. Both complete serial checks pass on the tree of this part.

## Part 5: `coalton-linalg/lsqr`

The newest package, Batch 20, had been reviewed by no one but its author. Two reviewers
read it, one the code and one the mathematics. The recurrences, both identities for the
estimates, the reflections, the minimum norm and the end of the sequences were re-derived
and checked at eighty digits and in rational arithmetic, and held. What did not:

| Finding | Before | Now |
| --- | --- | --- |
| A right-hand side with entries more than 2^1022 apart | The identity, `b = (2^100, 2^-980)`, tolerances of zero: `(2^100, 0)`, `LsqrResidualMet`, a residual of 0. The residual was formed from the right-hand side as divided, which had lost the small entry | `LsqrExhausted`, a residual of 2^-980 and a gradient of 2^-980, both exact. The residual is formed from b as given, entry by entry, and measured at its own power of two; a second call on it finds the rest |
| A norm whose sum of squares is above the largest number | `LsqrNonFinite` on problems with nothing non-finite in them: a column 2^-512 of the other, in one case with the exact answer already written into x; a row 2^512 above the other; `A = [1; 0]` with `b = (2^-512, 1)`, of condition number 1 | Measured again with the entries lowered, as a sum too small was already measured again with them lifted. All are solved. The iterate's norm is formed only under a damping |
| A gradient brought to the caller's units one factor at a time | Operator at 2^-500, right-hand side at 2^500: a gradient of 2^-601 reported as 0 and `LsqrGradientMet`, where the same problem at scale 1 is `LsqrLimit` with 2^-601. The running estimates late in a run were zeroed the same way, and bought confirmations | The two powers together where their product is a number the format holds. One problem at three scales gives one history, entry for entry, and the same confirmations |
| The gradient's two parts in one fixed unit | Not reached before, but the new residual made it reachable: an image and a damped term of very different sizes | Formed in the units of the larger, the other brought to them |
| A right-hand side the transpose sends to zero, under a large damping | `[1; 0]`, `b = (0, 3)`, damping 8 is `LsqrGradientMet`; all three times 2^500 was `LsqrOutOfRange` | `LsqrGradientMet` at every scale: zero solves it whatever the damping, and the ceiling is not asked |

**Three sentences that were false, with the code as it was.**

- The operator norm estimate was said to stay below the Frobenius norm "while the method
  has taken no more steps than it would in exact arithmetic". It exceeds it as soon as a
  singular value has converged and is counted again: at the fourth step of forty on a
  matrix with one singular value of 1000 and the rest between 1 and 2, and by a factor of
  3.46 at the fortieth. The test that held the claim used matrices with no dominant
  value. The contract says what is true, and a test holds the counterexample.
- "Neither the residual without the damped term nor the gradient is monotone." The
  residual is: the objective does not rise and the norm does not fall. The monotone test
  now measures it, on all 466 iterates.
- "A wrong estimate does not cost two products a step." A tolerance below what rounding
  lets the gradient reach does: 35 confirmations in forty steps at 5 × 10⁻¹⁵ on a matrix
  whose gradient stops at 7 × 10⁻¹⁵. That is measured, stated, and not changed.

### Verification of part 5

Both complete serial checks pass on the tree of this part.

- **Each fix by a test seen to fail first**, every norm compared with the exact one as
  squares in rational arithmetic, so that no small one is lost on the way: the five rows
  of the table at twenty inputs, among them the second call that finds what the first
  could not see, and a damped problem whose iterate has an entry of 2^520.
- **The suite of Batch 20 passes as it was**: 120 problems against exact arithmetic with
  the same worst distance, 4.507 × 10⁻¹¹, every status, the scale tests, zero bytes a
  call. The recomputation in the test language was brought to the changed iteration and
  agrees bit for bit throughout.
- **Planted source defects.** Ten edits to the source, each applied alone with the package's suite run, each putting one of the old behaviours back or breaking one part of the new residual and gradient. Every planted source compiled and every one failed a test at once. The restored source passes. The thirty-two edits of Batch 20 were not
  made again: the lines most of them changed are as they were.
- The aliasing roster and the scale sweep pass with their counts unchanged.

### Left as it was, in part 5

- **A tolerance below what rounding allows costs a confirmation at nearly every step**, as
  above. Nothing recognizes a stall.
- **The powers of two are still the right-hand side's and the first image's.** Where b is
  all but orthogonal to the operator's range the iterate is held far below its unit and
  can lose digits, and a problem can be `LsqrOutOfRange` that another choice would have
  served: `[1; 0]`, `b = (2^-501, 1)` under a damping of 1.
- **The last step allowed makes a transpose product the iterate does not need**, and a
  value in it that is not finite ends the call with a finished iterate.
- **A tolerance of infinity meets a norm of infinity.**
- Three of eight ways to `LsqrNonFinite` have a test; the others are reached by reading.

## Part 6: the two rosters

The rosters are guards: each claims that every exported call of a family is answered
for. A reviewer attacked the guards themselves.

- **Any condition counted as a refusal.** A pair declared refused passed if the aliased
  call died of anything at all, and a run declared refused at the top of the range passed
  for any error there. A refusal is now a status that is not the one needed, or the
  library's typed error; any other condition is not caught and ends the roster's run.
  Every pair and every run that was a refusal still is one: 18 pairs, and the same
  ranges on all 37 entries of the sweep.
- **An output that could not fail.** Six norms were held only to a bound that the status
  before them had already checked. Five are now held to their scaled values, bit for bit:
  the residual of conjugate gradients in two entries, and the residual and the gradient
  of least squares on an operator in three. They hold at every scale their entries
  range. The sixth, in the preconditioned entry, has no scaled value at 2^1020, where the
  iterates are declared only close, and keeps its bound.
- **`:lost` took any finite number**, at the top of the range as well as the bottom. It is
  now zero or a number below the normal range, where the right value is below it; above
  the range a finite number is `:wrong`.
- **The sources were read two directories deep**, where the page said every depth. Every
  depth is read.

Controls, in the rosters' own files as theirs are: the typed error is a refusal in each
roster; a plain error in a run or in an aliased call is not caught; five numbers on
either side of the range are classed lost or wrong. Both rosters pass with their counts
unchanged, which is itself the finding that no entry had rested on the looser rules.

**Left as it is, and now written on [the rosters' page](rosters.md#membership-is-read-from-the-sources):**
a call that takes its numbers inside a stored matrix or a handle and returns no status
is not a member of the sweep (`set-matrix!` of damped least squares, `transpose` and
`permute`); an array parameter is recognized by how its type is written; a pair is of
the entry's own fixture; the one exemption rests on today's caller.

## What is still open

The style guide's rule that a private single constructor is named with a leading `%`
(every handle here is `Mk...`), `if` forms whose branches act or bind, three central
abbreviations, a parameter order, and constants defined in two places: a reviewer listed
them against the style rules. They are the same in every package and are to be decided
once, with the style guide's author, and not package by package in a batch of fixes.
