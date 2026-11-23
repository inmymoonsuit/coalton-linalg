# Kept rows and the certificate

State: implemented in Batches 18a, the kept rows and the certificate, 18b, the step, and
18c and 18d, the refresh of kept vectors, and 19, the dependent sets. The package is
`coalton-linalg/sparse/rank-deficient`. [The design](design-rank-deficient-step.md) has
the mathematics and the measurements that shaped it. The batch records,
[18a](batch18a.md), [18b](batch18b.md), [18c](batch18c.md), [18d](batch18d.md) and
[19](batch19.md), have the evidence.

## What it is for

A matrix with redundant rows is rank deficient, and a solver that takes many steps on
one such matrix, its values changing a little each time, needs three things at every
step, cheaply and **holding as stated or refused by name**:

- the rows the others depend on, and for each other row the combination that gives it;
- a bound on how well each combination holds for the values of this step;
- a proof that no singular value lies in a band about the floor the solver cuts at.

[Dependent rows](dependent-rows.md) finds the first once. This package keeps what was
found, proves the second again for new values without finding anything, issues the
third, and then takes the step: the minimum-norm least squares step of the matrix with
its small singular values cut off, with a proved bound on how far the step written is
from it.

## Interface

- `analyze(matrix) -> KeptRows`, `analyze-with-limit(matrix, maximum-entries)`: once for
  a pattern. No value is read.
- `find-kept-rows!(kept, matrix, threshold, drop) -> FindStatus`: the search of
  [dependent rows](dependent-rows.md#interface) at a distance threshold, and then the
  vector of every row left out. In each vector an entry below `drop` times the vector's
  largest is set to zero, and the bound is proved for the vector as kept. **This is the
  expensive call**, a millisecond or so a row left out, and it allocates. It is made
  once for a kept set that will serve many matrices.
- `reprove!(kept, matrix) -> ReproveStatus`: for a matrix of the same pattern with other
  values, every kept vector's bound proved again against them. Nothing is searched, no
  vector changes, nothing is allocated. Another pattern or shape raises the typed error
  of the Gram refill.
- `refresh!(kept, matrix, gap-target) -> RefreshStatus`: `reprove!`, and then every kept
  vector whose bound is above its share of the target replaced by the row's multipliers
  for these values. `refreshed-count(kept)` is how many the latest call replaced.
  `refresh-to!(kept, matrix, gap-target, fraction)` is the same with the fraction of its
  share a replaced vector is corrected to; `refresh!` uses `default-refresh-fraction`, a
  sixteenth.
- `certify!(kept, lower-floor, upper-floor) -> CertificateStatus`: that the gap bound is
  below the lower floor and every singular value of the kept rows above the upper, both
  floors in the matrix's own units.
- `proved-lower-bound(kept)`: what the standing certificate proved of the kept rows, at
  least the upper floor asked for and usually far above it.
- `residual-bound(kept, position)`, `gap-bound(kept)`: the bounds as last proved.
- `kept-count`, `left-out-count`, `left-out-rows`; `multiplier-start`,
  `multiplier-count`, `multiplier-rows`, `multiplier-values` for the vectors as kept;
  `kept-found?`, `certified?`, `row-count`, `column-count`.
- `step!(kept, r, out) -> StepStatus`: the step against a residual, for the matrix the
  standing certificate is of. Nothing is allocated.
- `step-solve-bound(kept)`, `step-bound(kept)`: the two bounds of the standing step, −1
  when none stands; `step-norm-term`, `step-outside-term`: the two terms that the second
  adds to the first. `step-taken?`, `step-corrections`, `correction-limit`.
- `system-iterations`, `system-factor-count`, `system-iteration-limit`: what the system
  of the rows left out cost the standing step, and since the search.
- `class-supports!(kept, cut) -> SupportStatus`: each multiplier classed as in its row's
  dependency or out of it, and every row labelled with its dependent set;
  `dependent-sets`, `dependent-set-count`, `support-flags`, `multiplier-error-bound`,
  `supports-classed?`, `undecided-vector`, `undecided-entry`.
- `solve-growth`, `solve-backward-bound`, `predicted-rate`: diagnostics of the factor the
  steps solve with, −1 when there is none.
- `largest-singular-bounds!(matrix, steps, work-columns, work-rows, out)`: a lower and an
  upper bound on the largest singular value of any stored matrix, which is what a floor
  given as a fraction of it needs. The three arrays must be distinct.

| FindStatus | Meaning |
| --- | --- |
| FindFound | A kept set and the vectors of the rows left out stand |
| FindNonFinite | An intermediate result of the search, or a bound, was not finite; nothing stands |
| FindUnderflow | The search's `RowsUnderflow`; nothing stands |
| FindBadThreshold | The search's `RowsBadThreshold`; nothing stands |
| FindBadDrop | The drop tolerance was not above 0 and below 1; nothing stands |
| FindVectorRefused | The vector of a row left out could not be computed; nothing stands |

| ReproveStatus | Meaning |
| --- | --- |
| Reproved | Every bound stands for the values given |
| ReproveNothingKept | No kept set stands |
| ReproveNonFinite | A bound was not finite; the vectors stand, the bounds and any certificate do not. A later reproof against values that can be bounded restores them |

| CertificateStatus | Meaning |
| --- | --- |
| Certified | The statement below is proved |
| CertificateGapAboveFloor | The gap bound is not below the lower floor: a singular value beyond the kept count may reach it |
| CertificateKeptNotAboveFloor | The kept rows were not proved above the upper floor, by a factorization accurate at that floor |
| CertificateUnresolved | The proof's own error exceeds the upper floor: the floor is below what this method resolves |
| CertificateBadFloor | A floor was not positive and finite, or the lower exceeded the upper |
| CertificateNothingKept | No kept set stands, or its bounds do not |

A certificate is of the values the bounds were last proved for, and does not survive a
search or a reproof, refused or not. A refused search leaves no kept set. One handle
serves one call at a time.

`ReproveNonFinite` needs a kept vector with entries near 2^510 and values that give them
something to multiply; `FindVectorRefused` and a bound that is not finite at the search
itself have not been produced by any input tried.

## The vectors as kept

For a row t left out the search gives a vector with 1 at t and minus its multipliers at
kept rows. **Most of what is computed there is rounding.** On the patterns this was
measured on, a row depends on three or four others and its computed vector has an entry
at nearly every kept row the factor reaches, all but those three or four below 10⁻¹³
(the [design](design-rank-deficient-step.md#what-was-measured) has the counts). So an
entry below `drop` times the vector's largest is set to zero. The tolerance is relative
to each vector by itself, because the multipliers of one row can combine rows of
different units; 10⁻¹² is the value the measurements suggest.

**Dropping costs no proof anything.** The bound returned is proved for the vector as
kept, whatever was dropped from it, so every statement below is about the vectors a
caller reads through `multiplier-rows` and `multiplier-values`. A larger tolerance gives
shorter vectors and larger bounds, never a false one; the tests cut vectors at 10⁻³,
at a half and at 0.99 and hold every bound against exact arithmetic.

`multiplier-values` holds the vector's own entries at kept rows. A row's multipliers,
the coefficients that give it from the kept rows, are those entries with the sign
changed.

## What is proved

Write `J_S` for the kept rows and r for their number.

**Each kept vector is within its bound.** For the vector v at a position,
`‖Jᵀ v‖₂ ≤ residual-bound`, the product exact, for the values the bounds were last proved
for. The bound is formed as [dependent rows](dependent-rows.md#what-is-proved) forms
its own, with one difference: only the rows the vector touches are read, and only the
columns those rows reach are summed, the others being exactly zero. So a reproof costs
in proportion to the vectors' entries and their rows' lengths, not to the matrix.

**At most r singular values exceed the gap bound.** `gap-bound` is the root of the sum of
the squared bounds, rounded upward, and bounds singular value r + 1 of the matrix by the
argument dependent rows gives.

**`Certified` means: `gap-bound < lower-floor`, and every singular value of `J_S` exceeds
`proved-lower-bound ≥ upper-floor`.** Since deleting rows cannot raise a singular value,
the r largest of the matrix exceed that bound too. So exactly r singular values are above
the upper floor and the rest are below the lower: none lies between.

A row with no entry, or of stored zeros, is left out with a vector of its own entry alone
and a bound at the bottom of the format. For a matrix of zeros no row is kept: the
statement about the kept rows is then empty, and the certificate is the gap bound's alone.

The second half is [the rank certificate from the inertia](quasidefinite.md#a-rank-certificate-from-the-inertia),
applied to the kept rows with t = s = τ: one signed factorization of
`[[τ I, J_Sᵀ], [J_S, τ I]]`, its pivot counts and its measured backward bound, with the
margins that page requires. τ is tried on a short ladder, from the bound that served the
last call or from an eighth of the power of two at the matrix's largest entry, down by
fours, and last at the floor itself; the first that is proved is `proved-lower-bound`.
So the bound proved is a lower bound and not an estimate: it can be a factor of four, or
more at the top of the ladder, below the smallest singular value. A caller that keeps
certifying one kept set pays for one factorization a call.

**A bound is proved only with room for the factorization's own error.** The matrix
factored is shifted by more than the bound asked about, by at least the factorization's
measured backward bound, so that kept rows whose smallest singular value equals a bound
to rounding are never certified at it.

**The two refusals of the kept rows are different statements.**
`CertificateKeptNotAboveFloor` is returned when the factorization at the floor was
accurate, its backward bound under an eighth of the floor, and did not give the counts:
the rows are at or below the floor, or within the proof's small error of it.
`CertificateUnresolved` is returned when the proof's own error is not small beside the
floor, and says nothing about the rows. Two rows 4 × 10⁻⁸ apart have a smallest singular
value near 2 × 10⁻⁸: asked for a floor of 10⁻⁶ the answer is the first, and for a floor of
10⁻¹⁰, which is true of them, the second. The proof resolves about a ten-millionth of
the matrix's size.

**The search's threshold and the floors are different numbers.** The search is run at a
tolerance well inside the gap between the singular values that are kept and those that
are not. Run at the floor itself it works at its own resolution and can keep a dependent
row, and the certificate then refuses. Certifying at a floor below the search's
threshold is safe: a row left out whose distance lies between the two makes the gap
bound at least the floor.

**Everything is on the matrix divided by a power of two**, the one at its largest entry.
The rows left out, the vectors and the statuses at any scale the format holds are those
at ordinary scale, and the bounds are those times the scale;
[the scale-sweep roster](rosters.md#what-is-declared) holds that from 2^-800 to 2^1020.

## Refreshing the kept vectors

A kept vector says how one row is made from the kept rows, and it stays true only while
every row in that combination keeps its values. Rows that say two things coincide do;
a row that is a distance does not, and a vector through one is out of date at the next
step. `reprove!` then shows a large bound and the certificate refuses on the gap.

`refresh!(kept, J, gap-target)` does what `reprove!` does and then **replaces each vector
that no longer serves**, from the current values:

1. The values are taken and every bound proved, as by `reprove!`.
2. A vector is out of date when its bound is above its share of the target,
   `gap-target / (2 (d + 1)^½)` with d the rows left out. With every bound under its share
   the gap bound is under half the target. **Give a target below the lower floor the
   certificate will be asked at**: at the floor itself a vector would be replaced only
   once the certificate had already refused. Half the floor is what the measurements here
   use, which replaces a vector at about an eightieth of the floor when 400 rows are left
   out.
3. For each such row a, its multipliers are found again, by one of two routes.
   **On the rows its vector already holds**, where there are at most 16 of them: the least
   squares multipliers against those rows alone, by their normal equations and one
   correction. This serves whenever the row is still a combination of the same rows, and
   costs microseconds. **On all the kept rows** otherwise, `w = (J_S J_Sᵀ)⁻¹ J_S aᵀ` by
   their augmented matrix: one factorization, made at the first vector that needs it, and
   for each a solve and corrections.
   Either is taken when its residual is under `fraction` of the vector's share, a
   sixteenth unless `refresh-to!` is given another. A fraction of zero takes the first
   route only at rounding, and corrects the second until a correction no longer halves
   the residual.
4. The new vector is stored as a search stores one, with the search's own drop tolerance,
   and **its bound is proved for the vector as stored and these values**.

**Nothing is proved by the solves, and nothing needs to be.** A solve that went badly
gives a vector with a wide bound, and the certificate refuses; it cannot give a false
one. For that reason the second route's shift is taken from the bound a certificate
proved at an *earlier* call: it is a shift and nothing more, no bound of this call rests
on it, and the next step makes its own factor from the bound proved at its own call.
The first route needs no shift, and so no certificate to have been issued.

**How far to correct is a trade the caller makes.** The first route leaves a residual
near rounding where the row is exactly a combination of its rows, and just under the
fraction where it is only nearly one. The second usually overshoots the fraction by
decades. A smaller fraction buys a smaller gap bound, and with it a smaller bound on
the step, for more corrections. **A row that is not, at these values, a combination of
a few kept rows to within its share takes whatever rows the second route gives it,
which can be most of them**: such vectors are long, and the system of the rows left out
is then dense.

The kept set does not change. A row that has left the span of the kept rows gets the
best combination there is, its bound is then its distance from the span, and the
certificate refuses on the gap: the rank has changed, and that is the caller's signal to
search again. No certificate and no step stand after a refresh.

| RefreshStatus | Meaning |
| --- | --- |
| Refreshed | Every bound stands for the values given; `refreshed-count` vectors were replaced, possibly none |
| RefreshNothingKept | No kept set stands |
| RefreshBadTarget | The target was not positive and finite, or the fraction not from 0 to 1; nothing changed |
| RefreshNoShift | A vector needed the second route, and no certificate has been issued for this kept set to take its shift from; the bounds stand, and the vectors the first route served were replaced |
| RefreshFactorRefused | A vector needed the second route, and the kept rows' solve could not be factored; the bounds stand, and the vectors the first route served were replaced |
| RefreshNonFinite | A bound was not finite; the vectors stand, the bounds do not |

A vector whose solve refuses, or whose multipliers are not finite, is left as it was.
Another pattern or shape raises the typed error of the Gram refill.

**The system of the rows left out is not factored again by a refresh.** Where no vector
changed the rows it holds, the factor the search made stays, now a factor of earlier
vectors, and the step's solves of that system become conjugate gradients on the current
system preconditioned by it. With a factor of the vectors as they stand a solve is one
solve, exactly as before a refresh; with an earlier one it is a few iterations. A solve
makes at most `system-iteration-limit` iterations, 8, to a residual of 2⁻⁴⁰ of its
right-hand side; **if it has not converged by then the factor is made again**, from the
vectors as they stand, and the solve is one solve. Where a vector did change its rows
the system is analyzed and factored again by the refresh. `system-iterations` is the
iterations of the standing step and `system-factor-count` the factorizations since the
search. Nothing here is proved or needs to be: what the step returns for this system is
[a proved residual](#how-it-is-computed).

**What it costs.** Where no vector changes its rows and the arrays have reached their
size, a refresh allocates nothing. [The record](batch18d.md#what-it-costs) has the times
along a consumer's drags: under a tenth of a millisecond where one vector keeps its rows,
and 0.16 to 0.30 ms for each vector that needs the second route.

## The step

`step!(kept, r, out)` writes the step a solver takes against a residual r on the matrix
the standing certificate is of:

    out = −B⁺ r ,        B = C J_S ,

where `J_S` is the kept rows and C puts each row left out back as its kept combination of
them: the identity on the kept rows over W, whose row for a row left out is that row's
multipliers. `B⁺ r` is the minimum-norm least squares solution for B. The sign is that
of [the damped step](damped.md), so that a loop may take either.

**B is not the matrix, and the step says how far that takes it.** B has rank r exactly
and differs from the matrix by at most the gap bound. Two bounds stand after a step:

- `step-solve-bound`: on `‖out + B⁺ r‖₂`, what the two solves below left.
- `step-bound`: on the distance from out to **the step of the matrix truncated at the kept
  count**, `−J_r⁺ r` with `J_r` the nearest matrix of rank r, which is what a singular
  value decomposition cut at the certified floor gives. It includes the first.
- `step-norm-term` and `step-outside-term`: what the second adds to the first, in its two
  terms, so that a caller can see which one a wide bound comes from.

With γ the gap bound, τ the bound the certificate proved, x the exact step for B and
`x*` the truncated one, the second adds to the first

    ‖x − x*‖₂  ≤  2γ ( 2 X / τ + ρ / τ² ) / (1 − 4γ/τ) ,

where X bounds `‖x‖₂` and ρ bounds the residual outside the range, `‖r − J x‖₂ + γ X`;
[the design](design-rank-deficient-step.md#how-far-x-is-from-x) derives it. The term in
X is `step-norm-term` and the term in ρ `step-outside-term`; the second carries the
square of 1/τ and is the larger wherever the residual is far from the range. Each quantity
in it is proved: the norms from computed vectors with their rounding allowed for, γ and τ
by the certificate. **Where `4γ/τ` is not below 1 the inequality has no solution and the
call is `StepBoundNotReached`.** A certificate holds whenever γ is below τ, so vectors
that are far out of date can be certified and still be refused here.

| StepStatus | Meaning |
| --- | --- |
| StepTaken | The step was written and its two bounds stand |
| StepNoCertificate | No certificate stands for the values held |
| StepBadSize | The residual has not the row count, or the output not the column count |
| StepNonFiniteRhs | An entry of the residual was not finite |
| StepFactorRefused | The system of the rows left out, or the kept rows' solve, could not be factored |
| StepBoundNotReached | `4γ/τ ≥ 1`; or a residual of a solve was not finite; or the refinement was still at least halving its bound at its last correction |
| StepNonFiniteResult | The step cannot be held in the caller's units |

Every status but the first writes nothing to `out` and leaves no step standing. A step
needs a standing certificate and does not survive a search, a reproof or another
certificate. The residual and the output may be one array where the matrix is square.

### How it is computed

1. `y = C⁺ r`, a least squares problem whose normal matrix `I + Wᵀ W` is solved through
   the smaller `I + W Wᵀ`, of the order of the rows left out and as sparse as W. It is
   analyzed and factored once, by the search; its eigenvalues are at least 1. It is
   corrected while its proved residual at least halves, and the candidate with the
   smallest is kept.
2. `x = J_S⁺ y` by the kept rows' augmented matrix, on the analysis the certificate uses,
   with its own two shifts, refined against the same matrix with no bottom block.

Nothing is solved on the kept rows' Gram matrix, which a row touching every column makes
dense. [The design](design-rank-deficient-step.md#the-solve-on-the-kept-rows) has the
shifts, the contraction of a correction by blocks, and the measurement they were chosen
on.

**What the solve's bound rests on is the residuals and nothing else.** For any pair
`(x̂, λ̂)`, with `ρ₁ = τ x̂ + J_Sᵀ λ̂` and `ρ₂ = J_S x̂ − y`,

    x̂ − J_S⁺ y  =  J_S⁺ ρ₂  +  (I − J_S⁺ J_S) ρ₁ / τ        exactly,

the two parts orthogonal, so the distance is at most `(‖ρ₁‖² + ‖ρ₂‖²)^½ / τ`. Each residual
is formed with a bound on its own rounding, and **τ is the one the standing certificate
proved, for these kept rows and these values**. The error of y is bounded the same way,
by `‖Cᵀ(r − C ŷ)‖`, and reaches x divided by τ. No rate of convergence enters.

**The bound on y widens with the multipliers.** That residual's own rounding grows with
the square of the largest multiplier, and the bound with it: over the step's norm it was
9 × 10⁻⁹ at a multiplier of 2¹⁰, 4 × 10⁻⁵ at 2¹⁶ and 10⁻² at 2²⁰. A certificate limits
the multipliers to about the matrix's size over the bound it proved; where the rows are
of comparable size they are of order 1.

**The refinement stops by a rule on that bound.** It is formed after the first solve and
after each correction, and the iterate with the smallest is kept.

- A correction that leaves the bound above half the smallest so far ends the refinement.
  The kept iterate is the step, **with its bound, however wide**: the caller compares it
  with its own gate.
- At most `correction-limit` corrections are made, 8. If the last still at least halved
  the bound the refinement had not finished and the call is `StepBoundNotReached`. A
  correction contracts by 2⁻¹⁰ or better while it is above rounding, so six exhaust the
  precision; two are left for corrections that halve the bound by chance at the bottom.

`step-corrections` is how many the standing step made. The factorization for the solve is
made by the first step after a certificate and serves every later step until the next
certificate, reproof or search.

**Three diagnostics, on which nothing is refused**: `solve-growth` and
`solve-backward-bound` of that factor, and `predicted-rate`, the contraction that
backward bound predicts, `(2⁻¹⁰ + η) / (1 − η)` with `η = 2β/τ`. The prediction
overstates by decades, since the backward bound is a bound through a trace; the bounds
of a step do not use it.

**The residual is divided by the power of two at its largest entry**, as the matrix is by
its own, and the step is multiplied back by the two. A step at any scale the format holds
is the step at ordinary scale, bit for bit, and the bounds are those times the scale.

## Dependent sets

Which rows depend on which is read from the kept vectors: the vector of a row left out
holds a multiplier at each kept row that row is made from. **A multiplier computed in
floating point is never exactly zero or not**, so the question is put with a cut, and
answered only where it can be proved.

`class-supports!(kept, cut)` works on the matrix the standing certificate is of. For the
vector of each row left out:

- δ bounds the 2-norm of the difference between the stored multipliers ŵ and the exact
  ones w, the multipliers of the row against the kept rows for the matrix as stored. Since
  the kept rows are independent, `ŵ − w = (ŵ J_S − a) J_S⁺`, the vector's residual times
  the pseudo-inverse of the kept rows, so **δ is the vector's proved residual bound over
  the bound the certificate proved on the kept rows**. `multiplier-error-bound` returns it.
- The cut is relative to the vector's own largest entry, its own 1 included, which lies
  between `top_lo = max(1, top − δ)` and `top_hi = top + δ` for the exact vector. So the
  units of the rows cancel within a vector.
- A stored multiplier is **in** when `|ŵ| − δ > cut · top_hi` and **out** when
  `|ŵ| + δ < cut · top_lo`. A kept row the vector does not hold, having never been reached
  or having been dropped by the search, is out when `δ < cut · top_lo`.
- **Anything else is a refusal**, `SupportsUndecided`, with the vector and the entry
  readable: a multiplier wrongly classed would join two sets or part one.

Then every row is given a label. Rows are joined through every multiplier classed in;
**each set that holds a row left out is a dependent set**, a row left out with nothing
classed in being a set of one, as a row of zeros is; and a kept row that no vector holds
is in no dependency. The sets do not depend on which rows the search happened to keep:
they are the connected components of the rows as a matroid.

- `dependent-sets(kept)`: a borrowed array as long as the row count. Rows with one label
  depend on one another; the labels run from 0 to `dependent-set-count` less one, and
  **a row in no dependency has the row count as its label**.
- `support-flags(kept)`: beside each entry of `multiplier-values`, 1 where it is classed
  in and 0 where out.
- `undecided-vector`, `undecided-entry`: after a refusal, the position of the vector and
  which of its entries; an entry equal to the vector's `multiplier-count` means a kept row
  it does not hold could not be classed out, which is δ being too large for the cut.

| SupportStatus | Meaning |
| --- | --- |
| SupportsClassed | Every multiplier is classed and every row has its label |
| SupportsNoCertificate | No certificate stands for the values held; nothing was classed |
| SupportsBadCut | The cut was not above 0 and below 1; nothing was classed |
| SupportsUndecided | A multiplier could not be put on either side of the cut; nothing stands |

Classed supports stand until the next search, reproof, refresh or certificate. The call
allocates only when the kept vectors have grown since it last ran.

"Out" means below the cut for the matrix as stored. A multiplier that is exactly zero in
a model and 10⁻¹⁷ in the stored matrix is classed out only because 10⁻¹⁷ is below the cut;
no computation on rounded data proves more.

## The right null space

**The package applied to the transpose gives a basis of the right null space**, with
nothing added. Let `K = analyze(Jᵀ)` and search it. Its rows are J's columns, so:

- the rows it leaves out are the columns of J that depend on the others, n less the rank
  of them;
- the kept vector of each is a vector z over J's columns with 1 at that column, its
  stored values at kept columns and zero elsewhere, and its `residual-bound` is a proved
  bound on `‖J z‖₂`;
- the vectors are independent, each having 1 where every other has 0, so they are a
  basis of what J sends to within those bounds of zero;
- `certify!(K, lower-floor, upper-floor)` is the same certificate for the singular
  values, which J and its transpose share: with it, exactly as many singular values as
  kept columns are above the upper floor, and the basis has the dimension of the null
  space at that floor.

The basis is sparse where the multipliers are and **not orthonormal**. An orthonormal one,
for a caller that projects onto it, is designed in
[the design](design-rank-deficient-step.md#a-basis-of-the-right-null-space) and not built.
`reprove!` and `refresh!` keep the basis for new values as they keep any kept vectors.
Two handles, one on J and one on its transpose, share nothing and may be used at
different times.

## Bounds on the largest singular value

`largest-singular-bounds!` writes two numbers with `out[0] ≤ σ₁ ≤ out[1]`.

- The upper is `(‖J‖₁ ‖J‖_∞)^½`, rounded upward.
- The lower is `‖J v‖₂ / ‖v‖₂` for the vector the power method reaches in `steps` steps
  from the vector of ones. That quotient is a lower bound whatever v is; what is returned
  is lowered by a bound on the rounding of that one product and of the two lengths, so
  that it holds for the exact quantities. It is 0 for a matrix of zeros.

Thirty steps gave a pair within a factor of 2.2 on every matrix tested. A caller that
needs a floor as a fraction of the largest singular value multiplies the fraction by
each: the lower floor by `out[0]` and the upper by `out[1]`.

## Using it

```
kept = analyze(J)                                   once for the pattern
find-kept-rows!(kept, J, tolerance, 1e-12)          once for a kept set
each later matrix of the pattern:
    refresh!(kept, J, lower-floor / 2)              or reprove!(kept, J) where no row changes
    certify!(kept, lower-floor, upper-floor)        Certified: the kept set serves
    gap above floor after a refresh: the rank has changed; search again
    step!(kept, r, p)                               StepTaken: p is the step
    step-bound(kept)                                compare with the solver's own gate
    class-supports!(kept, 1e-6)                     when the dependent sets are wanted
```

[The first example](../examples/kept-rows.ct) finds a redundant row, proves its bound
again for other values, reads what the certificate proved, refreshes the vector, and
labels the dependent sets.
[The second](../examples/rank-deficient-step.ct) takes the step on two equal rows and
reads its bound.

## Representation choices

The certificate is proved on the kept rows' augmented matrix and not on their Gram
matrix. Both proofs were measured, and prove the same bound; the Gram matrix of a
pattern with a row that touches every column is dense, and a factorization of it at
every call cost a thousand times the other.

The vectors are kept one after another in two arrays that grow by doubling during the
search. The pattern is held a second time by rows, which is what lets a bound read only
the rows a vector touches.

The system of the rows left out is factored by the search and not by the step, because
W changes only when the vectors do.

The handle owns its own search. A caller that already holds a
[dependent rows](dependent-rows.md) analysis of the same pattern has two; the pattern
analysis is the smaller cost of either.
