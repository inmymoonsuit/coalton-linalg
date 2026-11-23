# Kept rows and the certificate

State: the first of two parts, implemented in Batch 18a. The package is
`coalton-linalg/sparse/rank-deficient`. [The design](design-rank-deficient-step.md) has
the mathematics, the measurements that shaped it and what the second part, the step
itself, will add. [The batch record](batch18a.md) has the evidence.

## What it is for

A matrix with redundant rows is rank deficient, and a solver that takes many steps on
one such matrix, its values changing a little each time, needs three things at every
step, cheaply and **holding as stated or refused by name**:

- the rows the others depend on, and for each other row the combination that gives it;
- a bound on how well each combination holds for the values of this step;
- a proof that no singular value lies in a band about the floor the solver cuts at.

[Dependent rows](dependent-rows.md) finds the first once. This package keeps what was
found, proves the second again for new values without finding anything, and issues the
third. It solves nothing yet; the step built on it is the second part.

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
- `certify!(kept, lower-floor, upper-floor) -> CertificateStatus`: that the gap bound is
  below the lower floor and every singular value of the kept rows above the upper, both
  floors in the matrix's own units.
- `proved-lower-bound(kept)`: what the standing certificate proved of the kept rows, at
  least the upper floor asked for and usually far above it.
- `residual-bound(kept, position)`, `gap-bound(kept)`: the bounds as last proved.
- `kept-count`, `left-out-count`, `left-out-rows`; `multiplier-start`,
  `multiplier-count`, `multiplier-rows`, `multiplier-values` for the vectors as kept;
  `kept-found?`, `certified?`, `row-count`, `column-count`.
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
    reprove!(kept, J)
    certify!(kept, lower-floor, upper-floor)        Certified: the kept set serves
    gap above floor: the vectors are out of date for these values; search again
```

[The example](../examples/kept-rows.ct) finds a redundant row, proves its bound again for
other values, and reads what the certificate proved.

## Representation choices

The certificate is proved on the kept rows' augmented matrix and not on their Gram
matrix. Both proofs were measured, and prove the same bound; the Gram matrix of a
pattern with a row that touches every column is dense, and a factorization of it at
every call cost a thousand times the other.

The vectors are kept one after another in two arrays that grow by doubling during the
search. The pattern is held a second time by rows, which is what lets a bound read only
the rows a vector touches.

The handle owns its own search. A caller that already holds a
[dependent rows](dependent-rows.md) analysis of the same pattern has two; the pattern
analysis is the smaller cost of either.
