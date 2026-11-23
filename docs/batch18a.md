# Batch 18a: kept rows and the certificate

State: completed; verified at safety 1 and safety 0 on 2026-11-23, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract is
[rank-deficient.md](rank-deficient.md) and the design
[design-rank-deficient-step.md](design-rank-deficient-step.md).

## What this is

The first of the two parts of a step on a rank-deficient matrix, the second item of a
consumer's list. That consumer's loop takes many steps on one pattern whose rows are
redundant, the values moving a little each time, and at each it must know that the
rank it is using is the rank: that no singular value lies in a band about the floor it
cuts at. The design was reviewed three times before any code, and measured on the
consumer's own patterns between the first review and the second. This batch builds what
the step stands on; the step itself is 18b.

## What was built

`coalton-linalg/sparse/rank-deficient`:

- `find-kept-rows!` runs [the search for dependent rows](dependent-rows.md) once and keeps
  the vector of every row left out. An entry below a drop tolerance times the vector's
  largest is set to zero, and the bound is then proved for the vector as kept.
- `reprove!` proves every kept vector's bound again for other values of the pattern. It
  reads only the rows the vectors touch, searches nothing and allocates nothing.
- `certify!` proves that the gap bound is below a lower floor and every singular value of
  the kept rows above an upper floor, or refuses by one of five names.
- `largest-singular-bounds!` bounds the largest singular value of any stored matrix from
  both sides, for a caller whose floors are fractions of it.

Three things in it came from the measurement and not from the first design.

**The vectors are kept sparse.** A computed vector has an entry at nearly every kept row
the factor reaches, and on the consumer's patterns all but three or four of them are
below 10⁻¹³. Kept dense they would be millions of numbers; kept above a relative
tolerance they are a few a row, and the bound, proved for what is kept, moved by less
than 2%.

**The search's threshold and the certificate's floors are separate arguments.** A search
run at the floor itself works at its own resolution and can keep a dependent row; run
well inside the gap it does not.

**The proof is on the kept rows' augmented matrix.** Both proofs the design named were
measured and prove the same bound. The one on the Gram matrix factors a matrix that a
row touching every column makes dense.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`.

- **Bounds against exact arithmetic.** Forty-two seeded matrices of small integers on
  seven shapes, from 3 by 4 to 13 by 10, each with rows that are combinations of others:
  the kept count is the exact rank, and each of 142 vectors' bounds is at least the exact
  norm of the exact product of the stored matrix with the vector as kept, at four drop
  tolerances from 10⁻¹² to 0.99. Each bound is also no looser than
  [the bound over the whole matrix](dependent-rows.md#interface), and the gap bound is at
  least the root of the sum of the exact squares.
- **Certificates against exact arithmetic.** On the same matrices, 546 certificates
  asked at floors from 10⁻⁹ to 10³: 399 issued, and for each the Gram matrix of the kept
  rows less the square of the bound proved is positive definite exactly; 147 refused,
  none at a floor a fifth of a bound proved before.
- **The proof's margin.** Twenty pairs of rows that are a rotation times diag(1, 1/16),
  so that their smallest singular value is 1/16 to rounding, which is the first bound the
  ladder tries. For eight of them it is exactly at or below 1/16. Every bound proved is
  true exactly; none is 1/16.
- **Reproof.** Fifteen matrices with every value moved by a part in a million: the
  vectors are unchanged bit for bit, the new bounds hold exactly for the moved values, the
  gap bound grows a thousandfold and more, and the certificate that held at 10⁻⁹ is
  refused on the gap. Back at the first values the first bounds return bit for bit.
- **Across the range.** One matrix at six scales from 2^-900 to 2^900 gives the same rows
  and the same vectors bit for bit, and bounds and a proved lower bound that are those at
  scale 1 times the scale. [The scale-sweep roster](rosters.md#what-is-declared) has two
  entries for the package: the kept rows are exact from 2^-800 up with a bound looser at
  2^-1000 and no input at 2^-1020, and the bounds on the largest singular value are exact
  everywhere.
- **The largest singular value.** Thirty-six seeded matrices on six shapes: the square of
  the upper bound less `Jᵀ J` has no negative eigenvalue, and `Jᵀ J` less the square of
  the lower has one that is not negative, both exactly. The widest pair is a factor 2.18
  apart after thirty steps.
- **The drop.** Two rows that are 1/1024 of another and that plus a third: cut at a half,
  the first vector keeps its one entry of 1/1024 and the second loses its own, and the
  second bound is the length of what was lost, the root of 5, to twelve digits.
- **Refusals and edges.** Each value of each status that an input reaches is produced,
  with the handle's state checked after it: before any search; five drop tolerances and
  three thresholds that cannot be used, each after a standing set; the search's
  underflow; six pairs of floors; a floor above the kept rows; a floor below what the
  proof resolves, against rows 4 × 10⁻⁸ apart; a floor off the ladder, proved as itself;
  a matrix of another pattern or shape; a matrix of zeros; a row with no entry and a row
  of stored zeros; one entry; and a reproof whose bound cannot be held, a kept vector
  with an entry of 2^510 against values of ordinary size, after which the vectors stand,
  the bounds do not, and a reproof against the first values restores them.
- **Allocation.** A compiled caller that proves the bounds again and certifies: zero
  bytes over 10,000 calls after warm-up.
- **Aliasing.** `largest-singular-bounds!` joins [the aliasing roster](rosters.md): each
  of its three pairs of arrays is refused.
- **Planted source defects.** Fifty-eight edits to the new source, each applied
  alone, with the package's suite and both rosters run. Every planted source compiled.
  Fifty-seven failed a test on the first run and one was missed: with every
  factorization counted as accurate at the floor, the one test that tells
  `CertificateUnresolved` from `CertificateKeptNotAboveFloor` asked at a single floor,
  where the planted source happened to answer with the right name. That test now asks at
  four floors, three of which the planted source answers wrongly; run again it is caught,
  and the restored source passes. Five of the tests above were written while the edits
  were being chosen, because an edit had nothing to fail: the relative drop, the floor
  off the ladder, the proof's margin, the state after a refused search, and the reproof
  that cannot hold a bound.

| Planted defect | Number | Caught |
| --- | --- | --- |
| The pattern by rows reads an entry from the wrong place; a new handle reads as having proved a bound; a count or a vector is read where none stands | 4 | yes |
| A vector's bound: one below the normal range is not raised; the row left out is missing from its own product; the rounding of a sum is not allowed for; magnitudes are summed with their signs; the workspace is not cleared between vectors, two ways; the values are not divided | 7 | yes |
| The gap bound is of half the squares; a bound that is not finite is taken as proved | 2 | yes |
| The search: any drop tolerance, or a tolerance of one, is taken; a refused search leaves the earlier set, or its bounds; a search leaves the earlier certificate, or the bound it proved; each of the search's own refusals is reported under another name | 8 | yes |
| The vectors: the rows left out are counted among the kept; the drop is relative to the whole vector with its own 1, or is an absolute size, or is not made; the vectors already kept are lost when the arrays grow | 5 | yes |
| The reproof: made with no kept set; leaves the certificate, or the bound it proved; does not check the pattern; proves the bounds for the values already held; leaves the earlier bounds after one that cannot be held | 6 | yes |
| A bound as read: a residual bound or the gap bound is not multiplied back; a residual bound is read past the last | 3 | yes |
| The largest singular value: a work vector too long, or the output shared with one, is taken; the upper bound is the 1-norm alone; the power method never leaves the vector of ones; the lower bound is not divided by the vector's length; either bound is not multiplied back | 7 | yes |
| The proof: the pivots' signs are not counted; the bottom block has the other sign; a bound is proved with no margin for the factorization's error; no factorization counts as accurate at the floor | 4 | yes |
| The proof: every factorization counts as accurate at the floor | 1 | on the second run, after the test was strengthened |
| The certificate: asked of bounds that do not stand; a refusal leaves the earlier one, or the bound it proved; a lower floor of zero, or one above the upper, is taken; the gap bound is not compared with the floor, or is compared in the divided units; the kept rows' values are read from the wrong places; the upper floor is not divided; the floor itself is never tried; the bound proved is not multiplied back | 11 | yes |

## What it costs

Timed once on the consumer's seven patterns, read in place and not part of the repository,
on a machine otherwise idle, safety 1, one process. The search is at 10⁻⁶ of the upper
bound on the largest singular value, the drop at 10⁻¹², the floors at 2 × 10⁻⁸ of the
lower bound and 10⁻⁶ of the upper. Each reproof time is the mean of twenty calls.

| Pattern | Rows by columns | Left out | Entries a vector | Search with vectors | Reproof | Certificate, first | Certificate, later |
| --- | --- | --- | --- | --- | --- | --- | --- |
| wheel, well constrained | 1,604 by 1,602 | 2 | 204 | 0.24 s | 0.012 ms | 1 ms | under 1 ms |
| wheel, redundant | 2,004 by 1,602 | 402 | 4.0 | 0.86 s | 0.041 ms | 1 ms | under 1 ms |
| wheel, redundant and loose | 1,803 by 1,602 | 399 | 3.0 | 0.77 s | 0.036 ms | 1 ms | under 1 ms |
| lattice 15, redundant | 4,414 by 3,600 | 814 | 4.4 | 0.57 s | 0.088 ms | 2 ms | 1 ms |
| lattice 15, redundant and loose | 4,384 by 3,600 | 812 | 4.0 | 0.57 s | 0.081 ms | 2 ms | 1 ms |
| lattice 21, redundant | 8,698 by 7,056 | 1,642 | 4.2 | 2.50 s | 0.18 ms | 6 ms | 2 ms |
| lattice 21, redundant and loose | 8,656 by 7,056 | 1,640 | 4.0 | 1.99 s | 0.17 ms | 5 ms | 2 ms |

**A reproof costs about a tenth of a microsecond for each row left out** where a vector
has four entries, and six where it has two hundred. The design's measurement proved
each bound over the whole matrix and took 0.004 to 0.09 s for the same work; finding the
vectors afresh took 0.4 to 2.5 s. So a step that keeps its vectors pays for them once.

Every certificate was issued, at both pairs of floors. The gap bounds are 10⁻¹⁵ to
5 × 10⁻¹⁴ of the scale, and 3 × 10⁻¹³ to 3 × 10⁻¹¹ after every value was moved by a part
in 10¹²: the bound follows the values, as it should. The bounds proved on the kept rows
are 3.3 × 10⁻⁶ of the scale on the wheels that keep every column, and 4.2 × 10⁻⁵ to
1.7 × 10⁻⁴ elsewhere, each within the ladder's factor of four of what the design's
finer search measured. The first certificate after a search walks the ladder down from
an eighth; each later one starts from the remembered bound.

The analysis took 0.01 to 0.08 s and the two bounds on the largest singular value, thirty
steps, 0.4 to 1.5 ms; the pair was a factor 1.13 to 1.94 apart.

## Limits

- **Nothing is solved.** The step, its bound to the truncated answer, the kept rows alone
  as an option and the refresh of one vector from current values are the next batch.
- **The bound proved on the kept rows is a lower bound from a ladder**, a factor of four
  apart, and can be that far below the smallest singular value, or further at the top.
- **The proof resolves about a ten-millionth of the matrix's size.** Below that it
  answers `CertificateUnresolved`. The squared proof of the design, which resolves less,
  is not built.
- **Out-of-date vectors are refused, not repaired.** The caller searches again; a refresh
  of single vectors is the next batch.
- **`FindVectorRefused`, and a bound that is not finite at the search itself, are not
  produced by any test**: no input tried reaches them.
- **The handle holds a second analysis of the pattern**, beside one a caller may already
  hold for [dependent rows](dependent-rows.md).
- **The exact tests are small**, up to 13 by 10. The consumer's patterns, up to 8,698
  rows, are timed and not checked against exact arithmetic.
