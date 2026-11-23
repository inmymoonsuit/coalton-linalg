# Batch 19: dependent sets, and the right null space

State: completed; verified at safety 1 and safety 0 on 2026-11-23, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract is
[rank-deficient.md](rank-deficient.md#dependent-sets) and the design
[design-rank-deficient-step.md](design-rank-deficient-step.md).

## What this is

The last item of the design a consumer reviewed: certified supports of the left null
vectors, and a basis of the right null space. They serve readers that say which of a
problem's conditions conflict and in which directions it is still free; they are not on
the path of a step. The consumer was asked what those readers take before anything was
written, and answered from its own code:

- **Labels are enough.** Its reader merges supports into sets and reads nothing else, so
  what is returned is a label for each row, not lists and not multipliers.
- **One cut, relative to each vector.** Its reader takes a row as part of a dependency
  when its weight is above a millionth of that null vector's largest, the vector's own 1
  included. So the cut is that fraction, and the proved bound is a band on either side.
- **The orthonormal basis can wait.** Its reader of freedoms wants one, and is not on the
  new route yet. The sparse basis is built and documented; the other stays designed.

## What was built

`class-supports!(kept, cut)` in `coalton-linalg/sparse/rank-deficient`, with
`dependent-sets`, `dependent-set-count`, `support-flags`, `multiplier-error-bound`,
`supports-classed?`, `undecided-vector` and `undecided-entry`.

- Each stored multiplier is classed in or out against the cut times its vector's largest
  entry, by a proved bound δ on its distance from the exact one; a kept row a vector does
  not hold is classed out where δ is under the cut; anything undecided is a refusal that
  names the vector and the entry.
- Rows are joined through the multipliers classed in. Every set holding a row left out
  gets a label; a kept row no vector holds gets the row count, which says it is in no
  dependency.

**The bound changed while building, for the better.** The design had
`δ = ‖g‖ / τ²`, the vector's residual carried through the kept rows twice. On the
consumer's redundant lattice of 21 that left a kept row undecided at a cut of a millionth,
δ being 4 × 10⁻⁶. The kept rows being independent, the stored multipliers are themselves
times `J_S J_S⁺`, so their distance from the exact ones is the vector's residual times
`J_S⁺`: **δ is the vector's proved residual bound over τ, once.** On that lattice it is
2 × 10⁻¹⁰, and the pattern is classed. The exact tests below hold the new bound against
rational arithmetic.

**The right null space needed no code.** The package applied to a matrix's transpose
leaves out the columns that depend on the others and keeps, for each, a vector with 1 at
that column, with a proved bound on the matrix's product with it: a basis, with its
certificate. It is now a section of the contract and has tests of its own.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`.

- **The sets, against exact arithmetic.** Forty-five seeded matrices of small integers on
  nine shapes. The exact multipliers of each row left out are found in rational
  arithmetic, rows are joined wherever one is not zero, and the sets that result are the
  sets the library labels, on every matrix: 76 sets in all, 29 matrices with more than
  one, and 91 rows in no dependency, each carrying the row count.
- **The bound, against exact arithmetic.** On the same matrices, for every vector, the
  exact 2-norm of its stored multipliers less the exact ones is at most
  `multiplier-error-bound`.
- **At the cut.** A row that is one kept row plus a multiple of another: with a multiple of
  a thousandth, in, and all three rows one set; of a billionth, out, the other kept row
  in no dependency; of exactly a millionth, the cut, refused, naming that entry. A row
  whose multipliers are a thousandth and a hundred-millionth: the second is out beside
  the row's own 1, where beside the first alone it would be in.
- **Where the bound is too large for the cut.** Two kept rows a part in a million apart:
  classed at a cut of a millionth, and refused at 10⁻¹² and 10⁻¹⁵ with the entry named as
  the vector's count, a kept row it does not hold.
- **Refusals and state.** No certificate, before a search and after one; five cuts that
  cannot be used; a classing withdrawn by a certificate, a reproof, a refresh and a
  search, and by a refusal. A row of zeros and a row with no entry are each a set of one.
- **Across the range.** The same matrix at five scales from 2^-900 to 2^900 gives the same
  labels and the same bounds; [the scale-sweep roster](rosters.md#what-is-declared) has an
  entry for the call.
- **The right null space.** Twenty-five matrices on five shapes whose columns depend on
  one another: the handle on the transpose keeps as many columns as the rank and leaves
  out the rest; each of its 61 vectors has 1 at its own column and 0 at every other
  column left out; and each one's bound is at least the exact norm of the exact product
  of the matrix with it.
- **Allocation.** A compiled caller that proves the bounds again, certifies and classes:
  zero bytes over 10,000 calls after warm-up.
- **Planted source defects.** Fifteen edits to the new source, each applied alone,
  with the package's suite and both rosters run. Every planted source compiled and every
  one failed a test, on the first run; the restored source passes. One test was added
  while the edits were being chosen, because an edit had nothing to fail: the row whose
  multipliers are a thousandth and a hundred-millionth, which tells a cut against the
  vector's largest entry, its own 1 counted, from a cut against its largest multiplier.

| Planted defect | Number | Caught |
| --- | --- | --- |
| Refusals and state: supports are classed with no certificate; any cut is taken; a cut of one is taken; a refusal leaves the earlier classing standing; what withdraws a step leaves the classing standing; the count of sets is read when no classing stands; the vector that could not be classed is not named | 7 | yes |
| The bound: the multipliers are taken to be exact; the bound is the vector's residual bound, not that over the bound proved on the kept rows | 2 | yes |
| The classing: the cut is against the largest multiplier, the row's own 1 not counted; a multiplier below the cut is never classed out; a kept row the vector does not hold is classed out whatever the bound; every stored multiplier is flagged in | 4 | yes |
| The sets: rows are joined through multipliers classed out as well; only the representative of a set is labelled | 2 | yes |

## On the consumer's patterns

Seven patterns a consumer supplied, each run once in a process of its own at safety 1 on
a machine otherwise idle. The search's threshold is a millionth of the matrix's largest
singular value as bounded from above, the certificate's floors are that and 2 × 10⁻⁸ of
the lower bound, and the cut is a millionth. Every pattern is certified and classed.

**The rows**, which is what a reader of conflicts asks for:

| Pattern | Rows by columns | Left out | Dependent sets | Rows in no dependency | Largest δ | Classing |
| --- | --- | --- | --- | --- | --- | --- |
| wheel, well constrained | 1,604 by 1,602 | 2 | 1 | 1,198 | 1.5 × 10⁻⁹ | 0.007 ms |
| wheel, redundant | 2,004 by 1,602 | 402 | 1 | 999 | 2.3 × 10⁻¹⁰ | 0.027 ms |
| wheel, redundant and loose | 1,803 by 1,602 | 399 | 1 | 1,003 | 1.1 × 10⁻¹² | 0.023 ms |
| lattice 15, redundant | 4,414 by 3,600 | 814 | 2 | 1,858 | 1.1 × 10⁻¹⁰ | 0.047 ms |
| lattice 15, redundant and loose | 4,384 by 3,600 | 812 | 28 | 1,920 | 8.8 × 10⁻¹² | 0.044 ms |
| lattice 21, redundant | 8,698 by 7,056 | 1,642 | 2 | 3,610 | 1.6 × 10⁻¹⁰ | 0.109 ms |
| lattice 21, redundant and loose | 8,656 by 7,056 | 1,640 | 40 | 3,696 | 8.8 × 10⁻¹² | 0.101 ms |

The classing is a pass over the stored multipliers and costs a tenth of a millisecond at
most, against a search of 0.2 to 3 seconds. The largest δ is 1.5 × 10⁻⁹, under a
six-hundredth of the cut. By the design's first bound δ reached 4.5 × 10⁻⁴ on the first
wheel, 6.8 × 10⁻⁵ on the second, and 2.6 × 10⁻⁶ and 3.8 × 10⁻⁶ on the two redundant
lattices, above the cut on four of the seven; the lattice of 21 was run that way while
building, and refused.

**The columns**, the same package on each pattern's transpose, which is a basis of the
right null space where there is one:

| Pattern | Search | Columns left out | Gap bound, of the scale | Sets among the columns | Columns in none | Largest δ |
| --- | --- | --- | --- | --- | --- | --- |
| wheel, well constrained | 1 ms | 0 | | | | |
| wheel, redundant | 1 ms | 0 | | | | |
| wheel, redundant and loose | 24 ms | 198 | 6.5 × 10⁻¹⁵ | 198 | 24 | 5.5 × 10⁻¹¹ |
| lattice 15, redundant | 5 ms | 0 | | | | |
| lattice 15, redundant and loose | 14 ms | 28 | 3.1 × 10⁻¹⁴ | 2 | 180 | 3.5 × 10⁻¹¹ |
| lattice 21, redundant | 11 ms | 0 | | | | |
| lattice 21, redundant and loose | 36 ms | 40 | 4.4 × 10⁻¹⁴ | 2 | 252 | 4.3 × 10⁻¹¹ |

The four patterns that keep every column have no right null space and the search says
so in milliseconds. On the three loose ones the columns left out are as many as the rows'
search found the rank short by (1,602 less 1,404, 3,600 less 3,572, 7,056 less 7,016),
the two searches sharing nothing. The gap bound is the bound on the matrix's product
with the whole basis, and is at rounding. The search of the columns is from some thirty
to some hundreds of times quicker than the search of the rows; why was not looked into.
The columns classed too, though nothing was asked of that: on the loose wheel each
free direction stands alone, 198 sets of columns; on the loose lattices the 28 and the
40 fall into two sets.

## Limits

- **Corrected 2026-11-23: the sets are those of the kept rows at hand.** This record and
  the contract said they did not depend on the rows kept. That holds for the exact test
  above, which joins rows wherever a multiplier is not zero, and not for a cut: the
  "at the cut" case above is one circuit of three rows labelled as a set of two and a row
  in no dependency. An independent review found it;
  [the contract](rank-deficient.md#dependent-sets) now says what is true and the tests
  hold the same rows in two orders. The same review found that a vector holding every
  kept row was refused where its bound was above the cut, for a kept row that does not
  exist; that is fixed, and the "bound too large for the cut" case above now classes.
- **A multiplier within δ of the cut is a refusal**, and so is a δ above the cut where a
  kept row is not held. Nothing is classed on a guess.
- **δ carries one over the bound proved on the kept rows**, which the ladder leaves up to
  four times below their smallest singular value.
- **The orthonormal basis of the right null space is not built.** The sparse one is not
  orthonormal, and a caller that projects onto it must allow for that.
- **The handle on the transpose is a second search**, of the columns; nothing is shared
  with the handle on the matrix.
- **`multiplier-error-bound` is for the stored vector**, entries dropped by the search
  counted as zeros.
- **The exact tests are small**, up to 13 by 10. The consumer's patterns are classed and
  their figures reported, not checked against exact arithmetic.
