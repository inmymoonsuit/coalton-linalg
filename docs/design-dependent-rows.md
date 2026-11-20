# Design: dependent rows of a sparse matrix, with certificates

State: **a design for review; nothing here is built.** It answers the roadmap item
"dependent rows with a certificate" and is written to be read by the consumer that asked
for it before any code is written. Where it says "proved" it means a statement that
follows from quantities the library computes together with a bound on their rounding;
where it says "found" it means a choice that is then checked, and that could have been
made differently.

## The question

A matrix A with m rows and n columns, sparse, has rows that are constraints. When some
constraints are redundant or in conflict, the rows are linearly dependent, and a caller
wants three things:

1. **A set S of rows to keep** that is independent, with every other row dependent on it.
2. **For each row left out, the kept rows it depends on**, with the multipliers. That set
   is what a user is shown: constraint k is redundant with constraints i, j and l.
3. **To know what of this is proved** and at what tolerance.

"Dependent" has no exact meaning for floating-point data: any matrix is within rounding of
a full-rank one. The design therefore fixes, first, which statements can be proved.

## What can be proved

Write `A_S` for the kept rows, `a_k` for a row left out, σ_i for singular values in
descending order, r for the number of kept rows.

**P1, independence of the kept rows.** `σ_min(A_S) > τ`. This is the
[rank certificate](quasidefinite.md#a-rank-certificate-from-the-inertia) the library
already has, applied to `A_S`: one signed factorization of `[[t I, A_Sᵀ], [A_S, s I]]`
with `s t = τ²`, the pivot counts, and `factor-backward-bound`. It is a proof.

**P2, dependence of a row left out.** `‖a_k − yᵀ A_S‖₂ ≤ ε_k` for multipliers y that are
exhibited. This needs no factorization to verify: the residual is formed from the stored
entries, and its rounding is bounded by `γ (|a_k| + |y|ᵀ |A_S|)`, so an upper bound on its
norm is a proof that a_k lies within that distance of the row space of `A_S`. It is a
proof whatever method produced y, which is what lets the method that *finds* y be a
heuristic.

**P3, the rank with a gap.** From P1 and P2 together, with no further computation:

- `σ_r(A) ≥ σ_min(A_S) > τ`, because adding rows to a matrix cannot lower its singular
  values, counted from the largest.
- `σ_{r+1}(A) ≤ √(Σ_k ε_k²)`, because replacing each row left out by its projection gives
  a matrix of rank r within that Frobenius distance of A.

So where `√(Σ ε_k²) < τ`, **A has numerical rank exactly r with a proved gap**: r singular
values above τ and the rest below the bound. Where the gap does not open, the report
says so, and the two bounds are still each true.

**P4, the multipliers are the unique ones for S.** `A_S` has independent rows, so the
combination of them closest to a_k is unique: the dependency set of k *relative to S* is
determined. A different S gives a different, equally valid set. The design reports the
one for the S it found, and does not claim it is the smallest over all choices of S.

**P5, a one-sided count without the rows.** The inertia of `[[t I, Aᵀ], [A, s I]]` for the
whole matrix counts the singular values of A on each side of `√(s t)`, and one
factorization proves one side of that count. With the margins added, as the certificate
adds them, c negative pivots prove that **at least c singular values exceed `√(s t)`**;
with the margins subtracted, p positive pivots beyond n prove that at least p are at most
`√(s t)`. A bracket takes two factorizations at two thresholds. This is weaker than it
looks, for the reason the certificate's contract gives: that matrix is not quasi-definite,
and with singular values on both sides of the threshold its trailing complement is
indefinite, so a small pivot and a backward bound too large for any margin are ordinary
outcomes, each meaning *not proved*. P5 is offered as a cheap first look, and P3 is the
statement the design relies on.

**What is not proved:** that S is the best conditioned choice, that a reported dependency
set is minimal over all bases, or anything at a tolerance below what the arithmetic
resolves, which differs between the two methods below.

## Finding S and the multipliers: two methods

Both use the elimination order and the structure the library already computes for
`A Aᵀ`, and both are followed by the same proofs.

### Method G: skip small pivots in the factorization of A Aᵀ

Factor `A Aᵀ` in the analyzed order. In exact arithmetic the pivot of row k is the
squared distance from a_k to the span of the rows kept before it. When a pivot falls at or
below a threshold θ², leave row k out: give its column of the factor no entries and carry
on. What results is the factor of `A_S A_Sᵀ` for the rows kept, in the structure already
analyzed. Row k of the factor, already computed when its pivot is examined, holds the
components of a_k along the kept rows before it, so the multipliers are one triangular
solve with the kept part of the factor. They can be refined with
[`refine!`](refinement.md) against `A_S A_Sᵀ` applied as two products.

- **Cost.** One numeric factorization of the kind the library does now, plus a triangular
  solve for each row left out. The kernel change is one test and one branch, and keeping
  the row of the factor that is now discarded.
- **Resolution.** A squared distance is computed as a difference of sums of squares, so it
  is lost below about u times the row's squared norm: **distances below √u, about
  1.5 × 10⁻⁸ of the row's size, cannot be told from zero.** A row whose true distance lies
  between the caller's tolerance and that floor is left out by the factorization although
  it is independent at the caller's scale. P2 then reports its true distance, which the
  residual resolves to rounding, so the error is visible, but the method cannot keep such
  a row: `A_S A_Sᵀ` would be singular to working precision.
- **Fits** a caller whose tolerance for dependence is no finer than about 10⁻⁷ relative.

### Method Q: a triangular factor of Aᵀ by rotations, without Q

Compute the upper triangular R with `Aᵀ = Q R` by Givens rotations, processing the rows of
Aᵀ one at a time and never storing Q. R has the structure of the Cholesky factor of
`A Aᵀ`, which the analysis already gives, so it is stored in the factor's existing arrays,
one row of R where a column of the factor is now. A row of Aᵀ is scattered into a work
vector; its leading entry is rotated against the row of R that owns that position, which
zeroes it and fills the work vector only within that row's structure; the next position
is the smaller of that row's first later index and the incoming row's next index, which
is a walk up the elimination tree; a row that reaches a position no row of R owns is
installed there.

A diagonal entry of R means something only when every row has been processed: an entry
installed early is one row's contribution, and later rotations enlarge it. With all rows
in, the diagonal entry at position k is the distance from row k of A to the span of the
rows of A before it in the order, accumulated as a root of a sum of squares and not as a
difference of squares.

Rank is then detected by **one sweep over the positions in order**. Where the diagonal at
position k is at most θ in magnitude, row k of A is left out: the diagonal is set to
zero, and row k of R is taken out and processed again as an incoming row, which carries
its remaining entries into the rows after it along the same path in the tree, so nothing
is filled outside the analyzed structure. The rows before k are untouched, and by the
time position k is examined every earlier decision has been applied, so its diagonal is
final: the distance from row k of A to the span of the rows *kept* before it. At the end
**the empty positions are the rows of A left out**, and column k of R, above the diagonal,
holds the components of a left-out row along the kept ones. Its multipliers are one
triangular solve with the kept part of R.

- **Cost, not yet measured.** Aᵀ has one row for each column of A, and each walks a path in
  the tree, at four multiplications for each entry of each row of R it meets. The total
  depends on the order the rows are taken in as well as on the structure, and can exceed
  the factorization of `A Aᵀ` by a factor that grows with the number of columns over the
  number of rows. **The design does not claim a figure.** A count of rotations from the
  pattern alone, on the consumer's own matrices, is the first thing to be done and decides
  whether this method is affordable as a routine step or only as a diagnosis after a
  failure. A new numeric kernel of moderate size; no new analysis.
- **Resolution.** Rounding, about u times the size of A: **distances down to 10⁻¹⁵ of the
  row's size are resolved.** Each diagonal set to zero is a perturbation of A of at most θ
  in norm, which the proofs do not need bounded in advance, since P1 and P2 are checked
  afterwards on A itself.
- **Fits** a caller whose tolerance is finer than Method G resolves, and gives R, which is
  also a factor for least squares on A without forming `A Aᵀ`.

### What the two share, and the recommendation

Neither method pivots for rank: both take the rows in the fill-reducing order and keep
the first independent ones. That is what keeps them sparse, and it is why the kept set
can be badly conditioned without either method noticing: every kept distance exceeds θ
and the smallest singular value of `A_S` can still be far below θ. **P1 is what catches
that**, and is the reason the certificate is a separate step and not a by-product. If P1
fails at the caller's τ, the report is "S found, independence not proved at τ", with the
largest τ that is proved; a caller may then raise θ and repeat.

**Recommendation.** P1, P2 and P3 are built first whichever method finds the set, since
they are what turns a found set into a result and they cost little: P2 is new and small,
and P1 and P3 are existing calls and arithmetic. For the finding, Method Q if the
consumer's tolerance is finer than about 10⁻⁷ and its rotation count is affordable, and
Method G otherwise: G is less work and reuses more, and at a fine tolerance it cannot do
the job, so building it first would then be building something to replace. The consumer's
tolerance and the rotation count on its matrices decide this; the first is a question
below and the second is a measurement to make before any kernel is written.

## Proposed interface

Names are provisional. Following the library's conventions: a handle made once for a
pattern, repeated numeric calls that allocate nothing, outcomes as values.

- `analyze-rows(matrix) -> RowAnalysis`: the plan for `A Aᵀ`, its analysis, the work
  vectors and the storage for R. Once for a pattern.
- `factor-rows!(analysis, values, threshold) -> RowStatus`: Method Q on these values.
  Afterwards:
  - `kept-count`, `kept-rows` (a borrowed index array, the rows of S), `dropped-rows`;
  - `row-distance(k)`: the diagonal of R for a kept row, its distance to the rows kept
    before it;
  - `multipliers!(analysis, k, out) -> SolveStatus`: for a row left out, its multipliers
    on the kept rows, written into a vector indexed by row of A, zero elsewhere.
- `dependence-residual!(matrix, values, k, multipliers, work) -> F64`: P2, an upper bound
  on `‖a_k − yᵀ A‖₂` with the rounding of its own computation included. It takes any
  multipliers, so a caller can zero the small ones and ask again for the reduced set.
- P1 and P5 are the existing certificate. **P1 needs no second pattern.** With the values
  of the rows left out written as zeros, `[[t I, Aᵀ], [A, s I]]` is, after a permutation,
  the same matrix for `A_S` beside a block `s I` for the rows left out, so the plan and
  the analysis already made for A serve: the count to require is r negative pivots and
  `n + m − r` positive, and the argument from the margins is unchanged. This rests on the
  factorization accepting stored entries that are zero, which the tests of the batch must
  establish before it is relied on.

## What the tests would establish

- **Exact truth.** Integer matrices with planted dependencies: the rank, a basis and the
  multipliers are computed in rational arithmetic, and the kept set must be a basis and
  the multipliers must agree.
- **The bounds, exactly.** P2's returned bound against the exact distance in rationals,
  on every row left out. P3's two bounds against exact singular values squared, through
  the exact inertia of `A Aᵀ − λ I` at the bounds.
- **Resolution, as a roster.** Rows at distances 10⁻⁴ down to 10⁻¹⁵ of their size from the
  span of others, with a threshold on each side: kept above, left out below, the distance
  reported to rounding.
- **The case that fools the method.** A matrix whose rows are each far from the span of
  those before them and whose smallest singular value is tiny, in the manner of the Kahan
  matrix: every row is kept, and P1 must refuse to certify it.
- **Structure.** Rows of zeros, repeated rows, a row dependent on rows that come after it
  in the order, a dense column, more rows than columns.
- Planted source controls, allocation, both policies, as for every batch.

## Questions for the consumer

1. **What is the dependence tolerance, relative to a row's size?** Finer than about 10⁻⁷
   means Method Q; coarser allows Method G.
2. **Is the answer a proof or a diagnosis?** If the kept set must be independent at a
   stated τ for a later proof to stand on, P1 is required at every call. If the set is
   shown to a user as "these constraints conflict", P2 alone may be enough, and P1 is run
   only when asked.
3. **Does a dependency set have to be small?** The set relative to S is unique but depends
   on the order. If a user is better served by the smallest set, that is a different
   problem (a shortest dependent combination), not solved by either method, and it should
   be said now.
4. **How large, and how many left out?** Each row left out costs a triangular solve and a
   residual. For Method Q the shape matters more than the size: how many columns for each
   row.
5. **Is a row of A a constraint, always?** The methods factor Aᵀ; a caller with constraints
   as columns passes the transpose.

## Provenance, in advance

Building a triangular factor row by row with rotations, on the structure of the normal
equations' factor, is the method of George and Heath (1980), and Heath (1982) extended it
to rank-deficient problems; skipping small pivots in a factorization of the normal
equations is common practice in sparse least squares. All are known here from
descriptions. The sweep set out above is derived in this document from what the diagonal
of R measures; how closely it coincides with the 1982 treatment has not been checked
against the paper, and the comparison is to be made and recorded in
[provenance.md](provenance.md) before anything is built. The statements P1 to P5 and the
choice to prove them afterwards, whatever found the set, are this library's. No
implementation will be consulted.
