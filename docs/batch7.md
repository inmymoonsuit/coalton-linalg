# Batch 7: signed pivots and the augmented matrix

State: completed; verified at safety 1 and safety 0 on 2026-11-20, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract and the derivation are in
[quasidefinite.md](quasidefinite.md).

## What was delivered, and why

A consumer's minimum-norm step on a wide matrix A went through A Aᵀ and the positive
definite factorization. Where A has one dense column, A Aᵀ is dense: on a public hub
system of 354 rows the product's factor holds 40,908 entries and takes 3.08 ms, against a
matrix of 1,200 entries. This batch gives a route whose cost follows A:

- `coalton-linalg/sparse/augmented`: a plan for `K = [[top I, Aᵀ], [A, −bottom I]]` on a
  fixed pattern, analyzed once and refilled without allocating.
- `factor-signed!` in `coalton-linalg/sparse/spd-ldlt`: the library's L D Lᵀ with pivots of
  either sign, a refusal only at an exactly zero pivot, and no numerical pivoting.
- `factor-growth`: the exact largest entry of `|L| |D| |Lᵀ|` over the largest magnitude
  factored, measured by every signed factorization.
- `positive-pivot-count` and `negative-pivot-count`.
- `analyze-quasidefinite`: minimum-degree analysis that holds the trailing block back, each
  eliminated leading index releasing one held neighbour.
- An artifact driver, `scripts/check.ps1 -Task Augmented`, and a compiled example.

**The route is not accurate on every matrix, and the batch does not claim it is.** A
factorization without numerical pivoting is accurate only when its factors stay small.
What the batch delivers is the factorization, an order that keeps the factors small where
a dense column is the difficulty, and an exact figure that says, for each factorization,
whether they stayed small. The caller chooses between this route and the product from that
figure.

## How the design was reached

Three findings shaped it, each a measurement made with NumPy on the public systems before
the library code was written, and repeated by the library afterwards.

1. **Plain minimum degree is unusable for a small trailing block.** A row of A with few
   entries is eliminated before any of its columns, its pivot is exactly −bottom, and
   growth is two over bottom relative to the row's size. The table below has it on every
   system.
2. **Holding every row until one of its columns has gone is not enough.** Two rows that
   share their only eliminated column leave the second with a pivot of −bottom again:
   growth 1.5 × 10¹⁵ at a bottom of 2⁻⁵⁰ on a small seeded case. Releasing one row a
   column removes that.
3. **No rule that reads only the pattern is enough.** On the grid system growth is between
   1.3 and 1.9 over the relative bottom under the held order and 3.5 under the plain one.
   Its rows
   are independent but nearly dependent on the columns eliminated before them. This is why
   the growth figure is part of the interface and not a diagnostic.

An earlier draft reported the largest multiplier instead of growth. It was dropped: on the
hub system the largest multiplier is 2.5 × 10⁴ while growth is 244 and the solve is
accurate to rounding, so the multiplier would have refused a good factor.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`, with every earlier roster unchanged. For this batch:

- **Reconstruction, entry by entry.** 234 augmented matrices and 48 symmetric matrices
  with a strictly dominant diagonal of mixed signs: named cases, seeded sparse sources at
  several shapes and densities with trailing blocks from 1 down to 2⁻⁵⁰, a hub, a chain,
  shared-column cases, stored zeros and shifts. For every stored entry the exact
  difference between the permuted matrix and L D Lᵀ, in rational arithmetic, is within
  `(p + 2) u |L| |D| |Lᵀ|`, p being the products that entry received. The largest was
  0.424 of that allowance.
- **Growth is exact.** In every case the reported figure equals the exact one, computed in
  rational arithmetic from the exported factor, to within 2⁻⁴⁰ relative. The largest growth
  in the roster is 9.8 × 10¹⁵, and the bounds hold there as they do at 1.
- **Solves, row by row.** For each case one solve with a seeded right-hand side: each row
  of the exact residual is within `(3w + 5) u |L| |D| |Lᵀ| |x̂|`. The largest was 0.145 of
  that allowance.
- **Inertia.** The pivot counts equal the exact inertia from a separate rational
  elimination, and for a full-row-rank source they are n and m.
- **The held order.** In every held case no row is eliminated before a column it has an
  entry in. On a hub of twelve rows the dense column is in the last two positions and the
  factor has no fill at all. With `leading` 0 or n the permutation equals plain minimum
  degree. Over the seeded roster the plain order's growth at a trailing block of 2⁻³⁰
  reaches 1.5 × 10¹⁰ where the held order's stays at or below 25.8 at every block size.
- **One kernel.** On a positive definite matrix `factor-signed!` and `factor!` give the
  same factor bit for bit; negated, the first succeeds with every pivot negative and the
  second refuses. A shift applied by the call equals the matrix stored with the shift
  subtracted, bit for bit, and growth is measured against the shifted diagonal.
- **Dependent rows.** A 3-by-2 source of rank 2 factors with counts (2, 3) at trailing
  blocks of 2⁻¹⁰, 2⁻⁴⁰ and 2⁻⁵⁰, and is refused as `SignedZeroPivot` at 2⁻⁵² and 2⁻⁶⁰,
  where the block is lost to rounding. The contract says what a caller may conclude.
- **Lifetime.** A refusal withdraws the factor, the counts, the growth and the stored
  estimate; `factor!` after a signed factor reports no growth; a refactorization on one
  analysis reports the growth a fresh analysis reports; the inverse-norm estimate works
  from a signed factor and lies within the estimator's bound of the exact value.
- **The plan.** Exact refills, a source column and a source row with no entry, stored
  zeros kept, both scalars placed and signed, every refusal by its constructor with the
  output untouched, the pattern comparison of the checked refill on dimensions, pointers
  and row indices separately, and the entry budget.
- **Interface.** Three external callers are rejected by the compiler: constructing a plan
  directly, reaching its source pattern, and using a `FactorStatus` where a
  `SignedFactorStatus` is required. A public caller is accepted.
- **Allocation.** A compiled Coalton caller making refill, signed factorization and solve:
  zero bytes over 10,000 calls after 10,000 warm-up calls, at both policies.
- **The driver.** Its regression test checks every written file of a three-row system
  against exact values in both orders, a zero pivot named with its column, a typed
  refusal file, a skipped system, and the removal of stale artifacts.

**Planted source controls.** Each of these edits was applied alone to the source named and
the focused tests run. Every planted source compiled, every one failed the tests, and the
restored source passed. Three were missed by the first version of the tests, which then
gained the cases that catch them: growth against the unshifted diagonal, growth carried
over from an earlier factorization, and a checked refill that compared pointers only.

| Source | Planted defect | Caught |
| --- | --- | --- |
| spd-ldlt | The signed entry refuses every pivot that is not positive | yes |
| spd-ldlt | The positive definite entry accepts a negative pivot | yes |
| spd-ldlt | A negative pivot is counted as positive | yes |
| spd-ldlt | A refusal keeps the pivot counts | yes |
| spd-ldlt | Growth leaves out the pivots themselves | yes |
| spd-ldlt | Growth adds the multipliers unsquared | yes |
| spd-ldlt | Growth is not divided by the size of the matrix | yes |
| spd-ldlt | Growth measures the diagonal before the shift | yes |
| spd-ldlt | A new signed factorization starts from the previous bounds | yes |
| spd-ldlt | The positive definite entry reports a growth | yes |
| spd-ldlt | A leading dimension above the size is accepted | yes |
| ordering | Held vertices are never released | yes |
| ordering | An eliminated column releases every row it held | yes |
| ordering | Every vertex is held from the start | yes |
| ordering | A held vertex's key is raised by one only | yes |
| augmented | The trailing diagonal keeps the sign it was given | yes |
| augmented | The leading diagonal is always one | yes |
| augmented | The source rows are not moved below the leading block | yes |
| augmented | The budget ignores the source entries | yes |
| augmented | Non-finite diagonal scalars are accepted | yes |
| augmented | The output may be the input array | yes |
| augmented | The checked refill does not compare the row indices | yes |

## Measurements

Six public wide systems from a geometric constraint solver, at safety 1, one process,
times the best of 20 after warm-up. Their entries are near 10⁻³ and their smallest
singular values near 5 × 10⁻⁶; top is 1 and bottom is 10⁻¹⁶ of the largest squared row
norm, about 2 × 10⁻²², which puts the regularized solution within 2 × 10⁻¹¹ of the
minimum-norm one on five of them. Entries count the factor with its diagonal. The last
three columns are the positive definite route on A Aᵀ, analyzed and factored by `factor!`
in the same process.

| System | Rows × columns | Order | Entries | Growth | Analysis ms | Factor ms | Solve ms | A Aᵀ entries | A Aᵀ analysis ms | A Aᵀ factor ms |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Hub, 100 spokes | 354 × 402 | held | 2,362 | 244 | 0.65 | 0.015 | 0.006 | 40,908 | 58.2 | 3.08 |
| | | plain | 2,355 | 2.0 × 10¹⁶ | 0.39 | 0.015 | 0.006 | | | |
| Chain, 25 | 352 × 400 | held | 1,958 | 13.0 | 0.56 | 0.012 | 0.005 | 1,202 | 0.30 | 0.006 |
| | | plain | 1,950 | 2.0 × 10¹⁶ | 0.54 | 0.012 | 0.005 | | | |
| Chain, 100 | 1402 × 1600 | held | 7,808 | 13.0 | 2.36 | 0.053 | 0.022 | 4,802 | 1.26 | 0.027 |
| | | plain | 7,800 | 2.0 × 10¹⁶ | 2.20 | 0.056 | 0.024 | | | |
| Frame, medium | 602 × 950 | held | 3,754 | 6.1 | 1.12 | 0.022 | 0.010 | 1,708 | 0.47 | 0.009 |
| | | plain | 3,654 | 1.0 × 10¹⁶ | 1.09 | 0.022 | 0.011 | | | |
| Frame, large | 2052 × 2400 | held | 11,804 | 15.4 | 3.32 | 0.089 | 0.033 | 6,859 | 1.27 | 0.040 |
| | | plain | 11,254 | 2.0 × 10¹⁶ | 2.13 | 0.084 | 0.034 | | | |
| Grid, 15 × 15 | 1277 × 1290 | held | 21,361 | 1.9 × 10¹⁶ | 10.5 | 0.358 | 0.041 | 45,614 | 26.8 | 1.22 |
| | | plain | 13,354 | 3.5 × 10¹⁶ | 3.83 | 0.152 | 0.027 | | | |

Refill takes between 0.002 and 0.011 ms and the measured calls allocate nothing.

**What the table says.**

- **On the hub the augmented route is the point of the batch**: a factor of 2,362 entries
  against 40,908, factored in 0.015 ms against 3.08 ms and analyzed in 0.65 ms against
  58 ms, with growth 244 at every bottom tried.
- **Elsewhere the product is cheaper.** Where A Aᵀ stays sparse its factor is about half
  the size of the augmented one and about twice as fast. The augmented route is for the
  dense-column case, not a replacement.
- **On the grid the route does not work at a small bottom, and says so.** Growth is 10¹⁶
  under both orders; the product is the route there, at 1.22 ms.

**Accuracy against a singular value reference**, computed afterwards with NumPy from the
driver's written solutions, as the largest error in x relative to its largest entry:

| System | Order | Relative bottom | Growth | Error against the regularized solution | Error against the minimum-norm solution |
| --- | --- | ---: | ---: | ---: | ---: |
| Hub | held | 10⁻¹⁶ | 244 | 3.8 × 10⁻¹⁴ | 5.0 × 10⁻¹² |
| Chain, 100 | held | 10⁻¹⁶ | 13.0 | 6.0 × 10⁻¹⁴ | 6.9 × 10⁻¹² |
| Chain, 100 | plain | 10⁻¹² | 2.0 × 10¹² | 3.9 × 10⁻⁷ | 3.9 × 10⁻⁷ |
| Chain, 100 | plain | 10⁻¹⁶ | 2.0 × 10¹⁶ | 1.3 × 10⁻¹ | 1.3 × 10⁻¹ |
| Frame, large | held | 10⁻¹⁶ | 15.4 | 2.0 × 10⁻¹⁴ | 1.6 × 10⁻¹¹ |
| Frame, large | plain | 10⁻¹⁶ | 2.0 × 10¹⁶ | 2.2 × 10⁻² | 2.2 × 10⁻² |
| Grid | held | 10⁻⁴ | 1.3 × 10⁴ | 1.9 × 10⁻¹⁴ | 1.0 |
| Grid | held | 10⁻¹² | 1.7 × 10¹² | 7.1 × 10⁻⁷ | 7.3 × 10⁻⁶ |
| Grid | held | 10⁻¹⁶ | 1.9 × 10¹⁶ | 4.4 × 10⁻³ | 4.4 × 10⁻³ |

Where growth is small the solve is accurate to rounding and the remaining distance to the
minimum-norm solution is the regularization. Where growth is large the error follows it.
The first grid row shows the other limit: at a bottom large enough for the factorization
to be accurate, the regularized solution is not near the minimum-norm one, because the
grid's smallest singular value is 5 × 10⁻⁷.

One exception is worth knowing: on the hub the plain order's x is accurate to 3.8 × 10⁻¹⁴
at growth 2 × 10¹⁶, while the residual of the whole augmented system is 8 × 10⁻⁷. Growth
bounds the error; it does not predict it. A caller that gates on growth will sometimes
refuse a usable factor, never accept a bad one on the strength of the figure.

## Limits

- No numerical pivoting. A symmetric indefinite matrix that is not quasi-definite can be
  refused at a zero pivot or factor with unbounded growth; the dense
  [indefinite factorization](indefinite-ldlt.md) pivots, this one does not.
- The solution is that of the regularized system. Correcting it to the minimum-norm
  solution is the caller's, for example by refinement with repeated solves.
- The held order's guarantee is structural: no row before a column of its own, no two
  rows released by one column. It is not a bound on growth.
- The measured growths are on six systems of one kind. They say the hub case works and
  the grid case does not; they do not classify patterns in general.
- The ordering is the exact minimum-degree analysis already in the library, whose cost
  grows faster than the factor on large patterns: 10.5 ms on the grid against a 0.36 ms
  factorization. A faster analysis is the next item on the roadmap, to be decided by a
  consumer's measurement of where its time goes.
- One host and finite rosters establish neither universal behaviour nor portability beyond
  the supported runtime.
