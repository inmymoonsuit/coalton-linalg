# Roadmap

## Completed: batches 1 and 2

The standalone library now provides sparse SPD LDLᵀ, stable vector norms,
validated rectangular CSC/COO construction, explicit duplicate handling,
transpose, row/column permutations, and both matrix-vector product directions.
The SPD analysis consumes the checked sparse structure through analyze-csc.

Validation includes the original 20 matrices, 336 seeded exact-input cases,
252 permutation/scaling relations, 96 repeated-refactor cases, floating-point
boundary observations, independent primitive tests, and planted checker defects.
The toolchain is content-pinned and was rebuilt from an empty cache.
See [batch record](docs/batches-1-2.md) for measured scope and limitations.

## Initial review snapshot

The [Coalton design review](docs/design.md#coalton-interface-review-2026-11-19)
led to private Symbolic/NumericFactor handles, distinct named factor/solve
statuses, and explicit factor!/solve! names. External compile probes cover
construction, destructuring, internal field access, and status-type confusion.
That snapshot's cached results were removed by the revision below. Private
handles still do not enforce read-only borrowing of inspection arrays.

The repository now supplies dependency registration, a source acquisition
manifest, setup instructions, and the maintainer-selected MIT license.
Fresh public-source acquisition and complete serial checks at both safety
policies passed in an isolated source export. [The review record](docs/initial-review.md)
states the measured evidence and limitations.

## Pre-publication revision

[The focused revision](docs/revision-2.md) simplifies factor statuses and transparent
handles, makes latest-diagnostic lifetime explicit, reuses Coalton's standard
numeric operations, supplies typed errors, and bounds the current analysis bitmap.
The complete serial checks passed at both safety policies with unchanged numerical
artifacts. The exact-singular regression demonstrates the computed-pivot contract.

## Batches 3 through 8

Batches 3 through 8 are completed, Batch 3 in four steps 3a through 3d, with one addition,
[6a](docs/batch6a.md): a Cholesky-form numeric entry for positivity proofs. Batches 7 and
8 and what follows them are set by a consuming application's measurements; what follows
is listed under *Next*. [The batch record](docs/batch3a.md)
reports independent review, both safety gates and measurements through 65,536 rows. [Detailed contracts and acceptance
gates](docs/planned-batches.md) incorporate application feedback while preserving
this library's generic mathematical scope. Each step receives independent review
and its own measurements. The compact analysis replaces the old 4096-row cap
with an explicit factor-entry budget and overflow-safe representability limits.

### Batch 3: scalable sparse analysis and repeated operations

1. **Completed (3a):** compact graph analysis preserves deterministic exact
   minimum-degree ordering, uses column counts for allocation, and enforces explicit
   total-factor-entry budgets. Independent checks cover tree/postorder identities;
   retaining those data awaits a traversal consumer. The old 4096 cap is removed
   after measured cases through 65,536. Numeric artifacts are unchanged.
2. **Completed (3b):** matvec! and transpose-matvec! overwrite caller storage,
   share the existing arithmetic, and check aliases before mutation. Independent
   exact products, source controls and six allocation measurements per policy pass.
   [The batch record](docs/batch3b.md) specifies refusal and allocation scope.
3. **Completed (3c):** fixed-pattern plans for A^T A and A A^T write the lower
   triangle with every diagonal, retain structural zeros and refill caller storage
   without allocating. Exact rational references, planted source controls and
   eight allocation measurements per policy pass. [The batch record](docs/batch3c.md)
   holds the measurements; [the contract](docs/gram.md) the calling rules.
4. **Completed (3d):** a lower estimate of the 1-norm of the inverse of the factored,
   shifted matrix in at most 11 solves, the matching matrix norm, and their reciprocal
   product. It is neither a singular-value ratio nor a rank certificate. Exact rational
   references, planted source controls and an allocation measurement per policy pass.
   [The batch record](docs/batch3d.md) holds the measurements; [the contract](docs/condition.md)
   the meaning.

### Batch 4: dense symmetric indefinite LDL^T and inertia

**Completed.** `coalton-linalg/dense/indefinite-ldlt` factors a dense symmetric matrix
with 1-by-1 and 2-by-2 pivots under bounded multipliers, solves, and reports the inertia
of the computed block diagonal with both eigenvalues of each 2-by-2 block. An exactly zero
pivot completes with a zero count and its solve is refused by name. 365 exact-reference
factorizations, shifted eigenvalue counts, augmented systems, thirteen planted source
controls and an allocation measurement per policy pass. [The batch record](docs/batch4.md)
holds the measurements; [the contract](docs/indefinite-ldlt.md) the derivation and what a
count near a threshold does not say. Dense delivery does not establish large sparse scaling.

### Batch 5: orthogonal decompositions

**Completed.** `coalton-linalg/dense/orthogonal` gives column-pivoted Householder QR of
any shape, a rank the caller chooses from the diagonal of R against an explicit strict
threshold, the complete orthogonal decomposition for that rank, minimum-norm least
squares with its residual norm, an orthonormal nullspace basis, and Q applied without
being formed. 215 exact-reference decompositions, the Kahan matrix, seventeen planted
source controls and an allocation measurement per policy pass. [The batch
record](docs/batch5.md) holds the measurements; [the contract](docs/orthogonal.md) the
derivation and what the diagonal of R does not say.

### Batch 6: dense symmetric eigenvalues and SVD

**Completed.** `coalton-linalg/dense/jacobi` computes all eigenvalues and eigenvectors of a
dense symmetric matrix by cyclic Jacobi rotations, and all singular values and vectors of
a dense matrix of any shape by one-sided Jacobi rotations on its columns, with a complete
right factor whose trailing columns are a nullspace basis, a rank count against the
caller's threshold, and a sweep limit whose exhaustion is a named outcome. Every computed
value in 280 test matrices is bracketed by an exact inertia count; twenty planted source
controls and two allocation measurements per policy pass. [The batch
record](docs/batch6.md) holds the measurements, including how slow these routines are
against a tridiagonalizing reference; [the contract](docs/jacobi.md) the derivation and
what is promised against what is only observed.

### Batch 7: signed pivots and the augmented matrix

**Completed.** For a wide matrix with a dense column, where A Aᵀ is itself dense:
`coalton-linalg/sparse/augmented` lays out `[[top I, Aᵀ], [A, −bottom I]]` on a fixed
pattern, `factor-signed!` factors it with pivots of either sign and no numerical pivoting,
`analyze-quasidefinite` orders it so that no row precedes a column of its own, and
`factor-growth` reports exactly how much accuracy each factorization kept, because no
order fixed in advance is accurate on every matrix. On a public hub system the factor has
2,362 entries against 40,908 for the product and takes 0.015 ms against 3.08 ms; on a grid
system growth is large and the product remains the route. 282 exact-reference
factorizations with entrywise bounds, twenty-two planted source controls and an allocation
measurement per policy pass. [The batch record](docs/batch7.md) holds the measurements;
[the contract](docs/quasidefinite.md) the bounds and what growth does and does not say.

### Batch 8: the same elimination order, computed faster

**Completed.** Analysis returns the permutation and structure it returned before, for
every input, from a representation that does not store the elimination graph: eliminated
vertices are kept as neighbour sets that absorb one another, a hub keeps a row of flags,
and columns are sorted by one pass through the rows. On thirteen public patterns
`analyze` takes between 1.7 and 25 times less, 5.3 times less on the loose grid a
consumer measured as its cost. An independent set-based elimination agrees on 1,188
analyses and agreed with the earlier code first; 392 driver output files are unchanged
byte for byte; twenty-three planted source controls pass. [The batch
record](docs/batch8.md) holds the measurements and what was tried and removed.

### Next

Each decided by a consumer's measurement and not before it:

1. The fixed-pattern Gram analysis, which a consumer measured as large as the ordering
   on long chains.
2. A further gain in ordering time, which now needs an approximate degree or merged
   vertices and so a decision about changing the order.
3. Dense fronts in the numeric factor, only if factorization comes to dominate.

### Later, when measured need justifies it

Lanczos, LSMR, conjugate gradients for SPD systems, sparse QR, sparse indefinite
factorization, supernodal methods and numeric types beyond F64. Application
Jacobians and dual-number residual evaluation remain application code.

## Release work

Configure CI workers with the pinned dependencies;
and expand supported hosts only after measured validation. The current supported
and tested host is Windows/x86-64 SBCL. The initial numeric type remains F64.

Applications own geometric Jacobians, nonlinear iteration, rank thresholds,
units, and integration. No SolveSpace promotion decision is made here.