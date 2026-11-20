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

## Batches 3 through 10

Batches 3 through 10 are completed, Batch 3 in four steps 3a through 3d, with one addition,
[6a](docs/batch6a.md): a Cholesky-form numeric entry for positivity proofs. Batches 7 and
8 were set by a consuming application's measurements; what follows is listed under
*Next*. [The batch record](docs/batch3a.md)
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

**Completed (7a):** `factor-backward-bound` gives the signed factorization's backward error
as one number, with what underflow can add, and the contract sets out how the pivot
counts and that number certify that every singular value of a sparse matrix exceeds a
threshold. Over 343 cases checked against exact inertia no certificate was issued for a
false statement. [The record](docs/batch7a.md).

### Batch 8: the same elimination order, computed faster

**Completed.** Analysis returns the permutation and structure it returned before, for
every input, from a representation that does not store the elimination graph: eliminated
vertices are kept as neighbour sets that absorb one another, a hub keeps a row of flags,
and columns are sorted by one pass through the rows. On thirteen public patterns
`analyze` takes between 1.7 and 25 times less, 5.3 times less on the loose grid a
consumer measured as its cost. An independent set-based elimination agrees on 1,188
analyses and agreed with the earlier code first; 392 driver output files are unchanged
byte for byte; twenty-three planted source controls pass. The fixed-pattern Gram analysis
was then changed the same way, without sorts and with an unchanged pattern, and takes
between 1.7 and 13 times less on the same systems. [The batch record](docs/batch8.md)
holds the measurements and what was tried and removed.

### Batch 9: scaling by powers of two

**Completed.** `coalton-linalg/sparse/scaling` equilibrates the rows and columns of a
sparse matrix so that each has its largest magnitude in [1/2, 1), scales a symmetric
positive definite matrix by its diagonal into [1/2, 2), and carries vectors through the
factors. Every factor is a power of two, so each scaled value is the exact product. 246
matrices checked exactly, twenty planted source controls and an allocation measurement
per policy pass. [The batch record](docs/batch9.md) holds the controls; [the
contract](docs/scaling.md) says which solutions each kind of scaling leaves unchanged.

### Batch 10: iterative refinement, and a solve on the Cholesky-form factor

**Completed.** `coalton-linalg/refinement` refines any approximate solve against an
operator the caller applies: the library forms the residual, stops on an absolute
tolerance in the infinity norm, on a caller's stall factor or at a step limit, returns
the best iterate seen and keeps the history. `solve-cholesky!` solves with the proof
factor, so a proof factorization and a refinement can replace a second factorization.
Sixty refinements reproduced bit for bit by an independent recomputation, both uses
checked against exact rational solutions, twenty-two planted source controls and two
allocation measurements per policy pass. [The batch record](docs/batch10.md) holds the
controls; [the contract](docs/refinement.md) says what refinement converges to.

### Next

From 2026-11-20 the library develops ahead of any one consumer's measurements, toward
what solving in CAD and engineering work is likely to need. The order is a consumer's
ranking of what it would use first:

1. Dependent rows with a certificate: which rows of a sparse matrix are dependent, and
   what can be proved about it. A written design comes before any code.
2. Sparse LU with partial pivoting for a square nonsymmetric matrix.
3. Iterative least squares and minimum norm on an operator, with damping, and conjugate
   gradients for the definite case.
4. A few eigenvalues of a sparse symmetric matrix: Lanczos, shift-invert through the
   signed factorization, and the inertia as a count that none was missed.
5. Continuous integration on the pinned toolchain.

**Measured since, and settled:** a consumer timed its whole first step on its loose grids
after Batch 8. Analysis fell from 141 to 17 ms and from 272 to 47 ms and is now a sixth
of that step, so no further ordering work and no reuse of an analysis across patterns is
warranted. The numeric factorization is now about 97% of its steady step. Dense fronts
in the numeric factor are therefore justified by measurement; the expected gain without
vector instructions or a tuned dense kernel is a factor of about one and a half to two,
and where they sit in this list is a priority decision not yet taken.

### Later

Sparse QR, sparse indefinite factorization with pivoting, supernodal methods and numeric
types beyond F64. Application
Jacobians and dual-number residual evaluation remain application code.

## Release work

Configure CI workers with the pinned dependencies;
and expand supported hosts only after measured validation. The current supported
and tested host is Windows/x86-64 SBCL. The initial numeric type remains F64.

Applications own geometric Jacobians, nonlinear iteration, rank thresholds,
units, and integration. No SolveSpace promotion decision is made here.