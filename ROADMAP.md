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

## Batches 3 through 6

Batch 3 is completed in four steps, 3a through 3d; Batches 4 through 6 are planned. [The batch record](docs/batch3a.md)
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

Add stable 1-by-1 and 2-by-2 pivoting, solves and block-correct inertia. Validate
shifted eigenvalue counts and augmented systems with exact/reference cases.
Numerical threshold reliability, variable nullity and singular-solve behavior
have distinct contracts. Dense delivery does not establish large sparse scaling.

### Batch 5: orthogonal decompositions

Build Householder QR and full-rank least squares, then column pivoting and complete
orthogonal decomposition for rank-deficient minimum-norm least squares. Caller
rank tolerances remain explicit. This extends and moves the earlier Batch 4 QR
plan; independent verification is required before Batch 6 exists.

### Batch 6: dense symmetric eigenvalues and SVD

Develop symmetric Jacobi and direct one-sided Jacobi SVD, with bounded work,
convergence reports and nullspace/subspace bases. Test accuracy under explicit
assumptions; use these routines as additional cross-checks while retaining
independent references for rank decisions.

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