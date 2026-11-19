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
Cached-result and workspace contracts remain explicit; private handles do not
enforce read-only array borrowing.

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

## Batch 3: scalable symbolic factorization

Measure ordering/fill/memory on larger sparse families, replace the unconditional
n² adjacency bitmap with sparse graph storage, and preserve the deterministic
ordering contract. Fill can itself become dense; avoid an unconditional linear
memory promise. Deliver a reviewed first release candidate with measured limits.

## Batch 4: QR and full-column-rank least squares

Build independently derived Householder transformations, bounded dense QR, and
full-column-rank rectangular least-squares solving. Verify reconstruction,
orthogonality, and optimality through independent references. Rank-deficient
minimum-norm solving and sparse QR require subsequent explicit contracts.

## Release work

Configure CI workers with the pinned dependencies;
and expand supported hosts only after measured validation. The current supported
and tested host is Windows/x86-64 SBCL. The initial numeric type remains F64.

Applications own geometric Jacobians, nonlinear iteration, rank thresholds,
units, and integration. No SolveSpace promotion decision is made here.