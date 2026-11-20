# Planned numerical batches

State: Batches 3 (3a through 3d), 4, 5 and 6 completed.
[Batch 6 evidence](batch6.md) records the Jacobi eigenvalue and singular value routines.
[Batch 5 evidence](batch5.md) records the orthogonal decompositions.
[Batch 4 evidence](batch4.md) records the dense indefinite factorization.
[Batch 3d evidence](batch3d.md) records the condition estimate.
[Batch 3c evidence](batch3c.md) records the Gram products.
[Batch 3b evidence](batch3b.md) records buffered-product delivery.
[Batch 3a evidence](batch3a.md) records the implemented contract and measurements.
Updated 2026-11-19 after application review.
[ROADMAP.md](../ROADMAP.md) owns the sequence; this page records contracts and
acceptance gates. Proposed operation names are not yet public-interface promises.
Application measurements motivate the work; they are not library benchmark results.
The [current interface](api.md) documents the new compact analysis and entry budgets.

## Batch 3: sparse analysis and repeated operations

Deliver the following as separately reviewed steps, preserving existing mathematical
semantics. Begin with the graph rewrite; buffer products can be reviewed separately
before the Gram refill machinery. Condition estimation follows the factor contract.

### 3a. Scalable symbolic analysis

Implementation note: column counts drive factor allocation. Tree and postorder
invariants are checked independently from exported L; production tree storage
waits for a traversal consumer and does not change the minimum-degree permutation.

Replace the unconditional square bitmap with compact graph storage. Preserve exact
minimum-degree selection and smallest-original-index tie breaking. Storage should
follow vertices, input structure and accumulated fill, rather than allocating an
n-by-n graph regardless of sparsity. Fill itself can be quadratic; neither linear
time nor linear memory is promised for every input.

Derive elimination-tree, column-count and traversal data where they support factor
allocation and traversal. A tree postorder may schedule internal work, but must not
silently change the documented final permutation. Any proposal to change ordering
or floating accumulation order needs an explicit contract review.

Replace the fixed 4096 limit only after measuring larger sparse families. Retain
overflow-safe index/size checks and explicit resource refusals before growth; lifting
a dimension cap does not mean attempting unlimited allocation. Avoid keeping both
a full explicit fill graph and redundant factor structure without measuring why.

Acceptance: independent set-graph ordering/fill checks; deterministic ties and
permutations; empty/disconnected, chain, grid, irregular and dense-fill families;
measured peak memory and analysis time, including sizes beyond 4096. For the old
roster, require equivalent structure and unchanged numeric artifacts when the
arithmetic order is unchanged. No fixed benchmark promise before the brief supplies
sizes and resource budgets.

### 3b. Products into caller storage

Completed: the three-argument Void operations use no scratch. Preflight refusals
preserve output; arithmetic refusal leaves it unspecified. See the current
[caller-buffer contract](primitives.md#caller-provided-output).

Add matvec! and transpose-matvec! using caller output and, where needed, explicitly
prepared scratch. Keep allocating convenience operations on the same mathematical
implementation. State overwrite semantics, dimension checks, nonfinite handling,
input/output alias rules and output state on refusal before freezing signatures.
Do not promise both aliasing and unchanged-on-failure output without enough scratch
to implement those promises.

Acceptance: independent dense/exact products in both directions, rectangular and
empty shapes, structural zeros, allowed and forbidden aliases, all refusal paths,
and measured zero allocation after preparation through a real Coalton caller.
Centralize the standard finite? inline hint if it simplifies the shared internals;
retain actual allocation checks at every new repeated operation.

### 3c. Fixed-pattern Gram products

Completed: see [the contract](gram.md) and [the batch record](batch3c.md).

Provide analysis plus repeated numeric refill for A^T A and A A^T. The plan is tied
to a checked canonical input pattern; values may change, structure may not. Specify
owned pattern storage, output lifetime and caller/workspace buffers explicitly.
Use one documented symmetric-triangle convention compatible with SPD analysis;
include every diagonal, including structural zeros for empty columns or rows.
Retain potential product entries when cancellation makes their current value zero.

Share the symbolic/refill mechanism between orientations. Measure symbolic product
storage and arithmetic work separately: caching every multiplication contribution
can cost much more than storing the resulting Gram pattern. Gram matrices can be
dense even when A is sparse. Make growth limits and failed-refill validity explicit.

Acceptance: exact small products, independently generated structural patterns,
repeated values and cancellation, both orientations, changed-pattern rejection,
empty dimensions, density growth, deterministic accumulation, and warmed allocation
measurements. Verify construction rounding; an exactly rank-deficient generating
matrix does not establish exact singularity of its rounded Gram output.

### 3d. Condition estimate using a prepared factor

Implementation note: delivered as `estimate-inverse-norm1!`, `shifted-norm1` and
`reciprocal-condition1`; [the contract](condition.md) is the current statement.
The estimate is tied to one factor generation; the matching norm is computed from
values the caller passes, because the library retains none.

Start with a named reciprocal condition estimate in a specified norm, such as
rcond_1, using factor solves and a caller-supplied or retained matching matrix norm.
Estimate the actual factored matrix, including any shift. Require matching numeric
generations, valid-factor checks, bounded iterations/workspace, failure/convergence
status and a documented cost. Do not form an explicit inverse merely for this task.
The factor diagonal ratio is not a matrix condition number. [Factor-based condition
estimation conventions](https://netlib.org/lapack/lug/node38.html).

This is not automatically the smallest/largest singular-value ratio. For a full-rank
matrix and its nonsingular, unshifted Gram matrix in exact arithmetic,
rcond_2(Gram) = (sigma_min(A)/sigma_max(A))^2. A 1-norm estimate cannot simply replace
the 2-norm quantity; a shifted factor describes a different matrix. Specify whether
zero singular values or only the nonzero spectrum are intended. If the application
requires a spectral ratio, keep that requirement visible for a later appropriate
estimator or decomposition; do not relabel the cheaper diagnostic.

Acceptance: diagonal/exact small references, scaled and nearly singular families,
comparison against independently computed condition numbers in the same norm,
estimator failure cases and counted solve calls. An estimate is not a certified
bound unless an error analysis and tests support that stronger contract.

## Batch 4: dense symmetric indefinite factorization and inertia

Implementation note: delivered as `coalton-linalg/dense/indefinite-ldlt`;
[the contract](indefinite-ldlt.md) is the current statement and carries the derivation.
The pivot search is the bounded variant that follows growing off-diagonal maxima to a
mutual maximum. Zero pivots are counted only when exactly zero in the computed arithmetic.

Build a separate symmetric-indefinite LDL^T interface, initially dense F64, using
permutations and stable 1-by-1 / 2-by-2 pivot selection of the Bunch-Kaufman family.
Keep the existing sparse SPD contract unchanged. Introduce only the dense storage
and orthogonal/numeric helpers required by concrete consumers, then reuse them.
Record a prose derivation before implementation. [Block-diagonal factor convention](https://netlib.org/lapack/lug/node38.html).

Return inertia as numbers of positive, negative and zero eigenvalues of the block
D, counting both eigenvalues of each 2-by-2 block. Counting pivot blocks or signs
of their diagonal entries is wrong: [[0,1],[1,0]] contributes one positive and one
negative value. Separate a completed singular factor/inertia result from whether
that factor can solve a system; report breakdown rather than invent missing counts.

By congruence, exact inertia of G-tI counts eigenvalues of G strictly below, equal
to and above t. Floating factors define an approximation whose reconstruction
and backward error must be checked; computed block signs alone do not establish
that it is close enough to preserve threshold counts. Near-threshold counts need
a documented reliability/error contract and caller threshold, not an exact rank
certificate. Separate numerical pivot selection for stability from the
application's definition of numerical rank.

For J of shape m-by-n, variable nullity is n-rank(J). Inertia of J J^T concerns an
m-dimensional spectrum and does not directly equal that nullity. For example,
J=[1,0] has a one-dimensional nullspace but J J^T=[1] has no zero eigenvalue.
Translate dimensions and equality-at-threshold policy explicitly. A nonnegative
singular-value cutoff tau corresponds to a Gram cutoff tau^2 in exact arithmetic,
with overflow/underflow and Gram-rounding limitations accounted for.

The augmented matrix [[I,J^T],[J,0]] with right-hand side [0,b] expresses minimum
Euclidean norm subject to Jx=b. With r=rank(J), its exact inertia is (n,r,m-r):
it is congruent to diag(I_n,-J J^T). It is nonsingular when J has full row rank.
With dependent rows it is singular; compatible b still has a unique minimum-norm
x but nonunique multipliers, while incompatible b makes the equality problem
infeasible. A general indefinite solve is not automatically rank-deficient least
squares. This follows from the [KKT invertibility conditions](https://web.stanford.edu/~boyd/vmls/vmls.pdf),
section 16.2, specialized to an identity objective matrix. Dense delivery verifies
the mathematics; it does not yet remove a large sparse augmented-system bottleneck.

Acceptance: reconstruction with permutations/block D, residuals and scaling;
2-by-2 cases whose inertia differs from diagonal-sign counts; exactly constructed
congruences with known inertia; singular blocks and clustered thresholds; compatible
and incompatible augmented systems; independent references and planted controls.

## Batch 5: orthogonal decompositions and minimum-norm least squares

Implementation note: delivered as `coalton-linalg/dense/orthogonal`;
[the contract](orthogonal.md) is the current statement and carries the derivation. Column
pivoting is always on; plain QR, full-rank least squares and the full-row-rank wide case
are selected ranks of the one factorization, not separate entry points. The rank
criterion is one absolute, strict threshold on the diagonal of R, supplied by the caller.

Stages: stable Householder transformations and QR; full-column-rank least squares
(and transpose/LQ treatment of full-row-rank wide systems); column-pivoted QR;
complete orthogonal decomposition (COD) for rank-deficient minimum-norm solutions.
State the numerical rank criterion, scaling, equality boundary and caller tolerance.
Column pivoting is not guaranteed to reveal numerical rank reliably for every
matrix; diagonal entries of R are not singular values. A basic pivoted-QR solution
need not have minimum norm. [QR rank conventions](https://netlib.org/lapack/lug/node42.html)
and [why the complete factorization is needed](https://netlib.org/lapack/lug/node43.html).

Keep reflector application and transpose application available without forcing a
full explicit Q. Define permutation direction and distinguish residual minimization
from minimum norm among residual minimizers. For rank truncation, document the
nearby problem actually solved. Accept arbitrary rectangular shapes, including
rank-zero and empty cases, with explicit ownership and failure semantics.

Acceptance: reconstruction, orthogonality, scaled residual optimality and solution
orthogonality to the selected nullspace; tall/wide and consistent/inconsistent
systems; known-rank exact inputs, tolerance transitions and adversarial rank-estimation
families such as the [Kahan matrix](https://math.nist.gov/MatrixMarket/deli/Kahan/information.html);
independent high-precision or exact checks where applicable. Do not wait for Batch 6 to validate this batch.

## Batch 6: dense symmetric eigenproblems and Jacobi SVD

Implementation note: delivered as `coalton-linalg/dense/jacobi`;
[the contract](jacobi.md) is the current statement and carries the derivation. The right
singular factor is always complete, n-by-n, so the nullspace basis needs no separate
completion; the price is that a wide matrix costs far more than its transpose. Absolute
accuracy is tested by exact inertia brackets; relative accuracy of small values is
observed on two graded matrices and not promised.

Develop symmetric Jacobi eigenvalue/eigenvector operations and a direct one-sided
Jacobi SVD. Do not obtain the SVD by first forming a Gram matrix, which squares the
conditioning and can lose information about small singular values. Start with
small-to-moderate dense sizes and measured size/work limits, explicit scaling,
stopping rules, sweep limits and nonconvergence results.

Jacobi's accuracy advantages require hypotheses and a careful implementation;
relative accuracy of arbitrarily small eigenvalues of every symmetric matrix is
not guaranteed. [Demmel and Veselic's accuracy result](https://epubs.siam.org/doi/abs/10.1137/0613074)
distinguishes symmetric positive definite eigenproblems and one-sided SVD results.
Preserve absolute-error and residual reporting near zero, where relative error
can be undefined or uninformative.

Provide nullspace/subspace bases under a caller-selected rank policy. Specify
economy versus complete factors: economy factors can omit nullspace directions
caused by rectangular shape, so a full nullspace operation needs an explicit or
implicit completion. Repeated or clustered values make individual basis vectors
nonunique: compare residuals,
orthogonality and subspace projectors, not raw vectors or their signs. These
routines add a library cross-check; they do not become the sole reference oracle
for every rank decision. Retain independently derived exact/high-precision tests,
known spectra, planted defects and externally held-back cases.

## Later work and scope

Only after measurements justify it: Lanczos for extreme eigenvalues; LSMR for
large rectangular least squares; conjugate gradients for SPD systems; sparse QR,
sparse indefinite factorization, supernodal methods and numeric types beyond F64.
Do not send indefinite or general rectangular systems to ordinary CG. Any claim
that augmented systems scale to large assemblies needs a sparse/iterative path
and separate measurements; dense Batch 4 alone cannot establish it.

Geometry, residual functions, dual-number Jacobians, nonlinear iteration, model
rank policy, solver routing and application integration stay outside this library.
Existing [clean-room provenance](provenance.md) and serial verification rules apply.
Source references above support mathematical contracts, not implementation reuse.

## Input requested from the next brief

Start with Batch 3, not one monolithic brief for all future algorithms. Supply:

- Generic matrix families with dimensions, input/factor/Gram nonzeros, repeated-value
  workload, storage convention and reproducible seeds; identify redistribution rights.
- Desired measurements and explicit budgets for analysis time, peak memory, refill
  time and allocations. Distinguish measured requirements from proposed targets.
- Exact meaning of each condition/rank quantity: matrix, norm, shift, tolerance,
  scaling and threshold equality. Mark uncertain classifications explicitly.
- Expected properties established from serialized F64 values, with independent
  exact/high-precision evidence where claimed. Include positive and negative controls.
- A separately maintained held-back suite and a recorded outcome format. Its owner
  may run it externally; no hidden case-specific acceptance rules enter the library.

Later dense briefs should include adversarial pivots, singular/clustered spectra,
rank-deficient and incompatible right-hand sides, and basis-invariant comparisons.