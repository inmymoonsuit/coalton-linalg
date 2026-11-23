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

## Batches 3 through 14

Batches 3 through 14 are completed, Batch 3 in four steps 3a through 3d, with one addition,
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
A held-back review later found it inexact where an intermediate product left the normal
range; [Batch 9a](docs/batch9a.md) corrects that and names the rows and columns that no
factor reaches.

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

### Batch 11: dependent rows, and a factorization that leaves small pivots out

**Completed.** `factor-skipping!` factors a symmetric matrix in the analyzed order and
leaves out each index whose pivot is at or below a floor, so that the factor of the kept
principal submatrix stands; `solve-skipping!` solves with it and `skipped-null-vector!`
gives an index left out its combination of the kept ones.
`coalton-linalg/sparse/dependent-rows` applies that to `A Aᵀ`: the rows to keep, a left
null vector for each row left out, refined against A, a bound on what each leaves that
holds for the exact product, and from those a bound on the singular values beyond the
kept count. That the kept rows are independent is proved by the existing certificate on
the whole matrix's own plan. [The design](docs/design-dependent-rows.md) came first and
was read by the consumer that asked for it; [the batch record](docs/batch11.md) holds the
controls, and [the contract](docs/dependent-rows.md) says what is proved and where the
method stops resolving.

### Batch 12: sparse LU with partial pivoting

**Completed.** `coalton-linalg/sparse/lu` factors a square sparse matrix as `P A Q = L U`,
choosing each pivot from the values: the largest magnitude in the column, or the
matrix's diagonal entry where it is within a threshold of that. Each column costs what
it changes, found by a depth-first search over the structure of L. The factors' storage
grows up to a budget. It solves with the matrix and with its transpose, measures its
growth, and bounds its backward error by one number that holds with underflow.
[The batch record](docs/batch12.md) holds the controls; [the contract](docs/lu.md) says
what is claimed and what is not built.

### Batch 13: conjugate gradients

**Completed.** `coalton-linalg/conjugate-gradient` solves a symmetric positive definite
system through nothing but products with its operator, with or without a
preconditioner, which may be any solve this library has. Convergence is reported on the
residual formed from the operator and not on the method's running one, and a direction
or a preconditioned residual that shows the operator or the preconditioner not positive
definite ends the call by name. [The batch record](docs/batch13.md) holds the controls;
[the contract](docs/conjugate-gradient.md) says what a converged call does and does not
establish.

### Batch 14: an outside audit at the extremes of the range

**Completed.** An audit contributed as a pull request found eight defects at the ends of
the floating-point range and at an aliasing the tests had not offered: a right-hand side
overwritten when it shares the iterate's array, convergence declared on a residual whose
square underflowed, a dense two-by-two pivot solve that returned zeros at large scale, a
small singular value reported as zero, a reflector wrong for a subnormal column, an
unchecked overflowing multiplier in sparse LU, and two smaller ones. Each came with a
test that fails on the earlier sources. Integrating it showed that the pinned compiler
miscompiles a function holding a cell when inlined from another package; the code was
arranged to avoid that and the build now fails on a compiler warning.
[The record](docs/batch14.md) has each finding and how it was confirmed.

### Batch 15: the audit's two questions asked of every call

**Completed.** The audit's findings were answers to questions nobody had asked of the
routines. Two of them are now asked of every exported call at once, by rosters that read
their members from the sources: what a call does when two of its arrays are one array,
and what it does with the same problem at every scale of the format. The first found
nothing. The second found one silent wrong answer, in the search for dependent rows at
very small scales, which is now refused by name, and measured the range of every
routine; conjugate gradients and dependent rows reach less than the whole range, and
their contracts now say how much. [The rosters' page](docs/rosters.md) is the contract
and [the record](docs/batch15.md) the evidence.

### Batch 16: the range of conjugate gradients and of dependent rows

**Completed.** The two limits Batch 15 measured are removed. Conjugate gradients carry the
iteration on the residual divided by a power of two, and the search for dependent rows is
made on the matrix divided by one, so a right-hand side, an operator or a matrix of any
size the format holds is solved as its ordinary-sized counterpart is, bit for bit. A row
too small to square beside the largest entry, under a threshold of zero, is refused by a
status of its own, `RowsUnderflow`. [The record](docs/batch16.md) has the measurements
before and after.

### Batch 17: damped least squares

**Completed.** `coalton-linalg/sparse/damped`: for a loop that tries one damping after
another on the same matrix, the step `−(Jᵀ J + μ D)⁻¹ Jᵀ r` with one analysis for the
pattern, the values of J set once an iterate, and a factorization and a step for each μ,
nothing analyzed again and nothing allocated. D is a positive diagonal the caller gives,
or the columns' own squared lengths with a floor. The matrix and the residual are each
divided by a power of two, so the step at any scale is the step at ordinary scale.
[The record](docs/batch17.md) has the exact references and the controls;
[the contract](docs/damped.md) says what the package leaves to its caller.

### Batch 18a: kept rows and the certificate

**Completed.** `coalton-linalg/sparse/rank-deficient`, the first of the two parts of a
step on a rank-deficient matrix. A search is made once; the vector of each row left out
is kept, with entries at rounding level dropped and the bound proved for the vector as
kept; for each later matrix of the pattern the bounds are proved again by a pass over
the rows the vectors touch, with no search and no allocation; and a certificate is
issued, or refused by name, that the gap bound is below one floor and every singular
value of the kept rows above another. The same package bounds the largest singular
value from both sides. It solves nothing yet. [The contract](docs/rank-deficient.md)
says what is proved and [the record](docs/batch18a.md) has the exact references.

### Batch 18b: the step on a rank-deficient matrix

**Completed.** `step!` in `coalton-linalg/sparse/rank-deficient`: for the matrix a
certificate stands on, the minimum-norm least squares step of the matrix whose rows left
out are replaced by their kept combinations, with two proved bounds, on what the solves
left and on the distance to the step of the matrix truncated at the kept count. The
solve on the kept rows is through their augmented matrix, on the certificate's own
analysis, refined, and its bound is read from the residuals and the bound the
certificate proved, with no rate of convergence assumed. [The contract](docs/rank-deficient.md#the-step)
says what is proved and [the record](docs/batch18b.md) has the exact references.

### Next

From 2026-11-22 the order is that of a consumer that is moving a nonlinear solver onto
the library, agreed with it in advance, plan before code:

1. A step on a rank-deficient matrix: the minimum-norm least squares step that a
   truncated singular value decomposition gives at a relative floor, with a certificate
   that no singular value lies in a band about the floor, or a refusal by name. Built on
   [dependent rows](docs/dependent-rows.md). **[The design](docs/design-rank-deficient-step.md)
   is reviewed and measured. The certificate and the kept vectors are Batch 18a and the
   step with its bounds Batch 18b. Left: the refresh of one kept vector from current
   values, for rows whose entries change at every call. The kept rows alone, as an
   option, are not built; the consumer will not call it.** The same note covers
   certified bases: the supports of left null vectors, and a right null-space basis with
   a proved residual.
2. Those bases, as their own batch after the step.
3. Least squares with bounds on linear functionals of the unknown, by an active set whose
   rows' rank is proved or refused.
4. Supernodal or multifrontal factorization, only after measurements say where real
   problems spend their time.

Kept on the roadmap for the library's own users, though that consumer has no use for it:

5. Iterative least squares and minimum norm on an operator, with damping. Conjugate
   gradients for the definite case is Batch 13.
6. A few eigenvalues of a sparse symmetric matrix: Lanczos, shift-invert through the
   signed factorization, and the inertia as a count that none was missed.
7. Continuous integration on the pinned toolchain.

Dependent rows, first in the ranking before this one, is Batch 11, and sparse LU, second, is Batch 12.
Two things are left open by those batches and wait for a consumer's use:

- For dependent rows, a second method: a triangular factor of Aᵀ by rotations, which
  resolves distances the first cannot. It is designed and not built, and is wanted
  only where a tolerance is finer than about 10⁻⁷ of a row's size.
- For sparse LU, a refactorization on a fixed pivot sequence and a condition estimate.

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