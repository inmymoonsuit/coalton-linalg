# Batch 6: symmetric eigenvalues and singular values by Jacobi rotations

State: completed; verified at safety 1 and safety 0 on 2026-11-19, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract and the derivation are in
[jacobi.md](jacobi.md).

## What was delivered

A new package, `coalton-linalg/dense/jacobi`:

- `eigen!`: all eigenvalues, ascending, and eigenvectors of a dense symmetric matrix by
  cyclic Jacobi rotations.
- `svd!`: all singular values, descending, and both sets of vectors of a dense matrix of
  any shape by one-sided Jacobi rotations on its columns. AᵀA is not formed. The right
  factor is always a complete n-by-n orthogonal matrix, so its trailing columns are a
  nullspace basis of the right dimension for tall, square and wide matrices alike.
- `rank-above`: a count of singular values against one absolute, strict threshold.
- A sweep limit as an argument, with a default of 60, and `JacobiNotConverged` when it is
  reached.

A compiled example, [singular-rank.ct](../examples/singular-rank.ct), reads the rank of a
system with a repeated equation from its singular values. The prose derivation in the
contract page was written before the implementation.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`, with every earlier roster unchanged. For this batch:

- **Every eigenvalue and singular value is bracketed exactly.** For a computed value λ and
  a tolerance δ, an exact symmetric elimination in rational arithmetic counts the
  eigenvalues of the stored matrix below λ − δ and at or below λ + δ (Sylvester's law of
  inertia). Value k must have at most k below and at least k + 1 at or below. This uses no
  eigenvalue routine. For singular values the counted matrix is the exact AᵀA and the
  shifts are (σ ± δ)². δ is 16 units of roundoff times n times the size of the matrix.
- **Eigenvalues, 136 matrices** up to order 8: named cases, 96 seeded integer matrices, 27
  positive semidefinite matrices of every rank, and one graded positive definite matrix.
  V is orthogonal to 23.1 units of roundoff against an asserted 8(n + 1); A V − V Λ is at
  most 2.5 units of roundoff times n max|A| against an asserted 16; the most sweeps
  taken was 8. The graded matrix's five eigenvalues, spanning 80 binary orders, are each
  bracketed to a relative 2⁻⁴⁰.
- **Singular values, 144 matrices** in sixteen shapes up to 7-by-5 and 5-by-8 at every rank
  from zero to full, with empty shapes: V is orthogonal to 18.1 units; A V − U Σ is at most
  1.5 units of roundoff times n times the largest column norm; the kept columns of U are
  orthonormal; a zero value has a zero left vector; the rank at 10⁻⁹ of the largest value
  equals the exact rank; the most sweeps taken was 9. A matrix with graded columns has its
  three singular values, spanning 40 binary orders, each bracketed to a relative 2⁻³⁰.
- **Identities.** Scaling by 2⁻⁶⁰⁰, 2⁻³⁰, 2²⁸ and 2⁶⁰⁰ scales the values exactly and leaves
  every vector unchanged bit for bit, for both routines. A repeated call repeats every bit.
  The transpose has the same singular values. For a positive semidefinite matrix the
  singular values are the eigenvalues.
- **The two rank-revealing routines agree.** On five rank-deficient shapes the minimum-norm
  solution assembled from the singular vectors equals the one from
  [the orthogonal decomposition](orthogonal.md) to 10⁻¹⁰, and the two ranks are equal.
- **Lifetime and refusals.** A sweep limit of 0 or 1 gives `JacobiNotConverged` and leaves
  nothing standing; the reported sweep count, given back as the limit, converges again;
  every refusal withdraws an earlier result; overflowing results are refusals; thresholds
  are strict; storage limits are typed errors, including the right factor's n² entries.
- **Encapsulation.** Five external compile probes that must be rejected and one public
  control that must compile and run.
- **Allocation.** Compiled Coalton callers recomputing the eigenvalues of a 10-by-10 matrix
  and the singular values of an 8-by-10 matrix of rank 5: zero bytes over 2,000 calls after
  10,000 warm-up calls, at both policies.

**Planted source controls.** Each of these edits to `jacobi.ct` was applied alone and the
focused tests run; every one failed, and the restored source passed.

| Planted defect | Caught |
| --- | --- |
| The larger root of the rotation equation is taken | yes |
| The rotation ignores the sign of ζ | yes |
| Eigenvalues: the rotated diagonal entry is not set | yes |
| Eigenvalues: the vectors are not rotated | yes |
| Eigenvalues: the two rows are not rotated with the two columns | yes |
| Eigenvalues: values are returned unsorted | yes |
| Eigenvalues: sorting moves values without their vectors | yes |
| Eigenvalues: the convergence test is 10⁻⁸ of the two diagonal entries | yes |
| Eigenvalues: reaching the sweep limit is reported as success | yes |
| Eigenvalues: the matrix is never scaled | yes |
| Eigenvalues: values are not scaled back | yes |
| Singular values: left vectors are not normalized | yes |
| Singular values: the right factor is not rotated | yes |
| Singular values: values are not scaled back | yes |
| Singular values: the matrix is never scaled | yes |
| Singular values: `rank-above` counts a value equal to the threshold | yes |
| Singular values: sorting moves values without their right vectors | yes |
| Singular values: tiny columns are rotated like any other | yes |
| Singular values: reaching the sweep limit is reported as success | yes |
| Singular values: the convergence test is 10⁻⁸ of the two column norms | yes |

## Two things the controls found

- **Three planted edits failed for the wrong reason.** They left a variable unused, the
  compiler refused the file, and the runner counted the missing pass line as a catch. The
  logs of every earlier batch's controls were then searched for a compile failure, and none
  was found. The three were re-planted so that the source compiles.
- **A safeguard nothing needed was removed.** The eigenvalue routine first skipped any entry
  below 2⁻¹⁰⁴ of the largest. Re-planted properly, the control that disabled this floor was
  missed: every test passed without it, sweep bounds included. The floor was taken out, and
  the routine now has the one relative test its contract states; positive semidefinite
  matrices of order 200 and rank 66 converge in 16 sweeps without it. The same floor in the
  singular value routine is needed: without it a 2-by-4 matrix does not converge, because
  its two surplus columns cannot be made orthogonal to anything and only shrink.

## Measurements

Observations, not targets. Nine symmetric matrices generated outside the library: `random`
is a symmetrized standard normal matrix; `lowrank` is BᵀB with n/3 rows; `spread` is
Q diag(1 … 10⁻¹²) Qᵀ with geometric spacing. Twelve rectangular matrices from the
[Batch 5 record](batch5.md): `full` with singular values from 1 to 10⁻³, `half` with half as
many from 1 to 10⁻⁶. The library's outputs were written to files and checked by a separate
NumPy program.

- Value error: the largest difference from the reference values, in units of roundoff
  times n times the largest entry (eigenvalues) or the largest column norm (singular values).
- Residual: max|A V − V Λ| or max|A V − U Σ| in the same units. Orthogonality:
  max|VᵀV − I| in units of roundoff.
- Times are the best of 3 after one warm-up call, in milliseconds; bytes are over 3 warmed
  calls. One process at a time under a 1024 MB heap cap.
- Reference: `numpy.linalg.eigh` and `scipy.linalg.svd` with the `gesvd` driver, economy
  factors, best of 3, called from Python over an optimized, possibly multithreaded BLAS.
  Reference points for scale, not like-for-like comparisons.

| Symmetric matrix | n | Sweeps | Value error | Residual | Orthogonality | `eigen!`, s1 | `eigen!`, s0 | Bytes | Reference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| random | 50 | 9 | 1.64 | 0.79 | 102 | 1.73 | 1.46 | 0 | 0.154 |
| random | 100 | 10 | 1.92 | 0.44 | 228 | 13.64 | 11.80 | 0 | 1.310 |
| random | 200 | 10 | 1.47 | 0.43 | 448 | 116.91 | 103.60 | 0 | 4.233 |
| lowrank | 50 | 14 | 0.58 | 0.71 | 114 | 2.26 | 1.91 | 0 | 0.125 |
| lowrank | 100 | 15 | 0.51 | 0.19 | 232 | 17.69 | 15.67 | 0 | 1.243 |
| lowrank | 200 | 16 | 0.66 | 0.28 | 532 | 154.17 | 136.52 | 0 | 3.228 |
| spread | 50 | 17 | 1.00 | 0.65 | 114 | 2.60 | 2.13 | 0 | 0.122 |
| spread | 100 | 21 | 0.49 | 0.31 | 270 | 20.52 | 17.44 | 0 | 2.979 |
| spread | 200 | 23 | 1.91 | 0.55 | 512 | 172.09 | 152.59 | 0 | 4.451 |

| Shape | Rank | Sweeps | Value error | Residual | Orthogonality | Null A v | `svd!`, s1 | `svd!`, s0 | Bytes | Reference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 200 × 50 | 50 | 13 | 2.75 | 0.15 | 120 | — | 6.70 | 6.04 | 0 | 0.653 |
| 200 × 50 | 25 | 20 | 1.33 | 0.11 | 124 | 9.9e-17 | 10.05 | 9.10 | 0 | 0.569 |
| 400 × 100 | 100 | 14 | 2.87 | 0.12 | 236 | — | 54.59 | 49.40 | 0 | 17.940 |
| 400 × 100 | 50 | 24 | 3.56 | 0.11 | 308 | 1.3e-16 | 91.37 | 85.91 | 0 | 6.205 |
| 100 × 100 | 100 | 16 | 2.55 | 0.17 | 268 | — | 21.05 | 19.61 | 0 | 2.202 |
| 100 × 100 | 50 | 25 | 3.22 | 0.16 | 318 | 2.1e-16 | 32.74 | 30.26 | 0 | 1.741 |
| 200 × 200 | 200 | 17 | 3.68 | 0.17 | 526 | — | 168.90 | 153.04 | 0 | 12.619 |
| 200 × 200 | 100 | 30 | 4.16 | 0.11 | 800 | 2.4e-16 | 292.95 | 267.70 | 0 | 9.461 |
| 50 × 200 | 50 | 14 | 4.25 | 0.21 | 652 | 2.5e-16 | 87.32 | 79.80 | 0 | 0.779 |
| 50 × 200 | 25 | 17 | 4.43 | 0.20 | 852 | 2.3e-16 | 112.30 | 100.45 | 0 | 0.755 |
| 100 × 400 | 100 | 14 | 4.79 | 0.14 | 1376 | 2.9e-16 | 710.49 | 640.08 | 0 | 5.946 |
| 100 × 400 | 50 | 22 | 5.38 | 0.24 | 2298 | 2.9e-16 | 1074.54 | 952.86 | 0 | 5.935 |

- **Values.** Every eigenvalue and singular value is within 5.4 units of the reference.
  Relative to the reference, the kept singular values agree to 7 × 10⁻¹⁴ on the `full`
  matrices and 6 × 10⁻¹² on the `half` ones, whose condition number is a thousand times
  larger; at that level the comparison also measures the reference.
- **Rank.** The count above 10⁻¹⁰ of the largest value equals the designed rank on all
  twelve. For the wide matrices the n − m surplus values are at most 4.6 × 10⁻³³.
- **Nullspace.** A times every trailing column of V is at most 2.9 × 10⁻¹⁶.
- **Orthogonality** of V degrades with the dimension and the sweep count, to about six
  units of roundoff a column at worst here.
- **Sweeps.** The random symmetric matrices and the `full` rectangular ones take 9 to 17
  sweeps. The spread spectrum takes up to 23 and the rank-deficient matrices up to 30,
  because columns at rounding level are orthogonalized among themselves as well. The
  default limit is 60.
- **Speed.** These routines are slow against the reference: 7 to 48 times for eigenvalues
  at orders 50 to 200, 3 to 31 times for tall and square singular value problems, and more
  than a hundred times for wide ones, whose cost grows with the square of the number of
  columns. A caller needing only the singular values of a wide matrix should pass its
  transpose. Safety 0 is 6% to 18% faster than safety 1, and the two policies produce
  identical output files.

## After review

**A held-back set found a wrong answer returned as a success, and it is fixed.** The
reviewer ran 53 runs through the driver below, and 52 passed.
The one failure was the symmetric matrix with diagonal 10³⁰⁸ and −10³⁰⁸ and off-diagonal
10³⁰⁷. The difference of its diagonal entries overflowed, the rotation computed from the
infinite ratio was the identity, and the routine still set the off-diagonal entry to zero
and counted the rotation. It returned ±10³⁰⁸ with unit vectors as `JacobiSuccess`, where the
eigenvalues are ±1.005 × 10³⁰⁸: each value 0.5% off, and nothing reported.

The eigenvalue routine now scales its matrix by an exact power of two first, as the
singular value routine already did, so that difference cannot overflow. The failing matrix,
a second with equal coupling whose eigenvalues ±1.414 × 10³⁰⁸ are representable, and a
third near the smallest normal number are in the exact roster, each value bracketed by the
exact inertia count. Two planted controls, one disabling the scaling and one the scaling
back, are caught. The rosters that passed before the fix did not hold a matrix whose
diagonal difference overflows; the scaling tests went only to 2⁶⁰⁰.

The same set measured, from the files alone: VᵀV within 464 units of roundoff at n = 200 and
within 3024 at 402 columns after 32 sweeps; A V within 4 and 5 units of the size of A; every
rank line equal to the strict count, and equal to the proved exact rank on every integer
system; and on a built matrix with singular values 0.04 decades above and 0.05 below a
relative threshold, a count of 11 from the singular values where the diagonal of R in
[the orthogonal decomposition](orthogonal.md) gave 12.

The reviewer then ran the fixed routine on the first set again and on 13 new symmetric
systems at the ends of the range: 66 of 66 runs passed. The failing matrix returns values
within one unit of roundoff of the exact ones. The new systems include diagonal entries of
mixed sign near 10³⁰⁸, a matrix whose larger eigenvalue is not representable and is refused
as `JacobiNonFinite`, the smallest subnormal in every entry, a diagonal matrix 616 decimal
orders across, and an integer matrix scaled by 2⁹⁰⁰, 2¹⁰⁰⁰, 2⁻⁹⁰⁰ and 2⁻¹⁰⁶⁰; the first three
scalings give exactly scaled values and vectors identical bit for bit. At 2⁻¹⁰⁶⁰ the
eigenvalues are subnormal and each returned value is the exact one rounded to the nearest
subnormal, which the contract now states.

An artifact driver, `scripts/check.ps1 -Task Jacobi`, was added at the reviewer's request so
that a held-back set can be checked from files alone. It writes the values, the vectors and
the sweep count of each run by sweep limit, and for singular value problems the rank at
each named threshold; [the contract page](jacobi.md#artifact-driver) has the protocol. Its
regressions cover exact eigenvalue and singular value files, a named non-convergence at a
sweep limit, strict rank lines, a named overflow refusal, one name used for both kinds, and
removal of stale artifacts.

## Limits

- Small and moderate sizes only, and an order of magnitude or more slower than
  tridiagonalizing and bidiagonalizing methods. The reason to have them is that they are
  simple, accurate, and independent of the factorizations they cross-check.
- The singular values decide rank better than the diagonal of R does, and the decision is
  still the caller's threshold.
- Absolute accuracy relative to the size of the matrix is tested; relative accuracy of
  small values is observed on two graded matrices and not promised.
- Eigenvectors and singular vectors of repeated or clustered values are not unique.
- One host and finite rosters establish neither universal accuracy nor portability beyond
  the supported runtime.
