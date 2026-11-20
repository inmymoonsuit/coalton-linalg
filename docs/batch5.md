# Batch 5: orthogonal decompositions and minimum-norm least squares

State: completed; verified at safety 1 and safety 0 on 2026-11-19, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract and the derivation are in
[orthogonal.md](orthogonal.md).

## What was delivered

A new package, `coalton-linalg/dense/orthogonal`, for a dense matrix of any shape:

- `factor!`: column-pivoted Householder QR, `A P = Q R`.
- `pivot-magnitude` and `rank-above`: the diagonal of R and a count against one absolute,
  strict threshold, so that the rank is the caller's decision.
- `set-rank!`: the complete orthogonal decomposition for the chosen rank, without a new
  factorization.
- `solve!`: the minimum-norm least-squares solution of the rank-truncated problem, with
  its residual norm.
- `nullspace-vector!`: an orthonormal basis of the selected nullspace.
- `apply-q!` and `apply-qt!`: Q and its transpose applied to a vector without forming Q.

Storage is refused by name above an entry budget, through a small private `dense/limits`
package that the dense indefinite factorization now shares. A compiled example,
[dependent-constraints.ct](../examples/dependent-constraints.ct), solves a system with a
repeated equation.

The prose derivation in the contract page was written before the implementation. One
scope decision is recorded there: pivoting is always on, and plain QR, full-rank least
squares and the full-row-rank case are selected ranks, not separate entry points.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`, with every earlier roster unchanged. For this batch:

- **Exact references, 215 decompositions**, 167 of them rank deficient: 14 named cases
  including empty shapes, and 201 seeded integer matrices in sixteen shapes up to
  7-by-5 and 5-by-8, each at every rank from zero to full. For each:
  - Q is read back one column at a time through `apply-q!` and multiplied out in exact
    rational arithmetic: QᵀQ is within 7.9 units of roundoff of the identity, against an
    asserted 8(m + 1), and Q R is within 5.7 units of roundoff times the largest column
    norm of A P, against an asserted 8(m + n + 1);
  - `apply-qt!` agrees with the transpose of what `apply-q!` applies;
  - the diagonal of R does not grow, and each diagonal entry dominates the rest of its
    row block, as pivoting by largest remaining norm requires;
  - the rank from `rank-above` at 10⁻⁹ of the largest diagonal entry equals the exact rank
    from an exact elimination;
  - the solution is within 2⁻³⁶ of the exact minimum-norm least-squares solution, computed
    from an exact rank factorization; the largest error was 127 units of roundoff;
  - the reported residual norm matches the exact one;
  - every nullspace vector is orthonormal to the others, is annihilated by A, and is
    orthogonal to the solution, each to 2⁻⁴⁰.
- **Rank decisions.** Thresholds at, between and beyond the diagonal entries of a matrix
  with entries 4, 2, 1, including equality; a rank below the exact one solving the
  truncated problem exactly; a rank above `step-count` and a rank reaching an exact zero
  pivot each refused by name.
- **The Kahan matrix of order 30.** No column is exchanged, the last diagonal entry of R is
  0.130, and the smallest singular value is at most 4.2 × 10⁻⁵ by the exact inverse of the
  stored matrix: a factor above 3000. Reconstruction and orthogonality hold as elsewhere.
- **Identities.** Scaling by 2⁻⁶⁰⁰, 2⁻³¹, 2²⁹ and 2⁶⁰⁰ leaves the permutation and the
  reflectors unchanged bit for bit and scales R and the solution exactly. A repeated
  factorization repeats every bit; selecting another rank and returning gives the same
  solution; Qᵀ undoes Q.
- **Lifetime and refusals.** Nothing before a factorization; a factor without a rank does
  not solve; every solve refusal leaves the output alone and withdraws the residual; a new
  factorization withdraws the rank; three kinds of refused factorization withdraw
  everything; an overflowing column norm and an overflowing solution are refusals; storage
  limits are typed errors.
- **Encapsulation.** Five external compile probes that must be rejected and one public
  control that must compile and run.
- **Allocation.** A compiled Coalton caller that factors a 9-by-12 matrix of rank 6, picks
  the rank by a relative threshold, completes, solves and reads a nullspace vector: zero
  bytes over 10,000 calls after 10,000 warm-up calls, at both policies.

**Planted source controls.** Each of these edits to `orthogonal.ct` was applied alone and
the focused tests run; every one failed, and the restored source passed.

| Planted defect | Caught |
| --- | --- |
| β takes the sign of α | yes |
| Norms are never scaled | yes |
| No column is exchanged | yes |
| The fused measurement trusts its unscaled sum at any size | yes |
| A skipped reflection leaves the remaining norms measured over one row too many | yes |
| The permutation is not recorded | yes |
| `rank-above` counts an entry equal to the threshold | yes |
| Completion leaves the rows above a reflection untouched | yes |
| The free unknowns are not set to zero | yes |
| The solution skips the reflections from the right | yes |
| The solution is not returned to the original column order | yes |
| The residual counts every entry of Qᵀb | yes |
| Q is applied in the order of its transpose | yes |
| A new factorization keeps the selected rank | yes |
| A rank reaching an exact zero pivot is accepted | yes |
| A refused solve keeps the earlier residual | yes |
| The nullspace vector starts from the wrong unit vector | yes |

## Two things found while building it

- **The Kahan matrix needs its perturbation.** As first written, without the customary
  small diagonal perturbation, rounding broke the exact ties among its column norms,
  columns were exchanged, and the factorization revealed the small singular value: at
  order 24 the last diagonal entry was 6.0 × 10⁻⁴ against a bound of 4.1 × 10⁻⁴. The test
  now perturbs the diagonal by a few rounding units and asserts that no column moves.
- **A double returned from a helper was boxed.** The first restructuring for speed
  returned each column's new norm as an F64 from the function that applies a reflection.
  The allocation check failed at 16 bytes a call. Writing the norm into the norms array
  inside that function removed it. Other F64-returning helpers in the same file, called
  from more than one place, did not allocate; the cause was not investigated further.

## Measurements

Observations, not targets. Eighteen matrices generated outside the library as U Σ Vᵀ with
orthogonal U and V from NumPy: `full` has singular values spaced geometrically from 1 to
10⁻³; `half` has half as many, from 1 to 10⁻⁶, and the rest zero up to rounding. Rank is
chosen at 10⁻¹⁰ of the largest diagonal entry of R. The library's outputs were written to
files and checked by a separate NumPy program.

- Solution error: the largest entry of the difference from the solution through a singular
  value decomposition truncated at the same relative 10⁻¹⁰, over the largest entry of that
  solution.
- Last / σ: the last kept diagonal entry of R over the singular value of the same index.
- Times are the best of 20 (5 for the larger shapes) after one warm-up call, in
  milliseconds. Bytes are over 20 warmed factor, rank and solve triples. One process at a
  time under a 1024 MB heap cap.
- Reference: `scipy.linalg.qr` with pivoting, R only, and `scipy.linalg.lstsq` with the
  `gelsy` driver, SciPy 1.18.1, best of 5, called from Python over an optimized, possibly
  multithreaded BLAS. Reference points for scale, not like-for-like comparisons.

| Shape | Rank | Solution error | Last / σ | `factor!`, s1 | `factor!`, s0 | `set-rank!`, s1 | `solve!`, s1 | Bytes | Reference QR | Reference solve |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 200 × 50 | 50 | 1.41e-13 | 2.351 | 0.3376 | 0.3327 | 0.0019 | 0.0142 | 0 | 0.188 | 0.256 |
| 200 × 50 | 25 | 3.29e-10 | 0.835 | 0.3266 | 0.3368 | 0.0136 | 0.0142 | 0 | 0.198 | 0.253 |
| 400 × 100 | 100 | 1.41e-13 | 2.509 | 2.3999 | 2.3257 | 0.0075 | 0.0563 | 0 | 1.532 | 1.686 |
| 400 × 100 | 50 | 4.71e-10 | 0.714 | 2.3469 | 2.3852 | 0.1162 | 0.0572 | 0 | 1.779 | 1.769 |
| 800 × 200 | 200 | 2.05e-13 | 2.528 | 16.5778 | 16.7356 | 0.0266 | 0.2077 | 0 | 7.367 | 8.180 |
| 800 × 200 | 100 | 1.17e-10 | 0.771 | 17.2425 | 17.4006 | 0.7881 | 0.2144 | 0 | 8.035 | 8.970 |
| 100 × 100 | 100 | 4.55e-14 | 2.621 | 0.5854 | 0.5869 | 0.0072 | 0.0115 | 0 | 0.276 | 0.341 |
| 100 × 100 | 50 | 5.51e-11 | 0.936 | 0.5850 | 0.5735 | 0.1097 | 0.0131 | 0 | 0.261 | 0.348 |
| 200 × 200 | 200 | 6.86e-14 | 2.961 | 3.8151 | 3.8127 | 0.0266 | 0.0426 | 0 | 1.197 | 1.429 |
| 200 × 200 | 100 | 4.01e-10 | 0.695 | 3.8066 | 3.8201 | 0.8007 | 0.0495 | 0 | 1.202 | 1.875 |
| 400 × 400 | 400 | 6.24e-14 | 2.889 | 26.1317 | 26.2716 | 0.1126 | 0.1651 | 0 | 6.108 | 6.827 |
| 400 × 400 | 200 | 2.83e-10 | 0.602 | 26.3616 | 26.5231 | 5.9976 | 0.1856 | 0 | 5.939 | 8.566 |
| 50 × 200 | 50 | 1.02e-13 | 0.419 | 0.5562 | 0.5419 | 0.3008 | 0.0148 | 0 | 0.261 | 0.429 |
| 50 × 200 | 25 | 7.43e-11 | 0.384 | 0.5314 | 0.5354 | 0.0903 | 0.0092 | 0 | 0.203 | 0.304 |
| 100 × 400 | 100 | 1.52e-13 | 0.389 | 3.3629 | 3.2781 | 2.2223 | 0.0563 | 0 | 1.805 | 3.332 |
| 100 × 400 | 50 | 3.09e-10 | 0.510 | 3.3430 | 3.3213 | 0.6619 | 0.0348 | 0 | 1.751 | 2.485 |
| 200 × 800 | 200 | 1.40e-13 | 0.311 | 21.9052 | 22.0046 | 17.5811 | 0.2140 | 0 | 7.696 | 13.948 |
| 200 × 800 | 100 | 1.74e-10 | 0.276 | 23.0258 | 23.0210 | 5.2162 | 0.1346 | 0 | 7.635 | 10.459 |

- **Rank.** On all eighteen the rank chosen from R at the relative threshold equals the
  designed rank.
- **Accuracy.** Solutions agree with the reference to about 2 × 10⁻¹³ where the condition
  number is 10³ and to about 5 × 10⁻¹⁰ where it is 10⁶, the size conditioning allows. Both
  sampled nullspace vectors of every rank-deficient or wide matrix have unit norm to
  7 × 10⁻¹⁶ and are mapped by A to at most 6 × 10⁻¹⁷.
- **The diagonal of R against the singular values.** The last kept diagonal entry was
  between 0.28 and 2.96 times the singular value of the same index: an indicator within
  a small factor on these matrices, on either side, and not that singular value.
- **Speed.** `factor!` takes between 1.3 and 4.4 times as long as the reference pivoted QR. Measuring
  the remaining column norms inside the pass that updates each column, and storing the
  kept rows of R row by row for completion, made `factor!` 1.7 to 2.7 times faster than
  the first working version, with the same tests and controls passing before and after.
  `set-rank!` costs about r²(n − r) multiply-adds and is the larger part of the work for
  a wide full-row-rank matrix.
- **Policy.** Safety 0 and safety 1 times are within 4% of each other, and the two
  policies produce identical output files.

## Limits

- The rank is the caller's decision. The diagonal of R can overstate the smallest singular
  value by a large factor, as the Kahan matrix shows; column pivoting is not a guarantee.
- A selected rank below the exact one solves a nearby problem; the residual reported is
  that problem's.
- Dense storage and unblocked, sequential arithmetic: small and moderate sizes.
- The result depends on column scaling. No equilibration is applied.
- One host and finite rosters establish neither universal accuracy nor portability beyond
  the supported runtime.
