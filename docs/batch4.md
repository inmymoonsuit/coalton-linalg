# Batch 4: dense symmetric indefinite LDLᵀ and inertia

State: completed; verified at safety 1 and safety 0 on 2026-11-19, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract and the derivation are in
[indefinite-ldlt.md](indefinite-ldlt.md).

## What was delivered

A new package, `coalton-linalg/dense/indefinite-ldlt`: `make-factor`, `factor!`, `solve!`,
inertia counts and inspection views for `P (A - shift I) Pᵀ = L D Lᵀ` of a dense symmetric
matrix, with 1-by-1 and 2-by-2 pivots chosen by a column search that bounds every
multiplier by 1/(1 - α) ≈ 2.78. The existing sparse SPD contract is unchanged.

- The inertia counts both eigenvalues of each 2-by-2 block.
- A singular matrix completes as `FactorSuccess` with a nonzero `zero-count`; `solve!`
  then returns `SolveSingular`.
- A 2-by-2 block is inverted in a form divided by its off-diagonal entry, so no square of
  an entry is formed.
- A compiled example, [minimum-norm.ct](../examples/minimum-norm.ct), solves a
  minimum-norm problem through the augmented system and reads its inertia first.

The prose derivation in the contract page was written before the implementation.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`, with every earlier roster unchanged. For this batch:

- **Exact references, 365 factorizations** holding 1,234 1-by-1 and 368 2-by-2 pivots:
  13 named cases, 320 seeded integer matrices of order 1 through 10 (half with a zero
  diagonal, each at shift 0 and shift 1/2), and 32 congruences `S diag(±1) Sᵀ`. For each:
  - the exported factor is multiplied out in exact rational arithmetic and compared with
    the exact input; the largest error was 0.96 units of roundoff times
    (max|A| + max(|L||D||L|ᵀ)), against an asserted bound of 4(n + 1);
  - the permutation, the block layout, the unit diagonal and zero upper triangle of L,
    the 2-by-2 acceptance conditions and the multiplier bound are checked directly;
  - the reported inertia equals the inertia of the exported block diagonal, computed
    from the exact trace and determinant of each 2-by-2 block;
  - where the exact matrix is nonsingular (345 of 365), the reported inertia equals the
    exact inertia from an independent exact symmetric elimination with its own pivot
    rule, and a solve's exact residual is within 8(n + 1) units; the largest was 0.56.
- **Inertia cases.** Three 2-by-2 matrices whose inertia differs from the count of blocks
  or of diagonal signs; three exactly singular matrices with exact zero counts; the
  second-difference matrix at seven shifts against its known spectrum; 24 augmented
  matrices of full-row-rank J with inertia (n, m, 0), whose solve returns the minimum-norm
  solution to 2⁻⁴⁰ against the exact one; one augmented matrix with dependent rows, inertia
  (3, 1, 1), whose solve is refused as singular with its output untouched.
- **Identities.** Scaling by 2⁻⁶⁰⁰, 2⁻³¹, 2²⁹ and 2⁶⁰⁰ leaves the permutation, the blocks
  and L unchanged bit for bit and scales D and the solution exactly. A shift applied by
  `factor!` equals the same matrix stored with the shift subtracted. NaN in the strict
  upper triangle changes nothing. A repeated call repeats every bit.
- **Lifetime and refusals.** No factor before the first call; every solve refusal leaves
  the output and the factor alone; five kinds of refused factorization each clear a
  standing factor and its inertia; overflow during elimination, an overflowing shifted
  diagonal and an overflowing solution are refusals; the dimension bound is a typed error.
- **Encapsulation.** Five external compile probes that must be rejected, read for their
  reasons, and one public control that must compile and run.
- **Allocation.** A compiled Coalton caller that refactors a 12-by-12 matrix with 2-by-2
  pivots, reads the inertia and solves: zero bytes over 10,000 calls after 10,000 warm-up
  calls, at both policies.

**Planted source controls.** Each of these edits to `indefinite-ldlt.ct` was applied alone
and the focused tests run; every one failed, and the restored source passed.

| Planted defect | Caught |
| --- | --- |
| A mutual maximum is taken as a 1-by-1 pivot | yes |
| A 2-by-2 block counts as two positive eigenvalues | yes |
| The exchange skips the entries between its two indices | yes |
| The 2-by-2 block is inverted through its unscaled determinant | yes, after a test change |
| The shift is ignored | yes |
| A zero pivot is not counted | yes |
| Solve does not undo the permutation | yes |
| The pivot threshold is a thousandth | yes |
| A refusal keeps the inertia counted so far | yes |
| The strict upper triangle of the input is read | yes |
| The search stops at the first column whatever its diagonal | yes |
| The backward substitution skips a column's first contribution | yes |
| Solve inverts the 2-by-2 block through its unscaled determinant | yes |

The fourth control was missed on its first run: the scaling test used 2^±400, where the
square of an entry is still representable. The test now uses 2^±600 and catches it. Writing
the controls for a kept partial inertia and for the threshold also added an assertion each
that the first draft of the tests lacked.

## Measurements

Observations, not targets. Twenty seeded matrices generated outside the library: `definite`
is BBᵀ/n + I; `random` has standard normal entries, symmetrized; `zerodiag` is the same with
a zero diagonal; `augmented` is `[[I, Jᵀ], [J, 0]]` with a standard normal J of n/3 rows.
The library's factors were written to files and checked by a separate NumPy program.

- Multiplier: the largest magnitude in L below the diagonal.
- Reconstruction: max|PAPᵀ - LDLᵀ| in units of roundoff (2⁻⁵³) times
  (max|A| + max(|L||D||L|ᵀ)).
- Residual: max|Ax - b| in units of roundoff times (‖A‖∞‖x‖∞ + ‖b‖∞).
- Times are the best of 20 (5 above n = 300) after one warm-up call, in milliseconds.
  Bytes are over 20 warmed factor-and-solve pairs. One process at a time under a 1024 MB
  heap cap.
- Reference: `scipy.linalg.ldl` of SciPy 1.18.1, best of 5, called from Python. It runs a
  blocked routine over an optimized, possibly multithreaded BLAS and then assembles its
  own output arrays. It is a reference point for scale, not a like-for-like comparison.

| Family | n | 2-by-2 pivots | Multiplier | Reconstruction | Residual | `factor!`, s1 | `factor!`, s0 | `solve!`, s1 | Bytes | Reference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| definite | 50 | 0 | 0.246 | 1.57 | 0.70 | 0.0226 | 0.0212 | 0.0018 | 0 | 0.060 |
| definite | 100 | 0 | 0.160 | 3.51 | 0.85 | 0.1564 | 0.1557 | 0.0065 | 0 | 0.107 |
| definite | 200 | 0 | 0.148 | 5.37 | 1.30 | 1.1826 | 1.1636 | 0.0247 | 0 | 0.311 |
| definite | 400 | 0 | 0.103 | 8.28 | 1.58 | 8.7358 | 8.6083 | 0.0930 | 0 | 1.882 |
| definite | 800 | 0 | 0.084 | 13.60 | 1.63 | 75.7755 | 72.5139 | 0.4208 | 0 | 10.995 |
| random | 50 | 8 | 1.556 | 0.53 | 0.82 | 0.0245 | 0.0218 | 0.0018 | 0 | 0.191 |
| random | 100 | 21 | 1.559 | 0.36 | 1.22 | 0.1581 | 0.1450 | 0.0065 | 0 | 0.334 |
| random | 200 | 53 | 1.641 | 0.32 | 1.96 | 1.1413 | 0.9889 | 0.0251 | 0 | 0.882 |
| random | 400 | 96 | 1.634 | 0.26 | 2.77 | 8.2976 | 7.5141 | 0.0906 | 0 | 3.431 |
| random | 800 | 230 | 1.552 | 0.24 | 6.44 | 69.4895 | 59.5075 | 0.4212 | 0 | 14.055 |
| zerodiag | 50 | 15 | 1.551 | 0.44 | 0.93 | 0.0238 | 0.0195 | 0.0017 | 0 | 0.152 |
| zerodiag | 100 | 28 | 1.556 | 0.31 | 2.03 | 0.1513 | 0.1353 | 0.0065 | 0 | 0.339 |
| zerodiag | 200 | 61 | 1.716 | 0.35 | 1.82 | 1.0996 | 0.9453 | 0.0248 | 0 | 0.778 |
| zerodiag | 400 | 112 | 1.574 | 0.32 | 3.13 | 8.2669 | 7.3864 | 0.0913 | 0 | 3.220 |
| zerodiag | 800 | 235 | 1.700 | 0.27 | 5.64 | 67.3888 | 62.4792 | 0.3569 | 0 | 12.855 |
| augmented | 50 | 5 | 1.521 | 0.34 | 0.70 | 0.0235 | 0.0204 | 0.0019 | 0 | 0.107 |
| augmented | 100 | 11 | 1.532 | 0.83 | 1.25 | 0.1532 | 0.1402 | 0.0066 | 0 | 0.180 |
| augmented | 200 | 13 | 1.551 | 1.40 | 1.39 | 1.1290 | 1.0722 | 0.0237 | 0 | 0.419 |
| augmented | 400 | 25 | 1.555 | 2.11 | 1.67 | 8.3305 | 7.8493 | 0.0916 | 0 | 2.211 |
| augmented | 800 | 36 | 1.560 | 2.22 | 1.72 | 70.6815 | 68.5071 | 0.3678 | 0 | 11.187 |

- **Inertia.** On all twenty the reported inertia equals the count of eigenvalue signs
  from NumPy's symmetric eigenvalue routine; the augmented matrices give (n - m, m, 0)
  for their n - m unknowns and m constraints.
- **Accuracy.** Reconstruction stays within 14 units and residuals within 7. The
  reconstruction figure grows with n on the definite family, to 13.6 at n = 800.
- **Speed.** At n = 800 a factorization's n³/3 multiply-adds run at about 2.4 thousand
  million a second. It is faster than the reference call at n = 50, and at n = 100 except
  on the definite family; it is five to seven times slower at n = 800, where the
  reference's blocked arithmetic pays off.
- **Policy.** Safety 0 is between under 1% and 18% faster than safety 1 on `factor!`. The
  two policies' factor and solution files are identical byte for byte.

## After review

A reviewer of the first delivery found no defect and made three points, each taken up:

- The dimension bound alone refused nothing a host would not already fail to allocate.
  `make-factor` now refuses n² above `default-entry-limit` (16,777,216, n = 4096) by name,
  and `make-factor-with-limit` takes the budget from the caller, as sparse analysis does.
  **A caller of the first delivery that made a factor above n = 4096 now needs the explicit
  budget.**
- A comment and the contract gave the scaled determinant's range for same-sign diagonals
  only. Both now give the range that holds for either sign. No code changed.
- An artifact driver, `scripts/check.ps1 -Task Indefinite`, writes each system's factor,
  inertia and solution at named shifts for an independent checker; the contract page has
  the file protocol. Its regressions cover an exact factor file, shifts from the meta
  file, a named singular solve, a named factor refusal, and removal of stale artifacts.

## Limits

- Dense storage and cubic work: this is for small and moderate n. At 70 ms for n = 800 and
  cubic growth, n = 3000 would take about four seconds; that figure is extrapolated, not
  measured. It does not remove a large sparse augmented-system bottleneck.
- The inertia is that of the computed D. Near a threshold the caller's shift and the
  measured reconstruction error decide what a count means; no rank decision is made here.
- The column search bounds multipliers; no bound on entry growth is claimed for it, and
  none was approached on these rosters.
- A singular factor does not solve consistent singular systems.
- One host and finite rosters establish neither universal accuracy nor portability beyond
  the supported runtime.
