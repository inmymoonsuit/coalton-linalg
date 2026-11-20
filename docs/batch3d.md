# Batch 3d: condition estimate from a prepared factor

State: completed; verified at safety 1 and safety 0 on 2026-11-19, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract is [condition.md](condition.md).

## What was delivered

Three operations in `coalton-linalg/sparse/spd-ldlt`, for the matrix `A - shift I` that
`factor!` has just factored:

- `estimate-inverse-norm1!` gives a lower estimate of the 1-norm of its inverse in at most
  11 solves with the factor, and stores it on the factor until the next `factor!`.
- `shifted-norm1` gives its 1-norm from the lower values the caller factored.
- `reciprocal-condition1` combines the two, with fixed conventions for zero and for a
  missing estimate.

Outcomes are a new `ConditionStatus`. Nothing is allocated on success; two vectors of
length n were added to the factor's owned workspace at analysis.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`, with every earlier roster unchanged. For this batch:

- **Exact references.** Each test matrix is inverted in exact rational arithmetic on its
  stored F64 values by a dense elimination written for the tests. The estimate must not
  exceed the exact norm beyond a stated rounding allowance and must reach a pinned
  fraction of it; `shifted-norm1` must equal the exact norm.
  - Four diagonal matrices: the estimate is exactly the norm.
  - `[[1, 1-e], [1-e, 1]]` for e = 2⁻¹⁰, 2⁻²⁰, 2⁻³⁰: the ascent alone returns about 1/2
    against 1/e; the estimate is the norm to six digits.
  - A singular integer Gram matrix plus e times the identity, same three e: the norm to
    six digits.
  - A tridiagonal matrix, a shifted 3-by-3, and 48 seeded integer matrices BᵀB + I of
    order 2 through 7: the estimate is between 0.571 and 1.000 of the norm.
- **Identities.** Scaling a matrix by 2⁻²⁰⁰, 2⁻³, 2⁷ and 2²⁰⁰ scales the estimate by exactly
  the inverse power. A shift applied by `factor!` gives the same estimate, bit for bit, as
  the same matrix stored with the shift already subtracted.
- **Lifetime and refusals.** No estimate before a factor, after a refused factor, or after
  a refactor; a repeated estimate is unchanged; a 1-by-1 pivot of 2⁻¹⁰⁶⁰ factors and then
  refuses the estimate as non-finite after one counted solve; the empty matrix succeeds
  with 0 and no solve. `shifted-norm1` refuses a wrong count, a non-finite shift or value,
  and an overflowing column sum, each by its typed error through a compiled caller.
- **Allocation.** A compiled Coalton caller that refactors, takes the norm, estimates and
  forms the reciprocal: zero bytes over 10,000 calls after 10,000 warm-up calls, at both
  policies.

**Planted source controls.** Each of these edits to `spd-ldlt.ct` was applied alone and the
focused tests run; every one failed, and the restored source passed.

| Planted defect | Caught |
| --- | --- |
| The alternating probe never raises the estimate | yes |
| The ascent stops after its first step | yes |
| `factor!` leaves the previous estimate standing | yes |
| The alternating bound divides by n, not by the probe's own norm | yes |
| `shifted-norm1` omits the mirrored entry | yes |
| `shifted-norm1` ignores the shift | yes |
| The estimate does not check that a factor stands | yes |

## Measurements

Observations, not targets. The thirteen sparse systems of the [Batch 3c record](batch3c.md),
each at shift 0 and at the positive shift its supplier marked as still factorable, about half
its smallest eigenvalue: 26 factors. The reference is the dense shifted matrix inverted in
double precision by NumPy, read for its largest absolute column sum; it was used as a
comparison only, after the implementation was written. Times are the best of 20 after one
warm-up call, in milliseconds, at shift 0. Bytes are over 1000 warmed estimates. One process
at a time under a 1024 MB heap cap.

| System | n | Condition, 1-norm | Estimate / reference | Solves | Estimate, s1 | Estimate, s0 | One solve, s1 | `factor!`, s1 | Bytes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| chain25 | 400 | 5.3e4 | 1.0000000000 | 5 | 0.0179 | 0.0147 | 0.0030 | 0.0052 | 0 |
| chain100 | 1,600 | 8.2e5 | 1.0000000000 | 5 | 0.0683 | 0.0583 | 0.0115 | 0.0213 | 0 |
| chain200 | 3,200 | 3.2e6 | 1.0000000000 | 5 | 0.1611 | 0.1185 | 0.0245 | 0.0460 | 0 |
| grid15x15 | 1,290 | 5.3e4 | 1.0000000000 | 7 | 0.1333 | 0.0954 | 0.0175 | 0.1091 | 0 |
| grid21x21 | 2,562 | 1.5e5 | 1.0000000000 | 7 | 0.3118 | 0.2065 | 0.0412 | 0.3194 | 0 |
| wheel100 | 402 | 2.0e6 | 1.0000000000 | 5 | 0.0205 | 0.0163 | 0.0036 | 0.0095 | 0 |
| wheel400 | 1,602 | 1.2e8 | 1.0000000000 | 5 | 0.0840 | 0.0698 | 0.0145 | 0.0398 | 0 |
| loosechain25 | 352 | 2.1e5 | 1.0000000000 | 5 | 0.0173 | 0.0134 | 0.0029 | 0.0060 | 0 |
| loosechain100 | 1,402 | 2.1e5 | 1.0000000000 | 5 | 0.0685 | 0.0536 | 0.0118 | 0.0247 | 0 |
| medium | 602 | 1.0e5 | 1.0000000000 | 5 | 0.0246 | 0.0194 | 0.0041 | 0.0072 | 0 |
| large | 2,052 | 4.1e5 | 1.0000000000 | 5 | 0.0996 | 0.0797 | 0.0170 | 0.0348 | 0 |
| loosegrid15x15 | 1,277 | 4.2e7 | 1.0000000000 | 5 | 0.3714 | 0.2418 | 0.0741 | 1.3250 | 0 |
| loosewheel100 | 354 | 3.8e6 | 1.0000000000 | 5 | 0.3043 | 0.1931 | 0.0602 | 3.0156 | 0 |

- **Agreement.** On all 26 factors the estimate agrees with the reference to ten digits;
  the smallest ratio is 0.9999999999, on chain200 at its positive shift. The shifted
  factors have about twice the condition of the unshifted ones and the same solve counts.
  `shifted-norm1` agrees with the reference norm to 2.2e-14 relative or better; the
  difference is summation order.
- **Cost.** An estimate takes 5 or 7 solves here and costs 5 to 8 times one measured solve.
  Against the factorization it follows, it ranges from a tenth of its time (loosewheel100)
  to 3.4 times it (chain25): where the factor is nearly as sparse as the matrix, a
  factorization costs about two solves.
- These systems all reach the norm. The small seeded roster above shows that is not
  general: there the estimate fell to 0.571 of the norm.

## Limits

- The estimate can fall short of the norm, so the reciprocal condition is an optimistic
  figure. No bound on the shortfall is claimed.
- It is the 1-norm of the factored, shifted matrix. A singular-value ratio, a rank decision
  and the conditioning of a matrix the caller formed A from are different quantities;
  [the contract](condition.md) relates them where a relation exists.
- One host and finite rosters establish neither universal accuracy nor portability beyond
  the supported runtime.
