# Batch 10: iterative refinement, and a solve on the Cholesky-form factor

State: completed; verified at safety 1 and safety 0 on 2026-11-20, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contracts are in
[refinement.md](refinement.md) and [cholesky-form.md](cholesky-form.md).

## What was delivered, and why

- `coalton-linalg/refinement`: `refine!`, which improves a solution through any
  approximate solve and an operator the caller applies, with `make-refiner` for its
  workspace and accessors for the step count, the best iterate's index and the residual
  history.
- `solve-cholesky!` in `coalton-linalg/sparse/spd-ldlt`: a solve with the Cholesky-form
  factor, of the shifted matrix that factor is of.

A consumer of this library was refining by hand, three fixed steps against an
unregularized operator, and asked for the general form. Its measurements then showed why
the general form matters more than the convenience. With the elimination order no longer
its cost, the numeric factorization is about 97% of its steady step, and it factors twice
at each iteration: once at a shift, in Cholesky's arithmetic, to prove a matrix positive
definite, and once unshifted to solve with it. The two matrices differ by the shift.
**Solving with the proof factor and refining against the unshifted matrix needs one
factorization where there were two.** That needed a solve on the proof factor, which by
design had none, and the refinement call.

The shape of the call was fixed with that consumer before any code:

- The operator is "apply A", and the library forms the residual, so that every caller
  computes it the same way.
- The solve is a function too, so one call serves the solving factor, the signed factor,
  the proof factor and a dense one.
- The norm is the infinity norm against an absolute tolerance, and the history is in the
  same norm.
- The stall factor is the caller's, beside the step limit, and a stall is its own
  outcome so a caller can count how often it had to fall back.
- The vector returned is the best iterate seen, with its index, not the last.

## A defect the tests found in the first version

The residual's norm was a running maximum written so that a NaN "is kept": an entry
replaces the maximum unless it is at most the maximum, and a NaN is not at most anything.
That keeps a NaN only until the next finite entry, which is not at most the NaN either
and replaces it. A right-hand side with one NaN in the middle was therefore refined as if
it were finite. The test for a non-finite right-hand side caught it on the first run.
The residual now notes an entry that is not finite for itself.

The same comparison stands in two places in the orthogonal decomposition. There a lost
NaN is recovered, because the sum of squares that follows contains it, except for a
segment in which every entry after the NaN is exactly zero. No call through the public
interface was found that produces such a segment, since a reflector applied with an
infinite scalar leaves no finite entry behind it. It is corrected in a separate change
with a test of the helper itself.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`, with every earlier roster unchanged. For this batch:

- **The iteration, bit for bit.** Sixty systems: positive definite matrices of orders 1
  to 8, solved through a factor of the matrix plus 0, 1/64, 1/4 and 2 times the identity.
  Each call converges, in at most 14 steps. The iteration is recomputed independently in
  the tests with the library's own product and solve, and every recorded norm and the
  returned vector agree with it bit for bit. The norms fall at every step. The returned
  vector agrees with the exact rational solution to 2⁻³⁰ of its size.
- **Outcomes.** A correction of zero stalls and returns the starting vector untouched; a
  good step followed by a bad one returns the first iterate, at the limit with the stall
  test off and as a stall with it on; a solve that brings the norm down by a quarter each
  step passes a stall factor of 1 and stalls at 1/2, keeping the one step that improved;
  a limit of zero judges the starting vector as it stands; a solve that refuses on its
  second call leaves the first iterate; a NaN in the right-hand side and an infinite
  correction are outcomes that return the starting vector; wrong lengths, a negative or
  NaN tolerance and a zero or NaN stall factor are refused before anything is read.
- **A value the operator never reads.** For `[[1, 0], [0, 0]]` and a correction whose
  second entry is infinite, the residual is finite and zero. The call reports
  `RefineNonFinite` and returns the starting vector.
- **Both uses.** On sixteen wide matrices: `J Jᵀ y = b` from the Cholesky-form factor of
  `J Jᵀ − I/8`, refined against J (Jᵀ y) applied as two products, converges in at most 8
  steps to the exact rational solution; and the minimum-norm x of `J x = b` from the
  signed factor of `[[I, Jᵀ], [J, −2⁻⁶ I]]`, refined against `[[I, Jᵀ], [J, 0]]`, converges
  in at most 5 steps to the exact minimum-norm solution.
- **The proof solve.** On the 107 matrices of the Cholesky-form roster, each row of the
  exact residual of `(A − shift I) x = b` is within `(3w + 5) u |G| |Gᵀ| |x|`; the largest
  was 0.338 of that. It refuses a solving factor, wrong lengths and a non-finite
  right-hand side with the output untouched, agrees with itself in place bit for bit,
  and solves with the shift it was factored at.
- **Allocation.** A compiled caller refining from the same start each round, and one
  making the proof factorization and two proof solves: zero bytes over 10,000 rounds
  after 10,000 warm-up rounds, at both policies.

**Planted source controls.** Each of these edits was applied alone and the focused tests
run. Every planted source compiled, every one failed the tests, and the restored source
passed.

| Source | Planted defect | Caught |
| --- | --- | --- |
| refinement | The residual is A x − b | yes |
| refinement | An entry that is not finite is not noted | yes |
| refinement | The best iterate is not saved | yes |
| refinement | The best iterate's step is not recorded | yes |
| refinement | A step that leaves the norm equal is not a stall | yes |
| refinement | The stall factor is taken as 1 | yes |
| refinement | A stalled call returns its last iterate | yes |
| refinement | A call at its limit returns its last iterate | yes |
| refinement | A met tolerance does not stop the steps | yes |
| refinement | A refusing solve is taken as a solve | yes |
| refinement | An update that is not finite goes unnoticed where the operator does not read it | yes |
| refinement | A tolerance that is not a number is accepted | yes |
| refinement | A stall factor of zero is accepted | yes |
| refinement | The vector lengths are not checked | yes |
| refinement | The step count is not recorded | yes |
| refinement | The residual count is not recorded | yes |
| spd-ldlt | The proof solve's forward pass does not divide by the diagonal | yes |
| spd-ldlt | The proof solve's backward pass divides by the pivot, not its root | yes |
| spd-ldlt | The proof solve's backward pass leaves out the column | yes |
| spd-ldlt | The proof solve accepts a solving factor | yes |
| spd-ldlt | The proof solve accepts a right-hand side that is not finite | yes |
| spd-ldlt | The proof solve does not undo the permutation | yes |

## Limits

- The residual is formed in working precision. Refinement removes the error of the
  approximate solve and converges to the conditioning of the system, not past it.
- A shifted factor converges only when the shift is small beside the smallest eigenvalue:
  the error is multiplied by at most `σ / (λ_min − σ)` a step, which exceeds 1 when the
  smallest eigenvalue is below twice the shift. A proof that passes with little margin
  is such a case. The stall outcome is how a caller learns of it.
- How many factorizations this saves a consumer is the consumer's to measure, and is not
  measured here.
- One host and finite rosters establish neither universal behaviour nor portability beyond
  the supported runtime.
