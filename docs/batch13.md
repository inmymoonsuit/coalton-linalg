# Batch 13: conjugate gradients

State: completed; verified at safety 1 and safety 0 on 2026-11-20, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract is in
[conjugate-gradient.md](conjugate-gradient.md).

## What was delivered, and why

A new package, `coalton-linalg/conjugate-gradient`:

- `cg!`: conjugate gradients for a symmetric positive definite operator the caller applies,
  to an absolute tolerance on the 2-norm of the residual.
- `pcg!`: the same with a preconditioner the caller supplies as a function.
- `make-cg` for the workspace, and accessors for the step count, the residual of the
  vector returned, the running residual's history and the number of times the running
  residual was replaced.

Until now every solve in this library went through a factorization. A system whose
matrix is a product never formed, or is too large to factor, or is known only as a
routine that multiplies a vector, had no entry. This is the first of the iterative
methods on the roadmap and the one the others are built from. It is built ahead of a
consumer's measurement.

**The shape follows `refine!`**, which a consumer fixed before that was written: the
operator and the preconditioner are functions passed as arguments, the workspace is only
storage, outcomes are values, and the history is read by accessors.

## What was decided, and why

**Convergence is not reported on the method's own residual.** The running residual is a
recurrence, and a recurrence drifts from what it stands for. A call that says "converged"
is saying something about the vector it returns, so the residual of that vector is
formed from the operator before the word is used, and if it is not within the tolerance
the method starts again from it. One product is the price. The tests corrupt a single
product by a thousandth: the running residual then meets a tolerance of 2⁻⁴⁰ of the
right-hand side while the residual of the vector does not, and the call goes on, taking
more steps than the honest run, until the latter is met.

**The vector returned is the last iterate.** `refine!` returns the best iterate seen,
measured by residual. Here the residual is not what improves monotonically; the error in
the operator's norm is, and by that measure the last iterate is the best. The tests
measure it.

**The norm is the 2-norm**, the one the method's own quantities give for nothing, and it
is formed as a plain sum of squares. A system with entries near 10¹⁵⁴ overflows it and is
reported as not finite; the contract says to scale such a system.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`, with every earlier roster unchanged. For this batch:

- **Ninety-six solves, each held three ways.** Positive definite matrices of orders 1 to
  12, each solved from zero, from another starting vector, preconditioned by the
  diagonal, and preconditioned by a factor of the matrix itself, to a tolerance of 2⁻⁴⁰
  of the right-hand side's norm.
  1. *The residual, exactly.* `b − A x` for the vector returned is formed in rational
     arithmetic; the norm reported agrees with it to the rounding of one product and is
     within the tolerance.
  2. *What the method is.* The error of every iterate is measured exactly in the norm of
     the operator, `(x − x*)ᵀ A (x − x*)` with x* the exact rational solution, and never
     rises from one iterate to the next by more than a part in 2³⁰.
  3. *The iteration, bit for bit.* It is recomputed independently in the tests with the
     library's own product, and every running norm and the vector returned agree with it
     bit for bit.

  The plain method took at most 13 steps, at order 12, and never more than the order
  plus three; preconditioned by a factor of the matrix it took at most two.
- **Outcomes.** Sizes and tolerances are refused before anything is read. A right-hand
  side or a start that is not finite, a product that overflows, a curvature that is not
  a number, a step length that overflows and a running residual whose norm overflows
  while everything before it is finite are each `CgNonFinite`, the step length before the
  vector is touched. An operator with a direction of zero curvature, one of negative
  curvature and the zero operator are `CgNotPositive` with the start returned. A
  preconditioner that is negative, one that is zero, and one that is positive on its
  first call and negative on its second are `CgNotPositive`; one that refuses at once and
  one that refuses on its second call are `CgPreconditionerRefused`. A limit of zero
  judges the start and leaves it untouched; a limit met reports the residual of what is
  returned; convergence on the last step allowed is convergence. A tolerance of zero is
  met where the arithmetic is exact. A refused call withdraws the earlier report.
- **One wrong product.** As above; and with the call stopped at its limit two steps after
  the wrong product, the report agrees with the exact residual of the vector to a part
  in 10⁹ while the running residual differs from it by more than a part in 10⁶.
- **Lifetime.** One workspace through plain and preconditioned calls gives the same vector
  for the same call, bit for bit.
- **Encapsulation.** A compiled external caller reaches the workspace through its public
  operations; three probes that reach into it or pass its status for another's are
  rejected by the compiler.
- **Allocation.** Compiled callers of `cg!` and of `pcg!` with a factor's solve as
  preconditioner: zero bytes over 10,000 rounds after 10,000 warm-up rounds, at both
  policies.

**Planted source controls.** Each of these edits was applied alone to the final source
and the focused tests run. Every planted source compiled, every one failed the tests, and
the restored source passed.

| Planted defect | Caught |
| --- | --- |
| The running residual's word is taken for convergence | yes |
| The residual is A x − b | yes |
| The step is taken backwards | yes |
| The running residual is not updated by the step | yes |
| The new direction's weight is inverted | yes |
| The new direction forgets the old one | yes |
| A direction without positive curvature is used | yes |
| A first preconditioned residual without a positive product is used | yes |
| A later preconditioned residual without a positive product is used | yes |
| A first refusal of the preconditioner is ignored | yes |
| A later refusal of the preconditioner is ignored | yes |
| A first residual that is not finite is accepted | yes |
| A curvature that is not a number is called not positive | yes |
| A step length that is not finite is applied | yes |
| A running residual that is not finite is accepted | yes |
| The vector lengths are not checked | yes |
| A tolerance that is not a number is accepted | yes |
| A negative tolerance is accepted | yes |
| A refused call leaves the earlier report standing | yes |
| A refused call leaves the earlier step count standing | yes |
| The step count is not recorded | yes |
| The residual replacements are not counted | yes |
| A running norm is not recorded | yes |
| At the limit the running residual is reported | yes |
| Convergence at the last step allowed is called the limit | yes |
| A converged call reports no residual | yes |
| The plain method's preconditioner is not the identity | yes |
| The first direction is not the preconditioned residual | yes |

One control of the first run was not one: it doubled the plain method's preconditioned
residual, which changes nothing, not even a bit, since every quantity that follows is
scaled by a power of two and the step is their ratio. It was replaced.

## Limits

- **The matrices are small and well conditioned**, so that exact solutions and exact
  errors can be computed. How the method behaves on a large or badly conditioned system,
  how many steps it takes and where it stalls, is not measured here; the contract promises
  no rate.
- **A residual within tolerance bounds the error only through the inverse's norm**, which
  this package does not estimate.
- **The operator's symmetry is assumed**, not checked, and `CgNotPositive` is a finding
  about one direction and not a test of the operator.
- The bit-for-bit recomputation shares the method's formulas with the library by
  necessity: it shows that the library computes those formulas and nothing else. That the
  formulas are the right ones is what the exact residuals and the exact errors show.
- One host and finite rosters establish neither universal behaviour nor portability beyond
  the supported runtime.
