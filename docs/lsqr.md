# Least squares on an operator

Package: `coalton-linalg/lsqr`. Source: [lsqr.ct](../src/lsqr.ct).
State: completed as Batch 20; verified at safety 1 and safety 0. The [batch
record](batch20.md) holds the controls and [the design](design-iterative-least-squares.md)
the derivation.

## What it is for

A matrix A with m rows and n columns that is available only as two functions, one that
writes `A v` and one that writes `Aᵀ u`: a matrix too large to factor, a product never
formed, a factor that fills. The call finds

    the x that minimizes   ‖b − A x‖₂² + λ² ‖x‖₂²

for a damping λ ≥ 0 the caller gives, and with λ = 0, where the minimizer is not unique,
the one of least norm. A may be of any shape and any rank.
[Conjugate gradients](conjugate-gradient.md) answer the same question through
`Aᵀ A + λ² I`, whose condition number is the square of the problem's; this method works on
A and Aᵀ and never forms that product. The factorizations of this library are the tool
where they fit; this is the tool where they do not.

## Interface

- `make-lsqr(rows, columns, step-limit) -> LsqrWorkspace` allocates the workspace once:
  three vectors of each length and two numbers a step.
- `lsqr!(workspace, apply, apply-transpose, b, x, damping, residual-tolerance,
  gradient-tolerance) -> LsqrStatus`.
  - `apply` is a function of two arrays, `(v, out)`, that **writes A v into out**, v of
    length `columns` and out of length `rows`. `apply-transpose` writes `Aᵀ u` into an
    array of length `columns`.
  - `b` is the right-hand side. It is read once, before x is written, so **b and x may be
    one array** where the matrix is square.
  - `x` is written only. **There is no starting vector**: the method starts from zero,
    which is what makes its answer the one of least norm.
  - `damping` is λ, not its square.
  - The tolerances are absolute, each on the **2-norm** of a quantity formed from the
    operators for the vector returned: the residual `b − A x`, and the gradient
    `Aᵀ(b − A x) − λ² x`. Either may be zero, which that quantity then meets only by being
    zero.
- `step-count(workspace)`: steps the latest call took. Each applies the operator once and
  its transpose once.
- `confirmation-count(workspace)`: how often the latest call formed the residual and the
  gradient from the operators, at one product with each.
- `residual-norm(workspace)`, `gradient-norm(workspace)`: the two norms for the vector
  returned, **formed from the operators**, after a call that ended `LsqrResidualMet`,
  `LsqrGradientMet`, `LsqrExhausted` or `LsqrLimit`; −1 otherwise.
- `recurrence-residual-norms(workspace)`, `recurrence-gradient-norms(workspace)`: the
  method's running estimates of the two after each step, entry 0 for the zero vector.
- `operator-norm-estimate(workspace)`: the Frobenius norm of the bidiagonal matrix the call
  reached; see [below](#what-is-and-is-not-claimed) for what it is worth.
- `rows(workspace)`, `columns(workspace)`, `step-limit(workspace)`.

| LsqrStatus | Meaning |
| --- | --- |
| LsqrResidualMet | `‖b − A x‖`, formed from the operator, is within the residual tolerance |
| LsqrGradientMet | The gradient, formed from the operators, is within its tolerance, and the residual is not within its own |
| LsqrExhausted | The method's sequences ended, so that in exact arithmetic x is the solution, and neither norm as formed is within its tolerance. A diagonal entry of the reduction that has underflowed to zero ends them as a vanished norm does |
| LsqrLimit | The step limit was reached with neither tolerance met; the two norms say how far it got |
| LsqrNonFinite | b, an operator or an update produced a value that was not finite |
| LsqrOutOfRange | The right-hand side's size divided by the size of the transpose's first image is not a number the format holds, or the damping is more than 2^500 times that size. Never where that image is zero: zero then solves the problem whatever the damping |
| LsqrBadSize | b or x does not have the workspace's length for it |
| LsqrBadTolerance | A tolerance is negative or not a number |
| LsqrBadDamping | The damping is negative, infinite or not a number |

The three `Bad` statuses are returned before anything is read or written, and
`LsqrOutOfRange` with x untouched. A call allocates nothing; the two functions are the
caller's and are made once, outside the loop that calls. An exception raised inside either
passes through. One workspace serves one call at a time.

**What x holds on the other outcomes.** After the first four, the last iterate. After
`LsqrNonFinite`: untouched where the right-hand side or the transpose's first image was
not finite; otherwise the iterate at which the value arose, which is the last one
completed where an operator's image was at fault.

## What is and is not claimed

- **The two norms reported are of the vector returned.** The method carries estimates of
  both forward by recurrences, which is what makes a step cost two products, and those
  estimates drift from the truth. When an estimate is within its tolerance the residual
  and the gradient are formed from the operators, at two more products, and the
  tolerances are tested on those. `LsqrResidualMet` and `LsqrGradientMet` therefore mean
  what they say of the vector returned, to the rounding of two products. At the limit and
  at the end of the sequences the report is likewise of what was formed.
- **A confirmation that fails does not restart the method.** The iteration goes on from
  where it was, and the estimate that triggered must then fall by the factor it was out
  by before it triggers again; one that was exactly zero is switched off. So a wrong
  estimate need not cost two products a step. In the tests one wrong product, after which
  the estimates are wrong for good, costs two confirmations in forty steps; a product
  wrong by less, after which the estimate is within the tolerance at 22 of forty steps,
  costs four; and an estimate that is zero at three steps running costs one.
  **A tolerance below what rounding lets the norm reach does cost them.** The estimate
  then goes on falling while the norm formed stays where it is, and each fall triggers
  again: on a four by three matrix whose gradient stops at 7 × 10⁻¹⁵, a tolerance of
  5 × 10⁻¹⁵ costs 35 confirmations in forty steps, and one of 10⁻¹⁵ costs 32. Nothing in
  the call recognizes a stall; the step limit is what bounds it.
- **Which tolerance to give.** A system that has a solution is asked for its residual. A
  least squares problem, whose residual does not go to zero, is given a residual
  tolerance of zero and stopped on its gradient, which is zero exactly at a minimizer;
  so is every damped problem. The gradient has the size of the operator times the
  right-hand side, and a tolerance for it is a fraction of that product.
- **A gradient within tolerance is not an error within tolerance.** The error can exceed
  the gradient by the square of the reciprocal of the smallest singular value that is not
  zero, or by 1/λ² with damping. Nothing here bounds it.
- **The minimum norm is a property in exact arithmetic.** Every iterate is a combination
  of images of the transpose, and among the minimizers exactly one is in the transpose's
  range: the shortest. In floating point each step adds a rounding's worth outside that
  range, which the operator cannot see and nothing removes. On the rank-deficient
  matrices of the tests that part is at most 3 × 10⁻¹⁵ of the vector.
- **Three quantities are monotone in exact arithmetic**: from one iterate to the next the
  damped objective `‖b − A x‖² + λ²‖x‖²` does not rise, the norm of x does not fall, and the
  distance from x to the minimizer does not rise. The tests measure each exactly on every
  iterate of 120 small problems and find none broken by more than a part in 2³⁰, a
  distance already at rounding level not being compared. **The residual without the damped
  term does not rise either**, since the objective does not rise and the norm does not
  fall, and the tests measure it with the other three. **The gradient is not monotone.**
  An earlier version of this page said the residual was not; an independent review
  showed that it is.
- **In exact arithmetic the method ends in at most as many steps as the rank of A.** In
  floating point the vectors it builds lose their orthogonality; it then takes more
  steps than that, and on a badly conditioned matrix many more. The count grows with the
  condition number of A, where conjugate gradients on a definite matrix grow with its
  root. **No rate is promised**; the call stops on what it measures, and the limit is the
  caller's.
- **`operator-norm-estimate` is an estimate and proves nothing, in either direction.** It
  is below the Frobenius norm of A only until a singular value has converged. The vectors
  then lose their orthogonality along it and it is counted again, and a dominant singular
  value converges within a few steps: on a matrix of order 40 with one singular value of
  1000 and the rest between 1 and 2 the estimate is above that norm at the fourth step,
  1.41 of it at the sixth and 3.46 of it at the fortieth. A matrix with no dominant value
  shows nothing of the kind until past its rank, and an earlier version of this page gave
  that as the rule. It is offered for a first idea of the operator's size, not as a bound,
  and not for a tolerance that has to mean something.
- **An operator pair that are not transposes of one another is not detected.** The
  estimates are then wrong, and the confirmations, made with the same pair, decide.
- **The size of the problem does not limit it.** The right-hand side is divided by the
  power of two at its largest entry and every image of either operator by the power of
  two at the largest entry of the transpose's first image. The iterate is multiplied back
  by one power of two when the call returns. Those divisions change no significand, so
  the same problem with b multiplied by a power of two, or A and λ by another, gives the
  same iterates bit for bit wherever the operators' own products stay in range.
- **Sizes inside one problem are another matter, and three things are said of them.**
  - **An entry of b more than 2^1022 below its largest is not seen by the iteration.**
    It is seen by the report: the residual is formed from b as it was given, entry by
    entry, and measured at its own size. For the identity and `b = (2^100, 2^-980)` at
    tolerances of zero the call returns `(2^100, 0)` with `LsqrExhausted` and a residual
    of 2^-980, where before Batch 21 it reported `LsqrResidualMet` and a residual of 0. A
    second call on that residual finds the rest.
  - **Where b is all but orthogonal to the operator's range the things iterated are not
    near 1.** The first image is then far below the operator's size, and a left or right
    vector can be 2^512 or more before it is made a unit. A norm whose sum of squares is
    above the largest number is measured again with the entries lowered, as one too small
    is measured again with them lifted, and the call goes on. The iterate is then held far
    below its unit and can lose digits; the norms formed say what was reached.
  - **An iterate with entries above 2^512** has a sum of squares that is not a number. It
    is read only under a damping, and is then measured again likewise.
  [The scale-sweep roster](rosters.md#what-is-declared) pins it: identical from 2^-800 to
  2^1020 with the right-hand side scaled and with the matrix and the damping scaled, and
  to 2^1000 with the matrix and the right-hand side scaled together, past which the
  right-hand side used there cannot be held.
- **`LsqrOutOfRange` is where the two powers of two do not serve.** The right-hand side's
  size over that of the transpose's first image is above the largest number or below the
  smallest normal one; or the damping is more than 2^500 times that image's size. That
  size is the operator's own, and the ratio the solution's scale, where b has a fair part
  in the operator's range; where it has not, both are off by as much as b is orthogonal
  to it, and a problem can be refused that another pair of powers would have served:
  `A = [1; 0]`, `b = (2^-501, 1)` under a damping of 1 is. A right-hand side the
  transpose sends to zero exactly is not refused: zero solves it. A
  damped solution is of the order of the operator over the damping's square, and beyond
  that ratio the iterate could not be held in the iteration's units. Up to it the problem
  is solved: the tests take a damping 2^480 times the operator.
- **A norm too small or too large to hold.** A gradient above the largest number is
  reported as infinity and meets no tolerance. A norm that is not zero and whose value in
  the caller's units rounds to zero is reported as zero and **does not meet a tolerance of
  zero**; it meets any other. A norm is brought to the caller's units by the two powers
  of two together where their product is a number the format holds, so that a small one
  and a large one do not lose between them a value that is in range: before Batch 21 a
  gradient of 2^-601 was reported as zero, and called met, where the operator was at
  2^-500 and the right-hand side at 2^500.
- **In the steps the operators are given vectors of unit size**, and they are assumed
  linear, as the method assumes. In a confirmation the transpose is given the residual at
  unit size, and the operator the iterate in the iteration's units, which is of the
  solution's size over the right-hand side's and grows with the condition number. One
  that is not linear, an absolute threshold inside it for instance, sees arguments of
  other sizes than the caller's.
- **If an operator signals a condition the call does not return**, and x then holds the
  iterate reached in the iteration's units, not the caller's.
- **The last step allowed makes a product it does not need for the iterate.** Its
  transpose product serves the last entry of the gradient's history and tells the end of
  the sequences from the limit. A value that is not finite in it ends the call as
  `LsqrNonFinite` with an iterate that was finished.

## The method

The bidiagonalization of Golub and Kahan started from b, with the least squares problem
solved on the bidiagonal matrix as it grows: the method Paige and Saunders published in
1982. [The design note](design-iterative-least-squares.md#the-method) derives it. In
short, with unit vectors u and v and positive α and β,

    β₁ u₁ = b,   α₁ v₁ = Aᵀ u₁
    repeat for k = 1, 2, …:
        β_{k+1} u_{k+1} = A v_k − α_k u_k
        α_{k+1} v_{k+1} = Aᵀ u_{k+1} − β_{k+1} v_k
        two reflections reduce the new column of [B; λ I], giving ρ_k, θ_{k+1}, φ_k
        x ← x + (φ_k / ρ_k) w,      w ← v_{k+1} − (θ_{k+1} / ρ_k) w

Each α and β is the norm that makes its vector a unit vector. The estimates are
`α_{k+1} β_{k+1} |φ_k| / ρ_k` for the gradient and, for the residual, the one number the
reduction could not remove where there is no damping, and with damping the root of what
it could not remove less `λ² ‖x‖²`, which is a difference of squares and a poor estimate
where the damped term is most of it. Elsewhere no squared norm is formed: each of the
reflections' lengths is the larger magnitude times `√(1 + q²)`, q the ratio of the
smaller to it, `λ² x` in the gradient is `λ (λ x)`, and the norm of a vector whose squares
are too small to hold is taken again with its entries lifted.

A vector u or v whose norm is exactly zero ends the sequences: nothing remains to add. So
do a right-hand side of zero, for which zero is returned with no product made, and one
the transpose sends to zero, for which zero is a minimizer.

## Using it

**With a stored matrix.** `apply` is [`matvec!`](primitives.md) and `apply-transpose` is
`transpose-matvec!`; [the example](../examples/operator-least-squares.ct) fits a line and
finds the shortest solution of an underdetermined system that way.

**With a product never formed.** For `A = J S`, `apply` is two products through a work
vector, `v ↦ J (S v)`, and `apply-transpose` is `u ↦ Sᵀ (Jᵀ u)`.

**With a preconditioner.** A right preconditioner is a change of variable: pass
`v ↦ A (M⁻¹ v)` and `u ↦ M⁻ᵀ (Aᵀ u)`, and apply M⁻¹ to the result. The tolerances, the
damping and the minimum norm are then those of the new variable. Dividing each column by
its norm, as [column scaling](scaling.md) gives it, is such an M and costs nothing.

**With a relative tolerance.** The tolerances are absolute. For the residual, a fraction
of the norm of b, which [`norm2`](primitives.md) gives. For the gradient, a fraction of
the operator's norm times that; where the operator's norm is not known,
`operator-norm-estimate` after a short first call is an estimate of it.

**To go further than one call reaches.** Without damping, a second call on the residual
of the first, `b − A x`, gives a correction to add to x; both are in the transpose's
range, so the sum still has the least norm. With damping that is another problem and the
sum is not its answer.

## Representation choices

The operators are functions passed as arguments, as in
[`cg!`](conjugate-gradient.md) and [`refine!`](refinement.md), so the workspace is only
storage and serves any operator of its shape. The right-hand side is copied, divided,
into the workspace, which is what lets it share an array with the solution and what the
confirmations read. The vector returned is the last iterate: in the three measures in
which the iterates improve, the last is the best.

The vectors u and v are not kept, so they cannot be made orthogonal again; that is the
cost of storing six vectors and not the basis, and it is why the step count is not
bounded by the rank.
