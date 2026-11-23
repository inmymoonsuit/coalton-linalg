# Design: least squares on an operator, by bidiagonalization

State: **built as Batch 20**; [the contract](lsqr.md) is the description of what exists
and [the record](batch20.md) says where building departed from this page. The text below
is the design as it was written, before any code, kept for its derivation. It answers the
roadmap item "iterative least squares and minimum norm on an operator, with damping". No
consumer asked for it; it is for the library's own users, and it changes no existing
call. Where it says "in exact arithmetic" it means a property of the method
that rounding weakens and that the tests would measure, not one the call proves.

## The question

A matrix A with m rows and n columns is available only as two functions: one that
writes `A v` and one that writes `Aᵀ u`. It may be too large to factor, a product never
formed, or have a factor that fills. A caller wants

    minimize over x:   ‖b − A x‖₂² + λ² ‖x‖₂²

for a damping λ ≥ 0 it supplies, and with λ = 0, where the minimizer is not unique, the
one of least norm. [Conjugate gradients](conjugate-gradient.md) answer this through
`Aᵀ A + λ² I`, whose condition number is the square of the problem's. The method here
works on A and Aᵀ and never forms that product, in the operator or in the recurrence.

Three things are wanted of the answer, in this library's usual order:

1. **What the call reports is formed from the operators**, not read from the method's
   running quantities: the residual `b − A x` and the gradient `Aᵀ(b − A x) − λ² x` of
   the vector returned.
2. **Every outcome has a name**, and a refusal is made before anything is written.
3. **The size of the problem does not limit it**: the same problem at another power of
   two is solved bit for bit the same.

## The method

The bidiagonalization of Golub and Kahan, started from b, with the least squares problem
solved on the bidiagonal matrix as it grows: the method Paige and Saunders published in
1982. It is set out here from the mathematics, because the code is to be written from
this page.

### The two sequences

Unit vectors u₁, u₂, … of length m and v₁, v₂, … of length n, and positive numbers α_k,
β_k, defined by

    β₁ u₁ = b,                      α₁ v₁ = Aᵀ u₁,
    β_{k+1} u_{k+1} = A v_k − α_k u_k,
    α_{k+1} v_{k+1} = Aᵀ u_{k+1} − β_{k+1} v_k,

each α and β being the norm that makes its vector a unit vector. In exact arithmetic the
u are orthonormal and the v are orthonormal, and with `U_{k+1} = [u₁ … u_{k+1}]`,
`V_k = [v₁ … v_k]` and B_k the k+1 by k matrix with α₁ … α_k on its diagonal and
β₂ … β_{k+1} below it,

    A V_k = U_{k+1} B_k,            Aᵀ U_{k+1} = V_k B_kᵀ + α_{k+1} v_{k+1} e_{k+1}ᵀ.

Each step costs one product with A and one with Aᵀ.

### The small problem

Take `x_k = V_k y`. Then `b − A x_k = U_{k+1} (β₁ e₁ − B_k y)`, and since the columns of
U and of V are orthonormal the caller's problem restricted to the span of V_k is

    minimize over y:   ‖β₁ e₁ − B_k y‖² + λ² ‖y‖².

That is a least squares problem with the matrix `[B_k; λ I]`, and it is solved by
reducing that matrix to upper bidiagonal form with plane reflections, two a column,
which needs only the column just added. Before step k the reduction has left a diagonal
entry ρ̄_k not yet combined with what is below it and a right-hand side entry φ̄_k, with
`ρ̄₁ = α₁` and `φ̄₁ = β₁`. Step k:

    against the damping:     ρ̃ = hypot(ρ̄_k, λ),   c₁ = ρ̄_k / ρ̃,   s₁ = λ / ρ̃,
                             ψ_k = s₁ φ̄_k,          φ̃ = c₁ φ̄_k
    against the next row:    ρ_k = hypot(ρ̃, β_{k+1}),   c = ρ̃ / ρ_k,   s = β_{k+1} / ρ_k,
                             θ_{k+1} = s α_{k+1},   ρ̄_{k+1} = −c α_{k+1},
                             φ_k = c φ̃,            φ̄_{k+1} = s φ̃

The reduced matrix R_k has ρ₁ … ρ_k on its diagonal and θ₂ … θ_k above it, the reduced
right-hand side is (φ₁ … φ_k), and what cannot be removed is φ̄_{k+1} and the ψ. So the
small problem's solution is `R_k y_k = f_k`, and

    ‖b − A x_k‖² + λ² ‖x_k‖² = φ̄_{k+1}² + ψ₁² + … + ψ_k².

### The iterate without the basis

`x_k = V_k R_k⁻¹ f_k`. The columns of `V_k R_k⁻¹` satisfy a two-term recurrence because R
is bidiagonal; writing w_k for the k-th of them times ρ_k,

    w₁ = v₁,        x_k = x_{k−1} + (φ_k / ρ_k) w_k,        w_{k+1} = v_{k+1} − (θ_{k+1} / ρ_k) w_k,

so nothing of V or U is kept but the latest vector of each.

### What the scalars say of the iterate

- **The gradient.** With `t = β₁ e₁ − B_k y_k`, the relations above give
  `Aᵀ(b − A x_k) − λ² x_k = V_k (B_kᵀ t − λ² y_k) + α_{k+1} t_{k+1} v_{k+1}`. The first
  term is the small problem's own optimality condition and is zero. The last row of B_k
  holds only β_{k+1}, so `t_{k+1} = −β_{k+1} (y_k)_k`, and the last equation of
  `R_k y_k = f_k` gives `(y_k)_k = φ_k / ρ_k`. So

      ‖Aᵀ(b − A x_k) − λ² x_k‖ = α_{k+1} β_{k+1} |φ_k| / ρ_k.

- **The residual.** From the identity above, `‖b − A x_k‖² = φ̄_{k+1}² + Σψ² − λ² ‖x_k‖²`,
  the norm of x_k being accumulated as the iterate is updated. Without damping it is
  `|φ̄_{k+1}|` alone. With damping it is a difference, and a poor estimate where the
  damped term is most of the sum.

Both are exact in exact arithmetic and drift from the truth in floating point, as the u
and the v lose their orthogonality. They are used here for one thing only: to decide
when it is worth forming the true quantities.

### Where the sequences end

If `β_{k+1} = 0` the span of the v is invariant, `φ̄_{k+1} = 0` and x_k solves the problem
exactly. If `α_{k+1} = 0` the gradient is zero and x_k is a minimizer. Either ends the
iteration; so do `b = 0`, where x = 0, and `Aᵀ b = 0`, where x = 0 is a minimizer. In
exact arithmetic one of these happens within min(m, n) steps. In floating point a norm
is rarely exactly zero and the method is stopped by its tolerances or its limit.

### The minimum norm

Every x_k is a combination of the v, and every v is in the range of Aᵀ. Among the
minimizers of `‖b − A x‖` exactly one is in the range of Aᵀ, and it is the one of least
norm. So with λ = 0 the method converges to the minimum-norm least squares solution,
whatever the rank and the shape of A, **because it starts from zero**: this is why the
call takes no starting vector. In floating point each v acquires a component outside
that range by rounding, which the operator cannot see and nothing removes.

### Three monotone quantities

In exact arithmetic, from one iterate to the next: the damped objective
`‖b − A x_k‖² + λ²‖x_k‖²` does not rise, since each iterate minimizes it over a larger
span; the norm of x_k does not fall; and the distance from x_k to the solution does not
rise. The last two are properties of conjugate gradients started from zero, which this
method is, in exact arithmetic, on the normal equations.

## Interface

Package `coalton-linalg/lsqr`, one new source file, in the shape of
[conjugate gradients](conjugate-gradient.md#interface).

- `make-lsqr(rows, columns, step-limit) -> LsqrWorkspace` allocates once: three vectors
  of each length and two histories.
- `lsqr!(workspace, apply, apply-transpose, b, x, damping, residual-tolerance,
  gradient-tolerance) -> LsqrStatus`.
  - `apply` is a function of two arrays `(v, out)` that writes `A v` into out, v of
    length `columns` and out of length `rows`; `apply-transpose` writes `Aᵀ u`.
  - `b` is read once, before x is written, so **b and x may be one array** where the
    matrix is square. `x` is written only: there is no starting vector.
  - `damping` is λ, not its square.
  - The tolerances are absolute, each on the 2-norm of a quantity formed from the
    operators for the vector returned: `b − A x`, and `Aᵀ(b − A x) − λ² x`. Either may be
    zero, which that quantity then meets only by being zero.
- `step-count`, `confirmation-count`: steps taken, and how often the true quantities
  were formed. A step is two products and a confirmation two more.
- `residual-norm`, `gradient-norm`: the two norms for the vector returned, formed from
  the operators; −1 where no vector was returned.
- `recurrence-residual-norms`, `recurrence-gradient-norms`: the running estimates after
  each step, as they were.
- `operator-norm-estimate`: the Frobenius norm of the bidiagonal matrix reached. In exact
  arithmetic it is at most the Frobenius norm of A and it grows with the steps; it is an
  estimate from below and proves nothing.
- `rows`, `columns`, `step-limit`.

| LsqrStatus | Meaning |
| --- | --- |
| LsqrResidualMet | `‖b − A x‖`, formed from the operator, is within the residual tolerance |
| LsqrGradientMet | The gradient, formed from the operators, is within its tolerance, and the residual is not within its own |
| LsqrExhausted | The sequences ended, so that in exact arithmetic x is the solution, and neither norm as formed is within its tolerance |
| LsqrLimit | The step limit was reached with neither tolerance met |
| LsqrNonFinite | b, an operator or an update produced a value that was not finite |
| LsqrOutOfRange | The right-hand side's size divided by the operator's is not a number the format holds, so neither is the solution's scale |
| LsqrBadSize | b or x does not have the workspace's length for it |
| LsqrBadTolerance | A tolerance is negative or not a number |
| LsqrBadDamping | The damping is negative, infinite or not a number |

The three `Bad` statuses and `LsqrOutOfRange` are returned with x untouched. After the
first four x is the last iterate and both norms are reported. After `LsqrNonFinite` x is
the last iterate completed, or untouched where the call ended before its first step.

## Decisions

**D1. One method, this one.** Conjugate gradients on the normal equations written out
directly are the same iterates in exact arithmetic; the paper that introduced this form
compared the two and found this one the more reliable on badly conditioned problems,
which is not measured here. A variant that makes the gradient fall monotonically exists
and is left for a caller who needs it.

**D2. Confirmation instead of trust.** When a running estimate is within its tolerance
the residual and the gradient are formed from the operators, two products, and the
tolerances are tested on those. If neither is met the iteration goes on, **from where it
was**: the sequences are still what they were and only the estimate was wrong. So that
a wrong estimate does not buy a confirmation at every step, the estimate that triggered
must then fall by the factor it was out by: its trigger becomes the tolerance times the
estimate over the true value. A trigger that would be zero is switched off. At the limit
and at the end of the sequences the true quantities are formed whatever the estimates
say. Conjugate gradients here restart from the residual formed. This method does not: a
damped problem restarted from a nonzero vector has a right-hand side under its damping
rows, which the two reflections a step do not carry, and one rule is kept for both.

**D3. No stall status.** A tolerance finer than rounding allows ends at the limit, as in
conjugate gradients, with what was reached reported. A rule that named a stall would be
a judgement on two confirmations, and the limit is the caller's.

**D4. Absolute tolerances, the caller's.** The usual relative tests divide by an estimate
of the operator's norm that the method makes as it goes. Here the estimate is offered
for reading and the tolerances are numbers the caller gives, as everywhere in this
library.

**D5. Everything iterated is near 1.** The right-hand side is divided by the power of
two at its largest entry and each image of the operator by the power of two at the
largest entry of `Aᵀ u₁`. The u, v and w are then of unit size, the α, β, ρ, θ and φ are
near 1, and the iterate is accumulated in those units and multiplied back, by one power
of two, when the call returns. All of those divisions are exact, so the same problem with
b multiplied by a power of two, or A and λ by another, gives the same iterates bit for
bit wherever the operator's own products stayed in range. The operators are given
vectors of unit size, and are assumed linear.

**D6. The squared norms are not formed.** Each `hypot` is computed as the larger
magnitude times `√(1 + q²)`, q the ratio of the smaller to it, and λ² x in the gradient as
λ (λ x). A norm of a vector that is all tiny is measured again with its entries lifted,
as in conjugate gradients.

**D7. No reorthogonalization.** The u and v are not kept, so they cannot be made
orthogonal again. On a badly conditioned problem the method then takes more steps than
min(m, n). This is the cost of storing six vectors and not the basis.

**D8. Preconditioning is composition.** A right preconditioner is a change of variable:
the caller passes `v ↦ A (M⁻¹ v)` and `u ↦ M⁻ᵀ (Aᵀ u)` and applies M⁻¹ to the result. The
tolerances, the damping and the minimum norm are then in the new variable. No argument
is added for it; the contract page would show it with
[column scaling](scaling.md).

## What would be claimed, and what not

- **The two norms reported are of the vector returned**, formed from the operators, to
  the rounding of two products. `LsqrResidualMet` and `LsqrGradientMet` mean what they
  say of them.
- **A gradient within tolerance is not an error within tolerance.** The error can exceed
  the gradient by the square of the reciprocal of the smallest singular value, or by
  1/λ² with damping. Nothing here bounds it.
- **No rate and no step count is promised.** The count grows with the condition number
  of A, not with its square root as for conjugate gradients on a definite matrix.
- **The minimum norm, the three monotone quantities and the end within min(m, n) steps
  are properties in exact arithmetic.** The tests would measure each exactly on small
  matrices and report by how much rounding breaks it.
- **An operator pair that are not transposes of one another is not detected.** The
  estimates are then wrong and the confirmations, which use the same pair, decide.

## Verification planned

1. **Against exact arithmetic.** Seeded matrices of small integers on shapes tall, square
   and wide, of full and of deficient rank, with rows and columns of zeros, at several
   dampings from zero up. The exact minimizer, of minimum norm where λ = 0, is found in
   rational arithmetic. Run to a fine tolerance, the vector returned is near it; the
   norms reported are the exact norms of the exact residual and gradient of the vector
   returned; a status that says a tolerance is met is true of those exact norms.
2. **The minimum norm.** On the rank-deficient cases, the exact component of x in the
   null space of A, against the norm of x.
3. **The monotone quantities**, measured exactly for every iterate, each iterate being
   the result of a call with that step limit.
4. **An independent recomputation** of the recurrences in the test language reproduces
   every iterate bit for bit, and the histories.
5. **Every status**, each by a case made for it, with x untouched where that is said. A
   badly conditioned case run below what rounding allows must show confirmations that
   failed and must not return a status that is false of the exact norms.
6. **Scale.** The same iterates with b multiplied by powers of two across the range, and
   with A and λ; entries in [the scale-sweep roster](rosters.md) that pin how far.
7. **Aliasing.** b and x one array, in [the aliasing roster](rosters.md).
8. **Allocation.** Zero bytes a call after warm-up, through a compiled caller.
9. **Planted source defects**, each applied alone.
10. **One measured problem** of a size that is not small, sparse, through the existing
    in-place products: steps, confirmations and time, beside conjugate gradients on its
    normal equations.

## Not in this design

- A starting vector, for the reasons under *The minimum norm* and D2. An undamped
  solution can be improved by a second call on its residual, added to the first, which
  the contract page would show.
- Reorthogonalization, a condition estimate, an error bound.
- The variant that minimizes the gradient, and the method for consistent systems that
  minimizes the error.
- Bounds on the unknowns; that is another roadmap item.

## Sources

The bidiagonalization is that of Golub and Kahan (1965) and its use for least squares,
with damping carried by a second reflection a step, that of Paige and Saunders (1982).
The recurrences, the two scalar identities and the properties above are derived on this
page and the code would be written from it. The confirmation rule, the scaling and the
statuses are this library's. No implementation and no listing would be opened.
