# Conjugate gradients

Package: `coalton-linalg/conjugate-gradient`. Source:
[conjugate-gradient.ct](../src/conjugate-gradient.ct).
State: completed as Batch 13; verified at safety 1 and safety 0. The [batch
record](batch13.md) holds the controls.

## What it is for

A system `A x = b` whose matrix is symmetric positive definite and is available only as
something that multiplies a vector: a matrix too large to factor, a product of matrices
that is never formed, such as `J (Jᵀ v)`, or a matrix-free operator. The method needs
nothing of the matrix but its product with a vector, four vectors of storage, and, to be
fast, an approximate solve to precondition it with. The factorizations of this library
are the tool where they fit in memory and time; this is the tool where they do not, and
the two combine, since any of their solves is a preconditioner.

## Interface

- `make-cg(dimension, step-limit) -> CgWorkspace` allocates the workspace once: four
  vectors of the dimension and one number a step.
- `cg!(workspace, apply, b, x, tolerance) -> CgStatus`.
  - `apply` is a function of two arrays, `(v, out)`, that **writes A v into out**.
  - `b` is the right-hand side. `x` holds the starting vector on entry, commonly zeros,
    and the last iterate on return.
  - `tolerance` is absolute and is compared with the **2-norm** of the residual
    `b − A x`.
- `pcg!(workspace, apply, precondition, b, x, tolerance) -> CgStatus`: the same with a
  preconditioner, a function of two arrays, `(r, z)`, that writes into z an approximate
  solution of `A z = r` and **returns whether it could**. It must itself act as a
  symmetric positive definite matrix on r: `solve!` with a factor of a nearby matrix, a
  division by the diagonal, any fixed combination of such.
- `step-count(workspace)`: steps the latest call took. Each applies the operator once.
- `residual-norm(workspace)`: the 2-norm of `b − A x` for the vector returned, **formed
  from the operator**, after a call that ended `CgConverged` or `CgLimit`; −1 otherwise.
- `recurrence-norms(workspace)`: the 2-norm of the running residual after each step, entry
  0 for the starting vector.
- `restart-count(workspace)`: how often the running residual was replaced by one formed
  from the operator.
- `dimension(workspace)`, `step-limit(workspace)`.

| CgStatus | Meaning |
| --- | --- |
| CgConverged | The residual of the vector returned, formed from the operator, is within the tolerance |
| CgLimit | The step limit was reached with the tolerance not met; `residual-norm` says how far it got |
| CgNotPositive | A direction with no positive curvature, or a preconditioned residual with no positive product with the residual, was met |
| CgNonFinite | The operator, the preconditioner or an update produced a value that was not finite |
| CgPreconditionerRefused | The preconditioner returned false |
| CgBadSize | b or x does not have the workspace's dimension |
| CgBadTolerance | The tolerance is negative or not a number |

`CgBadSize` and `CgBadTolerance` are returned before anything is read or written. A call
allocates nothing; the two functions are the caller's and are made once, outside the loop
that calls. An exception raised inside either passes through. One workspace serves one
call at a time.

**What x holds on the other outcomes.** `CgNotPositive` and `CgPreconditionerRefused` leave
the iterate reached, the starting vector if no step was completed; it is a valid iterate
and no worse than the start in the sense of the next section. A step length that is not
finite is reported before x is touched. After any other `CgNonFinite`, x is the iterate
at which the value arose and may itself hold one.

## What is and is not claimed

- **Convergence is reported on the residual of the vector, not on the method's own
  running one.** The method carries a residual forward by a recurrence, which is what
  makes a step cost one product, and that recurrence drifts from `b − A x`. When the
  running residual meets the tolerance the residual is formed from the operator, at the
  cost of one more product, and the test is made on that. If it fails, the method starts
  again from the residual just formed and carries on. `CgConverged` therefore means what
  it says of the vector returned, to the rounding of one product. At the limit the report
  is likewise the residual formed from the operator.
- **Each iterate is at least as near the solution as the one before, in the norm of A.**
  That is the property that defines the method: the k-th iterate minimizes
  `(x − x*)ᵀ A (x − x*)` over the starting vector plus the span of the first k directions,
  and the spans are nested. The residual's 2-norm has no such property and can rise from
  one step to the next; `recurrence-norms` shows it as it was. The tests measure the error
  of every iterate in the norm of A, exactly, and find it never rises by more than
  rounding.
- **In exact arithmetic the method ends in at most as many steps as the dimension.** In
  floating point it does not: the directions lose their conjugacy, and a badly
  conditioned system takes more steps than its dimension or stalls above the tolerance.
  The number of steps to a given accuracy grows with the root of the condition number of
  the preconditioned matrix. **No rate is promised**; the call stops on what it measures,
  and the limit is the caller's.
- **A residual within tolerance is not an error within tolerance.** The error can exceed
  the residual by the norm of the inverse. A caller who needs the error bounded needs a
  bound on the smallest eigenvalue, which this package does not supply.
- **`CgNotPositive` is a finding about the direction met, not a test of the operator.** A
  call can converge on an indefinite operator without meeting such a direction, and an
  operator that is not symmetric is not detected at all.
- **The tolerance is the caller's**, absolute, in the 2-norm, which is formed as the root
  of a plain sum of squares and so overflows for vectors of magnitude near 10¹⁵⁴; such a
  system is scaled first.

## The method

With M the preconditioner as a matrix, `z = M⁻¹ r`:

    r = b − A x,   z = M⁻¹ r,   p = z,   ρ = rᵀ z
    repeat:
        q = A p
        α = ρ / (pᵀ q)
        x ← x + α p
        r ← r − α q
        z = M⁻¹ r,   ρ' = rᵀ z
        p ← z + (ρ' / ρ) p,   ρ ← ρ'

The step length α is the one that minimizes the error in the norm of A along p, and the
new direction is the preconditioned residual made conjugate to the old one, `pᵀ A p_old =
0`; that the one correction makes it conjugate to all the earlier ones is the theorem
the method rests on. `pᵀ q` is the curvature of the quadratic along p and is positive for
every nonzero p when A is positive definite; ρ is positive for every nonzero r when M is.
The call checks both at every step and stops with `CgNotPositive` where one fails, before
the step is taken.

The plain method is the same with M the identity.

## Using it

**With a product that is never formed.** For `A = J Jᵀ`, `apply` is two products through a
work vector, `v ↦ J (Jᵀ v)`.

**With a factor as preconditioner.** A factor of a matrix near A, a shifted one, an
earlier one in a sequence, one with small entries dropped, gives
`precondition = r ↦ solve!(factor, r, ·) == SolveSuccess`. With a factor of A itself the
method ends in one step, which is how the tests check the preconditioned form.

**With the diagonal.** `precondition` divides by the diagonal of A. It costs nothing and
helps where the difficulty is one of scale.

## Representation choices

The operator and the preconditioner are functions passed as arguments, as in
[`refine!`](refinement.md), so the workspace is only storage and serves any system of its
size. The history is kept in the workspace and read by accessors. The 2-norm is used, and
not the infinity norm of `refine!`, because the method's own quantities are inner
products and the running norm then costs nothing more.

The vector returned is the last iterate and not the best seen: the measure in which the
iterates improve, the norm of A, is not one the call can evaluate, and in that measure
the last is the best.
