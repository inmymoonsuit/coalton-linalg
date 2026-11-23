# Iterative refinement

Package: `coalton-linalg/refinement`. Source: [refinement.ct](../src/refinement.ct).
State: completed as Batch 10; verified at safety 1 and safety 0. The [batch
record](batch10.md) holds the controls.

## What it is for

A solve is often available for a matrix that is close to the one wanted and not equal to
it: a factor of a regularized matrix, a factor computed at a shift, a factor that lost
digits to growth. Refinement turns that solve into a solution of the system that was
meant, using only products with the true matrix:

    r = b − A x,      solve the nearby system for c from r,      x ← x + c,

repeated until the residual is small enough or stops falling.

Two uses in this library:

- **One factorization where there were two.** A proof of positive definiteness factors
  `A − σ I` through [`factor-cholesky!`](cholesky-form.md); the step that follows wants a
  solve with A. `solve-cholesky!` solves with the proof factor, and refinement against A
  corrects the shift, so the second factorization is not needed.
- **A regularized augmented solve corrected to the unregularized one.** The
  [signed factorization](quasidefinite.md) of `[[top I, Aᵀ], [A, −bottom I]]` solves a
  regularized problem; refinement against `[[I, Aᵀ], [A, 0]]` returns the minimum-norm
  solution itself.

## Interface

- `make-refiner(dimension, step-limit) -> Refiner` allocates the workspace once: three
  working vectors, one vector for a shared right-hand side, and the norm history.
- `refine!(refiner, apply, solve, b, x, tolerance, stall) -> RefineStatus`.
  - `apply` is a function of two arrays, `(v, out)`, that **writes A v into out**. The
    library forms `r = b − A x` itself, in working precision, so every caller's residual
    is the same residual.
  - `solve` is a function of two arrays, `(r, out)`, that writes an approximate solution
    of `A c = r` into out and **returns whether it could**. Any solve will do: `solve!`,
    `solve-cholesky!`, a dense factor's, or the caller's own.
  - `b` is the right-hand side. `x` holds the starting vector on entry, commonly zeros or
    a first solve, and **the best iterate seen on return**, in every outcome. They may
    be the same array; residuals still use the right-hand side as it was on entry.
  - `tolerance` is absolute and is compared with the **infinity norm** of the residual.
  - `stall` is the caller's factor: a step that leaves the residual norm at or above
    `stall` times the previous one ends the call. 1 demands only that the norm fall; 1/2
    that it at least halve; a very large number turns the test off.
- `step-count(refiner)`: corrections the latest call applied.
- `best-step(refiner)`: which iterate was returned, 0 for the starting vector.
- `residual-norms(refiner)`, `residual-count(refiner)`: the infinity norm of the residual of
  each iterate of the latest call, entry k after k corrections, and how many were recorded.
  They are in the norm the tolerance is in, so they can be compared with a caller's own
  test directly.
- `dimension(refiner)`, `step-limit(refiner)`.

| RefineStatus | Meaning |
| --- | --- |
| RefineConverged | A residual norm reached the tolerance |
| RefineStalled | A step did not bring the norm below `stall` times the one before; the tolerance was not met |
| RefineLimit | The step limit was reached with the tolerance not met |
| RefineNonFinite | The operator, the solve or an update produced a value that was not finite |
| RefineSolveRefused | The solve function returned false |
| RefineBadSize | b or x does not have the refiner's dimension |
| RefineBadTolerance | The tolerance is negative or not a number, or the stall factor is not positive |

`RefineBadSize` and `RefineBadTolerance` are returned before anything is read or written.
A call allocates nothing; the two functions are the caller's and are made once, outside
the loop that calls `refine!`. An exception raised inside `apply` or `solve` passes
through, and x then holds the iterate reached, not necessarily the best. One refiner
serves one call at a time.

## What is and is not claimed

- **The vector returned is the best seen, not the last.** If a step makes the residual
  grow, the earlier iterate is restored. Ties go to the earlier iterate.
- **Refinement converges to the conditioning of the solve, not past it.** The residual is
  formed in the same precision as everything else. With an operator applied in working
  precision the residual cannot be driven below about `u (|A| |x| + |b|)`, and the
  solution cannot be made more accurate than the conditioning of A allows at that
  residual. What refinement removes is the error of the *approximate solve*: a shift, a
  regularization, a factorization's growth. A tolerance below the residual floor ends in
  `RefineStalled` or `RefineLimit`, with the best iterate returned.
- **Convergence needs the solve to be close enough.** If the solve is of `B` in place of
  A, each step multiplies the error by `I − B⁻¹ A`. For `B = A − σ I` with A positive
  definite that is a factor of at most `σ / (λ_min − σ)`: a shift small beside the
  smallest eigenvalue converges in a step or two, and one near it converges slowly or not
  at all. The stall test is how a caller finds out; it then has the best iterate and can
  fall back.
- **No rate is promised and none is assumed.** The call stops on what it measures.
- **The tolerance is the caller's.** It is absolute; a caller wanting a relative test
  multiplies by its own scale.

## Using it

**With a proof factor.** After `factor-cholesky!(symbolic, values, σ)`:

    apply = v ↦ A v                         (for A = J Jᵀ: two products, J (Jᵀ v), A never formed)
    solve = r ↦ solve-cholesky!(factor, r, ·) == SolveSuccess
    refine!(refiner, apply, solve, b, x, tolerance, stall)

with x zeros, or the result of one `solve-cholesky!`.

**With a regularized augmented factor.** After `factor-signed!` of
`[[I, Aᵀ], [A, −bottom I]]`, `apply` is `(x, y) ↦ (x + Aᵀ y, A x)`, the right-hand side is
`(0, b)`, and the leading part of the refined vector is the minimum-norm solution of
`A x = b`.

The tests do both and compare the results with exact rational solutions.

## Derivation

Let the solve be exact for a matrix B, and let `x*` solve `A x* = b`. Then
`c = B⁻¹ (b − A x) = B⁻¹ A (x* − x)`, so the new error is `(I − B⁻¹ A)(x* − x)`. The
iteration converges for every starting vector exactly when the spectral radius of
`I − B⁻¹ A` is below 1, and the residual, `A (x* − x)`, obeys the same recurrence with
`I − A B⁻¹`. For `B = A − σ I` and symmetric positive definite A with eigenvalues λ, the
eigenvalues of `I − B⁻¹ A` are `−σ / (λ − σ)`.

For the augmented pair, `K_δ = [[I, Aᵀ], [A, −δ I]]` and `K₀ = [[I, Aᵀ], [A, 0]]` differ
only in the trailing block, and `I − K_δ⁻¹ K₀ = K_δ⁻¹ (K_δ − K₀)` has the eigenvalues
`δ / (δ + μ)` for the eigenvalues μ of `A Aᵀ`, and zero. So the iteration converges for
every `δ > 0` when A has full row rank, quickly when δ is small beside the smallest
eigenvalue of `A Aᵀ`.

A computed solve is not exact for any one B, and the argument is then about what is
measured: each residual is recorded, and the call acts on the record.

## Representation choices

The operator and the solve are functions, so that one call serves every factor in the
library and matrices that are never assembled. They are passed as arguments and not
stored in the refiner, so a refiner is only workspace and can be reused for another
system of the same size. The history is kept in the refiner and read by accessors, in
the manner of the factorizations' latest-call diagnostics.
