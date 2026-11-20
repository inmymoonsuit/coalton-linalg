# coalton-linalg

An independent numerical linear algebra library written in Coalton, for 64-bit SBCL, in
double precision. It offers a sparse symmetric positive definite factorization with the
sparse operations around it, and dense factorizations for the problems that one cannot
take: indefinite systems, rank-deficient least squares, eigenvalues and singular values.

The library takes matrices and vectors and returns factors, solutions and counts.
Thresholds, rank decisions and iteration belong to the caller: nothing here decides a
rank, perturbs a pivot or retries. No application code is loaded or required.

## What it contains

| Package | What it does | Contract | Record |
| --- | --- | --- | --- |
| `coalton-linalg/errors` | Typed exceptions for malformed input and exhausted budgets | [errors](docs/errors.md) | |
| `coalton-linalg/numeric` | Euclidean and infinity norms that do not overflow early | [primitives](docs/primitives.md) | [1 and 2](docs/batches-1-2.md) |
| `coalton-linalg/refinement` | Iterative refinement of any approximate solve against an operator the caller applies: stopping rule, stall test, best iterate and residual history | [refinement](docs/refinement.md) | [10](docs/batch10.md) |
| `coalton-linalg/conjugate-gradient` | Conjugate gradients, plain and preconditioned, for a symmetric positive definite operator the caller applies; convergence reported on the residual formed from the operator | [conjugate gradients](docs/conjugate-gradient.md) | [13](docs/batch13.md) |
| `coalton-linalg/sparse/csc` | Checked CSC and COO construction, transpose, permutations, products into caller storage | [primitives](docs/primitives.md) | [3b](docs/batch3b.md) |
| `coalton-linalg/sparse/gram` | AᵀA and AAᵀ on a fixed pattern: analyze once, refill values | [gram](docs/gram.md) | [3c](docs/batch3c.md) |
| `coalton-linalg/sparse/spd-ldlt` | Sparse `P (A − shift I) Pᵀ = L D Lᵀ` with exact minimum-degree analysis done once, repeated factor and solve, a 1-norm condition estimate, a Cholesky-form entry for positivity proofs, a signed-pivot entry that measures its own growth, and an entry that leaves small pivots out and factors what is kept | [factorization](docs/spd-ldlt.md), [condition](docs/condition.md), [Cholesky form](docs/cholesky-form.md), [signed pivots](docs/quasidefinite.md), [skipped pivots](docs/skipping.md) | [3a](docs/batch3a.md), [3d](docs/batch3d.md), [6a](docs/batch6a.md), [7](docs/batch7.md), [7a](docs/batch7a.md), [8](docs/batch8.md), [10](docs/batch10.md), [11](docs/batch11.md) |
| `coalton-linalg/sparse/dependent-rows` | Which rows of a sparse matrix depend on the others: a kept set, a left null vector for each row left out, and proved bounds on what each leaves and on the singular values beyond the kept count | [dependent rows](docs/dependent-rows.md) | [11](docs/batch11.md) |
| `coalton-linalg/sparse/lu` | Sparse `P A Q = L U` for a square nonsymmetric matrix, with partial pivoting chosen from the values, solves with the matrix and its transpose, measured growth and a bound on the backward error | [LU](docs/lu.md) | [12](docs/batch12.md) |
| `coalton-linalg/sparse/scaling` | Row and column equilibration, and diagonal scaling of a symmetric positive definite matrix, by powers of two: exact, with the vector maps | [scaling](docs/scaling.md) | [9](docs/batch9.md) |
| `coalton-linalg/sparse/augmented` | The augmented matrix `[[top I, Aᵀ], [A, −bottom I]]` of a sparse A on a fixed pattern: analyze once, refill values | [signed pivots](docs/quasidefinite.md) | [7](docs/batch7.md) |
| `coalton-linalg/dense/indefinite-ldlt` | Dense symmetric indefinite LDLᵀ with 1-by-1 and 2-by-2 pivots, solves, and inertia | [indefinite](docs/indefinite-ldlt.md) | [4](docs/batch4.md) |
| `coalton-linalg/dense/orthogonal` | Dense column-pivoted QR, complete orthogonal decomposition, minimum-norm least squares and a nullspace basis at a rank the caller chooses | [orthogonal](docs/orthogonal.md) | [5](docs/batch5.md) |
| `coalton-linalg/dense/jacobi` | Dense symmetric eigenvalues, and singular values with a complete right factor, by Jacobi rotations | [jacobi](docs/jacobi.md) | [6](docs/batch6.md) |

The [interface reference](docs/api.md) indexes every exported symbol. The symmetric sparse
factorization does no numerical pivoting: its order is fixed at analysis, which is safe for
a positive definite matrix and, for the signed-pivot entry, is checked by a growth figure
the caller reads after every factorization; the dense routines
are unblocked and sequential, for small and moderate sizes. The [roadmap](ROADMAP.md) says
what is and is not planned.

Conventions the packages share:

- **Allocate once, repeat without allocating.** A handle is made for a size or a pattern
  and reused. A successful repeated call allocates nothing after warm-up, and each record
  measures that through a compiled Coalton caller at both safety policies.
- **Outcomes are values, misuse is an exception.** A refused pivot, a singular factor or an
  exhausted sweep limit is a constructor of a small enum the caller matches. Malformed
  structure, wrong dimensions and exceeded storage budgets raise a typed `LinalgError`.
- **A refusal leaves nothing half-written.** A solve writes its output only after every
  value is finite, and a refused factorization withdraws whatever stood before. The two
  operations that work in place say so in their contracts.
- **Handles are private nominal types.** Their inspection arrays are borrowed views,
  read-only by contract and rewritten by the next call.
- **Storage is refused by name**, against an explicit entry budget, before it is allocated.

## How it is checked

Each record states what was run. The pattern is the same throughout:

- **Exact references.** Results are compared with exact rational arithmetic on the stored
  F64 values: products multiplied out, inertia by exact elimination, minimum-norm solutions
  from an exact rank factorization, eigenvalues bracketed by exact counts. No floating
  routine is trusted as the reference for another.
- **Planted defects.** Each batch from 3c on has a list of one-line source faults, each
  applied alone and watched to fail the tests. Where a fault was missed, the record says
  so and what changed.
- **Compile probes.** External Coalton callers that must be rejected, for private
  construction, destructuring and field access, beside one that must compile and run.
- **Outside comparison, afterward.** NumPy and SciPy check generated matrices and give a
  timing reference, after the implementation exists. The records give the factors by which
  this library is slower.

[Provenance](docs/provenance.md) names the published method behind each algorithm and says
what was and was not consulted while writing it.

## Reviewing it

Start with the [interface reference](docs/api.md) and the six compiled examples under
`examples/`. For a reader who knows Coalton better than this library's authors do, these are
the places where a second opinion is most wanted:

- **Handle and status representation.** `(repr :transparent)` over a private struct for
  handles, `(repr :enum)` for outcomes, and latest-call accessors in place of result
  objects. The reasoning is in [the factorization contract](docs/spd-ldlt.md#representation-choices)
  and [the design review](docs/design.md#coalton-interface-review-2026-11-19).
- **Hot loops.** `for` with `:while`, `cell`s as accumulators, `(inline (finite? x))` on
  every checked value, and the block-compilation settings in the system definition.
  [The dense factorization](src/dense/indefinite-ldlt.ct) is the shortest complete example.
- **Where a double is boxed.** One helper returning an F64 allocated 16 bytes a call until
  it wrote into an array instead, while others in the same file did not;
  [the Batch 5 record](docs/batch5.md#two-things-found-while-building-it) has the case.
- **`lisp` forms.** There are few: constants (powers of two and the representable array
  size) and two calls to `cl:sort` during sparse analysis.
- **The safety policy** is a reader feature set by the check runner, and every check runs
  at both settings.

## Load and use

Register this directory with ASDF, with Coalton and its dependencies installed.
Start a fresh Lisp process with COALTON_ENV=release before loading Coalton;
development-mode performance is unverified. See the toolchain notes below.

```lisp
(asdf:load-system "coalton-linalg")
(asdf:load-system "coalton-linalg/examples")
(coalton:coalton (coalton-linalg/examples:solve-2x2)) ; => #(1.0d0 2.0d0)
```

The public packages are `coalton-linalg/errors`, `coalton-linalg/numeric`,
`coalton-linalg/refinement`, `coalton-linalg/conjugate-gradient`,
`coalton-linalg/sparse/csc`, `coalton-linalg/sparse/gram`, `coalton-linalg/sparse/scaling`,
`coalton-linalg/sparse/augmented`,
`coalton-linalg/sparse/spd-ldlt`, `coalton-linalg/sparse/dependent-rows`, `coalton-linalg/sparse/lu`,
`coalton-linalg/dense/indefinite-ldlt`, `coalton-linalg/dense/orthogonal` and
`coalton-linalg/dense/jacobi`. The [interface reference](docs/api.md) indexes
all exported symbols and links to their contracts. The compiled examples are
[an SPD solve](examples/spd-solve.ct), [sparse construction and products](examples/sparse-operations.ct),
[normal equations](examples/normal-equations.ct), [a minimum-norm solve through an augmented
system](examples/minimum-norm.ct), [dependent constraints](examples/dependent-constraints.ct),
[a rank from singular values](examples/singular-rank.ct),
[a sparse minimum-norm solve with signed pivots](examples/augmented-minimum-norm.ct),
[a solve through a scaled matrix](examples/scaled-solve.ct),
[redundant constraints found and bounded](examples/redundant-rows.ct),
[a nonsymmetric solve with row exchanges](examples/nonsymmetric-solve.ct) and
[a solve with a product that is never formed](examples/matrix-free-solve.ct). The
[mathematical glossary](CONTEXT.md) explains shared terms.

The main system loads only the library. `/tests`, `/examples`, `/tools`, and
`/tools-tests` are separate ASDF systems. Python is only used for independent
verification, never by the library.

## Run the complete checks

Copy `config.local.example.psd1` to `config.local.psd1` and point it at existing
SBCL, Python, CoaltonRoot, and QuicklispRoot. The [setup recipe](docs/setup.md)
specifies how to obtain the pinned sources. Local configuration is ignored by
Git. The check runner registers only those roots; it needs no private setup file
and downloads nothing.

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\check.ps1
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\check.ps1 -Safety 0
```

Run these **serially**; each takes about a minute. They rebuild the library for the
chosen safety policy and run every roster the batch records describe: interface and
allocation regressions, compiled Coalton callers and examples, the exact-reference tests of
each package, the external compile probes, the artifact drivers' own regressions, and the
independent Python verifiers with their planted controls. Generated artifacts stay under
`build/`. See the [initial review evidence](docs/initial-review.md) for fresh-source setup
verification and [pinned toolchain details](docs/toolchain.md) for the measurement protocol.
The planted source defects are run by hand for each batch, one process at a time, and are
recorded in the batch records, not rerun by this command.

## Repository map

| Path | Purpose |
| --- | --- |
| `src/errors.ct` | Typed validation and primitive arithmetic errors |
| `src/numeric.ct` | Stable norms |
| `src/refinement.ct` | Iterative refinement against an operator |
| `src/conjugate-gradient.ct` | Conjugate gradients against an operator, plain and preconditioned |
| `src/sparse/csc.ct` | Owned sparse matrices and operations |
| `src/sparse/gram.ct` | Fixed-pattern AᵀA and AAᵀ products |
| `src/sparse/augmented.ct` | Fixed-pattern augmented matrix of a sparse A |
| `src/sparse/scaling.ct` | Exact scaling of rows, columns and symmetric diagonals by powers of two |
| `src/sparse/ordering.ct` | Private minimum-degree ordering, with held indices |
| `src/sparse/limits.ct`, `src/dense/limits.ct` | Private dimension bounds and entry budgets |
| `src/sparse/spd-ldlt.ct` | SPD, signed-pivot and pivot-skipping factorization, solves and condition estimate |
| `src/sparse/dependent-rows.ct` | Dependent rows of a sparse matrix, their null vectors and proved bounds |
| `src/sparse/lu.ct` | Sparse LU with partial pivoting, and its solves |
| `src/dense/indefinite-ldlt.ct` | Dense symmetric indefinite factorization and inertia |
| `src/dense/orthogonal.ct` | Dense pivoted QR, complete orthogonal decomposition, least squares |
| `src/dense/jacobi.ct` | Dense symmetric eigenvalues and singular values by Jacobi rotations |
| `tests/` | Interface/runtime tests and small generic fixtures |
| `tools/` | Optional artifact drivers and independent Python verifiers |
| `examples/` | Executable caller examples, exercised by tests |
| `scripts/` | Configurable serial checks and inspection |
| `docs/` | Contracts, derivation, provenance, verification, and design sources |
| `toolchain.lock.json` | Installed dependency content and source-resolution pin |
| `toolchain.sources.json` | Versioned public acquisition coordinates and archive hashes |

The original prototype and its CAD/held-out results remain separate and unchanged.
This repository carries reusable code and generic tests, not those application
fixtures or historical benchmark dumps.

Licensed under [MIT](LICENSE). [Provenance](docs/provenance.md) records the
independent origin of the implementation. Dependency licenses remain their own;
dependency implementations are not bundled in this repository.