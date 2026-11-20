# coalton-linalg

An independent numerical linear algebra library in Coalton, initially targeting
64-bit SBCL. Its first implemented module is sparse, symmetric positive-definite
LDLᵀ factorization with reusable symbolic analysis. It also provides
validated rectangular sparse matrices, matrix operations, and stable vector norms.

The library accepts matrices and vectors. Applications own geometry, nonlinear
iteration, rank policy, units, and integration. No Blender or CAD dependency is
loaded or required.

## Current scope

- Owned rectangular CSC/COO construction with explicit duplicate handling.
- Transpose, row/column permutations, A*x and transpose(A)*x.
- Stable Euclidean/infinity norms using Coalton standard numeric primitives.
- Lower compressed sparse column (CSC) factorization input, double-precision values.
- Deterministic exact minimum-degree ordering and symbolic factor structure.
- Factorization of `P (A - shift I) Pᵀ = L D Lᵀ`, without numerical pivoting.
- Repeated factorization and solves using owned reusable workspace.
- Explicit refusals, borrowed-factor ownership, and independent numerical checks.

This is an early library, not a general indefinite solver. Symbolic analysis uses
compact adjacency tables and an indexed minimum-degree heap. Storage follows the
input and factor fill; fill can itself be dense. Analysis defaults to a budget of
1,048,576 total factor entries (including the implicit diagonal), with explicit
entry points for a caller-selected budget. See [factorization contract](docs/spd-ldlt.md),
[primitive contracts](docs/primitives.md), and [roadmap](ROADMAP.md).

## Review this snapshot

Start with the [setup recipe](docs/setup.md), [public interface](docs/api.md),
and the two executable examples below. The mathematical implementation lives under src/: numeric primitives,
rectangular CSC operations, sparse ordering, and SPD LDLᵀ.

Review priorities are Coalton idioms and type design; private handles and borrowed
array semantics; and the representation/calling choices behind repeated
factor!/solve! calls. The [design review](docs/design.md#coalton-interface-review-2026-11-19)
records the earlier trade-offs; the [pre-publication revision](docs/revision-2.md)
settles the factor interface. [Batch 3a](docs/batch3a.md) preserves its numerical
results and adds compact analysis with explicit resource budgets. F64 remains the only numeric type.

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
`coalton-linalg/sparse/csc`, and `coalton-linalg/sparse/spd-ldlt`. The [interface reference](docs/api.md) indexes
all exported symbols and links to their contracts. [A compiled SPD example](examples/spd-solve.ct)
uses the local nickname `ldlt`; [a sparse example](examples/sparse-operations.ct)
shows construction and products. The [mathematical glossary](CONTEXT.md) explains
shared terms.

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

Run these **serially**. They cover interface and allocation regressions, a Coalton
caller and examples, fixture tooling, checker controls, all 20 original matrices,
336 seeded matrices, 96 repeated refactors, independent primitive tests, and 161
structural cases (13 above the former 4096 cap).
Generated artifacts stay under `build/`. The scripts rebuild the library for
each safety policy. See the [compact-analysis evidence](docs/batch3a.md) for the latest full runs, the
[initial review evidence](docs/initial-review.md) for fresh-source setup verification, and
[pinned toolchain details](docs/toolchain.md) for the measurement protocol.

## Repository map

| Path | Purpose |
| --- | --- |
| `src/errors.ct` | Typed validation and primitive arithmetic errors |
| `src/numeric.ct` | Stable norms |
| `src/sparse/csc.ct` | Owned sparse matrices and operations |
| `src/sparse/ordering.ct` | Private compact graph ordering |
| `src/sparse/spd-ldlt.ct` | SPD factorization and solves |
| `tests/` | Interface/runtime tests and small generic fixtures |
| `tools/` | Optional artifact driver and independent Python verifier |
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