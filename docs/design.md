# Layout decisions and upstream guidance

Consulted 2026-11-19. These are documentation sources, not sources of numerical
implementation code.

| Source | Guidance used here |
| --- | --- |
| [Coalton: Introduction and Getting Started](https://coalton-lang.github.io/manual/topics/introduction/) | Use the normal Common Lisp toolchain, an ASDF system, coalton-asdf, and .ct source components. |
| [ASDF: Predefined operations](https://asdf.common-lisp.dev/asdf/Predefined-operations-of-ASDF.html) | Put tests in a separate system and attach them to test-op; signal failures rather than merely print them. |
| [ASDF: Components](https://asdf.common-lisp.dev/asdf/Components.html) | Name secondary systems primary/name so ASDF can discover them from the primary .asd file. |
| [Eigen: Solving Sparse Linear Systems](https://libeigen.gitlab.io/eigen/docs-3.4/group__TopicSparseSystems.html) | Keep structural analysis separate from numeric factorization and repeated solves when the pattern is unchanged. |

Our layout choices follow from those conventions and this library's needs:
`src/` contains the mathematical implementation; `tests/` checks its interface;
`tools/` produces and independently verifies artifacts; `examples/` demonstrates
a caller. Loading the primary system has no file-output or test side effects.

The first factorization remains one coherent module. Splitting storage, ordering,
and every helper into public packages now would enlarge the interface before
another implemented operation needs them. Extract shared modules when actual
callers establish a shared contract.

The main interface includes storage ordering, mutation, lifetime, concurrency,
numerical assumptions, failure semantics, and complexity—not just function types.
The independent verifier crosses the artifact interface and imports no Coalton
implementation. No application-specific fixture name changes an acceptance rule.

The current SBCL implementation assumptions are explicit. Configurable paths and
a fallback timer are useful groundwork, but do not establish cross-platform
numerical or performance support. A CI configuration will follow a reproducible,
fully specified dependency setup; no untested workflow is presented as working CI.

## Coalton interface review, 2026-11-19

The constructor/status recommendations below were subsequently implemented for
the initial review snapshot; they describe the reasoning for that change.
The current signatures and guarantees are in [the interface reference](api.md).
No broader numeric/storage generalization was introduced.

Scope: source and contract review of the three mathematical modules, with an
independent reviewing agent. This review changes documentation only. It does
not rerun numerical tests or establish new speed/allocation measurements;
[batch results](batches-1-2.md) remain the latest measured evidence.

### Strengths already used

Coalton offers first-class functions, static types, algebraic data types, and
type classes; the host Lisp supplies interactive development and macros.
These are tools for expressing and checking abstractions, not evidence that
an arbitrary function call is faster than a C++ or Rust call.
See the [language guide](https://coalton-lang.github.io/manual/topics/language-guide/).

The code already uses typed functions, a private Matrix constructor, an explicit
DuplicatePolicy type, shared finite checks, and a checked CSC bridge that reuses
symbolic analysis. The analyse/factor/solve split separates reusable structure
from changing values. These choices are worth retaining.

Typed mutable arrays and local accumulation are appropriate here. The pinned
compiler documents elimination of non-escaping cells in release mode, and
different calling conventions for direct and higher-order calls. Replacing
loops with lists, callbacks, or closures is not inherently an improvement.
The existing dynamic-shift caller and identity/allocation checks are especially
important when changing a function's ABI.
See [performance guidance](https://coalton-lang.github.io/manual/topics/performance-and-optimization/)
and [allocation evidence](runtime-allocation.md).

### Recommended interface work

1. **Hide workspace construction before stabilizing the interface.**
   Symbolic and NumericFactor use exported define-struct names, exposing both
   constructors and positional destructuring of workspace/state arrays.
   Matrix already uses a separate private constructor. Give the workspace
   types the same separation, while retaining Pattern as an input descriptor
   checked by analyse. This is an encapsulation gap, not a newly observed
   numerical failure. Validate with a real external Coalton caller that cannot
   fabricate/destructure the workspace, plus existing lifecycle/identity tests.
   [Coalton's struct definition](https://coalton-lang.github.io/manual/operators/define-struct/)
   confirms the type and constructor share a name.

2. **Use distinct named factor and solve statuses.**
   The original result-code and solve returned UFix with different meanings for the
   same numbers. Typed FactorStatus/SolveStatus values would express those
   choices and support pattern matching. They can be introduced without
   exposing ResultStorage or replacing the cached FactorResult representation.
   Decide deliberately whether numeric access remains for host interoperability;
   do not grow a second status mapping independently. Check all outcomes,
   external typed callers, cached identity, and allocations at both safety
   policies before changing this interface. No zero-cost claim is made for
   this proposed design.

3. **Keep ownership limits visible.**
   Hiding constructors does not make returned LispArray values read-only.
   Inspection is borrowed by contract, and refactorization changes an existing
   factor handle. Keep those obligations next to the functions; a future checked
   host-language adapter or copy-out operation should answer a concrete caller
   need. Do not claim compile-time lifetime protection.

Malformed-construction errors currently use ordinary Lisp errors with text.
Give them structured condition categories if a consumer needs programmatic
recovery; do not encourage matching error strings. Revisit mutating names
(factor/solve versus names ending in !) when settling the interface, not as an
isolated rename.

### Abstractions to defer

Do not generalize F64 to an arbitrary Num constraint yet: finiteness, square
roots, rounding, and positive-pivot behavior require more than arithmetic
operators. Add a type class when a second implemented numeric/storage type
establishes its laws. Coalton's [specialization mechanism](https://coalton-lang.github.io/manual/operators/specialize/)
can optimize known concrete types, but optimization must not change behavior.

Do not add macros merely to shorten loops. A macro earns its place when it
centralizes an actual repeated language pattern. During Batch 3, private
functions separating graph construction, ordering, and scatter preparation may
improve analyse's readability; they do not need to become public packages.
The Boolean-directed product helper is a measurement candidate, not a proven
performance problem.

### Documentation roles

The [mathematical glossary](../CONTEXT.md) defines shared terms. The
[public symbol index](api.md) locates types and functions. The existing per-module
contracts remain authoritative for behavior, and executable examples show valid
calling sequences. Keep each fact in its owning document and link to it.

This division follows established numerical-library documentation practice:
[LAPACK documents conventions and routine specifications](https://netlib.org/lapack/lug/node15.html);
[SciPy's sparse guide explains storage concepts](https://docs.scipy.org/doc/scipy/tutorial/sparse.html)
alongside its separate reference. There is no required glossary filename or
reason to build a documentation generator for this small interface yet.
Only documentation was consulted, not numerical implementations.

Language details were checked against the local content-pinned Coalton source
and manual at b6e284bdd91ef9ba0fb303e0798c9818b9b3e268. Online manual pages may
advance beyond that revision.
