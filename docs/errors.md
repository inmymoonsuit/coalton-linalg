# Typed errors

Package: `coalton-linalg/errors`, source: [errors.ct](../src/errors.ct).

`LinalgError` is a Coalton exception family. Catch a constructor for a particular
failure, or catch the family and match its constructors. All cases are Common Lisp
`error` conditions as well. Message strings supply human context and are not stable
codes; callers should match constructors, not text.

| Constructor | Meaning and payload |
| --- | --- |
| `InvalidStructure(message)` | Malformed CSC/COO arrays, invalid indices or ordering, or invalid lower/diagonal SPD structure. |
| `InvalidDimensions(message)` | A matrix/vector size mismatch or nonsquare SPD matrix. |
| `InvalidPermutation(message)` | Wrong permutation length, out-of-range or repeated index. |
| `DuplicateEntry(row,column)` | A repeated coordinate in canonical CSC or COO with RejectDuplicates; zero-based caller coordinates. |
| `NonFiniteInput(message)` | Nonfinite entries passed to matrix construction, products or norms. |
| `NonFiniteResult(message)` | Nonfinite primitive result or duplicate accumulation. |
| `AnalysisLimitExceeded(requested,maximum)` | Analysis dimension exceeds its representability bound; both payloads UFix. |
| `AnalysisStorageExceeded(required,budget)` | Analysis total factor entries, including the diagonal, exceed the entry budget. Required is an Integer lower bound; budget is the effective UFix budget (capped for representability). |

SPD Pattern validation reports InvalidStructure for repeated rows as well as other
structural violations. General CSC/COO constructors specifically identify repeated
coordinates as DuplicateEntry. Where several checks fail, callers should not assume
an ordering except the documented SPD dimension/entry-budget validation order. No partially constructed
matrix or analysis is returned on an exception.

Factorization and solve have their own FactorStatus and SolveStatus values. Their
expected numerical/argument refusals remain ordinary return values and do not throw
LinalgError. A saved status describes one call; current factor readiness and the
latest failure column have the lifetimes in [the factorization contract](spd-ldlt.md).

Raw Lisp entry points still require the declared argument types and specialized
arrays. This exception family does not turn them into a dynamically checked FFI.
Floating-point traps must have the supported disabled policy; enabling traps can
produce host arithmetic exceptions before explicit finiteness checks.
