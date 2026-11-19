# Public interface reference

This is the symbol index for the four library packages, excluding examples and
tooling. The linked contracts specify storage, errors, ownership, numerical
assumptions, and cost. See the [glossary](../CONTEXT.md) for mathematical terms.

Below, `V` abbreviates `(coalton/lisparray:LispArray F64)`, and `I` abbreviates
`(coalton/lisparray:LispArray UFix)`. These are documentation abbreviations,
not exported type aliases. `*` separates arguments; `->` introduces the result.
Dimensions and indices use UFix; matrix indices are zero-based. Factor and solve
outcomes have distinct named types.

## coalton-linalg/errors

[Typed error contract](errors.md) � [Source](../src/errors.ct)

| Export | Type or role |
| --- | --- |
| `LinalgError` | Exception family for construction, validation and primitive arithmetic |
| `InvalidStructure`, `InvalidDimensions`, `InvalidPermutation` | String -> LinalgError |
| `DuplicateEntry` | UFix * UFix -> LinalgError; row, column |
| `NonFiniteInput`, `NonFiniteResult` | String -> LinalgError |
| `AnalysisLimitExceeded` | UFix * UFix -> LinalgError; requested dimension, maximum |

## coalton-linalg/numeric

[Contract and derivation](primitives.md#numeric) · [Source](../src/numeric.ct)

| Export | Type |
| --- | --- |
| `norm2`, `norm-inf` | V -> F64 |

Use `coalton/math:finite?` for finiteness; this library does not duplicate that export.

## coalton-linalg/sparse/csc

[Contract](primitives.md#canonical-rectangular-csc) · [Source](../src/sparse/csc.ct)

| Export | Type or role |
| --- | --- |
| `Matrix` | Validated matrix type; construction is through from-csc/from-coo |
| `DuplicatePolicy` | Choice of duplicate-coordinate handling |
| `RejectDuplicates`, `SumDuplicates` | DuplicatePolicy values |
| `from-csc` | UFix * UFix * I * I * V -> Matrix |
| `from-coo` | UFix * UFix * I * I * V * DuplicatePolicy -> Matrix |
| `row-count`, `column-count`, `entry-count` | Matrix -> UFix |
| `column-pointers`, `row-indices` | Matrix -> I |
| `entries` | Matrix -> V |
| `transpose` | Matrix -> Matrix |
| `permute` | Matrix * I * I -> Matrix |
| `matvec`, `transpose-matvec` | Matrix * V -> V |

Constructor argument order is in the contract. Inspection arrays are borrowed
read-only views by agreement; they are ordinary mutable arrays at runtime.
Constructors and transforming operations return owned storage.

## coalton-linalg/sparse/spd-ldlt

[Contract and outcome meanings](spd-ldlt.md) · [Source](../src/sparse/spd-ldlt.ct)

| Export | Type or role |
| --- | --- |
| `maximum-analysis-size` | UFix; current analysis dimension cap, 4096 |
| `Pattern` | Lower CSC input; constructor UFix * I * I -> Pattern |
| `Symbolic` | Analysis and reusable workspace; obtain from analyze/analyze-csc |
| `NumericFactor` | Borrowed factor; obtain from symbolic-factor |
| `analyze` | Pattern -> Symbolic |
| `analyze-csc` | coalton-linalg/sparse/csc:Matrix -> Symbolic |
| `FactorStatus` | Factorization outcome type |
| `FactorSuccess`, `FactorNonPositive`, `FactorNonFinite`, `FactorBadValues`, `FactorBadShift` | FactorStatus values |
| `SolveStatus` | Solve outcome type |
| `SolveSuccess`, `SolveInvalidFactor`, `SolveBadSize`, `SolveNonFiniteRhs`, `SolveNonFiniteResult` | SolveStatus values |
| `factor!` | Symbolic * V * F64 -> FactorStatus |
| `refused-column` | Symbolic -> IFix; latest factor! column, otherwise -1 |
| `symbolic-factor` | Symbolic -> NumericFactor |
| `factor-permutation`, `factor-columns`, `factor-rows` | NumericFactor -> I |
| `factor-values`, `factor-diagonal` | NumericFactor -> V |
| `factor-ready?` | NumericFactor -> Boolean |
| `solve!` | NumericFactor * V * V -> SolveStatus |

Every factor call overwrites or invalidates the owning workspace's previous
factor. Solve mutates scratch and, on success, the supplied output. Handles and
inspection arrays are borrowed; calls sharing a workspace are not concurrent.

Symbolic and NumericFactor have private constructors and hide their workspace
fields from ordinary external Coalton callers. The public inspection arrays
remain mutable at runtime and read-only by contract. An earlier FactorSuccess
describes its original call; the separately borrowed factor can later become invalid.
See the [design review](design.md#coalton-interface-review-2026-11-19).

## Choosing the entry point

Use [the SPD example](../examples/spd-solve.ct) for analyze/factor!/solve! and
[the sparse example](../examples/sparse-operations.ct) for matrix construction
and products. Both are compiled and exercised by the library tests.

A general CSC matrix represents exactly its stored entries; it has no implicit
symmetry. The SPD interface specifically interprets its lower triangle as a
symmetric matrix. The two contracts must not be interchanged.

Coalton checks declared argument types, not matrix dimensions, finiteness, array
ownership, or future handle validity. Those properties have runtime checks or
calling obligations in the contracts. Raw Common Lisp entry points require the
documented specialized arrays and are not a dynamically checked external API.