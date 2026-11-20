# Public interface reference

This is the symbol index for the six public library packages, excluding examples,
tooling and the private sparse/ordering package. The linked contracts specify storage, errors, ownership, numerical
assumptions, and cost. See the [glossary](../CONTEXT.md) for mathematical terms.

Below, `V` abbreviates `(coalton/lisparray:LispArray F64)`, and `I` abbreviates
`(coalton/lisparray:LispArray UFix)`. These are documentation abbreviations,
not exported type aliases. `*` separates arguments; `->` introduces the result.
Dimensions and indices use UFix; matrix indices are zero-based. Factor and solve
outcomes have distinct named types.

## coalton-linalg/errors

[Typed error contract](errors.md) · [Source](../src/errors.ct)

| Export | Type or role |
| --- | --- |
| `LinalgError` | Exception family for construction, validation and primitive arithmetic |
| `InvalidStructure`, `InvalidDimensions`, `InvalidPermutation`, `InvalidAliasing` | String -> LinalgError |
| `DuplicateEntry` | UFix * UFix -> LinalgError; row, column |
| `NonFiniteInput`, `NonFiniteResult` | String -> LinalgError |
| `AnalysisLimitExceeded` | UFix * UFix -> LinalgError; requested dimension, representability maximum |
| `AnalysisStorageExceeded` | Integer * UFix -> LinalgError; required analyzed output entries (factor or Gram), budget |

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
| `matvec!`, `transpose-matvec!` | Matrix * V * V -> Void; caller-provided output |

Constructor argument order is in the contract. Inspection arrays are borrowed
read-only views by agreement; they are ordinary mutable arrays at runtime.
Constructors and transforming operations return owned storage.

## coalton-linalg/sparse/gram

[Contract](gram.md) · [Source](../src/sparse/gram.ct)

| Export | Type or role |
| --- | --- |
| `GramPlan` | Fixed-pattern product plan with private workspace; obtain from the analyze operations |
| `analyze-ata`, `analyze-aat` | coalton-linalg/sparse/csc:Matrix -> GramPlan; lower A^T A or A A^T |
| `analyze-ata-with-limit`, `analyze-aat-with-limit` | Matrix * UFix -> GramPlan; explicit stored-entry budget |
| `dimension`, `entry-count` | GramPlan -> UFix |
| `column-pointers`, `row-indices` | GramPlan -> I; lower CSC pattern with every diagonal |
| `refill!` | GramPlan * V * V -> Void; source values in analyzed order, caller output |
| `refill-csc!` | GramPlan * Matrix * V -> Void; checks the source pattern first |
| `maximum-analysis-size`, `default-analysis-entry-limit` | UFix; the bounds shared with SPD analysis |

A plan owns scratch, so refills sharing one plan are not reentrant. The output
values belong to the caller and feed `factor!` with a `Pattern` built from the
plan's dimension, column pointers and row indices.

## coalton-linalg/sparse/spd-ldlt

[Contract and outcome meanings](spd-ldlt.md) · [Condition estimate](condition.md) · [Source](../src/sparse/spd-ldlt.ct)

| Export | Type or role |
| --- | --- |
| `maximum-analysis-size` | UFix; representability bound for dimensions and effective entry budgets |
| `default-analysis-entry-limit` | UFix; 1,048,576 total factor entries including the diagonal |
| `Pattern` | Lower CSC input; constructor UFix * I * I -> Pattern |
| `Symbolic` | Analysis and reusable workspace; obtain from analyze/analyze-csc |
| `NumericFactor` | Borrowed factor; obtain from symbolic-factor |
| `analyze` | Pattern -> Symbolic |
| `analyze-csc` | coalton-linalg/sparse/csc:Matrix -> Symbolic |
| `analyze-with-limit` | Pattern * UFix -> Symbolic |
| `analyze-csc-with-limit` | coalton-linalg/sparse/csc:Matrix * UFix -> Symbolic |
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
| `ConditionStatus` | Condition-estimate outcome type |
| `ConditionSuccess`, `ConditionInvalidFactor`, `ConditionNonFiniteResult` | ConditionStatus values |
| `estimate-inverse-norm1!` | NumericFactor -> ConditionStatus; lower estimate of the inverse's 1-norm |
| `inverse-norm1-estimate` | NumericFactor -> F64; latest estimate for the standing factor, otherwise -1 |
| `estimate-solve-count` | NumericFactor -> UFix; solves made by the latest estimate call |
| `maximum-estimate-solves` | UFix; 11, the bound on that count |
| `shifted-norm1` | Symbolic * V * F64 -> F64; 1-norm of A - shift I from the analyzed lower values |
| `reciprocal-condition1` | F64 * F64 -> F64; 1 / (norm * inverse-norm), or 0 when that is not positive and finite |

Every factor call overwrites or invalidates the owning workspace's previous
factor. Solve mutates scratch and, on success, the supplied output. Handles and
inspection arrays are borrowed; calls sharing a workspace are not concurrent.

Symbolic and NumericFactor have private constructors and hide their workspace
fields from ordinary external Coalton callers. The public inspection arrays
remain mutable at runtime and read-only by contract. An earlier FactorSuccess
describes its original call; the separately borrowed factor can later become invalid.
See the [design review](design.md#coalton-interface-review-2026-11-19).

## coalton-linalg/dense/indefinite-ldlt

[Contract and derivation](indefinite-ldlt.md) · [Source](../src/dense/indefinite-ldlt.ct)

A dense symmetric matrix is a V of n² elements in column-major order, entry (i,j) at
i + j n, of which only the lower triangle is read.

| Export | Type or role |
| --- | --- |
| `Factor` | Reusable dense factor with private workspace; obtain from make-factor |
| `make-factor` | UFix -> Factor; storage for n-by-n problems under default-entry-limit |
| `make-factor-with-limit` | UFix * UFix -> Factor; explicit budget of n² entries |
| `maximum-dimension`, `default-entry-limit` | UFix; 32,768 and 16,777,216 |
| `pivot-threshold` | F64; (1 + √17)/8, the 1-by-1 acceptance ratio |
| `dimension` | Factor -> UFix |
| `FactorStatus` | Factorization outcome type |
| `FactorSuccess`, `FactorNonFinite`, `FactorBadSize`, `FactorBadShift` | FactorStatus values |
| `SolveStatus` | Solve outcome type |
| `SolveSuccess`, `SolveInvalidFactor`, `SolveBadSize`, `SolveSingular`, `SolveNonFiniteRhs`, `SolveNonFiniteResult` | SolveStatus values |
| `factor!` | Factor * V * F64 -> FactorStatus; P (A - shift I) Pᵀ = L D Lᵀ |
| `solve!` | Factor * V * V -> SolveStatus |
| `factor-ready?` | Factor -> Boolean |
| `positive-count`, `negative-count`, `zero-count` | Factor -> UFix; inertia of the computed block diagonal |
| `factor-permutation`, `factor-block-sizes` | Factor -> I |
| `factor-lower`, `factor-diagonal`, `factor-subdiagonal` | Factor -> V |

`FactorStatus` and `SolveStatus` here are distinct types from the sparse package's. A
singular matrix is `FactorSuccess` with a nonzero `zero-count`. Inspection arrays are
borrowed, read-only by contract, and rewritten by the next `factor!`.

## Choosing the entry point

Use [the SPD example](../examples/spd-solve.ct) for analyze/factor!/solve! and
[the sparse example](../examples/sparse-operations.ct) for matrix construction
and products. [The normal-equations example](../examples/normal-equations.ct) joins Gram
products to the SPD factorization, and [the minimum-norm example](../examples/minimum-norm.ct)
uses the dense indefinite factorization. All are compiled and exercised by the library tests.

A general CSC matrix represents exactly its stored entries; it has no implicit
symmetry. The SPD interface specifically interprets its lower triangle as a
symmetric matrix. The two contracts must not be interchanged.

Coalton checks declared argument types, not matrix dimensions, finiteness, array
ownership, or future handle validity. Those properties have runtime checks or
calling obligations in the contracts. Raw Common Lisp entry points require the
documented specialized arrays and are not a dynamically checked external API.