# Public interface reference

This is the symbol index for the eight public library packages, excluding examples,
tooling and the private sparse/ordering, sparse/limits and dense/limits packages. The linked contracts specify storage, errors, ownership, numerical
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

## coalton-linalg/refinement

[Contract](refinement.md) · [Source](../src/refinement.ct)

| Export | Type or role |
| --- | --- |
| `Refiner` | Workspace and latest-call history; obtain from make-refiner |
| `make-refiner` | UFix * UFix -> Refiner; dimension, step limit |
| `dimension`, `step-limit` | Refiner -> UFix |
| `RefineStatus` | Outcome type of refine! |
| `RefineConverged`, `RefineStalled`, `RefineLimit`, `RefineNonFinite`, `RefineSolveRefused`, `RefineBadSize`, `RefineBadTolerance` | RefineStatus values |
| `refine!` | Refiner * (V * V -> Void) * (V * V -> Boolean) * V * V * F64 * F64 -> RefineStatus; apply, solve, right-hand side, solution, tolerance, stall factor |
| `step-count`, `best-step`, `residual-count` | Refiner -> UFix; of the latest call |
| `residual-norms` | Refiner -> V; infinity norms of the latest call's residuals |

The solution array holds the best iterate seen on every return. It may also be the
right-hand side on entry. A call allocates nothing.

## coalton-linalg/conjugate-gradient

[Contract](conjugate-gradient.md) · [Source](../src/conjugate-gradient.ct)

| Export | Type or role |
| --- | --- |
| `CgWorkspace` | Workspace and latest-call history; obtain from make-cg |
| `make-cg` | UFix * UFix -> CgWorkspace; dimension, step limit |
| `dimension`, `step-limit` | CgWorkspace -> UFix |
| `CgStatus` | Outcome type of cg! and pcg! |
| `CgConverged`, `CgLimit`, `CgNotPositive`, `CgNonFinite`, `CgPreconditionerRefused`, `CgBadSize`, `CgBadTolerance` | CgStatus values |
| `cg!` | CgWorkspace * (V * V -> Void) * V * V * F64 -> CgStatus; apply, right-hand side, solution, tolerance |
| `pcg!` | CgWorkspace * (V * V -> Void) * (V * V -> Boolean) * V * V * F64 -> CgStatus; apply, precondition, right-hand side, solution, tolerance |
| `step-count`, `restart-count` | CgWorkspace -> UFix; of the latest call |
| `residual-norm` | CgWorkspace -> F64; 2-norm of b − A x for the vector returned, formed from the operator, otherwise -1 |
| `recurrence-norms` | CgWorkspace -> V; 2-norms of the latest call's running residual |

The operator must be symmetric positive definite and the tolerance is absolute, in the
2-norm. The right-hand side and solution may share one array. A call allocates nothing.

## coalton-linalg/lsqr

[Contract](lsqr.md) · [Source](../src/lsqr.ct)

| Export | Type or role |
| --- | --- |
| `LsqrWorkspace` | Workspace and latest-call history; obtain from make-lsqr |
| `make-lsqr` | UFix * UFix * UFix -> LsqrWorkspace; rows, columns, step limit |
| `rows`, `columns`, `step-limit` | LsqrWorkspace -> UFix |
| `LsqrStatus` | Outcome type of lsqr! |
| `LsqrResidualMet`, `LsqrGradientMet`, `LsqrExhausted`, `LsqrLimit`, `LsqrNonFinite`, `LsqrOutOfRange`, `LsqrBadSize`, `LsqrBadTolerance`, `LsqrBadDamping` | LsqrStatus values |
| `lsqr!` | LsqrWorkspace * (V * V -> Void) * (V * V -> Void) * V * V * F64 * F64 * F64 -> LsqrStatus; apply, apply-transpose, right-hand side, solution, damping, residual tolerance, gradient tolerance |
| `step-count`, `confirmation-count` | LsqrWorkspace -> UFix; of the latest call |
| `residual-norm`, `gradient-norm` | LsqrWorkspace -> F64; 2-norms of b − A x and of Aᵀ(b − A x) − damping² x for the vector returned, formed from the operators, otherwise -1 |
| `recurrence-residual-norms`, `recurrence-gradient-norms` | LsqrWorkspace -> V; the latest call's running estimates of the two |
| `operator-norm-estimate` | LsqrWorkspace -> F64; the Frobenius norm of the bidiagonal matrix reached, an estimate only; -1 where no product was made |

The operator may be of any shape and rank. The solution array is written only: the method
starts from zero, and without damping returns the minimizer of least norm. The tolerances
are absolute, in the 2-norm. The right-hand side and solution may share one array. A call
allocates nothing.

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
| `check-csc!` | GramPlan * Matrix -> Void; the source pattern check alone |
| `maximum-analysis-size`, `default-analysis-entry-limit` | UFix; the bounds shared with SPD analysis |

A plan owns scratch, so refills sharing one plan are not reentrant. The output
values belong to the caller and feed `factor!` with a `Pattern` built from the
plan's dimension, column pointers and row indices.

## coalton-linalg/sparse/spd-ldlt

[Contract and outcome meanings](spd-ldlt.md) · [Condition estimate](condition.md) · [Cholesky form](cholesky-form.md) · [Signed pivots](quasidefinite.md) · [Skipped pivots](skipping.md) · [Source](../src/sparse/spd-ldlt.ct)

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
| `factor-cholesky!` | Symbolic * V * F64 -> FactorStatus; Cholesky-form factor of the shifted matrix, for a proof |
| `cholesky-ready?` | NumericFactor -> Boolean |
| `solve-cholesky!` | NumericFactor * V * V -> SolveStatus; solves (A − shift I) x = b with the Cholesky-form factor |
| `analyze-quasidefinite` | Pattern * UFix -> Symbolic; indices from the second argument on are held back |
| `analyze-quasidefinite-with-limit` | Pattern * UFix * UFix -> Symbolic; the last argument is the entry budget |
| `SignedFactorStatus` | Signed factorization outcome type |
| `SignedSuccess`, `SignedZeroPivot`, `SignedNonFinite`, `SignedBadValues`, `SignedBadShift` | SignedFactorStatus values |
| `factor-signed!` | Symbolic * V * F64 -> SignedFactorStatus; pivots of either sign, no numerical pivoting |
| `factor-growth` | NumericFactor -> F64; measured by the latest successful factor-signed!, otherwise -1 |
| `factor-backward-bound` | NumericFactor -> F64; upper bound on the Frobenius norm of P K Pᵀ − L D Lᵀ for that factor, otherwise -1 |
| `positive-pivot-count`, `negative-pivot-count` | NumericFactor -> UFix; signs of D in the standing factor |
| `SkipFactorStatus` | Pivot-skipping factorization outcome type |
| `SkipSuccess`, `SkipNonFinite`, `SkipBadValues`, `SkipBadShift`, `SkipBadFloor` | SkipFactorStatus values |
| `factor-skipping!` | Symbolic * V * F64 * F64 -> SkipFactorStatus; shift, then the pivot floor: an index whose pivot is at or below it is left out |
| `skipping-ready?` | NumericFactor -> Boolean |
| `skipped-count` | NumericFactor -> UFix; indices the standing pivot-skipping factor left out |
| `solve-skipping!` | NumericFactor * V * V -> SolveStatus; solves with the kept principal submatrix, zero at the indices left out |
| `NullVectorStatus` | Outcome type of skipped-null-vector! |
| `NullVectorSuccess`, `NullVectorNoFactor`, `NullVectorNotSkipped`, `NullVectorBadSize`, `NullVectorNonFinite` | NullVectorStatus values |
| `skipped-null-vector!` | Symbolic * UFix * V -> NullVectorStatus; for an index left out, its combination of the kept indices before it |
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

## coalton-linalg/sparse/damped

[Contract](damped.md) · [Source](../src/sparse/damped.ct)

| Export | Type or role |
| --- | --- |
| `DampedPlan` | Analysis, factor and workspace for the damped normal equations of one pattern; obtain from analyze |
| `analyze` | coalton-linalg/sparse/csc:Matrix -> DampedPlan |
| `analyze-with-limit` | coalton-linalg/sparse/csc:Matrix * UFix -> DampedPlan; one entry budget for Jᵀ J and for its factor |
| `default-analysis-entry-limit` | UFix; 1,048,576 |
| `row-count`, `column-count` | DampedPlan -> UFix |
| `set-matrix!` | DampedPlan * coalton-linalg/sparse/csc:Matrix -> Void; the values of J for the analyzed pattern |
| `matrix-set?`, `damping-set?`, `factor-ready?` | DampedPlan -> Boolean |
| `DampingStatus` | Outcome type of set-damping! and set-column-damping! |
| `DampingSet`, `DampingBadSize`, `DampingNotPositive`, `DampingNoMatrix` | DampingStatus values |
| `set-damping!` | DampedPlan * V -> DampingStatus; a positive diagonal D in the units of the matrix squared |
| `set-column-damping!` | DampedPlan * F64 -> DampingStatus; D = diag(Jᵀ J) + floor, following the matrix |
| `DampedStatus` | Outcome type of factor! |
| `DampedSuccess`, `DampedNotPositive`, `DampedNonFinite`, `DampedBadDamping`, `DampedNotReady` | DampedStatus values |
| `factor!` | DampedPlan * F64 -> DampedStatus; the factor of Jᵀ J + mu D |
| `refused-column` | DampedPlan -> IFix; the column whose pivot the latest factor! refused, or -1 |
| `StepStatus` | Outcome type of step! |
| `StepSuccess`, `StepNoFactor`, `StepBadSize`, `StepNonFiniteRhs`, `StepNonFiniteResult` | StepStatus values |
| `step!` | DampedPlan * V * V -> StepStatus; the step -(Jᵀ J + mu D)⁻¹ Jᵀ r |

The plan's constructor and storage are private. One plan serves one call at a time.

## coalton-linalg/sparse/dependent-rows

[Contract](dependent-rows.md) · [Source](../src/sparse/dependent-rows.ct)

| Export | Type or role |
| --- | --- |
| `RowAnalysis` | Analysis, factor and workspace for the rows of one pattern; obtain from analyze |
| `analyze` | coalton-linalg/sparse/csc:Matrix -> RowAnalysis |
| `analyze-with-limit` | coalton-linalg/sparse/csc:Matrix * UFix -> RowAnalysis; one entry budget for A Aᵀ and for its factor |
| `default-analysis-entry-limit` | UFix; 1,048,576 |
| `row-count`, `column-count` | RowAnalysis -> UFix |
| `RowStatus` | Outcome type of find-dependent-rows! |
| `RowsFound`, `RowsNonFinite`, `RowsUnderflow`, `RowsBadThreshold` | RowStatus values |
| `find-dependent-rows!` | RowAnalysis * coalton-linalg/sparse/csc:Matrix * F64 -> RowStatus; the threshold is a distance |
| `rows-found?` | RowAnalysis -> Boolean |
| `check-pattern!` | RowAnalysis * coalton-linalg/sparse/csc:Matrix -> Void; the typed error of the Gram refill unless the matrix has the analyzed pattern |
| `kept-count`, `left-out-count` | RowAnalysis -> UFix |
| `left-out-rows` | RowAnalysis -> I; borrowed, its first left-out-count entries are the rows left out |
| `kept-values!` | RowAnalysis * V -> Boolean; the matrix's values with zeros for the rows left out |
| `NullStatus` | Outcome type of null-vector! |
| `NullFound`, `NullNoRows`, `NullBadPosition`, `NullBadSize`, `NullNonFinite` | NullStatus values |
| `null-vector!` | RowAnalysis * UFix * V -> NullStatus; the left null vector of the row at a position of left-out-rows |
| `refinement-step-limit` | UFix; 4, the most corrections null-vector! applies |
| `residual-bound` | RowAnalysis * UFix -> F64; proved bound on the 2-norm of Aᵀ v for that vector, otherwise -1 |
| `gap-bound` | RowAnalysis -> F64; upper bound on the singular value after the kept count, otherwise -1 |
| `transpose-product-bound` | coalton-linalg/sparse/csc:Matrix * V -> F64; proved bound on the 2-norm of Aᵀ v |

A search withdraws the earlier set and its bounds. After the analysis only
`transpose-product-bound` allocates. One analysis serves one call at a time.

## coalton-linalg/sparse/rank-deficient

[Contract](rank-deficient.md) · [Source](../src/sparse/rank-deficient.ct)

| Export | Type or role |
| --- | --- |
| `KeptRows` | The kept rows of one pattern, the vectors of the rows left out, and what the certificate factors; obtain from analyze |
| `analyze` | coalton-linalg/sparse/csc:Matrix -> KeptRows |
| `analyze-with-limit` | coalton-linalg/sparse/csc:Matrix * UFix -> KeptRows; one entry budget for every product and factor |
| `default-analysis-entry-limit` | UFix; 1,048,576 |
| `row-count`, `column-count` | KeptRows -> UFix |
| `FindStatus` | Outcome type of find-kept-rows! |
| `FindFound`, `FindNonFinite`, `FindUnderflow`, `FindBadThreshold`, `FindBadDrop`, `FindVectorRefused` | FindStatus values |
| `find-kept-rows!` | KeptRows * coalton-linalg/sparse/csc:Matrix * F64 * F64 -> FindStatus; a distance threshold, then a drop tolerance relative to each vector's largest entry |
| `kept-found?` | KeptRows -> Boolean |
| `kept-count`, `left-out-count` | KeptRows -> UFix |
| `left-out-rows` | KeptRows -> I; borrowed, its first left-out-count entries are the rows left out |
| `multiplier-start`, `multiplier-count` | KeptRows * UFix -> UFix; where the kept vector at a position starts, and how many entries it holds at kept rows |
| `multiplier-rows` | KeptRows -> I; borrowed, the kept row of each entry of every kept vector |
| `multiplier-values` | KeptRows -> V; borrowed, the vector's value at each of those rows |
| `ReproveStatus` | Outcome type of reprove! |
| `Reproved`, `ReproveNothingKept`, `ReproveNonFinite` | ReproveStatus values |
| `reprove!` | KeptRows * coalton-linalg/sparse/csc:Matrix -> ReproveStatus; every kept vector's bound for other values of the pattern |
| `RefreshStatus` | Outcome type of refresh! |
| `Refreshed`, `RefreshNothingKept`, `RefreshBadTarget`, `RefreshNoShift`, `RefreshFactorRefused`, `RefreshNonFinite` | RefreshStatus values |
| `refresh!` | KeptRows * coalton-linalg/sparse/csc:Matrix * F64 -> RefreshStatus; reprove!, then each vector above its share of the gap target replaced from the values given |
| `refresh-to!` | KeptRows * coalton-linalg/sparse/csc:Matrix * F64 * F64 -> RefreshStatus; refresh! with the fraction of its share a replaced vector is corrected to, from 0 to 1 |
| `default-refresh-fraction` | F64; 0.0625 |
| `refreshed-count` | KeptRows -> UFix; vectors the latest refresh! replaced |
| `residual-bound` | KeptRows * UFix -> F64; proved bound on the 2-norm of Jᵀ v for the kept vector, otherwise -1 |
| `gap-bound` | KeptRows -> F64; upper bound on the singular value after the kept count, otherwise -1 |
| `CertificateStatus` | Outcome type of certify! |
| `Certified`, `CertificateGapAboveFloor`, `CertificateKeptNotAboveFloor`, `CertificateUnresolved`, `CertificateBadFloor`, `CertificateNothingKept` | CertificateStatus values |
| `certify!` | KeptRows * F64 * F64 -> CertificateStatus; a lower and an upper floor, in the matrix's units |
| `certified?` | KeptRows -> Boolean |
| `proved-lower-bound` | KeptRows -> F64; what the standing certificate proved every singular value of the kept rows exceeds, otherwise -1 |
| `largest-singular-bounds!` | coalton-linalg/sparse/csc:Matrix * UFix * V * V * V -> Void; a lower and an upper bound on the largest singular value, into the third array |
| `StepStatus` | Outcome type of step! |
| `StepTaken`, `StepNoCertificate`, `StepBadSize`, `StepNonFiniteRhs`, `StepFactorRefused`, `StepBoundNotReached`, `StepNonFiniteResult` | StepStatus values |
| `step!` | KeptRows * V * V -> StepStatus; the step −B⁺ r against a residual, for the certified matrix |
| `step-taken?` | KeptRows -> Boolean |
| `step-solve-bound`, `step-bound` | KeptRows -> F64; proved bounds on the standing step's distance from the exact step for B, and from the truncated matrix's step; otherwise -1 |
| `step-norm-term`, `step-outside-term` | KeptRows -> F64; the two terms step-bound adds to step-solve-bound, in the step's norm and in the residual outside the range; otherwise -1 |
| `step-corrections` | KeptRows -> UFix; corrections the standing step's refinement made |
| `system-iterations` | KeptRows -> UFix; preconditioned iterations the standing step made on the system of the rows left out |
| `system-factor-count` | KeptRows -> UFix; factorizations of that system since the search, the search's included |
| `system-iteration-limit` | UFix; 8 |
| `SupportStatus` | Outcome type of class-supports! |
| `SupportsClassed`, `SupportsNoCertificate`, `SupportsBadCut`, `SupportsUndecided` | SupportStatus values |
| `class-supports!` | KeptRows * F64 -> SupportStatus; a cut relative to each vector's largest entry |
| `supports-classed?` | KeptRows -> Boolean |
| `dependent-sets` | KeptRows -> I; borrowed, each row's label; the row count for a row in no dependency |
| `dependent-set-count` | KeptRows -> UFix |
| `support-flags` | KeptRows -> I; borrowed, beside multiplier-values: 1 classed in, 0 classed out |
| `multiplier-error-bound` | KeptRows * UFix -> F64; proved bound on the 2-norm of a vector's stored multipliers less the exact ones, otherwise -1 |
| `undecided-vector`, `undecided-entry` | KeptRows -> UFix; what the latest classing could not decide |
| `correction-limit` | UFix; 8 |
| `solve-growth`, `solve-backward-bound`, `predicted-rate` | KeptRows -> F64; diagnostics of the factor the steps solve with, otherwise -1 |

The handle's constructor and storage are private. A search replaces the kept set, the
vectors and the two borrowed vector arrays; a search or a reproof withdraws the
certificate, and each of them and a certificate withdraws the step; a refresh does both
and replaces the two borrowed vector arrays. `analyze` and `find-kept-rows!` allocate;
`refresh!` does while its arrays grow or where a vector changes the rows it holds,
`class-supports!` does at its first call and when the vectors have grown, and nothing
else does. One handle serves one call at a time.

## coalton-linalg/sparse/lu

[Contract](lu.md) · [Source](../src/sparse/lu.ct)

| Export | Type or role |
| --- | --- |
| `Lu` | Factors and workspace for square matrices of one order; obtain from make-lu |
| `make-lu` | coalton-linalg/sparse/csc:Matrix -> Lu; column-order of the pattern, default budget |
| `make-lu-with-order` | UFix * I * UFix -> Lu; order, the column eliminated at each step, entry budget |
| `column-order` | coalton-linalg/sparse/csc:Matrix * UFix -> I; minimum degree on the pattern of AᵀA |
| `default-entry-limit` | UFix; 1,048,576 entries of L and U together, off their diagonals |
| `LuStatus` | Outcome type of factor! |
| `LuSuccess`, `LuSingular`, `LuNonFinite`, `LuStorageExceeded`, `LuBadMatrix`, `LuBadThreshold` | LuStatus values |
| `factor!` | Lu * coalton-linalg/sparse/csc:Matrix * F64 -> LuStatus; P A Q = L U, the threshold in (0, 1] |
| `LuSolveStatus` | Outcome type of the solves |
| `LuSolveSuccess`, `LuSolveNoFactor`, `LuSolveBadSize`, `LuSolveNonFiniteRhs`, `LuSolveNonFiniteResult` | LuSolveStatus values |
| `solve!` | Lu * V * V -> LuSolveStatus; A x = b |
| `solve-transpose!` | Lu * V * V -> LuSolveStatus; Aᵀ x = b |
| `lu-ready?` | Lu -> Boolean |
| `dimension`, `entry-limit` | Lu -> UFix |
| `singular-column` | Lu -> UFix; the column at which the latest factor! found no pivot, otherwise the order |
| `lower-entry-count`, `upper-entry-count` | Lu -> UFix; entries of the standing factors off their diagonals |
| `lu-growth` | Lu -> F64; largest magnitude in U over the largest in A, otherwise -1 |
| `lu-backward-bound` | Lu -> F64; upper bound on the Frobenius norm of P A Q − L U, otherwise -1 |
| `row-permutation`, `column-permutation` | Lu -> I; the row and the column of the matrix taken at each step |
| `lower-pointers`, `lower-rows`, `upper-pointers`, `upper-rows` | Lu -> I; the factors by column |
| `lower-values`, `upper-values`, `upper-diagonal` | Lu -> V |

The factors' storage grows up to the budget and the row and value views are replaced
when it does. A refused factor! leaves no factors. One handle serves one call at a time.

## coalton-linalg/sparse/scaling

[Contract](scaling.md) · [Source](../src/sparse/scaling.ct)

| Export | Type or role |
| --- | --- |
| `equilibrate!` | Matrix * V * V * V -> Void; values, then row and column scale vectors to write, all powers of two |
| `scale!` | Matrix * V * V * V * V -> Void; values, row scales, column scales, caller output |
| `scale-vector!` | V * V * V -> Void; scales, input, output, which may be the input |
| `row-maxima!`, `column-maxima!` | Matrix * V * V -> Void; values, caller output |
| `DiagonalStatus` | Outcome type of equilibrate-symmetric! |
| `DiagonalScaled`, `DiagonalNotPositive` | DiagonalStatus values |
| `equilibrate-symmetric!` | Matrix * V * V -> DiagonalStatus; lower triangle with a positive diagonal, values, scale vector to write |
| `scale-symmetric!` | Matrix * V * V * V -> Void; values, scales, caller output |

The matrix supplies the pattern and the values array the numbers, in the matrix's stored
order; it may be the matrix's own entries. No call allocates.

## coalton-linalg/sparse/augmented

[Contract](quasidefinite.md) · [Source](../src/sparse/augmented.ct)

| Export | Type or role |
| --- | --- |
| `AugmentedPlan` | Fixed-pattern plan for `[[top I, Aᵀ], [A, −bottom I]]`; obtain from the analyze operations |
| `analyze` | coalton-linalg/sparse/csc:Matrix -> AugmentedPlan |
| `analyze-with-limit` | Matrix * UFix -> AugmentedPlan; explicit stored-entry budget |
| `dimension`, `leading-dimension`, `entry-count` | AugmentedPlan -> UFix; n + m, n, and the stored lower entries |
| `column-pointers`, `row-indices` | AugmentedPlan -> I; lower CSC pattern with every diagonal |
| `refill!` | AugmentedPlan * V * F64 * F64 * V -> Void; source values, top, bottom, caller output |
| `refill-csc!` | AugmentedPlan * Matrix * F64 * F64 * V -> Void; checks the source pattern first |
| `maximum-analysis-size`, `default-analysis-entry-limit` | UFix; the bounds shared with SPD analysis |

The output values belong to the caller and feed `factor-signed!` with a `Pattern` built
from the plan's dimension, column pointers and row indices, analyzed by
`analyze-quasidefinite` with the plan's leading dimension.

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

`FactorStatus` and `SolveStatus` here are distinct types from the sparse package's. An
exactly zero pivot is `FactorSuccess` with a nonzero `zero-count`. Inspection arrays are
borrowed, read-only by contract, and rewritten by the next `factor!`.

## coalton-linalg/dense/orthogonal

[Contract and derivation](orthogonal.md) · [Source](../src/dense/orthogonal.ct)

A dense m-by-n matrix is a V of m n elements in column-major order, entry (i,j) at i + j m.

| Export | Type or role |
| --- | --- |
| `Decomposition` | Reusable decomposition of one shape with private workspace; obtain from make-decomposition |
| `make-decomposition` | UFix * UFix -> Decomposition; rows, columns, under default-entry-limit |
| `make-decomposition-with-limit` | UFix * UFix * UFix -> Decomposition; explicit budget of m n entries |
| `maximum-dimension`, `default-entry-limit` | UFix; 32,768 and 16,777,216 |
| `row-count`, `column-count`, `step-count` | Decomposition -> UFix; step-count is the smaller dimension |
| `FactorStatus` | `FactorSuccess`, `FactorNonFinite`, `FactorBadSize` |
| `RankStatus` | `RankSuccess`, `RankInvalidFactor`, `RankTooLarge`, `RankSingular`, `RankNonFinite` |
| `SolveStatus` | `SolveSuccess`, `SolveInvalidFactor`, `SolveNoRank`, `SolveBadSize`, `SolveBadIndex`, `SolveNonFiniteRhs`, `SolveNonFiniteResult` |
| `factor!` | Decomposition * V -> FactorStatus; A P = Q R with column pivoting |
| `factor-ready?`, `rank-set?` | Decomposition -> Boolean |
| `pivot-magnitude` | Decomposition * UFix -> F64; magnitude of a diagonal entry of R |
| `rank-above` | Decomposition * F64 -> UFix; leading diagonal entries strictly above an absolute threshold |
| `set-rank!` | Decomposition * UFix -> RankStatus; complete the decomposition for a chosen rank |
| `selected-rank` | Decomposition -> UFix |
| `solve!` | Decomposition * V * V -> SolveStatus; minimum-norm least squares for the selected rank |
| `solve-residual-norm` | Decomposition -> F64; latest successful solve, otherwise -1 |
| `nullspace-vector!` | Decomposition * UFix * V -> SolveStatus; one orthonormal basis vector |
| `apply-q!`, `apply-qt!` | Decomposition * V -> SolveStatus; in place on a vector of row-count entries |
| `factor-permutation` | Decomposition -> I |
| `factor-work`, `factor-tau` | Decomposition -> V; R with reflector tails below it, and the reflector scalars |

The three status types are distinct from each other and from the other packages' types of
the same names. The rank is always the caller's choice; nothing is solved until
`set-rank!` succeeds.

## coalton-linalg/dense/jacobi

[Contract and derivation](jacobi.md) · [Source](../src/dense/jacobi.ct)

| Export | Type or role |
| --- | --- |
| `JacobiStatus` | `JacobiSuccess`, `JacobiNotConverged`, `JacobiNonFinite`, `JacobiBadSize` |
| `maximum-dimension`, `default-entry-limit`, `default-sweep-limit` | UFix; 32,768, 16,777,216 and 60 |
| `SymmetricEigen` | Reusable eigenvalue computation of one dimension; obtain from make-symmetric-eigen |
| `make-symmetric-eigen`, `make-symmetric-eigen-with-limit` | UFix -> SymmetricEigen; UFix * UFix -> SymmetricEigen with an entry budget |
| `eigen!` | SymmetricEigen * V -> JacobiStatus; lower triangle of a column-major n-by-n array |
| `eigen-with-sweep-limit!` | SymmetricEigen * V * UFix -> JacobiStatus |
| `eigenvalues`, `eigenvectors` | SymmetricEigen -> V; ascending values, and their unit vectors as columns |
| `eigen-ready?`, `eigen-dimension`, `eigen-sweeps` | SymmetricEigen -> Boolean, UFix, UFix |
| `Svd` | Reusable singular value decomposition of one shape; obtain from make-svd |
| `make-svd`, `make-svd-with-limit` | UFix * UFix -> Svd; UFix * UFix * UFix -> Svd with an entry budget |
| `svd!` | Svd * V -> JacobiStatus; a column-major m-by-n array |
| `svd-with-sweep-limit!` | Svd * V * UFix -> JacobiStatus |
| `singular-values` | Svd -> V; n values, descending |
| `left-vectors`, `right-vectors` | Svd -> V; m-by-n and a complete n-by-n, column-major |
| `rank-above` | Svd * F64 -> UFix; values strictly above an absolute threshold |
| `svd-ready?`, `svd-rows`, `svd-columns`, `svd-sweeps` | Svd -> Boolean, UFix, UFix, UFix |

With r from `rank-above`, columns r to n - 1 of the right vectors are an orthonormal
nullspace basis for any shape. A wide matrix costs far more than its transpose.

## Choosing the entry point

Use [the SPD example](../examples/spd-solve.ct) for analyze/factor!/solve! and
[the sparse example](../examples/sparse-operations.ct) for matrix construction
and products. [The normal-equations example](../examples/normal-equations.ct) joins Gram
products to the SPD factorization, and [the minimum-norm example](../examples/minimum-norm.ct)
uses the dense indefinite factorization. [The dependent-constraints
example](../examples/dependent-constraints.ct) picks a rank and solves through the
orthogonal decomposition, and [the singular-rank example](../examples/singular-rank.ct) reads
the same rank from singular values. All are compiled and exercised by the library tests.

For a rank decision there are three tools, in rising cost and reliability: the inertia of
[the dense indefinite factorization](indefinite-ldlt.md) either side of a shift, the
diagonal of R in [the orthogonal decomposition](orthogonal.md), and the singular values in
[the Jacobi package](jacobi.md). Each leaves the threshold to the caller.

A general CSC matrix represents exactly its stored entries; it has no implicit
symmetry. The SPD interface specifically interprets its lower triangle as a
symmetric matrix. The two contracts must not be interchanged.

Coalton checks declared argument types, not matrix dimensions, finiteness, array
ownership, or future handle validity. Those properties have runtime checks or
calling obligations in the contracts. Raw Common Lisp entry points require the
documented specialized arrays and are not a dynamically checked external API.