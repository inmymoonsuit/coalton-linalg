# Sparse SPD LDLᵀ

Package: `coalton-linalg/sparse/spd-ldlt`. Source:
[src/sparse/spd-ldlt.ct](../src/sparse/spd-ldlt.ct).

## Interface

`Pattern(size, columns, rows)` describes the lower triangle in zero-based CSC.
There are n+1 column pointers, starting at 0 and ending at the row-array length.
Each column starts with its diagonal; rows then strictly increase below it.
Duplicates are invalid. Explicit numeric zeros remain structural entries.
The empty n=0 matrix is supported.

- `analyze(pattern) -> Symbolic` validates the structure and owns its derived
  data. Malformed structure raises InvalidStructure. Caller pattern arrays are not retained.
- `analyze-csc(matrix) -> Symbolic` accepts a checked square lower CSC matrix
  with every diagonal present; values remain explicit arguments to factor!.
  A nonsquare matrix raises InvalidDimensions.
- `factor!(symbolic, values, shift) -> FactorStatus` computes
  `P (A - shift I) Pᵀ = L D Lᵀ`. Values follow the original CSC order; a changed
  structure needs a new analysis. Values are not retained.
- `solve!(factor, b, out) -> SolveStatus` solves in original caller order. b and
  out may be the same array. On any refusal out remains unchanged.

`factor!` directly returns one of:

| FactorStatus | Meaning |
| --- | --- |
| FactorSuccess | Success; symbolic-factor returns the now-ready borrowed factor |
| FactorNonPositive | Computed pivot is nonpositive |
| FactorNonFinite | Nonfinite input or arithmetic |
| FactorBadValues | Wrong value count |
| FactorBadShift | Nonfinite shift |

`refused-column(symbolic)` is the latest factor! call's elimination column in
original numbering for FactorNonPositive and FactorNonFinite. It is not necessarily
the location of a nonfinite input entry. It is -1 before any call and after every
other outcome, including BadValues, BadShift and empty-matrix success. Read/copy
it before the next factor! call if a diagnostic must be retained. Solve calls do
not change it. `symbolic-factor(symbolic)` always returns the same owned handle;
`factor-ready?` reports whether that handle currently holds a successful factor.

SolveStatus is a distinct type with these outcomes:

| SolveStatus | Meaning |
| --- | --- |
| SolveSuccess | Output contains the solution |
| SolveInvalidFactor | Factor has not succeeded or was invalidated |
| SolveBadSize | Input/output vector-size mismatch |
| SolveNonFiniteRhs | Nonfinite right-hand side |
| SolveNonFiniteResult | Nonfinite computed solution |

Use Coalton matching or equality with the named values. There is no numeric
encoding contract. Raw Lisp callers can obtain a status value through, for
example, `(coalton:coalton coalton-linalg/sparse/spd-ldlt:SolveSuccess)`.

Inspection functions return permutation, CSC column pointers/rows, strict-lower
values of L, diagonal D, and factor-ready?. The diagonal of L is implicitly one.
The fixture tool's nnzL count includes this unit diagonal: n plus the stored
strict-lower entries, including stored zeros. It is not the length of factor-values.
Permutation[k] is the original index at elimination position k, so
`(P x)[k] = x[permutation[k]]`.

## Ownership and calls

Symbolic and NumericFactor are public types with private constructors. Analysis
creates their handles; callers using the public interface cannot construct or destructure them or
use generated record access to their workspace fields. Public inspection arrays
remain borrowed, ordinary mutable arrays, not statically protected views.

One Symbolic owns one mutable factor and scratch workspace. Every factor! call,
including a refusal, overwrites or invalidates earlier contents. A saved
FactorStatus describes its original call and carries no factor or column snapshot.
The factor handle is a borrowed view: a later successful refactor changes its data,
and a later refusal invalidates it even though an earlier status remains FactorSuccess.
factor-ready? and solve! inspect the factor's current validity.
Inspection arrays are borrowed, read-only views. Do not mutate them, pass them as
output buffers, or alias them with input values.

Calls sharing one Symbolic/factor are not reentrant or concurrent. Use separate
analyses for concurrent work. Repeated factor!/solve! calls reuse storage.
The tested allocation contract is after runtime/caller warmup; construction,
analysis, inspection by copying, file output, and cold startup are outside it.

Coalton callers receive static checking. Raw Common Lisp callers must supply the
declared simple specialized arrays: F64 uses SIMPLE-ARRAY DOUBLE-FLOAT and DOUBLE-FLOAT shifts; UFix uses
SIMPLE-ARRAY (unsigned-byte 62) on the tested 64-bit SBCL. These raw entry points are
not a dynamically checked foreign-language API.

## Numerical assumptions and limits

No pivoting, jitter, automatic damping, rank tolerance, or zero pruning is added.
A pivot must be finite and strictly positive. This floating-point decision is not
an exact rank or positive-definiteness certificate. The caller chooses its shift.
The generic regression matrix `[[3,1,4],[1,1,2],[4,2,6]]` is exactly rank two:
its null vector is `[-1,-1,1]` and its leading 2x2 minor is 2. Nevertheless, the
supported runtime accepts it with a final computed pivot of 2^-52. The regression
proves the input rank using exact arithmetic on the actual F64 entries; it does
not infer singularity from the way a rounded matrix was generated.

The tested Coalton setup disables floating-point traps. Explicit nonfinite checks
provide the documented refusals; a caller enabling traps changes that behavior.
Underflow is not diagnosed separately. Bitwise determinism is tested for repeated
input on the same runtime and policy; cross-platform identity is not promised.

Analysis accepts dimensions at most `maximum-analysis-size`, currently 4096.
Larger dimensions raise `AnalysisLimitExceeded(requested, maximum)` before pattern
validation, dimension arithmetic or workspace allocation. Comparing the dimension
first prevents oversized UFix inputs from overflowing n+1 or n*n. This engineering
cap bounds the adjacency bitmap to 16 MiB; it is not an accuracy threshold, a total
memory budget or a performance guarantee. Analysis of a nonsquare CSC matrix first
reports InvalidDimensions. See [typed errors](errors.md).

Analysis uses exact minimum degree, not AMD. Its temporary adjacency bitmap takes
n² bytes, with further ordering/frontier storage. A simple dense-graph worst-case
time bound is O(n³). Numeric storage follows input and factor structure; factor
work depends on fill, while solve is O(n + nnz(L)). Factor fill can itself become
dense. There is no large-scale sparse-analysis claim.

## Derivation

Initially connect every structural off-diagonal entry, including zeros. Select
the live vertex of minimum degree, with smallest original index breaking ties.
Connect each pair of surviving neighbors and eliminate the vertex. Those
neighbors form its strict-lower factor column; translate to elimination order.

For numeric column k, scatter the permuted A column and subtract the shift at k.
Each previous column j containing row k contributes
`L[k,j] D[j] L[k:n,j]`; subtract it over the precomputed structure. The diagonal
work entry becomes D[k]; divide the remaining column entries by that pivot.
Elimination cliques ensure these updates stay within the symbolic frontier.

For solving, gather Pb, solve Ly=Pb, divide Dz=y, solve Lᵀw=z, then scatter x=Pᵀw.
Only commit out after every computed result is finite.

## Representation choices

Private nominal handle types use `(repr :transparent)` around private records.
They remain distinct Coalton types with private constructors, but their runtime
value is the underlying record reference, without a separate wrapper. A public
alias would still permit type-directed dot-field access; external compile probes
cover that distinction.

The pinned compiler excludes transparent types from product-return unpacking.
The public handle-returning accessor therefore uses the automatic ABI, with no
Lisp-ABI override. FactorStatus and SolveStatus use symbol-backed enums. The
[external caller](../tests/callers.ct) exercises direct and first-class handle
access plus computed F64 shifts; [interface tests](../tests/interface.lisp) check
identity, latest diagnostics, factor lifetimes and allocations. These are measured
contracts for the pinned compiler, not promises about every future compiler.
