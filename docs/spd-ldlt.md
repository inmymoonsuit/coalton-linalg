# Sparse SPD LDLᵀ

Package: `coalton-linalg/sparse/spd-ldlt`. Source:
[src/sparse/spd-ldlt.ct](../src/sparse/spd-ldlt.ct).

This page is the positive definite contract. The same package has two further numeric
entries over the same analysis, each with its own page:
[a Cholesky-form entry for proofs](cholesky-form.md) and
[a signed-pivot entry for quasi-definite matrices](quasidefinite.md).

## Interface

`Pattern(size, columns, rows)` describes the lower triangle in zero-based CSC.
There are n+1 column pointers, starting at 0 and ending at the row-array length.
Each column starts with its diagonal; rows then strictly increase below it.
Duplicates are invalid. Explicit numeric zeros remain structural entries.
The empty n=0 matrix is supported.

- `analyze(pattern) -> Symbolic` validates the structure and owns its derived
  data, using default-analysis-entry-limit. Malformed structure raises
  InvalidStructure. Caller pattern arrays are not retained.
- `analyze-csc(matrix) -> Symbolic` accepts a checked square lower CSC matrix
  with every diagonal present; values remain explicit arguments to factor!.
  A nonsquare matrix raises InvalidDimensions.
- `analyze-with-limit(pattern, maximum-entries)` and
  `analyze-csc-with-limit(matrix, maximum-entries)` accept an explicit total nnz(L)
  budget including the implicit diagonal. Their result and ownership are identical
  to the default-budget operations.
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

[A Cholesky-form entry](cholesky-form.md), `factor-cholesky!`, computes G with G Gᵀ in
place of L D Lᵀ over the same analysis, for proofs of positive definiteness; it overwrites
the factor and no solve accepts it.

After a successful factor, [the condition estimate](condition.md) gives a lower estimate of
the 1-norm of the inverse of the shifted matrix by solves with the same factor.

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

Analysis uses exact minimum degree, not AMD. Its indexed binary heap compares
(live degree, original index); compact mutable adjacency tables hold the live graph.
Saved elimination neighbors form the factor structure. Column counts allocate the
factor arrays; there is no unconditional n-by-n graph and no full-vertex scan per
pivot. With f = nnz(L) including the diagonal, live analysis storage is O(n + f),
and retained numeric/symbolic storage is O(nnz(A) + f + n). These are space bounds,
not total-allocation or speed guarantees. Hash-table capacity and runtime overhead
matter; garbage collection can keep more than the live data resident.

Heap work is O((n + f) log n). Pairwise clique formation examines
O(sum_k degree_k^2) candidate edges, with expected constant-time hash operations.
Sorting each saved column contributes O(sum_k degree_k log degree_k). Thus dense
fill can still give cubic analysis time and quadratic storage. Factor work depends
on fill; solve is O(n + f). No large sparse indefinite/QR capability is implied.

The default entry budget is `default-analysis-entry-limit = 1048576`. It counts
**all factor entries, including n implicit unit diagonals and structural zeros**.
It is an engineering resource limit, not a rank tolerance, a byte budget or a
promise that the host can allocate that much storage. Explicit-limit entry points
accept any UFix budget, including zero for an empty matrix. The effective budget
is the smaller of that request and maximum-analysis-size, so every resulting
array length remains representable; refusals report this effective budget. The monotone count is
input nnz plus distinct fill edges; on success it equals final nnz(L). Before
inserting an edge that would exceed the budget, analysis throws
`AnalysisStorageExceeded(required, budget)`. Required is the first observed lower
bound exceeding the budget, not a prediction of the completed factor's size.
The required-count payload uses Integer; the effective budget uses UFix.
No partially analyzed handle is returned.

Validation order is explicit: CSC input must first be square; dimension must then
be at most `maximum-analysis-size`, the representable dimension bound
min(array-dimension-limit - 2, floor(most-positive-fixnum / 2)) on supported SBCL.
The same bound caps the effective entry budget before growth.
This check precedes n+1 and heap child arithmetic. Exceeding it raises
`AnalysisLimitExceeded(requested, maximum)`. Next n must fit the entry budget;
then the complete pattern is validated, then input nnz is checked against the
budget, before graph allocation. The 4096 engineering dimension cap is removed.
Normal host allocation failures remain possible inside any entry budget.
See [typed errors](errors.md).

The elimination tree is recoverable from the smallest strict-lower row of each
factor column (an empty column is a root). A tree postorder does not replace the
public minimum-degree permutation. This implementation does not retain a duplicate
tree or postorder without a traversal consumer; independent structural checks
verify the tree implied by the factor. A later tree-based factor traversal may
retain that data with its own measured contract.

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
