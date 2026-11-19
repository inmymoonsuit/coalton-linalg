# Sparse SPD LDLᵀ

Package: `coalton-linalg/sparse/spd-ldlt`. Source:
[src/sparse/spd-ldlt.ct](../src/sparse/spd-ldlt.ct).

## Interface

`Pattern(size, columns, rows)` describes the lower triangle in zero-based CSC.
There are n+1 column pointers, starting at 0 and ending at the row-array length.
Each column starts with its diagonal; rows then strictly increase below it.
Duplicates are invalid. Explicit numeric zeros remain structural entries.
The empty n=0 matrix is supported.

- `analyse(pattern) -> Symbolic` validates the structure and owns its derived
  data. Malformed structure signals an error. Caller pattern arrays are not retained.
- `analyse-csc(matrix) -> Symbolic` accepts a checked square lower CSC matrix
  with every diagonal present; values remain explicit arguments to factor!.
- `factor!(symbolic, values, shift) -> FactorResult` computes
  `P (A - shift I) Pᵀ = L D Lᵀ`. Values follow the original CSC order; a changed
  structure needs a new analysis. Values are not retained.
- `solve!(factor, b, out) -> SolveStatus` solves in original caller order. b and
  out may be the same array. On any refusal out remains unchanged.

`factor-status(result) -> FactorStatus` returns one of:

| FactorStatus | Meaning |
| --- | --- |
| FactorSuccess | Success; result-factor returns the borrowed factor |
| FactorNonPositive | Computed pivot is nonpositive |
| FactorNonFinite | Nonfinite input or arithmetic |
| FactorBadValues | Wrong value count |
| FactorBadShift | Nonfinite shift |

For FactorNonPositive and FactorNonFinite, result-column is the current
elimination column in original numbering. It is not necessarily the location
of a nonfinite input entry. Other outcomes return -1. Calling result-factor
on a refusal signals an error.

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
including a refusal, overwrites or invalidates earlier contents. A FactorResult retains
its original outcome. Its successful factor handle is a borrowed view, not a
snapshot: a later successful refactor changes its data, and a later refusal
invalidates it even though the earlier result still reports FactorSuccess.
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

The tested Coalton setup disables floating-point traps. Explicit nonfinite checks
provide the documented refusals; a caller enabling traps changes that behavior.
Underflow is not diagnosed separately. Bitwise determinism is tested for repeated
input on the same runtime and policy; cross-platform identity is not promised.

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

Private nominal handle types wrap private records. A public alias of a private
record would still permit Coalton's type-directed dot-field access; the external
compile probes cover that distinction. Handle wrappers are created during
analysis, and hot calls unwrap them to use the existing workspace.

Cached result and factor identity are part of the calling contract. The opaque
FactorResult and Lisp ABI on identity-returning accessors preserve that contract;
factor! retains the compiler's automatic calling convention for its F64 shift.
The [external caller](../tests/callers.ct) supplies computed, changing shifts,
and the [interface tests](../tests/interface.lisp) check cached identity and
allocations. These choices are deliberate review points for a language maintainer.
