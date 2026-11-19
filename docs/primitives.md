# Numeric and sparse matrix primitives

These modules use F64 values and UFix indices on the supported 64-bit SBCL
toolchain. Coalton checks the calling types. Raw Lisp callers must supply simple
specialized arrays: DOUBLE-FLOAT for values, (UNSIGNED-BYTE 62) for indices.
Malformed contents signal Common Lisp errors. Floating-point traps must have
the documented disabled policy.

## Numeric

Package: `coalton-linalg/numeric`, source: [numeric.ct](../src/numeric.ct).

- `finite?(x)` is false for NaN and either infinity.
- `norm2(vector)` computes the Euclidean norm with scaled sum-of-squares
  accumulation. Empty/all-zero vectors return positive zero. Nonfinite input or
  a nonfinite computed result signals an error.
- `norm-inf(vector)` returns the greatest absolute component, or positive zero
  for an empty vector. Nonfinite input signals an error.

Neither norm mutates input. Norm2 avoids unnecessarily squaring the largest or
smallest magnitude directly. It does not promise correct rounding or compensated
accumulation of every term.

Derivation: represent a partial sum of squares as scale² times s. On seeing a
larger absolute value a, replace s by 1 + s*(scale/a)² and scale by a. Otherwise
add (a/scale)². Zero components leave this representation unchanged. Return
scale*sqrt(s). This is independently derived from that identity.

## Canonical rectangular CSC

Package: `coalton-linalg/sparse/csc`, source: [csc.ct](../src/sparse/csc.ct).

`Matrix` has a private constructor. Public constructors validate inputs and own
their storage. A matrix has m rows and n columns, including zero dimensions.
Its column-pointer array has length n+1, begins at zero, is nondecreasing, and
ends at the number of stored entries. Row indices are in range and strictly
increase within each column. Empty columns are valid. All stored values are
finite; explicit zeros remain stored.

| Function | Contract |
| --- | --- |
| `from-csc(m,n,pointers,rows,entries)` | Validate canonical CSC and copy all three arrays. |
| `from-coo(m,n,rows,columns,entries,policy)` | Validate unsorted coordinate entries and construct canonical owned CSC. |
| `row-count`, `column-count`, `entry-count` | Dimensions and number of stored entries, including zeros. |
| `column-pointers`, `row-indices`, `entries` | Borrowed read-only inspection arrays. |
| `transpose(A)` | A freshly owned, canonical transpose, preserving stored zeros. |
| `permute(A,row-order,column-order)` | Fresh B where B[i,j]=A[row-order[i],column-order[j]]. |
| `matvec(A,x)` | A freshly allocated vector A*x. |
| `transpose-matvec(A,x)` | A freshly allocated vector transpose(A)*x, without constructing the transpose. |

Coordinate arrays have equal length. Out-of-range coordinates and nonfinite
values are errors. `RejectDuplicates` signals an error if a coordinate repeats.
`SumDuplicates` accumulates each coordinate's entries **in original input order**,
using ordinary F64 addition. A nonfinite intermediate sum is an error even if
later cancellation would yield a finite exact result. Cancellation to zero
retains that coordinate. Different ordering of duplicate inputs can round
differently; the operation makes no order-independent summation promise.

Permutation arrays must be bijections of the corresponding dimensions. Both
are **new index to old index**, including non-involutive permutations. A row or
column permutation never changes the number of stored entries.

All operations leave their inputs unchanged. Constructors, transpose and
permute return independent storage. Inspection arrays must not be changed by
callers; the borrowed-view convention is contractual, not a runtime copy or
write-protection mechanism. There is no public matrix mutation interface in this
batch, and no zero-allocation claim for these allocating operations.

Products validate every vector component, including those unused by an empty
row/column. Input length mismatch or nonfinite input/output signals an error.
They use deterministic ordinary F64 accumulation in CSC traversal order. They
do not offer compensated products or an exact accumulator; a nonfinite
intermediate is not recovered through later mathematical cancellation.
No partial caller-owned output exists to expose on refusal.

## Cost and construction

Canonical CSC construction takes O(n + nnz) validation/copy time and owned
storage; it does not scan all m rows.
COO construction uses a stable bottom-up merge of coordinate indices by column
then row: O(nnz log nnz + n) time and O(nnz + n) temporary/owned storage.
A scan folds each equal-coordinate run under the selected duplicate policy.

Transpose counts entries per destination column, prefix-sums counts, then scatters
while visiting original columns in ascending order. This yields sorted destination
rows directly: O(m+n+nnz) time/storage. Permute validates both bijections and
uses the COO constructor to recover canonical order: O(m+n+nnz log nnz).
Products cost O(m+n+nnz), including input validation and result initialization.
Norms take O(vector length) time with constant scalar workspace.

## Consuming the structure in SPD analysis

`coalton-linalg/sparse/spd-ldlt:analyse-csc(matrix)` requires a square matrix whose
stored pattern is lower triangular with every diagonal present. It calls the
existing symbolic analysis, which owns its derived structure. Numeric values
remain explicit arguments to `factor!`; use `csc:entries` initially, or a new
values array in exactly the same order for refactoring.

A general CSC matrix means exactly its stored entries. It has **no implicit
symmetry**: matvec on lower-triangle storage multiplies that lower triangle.
The SPD factorization interprets its specified lower triangle as a symmetric
matrix. The bridge does not change general matrix-product semantics.

See [sparse-operations.ct](../examples/sparse-operations.ct) for a compiled Coalton
caller and [primitives.lisp](../tests/primitives.lisp) for independent tests.