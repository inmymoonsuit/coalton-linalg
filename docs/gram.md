# Fixed-pattern sparse Gram products

Package: coalton-linalg/sparse/gram. Source: [gram.ct](../src/sparse/gram.ct).
State: completed as Batch 3c; verified at safety 1 and safety 0. The
[batch record](batch3c.md) holds the measurements and the controls.

## Analysis and structure

`analyze-ata(A)` prepares A^T A; `analyze-aat(A)` prepares A A^T. The input is a
checked general rectangular CSC Matrix. Values are not inspected, copied or
retained by analysis. The plan owns a copy of both source dimensions and its CSC
pattern, plus numeric-refill workspace. Changing numeric values needs no new
analysis; changing the source pattern or either dimension does.

The output is **only the lower triangle of a symmetric Gram matrix**, in canonical
zero-based CSC. Every diagonal is present, including zero diagonals for empty
source columns (ATA) or rows (AAT). Potential entries remain structural even when
stored zeros or cancellation make the present value zero. General CSC matvec has
no implicit symmetry and must not be mistaken for a symmetric product on this
lower storage.

A private transparent GramPlan handle exposes `dimension`, `column-pointers`,
`row-indices` and `entry-count`. Inspection arrays are borrowed read-only views.
Analysis returns a new independent plan; repeated inspectors return its existing
arrays. Numeric output is caller-owned, with no borrowed numeric result or cached
result-validity flag.

`analyze-ata-with-limit(A,budget)` and `analyze-aat-with-limit(A,budget)` select the
maximum number of stored Gram entries, including every explicit diagonal.
Defaults use default-analysis-entry-limit=1,048,576. The effective budget is capped
by maximum-analysis-size, the same representability bound used by SPD analysis.
An oversized output dimension raises AnalysisLimitExceeded before arithmetic or
allocation based on that dimension. A diagonal lower bound or counted output
pattern over budget raises AnalysisStorageExceeded(required,effective-budget).

This budget counts **Gram pattern entries**, not input entries, bytes or arithmetic
operations. A tall one-column source may contain far more entries than its scalar
ATA result and still succeed with budget 1. Source-pattern ownership and incidence
workspace require their own O(source-columns + nnz(A)) storage. Allocation failures
can occur inside any output budget; a small output pattern can still require much
input storage and repeated work.

## Numeric refill

- `refill!(plan,values,out) -> Void` accepts values in exactly the analyzed source
  CSC order. It validates the lengths of values and out, rejects their object
  identity with InvalidAliasing, then validates every finite input before writing.
  Equal length cannot detect incorrectly reordered numeric values.
- `refill-csc!(plan,A,out) -> Void` first compares **both** source dimensions,
  column pointers and row indices with the owned source pattern. Dimension changes
  raise InvalidDimensions; any pattern difference raises InvalidStructure. It then
  applies the ordinary refill checks to the matrix's entries.

Output is a caller-owned F64 array of entry-count(plan) elements. It must not be
any library's borrowed inspection array. Runtime alias detection compares it with
the actual input values array, including when empty. Existing output contents need
not be finite: every slot is overwritten. Never use csc:entries of a constructed
output Matrix as the destination; those arrays are read-only by contract.

A product or sum too large to hold is NonFiniteResult. One too small is not
refused: a refill returns what its products round to, so source entries below
about 10⁻¹⁵⁴ give squares that are zero or have lost digits, with no signal.
A factorization of such a result refuses a pivot; a caller with entries that
small scales them first. [The scale-sweep roster](rosters.md#what-is-declared)
records where this begins.

All preflight failures leave output, input and plan structure unchanged. Once
arithmetic starts, NonFiniteResult leaves all output contents unspecified; input
and structure remain unchanged. The plan remains reusable after every refusal:
the next call reinitializes output and its scatter workspace. Numeric results
remain with the caller and cannot be changed by later refills into other buffers.

A plan owns mutable scratch, so calls sharing it are not reentrant or concurrent,
even with distinct outputs. Separate plans permit independent calls. Successful
refills allocate no new arrays; warmed zero-allocation claims require measurement
on the pinned runtime. Typed exceptions may allocate.

## Arithmetic and storage

For ATA, G[i,j] = sum_k A[k,i] A[k,j]. For AAT,
G[i,j] = sum_k A[i,k] A[j,k]. Each stored output accumulates ordinary F64 products
in increasing original inner index k. There is no compensation, pruning, damping,
rank threshold or definiteness decision. A zero or singular product can succeed.
Rounded Gram entries need not retain exact algebraic rank or positive
semidefiniteness. Input finiteness and final output finiteness are checked;
overflow is not recovered by later exact cancellation. Underflow follows F64.

Internally normalize to G=B^T B with B=A or A^T, retaining original numeric-value
positions. Incidence groups contain only represented inner indices. An unused
inner dimension is never expanded into a row-pointer array: for example, ATA of
one entry at row maxUFix-1 in a maxUFix-by-1 matrix still has scalar structure.
AAT of that matrix refuses its unrepresentable output dimension.

Column candidates are marked and counted before output-pattern allocation. The
pattern is then filled by row: row i is appended to every column at or before it
that it shares an inner index with, so each column receives its rows in increasing
order, diagonal first, and nothing is sorted. The marker stamps need no reset
between the passes: the first stamps an index only with a column at or before it,
and the second asks whether a column carries the stamp of a row after it.
Numeric refill maps each current output column's rows into preallocated scatter
positions and accumulates directly; no list of all individual multiplication
contributions is cached.

Each stored entry of B keeps two indices into the sorted incidence: its own
position and the exclusive end of its inner-index group. Within a group the
owning columns ascend, because B is canonical CSC with one entry per column and
inner index. So the range from an entry's own position to its group's end holds
exactly the columns at or after its own, which are the lower-triangle partners.
The counting pass and refill scan that range and nothing before it. The fill
uses the other half, from the group's beginning to the entry's own position,
which analysis keeps only while it runs.

The incidence is ordered by inner index and then by position without comparing:
positions start in increasing order and are distributed by one base-2048 digit of
the inner index at a time, least significant first, each pass keeping the order it
was given among equal digits. The number of passes follows the largest inner index
present, at most six, and no array is indexed by an inner index.

With d the output dimension, s=nnz(A), g=nnz(G), live and retained storage is
O(source-columns + s + d + g). Source ordering costs O(s) for each digit of the
largest inner index; symbolic candidate traversal and refill work are sum_k degree_k (degree_k + 1) / 2, where degree_k is
a source row degree for ATA or column degree for AAT: one product for each
lower-triangle pair sharing inner index k. Refill also validates values and
initializes output. The checked-Matrix variant additionally compares the entire
source CSC pattern. Dense Gram growth is possible for sparse input.

## Factorization integration

Prepare a Gram plan and then one SPD analysis from its lower CSC structure:
`ldlt:Pattern(gram:dimension(plan), gram:column-pointers(plan), gram:row-indices(plan))`.
Allocate the Gram output values array once. Refill it each iteration, then pass it
to factor! with the caller's chosen shift. The product step makes no SPD or rank
promise, so the caller still checks the factorization status.

A compiled example is [normal-equations.ct](../examples/normal-equations.ct): one
plan and one symbolic analysis, then refill-csc!, factor! and solve!.

## Which refill a repeated loop should call

`refill!` trusts that its values are in the analyzed source order and costs only
the products. `refill-csc!` first compares the whole source pattern, O(nnz(A)),
on every call. Use `refill-csc!` where a Matrix arrives from elsewhere and its
pattern is not already known to be the analyzed one; use `refill!` inside an
iteration that rewrites values in place.

## Artifact driver

The optional tools system writes reviewable outputs for a folder of systems:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/check.ps1 -Task Gram `
  -InputFolder <folder> -OutputFolder build/<name> -AllocationRepeats 1000
```

Each `<name>.J.txt` holds `rows columns nnz`, then zero-based `row column value`
triplets in any order; duplicates are refused. `<name>.meta.txt` must carry a line
`gram: Jt J` or `gram: J Jt` naming the product. A system without one is a driver
error, never a default.

For each system the driver writes `<name>.G.out.txt` (`n nnz`, then
`row column value` for the lower triangle with every diagonal, by column then row,
17 significant digits), and `<name>.G2.out.txt` from a second refill with every
input multiplied by -1/2, whose values must be exactly a quarter of the first's.
`gram-table.txt` and standard output carry one line a system:
`name orientation n nnzG analyse_ms refill_ms refill_bytes`, the times being the
best of the requested repeats and the bytes measured over at least 1000 warmed
refills. A typed library refusal writes only `<name>.refused.out.txt` holding the
LinalgError constructor's name. Reruns replace every artifact of a system.
