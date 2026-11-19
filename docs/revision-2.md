# Pre-publication revision

State: implemented and verified locally, 2026-11-19. Comparison base: `78f2daa`.
The [initial review snapshot](initial-review.md) is a historical record; the
[current API](api.md) and [machine-readable results](revision-2-results.json)
describe this revision. Nothing has been pushed.

## Resulting interface

- `factor!` returns FactorStatus directly. FactorResult, factor-status,
  result-factor and result-column are removed. `symbolic-factor` supplies the
  same borrowed factor handle; `refused-column` describes the latest call and
  resets on every call. Saved statuses do not certify current factor readiness.
- Symbolic and NumericFactor have transparent private constructors. Their runtime
  values are the underlying records, with no separate wrappers or Lisp-ABI
  overrides. Real external callers check identity and encapsulation.
- `analyze` and `analyze-csc` use American spelling. The numeric package exports
  only norm2 and norm-inf; callers use Coalton's standard finite? themselves.
- A small [typed exception family](errors.md) identifies construction/primitive
  failures. Factor and solve retain their distinct ordinary status returns.
- Analysis refuses dimensions above the exported maximum-analysis-size, 4096,
  before size arithmetic or allocation. The 16 MiB bound covers the adjacency
  bitmap only; total storage and factor fill can be larger. It is an engineering
  limit pending sparse graph analysis, not a numerical threshold.

These are deliberate pre-publication API changes. Examples, tests and optional
fixture tools use the new signatures. Arithmetic, ordering, shift semantics and
positive-pivot policy remain unchanged.

## Measured verification

The complete check runner passed serially at safety 1, then safety 0 on the pinned
Windows/x86-64 toolchain. Existing installed dependencies and cache were reused;
this revision did not repeat the earlier public-source acquisition exercise.
Content verification covered 18 components and 1,314 files, with 28 loaded
ASDF dependency systems resolved to their pinned sources.

Both policies passed all 20 base matrices and 336 seeded matrices, 252
metamorphic relations, 96 repeated refactors, primitive tests, executable examples,
checker controls and tool isolation. All 2,090 numerical output artifacts (74 base,
2,016 stress) were byte-identical across policies and to the previous snapshot;
a second agent independently compared both rosters and bytes. This is bounded
regression evidence, not a proof for arbitrary matrices.

New interface checks cover 23 primitive error cases, duplicate coordinates in
both input formats, Lisp condition interoperability, typed SPD validation and
resource-limit payloads, the actual 4096-row diagonal boundary, rejection of 4097
and maximum-UFix dimensions, latest diagnostic resets and independent workspaces.
Eight encapsulation probes and four status-type probes reject invalid external
callers; valid controls compile and run. All ten status outcomes are exercised.

Fourteen warmed allocation paths measured zero bytes for 10,000 calls each at
both safety policies, as did the existing computed-shift Coalton caller. Three
handle routes preserve identity: raw Lisp, direct cross-file Coalton and an
accessor supplied at runtime as a first-class function. Construction and cold
startup are outside the allocation claim. The runner's existing finalizer
isolation and restoration control remain unchanged.

## Failures observed before correction

With only the exception definitions and tests added, all 22 failing primitive
operations produced the old untyped error category; a valid operation control
passed. The old analysis also returned the wrong category for maximum-UFix with
empty arrays, safely demonstrating the missing early limit gate. The same checks
pass after implementation.

Ordinary calls to the standard generic finite? introduced allocation: the 3x3
factor-success test measured 1,769,472 bytes over 10,000 calls. A separate probe
performed 100,000 finite checks: ordinary calls allocated 1,572,864 bytes and
call-site inline hints allocated zero. Explicit `(inline (math:finite? ...))`
hints (or the imported name in SPD) retain the standard implementation while
exposing its F64 operations to the compiler. The full allocation checks then
passed with their original zero-byte threshold. No global specialization or
upstream dependency modification was needed.

## Exact-input numerical boundary

The external review supplied three floating Gram candidates. Exact rational
inspection of their serialized binary F64 values found all three full rank and
positive definite, so none was imported as an exactly singular example.

The independent generic witness is `[[3,1,4],[1,1,2],[4,2,6]]`. Its entries are
exactly representable; its null vector is [-1,-1,1], its leading 2x2 minor is 2,
and its determinant is zero, proving rank two. The baseline and revised library
accept it with final computed pivot 2^-52. The retained test proves rank from the
actual loaded F64 entries using exact rational arithmetic. Changing only the last
diagonal to 7 gives a rank-three control with determinant 2 and a checked solve.
No hidden rank tolerance was introduced.

## Publication and remaining scope

The unpublished initial commit's author and committer emails were corrected to
the maintainer-supplied GitHub noreply identity; this repository uses it for
future commits. The license remains MIT. No dependency implementations, external
review files, held-out matrices, local configuration or generated logs were added.

Sparse graph scalability remains Batch 3. Buffered matrix-vector operations,
permutation optimization, broader numeric types and additional hosts remain
separate work. Matrix representation was retained: owned canonical matrix data
and mutable workspace handles have different responsibilities. Independent
source/test reviewers found no remaining substantive issue in this batch.
