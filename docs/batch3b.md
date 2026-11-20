# Batch 3b: products into caller storage

State: completed; independently reviewed and verified at safety 1 and safety 0.
Measured 2026-11-19 on the pinned Windows/x86-64 SBCL/Coalton toolchain.

The CSC package now exposes matvec! and transpose-matvec!, returning Void and
writing into a caller-provided output vector. They share the original arithmetic
loop with allocating products. Input validation happens before output allocation
in convenience operations and before any output write in buffered operations.
See [the calling contract](primitives.md#caller-provided-output).

The output must differ by object identity from the input and matrix values.
Read-only input aliases with matrix values remain permitted. Dimensions, aliases
and all input values are checked before mutation, including unused inputs and
empty-array aliases. Prior output values are ignored and overwritten. Preflight
refusals preserve every array; NonFiniteResult leaves output unspecified and
preserves matrix/input. Successful calls need no workspace beyond caller output;
exceptions can allocate. A scratch/transaction interface was deliberately not
introduced without a consumer needing it.

The new InvalidAliasing constructor extends LinalgError. Callers making an
exhaustive match on that exception family must include it. Existing allocating
product arithmetic, validation order and input-preservation contracts are retained.

## Verification

- 104 exact cases, 416 products: empty, rectangular, stored-zero, cancellation and
  seeded families. Independent rational dot products supply expected results.
  Dyadic values keep every partial sum representable; cases span scales from
  2^-200 to 2^200. Each seeded dimension reaches 19; the roster contains no 19x19 case.
- Two separate 2^54,1,-2^54 cases check ordinary CSC accumulation order, rather than
  incorrectly requiring exact summation. Allocating and buffered result bits agree
  on all tested cases.
- 14 successful ownership/overwrite cases, 6 dimension refusals, 12 alias refusals,
  14 nonfinite-input refusals and 6 arithmetic refusals use real public entry points.
  Preflight snapshots compare exact bits, including NaN payloads and signed zero.
- Six actual compiled Coalton loops each warm 10,000 calls then measure 10,000 calls:
  both orientations on rectangular, empty, and permitted read-alias cases.
  Every measurement reported zero bytes allocated at both safety policies.
  The disposable runner pauses SBCL finalization during the measurement window
  and restores it afterward; see [the measurement policy](toolchain.md).
- Six planted oracle defects first pass a correct baseline, then must raise the
  checker's own mismatch condition; unrelated errors cannot count as detection.
- Three independently compiled source mutations in an isolated export were watched
  fail the intended public behavior assertions at safety 1: removing output zeroing,
  bypassing alias checks, and writing before input validation. The unmodified
  source passed the same tests. Main working sources were never mutated.

The source mutation outcomes were, respectively, success replaced by typed
NonFiniteResult; expected InvalidAliasing replaced by success; and changed output
on a preflight refusal. Their source hashes are in the batch measurement record; detailed logs remain
in ignored build/batch3b. These are behavioral red controls, not merely missing exports.

Both required serial full checks passed: scripts/check.ps1 -Safety 1, then
-Safety 0.
They also exercise the existing factor, solve, sparse-analysis, primitive,
encapsulation, fixture and toolchain gates. No Blender or application test is
part of this library batch. [The machine-readable record](batch3b-results.json)
identifies the measured source and outcomes.
