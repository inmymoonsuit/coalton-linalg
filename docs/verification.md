Latest interface and setup results: [initial review snapshot](initial-review.md).
The migration measurements below remain a historical record.

# Verification and supported toolchain

Current batch results are in [batches-1-2.md](batches-1-2.md).
The migration measurements below remain a frozen historical record.

## Reproduce

The complete check command is:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\check.ps1
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\check.ps1 -Safety 0
```

Run serially from this repository, with the ignored local configuration pointing
to installed dependencies. The checked-in example configuration contains no
desktop-specific paths. Command-line settings override that configuration.

Measured 2026-11-19 on Windows x86-64: SBCL 2.6.10, Coalton commit
`b6e284bdd91ef9ba0fb303e0798c9818b9b3e268`, release mode, Python 3.13.13.
The Python verifier uses only the standard library. The initial migration recorded
versions only; Batch 1 added [toolchain.lock.json](../toolchain.lock.json) and
[source/resolution checks](toolchain.md), verified through a fresh-cache build.

Ordinary ASDF loading requires registered dependencies. Use a fresh process with
`COALTON_ENV=release` set before loading Coalton; development-mode performance is
not verified. The prepared setup also enables specialization and emitted type
annotations. Do not reuse dependency caches across incompatible compiler modes.
The script forces a rebuild of our library when changing its safety feature;
unchanged dependency caches are reused. A new empty cache rebuilds the installed
sources. No dependency installation is performed.

`asdf:test-system "coalton-linalg"` runs the library interface/runtime/example
tests. The PowerShell command additionally runs optional tooling tests, emits
generic matrix artifacts, tests the independent verifier, verifies all artifacts,
and regenerates/compares fixture contents. It also runs 336 seeded cases,
independent primitive tests, and content/source-resolution checks.
Failures signal errors and propagate
to a nonzero process exit.

## Recorded migration result

Both safety policies passed. [verification-results.json](verification-results.json)
records the counts; generated logs and artifacts are under ignored `build/`.

- 20 generic matrices: 14 successful factors/solves and 6 expected refusals.
- 75 fixture files reproduced from the independent generator.
- 74 factor/solution/refusal artifacts byte-identical across both policies and
  the original prototype's generic results.
- Maximum independent normwise backward error: 2.18417455838898e-17.
- The independent verifier rejected all 10 planted artifact defects.
- Interface tests covered malformed CSC, input ownership, refactoring lifetime,
  refusal invalidation, structural zeros, nonfinite input/arithmetic, aliasing,
  untouched output on refusal, empty systems, and repeated-call determinism.
- The compiled Coalton example returned [1,2]; a cross-package caller exercised
  10,000 computed, changing F64 shifts.
- Primary-system loading was checked to exclude tools, tests, and examples.
- The parser's missing-order-check mutant was detected and the original restored.
  Success/refusal/success reruns replaced stale owned artifacts correctly.

The algorithms are unchanged after namespace/safety-feature normalization.
[migration.json](migration.json) records original source hashes and this migration's
destination hashes; it is a frozen record, not a claim about future source versions.

## Allocation evidence and limits

The current interface measurement runs full GC, performs a fixed 10,000-call
warmup, then measures a separate 10,000 calls without another GC between them.
All 11 success/refusal paths and the changing-F64 caller measured zero bytes
after warmup at both policies. A deliberately allocating cons-cell caller fails
the same measurement, demonstrating that the warmup does not disable detection.

Initial measurements with the previous 64-call warmup reported roughly
131–164 KB during the first factor-success measurement. Repeated isolated
measurements then reached zero (one recorded series: 32,752, 0, 0 bytes).
At migration, the mechanism causing that transient had not been established.
Batch 1 subsequently isolated a background finalizer-thread contribution using
[controlled fresh-process probes](runtime-allocation.md). The inspected
factor wrapper and unboxed worker showed no allocation calls; this does not
prove that an entire cold process allocates nothing. SBCL describes
[GET-BYTES-CONSED](https://www.sbcl.org/manual/#Introspection-and-Tuning)
as a program-wide counter. Report only the measured steady-state result.
The controlled first factor batch reported zero bytes with background work
separated. This does not certify every possible call/input or cold-start allocation.
The canonical disposable runner now pauses/restores its own finalizer thread
during bounded checks; raw ASDF calls do not. See [the protocol](toolchain.md).

Local diagnostic logs include `build/compiled-caller-red.log`,
`build/repeat-allocation.log`, and `build/parser-asdf-red.log`.
The parser log records the ASDF package-context defect before its fix.
The permanent test reads its mutation source with an explicit package binding.

## Numerical oracle

The Python tool imports no library implementation. It reconstructs P A Pᵀ and
L D Lᵀ, including positions absent from the input pattern, and computes residuals
using separately accumulated products. It checks finite positive D, permutation,
structural zeros, fill structure, and shifted factor/refusal artifacts.

For matrices up to 32 rows it additionally uses exact rational dense
factorization and independently constructed graph order. All 20 checked-in
cases are small enough for that oracle. Exact arithmetic is a comparison tool,
not an assertion that floating-point pivot decisions must agree arbitrarily near
a numerical boundary.

For larger supplied fixtures, refusal validation can rely on declared metadata
and a column-range check; it does not independently recompute the first failing
pivot. The small-suite result must not be described as certifying arbitrary
large refusals. Fill targets are optional external benchmark criteria, not
mathematical correctness guarantees. No fixture has a name-based exemption.

## Optional external fixtures and inspection

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\check.ps1 -Task Fixtures -InputFolder <fixture-directory> -OutputFolder .\build\external
python .\tools\verify.py verify <fixture-directory> .\build\external --require-shifts
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\check.ps1 -Task Inspect
```

The fixture task emits artifacts and internal checks; run the independent Python
command separately for external fixtures. Input files are read only. Output
folders used by the wrapper must be under this repository's `build/`.
Reruns replace only the driver's named outputs for selected fixtures and its table.

Each `name.A.txt` has an n/nnz header and lower-triangle row/column/value triplets,
sorted by column then row, including diagonals. `name.b.txt` is the RHS.
Optional `.meta.txt` carries shift/reference metadata and `.x.txt` a reference
solution used only by the verifier. The fixture protocol requires n>=1;
the direct library interface also supports n=0.

Outputs are `.perm.out.txt`, `.L.out.txt` (strict-lower triplets),
`.D.out.txt`, and `.x.out.txt`, or `.refused.out.txt`.
Shifted variants use `.shift-kept` or `.shift-refused` before the suffix.
Indices are zero-based and refusal columns use original caller numbering.

The driver contains dense reconstruction scratch (O(n²) doubles) and is not a
large-matrix production client. Timings in routine checks are smoke measurements,
not performance claims. Windows uses QueryPerformanceCounter; other hosts fall
back to the Common Lisp internal clock and are untested. A reported zero duration
means below clock resolution, not zero execution cost. Comparable benchmarks
need sufficient batching, a documented clock, and a controlled workload.