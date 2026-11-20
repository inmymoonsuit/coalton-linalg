# Batch 3a: compact symbolic analysis

State: completed; independently reviewed and verified at safety 1 and safety 0.
Measured 2026-11-19 on Windows/x86-64 SBCL 2.6.10 and the pinned Coalton compiler.
The machine-readable [measurement record](batch3a-results.json) accompanies this page.

## Change and review

Analysis uses live adjacency hash tables and an indexed heap ordered by exact
minimum degree, then original vertex index. It no longer allocates a square bitmap
or scans every vertex at each elimination. Column counts allocate the factor
structure; numeric factorization and solve algorithms are unchanged.

The independent design review rejected sorted-array graph replacement because hub
updates could copy quadratic total storage even on sparse families. It also caught
the heap invariant: update one degree key and repair immediately, rather than
changing all keys before repairs. The implementation uses Coalton's existing hash
tables; no numerical implementation was consulted or copied.

A second reviewer examined resource errors and the measurement harness. Source
review verifies the fill-budget guard precedes insertion. Public tests verify
refusals and payloads; they do not inspect private insertion timing. Tree/ancestor
and deterministic postorder identities are checked from exported L. Production
tree storage is deferred until a factor traversal needs it, avoiding unused state.

## Calling-contract changes

`analyze` and `analyze-csc` now default to at most 1,048,576 **total factor entries**,
including the implicit unit diagonal and structural zeros. New explicit-limit
entry points accept a caller budget, capped by the documented representability
bound. AnalysisStorageExceeded reports the first known required count and effective
budget, before excessive fill is inserted. See [the full contract](spd-ldlt.md).

This deliberately changes resource refusal behavior: dimensions above 4096 can now
succeed, while fill-heavy matrices below 4096 that previously passed analysis can
exceed the new default budget. Use an explicit entry budget when appropriate; it
is not a byte budget or an allocation guarantee. Exhaustive matches on LinalgError
must include the new AnalysisStorageExceeded constructor.

## Measurements

One analysis call per row below; these are observations, not performance targets.
Time and bytes-consed cover analyze only, excluding compilation, input parsing and
output. Bytes-consed is program-wide total allocation, not live or peak memory.
The last column is the Windows process peak working set, including compilation,
runtime, parsing, analysis and output. The process waits at both sampling points.
For these eight runs the final peak equaled the peak already reached by the loaded
compiler/runtime; subtracting the baseline would not give analysis memory.

| Family | n | nnz(L), with diagonal | Analysis ms | Allocated MiB | Process peak MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| chain-4097 | 4,097 | 8,193 | 2.550 | 3.31 | 150.50 |
| chain-16384 | 16,384 | 32,767 | 10.995 | 13.26 | 150.67 |
| chain-65536 | 65,536 | 131,071 | 45.536 | 53.02 | 150.64 |
| star-16384 | 16,384 | 32,767 | 11.374 | 13.70 | 150.73 |
| grid-65x65 | 4,225 | 73,991 | 49.959 | 16.22 | 150.47 |
| bipartite-fill-4480 | 4,480 | 110,880 | 53.283 | 23.02 | 150.41 |
| hub-cycles-2048 | 6,145 | 16,385 | 4.781 | 5.35 | 150.60 |
| k2r-8192 | 8,194 | 24,579 | 6.678 | 7.65 | 150.71 |

For n=65,536 the previous unconditional bitmap alone would require 4 GiB. That is a
calculated comparison, not an old implementation benchmark: its 4096 cap prevented
that run. The new chain run completed with the process peak shown above. Sparse
numeric storage and compact analysis do not prevent quadratic fill or cubic
clique-formation work on other graph families.

## Verification and reproduction

The same 4097-row chain regression failed on the published baseline with
AnalysisLimitExceeded(4097,4096), then passed after the rewrite. The retained ignored
logs distinguish early compiler integration errors from successful validation.

The structural roster contains all labeled graphs up to four vertices, 64 seeded
cases and generic chain, star, grid, wheel, disconnected, bipartite and shared-hub
families: 161 cases total, 13 above 4096. A set-graph oracle verifies exact ordering and
fill; analytic large-family formulas are cross-checked against it on 105 small
instances. Four planted defects are rejected: wrong ties, missing fill, extra fill
and malformed column pointers. Every structure is checked for forest/ancestor and
postorder identities independently of production data structures.

Both complete policy gates passed. Safety 1 first checked 160 cases and then the
added 65,536-row chain, followed by verification of all 161 case artifacts; safety 0 checked
all 161 cases in its full run. All 2,090 numeric output files for each policy are
byte-identical to the published baseline. The 483 structural output files are
byte-identical between policies. Existing 356 matrix cases, 252 metamorphic relations,
96 repeated refactors, typed refusals, handle identity, examples, toolchain checks
and all 14 warmed allocation paths passed; the computed-shift caller also remained
allocation-free. Expected negative compilation probes are reported separately
from successful library compilation.

The current serial commands are:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/check.ps1 -Safety 1
powershell -NoProfile -ExecutionPolicy Bypass -File scripts/check.ps1 -Safety 0
```

The first run included 160 structural cases; chain-65536 was added afterward, measured
and independently checked separately, and the complete 161-case artifact roster was
then verified. The current runner generates all 161 cases. Peak-memory reproduction,
from a PowerShell session after generation, uses:

```powershell
& ./scripts/measure-analysis.ps1 -InputFolder ./build/sparse-analysis-input `
  -OutputFolder ./build/analysis-memory -Case @('chain-4097','chain-16384',`
  'chain-65536','star-16384','grid-65x65','bipartite-fill-4480',`
  'hub-cycles-2048','k2r-8192') -Safety 1
```

The process harness has a per-case deadline, writes failure logs, and removes a
stale success report before starting. Each measured case also goes through the
same independent structural verifier. A single host and a finite family roster
establish neither universal speed nor portability beyond the supported runtime.
