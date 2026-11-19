# Batches 1 and 2

State: completed, 2026-11-19.

Batch 1 established seeded independent numerical validation, separate tests at
floating-point boundaries, controlled allocation diagnosis, and an offline
toolchain lock. Batch 2 added reusable numeric and sparse matrix operations.

## Delivered

- 336 reproducible small exact-input cases across seven matrix families and
  four variants; 252 permutation/scaling relations.
- 96 repeated refactors plus explicit rounding/subnormal/shift tests and
  24 rounded-Gram observations.
- Stable Euclidean/infinity norms and a shared finite predicate.
- Owned canonical rectangular CSC, COO construction with explicit duplicate
  handling, transpose, permutations, A*x and transpose(A)*x.
- analyse-csc integration with the existing SPD factorization.
- Independent dense-reference primitive tests and two compiled Coalton examples.
- A content pin for 17 installed components/1,311 files, with actual source
  resolution checks for 28 loaded dependency systems.

The implementation and independent tests had separate owners. Review challenged
duplicate accumulation order, preservation of stored zeros, permutation direction,
ownership, empty dimensions, unused nonfinite vector entries, numerical-boundary
acceptance, and allocation attribution before or during implementation.

The shared-concern survey found a finite predicate in SPD LDL and lower-CSC
validation specific to that algorithm. The predicate now lives in numeric.
SPD-specific lower/diagonal checks remain with the factorization. The generic
matrix constructor permits rectangular storage. Existing Pattern callers remain
supported, and analyse-csc consumes the new structure through the existing analysis.

## Measured completion

The full serial check passed at speed3/safety1 from a new empty compiler cache,
then at speed3/safety0 using those freshly built dependencies. That cache contained
336 FASLs after the source build. The original 20 cases and all 336 seeded cases
passed at both policies. Their 74 and 2,016 respective artifact files were
byte-identical across policies.

For the 336 generated cases (at most 12 rows), maximum normwise backward error
was 1.1696176720330866e-16; maximum relative forward error was
7.450580596923828e-9. These measure different properties and are checked against
separate family-specific rules. They do not establish accuracy for arbitrary
large or ill-conditioned systems.

Review also strengthened the diagonal-scaling check to compare each pivot
componentwise; a planted tiny-pivot defect that the former global comparison
accepted is now rejected. Both artifact sets passed the strengthened check.

[Machine-readable results](batches-1-2-results.json) record the remaining counts.
[The numerical validation design](batch1-validation.md),
[primitive contracts](primitives.md), [toolchain protocol](toolchain.md), and
[allocation experiment](runtime-allocation.md) specify assumptions and limits.

The startup allocation anomaly was traced to background-runtime contamination of
the process-wide counter through controlled finalizer-thread experiments.
The disposable runners isolate that thread during bounded checks and prove
restoration on a planted failure. The library itself never changes thread state.

No Blender integration, dependency installation, commits, or publication occurred.
The original prototype remains unchanged. Scalable symbolic ordering and QR are
the next batches; CI setup and license selection remain release work.