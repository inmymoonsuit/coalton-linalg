# Changelog

Changes a caller can meet. Each release's evidence is in the batch record it links.

## 0.2.0 — 2026-11-19

Added:

- `coalton-linalg/sparse/gram`: fixed-pattern AᵀA and AAᵀ. Analyze a checked CSC matrix once,
  then `refill!` caller-owned values without allocating. The output is the lower triangle with
  every diagonal, the pattern SPD analysis takes. [Contract](docs/gram.md) ·
  [record](docs/batch3c.md).
- `matvec!` and `transpose-matvec!` in `coalton-linalg/sparse/csc`, writing into caller storage.
  [Record](docs/batch3b.md).
- `analyze-with-limit` and `analyze-csc-with-limit` in `coalton-linalg/sparse/spd-ldlt`, taking
  an explicit budget of total factor entries. [Record](docs/batch3a.md).
- A compiled normal-equations example, `examples/normal-equations.ct`.

Changed, and able to break a caller:

- **SPD analysis no longer refuses dimensions above 4096.** It refuses instead when the factor
  would exceed an entry budget, 1,048,576 total factor entries by default. A matrix below 4096
  whose factor fills beyond that budget, which 0.1.0 analyzed, now needs an explicit larger budget.
- **`LinalgError` gained two constructors**, `InvalidAliasing` and `AnalysisStorageExceeded`.
  An exhaustive match on the exception family must add them.
- `maximum-analysis-size` is now the representability bound for dimensions and budgets, not 4096.

Unchanged: the elimination order, every numeric factor and solve result on the published test
rosters (byte-identical output files), and the positive-pivot contract.

## 0.1.0 — 2026-11-19

The first published snapshot: validated CSC and COO construction, transpose, permutation,
allocating matrix-vector products, stable norms, and sparse SPD LDLᵀ with reusable symbolic
analysis, typed factor and solve statuses, and typed construction errors.
