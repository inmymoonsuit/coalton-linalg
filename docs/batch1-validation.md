# Batch 1: numerical validation

This page describes the tests and their acceptance rules. Actual run results belong
in generated JSON reports; implementing a test is not evidence that it has passed.

## Seeded exact-input lane

`tools/stress.py` imports the existing independent Python verifier, never the Coalton
implementation. It generates 336 cases by default: 12 seeds, seven mathematical
families, and four variants of each family. Every matrix has at most 12 rows.

| Family | Construction and purpose |
| --- | --- |
| Diagonally dominant | Sparse signed dyadic off-diagonals; each diagonal exceeds its off-diagonal absolute row sum by at least 2. |
| Disconnected | The same dominance construction confined to separate blocks, including isolated vertices. |
| Known LDL product | Form a symmetric matrix with exact rational arithmetic from a small unit-lower matrix and positive dyadic D. |
| Cancellation | An exactly cancelling input off-diagonal is explicitly stored as a structural zero; other factor entries remain nonzero. |
| Near-singular | A 2×2 matrix with a positive, representable margin from 2^-20 through 2^-31. |
| Zero diagonal | An isolated exactly zero diagonal gives a predictable singular refusal. |
| Indefinite | A well-separated indefinite 2×2 block among positive diagonal entries gives a predictable refusal. |

Each family has a base case, a deterministic simultaneous row/column permutation,
and scales of 2^-300 and 2^300. The RHS is constructed independently from a known
dyadic solution and transformed with the matrix. The default expected roster is
240 successes and 96 refusals.

The generator uses a specified 32-bit recurrence for choices; its manifest records
seeds, transforms, construction metadata, expected outcomes, forward-error test
budgets, and normalized-text hashes of the inputs. It refuses a destination
containing unrelated matrix names, rather than silently retaining extra cases.

Each emitted factor must satisfy the existing independent reconstruction,
permutation/structure, shifted-factor, refusal-column, and residual checks.
Additional checks require the known solution within a stated budget, invariant
success/refusal outcomes, permuted/scaled solution correspondence, scaled D,
and unchanged L under power-of-two scaling. They also check the input transformations
themselves. There are 252 nontrivial base/variant relations in the default run.

The forward-error budget is conservative and family-specific; it is not a
condition-number estimate or a general forward-error guarantee. Backward error
and forward error are reported separately. All cases stay far enough from
overflow/underflow in this pipeline that a nonfinite intermediate is a defect,
not an expected large-scale test outcome.

## Floating-point boundary lane

`tests/property-tests.lisp` runs through the public interface and has separate rules:

- Reuse one Symbolic for 96 changed numeric systems, alternating shifts/signs,
  explicit zeros, aliasing, and deliberate refusal followed by recovery.
  Independently multiply exported L/D factors and solve known RHS values.
- Show that a positive mathematical margin can disappear during input rounding:
  1 + 2^-54 rounds to 1 in binary64. The resulting stored singular matrix must
  be refused; its pre-rounding mathematical intent is not an input to the library.
- Exercise a scalar diagonal with shifts immediately below, at, and above it,
  and a least-positive subnormal pivot with a zero RHS.
- Observe 24 nearly singular Gram matrices assembled in floating point.
  Successful factorization or a nonpositive-pivot refusal is allowed here.
  A success must reconstruct the actual rounded matrix; a refusal must invalidate
  the factor and name a valid original column. Counts are reported without
  labeling them exact rank decisions.

The distinction is deliberate: the exact-input lane demands independent agreement
where it is meaningful; the boundary lane checks the numerical contract without
pretending exact arithmetic dictates every rounded pivot decision.

## Commands and controls

Run these serially, with an already installed Python/SBCL toolchain:

```powershell
python .\tools\stress.py generate .\build\stress-input
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\check.ps1 -Task Fixtures -InputFolder .\build\stress-input -OutputFolder .\build\stress-safety1 -Repeats 1
python .\tools\stress.py verify .\build\stress-input .\build\stress-safety1 --json .\build\stress-safety1\verification.json
python .\tools\stress.py selftest .\build\stress-selftest
```

Repeat the fixture run with `-Safety 0` and a separate output folder. The parent
check runner may automate these steps; the independent CLI remains available.

The stress verifier's self-test uses 28 exact oracle-output cases and six planted
controls. Five corrupt artifacts or inputs: a wrong permuted solution, wrong scaled
solution, wrong scaled diagonal, modified transformed input, and missing output.
The sixth compares expected pivots [1, 2^-100] with corrupted pivots [1, 2^-99]:
the former whole-vector comparison accepts the doubled tiny pivot, and the new
componentwise diagonal check must reject it. The evidence explicitly records both
the former acceptance and the new rejection. Solution comparisons retain their
normwise criterion.

Every diagonal scaling relation now checks each finite positive pivot relative to
that same pivot, with no absolute floor or scale taken from larger pivots. Thus a
large diagonal cannot hide a relative error in a small one.

The self-test restores inputs/outputs afterward and regenerates the seed set in a
separate directory to compare both file roster and contents. These controls
supplement the original verifier's ten artifact-corruption controls.

No third-party numerical implementation or published algorithm listing was read
or copied. The constructions follow diagonal dominance, matrix multiplication,
coordinate permutation/scaling, and the library's stated equations.
