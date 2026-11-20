# Batch 9a: scaling at the edges of the floating-point range

State: completed; verified at safety 1 and safety 0 on 2026-11-20, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. A correction to [Batch 9](batch9.md); the
contract is [scaling.md](scaling.md).

## What was wrong, and who found it

A held-back review of Batches 9 to 13 by the consumer that reviews this library ran
33,111 checks against exact rational arithmetic on matrices of its own. Refinement,
dependent rows, sparse LU and conjugate gradients passed every one. Three failed, all in
scaling, from two causes; reading the source against them showed two more of the same
kind.

1. **`scale!` was not exact where an intermediate left the range.** It formed `(r a) c`.
   For `[[1, 1e-310]]`, `equilibrate!` gives r = 1/2 and c = 2¹⁰²³; `r a` is subnormal and
   loses a bit before c restores the magnitude, and the value written was
   0.004494232837155998 where the exact product is 0.004494232837155776. The result is a
   normal number, so the contract's exception for values below the normal range did not
   cover it. *Found by the review.*
2. **The contract promised of every column what cannot hold of all.** A column whose
   entries fall below the normal range under the row factors cannot be brought to 1/2 by
   any factor there is, and the contract named that exception for rows only. *Found by
   the review.*
3. **`scale-symmetric!` had the same defect as `scale!`**, `(d a) d` in place of
   `(r a) c`. *Found here, on reading the source.*
4. **A column's factor could be wrong where its row-scaled entries were subnormal.** The
   column pass formed each row-scaled magnitude in floating point. One just under a power
   of two was rounded up to it, and the factor came out half of what the rule gives, so
   the column's largest scaled magnitude was 1/2 − 2⁻⁵³. One that rounded to zero made a
   column with a nonzero look empty: it kept the factor 1, and `scale!` then wrote zero
   for its only entry, which the largest factor would have kept exactly. The tests of
   Batch 9 stated that last behaviour as the expected one. *Found here.*

## What was changed

- **One product for both scaling calls.** The two factors are multiplied with each other
  first, and the value once by the result; where the factors' own product leaves the
  range they are applied in turn, which then moves the value one way only. It is exact
  for powers of two whenever the result is zero or normal, whatever the three sizes.
  [The contract](scaling.md#derivation) has the argument.

  An ordering of the two multiplications was considered first and is not enough: for
  a = 2.5, r = 2⁻¹⁰⁷⁴ and c = 2¹⁰²³ one order rounds and the other overflows. The reviewer
  showed that no order is exact for every pair a caller may choose.
- **The column pass works on magnitudes lifted by 2¹⁰²³**, so that none is rounded or lost
  where it matters, and gives a column beyond reach the largest factor.
- **The contract names the exception for rows and columns alike**, says how to find such a
  row or column, and warns that a threshold stated against 1 treats it as negligible.

**What a caller will see change.** For a matrix whose entries span less than about 2¹⁰²²
within every row, nothing: the factors and the scaled values are the same bits. Beyond
that, a column that used to keep the factor 1 or a factor one step short now gets the
largest factor, and values that used to be rounded are exact. `scale!` and
`scale-symmetric!` with factors that are not powers of two, which the contract makes no
claim about, may differ in the last bit from before, since the multiplications are
grouped differently.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`, with every earlier roster unchanged, the 246 cases of Batch 9 among them.
The new tests were written before the fix and failed on the source as it was.

- **The product, on factors a caller chose.** Eleven cases with a normal or zero result:
  the one where neither order is exact and its mirror; each order lossy alone, and
  overflowing alone; factors whose own product overflows, and underflows to zero; a
  subnormal value brought up; a subnormal factor; a zero factor. Each written value is the
  exact rational product. The symmetric form with 2.5 between factors of 2¹⁰²³ and 2⁻¹⁰⁷⁴
  likewise.
- **The whole range, exactly.** 150 matrices: the reviewer's two; a column whose
  row-scaled largest entry is just under 2⁻¹⁰²³; one whose only entry is lost even when
  lifted, the smallest number beside the largest; a row of subnormals; the smallest and the
  largest numbers in both rows; and 144 seeded
  matrices up to 8 by 8 whose entries range over the whole of the floating-point range,
  subnormals included. In each: every factor is a power of two; every scaled value whose
  exact product is zero or normal is that product, 1,496 values in all; every row and
  column with a nonzero has its exact largest scaled magnitude below 1, and at least 1/2
  unless its factor is 2¹⁰²³; a row or column with no nonzero has the factor 1.
- **Allocation.** The compiled caller of Batch 9, unchanged: zero bytes over 10,000 rounds
  after 10,000 warm-up rounds, at both policies.

**Planted source controls.** All twenty of Batch 9's, with the sites that moved brought up
to date, and ten for this correction. Each edit was applied alone to the final source and
the scaling tests run. Every planted source compiled, every one failed the tests, and the
restored source passed. The ten:

| Planted defect | Caught |
| --- | --- |
| The two factors are always taken in turn, as before the correction | yes |
| The two factors' own product is always used | yes |
| Factors whose product overflows are multiplied together all the same | yes |
| Factors whose product underflows to zero are multiplied together all the same | yes |
| The column pass forms the row-scaled entry before lifting it | yes |
| A column whose entries vanish under the row factors is treated as empty | yes |
| A column beyond reach is not held to the largest factor | yes |
| The column factor keeps the lift | yes |
| The symmetric form takes its factors in turn, as before the correction | yes |
| `scale!` takes its factors in turn, row first, as before the correction | yes |

## Limits

- **A product below the normal range is still rounded**, and can be rounded twice where the
  factors are applied in turn. Nothing is claimed of it but that it is finite.
- **A row or column beyond reach stays beyond reach.** The largest factor is the most a
  floating-point number allows; such a matrix wants a different formulation, not a better
  scaling.
- `equilibrate-symmetric!` is unchanged. Its own exception, a subnormal diagonal held to
  2⁵¹¹, was not among the reviewer's cases and was not revisited beyond the product it
  shares with `scale-symmetric!`.
- The reviewer's write-up lists what its review did not cover: allocation, the kept-rows
  certificate of dependent rows, the growth of LU storage, a preconditioner that is not
  positive definite, and dependent rows with hundreds of rows left out. Those remain as
  this library's own tests leave them.
- One host and finite rosters establish neither universal behaviour nor portability beyond
  the supported runtime.
