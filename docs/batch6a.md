# Batch 6a: Cholesky-form entry for positivity proofs

State: completed; verified at safety 1 and safety 0 on 2026-11-20, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract is
[cholesky-form.md](cholesky-form.md).

## What was delivered, and why

`factor-cholesky!` and `cholesky-ready?` in `coalton-linalg/sparse/spd-ldlt`: a second
numeric entry over the existing analysis that computes the Cholesky factor G with
`P (A − shift I) Pᵀ = G Gᵀ` in Cholesky's own arithmetic.

A consumer proves a matrix positive definite by running a floating-point Cholesky on a
shifted copy, and the theorem it uses is stated for that arithmetic. `factor!` computes
L D Lᵀ with products of three numbers, so the theorem did not apply to it as written.
There were two ways to close that: a rounding analysis of the L D Lᵀ kernel, or an entry
that does what the theorem assumes. The second was chosen because it leaves nothing to
trust but the code's correspondence to six stated steps.

The entry is for the proof only. It overwrites the Symbolic's one factor and no solve
accepts the result.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`, with every earlier roster unchanged. For this batch:

- **The bound, entry by entry.** 107 positive definite matrices up to order 10: named
  cases, 81 seeded sparse matrices at three densities, 18 dense Gram matrices with
  irrational factors, an arrow, a tridiagonal, a stored zero, and shifts. For each of the
  1,567 stored factor entries examined, the exact difference between the matrix and G Gᵀ, in rational
  arithmetic, is within `(p + 2) u Σ |g||g|`, p being the products that entry received.
  The largest was 0.875 of that allowance.
- **Structure.** The permutation, column pointers and row indices equal those of `factor!`
  on the same pattern; the row counts tallied from the row indices equal the counts read
  from the dense structure; nothing nonzero lies outside the structure.
- **Exactness.** `[[4,2],[2,5]]` gives G = `[[2,0],[1,2]]` exactly, and a shift applied by
  the call equals the matrix stored with the shift subtracted, bit for bit.
- **Lifetime.** After the call `cholesky-ready?` is true, `factor-ready?` false, `solve!`
  refuses without writing and the condition estimate refuses; a refused call withdraws a
  standing solving factor and its estimate; `factor!` restores solving.
- **Refusals.** Four matrices with a pivot that is negative or exactly zero are refused as
  `FactorNonPositive` at the named column; wrong value count, non-finite shift and
  non-finite value by their own names.
- **Allocation.** A compiled Coalton caller making the proof call at a shift and then the
  solving call: zero bytes over 10,000 calls after 10,000 warm-up calls, at both policies.

**Planted source controls.** Each of these edits inside `factor-cholesky!` was applied
alone and the focused tests run. Every planted source compiled, every one failed the
tests, and the restored source passed.

| Planted defect | Caught |
| --- | --- |
| The multiplier carries a pivot, as in the L D Lᵀ form | yes |
| The diagonal is the pivot, not its square root | yes |
| A zero pivot is accepted | yes |
| The shift is ignored | yes |
| The proof factor is marked ready for solving | yes |
| The update skips the diagonal term | yes |
| A refused proof call leaves the previous factor ready | yes |
| The stored estimate survives the proof call | yes |
| A refusal names no column | yes |
| The entries of a column are divided by the pivot, not by its root | yes |

## Limits

- The library states and tests the arithmetic. Whether a theorem's hypotheses hold for a
  given matrix and shift, including its terms for underflow, is the caller's to establish.
- The row count t′ⱼ is a count of row j's entries left of the diagonal. The contract page
  says so because the column length is the easy number to reach for and is the wrong one.
- One host and finite rosters establish neither universal behaviour nor portability beyond
  the supported runtime.
