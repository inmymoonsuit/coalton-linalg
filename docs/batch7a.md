# Batch 7a: a backward bound, and a rank certificate from the inertia

State: completed; verified at safety 1 and safety 0 on 2026-11-20, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract is in
[quasidefinite.md](quasidefinite.md).

## What was delivered, and why

`factor-backward-bound` in `coalton-linalg/sparse/spd-ldlt`: for the factor of the latest
successful `factor-signed!`, one number that bounds the Frobenius norm, and so the
2-norm, of `P K̃ Pᵀ − L D Lᵀ`.

A consumer of the signed factorization uses the pivot counts as a proof: that a matrix
`[[t I, Aᵀ], [A, s I]]` has exactly as many negative eigenvalues as A has rows says that
every singular value of A exceeds √(s t), and the computed pivots prove it once two
diagonal margins cover the factorization's backward error. On a system with a dense
column that proof replaced a dense factorization of A Aᵀ, and the consumer measured its
whole step at 1.80 ms against 876 ms.

The proof needs the backward error as a number with every constant accounted for. The
consumer was assembling it from the growth figure, a tally of the factor's longest row and
its own inflation for rounding, and its version had no term for underflow. **One
constant belongs in one place**, so the library now computes the bound in the
factorization that knows all its parts, and the contract states the certificate with its
conditions.

**A correction to Batch 7 came with it.** The entrywise bound of that batch was stated
without its condition. It is the standard result, which assumes no underflow; near the
bottom of the range it is false, and a test here exhibits a two-by-two matrix whose
second pivot is off by four million times that bound. The contract now says "barring
underflow" where it states the entrywise form, and the new figure includes the term.

## What the number is

    γ₍r+2₎ · trace(|L| |D| |Lᵀ|)  +  η · N (r + 2) (2 + 2 max|L|),   rounded upward,

with r the longest row of L left of the diagonal, N the order, u = 2⁻⁵³,
`γₘ = m u / (1 − m u)` and η = 2⁻¹⁰⁷⁴. Three choices, each from a measurement or a count:

- **The trace, not N times the largest entry.** `|L| |D| |Lᵀ|` is positive semidefinite, so
  its Frobenius norm is at most its trace, and the kernel was already accumulating its
  diagonal for the growth figure. The trace is never larger than N times the largest
  diagonal entry and is usually several times smaller.
- **The longest row, not the longest column.** The number of products an entry receives
  is bounded by its row's entries left of the diagonal.
- **No pivot term in the underflow part.** A first version carried the largest pivot as
  well. Working out how a planted control for it could be caught showed that none could:
  that contribution is below the relative term for the same pivot by nine hundred orders
  of magnitude. It was removed.

It costs one comparison for each entry of L and one pass over the order, in the signed
factorization only, and the measured calls still allocate nothing.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`, with every earlier roster unchanged. For this batch:

- **The bound, exactly.** In each of the 282 signed factorizations of the Batch 7 roster,
  and in its shifted and named cases, the exact Frobenius norm of
  `P K̃ Pᵀ − L D Lᵀ`, computed in rational arithmetic from the exported factor, is at most
  the value returned. The largest was 0.127 of it.
- **The formula, exactly.** In each, the value returned is at least the formula above
  evaluated in rational arithmetic on the exported factor, and exceeds it by no more than
  one part in 2³⁰ and eight η: it is the stated expression rounded up, not a looser one.
- **Underflow.** A matrix `[[a, b], [b, c]]` with a = 2⁻¹⁰⁶⁰, b = (2¹⁴ + 1) η and
  c = 5 · 2⁻¹⁰⁶⁰ has a multiplier of exactly 1 + 2⁻¹⁴, and its product with b is not a
  multiple of η. The second pivot is off by exactly 2⁻¹⁴ η, more than a thousand times
  the relative bound, and within the value returned, which is 27 η.
- **Lifetime.** −1 before any factorization, after `factor!`, after `factor-cholesky!` and
  after every refusal; a positive number after a signed factorization.
- **The certificate.** The procedure of the contract page, written out in the tests, on
  342 cases: named cases with a known answer, dependent rows at four thresholds, a nearly
  dependent pair on both sides of its smallest singular value, a threshold exactly at a
  singular value, and 54 seeded sources at six thresholds each, from 2⁻⁴⁰ to 40. The truth is the
  exact inertia of `A Aᵀ − s t I` in rational arithmetic. **No certificate was issued for a
  false statement**: 201 were issued, all true; 141 were withheld, all false; none was
  withheld though true.
- **Allocation.** The compiled caller of Batch 7, refill, signed factorization and two
  solves: zero bytes over 10,000 calls after 10,000 warm-up calls, at both policies.

**Planted source controls.** Each of these edits to `spd-ldlt.ct` was applied alone and
the focused tests run. Every planted source compiled, every one failed the tests, and the
restored source passed.

| Planted defect | Caught |
| --- | --- |
| The bound leaves out the underflow term | yes |
| The bound takes the largest diagonal entry in place of the trace | yes |
| The bound takes no account of how long a row of L is | yes |
| The bound measures the longest column in place of the longest row | yes |
| The bound is not rounded up | yes |
| The underflow term ignores the size of the multipliers | yes |
| The underflow term is counted once, not once for each entry of a row | yes |
| The positive definite entry leaves a signed factor's bound standing | yes |
| The Cholesky-form entry leaves a signed factor's bound standing | yes |
| The Cholesky-form entry leaves a signed factor's growth standing | yes |
| The accessor returns the growth | yes |

The tenth control is about Batch 7's growth figure, not this batch's bound: nothing had
tested that the Cholesky-form entry withdraws it. The test was added here.

## A consumer's figures

Reported by the application that uses the certificate, on its own fixtures, with the
bound assembled on its side before this batch:

- A stored under-defined hub of 100 spokes, 352 rows over 402 columns, at t the largest
  entry: growth 243.8, and one factorization for each proof once the first margin is
  learnt.
- The same system with one row repeated is rank deficient and is not certified; a planted
  proof that ignored the pivot counts certified it, and was caught.

## Limits

- The bound is on the factorization's backward error, not on a solve's. A solve's residual
  has the entrywise bound of Batch 7, barring underflow.
- The certificate is sufficient, not necessary: a true statement within the margin of
  false is not certified. On this roster none was withheld; a roster with singular values
  within rounding of the threshold would show some.
- The certificate proves a statement about the matrix whose values were passed, as
  floating-point numbers. Whether those values represent what the caller meant is the
  caller's.
- Overflow is not bounded: a non-finite intermediate is refused by the factorization, and
  a bound too large to represent is returned as infinite.
- One host and finite rosters establish neither universal behaviour nor portability beyond
  the supported runtime.
