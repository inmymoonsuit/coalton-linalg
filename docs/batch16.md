# Batch 16: the range of conjugate gradients and of dependent rows

State: completed; verified at safety 1 and safety 0 on 2026-11-22, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contracts it touches are
[conjugate-gradient.md](conjugate-gradient.md), [dependent-rows.md](dependent-rows.md),
[gram.md](gram.md) and [rosters.md](rosters.md).

## What this is

[Batch 15](batch15.md) asked every routine what it does with the same problem at every
scale, and measured three places where the answer was a refusal although the result could
be held. They were one cause seen three times: a method that squares its data reaches
half the exponent range. This batch removes it from two of them and gives the refusal
that is left a name of its own.

The remedy is the one the references give for it, and is the same in both: divide the
data by a power of two before the method squares it, and multiply back what has units.
Dividing by a power of two changes no significand, so nothing is rounded differently.

## What changed

- **Conjugate gradients.** `cg!` and `pcg!` carry the iteration on the residual divided by
  the power of two at the first residual's largest entry. The preconditioned residual, the
  direction and its image are divided with it, each being linear in the residual; the
  iterate and the norms reported are not. Inner products are of numbers near 1 whatever
  the size of the right-hand side, and the operator is applied to a direction near 1
  whatever its own size.
- **Dependent rows.** `find-dependent-rows!` makes its search on the matrix divided by the
  power of two at its largest entry, the threshold divided with it. The analysis keeps
  the values as given beside the divided ones, so `kept-values!` hands back what the
  caller passed; the bounds are formed on the divided matrix and multiplied back when
  read, raised by the smallest positive number where the product falls below the normal
  range. `transpose-product-bound` divides its matrix and its vector the same way.
- **`RowsUnderflow`**, a new value of `RowStatus`: under a threshold of zero, a row
  holding a nonzero entry that is too small beside the matrix's largest entry for its
  square to be held. Version 0.13.4 reported it as `RowsNonFinite`.
- **`gram:check-csc!`**, the source pattern check of `refill-csc!` by itself, so that the
  search can hold a matrix to its plan and then refill from divided values.
- One helper, `binade!`, in the internal float package: the power of two at a magnitude
  and its reciprocal, written into an array. Both routines use it. Its loops end whatever
  they are given, so a value that is not finite needs no separate test in its callers: it
  stays not finite when divided, and the check that already followed refuses it.

## What it measured

[The scale-sweep roster](rosters.md#the-scale-sweep-roster), exponents of two, before and
after. "Holds" is every output bit for bit at its scaled value.

| Entry | 0.13.4 | Now |
| --- | --- | --- |
| `cg!`, right-hand side scaled | holds -800 to 400; `CgNonFinite` from 510 | holds -800 to 1020 |
| `cg!`, operator and right-hand side scaled | holds -512 to 300; `CgNotPositive` at -520 and below, `CgNonFinite` from 400 | holds -800 to 1020 |
| `pcg!`, operator, preconditioner and right-hand side scaled | holds -800 to 400; `CgNonFinite` from 510 | holds -800 to 1000; within tolerance at 1020 |
| Dependent rows, threshold 2^-20 of the entries | holds -400 to 400; refused outside | holds -800 to 1000; succeeds at -1000 and 1020 with a bound that is looser there |
| Dependent rows, threshold zero | holds -512 to 510; refused outside | holds over the whole sweep |

Below -800 the conjugate gradient entries have no input: their tolerance, 2^-40 of the
right-hand side, is not a normal number there. At 2^1020 the dependent-rows bound on a
transposed product is infinite because its value, 19 times 2^1020, cannot be held.

**On data that worked before, nothing moved.** The ninety-six conjugate gradient solves
that the tests recompute step by step agree bit for bit, as they did, and so do the 118
exact dependent-rows cases, the 364 certificates and the allocation measurements, which
are still zero bytes for a compiled caller.

## Two things the planted defects changed

The first run of the planted defects missed two, both in conjugate gradients, and each
was a finding about the code and not about the plant.

- **A refusal that rejected a solvable problem.** The first version ended `CgNonFinite`
  when a running residual's norm, multiplied back, was not finite. No test reached it, so
  the plant that removed it changed nothing. Looking for a test that would reach it gave
  a system it was wrong for: the operator diag(1, 1000) under the right-hand side
  (1.5e308, 5e307) triples its residual in the first step, to an entry of -4.5e308 that no
  double holds, and is solved exactly in the second. The refusal is removed. The norm of
  such a step is recorded as infinity and the call goes on, and that system is a test.
  The same search showed the iterate's update multiplying a direction back to a number
  that cannot be held before applying a step length that would have brought it down; the
  two factors, the step length and the power of two, are now multiplied first.
- **A defect the method hides.** Choosing a new power of two when a residual is formed
  again, while the step is still multiplied back by the old one, makes every later step
  the wrong length. The call still converges: the confirmation on a residual formed from
  the operator catches each bad result and starts again. Every assertion in the suite
  held. What gives it away is the cost, so the test of one wrong product now also asserts
  two replacements and at most sixteen steps for an order of eight.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`.

- **Red first.** Each limit above was measured on the earlier sources by the roster of
  Batch 15 before any of this was written; the table's first column is that run. The test
  that expected `CgNonFinite` from an operator and a right-hand side of 1e200, and the
  one that expected an infinite bound from two products of 1e310 that cancel, failed on
  the new code and were rewritten to say what is now true, each keeping a case that
  still overflows. The assertion that a threshold of 1e200 is refused was rewritten with
  them without being run first.
- **Direct tests.** A matrix at 2^-1000, 2^-600, 2^600 and 2^1000 gives the rows left
  out, the first vector and the values of the matrix at scale 1, bit for bit, under a
  threshold of zero and under one scaled with it. One row 2^600 times smaller than the
  rest is `RowsUnderflow` under a threshold of zero and is left out under one that can be
  held. The product bound is within its margin of the exact norm with either factor at
  2^600, 2^-600, and at opposite ends.
- **Planted source defects.** Twenty-four edits to the three changed sources, each
  applied alone, with the dependent-rows and conjugate gradient suites and both rosters
  run. Every planted source compiled, every one failed a test, and the restored sources
  passed.

| Planted defect | Number | Caught |
| --- | --- | --- |
| Conjugate gradients: the residual is not divided; the step, a running norm, or the norm of a residual formed from the operator is not multiplied back | 4 | yes |
| Conjugate gradients: a residual formed again chooses a new power of two; a step that is not finite in the iterate's units is applied; the direction is multiplied back before the step length is applied | 3 | yes |
| Dependent rows: the threshold, or the matrix, is not divided; the values handed back are the divided ones | 3 | yes |
| Dependent rows: a residual bound, or the gap bound, is not multiplied back; a bound multiplied back below the normal range is not raised | 3 | yes |
| The standalone bound: not raised below the normal range; always multiplied by the two powers' product; the matrix, or the vector, not divided | 4 | yes |
| Dependent rows: a row too small to square is left out under a threshold of zero; a positive threshold too small to square is searched with; a threshold beyond every row is refused; a threshold that is not finite is taken | 4 | yes |
| The power of two: not raised with its magnitude; its reciprocal not raised as it falls; it falls below the least normal number | 3 | yes |

  This is the second run. The first, of twenty-three, missed the two described above,
  and the code and the tests were changed for each before this one.

## Limits

- **What a square cannot hold is now a matter of sizes inside one matrix.** A row more
  than about 2^511 times smaller than the largest entry, or a positive threshold that
  small, is refused by name. A caller with such rows scales them first.
- **A Gram refill by itself still returns zeros for squares that underflow.** The search
  no longer meets it, since it refills from divided values, but a caller of `refill!` does.
- **The operator and the preconditioner of conjugate gradients are called with divided
  vectors.** For a linear one that changes nothing it computes. One that is not linear
  was never the operator the method assumes, and now sees arguments of another size.
- **Conjugate gradients on an operator of condition 1e400 are not rescued by this.** With
  diag(1e-200, 1e200) the running residual's square leaves the range in the divided
  arithmetic too, and the call ends `CgNonFinite`.
- **The factorization's refusal inside the dependent-rows search still has no test that
  reaches it**; with the search made on a matrix whose largest entry is near 1 it is
  harder to reach than before.
- **`binade!` takes up to 1023 steps of a loop for a magnitude at the end of the range**,
  once a call. It is not measured against the solve it precedes.
