# Batch 20: least squares on an operator

State: completed; verified at safety 1 and safety 0 on 2026-11-23, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract is [lsqr.md](lsqr.md) and the
design [design-iterative-least-squares.md](design-iterative-least-squares.md).

## What this is

The roadmap's item for the library's own users: least squares, and the minimum norm,
on a matrix that is available only as its product and its transpose's, with damping. No
consumer asked for it and it changes no existing call. The design was written first, with
the recurrences derived on the page, and the code was written from it.

## What was built

`coalton-linalg/lsqr`, one source file: `make-lsqr`, `lsqr!` and ten accessors.

- **The method** is the bidiagonalization of Golub and Kahan started from the right-hand
  side, with the damped least squares problem solved on the bidiagonal matrix as it
  grows: two products a step, six vectors of storage, and nothing of `Aᵀ A` formed.
- **What a status says is of norms formed from the operators.** The method's running
  estimates of the residual and of the gradient only decide when to form them. A
  confirmation that fails does not restart anything; the estimate that triggered must
  then fall by the factor it was out by, and one that was exactly zero is switched off.
- **Everything iterated is near 1.** The right-hand side and the operator's images are
  each divided by a power of two, the iterate is accumulated in those units and
  multiplied back once, and the same problem at another power of two is solved bit for
  bit the same.
- **It starts from zero and takes no starting vector**, which is what makes its answer,
  without damping, the minimizer of least norm whatever the shape and rank.

## Where building departed from the design

- **The operator norm estimate is not a bound.** The design had it at most the Frobenius
  norm of the operator in exact arithmetic. The first run of the exact tests showed 90
  for a matrix whose Frobenius norm is 71, after seven steps on six columns: past the
  steps exact arithmetic would take, the lost orthogonality lets the same directions be
  counted again. The contract now says so, and the test holds the estimate to the norm
  only for as many steps as the rank.
- **A ceiling on the damping.** A damped solution is of the order of the operator over
  the damping's square, so beyond a damping 2^500 times the operator the iterate could
  not be held in the iteration's units. That is refused as `LsqrOutOfRange`, which the
  design had only for the ratio of the right-hand side to the operator.
- **Without damping the residual's estimate is a magnitude.** The design's one formula,
  a root of a difference of squares, lost a residual below 2^-511 of the right-hand side
  to underflow. Without damping the estimate is now the one number the reduction left.
- **A diagonal entry of the reduction that vanishes ends the sequences**, as a vanished
  norm does.
- **The step count is bounded by the rank**, in exact arithmetic, not by the smaller
  dimension as the design said; the tests hold it to the rank and four more.
- **After `LsqrNonFinite` x is the iterate at which the value arose**, the last one
  completed where an operator's image was at fault, and untouched where the call ended
  before its first step.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`.

- **Against exact arithmetic.** 120 problems of small integers on twelve shapes, tall,
  square and wide, of full rank and not, two with a row and a column of zeros, at
  dampings 0, 1/4, 1 and 8, and with a right-hand side in the operator's range where
  there is no damping. For each: the two norms reported are the exact norms of the exact
  residual and gradient of the vector returned, to rounding; the status names a tolerance
  that those norms meet; and the vector is within 5 × 10⁻¹¹ of the exact minimizer, of
  minimum norm where there is no damping, against a tolerance of 2⁻³⁰ on the gradient.
  The most steps taken is 7, and none takes more than its rank and four.
- **The minimum norm.** Where there is no damping, the exact part of the vector outside
  the range of the transpose is at most 3 × 10⁻¹⁵ of the vector.
- **The three monotone quantities**, exactly, on all 466 iterates: the damped objective
  does not rise, the norm does not fall and the distance to the minimizer does not rise,
  by more than a part in 2³⁰, a distance already at rounding level not being compared.
  Fifty iterates were also obtained as the result of a call with that step limit, and
  are the same bit for bit.
- **An independent recomputation.** The iteration written again in the test language,
  with the library's two products, reproduces on every one of those problems and on
  the cases below that use a stored matrix: the outcome, the vector, the steps, the
  confirmations, both histories, both reports and the estimate, bit for bit.
- **Every status**, each by a case made for it: sizes, tolerances and dampings refused
  before anything is read or written; a right-hand side and a first image that are not
  finite, x untouched; an image that is not finite in the first step and in the second,
  x the last iterate completed; sizes the format cannot hold together, three ways, and a
  damping above the ceiling, x untouched; a right-hand side of zero, with no product
  made; one the transpose sends to zero; sequences that end where the arithmetic is
  exact, and where it is not, on five one-column problems of which more than one ends at
  a gradient a rounding away from zero; a limit of zero and of one. A report is withdrawn
  by a refused call.
- **A confirmation that fails.** One wrong product, in the second step, under a gradient
  tolerance and under a residual tolerance: the estimates then meet the tolerance and the
  true norms do not, by a factor of more than a thousand. The call reports neither
  tolerance met, reports the exact norms of what it returns, and makes two confirmations
  in forty steps, not one a step. With the product wrong by a part in 2¹⁴ the gradient's
  estimate is within the tolerance at 22 of forty steps, and four confirmations are made.
  An estimate that is exactly zero at three steps running costs one confirmation.
- **Vectors and norms far from unit size.** An operator of diag(1, 2^-200) under
  (1, 2^-600): a left vector of 2^-600 before it is normalized and a gradient of 2^-800,
  each measured as that number and not as the zero its square rounds to. A gradient too
  small to hold in the caller's units does not meet a tolerance of zero and meets any
  other. A damping 2^480 times the operator is solved to a part in 2⁴⁰ of the exact
  minimizer.
- **Scale.** Four problems, each at eight combinations of powers of two for the
  right-hand side and for the operator with its damping, from 2^-900 to 2^900: the same
  outcome, steps and confirmations, and the vector, both reports, the estimate and both
  histories each the scaled value bit for bit. [The scale-sweep
  roster](rosters.md#what-is-declared) has three entries for the call, exact from 2^-800
  to the top of the range.
- **Aliasing.** The right-hand side and the solution as one array give the same vector,
  here and in [the aliasing roster](rosters.md#the-aliasing-roster).
- **Allocation.** A compiled caller that solves a damped problem: zero bytes over 10,000
  calls after warm-up.
- **Planted source defects.** Thirty-two edits to the new source, each applied alone,
  with the package's suite and both rosters run. Every planted source compiled. Thirty-one
  failed a test on the first run. One did not: with the gradient's estimate triggering as
  before after a confirmation that failed, nothing failed, because in the wrong-product
  case the estimate met the tolerance once and then rose above it. A product wrong by a
  part in 2¹⁴ was added, after which 22 estimates in forty steps are within the tolerance
  and four confirmations are made; with that test the edit fails, and the restored source
  passes. Five tests were added while the edits were being chosen, because an edit had
  nothing to fail: the damping's ceiling and a damping just under it, the vectors far
  below unit size, a wrong product under a residual tolerance, the one-column problems
  whose sequences end inside the iteration, and the estimate of zero.

| Planted defect | Number | Caught |
| --- | --- | --- |
| What a status says: a tolerance is called met on its estimate, the gradient's and the residual's; an estimate triggers as before after a confirmation that failed, the gradient's and the residual's; a trigger of zero is not switched off; a norm that rounded to zero meets a tolerance of zero; the end of the sequences is not seen; it is reported as the limit; the confirmation at the end is not counted | 9 | yes, one after a test was added |
| The method's arithmetic: the damping is left out of the gradient formed; it is not divided as the operator is; the reflection against it is left out; the direction does not carry the one before; the gradient's estimate leaves out the entry below the diagonal; the residual's leaves out the damped term; the operator norm estimate is not in the caller's units; it leaves out the entries below the diagonal; a norm whose square is too small to hold is zero; a reflection's length is the root of a sum of squares | 10 | yes |
| The iterate: it is not multiplied back into the caller's units; it does not start from zero; a zero right-hand side leaves x as it was; after a value that is not finite it is left in the iteration's units | 4 | yes |
| Refusals: an image that is not finite is not noticed; the transpose's first image is not checked; the right-hand side is not checked; one wrong length is taken; one bad tolerance is taken; an infinite damping is not a bad damping; the sizes' ratio is not checked; the damping has no ceiling; a refused call leaves the earlier report standing | 9 | yes |

## What it costs

Two sparse problems that are not small, each solved here and by
[conjugate gradients](conjugate-gradient.md) on the equations the same vector satisfies,
to the same tolerance on the same quantity. The shortest of five runs, one process,
safety 1, on a machine otherwise idle; the products are the library's in-place ones.

**Tall, 39,999 by 20,000**: the identity over half the first differences, every singular
value between 1 and √2, a right-hand side of fixed pseudo-random numbers, stopped when the
gradient is 10⁻¹⁰ of its first value.

| | Steps | Time | Gradient reached, of its first value |
| --- | --- | --- | --- |
| This method | 14, one confirmation | 5.7 ms | 2.6 × 10⁻¹¹ |
| Conjugate gradients on `Aᵀ A x = Aᵀ b` | 14, one replacement | 6.0 ms | 2.6 × 10⁻¹¹ |
| This method, damping 1 | 11, one confirmation | 4.6 ms | 1.5 × 10⁻¹¹ |

The two vectors differ by 2 × 10⁻¹⁶ of their norm.

**Wide, 1,999 by 2,000**: the first differences alone, condition number about 1,270, a
right-hand side in their range, the shortest solution wanted, stopped when the residual is
10⁻¹⁰ of the right-hand side.

| | Steps | Time | Residual reached, of the right-hand side |
| --- | --- | --- | --- |
| This method | 1,999, one confirmation | 50 ms | 7 × 10⁻¹⁴ |
| Conjugate gradients on `A Aᵀ y = b`, then `x = Aᵀ y` | 1,999, one replacement | 57 ms | 1.1 × 10⁻¹³ |

The two vectors differ by 10⁻¹³ of their norm. The part of this method's outside the
range of the transpose, which for this matrix is its component along the vector of ones,
is 2 × 10⁻¹⁷ of its norm.

**What this shows and what it does not.** The method takes the steps conjugate gradients
take, as it should, being the same iteration in exact arithmetic, and no more time: a
step is two products either way. On the wide problem both end at the rank, 1,999 steps,
not sooner. Neither problem is conditioned badly enough to separate the two methods in
accuracy, and none was looked for that would.

## Limits

- **No bound on the error.** A gradient within tolerance says nothing of the distance to
  the minimizer without the smallest singular value, which this package does not
  estimate.
- **No reorthogonalization.** The step count is not bounded by the rank in floating
  point, and on a badly conditioned operator it is far above it.
- **No starting vector**, and so no restart from an earlier answer but by a second call
  on its residual, without damping.
- **The operator norm estimate is not a bound** in either direction once the method has
  run past the steps exact arithmetic would take.
- **With damping the residual's running estimate is a difference of squares.** It is a
  trigger only, and a damped problem is stopped on its gradient.
- **The exact tests are small**, up to 8 by 6. The two larger problems are timed and
  compared with another method, not checked against exact arithmetic.
- **That this method loses less to rounding than conjugate gradients on the normal
  equations** is the published comparison's finding and is not measured here beyond the
  two problems above.
