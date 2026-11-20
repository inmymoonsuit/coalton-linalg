# Condition estimate from a prepared factor

Package: coalton-linalg/sparse/spd-ldlt. Source: [spd-ldlt.ct](../src/sparse/spd-ldlt.ct).
State: completed as Batch 3d; verified at safety 1 and safety 0. The
[batch record](batch3d.md) holds the measurements and the controls.

## What is estimated

After a successful `factor!(symbolic, values, shift)`, write M for the matrix the
factor represents: `A - shift I` in the caller's original ordering, as the computed
L and D reproduce it. The quantity is the **1-norm of the inverse of M**, the
largest absolute column sum of M⁻¹. M is symmetric, so this equals its infinity
norm.

- `estimate-inverse-norm1!(factor) -> ConditionStatus` computes a lower estimate of
  that norm by solves with the factor. It never forms an inverse.
- `inverse-norm1-estimate(factor) -> F64` is the latest successful estimate for the
  factor as it now stands, or -1 when none stands.
- `estimate-solve-count(factor) -> UFix` is the number of solves the latest
  estimate call made, counting one that was refused.
- `maximum-estimate-solves` is the bound on that count, 11.
- `shifted-norm1(symbolic, values, shift) -> F64` is the exact-in-structure 1-norm
  of `A - shift I` for values in the analyzed lower CSC order: each stored lower
  entry counts toward its own column and, off the diagonal, toward its mirror's.
- `reciprocal-condition1(norm, inverse-norm) -> F64` is `1 / (norm * inverse-norm)`
  when that product is positive and finite. Two zeros, the empty matrix, give 1.
  Anything else gives 0, including the -1 of a missing estimate.

Typical use, after `factor!` reports `FactorSuccess`:

```lisp
(let norm = (ldlt:shifted-norm1 symbolic values shift))
(let factor = (ldlt:symbolic-factor symbolic))
(when (== (ldlt:estimate-inverse-norm1! factor) ldlt:ConditionSuccess)
  (ldlt:reciprocal-condition1 norm (ldlt:inverse-norm1-estimate factor)))
```

## What the number means, and what it does not

**It is a lower estimate, not a bound that has been certified.** Each solve
applies M⁻¹ to a vector b and yields the true lower bound |M⁻¹b|₁ / |b|₁ on the
norm, up to the rounding of that solve. The estimate is the largest such value
seen. It can fall short of the norm, so the reciprocal condition can come out
**larger than the true one: a matrix can look better conditioned than it is, and
never worse** beyond rounding. On the tested rosters the estimate was between 0.57
and 1.00 of the exact norm; that is an observation about those matrices.

**It describes the shifted matrix that was factored**, not A, and not the matrix a
caller formed A from. For a Gram matrix G = BᵀB of a full-rank B, in exact
arithmetic the 2-norm condition of G is the square of B's.

**It is the 1-norm quantity, not a ratio of singular values.** For a symmetric
n-by-n M the two are related by cond₂(M) ≤ cond₁(M) ≤ n cond₂(M), which follows
from |M|₂ ≤ |M|₁ ≤ √n |M|₂ applied to M and to its inverse. A caller that needs
the spectral ratio itself needs a different routine; this one is not relabelled
as that.

**It is not a rank decision.** A matrix that is exactly singular can factor with
a tiny positive rounded pivot ([the factor contract](spd-ldlt.md) has the case);
the estimate is then very large, and deciding what counts as too large is the
caller's policy.

## Outcomes and lifetime

| ConditionStatus | Meaning |
| --- | --- |
| ConditionSuccess | The estimate is stored on the factor |
| ConditionInvalidFactor | No successful factor stands; no solve was made |
| ConditionNonFiniteResult | A solve or a sum was not finite; no estimate is stored |

The estimate belongs to one factor generation. Every `factor!` call on the owning
Symbolic clears it, whatever that call returns, and every estimate call replaces
it, so `inverse-norm1-estimate` is -1 from a `factor!` until the next successful
estimate. The empty matrix succeeds with an estimate of 0 and no solve.

`shifted-norm1` takes the values and the shift as arguments because the library
keeps neither. **Nothing ties the values passed to it to the values that were
factored**: the caller passes the same array and the same shift. It raises
InvalidDimensions for a wrong value count, NonFiniteInput for a non-finite shift or
value, and NonFiniteResult when a column sum overflows. It does not change the
factor or its validity.

All three operations use workspace owned by the factor, including the scratch
`solve!` uses. They are not reentrant with each other or with `solve!` and
`factor!` on the same Symbolic. A successful estimate allocates nothing after
warm-up on the measured runtime; typed exceptions may allocate.

## Cost

At most 11 solves, each O(n + nnz(L)), plus O(n) work between them. The observed
count on the measured systems is in the batch record. `shifted-norm1` is
O(n + nnz(A)).

## Method

The method is Hager's 1-norm estimator with Higham's refinements, a published algorithm,
written here from its description; [provenance](provenance.md#batch-3d-condition-estimate)
says what was and was not consulted. Its derivation:

Write f(x) = |M⁻¹x|₁. It is convex, so over the set |x|₁ ≤ 1 it takes its maximum
at a vertex ±eⱼ, and that maximum is the norm sought. At a point x let y = M⁻¹x
and let ξ hold the signs of y. Then f(x) = ξᵀy, and z = M⁻¹ξ (one more solve,
since M is symmetric) is a subgradient: f(w) ≥ f(x) + zᵀ(w - x) for every w.
Taking w = ±eⱼ shows that a vertex can only be promised an increase when
|zⱼ| > zᵀx.

1. Start at the uniform vector with every entry 1/n.
2. Solve for y and record |y|₁. Stop the ascent if it did not grow.
3. Solve for z. If the largest |zⱼ| is at most zᵀx, or j is the vertex just
   left, stop the ascent; otherwise move x to eⱼ and repeat from step 2, for at
   most five rounds.
4. For n > 1, solve once more for the vector b with alternating signs whose
   magnitudes grow evenly from 1 to 2. Its 1-norm is 3n/2, so 2|M⁻¹b|₁ / (3n) is
   also a lower bound; keep it if it is larger.

Step 4 is there because step 3 can stop at once. For two nearly parallel rows,
`[[1, 1-e], [1-e, 1]]`, the uniform start is almost the direction the inverse
amplifies least, the subgradient there promises nothing, and the ascent returns
about 1/2 against a norm near 1/e. The alternating vector is close to the
direction that matters and recovers it. The tests hold that family and a planted
control removes step 4.

The estimate is deterministic: no random start is used.
