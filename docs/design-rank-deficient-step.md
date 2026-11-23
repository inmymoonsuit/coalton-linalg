# Design: a step on a rank-deficient matrix, and certified bases

State: a design for review. Nothing here is built. It is written for the consumer that
asked for it, whose terms are answered where they apply, and it follows
[dependent rows](dependent-rows.md), whose proofs it reuses.

## What is asked for

A nonlinear solver's matrix J is often rank deficient: a model with redundant conditions
has rows that depend on other rows. Its step is then the one a truncated singular value
decomposition gives. With singular values `σ₁ ≥ σ₂ ≥ …`, a relative floor φ and a
reference scale s, keep the singular values above `φ s` and solve in their subspace:

    x* = Σ over σᵢ > φ s of  vᵢ (uᵢᵀ r) / σᵢ .

That is the minimum-norm least squares solution of the truncated matrix, whether or not r
is in its range. s is J's own largest singular value unless the caller gives another, as
it does when one floor is to cut several coupled matrices alike.

Three things are wanted, each **holding as stated or refused by name**:

1. The step `x*`, with a bound on how far the vector returned is from it, and a
   certificate that no singular value of J lies in a band about `φ s`.
2. The supports of the left null vectors: for each dependent row, which of the rows it
   depends on have a multiplier that is clearly not zero and which clearly zero.
3. A basis of the right null space with a proved residual.

The singular values themselves are never computed. For a sparse J of some thousands of
rows they cost far more than the step.

## What dependent rows already gives

[The search](dependent-rows.md) splits the m rows into a kept set S of r rows and d rows
left out, and for each row t left out gives a vector `v_t` with 1 at t, minus its
multipliers at kept rows, and a proved bound `ε_t` on `‖Jᵀ v_t‖₂`. Order the rows so
that the kept ones come first, write `J_S` for them and `w_t` for the multipliers of row
t as a row vector. Then

    J = C J_S + E,      C = [ I ]      E = [  0  ]
                            [ W ],         [ E_T ],

W the d-by-r matrix of the rows `w_t`, and row t of `E_T` the residual `Jᵀ v_t`, of norm
at most `ε_t`. So

    ‖E‖₂ ≤ ‖E‖_F ≤ γ := ( Σ ε_t² )^½ ,

and γ is the [gap bound](dependent-rows.md#what-is-proved) the library already returns.

## The inequalities

All four are about exact singular values of the matrix as stored.

**(a) `σ_{r+1}(J) ≤ γ`.** `C J_S` has rank at most r, and the singular values of a sum
satisfy `σ_{i+j−1}(A + B) ≤ σᵢ(A) + σⱼ(B)` (Weyl's inequality for singular values). With
`A = C J_S`, `i = r + 1`, `B = E`, `j = 1`: `σ_{r+1}(J) ≤ 0 + ‖E‖₂ ≤ γ`.

**(b) `σ_r(J) ≥ σ_r(J_S)`.** `J_S` is J with rows deleted, so
`Jᵀ J = J_Sᵀ J_S + J_Tᵀ J_T`, a positive semidefinite matrix added to `J_Sᵀ J_S`. Adding
one cannot lower any eigenvalue (Weyl's monotonicity theorem), and the singular values
are the roots of those eigenvalues. No correction for C is needed.

**(c) A lower bound τ on `σ_r(J_S)`.** `σ_r(J_S)² = λ_min(J_S J_Sᵀ)`. The
[Cholesky-form factorization](cholesky-form.md) of `J_S J_Sᵀ − c I` that runs to
completion proves, by the theorem that page states, that the matrix factored is positive
definite once c is lowered by the page's rounding allowance, and so `λ_min > τ²` for the
τ that allowance leaves. This is the proof the consumer already uses for a full-rank
Gram matrix; here it is applied to the kept rows. The kept rows' Gram matrix is formed
on the divided matrix, as the search is.

**(d) Bounds on `σ₁(J)`**, needed only when the caller gives no reference scale. From
below, `‖J e‖₂ / ‖e‖₂` for any vector, and a few steps of the power method on `Jᵀ J`
make it tight; the value is a lower bound whatever the vector, and the rounding of one
product is bounded as [the product bound](dependent-rows.md#what-is-proved) is. From
above, `σ₁ ≤ (‖J‖₁ ‖J‖_∞)^½`, which costs one pass, or the same Cholesky proof on
`t² I − Jᵀ J`, which costs a factorization and is as tight as the t tried.

## The certificate

With a reference scale known to lie in `[s_lo, s_hi]`, both ends the caller's s when one
is given:

    γ < φ s_lo        and        τ > φ s_hi .

Then by (a) and (b) exactly r singular values of J are above every possible `φ s`, the
truncation keeps exactly r, and none lies strictly between γ and τ, an interval that
contains every possible `φ s`. **There is no further factor β**: the band is the gap
between the two proved bounds, and each is proved against exact arithmetic. What is not
covered is another computation's own rounding. A singular value decomposition computed
elsewhere agrees on which values are above the floor whenever its error is smaller than
the distance from its floor to γ and to τ, and this design proves nothing about that
error.

## The step

    x = J_S⁺ ( C⁺ r ).

- `y = C⁺ r` is a least squares problem in C, whose normal matrix is `I + Wᵀ W`. By the
  inversion lemma `(I + Wᵀ W)⁻¹ = I − Wᵀ (I_d + W Wᵀ)⁻¹ W`, so it costs one dense
  symmetric positive definite system of order d, the number of rows left out, with every
  eigenvalue at least 1.
- `x = J_Sᵀ (J_S J_Sᵀ)⁻¹ y` is the minimum-norm solve on the kept rows, with the factor
  the search already holds and one or two steps of [refinement](refinement.md).

Because C has full column rank and `J_S` full row rank, `(C J_S)⁺ = J_S⁺ C⁺`: x is
exactly the minimum-norm least squares solution for the rank-r matrix `B = C J_S`.

**It is not the answer of solving the kept rows alone.** Take J of two equal rows,
`[1; 1]`, and `r = [1, 3]`. One row is kept, `W = [1]`. The least squares answer is 2,
the mean. Solving the kept row gives 1: it never sees the residual on the row left out.
`C⁺ r = (1 + 3) / 2 = 2` does, and `J_S⁺ · 2 = 2`.

## How far x is from x*

`x*` is the same solve for `J_r`, the best rank-r approximation of J, and x is for B. The
two matrices have the same rank and are close: `‖J − B‖₂ ≤ γ` and
`‖J − J_r‖₂ = σ_{r+1}(J) ≤ γ`, so `‖Δ‖₂ ≤ 2γ` for `Δ = B − J_r`.

The pseudo-inverses of two matrices differ by (Wedin, 1973; the decomposition is an
identity and can be checked by multiplying out)

    B⁺ − A⁺ = − B⁺ Δ A⁺  +  B⁺ B⁺ᵀ Δᵀ (I − A A⁺)  +  (I − B⁺ B) Δᵀ A⁺ᵀ A⁺ ,     Δ = B − A.

Apply it to r with `A = J_r`, write `x* = A⁺ r` and `ρ = (I − A A⁺) r` for the part of r
outside the range, and take norms, using `‖B⁺‖₂ = 1/σ_r(B)`, `‖A⁺‖₂ = 1/σ_r(J)` and
`‖I − B⁺ B‖₂ ≤ 1`:

    ‖x − x*‖₂  ≤  2γ ( ‖x*‖₂ / σ_r(B)  +  ‖ρ‖₂ / σ_r(B)²  +  ‖x*‖₂ / σ_r(J) ).

Both `σ_r(B) ≥ σ_min(C) σ_r(J_S) ≥ τ`, since `Cᵀ C = I + Wᵀ W` has no eigenvalue below
1, and `σ_r(J) ≥ τ` by (b). So

    ‖x − x*‖₂  ≤  2γ ( 2 ‖x*‖₂ / τ  +  ‖ρ‖₂ / τ² ).

**The second term is the one for r not in the range**, and carries the square of 1/τ, as
least squares does. Neither `‖x*‖` nor `‖ρ‖` is known, and both are bounded from what
is. `x*` minimizes the residual of `J_r`, so `‖ρ‖ ≤ ‖r − J_r x‖ ≤ ‖r − J x‖ + γ ‖x‖`,
from the computed x. And `‖x*‖ ≤ ‖x‖ + ‖x − x*‖`, which puts the bound on both sides
of one scalar inequality; it is solved once, and has a solution when `4γ/τ < 1`. Where it
has none the call refuses. To this is added the rounding of the two solves, bounded by
their refinements' final residuals over τ.

## Supports of left null vectors

For a row t left out, its exact multipliers are `w_t = a_t J_S⁺`, unique because the kept
rows are independent. The computed ones `ŵ_t` differ by

    ŵ_t − w_t = g (J_S J_Sᵀ)⁻¹ ,      g = (a_t − ŵ_t J_S) J_Sᵀ ,

so `‖ŵ_t − w_t‖₂ ≤ ‖g‖₂ / τ²`, with `‖g‖` bounded from the computed product as a product
bound is. Call that δ. Against two thresholds the caller gives, `lo < hi`, an entry is

- **clearly not zero** when `|ŵ| − δ > hi`,
- **clearly zero** when `|ŵ| + δ < lo`,
- **exactly zero** when the kept row cannot be reached from t in the factor's graph, which
  is a property of the pattern and needs no bound,

and otherwise the call refuses, naming the first entry it could not class. "Zero" means
below `lo` for the matrix as stored. A multiplier that is exactly zero in the model and
1e-17 in the stored matrix is classed zero only because 1e-17 is below `lo`; no
computation on rounded data proves more.

## A basis of the right null space

Dependent rows applied to `Jᵀ` finds the columns of J that depend on the others, n less
that search's own kept count, which is r when both certificates hold, and for each a
vector z with 1 at that column,
minus its multipliers at kept columns, zero elsewhere, and a proved bound on `‖J z‖₂`.
The n − r vectors are independent, since each has a 1 where the others have 0, and form

    Z = [ −Mᵀ ]     kept columns
        [  I  ]     columns left out,

a basis that is sparse where the multipliers are, and **not orthonormal**. An orthonormal
one is Q from a Householder factorization `Z = Q R`: `σ_min(Z) ≥ 1` for the same reason
as for C, so `‖R⁻¹‖₂ ≤ 1` and `‖J Q‖_F ≤ ‖J Z‖_F`, to which the factorization's own
rounding is added. Q is dense, n by n − r.

The same floor and reference scale apply: the basis is certified when the certificate
for `Jᵀ` holds, which is the certificate above with the transposed search's γ and τ.

## Refusals

By name, and distinct:

| Refusal | When |
| --- | --- |
| Band not separated | `γ ≥ φ s_lo` or `τ ≤ φ s_hi`: a singular value may lie at the floor |
| Below the method's resolution | `φ s` is below what the search can resolve, about 1.5e-8 of a row's size: no certificate could be issued at this floor by this method |
| Bad reference scale | s not finite or not positive |
| Bound not reached | `4γ/τ ≥ 1`, or a solve whose refinement did not converge: the step is returned by no call |
| Support undecided | an entry between the two thresholds, with its row and position |

The existing refusals of the search and of the factorizations pass through under their
own names.

## Cost, and what will be timed

With d rows left out: d null vectors, each one sparse triangular solve and up to four
refinement steps; W, d by r, held densely or by rows; one dense factorization of order d;
then two sparse solves a step. For d in the hundreds the dense system is cubic in d and
is expected to be the smaller part; that is a prediction and will be measured, on the
patterns the consumer names: a lattice of n by n rectangles with every shared edge
carrying a coincident pair and a duplicated horizontal or vertical condition, at n = 15
and 21, and a wheel with a redundant distance on every spoke.

## Questions for the reviewer

1. **Which basis?** The sparse Z with a bound for each vector, or the dense orthonormal Q
   with one bound for all, or both. Q costs `n (n − r)²` operations and `n (n − r)`
   numbers.
2. **The bound on `‖x − x*‖`** contains `‖ρ‖/τ²`. For an r far outside the range and a τ
   near the floor it can exceed `‖x‖` itself. Is a step returned with such a bound of
   use, or should a bound above a caller's fraction of `‖x‖` be a refusal?
3. **The reference scale as an interval.** When the caller gives s it is taken as exact.
   When it does not, `[s_lo, s_hi]` from (d) widens the band that must be clear by
   `s_hi / s_lo`. Is a default built on the cheap upper bound acceptable, at the price of
   refusing some floors a tighter bound would pass, or is the factorization for the tight
   one worth its cost here?
