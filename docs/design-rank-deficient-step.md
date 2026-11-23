# Design: a step on a rank-deficient matrix, and certified bases

State: a design, reviewed by the consumer that asked for it and agreed with the terms
recorded at its end. Nothing here is built, and one measurement is owed before any of it
is. It follows [dependent rows](dependent-rows.md), whose proofs it reuses.

## What is asked for

A nonlinear solver's matrix J is often rank deficient: a model with redundant conditions
has rows that depend on other rows. Its step is then the one a truncated singular value
decomposition gives. With singular values `σ₁ ≥ σ₂ ≥ …`, a relative floor φ and a
reference scale s, keep the singular values above `φ s` and solve in their subspace:

    x* = Σ over σᵢ > φ s of  vᵢ (uᵢᵀ r) / σᵢ .

That is the minimum-norm least squares solution of the truncated matrix, whether or not r
is in its range. s is J's own largest singular value unless the caller gives another, as
it does when one floor is to cut several coupled matrices alike. **The caller gives it as
an interval `[s_lo, s_hi]`**: a largest singular value over several matrices is known to
it only through bounds, and a scale taken as exact would claim more than is known.

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
at most `ε_t`. **W is exactly the stored doubles of the computed multipliers**, the ones
the null vectors carry. Then `J = C J_S + E` is an identity about the stored matrices,
with E the exact residuals of those vectors, which are what the bounds bound; and C is
the matrix actually used when the step is solved for. So

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

**(c) A lower bound τ on `σ_r(J_S)`.** Two proofs are available, and which is used is
what the measurement below decides.

- *On the squared matrix.* `σ_r(J_S)² = λ_min(J_S J_Sᵀ)`. The
  [Cholesky-form factorization](cholesky-form.md) of `J_S J_Sᵀ − c I` that runs to
  completion proves, by the theorem that page states, that the matrix factored is positive
  definite once c is lowered by the page's rounding allowance, and so `λ_min > τ²` for the
  τ that allowance leaves. It costs one factorization on the pattern the search already
  analyzed. Its limit is that the allowance is a multiple of the unit roundoff times the
  size of `J_S J_Sᵀ`, so it cannot prove a τ much below `(n u)^½ σ₁`: about 5e-7 of
  `σ₁` at a thousand rows.
- *On the matrix itself.* [The rank certificate from the inertia](quasidefinite.md#a-rank-certificate-from-the-inertia)
  proves that every singular value of a matrix exceeds `(s t)^½` from one signed
  factorization of `[[t I, J_Sᵀ], [J_S, s I]]` and its measured backward bound, without
  forming `J_S J_Sᵀ`. Its resolution is that backward bound, a multiple of the unit
  roundoff times the matrix's own size and the growth measured, so it can reach a τ far
  below the first proof's where growth allows. It costs an analysis and a factorization of
  the larger matrix.

**(d) Bounds on `σ₁(J)`**, needed only when the caller gives no reference scale. From
below, `‖J e‖₂ / ‖e‖₂` for any vector, and a few steps of the power method on `Jᵀ J`
make it tight; the value is a lower bound whatever the vector, and the rounding of one
product is bounded as [the product bound](dependent-rows.md#what-is-proved) is. From
above, `σ₁ ≤ (‖J‖₁ ‖J‖_∞)^½`, which costs one pass, or the same Cholesky proof on
`t² I − Jᵀ J`, which costs a factorization and is as tight as the t tried.

## The certificate

With the reference scale in `[s_lo, s_hi]`, the caller's interval or the one (d) gives
for J itself:

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
their refinements' final residuals over τ. **Each of those is a proved residual**, formed
with its [product bound](dependent-rows.md#what-is-proved) as the null vectors' are, and
not a residual computed in floating point and taken at its word.

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

**The supports give the dependent sets.** A caller that wants the groups of rows that
depend on one another merges rows whose supports overlap. Each vector here is the
fundamental circuit of its row with respect to the kept set, and a reduced echelon basis
of the left null space is the set of fundamental circuits with respect to another
independent set, the rows that are not its pivots. Merging by overlap gives, from
either, the connected components of the rows as a matroid, which do not depend on the
independent set chosen. That was checked in exact rational arithmetic on four hundred
small matrices, a hundred and twenty of them with two or more components, before being
written here. What a floating-point computation adds is only the classing of entries
above: a multiplier wrongly classed joins or parts two groups, which is why an entry
that cannot be classed is a refusal.

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
rounding is added. Q is dense, n by n − r. **Both are offered**: Z always, with a bound
for each vector, being cheap and certified as it stands; Q on request, with one bound
for all, for a caller that projects onto the basis and steps along unit directions. Its
`n (n − r)²` operations are less than such a caller spends re-solving once a direction.

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
own names. **A bound that is reached is returned however wide it is.** How wide a bound
is acceptable is the caller's to say, and a returned bound holds as stated; the library
takes no argument for it.

## Cost, and what will be timed

With d rows left out: d null vectors, each one sparse triangular solve and up to four
refinement steps; W, d by r, held densely or by rows; one dense factorization of order d;
then two sparse solves a step. For d in the hundreds the dense system is cubic in d and
is expected to be the smaller part; that is a prediction and will be measured, on the
patterns the consumer names: a lattice of n by n rectangles with every shared edge
carrying a coincident pair and a duplicated horizontal or vertical condition, at n = 15
and 21, and a wheel with a redundant distance on every spoke.

## The measurement owed before code

The floor a nonlinear solver cuts its step at can be as low as 2e-8 of `σ₁`, and two
limits of this design sit near that figure: the search resolves about 1.5e-8 of a row,
and the first proof of (c) cannot certify a τ much below 5e-7 of `σ₁` at a thousand rows.
For rows that are redundant exactly, γ is near the unit roundoff times the matrix and is
not the concern. τ is. A refusal there is safe, since it costs time and not agreement,
but how often it would happen decides which proof (c) is built on.

So before any of this is built, on the two timing patterns and on one well-constrained
matrix, at floors of 2e-8 and of 1e-6:

- `γ / s` and `τ / s`, with τ from each of the two proofs;
- how often each refusal would fire, and for "band not separated" which side fails;
- how often the width of the default scale interval alone is what fails it.

If τ from the squared proof is the binding side at 2e-8 on ordinary matrices, the step is
built on the proof on the matrix itself, and this note is changed to say so.

## Decided in review

- **The reference scale is an interval the caller gives**; J's own, from (d), is the
  default. The default uses the cheap upper bound and the power method's lower bound. If
  the measurement shows its width refusing real matrices, the tight bound is a retry the
  caller can make, not a cost every call pays.
- **Both bases**: the sparse one always, the orthonormal one on request.
- **A bound that is reached is returned**, however wide; "bound not reached" stays for
  `4γ/τ ≥ 1` and for a solve whose refinement did not converge.
- **The supports are per entry**, in three classes against two thresholds the caller
  gives, and are what the dependent sets are formed from.
- **No factor β**: the band is the gap between two proved bounds.
