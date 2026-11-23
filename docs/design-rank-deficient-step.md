# Design: a step on a rank-deficient matrix, and certified bases

State: a design, reviewed four times by the consumer that asked for it. The measurement its
first review asked for is done and is reported below; it changed the design in four
places, each marked. The certificate, the kept vectors and the step are built, as
[kept rows, the certificate and the step](rank-deficient.md); the refresh of one vector
and the bases are not. It follows
[dependent rows](dependent-rows.md), whose proofs it reuses.

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
at most `ε_t`. **W is exactly the stored doubles of the multipliers that are used**, and
the bounds are proved for the vectors that carry them. Then `J = C J_S + E` is an
identity about the stored matrices, with E the exact residuals of those vectors; and C
is the matrix actually used when the step is solved for.

**The multipliers that are used are not all the ones computed** (changed by the
measurement). A dependent row depends on a few rows, three to four on the patterns
measured, and its computed vector has an entry at nearly every kept row the factor
reaches: up to 5.7 million entries where 6,500 are real, the rest below 10⁻¹³ with
nothing between there and 0.1. A multiplier below a drop tolerance is set to zero, and
the bound of the reduced vector is proved again by the product bound. That costs no
proof anything: the decomposition is about whatever W is used, and E is whatever that W
leaves. On the patterns the bound moved by under two percent. So

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

**(c) A lower bound τ on `σ_r(J_S)`.** Two proofs are available. **The step is built on
the second** (decided by the measurement): both proved the same τ to three figures on
every pattern, and they differ in cost, by a factor of a thousand where the kept rows'
Gram matrix is dense.

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
between the two proved bounds, and each is proved against exact arithmetic.

**The search's threshold is not the floor** (changed by the measurement). The search is
run at a rank tolerance well inside the gap, and the certificate is then asked for
whatever floor the caller needs. They are two arguments. Run at a floor of 2e-8 of the
scale the search is at its own resolution, and on two of the seven patterns it kept two
rows that depend on the others; both proofs of (c) then fail, as they must, and the
certificate refuses. Run at 1e-6 it kept the right rows on all seven and the certificate
held at 2e-8 and at 1e-6 alike. Searching at one tolerance and certifying at a lower
floor is safe by construction: a row left out whose distance lies between the two makes
γ at least `φ s`, and the certificate refuses. What is not
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
- `x = J_Sᵀ (J_S J_Sᵀ)⁻¹ y` is the minimum-norm solve on the kept rows. It is made on
  their augmented matrix and never on `J_S J_Sᵀ`, which a row touching every column
  makes dense; [the solve on the kept rows](#the-solve-on-the-kept-rows) says how, and
  what is proved of it.

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
| Bound not reached | `4γ/τ ≥ 1`; or a residual that is not finite; or a refinement whose proved bound was still at least halving at its last allowed correction: the step is returned by no call |
| Support undecided | an entry between the two thresholds, with its row and position |
| Kept vectors out of date | a kept vector's bound, proved again for the current values and refined, still does not come under the floor: the kept set is no longer right, and the caller searches again |

The existing refusals of the search and of the factorizations pass through under their
own names. **A bound that is reached is returned however wide it is.** How wide a bound
is acceptable is the caller's to say, and a returned bound holds as stated; the library
takes no argument for it.

## Keeping the vectors from one call to the next

(Added in the second review.) Finding a null vector costs 0.5 to 1.8 milliseconds, by
refinement, and a matrix with 1,640 redundant rows has 1,640 of them: 2 to 2.5 seconds a
call, where the solver calling this takes milliseconds a step. So the vectors are found
once and kept.

- A vector found for one matrix is kept with its row, its multipliers and its bound.
- At each later call, on a matrix of the same pattern, **its bound is proved again
  against the current values**, by the product bound, and nothing else is recomputed.
  Where the rows of a dependency do not change with the values, as rows that state that
  two things coincide or are level do not, the exact multipliers are constant, the kept
  vector is still the right one, and only the product is paid.
- A vector whose bound, proved again, fails the certificate is refined from where it
  stands, and only that one. A row whose own entries change at each call, a distance for
  instance, is the case.
- The kept set is the search's and is kept with them. If a re-proved bound cannot be
  brought under the floor the kept set is no longer right for these values: the call
  refuses by name, and the caller searches again.

γ is then the root of the sum of the bounds as proved for the current values, and τ comes
from one signed factorization of the current kept rows. Nothing a certificate rests on
is carried over from an earlier matrix; only the work of finding the vectors is.

W changes only where a vector was refined, and the system of order d below is factored
again only then.

## The system of order d

(Measured in the second review.) `I + W Wᵀ` is as sparse as W, and W, once the entries
that are rounding are dropped, has three to four entries a row. On the lattices the
product has about 2,000 entries at d = 812 and 4,000 at d = 1,640, and its sparse
factorization takes a tenth of a millisecond. The same matrix handed to the dense
symmetric factorization takes 0.08 s at 812 and about 1.3 s at 1,640. On the wheels one
dependency touches every spoke, the product is dense at d = 400, and either
factorization takes about 0.02 s.

So it is factored sparsely, by the [Gram product](gram.md) of W and the
[positive-pivot factorization](spd-ldlt.md), and it is not the cubic cost a dense W would
have made it.

## The solve on the kept rows

(Added in the fourth review, after a measurement made for it.) `x = J_S⁺ y` is the
solution of

    K [x; λ] = [0; y],        K = [[α I, J_Sᵀ], [J_S, 0]],

for any α > 0: the first row puts x in the row space and the second makes it solve. K is
as sparse as `J_S`. It is solved with a factorization of a nearby matrix and
[refined](refinement.md) against K itself:

    M = [[α I, J_Sᵀ], [J_S, −ε I]],        u ← u + M⁻¹ (b − K u),    u = [x; λ],  b = [0; y].

**Two shifts on one analysis.** The certificate factors `[[τ I, J_Sᵀ], [J_S, τ I]]` and
the solve factors M. They have one pattern and differ in the two diagonal blocks, so the
symbolic analysis made when the kept set was found serves both, and the solve is a second
numeric factorization at each call. The solve uses

    α = τ,        ε = τ · 2⁻¹⁰,        so   α ε = τ² · 2⁻¹⁰ ,

with **τ the bound the certificate has just proved for these kept rows and these values**.
Nothing is carried from an earlier call: a τ remembered from one is only where the
certificate's ladder starts, and is proved again before anything uses it.

**What one correction does, by blocks.** In exact arithmetic the error goes to `G e` with
`G = M⁻¹ (M − K)`, and `M − K = diag(0, −ε I)`. So G's first block column is zero, and
with `A = J_S J_Sᵀ`:

    e_λ ← G₂₂ e_λ ,     G₂₂ = α ε (A + α ε I)⁻¹ ,          ‖G₂₂‖₂ = α ε / (σ_r² + α ε) ,
    e_x ← G₁₂ e_λ ,     G₁₂ = −ε J_Sᵀ (A + α ε I)⁻¹ ,      ‖G₁₂‖₂ = max ε σ / (σ² + α ε) ≤ ½ (ε/α)^½ .

**G is not normal, and no single rate holds on the stacked vector**: `α ε / (σ_r² + α ε)`
is its spectral radius and not a norm of it. The statement is the pair: the λ error
contracts in the 2-norm by the first figure, and the x error after a correction is at
most the second times the λ error before it. With the shifts above and `σ_r > τ`, the
first is below `2⁻¹⁰` and so is the second. Both were checked on four hundred random
matrices against the measured ratio of errors at each step, not against the spectral
radius: on the 284 steps that were above rounding, measured over bound was at most 1.000.

**The computed factor.** The factor is of `M + E` with `‖E‖₂ ≤ β`, β its
[measured backward bound](quasidefinite.md#accuracy-growth-and-the-elimination-order).
The iteration matrix is then `(M + E)⁻¹ (M − K + E)`, and with `η = β ‖M⁻¹‖₂`,
`‖M⁻¹‖₂ ≤ 2/τ` for these shifts, a rate of at most `(ρ₀ + η) / (1 − η)` follows where
`ρ₀ = 2⁻¹⁰`. **That sum is a diagnostic and decides nothing.** β is a bound through the
trace of `|L| |D| |Lᵀ|` and grows with the order; on the largest pattern measured it
makes the sum 0.9 where the contraction measured is 2 × 10⁻⁴. A refusal on it would
refuse solves that reach rounding level in two corrections.

**What is proved is read from the residuals, and needs no rate.** For any pair
`(x̂, λ̂)`, however it was reached, write `ρ₁ = α x̂ + J_Sᵀ λ̂` and `ρ₂ = J_S x̂ − y`, and
`P = J_S⁺ J_S` for the projection on the row space. Then `P x̂ − J_S⁺ y = J_S⁺ ρ₂`, and
`(I − P) ρ₁ = α (I − P) x̂` because `J_Sᵀ λ̂` lies in the row space. So

    x̂ − J_S⁺ y  =  J_S⁺ ρ₂  +  (I − P) ρ₁ / α        exactly,

the two parts orthogonal, and

    ‖x̂ − J_S⁺ y‖₂  ≤  ( (‖ρ₂‖₂ / τ)² + (‖ρ₁‖₂ / α)² )^½ .

Each residual is formed with [its product bound](dependent-rows.md#what-is-proved), so
what enters is a proved upper bound on its norm for the exact product, and τ is proved.
This is the bound the step returns for its solve. The same argument serves the system of
order d: `Cᵀ C = I + Wᵀ W` has no eigenvalue below 1, so for a computed ŷ with
`s = Cᵀ r − (I + Wᵀ W) ŷ`, `‖ŷ − C⁺ r‖₂ ≤ ‖s‖₂`, and that error reaches x divided by τ.

**The stopping rule.** After the first solve and after each correction the proved bound
is formed, and the iterate with the smallest bound so far is the one kept.

- A correction that leaves the bound above half the smallest so far ends the refinement:
  it has stalled, at rounding or above it. **The kept iterate is returned with its
  bound, whatever the bound is**; the caller compares it with its own gate.
- At most eight corrections are made. If the eighth still at least halved the bound, the
  refinement had not finished, and the call refuses with "bound not reached". Eight and
  not six: a correction contracts by 2⁻¹⁰ or better while it is above rounding, so six
  exhaust the precision, and on small matrices whose smallest singular value is near τ
  five were seen; two are left for corrections that halve the bound by chance at the
  bottom, so that a finished refinement is not refused.
- A residual or a bound that is not finite refuses with "bound not reached".

**What was measured**, on the seven patterns, for a right-hand side whose minimum-norm
solution is known, at α of 1, τ^½ and τ and `α ε = τ² 2⁻ᵏ`:

| | k = 0 | k = 10 | k = 20 | k = 30, 40 |
| --- | --- | --- | --- | --- |
| Contraction a correction | 0.05 to 0.17 | 5e-5 to 2e-4 | 1e-7 to 2e-5; on the lattices the first solve is off by 3e-4 to 0.3 of the answer | wheels as at 20; lattices 3e-3 to 0.05, or no factor |
| Growth, lattices | 8e3 to 3e4 | 4e6 to 2e7 | 4e9 to 2e10 | 5e12, or a zero pivot |
| Growth, wheels | 2e3 to 3e4 | the same | the same | the same |
| Proved bound over ‖x‖ after one solve | 9e-5 to 2e-3 | 9e-8 to 4e-6 | 9e-11 to 3e-3 | |
| after one correction | 3e-9 to 1e-6 | 1e-12 to 3e-11 | 6e-14 to 2e-9 | |
| after two | 2e-10 to 1e-7 | 5e-14 to 1e-12, its floor | 3e-14 to 1e-12 | |

**α changes no iterate**: the three values gave the same errors to every digit printed,
and differ only in what growth and β read. So α is τ, the certificate's own top block.
**k is fixed at 10.** At 0 the contraction is that of `τ² / (σ_r² + τ²)` with `σ_r` about
three times τ, and a dozen corrections would be needed; at 20 the lattices' growth is
10¹⁰ and the first correction is poor; at 30 they give a zero pivot. At 10 every pattern
contracts by four decades a correction and the true error is about 10⁻¹⁵ of the step
after two. The second factorization is to be timed with the build. A matrix whose
factor is poor at 10
cannot give a wrong step, only a refusal: its bound does not fall, and the bound is what
is returned.

The diagnostic sum at k = 10, `η = 2β/τ`: 2e-4 on the loose wheel, 0.006 and 0.02 on the
loose lattices, 0.09 on the other wheels, 0.11 on the lattice of 15 and 0.48 on the
lattice of 21.

## An option: the kept rows alone

(Put by the reviewer; checked and measured here.) Let `x_S = J_S⁺ r_S`, the minimum-norm
solve on the kept rows with nothing from the rows left out. Since `B⁺ = J_S⁺ C⁺` and
`C⁺ C = I`, `x_S = B⁺ (C r_S)`, so

    x − x_S = B⁺ [ 0 ; r_T − W r_S ] ,        ‖x − x_S‖₂ ≤ ‖r_T − W r_S‖₂ / τ .

With `e = r_T − J_T x_S`, one product that needs no W, and `J_T = W J_S + E_T`,
`J_S x_S = r_S`:

    ‖r_T − W r_S‖₂ ≤ ‖e‖₂ + γ ‖x_S‖₂ ,   so   ‖x − x_S‖₂ ≤ ( ‖e‖₂ + γ ‖x_S‖₂ ) / τ .

The bound from x to x* adds to it. The identity and both inequalities were checked
numerically on three hundred small cases before being written here.

**It is as good as the residual on the rows left out is consistent.** Measured on the
seven patterns, with a residual in the range of J, and with noise of 10⁻⁹ and of 10⁻³ of
its size added to every row. Figures are over the norm of the step:

| Residual | The step's own bound to x* | Kept rows alone: bound to the step | Kept rows alone: measured |
| --- | --- | --- | --- |
| in the range | 5e-10 to 6e-8 | 1e-10 to 2e-8 | 1e-15 to 1e-11 |
| noise of 1e-9 | 5e-10 to 6e-8 | 4e-7 to 3e-5 | 1e-9 to 2e-7 |
| noise of 1e-3 | 6e-10 to 2e-7 | 0.4 to 26 | 1e-3 to 0.2 |

The step's own bound does not depend on how far r is from the range, because it solves
the least squares problem; the kept rows alone ignore the part of r on the rows left out,
and their bound says so. And with W sparse **the option saves little**: on these
patterns the full solve takes 1.4 to 3.5 milliseconds and the kept rows alone 0.3 to 1.5,
in a prototype whose products are not the library's.

It is offered as an option, with its bound returned, for a caller whose residual on the
rows left out is known to be consistent. **Whether a step that differs from the
truncated answer by a stated bound is acceptable is the caller's to decide**, and the
full step stays the default.

## What was measured

Seven patterns from the consumer: a wheel of 400 spokes and lattices of 15 by 15 and 21
by 21 rectangles, each well constrained, with redundant conditions, and with redundant
conditions and too few; 1,604 to 8,698 rows. A scratch driver over the library's
existing calls, on a machine otherwise idle unless said.

**The certificate**, with the search at 1e-6 of the scale's upper bound:

| Pattern | Rows by columns | Kept, left out | γ over the scale | τ over the scale |
| --- | --- | --- | --- | --- |
| wheel, well constrained | 1,604 by 1,602 | 1,602, 2 | 4.3e-14 | 1.2e-5 |
| wheel, redundant | 2,004 by 1,602 | 1,602, 402 | 3.8e-13 | 1.2e-5 |
| wheel, redundant and loose | 1,803 by 1,602 | 1,404, 399 | 3.8e-13 | 1.8e-4 |
| lattice 15, redundant | 4,414 by 3,600 | 3,600, 814 | 7.1e-14 | 1.3e-4 |
| lattice 15, redundant and loose | 4,384 by 3,600 | 3,572, 812 | 6.7e-14 | 4.9e-4 |
| lattice 21, redundant | 8,698 by 7,056 | 7,056, 1,642 | 9.9e-14 | 9.3e-5 |
| lattice 21, redundant and loose | 8,656 by 7,056 | 7,016, 1,640 | 9.6e-14 | 4.1e-4 |

The kept counts are the ranks the consumer's singular value decompositions give where it
has them, and the column counts where the columns are independent. The band from γ to τ
contains 2e-8 and 1e-6 on every pattern, by eight decades in all. Both proofs of (c)
gave the τ shown. Their own resolution, read from what they need to succeed, is about
5e-7 to 3e-6 of the scale for the squared proof and 2e-8 to 3e-7 for the other: neither
would certify kept rows that really sat near 2e-8, and none of these do. The default
scale interval was 1.13 to 1.94 wide and refused nothing.

**The parts of a call:**

| | Wheels, d about 400 | Lattice 15, d about 812 | Lattice 21, d about 1,640 |
| --- | --- | --- | --- |
| The search | 0.19 to 0.24 s | under 0.01 s | under 0.01 s |
| Every null vector found afresh | 0.63 to 0.73 s | 0.44 to 0.54 s | 1.9 to 2.5 s |
| Every bound proved again | 0.004 to 0.009 s | 0.023 s | 0.09 s |
| Multipliers computed, and real | 15,000 to 75,000; 3 to 4 a row | 1.4 million; 4.0 to 4.4 a row | 5.7 million; 4.0 to 4.2 a row |
| `I + W Wᵀ`, sparse | 0.02 s | 0.0001 s | 0.0002 s |
| The same, dense | 0.02 s | 0.08 s | 1.3 s, beside another job |
| τ, one signed factorization | under 0.001 s | under 0.001 s | under 0.001 s |

The search is slow on the wheels because their `J Jᵀ` is dense. "Every bound proved
again" used the library's product bound once a vector, over the whole matrix; a product
over only the rows a vector touches would cost in proportion to its four entries and is
not built. The squared proof of (c) took 0.18 s a factorization on the wheels. The two
rows left out of the well-constrained wheel are the exception to "three to four": each
depends on some two hundred rows.

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

In the second review, after the measurement:

- **The search's threshold and the certified floor are two arguments.**
- **The step is built on the proof on the matrix itself**, with t = τ; with t the
  largest entry it certified nothing on the lattices. The squared proof is the stated
  alternative. A new analysis when the kept set changes is accepted: the kept set is a
  matter of structure and changes only when the rank does.
- **Multipliers below a drop tolerance are set to zero** and the bounds proved again.
- **The vectors are kept from call to call** and their bounds proved again each time.
- **The system of order d is factored sparsely.**
- **The kept rows alone are an option**, with the bound stated; whether a caller may use
  it is the caller's decision.

In the third, before the build:

- **The drop tolerance is an argument of its own**, relative to each vector's largest
  multiplier and not to the matrix, with 10⁻¹² the value the measurement suggests. It is
  not tied to the thresholds of the supports: a support's class is read from a multiplier
  as it is used, after the drop.
- **Kept vectors that are out of date are a refusal by name.** The library does not
  search again by itself; the caller decides when a search is worth its cost.
- **The bound of one vector reads only the rows that vector touches**, so that proving
  every bound again costs in proportion to the vectors and not to the matrix.

In the fourth, on the solve, after its measurement:

- **The solve has its own shifts**, `α = τ` and `ε = τ · 2⁻¹⁰`, on the certificate's
  analysis. k is fixed at 10 and not searched for at each call.
- **The contraction is stated by blocks**, because G is not normal.
- **The bound returned is a posteriori**, from proved residuals and the τ proved at
  this call. A rate predicted from the factor's backward bound is reported and refuses
  nothing.
- **A refinement that stalls returns its bound**; one still halving at its eighth
  correction, or a residual that is not finite, is "bound not reached".
- **The kept rows alone are not built.** The consumer will not call the option.

**Built so far:** the search with kept vectors, the reproof, the certificate, the bounds
on the largest singular value, and the step with its two bounds, as
[kept rows, the certificate and the step](rank-deficient.md). The refresh of one vector,
the option of the kept rows alone, the supports and the bases are not built.
