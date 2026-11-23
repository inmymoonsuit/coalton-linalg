# Provenance and licensing status

The initial code was independently developed in an earlier sparse-factorization
prototype for this project. This repository copies that work; the original
prototype, supplied CAD fixtures, held-out fixtures, reference results, historical
logs, and benchmark dumps were left in place.

At the initial migration, the algorithm source was unchanged apart from package
and safety-feature names. Batches 1 and 2 subsequently moved the shared finite
predicate, added the checked CSC entry, and implemented new numeric/CSC primitives.
A migration manifest records source hashes and destination paths. Test harness and
driver changes make the work loadable through ASDF and remove machine paths and
application-specific exemptions.

## Mathematical origin

- Ordering and symbolic structure were derived from elimination-graph adjacency:
  choose the smallest live degree, break ties by original index, connect surviving
  neighbors, and eliminate the vertex.
- Numeric LDLᵀ was derived by equating columns of the matrix product.
- Solve was derived from forward substitution, diagonal division, and backward
  substitution.
- A separate reviewing agent derived the Python verifier from matrix products,
  residuals, exact rational arithmetic, and set-based graph elimination.

No SuiteSparse, Eigen, faer, SuperLU, MUMPS, or SciPy implementation, or textbook/
paper algorithm listing, was used as an implementation source. API documentation
is reference material only. Coalton's language/array/cell/ABI documentation was
consulted for host-language contracts. Build-system syntax follows the official
documentation listed in [design.md](design.md).

The small checked-in fixtures are generated mathematical cases from our independent
verifier. They contain no CAD files or vendor-derived geometry.

## License

The maintainer selected [MIT](../LICENSE) on 2026-11-19 for this independent
library. The copyright notice uses the maintainer's public GitHub identity,
inmymoonsuit, and coalton-linalg contributors. This choice does not derive from
the dependency or any consuming application.

No third-party numerical implementation is bundled. Dependency licenses remain
with their respective projects. The standard MIT text follows
[the published license](https://choosealicense.com/licenses/mit/).
Record explicit attribution and notices for any future permitted source adoption.
## Batches 1 and 2

Stable norms were derived from the scaled sum-of-squares identity. CSC construction
uses stable coordinate sorting and ordered duplicate folding; transpose uses counts,
prefix sums, and scatter; permutations use explicit coordinate substitution.
Products follow their defining sums. These new algorithms are written in Coalton.
Independent reviewers supplied dense/exact tests and challenged the calling contract.
The SBCL runtime source was consulted for measurement semantics, not numerical algorithms.

## Pre-publication revision

Standard Coalton finiteness and absolute-value operations replace local duplicates;
their implementations are dependencies, not copied source. Typed exceptions use
Coalton's documented language facility.

The exact-singular regression is independently derived from the integer Gram
identity for rows [1,1,2], [1,0,1], [1,0,1]. It stores the resulting exact integer
entries directly and proves rank using rational arithmetic on the F64 inputs.
The external review's floating Gram candidates were inspected but not imported:
all three were full rank as serialized. No external generator or factorization
code was copied into the library or verifier.

## Batch 3a: compact symbolic analysis

The graph-elimination definition is unchanged. The implementation uses Coalton's
mutable hash tables for adjacency and an independently written indexed binary
heap keyed by degree and original index. Standard-library container contracts
were read; no numerical-library implementation or algorithm listing was used.
The entry-budget identity is n plus original edges plus distinct fill edges.
Independent set elimination, analytic graph families and planted artifact defects
check the resulting structure. This remains an independently derived implementation.

## Batch 3b: caller-buffer products

Both output-buffer product directions reuse the original independently derived
CSC sum kernel. Identity checks use Coalton's existing pointer-equality primitive;
no new array representation or numerical dependency is introduced. Independent
rational dot products and adversarial calling-contract tests verify the new entry
points without relying only on the allocating versions as a reference.

## The clean-room question, decided

Several sections below, written from Batch 3d onward, say that whether following a published
method closely from its description meets this repository's clean-room rule was an open
decision of the maintainer's. **He decided it on 2026-11-22: an implementation written
clean-room from a published method's description is acceptable, and textbook references
and pseudocode may be read as references.** What is unchanged: code is never copied, in
any form; no numerical-library implementation and no listing in a real programming
language is consulted; and each batch names what was read. Every batch before this
ruling was written without opening a pseudocode box, as its section says.

## Batch 20: least squares on an operator

The bidiagonalization is that of Golub and Kahan (1965), and its use for least squares,
with damping carried by a second reflection a step, is the method Paige and Saunders
published in 1982. Both are known here from the literature. The recurrences, the two
identities that give the gradient and the residual from the scalars, and the properties
the contract states were derived in [the design](design-iterative-least-squares.md)
before any code, from the two relations that define the sequences, and the code was
written from that page. No paper was opened while it was written, and no implementation
and no listing at any time. Forming both norms from the operators before a tolerance is
called met, the rule for the estimate after a confirmation that failed, the divisions by
powers of two, the lengths formed without a squared norm, and the statuses were chosen
here. That the method is more reliable than conjugate gradients on the normal equations
is the published comparison's finding, cited as such and not measured here. The
reference values in the tests are exact minimizers, residuals and gradients in rational
arithmetic, and the iteration is recomputed independently and compared bit for bit.

## Batch 19: dependent sets, and the right null space

No source was opened. That the vectors of the rows left out are the fundamental circuits
of those rows, and that merging rows whose circuits overlap gives the components of the
rows as a matroid whatever independent set was kept, is textbook matroid theory, stated
in the design and checked there in exact arithmetic on four hundred matrices before any
code. Joining sets through a forest of parents with the paths shortened is the classical
union-find structure, written from its description. The bound on a stored multiplier's
distance from the exact one was worked out here, in its second form while building: the
first, in the design, carried the square of one over the bound proved on the kept rows
and lost an entry on a consumer's pattern. The cut relative to each vector's largest
entry is the consumer's own rule, taken as a statement of what to compute. The reference
values in the tests are exact multipliers and exact components in rational arithmetic.

## Batch 18d: a refresh on a vector's own rows, and a factor kept

No source was opened. Conjugate gradients with a preconditioner is the classical
iteration, written here from the plain one this library already has in Batch 13, whose
section names its sources, with the preconditioner's solve put where the residual is
read; using a factorization of a nearby matrix as the preconditioner is as old as the
method. The least squares fit on a vector's own rows is the normal equations of a
handful of rows with a Cholesky factor in place and one step of refinement. That
neither needs a proof, because the product bound proves the vector and the proved
residual the system's solve, is this design's own. The reference values in the tests
are exact, in rational arithmetic.

## Batch 18c: kept vectors refreshed

No source was opened. The multipliers of a row against the kept rows are a least squares
problem, solved through the same augmented matrix and refinement as the step of Batch
18b, whose sources that section names. That nothing need be proved of the solve, because
the product bound of Batch 11 proves the vector it produces, is this design's own; so is
refreshing at a share of the floor and not at the floor, which the consumer asked for.
The reference values in the tests are exact products and exact minimum-norm least
squares solutions in rational arithmetic.

## Batch 18b: the step on a rank-deficient matrix

No source was opened for this batch; the mathematics is in
[the design](design-rank-deficient-step.md), written and reviewed before the code. Three
things in it are known from the literature as statements. That the minimum-norm solution
of a full-row-rank system is the first block of the solution of the augmented system,
and that a system factored with a small regularizing block can be refined against the
system without it, are classical in least squares (Björck's work on the augmented
system); the form here, with both diagonal blocks tied to the bound a certificate
proved, was chosen from a measurement. The decomposition of the difference of two
pseudo-inverses is Wedin's (1973); it is an identity, and it and the bound taken from
it were checked by multiplying out numerically before they were used. The inversion
lemma is textbook. **What the step's bound rests on was worked out here and checked by
the consumer that asked for the step**: that for any pair the error of the kept rows'
solve is an exact expression in the two residuals, in two orthogonal parts, so that no
rate of convergence is needed; and that the iteration matrix of the refinement is not
normal, so that its contraction is stated by blocks, which that consumer pointed out.
No listing or pseudocode was read and no implementation. The reference values in the
tests are exact minimum-norm least squares solutions in rational arithmetic, by the
tests' own elimination.

## Batch 18a: kept rows and the certificate

Nothing here is a method taken from a source, and no source was opened for it: it is
this library's own earlier pieces put together under
[a design](design-rank-deficient-step.md) written and reviewed before the code. The
search and the vector of a row left out are Batch 11's. The bound on one vector's product
is Batch 11's, the standard rounding analysis of an inner product as that batch records
it, applied to only the rows a vector touches. The certificate on the kept rows is the
one from the inertia that Batch 7a's section records, with both diagonal blocks set to the
bound asked about. The bounds on the largest singular value are two textbook facts,
that the square of the 2-norm is at most the product of the 1-norm and the
infinity-norm and that `‖J v‖ / ‖v‖` is a lower bound for any v, with the power method
choosing v; the rounding allowance was worked out here. No listing or pseudocode was
read and no implementation. The structure follows from the contract and from a
measurement on a consumer's patterns: the vectors are kept sparse because most of a
computed vector is rounding, the bound is proved for the vector as kept and not as
computed, the proof is on the augmented matrix because the kept rows' Gram matrix can
be dense, and the ladder of bounds tried is there so that a steady caller pays for one
factorization. The reference values in the tests are exact rational products and exact
inertias from the stored doubles, by the tests' own elimination.

## Batch 17: damped least squares

The damped step is the normal equations with a multiple of a positive diagonal added, the
method of Levenberg (1944) and, with the diagonal taken from the columns of the matrix,
of Marquardt (1963). It was written from that statement of it, which is one equation, on
top of this library's own Gram product and factorization. No listing or pseudocode was
read and no implementation. The structure differs from a textbook's because of this
library's contract: the matrix and the residual are divided by powers of two before
anything is squared, the damping is held in both the caller's units and the divided
ones, and every refusal is a named status that leaves the output unwritten. The
reference values in the tests are steps in exact rational arithmetic by a dense
elimination written for the tests.

## Batch 16: the range of conjugate gradients and of dependent rows

Both changes are one idea applied twice: divide the data by a power of two before a method
squares it, and multiply back what has units. Its sources are prose in Higham, *Accuracy
and Stability of Numerical Algorithms*, second edition: section 27.8, on intermediate
quantities that overflow or underflow in the computation of representable ones, which
notes that for about half of all machine numbers the square leaves the range and that
scaling by a power of the machine base adds no rounding error; and section 9.8, which
states Bauer's result that when rows and columns are scaled by powers of the base the
rounding errors of an elimination scale exactly as the matrix does. That conjugate
gradients do not depend on the length of a search direction follows from the iteration as
the contract page gives it: the step length is the ratio that the direction's length
cancels from. No listing or pseudocode was read for either, and no implementation. The
structure is this library's own: one power of two for a whole call, chosen from the first
residual or from the matrix's largest entry, held in the workspace, with the division
done in the pass that already copies the data.

## Batch 14: an outside audit

The changes of this batch were written by Robert Smith and contributed to this
repository as pull request 1; the commits carry his name. They are his own work on this
library's sources, contributed under the repository's licence. What was added here is
the placement of three helpers in the file that uses them, a guard in the build script,
the planted controls, and the record. No numerical-library implementation and no
algorithm listing was opened here in integrating it.

## Batch 13: conjugate gradients

The method of conjugate gradients is that of Hestenes and Stiefel (1952), and its
preconditioned form is the same method applied to a transformed system; both are
textbook material and are written here from the recurrences as they are derived on
[the contract page](conjugate-gradient.md#the-method): the step that minimizes the
error in the operator's norm along a direction, and the one correction that keeps the
directions conjugate. Replacing the running residual by one formed from the operator is
a known remedy for its drift; where it is done here, at the moment the running residual
claims convergence and at the limit, and the restart that follows a failed confirmation,
were chosen here. The checks for a direction without positive curvature and for a
preconditioned residual without a positive product follow from what those quantities
are. No numerical-library implementation and no algorithm listing was opened. The tests
recompute the iteration independently and measure residuals and errors in exact
rational arithmetic.

## Batch 12: sparse LU with partial pivoting

Gaussian elimination with partial pivoting, and its rounding analysis, are classical. The
organization used here, a column at a time with the set of earlier columns that a column
needs found by a depth-first search over the structure of L and applied in the order the
search finishes them, is the method of Gilbert and Peierls (1988), known here from its
description: that the solution of a sparse triangular system has its nonzeros at the
nodes reachable from the right-hand side's, and that a reverse finishing order is a
valid order of computation. The search without recursion, the storage, the pivot rule's
tie-breaking, the treatment of a refused call and the one-number backward bound with its
underflow term were worked out here. Preferring the diagonal within a threshold is common
practice in sparse elimination. Ordering the columns by the pattern of `Aᵀ A` is likewise
common practice and uses this library's own ordering. No numerical-library
implementation and no algorithm listing was opened. The tests' reference is an exact
rational product of the exported factors.

## Batch 11: dependent rows and skipped pivots

Leaving out an index whose pivot is negligible, in a factorization of a positive
semidefinite matrix, is long-standing practice and follows from what the pivot is: the
Schur complement of the index on those kept before it. The null vector of an index left
out is one triangular solve with the row of the factor that index leaves behind; that
and its error bound are derived on [the contract page](skipping.md#the-null-vector).
Refining the multipliers by a residual formed with the matrix itself and a correction
from the factor of the normal equations is the method of corrected seminormal equations,
known here from its description in the literature on least squares (Björck, 1987); it is
expressed through this library's own `refine!`. The bound on a transposed product is the
standard rounding analysis of an inner product, with the underflow and summation terms
worked out here. That a matrix within a Frobenius distance of one of rank r has its
singular values beyond r below that distance, and that adding rows does not lower
singular values, are textbook facts. Running the rank certificate on the whole matrix's
plan with zeros for the rows left out, and the count it then requires, were worked out
here. No numerical-library implementation and no algorithm listing was opened. The
second method of [the design](design-dependent-rows.md), by rotations, is not built; the
comparison with the literature that the design asks for remains to be made before it is.
The tests' references are exact rational eliminations, inversions and inertias written
for the tests.

## Batch 10: iterative refinement and the Cholesky-form solve

Iterative refinement is classical, from Wilkinson's work on rounding errors in the 1960s,
and its convergence as a stationary iteration with matrix `I − B⁻¹ A` is a line of
algebra, written on [the contract page](refinement.md#derivation). What is specific
here was decided with the consumer that uses it: the library forms the residual, the
norm is the infinity norm against an absolute tolerance, the stall factor is the
caller's, and the vector returned is the best iterate and not the last. The two
eigenvalue computations on the contract page, for a shifted factor and for the
regularized augmented matrix, were worked out here. The triangular solves with the
Cholesky-form factor follow from `G Gᵀ`. No numerical-library implementation and no
algorithm listing was opened. The tests recompute the iteration independently and
compare with exact rational solutions.

## Batch 9a: scaling at the edges of the range

The defects were found by an independent held-back review and by reading the source
against its cases. Multiplying the two factors first, and applying them in turn only
where their own product leaves the range, was worked out here, after the reviewer showed
by example that no fixed order of the two multiplications is exact; the reviewer's own
suggestion, adding the factors' exponents and scaling once, is the same idea, reached
here without decoding a number. Lifting the column pass by the largest power of two was
worked out here. No numerical-library implementation and no algorithm listing was
opened.

## Batch 9: scaling by powers of two

Scaling rows and then columns by their largest magnitudes is classical equilibration,
and rounding the factors to powers of two so that the scaling is exact is common practice;
both were known from descriptions. That one pass each way suffices for the property the
contract states is argued on [the contract page](scaling.md#derivation), and the way the
power of two is found, by ten exact multiplications or divisions and no decoding of the
number, was worked out here after the first version, which decoded it, was measured to
allocate. Scaling a positive definite matrix by its diagonal is Jacobi scaling, and the
statement that it is near the best diagonal scaling is van der Sluis's (1969), cited from
its statement and not used by the code. No numerical-library implementation and no
algorithm listing was opened. The tests compare every scaled value with the exact
rational product.

## Batch 8: faster minimum-degree analysis

**The representation is published; the code was written from its description.** Keeping
eliminated vertices as sets of neighbours that absorb one another, in place of the
elimination graph, is the quotient-graph idea of George and Liu, set out in their
accounts of minimum-degree ordering (1980 and 1989). It was known from descriptions. The
argument that it reaches the same graph, and that a vertex's two lists fit in the space
of its initial degree, is written out on [the contract page](spd-ldlt.md#derivation).
The degree here stays exact; the approximate degree of Amestoy, Davis and Duff (1996),
which is what makes the well-known ordering codes fast, was deliberately not used, because
it changes the order.

The row of flags for a vertex of large initial degree was worked out here from a
measurement: the first version made the hub pattern slower, and [the batch
record](batch8.md) has the figures. Sorting every column by one pass through the rows is
the ordinary transpose-twice argument. No ordering code, listing or pseudocode box was
opened while writing any of it. The tests' reference is a set-based elimination written
for the tests, run against the earlier code first; the complete check's independent
verifier is the Python program of Batch 3a.

Whether following a published representation this closely from its description meets
this repository's clean-room rule is the same open maintainer decision as for Batch 3d.

The Gram analysis's grouping by digits is the classical least-significant-digit
distribution sort, written from its definition, and its fill by row follows from the
order in which a group lists its owners; neither was taken from a source.

## Batch 7a: backward bound and rank certificate

The norm bound is the entrywise bound of Batch 7 summed, with one observation made here:
`|L| |D| |Lᵀ|` is positive semidefinite, so its Frobenius norm is at most its trace, which
the kernel was already accumulating. The underflow term is the standard model of rounding
with a gradual-underflow term, counted for this kernel on [the contract
page](quasidefinite.md#derivation).

**The certificate's argument is a consumer's.** That the inertia of
`[[t I, Aᵀ], [A, s I]]` counts the singular values of A above √(s t), and that the computed
pivots prove it once two diagonal margins cover the backward error, was proposed by the
author of an application that uses this library, and checked here before it was written
down; the conditions about underflow, the rounded margins and the rounding of the growth
figure were added in that check. Sylvester's law of inertia and the monotonicity of
eigenvalues under a positive semidefinite change are classical. The tests' truth is an
exact rational inertia of `A Aᵀ − s t I`, written for the tests.

## Batch 7: signed pivots and the augmented matrix

**The setting is published; the ordering rule and the growth figure were worked out here.**
That a matrix `[[H, Aᵀ], [A, −G]]` with H and G positive definite has an L D Lᵀ
factorization with diagonal D under every symmetric permutation is Vanderbei's result on
symmetric quasi-definite matrices (1995), and that such a factorization can be unstable,
with accuracy governed by the size of the factors, is the analysis of Gill, Saunders and
Shinnerl (1996). Both were known from their statements; the short proofs on [the contract
page](quasidefinite.md) were written out here. Regularizing a minimum-norm problem as an
augmented system with two small diagonal blocks is common practice in sparse least squares
and optimization.

The signed entry is the library's own L D Lᵀ kernel with one test changed. The growth
figure follows from an observation made here, that `|L| |D| |Lᵀ|` is a Gram matrix and so
has its largest entry on its diagonal, which makes the exact figure cost one multiply-add
an entry. The rounding bounds are the standard componentwise analysis of this arithmetic,
with constants counted for this kernel and checked exactly by the tests.

The held ordering came from experiment. Plain minimum degree was measured to give growth
of about two over the trailing diagonal; holding rows until a column of theirs had gone
was tried, failed on two rows sharing their only column, and became the rule that a column
releases one row. Orderings that place a constraint after a variable it touches are known
in the literature on saddle-point matrices, from descriptions only; no ordering code,
listing or pseudocode box was opened. The experiments used NumPy on public fixtures, with
orderings written for the experiment, and are summarized in [the batch record](batch7.md).
The tests' references are exact rational products of the exported factor, exact inertia
by a separate elimination, and exact residuals, all written for the tests. NumPy was used
after the implementation as a singular-value reference for the measured solutions.

Whether following published results this closely from their statements meets this
repository's clean-room rule is the same open maintainer decision as for Batch 3d.

## Batch 6a: Cholesky-form entry

The Cholesky factorization is classical and its column form follows from equating the
entries of G Gᵀ with the matrix. The entry was written from that identity, over the
elimination order and structure the library already had. Its purpose is to make a
published theorem applicable, Rump's verification of positive definiteness (BIT 46, 2006),
which a consumer cites; the theorem was not used to write the code, and the componentwise
bound the tests check is the standard rounding analysis of this arithmetic, stated on the
contract page with its constants. No numerical-library implementation and no algorithm
listing or pseudocode box was opened while writing it. The tests' reference is the exact
rational product of the exported factor.

## Batch 6: Jacobi eigenvalues and singular values

**These are classical published methods.** The rotation that annihilates one off-diagonal
entry and its use in cyclic sweeps is Jacobi's method (1846); the choice of the smaller
root and the update of the two diagonal entries by t b follow Rutishauser's formulation
(1966); rotating pairs of columns to orthogonality is the one-sided method of Hestenes
(1958); stopping on a test relative to the two diagonal entries or the two column norms is
the criterion analysed by Demmel and Veselic (1992), whom the planned-batch page cites. The
code was written from the derivation on [the contract page](jacobi.md), which was written
first: the rotation equation, its smaller root, the invariant kept by accumulating the
rotations, and the reason the one-sided form needs no Gram matrix are worked out there. The
methods were known from their descriptions. No numerical-library implementation and no
algorithm listing or pseudocode box was opened while writing it; the same open maintainer
decision as for Batches 3d, 4 and 5 applies.

The tests' reference is this library's own choice: every computed value is bracketed by an
exact count of eigenvalues below two rational shifts, by symmetric elimination in rational
arithmetic and Sylvester's law of inertia, so no eigenvalue routine is trusted as a
reference. NumPy and SciPy were used after the implementation, as a comparison and a timing
reference only.

## Batch 5: orthogonal decompositions

**These are standard published methods.** Householder reflections, QR with column
pivoting by largest remaining norm (Businger and Golub, 1965), and the reduction of the
kept upper trapezoid to a triangle by reflections from the right, which gives the complete
orthogonal decomposition and the minimum-norm solution, are textbook material; the
planned-batch page asked for them by name and links prose descriptions. The code was
written from the derivation on [the contract page](orthogonal.md), which was written
first: the reflector's formulas are verified there from its definition, and the solve and
the nullspace follow from the factorization's block form. The methods were known from
their descriptions. No numerical-library implementation and no algorithm listing or
pseudocode box was opened while writing it; the same open maintainer decision as for
Batches 3d and 4 applies.

Two choices are this library's own. Column norms are measured afresh from the current
entries in the pass that updates each column, where the classical method updates old norms
by a formula and guards the update. Sums of squares are scaled by fixed powers of two
chosen by comparison. The Kahan matrix and its small diagonal perturbation are a published
test case, used from its description. The tests' references are an exact rational product
of Q read back column by column, and an exact minimum-norm solution from an exact rank
factorization, both written for the tests. NumPy and SciPy were used after the
implementation, as a comparison and a timing reference only.

## Batch 4: dense symmetric indefinite factorization

**This factorization follows published methods.** Block-diagonal pivoting with 1-by-1 and
2-by-2 pivots and the threshold (1 + sqrt 17)/8 is the method of Bunch and Parlett (1971)
and Bunch and Kaufman (1977); the column search that follows growing off-diagonal maxima
until it meets a mutual maximum is the bounded variant of Ashcraft, Grimes and Lewis
(1998), often called rook pivoting. The planned-batch page asked for this family by name.
The code was written from the derivation on [the contract page](indefinite-ldlt.md), which
was written first: the multiplier bounds, the threshold equation, the termination argument
and the lower-triangle exchange are worked out there from the elimination identities. The
methods were known from their descriptions. No numerical-library implementation and no
algorithm listing or pseudocode box was opened while writing it. Whether following
published methods this closely from their descriptions meets this repository's clean-room
rule is the same open maintainer decision as for Batch 3d.

Inverting the 2-by-2 block in a form divided by its off-diagonal entry was worked out here
from the bound on the scaled determinant. The tests' references are an exact rational
product of the exported factor, an exact symmetric elimination with a different pivot rule
for the inertia, and an exact rational inverse, all written for the tests. NumPy and SciPy
were used after the implementation, as an eigenvalue-sign comparison and a timing
reference only.

## Batch 3d: condition estimate

The estimate was derived from convexity: |inverse(M) x|_1 is convex in x, its maximum over
the unit 1-norm ball is the norm sought and is attained at a vertex, and one further solve
with the sign vector gives a subgradient that says whether any vertex is promised an
increase. [The contract](condition.md) carries the derivation.

**This routine follows a published algorithm.** It is Hager's 1-norm estimator (1984) with
Higham's refinements (1988), the method widely used for factor-based condition estimates.
Every step has the published shape: the uniform start, the move to the unit vector the
subgradient favours, the stop when no vertex promises an increase, the stop on a repeated
vertex, the cap of five rounds, and the final alternating-sign vector with magnitudes from 1
to 2 and its bound 2|y|_1/(3n). It was written from recall of the method's description and
from the derivation on the contract page. No numerical-library implementation and no
algorithm listing or pseudocode box was opened while writing it. Whether following a
published method this closely from its description meets this repository's clean-room rule
is the maintainer's decision and is open; the alternative is an estimator of a different
construction under the same tests. The reference values in the tests are exact rational
inverses computed by a separate dense elimination written for the tests.

## Batch 3c: fixed-pattern Gram products

The product pattern and refill follow the defining sums G[i,j] = sum_k B[k,i] B[k,j]
with B = A or its transpose. Entries are grouped by inner index with a stable sort;
the lower-triangle scan starts at an entry's own sorted position, which follows
from canonical CSC ordering alone. No numerical-library implementation or algorithm
listing was consulted. Independent exact rational products check every value. On
thirteen application systems the stored values also agree bit for bit with an
independently produced reference product wherever that reference stores an entry;
the reference was used as a comparison only, after the implementation was written.
