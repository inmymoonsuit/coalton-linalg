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
