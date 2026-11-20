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
