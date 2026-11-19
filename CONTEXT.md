# Numerical linear algebra

Shared mathematical vocabulary for coalton-linalg. Function signatures and
calling conventions live in the [interface reference](docs/api.md); this glossary
defines concepts, not additional behavior.

## Language

**Sparse matrix**:
A matrix represented using selected stored entries, with omitted entries
understood to be zero.

**Sparsity pattern**:
The set of stored matrix coordinates, independently of their numerical values.
A stored zero belongs to the pattern.

**Stored entry**:
An explicitly represented coordinate and its value, including a possible zero.

**nnz**:
An entry count under a stated storage convention: CSC counts stored entries,
including zeros. Factor counts may also include implicitly represented entries,
such as a unit diagonal, and must specify that convention.

**CSC (compressed sparse column)**:
A column-oriented sparse representation using column offsets, row indices,
and values.

**COO (coordinate format)**:
A sparse representation listing row, column, and value for each supplied entry.
Repeated coordinates are duplicates, not additional matrix positions.

**Canonical CSC**:
CSC with strictly increasing row indices within each column and no duplicate
coordinates; canonical does not mean that stored zeros have been removed.

**Symmetric positive definite (SPD)**:
A real symmetric matrix A for which xᵀAx > 0 for every nonzero real vector x.
A floating-point factorization's acceptance is not an exact SPD certificate.

**LDLᵀ factorization**:
A decomposition into a unit lower-triangular matrix L, diagonal matrix D, and Lᵀ,
possibly after a symmetric permutation. For an SPD matrix in exact arithmetic,
the diagonal entries of D are positive.

**Symbolic analysis**:
Preparation of an elimination ordering and factor structure from the sparsity
pattern. It is distinct from computing numerical factor values.

**Numeric factorization**:
Computation of factor values for a particular matrix. Refactorization repeats
this computation for changed values with the same analyzed pattern.

**Solve**:
Computation of a solution to a linear system for a supplied right-hand side.
It is distinct from a nonlinear solver's iteration and convergence decisions.

**Ordering / permutation**:
A bijective reindexing of variables; an elimination ordering specifies the order
in which variables are eliminated. This differs from numerical pivot selection.

**Pivot**:
The current diagonal quantity used in an elimination step. A pivot test in
floating-point arithmetic does not by itself determine exact rank.

**Fill**:
Factor structure at positions absent from the input pattern; symbolic fill can
remain present even when cancellation makes its computed value zero.

**Shift**:
A scalar multiple of the identity added to or subtracted from a matrix.
The factorization contract specifies the sign; applications choose the scalar.

**Residual**:
The discrepancy b - A*x̂ for a computed solution x̂. A small residual alone
does not establish closeness to the exact solution.

**Forward error**:
The difference between the computed solution and a reference exact solution,
usually measured with a stated norm and normalization.

**Backward error**:
A measure of the perturbation to the input problem needed to make the computed
answer exact. Numerical checks must state the perturbation model and normalization.

**Conditioning**:
Sensitivity of the mathematical answer to changes in the problem's data.
It is a property of the problem, distinct from an algorithm's numerical stability.

**Refusal**:
An operation's reported inability to produce a result under its contract.
A numerical refusal is not automatically a geometric diagnosis or rank decision.