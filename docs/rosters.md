# The rosters

Two questions are asked of every call in the library at once, by test, rather
than of each routine when someone remembers to. Both came from an outside audit
that found defects the routines' own tests passed ([record](batch14.md)); the
[contributor notes](../CONTRIBUTING.md#four-questions-for-every-numeric-routine)
have the four questions that audit amounts to. This page is the contract of the
two that are mechanized, and how a new call joins them.

## Membership is read from the sources

Neither roster lists its members by hand. Each reads every `src/**/*.ct` file,
takes the names its package exports, and keeps those whose declared type fits a
rule:

| Roster | A member is an exported function whose type |
| --- | --- |
| Aliasing | has two or more arrays among its own parameters |
| Scale sweep | mentions a floating-point number anywhere, or returns a status |

The roster then fails for a member it has no entry for, and for an entry that
names no member. So a call added to the library joins by being declared, and
the suite is red until someone has answered for it. Arrays passed inside a
handle, and functions taken as arguments, are not counted as array parameters.

**What the rules do not see**, found by an independent review and left as it is:

- **A call that takes its numbers inside a stored matrix or a handle and returns no
  status is not a member of the scale sweep.** `set-matrix!` of
  [damped least squares](damped.md), which divides the matrix and squares it, and
  `transpose` and `permute` of a stored matrix, are three. The damped entries of the
  sweep call the first on the way to the calls they answer for; nothing makes them.
- **An array parameter is recognized by how its type is written**: the array package's
  usual nickname, not parenthesized. Every source writes it so today.
- **A pair is two positions of one element type and one length in the entry's own
  arrays.** A pair that exists only for other shapes is in no declaration, and the roster
  does not ask that each be reached by some fixture. Three positions as one array are
  never tried.
- **The one exemption, `minimum-degree`, rests on its caller**: the package is private and
  its one caller passes the two arrays of a validated pattern, which cannot be one
  array. That is a statement about today's caller, and nothing fails when another is
  added.

## The aliasing roster

[`tests/aliasing.lisp`](../tests/aliasing.lisp). The contract every member is
held to:

> A call given one array in two argument positions either does exactly what it
> does when given two arrays of equal contents, or refuses before it writes
> anything.

For each member and each pair of its array arguments that could be one array —
the same element type and the same length — the roster makes the call twice:
once with two arrays of equal contents, once with one array in both positions.

- If the second call refuses, every argument must be bit for bit what it was. **A
  refusal is the library's typed error, or a status that is not success.** Any other
  condition, an index out of range at a checked safety level for one, is a defect: it is
  not caught, and it ends the roster's run.
- If it succeeds, the first must have succeeded too, every argument and every
  returned array must agree with the first call's bit for bit, and the shared
  array must hold what the first call wrote to whichever of the two it wrote.
  A call that writes both cannot succeed.

Each entry declares which of the two each pair does, and the roster fails when
a pair does the other. So accepting an alias that was refused, or refusing one
that was accepted, is a change to the interface that the entry has to be
changed to make.

The matrices behind the solves are chosen so that every factorization reorders
or pivots, and each asserts that it does. With a matrix factored in its own
order, a solve that read its right-hand side after beginning to write could
read the right entries by accident.

### What the entries declare

| Call | Pairs that may be one array | Pairs that are refused |
| --- | --- | --- |
| `cg!`, `pcg!`, `refine!`, `lsqr!` | `b` and `x` | |
| `solve!` of the sparse, dense symmetric, orthogonal and LU factorizations; `solve-cholesky!`, `solve-skipping!`, `solve-transpose!` | the right-hand side and the output | |
| `matvec!`, `transpose-matvec!` | | the vector and the output |
| `from-csc`, `from-coo`, `permute` | the two index arrays, which are only read | |
| `gram:refill!` | | the values and the output |
| `row-maxima!`, `column-maxima!` | | the values and the output |
| `equilibrate!` | | any two of the values and the two scale vectors |
| `scale!` | the values and either scale vector; the two scale vectors | the output and anything else |
| `scale-vector!` | the vector and the output; the scales and the vector | the scales and the output |
| `equilibrate-symmetric!` | | the values and the scales |
| `scale-symmetric!` | the values and the scales | the output and either |

`augmented:refill!` has no pair: its output is longer than its input by the
order of the augmented matrix. `ordering:minimum-degree` is exempt, with its
reason in the roster: the package is the library's own, and its callers pass
it the two arrays of a pattern they have validated.

### Adding a call

Add an entry to `alias-roster`: the member's name, a function that makes its
array arguments afresh, a function that makes the call on them and returns 0
for success, and the outcome of each pair. Run the suite; where a pair's
outcome was guessed wrong the failure prints what it measured. Decide whether
that is the contract you want before writing it down.

## The scale-sweep roster

[`tests/scale-sweep.lisp`](../tests/scale-sweep.lisp). Multiplying every input
of a problem by a power of two changes no significand, so every output has one
right value at every scale: its value at scale 1 times a known power of the
factor. A factor's `L` does not move, its `D` moves with the matrix, a solution
of a system whose two sides move together does not move, a Gram matrix moves
with the square. Each of thirty-seven entries runs one small problem at twenty-five scales from
2^-1020 to 2^1020 and states that power for each output. Every run is put in
one class:

| Class | The run |
| --- | --- |
| `:exact` | succeeded, and every output is its scaled value bit for bit |
| `:close` | succeeded, and every output is within the entry's tolerance of it, against the largest magnitude in its array |
| `:looser` | succeeded, and an upper bound came out above its scaled value, which a bound may |
| `:lost` | succeeded, and an output's right value is too small to be a normal number: what came back is zero or below the normal range. Only at the bottom of the range: where the right value is too large to hold, a finite number is `:wrong` |
| `:refused` | refused, by a status that is not the one the entry needs or by the library's typed error, and some output's right value cannot be held at that scale. A condition of any other kind is not a refusal and ends the run |
| `:limit` | refused although every output's right value can be held: the scale is outside what the routine's own arithmetic reaches |
| `:no-input` | was not made: an input of the problem, a tolerance usually, is not a normal number at that scale |
| `:wrong` | succeeded with an output that is neither its right value nor lost |

**`:wrong` fails wherever it appears.** Every other class is declared by the
entry, as the ranges of the sweep it holds over, and the roster fails when a
run leaves its declared class in either direction. A range of `:exact` that
shrinks is a regression. A range of `:limit` that shrinks is an improvement,
and the entry has to be changed to record it.

An entry names the calls it answers for, and the name is accepted only if the
entry's own code calls it. The entries call the library directly, not through
the compiled callers the other tests use.

### What is declared

Exponents are of two. "Everywhere" is the whole sweep, 2^-1020 to 2^1020.

| Entry | `:exact` or `:close` | Outside it |
| --- | --- | --- |
| `norm2`, `norm-inf`; the power of two at a magnitude | everywhere | |
| Stored matrices and their products | to 2^1000 | refused at 2^1020, where a product is too large to hold |
| Gram refill | -520 to 400 | refused above, where a square is too large; **`:lost` at -600 and below**: the squares underflow and zeros come back |
| Augmented refill; row, column and symmetric scaling | everywhere | |
| Sparse factorization with positive pivots, solve, condition estimate | everywhere | |
| Cholesky form and its solve; the factorization that leaves pivots out | everywhere | |
| Signed-pivot factorization, growth and bound | to 2^1000; the bound `:looser` at -1000 and below | refused at 2^1020, where a pivot is too large to hold |
| Dense symmetric factorization, with a two-by-two pivot at the large end and at the small; orthogonal factorization and everything read from it | everywhere | |
| Symmetric eigenvalues; singular values | everywhere | |
| Sparse LU, solves, growth and bound | everywhere; the bound `:looser` at both ends, infinite at 2^1020 | |
| Refinement | -800 up | no input below |
| `cg!`, right-hand side scaled; `cg!`, matrix and right-hand side scaled together; `pcg!`, matrix, preconditioner and right-hand side scaled together | -800 up | no input below |
| `lsqr!`, damped, right-hand side scaled; `lsqr!`, damped, matrix and damping scaled | -800 up | no input below |
| `lsqr!`, matrix and right-hand side scaled together | -800 to 1000 | no input below, nor at 2^1020, where its right-hand side cannot be held |
| Dependent rows, threshold 2^-20 of the entries | -1000 up; a bound `:looser` at -1000, and infinite at 2^1020 where it cannot be held | no input at -1020 |
| Dependent rows, threshold zero | everywhere | |
| Damped least squares: the step with the columns' own damping, and with the residual scaled alone | everywhere | |
| Damped least squares: the step with a given damping, which is scaled as the matrix squared | -510 to 510 | no input beyond, where the damping cannot be held |
| Kept rows, their vectors, the bounds proved again and the certificate, threshold and floors 2^-20 of the entries | -800 up; a bound `:looser` at -1000 | no input at -1020 |
| The step on certified kept rows, the matrix and the residual scaled together | -1000 up | no input at -1020 |
| The step on certified kept rows, the residual scaled | -800 up; a bound `:looser` at -1000 and below | |
| Kept vectors refreshed for other values | -800 up; a bound `:looser` at -1000 | no input at -1020 |
| Supports classed and dependent sets labelled | -1000 up | no input at -1020 |
| Bounds on the largest singular value | everywhere | |

One limit is left, in bold, a property of the method as it is built:

- **A Gram refill returns what its products round to.** Above the range that is
  a refusal. Below it there is none: the squares become zero, and a
  factorization of the result then refuses a pivot.

Two more were measured by the first run of this roster and have since been
removed, which is what a declared `:limit` is for: conjugate gradients ended
`CgNonFinite` from a right-hand side of 2^510, and `CgNotPositive` on a small
positive definite operator, and the search for dependent rows reached half the
exponent range. Both now work on their data divided by a power of two.
[The record](batch16.md) has the measurements before and after.

### What the sweep does not ask

- It multiplies a whole problem by one factor. It does not ask what happens
  when the entries of one matrix lie far apart; the [scaling](scaling.md)
  module's own range tests ask that.
- Its inputs are normal numbers. A scale at which an input would be subnormal is
  `:no-input`, and what a routine does with subnormal entries is its own tests'.
- The law is relative to scale 1. A defect that is the same at every scale is
  the other tests' to find.
- The problems are of order three and four. Growth with size is not swept.
- **A norm that a status already bounds is held to its scaled value where it has one.**
  The residual of conjugate gradients and the residual and gradient of least squares on
  an operator are outputs of the law in five entries: a value wrong by any factor fails
  there. In the preconditioned entry the iterates at 2^1020 are only close, the residual
  is then another number, and it is held to its bound alone, which the status restates.

### Adding a call

Add its name to the entry whose problem exercises it, or write an entry: a
function of the exponent that builds the problem with `sweep-number`,
`sweep-array`, `sweep-stored` or `sweep-column-major`, makes the calls, and
returns each output with its power. Give the declaration as `:measure` and run
the suite: it prints the ranges it found and what decided each run that was not
exact. Read them before writing them down. A `:limit` is a finding about the
routine, to fix or to state in its contract, not a line to copy.
