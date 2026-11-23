# Batch 15: the audit's two questions asked of every call

State: completed; verified at safety 1 and safety 0 on 2026-11-22, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract is [rosters.md](rosters.md);
the contracts it touches are [dependent-rows.md](dependent-rows.md),
[conjugate-gradient.md](conjugate-gradient.md) and [gram.md](gram.md).

## What this is

[The audit of Batch 14](batch14.md) found eight defects that every test then in the suite
passed. None was in a place the tests had looked and got wrong; each was the answer to a
question that had not been asked of that routine. Fixing the eight left the questions
unasked of everything else.

This batch asks two of them of every exported call at once. A roster reads its members
from the sources, so a call joins by being declared, and holds each to one stated
contract. The other two of [the four questions](../CONTRIBUTING.md#four-questions-for-every-numeric-routine)
are written into the contributor notes and stay with each routine's own tests.

## What was built

- **The aliasing roster**, [`tests/aliasing.lisp`](../tests/aliasing.lisp): 25 exported
  calls take two or more arrays. Each pair of arguments that could be one array is
  called both ways, 34 pairs in all, and must either do exactly what it does with two
  equal arrays or refuse before writing.
- **The scale-sweep roster**, [`tests/scale-sweep.lisp`](../tests/scale-sweep.lisp): 82
  exported calls take or return floating-point numbers, or return a status. 23 entries
  answer for all of them, each run at 25 scales from 2^-1020 to 2^1020, 575 runs, every
  output compared with its value at scale 1 times the power of the factor it must move
  by.
- **Membership from the sources.** Both read every source file's export list and
  declarations. A member without an entry fails the roster, and so does an entry that
  names no member or, in the sweep, an entry that answers for a call its code does not
  make.

## What the rosters found

**Aliasing: nothing.** Nineteen pairs behave as two equal arrays and fifteen are refused
before a write, each as its contract already said. The solves copy their right-hand side
before writing, and the two iterative methods save it, which was the audit's fix.

**The scale sweep: one silent wrong answer, and the range of every routine.**

| Finding | Measured | What was done |
| --- | --- | --- |
| `find-dependent-rows!` at small scales | With entries of size 2^-600 to 2^-1000 and a threshold 2^-20 of them, `RowsFound` with all four rows left out where one is dependent. The squares of the rows and of the threshold had underflowed to zero, and zero is at most zero. | Refused by name: `RowsBadThreshold` for a positive threshold whose square is not a normal number, `RowsNonFinite` under a threshold of zero for a row whose square is not. |
| `cg!` and `pcg!` at the large end | `CgNonFinite` from a right-hand side of 2^510, with a solution that can be held | Already stated in the contract as a limit of sums taken as they stand. Now measured and pinned. Not changed. |
| `cg!` with the operator scaled too | Holds from 2^-512 to 2^300. Below, `CgNotPositive` for a positive definite operator: the image of a direction underflows and its curvature is computed as zero. | Stated in the contract and pinned. Not changed. |
| Gram refill at small scales | Squares of entries of 2^-600 and below come back zero, with no refusal | Stated in the contract and pinned as `:lost`. Not changed. |
| Everything else | Bit for bit across the range, or within the entry's tolerance at one end of it, or refused only where an output cannot be held | Pinned |

Thirteen of the twenty-three entries hold over the whole sweep. Three more are refused only
at 2^1020, where an output is too large to hold, and refinement holds wherever its
tolerance is a number the format has. The table of every entry's declared ranges is on
[the rosters' page](rosters.md#what-is-declared).

The first three rows are one cause seen three times: a method that squares its data
reaches half the exponent range. The audit handled the small end of conjugate gradients
by carrying the exponent of a tiny inner product apart from its significand. The large
end, the scaled operator, and dependent rows are the same shape and are left as limits
here, on the principle that a change to a method's arithmetic is proposed before it is
made. [The roadmap](../ROADMAP.md) lists them.

## The refusal that replaced a silent answer

`find-dependent-rows!` leaves a row out when its squared distance from the rows kept
before it is at most the threshold squared. Both sides are squares, formed in `A Aᵀ` and
as `threshold * threshold`.

- A positive threshold below about 1.5 × 10⁻¹⁵⁴ has a square below the least normal
  number. It is now `RowsBadThreshold`, as a threshold whose square is not finite
  already was.
- Under a threshold whose square can be held, a row whose own square underflowed is
  shorter than the threshold. Leaving it out is what the contract says, so nothing is
  refused.
- Under a threshold of zero, a row holding a nonzero entry whose square is below the
  least normal number is `RowsNonFinite`. The search would otherwise leave it out as
  dependent.

`RowsNonFinite` is the nearer of the two statuses the type has, and its meaning in the
contract is widened to say what it now also reports. A status of its own would be the
clearer interface and would break every caller that matches the type exhaustively; that
is left for a decision.

One existing test changed its reason. A two-row matrix with entries 1e-160 and 1e150
under a threshold of zero was `RowsNonFinite` because the factorization's first
multiplier overflowed. It is `RowsNonFinite` now because the first row's square is too
small, which is met first. **No test here now reaches the factorization's own refusal
through the search**: with every row's square a normal number the first multipliers are
bounded, and a matrix that breaks that bound by rounding was not found.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`, with every earlier roster unchanged.

- **Red first.** The sweep's dependent-rows entry was run before the fix and reported
  the wrong count of rows left out at 2^-600, 2^-800 and 2^-1000; the same runs are
  `:limit` after it. The direct tests of the two refusals are in the dependent-rows
  suite.
- **In-file controls.** Each roster carries its own: a member the sources would add is
  reported until answered for; a call that reads what it has overwritten, and a refusal
  after a write, are reported by the aliasing verdict and a call that is right with one
  array is not; each class of the sweep is produced by a planted run, a value off its
  scaled value, one the format cannot hold, a refusal of each kind, a bound that came
  out larger and one that came out smaller; and an entry answering for a call it does
  not make is refused.
- **Planted source defects for the aliasing roster.** Twenty-four edits to the library's
  sources, each applied alone and the roster run. Every planted source compiled, every
  one failed the roster, and the restored sources passed.

| Planted defect | Number | Caught |
| --- | --- | --- |
| The iterate's array is read as the right-hand side, in conjugate gradients and in refinement | 2 | yes |
| A refusal of an alias is removed: the products, the Gram refill, and each of the twelve in the scaling module | 14 | yes |
| A solve clears its output before reading the right-hand side: the three sparse symmetric solves, both LU solves, the dense symmetric and the orthogonal | 7 | yes |
| A call taking two arrays is exported with no entry | 1 | yes |

- **Planted source defects for the scale sweep.** Eighteen edits, each applied alone and
  both rosters run. Every planted source compiled, every one failed the sweep, and the
  restored sources passed.

| Planted defect | Number | Caught |
| --- | --- | --- |
| The dependent-rows search takes a positive threshold whose square underflowed; under a threshold of zero it leaves out a row whose square underflowed | 2 | yes |
| Conjugate gradients take an inner product, or a sum of squares, that underflowed as it stands | 2 | yes |
| `norm2` squares its scale before taking the root | 1 | yes |
| A matrix-vector product that is not finite is returned | 1 | yes |
| Refinement compares a residual with the tolerance plus an absolute number | 1 | yes |
| A sparse pivot is judged against an absolute number: positive, signed, and in the Cholesky form | 3 | yes |
| Sparse LU never pivots on an entry below an absolute number | 1 | yes |
| The dense symmetric factorization does not move a large two-by-two block away from overflow | 1 | yes |
| The orthogonal factorization sums the squares of a large column, or of a small one, as they stand | 2 | yes |
| The symmetric eigenvalue routine rotates a large matrix, or a small one, as it stands; the singular value routine forms the moments of large columns as they stand | 3 | yes |
| A call taking floating-point arrays is exported with no entry | 1 | yes |

  **The first run of these missed three, and that is why three fixtures are what they
  are.** One plant was aimed at a tie-break and changed nothing; it was replaced. The
  other two were real: the dense symmetric matrix never formed the two-by-two pivot whose
  guard was planted against, and the symmetric eigenvalue matrix had a diagonal of one
  sign, whose differences cannot overflow. Each now has a matrix that reaches its guard,
  the dense one asserting that its first pivot is a block. **One planted defect was
  withdrawn as unreachable**: the same factorization's guard for a block of entries below
  twice the least normal number. The sweep gives a routine only normal numbers, and with
  those the guard changes at most the last bit of a product that is already subnormal;
  the roster passed with it removed. The audit's own test, on subnormal entries, is what
  holds it.

## Limits

- **The sweep multiplies a whole problem by one factor.** A matrix whose own entries lie
  far apart is not swept; the scaling module's range tests ask that of the scaling
  module, and nothing asks it of the factorizations.
- **The sweep's inputs are normal numbers.** What a routine does with subnormal entries
  is not swept, and a guard that only they reach is not held by it.
- **A defect that is the same at every scale is invisible to the sweep**, which compares
  each scale with scale 1. It is not a test of correctness.
- **The problems are of order three and four**, chosen so that every factorization
  pivots or reorders. Growth with size, and fill, are not swept.
- **That an entry answers for a call is checked by its code naming the call**, not by the
  call's being reached on every path, or with arguments that exercise it.
- **Arrays inside handles, and functions passed as arguments, are not array parameters**
  to the aliasing roster. A caller's operator that writes over its own input is the
  caller's.
- **The factorization's refusal inside the dependent-rows search has no test that
  reaches it**, as said above.
- **The limits measured are stated, not removed.** Conjugate gradients at the large end
  and with a scaled operator, dependent rows beyond half the range, and the Gram
  refill's silent underflow are as they were.
