# Batch 14: an outside audit at the extremes of the range

State: completed; verified at safety 1 and safety 0 on 2026-11-22, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contracts it touches are
[conjugate-gradient.md](conjugate-gradient.md), [refinement.md](refinement.md),
[indefinite-ldlt.md](indefinite-ldlt.md), [jacobi.md](jacobi.md),
[orthogonal.md](orthogonal.md), [lu.md](lu.md) and [quasidefinite.md](quasidefinite.md).

## What this is

Robert Smith audited the library and sent the result as pull request 1, three commits,
"Fix false convergence and factorization errors at extreme scales". Every change in it is
his. This record says what the audit found, how each finding was confirmed here, and the
one thing that had to be changed for his code to work on this library's pinned compiler.

## What the audit found

Each finding came with a test. **Each of those tests was run here against the sources as
they were before the audit, and each failed**, so every item below is a defect that was
present, shown by its own test:

| Where | What was wrong before | What the test saw on the old source |
| --- | --- | --- |
| `cg!`, `pcg!`, `refine!` | A right-hand side that is the same array as the iterate was overwritten by the first step, and every later residual was formed from the wrong vector | The solution of `2 x = 2` came back as 0.5 |
| `cg!`, `pcg!` | The residual's norm was the root of a sum of squares, and a square below the normal range is zero: a residual near 10⁻¹⁶⁰ had norm 0 and met any tolerance, a tolerance of zero included | A call with no steps allowed reported convergence on a residual of 2⁻⁵³⁰ |
| `pcg!` | After the last step allowed, the preconditioner was called again to prepare a step that would never be taken, and its refusal became the call's outcome | `CgPreconditionerRefused` where the answer is `CgLimit` with its residual |
| Dense indefinite LDLᵀ | The two-by-two pivot's scaled determinant was multiplied by the off-diagonal entry, which overflows for entries near the top of the range; the quotient then came out as zero, finite and wrong | A system whose solution is (2/17, 8/17) was solved as (0, 0), with success reported |
| Jacobi SVD | The whole matrix was scaled only above 2⁵⁰⁰ or below 2⁻⁵⁰⁰, so between those a small singular value's square could still underflow | A singular value of 2⁻⁵³⁸ beside one of 2⁻⁵⁰⁰ was reported as 0, and the rank with it |
| Householder QR | The reflector was formed from a column norm already rounded to the subnormal grid | For a column of the smallest positive number, Q applied to a unit vector was off by 0.29 |
| Sparse LU | With a threshold small enough, a pivot far below a candidate is accepted and the multiplier overflows; it was stored unchecked | `factor!` did not refuse a multiplier of 2¹⁰²⁴ |
| Augmented plan | A budget too small was refused only after the source pattern had been copied | Sixteen refusals allocated 4 MB |

The first seven are wrong answers or wrong outcomes; the last is a cost. Most of them are
in code that an independent held-back review had already passed, which is a fair measure
of what another reader with a different eye finds: all eight sit at the edges of the
floating-point range or at an aliasing the tests never offered.

## How it was fixed

The changes are the pull request's, and its commit messages and the contracts describe
them: the right-hand side is copied into storage the workspace already owns when, and
only when, it is the iterate's array; a sum of squares that is zero or subnormal is
formed again from lifted components; an inner product that underflowed is carried as a
value and a separate binary exponent until its ratio is taken; the two-by-two pivot's
divisor and right-hand side are scaled together; the SVD scales its inner products and
not its stored columns, and normalizes each final column at its own scale; the reflector
is formed before its head is rounded back to the input's scale; a multiplier that is not
finite is `LuNonFinite`; and the budget is checked before the copy.

## What had to change for this library's compiler

**The pull request did not pass here as sent.** Seven of its eight tests passed and the
tiny-residual test failed, on its own branch, at both safety levels.

The cause is not in the pull request's logic. Its three small helpers for scaled
arithmetic lived in a new package, `coalton-linalg/internal/float`, and were inlined into
the conjugate gradient file. Each holds a mutable cell. **The Coalton compiler this
library pins, revision `b6e284b`, miscompiles a function that holds a cell when it
inlines it from another package**: the cell's variable is left unbound at the point of
use. The compiler reports that as two warnings, "undefined variable", and the build goes
on. The running total inside the norm's recomputation was then lost, and the norm came
out as zero, which is the defect the helper exists to remove.

Three measurements settled it:

1. With the three helpers as ordinary functions, not inlined, there are no warnings and
   every numerical assertion passes, and the allocation measurement fails at about 49
   bytes a call, since a number returned from an ordinary function is boxed.
2. With the helpers inlined from within the same file as their caller, there are no
   warnings, every test passes, and nothing is allocated.
3. Removing only the explicit inlining at the two call sites changes nothing.

So the helpers were moved into the conjugate gradient file, their only user, unchanged.
The float package keeps the one constant two files share. The author presumably built
with a later compiler; that is for him to say, and the compiler revision is the first
thing to check when a result here differs from one elsewhere.

**A guard was added so that this cannot pass silently again.** The build script that
every check and every focused run goes through now fails if the compiler prints a
warning while building the library. Its control: run against the pull request's sources
as sent, it refuses the build, where before it built and ran.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`, with every earlier roster unchanged, the bit-for-bit recomputation of ninety-
six conjugate gradient solves among them: on ordinary data the audited code computes
exactly what the earlier code did.

- **The audit's eight tests**: red on the old sources, each with the observation in the
  table above, and green on the final tree.
- **Allocation**: the audit's own measurement on the banded inner products, and every
  earlier one, at zero bytes, at both policies.
- **Planted source controls on the conjugate gradient source**, which the audit reworked
  most: Batch 13's twenty-eight, with the eight whose sites moved brought up to date, and
  six for what the audit added. Each edit was applied alone and the conjugate gradient
  tests run. Every planted source compiled, every one failed the tests, and the restored
  source passed. The six:

| Planted defect | Caught |
| --- | --- |
| A right-hand side that is the iterate is not saved | yes |
| Refinement does not save a right-hand side that is the iterate | yes |
| A sum of squares that underflowed is taken as it stands | yes |
| An inner product that underflowed is taken as it stands | yes |
| A ratio of inner products ignores their separate exponents | yes |
| A direction is prepared after the last step allowed | yes |

## The pin, afterwards

The defect was already fixed upstream. On the day of the audit its author changed the
compiler so that a local variable that is assigned is no longer treated as a constant
when an inlined body is propagated, which is what had left the assignment with no
target. **The pinned compiler was then moved**, from `b6e284b` to `deb300be`, the head of
the upstream repository on 2026-11-22, which contains that change.

- The pull request's sources **as sent**, helpers in the separate package, were built
  with the new compiler under the warning guard: no warning, all eight of its tests
  pass, and nothing is allocated. With the old compiler that same tree is refused.
- Both complete checks pass on the new compiler at safety 1 and safety 0, with a cache
  that has held only this compiler's output and with the content lock rewritten for it:
  every Coalton file's hash, and the list of systems the build loads, which is one
  shorter, since the compiler no longer depends on `source-error`. The old compiler's
  checkout no longer matches the lock and is refused by the content check.
- The helpers stay where the integration put them. That arrangement is right under both
  compilers, and the guard stays whatever the compiler.

The old compiler's checkout is left in place beside the new one, since other work on the
same host is pinned to it.

## Limits

- **No planted controls were written for the dense changes**, in the indefinite
  factorization, the SVD and the reflector, or for the one-line changes in LU and the
  augmented plan. Their evidence is the audit's tests, each seen to fail without its fix.
  Sparse LU's forty controls of Batch 12 still apply to the source as it stands and were
  not run again.
- **The scaled-arithmetic helpers do not end on an argument that is zero or not finite.**
  Every caller excludes both before calling, by checks that are in the source for other
  reasons; nothing in the helpers themselves enforces it.
- **The compiler defect is known here only through this one case.** It is fixed in the
  compiler now pinned. The guard catches a warning; a miscompilation that produced none
  would be caught only by a test.
- **The new compiler is the upstream head of one day**, chosen because it holds the fix
  and verified by this library's own checks and nothing else. It is not a release.
- The tests added are at the edges the audit chose. That the same code is right
  everywhere between them rests on the earlier rosters passing unchanged.
- One host and finite rosters establish neither universal behaviour nor portability beyond
  the supported runtime.
