# Batch 3c: fixed-pattern Gram products

State: completed; verified at safety 1 and safety 0 on 2026-11-19, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The contract is [gram.md](gram.md).

## What was delivered

`coalton-linalg/sparse/gram` prepares a reusable plan for AᵀA or AAᵀ from a checked CSC
matrix and refills caller-owned values. The output is the lower triangle in canonical CSC with
every diagonal, the pattern `spd-ldlt` analyzes. Structural entries stay when their value
cancels to zero. Budgets and limits are shared with SPD analysis through a small private
`sparse/limits` package.

An optional driver, `scripts/check.ps1 -Task Gram`, writes each system's product, a rescaled
second refill and one measurement line for an independent reviewer. A compiled example,
[normal-equations.ct](../examples/normal-equations.ct), runs plan, analysis, refill, factor
and solve together.

## One change made after the first implementation

The first implementation scanned each inner-index group from its start and skipped the owners
before the current column. Within a group the owners ascend, so the scan now starts at the
entry's own sorted position. Work per inner index of degree d fell from d² to d(d+1)/2 in both
analysis and refill. The exact-reference tests and the allocation measurements passed before
and after the change with no edit to either.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`, with every earlier roster unchanged. For this batch:

- **Exact references.** 76 seeded and constructed cases in both orientations, 912 refills,
  each stored value compared with an independent exact rational product. Dyadic inputs keep
  the comparison exact.
- **Accumulation order.** Two orientations of a case whose rounded result depends on the order
  confirm the documented increasing-inner-index order.
- **Calling contract.** Lifetime, preflight preservation of output and input bits, both refill
  entry points, limits through the maximum UFix dimension, and eight arithmetic refusals each
  followed by a successful reuse of the same plan.
- **Encapsulation.** Three external compile probes that must be rejected and one public control.
- **Allocation.** Eight warmed paths, 10,000 calls each after 10,000 warm-up calls: zero bytes
  at both policies. Raw and checked refills, both orientations, non-empty and empty.
- **Driver.** Exact small products in both orientations, a source column with no entry, a
  typed refusal file, artifact replacement on rerun, and a refusal to guess an unnamed
  orientation. A guard fails if a new `LinalgError` constructor has no recorded spelling.

**Planted source controls.** Each of these edits to `gram.ct` was applied alone and the
focused tests run; every one failed all three test groups, and the restored source passed.

| Planted defect | Caught |
| --- | --- |
| Refill skips each entry's own product | yes |
| A group's end stops one entry short | yes |
| The mandatory diagonal is not collected | yes |
| AᵀA is analyzed as AAᵀ | yes |
| Collected rows are left unsorted | yes |
| Refill does not clear the output first | yes |

## Measurements

Observations, not targets. Thirteen sparse Jacobian-like systems supplied by the consuming
application's reviewer and kept outside this repository; the generic rosters above are the
library's own evidence. Times are the best of 20 after one warm-up call, in milliseconds.
Refill bytes are over 1000 warmed refills. One process at a time under a 1024 MB heap cap. The
desktop was not idle: 7% processor use and 96% of its commit limit in use by other programs at
the start.

| System | Product | n | Stored entries | Analyze, s1 | Refill, s1 | Refill, s0 | Refill bytes | `factor!`, s1 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| chain25 | AᵀA | 400 | 1,048 | 0.156 | 0.0064 | 0.0052 | 0 | 0.0052 |
| chain100 | AᵀA | 1,600 | 4,198 | 0.592 | 0.0268 | 0.0199 | 0 | 0.0218 |
| chain200 | AᵀA | 3,200 | 8,398 | 1.234 | 0.0559 | 0.0398 | 0 | 0.0460 |
| grid15x15 | AᵀA | 1,290 | 4,818 | 0.686 | 0.0308 | 0.0221 | 0 | 0.1081 |
| grid21x21 | AᵀA | 2,562 | 9,642 | 1.439 | 0.0618 | 0.0439 | 0 | 0.3397 |
| wheel100 | AᵀA | 402 | 1,697 | 0.273 | 0.0099 | 0.0077 | 0 | 0.0095 |
| wheel400 | AᵀA | 1,602 | 6,797 | 1.022 | 0.0450 | 0.0309 | 0 | 0.0396 |
| loosechain25 | AAᵀ | 352 | 1,152 | 0.148 | 0.0062 | 0.0047 | 0 | 0.0063 |
| loosechain100 | AAᵀ | 1,402 | 4,602 | 0.681 | 0.0255 | 0.0199 | 0 | 0.0245 |
| medium | AAᵀ | 602 | 1,608 | 0.283 | 0.0102 | 0.0076 | 0 | 0.0073 |
| large | AAᵀ | 2,052 | 6,309 | 1.059 | 0.0392 | 0.0287 | 0 | 0.0327 |
| loosegrid15x15 | AAᵀ | 1,277 | 10,191 | 0.832 | 0.0490 | 0.0359 | 0 | 1.3321 |
| loosewheel100 | AAᵀ | 354 | 35,705 | 2.192 | 0.1212 | 0.0880 | 0 | 2.9859 |

- **Refill against factorization.** The last column is `factor!` on the same system's reference
  pattern, measured in the same session by the fixture driver. A refill costs between 0.04 and
  1.4 times that factorization on these thirteen.
- **Agreement with an independent product.** The reviewer's files hold each system's product
  as computed by SciPy. Every entry stored in both agrees bit for bit. This library stores up
  to 850 further entries a system; each is an exact zero, a structural pair whose products
  cancel, which the reference had dropped.

## Observed beside this batch

SPD analysis of the two densest factors here takes longer under the compact ordering of Batch
3a than the bitmap it replaced did on the earlier prototype: 28.5 ms against 8.5 ms for
loosegrid15x15 (47,393 factor entries) and 53.6 ms against 7.6 ms for loosewheel100 (40,908).
The chains went the other way: 1.6 ms against 13.9 ms at 3,200 unknowns. The earlier figures
come from the prototype's own record, taken in a different session of the same machine. The
ordering looks up every pair of neighbours in a hash table at each elimination, which costs
most where fill is dense. Analysis runs once for a pattern. Nothing was changed for it in this
batch.

## Limits

- A Gram product squares the conditioning of its source. This package forms the product and
  decides nothing about rank or definiteness.
- Work follows the sum over inner indices of d(d+1)/2. A single dense row or column of the
  source makes its product dense, as loosewheel100 shows.
- One host and a finite roster establish neither universal speed nor portability beyond the
  supported runtime.
