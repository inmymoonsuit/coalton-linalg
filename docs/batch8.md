# Batch 8: the same elimination order, computed faster

State: completed; verified at safety 1 and safety 0 on 2026-11-20, Windows/x86-64,
SBCL 2.6.10 and the pinned Coalton compiler. The calling contract is unchanged and is in
[spd-ldlt.md](spd-ldlt.md).

## What changed, and what did not

**No call changes meaning and no result changes.** `analyze`, `analyze-csc`,
`analyze-quasidefinite` and their budgeted forms return the permutation, the factor
structure and the refusals they returned before, for every input. What changed is how the
order is found.

A consumer that analyses a pattern at every edit measured the analysis, not the
factorization, as its cost: 116 to 314 ms of minimum degree over three patterns on its
loose grid systems, against factorizations of about a millisecond. In this library's own
measurement the ordering was 70 to 88% of analysis, and most of that was the
representation:

- The elimination graph was stored whole, in one hash table a vertex, and eliminating a
  vertex of degree d looked up each of its d² pairs of neighbours.
- Each column was sorted twice with a general comparison sort, once in original numbering
  and once in elimination order.

Now:

- **The elimination graph is not stored.** An eliminated vertex leaves its neighbour
  set, which is its factor column and was being saved already, as an *element*. A live
  vertex holds two lists in space fixed at the start: the direct neighbours no element
  covers yet, and the elements it belongs to. Its degree is counted from those lists with
  stamps so that every vertex is counted once. An element is absorbed when one of its
  members is eliminated, since the new element contains it.
- **A hub keeps a row of flags.** A vertex adjacent to more than an eighth of the others
  belongs to hundreds of small elements that nothing absorbs, and re-reading them at each
  of its updates cost more than the hash tables had. Such a vertex keeps one byte a vertex
  saying which it is adjacent to now, and an update reads the new element against that
  row. The rows together take at most eight bytes for each stored neighbour.
- **Columns are sorted by one pass through the rows.** Filling the rows column by column
  leaves each row's columns in increasing order, and refilling the columns row by row
  leaves each column's rows in increasing order, in time linear in the factor.
- **The entry budget is a sum that was already there.** The diagonal, the columns
  emitted, and half the sum of the live degrees is the number of factor entries the
  eliminations so far commit to. It is checked after every elimination, so the verdict
  and the refusal's payload are what they were, and no column is allocated that the
  budget does not cover.

## How it was reached

The first version had no flag rows, and it made the hub pattern slower: 4.4 ms against
1.3 ms on a wheel of 400 spokes, while every other pattern improved. Flags for the hub's
direct neighbours alone brought that to 2.9 ms; the remaining cost was the hub's
elements, and the row became its whole adjacency. A shortcut for an ordinary vertex
already adjacent to every other live vertex was tried next and removed, because it changed
no measurement.

## Verification

Both complete serial checks pass on the final tree, `scripts/check.ps1 -Safety 1` then
`-Safety 0`, with every earlier roster unchanged. For this batch:

- **An independent elimination.** A reference written for the tests stores every fill
  edge in explicit neighbour sets and finds the least key by scanning every vertex; it
  shares no representation with the library. **It was first run against the code as it
  stood before this batch, and agreed**, so the reference is known to describe the old
  order. Against the new code it agrees on the permutation and on every column of the
  factor's structure in 1,188 analyses: paths, cycles, stars with the hub first and last,
  a wheel, two hubs over a path, complete graphs with and without one edge missing,
  grids, cliques joined at a vertex, an augmented pattern, 228 seeded graphs from 2 to 45
  vertices at densities from 4% to 90%, and six seeded graphs of 80 and 130 vertices; each
  with none, a third, half, all but one and all of its indices held.
- **The budget.** For each of those patterns a budget of exactly the factor's entries is
  accepted and one fewer is refused, naming that count and that budget.
- **The independent analysis verifier** of the complete check, a Python program with its
  own set-based elimination, passes on its 161 patterns, 13 of them above 4,096 vertices,
  including a star of 16,384 and a wheel of 4,097.
- **Outputs unchanged, byte for byte.** The fixture driver's 104 output files on thirteen
  public systems (permutations, factors, pivots, solutions, shifted variants) and the
  augmented driver's 288 output files on six more are identical to those written before
  the change.
- **Numeric rosters.** Every factorization test of the earlier batches runs on the new
  structure arrays and passes unchanged.

**Planted source controls.** Each of these edits was applied alone and the focused tests
run. Every planted source compiled, every one failed the tests, and the restored source
passed.

| Source | Planted defect | Caught |
| --- | --- | --- |
| ordering | A light vertex keeps direct neighbours that the new element covers | yes |
| ordering | Members of the new element are counted again through another element | yes |
| ordering | A vertex in two other elements is counted twice | yes |
| ordering | Absorbed elements stay in a light vertex's list | yes |
| ordering | The new element is not recorded for its members | yes |
| ordering | A light vertex's elements are not absorbed when it is eliminated | yes |
| ordering | A heavy vertex's elements are not absorbed when it is eliminated | yes |
| ordering | A heavy vertex's row keeps the eliminated vertex | yes |
| ordering | A heavy vertex's degree does not lose the eliminated vertex | yes |
| ordering | A heavy vertex counts itself among the new element's members | yes |
| ordering | A heavy vertex one neighbour short of all is taken as complete | yes |
| ordering | A heavy vertex's element list is never cleared of absorbed elements | yes |
| ordering | The eliminated vertex is not kept out of its own neighbours | yes |
| ordering | The budget leaves out the columns already emitted | yes |
| ordering | A budget of exactly the factor's entries is refused | yes |
| ordering | The refusal names one entry too many | yes |
| ordering | Among held neighbours of equal degree the later index is released | yes |
| ordering | The held neighbour released is the first met, whatever its degree | yes |
| ordering | A heavy vertex's row starts empty | yes |
| spd-ldlt | Row counts use original indices, not elimination positions | yes |
| spd-ldlt | A row entry records its own slot, not its place in the column | yes |
| spd-ldlt | A column's rows all land in its first place | yes |
| spd-ldlt | A column stores its own index in place of its rows | yes |

Two things the controls cannot see, and are not claimed: which vertices are treated as
hubs, and that the budget is refused early and not only at the end. Both change time and
memory, not results. The first is covered by the measurements below; the second by the
argument that a column's length is at most the live edges the previous check counted.

## Measurements

`analyze` on thirteen public patterns from a geometric constraint solver, safety 1, the
fixture driver, best of 20 after warm-up. "Before" was measured earlier the same day on
an idle desktop; "after" while another program was using about a fifth of the processor,
so the gains are if anything understated. Entries count the factor with its diagonal.

| Pattern | Order | Factor entries | Before, ms | After, ms | Ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| Chain, 25 | 400 | 1,098 | 0.275 | 0.133 | 2.1 |
| Chain, 100 | 1,600 | 4,398 | 0.752 | 0.430 | 1.7 |
| Chain, 200 | 3,200 | 8,798 | 1.628 | 0.957 | 1.7 |
| Wheel, 100 | 402 | 1,698 | 0.420 | 0.187 | 2.2 |
| Wheel, 400 | 1,602 | 6,798 | 1.154 | 0.482 | 2.4 |
| Grid, 15 × 15 | 1,290 | 9,418 | 3.056 | 1.010 | 3.0 |
| Grid, 21 × 21 | 2,562 | 21,631 | 8.046 | 2.459 | 3.3 |
| Frame, medium | 602 | 1,510 | 0.404 | 0.195 | 2.1 |
| Frame, large | 2,052 | 6,261 | 1.104 | 0.635 | 1.7 |
| Loose chain, 25 | 352 | 1,202 | 0.278 | 0.133 | 2.1 |
| Loose chain, 100 | 1,402 | 4,802 | 0.796 | 0.412 | 1.9 |
| Loose grid, 15 × 15 | 1,277 | 47,393 | 28.696 | 5.452 | 5.3 |
| Loose wheel, 100 (dense) | 354 | 40,908 | 59.230 | 2.391 | 24.8 |

Of the time after, the ordering is still the larger part: 5.0 of 5.5 ms on the loose grid
and 2.2 of 2.4 ms on the larger grid, measured separately; assembling the factor's
structure takes 0.02 to 0.5 ms.

The held order on the augmented patterns of [Batch 7](batch7.md), by the augmented driver
under the same conditions:

| System | Factor entries | Before, ms | After, ms |
| --- | ---: | ---: | ---: |
| Hub, 100 spokes | 2,362 | 0.58 | 0.30 |
| Chain, 25 | 1,958 | 0.57 | 0.27 |
| Chain, 100 | 7,808 | 1.54 | 0.84 |
| Frame, medium | 3,754 | 1.11 | 0.62 |
| Frame, large | 11,804 | 2.55 | 1.28 |
| Grid, 15 × 15 | 21,361 | 8.47 | 2.24 |

**What the tables say.** Sparse patterns with little fill gain about a factor of two,
which is the hash tables and the sorts. Patterns whose fill is dense gain most, because
dense fill is where elements absorb one another and nothing is read twice. The loose
grid is in between: its fill grows in fronts that overlap without nesting, and its
vertices read several large elements at each update.

## Addition: the fixed-pattern Gram analysis

The same consumer measured the Gram analyses as large as the ordering on its long chains.
[The Gram plan](gram.md) had the same two costs in another form, and the same answer
applied. **Its pattern, its refill and every refusal are unchanged.**

- Grouping the stored entries by inner index was a comparison sort through a closure. It
  is now a distribution by one base-2048 digit of the inner index at a time, least
  significant first, at most six passes. No array is indexed by an inner index, so a
  matrix with a vast unused height costs what it did.
- Each output column was collected and then sorted. The pattern is now filled by row: row
  i is appended to every column before it that it shares an inner index with, so each
  column receives its rows in increasing order and nothing is sorted. That uses the half
  of each group the counting pass does not, from the group's beginning to the entry's
  own position.

Verification, at both safety policies within the complete checks:

- The plan's exact-reference rosters pass unchanged: 76 seeded sources in both
  orientations with 912 refills against exact rational products, the lifetime, preflight,
  budget, overflow and encapsulation tests.
- **Scattered rows.** A^T A does not depend on what the rows are called. Twenty-four
  seeded sources are analysed twice, once with their twelve rows numbered 0 to 11 and
  once with the same rows scattered over indices that differ in every digit the
  distribution reads, up to two below the largest index. The patterns are equal and the
  refilled values are equal bit for bit. This test was run against the earlier code
  first, and passed.
- The Gram driver's 26 output files on the thirteen public systems are identical byte for
  byte to those written before.

Planted source controls, each applied alone, compiled and caught, with the restored source
passing:

| Planted defect | Caught |
| --- | --- |
| Only the lowest digit of an inner index is sorted on | yes |
| The placing pass reads the lowest digit whatever the counting pass read | yes |
| The sorted positions are read from the array a pass has just left | yes |
| A pass places equal digits from the end, reversing their order | yes |
| A group's beginning is never moved on | yes |
| A group's end stops one entry short | yes |
| The fill skips the first owner of each group | yes |
| The fill reads the owners after its own, not before | yes |
| The diagonal is not written | yes |
| The fill writes a row into a column once for each shared inner index | yes |
| The count takes a neighbour once for each shared inner index | yes |
| The count leaves out the diagonal | yes |
| A^T A is analysed as A A^T | yes |

One control was missed, and the source changed instead of the tests: resetting the marker
stamps between the two passes. It is not needed, because the first pass stamps an index
only with a column at or before it and the second asks about a row after it, so the reset
was removed.

Gram analysis on the public systems, the Gram driver, safety 1, best of 20. "Before" was
measured the day before on an idle desktop; "after" with another program using about half
the processor, so the gains are understated.

| System | Product | Order | Stored entries | Before, ms | After, ms |
| --- | --- | ---: | ---: | ---: | ---: |
| Chain, 25 | AᵀA | 400 | 1,048 | 0.156 | 0.077 |
| Chain, 100 | AᵀA | 1,600 | 4,198 | 0.592 | 0.127 |
| Chain, 200 | AᵀA | 3,200 | 8,398 | 1.234 | 0.266 |
| Wheel, 100 | AᵀA | 402 | 1,697 | 0.273 | 0.123 |
| Wheel, 400 | AᵀA | 1,602 | 6,797 | 1.022 | 0.187 |
| Grid, 15 × 15 | AᵀA | 1,290 | 4,818 | 0.686 | 0.148 |
| Grid, 21 × 21 | AᵀA | 2,562 | 9,642 | 1.439 | 0.311 |
| Frame, medium | AAᵀ | 602 | 1,608 | 0.283 | 0.159 |
| Frame, large | AAᵀ | 2,052 | 6,309 | 1.059 | 0.456 |
| Loose chain, 25 | AAᵀ | 352 | 1,152 | 0.148 | 0.087 |
| Loose chain, 100 | AAᵀ | 1,402 | 4,602 | 0.681 | 0.323 |
| Loose grid, 15 × 15 | AAᵀ | 1,277 | 10,191 | 0.832 | 0.193 |
| Loose wheel, 100 (dense) | AAᵀ | 354 | 35,705 | 2.192 | 0.169 |

## Limits

- The degree is exact. Its cost has no better worst-case bound than the earlier sum of
  squared degrees, and on a pattern like the loose grid the ordering still takes several
  times as long as the factorization. Going further means either an approximate degree,
  which changes the order and with it every fixed permutation a consumer holds, or
  treating vertices with identical neighbourhoods as one. Neither is done here.
- A vertex is treated as a hub by its initial degree only. One that becomes heavily
  connected through fill is not, and is updated through its lists.
- The measurements are on one host and one family of patterns, and the "after" column was
  taken on a shared desktop.
- One host and finite rosters establish neither universal behaviour nor portability beyond
  the supported runtime.
