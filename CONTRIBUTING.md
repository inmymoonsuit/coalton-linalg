# Contributing

Start with the relevant contract linked by the [interface reference](docs/api.md).
Changes to numeric behavior, ownership, or a refusal are interface changes and
must update that contract and its tests together. Keep the symbol index aligned
with package exports. Add shared mathematical terms to [CONTEXT.md](CONTEXT.md)
when they need a definition; keep implementation details in the contracts.

Use the existing configured toolchain and `scripts/check.ps1`. Run one build or
test job at a time on a shared desktop. A mathematical change needs independent
reconstruction/residual checks and a test that exposes the previous defect.
Keep allocation checks distinct from numerical correctness and timing.

## Four questions for every numeric routine

An outside audit found eight defects that every test then in the suite passed
([record](docs/batch14.md)). Each was the answer to one of four questions that
had not been asked of the routine. Ask all four of every routine you add or
change.

1. **Can an intermediate value overflow or underflow when the result is
   representable?** A sum of squares, a product of two scale factors, the
   determinant of a two-by-two block: each can leave the range of the format
   while the value the caller wants is well inside it. Scale before combining,
   or carry the exponent apart from the significand.
2. **What happens when two array arguments are one array?** The call either does
   exactly what it does with two arrays of equal contents, or refuses before it
   writes. A call that reads an input after it has begun to write its output
   needs a copy, or the refusal.
3. **Is work done for a step that will not be taken?** Preparing the next step
   after the last one allowed can cost an application of the caller's operator,
   and can turn a finished answer into a refusal when that preparation fails.
4. **Is anything allocated or written before a refusal?** Sizes, budgets and
   finiteness are checked before the first allocation and the first write, so a
   refused call leaves its outputs as they were and takes no storage.

The first two are asked by test, of every member of the family at once. The
[aliasing roster](tests/aliasing.lisp) and the [scale-sweep roster](tests/scale-sweep.lisp)
read their members from the sources, so an exported call joins by being
declared, and the suite fails until the roster has an entry for it or a written
reason for having none. [The rosters' page](docs/rosters.md) says what each
entry declares and how to add one. The last two have no roster: answer them in
the routine's own tests, with a refusal test that checks the outputs and the
allocation.

The primary ASDF system must remain usable without Python, fixtures, benchmarks,
or a CAD application. Add a new mathematical module when a concrete use needs it;
do not create empty packages or speculative backend abstractions.

The independently derived algorithms and tests have a clean-room provenance
requirement. Consult language documentation and mathematical prose, not other
numerical-library implementations or algorithm code listings. Record any permitted
new source and its license before adopting code. See [provenance](docs/provenance.md).

Generated logs, factor artifacts, compiled caches, and local configuration are not
source. Keep them under ignored locations. A concise evidence record should state
the revision/toolchain, command, counts, failures, and limitations.

Contributions are under the project's [MIT license](LICENSE). Use the
[setup recipe](docs/setup.md) and preserve its explicit dependency registration.
Hosted CI and additional platforms remain roadmap work; no automatic publishing
or external integration is configured.