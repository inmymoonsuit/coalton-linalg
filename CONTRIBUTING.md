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