# Initial review snapshot

State: verified locally, 2026-11-19.

This snapshot supplies an independent F64 linear algebra library: checked
rectangular CSC construction and operations, stable vector norms, and sparse
SPD LDLᵀ with reusable analysis and workspace. It is intended for mathematical
and Coalton interface review before further algorithms or application integration.

## Interface changes

The public factor/workspace types now have private nominal constructors around
private records. Plain aliases would not hide Coalton's type-directed field
access. Eight real external caller probes compiled before this change and were
rejected afterward; a valid public caller compiled and ran in both cases.
Four additional checks ensure the replacement constructors remain unexported.

FactorStatus and SolveStatus are distinct types with five named outcomes each.
External compiled callers exercise all ten outcomes, and four invalid
cross-family assignments are rejected after a valid matching control compiles.
An earlier successful result retains its outcome while its borrowed factor can
later become invalid; both the documentation and tests cover this distinction.

This is a deliberate pre-release interface change: factor/solve/result-code
became factor!/solve!/factor-status, solve returns a named status, and callers
cannot construct internal workspaces through the public interface. All examples,
tests, and optional tools use the new interface. The original prototype remains
separate and unchanged.

The arithmetic and elimination ordering are unchanged. Independent review found
one stale numeric comparison in a fixture callback; the compiler also rejected
it, and it was corrected before the final runs.

## Reproducible setup

The [setup recipe](setup.md) was exercised using a separate export of staged
source, a new detached checkout of the pinned public Coalton commit, and all
20 public Quicklisp downloads. Every download matched its recorded size/SHA256;
the extracted installation matched all 18 components and 1,314 files in the
expanded content lock. Existing installed SBCL and Python were used; no system
installation was modified. The SBCL MSI archive hash remains unestablished;
installed binary content is pinned.

Registration now resides in the repository. No private setup file or existing
compiled dependency cache was used in the final fresh-source run. Inherited
ASDF source/cache settings were overridden, and their designated cache remained
unused. A planted Quicklisp startup hook was refused before execution.

The first catch-all cache layout repeated long absolute source paths and the
fresh run exited with SBCL access-violation status -1073741819. Changing only
to a fresh short cache location passed. Explicit source-root mappings now
retain relative paths under separate cache namespaces: the final long-location
run passed with 338 FASLs and a longest completed FASL path of 240 characters.
Three planted old catch-all mappings were rejected by the real ASDF mapping
guard. These observations do not establish a general path threshold or the
crash's precise internal cause.

## Final validation

The final exported tree ran the complete check command serially at safety 1
from a nonexistent compiler cache, then at safety 0 using those freshly compiled
dependencies. Both runs passed:

- 20 baseline and 336 seeded matrix cases, 252 scaling/permutation relations,
  and 96 repeated refactor cases.
- Independent primitive checks, all outcome branches, ownership/refusal tests,
  compiled examples, and negative external caller checks.
- Eleven repeated factor/solve success/refusal paths each reported zero bytes
  for 10,000 calls after warmup; the computed-shift Coalton caller also reported
  zero bytes for 10,000 calls. The allocating control was detected.
- Dependency source resolution matched all 28 pinned loaded systems.
- 75 fixture files reproduced; the independent artifact verifier rejected ten
  planted defects, the stress checker six, and the toolchain checker six.
- 74 baseline and 2,016 seeded factor/solution/refusal artifacts were byte-identical
  across the final safety policies and the previous snapshot's safety-0 artifacts.

For the seeded cases, maximum normwise backward error was
1.1696176720330866e-16 and maximum relative forward error was
7.450580596923828e-9. The cases have at most 12 rows; these figures do not
establish accuracy for arbitrary large or ill-conditioned systems.
Allocation results apply to the documented bounded measurement protocol.

[Machine-readable evidence](initial-review-results.json) records the checks and
normalized-LF SHA256s of all 105 executable, test, fixture, and toolchain inputs
matched to the tested export. Documentation and the configuration template are
checked separately; the template's final change was a cache-ownership comment.

## Review and remaining scope

The [README](../README.md) identifies the small implementation surface and review
questions. The library is [MIT-licensed](../LICENSE), with its
[provenance](provenance.md) retained. The concise AGENTS.md is shared project
guidance; local configuration, dependencies, caches, logs, and generated artifacts
are excluded from Git.

Scalable symbolic graph storage, QR, hosted CI, and additional hosts remain
roadmap work. The current ordering still uses an n²-byte adjacency bitmap.
No Blender integration or SolveSpace promotion claim is part of this snapshot.