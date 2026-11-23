# Working on coalton-linalg

- This is an independent mathematical library. No Blender or CAD procedures apply.
- Read README.md and the relevant contract (docs/spd-ldlt.md or docs/primitives.md) before changing code; read ROADMAP.md before extending scope.
- Inspect current Git state before editing. Preserve other people's work.
- Keep mathematical algorithms in Coalton. Common Lisp handles hosting/tooling; Python is an independent test oracle.
- Preserve the clean-room provenance in docs/provenance.md. Do not consult or copy numerical-library implementations or published algorithm code listings.
- Keep the main ASDF system free of test, benchmark, Python, and application dependencies.
- Document types, storage, ownership, mutation, numerical assumptions, refusals, and complexity with each public interface.
- Keep analyze/factor!/solve! separate. Do not silently add damping, tolerances, pivoting, or geometry policy.
- Ask the four questions in CONTRIBUTING.md of every numeric routine you add or change: an intermediate that leaves the range, two arguments that are one array, work for a step not taken, anything allocated or written before a refusal. A new exported call fails the aliasing and scale-sweep rosters until it has an entry; add the entry, do not exempt it for convenience.
- Test through the public interface and independent mathematical identities. For a regression, observe its test fail without the fix; retain a reproducible control.
- Run scripts/check.ps1 serially at safety 1, then safety 0 for numeric or ABI changes. Use the existing configured toolchain; do not install dependencies without authorization.
- Keep generated files under build/, local machine paths in ignored config.local.psd1, and recorded results distinct from current measurements.
- Scope present portability to 64-bit SBCL. Measure before making allocation or speed claims.
- This project is MIT-licensed. Preserve the provenance record and record permitted third-party code at its point of use.