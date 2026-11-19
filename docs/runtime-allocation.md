# Allocation measurement: separating the caller from the runtime

The initial 131–164 KB readings did not establish that factorization allocated.
The counter includes the whole process, and garbage collection starts background
work. A long warmup merely moved that work outside the measurement; it did not
identify its source.

On 2026-11-19, three fresh Windows x86-64 SBCL 2.6.10 processes ran
`tools/runtime-probe.lisp`, with the pinned Coalton checkout in release mode and
library safety 1. Each created the same three-variable factorization workspace,
ran full GC, and measured without any factor-call warmup. Reporting and workspace
allocation were outside the measured intervals. A factor batch was 10,000 calls;
the idle control spun for 50 ms without calling the library.

| Fresh process | First measured operation | Bytes reported | First factor batch | Later factor batches |
| --- | --- | ---: | ---: | ---: |
| Factor first | Factor batch | 98,288 | 98,288 | 0 |
| Idle first | Idle control | 98,272 | 0 | 0 |
| Runtime finalizer thread paused | Factor batch | 0 | 0 | 0 |

All three processes detected the same deliberately allocating control:
131,008 bytes reported over 10,000 cons-cell calls. The counter is not disabled by
pausing the background thread. The paused case stopped only the disposable probe
process's finalizer thread and restored it with `unwind-protect`; neither the library nor ordinary ASDF test-system calls change that thread.
The canonical disposable check/fixture runners now apply this reviewed isolation
policy during bounded measurements and test failure-path restoration; see toolchain.md.

The 50 ms idle interval is a deliberately longer observation window, not a
duration-matched baseline: the first factor batches took 294–319 microseconds.
No baseline bytes are subtracted. The useful observation is that approximately
the same transient appeared without any factor calls, while pausing the runtime
thread removed it from the first factor batch.

These controls identify background-runtime contamination of the original
process-wide measurement. They do not identify the individual callback or
allocation site inside that thread. SBCL's versioned
[GC implementation](https://raw.githubusercontent.com/sbcl/sbcl/sbcl-2.6.10/src/code/gc.lisp)
confirms that `get-bytes-consed` reads process-wide usage and that post-GC work
notifies the finalizer thread asynchronously. This source reading supports the
experiment; the controls provide the causal distinction for this installation.

The result is **zero bytes reported for the first 10,000 calls on this prepared
factorization workspace when background work is separated**. It is not a claim
that loading Coalton, creating the symbolic analysis, or starting the process
allocates nothing. Nor is the region-accounted counter an exact per-object trace:
the cons-cell control's reported total itself shows that it should not be used
to count individual allocations. The compiled factor wrapper and worker were
also inspected without finding allocation calls; neither observation alone
certifies every possible input or future compiler version.

## Reproduce the diagnosis

Use the same environment and SBCL arguments as `scripts/check.ps1`. Run three
separate processes, setting `COALTON_LINALG_RUNTIME_PROBE` to `factor-first`,
`idle-first`, and `paused`, and use `--script tools/runtime-probe.lisp`.
The probe accepts only those modes and requires its allocating control to register
positive bytes. It is an optional diagnostic, not part of the numerical runtime.

The original run logs are the ignored files
`build/runtime-probe-factor-first.log`, `build/runtime-probe-idle-first.log`, and
`build/runtime-probe-paused.log`. Their figures above are a measured record;
reruns may distribute background allocation differently across intervals.
