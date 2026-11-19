# Reproducible offline checks

`toolchain.lock.json` records 18 installed components and 1,314 content-hashed
files: SBCL, Coalton, Quicklisp bootstrap/client/distribution metadata, and installed
dependency sources. It records the actual 28 loaded dependency systems, their
ASDF source locations relative to pinned roots, and compiler settings. Absolute
desktop paths are kept in ignored local configuration, not in the lock.

This is an installed-content lock. It is not a downloader or installer, and
contains no redistributed dependency implementations. It covers the tested
Windows x86-64 toolchain; other hosts need a separately reviewed environment.

The full `scripts/check.ps1` command checks source content before loading, then
checks the actual ASDF resolutions and release/specialization/type-annotation
settings after the numerical tests. Its checker plants six mismatches to ensure
changed/missing/extra sources and runtime metadata do not pass silently.
Failures propagate as nonzero process exits.

Set `CoaltonRoot` and `QuicklispRoot` in `config.local.psd1`, alongside the
existing SBCL, Python, and CacheRoot fields. The supplied example lists
these settings without machine-specific locations. The checker requires Python
3.11 or newer; this installation uses 3.13.13.

A cache is not source. Content and resolution checks do not certify previously
compiled FASLs. To rebuild dependencies from an empty cache, choose a new,
nonexistent directory under `build/`:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\check.ps1 -CacheRoot .\build\new-source-cache
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\check.ps1 -Safety 0 -CacheRoot .\build\new-source-cache
```

The second command reuses the freshly built dependency cache and recompiles our
library at its other safety policy. Do not call an existing directory a fresh
cache. Run all jobs serially.

The standalone script preloads the Lisp systems, then pauses **its disposable
process's** SBCL finalizer thread for bounded numerical/allocation checks,
restoring it with unwind-protect, including after an intentionally planted error.
The library and ordinary ASDF systems do not change runtime threads. This isolates
the known background contribution to SBCL's program-wide allocation counter;
[the experiment and limits](runtime-allocation.md) document why. The canonical
script is the allocation protocol; arbitrary interactive REPL background work
can contaminate a raw ASDF test-system measurement.

[The setup recipe](setup.md) and [acquisition manifest](../toolchain.sources.json)
describe obtaining the pinned inputs. The check runner uses repository-owned
registration with explicit roots and ignores inherited ASDF configuration.
Hosted CI worker setup remains separate release work.
The local full-check command is suitable as the eventual CI entry point once
a supported worker has the pinned dependencies.