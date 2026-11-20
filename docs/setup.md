# Standalone setup on Windows

The supported validation host is Windows x86-64 with SBCL 2.6.10 and Python
3.11 or newer. Python needs no third-party packages. The library itself uses
ASDF and Coalton; Python, PowerShell, and Quicklisp registration are check tooling.
Other operating systems and Lisp implementations are not covered by the recorded
results.

Two files describe different parts of reproducibility:

- [toolchain.sources.json](../toolchain.sources.json) gives exact acquisition
  coordinates, archive hashes where established, and installation destinations.
- [toolchain.lock.json](../toolchain.lock.json) pins installed file contents,
  compiler settings, and actual loaded ASDF systems. It is the acceptance gate
  after acquiring or copying dependencies.

A source archive hash does not establish a successful build, and a version string
does not establish identical installed content. The current lock includes the
tested Windows SBCL executable, core, and contrib files. An independently rebuilt
SBCL with the same version may legitimately fail that lock.

## Use existing dependencies

Copy [config.local.example.psd1](../config.local.example.psd1) to the ignored
`config.local.psd1`. Set `Sbcl`, `Python`, and absolute `CoaltonRoot` and
`QuicklispRoot` paths. The two source roots must match the installed-content lock.
`CacheRoot` may be omitted to use `build/cache`.

No private setup file or user Lisp init file is needed. The repository's
[register-dependencies.lisp](../scripts/register-dependencies.lisp) registers
dependencies. It expects a complete, isolated Quicklisp root without
`asdf-config/init.lisp` or `local-init/*.lisp` / `*.cl` startup hooks. Those hooks
are refused before they can run.

Run the checks serially, from this repository:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\check.ps1
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\check.ps1 -Safety 0
```

Use a new, nonexistent cache directory to rebuild dependency sources:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\check.ps1 -CacheRoot .\build\review-cache
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\check.ps1 -Safety 0 -CacheRoot .\build\review-cache
```

The first command fills the cache. The second reuses those dependencies and
recompiles the library at safety 0. Do not describe an existing cache as fresh.
The full command checks installed contents before loading, and checks actual
ASDF source resolution and compiler mode after the tests.
[Toolchain verification](toolchain.md) explains these checks and their limits.

## Acquire the pinned sources

These commands deliberately acquire dependencies. The normal check runner
does not download or install anything. Run the recipe in a checkout of this
repository, using an unused directory beneath `build/`.

Install the Windows x86-64 SBCL 2.6.10 binary referenced by the
[official platform table](https://www.sbcl.org/platform-table.html), or use an
already installed copy that passes the lock. The acquisition manifest records
its versioned MSI URL, but does not yet establish the MSI's own SHA256 from a
primary source or local archive. Verify the installed files with the lock before
using them; do not substitute a version-only check.

The Coalton checkout must reproduce the recorded line endings. Its upstream
attributes require LF for Lisp/Coalton sources, while the recorded Windows Git
checkout uses CRLF for other text files. A GitHub source archive without that
checkout transformation will not match the installed-content lock.

```powershell
$ErrorActionPreference = 'Stop'
$recipe = Get-Content -LiteralPath .\toolchain.sources.json -Raw | ConvertFrom-Json
$dependencyRoot = Join-Path (Get-Location).Path 'build\deps-review'
if (Test-Path -LiteralPath $dependencyRoot) {
    throw 'Choose a new dependency directory for this acquisition.'
}
$archiveRoot = Join-Path $dependencyRoot 'archives'
New-Item -ItemType Directory -Path $archiveRoot | Out-Null
$coaltonRoot = Join-Path $dependencyRoot 'coalton'
$quicklispRoot = Join-Path $dependencyRoot 'quicklisp'

git clone --no-checkout --config core.autocrlf=true --config core.eol=crlf $recipe.coalton.repository $coaltonRoot
if ($LASTEXITCODE -ne 0) { throw 'Coalton clone failed.' }
git -C $coaltonRoot checkout --detach $recipe.coalton.revision
if ($LASTEXITCODE -ne 0) { throw 'Pinned Coalton checkout failed.' }
```

Fetch the 20 versioned Quicklisp files listed in the manifest. All are hashed
before any archive is extracted or Lisp file is loaded. The 13 dependency
archive SHA256s were measured from the locally cached release archives; their
URLs, lengths, and legacy hashes came from the pinned release index. The client
tar SHA256 comes from the installed canonical client metadata. Every entry
records this distinction. Hashes identify the reviewed bytes; they are not
claims of upstream signatures.

```powershell
foreach ($item in $recipe.downloads) {
    $archivePath = Join-Path $archiveRoot $item.file
    Invoke-WebRequest -Uri $item.url -OutFile $archivePath
    $actualHash = (Get-FileHash -LiteralPath $archivePath -Algorithm SHA256).Hash
    if ((Get-Item -LiteralPath $archivePath).Length -ne $item.bytes -or
        $actualHash -ine $item.sha256) {
        throw "Downloaded content differs from the recorded artifact: $($item.id)"
    }
}
```

Keep the source archives' bytes unchanged when extracting. The paths in the
manifest are relative to `$dependencyRoot`. Windows' `tar` extracts the listed
tar/tgz archives without running their contents.

```powershell
foreach ($item in $recipe.downloads) {
    $archivePath = Join-Path $archiveRoot $item.file
    if ($item.PSObject.Properties['copy_to']) {
        $destinationPath = Join-Path $dependencyRoot $item.copy_to
        New-Item -ItemType Directory -Force -Path (Split-Path -Parent $destinationPath) | Out-Null
        Copy-Item -LiteralPath $archivePath -Destination $destinationPath
    } elseif ($item.PSObject.Properties['unpack_to']) {
        $destinationPath = Join-Path $dependencyRoot $item.unpack_to
        New-Item -ItemType Directory -Force -Path $destinationPath | Out-Null
        tar -xf $archivePath -C $destinationPath
        if ($LASTEXITCODE -ne 0) { throw "Extraction failed: $($item.id)" }
    }
}
```

The resulting layout is:

```text
build/deps-review/
  coalton/                       detached pinned Git checkout
  archives/                      verified downloaded archives and metadata
  quicklisp/
    setup.lisp
    client-info.sexp
    asdf.lisp                    retained; the runner uses SBCL's newer ASDF
    quicklisp/                   extracted client source
    dists/quicklisp/
      distinfo.txt
      releases.txt
      systems.txt
      software/
        acclimation-20241012-git/
        ...                      thirteen pinned release trees
```

No Quicklisp installed-system index, local-project registration, or dist update
is required: the repository registers the extracted ASDF sources directly.
Quicklisp's client is still loaded to preserve the currently reviewed runtime.

Check the complete installed content before compiling. Replace the SBCL path
if it was installed elsewhere:

```powershell
$sbclRoot = 'C:\Program Files\Steel Bank Common Lisp'
python .\tools\check-toolchain.py --sbcl-root $sbclRoot --coalton-root $coaltonRoot --quicklisp-root $quicklispRoot
if ($LASTEXITCODE -ne 0) { throw 'Installed toolchain content does not match.' }

powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\check.ps1 -Sbcl (Join-Path $sbclRoot 'sbcl.exe') -Python python -CoaltonRoot $coaltonRoot -QuicklispRoot $quicklispRoot -CacheRoot .\build\acquired-source-cache
if ($LASTEXITCODE -ne 0) { throw 'Fresh source validation failed.' }
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\check.ps1 -Safety 0 -Sbcl (Join-Path $sbclRoot 'sbcl.exe') -Python python -CoaltonRoot $coaltonRoot -QuicklispRoot $quicklispRoot -CacheRoot .\build\acquired-source-cache
if ($LASTEXITCODE -ne 0) { throw 'Safety-0 validation failed.' }
```

If network access is unavailable, copying the already pinned source trees and
bootstrap files into a new isolated directory can exercise registration,
compilation, and tests offline. That proves independence from the old source
paths and cache. It does not prove that the remote acquisition coordinates
currently deliver the pinned bytes. Keep those two results distinct.

## What the registration script controls

The PowerShell runner starts a new SBCL process with system and user init files
disabled. It supplies `COALTON_LINALG_COALTON_ROOT`,
`COALTON_LINALG_QUICKLISP_ROOT`, and the chosen cache as `XDG_CACHE_HOME`, and sets
`COALTON_ENV=release`. Registration enables Coalton specialization and type
annotations before compiling Coalton.

The script uses ASDF's documented
[source-registry configuration](https://asdf.common-lisp.dev/asdf/Configurations.html)
and [output translations](https://asdf.common-lisp.dev/asdf/Output-Configuration-DSL.html)
to ignore inherited source/cache settings. Source discovery uses this checkout,
the explicit Coalton tree, the named Quicklisp client and software trees, and
SBCL's implementation-provided contrib systems. Compiled outputs go to the chosen
cache, separated by implementation and by named `coalton`, `quicklisp`, and
`library` directories. Each mapping preserves the path relative to that source
root, without repeating its full absolute directory prefix. ASDF chooses the
longest matching source prefix, so a dependency nested inside the checkout still
uses its own namespace. A synthetic nested output probe checks all three mappings
before Quicklisp loads. The fallback for other sources remains separate.

Use a separate cache directory for each checkout and compiler mode. These short
namespaces do not include the original checkout path as a cache key.

The initial review exposed a limitation of the previous catch-all cache layout.
A fresh full run exited with SBCL access violation `-1073741819` after dependency
registration; completed compiled-file paths reached 288 characters. With the
same acquired sources and checkout, changing only to a fresh shorter cache path
completed the safety-1 checks and produced 338 FASLs; completed paths there
reached 268 characters. This does not establish a precise path-length threshold
or prove the crash's underlying cause. The explicit root mappings remove the
observed repetition of absolute source prefixes; unusually long cache roots
or source-relative names can still produce long paths.

The mappings shorten cache paths, not source paths. A consumer on this host and SBCL
version reported, on 2026-11-20, the same access violation with no message when the
Quicklisp root itself was 122 characters deep and none at 79, and that ASDF's recursive
scan of Quicklisp's software folder then reaches files beyond 260 characters. That was not
reproduced for this record. Keep the dependency roots short as well as the cache root.

After loading the named Quicklisp setup, the script restores ASDF's previous
search functions and clears Quicklisp local-project directories. It does not
call `ql:quickload`, register a neighboring checkout, or fetch missing systems.
Missing dependencies fail. Quicklisp may create its empty `local-projects`
directory during setup; source files are not changed by this registration.

The ordinary `coalton-linalg` ASDF system does not load the registration script.
Consumers may register Coalton and the library through normal ASDF facilities.
Use a fresh release-mode Coalton process for behavior comparable to the recorded
checks. Ordinary ASDF use is not a claim that another host or compiler has been
validated.