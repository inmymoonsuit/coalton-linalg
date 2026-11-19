param(
    [ValidateSet('Test','Fixtures','Inspect')][string]$Task = 'Test',
    [ValidateSet(0,1)][int]$Safety = 1,
    [string]$Sbcl = '',
    [string]$Python = '',
    [string]$CacheRoot = '',
    [string]$CoaltonRoot = '',
    [string]$QuicklispRoot = '',
    [string]$InputFolder = '',
    [string]$OutputFolder = '',
    [ValidateRange(1,1000000)][int]$Repeats = 20,
    [ValidateRange(1,1000000)][int]$AllocationRepeats = 100
)
$ErrorActionPreference = 'Stop'
$libraryRoot = Split-Path -Parent $PSScriptRoot
$configPath = Join-Path $libraryRoot 'config.local.psd1'
$config = @{}
if (Test-Path -LiteralPath $configPath) { $config = Import-PowerShellDataFile -LiteralPath $configPath }
foreach ($setting in @('Sbcl','Python','CacheRoot','CoaltonRoot','QuicklispRoot')) {
    if (-not $PSBoundParameters.ContainsKey($setting) -and $config.ContainsKey($setting)) {
        Set-Variable -Name $setting -Value $config[$setting]
    }
}
if (-not $Sbcl) { $Sbcl = 'sbcl' }
if (-not $Python) { $Python = 'python' }
$Sbcl = (Get-Command $Sbcl -CommandType Application -ErrorAction Stop).Source
if ($Task -eq 'Test') { $Python = (Get-Command $Python -CommandType Application -ErrorAction Stop).Source }
if (-not $CoaltonRoot -or -not $QuicklispRoot) {
    throw 'Set CoaltonRoot and QuicklispRoot to the pinned installed dependencies; see docs/setup.md.'
}
$CoaltonRoot = (Resolve-Path -LiteralPath $CoaltonRoot -ErrorAction Stop).Path
$QuicklispRoot = (Resolve-Path -LiteralPath $QuicklispRoot -ErrorAction Stop).Path
if (-not $CacheRoot) { $CacheRoot = Join-Path $libraryRoot 'build/cache' }
$CacheRoot = [System.IO.Path]::GetFullPath($CacheRoot)
if (-not $InputFolder) { $InputFolder = Join-Path $libraryRoot 'tests/fixtures' }
if (-not $OutputFolder) { $OutputFolder = Join-Path $libraryRoot "build/safety$Safety/fixtures" }
$InputFolder = [System.IO.Path]::GetFullPath($InputFolder)
$OutputFolder = [System.IO.Path]::GetFullPath($OutputFolder)
$buildRoot = [System.IO.Path]::GetFullPath((Join-Path $libraryRoot 'build'))
$buildPrefix = $buildRoot.TrimEnd('\','/') + [System.IO.Path]::DirectorySeparatorChar
if (-not $OutputFolder.StartsWith($buildPrefix, [System.StringComparison]::OrdinalIgnoreCase)) {
    throw 'OutputFolder must be a subdirectory of this library build directory.'
}
if ($Task -eq 'Fixtures' -and
    -not (Get-ChildItem -LiteralPath $InputFolder -Filter '*.A.txt' -File -ErrorAction Stop)) {
    throw "No matrix fixtures in $InputFolder"
}
[System.IO.Directory]::CreateDirectory($OutputFolder) | Out-Null
[System.IO.Directory]::CreateDirectory($CacheRoot) | Out-Null
$environment = @{
    XDG_CACHE_HOME = $CacheRoot.Replace('\','/')
    COALTON_ENV = 'release'
    COALTON_LINALG_COALTON_ROOT = $CoaltonRoot.Replace('\','/')
    COALTON_LINALG_QUICKLISP_ROOT = $QuicklispRoot.Replace('\','/')
    COALTON_LINALG_SAFETY = "$Safety"
    COALTON_LINALG_INPUT = $InputFolder.Replace('\','/')
    COALTON_LINALG_OUTPUT = $OutputFolder.Replace('\','/') + '/'
    COALTON_LINALG_METADATA = (Join-Path $buildRoot 'toolchain-metadata.tsv').Replace('\','/')
    COALTON_LINALG_STRESS_INPUT = (Join-Path $buildRoot 'stress-input').Replace('\','/') + '/'
    COALTON_LINALG_STRESS_OUTPUT = (Join-Path $buildRoot "safety$Safety/stress").Replace('\','/') + '/'
    COALTON_LINALG_REPEATS = "$Repeats"
    COALTON_LINALG_ALLOCATION_REPEATS = "$AllocationRepeats"
}
$previous = @{}
foreach ($name in $environment.Keys) {
    $previous[$name] = [Environment]::GetEnvironmentVariable($name, 'Process')
    [Environment]::SetEnvironmentVariable($name, $environment[$name], 'Process')
}
try {
    if ($Task -eq 'Test') {
        $toolchainCheck = Join-Path $libraryRoot 'tools/check-toolchain.py'
        & $Python $toolchainCheck --sbcl-root (Split-Path -Parent $Sbcl) --coalton-root $CoaltonRoot --quicklisp-root $QuicklispRoot
        if ($LASTEXITCODE -ne 0) { throw 'Installed toolchain differs from toolchain.lock.json.' }
        $stress = Join-Path $libraryRoot 'tools/stress.py'
        & $Python $stress generate $environment.COALTON_LINALG_STRESS_INPUT
        if ($LASTEXITCODE -ne 0) { throw 'Stress fixture generation failed.' }
        [System.IO.Directory]::CreateDirectory($environment.COALTON_LINALG_STRESS_OUTPUT) | Out-Null
    }
    $entry = switch ($Task) {
        Test { 'check.lisp' }
        Fixtures { 'fixtures.lisp' }
        Inspect { 'inspect.lisp' }
    }
    & $Sbcl --noinform --no-sysinit --no-userinit --disable-debugger --script (Join-Path $PSScriptRoot $entry)
    if ($LASTEXITCODE -ne 0) { throw "Lisp $Task failed with exit code $LASTEXITCODE." }
    if ($Task -eq 'Test') {
        $verifier = Join-Path $libraryRoot 'tools/verify.py'
        $selftest = Join-Path $buildRoot "safety$Safety/verifier-selftest"
        & $Python $verifier selftest $selftest
        if ($LASTEXITCODE -ne 0) { throw 'Independent verifier self-test failed.' }
        & $Python $verifier verify (Join-Path $libraryRoot 'tests/fixtures') $OutputFolder --require-shifts --json (Join-Path $OutputFolder 'verification.json')
        if ($LASTEXITCODE -ne 0) { throw 'Independent factor verification failed.' }
        $generated = Join-Path $selftest 'fixtures'
        $expectedFiles = @(Get-ChildItem -LiteralPath (Join-Path $libraryRoot 'tests/fixtures') -File | Sort-Object Name)
        $actualFiles = @(Get-ChildItem -LiteralPath $generated -File | Sort-Object Name)
        if (Compare-Object @($expectedFiles.Name) @($actualFiles.Name)) {
            throw 'Generated fixture roster differs from checked-in fixtures.'
        }
        foreach ($fixture in $expectedFiles) {
            $expectedText = [IO.File]::ReadAllText($fixture.FullName).Replace("`r`n", "`n")
            $actualText = [IO.File]::ReadAllText((Join-Path $generated $fixture.Name)).Replace("`r`n", "`n")
            if ($expectedText -cne $actualText) { throw "Fixture is not reproducible: $($fixture.Name)" }
        }
        Write-Output "FIXTURE-REPRODUCTION-PASS $($expectedFiles.Count) files"
        & $Python $toolchainCheck --sbcl-root (Split-Path -Parent $Sbcl) --coalton-root $CoaltonRoot --quicklisp-root $QuicklispRoot --metadata $environment.COALTON_LINALG_METADATA
        if ($LASTEXITCODE -ne 0) { throw 'Loaded toolchain differs from the pin.' }
        & $Python $toolchainCheck --selftest (Join-Path $buildRoot 'toolchain-selftest')
        if ($LASTEXITCODE -ne 0) { throw 'Toolchain checker control failed.' }
        & $Python $stress selftest (Join-Path $buildRoot 'stress-selftest')
        if ($LASTEXITCODE -ne 0) { throw 'Stress verifier control failed.' }
        & $Python $stress verify $environment.COALTON_LINALG_STRESS_INPUT $environment.COALTON_LINALG_STRESS_OUTPUT --json (Join-Path $environment.COALTON_LINALG_STRESS_OUTPUT 'verification.json')
        if ($LASTEXITCODE -ne 0) { throw 'Stress factor verification failed.' }
        Write-Output "CHECK-PASS safety=$Safety"
    }
}
finally {
    foreach ($name in $environment.Keys) {
        [Environment]::SetEnvironmentVariable($name, $previous[$name], 'Process')
    }
}