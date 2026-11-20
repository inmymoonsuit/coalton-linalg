param(
    [Parameter(Mandatory=$true)][string]$InputFolder,
    [Parameter(Mandatory=$true)][string]$OutputFolder,
    [Parameter(Mandatory=$true)][string[]]$Case,
    [ValidateSet(0,1)][int]$Safety = 1,
    [ValidateRange(10,3600)][int]$TimeoutSeconds = 180
)
# Whole-process Windows peak working set. Child handshakes keep the process alive
# at both samples; the final peak includes compilation, parsing and output.
$ErrorActionPreference = 'Stop'
$libraryRoot = Split-Path -Parent $PSScriptRoot
$config = Import-PowerShellDataFile (Join-Path $libraryRoot 'config.local.psd1')
$inputPath = [IO.Path]::GetFullPath($InputFolder)
$outputPath = [IO.Path]::GetFullPath($OutputFolder)
$buildPrefix = [IO.Path]::GetFullPath((Join-Path $libraryRoot 'build')).TrimEnd('\','/') + [IO.Path]::DirectorySeparatorChar
if (-not $outputPath.StartsWith($buildPrefix, [StringComparison]::OrdinalIgnoreCase)) {
    throw 'OutputFolder must be inside this library build directory.'
}
[IO.Directory]::CreateDirectory($outputPath) | Out-Null
$reportPath = Join-Path $outputPath 'memory.json'
if (Test-Path -LiteralPath $reportPath) { Remove-Item -LiteralPath $reportPath }
$cache = $config.CacheRoot
if (-not $cache) { $cache = Join-Path $libraryRoot 'build/cache' }
$cache = [IO.Path]::GetFullPath($cache)
[IO.Directory]::CreateDirectory($cache) | Out-Null
$records = @()
foreach ($name in $Case) {
    if ($name -cnotmatch '^[a-z0-9_-]+$') { throw 'Invalid generated case name.' }
    $pattern = Join-Path $inputPath ($name + '.pattern.txt')
    if (-not (Test-Path -LiteralPath $pattern -PathType Leaf)) { throw "Missing $pattern" }
    $caseOutput = Join-Path $outputPath $name
    [IO.Directory]::CreateDirectory($caseOutput) | Out-Null
    $start = New-Object Diagnostics.ProcessStartInfo
    $start.FileName = $config.Sbcl
    $scriptPath = Join-Path $PSScriptRoot 'measure-analysis.lisp'
    $start.Arguments = '--noinform --no-sysinit --no-userinit --disable-debugger --script "' + $scriptPath + '"'
    $start.WorkingDirectory = $libraryRoot
    $start.UseShellExecute = $false
    $start.CreateNoWindow = $true
    $start.RedirectStandardInput = $true
    $start.RedirectStandardOutput = $true
    $start.RedirectStandardError = $true
    $start.EnvironmentVariables['COALTON_ENV'] = 'release'
    $start.EnvironmentVariables['COALTON_LINALG_COALTON_ROOT'] = $config.CoaltonRoot
    $start.EnvironmentVariables['COALTON_LINALG_QUICKLISP_ROOT'] = $config.QuicklispRoot
    $start.EnvironmentVariables['COALTON_LINALG_SAFETY'] = "$Safety"
    $start.EnvironmentVariables['XDG_CACHE_HOME'] = $cache.Replace('\','/')
    $start.EnvironmentVariables['COALTON_LINALG_ANALYSIS_CASE'] = $pattern.Replace('\','/')
    $start.EnvironmentVariables['COALTON_LINALG_ANALYSIS_OUTPUT'] = $caseOutput.Replace('\','/') + '/'
    $process = New-Object Diagnostics.Process
    $process.StartInfo = $start
    $log = New-Object Collections.Generic.List[string]
    $baseline = $null
    $peak = $null
    $started = $false
    $stderr = $null
    $failure = $null
    $watch = [Diagnostics.Stopwatch]::StartNew()
    try {
        $started = $process.Start()
        if (-not $started) { throw 'Could not start SBCL.' }
        $stderr = $process.StandardError.ReadToEndAsync()
        while ($true) {
            $pending = $process.StandardOutput.ReadLineAsync()
            while (-not $pending.Wait(1000)) {
                if ($watch.Elapsed.TotalSeconds -gt $TimeoutSeconds) { throw 'SBCL exceeded the per-case deadline.' }
            }
            if ($watch.Elapsed.TotalSeconds -gt $TimeoutSeconds) { throw 'SBCL exceeded the per-case deadline.' }
            $line = $pending.Result
            if ($null -eq $line) { break }
            $log.Add($line)
            if ($line -eq 'MEASUREMENT-BASELINE') {
                $process.Refresh()
                $baseline = $process.PeakWorkingSet64
                $process.StandardInput.WriteLine('continue')
                $process.StandardInput.Flush()
            }
            elseif ($line -eq 'MEASUREMENT-DONE') {
                $process.Refresh()
                $peak = $process.PeakWorkingSet64
                $process.StandardInput.WriteLine('exit')
                $process.StandardInput.Flush()
            }
        }
        if (-not $process.WaitForExit(5000)) { throw 'SBCL did not exit after completing output.' }
        if (-not $stderr.Wait(5000)) { throw 'SBCL error stream did not close.' }
        if ($process.ExitCode -ne 0 -or $null -eq $baseline -or $null -eq $peak) {
            throw "Analysis measurement failed for $name; see its logs."
        }
        if ($baseline -le 0 -or $peak -lt $baseline) { throw 'Invalid process peak-memory samples.' }
        $records += [ordered]@{
            case = $name; safety = $Safety
            baseline_process_peak_working_set_bytes = $baseline
            final_process_peak_working_set_bytes = $peak
        }
        Write-Output "ANALYSIS-MEMORY $name baseline_peak=$baseline final_peak=$peak"
    }
    catch {
        $failure = $_.Exception.Message
        throw
    }
    finally {
        if ($started -and -not $process.HasExited) {
            $process.Kill()
            [void]$process.WaitForExit(5000)
        }
        [IO.File]::WriteAllLines((Join-Path $caseOutput 'process.log'), $log)
        $errorText = if ($null -ne $stderr -and $stderr.Wait(5000)) { $stderr.Result } else { 'Error stream unavailable.' }
        if ($failure) { $errorText += "`nHarness failure: $failure" }
        [IO.File]::WriteAllText((Join-Path $caseOutput 'stderr.log'), $errorText)
        $process.Dispose()
    }
}
$report = [ordered]@{
    format = 1
    scope = 'Windows whole-process peak working set; includes compiler/runtime, parsing, analysis and output. Baseline peak sampled after library load and full GC. Not phase-only peak or total allocation.'
    results = $records
}
[IO.File]::WriteAllText($reportPath, ($report | ConvertTo-Json -Depth 5))
