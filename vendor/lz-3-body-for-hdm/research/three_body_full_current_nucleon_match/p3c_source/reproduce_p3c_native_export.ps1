[CmdletBinding()]
param(
    [string]$OutputDirectory = ""
)

$ErrorActionPreference = "Stop"

function Get-RepoRoot {
    $root = & git -C $PSScriptRoot rev-parse --show-toplevel 2>$null
    if ($LASTEXITCODE -ne 0 -or [string]::IsNullOrWhiteSpace(($root -join ""))) {
        throw "Run this wrapper from a Git checkout; unable to resolve repository root."
    }
    return ($root -join "").Trim()
}

$repoRoot = Get-RepoRoot
$relativeSource = "v6/dmformfactor-V6.m"
$packagePath = Join-Path $repoRoot ($relativeSource -replace "/", [IO.Path]::DirectorySeparatorChar)
$expectedPackageSha256 = "dc0742f1c50b0a478adab62f2ddef6c82978ee034c1f177e0ed917deab5dd372"
if (-not (Test-Path -LiteralPath $packagePath -PathType Leaf)) {
    throw "Native AFH source is missing: $relativeSource"
}

$actualPackageSha256 = (Get-FileHash -LiteralPath $packagePath -Algorithm SHA256).Hash.ToLowerInvariant()
if ($actualPackageSha256 -ne $expectedPackageSha256) {
    throw "AFH source SHA-256 mismatch. Expected $expectedPackageSha256, got $actualPackageSha256"
}

$sourceDirectory = $PSScriptRoot
$wolframScript = Join-Path $sourceDirectory "export_p3c_afh_multipoles_native.wls"
$expectedWolframScriptSha256 = "dc2ae4ba6a7bdb03073e4c1f9ce498317112b29af6d890d915be745e610cc85a"
$actualWolframScriptSha256 = (Get-FileHash -LiteralPath $wolframScript -Algorithm SHA256).Hash.ToLowerInvariant()
if ($actualWolframScriptSha256 -ne $expectedWolframScriptSha256) {
    throw "Archived Wolfram exporter changed. Expected $expectedWolframScriptSha256, got $actualWolframScriptSha256"
}

if ([string]::IsNullOrWhiteSpace($OutputDirectory)) {
    $stamp = [DateTime]::UtcNow.ToString("yyyyMMddTHHmmssZ")
    $OutputDirectory = Join-Path $repoRoot "research/three_body_full_current_nucleon_match/results/2026-11-22/p3c_bound_current_amplitude/reproduction_reruns/$stamp"
}
if (-not [IO.Path]::IsPathRooted($OutputDirectory)) {
    $OutputDirectory = Join-Path $repoRoot $OutputDirectory
}
$OutputDirectory = [IO.Path]::GetFullPath($OutputDirectory)
if (Test-Path -LiteralPath $OutputDirectory) {
    throw "Refusing to overwrite an existing reproduction directory: $OutputDirectory"
}
$nativeOutputDirectory = Join-Path $OutputDirectory "native"
New-Item -ItemType Directory -Path $nativeOutputDirectory -Force | Out-Null

$previousAfhPath = $env:DMFormFactorPath
$previousOutputPath = $env:P3C_NATIVE_OUTDIR
try {
    # The archived .wls file's relative fallback walks up only four levels.
    # Set the package path explicitly so that fallback is never used.
    $env:DMFormFactorPath = $packagePath
    $env:P3C_NATIVE_OUTDIR = $nativeOutputDirectory

    $wolfram = Get-Command "wolframscript" -ErrorAction Stop
    # Use -file: wolframscript's -script argument parser splits paths at spaces.
    & $wolfram.Source -file $wolframScript
    if ($LASTEXITCODE -ne 0) {
        throw "wolframscript exited with code $LASTEXITCODE"
    }
}
finally {
    $env:DMFormFactorPath = $previousAfhPath
    $env:P3C_NATIVE_OUTDIR = $previousOutputPath
}

$expectedOutputs = @(
    "afh_v6_p3c_xe129_xe131_native_multipoles.csv",
    "afh_v6_p3c_xe129_xe131_native_obdme_terms.csv",
    "native_metadata.txt"
)
foreach ($name in $expectedOutputs) {
    $path = Join-Path $nativeOutputDirectory $name
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "Expected Native AFH output was not produced: $name"
    }
}

$os = Get-CimInstance Win32_OperatingSystem
$cpu = Get-CimInstance Win32_Processor | Select-Object -First 1
$gitCommit = (& git -C $repoRoot rev-parse HEAD).Trim()
$gitBranch = (& git -C $repoRoot rev-parse --abbrev-ref HEAD).Trim()
$wolframMetadata = Get-Content -LiteralPath (Join-Path $nativeOutputDirectory "native_metadata.txt") -Raw
$outputRecords = Get-ChildItem -LiteralPath $nativeOutputDirectory -File | ForEach-Object {
    [ordered]@{
        file = $_.Name
        bytes = $_.Length
        sha256 = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
    }
}

$provenance = [ordered]@{
    schema_version = "1.0"
    record = "p3c_native_afh_reproduction"
    run_utc = [DateTime]::UtcNow.ToString("o")
    repository = [ordered]@{
        root_alias = "."
        branch = $gitBranch
        commit = $gitCommit
    }
    source = [ordered]@{
        afh_package_path = $relativeSource
        afh_package_sha256 = $actualPackageSha256
        wolfram_exporter_path = "research/three_body_full_current_nucleon_match/results/2026-11-22/p3c_bound_current_amplitude/source/export_p3c_afh_multipoles_native.wls"
        wolfram_exporter_sha256 = $actualWolframScriptSha256
    }
    runtime = [ordered]@{
        os_caption = $os.Caption
        os_version = $os.Version
        powershell_version = $PSVersionTable.PSVersion.ToString()
        cpu_name = $cpu.Name
        logical_processors = $cpu.NumberOfLogicalProcessors
        physical_memory_bytes = [int64]$os.TotalVisibleMemorySize * 1024
        wolfram_executable = $wolfram.Source
        wolfram_native_metadata = $wolframMetadata
    }
    outputs = @($outputRecords)
    note = "This is a new reproduction run. The historical run_manifest.json and its hashed outputs are not modified."
}

$provenancePath = Join-Path $OutputDirectory "reproduction_provenance.json"
$provenance | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $provenancePath -Encoding UTF8
Write-Host "P3-C Native AFH reproduction completed."
Write-Host "Native output: $nativeOutputDirectory"
Write-Host "Sanitized provenance: $provenancePath"
