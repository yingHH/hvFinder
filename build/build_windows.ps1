# -*- powershell -*-
<#
.SYNOPSIS
    Rebuild the hvFinder desktop client (onefile exe) on Windows.

.DESCRIPTION
    - create the conda env from build/environment-win.yml when missing, else reuse it
    - build with hvFinder.spec (independent of the current directory)
    - verify dist/hvFinder.exe, print size and SHA256, write build/BUILD_INFO.txt

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File build\build_windows.ps1
#>
param(
    [string]$EnvName = "hvFinder"
)

$ErrorActionPreference = "Stop"

# repo root = parent of build/
$RepoRoot = Split-Path -Parent $PSScriptRoot
Set-Location $RepoRoot
Write-Host "== hvFinder build ==" -ForegroundColor Cyan
Write-Host "repo: $RepoRoot"

# 0) platform check: PyInstaller cannot cross-compile
if ($env:OS -ne "Windows_NT") {
    throw "This script must run on Windows; PyInstaller cannot produce a usable .exe on Linux/macOS."
}

# 1) conda executable
# Do not use the `conda` function (conda init injects a PowerShell function which
# re-invokes itself under -File and overflows the call depth); call the exe instead.
$condaExe = if ($Env:CONDA_EXE -like "*.exe") { $Env:CONDA_EXE } else { "conda.exe" }

# 2) environment: create it when missing, reuse it otherwise
$envYml = Join-Path $RepoRoot "build\environment-win.yml"
if (-not (Test-Path $envYml)) {
    throw "Missing $envYml; check the directory layout (build/ client/ server/)"
}

# read the env name from --json: the text output contains paths and a `*` marker
$envJson = & $condaExe env list --json
if ($LASTEXITCODE -ne 0) {
    throw "conda env list --json failed (exit $LASTEXITCODE); check that $condaExe works"
}
try {
    $envList = $envJson | Out-String | ConvertFrom-Json
} catch {
    throw "Cannot parse conda env list --json output: $($_.Exception.Message)"
}
$envNames = @($envList.envs | ForEach-Object { Split-Path -Leaf $_ })

if ($envNames -contains $EnvName) {
    Write-Host "Reusing existing conda environment $EnvName" -ForegroundColor Yellow
} else {
    Write-Host "conda environment $EnvName not found, creating it from build\environment-win.yml ..." -ForegroundColor Yellow
    & $condaExe env create -f $envYml -n $EnvName
    if ($LASTEXITCODE -ne 0) {
        throw "conda env create failed (exit $LASTEXITCODE)"
    }
}

# 3) dependency check (prints a ready-to-run fix command when something is missing)
& $condaExe run --no-capture-output -n $EnvName python -c "import PySide6, paramiko, yaml, cryptography, PyInstaller; print('deps OK')"
if ($LASTEXITCODE -ne 0) {
    throw "Environment $EnvName is missing build dependencies. Fix: conda env update -f `"$envYml`" --prune"
}

if (-not (Test-Path (Join-Path $RepoRoot "client\gui.py"))) {
    throw "Missing client\gui.py; check the directory layout (client/ server/)"
}

# 4) remove previous artifacts
if (Test-Path dist)  { Remove-Item -Recurse -Force dist }
if (Test-Path build\hvFinder) { Remove-Item -Recurse -Force build\hvFinder }

# 5) build
Write-Host "PyInstaller ..." -ForegroundColor Yellow
& $condaExe run --no-capture-output -n $EnvName python -m PyInstaller --clean --noconfirm hvFinder.spec
if ($LASTEXITCODE -ne 0) {
    throw "PyInstaller build failed (exit $LASTEXITCODE)"
}

$exe = Join-Path $RepoRoot "dist\hvFinder.exe"
if (-not (Test-Path $exe)) { throw "Build failed: $exe was not produced" }

# put settings.ini next to the exe so dist/ is self-contained
$settingsSrc = Join-Path $RepoRoot "client\settings.ini"
if (Test-Path $settingsSrc) {
    Copy-Item $settingsSrc (Join-Path $RepoRoot "dist\settings.ini") -Force
    Write-Host "Copied settings.ini into dist\" -ForegroundColor Green
} else {
    Write-Warning "client\settings.ini is missing; a default one will be created next to the exe on first start."
}

# 6) artifact details
# read the version with python -c: `python -V` writes to stderr on some versions
function Get-EnvVersion {
    param([string]$Module, [string]$Fallback = "unknown")
    if (-not $Module) {
        $out = & $condaExe run -n $EnvName python -c "import sys; print(sys.version.split()[0])"
    } else {
        $out = & $condaExe run -n $EnvName python -c "import $Module as m; print(m.__version__)"
    }
    if ($LASTEXITCODE -ne 0 -or -not $out) { return $Fallback }
    return ($out | Out-String).Trim()
}

$sizeMB = [math]::Round((Get-Item $exe).Length / 1MB, 1)
$sha    = (Get-FileHash $exe -Algorithm SHA256).Hash
$ver    = Get-EnvVersion -Module ""
$pyi    = Get-EnvVersion "PyInstaller"
$pyside = Get-EnvVersion "PySide6"

$info = @(
    "build_time   : $(Get-Date -Format o)",
    "host         : $env:COMPUTERNAME ($env:OS)",
    "python       : $ver",
    "pyinstaller  : $pyi",
    "pyside6      : $pyside",
    "exe          : dist/hvFinder.exe",
    "size_mb      : $sizeMB",
    "sha256       : $sha",
    "spec         : hvFinder.spec",
    "note         : settings.ini is read from the current working directory"
)
$info | Set-Content -Path "build\BUILD_INFO.txt" -Encoding UTF8

Write-Host ""
Write-Host "Build finished" -ForegroundColor Green
$info | ForEach-Object { Write-Host "  $_" }
Write-Host ""
Write-Host "Manual smoke check:" -ForegroundColor Yellow
Write-Host "  1) double-click dist\hvFinder.exe"
Write-Host "  2) the three tabs switch; Tab 1 connects over SSH"
Write-Host "  3) after closing, check that settings.ini holds the connection settings"
