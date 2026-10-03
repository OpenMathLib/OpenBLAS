#Requires -Version 7.4
# Builds and runs a CMake project against a Windows `make install` tree, with MinGW and MSVC.
param(
    [Parameter(Mandatory)][string]$Prefix,
    [Parameter(Mandatory)][ValidateSet('x64', 'arm64')][string]$Arch,
    [Parameter(Mandatory)][string]$MinGWCompiler,
    [Parameter(Mandatory)][string]$DefFile,
    [Parameter(Mandatory)][string]$BuildDir
)
$ErrorActionPreference = 'Stop'
$PSNativeCommandUseErrorActionPreference = $true

$Prefix = $Prefix -replace '\\', '/'
$env:PATH = "$($Prefix -replace '/', '\')\bin;$env:PATH"

function Test-Consumer([string]$name, [string]$implib, [string[]]$cmakeArgs) {
    $build = Join-Path $BuildDir $name
    cmake -S $PSScriptRoot -B $build "-DCMAKE_PREFIX_PATH=$Prefix" "-DEXPECTED_IMPLIB=$implib" @cmakeArgs
    cmake --build $build
    & "$build/consumer.exe"
}

Test-Consumer mingw libopenblas.dll.a @('-G', 'MinGW Makefiles', "-DCMAKE_C_COMPILER=$MinGWCompiler", '-DCMAKE_SH=CMAKE_SH-NOTFOUND')

$vcTools = @{ x64 = 'x86.x64'; arm64 = 'ARM64' }[$Arch]
$vs = & "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe" -latest -products * `
    -requires "Microsoft.VisualStudio.Component.VC.Tools.$vcTools" -property installationPath
$hostArch = if ($env:PROCESSOR_ARCHITECTURE -eq 'ARM64') { 'arm64' } else { 'amd64' }
& "$vs\Common7\Tools\Launch-VsDevShell.ps1" -Arch ($Arch -eq 'x64' ? 'amd64' : 'arm64') -HostArch $hostArch -SkipAutomaticLocation
Test-Consumer msvc-dlla libopenblas.dll.a @('-G', 'NMake Makefiles', '-DCMAKE_C_COMPILER=cl')

# release packages add an MSVC import library, which the config should prefer
lib /nologo "/machine:$Arch" "/def:$DefFile" /name:libopenblas.dll "/out:$Prefix/lib/libopenblas.lib"
Test-Consumer msvc-lib libopenblas.lib @('-G', 'NMake Makefiles', '-DCMAKE_C_COMPILER=cl')
