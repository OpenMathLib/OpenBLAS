This page describes how to natively build OpenBLAS library for windows on arm64 targets.

See [below](#xcomp) for how to cross-compile OpenBLAS for WoA on an x86_64 Windows host.

_We hope that the procedure can be simplified for/by LLVM 17+ and its flang-new compiler for the FORTRAN parts (LAPACK) of the code but we have no means of verifying that as of January 2024._


# Prerequisite

Following tools needs to be installed

## 1. Download and install clang for windows on arm

Find the latest LLVM build for WoA from [LLVM release page](https://releases.llvm.org/)

E.g: LLVM 12 build for WoA64 can be found [here](https://github.com/llvm/llvm-project/releases/download/llvmorg-12.0.0/LLVM-12.0.0-woa64.exe)

Run the LLVM installer and ensure that LLVM is added to environment PATH.

## 2. Download and install classic flang for windows on arm

Classic flang is the only available FORTRAN compiler for windows on arm for now and a pre-release build can be found [here](https://github.com/kaadam/flang/releases/tag/v0.1)

There is no installer for classic flang and the zip package can be extracted and the path needs to be added to environment PATH.

E.g: on PowerShell

`$env:Path += ";C:\flang_woa\bin"`

# Build

The following steps describe how to build the static library for OpenBLAS with and without LAPACK

## 1. Build OpenBLAS static library with BLAS and LAPACK routines with Make

Following command can be used to build OpenBLAS static library with BLAS and LAPACK routines

`make CC="clang-cl" HOSTCC="clang-cl" AR="llvm-ar" BUILD_WITHOUT_LAPACK=0 NOFORTRAN=0 DYNAMIC_ARCH=0 TARGET=ARMV8 ARCH=arm64 BINARY=64 USE_OPENMP=0 PARALLEL=1 RANLIB="llvm-ranlib" MAKE=make F_COMPILER=FLANG FC=FLANG FFLAGS_NOOPT="-march=armv8-a -cpp" FFLAGS="-march=armv8-a -cpp" NEED_PIC=0 HOSTARCH=arm64 libs netlib`

## 2. Build static library with BLAS routines using CMake

Classic flang has compatibility issues with cmake hence only BLAS routines can be compiled with CMake

`mkdir build`

`cd build`

`cmake ..  -G Ninja -DCMAKE_C_COMPILER=clang -DBUILD_WITHOUT_LAPACK=1 -DNOFORTRAN=1 -DDYNAMIC_ARCH=0 -DTARGET=ARMV8 -DARCH=arm64 -DBINARY=64 -DUSE_OPENMP=0 -DCMAKE_SYSTEM_PROCESSOR=ARM64 -DCMAKE_CROSSCOMPILING=1 -DCMAKE_SYSTEM_NAME=Windows`

`cmake --build . --config Release`

# Known issue

## `getarch.exe` execution error

If you notice that platform-specific headers by `getarch.exe` are not generated correctly, It could be due to a known debug runtime DLL issue for arm64 platforms. Please check out [link](https://linaro.atlassian.net/wiki/spaces/WOAR/pages/28677636097/Debug+run-time+DLL+issue#Workaround) for the workaround.

# <a id="xcomp"></a>Alternative: cross-compiling on Windows x86
## Prerequisites:
a working installation of LLVM and Ninja (the versions included in the latest Visual Studio 2022 will do, if you install its optional Llvm.Clang component) 

## Build
 1. Load into the appropriate cross compiling developer shell (in this case powershell, for cmd.exe it would be vcvarsamd64_arm64.bat): 
      'C:\Program Files\Microsoft Visual Studio\2022\Professional\Common7\Tools\Launch-VsDevShell.ps1' -Arch arm64 -HostArch amd64

 2. Change to the folder where you unpacked the OpenBLAS sources if you haven't already, and create a subfolder for building. 
    Then change into that folder. 

 3. Invoke the following CMake command

    ```
    cmake .. -G Ninja -DCMAKE_BUILD_TYPE=Release -DTARGET=ARMV8 -DCMAKE_CROSSCOMPILING=ON -DCMAKE_SYSTEM_NAME="Windows" -DARCH=arm     -DBINARY=64 -DCMAKE_SYSTEM_PROCESSOR=ARM64 -DCMAKE_C_COMPILER=clang-cl -DCMAKE_C_COMPILER_TARGET=arm64-pc-windows-msvc    -DCMAKE_ASM_COMPILER_TARGET=arm64-pc-windows-msvc.
    ``` 
    (Note the use of clang-cl instead of clang, and adding in the _COMPILER_TARGET entries. Without C_COMPILER_TARGET, clang-cl will be       trying to compile x64 machine code and not make it past the cmake generation step. Without ASM_COMPILER_TARGET the assembly files in openblas will fail to assemble because the assembler is trying to interpret them as x64 assembly.)

 4. Invoke Ninja for the actual compilation

