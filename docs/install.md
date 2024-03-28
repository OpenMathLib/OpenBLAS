# Install OpenBLAS

!!! note
    Lists of precompiled packages are not comprehensive, is not meant to validate nor endorse a particular third-party build over others, and may not always lead to the newest version


## Quick install

Precompiled packages have recently become available for a number of platforms through their normal installation procedures, so for users of desktop devices at least, the instructions below are mostly relevant when you want to try the most recent development snapshot from git. See your platform's relevant "Precompiled packages" section.

The [Conda-Forge](https://github.com/conda-forge) project maintains packages for the conda package manager at <https://github.com/conda-forge/openblas-feedstock>.

## Source
Download the latest [stable version](https://github.com/xianyi/OpenBLAS/releases) from release page.

## Platforms

### Linux

Just type `make` to compile the library.

Notes:

* OpenBLAS doesn't support g77. Please use gfortran or other Fortran compilers. e.g. `make FC=gfortran`.
* When building in an emulator (KVM,QEMU etc.) make sure that the combination of CPU features exposed to
  the virtual environment matches that of an existing CPU to allow detection of the cpu model to succeed. 
  (With qemu, this can be done by passing `-cpu host` or a supported model name at invocation)


#### Precompiled packages

##### Debian/Ubuntu/Mint/Kali
 OpenBLAS package is available in default repositories and can act as default BLAS in system

Example installation commands:
```bash
$ sudo apt update
$ apt search openblas
$ sudo apt install libopenblas-dev
$ sudo update-alternatives --config libblas.so.3
```
 Alternatively, if distributor's package proves unsatisfactory, you may try latest version of OpenBLAS, [Following guide in OpenBLAS FAQ](faq.md#debianlts)
 
##### openSuSE/SLE
 Recent OpenSUSE versions include OpenBLAS in default repositories and also permit OpenBLAS to act as replacement of system-wide BLAS.

 Example installation commands:
```bash
$ sudo zypper ref
$ zypper se openblas
$ sudo zypper in openblas-devel
$ sudo update-alternatives --config libblas.so.3
``` 
Should you be using older OpenSUSE or SLE that provides no OpenBLAS, you can attach optional or experimental openSUSE repository as a new package source to acquire recent build of OpenBLAS following [instructions on openSUSE software site](https://software.opensuse.org/package/openblas)

##### Fedora/CentOS/RHEL
Fedora provides OpenBLAS in default installation repositories.

To install it try following:
```bash
$ dnf search openblas
$ dnf install openblas-devel
```
For CentOS/RHEL/Scientific Linux packages are provided via [Fedora EPEL repository](https://fedoraproject.org/wiki/EPEL)

After adding repository and repository keys installation is pretty straightforward:
```bash
$ yum search openblas
$ yum install openblas-devel
```
No alternatives mechanism is provided for BLAS, and packages in system repositories are linked against NetLib BLAS or ATLAS BLAS libraries. You may wish to re-package RPMs to use OpenBLAS instead [as described here](https://fedoraproject.org/wiki/How_to_create_an_RPM_package)

##### Mageia
Mageia offers ATLAS and NetLIB LAPACK in base repositories.
You can build your own OpenBLAS replacement, and once installed in /opt
TODO: populate /usr/lib64 /usr/include accurately to replicate netlib with update-alternatives

##### Arch/Manjaro/Antergos
```bash
$ sudo pacman -S openblas
```

### Windows

The precompiled binaries available with each release (in <https://github.com/xianyi/OpenBLAS/releases>) are
created with MinGW using an option list of 
"NUM_THREADS=64 TARGET=GENERIC DYNAMIC_ARCH=1 DYNAMIC_OLDER=1 CONSISTENT_FPCSR=1" - they should work on
any x86 or x86_64 computer. The zip archive contains the include files, static and dll libraries as well
as configuration files for getting them found via CMAKE or pkgconfig - just create a suitable folder for
your OpenBLAS installation and unzip it there. (Note that you will need to edit the provided openblas.pc
and OpenBLASConfig.cmake to reflect the installation path on your computer, as distributed they have "win"
or "win64" reflecting the local paths on the system they were built on). Some programs will expect the DLL
name to be lapack.dll, blas.dll, or (in the case of the statistics package "R") even Rblas.dll to act as a
direct replacement for whatever other implementation of BLAS and LAPACK they use by default. Just copy the
openblas.dll to the desired name(s).
Note that the provided binaries are built with INTERFACE64=0, meaning they use standard 32bit integers for
array indexing and the like (as is the default for most if not all BLAS and LAPACK implementations). If the
documentation of whatever program you are using with OpenBLAS mentions 64bit integers (INTERFACE64=1) for
addressing huge matrix sizes, you will need to build OpenBLAS from source (or open an issue ticket to make
the demand for such a precompiled build known).

#### Precompiled packages

* <http://sourceforge.net/projects/openblas/files>
* <https://www.nuget.org/packages?q=openblas>

#### Visual Studio

As of OpenBLAS v0.2.15, we support MinGW and Visual Studio (using CMake to generate visual studio solution files &ndash; note that you will need at least version 3.11 of CMake for linking to work correctly) to build OpenBLAS on Windows.

Note that you need a Fortran compiler if you plan to build and use the LAPACK functions included with OpenBLAS. The sections below describe using either `flang` as an add-on to clang/LLVM or `gfortran` as part of MinGW for this purpose. If you want to use the Intel Fortran compiler `ifort` for this, be sure to also use the Intel C compiler `icc` for building the C parts, as the ABI imposed by `ifort` is incompatible with `msvc`.

##### 1. Native (MSVC) ABI

A fully-optimized OpenBLAS that can be statically or dynamically linked to your application can currently be built for the 64-bit architecture with the LLVM compiler infrastructure. We're going to use Miniconda3 to grab all of the tools we need, since some of them are in an experimental status. Before you begin, you'll need to have Microsoft Visual Studio 2015 or newer installed.

1. Install Miniconda3 for 64 bits using `winget install --id Anaconda.Miniconda3` or easily download from [conda.io](https://docs.conda.io/en/latest/miniconda.html).
2. Open the "Anaconda Command Prompt," now available in the Start Menu, or at `%USERPROFILE%\miniconda3\shell\condabin\conda-hook.ps1`.
3. In that command prompt window, use `cd` to change to the directory where you want to build OpenBLAS
4. Now install all of the tools we need:

   ```
   conda update -n base conda
   conda config --add channels conda-forge
   conda install -y cmake flang clangdev perl libflang ninja
   ```

5. Still in the Anaconda Command Prompt window, activate the MSVC environment for 64 bits with `vcvarsall x64`. On Windows 11 with Visual Studio 2022, this would be done by invoking:

   ```shell
   "c:\Program Files\Microsoft Visual Studio\2022\Preview\vc\Auxiliary\Build\vcvars64.bat"
   ```

   With VS2019, the command should be the same &ndash; except for the year number, obviously. For other/older versions of MSVC,
   the VS documentation or a quick search on the web should turn up the exact wording you need.

   Confirm that the environment is active by typing `link` &ndash; this should return a long list of possible options for the `link` command. If it just 
   returns "command not found" or similar, review and retype the call to vcvars64.bat.
   **NOTE:** if you are working from a Visual Studio Command prompt window instead (so that you do not have to do the vcvars call), you need to invoke 
   `conda activate` so that CONDA_PREFIX etc. get set up correctly before proceeding to step 6. Failing to do so will lead to link errors like 
   libflangmain.lib not getting found later in the build.

6. Now configure the project with CMake. Starting in the project directory, execute the following:

   ```
   set "LIB=%CONDA_PREFIX%\Library\lib;%LIB%"
   set "CPATH=%CONDA_PREFIX%\Library\include;%CPATH%"
   mkdir build
   cd build
   cmake .. -G "Ninja" -DCMAKE_CXX_COMPILER=clang-cl -DCMAKE_C_COMPILER=clang-cl -DCMAKE_Fortran_COMPILER=flang -DCMAKE_MT=mt -DBUILD_WITHOUT_LAPACK=no -DNOFORTRAN=0 -DDYNAMIC_ARCH=ON -DCMAKE_BUILD_TYPE=Release
   ```

   You may want to add further options in the `cmake` command here &ndash; for instance, the default only produces a static .lib version of the library. If you would rather have a DLL, add -DBUILD_SHARED_LIBS=ON above. Note that this step only creates some command files and directories, the actual build happens next.


7. Build the project:

   ```
   cmake --build . --config Release
   ```
   This step will create the OpenBLAS library in the "lib" directory, and various build-time tests in the `test`, `ctest` and `openblas_utest` directories. However it will not separate the header files you might need for building your own programs from those used internally. To put all relevant files in a more convenient arrangement, run the next step.

8. Install all relevant files created by the build

   ```
   cmake --install . --prefix c:\opt -v
   ```
   This will copy all files that are needed for building and running your own programs with OpenBLAS to the given location, creating appropriate subdirectories for the individual kinds of files. In the case of "C:\opt" as given above, this would be C:\opt\include\openblas for the header files, 
   C:\opt\bin for the libopenblas.dll and C:\opt\lib for the static library. C:\opt\share holds various support files that enable other cmake-based build scripts to find OpenBLAS automatically.

###### Visual studio 2017+ (C++2017 standard)

In newer visual studio versions, Microsoft has changed [how it handles complex types](https://docs.microsoft.com/en-us/cpp/c-runtime-library/complex-math-support?view=msvc-170#types-used-in-complex-math). Even when using a precompiled version of OpenBLAS, you might need to define `LAPACK_COMPLEX_CUSTOM` in order to define complex types properly for MSVC. For example, some variant of the following might help:

```
#if defined(_MSC_VER)
    #include <complex.h>
    #define LAPACK_COMPLEX_CUSTOM
    #define lapack_complex_float _Fcomplex
    #define lapack_complex_double _Dcomplex
#endif
```

For reference, see https://github.com/xianyi/OpenBLAS/issues/3661, https://github.com/Reference-LAPACK/lapack/issues/683, and https://stackoverflow.com/questions/47520244/using-openblas-lapacke-in-visual-studio.

###### CMake and Visual Studio

To build OpenBLAS for the 32-bit architecture, you'll need to use the builtin Visual Studio compilers.

!!! note
    This method may produce binaries which demonstrate significantly lower performance than those built with the other methods. (The Visual Studio compiler does not support the dialect of assembly used in the cpu-specific optimized files, so only the "generic" TARGET which is
    written in pure C will get built. For the same reason it is not possible (and not necessary) to use -DDYNAMIC_ARCH=ON in a Visual Studio build) You may consider building for the 32-bit architecture using the GNU (MinGW) ABI.

####### 1. Install CMake at Windows

####### 2. Use CMake to generate Visual Studio solution files

```
# Do this from Powershell so cmake can find visual studio
cmake -G "Visual Studio 14 Win64" -DCMAKE_BUILD_TYPE=Release .
```

###### Build the solution at Visual Studio

Note that this step depends on perl, so you'll need to install perl for windows, and put perl on your path so VS can start perl (http://stackoverflow.com/questions/3051049/active-perl-installation-on-windows-operating-system).

Step 2 will build the OpenBLAS solution, open it in VS, and build the projects. Note that the dependencies do not seem to be automatically configured: if you try to build libopenblas directly, it will fail with a message saying that some .obj files aren't found, but if you build the projects libopenblas depends on before building libopenblas, the build will succeed.

###### Build OpenBLAS for Universal Windows Platform

OpenBLAS can be built for use on the [Universal Windows Platform](https://en.wikipedia.org/wiki/Universal_Windows_Platform) using a two step process since commit [c66b842](https://github.com/xianyi/OpenBLAS/commit/c66b842d66c5516e52804bf5a0544d18b1da1b44).

####### 1. Follow steps 1 and 2 above to build the Visual Studio solution files for Windows. This builds the helper executables which are required when building the OpenBLAS Visual Studio solution files for UWP in step 2.

####### 2. Remove the generated CMakeCache.txt and CMakeFiles directory from the OpenBLAS source directory and re-run CMake with the following options:

```
# do this to build UWP compatible solution files
cmake -G "Visual Studio 14 Win64" -DCMAKE_SYSTEM_NAME=WindowsStore -DCMAKE_SYSTEM_VERSION="10.0" -DCMAKE_SYSTEM_PROCESSOR=AMD64 -DVS_WINRT_COMPONENT=TRUE -DCMAKE_BUILD_TYPE=Release .
```

####### Build the solution with Visual Studio

This will build the OpenBLAS binaries with the required settings for use with UWP.

##### 2. GNU (MinGW) ABI

The resulting library can be used in Visual Studio, but it can only be linked dynamically. This configuration has not been thoroughly tested and should be considered experimental.

###### Incompatible x86 calling conventions

Due to incompatibilities between the calling conventions of MinGW and Visual Studio you will need to make the following modifications ( **32-bit only** ):

1. Use the newer GCC 4.7.0. The older GCC (<4.7.0) has an ABI incompatibility for returning aggregate structures larger than 8 bytes with MSVC. 


###### Build OpenBLAS on Windows OS
1. Install the MinGW (GCC) compiler suite, either 32-bit (http://www.mingw.org/) or 64-bit (http://mingw-w64.sourceforge.net/). Be sure to install its gfortran package as well (unless you really want to build the BLAS part of OpenBLAS only) and check that gcc and gfortran are the same version &ndash; mixing compilers from different sources or release versions can lead to strange error messages in the linking stage. In addition, please install MSYS with MinGW.
1. Build OpenBLAS in the MSYS shell. Usually, you can just type "make". OpenBLAS will detect the compiler and CPU automatically. 
1. After the build is complete, OpenBLAS will generate the static library "libopenblas.a" and the shared dll library "libopenblas.dll" in the folder. You can type "make PREFIX=/your/installation/path install" to install the library to a certain location.

!!! note
    We suggest using official MinGW or MinGW-w64 compilers. A user reported that s/he met `Unhandled exception` by other compiler suite. https://groups.google.com/forum/#!topic/openblas-users/me2S4LkE55w

Note also that older versions of the alternative builds of mingw-w64 available through http://www.msys2.org may contain a defect that leads to a compilation failure accompanied by the error message
```
<command-line>:0:4: error: expected identifier or '(' before numeric constant
```
If you encounter this, please upgrade your msys2 setup or see https://github.com/xianyi/OpenBLAS/issues/1503 for a workaround.

###### Generate import library (before 0.2.10 version)

1. First, you will need to have the `lib.exe` tool in the Visual Studio command prompt. 
1. Open the command prompt and type `cd OPENBLAS_TOP_DIR/exports`, where OPENBLAS_TOP_DIR is the main folder of your OpenBLAS installation.
1. For a 32-bit library, type `lib /machine:i386 /def:libopenblas.def`. For 64-bit, type `lib /machine:X64 /def:libopenblas.def`.
1. This will generate the import library "libopenblas.lib" and the export library "libopenblas.exp" in OPENBLAS_TOP_DIR/exports. Although these two files have the same name, they are totally different.

###### Generate import library (0.2.10 and after version)
1. OpenBLAS already generated the import library "libopenblas.dll.a" for "libopenblas.dll".

###### generate windows native PDB files from gcc/gfortran build
Tool to do so is available at https://github.com/rainers/cv2pdb

###### Use OpenBLAS .dll library in Visual Studio
1. Copy the import library (before 0.2.10: "OPENBLAS_TOP_DIR/exports/libopenblas.lib", 0.2.10 and after: "OPENBLAS_TOP_DIR/libopenblas.dll.a") and .dll library "libopenblas.dll" into the same folder(The folder of your project that is going to use the BLAS library. You may need to add the libopenblas.dll.a to the linker input list: properties->Linker->Input).
1. Please follow the documentation about using third-party .dll libraries in MS Visual Studio 2008 or 2010.  Make sure to link against a library for the correct architecture.  For example, you may receive an error such as "The application was unable to start correctly (0xc000007b)" which typically indicates a mismatch between 32/64-bit libraries.

!!! note
    If you need CBLAS, you should include cblas.h in /your/installation/path/include in Visual Studio.  Please read [this page](http://github.com/xianyi/OpenBLAS/issues/95).

###### Limitations
* Both static and dynamic linking are supported with MinGW.  With Visual Studio, however, only dynamic linking is supported and so you should use the import library.
* Debugging from Visual Studio does not work because MinGW and Visual Studio have incompatible formats for debug information (PDB vs. DWARF/STABS).  You should either debug with GDB on the command-line or with a visual frontend, for instance [Eclipse](http://www.eclipse.org/cdt/) or [Qt Creator](http://qt.nokia.com/products/developer-tools/).


#### Windows on Arm

##### Prerequisites

Following tools needs to be installed

###### 1. Download and install clang for windows on arm

Find the latest LLVM build for WoA from [LLVM release page](https://releases.llvm.org/)

E.g: LLVM 12 build for WoA64 can be found [here](https://github.com/llvm/llvm-project/releases/download/llvmorg-12.0.0/LLVM-12.0.0-woa64.exe)

Run the LLVM installer and ensure that LLVM is added to environment PATH.

###### 2. Download and install classic flang for windows on arm

Classic flang is the only available FORTRAN compiler for windows on arm for now and a pre-release build can be found [here](https://github.com/kaadam/flang/releases/tag/v0.1)

There is no installer for classic flang and the zip package can be extracted and the path needs to be added to environment PATH.

E.g: on PowerShell

```
$env:Path += ";C:\flang_woa\bin"
```

##### Build

The following steps describe how to build the static library for OpenBLAS with and without LAPACK

###### 1. Build OpenBLAS static library with BLAS and LAPACK routines with Make

Following command can be used to build OpenBLAS static library with BLAS and LAPACK routines

```bash
$ make CC="clang-cl" HOSTCC="clang-cl" AR="llvm-ar" BUILD_WITHOUT_LAPACK=0 NOFORTRAN=0 DYNAMIC_ARCH=0 TARGET=ARMV8 ARCH=arm64 BINARY=64 USE_OPENMP=0 PARALLEL=1 RANLIB="llvm-ranlib" MAKE=make F_COMPILER=FLANG FC=FLANG FFLAGS_NOOPT="-march=armv8-a -cpp" FFLAGS="-march=armv8-a -cpp" NEED_PIC=0 HOSTARCH=arm64 libs netlib
```

###### 2. Build static library with BLAS routines using CMake

Classic flang has compatibility issues with cmake hence only BLAS routines can be compiled with CMake

```bash
$ mkdir build
$ cd build
$ cmake ..  -G Ninja -DCMAKE_C_COMPILER=clang -DBUILD_WITHOUT_LAPACK=1 -DNOFORTRAN=1 -DDYNAMIC_ARCH=0 -DTARGET=ARMV8 -DARCH=arm64 -DBINARY=64 -DUSE_OPENMP=0 -DCMAKE_SYSTEM_PROCESSOR=ARM64 -DCMAKE_CROSSCOMPILING=1 -DCMAKE_SYSTEM_NAME=Windows
$ cmake --build . --config Release
```

###### `getarch.exe` execution error

If you notice that platform-specific headers by `getarch.exe` are not generated correctly, It could be due to a known debug runtime DLL issue for arm64 platforms. Please check out [link](https://linaro.atlassian.net/wiki/spaces/WOAR/pages/28677636097/Debug+run-time+DLL+issue#Workaround) for the workaround.

#### MinGW import library

Microsoft Windows has this thing called "import libraries". You don't need it in MinGW because the `ld` linker from GNU Binutils is smart, but you may still want it for whatever reason.

##### Make the `.def`

Import libraries are compiled from a list of what symbols to use, `.def`. This should be already in your `exports` directory: `cd OPENBLAS_TOP_DIR/exports`.

##### Making a MinGW import library

MinGW import libraries have the suffix `.a`, same as static libraries. (It's actually more common to do `.dll.a`...)

You need to first prepend `libopenblas.def` with a line `LIBRARY libopenblas.dll`:

    cat <(echo "LIBRARY libopenblas.dll") libopenblas.def > libopenblas.def.1
    mv libopenblas.def.1 libopenblas.def

Now it probably looks like:

    LIBRARY libopenblas.dll
    EXPORTS
	   caxpy=caxpy_  @1
	   caxpy_=caxpy_  @2
           ...

Then, generate the import library: `dlltool -d libopenblas.def -l libopenblas.a`

Again, there is basically **no point** in making an import library for use in MinGW. It actually slows down linking.

##### Making a MSVC import library

Unlike MinGW, MSVC absolutely requires an import library. Now the C ABI of MSVC and MinGW are actually identical, so linking is actually okay. (Any incompatibility in the C ABI would be a bug.)

The import libraries of MSVC have the suffix `.lib`. They are generated from a `.def` file using MSVC's `lib.exe`. See [the MSVC instructions](use_visual_studio.md#generate-import-library-before-0210-version).

##### Notes

* Always remember that MinGW is **not the same** as MSYS2 or Cygwin. MSYS2 and Cygwin are full POSIX environments with a lot of magic such as `fork()` and its own `malloc()`. MinGW, which builds on the normal Microsoft C Runtime, has none of that. Be clear about which one you are building for.

### Mac OSX

If your CPU is Sandy Bridge, please use Clang version 3.1 and above. The Clang 3.0 will generate the wrong AVX binary code of OpenBLAS.

#### Precompiled packages

<https://www.macports.org/ports.php?by=name&substr=openblas>

`brew install openblas`

or using the conda package manager from
<https://github.com/conda-forge/miniforge#download>
(which also has packages for the new M1 cpu)

 `conda install openblas`

#### Build on Apple M1

On newer versions of Xcode and on arm64, you might need to compile with a newer macOS target (11.0) than the default (10.8) with `MACOSX_DEPLOYMENT_TARGET=11.0`, or switch your command-line tools to use an older SDK (e.g., [13.1](https://developer.apple.com/download/all/?q=Xcode%2013)).

* without Fortran compiler （cannot build LAPACK)
  ```bash
  $ make CC=cc NOFORTRAN=1
  ```
* with Fortran compiler (you could `brew install gfortran`) https://github.com/xianyi/OpenBLAS/issues/3032
  ```bash
  $ export MACOSX_DEPLOYMENT_TARGET=11.0
  $ make CC=cc FC=gfortran
  ```

### Android

#### Prerequisites

In addition to the Android NDK, you will need both Perl and a C compiler on the build host as these are currently
required by the OpenBLAS build environment.


#### Building with android NDK using clang compiler
Around version 11 Android NDKs stopped supporting gcc, so you would need to use clang to compile OpenBLAS. clang is supported from OpenBLAS 0.2.20 version onwards. See below sections on how to build with clang for ARMV7 and ARMV8 targets. The same basic principles as described below for ARMV8 should also apply to building an x86 or x86_64 version (substitute something like NEHALEM for the target instead of ARMV8 and replace all the aarch64 in the toolchain paths obviously)
"Historic" notes:
Since version 19 the default toolchain is provided as a standalone toolchain, so building one yourself following [building a standalone toolchain](http://developer.android.com/ndk/guides/standalone_toolchain.html) should no longer be necessary.
If you want to use static linking with an old NDK version older than about r17, you need to choose an API level below 23 currently due to NDK bug 272 (https://github.com/android-ndk/ndk/issues/272 , the libc.a lacks a definition of stderr) that will probably be fixed in r17 of the NDK.

#### Build ARMV7 with clang
```
## Set path to ndk-bundle
export NDK_BUNDLE_DIR=/path/to/ndk-bundle

## Set the PATH to contain paths to clang and arm-linux-androideabi-* utilities
export PATH=${NDK_BUNDLE_DIR}/toolchains/arm-linux-androideabi-4.9/prebuilt/linux-x86_64/bin:${NDK_BUNDLE_DIR}/toolchains/llvm/prebuilt/linux-x86_64/bin:$PATH

## Set LDFLAGS so that the linker finds the appropriate libgcc
export LDFLAGS="-L${NDK_BUNDLE_DIR}/toolchains/arm-linux-androideabi-4.9/prebuilt/linux-x86_64/lib/gcc/arm-linux-androideabi/4.9.x"

## Set the clang cross compile flags
export CLANG_FLAGS="-target arm-linux-androideabi -marm -mfpu=vfp -mfloat-abi=softfp --sysroot ${NDK_BUNDLE_DIR}/platforms/android-23/arch-arm -gcc-toolchain ${NDK_BUNDLE_DIR}/toolchains/arm-linux-androideabi-4.9/prebuilt/linux-x86_64/"

#OpenBLAS Compile
make TARGET=ARMV7 ONLY_CBLAS=1 AR=ar CC="clang ${CLANG_FLAGS}" HOSTCC=gcc ARM_SOFTFP_ABI=1 -j4
```
On a Mac, it may also be necessary to give the complete path to the `ar` utility in the make command above, like so:
```
AR=${NDK_BUNDLE_DIR}/toolchains/arm-linux-androideabi-4.9/prebuilt/darwin-x86_64/bin/arm-linux-androideabi-gcc-ar
```
otherwise you may get a linker error complaining about a "malformed archive header name at 8" when the native OSX ar command was invoked instead.
 
#### Build ARMV8 with clang
```
## Set path to ndk-bundle
export NDK_BUNDLE_DIR=/path/to/ndk-bundle/

## Export PATH to contain directories of clang and aarch64-linux-android-* utilities
export PATH=${NDK_BUNDLE_DIR}/toolchains/aarch64-linux-android-4.9/prebuilt/linux-x86_64/bin/:${NDK_BUNDLE_DIR}/toolchains/llvm/prebuilt/linux-x86_64/bin:$PATH

## Setup LDFLAGS so that loader can find libgcc and pass -lm for sqrt
export LDFLAGS="-L${NDK_BUNDLE_DIR}/toolchains/aarch64-linux-android-4.9/prebuilt/linux-x86_64/lib/gcc/aarch64-linux-android/4.9.x -lm"

## Setup the clang cross compile options
export CLANG_FLAGS="-target aarch64-linux-android --sysroot ${NDK_BUNDLE_DIR}/platforms/android-23/arch-arm64 -gcc-toolchain ${NDK_BUNDLE_DIR}/toolchains/aarch64-linux-android-4.9/prebuilt/linux-x86_64/"

## Compile
make TARGET=ARMV8 ONLY_CBLAS=1 AR=ar CC="clang ${CLANG_FLAGS}" HOSTCC=gcc -j4
```
Note: Using TARGET=CORTEXA57 in place of ARMV8 will pick up better optimized routines. Implementations for CORTEXA57 target is compatible with all other armv8 targets.

Note: For NDK 23b, something as simple as 
```
export PATH=/opt/android-ndk-r23b/toolchains/llvm/prebuilt/linux-x86_64/bin/:$PATH
make HOSTCC=gcc CC=/opt/android-ndk-r23b/toolchains/llvm/prebuilt/linux-x86_64/bin/aarch64-linux-android31-clang ONLY_CBLAS=1 TARGET=ARMV8
```
appears to be sufficient on Linux.

#### Alternative script which was tested on OSX with NDK(21.3.6528147)
This script will build openblas for 3 architecture (ARMV7,ARMV8,X86) and put them with `sudo make install` to `/opt/OpenBLAS/lib`
```
export NDK=YOUR_PATH_TO_SDK/Android/sdk/ndk/21.3.6528147
export TOOLCHAIN=$NDK/toolchains/llvm/prebuilt/darwin-x86_64

make clean
make \
    TARGET=ARMV7 \
    ONLY_CBLAS=1 \
    CC="$TOOLCHAIN"/bin/armv7a-linux-androideabi21-clang \
    AR="$TOOLCHAIN"/bin/arm-linux-androideabi-ar \
    HOSTCC=gcc \
    ARM_SOFTFP_ABI=1 \
    -j4
sudo make install

make clean
make \
    TARGET=CORTEXA57 \
    ONLY_CBLAS=1 \
    CC=$TOOLCHAIN/bin/aarch64-linux-android21-clang \
    AR=$TOOLCHAIN/bin/aarch64-linux-android-ar \
    HOSTCC=gcc \
    -j4
sudo make install

make clean
make \
    TARGET=ATOM \
    ONLY_CBLAS=1 \
    CC="$TOOLCHAIN"/bin/i686-linux-android21-clang \
    AR="$TOOLCHAIN"/bin/i686-linux-android-ar \
    HOSTCC=gcc \
    ARM_SOFTFP_ABI=1 \
    -j4
sudo make install

## This will build for x86_64 
make clean
make \
    TARGET=ATOM BINARY=64\
    ONLY_CBLAS=1 \
    CC="$TOOLCHAIN"/bin/x86_64-linux-android21-clang \
    AR="$TOOLCHAIN"/bin/x86_64-linux-android-ar \
    HOSTCC=gcc \
    ARM_SOFTFP_ABI=1 \
    -j4
sudo make install
```
Also you can find full list of target architectures in [TargetsList.txt](https://github.com/xianyi/OpenBLAS/blob/develop/TargetList.txt)

***
anything below this line should be irrelevant nowadays unless you need to perform software archeology
***
#### Building OpenBLAS with very old gcc-based versions of the NDK, without Fortran

The prebuilt Android NDK toolchains do not include Fortran, hence parts like LAPACK cannot be built. You can still build OpenBLAS without it. For instructions on how to build OpenBLAS with Fortran, see the [next section](#building-openblas-with-fortran).

To use easily the prebuilt toolchains, follow [building a standalone toolchain](http://developer.android.com/ndk/guides/standalone_toolchain.html) for your desired architecture.
This would be `arm-linux-androideabi-gcc-4.9` for ARMV7 and `aarch64-linux-android-gcc-4.9` for ARMV8.

You can build OpenBLAS (0.2.19 and earlier) with:
```
## Add the toolchain to your path
export PATH=/path/to/standalone-toolchain/bin:$PATH

## Build without Fortran for ARMV7
make TARGET=ARMV7 HOSTCC=gcc CC=arm-linux-androideabi-gcc NOFORTRAN=1 libs
## Build without Fortran for ARMV8
make TARGET=ARMV8 BINARY=64 HOSTCC=gcc CC=aarch64-linux-android-gcc NOFORTRAN=1 libs
```

Since we are cross-compiling, we make the `libs` recipe, not `all`. Otherwise you will get errors when trying to link/run tests as versions up to and including 0.2.19 cannot build a shared library for Android. 

From 0.2.20 on, you should leave off the "libs" to get a full build, and you may want to use the softfp ABI instead of the deprecated hardfp one on ARMV7 so you would use 
```
## Add the toolchain to your path
export PATH=/path/to/standalone-toolchain/bin:$PATH

## Build without Fortran for ARMV7
make TARGET=ARMV7 ARM_SOFTFP_ABI=1 HOSTCC=gcc CC=arm-linux-androideabi-gcc NOFORTRAN=1
## Build without Fortran for ARMV8
make TARGET=ARMV8 BINARY=64 HOSTCC=gcc CC=aarch64-linux-android-gcc NOFORTRAN=1
```

If you get an error about stdio.h not being found, you need to specify your sysroot in the CFLAGS argument to `make` like
```CFLAGS=--sysroot=$NDK/platforms/android-16/arch-arm```
When you are done, install OpenBLAS into the desired directory. Be sure to also use all command line options
here that you specified for building, otherwise errors may occur as it tries to install things you did not build:
```
make PREFIX=/path/to/install-dir TARGET=... install
```

#### Building OpenBLAS with Fortran

Instructions on how to build the GNU toolchains with Fortran can be found [here](https://github.com/buffer51/android-gfortran). The [Releases section](https://github.com/buffer51/android-gfortran/releases) provides prebuilt versions, use the standalone one.

You can build OpenBLAS with:
```
## Add the toolchain to your path
export PATH=/path/to/standalone-toolchain-with-fortran/bin:$PATH

## Build with Fortran for ARMV7
make TARGET=ARMV7 HOSTCC=gcc CC=arm-linux-androideabi-gcc FC=arm-linux-androideabi-gfortran libs
## Build with LAPACK for ARMV8
make TARGET=ARMV8 BINARY=64 HOSTCC=gcc CC=aarch64-linux-android-gcc FC=aarch64-linux-android-gfortran libs
```

As mentioned above you can leave off the `libs` argument here when building 0.2.20 and later, and you may want to add ARM_SOFTFP_ABI=1 when building for ARMV7.

#### Linking OpenBLAS (0.2.19 and earlier) for ARMV7

If you are using `ndk-build`, you need to set the ABI to hard floating points in your Application.mk:
```
APP_ABI := armeabi-v7a-hard
```

This will set the appropriate flags for you. If you are not using `ndk-build`, you will want to add the following flags:
```
TARGET_CFLAGS += -mhard-float -D_NDK_MATH_NO_SOFTFP=1
TARGET_LDFLAGS += -Wl,--no-warn-mismatch -lm_hard
```

From 0.2.20 on, it is also possible to build for the softfp ABI by specifying ARM_SOFTFP_ABI=1 during the build.
In that case, also make sure that all your dependencies are compiled with -mfloat-abi=softfp as well, as mixing
"hard" and "soft" floating point ABIs in a program will make it crash.

### iPhone/iOS

As none of the current developers uses iOS, the following instructions are what was found to work in our Azure CI setup, but as far as we know this builds a fully working OpenBLAS for this platform.

Go to the directory where you unpacked OpenBLAS,and enter the following commands:
```
     CC=/Applications/Xcode_12.4.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang

CFLAGS= -O2 -Wno-macro-redefined -isysroot /Applications/Xcode_12.4.app/Contents/Developer/Platforms/iPhoneOS.platform/Developer/SDKs/iPhoneOS14.4.sdk -arch arm64 -miphoneos-version-min=10.0

make TARGET=ARMV8 DYNAMIC_ARCH=1 NUM_THREADS=32 HOSTCC=clang NOFORTRAN=1
```
Adjust MIN_IOS_VERSION as necessary for your installation, e.g. change the version number
to the minimum iOS version you want to target and execute this file to build the library.

### MIPS

For mips targets you will need latest toolchains
P5600 - MTI GNU/Linux Toolchain
I6400, P6600 - IMG GNU/Linux Toolchain

The download link is below
(http://codescape-mips-sdk.imgtec.com/components/toolchain/2016.05-03/downloads.html)

You can use following commandlines for builds


    IMG_TOOLCHAIN_DIR={full IMG GNU/Linux Toolchain path including "bin" directory -- for example, /opt/linux_toolchain/bin}
    IMG_GCC_PREFIX=mips-img-linux-gnu
    IMG_TOOLCHAIN=${IMG_TOOLCHAIN_DIR}/${IMG_GCC_PREFIX}

    I6400 Build (n32):
    make BINARY=32 BINARY32=1 CC=$IMG_TOOLCHAIN-gcc AR=$IMG_TOOLCHAIN-ar FC="$IMG_TOOLCHAIN-gfortran -EL -mabi=n32" RANLIB=$IMG_TOOLCHAIN-ranlib HOSTCC=gcc CFLAGS="-EL" FFLAGS=$CFLAGS LDFLAGS=$CFLAGS TARGET=I6400

    I6400 Build (n64):
    make BINARY=64 BINARY64=1 CC=$IMG_TOOLCHAIN-gcc AR=$IMG_TOOLCHAIN-ar FC="$IMG_TOOLCHAIN-gfortran -EL" RANLIB=$IMG_TOOLCHAIN-ranlib HOSTCC=gcc CFLAGS="-EL" FFLAGS=$CFLAGS LDFLAGS=$CFLAGS TARGET=I6400

    P6600 Build (n32):
    make BINARY=32 BINARY32=1 CC=$IMG_TOOLCHAIN-gcc AR=$IMG_TOOLCHAIN-ar FC="$IMG_TOOLCHAIN-gfortran -EL -mabi=n32" RANLIB=$IMG_TOOLCHAIN-ranlib HOSTCC=gcc CFLAGS="-EL" FFLAGS=$CFLAGS LDFLAGS=$CFLAGS TARGET=P6600

    P6600 Build (n64):
    make BINARY=64 BINARY64=1 CC=$IMG_TOOLCHAIN-gcc AR=$IMG_TOOLCHAIN-ar FC="$IMG_TOOLCHAIN-gfortran -EL" RANLIB=$IMG_TOOLCHAIN-ranlib HOSTCC=gcc CFLAGS="-EL" FFLAGS="$CFLAGS" LDFLAGS="$CFLAGS" TARGET=P6600

    MTI_TOOLCHAIN_DIR={full MTI GNU/Linux Toolchain path including "bin" directory -- for example, /opt/linux_toolchain/bin}
    MTI_GCC_PREFIX=mips-mti-linux-gnu
    MTI_TOOLCHAIN=${IMG_TOOLCHAIN_DIR}/${IMG_GCC_PREFIX}

    P5600 Build:

    make BINARY=32 BINARY32=1 CC=$MTI_TOOLCHAIN-gcc AR=$MTI_TOOLCHAIN-ar FC="$MTI_TOOLCHAIN-gfortran -EL"    RANLIB=$MTI_TOOLCHAIN-ranlib HOSTCC=gcc CFLAGS="-EL" FFLAGS=$CFLAGS LDFLAGS=$CFLAGS TARGET=P5600

### FreeBSD

You will need to install the following tools from the FreeBSD ports tree:
* lang/gcc [1]
* lang/perl5.12
* ftp/curl
* devel/gmake
* devel/patch

To compile run the command:

    $ gmake CC=gcc46 FC=gfortran46

Note that you need to build with GNU make and manually specify the compiler, otherwhise gcc 4.2 from the base system would be used.

[1]: [Removal of Fortran from the FreeBSD base system](http://www.bsdunix.ch/serendipity/index.php?/archives/345-Removal-of-Fortran-from-the-FreeBSD-base-system.html)


```
pkg install openblas
```

see <https://www.freebsd.org/ports/index.html>

### Cortex-M

Cortex-M is a widely used microcontroller that is present in a variety of industrial and consumer electronics.
A common variant of the Cortex-M is the STM32F4xx series. Here, we will give instructions for building for
the STM32F4xx.

First, install the embedded arm gcc compiler from the arm website. Then, create the following toolchain file and build as follows.

```cmake
# cmake .. -G Ninja -DCMAKE_C_COMPILER=arm-none-eabi-gcc -DCMAKE_TOOLCHAIN_FILE:PATH="toolchain.cmake" -DNOFORTRAN=1 -DTARGET=ARMV5 -DEMBEDDED=1

set(CMAKE_SYSTEM_NAME Generic)
set(CMAKE_SYSTEM_PROCESSOR arm)

set(CMAKE_C_COMPILER "arm-none-eabi-gcc.exe")
set(CMAKE_CXX_COMPILER "arm-none-eabi-g++.exe")

set(CMAKE_EXE_LINKER_FLAGS "--specs=nosys.specs" CACHE INTERNAL "")

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)
```

In your embedded application, the following functions need to be provided for OpenBLAS to work correctly:

```C
void free(void* ptr);
void* malloc(size_t size);
```

!!! note
    If you are developing for an embedded platform, it is your responsibility to make sure that the device has sufficient memory for malloc calls. [Libmemory][2] provides one implementation of malloc for embedded platforms.

[2]: https://github.com/embeddedartistry/libmemory
