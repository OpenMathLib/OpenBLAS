As of OpenBLAS v0.2.15, we support MinGW and Visual Studio (using CMake to generate visual studio solution files &ndash; note that you will need at least version 3.11 of CMake for linking to work correctly) to build OpenBLAS on Windows.

Note that you need a Fortran compiler if you plan to build and use the LAPACK functions included with OpenBLAS. The sections below describe using either `flang` as an add-on to clang/LLVM or `gfortran` as part of MinGW for this purpose. If you want to use the Intel Fortran compiler `ifort` for this, be sure to also use the Intel C compiler `icc` for building the C parts, as the ABI imposed by `ifort` is incompatible with `msvc`.

## 1. Native (MSVC) ABI

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

### Visual studio 2017+ (C++2017 standard)

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

### CMake and Visual Studio

To build OpenBLAS for the 32-bit architecture, you'll need to use the builtin Visual Studio compilers.

**[Notice]** This method may produce binaries which demonstrate significantly lower performance than those built with the other methods. (The Visual Studio compiler does not support the dialect of assembly used in the cpu-specific optimized files, so only the "generic" TARGET which is
written in pure C will get built. For the same reason it is not possible (and not necessary) to use -DDYNAMIC_ARCH=ON in a Visual Studio build) You may consider building for the 32-bit architecture using the GNU (MinGW) ABI.

#### 1. Install CMake at Windows

#### 2. Use CMake to generate Visual Studio solution files

```
# Do this from Powershell so cmake can find visual studio
cmake -G "Visual Studio 14 Win64" -DCMAKE_BUILD_TYPE=Release .
```

### Build the solution at Visual Studio

Note that this step depends on perl, so you'll need to install perl for windows, and put perl on your path so VS can start perl (http://stackoverflow.com/questions/3051049/active-perl-installation-on-windows-operating-system).

Step 2 will build the OpenBLAS solution, open it in VS, and build the projects. Note that the dependencies do not seem to be automatically configured: if you try to build libopenblas directly, it will fail with a message saying that some .obj files aren't found, but if you build the projects libopenblas depends on before building libopenblas, the build will succeed.

### Build OpenBLAS for Universal Windows Platform

OpenBLAS can be built for use on the [Universal Windows Platform](https://en.wikipedia.org/wiki/Universal_Windows_Platform) using a two step process since commit [c66b842](https://github.com/xianyi/OpenBLAS/commit/c66b842d66c5516e52804bf5a0544d18b1da1b44).

#### 1. Follow steps 1 and 2 above to build the Visual Studio solution files for Windows. This builds the helper executables which are required when building the OpenBLAS Visual Studio solution files for UWP in step 2.

#### 2. Remove the generated CMakeCache.txt and CMakeFiles directory from the OpenBLAS source directory and re-run CMake with the following options:

```
# do this to build UWP compatible solution files
cmake -G "Visual Studio 14 Win64" -DCMAKE_SYSTEM_NAME=WindowsStore -DCMAKE_SYSTEM_VERSION="10.0" -DCMAKE_SYSTEM_PROCESSOR=AMD64 -DVS_WINRT_COMPONENT=TRUE -DCMAKE_BUILD_TYPE=Release .
```

#### Build the solution with Visual Studio

This will build the OpenBLAS binaries with the required settings for use with UWP.

## 2. GNU (MinGW) ABI

The resulting library can be used in Visual Studio, but it can only be linked dynamically. This configuration has not been thoroughly tested and should be considered experimental.

### Incompatible x86 calling conventions

Due to incompatibilities between the calling conventions of MinGW and Visual Studio you will need to make the following modifications ( **32-bit only** ):

1. Use the newer GCC 4.7.0. The older GCC (<4.7.0) has an ABI incompatibility for returning aggregate structures larger than 8 bytes with MSVC. 


### Build OpenBLAS on Windows OS
1. Install the MinGW (GCC) compiler suite, either 32-bit (http://www.mingw.org/) or 64-bit (http://mingw-w64.sourceforge.net/). Be sure to install its gfortran package as well (unless you really want to build the BLAS part of OpenBLAS only) and check that gcc and gfortran are the same version &ndash; mixing compilers from different sources or release versions can lead to strange error messages in the linking stage. In addition, please install MSYS with MinGW.
1. Build OpenBLAS in the MSYS shell. Usually, you can just type "make". OpenBLAS will detect the compiler and CPU automatically. 
1. After the build is complete, OpenBLAS will generate the static library "libopenblas.a" and the shared dll library "libopenblas.dll" in the folder. You can type "make PREFIX=/your/installation/path install" to install the library to a certain location.

**[Notice]** We suggest using official MinGW or MinGW-w64 compilers. A user reported that s/he met `Unhandled exception` by other compiler suite. https://groups.google.com/forum/#!topic/openblas-users/me2S4LkE55w

Note also that older versions of the alternative builds of mingw-w64 available through http://www.msys2.org may contain a defect that leads to a compilation failure accompanied by the error message
```
<command-line>:0:4: error: expected identifier or '(' before numeric constant
```
If you encounter this, please upgrade your msys2 setup or see https://github.com/xianyi/OpenBLAS/issues/1503 for a workaround.

### Generate import library (before 0.2.10 version)

1. First, you will need to have the `lib.exe` tool in the Visual Studio command prompt. 
1. Open the command prompt and type `cd OPENBLAS_TOP_DIR/exports`, where OPENBLAS_TOP_DIR is the main folder of your OpenBLAS installation.
1. For a 32-bit library, type `lib /machine:i386 /def:libopenblas.def`. For 64-bit, type `lib /machine:X64 /def:libopenblas.def`.
1. This will generate the import library "libopenblas.lib" and the export library "libopenblas.exp" in OPENBLAS_TOP_DIR/exports. Although these two files have the same name, they are totally different.

### Generate import library (0.2.10 and after version)
1. OpenBLAS already generated the import library "libopenblas.dll.a" for "libopenblas.dll".

### generate windows native PDB files from gcc/gfortran build
Tool to do so is available at https://github.com/rainers/cv2pdb

### Use OpenBLAS .dll library in Visual Studio
1. Copy the import library (before 0.2.10: "OPENBLAS_TOP_DIR/exports/libopenblas.lib", 0.2.10 and after: "OPENBLAS_TOP_DIR/libopenblas.dll.a") and .dll library "libopenblas.dll" into the same folder(The folder of your project that is going to use the BLAS library. You may need to add the libopenblas.dll.a to the linker input list: properties->Linker->Input).
1. Please follow the documentation about using third-party .dll libraries in MS Visual Studio 2008 or 2010.  Make sure to link against a library for the correct architecture.  For example, you may receive an error such as "The application was unable to start correctly (0xc000007b)" which typically indicates a mismatch between 32/64-bit libraries.

**[Notice]** If you need CBLAS, you should include cblas.h in /your/installation/path/include in Visual Studio.  Please read [this page](http://github.com/xianyi/OpenBLAS/issues/95).

### Limitations
* Both static and dynamic linking are supported with MinGW.  With Visual Studio, however, only dynamic linking is supported and so you should use the import library.
* Debugging from Visual Studio does not work because MinGW and Visual Studio have incompatible formats for debug information (PDB vs. DWARF/STABS).  You should either debug with GDB on the command-line or with a visual frontend, for instance [Eclipse](http://www.eclipse.org/cdt/) or [Qt Creator](http://qt.nokia.com/products/developer-tools/).
