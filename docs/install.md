## Quick Installation
[Precompiled packages](https://github.com/xianyi/OpenBLAS/wiki/Precompiled-installation-packages) have recently become available for a number of platforms through their normal installation procedures, so for users of desktop devices at least, the instructions below are mostly relevant when you want to try the most recent development snapshot from git.

### Linux

Just type `make` to compile the library.

Notes
* OpenBLAS doesn't support g77. Please use gfortran or other Fortran compilers. e.g. `make FC=gfortran`.
* When building in an emulator (KVM,QEMU etc.) make sure that the combination of CPU features exposed to
  the virtual environment matches that of an existing CPU to allow detection of the cpu model to succeed. 
  (With qemu, this can be done by passing `-cpu host` or a supported model name at invocation)

### Windows

See [[How-to-use-OpenBLAS-in-Microsoft-Visual-Studio]].

The precompiled binaries available with each release (in https://github.com/xianyi/OpenBLAS/releases) are
created with MinGW (as described in Section 2 of the wiki page mentioned above) using an option list of 
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

### Mac OSX

If your CPU is Sandy Bridge, please use Clang version 3.1 and above. The Clang 3.0 will generate the wrong AVX binary code of OpenBLAS.

#### Build on Apple M1

* without Fortran compiler （cannot build LAPACK)
```
    $ make CC=cc NOFORTRAN=1
```
* with Fortran compiler (you could `brew install gfortran`) https://github.com/xianyi/OpenBLAS/issues/3032
```
    $ export MACOSX_DEPLOYMENT_TARGET=11.0
    $ make CC=cc FC=gfortran
```
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

### Android 
See [this page](https://github.com/xianyi/OpenBLAS/wiki/How-to-build-OpenBLAS-for-Android)

### MIPS
See [this page](https://github.com/xianyi/OpenBLAS/wiki/Faq#mips)