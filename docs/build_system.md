This page describes the Make-based build, which is the default/authoritative
build method. Note that the OpenBLAS repository also supports building with
CMake (not described here) - that generally works and is tested, however there
may be small differences between the Make and CMake builds.

!!! warning
    This page is made by someone who is not the developer and should not be considered as an official documentation of the build system. For getting the full picture, it is best to read the Makefiles and understand them yourself.

## Makefile dep graph

```
Makefile                                                        
|                                                               
|-----  Makefile.system # !!! this is included by many of the Makefiles in the subdirectories !!!
|       |
|       |=====  Makefile.prebuild # This is triggered (not included) once by Makefile.system 
|       |       |                 # and runs before any of the actual library code is built.
|       |       |                 # (builds and runs the "getarch" tool for cpu identification,
|       |       |                 # runs the compiler detection scripts c_check and f_check) 
|       |       |
|       |       -----  (Makefile.conf) [ either this or Makefile_kernel.conf is generated ] 
|       |       |                            { Makefile.system#L243 }
|       |       -----  (Makefile_kernel.conf) [ temporary Makefile.conf during DYNAMIC_ARCH builds ]
|       |
|       |-----  Makefile.rule # defaults for build options that can be given on the make command line
|       |
|       |-----  Makefile.$(ARCH) # architecture-specific compiler options and OpenBLAS buffer size values
|
|~~~~~ exports/
|
|~~~~~ test/
|
|~~~~~ utest/  
|
|~~~~~ ctest/
|
|~~~~~ cpp_thread_test/
|
|~~~~~ kernel/
|
|~~~~~ ${SUBDIRS}
|
|~~~~~ ${BLASDIRS}
|
|~~~~~ ${NETLIB_LAPACK_DIR}{,/timing,/testing/{EIG,LIN}}
|
|~~~~~ relapack/
```

## Important Variables

Most of the tunable variables are found in [Makefile.rule](https://github.com/xianyi/OpenBLAS/blob/develop/Makefile.rule), along with their detailed descriptions.<br/>
Most of the variables are detected automatically in [Makefile.prebuild](https://github.com/xianyi/OpenBLAS/blob/develop/Makefile.prebuild), if they are not set in the environment.

### CPU related
```
ARCH         - Target architecture (eg. x86_64)
TARGET       - Target CPU architecture, in case of DYNAMIC_ARCH=1 means library will not be usable on less capable CPUs
TARGET_CORE  - TARGET_CORE will override TARGET internally during each cpu-specific cycle of the build for DYNAMIC_ARCH
DYNAMIC_ARCH - For building library for multiple TARGETs (does not lose any optimizations, but increases library size)
DYNAMIC_LIST - optional user-provided subset of the DYNAMIC_CORE list in Makefile.system
```

### Toolchain related
```
CC                 - TARGET C compiler used for compilation (can be cross-toolchains)
FC                 - TARGET Fortran compiler used for compilation (can be cross-toolchains, set NOFORTRAN=1 if used cross-toolchain has no fortran compiler)
AR, AS, LD, RANLIB - TARGET toolchain helpers used for compilation (can be cross-toolchains)

HOSTCC             - compiler of build machine, needed to create proper config files for target architecture
HOST_CFLAGS        - flags for build machine compiler
```

### Library related
```
BINARY          - 32/64 bit library

BUILD_SHARED    - Create shared library
BUILD_STATIC    - Create static library

QUAD_PRECISION  - enable support for IEEE quad precision [ largely unimplemented leftover from GotoBLAS, do not use ]
EXPRECISION     - Obsolete option to use float80 of SSE on BSD-like systems
INTERFACE64     - Build with 64bit integer representations to support large array index values [ incompatible with standard API ]

BUILD_SINGLE    - build the single-precision real functions of BLAS [and optionally LAPACK] 
BUILD_DOUBLE    - build the double-precision real functions
BUILD_COMPLEX   - build the single-precision complex functions
BUILD_COMPLEX16 - build the double-precision complex functions
(all four types are included in the build by default when none was specifically selected)

BUILD_BFLOAT16  - build the "half precision brainfloat" real functions 
 
USE_THREAD      - Use a multithreading backend (default to pthread)
USE_LOCKING     - implement locking for thread safety even when USE_THREAD is not set (so that the singlethreaded library can
                  safely be called from multithreaded programs)
USE_OPENMP      - Use OpenMP as multithreading backend
NUM_THREADS     - define this to the maximum number of parallel threads you expect to need (defaults to the number of cores in the build cpu)
NUM_PARALLEL    - define this to the number of OpenMP instances that your code may use for parallel calls into OpenBLAS (default 1,see below)

```


OpenBLAS uses a fixed set of memory buffers internally, used for communicating
and compiling partial results from individual threads. For efficiency, the
management array structure for these buffers is sized at build time - this
makes it necessary to know in advance how many threads need to be supported on
the target system(s).

With OpenMP, there is an additional level of complexity as there may be calls
originating from a parallel region in the calling program. If OpenBLAS gets
called from a single parallel region, it runs single-threaded automatically to
avoid overloading the system by fanning out its own set of threads. In the case
that an OpenMP program makes multiple calls from independent regions or
instances in parallel, this default serialization is not sufficient as the
additional caller(s) would compete for the original set of buffers already in
use by the first call. So if multiple OpenMP runtimes call into OpenBLAS at the
same time, then only one of them will be able to make progress while all the
rest of them spin-wait for the one available buffer. Setting `NUM_PARALLEL` to
the upper bound on the number of OpenMP runtimes that you can have in a process
ensures that there are a sufficient number of buffer sets available.
