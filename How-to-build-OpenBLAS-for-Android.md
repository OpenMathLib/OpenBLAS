## Prerequisites

In addition to the Android NDK, you will need a C compiler on the build host as this is currently
required by the OpenBLAS build environment.


# Building with android NDK using clang compiler
Around version 11 Android NDKs stopped supporting gcc, so you would need to use clang to compile OpenBLAS. clang is supported from OpenBLAS 0.2.20 version onwards. See below sections on how to build with clang for ARMV7 and ARMV8 targets. The same basic principles as described below for ARMV8 should also apply to building an x86 or x86_64 version (substitute something like NEHALEM for the target instead of ARMV8 and replace all the aarch64 in the toolchain paths obviously)
"Historic" notes:
Since version 19 the default toolchain is provided as a standalone toolchain, so building one yourself following [building a standalone toolchain](http://developer.android.com/ndk/guides/standalone_toolchain.html) should no longer be necessary.
If you want to use static linking with an old NDK version older than about r17, you need to choose an API level below 23 currently due to NDK bug 272 (https://github.com/android-ndk/ndk/issues/272 , the libc.a lacks a definition of stderr) that will probably be fixed in r17 of the NDK.

## Build ARMV7 with clang
```
# Set path to ndk-bundle
export NDK_BUNDLE_DIR=/path/to/ndk-bundle

# Set the PATH to contain paths to clang and arm-linux-androideabi-* utilities
export PATH=${NDK_BUNDLE_DIR}/toolchains/arm-linux-androideabi-4.9/prebuilt/linux-x86_64/bin:${NDK_BUNDLE_DIR}/toolchains/llvm/prebuilt/linux-x86_64/bin:$PATH

# Set LDFLAGS so that the linker finds the appropriate libgcc
export LDFLAGS="-L${NDK_BUNDLE_DIR}/toolchains/arm-linux-androideabi-4.9/prebuilt/linux-x86_64/lib/gcc/arm-linux-androideabi/4.9.x"

# Set the clang cross compile flags
export CLANG_FLAGS="-target arm-linux-androideabi -marm -mfpu=vfp -mfloat-abi=softfp --sysroot ${NDK_BUNDLE_DIR}/platforms/android-23/arch-arm -gcc-toolchain ${NDK_BUNDLE_DIR}/toolchains/arm-linux-androideabi-4.9/prebuilt/linux-x86_64/"

#OpenBLAS Compile
make TARGET=ARMV7 ONLY_CBLAS=1 AR=ar CC="clang ${CLANG_FLAGS}" HOSTCC=gcc ARM_SOFTFP_ABI=1 -j4
```
On a Mac, it may also be necessary to give the complete path to the `ar` utility in the make command above, like so:
```
AR=${NDK_BUNDLE_DIR}/toolchains/arm-linux-androideabi-4.9/prebuilt/darwin-x86_64/bin/arm-linux-androideabi-gcc-ar
```
otherwise you may get a linker error complaining about a "malformed archive header name at 8" when the native OSX ar command was invoked instead. Note that with recent NDK versions, the AR tool may be named `llvm-ar` rather than what is assumed above. 
 
## Build ARMV8 with clang
```
# Set path to ndk-bundle
export NDK_BUNDLE_DIR=/path/to/ndk-bundle/

# Export PATH to contain directories of clang and aarch64-linux-android-* utilities
export PATH=${NDK_BUNDLE_DIR}/toolchains/aarch64-linux-android-4.9/prebuilt/linux-x86_64/bin/:${NDK_BUNDLE_DIR}/toolchains/llvm/prebuilt/linux-x86_64/bin:$PATH

# Setup LDFLAGS so that loader can find libgcc and pass -lm for sqrt
export LDFLAGS="-L${NDK_BUNDLE_DIR}/toolchains/aarch64-linux-android-4.9/prebuilt/linux-x86_64/lib/gcc/aarch64-linux-android/4.9.x -lm"

# Setup the clang cross compile options
export CLANG_FLAGS="-target aarch64-linux-android --sysroot ${NDK_BUNDLE_DIR}/platforms/android-23/arch-arm64 -gcc-toolchain ${NDK_BUNDLE_DIR}/toolchains/aarch64-linux-android-4.9/prebuilt/linux-x86_64/"

# Compile
make TARGET=ARMV8 ONLY_CBLAS=1 AR=ar CC="clang ${CLANG_FLAGS}" HOSTCC=gcc -j4
```
Note: Using TARGET=CORTEXA57 in place of ARMV8 will pick up better optimized routines. Implementations for CORTEXA57 target is compatible with all other armv8 targets.

Note: For NDK 23b, something as simple as 
```
export PATH=/opt/android-ndk-r23b/toolchains/llvm/prebuilt/linux-x86_64/bin/:$PATH
make HOSTCC=gcc CC=/opt/android-ndk-r23b/toolchains/llvm/prebuilt/linux-x86_64/bin/aarch64-linux-android31-clang ONLY_CBLAS=1 TARGET=ARMV8
```
appears to be sufficient on Linux. On OSX, setting AR to the ar provided in the "bin" path of the NDK (probably llvm-ar) is also necessary.

## Alternative script which was tested on OSX with NDK(21.3.6528147)
This script will build openblas for 3 architecture (ARMV7,ARMV8,X86) and put them with `sudo make install` to `/opt/OpenBLAS/lib`.
Of course you can also copy only the section that is of interest to you - also notice that the AR= line may need adapting to the
name of the ar tool provided in your $TOOLCHAIN/bin - for example llvm-ar in some recent NDK versions.
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

# This will build for x86_64 
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
## Building OpenBLAS with very old gcc-based versions of the NDK, without Fortran

The prebuilt Android NDK toolchains do not include Fortran, hence parts like LAPACK cannot be built. You can still build OpenBLAS without it. For instructions on how to build OpenBLAS with Fortran, see the [next section](#building-openblas-with-fortran).

To use easily the prebuilt toolchains, follow [building a standalone toolchain](http://developer.android.com/ndk/guides/standalone_toolchain.html) for your desired architecture.
This would be `arm-linux-androideabi-gcc-4.9` for ARMV7 and `aarch64-linux-android-gcc-4.9` for ARMV8.

You can build OpenBLAS (0.2.19 and earlier) with:
```
# Add the toolchain to your path
export PATH=/path/to/standalone-toolchain/bin:$PATH

# Build without Fortran for ARMV7
make TARGET=ARMV7 HOSTCC=gcc CC=arm-linux-androideabi-gcc NOFORTRAN=1 libs
# Build without Fortran for ARMV8
make TARGET=ARMV8 BINARY=64 HOSTCC=gcc CC=aarch64-linux-android-gcc NOFORTRAN=1 libs
```

Since we are cross-compiling, we make the `libs` recipe, not `all`. Otherwise you will get errors when trying to link/run tests as versions up to and including 0.2.19 cannot build a shared library for Android. 

From 0.2.20 on, you should leave off the "libs" to get a full build, and you may want to use the softfp ABI instead of the deprecated hardfp one on ARMV7 so you would use 
```
# Add the toolchain to your path
export PATH=/path/to/standalone-toolchain/bin:$PATH

# Build without Fortran for ARMV7
make TARGET=ARMV7 ARM_SOFTFP_ABI=1 HOSTCC=gcc CC=arm-linux-androideabi-gcc NOFORTRAN=1
# Build without Fortran for ARMV8
make TARGET=ARMV8 BINARY=64 HOSTCC=gcc CC=aarch64-linux-android-gcc NOFORTRAN=1
```

If you get an error about stdio.h not being found, you need to specify your sysroot in the CFLAGS argument to `make` like
```CFLAGS=--sysroot=$NDK/platforms/android-16/arch-arm```
When you are done, install OpenBLAS into the desired directory. Be sure to also use all command line options
here that you specified for building, otherwise errors may occur as it tries to install things you did not build:
```
make PREFIX=/path/to/install-dir TARGET=... install
```

## Building OpenBLAS with Fortran

Instructions on how to build the GNU toolchains with Fortran can be found [here](https://github.com/buffer51/android-gfortran). The [Releases section](https://github.com/buffer51/android-gfortran/releases) provides prebuilt versions, use the standalone one.

You can build OpenBLAS with:
```
# Add the toolchain to your path
export PATH=/path/to/standalone-toolchain-with-fortran/bin:$PATH

# Build with Fortran for ARMV7
make TARGET=ARMV7 HOSTCC=gcc CC=arm-linux-androideabi-gcc FC=arm-linux-androideabi-gfortran libs
# Build with LAPACK for ARMV8
make TARGET=ARMV8 BINARY=64 HOSTCC=gcc CC=aarch64-linux-android-gcc FC=aarch64-linux-android-gfortran libs
```

As mentioned above you can leave off the `libs` argument here when building 0.2.20 and later, and you may want to add ARM_SOFTFP_ABI=1 when building for ARMV7.

## Linking OpenBLAS (0.2.19 and earlier) for ARMV7

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
