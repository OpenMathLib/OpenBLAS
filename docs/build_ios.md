# How to build OpenBLAS for iPhone/iOS

As none of the current developers uses iOS, the following instructions are what was found to work in our Azure CI setup, but as far as we know this builds a fully working OpenBLAS for this platform.

Go to the directory where you unpacked OpenBLAS,and enter the following commands:
```
     CC=/Applications/Xcode_12.4.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/clang

CFLAGS= -O2 -Wno-macro-redefined -isysroot /Applications/Xcode_12.4.app/Contents/Developer/Platforms/iPhoneOS.platform/Developer/SDKs/iPhoneOS14.4.sdk -arch arm64 -miphoneos-version-min=10.0

make TARGET=ARMV8 DYNAMIC_ARCH=1 NUM_THREADS=32 HOSTCC=clang NOFORTRAN=1
```
Adjust MIN_IOS_VERSION as necessary for your installation, e.g. change the version number
to the minimum iOS version you want to target and execute this file to build the library.