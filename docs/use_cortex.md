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

Note: if you are developing for an embedded platform, it is your responsibility to make sure that the device
has sufficient memory for malloc calls. [Libmemory][1] provides one implementation of malloc for embedded
platforms.


[1]: https://github.com/embeddedartistry/libmemory

