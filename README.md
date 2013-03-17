# OpenBLAS

## Introduction
OpenBLAS is an optimized BLAS library based on GotoBLAS2 1.13 BSD version. OpenBLAS is an open source project supported by Lab of Parallel Software and Computational Science, ISCAS <http://www.rdcps.ac.cn>.

Please read the documents on OpenBLAS wiki pages <http://github.com/xianyi/OpenBLAS/wiki>.

## Installation
Download from project homepage. http://xianyi.github.com/OpenBLAS/

Or, check out codes from git://github.com/xianyi/OpenBLAS.git
### Normal compile
  * type "make" to detect the CPU automatically.
  or
  * type "make TARGET=xxx" to set target CPU, e.g. "make TARGET=NEHALEM". The full target list is in file TargetList.txt.

### Cross compile
Please set CC and FC with the cross toolchains. Then, set HOSTCC with your host C compiler. At last, set TARGET explicitly.

Examples:

On X86 box, compile this library for loongson3a CPU.

    make BINARY=64 CC=mips64el-unknown-linux-gnu-gcc FC=mips64el-unknown-linux-gnu-gfortran HOSTCC=gcc TARGET=LOONGSON3A

### Debug version

    make DEBUG=1

### Install to the directory (Optional)

Example:

    make install PREFIX=your_installation_directory

The default directory is /opt/OpenBLAS

## Support CPU & OS
Please read GotoBLAS_01Readme.txt

### Additional support CPU:

#### x86/x86-64:
- **Intel Xeon 56xx (Westmere)**: Used GotoBLAS2 Nehalem codes.
- **Intel Sandy Bridge**: Optimized Level-3 BLAS with AVX on x86-64.
- **AMD Bobcat**: Used GotoBLAS2 Barcelona codes.
- **AMD Bulldozer**: x86-64 S/DGEMM AVX kernels. (Thank Werner Saar)

#### MIPS64:
- **ICT Loongson 3A**: Optimized Level-3 BLAS and the part of Level-1,2.
- **ICT Loongson 3B**: Experimental

### Support OS:
- **GNU/Linux**
- **MingWin/Windows**: Please read <https://github.com/xianyi/OpenBLAS/wiki/How-to-use-OpenBLAS-in-Microsoft-Visual-Studio>.
- **Darwin/Mac OS X**: Experimental. Although GotoBLAS2 supports Darwin, we are the beginner on Mac OS X.
- **FreeBSD**: Supportted by community. We didn't test the library on this OS.

## Usages
Link with libopenblas.a or -lopenblas for shared library.

### Set the number of threads with environment variables. 

Examples:

    export OPENBLAS_NUM_THREADS=4

 or

    export GOTO_NUM_THREADS=4

 or 

    export OMP_NUM_THREADS=4

The priorities are OPENBLAS_NUM_THREADS > GOTO_NUM_THREADS > OMP_NUM_THREADS.

If you compile this lib with USE_OPENMP=1, you should set OMP_NUM_THREADS environment variable. OpenBLAS ignores OPENBLAS_NUM_THREADS and GOTO_NUM_THREADS with USE_OPENMP=1.

### Set the number of threads on runtime. 

We provided the below functions to controll the number of threads on runtime.

    void goto_set_num_threads(int num_threads);

    void openblas_set_num_threads(int num_threads);

If you compile this lib with USE_OPENMP=1, you should use the above functions, too.

## Report Bugs
Please add a issue in https://github.com/xianyi/OpenBLAS/issues

## Contact
OpenBLAS users mailing list: http://list.rdcps.ac.cn/mailman/listinfo/openblas

## ChangeLog
Please see Changelog.txt to obtain the differences between GotoBLAS2 1.13 BSD version.

## Troubleshooting
* Please read [Faq](https://github.com/xianyi/OpenBLAS/wiki/Faq) at first.
* Please use gcc version 4.6 and above to compile Sandy Bridge AVX kernels on Linux/MingW/BSD.
* Please use Clang version 3.1 and above to compile the library on Sandy Bridge microarchitecture. The Clang 3.0 will generate the wrong AVX binary code.
* The number of CPUs/Cores should less than or equal to 256. 
* On Linux, OpenBLAS sets the processor affinity by default. This may cause [the conflict with R parallel](https://stat.ethz.ch/pipermail/r-sig-hpc/2012-April/001348.html). You can build the library with NO_AFFINITY=1.
* On Loongson 3A. make test would be failed because of pthread_create error. The error code is EAGAIN. However, it will be OK when you run the same testcase on shell. 

## Specification of Git Branches
We used the git branching model in this article (http://nvie.com/posts/a-successful-git-branching-model/). 
Now, there are 4 branches in github.com.
  * The master branch. This a main branch to reflect a production-ready state.
  * The develop branch. This a main branch to reflect a state with the latest delivered development changes for the next release.
  * The loongson3a branch. This is a feature branch. We develop Loongson3A codes on this branch. We will merge this feature to develop branch in future.
  * The gh-pages branch. This is for web pages
