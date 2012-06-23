# OpenBLAS

## Introduction
OpenBLAS is an optimized BLAS library based on GotoBLAS2 1.13 BSD version. OpenBLAS is an open source project supported by Lab of Parallel Software and Computational Science, ISCAS <http://www.rdcps.ac.cn>.

Please read the documents on OpenBLAS wiki pages <http://github.com/xianyi/OpenBLAS/wiki>.

## Intallation
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

### Intall to the directory (Optional)

Example:

    make install PREFIX=your_installation_directory

The default directory is /opt/OpenBLAS

## Support CPU & OS
Please read GotoBLAS_01Readme.txt

### Additional support CPU:

#### x86/x86-64:
* Intel Xeon 56xx (Westmere). Used GotoBLAS2 Nehalem codes.
* Intel Sandy Bridge. Optimized Level-3 BLAS with AVX on x86-64.
* AMD Bobcat. Used GotoBLAS2 Barcelona codes.
#### MIPS64:
* ICT Loongson 3A. Optimized Level-3 BLAS and the part of Level-1,2.
* ICT Loongson 3B (Experimental)

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

### Set the number of threads with calling functions. 

Examples:

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
* Please use gcc version 4.6 and above to compile Sandy Bridge AVX kernels on Linux/MingW/BSD.
* Please use Clang version 3.1 and above to compile the library on Sandy Bridge microarchitecture. The Clang 3.0 will generate the wrong AVX binary code.
* The number of CPUs/Cores should less than or equal to 256. 
* On Loongson 3A. make test would be failed because of pthread_create error. The error code is EAGAIN. However, it will be OK when you run the same testcase on shell. 

## Specification of Git Branches
We used the git branching model in this article (http://nvie.com/posts/a-successful-git-branching-model/). 
Now, there are 4 branches in github.com.
  * The master branch. This a main branch to reflect a production-ready state.
  * The develop branch. This a main branch to reflect a state with the latest delivered development changes for the next release.
  * The loongson3a branch. This is a feature branch. We develop Loongson3A codes on this branch. We will merge this feature to develop branch in future.
  * The gh-pages branch. This is for web pages
