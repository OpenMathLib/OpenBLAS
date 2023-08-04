## Compile the library
### Normal compile
  * type `make` to detect the CPU automatically.
  or
  * type `make TARGET=xxx` to set target CPU, e.g. `make TARGET=NEHALEM`. The full target list is in file TargetList.txt.

### Cross compile
Please set `CC` and `FC` with the cross toolchains. Then, set `HOSTCC` with your host C compiler. At last, set `TARGET` explicitly.

Examples:

* On x86 box, compile the library for ARM Cortex-A9 linux.

Install only gnueabihf versions. Please check https://github.com/xianyi/OpenBLAS/issues/936#issuecomment-237596847

    make CC=arm-linux-gnueabihf-gcc FC=arm-linux-gnueabihf-gfortran HOSTCC=gcc TARGET=CORTEXA9

* On X86 box, compile this library for loongson3a CPU.

```
make BINARY=64 CC=mips64el-unknown-linux-gnu-gcc FC=mips64el-unknown-linux-gnu-gfortran HOSTCC=gcc TARGET=LOONGSON3A
```

* On X86 box, compile this library for loongson3a CPU with loongcc (based on Open64) compiler.

```
make CC=loongcc FC=loongf95 HOSTCC=gcc TARGET=LOONGSON3A CROSS=1 CROSS_SUFFIX=mips64el-st-linux-gnu-   NO_LAPACKE=1 NO_SHARED=1 BINARY=32
```

### Debug version

    make DEBUG=1

### Install to the directory (optional)

Example:

    make install PREFIX=your_installation_directory

The default directory is /opt/OpenBLAS. Note that any flags passed to `make` during build should also be passed to `make install` to circumvent any install errors, i.e. some headers not being copied over correctly.

For more information, please read [Installation Guide](install.md).

## Link the library

* Link shared library

```
gcc -o test test.c -I/your_path/OpenBLAS/include/ -L/your_path/OpenBLAS/lib -Wl,-rpath,/your_path/OpenBLAS/lib -lopenblas
```

The `-Wl,-rpath,/your_path/OpenBLAS/lib` option to linker can be omitted if you ran `ldconfig` to update linker cache, put `/your_path/OpenBLAS/lib` in `/etc/ld.so.conf` or a file in `/etc/ld.so.conf.d`, or installed OpenBLAS in a location part of `ld.so` default search path. Otherwise, linking at runtime will fail.

If the library is multithreaded, please add `-lpthread`. If the library contains LAPACK functions, please add `-lgfortran` or other Fortran libs, although if you only make calls to LAPACKE routines, i.e. your code has `#include "lapacke.h"` and makes calls to methods like `LAPACKE_dgeqrf`, `-lgfortran` is not needed.

* Link static library

```
gcc -o test test.c /your/path/libopenblas.a
```

You can download `test.c` from https://gist.github.com/xianyi/5780018 

## Code examples

### Call CBLAS interface
This example shows calling cblas_dgemm in C. https://gist.github.com/xianyi/6930656
```c
#include <cblas.h>
#include <stdio.h>

void main()
{
  int i=0;
  double A[6] = {1.0,2.0,1.0,-3.0,4.0,-1.0};         
  double B[6] = {1.0,2.0,1.0,-3.0,4.0,-1.0};  
  double C[9] = {.5,.5,.5,.5,.5,.5,.5,.5,.5}; 
  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,3,3,2,1,A, 3, B, 3,2,C,3);

  for(i=0; i<9; i++)
    printf("%lf ", C[i]);
  printf("\n");
}
```

```
gcc -o test_cblas_open test_cblas_dgemm.c -I /your_path/OpenBLAS/include/ -L/your_path/OpenBLAS/lib -lopenblas -lpthread -lgfortran
```

### Call BLAS Fortran interface

This example shows calling dgemm Fortran interface in C. https://gist.github.com/xianyi/5780018

```c
#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"
#include "time.h"

extern void dgemm_(char*, char*, int*, int*,int*, double*, double*, int*, double*, int*, double*, double*, int*);

int main(int argc, char* argv[])
{
  int i;
  printf("test!\n");
  if(argc<4){
    printf("Input Error\n");
    return 1;
  }

  int m = atoi(argv[1]);
  int n = atoi(argv[2]);
  int k = atoi(argv[3]);
  int sizeofa = m * k;
  int sizeofb = k * n;
  int sizeofc = m * n;
  char ta = 'N';
  char tb = 'N';
  double alpha = 1.2;
  double beta = 0.001;

  struct timeval start,finish;
  double duration;

  double* A = (double*)malloc(sizeof(double) * sizeofa);
  double* B = (double*)malloc(sizeof(double) * sizeofb);
  double* C = (double*)malloc(sizeof(double) * sizeofc);

  srand((unsigned)time(NULL));

  for (i=0; i<sizeofa; i++)
    A[i] = i%3+1;//(rand()%100)/10.0;

  for (i=0; i<sizeofb; i++)
    B[i] = i%3+1;//(rand()%100)/10.0;

  for (i=0; i<sizeofc; i++)
    C[i] = i%3+1;//(rand()%100)/10.0;
  //#if 0
  printf("m=%d,n=%d,k=%d,alpha=%lf,beta=%lf,sizeofc=%d\n",m,n,k,alpha,beta,sizeofc);
  gettimeofday(&start, NULL);
  dgemm_(&ta, &tb, &m, &n, &k, &alpha, A, &m, B, &k, &beta, C, &m);
  gettimeofday(&finish, NULL);

  duration = ((double)(finish.tv_sec-start.tv_sec)*1000000 + (double)(finish.tv_usec-start.tv_usec)) / 1000000;
  double gflops = 2.0 * m *n*k;
  gflops = gflops/duration*1.0e-6;

  FILE *fp;
  fp = fopen("timeDGEMM.txt", "a");
  fprintf(fp, "%dx%dx%d\t%lf s\t%lf MFLOPS\n", m, n, k, duration, gflops);
  fclose(fp);

  free(A);
  free(B);
  free(C);
  return 0;
}
```

```
gcc -o time_dgemm time_dgemm.c /your/path/libopenblas.a
./time_dgemm <m> <n> <k>
```

## Troubleshooting

* Please read [Faq](faq.md) at first.
* Please use gcc version 4.6 and above to compile Sandy Bridge AVX kernels on Linux/MingW/BSD.
* Please use Clang version 3.1 and above to compile the library on Sandy Bridge microarchitecture. The Clang 3.0 will generate the wrong AVX binary code.
* The number of CPUs/Cores should less than or equal to 256. On Linux x86_64(amd64), there is experimental support for up to 1024 CPUs/Cores and 128 numa nodes if you build the library with BIGNUMA=1.
* OpenBLAS does not set processor affinity by default. On Linux, you can enable processor affinity by commenting the line NO_AFFINITY=1 in Makefile.rule. But this may cause [the conflict with R parallel](https://stat.ethz.ch/pipermail/r-sig-hpc/2012-April/001348.html).
* On Loongson 3A. make test would be failed because of pthread_create error. The error code is EAGAIN. However, it will be OK when you run the same testcase on shell.

## BLAS reference manual

If you want to understand every BLAS function and definition, please read [Intel MKL reference manual](https://software.intel.com/en-us/intel-mkl/documentation) or [netlib.org](http://netlib.org/blas/)

Here are [OpenBLAS extension functions](extensions.md)
