Werner found a lot of bugs about lapack testing. In 0.2.9 version, we used a work-around,e.g. fallback to the old kernel, replacing the optimized kernel with reference kernel. We must fix them in future release.

* ~~Nehalem dgemm kernel. Fallback to core2 kernel.~~
* ~~Sandy bridge sgemm kernel. Fallback to  Nehalem kernels.~~
* ~~Sandy bridge cgemm kernel. Fallback to  Nehalem kernels.~~
* sbmv, zsbmv, smp bug (interface/sbmv.c zsbmv.c). Now, it is only sequential.
* zhbmv, smp bug. 
* ~~scal, zscal x86/x86_64 kernels. Fallback to C implementation.~~
* ~~potri, zpotri. Fallback to LAPACK reference implementation.~~
* ~~lauu2, zlauu2.  Fallback to LAPACK reference implementation.~~
* ~~trtri, ztrtri. Fallback to LAPACK reference implementation.~~
* ~~lauum, zlauum. Fallback to LAPACK reference implementation.~~
* ~~trti2, ztrti2. Fallback to LAPACK reference implementation.~~
* Disable SMP in ger.c and zger.c