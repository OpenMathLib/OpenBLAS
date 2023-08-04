# TODO

Werner found a lot of bugs about lapack testing. In 0.2.9 version, we used a work-around,e.g. fallback to the old kernel, replacing the optimized kernel with reference kernel. We must fix them in future release.

* <s>Nehalem dgemm kernel. Fallback to core2 kernel.</s>
* <s>Sandy bridge sgemm kernel. Fallback to  Nehalem kernels.</s>
* <s>Sandy bridge cgemm kernel. Fallback to  Nehalem kernels.</s>
* sbmv, zsbmv, smp bug (interface/sbmv.c zsbmv.c). Now, it is only sequential.
* zhbmv, smp bug.
* <s>scal, zscal x86/x86_64 kernels. Fallback to C implementation.</s>
* <s>potri, zpotri. Fallback to LAPACK reference implementation.</s>
* <s>lauu2, zlauu2.  Fallback to LAPACK reference implementation.</s>
* <s>trtri, ztrtri. Fallback to LAPACK reference implementation.</s>
* <s>lauum, zlauum. Fallback to LAPACK reference implementation.</s>
* <s>trti2, ztrti2. Fallback to LAPACK reference implementation.</s>
* Disable SMP in ger.c and zger.c
