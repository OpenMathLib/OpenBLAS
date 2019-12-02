/* %0 = "+r"(a_pointer), %1 = "+r"(b_pointer), %2 = "+r"(c_pointer), %3 = "+r"(ldc_in_bytes), %4 for k_count, %5 for c_store */
/* r11 = m(const), r12 = k << 5(const), r13 = k(const), r14 = b_head_pos(const), r15 = %1 + 3r12 */

#include "common.h"
#include <stdint.h>

//recommended settings: GEMM_Q=168, GEMM_P=384

/* m = 8 */ /* zmm8-zmm31 for accumulators, zmm3-zmm7 for temporary use, zmm0 for alpha, zmm1-zmm2 for control words of permutations */
#define KERNEL_k1m8n1 \
    "vmovupd (%0),%%zmm4; addq $64,%0;"\
    "vbroadcastsd (%1),%%zmm6; vfmadd231pd %%zmm4,%%zmm6,%%zmm8;"\
    "addq $8,%1;"
#define KERNEL_h_k1m8n2 \
    "vmovddup (%0),%%zmm4; vmovddup 8(%0),%%zmm5; prefetcht0 512(%0); addq $64,%0;"\
    "vbroadcastf32x4 (%1),%%zmm6; vfmadd231pd %%zmm4,%%zmm6,%%zmm8; vfmadd231pd %%zmm5,%%zmm6,%%zmm9;"
#define KERNEL_k1m8n2 KERNEL_h_k1m8n2 "addq $16,%1;"
#define KERNEL_h_k1m8n4 KERNEL_h_k1m8n2 "vbroadcastf32x4 16(%1),%%zmm7; vfmadd231pd %%zmm4,%%zmm7,%%zmm10; vfmadd231pd %%zmm5,%%zmm7,%%zmm11;"
#define KERNEL_k1m8n4 KERNEL_h_k1m8n4 "addq $32,%1;"
#define unit_kernel_k1m8n4(c1,c2,c3,c4, ...) \
    "vbroadcastf32x4  ("#__VA_ARGS__"),%%zmm6; vfmadd231pd %%zmm4,%%zmm6,"#c1"; vfmadd231pd %%zmm5,%%zmm6,"#c2";"\
    "vbroadcastf32x4 16("#__VA_ARGS__"),%%zmm7; vfmadd231pd %%zmm4,%%zmm7,"#c3"; vfmadd231pd %%zmm5,%%zmm7,"#c4";"
#define KERNEL_h_k1m8n8 KERNEL_h_k1m8n4 unit_kernel_k1m8n4(%%zmm12,%%zmm13,%%zmm14,%%zmm15,%1,%%r12,1)
#define KERNEL_k1m8n8 KERNEL_h_k1m8n8 "addq $32,%1;"
#define KERNEL_h_k1m8n12 KERNEL_h_k1m8n8 unit_kernel_k1m8n4(%%zmm16,%%zmm17,%%zmm18,%%zmm19,%1,%%r12,2)
#define KERNEL_k1m8n12 KERNEL_h_k1m8n12 "addq $32,%1;"
#define KERNEL_h_k1m8n16 KERNEL_k1m8n12 unit_kernel_k1m8n4(%%zmm20,%%zmm21,%%zmm22,%%zmm23,%%r15)
#define KERNEL_k1m8n16 KERNEL_h_k1m8n16 "addq $32,%%r15;"
#define KERNEL_h_k1m8n20 KERNEL_h_k1m8n16 unit_kernel_k1m8n4(%%zmm24,%%zmm25,%%zmm26,%%zmm27,%%r15,%%r12,1)
#define KERNEL_k1m8n20 KERNEL_h_k1m8n20 "addq $32,%%r15;"
#define KERNEL_h_k1m8n24 KERNEL_h_k1m8n20 unit_kernel_k1m8n4(%%zmm28,%%zmm29,%%zmm30,%%zmm31,%%r15,%%r12,2)
#define KERNEL_k1m8n24 KERNEL_h_k1m8n24 "addq $32,%%r15;"
#define INIT_m8n1 "vpxorq %%zmm8,%%zmm8,%%zmm8;"
#define INIT_m8n2 INIT_m8n1 "vpxorq %%zmm9,%%zmm9,%%zmm9;"
#define INIT_m8n4 INIT_m8n2 "vpxorq %%zmm10,%%zmm10,%%zmm10;vpxorq %%zmm11,%%zmm11,%%zmm11;"
#define unit_init_m8n4(c1,c2,c3,c4) \
    "vpxorq "#c1","#c1","#c1";vpxorq "#c2","#c2","#c2";vpxorq "#c3","#c3","#c3";vpxorq "#c4","#c4","#c4";"
#define INIT_m8n8 INIT_m8n4 unit_init_m8n4(%%zmm12,%%zmm13,%%zmm14,%%zmm15)
#define INIT_m8n12 INIT_m8n8 unit_init_m8n4(%%zmm16,%%zmm17,%%zmm18,%%zmm19)
#define INIT_m8n16 INIT_m8n12 unit_init_m8n4(%%zmm20,%%zmm21,%%zmm22,%%zmm23)
#define INIT_m8n20 INIT_m8n16 unit_init_m8n4(%%zmm24,%%zmm25,%%zmm26,%%zmm27)
#define INIT_m8n24 INIT_m8n20 unit_init_m8n4(%%zmm28,%%zmm29,%%zmm30,%%zmm31)
#define SAVE_h_m8n1 \
    "vunpcklpd %%zmm8,%%zmm8,%%zmm3; vunpckhpd %%zmm8,%%zmm8,%%zmm4; vmovapd %%zmm3,%%zmm5;"\
    "vpermt2pd %%zmm4,%%zmm1,%%zmm3; vpermt2pd %%zmm5,%%zmm2,%%zmm4;"\
    "vfmadd213pd (%2),%%zmm0,%%zmm3; vmovupd %%zmm3,(%2);"\
    "vfmadd213pd 64(%2),%%zmm0,%%zmm4; vmovupd %%zmm4,64(%2);"
#define unit_save_m8n2(c1,c2) \
    "vmovapd "#c1",%%zmm4; vpermt2pd "#c2",%%zmm1,"#c1"; vpermt2pd %%zmm4,%%zmm2,"#c2";"\
    "vunpcklpd "#c1","#c1",%%zmm4; vunpcklpd "#c2","#c2",%%zmm5; vunpckhpd "#c1","#c1",%%zmm6; vunpckhpd "#c2","#c2",%%zmm7;"\
    "vfmadd213pd (%5),%%zmm0,%%zmm4;vfmadd213pd 64(%5),%%zmm0,%%zmm5;vfmadd213pd (%5,%3,1),%%zmm0,%%zmm6;vfmadd213pd 64(%5,%3,1),%%zmm0,%%zmm7;"\
    "vmovupd %%zmm4,(%5); vmovupd %%zmm5,64(%5); vmovupd %%zmm6,(%5,%3,1); vmovupd %%zmm7,64(%5,%3,1); leaq (%5,%3,2),%5;"
#define SAVE_h_m8n2 "movq %2,%5;" unit_save_m8n2(%%zmm8,%%zmm9)
#define SAVE_h_m8n4  SAVE_h_m8n2  unit_save_m8n2(%%zmm10,%%zmm11)
#define SAVE_h_m8n8  SAVE_h_m8n4  unit_save_m8n2(%%zmm12,%%zmm13) unit_save_m8n2(%%zmm14,%%zmm15)
#define SAVE_h_m8n12 SAVE_h_m8n8  unit_save_m8n2(%%zmm16,%%zmm17) unit_save_m8n2(%%zmm18,%%zmm19)
#define SAVE_h_m8n16 SAVE_h_m8n12 unit_save_m8n2(%%zmm20,%%zmm21) unit_save_m8n2(%%zmm22,%%zmm23)
#define SAVE_h_m8n20 SAVE_h_m8n16 unit_save_m8n2(%%zmm24,%%zmm25) unit_save_m8n2(%%zmm26,%%zmm27)
#define SAVE_h_m8n24 SAVE_h_m8n20 unit_save_m8n2(%%zmm28,%%zmm29) unit_save_m8n2(%%zmm30,%%zmm31)
#define SAVE_m8(ndim) SAVE_h_m8n##ndim "addq $128,%2;"
#define COMPUTE_m8(ndim) \
    INIT_m8n##ndim\
    "movq %%r13,%4; movq %%r14,%1; leaq (%1,%%r12,2),%%r15; addq %%r12,%%r15; movq %2,%5;"\
    "cmpq $30,%4; jb "#ndim"008082f;"\
    #ndim"008081:\n\t"\
    KERNEL_k1m8n##ndim\
    KERNEL_k1m8n##ndim\
    KERNEL_k1m8n##ndim\
    "prefetcht1 (%5); prefetcht1 64(%5); prefetcht1 127(%5); addq %3,%5;"\
    KERNEL_k1m8n##ndim\
    KERNEL_k1m8n##ndim\
    KERNEL_k1m8n##ndim\
    "prefetcht1 (%8); addq $32,%8;"\
    "subq $6,%4; cmpq $30,%4; jnb "#ndim"008081b;"\
    "movq %2,%5;"\
    #ndim"008082:\n\t"\
    "testq %4,%4; jz "#ndim"008083f;"\
    "prefetcht0 (%5); prefetcht0 64(%5); prefetcht0 127(%5); addq %3,%5;"\
    KERNEL_k1m8n##ndim\
    "decq %4; jmp "#ndim"008082b;"\
    #ndim"008083:\n\t"\
    "prefetcht0 (%%r14); prefetcht0 64(%%r14);"\
    SAVE_m8(ndim)

/* m = 4 *//* ymm0 for alpha, ymm1-ymm3 for temporary use, ymm4-ymm15 for accumulators */
#define KERNEL_k1m4n1(b_addr) \
    "vmovupd (%0),%%ymm1; addq $32,%0;"\
    "vbroadcastsd ("#b_addr"),%%ymm2; vfmadd231pd %%ymm1,%%ymm2,%%ymm4;"\
    "addq $8,"#b_addr";"
#define KERNEL_h_k1m4n2(b_addr) \
    "vmovddup (%0),%%ymm1; vmovddup 8(%0),%%ymm2; addq $32,%0;"\
    "vbroadcastf128 ("#b_addr"),%%ymm3; vfmadd231pd %%ymm1,%%ymm3,%%ymm4; vfmadd231pd %%ymm2,%%ymm3,%%ymm5;"
#define KERNEL_k1m4n2(b_addr) KERNEL_h_k1m4n2(b_addr) "addq $16,"#b_addr";"
#define KERNEL_h_k1m4n4(b_addr) \
    KERNEL_h_k1m4n2(b_addr) "vbroadcastf128 16("#b_addr"),%%ymm3; vfmadd231pd %%ymm1,%%ymm3,%%ymm6; vfmadd231pd %%ymm2,%%ymm3,%%ymm7;"
#define KERNEL_k1m4n4(b_addr) KERNEL_h_k1m4n4(b_addr) "addq $32,"#b_addr";"
#define unit_kernel_k1m4n4(c1,c2,c3,c4,...) \
    "vbroadcastf128  ("#__VA_ARGS__"),%%ymm3; vfmadd231pd %%ymm1,%%ymm3,"#c1"; vfmadd231pd %%ymm2,%%ymm3,"#c2";"\
    "vbroadcastf128 16("#__VA_ARGS__"),%%ymm3; vfmadd231pd %%ymm1,%%ymm3,"#c3"; vfmadd231pd %%ymm2,%%ymm3,"#c4";"
#define KERNEL_h_k1m4n8(b_addr) KERNEL_h_k1m4n4(b_addr) unit_kernel_k1m4n4(%%ymm8,%%ymm9,%%ymm10,%%ymm11,b_addr,%%r12,1)
#define KERNEL_k1m4n8(b_addr) KERNEL_h_k1m4n8(b_addr) "addq $32,"#b_addr";"
#define KERNEL_h_k1m4n12(b_addr) KERNEL_h_k1m4n8(b_addr) unit_kernel_k1m4n4(%%ymm12,%%ymm13,%%ymm14,%%ymm15,b_addr,%%r12,2)
#define KERNEL_k1m4n12(b_addr) KERNEL_h_k1m4n12(b_addr) "addq $32,"#b_addr";"
#define INIT_m4n1 "vpxor %%ymm4,%%ymm4,%%ymm4;"
#define INIT_m4n2 INIT_m4n1 "vpxor %%ymm5,%%ymm5,%%ymm5;"
#define INIT_m4n4 INIT_m4n2 "vpxor %%ymm6,%%ymm6,%%ymm6;vpxor %%ymm7,%%ymm7,%%ymm7;"
#define unit_init_m4n4(c1,c2,c3,c4) \
    "vpxor "#c1","#c1","#c1";vpxor "#c2","#c2","#c2";vpxor "#c3","#c3","#c3";vpxor "#c4","#c4","#c4";"
#define INIT_m4n8  INIT_m4n4 unit_init_m4n4(%%ymm8,%%ymm9,%%ymm10,%%ymm11)
#define INIT_m4n12 INIT_m4n8 unit_init_m4n4(%%ymm12,%%ymm13,%%ymm14,%%ymm15)
#define SAVE_L_m4n1 \
    "vpermpd $216,%%ymm4,%%ymm3; vunpcklpd %%ymm3,%%ymm3,%%ymm1; vunpckhpd %%ymm3,%%ymm3,%%ymm2;"\
    "vfmadd213pd (%2),%%ymm0,%%ymm1; vfmadd213pd 32(%2),%%ymm0,%%ymm2; vmovupd %%ymm1,(%2); vmovupd %%ymm2,32(%2);"
#define unit_save_m4n2(c1,c2) \
    "vperm2f128 $2,"#c1","#c2",%%ymm2; vperm2f128 $19,"#c1","#c2","#c2"; vmovapd %%ymm2,"#c1";"\
    "vunpcklpd "#c1","#c1",%%ymm2; vunpcklpd "#c2","#c2",%%ymm3;"\
    "vfmadd213pd (%5),%%ymm0,%%ymm2; vfmadd213pd 32(%5),%%ymm0,%%ymm3; vmovupd %%ymm2,(%5); vmovupd %%ymm3,32(%5);"\
    "vunpckhpd "#c1","#c1",%%ymm2; vunpckhpd "#c2","#c2",%%ymm3;"\
    "vfmadd213pd (%5,%3,1),%%ymm0,%%ymm2; vfmadd213pd 32(%5,%3,1),%%ymm0,%%ymm3; vmovupd %%ymm2,(%5,%3,1); vmovupd %%ymm3,32(%5,%3,1);"\
    "leaq (%5,%3,2),%5;"
#define SAVE_L_m4n2 "movq %2,%5;" unit_save_m4n2(%%ymm4,%%ymm5)
#define SAVE_L_m4n4  SAVE_L_m4n2  unit_save_m4n2(%%ymm6,%%ymm7)
#define SAVE_L_m4n8  SAVE_L_m4n4  unit_save_m4n2(%%ymm8,%%ymm9)   unit_save_m4n2(%%ymm10,%%ymm11)
#define SAVE_L_m4n12 SAVE_L_m4n8  unit_save_m4n2(%%ymm12,%%ymm13) unit_save_m4n2(%%ymm14,%%ymm15)
#define SAVE_R_m4n4               unit_save_m4n2(%%ymm4,%%ymm5)   unit_save_m4n2(%%ymm6,%%ymm7)
#define SAVE_R_m4n8  SAVE_R_m4n4  unit_save_m4n2(%%ymm8,%%ymm9)   unit_save_m4n2(%%ymm10,%%ymm11)
#define SAVE_R_m4n12 SAVE_R_m4n8  unit_save_m4n2(%%ymm12,%%ymm13) unit_save_m4n2(%%ymm14,%%ymm15)
#define COMPUTE_L_m4(ndim,sim) \
    INIT_m4n##ndim\
    "movq %%r13,%4; movq %%r14,%1;"\
    #ndim""#sim"04042:\n\t"\
    "testq %4,%4; jz "#ndim""#sim"04043f;"\
    KERNEL_k1m4n##ndim(%1)\
    "decq %4; jmp "#ndim""#sim"04042b;"\
    #ndim""#sim"04043:\n\t"\
    SAVE_L_m4n##ndim "addq $64,%2;"
#define COMPUTE_R_m4(ndim,sim) \
    "subq %%r12,%0;"\
    INIT_m4n##ndim\
    "movq %%r13,%4; leaq (%%r14,%%r12,2),%%r15; addq %%r12,%%r15;"\
    #ndim""#sim"04042:\n\t"\
    "testq %4,%4; jz "#ndim""#sim"04043f;"\
    KERNEL_k1m4n##ndim(%%r15)\
    "decq %4; jmp "#ndim""#sim"04042b;"\
    #ndim""#sim"04043:\n\t"\
    SAVE_R_m4n##ndim
#define COMPUTE_m4_n1  COMPUTE_L_m4(1,33833)
#define COMPUTE_m4_n2  COMPUTE_L_m4(2,33833)
#define COMPUTE_m4_n4  COMPUTE_L_m4(4,33833)
#define COMPUTE_m4_n8  COMPUTE_L_m4(8,33833)
#define COMPUTE_m4_n12 COMPUTE_L_m4(12,33833)
#define COMPUTE_m4_n16 COMPUTE_L_m4(12,33733) COMPUTE_R_m4(4,33933)
#define COMPUTE_m4_n20 COMPUTE_L_m4(12,33633) COMPUTE_R_m4(8,33933)
#define COMPUTE_m4_n24 COMPUTE_L_m4(12,33533) COMPUTE_R_m4(12,33933)
#define COMPUTE_m4(ndim) COMPUTE_m4_n##ndim

/* m = 2 *//* vmm0 for alpha, vmm1-vmm3 for temporary use, vmm4-vmm15 for accumulators */
#define KERNEL_k1m2n1 \
    "vmovupd (%0),%%xmm1; addq $16,%0;"\
    "vmovddup (%1),%%xmm2; vfmadd231pd %%xmm1,%%xmm2,%%xmm4;"\
    "addq $8,%1;"
#define KERNEL_h_k1m2n2 \
    "vmovddup (%0),%%xmm1; vmovddup 8(%0),%%xmm2; addq $16,%0;"\
    "vmovupd (%1),%%xmm3; vfmadd231pd %%xmm1,%%xmm3,%%xmm4; vfmadd231pd %%xmm2,%%xmm3,%%xmm5;"
#define KERNEL_k1m2n2 KERNEL_h_k1m2n2 "addq $16,%1;"
#define unit_kernel_k1m2n4(c1,c2,...) \
    "vmovupd ("#__VA_ARGS__"),%%ymm3; vfmadd231pd %%ymm1,%%ymm3,"#c1"; vfmadd231pd %%ymm2,%%ymm3,"#c2";"
#define KERNEL_h_k1m2n4 \
    "vbroadcastsd (%0),%%ymm1; vbroadcastsd 8(%0),%%ymm2; addq $16,%0;"\
    unit_kernel_k1m2n4(%%ymm4,%%ymm5,%1)
#define KERNEL_k1m2n4 KERNEL_h_k1m2n4 "addq $32,%1;"
#define KERNEL_h_k1m2n8 KERNEL_h_k1m2n4 \
    unit_kernel_k1m2n4(%%ymm6,%%ymm7,%1,%%r12,1)
#define KERNEL_k1m2n8 KERNEL_h_k1m2n8 "addq $32,%1;"
#define KERNEL_h_k1m2n12 KERNEL_h_k1m2n8 \
    unit_kernel_k1m2n4(%%ymm8,%%ymm9,%1,%%r12,2)
#define KERNEL_k1m2n12 KERNEL_h_k1m2n12 "addq $32,%1;"
#define KERNEL_h_k1m2n16 KERNEL_k1m2n12 \
    unit_kernel_k1m2n4(%%ymm10,%%ymm11,%%r15)
#define KERNEL_k1m2n16 KERNEL_h_k1m2n16 "addq $32,%%r15;"
#define KERNEL_h_k1m2n20 KERNEL_h_k1m2n16 \
    unit_kernel_k1m2n4(%%ymm12,%%ymm13,%%r15,%%r12,1)
#define KERNEL_k1m2n20 KERNEL_h_k1m2n20 "addq $32,%%r15;"
#define KERNEL_h_k1m2n24 KERNEL_h_k1m2n20 \
    unit_kernel_k1m2n4(%%ymm14,%%ymm15,%%r15,%%r12,2)
#define KERNEL_k1m2n24 KERNEL_h_k1m2n24 "addq $32,%%r15;"
#define INIT_m2n1 "vpxor %%xmm4,%%xmm4,%%xmm4;"
#define INIT_m2n2 INIT_m2n1 "vpxor %%xmm5,%%xmm5,%%xmm5;"
#define unit_init_m2n4(c1,c2) "vpxor "#c1","#c1","#c1";vpxor "#c2","#c2","#c2";"
#define INIT_m2n4 unit_init_m2n4(%%ymm4,%%ymm5)
#define INIT_m2n8 INIT_m2n4 unit_init_m2n4(%%ymm6,%%ymm7)
#define INIT_m2n12 INIT_m2n8 unit_init_m2n4(%%ymm8,%%ymm9)
#define INIT_m2n16 INIT_m2n12 unit_init_m2n4(%%ymm10,%%ymm11)
#define INIT_m2n20 INIT_m2n16 unit_init_m2n4(%%ymm12,%%ymm13)
#define INIT_m2n24 INIT_m2n20 unit_init_m2n4(%%ymm14,%%ymm15)
#define SAVE_h_m2n1 \
    "vinsertf128 $1,%%xmm4,%%ymm4,%%ymm4; vpermilpd $12,%%ymm4,%%ymm4; vfmadd213pd (%2),%%ymm0,%%ymm4; vmovupd %%ymm4,(%2);"
#define SAVE_h_m2n2 \
    "vinsertf128 $1,%%xmm5,%%ymm4,%%ymm4; vunpcklpd %%ymm4,%%ymm4,%%ymm1; vunpckhpd %%ymm4,%%ymm4,%%ymm2;"\
    "vfmadd213pd (%2),%%ymm0,%%ymm1; vmovupd %%ymm1,(%2);"\
    "vfmadd213pd (%2,%3,1),%%ymm0,%%ymm2; vmovupd %%ymm2,(%2,%3,1);"
#define unit_save_m2n4(c1,c2) \
    "vperm2f128 $2,"#c1","#c2",%%ymm1; vunpcklpd %%ymm1,%%ymm1,%%ymm2; vunpckhpd %%ymm1,%%ymm1,%%ymm3;"\
    "vfmadd213pd (%5),%%ymm0,%%ymm2; vfmadd213pd (%5,%3,1),%%ymm0,%%ymm3; vmovupd %%ymm2,(%5); vmovupd %%ymm3,(%5,%3,1); leaq (%5,%3,2),%5;"\
    "vperm2f128 $19,"#c1","#c2",%%ymm1; vunpcklpd %%ymm1,%%ymm1,%%ymm2; vunpckhpd %%ymm1,%%ymm1,%%ymm3;"\
    "vfmadd213pd (%5),%%ymm0,%%ymm2; vfmadd213pd (%5,%3,1),%%ymm0,%%ymm3; vmovupd %%ymm2,(%5); vmovupd %%ymm3,(%5,%3,1); leaq (%5,%3,2),%5;"
#define SAVE_h_m2n4 "movq %2,%5;" unit_save_m2n4(%%ymm4,%%ymm5)
#define SAVE_h_m2n8 SAVE_h_m2n4 unit_save_m2n4(%%ymm6,%%ymm7)
#define SAVE_h_m2n12 SAVE_h_m2n8 unit_save_m2n4(%%ymm8,%%ymm9)
#define SAVE_h_m2n16 SAVE_h_m2n12 unit_save_m2n4(%%ymm10,%%ymm11)
#define SAVE_h_m2n20 SAVE_h_m2n16 unit_save_m2n4(%%ymm12,%%ymm13)
#define SAVE_h_m2n24 SAVE_h_m2n20 unit_save_m2n4(%%ymm14,%%ymm15)
#define SAVE_m2(ndim) SAVE_h_m2n##ndim "addq $32,%2;"
#define COMPUTE_m2(ndim) \
    INIT_m2n##ndim\
    "movq %%r13,%4; movq %%r14,%1; leaq (%1,%%r12,2),%%r15; addq %%r12,%%r15;"\
    #ndim"002022:\n\t"\
    "testq %4,%4; jz "#ndim"002023f;"\
    KERNEL_k1m2n##ndim\
    "decq %4; jmp "#ndim"002022b;"\
    #ndim"002023:\n\t"\
    SAVE_m2(ndim)

/* m = 1 *//* vmm0 for alpha, vmm1-vmm3 and vmm10-vmm15 for temporary use, vmm4-vmm9 for accumulators */
#define KERNEL_k1m1n1 \
    "vmovsd (%0),%%xmm1; addq $8,%0;"\
    "vfmadd231sd (%1),%%xmm1,%%xmm4; addq $8,%1;"
#define KERNEL_k1m1n2 \
    "vmovddup (%0),%%xmm1; addq $8,%0;"\
    "vfmadd231pd (%1),%%xmm1,%%xmm4; addq $16,%1;"
#define unit_kernel_k1m1n4(c1,...) \
    "vmovupd ("#__VA_ARGS__"),%%ymm2; vfmadd231pd %%ymm1,%%ymm2,"#c1";"
#define KERNEL_h_k1m1n4 \
    "vbroadcastsd (%0),%%ymm1; addq $8,%0;"\
    unit_kernel_k1m1n4(%%ymm4,%1)
#define KERNEL_k1m1n4 KERNEL_h_k1m1n4 "addq $32,%1;"
#define KERNEL_h_k1m1n8 KERNEL_h_k1m1n4 unit_kernel_k1m1n4(%%ymm5,%1,%%r12,1)
#define KERNEL_k1m1n8 KERNEL_h_k1m1n8 "addq $32,%1;"
#define KERNEL_h_k1m1n12 KERNEL_h_k1m1n8 unit_kernel_k1m1n4(%%ymm6,%1,%%r12,2)
#define KERNEL_k1m1n12 KERNEL_h_k1m1n12 "addq $32,%1;"
#define KERNEL_h_k1m1n16 KERNEL_k1m1n12 unit_kernel_k1m1n4(%%ymm7,%%r15)
#define KERNEL_k1m1n16 KERNEL_h_k1m1n16 "addq $32,%%r15;"
#define KERNEL_h_k1m1n20 KERNEL_h_k1m1n16 unit_kernel_k1m1n4(%%ymm8,%%r15,%%r12,1)
#define KERNEL_k1m1n20 KERNEL_h_k1m1n20 "addq $32,%%r15;"
#define KERNEL_h_k1m1n24 KERNEL_h_k1m1n20 unit_kernel_k1m1n4(%%ymm9,%%r15,%%r12,2)
#define KERNEL_k1m1n24 KERNEL_h_k1m1n24 "addq $32,%%r15;"
#define INIT_m1n1 INIT_m2n1
#define INIT_m1n2 INIT_m2n1
#define INIT_m1n4 "vpxor %%ymm4,%%ymm4,%%ymm4;"
#define INIT_m1n8 INIT_m1n4 "vpxor %%ymm5,%%ymm5,%%ymm5;"
#define INIT_m1n12 INIT_m1n8 "vpxor %%ymm6,%%ymm6,%%ymm6;"
#define INIT_m1n16 INIT_m1n12 "vpxor %%ymm7,%%ymm7,%%ymm7;"
#define INIT_m1n20 INIT_m1n16 "vpxor %%ymm8,%%ymm8,%%ymm8;"
#define INIT_m1n24 INIT_m1n20 "vpxor %%ymm9,%%ymm9,%%ymm9;"
#define SAVE_h_m1n1 \
    "vmovddup %%xmm4,%%xmm4; vfmadd213pd (%2),%%xmm0,%%xmm4; vmovupd %%xmm4,(%2);"
#define SAVE_h_m1n2 \
    "vunpcklpd %%xmm4,%%xmm4,%%xmm1; vunpckhpd %%xmm4,%%xmm4,%%xmm2;"\
    "vfmadd213pd (%2),%%xmm0,%%xmm1; vmovupd %%xmm1,(%2);"\
    "vfmadd213pd (%2,%3,1),%%xmm0,%%xmm2; vmovupd %%xmm2,(%2,%3,1);"
#define unit_save_m1n4(c1) \
    "vunpcklpd "#c1","#c1",%%ymm1; vunpckhpd "#c1","#c1",%%ymm2;"\
    "vmovupd (%5),%%xmm3; vinsertf128 $1,(%5,%3,2),%%ymm3,%%ymm3;"\
    "vfmadd213pd %%ymm3,%%ymm0,%%ymm1; vmovupd %%xmm1,(%5); vextractf128 $1,%%ymm1,(%5,%3,2); addq %3,%5;"\
    "vmovupd (%5),%%xmm3; vinsertf128 $1,(%5,%3,2),%%ymm3,%%ymm3;"\
    "vfmadd213pd %%ymm3,%%ymm0,%%ymm2; vmovupd %%xmm2,(%5); vextractf128 $1,%%ymm2,(%5,%3,2); addq %3,%5; leaq (%5,%3,2),%5;"
#define SAVE_h_m1n4 "movq %2,%5;" unit_save_m1n4(%%ymm4)
#define SAVE_h_m1n8 SAVE_h_m1n4 unit_save_m1n4(%%ymm5)
#define SAVE_h_m1n12 SAVE_h_m1n8 unit_save_m1n4(%%ymm6)
#define SAVE_h_m1n16 SAVE_h_m1n12 unit_save_m1n4(%%ymm7)
#define SAVE_h_m1n20 SAVE_h_m1n16 unit_save_m1n4(%%ymm8)
#define SAVE_h_m1n24 SAVE_h_m1n20 unit_save_m1n4(%%ymm9)
#define SAVE_m1(ndim) SAVE_h_m1n##ndim "addq $16,%2;"
#define COMPUTE_m1(ndim) \
    INIT_m1n##ndim\
    "movq %%r13,%4; movq %%r14,%1; leaq (%1,%%r12,2),%%r15; addq %%r12,%%r15;"\
    #ndim"001011:\n\t"\
    "testq %4,%4; jz "#ndim"001012f;"\
    KERNEL_k1m1n##ndim\
    "decq %4; jmp "#ndim"001011b;"\
    #ndim"001012:\n\t"\
    SAVE_m1(ndim)

#define COMPUTE(ndim) {\
    next_b = b_pointer + ndim * K;\
    __asm__ __volatile__(\
    "vbroadcastf32x4 (%6),%%zmm0; vmovups 16(%6),%%zmm1; vmovups 80(%6),%%zmm2;"\
    "movq %4,%%r13; movq %4,%%r12; salq $5,%%r12; movq %1,%%r14; movq %7,%%r11;"\
    "cmpq $8,%7;jb 33101"#ndim"f;"\
    "33109"#ndim":\n\t"\
    COMPUTE_m8(ndim)\
    "subq $8,%7;cmpq $8,%7;jnb 33109"#ndim"b;"\
    "33101"#ndim":\n\t"\
    "cmpq $4,%7;jb 33103"#ndim"f;"\
    COMPUTE_m4(ndim)\
    "subq $4,%7;"\
    "33103"#ndim":\n\t"\
    "cmpq $2,%7;jb 33104"#ndim"f;"\
    COMPUTE_m2(ndim)\
    "subq $2,%7;"\
    "33104"#ndim":\n\t"\
    "testq %7,%7;jz 33105"#ndim"f;"\
    COMPUTE_m1(ndim)\
    "33105"#ndim":\n\t"\
    "movq %%r13,%4; movq %%r14,%1; movq %%r11,%7;"\
    :"+r"(a_pointer),"+r"(b_pointer),"+r"(c_pointer),"+r"(ldc_in_bytes),"+r"(K),"+r"(ctemp),"+r"(const_val),"+r"(M),"+r"(next_b)\
    ::"r11","r12","r13","r14","r15","zmm0","zmm1","zmm2","zmm3","zmm4","zmm5","zmm6","zmm7","zmm8","zmm9","zmm10","zmm11","zmm12","zmm13","zmm14",\
    "zmm15","zmm16","zmm17","zmm18","zmm19","zmm20","zmm21","zmm22","zmm23","zmm24","zmm25","zmm26","zmm27","zmm28","zmm29","zmm30","zmm31",\
    "cc","memory");\
    a_pointer -= M * K; b_pointer += ndim * K; c_pointer += 2*(LDC * ndim - M);\
}
int __attribute__ ((noinline))
CNAME(BLASLONG m, BLASLONG n, BLASLONG k, double alphar, double alphai, double * __restrict__ A, double * __restrict__ B, double * __restrict__ C, BLASLONG LDC)
{
    if(m==0||n==0||k==0) return 0;
    int64_t ldc_in_bytes = (int64_t)LDC * sizeof(double) * 2;
    double constval[18]; constval[0] = alphar; constval[1] = alphai;
    double *const_val=constval;
    __asm__ __volatile__(
    "movq $0,16(%0);  movq $1,24(%0);  movq $8,32(%0); movq $9,40(%0);"
    "movq $2,48(%0);  movq $3,56(%0);  movq $10,64(%0);movq $11,72(%0);"
    "movq $12,80(%0); movq $13,88(%0); movq $4,96(%0); movq $5,104(%0);"
    "movq $14,112(%0);movq $15,120(%0);movq $6,128(%0);movq $7,136(%0);"
    ::"r"(const_val):"memory");//control words for permutations(vpermt2pd) before saving vectors of C elements
    int64_t M = (int64_t)m, K = (int64_t)k;
    BLASLONG n_count = n;
    double *a_pointer = A,*b_pointer = B,*c_pointer = C,*ctemp = C,*next_b = B;
    for(;n_count>23;n_count-=24) COMPUTE(24)
    for(;n_count>19;n_count-=20) COMPUTE(20)
    for(;n_count>15;n_count-=16) COMPUTE(16)
    for(;n_count>11;n_count-=12) COMPUTE(12)
    for(;n_count>7;n_count-=8) COMPUTE(8)
    for(;n_count>3;n_count-=4) COMPUTE(4)
    for(;n_count>1;n_count-=2) COMPUTE(2)
    if(n_count>0) COMPUTE(1)
    return 0;
}
