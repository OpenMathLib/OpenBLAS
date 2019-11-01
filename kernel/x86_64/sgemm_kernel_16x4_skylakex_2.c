/* %0 = "+r"(a_pointer), %1 = "+r"(b_pointer), %2 = "+r"(c_pointer), %3 = "+r"(ldc_in_bytes), %4 for k_count, %5 for c_store */
/* r12 = k << 4(const), r13 = k(const), r14 = b_head_pos(const), r15 = %1 + 3r12 */

#include "common.h"
#include <stdint.h>

/* m = 16 */ /* zmm8-zmm31 for accumulators, zmm1-zmm7 for temporary use, zmm0 for alpha */
#define KERNEL_k1m16n1 \
    "vmovups (%0),%%zmm4; addq $64,%0;"\
    "vbroadcastss (%1),%%zmm6; vfmadd231ps %%zmm4,%%zmm6,%%zmm8;"\
    "addq $4,%1;"
#define KERNEL_h_k1m16n2 \
    "vmovsldup (%0),%%zmm4; vmovshdup (%0),%%zmm5; prefetcht0 512(%0); addq $64,%0;"\
    "vbroadcastsd (%1),%%zmm6; vfmadd231ps %%zmm4,%%zmm6,%%zmm8; vfmadd231ps %%zmm5,%%zmm6,%%zmm9;"
#define KERNEL_k1m16n2 KERNEL_h_k1m16n2 "addq $8,%1;"
#define KERNEL_h_k1m16n4 KERNEL_h_k1m16n2 "vbroadcastsd 8(%1),%%zmm7; vfmadd231ps %%zmm4,%%zmm7,%%zmm10; vfmadd231ps %%zmm5,%%zmm7,%%zmm11;"
#define KERNEL_k1m16n4 KERNEL_h_k1m16n4 "addq $16,%1;"
#define unit_kernel_k1m16n4(c1,c2,c3,c4, ...) \
    "vbroadcastsd  ("#__VA_ARGS__"),%%zmm6; vfmadd231ps %%zmm4,%%zmm6,"#c1"; vfmadd231ps %%zmm5,%%zmm6,"#c2";"\
    "vbroadcastsd 8("#__VA_ARGS__"),%%zmm7; vfmadd231ps %%zmm4,%%zmm7,"#c3"; vfmadd231ps %%zmm5,%%zmm7,"#c4";"
#define KERNEL_h_k1m16n8 KERNEL_h_k1m16n4 unit_kernel_k1m16n4(%%zmm12,%%zmm13,%%zmm14,%%zmm15,%1,%%r12,1)
#define KERNEL_k1m16n8 KERNEL_h_k1m16n8 "addq $16,%1;"
#define KERNEL_h_k1m16n12 KERNEL_h_k1m16n8 unit_kernel_k1m16n4(%%zmm16,%%zmm17,%%zmm18,%%zmm19,%1,%%r12,2)
#define KERNEL_k1m16n12 KERNEL_h_k1m16n12 "addq $16,%1;"
#define KERNEL_h_k1m16n16 KERNEL_k1m16n12 unit_kernel_k1m16n4(%%zmm20,%%zmm21,%%zmm22,%%zmm23,%%r15)
#define KERNEL_k1m16n16 KERNEL_h_k1m16n16 "addq $16,%%r15;"
#define KERNEL_h_k1m16n20 KERNEL_h_k1m16n16 unit_kernel_k1m16n4(%%zmm24,%%zmm25,%%zmm26,%%zmm27,%%r15,%%r12,1)
#define KERNEL_k1m16n20 KERNEL_h_k1m16n20 "addq $16,%%r15;"
#define KERNEL_h_k1m16n24 KERNEL_h_k1m16n20 unit_kernel_k1m16n4(%%zmm28,%%zmm29,%%zmm30,%%zmm31,%%r15,%%r12,2)
#define KERNEL_k1m16n24 KERNEL_h_k1m16n24 "addq $16,%%r15;"
#define INIT_m16n1 "vpxorq %%zmm8,%%zmm8,%%zmm8;"
#define INIT_m16n2 INIT_m16n1 "vpxorq %%zmm9,%%zmm9,%%zmm9;"
#define INIT_m16n4 INIT_m16n2 "vpxorq %%zmm10,%%zmm10,%%zmm10;vpxorq %%zmm11,%%zmm11,%%zmm11;"
#define unit_init_m16n4(c1,c2,c3,c4) \
    "vpxorq "#c1","#c1","#c1";vpxorq "#c2","#c2","#c2";vpxorq "#c3","#c3","#c3";vpxorq "#c4","#c4","#c4";"
#define INIT_m16n8 INIT_m16n4 unit_init_m16n4(%%zmm12,%%zmm13,%%zmm14,%%zmm15)
#define INIT_m16n12 INIT_m16n8 unit_init_m16n4(%%zmm16,%%zmm17,%%zmm18,%%zmm19)
#define INIT_m16n16 INIT_m16n12 unit_init_m16n4(%%zmm20,%%zmm21,%%zmm22,%%zmm23)
#define INIT_m16n20 INIT_m16n16 unit_init_m16n4(%%zmm24,%%zmm25,%%zmm26,%%zmm27)
#define INIT_m16n24 INIT_m16n20 unit_init_m16n4(%%zmm28,%%zmm29,%%zmm30,%%zmm31)
#define SAVE_h_m16n1 "vfmadd213ps (%2),%%zmm0,%%zmm8; vmovups %%zmm8,(%2);"
#define unit_save_m16n2(c1,c2) \
    "vunpcklps "#c2","#c1",%%zmm6; vunpckhps "#c2","#c1",%%zmm7; vunpcklpd %%zmm7,%%zmm6,%%zmm4; vunpckhpd %%zmm7,%%zmm6,%%zmm5;"\
    "vfmadd213ps (%5),%%zmm0,%%zmm4; vfmadd213ps (%5,%3,1),%%zmm0,%%zmm5;"\
    "prefetcht1 127(%5); prefetcht1 127(%5,%3,1);"\
    "vmovups %%zmm4,(%5); vmovups %%zmm5,(%5,%3,1); leaq (%5,%3,2),%5;"
#define SAVE_h_m16n2 "movq %2,%5;" unit_save_m16n2(%%zmm8,%%zmm9)
#define SAVE_h_m16n4  SAVE_h_m16n2  unit_save_m16n2(%%zmm10,%%zmm11)
#define SAVE_h_m16n8  SAVE_h_m16n4  unit_save_m16n2(%%zmm12,%%zmm13) unit_save_m16n2(%%zmm14,%%zmm15)
#define SAVE_h_m16n12 SAVE_h_m16n8  unit_save_m16n2(%%zmm16,%%zmm17) unit_save_m16n2(%%zmm18,%%zmm19)
#define SAVE_h_m16n16 SAVE_h_m16n12 unit_save_m16n2(%%zmm20,%%zmm21) unit_save_m16n2(%%zmm22,%%zmm23)
#define SAVE_h_m16n20 SAVE_h_m16n16 unit_save_m16n2(%%zmm24,%%zmm25) unit_save_m16n2(%%zmm26,%%zmm27)
#define SAVE_h_m16n24 SAVE_h_m16n20 unit_save_m16n2(%%zmm28,%%zmm29) unit_save_m16n2(%%zmm30,%%zmm31)
#define SAVE_m16(ndim) SAVE_h_m16n##ndim "addq $64,%2;"
#define COMPUTE_m16(ndim) \
    INIT_m16n##ndim\
    "movq %%r13,%4; movq %%r14,%1; leaq (%1,%%r12,2),%%r15; addq %%r12,%%r15;"\
    "cmpq $4,%4; jb "#ndim"016162f;"\
    #ndim"016161:\n\t"\
    KERNEL_k1m16n##ndim\
    KERNEL_k1m16n##ndim\
    KERNEL_k1m16n##ndim\
    KERNEL_k1m16n##ndim\
    "subq $4,%4; cmpq $4,%4; jnb "#ndim"016161b;"\
    #ndim"016162:\n\t"\
    "testq %4,%4; jz "#ndim"016163f;"\
    KERNEL_k1m16n##ndim\
    "decq %4; jmp "#ndim"016162b;"\
    #ndim"016163:\n\t"\
    SAVE_m16(ndim)

/* m = 8 *//* ymm0 for alpha, ymm1-ymm3 for temporary use, ymm4-ymm15 for accumulators */
#define KERNEL_k1m8n1(b_addr) \
    "vmovups (%0),%%ymm1; addq $32,%0;"\
    "vbroadcastss ("#b_addr"),%%ymm2; vfmadd231ps %%ymm1,%%ymm2,%%ymm4;"\
    "addq $4,"#b_addr";"
#define KERNEL_h_k1m8n2(b_addr) \
    "vmovsldup (%0),%%ymm1; vmovshdup (%0),%%ymm2; addq $32,%0;"\
    "vbroadcastsd ("#b_addr"),%%ymm3; vfmadd231ps %%ymm1,%%ymm3,%%ymm4; vfmadd231ps %%ymm2,%%ymm3,%%ymm5;"
#define KERNEL_k1m8n2(b_addr) KERNEL_h_k1m8n2(b_addr) "addq $8,"#b_addr";"
#define KERNEL_h_k1m8n4(b_addr) \
    KERNEL_h_k1m8n2(b_addr) "vbroadcastsd 8("#b_addr"),%%ymm3; vfmadd231ps %%ymm1,%%ymm3,%%ymm6; vfmadd231ps %%ymm2,%%ymm3,%%ymm7;"
#define KERNEL_k1m8n4(b_addr) KERNEL_h_k1m8n4(b_addr) "addq $16,"#b_addr";"
#define unit_kernel_k1m8n4(c1,c2,c3,c4,...) \
    "vbroadcastsd  ("#__VA_ARGS__"),%%ymm3; vfmadd231ps %%ymm1,%%ymm3,"#c1"; vfmadd231ps %%ymm2,%%ymm3,"#c2";"\
    "vbroadcastsd 8("#__VA_ARGS__"),%%ymm3; vfmadd231ps %%ymm1,%%ymm3,"#c3"; vfmadd231ps %%ymm2,%%ymm3,"#c4";"
#define KERNEL_h_k1m8n8(b_addr) KERNEL_h_k1m8n4(b_addr) unit_kernel_k1m8n4(%%ymm8,%%ymm9,%%ymm10,%%ymm11,b_addr,%%r12,1)
#define KERNEL_k1m8n8(b_addr) KERNEL_h_k1m8n8(b_addr) "addq $16,"#b_addr";"
#define KERNEL_h_k1m8n12(b_addr) KERNEL_h_k1m8n8(b_addr) unit_kernel_k1m8n4(%%ymm12,%%ymm13,%%ymm14,%%ymm15,b_addr,%%r12,2)
#define KERNEL_k1m8n12(b_addr) KERNEL_h_k1m8n12(b_addr) "addq $16,"#b_addr";"
#define INIT_m8n1 "vpxor %%ymm4,%%ymm4,%%ymm4;"
#define INIT_m8n2 INIT_m8n1 "vpxor %%ymm5,%%ymm5,%%ymm5;"
#define INIT_m8n4 INIT_m8n2 "vpxor %%ymm6,%%ymm6,%%ymm6;vpxor %%ymm7,%%ymm7,%%ymm7;"
#define unit_init_m8n4(c1,c2,c3,c4) \
    "vpxor "#c1","#c1","#c1";vpxor "#c2","#c2","#c2";vpxor "#c3","#c3","#c3";vpxor "#c4","#c4","#c4";"
#define INIT_m8n8  INIT_m8n4 unit_init_m8n4(%%ymm8,%%ymm9,%%ymm10,%%ymm11)
#define INIT_m8n12 INIT_m8n8 unit_init_m8n4(%%ymm12,%%ymm13,%%ymm14,%%ymm15)
#define SAVE_L_m8n1 "vfmadd213ps (%2),%%ymm0,%%ymm4; vmovups %%ymm4,(%2);"
#define unit_save_m8n2(c1,c2) \
    "vunpcklps "#c2","#c1",%%ymm2; vunpckhps "#c2","#c1",%%ymm3;"\
    "vunpcklpd %%ymm3,%%ymm2,%%ymm1;vfmadd213ps (%5),     %%ymm0,%%ymm1;vmovups %%ymm1,(%5);"\
    "vunpckhpd %%ymm3,%%ymm2,%%ymm1;vfmadd213ps (%5,%3,1),%%ymm0,%%ymm1;vmovups %%ymm1,(%5,%3,1);"\
    "leaq (%5,%3,2),%5;"
#define SAVE_L_m8n2 "movq %2,%5;" unit_save_m8n2(%%ymm4,%%ymm5)
#define SAVE_L_m8n4  SAVE_L_m8n2  unit_save_m8n2(%%ymm6,%%ymm7)
#define SAVE_L_m8n8  SAVE_L_m8n4  unit_save_m8n2(%%ymm8,%%ymm9)   unit_save_m8n2(%%ymm10,%%ymm11)
#define SAVE_L_m8n12 SAVE_L_m8n8  unit_save_m8n2(%%ymm12,%%ymm13) unit_save_m8n2(%%ymm14,%%ymm15)
#define SAVE_R_m8n4               unit_save_m8n2(%%ymm4,%%ymm5)   unit_save_m8n2(%%ymm6,%%ymm7)
#define SAVE_R_m8n8  SAVE_R_m8n4  unit_save_m8n2(%%ymm8,%%ymm9)   unit_save_m8n2(%%ymm10,%%ymm11)
#define SAVE_R_m8n12 SAVE_R_m8n8  unit_save_m8n2(%%ymm12,%%ymm13) unit_save_m8n2(%%ymm14,%%ymm15)
#define COMPUTE_L_m8(ndim,sim) \
    INIT_m8n##ndim\
    "movq %%r13,%4; movq %%r14,%1;"\
    #ndim""#sim"882:\n\t"\
    "testq %4,%4; jz "#ndim""#sim"883f;"\
    KERNEL_k1m8n##ndim(%1)\
    "decq %4; jmp "#ndim""#sim"882b;"\
    #ndim""#sim"883:\n\t"\
    SAVE_L_m8n##ndim "addq $32,%2;"
#define COMPUTE_R_m8(ndim,sim) \
    "subq %%r12,%0; subq %%r12,%0;"\
    INIT_m8n##ndim\
    "movq %%r13,%4; leaq (%%r14,%%r12,2),%%r15; addq %%r12,%%r15;"\
    #ndim""#sim"882:\n\t"\
    "testq %4,%4; jz "#ndim""#sim"883f;"\
    KERNEL_k1m8n##ndim(%%r15)\
    "decq %4; jmp "#ndim""#sim"882b;"\
    #ndim""#sim"883:\n\t"\
    SAVE_R_m8n##ndim
#define COMPUTE_m8_n1  COMPUTE_L_m8(1,33833)
#define COMPUTE_m8_n2  COMPUTE_L_m8(2,33833)
#define COMPUTE_m8_n4  COMPUTE_L_m8(4,33833)
#define COMPUTE_m8_n8  COMPUTE_L_m8(8,33833)
#define COMPUTE_m8_n12 COMPUTE_L_m8(12,33833)
#define COMPUTE_m8_n16 COMPUTE_L_m8(12,33733) COMPUTE_R_m8(4,33933)
#define COMPUTE_m8_n20 COMPUTE_L_m8(12,33633) COMPUTE_R_m8(8,33933)
#define COMPUTE_m8_n24 COMPUTE_L_m8(12,33533) COMPUTE_R_m8(12,33933)
#define COMPUTE_m8(ndim) COMPUTE_m8_n##ndim

/* m = 4 *//* xmm0 for alpha, xmm1-xmm3 for temporary use, xmm4-xmm15 for accumulators */
#define KERNEL_k1m4n1(b_addr) \
    "vmovups (%0),%%xmm1; addq $16,%0;"\
    "vbroadcastss ("#b_addr"),%%xmm2; vfmadd231ps %%xmm1,%%xmm2,%%xmm4;"\
    "addq $4,"#b_addr";"
#define KERNEL_h_k1m4n2(b_addr) \
    "vmovsldup (%0),%%xmm1; vmovshdup (%0),%%xmm2; addq $16,%0;"\
    "vmovddup ("#b_addr"),%%xmm3; vfmadd231ps %%xmm1,%%xmm3,%%xmm4; vfmadd231ps %%xmm2,%%xmm3,%%xmm5;"
#define KERNEL_k1m4n2(b_addr) KERNEL_h_k1m4n2(b_addr) "addq $8,"#b_addr";"
#define KERNEL_h_k1m4n4(b_addr) \
    KERNEL_h_k1m4n2(b_addr) "vmovddup 8("#b_addr"),%%xmm3; vfmadd231ps %%xmm1,%%xmm3,%%xmm6; vfmadd231ps %%xmm2,%%xmm3,%%xmm7;"
#define KERNEL_k1m4n4(b_addr) KERNEL_h_k1m4n4(b_addr) "addq $16,"#b_addr";"
#define unit_kernel_k1m4n4(c1,c2,c3,c4,...) \
    "vmovddup  ("#__VA_ARGS__"),%%xmm3; vfmadd231ps %%xmm1,%%xmm3,"#c1"; vfmadd231ps %%xmm2,%%xmm3,"#c2";"\
    "vmovddup 8("#__VA_ARGS__"),%%xmm3; vfmadd231ps %%xmm1,%%xmm3,"#c3"; vfmadd231ps %%xmm2,%%xmm3,"#c4";"
#define KERNEL_h_k1m4n8(b_addr) KERNEL_h_k1m4n4(b_addr) unit_kernel_k1m4n4(%%xmm8,%%xmm9,%%xmm10,%%xmm11,b_addr,%%r12,1)
#define KERNEL_k1m4n8(b_addr) KERNEL_h_k1m4n8(b_addr) "addq $16,"#b_addr";"
#define KERNEL_h_k1m4n12(b_addr) KERNEL_h_k1m4n8(b_addr) unit_kernel_k1m4n4(%%xmm12,%%xmm13,%%xmm14,%%xmm15,b_addr,%%r12,2)
#define KERNEL_k1m4n12(b_addr) KERNEL_h_k1m4n12(b_addr) "addq $16,"#b_addr";"
#define INIT_m4n1 "vpxor %%xmm4,%%xmm4,%%xmm4;"
#define INIT_m4n2 INIT_m4n1 "vpxor %%xmm5,%%xmm5,%%xmm5;"
#define INIT_m4n4 INIT_m4n2 "vpxor %%xmm6,%%xmm6,%%xmm6;vpxor %%xmm7,%%xmm7,%%xmm7;"
#define unit_init_m4n4(c1,c2,c3,c4) \
    "vpxor "#c1","#c1","#c1";vpxor "#c2","#c2","#c2";vpxor "#c3","#c3","#c3";vpxor "#c4","#c4","#c4";"
#define INIT_m4n8  INIT_m4n4 unit_init_m4n4(%%xmm8,%%xmm9,%%xmm10,%%xmm11)
#define INIT_m4n12 INIT_m4n8 unit_init_m4n4(%%xmm12,%%xmm13,%%xmm14,%%xmm15)
#define SAVE_L_m4n1 "vfmadd213ps (%2),%%xmm0,%%xmm4; vmovups %%xmm4,(%2);"
#define unit_save_m4n2(c1,c2) \
    "vunpcklps "#c2","#c1",%%xmm2; vunpckhps "#c2","#c1",%%xmm3;"\
    "vunpcklpd %%xmm3,%%xmm2,%%xmm1;vfmadd213ps (%5),     %%xmm0,%%xmm1;vmovups %%xmm1,(%5);"\
    "vunpckhpd %%xmm3,%%xmm2,%%xmm1;vfmadd213ps (%5,%3,1),%%xmm0,%%xmm1;vmovups %%xmm1,(%5,%3,1);"\
    "leaq (%5,%3,2),%5;"
#define SAVE_L_m4n2 "movq %2,%5;" unit_save_m4n2(%%xmm4,%%xmm5)
#define SAVE_L_m4n4  SAVE_L_m4n2  unit_save_m4n2(%%xmm6,%%xmm7)
#define SAVE_L_m4n8  SAVE_L_m4n4  unit_save_m4n2(%%xmm8,%%xmm9)   unit_save_m4n2(%%xmm10,%%xmm11)
#define SAVE_L_m4n12 SAVE_L_m4n8  unit_save_m4n2(%%xmm12,%%xmm13) unit_save_m4n2(%%xmm14,%%xmm15)
#define SAVE_R_m4n4               unit_save_m4n2(%%xmm4,%%xmm5)   unit_save_m4n2(%%xmm6,%%xmm7)
#define SAVE_R_m4n8  SAVE_R_m4n4  unit_save_m4n2(%%xmm8,%%xmm9)   unit_save_m4n2(%%xmm10,%%xmm11)
#define SAVE_R_m4n12 SAVE_R_m4n8  unit_save_m4n2(%%xmm12,%%xmm13) unit_save_m4n2(%%xmm14,%%xmm15)
#define COMPUTE_L_m4(ndim,sim) \
    INIT_m4n##ndim\
    "movq %%r13,%4; movq %%r14,%1;"\
    #ndim""#sim"442:\n\t"\
    "testq %4,%4; jz "#ndim""#sim"443f;"\
    KERNEL_k1m4n##ndim(%1)\
    "decq %4; jmp "#ndim""#sim"442b;"\
    #ndim""#sim"443:\n\t"\
    SAVE_L_m4n##ndim "addq $16,%2;"
#define COMPUTE_R_m4(ndim,sim) \
    "subq %%r12,%0;"\
    INIT_m4n##ndim\
    "movq %%r13,%4; leaq (%%r14,%%r12,2),%%r15; addq %%r12,%%r15;"\
    #ndim""#sim"442:\n\t"\
    "testq %4,%4; jz "#ndim""#sim"443f;"\
    KERNEL_k1m4n##ndim(%%r15)\
    "decq %4; jmp "#ndim""#sim"442b;"\
    #ndim""#sim"443:\n\t"\
    SAVE_R_m4n##ndim
#define COMPUTE_m4_n1  COMPUTE_L_m4(1,55855)
#define COMPUTE_m4_n2  COMPUTE_L_m4(2,55855)
#define COMPUTE_m4_n4  COMPUTE_L_m4(4,55855)
#define COMPUTE_m4_n8  COMPUTE_L_m4(8,55855)
#define COMPUTE_m4_n12 COMPUTE_L_m4(12,55855)
#define COMPUTE_m4_n16 COMPUTE_L_m4(12,55755) COMPUTE_R_m4(4,55955)
#define COMPUTE_m4_n20 COMPUTE_L_m4(12,55655) COMPUTE_R_m4(8,55955)
#define COMPUTE_m4_n24 COMPUTE_L_m4(12,55555) COMPUTE_R_m4(12,55955)
#define COMPUTE_m4(ndim) COMPUTE_m4_n##ndim

/* m = 2 *//* xmm0 for alpha, xmm1-xmm3 and xmm10 for temporary use, xmm4-xmm9 for accumulators */
#define INIT_m2n1 "vpxor %%xmm4,%%xmm4,%%xmm4;"
#define KERNEL_k1m2n1(b_addr) \
    "vmovsd (%0),%%xmm1; addq $8,%0;"\
    "vbroadcastss ("#b_addr"),%%xmm2; vfmadd231ps %%xmm1,%%xmm2,%%xmm4;"\
    "addq $4,"#b_addr";"
#define SAVE_L_m2n1 "vmovsd (%2),%%xmm1; vfmadd213ps %%xmm1,%%xmm0,%%xmm4; vmovsd %%xmm4,(%2);"
#define INIT_m2n2 INIT_m2n1 "vpxor %%xmm5,%%xmm5,%%xmm5;"
#define KERNEL_k1m2n2(b_addr) \
    "vmovsd (%0),%%xmm1; addq $8,%0;"\
    "vbroadcastss  ("#b_addr"),%%xmm2; vfmadd231ps %%xmm1,%%xmm2,%%xmm4;"\
    "vbroadcastss 4("#b_addr"),%%xmm3; vfmadd231ps %%xmm1,%%xmm3,%%xmm5;"\
    "addq $8,"#b_addr";"
#define SAVE_L_m2n2 SAVE_L_m2n1 "vmovsd (%2,%3,1),%%xmm1; vfmadd213ps %%xmm1,%%xmm0,%%xmm5; vmovsd %%xmm5,(%2,%3,1);"
#define INIT_m2n4  INIT_m2n2
#define INIT_m2n8  INIT_m2n4 "vpxor %%xmm6,%%xmm6,%%xmm6; vpxor %%xmm7,%%xmm7,%%xmm7;"
#define INIT_m2n12 INIT_m2n8 "vpxor %%xmm8,%%xmm8,%%xmm8; vpxor %%xmm9,%%xmm9,%%xmm9;"
#define KERNEL_k1m2n4(b_addr) \
    "vmovups ("#b_addr"),%%xmm3; addq $16,"#b_addr";"\
    "vbroadcastss  (%0),%%xmm1; vfmadd231ps %%xmm3,%%xmm1,%%xmm4;"\
    "vbroadcastss 4(%0),%%xmm2; vfmadd231ps %%xmm3,%%xmm2,%%xmm5;"\
    "addq $8,%0;"
#define KERNEL_k1m2n8(b_addr) \
    "vmovups ("#b_addr"),%%xmm3; vmovups ("#b_addr",%%r12,1),%%xmm2; addq $16,"#b_addr";"\
    "vbroadcastss  (%0),%%xmm1; vfmadd231ps %%xmm3,%%xmm1,%%xmm4; vfmadd231ps %%xmm2,%%xmm1,%%xmm6;"\
    "vbroadcastss 4(%0),%%xmm1; vfmadd231ps %%xmm3,%%xmm1,%%xmm5; vfmadd231ps %%xmm2,%%xmm1,%%xmm7;"\
    "addq $8,%0;"
#define KERNEL_k1m2n12(b_addr) \
    "vmovups ("#b_addr"),%%xmm3; vmovups ("#b_addr",%%r12,1),%%xmm2; vmovups ("#b_addr",%%r12,2),%%xmm1; addq $16,"#b_addr";"\
    "vbroadcastss  (%0),%%xmm10; vfmadd231ps %%xmm3,%%xmm10,%%xmm4; vfmadd231ps %%xmm2,%%xmm10,%%xmm6; vfmadd231ps %%xmm1,%%xmm10,%%xmm8;"\
    "vbroadcastss 4(%0),%%xmm10; vfmadd231ps %%xmm3,%%xmm10,%%xmm5; vfmadd231ps %%xmm2,%%xmm10,%%xmm7; vfmadd231ps %%xmm1,%%xmm10,%%xmm9;"\
    "addq $8,%0;"
#define unit_save_m2n4(c1,c2) \
    "vunpcklps "#c2","#c1",%%xmm1; vunpckhps "#c2","#c1",%%xmm2;"\
    "vmovsd (%5),%%xmm3; vmovhpd (%5,%3,1),%%xmm3,%%xmm3; vfmadd213ps %%xmm3,%%xmm0,%%xmm1; vmovsd %%xmm1,(%5); vmovhpd %%xmm1,(%5,%3,1);"\
    "leaq (%5,%3,2),%5;"\
    "vmovsd (%5),%%xmm3; vmovhpd (%5,%3,1),%%xmm3,%%xmm3; vfmadd213ps %%xmm3,%%xmm0,%%xmm2; vmovsd %%xmm2,(%5); vmovhpd %%xmm2,(%5,%3,1);"\
    "leaq (%5,%3,2),%5;"
#define SAVE_L_m2n4  "movq %2,%5;" unit_save_m2n4(%%xmm4,%%xmm5)
#define SAVE_L_m2n8  SAVE_L_m2n4   unit_save_m2n4(%%xmm6,%%xmm7)
#define SAVE_L_m2n12 SAVE_L_m2n8   unit_save_m2n4(%%xmm8,%%xmm9)
#define SAVE_R_m2n4                unit_save_m2n4(%%xmm4,%%xmm5)
#define SAVE_R_m2n8  SAVE_R_m2n4   unit_save_m2n4(%%xmm6,%%xmm7)
#define SAVE_R_m2n12 SAVE_R_m2n8   unit_save_m2n4(%%xmm8,%%xmm9)
#define COMPUTE_L_m2(ndim,sim) \
    INIT_m2n##ndim\
    "movq %%r13,%4; movq %%r14,%1;"\
    #ndim""#sim"222:\n\t"\
    "testq %4,%4; jz "#ndim""#sim"223f;"\
    KERNEL_k1m2n##ndim(%1)\
    "decq %4; jmp "#ndim""#sim"222b;"\
    #ndim""#sim"223:\n\t"\
    SAVE_L_m2n##ndim "addq $8,%2;"
#define COMPUTE_R_m2(ndim,sim) \
    "salq $3,%%r13;subq %%r13,%0;sarq $3,%%r13;"\
    INIT_m2n##ndim\
    "movq %%r13,%4; leaq (%%r14,%%r12,2),%%r15; addq %%r12,%%r15;"\
    #ndim""#sim"222:\n\t"\
    "testq %4,%4; jz "#ndim""#sim"223f;"\
    KERNEL_k1m2n##ndim(%%r15)\
    "decq %4; jmp "#ndim""#sim"222b;"\
    #ndim""#sim"223:\n\t"\
    SAVE_R_m2n##ndim
#define COMPUTE_m2_n1  COMPUTE_L_m2(1,77877)
#define COMPUTE_m2_n2  COMPUTE_L_m2(2,77877)
#define COMPUTE_m2_n4  COMPUTE_L_m2(4,77877)
#define COMPUTE_m2_n8  COMPUTE_L_m2(8,77877)
#define COMPUTE_m2_n12 COMPUTE_L_m2(12,77877)
#define COMPUTE_m2_n16 COMPUTE_L_m2(12,77777) COMPUTE_R_m2(4,77977)
#define COMPUTE_m2_n20 COMPUTE_L_m2(12,77677) COMPUTE_R_m2(8,77977)
#define COMPUTE_m2_n24 COMPUTE_L_m2(12,77577) COMPUTE_R_m2(12,77977)
#define COMPUTE_m2(ndim) COMPUTE_m2_n##ndim

/* m = 1 *//* xmm0 for alpha, xmm1-xmm3 and xmm10 for temporary use, xmm4-xmm6 for accumulators */
#define INIT_m1n1 "vpxor %%xmm4,%%xmm4,%%xmm4;"
#define KERNEL_k1m1n1(b_addr) \
    "vmovss ("#b_addr"),%%xmm3; addq $4,"#b_addr";"\
    "vmovss (%0),%%xmm1; vfmadd231ss %%xmm3,%%xmm1,%%xmm4;"\
    "addq $4,%0;"
#define SAVE_L_m1n1 "vfmadd213ss (%2),%%xmm0,%%xmm4; vmovss %%xmm4,(%2);"
#define INIT_m1n2 INIT_m1n1
#define KERNEL_k1m1n2(b_addr) \
    "vmovsd ("#b_addr"),%%xmm3; addq $8,"#b_addr";"\
    "vbroadcastss  (%0),%%xmm1; vfmadd231ps %%xmm3,%%xmm1,%%xmm4;"\
    "addq $4,%0;"
#define SAVE_L_m1n2 \
    "vmovss (%2),%%xmm3; vinsertps $16,(%2,%3,1),%%xmm3,%%xmm3; vfmadd213ps %%xmm3,%%xmm0,%%xmm4;"\
    "vmovss %%xmm4,(%2); vextractps $1,%%xmm4,(%2,%3,1);"
#define INIT_m1n4  INIT_m1n2
#define INIT_m1n8  INIT_m1n4 "vpxor %%xmm5,%%xmm5,%%xmm5;"
#define INIT_m1n12 INIT_m1n8 "vpxor %%xmm6,%%xmm6,%%xmm6;"
#define KERNEL_k1m1n4(b_addr) \
    "vmovups ("#b_addr"),%%xmm3; addq $16,"#b_addr";"\
    "vbroadcastss  (%0),%%xmm1; vfmadd231ps %%xmm3,%%xmm1,%%xmm4;"\
    "addq $4,%0;"
#define KERNEL_k1m1n8(b_addr) \
    "vmovups ("#b_addr"),%%xmm3; vmovups ("#b_addr",%%r12,1),%%xmm2; addq $16,"#b_addr";"\
    "vbroadcastss  (%0),%%xmm1; vfmadd231ps %%xmm3,%%xmm1,%%xmm4; vfmadd231ps %%xmm2,%%xmm1,%%xmm5;"\
    "addq $4,%0;"
#define KERNEL_k1m1n12(b_addr) \
    "vmovups ("#b_addr"),%%xmm3; vmovups ("#b_addr",%%r12,1),%%xmm2; vmovups ("#b_addr",%%r12,2),%%xmm1; addq $16,"#b_addr";"\
    "vbroadcastss  (%0),%%xmm10; vfmadd231ps %%xmm3,%%xmm10,%%xmm4; vfmadd231ps %%xmm2,%%xmm10,%%xmm5; vfmadd231ps %%xmm1,%%xmm10,%%xmm6;"\
    "addq $4,%0;"
#define unit_save_m1n4(c1) \
    "vpxor %%xmm10,%%xmm10,%%xmm10; vmovsd "#c1",%%xmm10,%%xmm2; vmovhlps "#c1",%%xmm10,%%xmm1;"\
    "vmovss (%5),%%xmm3; vinsertps $16,(%5,%3,1),%%xmm3,%%xmm3; vfmadd213ps %%xmm3,%%xmm0,%%xmm2;"\
    "vmovss %%xmm2,(%5); vextractps $1,%%xmm2,(%5,%3,1); leaq (%5,%3,2),%5;"\
    "vmovss (%5),%%xmm3; vinsertps $16,(%5,%3,1),%%xmm3,%%xmm3; vfmadd213ps %%xmm3,%%xmm0,%%xmm1;"\
    "vmovss %%xmm1,(%5); vextractps $1,%%xmm1,(%5,%3,1); leaq (%5,%3,2),%5;"
#define SAVE_L_m1n4 "movq %2,%5;" unit_save_m1n4(%%xmm4)
#define SAVE_L_m1n8  SAVE_L_m1n4  unit_save_m1n4(%%xmm5)
#define SAVE_L_m1n12 SAVE_L_m1n8  unit_save_m1n4(%%xmm6)
#define SAVE_R_m1n4               unit_save_m1n4(%%xmm4)
#define SAVE_R_m1n8  SAVE_R_m1n4  unit_save_m1n4(%%xmm5)
#define SAVE_R_m1n12 SAVE_R_m1n8  unit_save_m1n4(%%xmm6)
#define COMPUTE_L_m1(ndim,sim) \
    INIT_m1n##ndim\
    "movq %%r13,%4; movq %%r14,%1;"\
    #ndim""#sim"112:\n\t"\
    "testq %4,%4; jz "#ndim""#sim"113f;"\
    KERNEL_k1m1n##ndim(%1)\
    "decq %4; jmp "#ndim""#sim"112b;"\
    #ndim""#sim"113:\n\t"\
    SAVE_L_m1n##ndim "addq $4,%2;"
#define COMPUTE_R_m1(ndim,sim) \
    "salq $2,%%r13;subq %%r13,%0;sarq $2,%%r13;"\
    INIT_m1n##ndim\
    "movq %%r13,%4; leaq (%%r14,%%r12,2),%%r15; addq %%r12,%%r15;"\
    #ndim""#sim"112:\n\t"\
    "testq %4,%4; jz "#ndim""#sim"113f;"\
    KERNEL_k1m1n##ndim(%%r15)\
    "decq %4; jmp "#ndim""#sim"112b;"\
    #ndim""#sim"113:\n\t"\
    SAVE_R_m1n##ndim
#define COMPUTE_m1_n1  COMPUTE_L_m1(1,99899)
#define COMPUTE_m1_n2  COMPUTE_L_m1(2,99899)
#define COMPUTE_m1_n4  COMPUTE_L_m1(4,99899)
#define COMPUTE_m1_n8  COMPUTE_L_m1(8,99899)
#define COMPUTE_m1_n12 COMPUTE_L_m1(12,99899)
#define COMPUTE_m1_n16 COMPUTE_L_m1(12,99799) COMPUTE_R_m1(4,99999)
#define COMPUTE_m1_n20 COMPUTE_L_m1(12,99699) COMPUTE_R_m1(8,99999)
#define COMPUTE_m1_n24 COMPUTE_L_m1(12,99599) COMPUTE_R_m1(12,99999)
#define COMPUTE_m1(ndim) COMPUTE_m1_n##ndim

/* %0 = "+r"(a_pointer), %1 = "+r"(b_pointer), %2 = "+r"(c_pointer), %3 = "+r"(ldc_in_bytes), %4 = "+r"(K), %5 = "+r"(ctemp) */
/* %6 = "+r"(&alpha), %7 = "+r"(M) */
/* r11 = m(const), r12 = k << 4(const), r13 = k(const), r14 = b_head_pos(const), r15 = %1 + 3r12 */

#define COMPUTE(ndim) {\
    __asm__ __volatile__(\
    "vbroadcastss (%6),%%zmm0;"\
    "movq %4,%%r13; movq %4,%%r12; salq $4,%%r12; movq %1,%%r14; movq %7,%%r11;"\
    "cmpq $16,%7;jb 33101"#ndim"f;"\
    "33109"#ndim":\n\t"\
    COMPUTE_m16(ndim)\
    "subq $16,%7;cmpq $16,%7;jnb 33109"#ndim"b;"\
    "33101"#ndim":\n\t"\
    "cmpq $8,%7;jb 33102"#ndim"f;"\
    COMPUTE_m8(ndim)\
    "subq $8,%7;"\
    "33102"#ndim":\n\t"\
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
    :"+r"(a_pointer),"+r"(b_pointer),"+r"(c_pointer),"+r"(ldc_in_bytes),"+r"(K),"+r"(ctemp),"+r"(alp),"+r"(M)\
    ::"r11","r12","r13","r14","r15","zmm0","zmm1","zmm2","zmm3","zmm4","zmm5","zmm6","zmm7","zmm8","zmm9","zmm10","zmm11","zmm12","zmm13","zmm14",\
    "zmm15","zmm16","zmm17","zmm18","zmm19","zmm20","zmm21","zmm22","zmm23","zmm24","zmm25","zmm26","zmm27","zmm28","zmm29","zmm30","zmm31",\
    "cc","memory");\
    a_pointer -= M * K; b_pointer += ndim * K;c_pointer += LDC * ndim - M;\
}
int __attribute__ ((noinline))
CNAME(BLASLONG m, BLASLONG n, BLASLONG k, float alpha, float * __restrict__ A, float * __restrict__ B, float * __restrict__ C, BLASLONG LDC)
{
    if(m==0||n==0||k==0||alpha==(float)0.0) return 0;
    int64_t ldc_in_bytes = (int64_t)LDC * sizeof(float);float ALPHA = alpha;
    int64_t M = (int64_t)m, K = (int64_t)k;
    BLASLONG n_count = n;
    float *a_pointer = A,*b_pointer = B,*c_pointer = C,*ctemp = C,*alp = &ALPHA;
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

#include <immintrin.h>
/* codes below are copied from the sgemm kernel written by Arjan van der Ven */

/*
 * "Direct sgemm" code. This code operates directly on the inputs and outputs
 * of the sgemm call, avoiding the copies, memory realignments and threading,
 * and only supports alpha = 1 and beta = 0.
 * This is a common case and provides value for relatively small matrixes.
 * For larger matrixes the "regular" sgemm code is superior, there the cost of
 * copying/shuffling the B matrix really pays off.
 */



#define DECLARE_RESULT_512(N,M) __m512 result##N##M = _mm512_setzero_ps()
#define BROADCAST_LOAD_A_512(N,M) __m512 Aval##M = _mm512_broadcastss_ps(_mm_load_ss(&A[k  + strideA * (i+M)]))
#define LOAD_B_512(N,M)  __m512 Bval##N = _mm512_loadu_ps(&B[strideB * k + j + (N*16)])
#define MATMUL_512(N,M)  result##N##M = _mm512_fmadd_ps(Aval##M, Bval##N , result##N##M)
#define STORE_512(N,M) _mm512_storeu_ps(&R[(i+M) * strideR + j+(N*16)], result##N##M)


#define DECLARE_RESULT_256(N,M) __m256 result##N##M = _mm256_setzero_ps()
#define BROADCAST_LOAD_A_256(N,M) __m256 Aval##M = _mm256_broadcastss_ps(_mm_load_ss(&A[k  + strideA * (i+M)]))
#define LOAD_B_256(N,M)  __m256 Bval##N = _mm256_loadu_ps(&B[strideB * k + j + (N*8)])
#define MATMUL_256(N,M)  result##N##M = _mm256_fmadd_ps(Aval##M, Bval##N , result##N##M)
#define STORE_256(N,M) _mm256_storeu_ps(&R[(i+M) * strideR + j+(N*8)], result##N##M)

#define DECLARE_RESULT_128(N,M) __m128 result##N##M = _mm_setzero_ps()
#define BROADCAST_LOAD_A_128(N,M) __m128 Aval##M = _mm_broadcastss_ps(_mm_load_ss(&A[k  + strideA * (i+M)]))
#define LOAD_B_128(N,M)  __m128 Bval##N = _mm_loadu_ps(&B[strideB * k + j + (N*4)])
#define MATMUL_128(N,M)  result##N##M = _mm_fmadd_ps(Aval##M, Bval##N , result##N##M)
#define STORE_128(N,M) _mm_storeu_ps(&R[(i+M) * strideR + j+(N*4)], result##N##M)

#define DECLARE_RESULT_SCALAR(N,M) float result##N##M = 0;
#define BROADCAST_LOAD_A_SCALAR(N,M) float Aval##M = A[k + strideA * (i + M)];
#define LOAD_B_SCALAR(N,M)  float Bval##N  = B[k * strideB + j + N];
#define MATMUL_SCALAR(N,M) result##N##M +=  Aval##M * Bval##N;
#define STORE_SCALAR(N,M)  R[(i+M) * strideR + j + N] = result##N##M;

int sgemm_kernel_direct_performant(BLASLONG M, BLASLONG N, BLASLONG K)
{
	int mnk = M * N * K;
	/* large matrixes -> not performant */
	if (mnk >= 28 * 512 * 512)
		return 0;

	/*
	 * if the B matrix is not a nice multiple if 4 we get many unaligned accesses,
	 * and the regular sgemm copy/realignment of data pays off much quicker
	 */
	if ((N & 3) != 0 && (mnk >= 8 * 512 * 512))
		return 0;

#ifdef SMP
	/* if we can run multithreaded, the threading changes the based threshold */
	if (mnk > 2 * 350 * 512 && num_cpu_avail(3)> 1)
		return 0;
#endif

	return 1;
}



void sgemm_kernel_direct (BLASLONG M, BLASLONG N, BLASLONG K, float * __restrict A, BLASLONG strideA, float * __restrict B, BLASLONG strideB , float * __restrict R, BLASLONG strideR)
{
	int i, j, k;

        int m4 = M & ~3;
	int m2 = M & ~1;

	int n64 = N & ~63;
	int n32 = N & ~31;
	int n16 = N & ~15;
	int n8 = N & ~7;
	int n4 = N & ~3;
	int n2 = N & ~1;

	i = 0;

	for (i = 0; i < m4; i+=4) {

		for (j = 0; j < n64; j+= 64) {
			k = 0;
			DECLARE_RESULT_512(0, 0);    DECLARE_RESULT_512(1, 0);    			DECLARE_RESULT_512(2, 0);    DECLARE_RESULT_512(3, 0);
			DECLARE_RESULT_512(0, 1);    DECLARE_RESULT_512(1, 1);    			DECLARE_RESULT_512(2, 1);    DECLARE_RESULT_512(3, 1);
			DECLARE_RESULT_512(0, 2);    DECLARE_RESULT_512(1, 2);    			DECLARE_RESULT_512(2, 2);    DECLARE_RESULT_512(3, 2);
			DECLARE_RESULT_512(0, 3);    DECLARE_RESULT_512(1, 3);    			DECLARE_RESULT_512(2, 3);    DECLARE_RESULT_512(3, 3);


			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_512(x, 0);
				BROADCAST_LOAD_A_512(x, 1);
				BROADCAST_LOAD_A_512(x, 2);
				BROADCAST_LOAD_A_512(x, 3);

				LOAD_B_512(0, x);		LOAD_B_512(1, x);			LOAD_B_512(2, x);		LOAD_B_512(3, x);

				MATMUL_512(0, 0);		MATMUL_512(1, 0);			MATMUL_512(2, 0);		MATMUL_512(3, 0);
				MATMUL_512(0, 1);		MATMUL_512(1, 1);			MATMUL_512(2, 1);		MATMUL_512(3, 1);
				MATMUL_512(0, 2);		MATMUL_512(1, 2);			MATMUL_512(2, 2);		MATMUL_512(3, 2);
				MATMUL_512(0, 3);		MATMUL_512(1, 3);			MATMUL_512(2, 3);		MATMUL_512(3, 3);
			}
			STORE_512(0, 0);		STORE_512(1, 0);			STORE_512(2, 0);		STORE_512(3, 0);
			STORE_512(0, 1);		STORE_512(1, 1);			STORE_512(2, 1);		STORE_512(3, 1);
			STORE_512(0, 2);		STORE_512(1, 2);			STORE_512(2, 2);		STORE_512(3, 2);
			STORE_512(0, 3);		STORE_512(1, 3);			STORE_512(2, 3);		STORE_512(3, 3);
		}

		for (; j < n32; j+= 32) {
			DECLARE_RESULT_512(0, 0);    DECLARE_RESULT_512(1, 0);
			DECLARE_RESULT_512(0, 1);    DECLARE_RESULT_512(1, 1);
			DECLARE_RESULT_512(0, 2);    DECLARE_RESULT_512(1, 2);
			DECLARE_RESULT_512(0, 3);    DECLARE_RESULT_512(1, 3);

			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_512(x, 0);
				BROADCAST_LOAD_A_512(x, 1);
				BROADCAST_LOAD_A_512(x, 2);
				BROADCAST_LOAD_A_512(x, 3);

				LOAD_B_512(0, x);		LOAD_B_512(1, x);

				MATMUL_512(0, 0);		MATMUL_512(1, 0);
				MATMUL_512(0, 1);		MATMUL_512(1, 1);
				MATMUL_512(0, 2);		MATMUL_512(1, 2);
				MATMUL_512(0, 3);		MATMUL_512(1, 3);
			}
			STORE_512(0, 0);		STORE_512(1, 0);
			STORE_512(0, 1);		STORE_512(1, 1);
			STORE_512(0, 2);		STORE_512(1, 2);
			STORE_512(0, 3);		STORE_512(1, 3);
		}

		for (; j < n16; j+= 16) {
			DECLARE_RESULT_512(0, 0);
			DECLARE_RESULT_512(0, 1);
			DECLARE_RESULT_512(0, 2);
			DECLARE_RESULT_512(0, 3);

		 	for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_512(x, 0);
				BROADCAST_LOAD_A_512(x, 1);
				BROADCAST_LOAD_A_512(x, 2);
				BROADCAST_LOAD_A_512(x, 3);

				LOAD_B_512(0, x);

				MATMUL_512(0, 0);
				MATMUL_512(0, 1);
				MATMUL_512(0, 2);
				MATMUL_512(0, 3);
			}
			STORE_512(0, 0);
			STORE_512(0, 1);
			STORE_512(0, 2);
			STORE_512(0, 3);
		}

		for (; j < n8; j+= 8) {
			DECLARE_RESULT_256(0, 0);
			DECLARE_RESULT_256(0, 1);
			DECLARE_RESULT_256(0, 2);
			DECLARE_RESULT_256(0, 3);

			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_256(x, 0);
				BROADCAST_LOAD_A_256(x, 1);
				BROADCAST_LOAD_A_256(x, 2);
				BROADCAST_LOAD_A_256(x, 3);

				LOAD_B_256(0, x);

				MATMUL_256(0, 0);
				MATMUL_256(0, 1);
				MATMUL_256(0, 2);
				MATMUL_256(0, 3);
			}
			STORE_256(0, 0);
			STORE_256(0, 1);
			STORE_256(0, 2);
			STORE_256(0, 3);
		}

		for (; j < n4; j+= 4) {
			DECLARE_RESULT_128(0, 0);
			DECLARE_RESULT_128(0, 1);
			DECLARE_RESULT_128(0, 2);
			DECLARE_RESULT_128(0, 3);

			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_128(x, 0);
				BROADCAST_LOAD_A_128(x, 1);
				BROADCAST_LOAD_A_128(x, 2);
				BROADCAST_LOAD_A_128(x, 3);

				LOAD_B_128(0, x);

				MATMUL_128(0, 0);
				MATMUL_128(0, 1);
				MATMUL_128(0, 2);
				MATMUL_128(0, 3);
			}
			STORE_128(0, 0);
			STORE_128(0, 1);
			STORE_128(0, 2);
			STORE_128(0, 3);
		}

		for (; j < n2; j+= 2) {
			DECLARE_RESULT_SCALAR(0, 0);	DECLARE_RESULT_SCALAR(1, 0);
			DECLARE_RESULT_SCALAR(0, 1);	DECLARE_RESULT_SCALAR(1, 1);
			DECLARE_RESULT_SCALAR(0, 2);	DECLARE_RESULT_SCALAR(1, 2);
			DECLARE_RESULT_SCALAR(0, 3);	DECLARE_RESULT_SCALAR(1, 3);

			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_SCALAR(x, 0);
				BROADCAST_LOAD_A_SCALAR(x, 1);
				BROADCAST_LOAD_A_SCALAR(x, 2);
				BROADCAST_LOAD_A_SCALAR(x, 3);

				LOAD_B_SCALAR(0, x);	LOAD_B_SCALAR(1, x);

				MATMUL_SCALAR(0, 0);	MATMUL_SCALAR(1, 0);
				MATMUL_SCALAR(0, 1);	MATMUL_SCALAR(1, 1);
				MATMUL_SCALAR(0, 2);	MATMUL_SCALAR(1, 2);
				MATMUL_SCALAR(0, 3);	MATMUL_SCALAR(1, 3);
			}
			STORE_SCALAR(0, 0);	STORE_SCALAR(1, 0);
			STORE_SCALAR(0, 1);	STORE_SCALAR(1, 1);
			STORE_SCALAR(0, 2);	STORE_SCALAR(1, 2);
			STORE_SCALAR(0, 3);	STORE_SCALAR(1, 3);
		}

		for (; j < N; j++) {
			DECLARE_RESULT_SCALAR(0, 0)
			DECLARE_RESULT_SCALAR(0, 1)
			DECLARE_RESULT_SCALAR(0, 2)
			DECLARE_RESULT_SCALAR(0, 3)

			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_SCALAR(0, 0);
				BROADCAST_LOAD_A_SCALAR(0, 1);
				BROADCAST_LOAD_A_SCALAR(0, 2);
				BROADCAST_LOAD_A_SCALAR(0, 3);

				LOAD_B_SCALAR(0, 0);

				MATMUL_SCALAR(0, 0);
				MATMUL_SCALAR(0, 1);
				MATMUL_SCALAR(0, 2);
				MATMUL_SCALAR(0, 3);
			}
			STORE_SCALAR(0, 0);
			STORE_SCALAR(0, 1);
			STORE_SCALAR(0, 2);
			STORE_SCALAR(0, 3);
		}
	}

	for (; i < m2; i+=2) {
		j = 0;

		for (; j < n64; j+= 64) {
			DECLARE_RESULT_512(0, 0);    DECLARE_RESULT_512(1, 0);    			DECLARE_RESULT_512(2, 0);    DECLARE_RESULT_512(3, 0);
			DECLARE_RESULT_512(0, 1);    DECLARE_RESULT_512(1, 1);    			DECLARE_RESULT_512(2, 1);    DECLARE_RESULT_512(3, 1);


			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_512(x, 0);
				BROADCAST_LOAD_A_512(x, 1);

				LOAD_B_512(0, x);		LOAD_B_512(1, x);			LOAD_B_512(2, x);		LOAD_B_512(3, x);

				MATMUL_512(0, 0);		MATMUL_512(1, 0);			MATMUL_512(2, 0);		MATMUL_512(3, 0);
				MATMUL_512(0, 1);		MATMUL_512(1, 1);			MATMUL_512(2, 1);		MATMUL_512(3, 1);
			}
			STORE_512(0, 0);		STORE_512(1, 0);			STORE_512(2, 0);		STORE_512(3, 0);
			STORE_512(0, 1);		STORE_512(1, 1);			STORE_512(2, 1);		STORE_512(3, 1);
		}

		for (; j < n32; j+= 32) {
			DECLARE_RESULT_512(0, 0);    DECLARE_RESULT_512(1, 0);
			DECLARE_RESULT_512(0, 1);    DECLARE_RESULT_512(1, 1);

			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_512(x, 0);
				BROADCAST_LOAD_A_512(x, 1);

				LOAD_B_512(0, x);		LOAD_B_512(1, x);

				MATMUL_512(0, 0);		MATMUL_512(1, 0);
				MATMUL_512(0, 1);		MATMUL_512(1, 1);
			}
			STORE_512(0, 0);		STORE_512(1, 0);
			STORE_512(0, 1);		STORE_512(1, 1);
		}


		for (; j < n16; j+= 16) {
			DECLARE_RESULT_512(0, 0);
			DECLARE_RESULT_512(0, 1);

			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_512(x, 0);
				BROADCAST_LOAD_A_512(x, 1);

				LOAD_B_512(0, x);

				MATMUL_512(0, 0);
				MATMUL_512(0, 1);
			}
			STORE_512(0, 0);
			STORE_512(0, 1);
		}

		for (; j < n8; j+= 8) {
			DECLARE_RESULT_256(0, 0);
			DECLARE_RESULT_256(0, 1);

			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_256(x, 0);
				BROADCAST_LOAD_A_256(x, 1);

				LOAD_B_256(0, x);

				MATMUL_256(0, 0);
				MATMUL_256(0, 1);
			}
			STORE_256(0, 0);
			STORE_256(0, 1);
		}

		for (; j < n4; j+= 4) {
			DECLARE_RESULT_128(0, 0);
			DECLARE_RESULT_128(0, 1);

			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_128(x, 0);
				BROADCAST_LOAD_A_128(x, 1);

				LOAD_B_128(0, x);

				MATMUL_128(0, 0);
				MATMUL_128(0, 1);
			}
			STORE_128(0, 0);
			STORE_128(0, 1);
		}
		for (; j < n2; j+= 2) {
			DECLARE_RESULT_SCALAR(0, 0);	DECLARE_RESULT_SCALAR(1, 0);
			DECLARE_RESULT_SCALAR(0, 1);	DECLARE_RESULT_SCALAR(1, 1);

			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_SCALAR(x, 0);
				BROADCAST_LOAD_A_SCALAR(x, 1);

				LOAD_B_SCALAR(0, x);	LOAD_B_SCALAR(1, x);

				MATMUL_SCALAR(0, 0);	MATMUL_SCALAR(1, 0);
				MATMUL_SCALAR(0, 1);	MATMUL_SCALAR(1, 1);
			}
			STORE_SCALAR(0, 0);	STORE_SCALAR(1, 0);
			STORE_SCALAR(0, 1);	STORE_SCALAR(1, 1);
		}

		for (; j < N; j++) {
			DECLARE_RESULT_SCALAR(0, 0);
			DECLARE_RESULT_SCALAR(0, 1);

			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_SCALAR(0, 0);
				BROADCAST_LOAD_A_SCALAR(0, 1);

				LOAD_B_SCALAR(0, 0);

				MATMUL_SCALAR(0, 0);
				MATMUL_SCALAR(0, 1);
			}
			STORE_SCALAR(0, 0);
			STORE_SCALAR(0, 1);
		}
	}

	for (; i < M; i+=1) {
		j = 0;
		for (; j < n64; j+= 64) {
			DECLARE_RESULT_512(0, 0);    DECLARE_RESULT_512(1, 0);    			DECLARE_RESULT_512(2, 0);    DECLARE_RESULT_512(3, 0);

			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_512(x, 0);
				LOAD_B_512(0, x);		LOAD_B_512(1, x);			LOAD_B_512(2, x);		LOAD_B_512(3, x);
				MATMUL_512(0, 0);		MATMUL_512(1, 0);			MATMUL_512(2, 0);		MATMUL_512(3, 0);
			}
			STORE_512(0, 0);		STORE_512(1, 0);			STORE_512(2, 0);		STORE_512(3, 0);
		}
		for (; j < n32; j+= 32) {
			DECLARE_RESULT_512(0, 0);    DECLARE_RESULT_512(1, 0);

			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_512(x, 0);
				LOAD_B_512(0, x);		LOAD_B_512(1, x);
				MATMUL_512(0, 0);		MATMUL_512(1, 0);
			}
			STORE_512(0, 0);		STORE_512(1, 0);
		}


		for (; j < n16; j+= 16) {
			DECLARE_RESULT_512(0, 0);

			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_512(x, 0);

				LOAD_B_512(0, x);

				MATMUL_512(0, 0);
			}
			STORE_512(0, 0);
		}

		for (; j < n8; j+= 8) {
			DECLARE_RESULT_256(0, 0);

			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_256(x, 0);
				LOAD_B_256(0, x);
				MATMUL_256(0, 0);
			}
			STORE_256(0, 0);
		}

		for (; j < n4; j+= 4) {
			DECLARE_RESULT_128(0, 0);

			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_128(x, 0);
				LOAD_B_128(0, x);
				MATMUL_128(0, 0);
			}
			STORE_128(0, 0);
		}

		for (; j < n2; j+= 2) {
			DECLARE_RESULT_SCALAR(0, 0);	DECLARE_RESULT_SCALAR(1, 0);

			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_SCALAR(x, 0);
				LOAD_B_SCALAR(0, 0);	LOAD_B_SCALAR(1, 0);
				MATMUL_SCALAR(0, 0);	MATMUL_SCALAR(1, 0);
			}
			STORE_SCALAR(0, 0);	STORE_SCALAR(1, 0);
		}

		for (; j < N; j++) {
			DECLARE_RESULT_SCALAR(0, 0);

			for (k = 0; k < K; k++) {
				BROADCAST_LOAD_A_SCALAR(0, 0);
				LOAD_B_SCALAR(0, 0);
				MATMUL_SCALAR(0, 0);
			}
			STORE_SCALAR(0, 0);
		}
	}
}
