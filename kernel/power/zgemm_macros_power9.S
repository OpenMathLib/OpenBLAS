/***************************************************************************
Copyright (c) 2013-2019, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

#if   defined(NN) || defined(NT) || defined(TN) || defined(TT)

    #define XSFADD_R1   xsadddp
    #define XSFADD_R2   xssubdp
    #define XSFADD_I1   xsadddp
    #define XSFADD_I2   xsadddp

#elif  defined(CN) || defined(CT) || defined(RN) || defined(RT)

    #define XSFADD_R1   xsadddp
    #define XSFADD_R2   xsadddp
    #define XSFADD_I1   xssubdp
    #define XSFADD_I2   xsadddp

#elif  defined(NC) || defined(TC) || defined(NR) || defined(TR)

    #define XSFADD_R1   xsadddp
    #define XSFADD_R2   xsadddp
    #define XSFADD_I1   xsadddp
    #define XSFADD_I2   xssubdp

#else       // CC || CR || RC || RR

    #define XSFADD_R1   xsadddp
    #define XSFADD_R2   xssubdp
    #define XSFADD_I1   xssubdp
    #define XSFADD_I2   xssubdp

#endif

.macro AGGREGATE_INTO_COMPLEX  FIRST_V, SECOND_V, OUTPUT_V
     AGGREGATE_INTO_COMPLEX_INNER \FIRST_V, \SECOND_V, \OUTPUT_V, vs0,vs1,vs2,vs3,vs4,vs5,vs6,vs7
.endm

.macro AGGREGATE_INTO_COMPLEX_INNER  FIRST_V, SECOND_V, OUTPUT_V ,TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,TEMP6,TEMP7,TEMP8
    xxlxor      \TEMP1, \TEMP1, \TEMP1
    xxlxor      \TEMP2, \TEMP2, \TEMP2
 
    xxswapd     \SECOND_V,  \SECOND_V           //   imagA*imagB, realA*imagB -> realA*imagB, imagA*imagB 

    XSFADD_I1   \TEMP2, \TEMP2, \FIRST_V        // realA*imagB
    XSFADD_I2   \TEMP2, \TEMP2, \SECOND_V       // imagA*realB

    xxswapd     \FIRST_V,   \FIRST_V            //imagA*realB, realA*realB -> realA*realB, imagA*realB   
    xxswapd     \SECOND_V,  \SECOND_V           //  reverse to original imagA*imagB, realA*imagB 

    XSFADD_R1   \TEMP1, \TEMP1, \FIRST_V        // realA*realB
    XSFADD_R2   \TEMP1, \TEMP1, \SECOND_V       // imagA*imagB

    xsmuldp     \TEMP3, \TEMP2, alpha_i     // imag*alpha_i
    xsmuldp     \TEMP4, \TEMP2, alpha_r     // imag*alpha_r 
    xsmuldp     \TEMP5, \TEMP1, alpha_r     // real*alpha_r 
    xsmuldp     \TEMP6, \TEMP1, alpha_i     // real*alpha_i

    xssubdp     \TEMP7, \TEMP5, \TEMP3      // real*alpha_r - imag*alpha_i
    xsadddp     \TEMP8, \TEMP6, \TEMP4      // real*alpha_i + imag*alpha_r
    xxpermdi    \OUTPUT_V,  \TEMP8, \TEMP7, 0   // merge real and imag part
.endm

/**********************************************************************************************
* Macros for N=2 and M=8
**********************************************************************************************/

#define unit_size 16
#define DISP32(ind,disp) (ind*unit_size*32+disp)
#define DISP16(ind,disp) (ind*unit_size*16+disp)
#define DISP8(ind,disp) (ind*unit_size*8+disp)
#define DISP4(ind,disp) (ind*unit_size*4+disp)
#define DISP2(ind,disp) (ind*unit_size*2+disp)
#define DISP1(ind,disp) (ind*unit_size+disp)

.macro Zero2x8
    xxlxor      vs32,   vs32,   vs32
    xxlxor      vs33,   vs33,   vs33
    xxlxor      vs34,   vs34,   vs34
    xxlxor      vs35,   vs35,   vs35
    xxlxor      vs36,   vs36,   vs36
    xxlxor      vs37,   vs37,   vs37
    xxlxor      vs38,   vs38,   vs38
    xxlxor      vs39,   vs39,   vs39
    xxlxor      vs40,   vs40,   vs40
    xxlxor      vs41,   vs41,   vs41
    xxlxor      vs42,   vs42,   vs42
    xxlxor      vs43,   vs43,   vs43
    xxlxor      vs44,   vs44,   vs44
    xxlxor      vs45,   vs45,   vs45
    xxlxor      vs46,   vs46,   vs46
    xxlxor      vs47,   vs47,   vs47
    xxlxor      vs48,   vs48,   vs48
    xxlxor      vs49,   vs49,   vs49
    xxlxor      vs50,   vs50,   vs50
    xxlxor      vs51,   vs51,   vs51 
    xxlxor      vs52,   vs52,   vs52
    xxlxor      vs53,   vs53,   vs53
    xxlxor      vs54,   vs54,   vs54
    xxlxor      vs55,   vs55,   vs55 
    xxlxor      vs56,   vs56,   vs56
    xxlxor      vs57,   vs57,   vs57
    xxlxor      vs58,   vs58,   vs58
    xxlxor      vs59,   vs59,   vs59 
    xxlxor      vs60,   vs60,   vs60
    xxlxor      vs61,   vs61,   vs61
    xxlxor      vs62,   vs62,   vs62
    xxlxor      vs63,   vs63,   vs63    
.endm

.macro LOAD2x8 Zero

    lxv     vs16,   0(BO)       // load real part from B
    lxv     vs17,   16(BO)      // load imag part from B
    lxv     vs18,   32(BO)      // load real part from B
    lxv     vs19,   48(BO)      // load imag part from B

    lxv     vs0,    0(AO)       // load real,imag from A
    lxv     vs1,    16(AO)      // load real,imag from A
    lxv     vs2,    32(AO)      // load real,imag from A
    lxv     vs3,    48(AO)      // load real,imag from A

    lxv     vs4,    64(AO)      // load real,imag from A
    lxv     vs5,    80(AO)      // load real,imag from A
    lxv     vs6,    96(AO)      // load real,imag from A
    lxv     vs7,    112(AO)     // load real,imag from A

.if \Zero==1
    Zero2x8 
.endif

.endm

.macro END2x8_NORMAL
   END2x8 AO,BO,128,64
.endm

.macro END2x8   AREG, BREG, OffsetA, OffsetB

.if \OffsetB != 0 
    addi        \BREG, \BREG, \OffsetB 
.endif
.if \OffsetA != 0 
    addi        \AREG, \AREG, \OffsetA 
.endif  

    xvmaddadp   vs32,   vs0,    vs16        // real*real, imag*real
    xvmaddadp   vs33,   vs0,    vs17        // real*imag, imag*imag
    xvmaddadp   vs34,   vs1,    vs16        // real*real, imag*real
    xvmaddadp   vs35,   vs1,    vs17        // real*imag, imag*imag
    xvmaddadp   vs36,   vs2,    vs16        // real*real, imag*real
    xvmaddadp   vs37,   vs2,    vs17        // real*imag, imag*imag
    xvmaddadp   vs38,   vs3,    vs16        // real*real, imag*real
    xvmaddadp   vs39,   vs3,    vs17        // real*imag, imag*imag
    xvmaddadp   vs40,   vs4,    vs16        // real*real, imag*real
    xvmaddadp   vs41,   vs4,    vs17        // real*imag, imag*imag
    xvmaddadp   vs42,   vs5,    vs16        // real*real, imag*real
    xvmaddadp   vs43,   vs5,    vs17        // real*imag, imag*imag
    xvmaddadp   vs44,   vs6,    vs16        // real*real, imag*real
    xvmaddadp   vs45,   vs6,    vs17        // real*imag, imag*imag
    xvmaddadp   vs46,   vs7,    vs16        // real*real, imag*real
    xvmaddadp   vs47,   vs7,    vs17        // real*imag, imag*imag

    xvmaddadp   vs48,   vs0,    vs18        // real*real, imag*real
    xvmaddadp   vs49,   vs0,    vs19        // real*imag, imag*imag
    xvmaddadp   vs50,   vs1,    vs18        // real*real, imag*real
    xvmaddadp   vs51,   vs1,    vs19        // real*imag, imag*imag
    xvmaddadp   vs52,   vs2,    vs18        // real*real, imag*real
    xvmaddadp   vs53,   vs2,    vs19        // real*imag, imag*imag
    xvmaddadp   vs54,   vs3,    vs18        // real*real, imag*real
    xvmaddadp   vs55,   vs3,    vs19        // real*imag, imag*imag
    xvmaddadp   vs56,   vs4,    vs18        // real*real, imag*real
    xvmaddadp   vs57,   vs4,    vs19        // real*imag, imag*imag
    xvmaddadp   vs58,   vs5,    vs18        // real*real, imag*real
    xvmaddadp   vs59,   vs5,    vs19        // real*imag, imag*imag
    xvmaddadp   vs60,   vs6,    vs18        // real*real, imag*real
    xvmaddadp   vs61,   vs6,    vs19        // real*imag, imag*imag
    xvmaddadp   vs62,   vs7,    vs18        // real*real, imag*real
    xvmaddadp   vs63,   vs7,    vs19        // real*imag, imag*imag

.endm

.macro KERNEL2x8_L      OffsetA,OffsetB, Index,IsLast  
  KERNEL2x8_2   AO,BO,   \OffsetA,\OffsetB, \Index,\IsLast ,0
.endm

.macro KERNEL2x8_E      OffsetA,OffsetB, Index,IsLast  
  KERNEL2x8_2   AO,BO,   \OffsetA,\OffsetB, \Index,\IsLast ,1
.endm

.macro KERNEL2x8_2   AREG,BREG,   OffsetA,OffsetB, Index,IsLast ,Complete

    lxv     vs8,     DISP16(\Index, 0 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs9,     DISP16(\Index,16 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs10,    DISP16(\Index,32 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs11,    DISP16(\Index,48 + \OffsetA)(\AREG)        // load real,imag from A 

    lxv     vs12,    DISP16(\Index, 64 + \OffsetA)(\AREG)       // load real,imag from A
    lxv     vs13,    DISP16(\Index,64+16 + \OffsetA)(\AREG)     // load real,imag from A
    lxv     vs14,    DISP16(\Index,64+32 + \OffsetA)(\AREG)     // load real,imag from A
    lxv     vs15,    DISP16(\Index,64+48 + \OffsetA)(\AREG)     // load real,imag from A

lxv     vs20,   DISP8(\Index,   0+\OffsetB)(\BREG)      // load real part from B
    lxv     vs21,   DISP8(\Index,16+\OffsetB)(\BREG)        // load imag part from B
    lxv     vs22,   DISP8(\Index,32+\OffsetB)(\BREG)        // load real part from B
    lxv     vs23,   DISP8(\Index,48+\OffsetB)(\BREG)        // load imag part from B

    xvmaddadp   vs32,   vs0,    vs16        // real*real, imag*real
    xvmaddadp   vs33,   vs0,    vs17        // real*imag, imag*imag
    xvmaddadp   vs34,   vs1,    vs16        // real*real, imag*real
    xvmaddadp   vs35,   vs1,    vs17        // real*imag, imag*imag
    xvmaddadp   vs36,   vs2,    vs16        // real*real, imag*real
    xvmaddadp   vs37,   vs2,    vs17        // real*imag, imag*imag
    xvmaddadp   vs38,   vs3,    vs16        // real*real, imag*real
    xvmaddadp   vs39,   vs3,    vs17        // real*imag, imag*imag
    xvmaddadp   vs40,   vs4,    vs16        // real*real, imag*real
    xvmaddadp   vs41,   vs4,    vs17        // real*imag, imag*imag
    xvmaddadp   vs42,   vs5,    vs16        // real*real, imag*real
    xvmaddadp   vs43,   vs5,    vs17        // real*imag, imag*imag
    xvmaddadp   vs44,   vs6,    vs16        // real*real, imag*real
    xvmaddadp   vs45,   vs6,    vs17        // real*imag, imag*imag
    xvmaddadp   vs46,   vs7,    vs16        // real*real, imag*real
    xvmaddadp   vs47,   vs7,    vs17        // real*imag, imag*imag

    xvmaddadp   vs48,   vs0,    vs18        // real*real, imag*real
    xvmaddadp   vs49,   vs0,    vs19        // real*imag, imag*imag
    xvmaddadp   vs50,   vs1,    vs18        // real*real, imag*real
    xvmaddadp   vs51,   vs1,    vs19        // real*imag, imag*imag
    xvmaddadp   vs52,   vs2,    vs18        // real*real, imag*real
    xvmaddadp   vs53,   vs2,    vs19        // real*imag, imag*imag
    xvmaddadp   vs54,   vs3,    vs18        // real*real, imag*real
    xvmaddadp   vs55,   vs3,    vs19        // real*imag, imag*imag
    xvmaddadp   vs56,   vs4,    vs18        // real*real, imag*real
    xvmaddadp   vs57,   vs4,    vs19        // real*imag, imag*imag
    xvmaddadp   vs58,   vs5,    vs18        // real*real, imag*real
    xvmaddadp   vs59,   vs5,    vs19        // real*imag, imag*imag
    xvmaddadp   vs60,   vs6,    vs18        // real*real, imag*real
    xvmaddadp   vs61,   vs6,    vs19        // real*imag, imag*imag
    xvmaddadp   vs62,   vs7,    vs18        // real*real, imag*real
    xvmaddadp   vs63,   vs7,    vs19        // real*imag, imag*imag

.if \Complete==0
    lxv     vs0,     DISP16(\Index,128+ + \OffsetA)(\AREG)      // load real,imag from A
    lxv     vs1,     DISP16(\Index,128+16 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs2,     DISP16(\Index,128+32 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs3,     DISP16(\Index,128+48 + \OffsetA)(\AREG)        // load real,imag from A

    lxv     vs4,     DISP16(\Index, 192 + \OffsetA)(\AREG)      // load real,imag from A
    lxv     vs5,     DISP16(\Index,192 +16 + \OffsetA)(\AREG)       // load real,imag from A
    lxv     vs6,     DISP16(\Index,192 +32 + \OffsetA)(\AREG)       // load real,imag from A
    lxv     vs7,     DISP16(\Index,192 +48 + \OffsetA)(\AREG)       // load real,imag from A

    lxv     vs16,   DISP8(\Index,   64+\OffsetB)(\BREG)     // load real part from B
    lxv     vs17,   DISP8(\Index,64+16+\OffsetB)(\BREG)     // load imag part from B
    lxv     vs18,   DISP8(\Index,64+32+\OffsetB)(\BREG)     // load real part from B
    lxv     vs19,   DISP8(\Index,64+48+\OffsetB)(\BREG)     // load imag part from B
.endif

.if \IsLast==1  
.if \Complete==1
    addi        \AREG, \AREG, DISP16(\Index,128+\OffsetA)
    addi        \BREG, \BREG,  DISP8(\Index,64+\OffsetB)
.else 
    addi        \AREG, \AREG, DISP16(\Index,256)
    addi        \BREG, \BREG,  DISP8(\Index,128)
.endif
.endif  

    xvmaddadp   vs32,   vs8,    vs20        // real*real, imag*real
    xvmaddadp   vs33,   vs8,    vs21        // real*imag, imag*imag
    xvmaddadp   vs34,   vs9,    vs20        // real*real, imag*real
    xvmaddadp   vs35,   vs9,    vs21        // real*imag, imag*imag
    xvmaddadp   vs36,   vs10,   vs20        // real*real, imag*real
    xvmaddadp   vs37,   vs10,   vs21        // real*imag, imag*imag
    xvmaddadp   vs38,   vs11,   vs20        // real*real, imag*real
    xvmaddadp   vs39,   vs11,   vs21        // real*imag, imag*imag
    xvmaddadp   vs40,   vs12,   vs20        // real*real, imag*real
    xvmaddadp   vs41,   vs12,   vs21        // real*imag, imag*imag
    xvmaddadp   vs42,   vs13,   vs20        // real*real, imag*real
    xvmaddadp   vs43,   vs13,   vs21        // real*imag, imag*imag
    xvmaddadp   vs44,   vs14,   vs20        // real*real, imag*real
    xvmaddadp   vs45,   vs14,   vs21        // real*imag, imag*imag
    xvmaddadp   vs46,   vs15,   vs20        // real*real, imag*real
    xvmaddadp   vs47,   vs15,   vs21        // real*imag, imag*imag

    xvmaddadp   vs48,   vs8,    vs22        // real*real, imag*real
    xvmaddadp   vs49,   vs8,    vs23        // real*imag, imag*imag
    xvmaddadp   vs50,   vs9,    vs22        // real*real, imag*real
    xvmaddadp   vs51,   vs9,    vs23        // real*imag, imag*imag
    xvmaddadp   vs52,   vs10,   vs22        // real*real, imag*real
    xvmaddadp   vs53,   vs10,   vs23        // real*imag, imag*imag
    xvmaddadp   vs54,   vs11,   vs22        // real*real, imag*real
    xvmaddadp   vs55,   vs11,   vs23        // real*imag, imag*imag
    xvmaddadp   vs56,   vs12,   vs22        // real*real, imag*real
    xvmaddadp   vs57,   vs12,   vs23        // real*imag, imag*imag
    xvmaddadp   vs58,   vs13,   vs22        // real*real, imag*real
    xvmaddadp   vs59,   vs13,   vs23        // real*imag, imag*imag
    xvmaddadp   vs60,   vs14,   vs22        // real*real, imag*real
    xvmaddadp   vs61,   vs14,   vs23        // real*imag, imag*imag
    xvmaddadp   vs62,   vs15,   vs22        // real*real, imag*real
    xvmaddadp   vs63,   vs15,   vs23        // real*imag, imag*imag

.endm

.macro KERNEL2x8 
  LOAD2x8 0
  END2x8  AO, BO, 128,64 
.endm

.macro SAVE2x8

    mr      T1, CO
    addi        T2, T1, 64

#ifndef TRMMKERNEL

    lxv     vs16,   0(T1)
    lxv     vs17,   16(T1)
    lxv     vs18,   32(T1)
    lxv     vs19,   48(T1)
    lxv     vs20,   0(T2)
    lxv     vs21,   16(T2)
    lxv     vs22,   32(T2)
    lxv     vs23,   48(T2)

#endif

   AGGREGATE_INTO_COMPLEX vs32,vs33,vs8
   AGGREGATE_INTO_COMPLEX vs34,vs35,vs9
   AGGREGATE_INTO_COMPLEX vs36,vs37,vs10
   AGGREGATE_INTO_COMPLEX vs38,vs39,vs11
   AGGREGATE_INTO_COMPLEX vs40,vs41,vs12
   AGGREGATE_INTO_COMPLEX vs42,vs43,vs13
   AGGREGATE_INTO_COMPLEX vs44,vs45,vs14
   AGGREGATE_INTO_COMPLEX vs46,vs47,vs15

#ifndef TRMMKERNEL

    xvadddp     vs8,    vs8,    vs16
    xvadddp     vs9,    vs9,    vs17
    xvadddp     vs10,   vs10,   vs18
    xvadddp     vs11,   vs11,   vs19
    xvadddp     vs12,   vs12,   vs20
    xvadddp     vs13,   vs13,   vs21
    xvadddp     vs14,   vs14,   vs22
    xvadddp     vs15,   vs15,   vs23

#endif

    stxv        vs8,    0(T1)
    stxv        vs9,    16(T1)
    stxv        vs10,   32(T1)
    stxv        vs11,   48(T1)
    stxv        vs12,   0(T2)
    stxv        vs13,   16(T2)
    stxv        vs14,   32(T2)
    stxv        vs15,   48(T2)

    add     T1, T1, LDC
    add     T2, T2, LDC

#ifndef TRMMKERNEL

    lxv     vs16,   0(T1)
    lxv     vs17,   16(T1)
    lxv     vs18,   32(T1)
    lxv     vs19,   48(T1)
    lxv     vs20,   0(T2)
    lxv     vs21,   16(T2)
    lxv     vs22,   32(T2)
    lxv     vs23,   48(T2)

#endif

   AGGREGATE_INTO_COMPLEX vs48,vs49,vs8
   AGGREGATE_INTO_COMPLEX vs50,vs51,vs9
   AGGREGATE_INTO_COMPLEX vs52,vs53,vs10
   AGGREGATE_INTO_COMPLEX vs54,vs55,vs11
   AGGREGATE_INTO_COMPLEX vs56,vs57,vs12
   AGGREGATE_INTO_COMPLEX vs58,vs59,vs13
   AGGREGATE_INTO_COMPLEX vs60,vs61,vs14
   AGGREGATE_INTO_COMPLEX vs62,vs63,vs15

#ifndef TRMMKERNEL

    xvadddp     vs8,    vs8,    vs16
    xvadddp     vs9,    vs9,    vs17
    xvadddp     vs10,   vs10,   vs18
    xvadddp     vs11,   vs11,   vs19
    xvadddp     vs12,   vs12,   vs20
    xvadddp     vs13,   vs13,   vs21
    xvadddp     vs14,   vs14,   vs22
    xvadddp     vs15,   vs15,   vs23

#endif

    stxv        vs8,    0(T1)
    stxv        vs9,    16(T1)
    stxv        vs10,   32(T1)
    stxv        vs11,   48(T1)
    stxv        vs12,   0(T2)
    stxv        vs13,   16(T2)
    stxv        vs14,   32(T2)
    stxv        vs15,   48(T2)
 
    addi        CO, CO, 128

.endm

/**********************************************************************************************
* Macros for N=2 and M=4
**********************************************************************************************/

.macro Zero2x4
    xxlxor      vs32,   vs32,   vs32
    xxlxor      vs33,   vs33,   vs33
    xxlxor      vs34,   vs34,   vs34
    xxlxor      vs35,   vs35,   vs35
    xxlxor      vs36,   vs36,   vs36
    xxlxor      vs37,   vs37,   vs37
    xxlxor      vs38,   vs38,   vs38
    xxlxor      vs39,   vs39,   vs39
    xxlxor      vs40,   vs40,   vs40
    xxlxor      vs41,   vs41,   vs41
    xxlxor      vs42,   vs42,   vs42
    xxlxor      vs43,   vs43,   vs43
    xxlxor      vs44,   vs44,   vs44
    xxlxor      vs45,   vs45,   vs45
    xxlxor      vs46,   vs46,   vs46
    xxlxor      vs47,   vs47,   vs47 
.endm

.macro LOAD2x4 Zero

    lxv     vs16,   0(BO)       // load real part from B
    lxv     vs17,   16(BO)      // load imag part from B
    lxv     vs18,   32(BO)      // load real part from B
    lxv     vs19,   48(BO)      // load imag part from B

    lxv     vs0,    0(AO)       // load real,imag from A
    lxv     vs1,    16(AO)      // load real,imag from A
    lxv     vs2,    32(AO)      // load real,imag from A
    lxv     vs3,    48(AO)      // load real,imag from A
 
.if \Zero==1
    Zero2x4 
.endif

.endm

.macro END2x4_NORMAL
   END2x4 AO,BO,64,64
.endm

.macro END2x4   AREG, BREG, OffsetA, OffsetB

.if \OffsetB != 0 
    addi        \BREG, \BREG, \OffsetB 
.endif
.if \OffsetA != 0 
    addi        \AREG, \AREG, \OffsetA 
.endif  

    xvmaddadp   vs32,   vs0,    vs16        // real*real, imag*real
    xvmaddadp   vs33,   vs0,    vs17        // real*imag, imag*imag
    xvmaddadp   vs34,   vs1,    vs16        // real*real, imag*real
    xvmaddadp   vs35,   vs1,    vs17        // real*imag, imag*imag
    xvmaddadp   vs36,   vs2,    vs16        // real*real, imag*real
    xvmaddadp   vs37,   vs2,    vs17        // real*imag, imag*imag
    xvmaddadp   vs38,   vs3,    vs16        // real*real, imag*real
    xvmaddadp   vs39,   vs3,    vs17        // real*imag, imag*imag

    xvmaddadp   vs40,   vs0,    vs18        // real*real, imag*real
    xvmaddadp   vs41,   vs0,    vs19        // real*imag, imag*imag
    xvmaddadp   vs42,   vs1,    vs18        // real*real, imag*real
    xvmaddadp   vs43,   vs1,    vs19        // real*imag, imag*imag
    xvmaddadp   vs44,   vs2,    vs18        // real*real, imag*real
    xvmaddadp   vs45,   vs2,    vs19        // real*imag, imag*imag
    xvmaddadp   vs46,   vs3,    vs18        // real*real, imag*real
    xvmaddadp   vs47,   vs3,    vs19        // real*imag, imag*imag

.endm

.macro KERNEL2x4_L      OffsetA,OffsetB, Index,IsLast  
  KERNEL2x4_2   AO,BO,   \OffsetA,\OffsetB, \Index,\IsLast ,0
.endm

.macro KERNEL2x4_E      OffsetA,OffsetB, Index,IsLast  
  KERNEL2x4_2   AO,BO,   \OffsetA,\OffsetB, \Index,\IsLast ,1
.endm

.macro KERNEL2x4_2   AREG,BREG,   OffsetA,OffsetB, Index,IsLast ,Complete

    lxv     vs8,     DISP8(\Index, 0 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs9,     DISP8(\Index,16 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs10,    DISP8(\Index,32 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs11,    DISP8(\Index,48 + \OffsetA)(\AREG)        // load real,imag from A

lxv     vs20,   DISP8(\Index,   0+\OffsetB)(\BREG)      // load real part from B
    lxv     vs21,   DISP8(\Index,16+\OffsetB)(\BREG)        // load imag part from B
    lxv     vs22,   DISP8(\Index,32+\OffsetB)(\BREG)        // load real part from B
    lxv     vs23,   DISP8(\Index,48+\OffsetB)(\BREG)        // load imag part from B

    xvmaddadp   vs32,   vs0,    vs16        // real*real, imag*real
    xvmaddadp   vs33,   vs0,    vs17        // real*imag, imag*imag
    xvmaddadp   vs34,   vs1,    vs16        // real*real, imag*real
    xvmaddadp   vs35,   vs1,    vs17        // real*imag, imag*imag
    xvmaddadp   vs36,   vs2,    vs16        // real*real, imag*real
    xvmaddadp   vs37,   vs2,    vs17        // real*imag, imag*imag
    xvmaddadp   vs38,   vs3,    vs16        // real*real, imag*real
    xvmaddadp   vs39,   vs3,    vs17        // real*imag, imag*imag

    xvmaddadp   vs40,   vs0,    vs18        // real*real, imag*real
    xvmaddadp   vs41,   vs0,    vs19        // real*imag, imag*imag
    xvmaddadp   vs42,   vs1,    vs18        // real*real, imag*real
    xvmaddadp   vs43,   vs1,    vs19        // real*imag, imag*imag
    xvmaddadp   vs44,   vs2,    vs18        // real*real, imag*real
    xvmaddadp   vs45,   vs2,    vs19        // real*imag, imag*imag
    xvmaddadp   vs46,   vs3,    vs18        // real*real, imag*real
    xvmaddadp   vs47,   vs3,    vs19        // real*imag, imag*imag

.if \Complete==0
    lxv     vs0,     DISP8(\Index,64+  \OffsetA)(\AREG)      // load real,imag from A
    lxv     vs1,     DISP8(\Index,64+16 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs2,     DISP8(\Index,64+32 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs3,     DISP8(\Index,64+48 + \OffsetA)(\AREG)        // load real,imag from A 

    lxv     vs16,   DISP8(\Index,   64+\OffsetB)(\BREG)     // load real part from B
    lxv     vs17,   DISP8(\Index,64+16+\OffsetB)(\BREG)     // load imag part from B
    lxv     vs18,   DISP8(\Index,64+32+\OffsetB)(\BREG)     // load real part from B
    lxv     vs19,   DISP8(\Index,64+48+\OffsetB)(\BREG)     // load imag part from B
.endif

.if \IsLast==1  
.if \Complete==1
    addi        \AREG, \AREG,  DISP8(\Index,64+\OffsetA)
    addi        \BREG, \BREG,  DISP8(\Index,64+\OffsetB)
.else 
    addi        \AREG, \AREG,  DISP8(\Index,128)
    addi        \BREG, \BREG,  DISP8(\Index,128)
.endif
.endif  

    xvmaddadp   vs32,   vs8,    vs20        // real*real, imag*real
    xvmaddadp   vs33,   vs8,    vs21        // real*imag, imag*imag
    xvmaddadp   vs34,   vs9,    vs20        // real*real, imag*real
    xvmaddadp   vs35,   vs9,    vs21        // real*imag, imag*imag
    xvmaddadp   vs36,   vs10,   vs20        // real*real, imag*real
    xvmaddadp   vs37,   vs10,   vs21        // real*imag, imag*imag
    xvmaddadp   vs38,   vs11,   vs20        // real*real, imag*real
    xvmaddadp   vs39,   vs11,   vs21        // real*imag, imag*imag
 
    xvmaddadp   vs40,   vs8,    vs22        // real*real, imag*real
    xvmaddadp   vs41,   vs8,    vs23        // real*imag, imag*imag
    xvmaddadp   vs42,   vs9,    vs22        // real*real, imag*real
    xvmaddadp   vs43,   vs9,    vs23        // real*imag, imag*imag
    xvmaddadp   vs44,   vs10,   vs22        // real*real, imag*real
    xvmaddadp   vs45,   vs10,   vs23        // real*imag, imag*imag
    xvmaddadp   vs46,   vs11,   vs22        // real*real, imag*real
    xvmaddadp   vs47,   vs11,   vs23        // real*imag, imag*imag

.endm

.macro KERNEL2x4 
  LOAD2x4 0
  END2x4  AO, BO, 64,64 
.endm

.macro SAVE2x4

    mr      T1, CO

#ifndef TRMMKERNEL

    lxv     vs16,   0(T1)
    lxv     vs17,   16(T1)
    lxv     vs18,   32(T1)
    lxv     vs19,   48(T1)

#endif

   AGGREGATE_INTO_COMPLEX vs32,vs33,vs8
   AGGREGATE_INTO_COMPLEX vs34,vs35,vs9
   AGGREGATE_INTO_COMPLEX vs36,vs37,vs10
   AGGREGATE_INTO_COMPLEX vs38,vs39,vs11

#ifndef TRMMKERNEL

    xvadddp     vs8,    vs8,    vs16
    xvadddp     vs9,    vs9,    vs17
    xvadddp     vs10,   vs10,   vs18
    xvadddp     vs11,   vs11,   vs19

#endif

    stxv        vs8,    0(T1)
    stxv        vs9,    16(T1)
    stxv        vs10,   32(T1)
    stxv        vs11,   48(T1)

    add     T1, T1, LDC

#ifndef TRMMKERNEL

    lxv     vs16,   0(T1)
    lxv     vs17,   16(T1)
    lxv     vs18,   32(T1)
    lxv     vs19,   48(T1)

#endif

   AGGREGATE_INTO_COMPLEX vs40,vs41,vs8
   AGGREGATE_INTO_COMPLEX vs42,vs43,vs9
   AGGREGATE_INTO_COMPLEX vs44,vs45,vs10
   AGGREGATE_INTO_COMPLEX vs46,vs47,vs11

#ifndef TRMMKERNEL

    xvadddp     vs8,    vs8,    vs16
    xvadddp     vs9,    vs9,    vs17
    xvadddp     vs10,   vs10,   vs18
    xvadddp     vs11,   vs11,   vs19

#endif

    stxv        vs8,    0(T1)
    stxv        vs9,    16(T1)
    stxv        vs10,   32(T1)
    stxv        vs11,   48(T1)
 
    addi        CO, CO, 64

.endm

/**********************************************************************************************
* Macros for N=2 and M=2
**********************************************************************************************/

.macro Zero2x2
    xxlxor      vs32,   vs32,   vs32
    xxlxor      vs33,   vs33,   vs33
    xxlxor      vs34,   vs34,   vs34
    xxlxor      vs35,   vs35,   vs35
    xxlxor      vs36,   vs36,   vs36
    xxlxor      vs37,   vs37,   vs37
    xxlxor      vs38,   vs38,   vs38
    xxlxor      vs39,   vs39,   vs39 
.endm

.macro LOAD2x2 Zero

    lxv     vs16,   0(BO)       // load real part from B
    lxv     vs17,   16(BO)      // load imag part from B
    lxv     vs18,   32(BO)      // load real part from B
    lxv     vs19,   48(BO)      // load imag part from B

    lxv     vs0,    0(AO)       // load real,imag from A
    lxv     vs1,    16(AO)      // load real,imag from A 
 
.if \Zero==1
    Zero2x2 
.endif

.endm

.macro END2x2_NORMAL
   END2x2 AO,BO,32,64
.endm

.macro END2x2   AREG, BREG, OffsetA, OffsetB

.if \OffsetB != 0 
    addi        \BREG, \BREG, \OffsetB 
.endif
.if \OffsetA != 0 
    addi        \AREG, \AREG, \OffsetA 
.endif  

    xvmaddadp   vs32,   vs0,    vs16        // real*real, imag*real
    xvmaddadp   vs33,   vs0,    vs17        // real*imag, imag*imag
    xvmaddadp   vs34,   vs1,    vs16        // real*real, imag*real
    xvmaddadp   vs35,   vs1,    vs17        // real*imag, imag*imag 

    xvmaddadp   vs36,   vs0,    vs18        // real*real, imag*real
    xvmaddadp   vs37,   vs0,    vs19        // real*imag, imag*imag
    xvmaddadp   vs38,   vs1,    vs18        // real*real, imag*real
    xvmaddadp   vs39,   vs1,    vs19        // real*imag, imag*imag 
  
.endm

.macro KERNEL2x2_L      OffsetA,OffsetB, Index,IsLast  
  KERNEL2x2_2   AO,BO,   \OffsetA,\OffsetB, \Index,\IsLast ,0
.endm

.macro KERNEL2x2_E      OffsetA,OffsetB, Index,IsLast  
  KERNEL2x2_2   AO,BO,   \OffsetA,\OffsetB, \Index,\IsLast ,1
.endm

.macro KERNEL2x2_2   AREG,BREG,   OffsetA,OffsetB, Index,IsLast ,Complete

    lxv     vs8,     DISP4(\Index, 0 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs9,     DISP4(\Index,16 + \OffsetA)(\AREG)        // load real,imag from A

lxv     vs20,   DISP8(\Index,   0+\OffsetB)(\BREG)      // load real part from B
    lxv     vs21,   DISP8(\Index,16+\OffsetB)(\BREG)        // load imag part from B
    lxv     vs22,   DISP8(\Index,32+\OffsetB)(\BREG)        // load real part from B
    lxv     vs23,   DISP8(\Index,48+\OffsetB)(\BREG)        // load imag part from B

    xvmaddadp   vs32,   vs0,    vs16        // real*real, imag*real
    xvmaddadp   vs33,   vs0,    vs17        // real*imag, imag*imag
    xvmaddadp   vs34,   vs1,    vs16        // real*real, imag*real
    xvmaddadp   vs35,   vs1,    vs17        // real*imag, imag*imag 

    xvmaddadp   vs36,   vs0,    vs18        // real*real, imag*real
    xvmaddadp   vs37,   vs0,    vs19        // real*imag, imag*imag
    xvmaddadp   vs38,   vs1,    vs18        // real*real, imag*real
    xvmaddadp   vs39,   vs1,    vs19        // real*imag, imag*imag 

.if \Complete==0
    lxv     vs0,     DISP4(\Index,32 + \OffsetA)(\AREG)      // load real,imag from A
    lxv     vs1,     DISP4(\Index,48+ \OffsetA)(\AREG)        // load real,imag from A 

    lxv     vs16,   DISP8(\Index,   64+\OffsetB)(\BREG)     // load real part from B
    lxv     vs17,   DISP8(\Index,64+16+\OffsetB)(\BREG)     // load imag part from B
    lxv     vs18,   DISP8(\Index,64+32+\OffsetB)(\BREG)     // load real part from B
    lxv     vs19,   DISP8(\Index,64+48+\OffsetB)(\BREG)     // load imag part from B
.endif

.if \IsLast==1  
.if \Complete==1
    addi        \AREG, \AREG,  DISP4(\Index,32+\OffsetA)
    addi        \BREG, \BREG,  DISP8(\Index,64+\OffsetB)
.else 
    addi        \AREG, \AREG,  DISP4(\Index,64)
    addi        \BREG, \BREG,  DISP8(\Index,128)
.endif
.endif  

    xvmaddadp   vs32,   vs8,    vs20        // real*real, imag*real
    xvmaddadp   vs33,   vs8,    vs21        // real*imag, imag*imag
    xvmaddadp   vs34,   vs9,    vs20        // real*real, imag*real
    xvmaddadp   vs35,   vs9,    vs21        // real*imag, imag*imag 
 
    xvmaddadp   vs36,   vs8,    vs22        // real*real, imag*real
    xvmaddadp   vs37,   vs8,    vs23        // real*imag, imag*imag
    xvmaddadp   vs38,   vs9,    vs22        // real*real, imag*real
    xvmaddadp   vs39,   vs9,    vs23        // real*imag, imag*imag 
     
.endm

.macro KERNEL2x2 
  LOAD2x2 0
  END2x2  AO, BO, 32,64 
.endm

.macro SAVE2x2

    mr      T1, CO

#ifndef TRMMKERNEL

    lxv     vs16,   0(T1)
    lxv     vs17,   16(T1)

#endif

   AGGREGATE_INTO_COMPLEX vs32,vs33,vs8
   AGGREGATE_INTO_COMPLEX vs34,vs35,vs9 

#ifndef TRMMKERNEL

    xvadddp     vs8,    vs8,    vs16
    xvadddp     vs9,    vs9,    vs17

#endif

    stxv        vs8,    0(T1)
    stxv        vs9,    16(T1)

    add     T1, T1, LDC

#ifndef TRMMKERNEL

    lxv     vs16,   0(T1)
    lxv     vs17,   16(T1)

#endif

   AGGREGATE_INTO_COMPLEX vs36,vs37,vs8
   AGGREGATE_INTO_COMPLEX vs38,vs39,vs9

#ifndef TRMMKERNEL

    xvadddp     vs8,    vs8,    vs16
    xvadddp     vs9,    vs9,    vs17

#endif

    stxv        vs8,    0(T1)
    stxv        vs9,    16(T1)
 
    addi        CO, CO, 32

.endm

/**********************************************************************************************
* Macros for N=2 and M=1
**********************************************************************************************/

.macro Zero2x1
    xxlxor      vs32,   vs32,   vs32
    xxlxor      vs33,   vs33,   vs33
    xxlxor      vs34,   vs34,   vs34
    xxlxor      vs35,   vs35,   vs35 
.endm

.macro LOAD2x1 Zero
    lxv     vs0,    0(AO)       // load real,imag from A 

    lxv     vs16,   0(BO)       // load real part from B
    lxv     vs17,   16(BO)      // load imag part from B
    lxv     vs18,   32(BO)      // load real part from B
    lxv     vs19,   48(BO)      // load imag part from B

.if \Zero==1
    Zero2x1 
.endif

.endm

.macro END2x1_NORMAL
   END2x1 AO,BO,16,64
.endm

.macro END2x1   AREG, BREG, OffsetA, OffsetB

.if \OffsetA != 0 
    addi        \AREG, \AREG, \OffsetA 
.endif  
.if \OffsetB != 0 
    addi        \BREG, \BREG, \OffsetB 
.endif

    xvmaddadp   vs32,   vs0,    vs16        // real*real, imag*real
    xvmaddadp   vs33,   vs0,    vs17        // real*imag, imag*imag 

    xvmaddadp   vs34,   vs0,    vs18        // real*real, imag*real
    xvmaddadp   vs35,   vs0,    vs19        // real*imag, imag*imag 
  
.endm

.macro KERNEL2x1_L      OffsetA,OffsetB, Index,IsLast  
  KERNEL2x1_2   AO,BO,   \OffsetA,\OffsetB, \Index,\IsLast ,0
.endm

.macro KERNEL2x1_E      OffsetA,OffsetB, Index,IsLast  
  KERNEL2x1_2   AO,BO,   \OffsetA,\OffsetB, \Index,\IsLast ,1
.endm

.macro KERNEL2x1_2   AREG,BREG,   OffsetA,OffsetB, Index,IsLast ,Complete

    lxv     vs8,     DISP2(\Index, 0 + \OffsetA)(\AREG)        // load real,imag from A

lxv     vs20,   DISP8(\Index,   0+\OffsetB)(\BREG)      // load real part from B
    lxv     vs21,   DISP8(\Index,16+\OffsetB)(\BREG)        // load imag part from B
    lxv     vs22,   DISP8(\Index,32+\OffsetB)(\BREG)        // load real part from B
    lxv     vs23,   DISP8(\Index,48+\OffsetB)(\BREG)        // load imag part from B

    xvmaddadp   vs32,   vs0,    vs16        // real*real, imag*real
    xvmaddadp   vs33,   vs0,    vs17        // real*imag, imag*imag 

    xvmaddadp   vs34,   vs0,    vs18        // real*real, imag*real
    xvmaddadp   vs35,   vs0,    vs19        // real*imag, imag*imag 

.if \Complete==0
    lxv     vs0,     DISP2(\Index,16 + \OffsetA)(\AREG)      // load real,imag from A 

    lxv     vs16,   DISP8(\Index,   64+\OffsetB)(\BREG)     // load real part from B
    lxv     vs17,   DISP8(\Index,64+16+\OffsetB)(\BREG)     // load imag part from B
    lxv     vs18,   DISP8(\Index,64+32+\OffsetB)(\BREG)     // load real part from B
    lxv     vs19,   DISP8(\Index,64+48+\OffsetB)(\BREG)     // load imag part from B
.endif

.if \IsLast==1  
.if \Complete==1
    addi        \AREG, \AREG,  DISP2(\Index,16+\OffsetA)
    addi        \BREG, \BREG,  DISP8(\Index,64+\OffsetB)
.else 
    addi        \AREG, \AREG,  DISP2(\Index,32)
    addi        \BREG, \BREG,  DISP8(\Index,128)
.endif
.endif  

    xvmaddadp   vs32,   vs8,    vs20        // real*real, imag*real
    xvmaddadp   vs33,   vs8,    vs21        // real*imag, imag*imag 
 
    xvmaddadp   vs34,   vs8,    vs22        // real*real, imag*real
    xvmaddadp   vs35,   vs8,    vs23        // real*imag, imag*imag  
     
.endm

.macro KERNEL2x1 
  LOAD2x1 0
  END2x1  AO, BO, 16,64 
.endm

.macro SAVE2x1

    mr      T1, CO
#ifndef TRMMKERNEL
    lxv     vs16,   0(T1)
#endif
   AGGREGATE_INTO_COMPLEX vs32,vs33,vs8

#ifndef TRMMKERNEL
    xvadddp     vs8,    vs8,    vs16
#endif

    stxv        vs8,    0(T1)

    add     T1, T1, LDC

#ifndef TRMMKERNEL
    lxv     vs16,   0(T1)
#endif

   AGGREGATE_INTO_COMPLEX vs34,vs35,vs8

#ifndef TRMMKERNEL
    xvadddp     vs8,    vs8,    vs16
#endif

    stxv        vs8,    0(T1)

    addi        CO, CO, 16

.endm

/**********************************************************************************************
* Macros for N=1 and M=8
**********************************************************************************************/
.macro Zero1x8
    xxlxor      vs32,   vs32,   vs32
    xxlxor      vs33,   vs33,   vs33
    xxlxor      vs34,   vs34,   vs34
    xxlxor      vs35,   vs35,   vs35
    xxlxor      vs36,   vs36,   vs36
    xxlxor      vs37,   vs37,   vs37
    xxlxor      vs38,   vs38,   vs38
    xxlxor      vs39,   vs39,   vs39
    xxlxor      vs40,   vs40,   vs40
    xxlxor      vs41,   vs41,   vs41
    xxlxor      vs42,   vs42,   vs42
    xxlxor      vs43,   vs43,   vs43
    xxlxor      vs44,   vs44,   vs44
    xxlxor      vs45,   vs45,   vs45
    xxlxor      vs46,   vs46,   vs46
    xxlxor      vs47,   vs47,   vs47     
.endm

.macro LOAD1x8 Zero

    lxv     vs16,   0(BO)       // load real part from B
    lxv     vs17,   16(BO)      // load imag part from B 

    lxv     vs0,    0(AO)       // load real,imag from A
    lxv     vs1,    16(AO)      // load real,imag from A
    lxv     vs2,    32(AO)      // load real,imag from A
    lxv     vs3,    48(AO)      // load real,imag from A

    lxv     vs4,    64(AO)      // load real,imag from A
    lxv     vs5,    80(AO)      // load real,imag from A
    lxv     vs6,    96(AO)      // load real,imag from A
    lxv     vs7,    112(AO)     // load real,imag from A

.if \Zero==1
    Zero1x8 
.endif

.endm

.macro END1x8_NORMAL
   END1x8 AO,BO,128,32
.endm

.macro END1x8   AREG, BREG, OffsetA, OffsetB

.if \OffsetB != 0 
    addi        \BREG, \BREG, \OffsetB 
.endif
.if \OffsetA != 0 
    addi        \AREG, \AREG, \OffsetA 
.endif  

    xvmaddadp   vs32,   vs0,    vs16        // real*real, imag*real
    xvmaddadp   vs33,   vs0,    vs17        // real*imag, imag*imag
    xvmaddadp   vs34,   vs1,    vs16        // real*real, imag*real
    xvmaddadp   vs35,   vs1,    vs17        // real*imag, imag*imag
    xvmaddadp   vs36,   vs2,    vs16        // real*real, imag*real
    xvmaddadp   vs37,   vs2,    vs17        // real*imag, imag*imag
    xvmaddadp   vs38,   vs3,    vs16        // real*real, imag*real
    xvmaddadp   vs39,   vs3,    vs17        // real*imag, imag*imag
    xvmaddadp   vs40,   vs4,    vs16        // real*real, imag*real
    xvmaddadp   vs41,   vs4,    vs17        // real*imag, imag*imag
    xvmaddadp   vs42,   vs5,    vs16        // real*real, imag*real
    xvmaddadp   vs43,   vs5,    vs17        // real*imag, imag*imag
    xvmaddadp   vs44,   vs6,    vs16        // real*real, imag*real
    xvmaddadp   vs45,   vs6,    vs17        // real*imag, imag*imag
    xvmaddadp   vs46,   vs7,    vs16        // real*real, imag*real
    xvmaddadp   vs47,   vs7,    vs17        // real*imag, imag*imag

.endm

.macro KERNEL1x8_L      OffsetA,OffsetB, Index,IsLast  
  KERNEL1x8_2   AO,BO,   \OffsetA,\OffsetB, \Index,\IsLast ,0
.endm

.macro KERNEL1x8_E      OffsetA,OffsetB, Index,IsLast  
  KERNEL1x8_2   AO,BO,   \OffsetA,\OffsetB, \Index,\IsLast ,1
.endm

.macro KERNEL1x8_2   AREG,BREG,   OffsetA,OffsetB, Index,IsLast ,Complete

    lxv     vs8,     DISP16(\Index, 0 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs9,     DISP16(\Index,16 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs10,    DISP16(\Index,32 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs11,    DISP16(\Index,48 + \OffsetA)(\AREG)        // load real,imag from A 

    lxv     vs12,    DISP16(\Index, 64 + \OffsetA)(\AREG)       // load real,imag from A
    lxv     vs13,    DISP16(\Index,64+16 + \OffsetA)(\AREG)     // load real,imag from A
    lxv     vs14,    DISP16(\Index,64+32 + \OffsetA)(\AREG)     // load real,imag from A
    lxv     vs15,    DISP16(\Index,64+48 + \OffsetA)(\AREG)     // load real,imag from A

    lxv     vs20,   DISP4(\Index,   0+\OffsetB)(\BREG)      // load real part from B
    lxv     vs21,   DISP4(\Index,16+\OffsetB)(\BREG)        // load imag part from B 

    xvmaddadp   vs32,   vs0,    vs16        // real*real, imag*real
    xvmaddadp   vs33,   vs0,    vs17        // real*imag, imag*imag
    xvmaddadp   vs34,   vs1,    vs16        // real*real, imag*real
    xvmaddadp   vs35,   vs1,    vs17        // real*imag, imag*imag
    xvmaddadp   vs36,   vs2,    vs16        // real*real, imag*real
    xvmaddadp   vs37,   vs2,    vs17        // real*imag, imag*imag
    xvmaddadp   vs38,   vs3,    vs16        // real*real, imag*real
    xvmaddadp   vs39,   vs3,    vs17        // real*imag, imag*imag
    xvmaddadp   vs40,   vs4,    vs16        // real*real, imag*real
    xvmaddadp   vs41,   vs4,    vs17        // real*imag, imag*imag
    xvmaddadp   vs42,   vs5,    vs16        // real*real, imag*real
    xvmaddadp   vs43,   vs5,    vs17        // real*imag, imag*imag
    xvmaddadp   vs44,   vs6,    vs16        // real*real, imag*real
    xvmaddadp   vs45,   vs6,    vs17        // real*imag, imag*imag
    xvmaddadp   vs46,   vs7,    vs16        // real*real, imag*real
    xvmaddadp   vs47,   vs7,    vs17        // real*imag, imag*imag

.if \Complete==0
    lxv     vs0,     DISP16(\Index,128+ + \OffsetA)(\AREG)      // load real,imag from A
    lxv     vs1,     DISP16(\Index,128+16 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs2,     DISP16(\Index,128+32 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs3,     DISP16(\Index,128+48 + \OffsetA)(\AREG)        // load real,imag from A

    lxv     vs4,     DISP16(\Index, 192 + \OffsetA)(\AREG)      // load real,imag from A
    lxv     vs5,     DISP16(\Index,192 +16 + \OffsetA)(\AREG)       // load real,imag from A
    lxv     vs6,     DISP16(\Index,192 +32 + \OffsetA)(\AREG)       // load real,imag from A
    lxv     vs7,     DISP16(\Index,192 +48 + \OffsetA)(\AREG)       // load real,imag from A

    lxv     vs16,   DISP4(\Index,   32+\OffsetB)(\BREG)     // load real part from B
    lxv     vs17,   DISP4(\Index,48+\OffsetB)(\BREG)     // load imag part from B 
.endif

.if \IsLast==1  
.if \Complete==1
    addi        \AREG, \AREG, DISP16(\Index,128+\OffsetA)
    addi        \BREG, \BREG, DISP4(\Index,32+\OffsetB)
.else 
    addi        \AREG, \AREG, DISP16(\Index,256)
    addi        \BREG, \BREG, DISP4(\Index,64)
.endif
.endif  

    xvmaddadp   vs32,   vs8,    vs20        // real*real, imag*real
    xvmaddadp   vs33,   vs8,    vs21        // real*imag, imag*imag
    xvmaddadp   vs34,   vs9,    vs20        // real*real, imag*real
    xvmaddadp   vs35,   vs9,    vs21        // real*imag, imag*imag
    xvmaddadp   vs36,   vs10,   vs20        // real*real, imag*real
    xvmaddadp   vs37,   vs10,   vs21        // real*imag, imag*imag
    xvmaddadp   vs38,   vs11,   vs20        // real*real, imag*real
    xvmaddadp   vs39,   vs11,   vs21        // real*imag, imag*imag
    xvmaddadp   vs40,   vs12,   vs20        // real*real, imag*real
    xvmaddadp   vs41,   vs12,   vs21        // real*imag, imag*imag
    xvmaddadp   vs42,   vs13,   vs20        // real*real, imag*real
    xvmaddadp   vs43,   vs13,   vs21        // real*imag, imag*imag
    xvmaddadp   vs44,   vs14,   vs20        // real*real, imag*real
    xvmaddadp   vs45,   vs14,   vs21        // real*imag, imag*imag
    xvmaddadp   vs46,   vs15,   vs20        // real*real, imag*real
    xvmaddadp   vs47,   vs15,   vs21        // real*imag, imag*imag

.endm

.macro KERNEL1x8 
  LOAD1x8 0
  END1x8  AO, BO, 128,32 
.endm

.macro SAVE1x8

     mr      T1, CO
    addi        T2, T1, 64

#ifndef TRMMKERNEL

    lxv     vs16,   0(T1)
    lxv     vs17,   16(T1)
    lxv     vs18,   32(T1)
    lxv     vs19,   48(T1)
    lxv     vs20,   0(T2)
    lxv     vs21,   16(T2)
    lxv     vs22,   32(T2)
    lxv     vs23,   48(T2)

#endif

   AGGREGATE_INTO_COMPLEX vs32,vs33,vs8
   AGGREGATE_INTO_COMPLEX vs34,vs35,vs9
   AGGREGATE_INTO_COMPLEX vs36,vs37,vs10
   AGGREGATE_INTO_COMPLEX vs38,vs39,vs11
   AGGREGATE_INTO_COMPLEX vs40,vs41,vs12
   AGGREGATE_INTO_COMPLEX vs42,vs43,vs13
   AGGREGATE_INTO_COMPLEX vs44,vs45,vs14
   AGGREGATE_INTO_COMPLEX vs46,vs47,vs15

#ifndef TRMMKERNEL

    xvadddp     vs8,    vs8,    vs16
    xvadddp     vs9,    vs9,    vs17
    xvadddp     vs10,   vs10,   vs18
    xvadddp     vs11,   vs11,   vs19
    xvadddp     vs12,   vs12,   vs20
    xvadddp     vs13,   vs13,   vs21
    xvadddp     vs14,   vs14,   vs22
    xvadddp     vs15,   vs15,   vs23

#endif

    stxv        vs8,    0(T1)
    stxv        vs9,    16(T1)
    stxv        vs10,   32(T1)
    stxv        vs11,   48(T1)
    stxv        vs12,   0(T2)
    stxv        vs13,   16(T2)
    stxv        vs14,   32(T2)
    stxv        vs15,   48(T2)

    addi        CO, CO, 128

.endm

/**********************************************************************************************
* Macros for N=1 and M=4
**********************************************************************************************/

.macro Zero1x4
    xxlxor      vs32,   vs32,   vs32
    xxlxor      vs33,   vs33,   vs33
    xxlxor      vs34,   vs34,   vs34
    xxlxor      vs35,   vs35,   vs35
    xxlxor      vs36,   vs36,   vs36
    xxlxor      vs37,   vs37,   vs37
    xxlxor      vs38,   vs38,   vs38
    xxlxor      vs39,   vs39,   vs39 
.endm

.macro LOAD1x4 Zero

    lxv     vs16,   0(BO)       // load real part from B
    lxv     vs17,   16(BO)      // load imag part from B 

    lxv     vs0,    0(AO)       // load real,imag from A
    lxv     vs1,    16(AO)      // load real,imag from A
    lxv     vs2,    32(AO)      // load real,imag from A
    lxv     vs3,    48(AO)      // load real,imag from A
 
.if \Zero==1
    Zero1x4 
.endif

.endm

.macro END1x4_NORMAL
   END1x4 AO,BO,64,32
.endm

.macro END1x4   AREG, BREG, OffsetA, OffsetB

.if \OffsetB != 0 
    addi        \BREG, \BREG, \OffsetB 
.endif
.if \OffsetA != 0 
    addi        \AREG, \AREG, \OffsetA 
.endif  

    xvmaddadp   vs32,   vs0,    vs16        // real*real, imag*real
    xvmaddadp   vs33,   vs0,    vs17        // real*imag, imag*imag
    xvmaddadp   vs34,   vs1,    vs16        // real*real, imag*real
    xvmaddadp   vs35,   vs1,    vs17        // real*imag, imag*imag
    xvmaddadp   vs36,   vs2,    vs16        // real*real, imag*real
    xvmaddadp   vs37,   vs2,    vs17        // real*imag, imag*imag
    xvmaddadp   vs38,   vs3,    vs16        // real*real, imag*real
    xvmaddadp   vs39,   vs3,    vs17        // real*imag, imag*imag

.endm

.macro KERNEL1x4_L      OffsetA,OffsetB, Index,IsLast  
  KERNEL1x4_2   AO,BO,   \OffsetA,\OffsetB, \Index,\IsLast ,0
.endm

.macro KERNEL1x4_E      OffsetA,OffsetB, Index,IsLast  
  KERNEL1x4_2   AO,BO,   \OffsetA,\OffsetB, \Index,\IsLast ,1
.endm

.macro KERNEL1x4_2   AREG,BREG,   OffsetA,OffsetB, Index,IsLast ,Complete

    lxv     vs8,     DISP8(\Index, 0 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs9,     DISP8(\Index,16 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs10,    DISP8(\Index,32 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs11,    DISP8(\Index,48 + \OffsetA)(\AREG)        // load real,imag from A

lxv     vs20,   DISP4(\Index,   0+\OffsetB)(\BREG)      // load real part from B
    lxv     vs21,   DISP4(\Index,16+\OffsetB)(\BREG)        // load imag part from B 

    xvmaddadp   vs32,   vs0,    vs16        // real*real, imag*real
    xvmaddadp   vs33,   vs0,    vs17        // real*imag, imag*imag
    xvmaddadp   vs34,   vs1,    vs16        // real*real, imag*real
    xvmaddadp   vs35,   vs1,    vs17        // real*imag, imag*imag
    xvmaddadp   vs36,   vs2,    vs16        // real*real, imag*real
    xvmaddadp   vs37,   vs2,    vs17        // real*imag, imag*imag
    xvmaddadp   vs38,   vs3,    vs16        // real*real, imag*real
    xvmaddadp   vs39,   vs3,    vs17        // real*imag, imag*imag

    xvmaddadp   vs40,   vs0,    vs18        // real*real, imag*real
    xvmaddadp   vs41,   vs0,    vs19        // real*imag, imag*imag
    xvmaddadp   vs42,   vs1,    vs18        // real*real, imag*real
    xvmaddadp   vs43,   vs1,    vs19        // real*imag, imag*imag
    xvmaddadp   vs44,   vs2,    vs18        // real*real, imag*real
    xvmaddadp   vs45,   vs2,    vs19        // real*imag, imag*imag
    xvmaddadp   vs46,   vs3,    vs18        // real*real, imag*real
    xvmaddadp   vs47,   vs3,    vs19        // real*imag, imag*imag

.if \Complete==0
    lxv     vs0,     DISP8(\Index,64+  \OffsetA)(\AREG)      // load real,imag from A
    lxv     vs1,     DISP8(\Index,64+16 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs2,     DISP8(\Index,64+32 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs3,     DISP8(\Index,64+48 + \OffsetA)(\AREG)        // load real,imag from A 

    lxv     vs16,   DISP4(\Index,   32+\OffsetB)(\BREG)     // load real part from B
    lxv     vs17,   DISP4(\Index,32+16+\OffsetB)(\BREG)     // load imag part from B 
.endif

.if \IsLast==1  
.if \Complete==1
    addi        \AREG, \AREG,  DISP8(\Index,64+\OffsetA)
    addi        \BREG, \BREG,  DISP4(\Index,32+\OffsetB)
.else 
    addi        \AREG, \AREG,  DISP8(\Index,128)
    addi        \BREG, \BREG,  DISP4(\Index,64)
.endif
.endif  

    xvmaddadp   vs32,   vs8,    vs20        // real*real, imag*real
    xvmaddadp   vs33,   vs8,    vs21        // real*imag, imag*imag
    xvmaddadp   vs34,   vs9,    vs20        // real*real, imag*real
    xvmaddadp   vs35,   vs9,    vs21        // real*imag, imag*imag
    xvmaddadp   vs36,   vs10,   vs20        // real*real, imag*real
    xvmaddadp   vs37,   vs10,   vs21        // real*imag, imag*imag
    xvmaddadp   vs38,   vs11,   vs20        // real*real, imag*real
    xvmaddadp   vs39,   vs11,   vs21        // real*imag, imag*imag
 
    xvmaddadp   vs40,   vs8,    vs22        // real*real, imag*real
    xvmaddadp   vs41,   vs8,    vs23        // real*imag, imag*imag
    xvmaddadp   vs42,   vs9,    vs22        // real*real, imag*real
    xvmaddadp   vs43,   vs9,    vs23        // real*imag, imag*imag
    xvmaddadp   vs44,   vs10,   vs22        // real*real, imag*real
    xvmaddadp   vs45,   vs10,   vs23        // real*imag, imag*imag
    xvmaddadp   vs46,   vs11,   vs22        // real*real, imag*real
    xvmaddadp   vs47,   vs11,   vs23        // real*imag, imag*imag

.endm

.macro KERNEL1x4 
  LOAD1x4 0
  END1x4  AO, BO, 64,32 
.endm

.macro SAVE1x4

    mr      T1, CO

#ifndef TRMMKERNEL

    lxv     vs16,   0(T1)
    lxv     vs17,   16(T1)
    lxv     vs18,   32(T1)
    lxv     vs19,   48(T1)

#endif

   AGGREGATE_INTO_COMPLEX vs32,vs33,vs8
   AGGREGATE_INTO_COMPLEX vs34,vs35,vs9
   AGGREGATE_INTO_COMPLEX vs36,vs37,vs10
   AGGREGATE_INTO_COMPLEX vs38,vs39,vs11

#ifndef TRMMKERNEL

    xvadddp     vs8,    vs8,    vs16
    xvadddp     vs9,    vs9,    vs17
    xvadddp     vs10,   vs10,   vs18
    xvadddp     vs11,   vs11,   vs19

#endif

    stxv        vs8,    0(T1)
    stxv        vs9,    16(T1)
    stxv        vs10,   32(T1)
    stxv        vs11,   48(T1) 
 
    addi        CO, CO, 64

.endm

/**********************************************************************************************
* Macros for N=1 and M=2
**********************************************************************************************/

.macro Zero1x2
    xxlxor      vs32,   vs32,   vs32
    xxlxor      vs33,   vs33,   vs33
    xxlxor      vs34,   vs34,   vs34
    xxlxor      vs35,   vs35,   vs35 
.endm

.macro LOAD1x2 Zero

    lxv     vs16,   0(BO)       // load real part from B
    lxv     vs17,   16(BO)      // load imag part from B 

    lxv     vs0,    0(AO)       // load real,imag from A
    lxv     vs1,    16(AO)      // load real,imag from A 
 
.if \Zero==1
    Zero1x2 
.endif

.endm

.macro END1x2_NORMAL
   END1x2 AO,BO,32,32
.endm

.macro END1x2   AREG, BREG, OffsetA, OffsetB

.if \OffsetB != 0 
    addi        \BREG, \BREG, \OffsetB 
.endif
.if \OffsetA != 0 
    addi        \AREG, \AREG, \OffsetA 
.endif  

    xvmaddadp   vs32,   vs0,    vs16        // real*real, imag*real
    xvmaddadp   vs33,   vs0,    vs17        // real*imag, imag*imag
    xvmaddadp   vs34,   vs1,    vs16        // real*real, imag*real
    xvmaddadp   vs35,   vs1,    vs17        // real*imag, imag*imag
  
.endm

.macro KERNEL1x2_L      OffsetA,OffsetB, Index,IsLast  
  KERNEL1x2_2   AO,BO,   \OffsetA,\OffsetB, \Index,\IsLast ,0
.endm

.macro KERNEL1x2_E      OffsetA,OffsetB, Index,IsLast  
  KERNEL1x2_2   AO,BO,   \OffsetA,\OffsetB, \Index,\IsLast ,1
.endm

.macro KERNEL1x2_2   AREG,BREG,   OffsetA,OffsetB, Index,IsLast ,Complete

    lxv     vs8,     DISP4(\Index, 0 + \OffsetA)(\AREG)        // load real,imag from A
    lxv     vs9,     DISP4(\Index,16 + \OffsetA)(\AREG)        // load real,imag from A

lxv     vs20,   DISP4(\Index,   0+\OffsetB)(\BREG)      // load real part from B
    lxv     vs21,   DISP4(\Index,16+\OffsetB)(\BREG)        // load imag part from B 

    xvmaddadp   vs32,   vs0,    vs16        // real*real, imag*real
    xvmaddadp   vs33,   vs0,    vs17        // real*imag, imag*imag
    xvmaddadp   vs34,   vs1,    vs16        // real*real, imag*real
    xvmaddadp   vs35,   vs1,    vs17        // real*imag, imag*imag  
.if \Complete==0
    lxv     vs0,     DISP4(\Index,32 + \OffsetA)(\AREG)      // load real,imag from A
    lxv     vs1,     DISP4(\Index,48+ \OffsetA)(\AREG)        // load real,imag from A 

    lxv     vs16,   DISP4(\Index,   32+\OffsetB)(\BREG)     // load real part from B
    lxv     vs17,   DISP4(\Index,32+16+\OffsetB)(\BREG)     // load imag part from B 
.endif

.if \IsLast==1  
.if \Complete==1
    addi        \AREG, \AREG,  DISP4(\Index,32+\OffsetA)
    addi        \BREG, \BREG,  DISP4(\Index,32+\OffsetB)
.else 
    addi        \AREG, \AREG,  DISP4(\Index,64)
    addi        \BREG, \BREG,  DISP4(\Index,64)
.endif
.endif  

    xvmaddadp   vs32,   vs8,    vs20        // real*real, imag*real
    xvmaddadp   vs33,   vs8,    vs21        // real*imag, imag*imag
    xvmaddadp   vs34,   vs9,    vs20        // real*real, imag*real
    xvmaddadp   vs35,   vs9,    vs21        // real*imag, imag*imag

.endm

.macro KERNEL1x2 
  LOAD1x2 0
  END1x2  AO, BO, 32,32 
.endm

.macro SAVE1x2

    mr      T1, CO

#ifndef TRMMKERNEL

    lxv     vs16,   0(T1)
    lxv     vs17,   16(T1)

#endif

   AGGREGATE_INTO_COMPLEX vs32,vs33,vs8
   AGGREGATE_INTO_COMPLEX vs34,vs35,vs9 

#ifndef TRMMKERNEL

    xvadddp     vs8,    vs8,    vs16
    xvadddp     vs9,    vs9,    vs17

#endif

    stxv        vs8,    0(T1)
    stxv        vs9,    16(T1)

addi        CO, CO, 32

.endm

/**********************************************************************************************
* Macros for N=1 and M=1
**********************************************************************************************/

.macro Zero1x1
    xxlxor      vs32,   vs32,   vs32
    xxlxor      vs33,   vs33,   vs33 
.endm

.macro LOAD1x1 Zero
    lxv     vs0,    0(AO)       // load real,imag from A 

    lxv     vs16,   0(BO)       // load real part from B
    lxv     vs17,   16(BO)      // load imag part from B 

.if \Zero==1
    Zero1x1 
.endif

.endm

.macro END1x1_NORMAL
   END1x1 AO,BO,16,32
.endm

.macro END1x1   AREG, BREG, OffsetA, OffsetB

.if \OffsetA != 0 
    addi        \AREG, \AREG, \OffsetA 
.endif  
.if \OffsetB != 0 
    addi        \BREG, \BREG, \OffsetB 
.endif

    xvmaddadp   vs32,   vs0,    vs16        // real*real, imag*real
    xvmaddadp   vs33,   vs0,    vs17        // real*imag, imag*imag 
 
  
.endm

.macro KERNEL1x1_L      OffsetA,OffsetB, Index,IsLast  
  KERNEL1x1_2   AO,BO,   \OffsetA,\OffsetB, \Index,\IsLast ,0
.endm

.macro KERNEL1x1_E      OffsetA,OffsetB, Index,IsLast  
  KERNEL1x1_2   AO,BO,   \OffsetA,\OffsetB, \Index,\IsLast ,1
.endm

.macro KERNEL1x1_2   AREG,BREG,   OffsetA,OffsetB, Index,IsLast ,Complete

    lxv     vs8,     DISP2(\Index, 0 + \OffsetA)(\AREG)        // load real,imag from A

    lxv     vs20,   DISP4(\Index,   0+\OffsetB)(\BREG)      // load real part from B
    lxv     vs21,   DISP4(\Index,16+\OffsetB)(\BREG)        // load imag part from B 

    xvmaddadp   vs32,   vs0,    vs16        // real*real, imag*real
    xvmaddadp   vs33,   vs0,    vs17        // real*imag, imag*imag  

.if \Complete==0
    lxv     vs0,     DISP2(\Index,16 + \OffsetA)(\AREG)      // load real,imag from A 

    lxv     vs16,   DISP4(\Index,   32+\OffsetB)(\BREG)     // load real part from B
    lxv     vs17,   DISP4(\Index,32+16+\OffsetB)(\BREG)     // load imag part from B 
.endif


.if \IsLast==1  
.if \Complete==1
    addi        \AREG, \AREG,  DISP2(\Index,16+\OffsetA)
    addi        \BREG, \BREG,  DISP4(\Index,32+\OffsetB)
.else 
    addi        \AREG, \AREG,  DISP2(\Index,32)
    addi        \BREG, \BREG,  DISP4(\Index,64)
.endif
.endif
  
    xvmaddadp   vs32,   vs8,    vs20        // real*real, imag*real
    xvmaddadp   vs33,   vs8,    vs21        // real*imag, imag*imag 
  
     
.endm

.macro KERNEL1x1 
  LOAD1x1 0
  END1x1  AO, BO, 16,32 

.endm  

.macro SAVE1x1

    mr      T1, CO
#ifndef TRMMKERNEL
    lxv     vs16,   0(T1)
#endif
   AGGREGATE_INTO_COMPLEX vs32,vs33,vs8

#ifndef TRMMKERNEL
    xvadddp     vs8,    vs8,    vs16
#endif

    stxv        vs8,    0(T1)

addi        CO, CO, 16

.endm


.macro ZCOPYB_2

        lxv          vs32,   0(BO)
        lxv          vs33,  16(BO)            
        addi            BO,     BO,     32
        xxspltd     vs40, vs32, 1
        xxspltd     vs41, vs32, 0     
        xxspltd     vs42, vs33, 1
        xxspltd     vs43, vs33, 0

        stxv         vs40,    0(BBO)
        stxv         vs41,   16(BBO)
        stxv         vs42,   32(BBO)
        stxv         vs43,   48(BBO)
        addi            BBO,    BBO,    64

.endm

.macro ZCOPYB_1

        lxv          vs32,   0(BO)              
        addi            BO,     BO,     16
        xxspltd     vs40, vs32, 1
        xxspltd     vs41, vs32, 0        
        stxv         vs40,    0(BBO)
        stxv         vs41,   16(BBO)

        addi            BBO,    BBO,    32

.endm

.macro ZCOPYB_8

        lxv          vs32,   0(BO)
        lxv          vs33,  16(BO)
        lxv          vs34,  32(BO)
        lxv          vs35,  48(BO) 

        lxv          vs36,   64+0(BO)
        lxv          vs37,  64+16(BO)
        lxv          vs38,  64+32(BO)
        lxv          vs39,  64+48(BO) 
        addi         BO, BO,    128
        xxspltd     vs40, vs32, 1
        xxspltd     vs41, vs32, 0
        xxspltd     vs42, vs33, 1
        xxspltd     vs43, vs33, 0
        xxspltd     vs44, vs34, 1
        xxspltd     vs45, vs34, 0
        xxspltd     vs46, vs35, 1
        xxspltd     vs47, vs35, 0    

        xxspltd     vs48, vs36, 1
        xxspltd     vs49, vs36, 0
        xxspltd     vs50, vs37, 1
        xxspltd     vs51, vs37, 0
        xxspltd     vs52, vs38, 1
        xxspltd     vs53, vs38, 0
        xxspltd     vs54, vs39, 1
        xxspltd     vs55, vs39, 0

        stxv         vs40,    0(BBO)
        stxv         vs41,   16(BBO)
        stxv         vs42,   32(BBO)
        stxv         vs43,   48(BBO) 

        stxv         vs44,    64+0(BBO)
        stxv         vs45,   64+16(BBO)
        stxv         vs46,   64+32(BBO)
        stxv         vs47,   64+48(BBO) 

        stxv         vs48,   128+ 0(BBO)
        stxv         vs49,   128+ 16(BBO)
        stxv         vs50,   128+ 32(BBO)
        stxv         vs51,   128+ 48(BBO) 

        stxv         vs52,   192 + 0(BBO)
        stxv         vs53,   192 + 16(BBO)
        stxv         vs54,   192+ 32(BBO)
        stxv         vs55,   192 + 48(BBO)
        addi            BBO,    BBO,    256

.endm

