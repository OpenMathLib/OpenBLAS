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


#define unit_size 16
#define DISP32(ind,disp) (ind*unit_size*32+disp)
#define DISP16(ind,disp) (ind*unit_size*16+disp)
#define DISP8(ind,disp) (ind*unit_size*8+disp)
#define DISP4(ind,disp) (ind*unit_size*4+disp)
#define DISP2(ind,disp) (ind*unit_size*2+disp)
#define DISP1(ind,disp) (ind*unit_size+disp)
#define DISPX(disp)  (disp)

/*	HELPERS FOR SAVE	*/

/* {r0,i0} and {r1,i1} into  {r0,r1} {i0,i1} */
.macro LOAD_COUPLE_AS_RR_II  VS_OUT1,VS_OUT2,VS_TEMP1,VS_TEMP2,REG,LOFFSET 
#ifndef TRMMKERNEL 
  lxv	\VS_TEMP1,	DISPX(\LOFFSET)(\REG)
  lxv	\VS_TEMP2,	DISPX(\LOFFSET+16)(\REG)
  xxmrgld  \VS_OUT1,\VS_TEMP1,\VS_TEMP2
  xxmrghd  \VS_OUT2,\VS_TEMP1,\VS_TEMP2	
#endif	
.endm

/*from 2 result {a0r*br,a0i*bi} and {a1r*br,a1i*bi} pack into {a0r*br,a1r*br} and {a0i*bi,a1i*bi}*/
.macro RESULT_INTO_REALREAL_IMAGEIMAGE VSIN1,VSIN2,VSOUT1,VSOUT2
	xxmrgld	\VSOUT1, \VSIN1,\VSIN2 /*  real*real from 2 results*/
	xxmrghd	\VSOUT2, \VSIN1,\VSIN2 /*  imag*imag from 2 results*/
.endm 

/*from 2 result {a0r*bi,a0i*br} and {a1r*bi,a1i*br} pack into {a0r*bi,a1r*bi} and {a0i*br,a1i*br}*/
.macro RESULT_INTO_REALIMAG_IMAGREAL VSIN1,VSIN2,VSOUT1,VSOUT2 
	xxmrgld	\VSOUT1, \VSIN1,\VSIN2 /*  real*imag */
	xxmrghd	\VSOUT2, \VSIN1,\VSIN2 /*  imag*real*/
.endm

/* {a0r*br op a0i*bi ,a1r*br op a1i*bi} ~ {r0,r1}; {a0r*bi op a0i*br ,a1r*bi op a1i*br} ~ {i0,i1}*/
.macro  AGGREGATE_REALS_IMAGES  VSINR_OUT1,VSINR,VSINI_OUT2,VSINI
#if	defined(NN) || defined(NT) || defined(TN) || defined(TT) 
	xvsubdp  \VSINR_OUT1,\VSINR_OUT1,\VSINR
	xvadddp  \VSINI_OUT2,\VSINI_OUT2,\VSINI  
#elif  defined(CN) || defined(CT) || defined(RN) || defined(RT) 
	xvadddp  \VSINR_OUT1,\VSINR_OUT1,\VSINR
	xvsubdp  \VSINI_OUT2,\VSINI_OUT2,\VSINI 
#elif  defined(NC) || defined(TC) || defined(NR) || defined(TR) 
	xvadddp  \VSINR_OUT1,\VSINR_OUT1,\VSINR
	xvsubdp  \VSINI_OUT2,\VSINI,\VSINI_OUT2  
#else	// CC || CR || RC || RR 
    /*we will assume {-alpha_r,-alpha_i} for this case */
    /*i1i2-r1r2 so we will negate alpha real instead to fix sign*/
	xvsubdp  \VSINR_OUT1,\VSINR,\VSINR_OUT1
    /*we will negate alpha image instead  instead to fix sign*/
	xvadddp  \VSINI_OUT2,\VSINI_OUT2,\VSINI 
#endif
.endm 

/* {i0,i1} * {alpha_i,alpha_i} - VSOUT1 ;VSOUT2 + {r0,r1}*{alpha_i,alpha_i} */
.macro MULT_APLHA_PART1  VSINRR,VSINII,VSOUT1,VSOUT2
#ifndef TRMMKERNEL  
	xvmsubadp \VSOUT1,\VSINII, alpha_i
	xvmaddadp  \VSOUT2,\VSINRR, alpha_i
#else 
	xvmuldp \VSOUT1,\VSINII, alpha_i 
	xvmuldp  \VSOUT2,\VSINRR, alpha_i
#endif 
.endm

/*   {r0,r1} * {alpha_r,alpha_r} -  VSOUT1 ;VSOUT2 + {i0,i1} * {alpha_r,alpha_r} */
.macro MULT_APLHA_PART2  VSINRR,VSINII,VSOUT1,VSOUT2 
	xvmsubadp  \VSOUT1,\VSINRR, alpha_r
	xvmaddadp \VSOUT2,\VSINII, alpha_r
.endm

/* unpack to store 2{r,r} {i,i} into  {r,i} {r,i} (big endian because of stxv) */
.macro UNPACK_FOR_STORE VSIN1,VSIN2,VSOUT1,VSOUT2 
	xxmrghd  \VSOUT1,\VSIN2,\VSIN1
	xxmrgld  \VSOUT2,\VSIN2,\VSIN1
.endm
.macro STORE_COUPLE REG,LOFFSET,VSIN1,VSIN2
	stxv	\VSIN1,	DISPX(\LOFFSET)(\REG)
	stxv	\VSIN2,	DISPX(\LOFFSET+16)(\REG)
.endm

.macro SAVE8 VSRes1,VSRes2,VSRes3,VSRes4,VSRes5,VSRes6,VSRes7,VSRes8,VSRes9,VSRes10,VSRes11,VSRes12,VSRes13,VSRes14,VSRes15,VSRes16,BASE_REG,LOFFSET
  RESULT_INTO_REALREAL_IMAGEIMAGE \VSRes1,\VSRes3,vs2,vs3
  LOAD_COUPLE_AS_RR_II	vs14,vs15,vs18,vs19,\BASE_REG,\LOFFSET
  RESULT_INTO_REALIMAG_IMAGREAL	\VSRes2,\VSRes4,vs4,vs5
  LOAD_COUPLE_AS_RR_II	vs16,vs17,vs20,vs21,\BASE_REG,(\LOFFSET+32)
  RESULT_INTO_REALREAL_IMAGEIMAGE \VSRes5,\VSRes7,vs6,vs7
  LOAD_COUPLE_AS_RR_II	vs24,vs25,vs18,vs19,\BASE_REG,(\LOFFSET +64)
  RESULT_INTO_REALIMAG_IMAGREAL	\VSRes6,\VSRes8,vs8,vs9 
  LOAD_COUPLE_AS_RR_II	vs26,vs27,vs20,vs21,\BASE_REG,(\LOFFSET+96)
  RESULT_INTO_REALREAL_IMAGEIMAGE \VSRes9,\VSRes11,vs10,vs11
  AGGREGATE_REALS_IMAGES	vs2,vs3,vs4,vs5
  RESULT_INTO_REALIMAG_IMAGREAL	\VSRes10,\VSRes12,vs12,vs13 
  AGGREGATE_REALS_IMAGES	vs6,vs7,vs8,vs9  
  RESULT_INTO_REALREAL_IMAGEIMAGE \VSRes13,\VSRes15,\VSRes1,\VSRes2
  MULT_APLHA_PART1	vs2,vs4, vs14,vs15
  RESULT_INTO_REALIMAG_IMAGREAL	\VSRes14,\VSRes16,\VSRes3,\VSRes4
  MULT_APLHA_PART1	vs6,vs8,vs16,vs17
  MULT_APLHA_PART2  vs2,vs4,vs14,vs15 
  AGGREGATE_REALS_IMAGES	vs10,vs11,vs12,vs13
  MULT_APLHA_PART2	vs6,vs8,vs16,vs17
  AGGREGATE_REALS_IMAGES	\VSRes1,\VSRes2,\VSRes3,\VSRes4	
  UNPACK_FOR_STORE	vs14,vs15,vs7,vs9
  MULT_APLHA_PART1	vs10,vs12, vs24,vs25
  UNPACK_FOR_STORE	vs16,vs17,vs3,vs5 
  MULT_APLHA_PART1	\VSRes1,\VSRes3, vs26,vs27
  STORE_COUPLE	\BASE_REG,\LOFFSET,vs7,vs9
  MULT_APLHA_PART2	vs10,vs12,vs24,vs25
  STORE_COUPLE	\BASE_REG,(\LOFFSET+32),vs3,vs5 
  MULT_APLHA_PART2	\VSRes1,\VSRes3, vs26,vs27
  UNPACK_FOR_STORE	vs24,vs25,vs10,vs12
  UNPACK_FOR_STORE	vs26,vs27,\VSRes1,\VSRes3
  STORE_COUPLE	\BASE_REG,(\LOFFSET +64),vs10,vs12
  STORE_COUPLE	\BASE_REG,(\LOFFSET+96),\VSRes1,\VSRes3
.endm

.macro SAVE4  VSRes1,VSRes2,VSRes3,VSRes4,VSRes5,VSRes6,VSRes7,VSRes8,BASE_REG,LOFFSET
  RESULT_INTO_REALREAL_IMAGEIMAGE \VSRes1,\VSRes3,vs2,vs3
  LOAD_COUPLE_AS_RR_II	vs14,vs15,vs18,vs19,\BASE_REG,\LOFFSET
  RESULT_INTO_REALIMAG_IMAGREAL	\VSRes2,\VSRes4,vs4,vs5
  LOAD_COUPLE_AS_RR_II	vs16,vs17,vs20,vs21,\BASE_REG,(\LOFFSET+32)
  RESULT_INTO_REALREAL_IMAGEIMAGE \VSRes5,\VSRes7,vs6,vs7
  RESULT_INTO_REALIMAG_IMAGREAL	\VSRes6,\VSRes8,vs8,vs9 
  AGGREGATE_REALS_IMAGES	vs2,vs3,vs4,vs5	
  AGGREGATE_REALS_IMAGES	vs6,vs7,vs8,vs9  
  MULT_APLHA_PART1	vs2,vs4, vs14,vs15
  MULT_APLHA_PART1	vs6,vs8, vs16,vs17
  MULT_APLHA_PART2	vs2,vs4, vs14,vs15 
  MULT_APLHA_PART2	vs6,vs8,vs16,vs17
  UNPACK_FOR_STORE	vs14,vs15,vs7,vs9
  UNPACK_FOR_STORE	vs16,vs17,vs3,vs5
  STORE_COUPLE	\BASE_REG,\LOFFSET,vs7,vs9
  STORE_COUPLE	\BASE_REG,(\LOFFSET+32),vs3,vs5
.endm


.macro SAVE2  VSRes1,VSRes2,VSRes3,VSRes4,BASE_REG,LOFFSET
  RESULT_INTO_REALREAL_IMAGEIMAGE \VSRes1,\VSRes3,vs2,vs3
  LOAD_COUPLE_AS_RR_II	vs14,vs15,vs18,vs19,\BASE_REG,\LOFFSET
  RESULT_INTO_REALIMAG_IMAGREAL	\VSRes2,\VSRes4,vs4,vs5	
  AGGREGATE_REALS_IMAGES	vs2,vs3,vs4,vs5	
  MULT_APLHA_PART1	vs2,vs4, vs14,vs15	
  MULT_APLHA_PART2	vs2,vs4, vs14,vs15  
  UNPACK_FOR_STORE	vs14,vs15,vs7,vs9	
  STORE_COUPLE	\BASE_REG,\LOFFSET,vs7,vs9  
.endm


.macro SAVE1  VSRes1,VSRes2,BASE_REG,LOFFSET
  RESULT_INTO_REALREAL_IMAGEIMAGE \VSRes1,\VSRes1,vs2,vs3
#ifndef TRMMKERNEL 
  lxv	vs18,	(\LOFFSET)(\BASE_REG) 
  xxmrgld  vs14,vs18,vs18
  xxmrghd  vs15,vs18,vs18	
#endif	
  RESULT_INTO_REALIMAG_IMAGREAL	\VSRes2,\VSRes2,vs4,vs5	
  AGGREGATE_REALS_IMAGES	vs2,vs3,vs4,vs5	
  MULT_APLHA_PART1	vs2,vs4, vs14,vs15	
  MULT_APLHA_PART2	vs2,vs4, vs14,vs15  
  UNPACK_FOR_STORE	vs14,vs15,vs7,vs9 
  xxmrghd  vs7,vs15,vs14	
  stxv	vs7,	(\LOFFSET)(\BASE_REG) 
.endm

/**********************************************************************************************
* Macros for N=2 and M=8
**********************************************************************************************/

.macro Zero2x8
	xxlxor	vs32,	vs32,	vs32
	xxlxor	vs33,	vs33,	vs33
	xxlxor	vs34,	vs34,	vs34
	xxlxor	vs35,	vs35,	vs35
	xxlxor	vs36,	vs36,	vs36
	xxlxor	vs37,	vs37,	vs37
	xxlxor	vs38,	vs38,	vs38
	xxlxor	vs39,	vs39,	vs39
	xxlxor	vs40,	vs40,	vs40
	xxlxor	vs41,	vs41,	vs41
	xxlxor	vs42,	vs42,	vs42
	xxlxor	vs43,	vs43,	vs43
	xxlxor	vs44,	vs44,	vs44
	xxlxor	vs45,	vs45,	vs45
	xxlxor	vs46,	vs46,	vs46
	xxlxor	vs47,	vs47,	vs47
	xxlxor	vs48,	vs48,	vs48
	xxlxor	vs49,	vs49,	vs49
	xxlxor	vs50,	vs50,	vs50
	xxlxor	vs51,	vs51,	vs51
	xxlxor	vs52,	vs52,	vs52
	xxlxor	vs53,	vs53,	vs53
	xxlxor	vs54,	vs54,	vs54
	xxlxor	vs55,	vs55,	vs55
	xxlxor	vs56,	vs56,	vs56
	xxlxor	vs57,	vs57,	vs57
	xxlxor	vs58,	vs58,	vs58
	xxlxor	vs59,	vs59,	vs59
	xxlxor	vs60,	vs60,	vs60
	xxlxor	vs61,	vs61,	vs61
	xxlxor	vs62,	vs62,	vs62
	xxlxor	vs63,	vs63,	vs63
.endm

.macro LOAD2x8 Zero

	lxv	vs16,	0(BO)	// load real imag from B
	lxv	vs18,	16(BO)	// load real,imag from B
	xxswapd	vs17, vs16
	xxswapd	vs19, vs18

	lxv	vs0,	0(AO)	// load real,imag from A
	lxv	vs1,	16(AO)	// load real,imag from A
	lxv	vs2,	32(AO)	// load real,imag from A
	lxv	vs3,	48(AO)	// load real,imag from A

	lxv	vs4,	64(AO)	// load real,imag from A
	lxv	vs5,	80(AO)	// load real,imag from A
	lxv	vs6,	96(AO)	// load real,imag from A
	lxv	vs7,	112(AO)	// load real,imag from A

.if \Zero==1
	Zero2x8
.endif

.endm

.macro END2x8_NORMAL
	END2x8 AO,BO,128,32
.endm

.macro END2x8	AREG, BREG, OffsetA, OffsetB

.if \OffsetB != 0
	addi	\BREG, \BREG, \OffsetB
.endif
.if \OffsetA != 0
	addi	\AREG, \AREG, \OffsetA
.endif

	xvmaddadp	vs32,	vs0,	vs16
	xvmaddadp	vs48,	vs0,	vs18

	xvmaddadp	vs34,	vs1,	vs16
	xvmaddadp	vs50,	vs1,	vs18

	xvmaddadp	vs36,	vs2,	vs16
	xvmaddadp	vs52,	vs2,	vs18

	xvmaddadp	vs38,	vs3,	vs16
	xvmaddadp	vs54,	vs3,	vs18

	xvmaddadp	vs40,	vs4,	vs16
	xvmaddadp	vs56,	vs4,	vs18

	xvmaddadp	vs42,	vs5,	vs16
	xvmaddadp	vs58,	vs5,	vs18

	xvmaddadp	vs44,	vs6,	vs16
	xvmaddadp	vs60,	vs6,	vs18

	xvmaddadp	vs46,	vs7,	vs16
	xvmaddadp	vs62,	vs7,	vs18


	xvmaddadp	vs33,	vs0,	vs17
	xvmaddadp	vs49,	vs0,	vs19

	xvmaddadp	vs35,	vs1,	vs17
	xvmaddadp	vs51,	vs1,	vs19

	xvmaddadp	vs37,	vs2,	vs17
	xvmaddadp	vs53,	vs2,	vs19

	xvmaddadp	vs39,	vs3,	vs17
	xvmaddadp	vs55,	vs3,	vs19

	xvmaddadp	vs41,	vs4,	vs17
	xvmaddadp	vs57,	vs4,	vs19

	xvmaddadp	vs43,	vs5,	vs17
	xvmaddadp	vs59,	vs5,	vs19

	xvmaddadp	vs45,	vs6,	vs17
	xvmaddadp	vs61,	vs6,	vs19

	xvmaddadp	vs47,	vs7,	vs17
	xvmaddadp	vs63,	vs7,	vs19

.endm

.macro KERNEL2x8_L	OffsetA,OffsetB, Index,IsLast
  KERNEL2x8_2	AO,BO,	\OffsetA,\OffsetB, \Index,\IsLast ,0
.endm

.macro KERNEL2x8_E	OffsetA,OffsetB, Index,IsLast
  KERNEL2x8_2	AO,BO,	\OffsetA,\OffsetB, \Index,\IsLast ,1
.endm


.macro KERNEL2x8_2	AREG,BREG,	OffsetA,OffsetB, Index,IsLast ,Complete

	lxv	vs20,	DISP4(\Index,	0+\OffsetB)(\BREG)	// load real,imag	from B
	lxv	vs22,	DISP4(\Index,16+\OffsetB)(\BREG)	// load real,imag  from B
	xvmaddadp	vs32,	vs0,	vs16
	xvmaddadp	vs48,	vs0,	vs18
	xvmaddadp	vs33,	vs0,	vs17
	xvmaddadp	vs49,	vs0,	vs19

	xxswapd	vs21, vs20
	xxswapd	vs23, vs22

	xvmaddadp	vs34,	vs1,	vs16
	xvmaddadp	vs50,	vs1,	vs18

	lxv	vs8,	DISP16(\Index, 0 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs9,	DISP16(\Index,16 + \OffsetA)(\AREG)	// load real,imag from A

	xvmaddadp	vs35,	vs1,	vs17
	xvmaddadp	vs51,	vs1,	vs19

	lxv	vs10,	DISP16(\Index,32 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs11,	DISP16(\Index,48 + \OffsetA)(\AREG)	// load real,imag from A

	xvmaddadp	vs36,	vs2,	vs16
	xvmaddadp	vs52,	vs2,	vs18

	lxv	vs12,	DISP16(\Index, 64 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs13,	DISP16(\Index,64+16 + \OffsetA)(\AREG)	// load real,imag from A

	xvmaddadp	vs37,	vs2,	vs17
	xvmaddadp	vs53,	vs2,	vs19

	lxv	vs14,	DISP16(\Index,64+32 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs15,	DISP16(\Index,64+48 + \OffsetA)(\AREG)	// load real,imag from A

.if \IsLast==1
.if \Complete==1 
	addi	\AREG, \AREG, DISP16(\Index,128+\OffsetA)
	addi	\BREG, \BREG,  DISP4(\Index,32+\OffsetB)
.endif
.endif


	xvmaddadp	vs38,	vs3,	vs16
	xvmaddadp	vs54,	vs3,	vs18

.if \Complete==0
	lxv	vs0,	DISP16(\Index,128+ + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs1,	DISP16(\Index,128+16 + \OffsetA)(\AREG)	// load real,imag from A
.endif


	xvmaddadp	vs39,	vs3,	vs17
	xvmaddadp	vs55,	vs3,	vs19

.if \Complete==0
	lxv	vs2,	DISP16(\Index,128+32 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs3,	DISP16(\Index,128+48 + \OffsetA)(\AREG)	// load real,imag from A
.endif
	xvmaddadp	vs40,	vs4,	vs16
	xvmaddadp	vs56,	vs4,	vs18

	xvmaddadp	vs41,	vs4,	vs17
	xvmaddadp	vs57,	vs4,	vs19

	xvmaddadp	vs42,	vs5,	vs16
	xvmaddadp	vs58,	vs5,	vs18
	xvmaddadp	vs43,	vs5,	vs17
	xvmaddadp	vs59,	vs5,	vs19

.if \Complete==0
	lxv	vs4,	DISP16(\Index, 192 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs5,	DISP16(\Index,192 +16 + \OffsetA)(\AREG)	// load real,imag from A
.endif

	xvmaddadp	vs44,	vs6,	vs16
	xvmaddadp	vs60,	vs6,	vs18
	xvmaddadp	vs45,	vs6,	vs17
	xvmaddadp	vs61,	vs6,	vs19

	xvmaddadp	vs46,	vs7,	vs16
	xvmaddadp	vs62,	vs7,	vs18
	xvmaddadp	vs47,	vs7,	vs17
	xvmaddadp	vs63,	vs7,	vs19

.if \Complete==0
	lxv	vs6,	DISP16(\Index,192 +32 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs7,	DISP16(\Index,192 +48 + \OffsetA)(\AREG)	// load real,imag from A
.endif

	xvmaddadp	vs32,	vs8,	vs20
	xvmaddadp	vs48,	vs8,	vs22
.if \Complete==0
	lxv	vs16,	DISP4(\Index, 32+\OffsetB)(\BREG)	// load real imag from B
	lxv	vs18,	DISP4(\Index, 48+\OffsetB)(\BREG)	// load real,imag from B
.endif
.if \Complete==0
.if \IsLast==1 
	addi	\AREG, \AREG, DISP16(\Index,256)
	addi	\BREG, \BREG,  DISP4(\Index,64)
.endif

.endif
	xvmaddadp	vs33,	vs8,	vs21
	xvmaddadp	vs49,	vs8,	vs23

.if \Complete==0
	xxswapd	vs17, vs16
	xxswapd	vs19, vs18
.endif

	xvmaddadp	vs34,	vs9,	vs20
	xvmaddadp	vs50,	vs9,	vs22
	xvmaddadp	vs35,	vs9,	vs21
	xvmaddadp	vs51,	vs9,	vs23

	xvmaddadp	vs36,	vs10,	vs20
	xvmaddadp	vs52,	vs10,	vs22
	xvmaddadp	vs37,	vs10,	vs21
	xvmaddadp	vs53,	vs10,	vs23

	xvmaddadp	vs38,	vs11,	vs20
	xvmaddadp	vs54,	vs11,	vs22
	xvmaddadp	vs39,	vs11,	vs21
	xvmaddadp	vs55,	vs11,	vs23

	xvmaddadp	vs40,	vs12,	vs20
	xvmaddadp	vs56,	vs12,	vs22
	xvmaddadp	vs41,	vs12,	vs21
	xvmaddadp	vs57,	vs12,	vs23

	xvmaddadp	vs42,	vs13,	vs20
	xvmaddadp	vs58,	vs13,	vs22
	xvmaddadp	vs43,	vs13,	vs21
	xvmaddadp	vs59,	vs13,	vs23

	xvmaddadp	vs44,	vs14,	vs20
	xvmaddadp	vs60,	vs14,	vs22
	xvmaddadp	vs45,	vs14,	vs21
	xvmaddadp	vs61,	vs14,	vs23

	xvmaddadp	vs46,	vs15,	vs20
	xvmaddadp	vs62,	vs15,	vs22
	xvmaddadp	vs47,	vs15,	vs21
	xvmaddadp	vs63,	vs15,	vs23

.endm

.macro KERNEL2x8
  LOAD2x8 0
  END2x8  AO, BO, 128,32
.endm

.macro SAVE2x8

	add	T1, CO ,LDC 
	SAVE8  vs32,vs33,vs34,vs35,vs36,vs37,vs38,vs39,vs40,vs41,vs42,vs43,vs44,vs45,vs46,vs47,CO,0
	SAVE8  vs48,vs49,vs50,vs51,vs52,vs53,vs54,vs55,vs56,vs57,vs58,vs59,vs60,vs61,vs62,vs63,T1,0  
	addi	CO, CO, 128

.endm

/**********************************************************************************************
* Macros for N=2 and M=4
**********************************************************************************************/

.macro Zero2x4
	xxlxor	vs32,	vs32,	vs32
	xxlxor	vs33,	vs33,	vs33
	xxlxor	vs34,	vs34,	vs34
	xxlxor	vs35,	vs35,	vs35
	xxlxor	vs36,	vs36,	vs36
	xxlxor	vs37,	vs37,	vs37
	xxlxor	vs38,	vs38,	vs38
	xxlxor	vs39,	vs39,	vs39
	xxlxor	vs40,	vs40,	vs40
	xxlxor	vs41,	vs41,	vs41
	xxlxor	vs42,	vs42,	vs42
	xxlxor	vs43,	vs43,	vs43
	xxlxor	vs44,	vs44,	vs44
	xxlxor	vs45,	vs45,	vs45
	xxlxor	vs46,	vs46,	vs46
	xxlxor	vs47,	vs47,	vs47
.endm

.macro LOAD2x4 Zero

	lxv	vs16,	0(BO)	// load real imag from B
	lxv	vs18,	16(BO)	// load real,imag from B
	xxswapd	vs17, vs16
	xxswapd	vs19, vs18

	lxv	vs0,	0(AO)	// load real,imag from A
	lxv	vs1,	16(AO)	// load real,imag from A
	lxv	vs2,	32(AO)	// load real,imag from A
	lxv	vs3,	48(AO)	// load real,imag from A

.if \Zero==1
	Zero2x4
.endif

.endm

.macro END2x4_NORMAL
	END2x4 AO,BO,64,32
.endm

.macro END2x4	AREG, BREG, OffsetA, OffsetB

.if \OffsetB != 0
	addi	\BREG, \BREG, \OffsetB
.endif
.if \OffsetA != 0
	addi	\AREG, \AREG, \OffsetA
.endif

	xvmaddadp	vs32,	vs0,	vs16
	xvmaddadp	vs33,	vs0,	vs17
	xvmaddadp	vs40,	vs0,	vs18
	xvmaddadp	vs41,	vs0,	vs19

	xvmaddadp	vs34,	vs1,	vs16
	xvmaddadp	vs35,	vs1,	vs17
	xvmaddadp	vs42,	vs1,	vs18
	xvmaddadp	vs43,	vs1,	vs19
	
	xvmaddadp	vs36,	vs2,	vs16
	xvmaddadp	vs37,	vs2,	vs17
	xvmaddadp	vs44,	vs2,	vs18
	xvmaddadp	vs45,	vs2,	vs19

	xvmaddadp	vs38,	vs3,	vs16
	xvmaddadp	vs39,	vs3,	vs17
	xvmaddadp	vs46,	vs3,	vs18
	xvmaddadp	vs47,	vs3,	vs19

.endm

.macro KERNEL2x4_L	OffsetA,OffsetB, Index,IsLast
  KERNEL2x4_2	AO,BO,	\OffsetA,\OffsetB, \Index,\IsLast ,0
.endm

.macro KERNEL2x4_E	OffsetA,OffsetB, Index,IsLast
  KERNEL2x4_2	AO,BO,	\OffsetA,\OffsetB, \Index,\IsLast ,1
.endm

.macro KERNEL2x4_2	AREG,BREG,	OffsetA,OffsetB, Index,IsLast ,Complete

	lxv	vs20,	DISP4(\Index,	0+\OffsetB)(\BREG)	// load real,imag	from B
	lxv	vs22,	DISP4(\Index,16+\OffsetB)(\BREG)	// load real,imag  from B
 
	xvmaddadp	vs32,	vs0,	vs16
	xvmaddadp	vs33,	vs0,	vs17
	xxswapd	vs21, vs20
	xxswapd	vs23, vs22	
	lxv	vs8,	DISP8(\Index, 0 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs9,	DISP8(\Index,16 + \OffsetA)(\AREG)	// load real,imag from A
	xvmaddadp	vs40,	vs0,	vs18
	xvmaddadp	vs41,	vs0,	vs19
	lxv	vs10,	DISP8(\Index,32 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs11,	DISP8(\Index,48 + \OffsetA)(\AREG)	// load real,imag from A
.if \IsLast==1
.if \Complete==1
	addi	\BREG, \BREG,  DISP4(\Index,32+\OffsetB) 
	addi	\AREG, \AREG,  DISP8(\Index,64+\OffsetA) 
.endif
.endif

	xvmaddadp	vs34,	vs1,	vs16
	xvmaddadp	vs35,	vs1,	vs17
	xvmaddadp	vs42,	vs1,	vs18
	xvmaddadp	vs43,	vs1,	vs19
	
	xvmaddadp	vs36,	vs2,	vs16
	xvmaddadp	vs37,	vs2,	vs17
.if \Complete==0
	lxv	vs0,	DISP8(\Index,64+  \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs1,	DISP8(\Index,64+16 + \OffsetA)(\AREG)	// load real,imag from A
.endif	
	xvmaddadp	vs44,	vs2,	vs18
	xvmaddadp	vs45,	vs2,	vs19
	
	xvmaddadp	vs38,	vs3,	vs16
	xvmaddadp	vs39,	vs3,	vs17
	xvmaddadp	vs46,	vs3,	vs18
	xvmaddadp	vs47,	vs3,	vs19


.if \Complete==0	
	lxv	vs2,	DISP8(\Index,64+32 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs3,	DISP8(\Index,64+48 + \OffsetA)(\AREG)	// load real,imag from A

.endif
	xvmaddadp	vs32,	vs8,	vs20
	xvmaddadp	vs33,	vs8,	vs21
.if \Complete==0
	lxv	vs16,	DISP4(\Index, 32+\OffsetB)(\BREG)	// load real imag from B
	lxv	vs18,	DISP4(\Index, 48+\OffsetB)(\BREG)	// load real,imag from B
.if \IsLast==1 
	addi	\AREG, \AREG,  DISP8(\Index,128)
	addi	\BREG, \BREG,  DISP4(\Index,64) 
.endif    
.endif

.if \Complete==0
	xxswapd	vs17, vs16
	xxswapd	vs19, vs18
.endif
 
	xvmaddadp	vs40,	vs8,	vs22
	xvmaddadp	vs41,	vs8,	vs23

	xvmaddadp	vs34,	vs9,	vs20
	xvmaddadp	vs35,	vs9,	vs21
	xvmaddadp	vs42,	vs9,	vs22
	xvmaddadp	vs43,	vs9,	vs23

	xvmaddadp	vs36,	vs10,	vs20
	xvmaddadp	vs37,	vs10,	vs21
	xvmaddadp	vs44,	vs10,	vs22
	xvmaddadp	vs45,	vs10,	vs23

	xvmaddadp	vs38,	vs11,	vs20
	xvmaddadp	vs39,	vs11,	vs21
	xvmaddadp	vs46,	vs11,	vs22
	xvmaddadp	vs47,	vs11,	vs23

.endm

.macro KERNEL2x4
  LOAD2x4 0
  END2x4  AO, BO, 64,32
.endm

.macro SAVE2x4 
	add	T1, CO ,LDC 
	SAVE4  vs32,vs33,vs34,vs35,vs36,vs37,vs38,vs39,CO,0
	SAVE4  vs40,vs41,vs42,vs43,vs44,vs45,vs46,vs47,T1,0  
	addi	CO, CO, 64

.endm

/**********************************************************************************************
* Macros for N=2 and M=2
**********************************************************************************************/

.macro Zero2x2
	xxlxor	vs32,	vs32,	vs32
	xxlxor	vs33,	vs33,	vs33
	xxlxor	vs34,	vs34,	vs34
	xxlxor	vs35,	vs35,	vs35
	xxlxor	vs36,	vs36,	vs36
	xxlxor	vs37,	vs37,	vs37
	xxlxor	vs38,	vs38,	vs38
	xxlxor	vs39,	vs39,	vs39
.endm

.macro LOAD2x2 Zero

	lxv	vs16,	0(BO)	// load real imag from B
	lxv	vs18,	16(BO)	// load real,imag from B
	xxswapd	vs17, vs16
	xxswapd	vs19, vs18

	lxv	vs0,	0(AO)	// load real,imag from A
	lxv	vs1,	16(AO)	// load real,imag from A


.if \Zero==1
	Zero2x2
.endif 
.endm

.macro END2x2_NORMAL
	END2x2 AO,BO,32,32
.endm

.macro END2x2	AREG, BREG, OffsetA, OffsetB

.if \OffsetB != 0
	addi	\BREG, \BREG, \OffsetB
.endif
.if \OffsetA != 0
	addi	\AREG, \AREG, \OffsetA
.endif

	xvmaddadp	vs32,	vs0,	vs16
	xvmaddadp	vs33,	vs0,	vs17
	xvmaddadp	vs36,	vs0,	vs18
	xvmaddadp	vs37,	vs0,	vs19

	xvmaddadp	vs34,	vs1,	vs16
	xvmaddadp	vs35,	vs1,	vs17 
	xvmaddadp	vs38,	vs1,	vs18
	xvmaddadp	vs39,	vs1,	vs19

.endm

.macro KERNEL2x2_L	OffsetA,OffsetB, Index,IsLast
  KERNEL2x2_2	AO,BO,	\OffsetA,\OffsetB, \Index,\IsLast ,0
.endm

.macro KERNEL2x2_E	OffsetA,OffsetB, Index,IsLast
  KERNEL2x2_2	AO,BO,	\OffsetA,\OffsetB, \Index,\IsLast ,1
.endm

.macro KERNEL2x2_2	AREG,BREG,	OffsetA,OffsetB, Index,IsLast ,Complete

	lxv	vs20,	DISP4(\Index,	0+\OffsetB)(\BREG)	// load real,imag	from B
	lxv	vs22,	DISP4(\Index,16+\OffsetB)(\BREG)	// load real,imag  from B
	xvmaddadp	vs32,	vs0,	vs16
	xvmaddadp	vs33,	vs0,	vs17
	xxswapd	vs21, vs20
	xxswapd	vs23, vs22

	lxv	vs8,	DISP4(\Index, 0 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs9,	DISP4(\Index,16 + \OffsetA)(\AREG)	// load real,imag from A
.if \IsLast==1
.if \Complete==1
	addi	\AREG, \AREG,  DISP4(\Index,32+\OffsetA)
	addi	\BREG, \BREG,  DISP4(\Index,32+\OffsetB) 
.endif
.endif 
	xvmaddadp	vs36,	vs0,	vs18
	xvmaddadp	vs37,	vs0,	vs19

	xvmaddadp	vs34,	vs1,	vs16
	xvmaddadp	vs35,	vs1,	vs17 
	xvmaddadp	vs38,	vs1,	vs18
	xvmaddadp	vs39,	vs1,	vs19

.if \Complete==0
	lxv	vs0,	DISP4(\Index,32 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs1,	DISP4(\Index,48+ \OffsetA)(\AREG)	// load real,imag from A
.endif
.if \Complete==0
	lxv	vs16,	DISP4(\Index, 32+\OffsetB)(\BREG)	// load real imag from B
	lxv	vs18,	DISP4(\Index, 48+\OffsetB)(\BREG)	// load real,imag from B
.if \IsLast==1 
	addi	\AREG, \AREG,  DISP4(\Index,64)
	addi	\BREG, \BREG,  DISP4(\Index,64)
.endif	
.endif

	xvmaddadp	vs32,	vs8,	vs20
	xvmaddadp	vs33,	vs8,	vs21

.if \Complete==0
	xxswapd	vs17, vs16
	xxswapd	vs19, vs18
.endif 
	xvmaddadp	vs36,	vs8,	vs22
	xvmaddadp	vs37,	vs8,	vs23

	xvmaddadp	vs34,	vs9,	vs20
	xvmaddadp	vs35,	vs9,	vs21 

	xvmaddadp	vs38,	vs9,	vs22
	xvmaddadp	vs39,	vs9,	vs23

.endm

.macro KERNEL2x2
  LOAD2x2 0
  END2x2  AO, BO, 32,32
.endm

.macro SAVE2x2 
	add	T1, CO ,LDC 
	SAVE2  vs32,vs33,vs34,vs35,CO,0
	SAVE2  vs36,vs37,vs38,vs39,T1,0 
	addi	CO, CO, 32 
.endm

/**********************************************************************************************
* Macros for N=2 and M=1
**********************************************************************************************/

.macro Zero2x1
	xxlxor	vs32,	vs32,	vs32
	xxlxor	vs33,	vs33,	vs33
	xxlxor	vs34,	vs34,	vs34
	xxlxor	vs35,	vs35,	vs35
.endm

.macro LOAD2x1 Zero
	lxv	vs0,	0(AO)	// load real,imag from A

	lxv	vs16,	0(BO)	// load real imag from B
	lxv	vs18,	16(BO)	// load real,imag from B

	xxswapd	vs17, vs16
	xxswapd	vs19, vs18
.if \Zero==1
	Zero2x1
.endif 
.endm

.macro END2x1_NORMAL
	END2x1 AO,BO,16,32
.endm

.macro END2x1	AREG, BREG, OffsetA, OffsetB

.if \OffsetA != 0
	addi	\AREG, \AREG, \OffsetA
.endif
.if \OffsetB != 0
	addi	\BREG, \BREG, \OffsetB
.endif

	xvmaddadp	vs32,	vs0,	vs16
	xvmaddadp	vs33,	vs0,	vs17

	xvmaddadp	vs34,	vs0,	vs18
	xvmaddadp	vs35,	vs0,	vs19

.endm

.macro KERNEL2x1_L	OffsetA,OffsetB, Index,IsLast
  KERNEL2x1_2	AO,BO,	\OffsetA,\OffsetB, \Index,\IsLast ,0
.endm

.macro KERNEL2x1_E	OffsetA,OffsetB, Index,IsLast
  KERNEL2x1_2	AO,BO,	\OffsetA,\OffsetB, \Index,\IsLast ,1
.endm

.macro KERNEL2x1_2	AREG,BREG,	OffsetA,OffsetB, Index,IsLast ,Complete

	lxv	vs20,	DISP4(\Index,	0+\OffsetB)(\BREG)	// load real,imag	from B
	lxv	vs22,	DISP4(\Index,16+\OffsetB)(\BREG)	// load real,imag  from B

	lxv	vs8,	DISP2(\Index, 0 + \OffsetA)(\AREG)	// load real,imag from A

	xxswapd	vs21, vs20
	xxswapd	vs23, vs22
.if \IsLast==1
.if \Complete==1
	addi	\AREG, \AREG,  DISP2(\Index,16+\OffsetA)
	addi	\BREG, \BREG,  DISP4(\Index,32+\OffsetB) 
.endif
.endif

	xvmaddadp	vs32,	vs0,	vs16
	xvmaddadp	vs33,	vs0,	vs17

	xvmaddadp	vs34,	vs0,	vs18
	xvmaddadp	vs35,	vs0,	vs19

.if \Complete==0
	lxv	vs0,	DISP2(\Index,16 + \OffsetA)(\AREG)	// load real,imag from A

.endif
.if \Complete==0
	lxv	vs16,	DISP4(\Index, 32+\OffsetB)(\BREG)	// load real imag from B
	lxv	vs18,	DISP4(\Index, 48+\OffsetB)(\BREG)	// load real,imag from B
.if \IsLast==1 
	addi	\AREG, \AREG,  DISP2(\Index,32)
	addi	\BREG, \BREG,  DISP4(\Index,64)
.endif	
.endif
 
.if \Complete==0
	xxswapd	vs17, vs16
	xxswapd	vs19, vs18
.endif

	xvmaddadp	vs32,	vs8,	vs20
	xvmaddadp	vs33,	vs8,	vs21

	xvmaddadp	vs34,	vs8,	vs22
	xvmaddadp	vs35,	vs8,	vs23

.endm

.macro KERNEL2x1
  LOAD2x1 0
  END2x1  AO, BO, 16,32
.endm

.macro SAVE2x1
	add	T1, CO ,LDC 
	SAVE1  vs32,vs33,CO,0
	SAVE1  vs34,vs35,T1,0  
	addi	CO, CO, 16 
.endm

/**********************************************************************************************
* Macros for N=1 and M=8
**********************************************************************************************/
.macro Zero1x8
	xxlxor	vs32,	vs32,	vs32
	xxlxor	vs33,	vs33,	vs33
	xxlxor	vs34,	vs34,	vs34
	xxlxor	vs35,	vs35,	vs35
	xxlxor	vs36,	vs36,	vs36
	xxlxor	vs37,	vs37,	vs37
	xxlxor	vs38,	vs38,	vs38
	xxlxor	vs39,	vs39,	vs39
	xxlxor	vs40,	vs40,	vs40
	xxlxor	vs41,	vs41,	vs41
	xxlxor	vs42,	vs42,	vs42
	xxlxor	vs43,	vs43,	vs43
	xxlxor	vs44,	vs44,	vs44
	xxlxor	vs45,	vs45,	vs45
	xxlxor	vs46,	vs46,	vs46
	xxlxor	vs47,	vs47,	vs47
.endm

.macro LOAD1x8 Zero

	lxv	vs16,	0(BO)	// load real imag from B
	xxswapd	vs17, vs16
	lxv	vs0,	0(AO)	// load real,imag from A
	lxv	vs1,	16(AO)	// load real,imag from A
	lxv	vs2,	32(AO)	// load real,imag from A
	lxv	vs3,	48(AO)	// load real,imag from A

	lxv	vs4,	64(AO)	// load real,imag from A
	lxv	vs5,	80(AO)	// load real,imag from A
	lxv	vs6,	96(AO)	// load real,imag from A
	lxv	vs7,	112(AO)	// load real,imag from A

.if \Zero==1
	Zero1x8
.endif

.endm

.macro END1x8_NORMAL
	END1x8 AO,BO,128,16
.endm

.macro END1x8	AREG, BREG, OffsetA, OffsetB

.if \OffsetB != 0
	addi	\BREG, \BREG, \OffsetB
.endif
.if \OffsetA != 0
	addi	\AREG, \AREG, \OffsetA
.endif

	xvmaddadp	vs32,	vs0,	vs16
	xvmaddadp	vs33,	vs0,	vs17
	xvmaddadp	vs34,	vs1,	vs16
	xvmaddadp	vs35,	vs1,	vs17
	xvmaddadp	vs36,	vs2,	vs16
	xvmaddadp	vs37,	vs2,	vs17
	xvmaddadp	vs38,	vs3,	vs16
	xvmaddadp	vs39,	vs3,	vs17
	xvmaddadp	vs40,	vs4,	vs16
	xvmaddadp	vs41,	vs4,	vs17
	xvmaddadp	vs42,	vs5,	vs16
	xvmaddadp	vs43,	vs5,	vs17
	xvmaddadp	vs44,	vs6,	vs16
	xvmaddadp	vs45,	vs6,	vs17
	xvmaddadp	vs46,	vs7,	vs16
	xvmaddadp	vs47,	vs7,	vs17

.endm

.macro KERNEL1x8_L	OffsetA,OffsetB, Index,IsLast
  KERNEL1x8_2	AO,BO,	\OffsetA,\OffsetB, \Index,\IsLast ,0
.endm

.macro KERNEL1x8_E	OffsetA,OffsetB, Index,IsLast
  KERNEL1x8_2	AO,BO,	\OffsetA,\OffsetB, \Index,\IsLast ,1
.endm

.macro KERNEL1x8_2	AREG,BREG,	OffsetA,OffsetB, Index,IsLast ,Complete

	lxv	vs20,	DISP2(\Index,	0+\OffsetB)(\BREG)	// load real,imag	from B
	xxswapd	vs21, vs20


	lxv	vs8,	DISP16(\Index, 0 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs9,	DISP16(\Index,16 + \OffsetA)(\AREG)	// load real,imag from A
	xvmaddadp	vs32,	vs0,	vs16
	xvmaddadp	vs33,	vs0,	vs17	
	lxv	vs10,	DISP16(\Index,32 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs11,	DISP16(\Index,48 + \OffsetA)(\AREG)	// load real,imag from A
	xvmaddadp	vs34,	vs1,	vs16
	xvmaddadp	vs35,	vs1,	vs17
	lxv	vs12,	DISP16(\Index, 64 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs13,	DISP16(\Index,64+16 + \OffsetA)(\AREG)	// load real,imag from A
	xvmaddadp	vs36,	vs2,	vs16
	xvmaddadp	vs37,	vs2,	vs17
	lxv	vs14,	DISP16(\Index,64+32 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs15,	DISP16(\Index,64+48 + \OffsetA)(\AREG)	// load real,imag from A

	xvmaddadp	vs38,	vs3,	vs16
	xvmaddadp	vs39,	vs3,	vs17
.if \Complete==0
	lxv	vs0,	DISP16(\Index,128+ + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs1,	DISP16(\Index,128+16 + \OffsetA)(\AREG)	// load real,imag from A
.endif	
	xvmaddadp	vs40,	vs4,	vs16
	xvmaddadp	vs41,	vs4,	vs17
.if \Complete==0	
	lxv	vs2,	DISP16(\Index,128+32 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs3,	DISP16(\Index,128+48 + \OffsetA)(\AREG)	// load real,imag from A
.endif	
	xvmaddadp	vs42,	vs5,	vs16
	xvmaddadp	vs43,	vs5,	vs17
	xvmaddadp	vs44,	vs6,	vs16
	xvmaddadp	vs45,	vs6,	vs17
.if \Complete==0
	lxv	vs4,	DISP16(\Index, 192 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs5,	DISP16(\Index,192 +16 + \OffsetA)(\AREG)	// load real,imag from A
.endif	
	xvmaddadp	vs46,	vs7,	vs16
	xvmaddadp	vs47,	vs7,	vs17


	xvmaddadp	vs32,	vs8,	vs20
	xvmaddadp	vs33,	vs8,	vs21
.if \Complete==0	
	lxv	vs6,	DISP16(\Index,192 +32 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs7,	DISP16(\Index,192 +48 + \OffsetA)(\AREG)	// load real,imag from A 
.endif
	xvmaddadp	vs34,	vs9,	vs20
	xvmaddadp	vs35,	vs9,	vs21
.if \Complete==0
	lxv	vs16,	DISP2(\Index, 16+\OffsetB)(\BREG)	// load real imag from B
	xxswapd	vs17,vs16
.endif
.if \IsLast==1
.if \Complete==1
	addi	\AREG, \AREG, DISP16(\Index,128+\OffsetA)
	addi	\BREG, \BREG, DISP2(\Index,16+\OffsetB)
.else
	addi	\AREG, \AREG, DISP16(\Index,256)
	addi	\BREG, \BREG, DISP2(\Index,32)
.endif
.endif
	xvmaddadp	vs36,	vs10,	vs20
	xvmaddadp	vs37,	vs10,	vs21

	xvmaddadp	vs38,	vs11,	vs20
	xvmaddadp	vs39,	vs11,	vs21

	xvmaddadp	vs40,	vs12,	vs20
	xvmaddadp	vs41,	vs12,	vs21
	xvmaddadp	vs42,	vs13,	vs20
	xvmaddadp	vs43,	vs13,	vs21
	xvmaddadp	vs44,	vs14,	vs20
	xvmaddadp	vs45,	vs14,	vs21
	xvmaddadp	vs46,	vs15,	vs20
	xvmaddadp	vs47,	vs15,	vs21

.endm

.macro KERNEL1x8
  LOAD1x8 0
  END1x8  AO, BO, 128,16
.endm

.macro SAVE1x8

	SAVE8  vs32,vs33,vs34,vs35,vs36,vs37,vs38,vs39,vs40,vs41,vs42,vs43,vs44,vs45,vs46,vs47,CO,0 
	addi	CO, CO, 128

.endm

/**********************************************************************************************
* Macros for N=1 and M=4
**********************************************************************************************/

.macro Zero1x4
	xxlxor	vs32,	vs32,	vs32
	xxlxor	vs33,	vs33,	vs33
	xxlxor	vs34,	vs34,	vs34
	xxlxor	vs35,	vs35,	vs35
	xxlxor	vs36,	vs36,	vs36
	xxlxor	vs37,	vs37,	vs37
	xxlxor	vs38,	vs38,	vs38
	xxlxor	vs39,	vs39,	vs39
.endm

.macro LOAD1x4 Zero

	lxv	vs16,	0(BO)	// load real imag from B
	xxswapd	vs17,vs16
	lxv	vs0,	0(AO)	// load real,imag from A
	lxv	vs1,	16(AO)	// load real,imag from A
	lxv	vs2,	32(AO)	// load real,imag from A
	lxv	vs3,	48(AO)	// load real,imag from A

.if \Zero==1
	Zero1x4
.endif

.endm

.macro END1x4_NORMAL
	END1x4 AO,BO,64,16
.endm

.macro END1x4	AREG, BREG, OffsetA, OffsetB

.if \OffsetB != 0
	addi	\BREG, \BREG, \OffsetB
.endif
.if \OffsetA != 0
	addi	\AREG, \AREG, \OffsetA
.endif

	xvmaddadp	vs32,	vs0,	vs16
	xvmaddadp	vs33,	vs0,	vs17
	xvmaddadp	vs34,	vs1,	vs16
	xvmaddadp	vs35,	vs1,	vs17
	xvmaddadp	vs36,	vs2,	vs16
	xvmaddadp	vs37,	vs2,	vs17
	xvmaddadp	vs38,	vs3,	vs16
	xvmaddadp	vs39,	vs3,	vs17

.endm

.macro KERNEL1x4_L	OffsetA,OffsetB, Index,IsLast
  KERNEL1x4_2	AO,BO,	\OffsetA,\OffsetB, \Index,\IsLast ,0
.endm

.macro KERNEL1x4_E	OffsetA,OffsetB, Index,IsLast
  KERNEL1x4_2	AO,BO,	\OffsetA,\OffsetB, \Index,\IsLast ,1
.endm

.macro KERNEL1x4_2	AREG,BREG,	OffsetA,OffsetB, Index,IsLast ,Complete

	lxv	vs20,	DISP2(\Index,	0+\OffsetB)(\BREG)	// load real,imag	from B
	xxswapd	vs21,vs20

	xvmaddadp	vs32,	vs0,	vs16
	xvmaddadp	vs33,	vs0,	vs17

	lxv	vs8,	DISP8(\Index, 0 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs9,	DISP8(\Index,16 + \OffsetA)(\AREG)	// load real,imag from A
	xvmaddadp	vs34,	vs1,	vs16
	xvmaddadp	vs35,	vs1,	vs17	
	lxv	vs10,	DISP8(\Index,32 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs11,	DISP8(\Index,48 + \OffsetA)(\AREG)	// load real,imag from A

	xvmaddadp	vs36,	vs2,	vs16
	xvmaddadp	vs37,	vs2,	vs17
	xvmaddadp	vs38,	vs3,	vs16
	xvmaddadp	vs39,	vs3,	vs17

	xvmaddadp	vs40,	vs0,	vs18
	xvmaddadp	vs41,	vs0,	vs19
	xvmaddadp	vs42,	vs1,	vs18
	xvmaddadp	vs43,	vs1,	vs19
	xvmaddadp	vs44,	vs2,	vs18
	xvmaddadp	vs45,	vs2,	vs19
	xvmaddadp	vs46,	vs3,	vs18
	xvmaddadp	vs47,	vs3,	vs19

.if \Complete==0
	lxv	vs0,	DISP8(\Index,64+  \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs1,	DISP8(\Index,64+16 + \OffsetA)(\AREG)	// load real,imag from A
.endif
.if \Complete==0	
	lxv	vs2,	DISP8(\Index,64+32 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs3,	DISP8(\Index,64+48 + \OffsetA)(\AREG)	// load real,imag from A
 
.endif
.if \Complete==0
	lxv	vs16,	DISP2(\Index, 16+\OffsetB)(\BREG)	// load real imag from B
	xxswapd	vs17,vs16
.endif
.if \IsLast==1
.if \Complete==1
	addi	\AREG, \AREG,  DISP8(\Index,64+\OffsetA)
	addi	\BREG, \BREG,  DISP2(\Index,16+\OffsetB)
.else
	addi	\AREG, \AREG,  DISP8(\Index,128)
	addi	\BREG, \BREG,  DISP2(\Index,32)
.endif
.endif

	xvmaddadp	vs32,	vs8,	vs20
	xvmaddadp	vs33,	vs8,	vs21
	xvmaddadp	vs34,	vs9,	vs20
	xvmaddadp	vs35,	vs9,	vs21
	xvmaddadp	vs36,	vs10,	vs20
	xvmaddadp	vs37,	vs10,	vs21
	xvmaddadp	vs38,	vs11,	vs20
	xvmaddadp	vs39,	vs11,	vs21

	xvmaddadp	vs40,	vs8,	vs22
	xvmaddadp	vs41,	vs8,	vs23
	xvmaddadp	vs42,	vs9,	vs22
	xvmaddadp	vs43,	vs9,	vs23
	xvmaddadp	vs44,	vs10,	vs22
	xvmaddadp	vs45,	vs10,	vs23
	xvmaddadp	vs46,	vs11,	vs22
	xvmaddadp	vs47,	vs11,	vs23

.endm

.macro KERNEL1x4
  LOAD1x4 0
  END1x4  AO, BO, 64,16
.endm

.macro SAVE1x4
	SAVE4  vs32,vs33,vs34,vs35,vs36,vs37,vs38,vs39,CO,0
	addi	CO, CO, 64

.endm

/**********************************************************************************************
* Macros for N=1 and M=2
**********************************************************************************************/

.macro Zero1x2
	xxlxor	vs32,	vs32,	vs32
	xxlxor	vs33,	vs33,	vs33
	xxlxor	vs34,	vs34,	vs34
	xxlxor	vs35,	vs35,	vs35
.endm

.macro LOAD1x2 Zero

	lxv	vs16,	0(BO)	// load real imag from B
	xxswapd	vs17,vs16
	lxv	vs0,	0(AO)	// load real,imag from A
	lxv	vs1,	16(AO)	// load real,imag from A
 
.if \Zero==1
	Zero1x2
.endif

.endm

.macro END1x2_NORMAL
	END1x2 AO,BO,32,16
.endm

.macro END1x2	AREG, BREG, OffsetA, OffsetB

.if \OffsetB != 0
	addi	\BREG, \BREG, \OffsetB
.endif
.if \OffsetA != 0
	addi	\AREG, \AREG, \OffsetA
.endif

	xvmaddadp	vs32,	vs0,	vs16
	xvmaddadp	vs33,	vs0,	vs17
	xvmaddadp	vs34,	vs1,	vs16
	xvmaddadp	vs35,	vs1,	vs17

.endm

.macro KERNEL1x2_L	OffsetA,OffsetB, Index,IsLast
  KERNEL1x2_2	AO,BO,	\OffsetA,\OffsetB, \Index,\IsLast ,0
.endm

.macro KERNEL1x2_E	OffsetA,OffsetB, Index,IsLast
  KERNEL1x2_2	AO,BO,	\OffsetA,\OffsetB, \Index,\IsLast ,1
.endm

.macro KERNEL1x2_2	AREG,BREG,	OffsetA,OffsetB, Index,IsLast ,Complete

	lxv	vs20,	DISP2(\Index,	0+\OffsetB)(\BREG)	// load real,imag	from B
	xxswapd	vs21,vs20

	xvmaddadp	vs32,	vs0,	vs16
	xvmaddadp	vs33,	vs0,	vs17

	lxv	vs8,	DISP4(\Index, 0 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs9,	DISP4(\Index,16 + \OffsetA)(\AREG)	// load real,imag from A

	xvmaddadp	vs34,	vs1,	vs16
	xvmaddadp	vs35,	vs1,	vs17
.if \Complete==0
	lxv	vs0,	DISP4(\Index,32 + \OffsetA)(\AREG)	// load real,imag from A
	lxv	vs1,	DISP4(\Index,48+ \OffsetA)(\AREG)	// load real,imag from A 
.endif
.if \Complete==0
	lxv	vs16,	DISP2(\Index, 16+\OffsetB)(\BREG)	// load real imag from B
	xxswapd	vs17,vs16
.endif
.if \IsLast==1
.if \Complete==1
	addi	\AREG, \AREG,  DISP4(\Index,32+\OffsetA)
	addi	\BREG, \BREG,  DISP2(\Index,16+\OffsetB)
.else
	addi	\AREG, \AREG,  DISP4(\Index,64)
	addi	\BREG, \BREG,  DISP2(\Index,32)
.endif
.endif

	xvmaddadp	vs32,	vs8,	vs20
	xvmaddadp	vs33,	vs8,	vs21
	xvmaddadp	vs34,	vs9,	vs20
	xvmaddadp	vs35,	vs9,	vs21

.endm

.macro KERNEL1x2
  LOAD1x2 0
  END1x2  AO, BO, 32,16
.endm

.macro SAVE1x2
	SAVE2  vs32,vs33,vs34,vs35,CO,0
	addi	CO, CO, 32 
.endm

/**********************************************************************************************
* Macros for N=1 and M=1
**********************************************************************************************/

.macro Zero1x1
	xxlxor	vs32,	vs32,	vs32
	xxlxor	vs33,	vs33,	vs33
.endm

.macro LOAD1x1 Zero
	lxv	vs0,	0(AO)	// load real,imag from A

	lxv	vs16,	0(BO)	// load real imag from B
	xxswapd vs17,  vs16
.if \Zero==1
	Zero1x1
.endif
 
.endm

.macro END1x1_NORMAL
	END1x1 AO,BO,16,16
.endm

.macro END1x1	AREG, BREG, OffsetA, OffsetB

.if \OffsetA != 0
	addi	\AREG, \AREG, \OffsetA
.endif
.if \OffsetB != 0
	addi	\BREG, \BREG, \OffsetB
.endif

	xvmaddadp	vs32,	vs0,	vs16
	xvmaddadp	vs33,	vs0,	vs17


.endm

.macro KERNEL1x1_L	OffsetA,OffsetB, Index,IsLast
  KERNEL1x1_2	AO,BO,	\OffsetA,\OffsetB, \Index,\IsLast ,0
.endm

.macro KERNEL1x1_E	OffsetA,OffsetB, Index,IsLast
  KERNEL1x1_2	AO,BO,	\OffsetA,\OffsetB, \Index,\IsLast ,1
.endm

.macro KERNEL1x1_2	AREG,BREG,	OffsetA,OffsetB, Index,IsLast ,Complete
	lxv	vs20,	DISP2(\Index,	0+\OffsetB)(\BREG)	// load real,imag	from B
	xxswapd vs21,  vs20

	lxv	vs8,	DISP2(\Index, 0 + \OffsetA)(\AREG)	// load real,imag from A
 
	xvmaddadp	vs32,	vs0,	vs16
	xvmaddadp	vs33,	vs0,	vs17

.if \Complete==0
	lxv	vs0,	DISP2(\Index,16 + \OffsetA)(\AREG)	// load real,imag from A
.endif
.if \Complete==0
	lxv	vs16,	DISP2(\Index, 16+\OffsetB)(\BREG)	// load real imag from B
	xxswapd vs17,  vs16	
.endif

.if \IsLast==1
.if \Complete==1
	addi	\AREG, \AREG,  DISP2(\Index,16+\OffsetA)
	addi	\BREG, \BREG,  DISP2(\Index,16+\OffsetB)
.else
	addi	\AREG, \AREG,  DISP2(\Index,32)
	addi	\BREG, \BREG,  DISP2(\Index,32)
.endif
.endif

	xvmaddadp	vs32,	vs8,	vs20
	xvmaddadp	vs33,	vs8,	vs21


.endm

.macro KERNEL1x1
  LOAD1x1 0
  END1x1  AO, BO, 16,16

.endm

.macro SAVE1x1 
	SAVE1  vs32,vs33,CO,0
	addi	CO, CO, 16 
.endm

