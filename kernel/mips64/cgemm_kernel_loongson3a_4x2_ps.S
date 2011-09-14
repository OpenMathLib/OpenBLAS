##define REALNAME gemm 
#define ASSEMBLER
#include "common.h"

#define FETCH	ld
#define	STACKSIZE	192
#define gsLQC1(base,fq,ft,offset) .word(0x32<<26|base<<21|ft<<16|0x1<<15|offset<<6|0x1<<5|fq)
#define gsSQC1(base,fq,ft,offset) .word(0x3A<<26|base<<21|ft<<16|0x1<<15|offset<<6|0x1<<5|fq)


#####	Parameter registers	####
#define M	$4
#define	N	$5
#define	K	$6
#define A	$8
#define B	$9
#define C	$10
#define LDC	$11

####	Pointer A, B, C	####
#define AO	$12
#define BO	$13

#define CO1	$14
#define CO2	$15

#define PREA	$18
#define PREB	$19

####	Used registers	####
#define A1	$f0
#define A2	$f1
#define A3	$f2
#define A4	$f3
#define A5	$f4
#define	A6	$f5
#define	A7	$f6
#define	A8	$f7

#define B1	$f8
#define B2	$f9
#define B3	$f10
#define B4	$f11
#define B5	$f12
#define	B6	$f13
#define	B7	$f14
#define	B8	$f15

#define C11	$f16
#define C12	$f17
#define C21	$f18
#define C22	$f19
#define C31	$f20
#define C32	$f21
#define C41	$f22
#define C42	$f23
#define C13	$f24
#define C14	$f25
#define C23	$f26
#define C24	$f27
#define C33	$f28
#define C34	$f29
#define C43	$f30
#define C44	$f31

#define I	$2
#define J	$3
#define L	$7

####	Alpha register	####
#define ALPHA	$f15

#define F31 31
#define F30 30
#define F29 29
#define F28 28
#define F27 27
#define F26 26
#define F25 25
#define F24 24 
#define F23 23
#define F22 22
#define F21 21
#define F20 20
#define F19 19
#define F18 18
#define F17 17
#define F16 16 
#define F15 15
#define F14 14
#define F13 13
#define F12 12
#define F11 11
#define F10 10
#define F9 9
#define F8 8
#define F7 7
#define F6 6
#define F5 5
#define F4 4 
#define F3 3 
#define F2 2 
#define F1 1 
#define F0 0

#define	R12	12
#define	R13	13

#define R14	14
#define R15	15
#define	R16	16
#define	R17	17

#if defined(TRMMKERNEL)
#define	OFFSET	$23
#define	KK		$24
#define TEMP	$25
#endif


	PROLOGUE

	daddiu	$sp,$sp,-STACKSIZE

	sd	$16,   0($sp)
	sd	$17,   8($sp)
	sd	$18,  16($sp)
	sd	$19,  24($sp)
	sd	$20,  32($sp)
	sd	$21,  40($sp)
	sd	$22,  48($sp)

	ST	$f24, 56($sp)
	ST	$f25, 64($sp)
	ST	$f26, 72($sp)
	ST	$f27, 80($sp)
	ST	$f28, 88($sp)

#if defined(TRMMKERNEL)
	sd	$23,  96($sp)
	sd	$24, 104($sp)
	sd	$25, 112($sp)

	LDARG	OFFSET, 160($sp)
#endif

#ifndef __64BIT__
	ST	$f20,120($sp)
	ST	$f21,128($sp)
	ST	$f22,136($sp)
	ST	$f23,144($sp)
#endif

	.align	4
.L2:
	dsra	J, N, 1				#	NR=2
	ST		$f15, 152($sp)

	dsll	LDC, LDC, ZBASE_SHIFT#	LDC*SIZE
	blez	J, .L1
	ST		$f16, 160($sp)

.L24:
	dsra	I, M, 2				#	MR=8
	move	AO, A				#	Reset A
	move	CO1, C

	daddu	CO2, C,   LDC
	blez	I, .L22
	daddu	C,   CO2, LDC

	.align	4
.L241:
	move	BO, B				#	Reset	B
	dsra	L, K, 2				#	UnRoll	K=64

	MTC		$0, C11				#	CLEAR REAULTS REGISTERS
	MOV		C12, C11		

	MOV		C21, C11
	MOV		C22, C11
	
	MOV		C31, C11
	MOV		C32, C11
	gsLQC1(R13, F9, F8, 0)		#	B1 B2

	MOV		C41, C11
	MOV		C42, C11
	gsLQC1(R12, F1, F0, 0)		#	A1 A2

	MOV		C13, C11
	MOV		C14, C11
	gsLQC1(R12, F3, F2, 1)		#	A3 A4

	MOV		C23, C11
	FETCH	$0, 0 * SIZE(CO1)
	MOV		C24, C11
	FETCH	$0, 4 * SIZE(CO1)
	
	MOV		C33, C11
	FETCH	$0, 0 * SIZE(CO2)
	MOV		C34, C11
	FETCH	$0, 4 * SIZE(CO2)
	
	MOV		C43, C11
	PLU		B3,	B1, B1

	MOV		C44, C11
	blez	L, .L242
	PLU		B4, B2, B2

.L2410:
	daddiu	L, L, -1
	gsLQC1(R13, F13, F12, 1)	#	B3 B4
	MADPS	C11, C11, A1, B1
	MADPS	C21, C21, A2, B1

	gsLQC1(R12, F5, F4, 2)		#	A5 A6
	MADPS	C12, C12, A1, B2
	MADPS	C22, C22, A2, B2

	gsLQC1(R12, F7, F6, 3)		#	A7 A8
	MADPS	C31, C31, A3, B1
	MADPS	C41, C41, A4, B1

	MADPS	C32, C32, A3, B2
	MADPS	C42, C42, A4, B2

	MADPS	C13, C13, A1, B3
	MADPS	C23, C23, A2, B3

	MADPS	C33, C33, A3, B3
	MADPS	C43, C43, A4, B3

	MADPS	C14, C14, A1, B4
	PLU		B7,	B5, B5

	MADPS	C24, C24, A2, B4
	PLU		B8, B6, B6

	MADPS	C34, C34, A3, B4
	MADPS	C44, C44, A4, B4

	gsLQC1(R13, F9, F8, 2)		#	B1 B2
	MADPS	C11, C11, A5, B5
	MADPS	C21, C21, A6, B5

	gsLQC1(R12, F1, F0, 4)		#	A1 A2
	MADPS	C12, C12, A5, B6
	MADPS	C22, C22, A6, B6

	gsLQC1(R12, F3, F2, 5)		#	A3 A4
	MADPS	C31, C31, A7, B5
	MADPS	C41, C41, A8, B5

	MADPS	C32, C32, A7, B6
	MADPS	C42, C42, A8, B6

	MADPS	C13, C13, A5, B7
	MADPS	C23, C23, A6, B7

	MADPS	C33, C33, A7, B7
	MADPS	C43, C43, A8, B7

	MADPS	C14, C14, A5, B8
	PLU		B3,	B1, B1

	MADPS	C24, C24, A6, B8
	PLU		B4, B2, B2

	MADPS	C34, C34, A7, B8
	MADPS	C44, C44, A8, B8
	
	gsLQC1(R13, F13, F12, 3)	#	B3 B4
	MADPS	C11, C11, A1, B1
	MADPS	C21, C21, A2, B1

	gsLQC1(R12, F5, F4, 6)		#	A5 A6
	MADPS	C12, C12, A1, B2
	MADPS	C22, C22, A2, B2

	gsLQC1(R12, F7, F6, 7)		#	A7 A8
	MADPS	C31, C31, A3, B1
	daddiu	BO, BO, 16 * SIZE	#	4KR*4NR
	MADPS	C41, C41, A4, B1

	MADPS	C32, C32, A3, B2
	MADPS	C42, C42, A4, B2
	daddiu	AO, AO, 32 * SIZE 	#	4KR*8MR

	MADPS	C13, C13, A1, B3
	MADPS	C23, C23, A2, B3

	MADPS	C33, C33, A3, B3
	MADPS	C43, C43, A4, B3

	MADPS	C14, C14, A1, B4
	PLU		B7,	B5, B5

	MADPS	C24, C24, A2, B4
	PLU		B8, B6, B6

	MADPS	C34, C34, A3, B4
	MADPS	C44, C44, A4, B4

	gsLQC1(R13, F9, F8, 0)		#	B1 B2
	MADPS	C11, C11, A5, B5
	MADPS	C21, C21, A6, B5

	gsLQC1(R12, F1, F0, 0)		#	A1 A2
	MADPS	C12, C12, A5, B6
	MADPS	C22, C22, A6, B6

	gsLQC1(R12, F3, F2, 1)		#	A3 A4
	MADPS	C31, C31, A7, B5
	MADPS	C41, C41, A8, B5

	MADPS	C32, C32, A7, B6
	MADPS	C42, C42, A8, B6

	MADPS	C13, C13, A5, B7
	MADPS	C23, C23, A6, B7

	MADPS	C33, C33, A7, B7
	MADPS	C43, C43, A8, B7

	MADPS	C14, C14, A5, B8
	PLU		B3,	B1, B1

	MADPS	C24, C24, A6, B8
	PLU		B4, B2, B2

	MADPS	C34, C34, A7, B8
	bgtz	L, .L2410
	MADPS	C44, C44, A8, B8


	.align	4
.L242:
	andi	L, K, 2
	blez	L, .L247
	NOP

	.align	4
.L247:
	andi	L, K, 1
	blez	L, .L240
	NOP


	.align	4
.L240:							#	Write Back
	daddiu	I, I, -1
	CVTU	A1, C11
	CVTU	A2, C21

	CVTU	A3, C31
	CVTU	A4, C41

	CVTU	A5, C13
	CVTU	A6, C23

	CVTU	A7, C33
	CVTU	A8, C43

	CVTU	B1, C12
	CVTU	B2, C22

	CVTU	B3, C32
	CVTU	B4, C42

	CVTU	B5, C14
	CVTU	B6, C24

	CVTU	B7, C34
	CVTU	B8, C44

#if   defined(NN) || defined(NT) || defined(TN) || defined(TT)
	/*	(a + bi) * (c + di) */
	SUB		C11, C11, A1		#	ac'+'bd
	SUB		C21, C21, A2
	LD		A1, 152($sp)		#	load alpha_r
#	LD		A1, 0 * SIZE(A)		#	load alpha_r
	SUB		C31, C31, A3
	LD		A2, 160($sp)		#	load alpha_i
#	LD		A2, 0 * SIZE(A)		#	load alpha_i

	SUB		C41, C41, A4
	ADD		C13, A5, C13		#	ad'+'cb
	ADD		C23, A6, C23
	ADD		C33, A7, C33
	ADD		C43, A8, C43
	SUB		C12, C12, B1
	SUB		C22, C22, B2
	SUB		C32, C32, B3
	SUB		C42, C42, B4
	ADD		C14, B5, C14
	ADD		C24, B6, C24
	ADD		C34, B7, C34
	ADD		C44, B8, C44

	LD		B1, 0 * SIZE(CO1)
	LD		B3, 2 * SIZE(CO1)
	LD		B5, 4 * SIZE(CO1)
	LD		B7, 6 * SIZE(CO1)
	LD		B2, 1 * SIZE(CO1)
	LD		B4, 3 * SIZE(CO1)
	LD		B6, 5 * SIZE(CO1)
	LD		B8, 7 * SIZE(CO1)

	MADD	B1, B1, C11, A1		#	A1 = alpha_r
	MADD	B3, B3, C21, A1
	MADD	B5, B5, C31, A1
	MADD	B7, B7, C41, A1
	MADD	B2, B2, C13, A1
	MADD	B4, B4, C23, A1
	MADD	B6, B6, C33, A1
	MADD	B8, B8, C43, A1
	NMSUB	B1, B1, C13, A2		#	A2 = alpha_i

	NMSUB	B3, B3, C23, A2
	NMSUB	B5, B5, C33, A2

	NMSUB	B7, B7, C43, A2
	MADD	B2, B2, C11, A2

	MADD	B4, B4, C21, A2
	MADD	B6, B6, C31, A2
	MADD	B8, B8, C41, A2

	LD		C13, 0 * SIZE(CO2)
	LD		C23, 2 * SIZE(CO2)
	LD		C33, 4 * SIZE(CO2)
	LD		C43, 6 * SIZE(CO2)
	LD		C11, 1 * SIZE(CO2)
	LD		C21, 3 * SIZE(CO2)
	LD		C31, 5 * SIZE(CO2)
	LD		C41, 7 * SIZE(CO2)

	MADD	C13, C13, C12, A1
	MADD	C23, C23, C22, A1

	MADD	C33, C33, C32, A1
	ST		B1, 0 * SIZE(CO1)

	MADD	C43, C43, C42, A1
	ST		B3, 2 * SIZE(CO1)

	MADD	C11, C11, C14, A1
	ST		B5, 4 * SIZE(CO1)

	MADD	C21, C21, C24, A1
	ST		B7, 6 * SIZE(CO1)

	MADD	C31, C31, C34, A1
	ST		B2, 1 * SIZE(CO1)

	MADD	C41, C41, C44, A1
	ST		B4, 3 * SIZE(CO1)

	NMSUB	C13, C13, C14, A2
	ST		B6, 5 * SIZE(CO1)

	NMSUB	C23, C23, C24, A2
	ST		B8, 7 * SIZE(CO1)

	NMSUB	C33, C33, C34, A2
	NMSUB	C43, C43, C44, A2

	MADD	C11, C11, C12, A2
	MADD	C21, C21, C22, A2
	
	MADD	C31, C31, C32, A2
	MADD	C41, C41, C42, A2

	ST		C13, 0 * SIZE(CO2)
	ST		C23, 2 * SIZE(CO2)
	ST		C33, 4 * SIZE(CO2)
	ST		C43, 6 * SIZE(CO2)
	ST		C11, 1 * SIZE(CO2)
	ST		C21, 3 * SIZE(CO2)
	ST		C31, 5 * SIZE(CO2)
	ST		C41, 7 * SIZE(CO2)
#endif

#if   defined(NR) || defined(NC) || defined(TR) || defined(TC)
	/*	(a + bi) * (c - di) */
	ADD		C11, A1, C11		#	ac'+'bd
	ADD		C21, A2, C21
#	LD		A1, 0 * SIZE(A)		#	load alpha_r
	LD		A1, 152($sp)		#	load alpha_r

	ADD		C31, A3, C31
	LD		A2, 160($sp)		#	load alpha_i
#	LD		A2, 0 * SIZE(A)		#	load alpha_r

	ADD		C41, A4, C41
	LD		B1, 0 * SIZE(CO1)

	SUB		C13, A5, C13		#	ad'+'cb
	LD		B3, 2 * SIZE(CO1)

	SUB		C23, A6, C23
	LD		B5, 4 * SIZE(CO1)

	SUB		C33, A7, C33
	LD		B7, 6 * SIZE(CO1)

	SUB		C43, A8, C43
	LD		B2, 1 * SIZE(CO1)

	ADD		C12, B1, C12
	LD		B4, 3 * SIZE(CO1)

	ADD		C22, B2, C22
	LD		B6, 5 * SIZE(CO1)

	ADD		C32, B3, C32
	LD		B8, 7 * SIZE(CO1)

	ADD		C42, B4, C42
	MADD	B1, B1, C11, A1		#	A1 = alpha_r

	SUB		C14, B5, C14
	MADD	B3, B3, C21, A1

	SUB		C24, B6, C24
	MADD	B5, B5, C31, A1

	SUB		C34, B7, C34
	MADD	B7, B7, C41, A1

	SUB		C44, B8, C44
	MADD	B2, B2, C13, A1

	MADD	B4, B4, C23, A1
	MADD	B6, B6, C33, A1

	MADD	B8, B8, C43, A1
	NMSUB	B1, B1, C13, A2		#	A2 = alpha_i

	NMSUB	B3, B3, C23, A2
	NMSUB	B5, B5, C33, A2
	LD		C13, 0 * SIZE(CO2)

	NMSUB	B7, B7, C43, A2
	MADD	B2, B2, C11, A2
	LD		C23, 2 * SIZE(CO2)

	MADD	B4, B4, C12, A2
	MADD	B6, B6, C13, A2
	LD		C33, 4 * SIZE(CO2)

	MADD	B8, B8, C14, A2
	LD		C43, 6 * SIZE(CO2)

	LD		C11, 1 * SIZE(CO2)
	LD		C21, 3 * SIZE(CO2)
	LD		C31, 5 * SIZE(CO2)
	MADD	C13, C13, C12, A1

	LD		C41, 7 * SIZE(CO2)
	MADD	C23, C23, C22, A1

	MADD	C33, C33, C32, A1
	ST		B1, 0 * SIZE(CO1)

	MADD	C43, C43, C42, A1
	ST		B3, 2 * SIZE(CO1)

	MADD	C11, C11, C14, A1
	ST		B5, 4 * SIZE(CO1)

	MADD	C21, C21, C24, A1
	ST		B7, 6 * SIZE(CO1)

	MADD	C31, C31, C34, A1
	ST		B2, 1 * SIZE(CO1)

	MADD	C41, C41, C44, A1
	ST		B4, 3 * SIZE(CO1)

	NMSUB	C13, C13, C14, A2
	ST		B6, 5 * SIZE(CO1)

	NMSUB	C23, C23, C24, A2
	ST		B8, 7 * SIZE(CO1)

	NMSUB	C33, C33, C34, A2
	NMSUB	C43, C43, C44, A2

	MADD	C11, C11, C12, A2
	MADD	C21, C21, C22, A2
	
	MADD	C31, C31, C32, A2
	MADD	C41, C41, C42, A2

	ST		C13, 0 * SIZE(CO2)
	ST		C23, 2 * SIZE(CO2)
	ST		C33, 4 * SIZE(CO2)
	ST		C43, 6 * SIZE(CO2)
	ST		C11, 1 * SIZE(CO2)
	ST		C21, 3 * SIZE(CO2)
	ST		C31, 5 * SIZE(CO2)
	ST		C41, 7 * SIZE(CO2)

#endif

#if	  defined(RN) || defined(RT) || defined(CN) || defined(CT)
	/*	(a - bi) * (c + di) */
	ADD		C11, A1, C11		#	ac'+'bd
	ADD		C21, A2, C21
#	LD		A1, 0 * SIZE(A)		#	load alpha_r
	LD		A1, 152($sp)		#	load alpha_r

	ADD		C31, A3, C31
#	LD		A2, 0 * SIZE(A)		#	load alpha_r
	LD		A2, 160($sp)		#	load alpha_i

	ADD		C41, A4, C41
	LD		B1, 0 * SIZE(CO1)

	SUB		C13, C13, A5		#	ad'+'cb
	LD		B3, 2 * SIZE(CO1)

	SUB		C23, C23, A6
	LD		B5, 4 * SIZE(CO1)

	SUB		C33, C33, A7
	LD		B7, 6 * SIZE(CO1)

	SUB		C43, C43, A8
	LD		B2, 1 * SIZE(CO1)

	ADD		C12, B1, C12
	LD		B4, 3 * SIZE(CO1)

	ADD		C22, B2, C22
	LD		B6, 5 * SIZE(CO1)

	ADD		C32, B3, C32
	LD		B8, 7 * SIZE(CO1)

	ADD		C42, B4, C42
	MADD	B1, B1, C11, A1		#	A1 = alpha_r

	SUB		C14, C14, B5
	MADD	B3, B3, C21, A1

	SUB		C24, C24, B6
	MADD	B5, B5, C31, A1

	SUB		C34, C34, B7
	MADD	B7, B7, C41, A1

	SUB		C44, C44, B8
	MADD	B2, B2, C13, A1

	MADD	B4, B4, C23, A1
	MADD	B6, B6, C33, A1

	MADD	B8, B8, C43, A1
	NMSUB	B1, B1, C13, A2		#	A2 = alpha_i

	NMSUB	B3, B3, C23, A2
	NMSUB	B5, B5, C33, A2
	LD		C13, 0 * SIZE(CO2)

	NMSUB	B7, B7, C43, A2
	MADD	B2, B2, C11, A2
	LD		C23, 2 * SIZE(CO2)

	MADD	B4, B4, C12, A2
	MADD	B6, B6, C13, A2
	LD		C33, 4 * SIZE(CO2)

	MADD	B8, B8, C14, A2
	LD		C43, 6 * SIZE(CO2)

	LD		C11, 1 * SIZE(CO2)
	LD		C21, 3 * SIZE(CO2)
	LD		C31, 5 * SIZE(CO2)
	MADD	C13, C13, C12, A1

	LD		C41, 7 * SIZE(CO2)
	MADD	C23, C23, C22, A1

	MADD	C33, C33, C32, A1
	ST		B1, 0 * SIZE(CO1)

	MADD	C43, C43, C42, A1
	ST		B3, 2 * SIZE(CO1)

	MADD	C11, C11, C14, A1
	ST		B5, 4 * SIZE(CO1)

	MADD	C21, C21, C24, A1
	ST		B7, 6 * SIZE(CO1)

	MADD	C31, C31, C34, A1
	ST		B2, 1 * SIZE(CO1)

	MADD	C41, C41, C44, A1
	ST		B4, 3 * SIZE(CO1)

	NMSUB	C13, C13, C14, A2
	ST		B6, 5 * SIZE(CO1)

	NMSUB	C23, C23, C24, A2
	ST		B8, 7 * SIZE(CO1)

	NMSUB	C33, C33, C34, A2
	NMSUB	C43, C43, C44, A2

	MADD	C11, C11, C12, A2
	MADD	C21, C21, C22, A2
	
	MADD	C31, C31, C32, A2
	MADD	C41, C41, C42, A2

	ST		C13, 0 * SIZE(CO2)
	ST		C23, 2 * SIZE(CO2)
	ST		C33, 4 * SIZE(CO2)
	ST		C43, 6 * SIZE(CO2)
	ST		C11, 1 * SIZE(CO2)
	ST		C21, 3 * SIZE(CO2)
	ST		C31, 5 * SIZE(CO2)
	ST		C41, 7 * SIZE(CO2)

#endif

#if   defined(RR) || defined(RC) || defined(CR) || defined(CC)
	/*	(a - bi) * (c - di) */
	SUB		C11, A1, C11		#	ac'+'bd
	SUB		C21, A2, C21
	LD		A1, 152($sp)		#	load alpha_r
#	LD		A1, 0 * SIZE(A)		#	load alpha_r

	SUB		C31, A3, C31
#	LD		A2, 0 * SIZE(A)		#	load alpha_i
	LD		A2, 160($sp)

	SUB		C41, A4, C41
	LD		B1, 0 * SIZE(CO1)

	ADD		C13, A5, C13		#	ad'+'cb
	LD		B3, 2 * SIZE(CO1)

	ADD		C23, A6, C23
	LD		B5, 4 * SIZE(CO1)

	ADD		C33, A7, C33
	LD		B7, 6 * SIZE(CO1)

	ADD		C43, A8, C43
	LD		B2, 1 * SIZE(CO1)

	SUB		C12, B1, C12
	LD		B4, 3 * SIZE(CO1)

	SUB		C22, B2, C22
	LD		B6, 5 * SIZE(CO1)

	SUB		C32, B3, C32
	LD		B8, 7 * SIZE(CO1)

	SUB		C42, B4, C42
	MADD	B1, B1, C11, A1		#	A1 = alpha_r

	ADD		C14, B5, C14
	MADD	B3, B3, C21, A1

	ADD		C24, B6, C24
	MADD	B5, B5, C31, A1

	ADD		C34, B7, C34
	MADD	B7, B7, C41, A1

	ADD		C44, B8, C44
	NMSUB	B2, B2, C13, A1

	NMSUB	B4, B4, C23, A1
	NMSUB	B6, B6, C33, A1

	NMSUB	B8, B8, C43, A1
	NMSUB	B1, B1, C13, A2		#	A2 = alpha_i

	NMSUB	B3, B3, C23, A2
	NMSUB	B5, B5, C33, A2
	LD		C13, 0 * SIZE(CO2)

	NMSUB	B7, B7, C43, A2
	MADD	B2, B2, C11, A2
	LD		C23, 2 * SIZE(CO2)

	MADD	B4, B4, C12, A2
	MADD	B6, B6, C13, A2
	LD		C33, 4 * SIZE(CO2)

	MADD	B8, B8, C14, A2
	LD		C43, 6 * SIZE(CO2)

	LD		C11, 1 * SIZE(CO2)
	LD		C21, 3 * SIZE(CO2)
	LD		C31, 5 * SIZE(CO2)
	MADD	C13, C13, C12, A1

	LD		C41, 7 * SIZE(CO2)
	MADD	C23, C23, C22, A1

	MADD	C33, C33, C32, A1
	ST		B1, 0 * SIZE(CO1)

	MADD	C43, C43, C42, A1
	ST		B3, 2 * SIZE(CO1)

	NMSUB	C11, C11, C14, A1
	ST		B5, 4 * SIZE(CO1)

	NMSUB	C21, C21, C24, A1
	ST		B7, 6 * SIZE(CO1)

	NMSUB	C31, C31, C34, A1
	ST		B2, 1 * SIZE(CO1)

	NMSUB	C41, C41, C44, A1
	ST		B4, 3 * SIZE(CO1)

	NMSUB	C13, C13, C14, A2
	ST		B6, 5 * SIZE(CO1)

	NMSUB	C23, C23, C24, A2
	ST		B8, 7 * SIZE(CO1)

	NMSUB	C33, C33, C34, A2
	NMSUB	C43, C43, C44, A2

	MADD	C11, C11, C12, A2
	MADD	C21, C21, C22, A2
	
	MADD	C31, C31, C32, A2
	MADD	C41, C41, C42, A2

	ST		C13, 0 * SIZE(CO2)
	ST		C23, 2 * SIZE(CO2)
	ST		C33, 4 * SIZE(CO2)
	ST		C43, 6 * SIZE(CO2)
	ST		C11, 1 * SIZE(CO2)
	ST		C21, 3 * SIZE(CO2)
	ST		C31, 5 * SIZE(CO2)
	ST		C41, 7 * SIZE(CO2)

#endif

	daddiu	CO1, CO1, 8 * SIZE
	bgtz	I, .L241
	daddiu	CO2, CO2, 8 * SIZE

	.align	4
.L22:
	andi	I, M, 2				#	MR=4
	blez	I, .L21
	NOP

	.align	4
.L21:
	andi	I, M, 1
	blez	I, .L20
	NOP

	.align	4
.L20:
	daddiu	J, J, -1
	move	B, BO
	bgtz	J, .L24
	NOP


	.align	4
.L1:
	andi	J, N, 1
	blez	J, .L999
	NOP

	.align	4
.L10:
	move	B, BO

.L999:
	ld	$16,   0($sp)
	ld	$17,   8($sp)
	ld	$18,  16($sp)
	ld	$19,  24($sp)
	ld	$20,  32($sp)
	ld	$21,  40($sp)
	ld	$22,  48($sp)

	LD	$f24, 56($sp)
	LD	$f25, 64($sp)
	LD	$f26, 72($sp)
	LD	$f27, 80($sp)
	LD	$f28, 88($sp)

#if defined(TRMMKERNEL)
	ld	$23,  96($sp)
	ld	$24, 104($sp)
	ld	$25, 112($sp)
#endif

#ifndef __64BIT__
	LD	$f20,120($sp)
	LD	$f21,128($sp)
	LD	$f22,136($sp)
	LD	$f23,144($sp)
#endif

	daddiu	$sp,$sp,STACKSIZE
	j	$31
	nop

	EPILOGUE
