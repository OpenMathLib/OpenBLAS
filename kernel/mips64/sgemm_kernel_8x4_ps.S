#define REALNAME ASMNAME
#define ASSEMBLER
#include "common.h"

#define gsLQC1(base,fq,ft,offset) .word(0x32<<26|base<<21|ft<<16|0x1<<15|offset<<6|0x1<<5|fq)
#define gsSQC1(base,fq,ft,offset) .word(0x3A<<26|base<<21|ft<<16|0x1<<15|offset<<6|0x1<<5|fq)


#define FETCH	ld
#define	STACKSIZE	192

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
#define CO3	$16
#define CO4	$17

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

	#.text
#.align	2
#	.globl	REALNAME
#	.set	nomips16
#	.ent	REALNAME
#	.type	REALNAME, @function
#REALNAME:
#	.frame	$fp,STACKSIZE,$31		# vars= 48, regs= 1/0, args= 0, gp= 0
#	.mask	0x40000000,-8
#	.fmask	0x00000000,0
#	.set	noreorder
#	.set	nomacro
	

	PROLOGUE

	daddiu	$sp,$sp,-STACKSIZE
	sd	$fp,184($sp)
	move	$fp,$sp

	sd	$16,   0($fp)
	sd	$17,   8($fp)
	sd	$18,  16($fp)
	sd	$19,  24($fp)
	sd	$20,  32($fp)
	sd	$21,  40($fp)
	sd	$22,  48($fp)

	ST	$f24, 56($fp)
	ST	$f25, 64($fp)
	ST	$f26, 72($fp)
	ST	$f27, 80($fp)
	ST	$f28, 88($fp)

#if defined(TRMMKERNEL)
	sd	$23,  96($fp)
	sd	$24, 104($fp)
	sd	$25, 112($fp)
#endif

#ifndef __64BIT__
	ST	$f20,120($fp)
	ST	$f21,128($fp)
	ST	$f22,136($fp)
	ST	$f23,144($fp)
#endif

	.align	4
.L4:
	dsra	J, N, 2				#	NR=4
	dsll	LDC, LDC, BASE_SHIFT#	LDC*SIZE

	ST		ALPHA, 152($fp)		#	Store alpha
	blez	J, .L2
	NOP


.L48:
	dsra	I, M, 3				#	MR=8
	dsll	PREA, K, BASE_SHIFT

	move	AO, A				#	Reset A
	move	CO1, C

	daddu	CO2, C,   LDC
	daddu	CO3, CO2, LDC

	daddu	CO4, CO3, LDC
	daddu	PREA, A, PREA

	blez	I, .L44
	daddu	C,   CO4, LDC

	.align	4
.L488:
	move	BO, B				#	Reset	B
	dsra	L, K, 2				#	UnRoll	K=8

	MTC		$0, C11				#	CLEAR REAULTS REGISTERS
	MOV		C12, C11		

	dsll	PREB, K, BASE_SHIFT
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
	
	daddu	PREB, B, PREB 
	MOV		C43, C11
	FETCH	$0, 0 * SIZE(CO3)

	MOV		C44, C11
	FETCH	$0, 4 * SIZE(CO3)

	PLU		B3,	B1, B1
	FETCH	$0, 0 * SIZE(CO4)

	PLU		B4, B2, B2
	blez	L, .L484
	FETCH	$0, 0 * SIZE(CO4)

.L4880:
	daddiu	L, L, -1
	MADPS	C11, C11, A1, B1
	MADPS	C21, C21, A2, B1
	gsLQC1(R13, F13, F12, 1)	#	B3 B4

	MADPS	C12, C12, A1, B2
	MADPS	C22, C22, A2, B2
	gsLQC1(R12, F5, F4, 2)		#	A5 A6

	MADPS	C31, C31, A3, B1
	MADPS	C41, C41, A4, B1
	gsLQC1(R12, F7, F6, 3)		#	A7 A8

	MADPS	C32, C32, A3, B2
	MADPS	C42, C42, A4, B2

	MADPS	C13, C13, A1, B3
	MADPS	C23, C23, A2, B3

	FETCH	$0, 0 * SIZE(PREA)
	MADPS	C33, C33, A3, B3
	MADPS	C43, C43, A4, B3

	FETCH	$0, 0 * SIZE(PREB)
	MADPS	C14, C14, A1, B4
	PLU		B7,	B5, B5

	FETCH	$0, 4 * SIZE(PREA)
	MADPS	C24, C24, A2, B4
	PLU		B8, B6, B6

	MADPS	C34, C34, A3, B4
	MADPS	C44, C44, A4, B4

	MADPS	C11, C11, A5, B5
	MADPS	C21, C21, A6, B5
	gsLQC1(R13, F9, F8, 2)		#	B1 B2

	MADPS	C12, C12, A5, B6
	MADPS	C22, C22, A6, B6
	gsLQC1(R12, F1, F0, 4)		#	A1 A2

	MADPS	C31, C31, A7, B5
	MADPS	C41, C41, A8, B5
	gsLQC1(R12, F3, F2, 5)		#	A3 A4

	MADPS	C32, C32, A7, B6
	MADPS	C42, C42, A8, B6

	FETCH	$0, 4 * SIZE(PREB)
	MADPS	C13, C13, A5, B7
	MADPS	C23, C23, A6, B7

	FETCH	$0, 8 * SIZE(PREA)
	MADPS	C33, C33, A7, B7
	MADPS	C43, C43, A8, B7

	FETCH	$0, 12 * SIZE(PREA)
	MADPS	C14, C14, A5, B8
	PLU		B3,	B1, B1

	MADPS	C24, C24, A6, B8
	PLU		B4, B2, B2

	MADPS	C34, C34, A7, B8
	MADPS	C44, C44, A8, B8

	MADPS	C11, C11, A1, B1
	MADPS	C21, C21, A2, B1
	gsLQC1(R13, F13, F12, 3)	#	B3 B4

	MADPS	C12, C12, A1, B2
	MADPS	C22, C22, A2, B2
	gsLQC1(R12, F5, F4, 6)		#	A5 A6

	MADPS	C31, C31, A3, B1
	MADPS	C41, C41, A4, B1
	gsLQC1(R12, F7, F6, 7)		#	A7 A8

	FETCH	$0, 16 * SIZE(PREA)
	MADPS	C32, C32, A3, B2
	MADPS	C42, C42, A4, B2

	MADPS	C13, C13, A1, B3
	daddiu	BO, BO, 16 * SIZE	#	4KR*4NR
	MADPS	C23, C23, A2, B3
	daddiu	AO, AO, 32 * SIZE 	#	4KR*8MR

	FETCH	$0, 20 * SIZE(PREA)
	MADPS	C33, C33, A3, B3
	MADPS	C43, C43, A4, B3

	FETCH	$0, 8 * SIZE(PREB)
	MADPS	C14, C14, A1, B4
	PLU		B7,	B5, B5

	MADPS	C24, C24, A2, B4
	PLU		B8, B6, B6

	MADPS	C34, C34, A3, B4
	MADPS	C44, C44, A4, B4

	MADPS	C11, C11, A5, B5
	MADPS	C21, C21, A6, B5
	gsLQC1(R13, F9, F8, 0)		#	B1 B2

	MADPS	C12, C12, A5, B6
	MADPS	C22, C22, A6, B6
	gsLQC1(R12, F1, F0, 0)		#	A1 A2

	MADPS	C31, C31, A7, B5
	MADPS	C41, C41, A8, B5
	gsLQC1(R12, F3, F2, 1)		#	A3 A4

	MADPS	C32, C32, A7, B6
	MADPS	C42, C42, A8, B6

	FETCH	$0, 12 * SIZE(PREB)
	MADPS	C13, C13, A5, B7
	MADPS	C23, C23, A6, B7

	FETCH	$0, 24 * SIZE(PREA)
	MADPS	C33, C33, A7, B7
	MADPS	C43, C43, A8, B7

	FETCH	$0, 28 * SIZE(PREA)
	MADPS	C14, C14, A5, B8
	PLU		B3,	B1, B1
	daddiu	PREB, PREB, 16 * SIZE	

	MADPS	C24, C24, A6, B8
	PLU		B4, B2, B2
	daddiu	PREA, PREA, 32 * SIZE
	
	MADPS	C34, C34, A7, B8
	bgtz	L, .L4880
	MADPS	C44, C44, A8, B8

	.align	4
.L484:
	andi	L, K, 4
	blez	L, .L482
	NOP

	.align	4
.L482:
	andi	L, K, 2
	blez	L, .L481
	NOP

	.align	4
.L481:				
	andi	L, K, 1
	blez	L, .L480
	NOP

	.align	4
.L480:							#	Write Back
	daddiu	I, I, -1
	CVTU	A1, C13				#	A1=C13.upper=c12
	CVTU	A2, C11				#	A2=C11.upper=c22

	CVTU	A3, C23				#	A3=C23.upper=c14
	LD		B1, 1 * SIZE(CO1)

	CVTU	A4, C21				#	A4=C21.upper=c24
	LD		B2, 1 * SIZE(CO2)

	CVTU	A5, C33				#	A5=C33.upper=c16
	LD		B3, 3 * SIZE(CO1)

	CVTU	A6, C31				#	A6=C31.upper=c26
	LD		B4, 3 * SIZE(CO2)

	CVTU	A7, C43				#	A7=C43.upper=c18
	LD		B5, 5 * SIZE(CO1)

	CVTU	A8, C41				#	A8=C41.upper=c28
	LD		B6, 5 * SIZE(CO2)

	MADD	A1, B1, A1, ALPHA		#	c12		
	LD		B7, 7 * SIZE(CO1)

	MADD	A2, B2, A2, ALPHA		#	c22
	LD		B1, 7 * SIZE(CO2)

	MADD	A3, B3, A3, ALPHA		#	c14
	LD		B2, 0 * SIZE(CO1)

	MADD	A4, B4, A4, ALPHA		#	c24
	LD		B3, 0 * SIZE(CO2)
	
	MADD	A5, B5, A5, ALPHA		#	c16
	LD		B4, 2 * SIZE(CO1)

	MADD	A6, B6, A6, ALPHA		#	c26
	LD		B5, 2 * SIZE(CO2)

	MADD	A7, B7, A7, ALPHA		#	c18
	LD		B6, 4 * SIZE(CO1)
	ST		A1, 1 * SIZE(CO1)

	MADD	A8, B1, A8, ALPHA		#	c28
	LD		B7, 4 * SIZE(CO2)
	ST		A2, 1 * SIZE(CO2)

	MADD	C11, B2, C11, ALPHA		#	c12
	LD		A1, 6 * SIZE(CO1)
	ST		A3, 3 * SIZE(CO1)

	MADD	C13, B3, C13, ALPHA		#	c22
	LD		A2, 6 * SIZE(CO2)
	ST		A4, 3 * SIZE(CO2)

	MADD	C21, B4, C21, ALPHA		#	c14
	ST		A5, 5 * SIZE(CO1)

	MADD	C23, B5, C23, ALPHA		#	c24
	ST		A6, 5 * SIZE(CO2)

	MADD	C31, B6, C31, ALPHA		#	c16
	ST		A7, 7 * SIZE(CO1)

	MADD	C33, B7, C33, ALPHA		#	c26
	ST		A8, 7 * SIZE(CO2)

	MADD	C41, A1, C41, ALPHA		#	c18
	ST		C11, 0 * SIZE(CO1)

	MADD	C43, A2, C43, ALPHA		#	c28
	ST		C13, 0 * SIZE(CO2)
	
	ST		C21, 2 * SIZE(CO1)
	ST		C23, 2 * SIZE(CO2)
	ST		C31, 4 * SIZE(CO1)
	ST		C33, 4 * SIZE(CO2)
	ST		C41, 6 * SIZE(CO1)

	CVTU	A1, C14				#	B1=C12.upper=c42
	ST		C43, 6 * SIZE(CO2)

	CVTU	A2, C12				#	B2=C14.upper=c32
	LD		B1, 1 * SIZE(CO3)

	CVTU	A3, C24				#	B3=C22.upper=c44
	LD		B2, 1 * SIZE(CO4)

	CVTU	A4, C22				#	B4=C24.upper=c34
	LD		B3, 3 * SIZE(CO3)

	CVTU	A5, C34				#	B5=C32.upper=c46
	LD		B4, 3 * SIZE(CO4)

	CVTU	A6, C32				#	B6=C24.upper=c36
	LD		B5, 5 * SIZE(CO3)

	CVTU	A7, C44				#	B7=C42.upper=c48
	LD		B6, 5 * SIZE(CO4)

	CVTU	A8, C42				#	A1=C44.upper=c38	
	LD		B7, 7 * SIZE(CO3)

	MADD	A1, B1, A1, ALPHA		#	c31
	LD		C11, 7 * SIZE(CO4)

	MADD	A2, B2, A2, ALPHA
	LD		C13, 0 * SIZE(CO3)
	
	MADD	A3, B3, A3, ALPHA
	LD		C21, 0 * SIZE(CO4)
	
	MADD	A4, B4, A4, ALPHA
	LD		C23, 2 * SIZE(CO3)

	MADD	A5, B5, A5, ALPHA
	LD		C31, 2 * SIZE(CO4)
	
	MADD	A6, B6, A6, ALPHA
	LD		C33, 4 * SIZE(CO3)
	
	MADD	A7, B7, A7, ALPHA
	LD		C41, 4 * SIZE(CO4)
	ST		A1, 1 * SIZE(CO3)

	MADD	A8, C11, A8, ALPHA
	LD		C43, 6 * SIZE(CO3)
	ST		A2, 1 * SIZE(CO4)

	MADD	C12, C13, C12, ALPHA
	LD		B1, 6 * SIZE(CO4)
	ST		A3, 3 * SIZE(CO3)

	MADD	C14, C21, C14, ALPHA
	ST		A4, 3 * SIZE(CO4)

	MADD	C22, C23, C22, ALPHA
	ST		A5, 5 * SIZE(CO3)

	MADD	C24, C31, C24, ALPHA
	ST		A6, 5 * SIZE(CO4)

	MADD	C32, C33, C32, ALPHA
	ST		A7, 7 * SIZE(CO3)

	MADD	C34, C41, C34, ALPHA
	ST		A8, 7 * SIZE(CO4)

	MADD	C42, C43, C42, ALPHA
	ST		C12, 0 * SIZE(CO3)

	MADD	C44, B1, C44, ALPHA
	ST		C14, 0 * SIZE(CO4)

	ST		C22, 2 * SIZE(CO3)
	daddiu	CO1, CO1, 8 * SIZE

	ST		C24, 2 * SIZE(CO4)
	daddiu	CO2, CO2, 8 * SIZE
	
	ST		C32, 4 * SIZE(CO3)
	ST		C34, 4 * SIZE(CO4)
	ST		C42, 6 * SIZE(CO3)
	ST		C44, 6 * SIZE(CO4)

	daddiu	CO3, CO3, 8 * SIZE
	bgtz	I, .L488
	daddiu	CO4, CO4, 8 * SIZE

.L44:

.L40:
	daddiu	J, J, -1
	move	B, BO

	bgtz	J, .L48
	NOP

	.align	4
.L2:							#	Nr=2
	andi	J, N, 2			
	blez	J, .L1
	NOP



	.align	4
.L1:
	andi	J, N, 1
	blez	J, .L999
	NOP



.L999:
	ld	$16,   0($fp)
	ld	$17,   8($fp)
	ld	$18,  16($fp)
	ld	$19,  24($fp)
	ld	$20,  32($fp)
	ld	$21,  40($fp)
	ld	$22,  48($fp)

	LD	$f24, 56($fp)
	LD	$f25, 64($fp)
	LD	$f26, 72($fp)
	LD	$f27, 80($fp)
	LD	$f28, 88($fp)

#if defined(TRMMKERNEL)
	ld	$23,  96($fp)
	ld	$24, 104($fp)
	ld	$25, 112($fp)
#endif

#ifndef __64BIT__
	LD	$f20,120($fp)
	LD	$f21,128($fp)
	LD	$f22,136($fp)
	LD	$f23,144($fp)
#endif

	move	$sp,$fp
	ld	$fp,184($sp)
	daddiu	$sp,$sp,STACKSIZE
	j	$31
	nop

	EPILOGUE
#	.set	macro
#	.set	reorder
#	.end	REALNAME
#	.size	REALNAME, .-REALNAME
#.ident	"GCC: (Debian 4.4.6-6) 4.4.6"
