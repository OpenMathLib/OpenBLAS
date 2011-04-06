#define gsLQC1(base,fq,ft,offset) .word(0x32<<26|base<<21|ft<<16|0x1<<15|offset<<6|0x1<<5|fq)
#define gsSQC1(base,fq,ft,offset) .word(0x3A<<26|base<<21|ft<<16|0x1<<15|offset<<6|0x1<<5|fq)
#define REALNAME ASMNAME 
#define PROLOGUE \
	.text ;\
	.set	mips64 ;\
	.align 5 ;\
	.globl	REALNAME ;\
	.ent	REALNAME ;\
	.type	REALNAME, @function ;\
REALNAME: ;\
	.set	noreorder ;\
	.set	nomacro

#define EPILOGUE \
	.set	macro ;\
	.set	reorder ;\
	.end	REALNAME
#define BASE_SHIFT 3
#define FETCH	ld

#define M	$4
#define	N	$5
#define	K	$6
#define A	$8
#define B	$9
#define C	$10
#define LDC	$11

#define AO	$12
#define BO	$13

#define I	$2
#define J	$3
#define L	$7

#define CO1	$14
#define CO2	$15
#define CO3	$16
#define CO4	$17

#define KCO	$18
#define MCO	$19
#define NCO	$20

#define SPANB	$21
#define SPANC	$22
#define PREB	$23
#define PREA	$24
#define SPANA	$25

#define ALPHA	$f15

#define R8	8
#define	R9	9
#define R14	14
#define R15	15
#define R16	16
#define R17 17

#define	t11	$f30
#define	t21	$f31
#define	t31	$f28
#define	t41	$f29

#define	t12	$f26
#define	t22	$f27
#define	t32	$f24
#define	t42	$f25

#define	t13	$f22
#define	t23	$f23
#define	t33	$f20
#define	t43	$f21

#define	t14	$f18
#define	t24	$f19
#define	t34	$f16
#define	t44	$f17

#define	c11	$f0
#define	c21	$f1
#define	c31	$f2
#define	c41	$f3

#define	c12	$f4
#define	c22	$f5
#define	c32	$f6
#define	c42	$f7

#define	c13	$f8
#define	c23	$f9
#define	c33	$f10
#define c43	$f11

#define	c14	$f12
#define	c24	$f13
#define	c34	$f14
#define	c44	$f0

#define	a0	$f0
#define	a1	$f1
#define	a2	$f2
#define	a3	$f3
#define	a4	$f4
#define	a5	$f5
#define	a6	$f6
#define	a7	$f7
#define	b0	$f8
#define	b1	$f9
#define	b2	$f10
#define b3	$f11
#define	b4	$f12
#define	b5	$f13
#define	b6	$f14
#define	b7	$f15

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

	PROLOGUE
	
	daddiu	$sp, $sp, -160
	sd	$16,   0($sp)
	sd	$17,   8($sp)
	sd	$18,  16($sp)
	sd	$19,  24($sp)
	sd	$20,  32($sp)
	sd	$21,  40($sp)
	sd	$22,  48($sp)
	sdc1	$f24, 56($sp)
	sdc1	$f25, 64($sp)
	sdc1	$f26, 72($sp)
	sdc1	$f27, 80($sp)
	sdc1	$f28, 88($sp)
	sd	$23,  96($sp)
	sd	$24, 104($sp)
	sd	$25, 112($sp)
	sdc1	$f20,120($sp)
	sdc1	$f21,128($sp)
	sdc1	$f22,136($sp)
	sdc1	$f23,144($sp)


	.align	5					#  BACKUP
.L0_N4:							#  Loop N
	sdc1	ALPHA,152($sp)		#  Backup	ALPHA
	move	MCO,M				#  Backup	M

	move	NCO,N				#  Backup	N
	move	KCO,K				#  Backup	K

	move	AO,A				#  Backup	A_addr
	move	BO,B				#  Backup	B_addr
	
	dsll	LDC,LDC,3			#  LDC*8Byte
	dsll	SPANB,KCO,5			#  SPANB=KC*NR(4)*8Byte=KC*2^5
	
	dsll	SPANA,KCO,5			#  SPANA = KCO*4mr*8Byte
	dsra	N,NCO,2				#  N=NCO/2
	beq		N,$0,.L0_N2			#  N=0,NCO<4
	dsll	SPANC,LDC,2			#  SPANC=LDC*4

.L0_N4_Lb:
	move	CO1,C				#  Set C	
	dsra	M,MCO,2				#  M=MCO/2
	
	move	A,AO				#  Reset A
	daddu	CO2,CO1,LDC

	daddu	CO3,CO2,LDC
	daddu	PREB,BO,SPANB		#  PreB point next panelB

	daddu	CO4,CO3,LDC
	beqz	M,.L14_M2
	daddu	PREA,AO,SPANA

.L10:
	dmtc1	$0,t11
	mov.d	t21,t11
	gsLQC1(R8,F1,F0,0)			#a0,a1

	mov.d	t31,t11
	mov.d	t41,t11
	gsLQC1(R9,F9,F8,0)			#b0,b1

	mov.d	t12,t11
	mov.d	t22,t11
	gsLQC1(R8,F3,F2,1)			#a2,a3
	
	mov.d	t32,t11
	mov.d	t42,t11
	gsLQC1(R9,F11,F10,1)		#b2,b3

	dsra	K,KCO,2				#  K=KCO/2
	mov.d	t13,t11
	
	mov.d	t23,t11
	mov.d	t33,t11
	
	mov.d	t43,t11
	mov.d	t14,t11
	
	mov.d	t24,t11
	mov.d	t34,t11
	
	mov.d	t44,t11
	beqz	K,.L15
	nop

.L11:							#  N=M=K=4
	gsLQC1(R8,F5,F4,2)			#  R8=A
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0

	gsLQC1(R9,F13,F12,2)		#  R9=B
	madd.d	t12,t12,a0,b1
	madd.d	t22,t22,a1,b1

	gsLQC1(R8,F7,F6,3)
	madd.d	t31,t31,a2,b0
	madd.d	t41,t41,a3,b0
	
	gsLQC1(R9,F15,F14,3)
	madd.d	t32,t32,a2,b1
	madd.d	t42,t42,a3,b1

	FETCH		$0,(PREB)
	madd.d	t13,t13,a0,b2
	madd.d	t23,t23,a1,b2
	
	madd.d	t14,t14,a0,b3
	madd.d	t24,t24,a1,b3
	
	FETCH		$0,(PREA)
	madd.d	t33,t33,a2,b2
	madd.d	t43,t43,a3,b2

	madd.d	t34,t34,a2,b3
	madd.d	t44,t44,a3,b3
								#load2 comp1
.L12:
	gsLQC1(R8,F1,F0,4)
	madd.d	t11,t11,a4,b4
	madd.d	t21,t21,a5,b4

	gsLQC1(R9,F9,F8,4)
	madd.d	t12,t12,a4,b5
	madd.d	t22,t22,a5,b5

	gsLQC1(R8,F3,F2,5)
	madd.d	t31,t31,a6,b4
	madd.d	t41,t41,a7,b4

	gsLQC1(R9,F11,F10,5)
	madd.d	t32,t32,a6,b5
	madd.d	t42,t42,a7,b5

	FETCH		$0,32(PREB)
	madd.d	t13,t13,a4,b6
	madd.d	t23,t23,a5,b6

	madd.d	t14,t14,a4,b7
	madd.d	t24,t24,a5,b7

	FETCH		$0,32(PREA)
	madd.d	t33,t33,a6,b6
	madd.d	t43,t43,a7,b6
	
	madd.d	t34,t34,a6,b7
	madd.d	t44,t44,a7,b7

.L13:
	gsLQC1(R8,F5,F4,6)	
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0

	gsLQC1(R9,F13,F12,6)
	madd.d	t12,t12,a0,b1
	madd.d	t22,t22,a1,b1

	gsLQC1(R8,F7,F6,7)
	madd.d	t31,t31,a2,b0
	madd.d	t41,t41,a3,b0

	gsLQC1(R9,F15,F14,7)
	madd.d	t32,t32,a2,b1
	madd.d	t42,t42,a3,b1
	daddu	A,A,128				#  A+=4(mr)*4(kr)*8Byte=128

	FETCH		$0,64(PREB)
	madd.d	t13,t13,a0,b2
	madd.d	t23,t23,a1,b2
	daddu	B,B,128

	madd.d	t14,t14,a0,b3
	madd.d	t24,t24,a1,b3

	FETCH		$0,64(PREA)
	madd.d	t33,t33,a2,b2
	madd.d	t43,t43,a3,b2

	madd.d	t34,t34,a2,b3
	madd.d	t44,t44,a3,b3
	
.L14:
	gsLQC1(R8,F1,F0,0)
	madd.d	t11,t11,a4,b4
	madd.d	t21,t21,a5,b4

	gsLQC1(R9,F9,F8,0)
	madd.d	t12,t12,a4,b5
	madd.d	t22,t22,a5,b5

	gsLQC1(R8,F3,F2,1)
	madd.d	t31,t31,a6,b4
	madd.d	t41,t41,a7,b4
	daddiu	K,K,-1

	gsLQC1(R9,F11,F10,1)
	madd.d	t32,t32,a6,b5
	madd.d	t42,t42,a7,b5

	FETCH		$0,96(PREB)
	madd.d	t13,t13,a4,b6
	madd.d	t23,t23,a5,b6

	madd.d	t14,t14,a4,b7
	madd.d	t24,t24,a5,b7

	FETCH		$0,96(PREA)
	madd.d	t33,t33,a6,b6
	madd.d	t43,t43,a7,b6
	daddu	PREB,PREB,128
	
	madd.d	t34,t34,a6,b7
	daddu	PREA,PREA,128
	bnez 	K,.L11
	madd.d	t44,t44,a7,b7

.L15:							#  N=4 M=4 K=2
	and		K,KCO,2				#  k = KCO&2
	beqz	K,.L18
	nop

.L16:			
	gsLQC1(R8,F5,F4,2)			#  R8=A
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0

	gsLQC1(R9,F13,F12,2)		#  R9=B
	madd.d	t12,t12,a0,b1
	madd.d	t22,t22,a1,b1

	gsLQC1(R8,F7,F6,3)
	madd.d	t31,t31,a2,b0
	madd.d	t41,t41,a3,b0

	gsLQC1(R9,F15,F14,3)
	madd.d	t32,t32,a2,b1
	madd.d	t42,t42,a3,b1
	daddu	A,A,64				#  A+=4(mr)*2(kr)*8Byte=64

	FETCH		$0,0(PREB)
	madd.d	t13,t13,a0,b2
	madd.d	t23,t23,a1,b2
	daddu	B,B,64

	madd.d	t14,t14,a0,b3
	madd.d	t24,t24,a1,b3

	FETCH		$0,0(PREA)
	madd.d	t33,t33,a2,b2
	madd.d	t43,t43,a3,b2

	madd.d	t34,t34,a2,b3
	madd.d	t44,t44,a3,b3
								
.L17:
	gsLQC1(R8,F1,F0,0)
	madd.d	t11,t11,a4,b4
	madd.d	t21,t21,a5,b4

	gsLQC1(R9,F9,F8,0)
	madd.d	t12,t12,a4,b5
	madd.d	t22,t22,a5,b5

	gsLQC1(R8,F3,F2,1)
	madd.d	t31,t31,a6,b4
	madd.d	t41,t41,a7,b4

	gsLQC1(R9,F11,F10,1)
	madd.d	t32,t32,a6,b5
	madd.d	t42,t42,a7,b5

	FETCH		$0,32(PREB)
	madd.d	t13,t13,a4,b6
	madd.d	t23,t23,a5,b6

	madd.d	t14,t14,a4,b7
	madd.d	t24,t24,a5,b7

	FETCH		$0,32(PREA)
	madd.d	t33,t33,a6,b6
	madd.d	t43,t43,a7,b6
	daddu	PREB,PREB,64
	
	madd.d	t34,t34,a6,b7
	madd.d	t44,t44,a7,b7
	daddu	PREA,PREA,64
	
.L18:							#  N=4, M=4, K=1
	and		K,KCO,1
	beqz	K,.L19				#  
	ldc1	ALPHA,152($sp)		#  Get ALPHA
	
	FETCH		$0,0(PREB)
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0
	daddu	A,A,32				#  A+=4(mr)*1(kr)*8Byte=32

	madd.d	t12,t12,a0,b1
	madd.d	t22,t22,a1,b1
	daddu	B,B,32

	FETCH		$0,0(PREA)
	madd.d	t31,t31,a2,b0
	madd.d	t41,t41,a3,b0
	daddu	PREB,PREB,32

	madd.d	t32,t32,a2,b1
	madd.d	t42,t42,a3,b1
	daddu	PREA,PREA,32

	madd.d	t13,t13,a0,b2
	madd.d	t23,t23,a1,b2

	madd.d	t14,t14,a0,b3
	madd.d	t24,t24,a1,b3

	madd.d	t33,t33,a2,b2
	madd.d	t43,t43,a3,b2

	madd.d	t34,t34,a2,b3
	madd.d	t44,t44,a3,b3

.L19:							#  Write Back
	ldc1	c11,0(CO1)			#  Fetch 16 C
	ldc1	c21,8(CO1)			
	ldc1	c31,16(CO1)
	ldc1	c41,24(CO1)

	ldc1	c12,0(CO2)
	madd.d	t11,c11,t11,ALPHA
	ldc1	c22,8(CO2)
	madd.d	t21,c21,t21,ALPHA
	ldc1	c32,16(CO2)
	madd.d	t31,c31,t31,ALPHA
	ldc1	c42,24(CO2)
	madd.d	t41,c41,t41,ALPHA

	ldc1	c13,0(CO3)
	madd.d	t12,c12,t12,ALPHA
	ldc1	c23,8(CO3)
	madd.d	t22,c22,t22,ALPHA
	ldc1	c33,16(CO3)
	madd.d	t32,c32,t32,ALPHA
	ldc1	c43,24(CO3)
	madd.d	t42,c42,t42,ALPHA

	ldc1	c14,0(CO4)
	madd.d	t13,c13,t13,ALPHA
	ldc1	c24,8(CO4)
	madd.d	t23,c23,t23,ALPHA
	ldc1	c34,16(CO4)
	madd.d	t33,c33,t33,ALPHA
	ldc1	c44,24(CO4)
	madd.d	t43,c43,t43,ALPHA

	sdc1	t11,0(CO1)
	madd.d	t14,c14,t14,ALPHA
	sdc1	t21,8(CO1)
	madd.d	t24,c24,t24,ALPHA
	sdc1	t31,16(CO1)
	madd.d	t34,c34,t34,ALPHA
	sdc1	t41,24(CO1)
	madd.d	t44,c44,t44,ALPHA
	daddiu	M,M,-1				#  M--

	sdc1	t12,0(CO2)
	sdc1	t22,8(CO2)
	sdc1	t32,16(CO2)
	sdc1	t42,24(CO2)

	sdc1	t13,0(CO3)
	sdc1	t23,8(CO3)
	sdc1	t33,16(CO3)
	sdc1	t43,24(CO3)

	FETCH	$0,32(CO1)
	FETCH	$0,32(CO2)
	FETCH	$0,32(CO3)
	FETCH	$0,32(CO4)

	sdc1	t14,0(CO4)
	daddu	CO1,CO1,32			#  COx += 4*8Byte
	sdc1	t24,8(CO4)
	daddu	CO2,CO2,32
	sdc1	t34,16(CO4)
	daddu	CO3,CO3,32
	sdc1	t44,24(CO4)
	move	B,BO				#  Reset B
	daddu	PREB,BO,SPANB
	bnez	M,.L10				#  M!=0
	daddu	CO4,CO4,32



.L14_M2:
	and		M,MCO,2				#  Remainder M = 2
	beqz	M,.L14_M1			
	nop

.L20:
	dmtc1	$0,t11
	mov.d	t21,t11
	gsLQC1(R8,F1,F0,0)			#a0,a1

	mov.d	t12,t11
	mov.d	t22,t11
	gsLQC1(R9,F9,F8,0)			#b0,b1
	
	dsra	K,KCO,2				#  K=KCO/2
	mov.d	t13,t11
	gsLQC1(R9,F11,F10,1)		#b2,b3
	
	mov.d	t23,t11
	mov.d	t14,t11
	
	mov.d	t24,t11
	beqz	K,.L25
	nop

.L21:							#  N=4 m=2,=K=4
	gsLQC1(R8,F5,F4,1)			#  R8=A
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0

	gsLQC1(R9,F13,F12,2)		#  R9=B
	madd.d	t12,t12,a0,b1
	madd.d	t22,t22,a1,b1

	gsLQC1(R9,F15,F14,3)
	madd.d	t13,t13,a0,b2
	madd.d	t23,t23,a1,b2
	
	gsLQC1(R8,F3,F2,2)
	madd.d	t14,t14,a0,b3
	madd.d	t24,t24,a1,b3
	
	gsLQC1(R9,F9,F8,4)
	madd.d	t11,t11,a4,b4
	madd.d	t21,t21,a5,b4

	gsLQC1(R9,F11,F10,5)
	madd.d	t12,t12,a4,b5
	madd.d	t22,t22,a5,b5

	gsLQC1(R8,F7,F6,3)	
	madd.d	t13,t13,a4,b6
	madd.d	t23,t23,a5,b6

	madd.d	t14,t14,a4,b7
	madd.d	t24,t24,a5,b7

	gsLQC1(R9,F13,F12,6)
	madd.d	t11,t11,a2,b0
	madd.d	t21,t21,a3,b0
	daddu	A,A,64				#  A+=2(mr)*4(kr)*8Byte=64

	gsLQC1(R9,F15,F14,7)
	madd.d	t12,t12,a2,b1
	madd.d	t22,t22,a3,b1
	daddiu	K,K,-1

	gsLQC1(R8,F1,F0,0)
	madd.d	t13,t13,a2,b2
	madd.d	t23,t23,a3,b2
	daddu	B,B,128				#  B+=4(nr)*4(kr)*8Byte=128

	madd.d	t14,t14,a2,b3
	madd.d	t24,t24,a3,b3

	gsLQC1(R9,F9,F8,0)
	madd.d	t11,t11,a6,b4
	madd.d	t21,t21,a7,b4

	gsLQC1(R9,F11,F10,1)
	madd.d	t12,t12,a6,b5
	madd.d	t22,t22,a7,b5

	madd.d	t13,t13,a6,b6
	madd.d	t23,t23,a7,b6

	madd.d	t14,t14,a6,b7
	bnez 	K,.L21
	madd.d	t24,t24,a7,b7

.L25:							#  N=4 M=2 K=2
	and		K,KCO,2				#  k = KCO&2
	beqz	K,.L28
	nop

.L26:			
	gsLQC1(R8,F5,F4,1)			#  R8=A
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0

	gsLQC1(R9,F13,F12,2)		#  R9=B
	madd.d	t12,t12,a0,b1
	madd.d	t22,t22,a1,b1
	daddu	A,A,32				#  A+=2(mr)*2(kr)*8Byte=32

	gsLQC1(R9,F15,F14,3)
	madd.d	t13,t13,a0,b2
	madd.d	t23,t23,a1,b2
	daddu	B,B,64

	madd.d	t14,t14,a0,b3
	madd.d	t24,t24,a1,b3

.L27:
	gsLQC1(R8,F1,F0,0)
	madd.d	t11,t11,a4,b4
	madd.d	t21,t21,a5,b4

	gsLQC1(R9,F9,F8,0)
	madd.d	t12,t12,a4,b5
	madd.d	t22,t22,a5,b5

	gsLQC1(R9,F11,F10,1)
	madd.d	t13,t13,a4,b6
	madd.d	t23,t23,a5,b6

	madd.d	t14,t14,a4,b7
	madd.d	t24,t24,a5,b7
	
.L28:							#  N=4, M=2, K=1
	and		K,KCO,1
	beqz	K,.L29				#  
	ldc1	ALPHA,152($sp)		#  Get ALPHA
	
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0
	daddu	A,A,16				#  A+=2(mr)*1(kr)*8Byte=16
	daddu	B,B,32

	madd.d	t12,t12,a0,b1
	madd.d	t22,t22,a1,b1

	madd.d	t13,t13,a0,b2
	madd.d	t23,t23,a1,b2

	madd.d	t14,t14,a0,b3
	madd.d	t24,t24,a1,b3

.L29:							#  Write Back
	ldc1	c11,0(CO1)			#  Fetch 16 C
	ldc1	c21,8(CO1)			

	ldc1	c12,0(CO2)
	ldc1	c22,8(CO2)
	
	ldc1	c13,0(CO3)
	madd.d	t11,c11,t11,ALPHA
	ldc1	c23,8(CO3)
	madd.d	t21,c21,t21,ALPHA

	ldc1	c14,0(CO4)
	madd.d	t12,c12,t12,ALPHA
	ldc1	c24,8(CO4)
	madd.d	t22,c22,t22,ALPHA

	sdc1	t11,0(CO1)
	madd.d	t13,c13,t13,ALPHA
	sdc1	t21,8(CO1)
	madd.d	t23,c23,t23,ALPHA

	sdc1	t12,0(CO2)
	madd.d	t14,c14,t14,ALPHA
	sdc1	t22,8(CO2)
	madd.d	t24,c24,t24,ALPHA

	sdc1	t13,0(CO3)
	move	B,BO				#  Reset B
	sdc1	t23,8(CO3)
	daddu	CO1,CO1,16			#  COx += 2*8Byte

	FETCH	$0,0(CO1)
	FETCH	$0,16(CO2)
	FETCH	$0,16(CO3)
	FETCH	$0,16(CO4)

	sdc1	t14,0(CO4)
	daddu	CO2,CO2,16
	sdc1	t24,8(CO4)
	daddu	CO3,CO3,16
	daddu	CO4,CO4,16



.L14_M1:
	and		M,MCO,1				#  Remainder M = 1
	beqz	M,.L0_N4_Loop		#  M = 0, finishing one panel B
	nop

.L30:
	ldc1	a0,0(A)
	dsra	K,KCO,2				#  K=KCO/2
	gsLQC1(R9,F9,F8,0)			#b0,b1
	dmtc1	$0,t11
	gsLQC1(R9,F11,F10,1)		#b2,b3
	mov.d	t12,t11
	mov.d	t13,t11
	beqz	K,.L25
	mov.d	t14,t11

.L31:							#  N=4 m=1,=K=4
	ldc1	a1,8(A)
	
	gsLQC1(R9,F13,F12,2)		#  R9=B
	madd.d	t11,t11,a0,b0
	madd.d	t12,t12,a0,b1
	
	gsLQC1(R9,F15,F14,3)
	madd.d	t13,t13,a0,b2
	madd.d	t14,t14,a0,b3
	
	ldc1	a2,16(A)
	
	gsLQC1(R9,F9,F8,4)
	madd.d	t11,t11,a1,b4
	madd.d	t12,t12,a1,b5
	
	gsLQC1(R9,F11,F10,5)
	madd.d	t13,t13,a1,b6
	madd.d	t14,t14,a1,b7

	ldc1	a3,24(A)
	daddiu	K,K,-1
	
	gsLQC1(R9,F13,F12,6)
	madd.d	t11,t11,a2,b0
	madd.d	t12,t12,a2,b1
	daddu	A,A,32				#  A+=1(mr)*4(kr)*8Byte=64
	
	gsLQC1(R9,F15,F14,7)
	madd.d	t13,t13,a2,b2
	madd.d	t14,t14,a2,b3
	daddu	B,B,128				#  B+=4(nr)*4(kr)*8Byte=128

	ldc1	a0,0(A)
	
	gsLQC1(R9,F9,F8,0)
	madd.d	t11,t11,a3,b4
	madd.d	t12,t12,a3,b5
	
	gsLQC1(R9,F11,F10,1)
	madd.d	t13,t13,a3,b6
	bnez 	K,.L31
	madd.d	t14,t14,a3,b7

.L35:							#  N=4 M=1 K=2
	and		K,KCO,2				#  k = KCO&2
	beqz	K,.L38
	nop

.L36:			
	ldc1	a1,8(A)	
	
	gsLQC1(R9,F13,F12,2)		#  R9=B
	madd.d	t11,t11,a0,b0
	madd.d	t12,t12,a0,b1
	daddu	A,A,16				#  A+=1(mr)*2(kr)*8Byte=32
	
	gsLQC1(R9,F15,F14,3)
	madd.d	t13,t13,a0,b2
	madd.d	t14,t14,a0,b3
	daddu	B,B,64


.L37:
	ldc1	a0,0(A)
	
	gsLQC1(R9,F9,F8,0)
	madd.d	t11,t11,a1,b4
	madd.d	t12,t12,a1,b5
	
	gsLQC1(R9,F11,F10,1)
	madd.d	t13,t13,a1,b6
	madd.d	t14,t14,a1,b7
	
.L38:							#  N=4, M=1, K=1
	and		K,KCO,1
	beqz	K,.L39				#  
	ldc1	ALPHA,152($sp)		#  Get ALPHA
	
	madd.d	t11,t11,a0,b0
	madd.d	t12,t12,a0,b1
	daddu	A,A,8				#  A+=1(mr)*1(kr)*8Byte=16
	daddu	B,B,32
	
	madd.d	t13,t13,a0,b2
	madd.d	t14,t14,a0,b3

.L39:							#  Write Back
	ldc1	c11,0(CO1)			#  Fetch 16 C
	ldc1	c12,0(CO2)
	ldc1	c13,0(CO3)
	ldc1	c14,0(CO4)
	
	madd.d	t11,c11,t11,ALPHA
	madd.d	t12,c12,t12,ALPHA
	madd.d	t13,c13,t13,ALPHA
	madd.d	t14,c14,t14,ALPHA

	sdc1	t11,0(CO1)
	sdc1	t12,0(CO2)
	sdc1	t13,0(CO3)
	sdc1	t14,0(CO4)


.L0_N4_Loop:
	daddu	BO,BO,SPANB			#  BO point to next panel B
	daddiu	N,N,-1				#  N--
	daddu	C,C,SPANC			#  C pointe to next panel C
	bnez	N,.L0_N4_Lb			#  N!=0
	move	B,BO				#  Set B



	.align	5					
.L0_N2:
	and		N,NCO,2				#  Remainder N = 2
	beqz	N,.L0_N1			#  N=0,NCO<2
	dsll	SPANC,LDC,1			#  SPANC=LDC*2

.L0_N2_Lb:
	move	CO1,C				#  Set C	
	dsra	M,MCO,2				#  M=MCO/2
	
	dsll	SPANB,KCO,4			#  SPANB=KC*NR(2)*8Byte=KC*16=KC*2^4
	move	A,AO				#  Reset A

	daddu	CO2,CO1,LDC
	beqz	M,.L12_M2
	daddu	PREA,AO,SPANA

.L40:						
	dmtc1	$0,t11
	mov.d	t21,t11
	gsLQC1(R8,F1,F0,0)			#a0,a1

	mov.d	t31,t11
	mov.d	t41,t11
	gsLQC1(R9,F9,F8,0)			#b0,b1

	dsra	K,KCO,2				#  K=KCO/2
	mov.d	t12,t11
	gsLQC1(R8,F3,F2,1)			#a2,a3
	
	mov.d	t22,t11
	mov.d	t32,t11
	
	mov.d	t42,t11
	beqz	K,.L45
	nop

.L41:							#  N=2,M=K=4
	gsLQC1(R8,F5,F4,2)			#  R8=A
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0

	gsLQC1(R9,F13,F12,1)		#  R9=B
	madd.d	t12,t12,a0,b1
	madd.d	t22,t22,a1,b1

	gsLQC1(R8,F7,F6,3)
	madd.d	t31,t31,a2,b0
	madd.d	t41,t41,a3,b0
	
	FETCH		$0,(PREA)
	madd.d	t32,t32,a2,b1
	madd.d	t42,t42,a3,b1

.L42:
	gsLQC1(R8,F1,F0,4)
	madd.d	t11,t11,a4,b4
	madd.d	t21,t21,a5,b4

	gsLQC1(R9,F11,F10,2)
	madd.d	t12,t12,a4,b5
	madd.d	t22,t22,a5,b5

	gsLQC1(R8,F3,F2,5)
	madd.d	t31,t31,a6,b4
	madd.d	t41,t41,a7,b4

	FETCH		$0,32(PREA)
	madd.d	t32,t32,a6,b5
	madd.d	t42,t42,a7,b5

.L43:
	gsLQC1(R8,F5,F4,6)	
	madd.d	t11,t11,a0,b2
	madd.d	t21,t21,a1,b2

	gsLQC1(R9,F15,F14,3)
	madd.d	t12,t12,a0,b3
	madd.d	t22,t22,a1,b3

	gsLQC1(R8,F7,F6,7)
	madd.d	t31,t31,a2,b2
	madd.d	t41,t41,a3,b2
	daddu	B,B,64				#  B+=2(nr)*4(kr)*8Byte=64

	FETCH		$0,64(PREA)
	madd.d	t32,t32,a2,b3
	madd.d	t42,t42,a3,b3
	daddu	A,A,128				#  A+=4(mr)*4(kr)*8Byte=128

.L44:
	gsLQC1(R8,F1,F0,0)
	madd.d	t11,t11,a4,b6
	madd.d	t21,t21,a5,b6
	daddiu	K,K,-1

	gsLQC1(R9,F9,F8,0)
	madd.d	t12,t12,a4,b7
	madd.d	t22,t22,a5,b7
	daddu	PREA,PREA,128

	gsLQC1(R8,F3,F2,1)
	madd.d	t31,t31,a6,b6
	madd.d	t41,t41,a7,b6

	FETCH		$0,-32(PREA)
	madd.d	t32,t32,a6,b7
	bnez 	K,.L41
	madd.d	t42,t42,a7,b7


.L45:							#  N=2 M=4 K=2
	and		K,KCO,2				#  k = KCO&2
	beqz	K,.L48
	nop

.L46:			
	gsLQC1(R8,F5,F4,2)			#  R8=A
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0

	gsLQC1(R9,F13,F12,1)		#  R9=B
	madd.d	t12,t12,a0,b1
	madd.d	t22,t22,a1,b1

	gsLQC1(R8,F7,F6,3)
	madd.d	t31,t31,a2,b0
	madd.d	t41,t41,a3,b0
	daddu	B,B,32				#  B+=2(nr)*2(kr)*8Byte=32

	FETCH		$0,0(PREA)
	madd.d	t32,t32,a2,b1
	madd.d	t42,t42,a3,b1
	daddu	A,A,64				#  A+=4(mr)*2(kr)*8Byte=64

.L47:
	gsLQC1(R8,F1,F0,0)
	madd.d	t11,t11,a4,b4
	madd.d	t21,t21,a5,b4

	gsLQC1(R9,F9,F8,0)
	madd.d	t12,t12,a4,b5
	madd.d	t22,t22,a5,b5

	gsLQC1(R8,F3,F2,1)
	madd.d	t31,t31,a6,b4
	madd.d	t41,t41,a7,b4

	FETCH		$0,32(PREA)
	madd.d	t32,t32,a6,b5
	madd.d	t42,t42,a7,b5
	daddu	PREA,PREA,64

	
.L48:							#  N=2, M=4, K=1
	and		K,KCO,1
	beqz	K,.L49				#  
	ldc1	ALPHA,152($sp)		#  Get ALPHA
	
	FETCH		$0,0(PREA)
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0
	daddu	A,A,32				#  A+=4(mr)*1(kr)*8Byte=32

	madd.d	t12,t12,a0,b1
	madd.d	t22,t22,a1,b1
	daddu	B,B,32
	daddu	PREA,PREA,32

	madd.d	t31,t31,a2,b0
	madd.d	t41,t41,a3,b0

	madd.d	t32,t32,a2,b1
	madd.d	t42,t42,a3,b1

.L49:							#  Write Back
	ldc1	c11,0(CO1)			#  Fetch 16 C
	ldc1	c21,8(CO1)			
	ldc1	c31,16(CO1)
	ldc1	c41,24(CO1)

	ldc1	c12,0(CO2)
	madd.d	t11,c11,t11,ALPHA
	ldc1	c22,8(CO2)
	madd.d	t21,c21,t21,ALPHA
	ldc1	c32,16(CO2)
	madd.d	t31,c31,t31,ALPHA
	ldc1	c42,24(CO2)
	madd.d	t41,c41,t41,ALPHA

	sdc1	t11,0(CO1)
	madd.d	t12,c12,t12,ALPHA
	sdc1	t21,8(CO1)
	madd.d	t22,c22,t22,ALPHA
	sdc1	t31,16(CO1)
	madd.d	t32,c32,t32,ALPHA
	sdc1	t41,24(CO1)
	madd.d	t42,c42,t42,ALPHA
	daddiu	M,M,-1				#  M--

	sdc1	t12,0(CO2)
	sdc1	t22,8(CO2)
	sdc1	t32,16(CO2)
	sdc1	t42,24(CO2)

	FETCH	$0,32(CO1)
	FETCH	$0,32(CO2)

	daddu	CO1,CO1,32			#  COx += 4*8Byte
	daddu	CO2,CO2,32
	bnez	M,.L40				#  M!=0
	move	B,BO				#  Reset B


.L12_M2:
	and		M,MCO,2				#  Remainder M = 2
	beqz	M,.L12_M1			
	nop

.L50:
	dsra	K,KCO,2				#  K=KCO/2
	dmtc1	$0,t11
	gsLQC1(R8,F1,F0,0)			#a0,a1

	mov.d	t21,t11
	mov.d	t12,t11
	gsLQC1(R9,F9,F8,0)			#b0,b1
	
	mov.d	t22,t11
	beqz	K,.L55
	nop

.L51:							#  N=2 m=2,=K=4
	gsLQC1(R8,F5,F4,1)			#  R8=A
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0

	gsLQC1(R9,F13,F12,1)		#  R9=B
	madd.d	t12,t12,a0,b1
	madd.d	t22,t22,a1,b1

	gsLQC1(R8,F3,F2,2)
	madd.d	t11,t11,a4,b4
	madd.d	t21,t21,a5,b4

	gsLQC1(R9,F11,F10,2)
	madd.d	t12,t12,a4,b5
	madd.d	t22,t22,a5,b5
	daddiu	K,K,-1

	gsLQC1(R8,F7,F6,3)	
	madd.d	t11,t11,a2,b2
	madd.d	t21,t21,a3,b2
	daddu	A,A,64				#  A+=2(mr)*4(kr)*8Byte=64

	gsLQC1(R9,F15,F14,3)
	madd.d	t12,t12,a2,b3
	madd.d	t22,t22,a3,b3
	daddu	B,B,64				#  B+=2(nr)*4(kr)*8Byte=128

	gsLQC1(R8,F1,F0,0)
	madd.d	t11,t11,a6,b6
	madd.d	t21,t21,a7,b6

	gsLQC1(R9,F9,F8,0)
	madd.d	t12,t12,a6,b7
	bnez 	K,.L51
	madd.d	t22,t22,a7,b7

.L55:							#  N=2 M=2 K=2
	and		K,KCO,2				#  k = KCO&2
	beqz	K,.L58
	nop

.L56:			
	gsLQC1(R8,F5,F4,1)			#  R8=A
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0
	daddu	A,A,32				#  A+=2(mr)*2(kr)*8Byte=32

	gsLQC1(R9,F13,F12,1)		#  R9=B
	madd.d	t12,t12,a0,b1
	madd.d	t22,t22,a1,b1
	daddu	B,B,32

.L57:
	gsLQC1(R8,F1,F0,0)
	madd.d	t11,t11,a4,b4
	madd.d	t21,t21,a5,b4

	gsLQC1(R9,F9,F8,0)
	madd.d	t12,t12,a4,b5
	madd.d	t22,t22,a5,b5

	
.L58:							#  N=2, M=2, K=1
	and		K,KCO,1
	beqz	K,.L59				#  
	ldc1	ALPHA,152($sp)		#  Get ALPHA
	
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0
	daddu	A,A,16				#  A+=2(mr)*1(kr)*8Byte=16
	daddu	B,B,16

	madd.d	t12,t12,a0,b1
	madd.d	t22,t22,a1,b1


.L59:							#  Write Back
	ldc1	c11,0(CO1)			#  Fetch 16 C
	ldc1	c21,8(CO1)			
	ldc1	c12,0(CO2)
	ldc1	c22,8(CO2)
	
	madd.d	t11,c11,t11,ALPHA
	madd.d	t21,c21,t21,ALPHA
	madd.d	t12,c12,t12,ALPHA
	madd.d	t22,c22,t22,ALPHA

	sdc1	t11,0(CO1)
	sdc1	t21,8(CO1)
	sdc1	t12,0(CO2)
	move	B,BO				#  Reset B
	sdc1	t22,8(CO2)
	daddu	CO1,CO1,16			#  COx += 2*8Byte
	daddu	CO2,CO2,16

	FETCH	$0,0(CO1)
	FETCH	$0,0(CO2)


.L12_M1:
	and		M,MCO,1				#  Remainder M = 1
	beqz	M,.L0_N2_Loop		#  M = 0, finishing one panel B
	nop

.L60:
	dsra	K,KCO,2				#  K=KCO/2
	dmtc1	$0,t11
	ldc1	a0,0(A)

	mov.d	t21,t11
	mov.d	t12,t11
	gsLQC1(R9,F9,F8,0)			#b0,b1
	
	mov.d	t22,t11
	beqz	K,.L65
	nop

.L61:							#  N=2 m=1,=K=4
	ldc1	a4,8(A)
	gsLQC1(R9,F13,F12,1)		#  R9=B
	ldc1	a2,16(A)
	madd.d	t11,t11,a0,b0
	madd.d	t12,t12,a0,b1

	gsLQC1(R9,F11,F10,2)
	madd.d	t11,t11,a4,b4
	madd.d	t12,t12,a4,b5
	daddiu	K,K,-1

	ldc1	a6,24(A)	
	madd.d	t11,t11,a2,b2
	daddu	A,A,32				#  A+=1(mr)*4(kr)*8Byte=64
	
	gsLQC1(R9,F15,F14,3)
	madd.d	t12,t12,a2,b3
	daddu	B,B,64				#  B+=2(nr)*4(kr)*8Byte=128

	ldc1	a0,0(A)
	gsLQC1(R9,F9,F8,0)
	madd.d	t11,t11,a6,b6
	bnez 	K,.L61
	madd.d	t12,t12,a6,b7

.L65:							#  N=2 M=1 K=2
	and		K,KCO,2				#  k = KCO&2
	beqz	K,.L68
	nop

.L66:			
	ldc1	a4,8(A)
	gsLQC1(R9,F13,F12,1)		#  R9=B
	madd.d	t11,t11,a0,b0
	madd.d	t12,t12,a0,b1
	daddu	A,A,16				#  A+=1(mr)*2(kr)*8Byte=32
	daddu	B,B,32

.L67:
	ldc1	a0,0(A)
	gsLQC1(R9,F9,F8,0)
	madd.d	t11,t11,a4,b4
	madd.d	t12,t12,a4,b5

	
.L68:							#  N=2, M=1, K=1
	and		K,KCO,1
	beqz	K,.L69				#  
	ldc1	ALPHA,152($sp)		#  Get ALPHA
	
	madd.d	t11,t11,a0,b0
	madd.d	t12,t12,a0,b1
	daddu	A,A,8				#  A+=1(mr)*1(kr)*8Byte=16
	daddu	B,B,16


.L69:							#  Write Back
	ldc1	c11,0(CO1)			#  Fetch 16 C
	ldc1	c12,0(CO2)
	
	madd.d	t11,c11,t11,ALPHA
	madd.d	t12,c12,t12,ALPHA

	sdc1	t11,0(CO1)
	move	B,BO				#  Reset B
	sdc1	t12,0(CO2)
	daddu	CO1,CO1,8			#  COx += 2*8Byte
	daddu	CO2,CO2,8

	FETCH	$0,0(CO1)
	FETCH	$0,0(CO2)


.L0_N2_Loop:
	daddu	BO,BO,SPANB			#  BO+=KC*2N
	move	B,BO				#  Set B
	daddu	C,C,SPANC			#  C+=LDC*2



	.align	5					
.L0_N1:
	and		N,NCO,1				#  Remainder N = 1
	beqz	N,.L999				#  N=0,NCO<1
	nop

	move	CO1,C				#  Set C	
	dsra	M,MCO,2				#  M=MCO/2
	
	move	A,AO				#  Reset A
	beqz	M,.L11_M2
	daddu	PREA,AO,SPANA


.L70:						
	dsra	K,KCO,2				#  K=KCO/2
	ldc1	b0,0(B)
	dmtc1	$0,t11
	gsLQC1(R8,F1,F0,0)			#a0,a1
	mov.d	t21,t11
	gsLQC1(R8,F3,F2,1)			#a2,a3
	mov.d	t31,t11
	beqz	K,.L75
	mov.d	t41,t11

.L71:							#  N=1,M=K=4
	ldc1	b4,8(B)
	gsLQC1(R8,F5,F4,2)			#  R8=A
	gsLQC1(R8,F7,F6,3)
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0

	FETCH		$0,(PREA)
	madd.d	t31,t31,a2,b0
	madd.d	t41,t41,a3,b0
	

.L72:
	ldc1	b2,16(B)
	gsLQC1(R8,F1,F0,4)
	gsLQC1(R8,F3,F2,5)
	madd.d	t11,t11,a4,b4
	madd.d	t21,t21,a5,b4

	FETCH		$0,32(PREA)
	madd.d	t31,t31,a6,b4
	madd.d	t41,t41,a7,b4


.L73:
	ldc1	b6,24(B)
	gsLQC1(R8,F5,F4,6)	
	gsLQC1(R8,F7,F6,7)
	madd.d	t11,t11,a0,b2
	madd.d	t21,t21,a1,b2
	daddu	B,B,32				#  B+=1(nr)*4(kr)*8Byte=64

	FETCH		$0,64(PREA)
	madd.d	t31,t31,a2,b2
	madd.d	t41,t41,a3,b2
	daddu	A,A,128				#  A+=4(mr)*4(kr)*8Byte=128

.L74:
	ldc1	b0,0(B)
	gsLQC1(R8,F1,F0,0)
	daddu	PREA,PREA,128
	gsLQC1(R8,F3,F2,1)
	madd.d	t11,t11,a4,b6
	madd.d	t21,t21,a5,b6
	daddiu	K,K,-1

	FETCH		$0,-32(PREA)
	madd.d	t31,t31,a6,b6
	bnez 	K,.L71
	madd.d	t41,t41,a7,b6



.L75:							#  N=2 M=4 K=2
	and		K,KCO,2				#  k = KCO&2
	beqz	K,.L78
	nop

.L76:			
	ldc1	b4,8(B)
	gsLQC1(R8,F5,F4,2)			#  R8=A
	gsLQC1(R8,F7,F6,3)
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0
	daddu	B,B,16				#  B+=1(nr)*2(kr)*8Byte=32

	FETCH		$0,0(PREA)
	madd.d	t31,t31,a2,b0
	madd.d	t41,t41,a3,b0
	daddu	A,A,64				#  A+=4(mr)*2(kr)*8Byte=64

.L77:
	ldc1	b0,0(B)
	gsLQC1(R8,F1,F0,0)
	gsLQC1(R8,F3,F2,1)
	madd.d	t11,t11,a4,b4
	madd.d	t21,t21,a5,b4

	FETCH		$0,32(PREA)
	madd.d	t31,t31,a6,b4
	madd.d	t41,t41,a7,b4
	daddu	PREA,PREA,64


	
.L78:							#  N=2, M=4, K=1
	and		K,KCO,1
	beqz	K,.L79				#  
	ldc1	ALPHA,152($sp)		#  Get ALPHA
	
	FETCH		$0,0(PREA)
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0
	daddu	A,A,32				#  A+=4(mr)*1(kr)*8Byte=32

	madd.d	t31,t31,a2,b0
	madd.d	t41,t41,a3,b0
	daddu	B,B,8
	daddu	PREA,PREA,32


.L79:							#  Write Back
	ldc1	c11,0(CO1)			#  Fetch 16 C
	ldc1	c21,8(CO1)			
	ldc1	c31,16(CO1)
	ldc1	c41,24(CO1)

	madd.d	t11,c11,t11,ALPHA
	madd.d	t21,c21,t21,ALPHA
	madd.d	t31,c31,t31,ALPHA
	madd.d	t41,c41,t41,ALPHA

	sdc1	t11,0(CO1)
	sdc1	t21,8(CO1)
	sdc1	t31,16(CO1)
	sdc1	t41,24(CO1)
	daddiu	M,M,-1				#  M--

	FETCH	$0,32(CO1)
	daddu	CO1,CO1,32			#  COx += 4*8Byte
	bnez	M,.L70				#  M!=0
	move	B,BO				#  Reset B



.L11_M2:
	and		M,MCO,2				#  Remainder M = 2
	beqz	M,.L11_M1			
	nop

.L80:						
	dsra	K,KCO,2				#  K=KCO/2
	ldc1	b0,0(B)
	dmtc1	$0,t11
	gsLQC1(R8,F1,F0,0)			#a0,a1
	mov.d	t21,t11
	beqz	K,.L85
	nop

.L81:							#  N=1,M=2,K=4
	ldc1	b4,8(B)
	gsLQC1(R8,F5,F4,1)			#  R8=A
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0

	ldc1	b2,16(B)
	gsLQC1(R8,F3,F2,2)
	madd.d	t11,t11,a4,b4
	madd.d	t21,t21,a5,b4

	ldc1	b6,24(B)
	daddu	B,B,32				#  B+=1(nr)*4(kr)*8Byte=32
	
	gsLQC1(R8,F7,F6,3)
	madd.d	t11,t11,a2,b2
	madd.d	t21,t21,a3,b2
	daddu	A,A,64				#  A+=2(mr)*4(kr)*8Byte=64

	ldc1	b0,0(B)
	daddiu	K,K,-1
	
	gsLQC1(R8,F1,F0,0)
	madd.d	t11,t11,a6,b6
	bnez 	K,.L81
	madd.d	t21,t21,a7,b6


.L85:							#  N=2 M=4 K=2
	and		K,KCO,2				#  k = KCO&2
	beqz	K,.L88
	nop

.L86:			
	ldc1	b4,8(B)
	daddu	B,B,16				#  B+=1(nr)*2(kr)*8Byte=16
	
	gsLQC1(R8,F5,F4,1)			#  R8=A
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0
	daddu	A,A,32				#  A+=2(mr)*2(kr)*8Byte=32

	ldc1	b0,0(B)
	gsLQC1(R8,F1,F0,0)
	madd.d	t11,t11,a4,b4
	madd.d	t21,t21,a5,b4

	
.L88:							#  N=2, M=4, K=1
	and		K,KCO,1
	beqz	K,.L89				#  
	ldc1	ALPHA,152($sp)		#  Get ALPHA
	
	madd.d	t11,t11,a0,b0
	madd.d	t21,t21,a1,b0
	daddu	A,A,16				#  A+=2(mr)*1(kr)*8Byte=16
	daddu	B,B,8


.L89:							#  Write Back
	ldc1	c11,0(CO1)			#  Fetch 16 C
	ldc1	c21,8(CO1)			

	madd.d	t11,c11,t11,ALPHA
	madd.d	t21,c21,t21,ALPHA

	sdc1	t11,0(CO1)
	sdc1	t21,8(CO1)

	FETCH	$0,16(CO1)
	daddu	CO1,CO1,16			#  COx += 2*8Byte
	move	B,BO				#  Reset B


.L11_M1:
	and		M,MCO,1				#  Remainder M = 1
	beqz	M,.L999				#  M = 0, End
	nop

.L90:						
	dsra	K,KCO,2				#  K=KCO/2
	ldc1	b0,0(B)
	ldc1	a0,0(A)
	beqz	K,.L95
	dmtc1	$0,t11

.L91:							#  N=1,M=1,K=4
	ldc1	b4,8(B)
	ldc1	a4,8(A)
	ldc1	b2,16(B)
	ldc1	a2,16(A)
	ldc1	b6,24(B)
	ldc1	a6,24(A)
	
	madd.d	t11,t11,a0,b0
	madd.d	t11,t11,a4,b4
	daddu	B,B,32				#  B+=1(nr)*4(kr)*8Byte=32
	daddu	A,A,32				#  A+=1(mr)*4(kr)*8Byte=32
	madd.d	t11,t11,a2,b2
	madd.d	t11,t11,a6,b6
	daddiu	K,K,-1
	
	ldc1	b0,0(B)
	bnez 	K,.L91
	ldc1	a0,0(A)


.L95:							#  N=2 M=4 K=2
	and		K,KCO,2				#  k = KCO&2
	beqz	K,.L98
	nop

.L96:			
	ldc1	b4,8(B)
	ldc1	a4,8(A)
	
	madd.d	t11,t11,a0,b0
	madd.d	t11,t11,a4,b4
	daddu	B,B,16				#  B+=1(nr)*2(kr)*8Byte=16
	daddu	A,A,16				#  A+=1(mr)*2(kr)*8Byte=32

	ldc1	b0,0(B)
	ldc1	a0,0(A)

	
.L98:							#  N=2, M=4, K=1
	and		K,KCO,1
	beqz	K,.L99				#  
	ldc1	ALPHA,152($sp)		#  Get ALPHA
	madd.d	t11,t11,a0,b0


.L99:							#  Write Back
	ldc1	c11,0(CO1)			#  Fetch 16 C
	madd.d	t11,c11,t11,ALPHA
	sdc1	t11,0(CO1)




.L999:							#  End
	ld	$16,   0($sp)
	ld	$17,   8($sp)
	ld	$18,  16($sp)
	ld	$19,  24($sp)
	ld	$20,  32($sp)
	ld	$21,  40($sp)
	ld	$22,  48($sp)
	ldc1	$f24, 56($sp)
	ldc1	$f25, 64($sp)
	ldc1	$f26, 72($sp)
	ldc1	$f27, 80($sp)
	ldc1	$f28, 88($sp)
	ld	$23,  96($sp)
	ld	$24, 104($sp)
	ld	$25, 112($sp)
	ldc1	$f20,120($sp)
	ldc1	$f21,128($sp)
	ldc1	$f22,136($sp)
	ldc1	$f23,144($sp)

	j	$31
	daddiu	$sp, $sp, 160

	EPILOGUE
