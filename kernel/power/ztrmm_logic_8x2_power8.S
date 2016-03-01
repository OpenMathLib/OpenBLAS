	srawi.		J,	N,	1
	ble		ZTRMM_L2_END

ZTRMM_L2_BEGIN:

	mr		CO,	C
	mr		AO,	A
	slwi		T1,	LDC	,	1
	add		C,	C,	T1

#if defined(LEFT)
	mr		KK,	OFFSET		// OFFSET -> KK
#endif

	srawi.		I,	M,	3
	ble		ZTRMM_L2x8_END

ZTRMM_L2x8_BEGIN:


#if  (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mr		BO,	B					// B -> BO
#else
	mr		BO,	B					// B -> BO
	slwi		T1,	KK,	5				// Number of values in B shifted
	slwi		T2,	KK,	7				// Number of values in A shifted
	add		BO,	BO,	T1				// Add values to BO
	add		AO,	AO,	T2				// Add values to AO
#endif

#if (defined(LEFT) &&  !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub		T1,	K,	KK				// K - KK -> TEMP1
#else
	mr		T1,	KK					// KK -> KTEMP
#ifdef LEFT
	addi		T1,	T1,	8				// KTEMP + Number of values in A -> KTEMP
#else
	addi		T1,	T1,	2				// KTEMP + Number of values in B -> KTEMP
#endif
#endif

	mr		KKK,	T1
	mr		K1,	T1
	srawi.		L,	K1,	3				// KTEMP / 8 -> L
	ble		ZTRMM_L2x8_SUB0
	cmpwi		cr0,	L,	1
	ble		ZTRMM_L2x8_SUB4

ZTRMM_L2x8_LOOP_START:

	dcbt		AO,	PRE
	LOAD2x8_1
	dcbt		AO,	PRE
	KERNEL2x8_I1
	dcbt		AO,	PRE
	KERNEL2x8_2
	dcbt		AO,	PRE
	KERNEL2x8_1
	dcbt		AO,	PRE
	KERNEL2x8_2

	dcbt		AO,	PRE
	KERNEL2x8_1
	dcbt		AO,	PRE
	KERNEL2x8_2
	dcbt		AO,	PRE
	KERNEL2x8_1
	dcbt		AO,	PRE
	KERNEL2x8_2

	addic.		L,	L,	-2
	ble		ZTRMM_L2x8_LOOP_END

	.align 5

ZTRMM_L2x8_LOOP:

	dcbt		AO,	PRE
	KERNEL2x8_1
	dcbt		AO,	PRE
	KERNEL2x8_2
	dcbt		AO,	PRE
	KERNEL2x8_1
	dcbt		AO,	PRE
	KERNEL2x8_2

	dcbt		AO,	PRE
	KERNEL2x8_1
	dcbt		AO,	PRE
	KERNEL2x8_2
	dcbt		AO,	PRE
	KERNEL2x8_1
	dcbt		AO,	PRE
	KERNEL2x8_2

	addic.		L,	L,	-1
	bgt		ZTRMM_L2x8_LOOP

ZTRMM_L2x8_LOOP_END:

	dcbt		AO,	PRE
	KERNEL2x8_1
	dcbt		AO,	PRE
	KERNEL2x8_2
	dcbt		AO,	PRE
	KERNEL2x8_1
	dcbt		AO,	PRE
	KERNEL2x8_2

	dcbt		AO,	PRE
	KERNEL2x8_1
	dcbt		AO,	PRE
	KERNEL2x8_2
	dcbt		AO,	PRE
	KERNEL2x8_1
	KERNEL2x8_E2

	b		ZTRMM_L2x8_SUB1

ZTRMM_L2x8_SUB4:

	dcbt		AO,	PRE
	KERNEL2x8_SUBI1
	dcbt		AO,	PRE
	KERNEL2x8_SUB1
	dcbt		AO,	PRE
	KERNEL2x8_SUB1
	dcbt		AO,	PRE
	KERNEL2x8_SUB1

	KERNEL2x8_SUB1
	KERNEL2x8_SUB1
	KERNEL2x8_SUB1
	KERNEL2x8_SUB1

	b		ZTRMM_L2x8_SUB1

ZTRMM_L2x8_SUB0:

	andi.		L,	K1,	7						// K1 & 7 -> L

	KERNEL2x8_SUBI1

	addic.		L,	L,	-1
	ble		ZTRMM_L2x8_SAVE
	b		ZTRMM_L2x8_SUB2

ZTRMM_L2x8_SUB1:

	andi.		L,	K1,	7						// K1 & 7 -> L
	ble		ZTRMM_L2x8_SAVE

ZTRMM_L2x8_SUB2:

	KERNEL2x8_SUB1

	addic.		L,	L,	-1
	bgt		ZTRMM_L2x8_SUB2

ZTRMM_L2x8_SAVE:

	SAVE2x8

#if  (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub		T1,	K,	KKK					// K - KKK -> TEMP1
	slwi		T2,	T1,	5			// TEMP1 * Number of values in B shifted -> TEMP2
	slwi		T1,	T1,	7			// TEMP1 * Number of values in A shifted -> TEMP1
	add		BO,	BO,	T2					// BO += TEMP2 * number of values in B shifted
	add		AO,	AO,	T1					// AO += TEMP1 * number of values in A shifted
#endif

#if defined(LEFT)
	addi		KK,	KK,	8				// KK += Number of values in A
#endif


	addic.		I,	I,	-1
	bgt		ZTRMM_L2x8_BEGIN

ZTRMM_L2x8_END:

ZTRMM_L2x4_BEGIN:
	andi.		T2,	M,	7
	ble		ZTRMM_L2x1_END

	andi.		T1,	M,	4
	ble		ZTRMM_L2x4_END

#if  (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mr		BO,	B					// B -> BO
#else
	mr		BO,	B					// B -> BO
	slwi		T1,	KK,	5				// Number of values in B shifted
	slwi		T2,	KK,	6				// Number of values in A shifted
	add		BO,	BO,	T1				// Add values to BO
	add		AO,	AO,	T2				// Add values to AO
#endif

#if (defined(LEFT) &&  !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub		T1,	K,	KK				// K - KK -> TEMP1
#else
	mr		T1,	KK					// KK -> KTEMP
#ifdef LEFT
	addi		T1,	T1,	4				// KTEMP + Number of values in A -> KTEMP
#else
	addi		T1,	T1,	2				// KTEMP + Number of values in B -> KTEMP
#endif
#endif

	mr		KKK,	T1
	mr		K1,	T1
	srawi.		L,	K1,	3				// KTEMP / 8 -> L
	ble		ZTRMM_L2x4_SUB0
	cmpwi		cr0,	L,	1
	ble		ZTRMM_L2x4_SUB4

ZTRMM_L2x4_LOOP_START:

	LOAD2x4_1
	KERNEL2x4_I1
	KERNEL2x4_2
	KERNEL2x4_1
	KERNEL2x4_2

	KERNEL2x4_1
	KERNEL2x4_2
	KERNEL2x4_1
	KERNEL2x4_2

	addic.		L,	L,	-2
	ble		ZTRMM_L2x4_LOOP_END

	.align 5

ZTRMM_L2x4_LOOP:

	KERNEL2x4_1
	KERNEL2x4_2
	KERNEL2x4_1
	KERNEL2x4_2

	KERNEL2x4_1
	KERNEL2x4_2
	KERNEL2x4_1
	KERNEL2x4_2

	addic.		L,	L,	-1
	bgt		ZTRMM_L2x4_LOOP

ZTRMM_L2x4_LOOP_END:

	KERNEL2x4_1
	KERNEL2x4_2
	KERNEL2x4_1
	KERNEL2x4_2

	KERNEL2x4_1
	KERNEL2x4_2
	KERNEL2x4_1
	KERNEL2x4_E2

	b		ZTRMM_L2x4_SUB1

ZTRMM_L2x4_SUB4:

	KERNEL2x4_SUBI1
	KERNEL2x4_SUB1
	KERNEL2x4_SUB1
	KERNEL2x4_SUB1

	KERNEL2x4_SUB1
	KERNEL2x4_SUB1
	KERNEL2x4_SUB1
	KERNEL2x4_SUB1

	b		ZTRMM_L2x4_SUB1

ZTRMM_L2x4_SUB0:

	andi.		L,	K1,	7						// K1 & 7 -> L

	KERNEL2x4_SUBI1

	addic.		L,	L,	-1
	ble		ZTRMM_L2x4_SAVE
	b		ZTRMM_L2x4_SUB2

ZTRMM_L2x4_SUB1:

	andi.		L,	K1,	7						// K1 & 7 -> L
	ble		ZTRMM_L2x4_SAVE

ZTRMM_L2x4_SUB2:

	KERNEL2x4_SUB1

	addic.		L,	L,	-1
	bgt		ZTRMM_L2x4_SUB2

ZTRMM_L2x4_SAVE:

	SAVE2x4

#if  (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub		T1,	K,	KKK					// K - KKK -> TEMP1
	slwi		T2,	T1,	5			// TEMP1 * Number of values in B shifted -> TEMP2
	slwi		T1,	T1,	6			// TEMP1 * Number of values in A shifted -> TEMP1
	add		BO,	BO,	T2					// BO += TEMP2 * number of values in B shifted
	add		AO,	AO,	T1					// AO += TEMP1 * number of values in A shifted
#endif

#if defined(LEFT)
	addi		KK,	KK,	4				// KK += Number of values in A
#endif


ZTRMM_L2x4_END:

ZTRMM_L2x2_BEGIN:

	andi.		T1,	M,	2
	ble		ZTRMM_L2x2_END

#if  (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mr		BO,	B					// B -> BO
#else
	mr		BO,	B					// B -> BO
	slwi		T1,	KK,	5				// Number of values in B shifted
	slwi		T2,	KK,	5				// Number of values in A shifted
	add		BO,	BO,	T1				// Add values to BO
	add		AO,	AO,	T2				// Add values to AO
#endif

#if (defined(LEFT) &&  !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub		T1,	K,	KK				// K - KK -> TEMP1
#else
	mr		T1,	KK					// KK -> KTEMP
#ifdef LEFT
	addi		T1,	T1,	2				// KTEMP + Number of values in A -> KTEMP
#else
	addi		T1,	T1,	2				// KTEMP + Number of values in B -> KTEMP
#endif
#endif

	mr		KKK,	T1
	mr		K1,	T1
	srawi.		L,	K1,	3				// KTEMP / 8 -> L
	ble		ZTRMM_L2x2_SUB0
	cmpwi		cr0,	L,	1
	ble		ZTRMM_L2x2_SUB4

ZTRMM_L2x2_LOOP_START:

	LOAD2x2_1
	KERNEL2x2_I1
	KERNEL2x2_2
	KERNEL2x2_1
	KERNEL2x2_2

	KERNEL2x2_1
	KERNEL2x2_2
	KERNEL2x2_1
	KERNEL2x2_2

	addic.		L,	L,	-2
	ble		ZTRMM_L2x2_LOOP_END

	.align 5

ZTRMM_L2x2_LOOP:

	KERNEL2x2_1
	KERNEL2x2_2
	KERNEL2x2_1
	KERNEL2x2_2

	KERNEL2x2_1
	KERNEL2x2_2
	KERNEL2x2_1
	KERNEL2x2_2

	addic.		L,	L,	-1
	bgt		ZTRMM_L2x2_LOOP

ZTRMM_L2x2_LOOP_END:

	KERNEL2x2_1
	KERNEL2x2_2
	KERNEL2x2_1
	KERNEL2x2_2

	KERNEL2x2_1
	KERNEL2x2_2
	KERNEL2x2_1
	KERNEL2x2_E2

	b		ZTRMM_L2x2_SUB1

ZTRMM_L2x2_SUB4:

	KERNEL2x2_SUBI1
	KERNEL2x2_SUB1
	KERNEL2x2_SUB1
	KERNEL2x2_SUB1

	KERNEL2x2_SUB1
	KERNEL2x2_SUB1
	KERNEL2x2_SUB1
	KERNEL2x2_SUB1

	b		ZTRMM_L2x2_SUB1

ZTRMM_L2x2_SUB0:

	andi.		L,	K1,	7						// K1 & 7 -> L

	KERNEL2x2_SUBI1

	addic.		L,	L,	-1
	ble		ZTRMM_L2x2_SAVE
	b		ZTRMM_L2x2_SUB2

ZTRMM_L2x2_SUB1:

	andi.		L,	K1,	7						// K1 & 7 -> L
	ble		ZTRMM_L2x2_SAVE

ZTRMM_L2x2_SUB2:

	KERNEL2x2_SUB1

	addic.		L,	L,	-1
	bgt		ZTRMM_L2x2_SUB2

ZTRMM_L2x2_SAVE:

	SAVE2x2

#if  (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub		T1,	K,	KKK					// K - KKK -> TEMP1
	slwi		T2,	T1,	5			// TEMP1 * Number of values in B shifted -> TEMP2
	slwi		T1,	T1,	5			// TEMP1 * Number of values in A shifted -> TEMP1
	add		BO,	BO,	T2					// BO += TEMP2 * number of values in B shifted
	add		AO,	AO,	T1					// AO += TEMP1 * number of values in A shifted
#endif

#if defined(LEFT)
	addi		KK,	KK,	2				// KK += Number of values in A
#endif


ZTRMM_L2x2_END:

ZTRMM_L2x1_BEGIN:

	andi.		T1,	M,	1
	ble		ZTRMM_L2x1_END

#if  (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mr		BO,	B					// B -> BO
#else
	mr		BO,	B					// B -> BO
	slwi		T1,	KK,	5				// Number of values in B shifted
	slwi		T2,	KK,	4				// Number of values in A shifted
	add		BO,	BO,	T1				// Add values to BO
	add		AO,	AO,	T2				// Add values to AO
#endif

#if (defined(LEFT) &&  !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub		T1,	K,	KK				// K - KK -> TEMP1
#else
	mr		T1,	KK					// KK -> KTEMP
#ifdef LEFT
	addi		T1,	T1,	1				// KTEMP + Number of values in A -> KTEMP
#else
	addi		T1,	T1,	2				// KTEMP + Number of values in B -> KTEMP
#endif
#endif

	mr		KKK,	T1
	mr		K1,	T1
	srawi.		L,	K1,	3				// KTEMP / 8 -> L
	ble		ZTRMM_L2x1_SUB0
	cmpwi		cr0,	L,	1
	ble		ZTRMM_L2x1_SUB4

ZTRMM_L2x1_LOOP_START:

	LOAD2x1_1
	KERNEL2x1_I1
	KERNEL2x1_2
	KERNEL2x1_1
	KERNEL2x1_2

	KERNEL2x1_1
	KERNEL2x1_2
	KERNEL2x1_1
	KERNEL2x1_2

	addic.		L,	L,	-2
	ble		ZTRMM_L2x1_LOOP_END

	.align 5

ZTRMM_L2x1_LOOP:

	KERNEL2x1_1
	KERNEL2x1_2
	KERNEL2x1_1
	KERNEL2x1_2

	KERNEL2x1_1
	KERNEL2x1_2
	KERNEL2x1_1
	KERNEL2x1_2

	addic.		L,	L,	-1
	bgt		ZTRMM_L2x1_LOOP

ZTRMM_L2x1_LOOP_END:

	KERNEL2x1_1
	KERNEL2x1_2
	KERNEL2x1_1
	KERNEL2x1_2

	KERNEL2x1_1
	KERNEL2x1_2
	KERNEL2x1_1
	KERNEL2x1_E2

	b		ZTRMM_L2x1_SUB1

ZTRMM_L2x1_SUB4:

	KERNEL2x1_SUBI1
	KERNEL2x1_SUB1
	KERNEL2x1_SUB1
	KERNEL2x1_SUB1

	KERNEL2x1_SUB1
	KERNEL2x1_SUB1
	KERNEL2x1_SUB1
	KERNEL2x1_SUB1

	b		ZTRMM_L2x1_SUB1

ZTRMM_L2x1_SUB0:

	andi.		L,	K1,	7						// K1 & 7 -> L

	KERNEL2x1_SUBI1

	addic.		L,	L,	-1
	ble		ZTRMM_L2x1_SAVE
	b		ZTRMM_L2x1_SUB2

ZTRMM_L2x1_SUB1:

	andi.		L,	K1,	7						// K1 & 7 -> L
	ble		ZTRMM_L2x1_SAVE

ZTRMM_L2x1_SUB2:

	KERNEL2x1_SUB1

	addic.		L,	L,	-1
	bgt		ZTRMM_L2x1_SUB2

ZTRMM_L2x1_SAVE:

	SAVE2x1

#if  (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub		T1,	K,	KKK					// K - KKK -> TEMP1
	slwi		T2,	T1,	5			// TEMP1 * Number of values in B shifted -> TEMP2
	slwi		T1,	T1,	4			// TEMP1 * Number of values in A shifted -> TEMP1
	add		BO,	BO,	T2					// BO += TEMP2 * number of values in B shifted
	add		AO,	AO,	T1					// AO += TEMP1 * number of values in A shifted
#endif

#if defined(LEFT)
	addi		KK,	KK,	1				// KK += Number of values in A
#endif


ZTRMM_L2x1_END:

	slwi		T1,	K,	5
	add		B,	B,	T1

#if !defined(LEFT)
	addi		KK,	KK,	2					// KK += Number of values in B
#endif


	addic.		J,	J,	-1
	bgt		ZTRMM_L2_BEGIN

	andi.		T2,	N,	1
	ble		L999

ZTRMM_L2_END:

	b		ZTRMM_L1_BEGIN

L999_H1:

	b		L999

ZTRMM_L1_BEGIN:

	andi.		T1,	N,	1
	ble		ZTRMM_L1_END
	mr		CO,	C
	mr		AO,	A

#if defined(LEFT)
	mr		KK,	OFFSET		// OFFSET -> KK
#endif

	srawi.		I,	M,	3
	ble		ZTRMM_L1x8_END

ZTRMM_L1x8_BEGIN:


#if  (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mr		BO,	B					// B -> BO
#else
	mr		BO,	B					// B -> BO
	slwi		T1,	KK,	4				// Number of values in B shifted
	slwi		T2,	KK,	7				// Number of values in A shifted
	add		BO,	BO,	T1				// Add values to BO
	add		AO,	AO,	T2				// Add values to AO
#endif

#if (defined(LEFT) &&  !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub		T1,	K,	KK				// K - KK -> TEMP1
#else
	mr		T1,	KK					// KK -> KTEMP
#ifdef LEFT
	addi		T1,	T1,	8				// KTEMP + Number of values in A -> KTEMP
#else
	addi		T1,	T1,	1				// KTEMP + Number of values in B -> KTEMP
#endif
#endif

	mr		KKK,	T1
	mr		K1,	T1
	srawi.		L,	K1,	3				// KTEMP / 8 -> L
	ble		ZTRMM_L1x8_SUB0
	cmpwi		cr0,	L,	1
	ble		ZTRMM_L1x8_SUB4

ZTRMM_L1x8_LOOP_START:

	dcbt		AO,	PRE
	LOAD1x8_1
	dcbt		AO,	PRE
	KERNEL1x8_I1
	dcbt		AO,	PRE
	KERNEL1x8_2
	dcbt		AO,	PRE
	KERNEL1x8_1
	dcbt		AO,	PRE
	KERNEL1x8_2

	dcbt		AO,	PRE
	KERNEL1x8_1
	dcbt		AO,	PRE
	KERNEL1x8_2
	dcbt		AO,	PRE
	KERNEL1x8_1
	dcbt		AO,	PRE
	KERNEL1x8_2

	addic.		L,	L,	-2
	ble		ZTRMM_L1x8_LOOP_END

	.align 5

ZTRMM_L1x8_LOOP:

	dcbt		AO,	PRE
	KERNEL1x8_1
	dcbt		AO,	PRE
	KERNEL1x8_2
	dcbt		AO,	PRE
	KERNEL1x8_1
	dcbt		AO,	PRE
	KERNEL1x8_2

	dcbt		AO,	PRE
	KERNEL1x8_1
	dcbt		AO,	PRE
	KERNEL1x8_2
	dcbt		AO,	PRE
	KERNEL1x8_1
	dcbt		AO,	PRE
	KERNEL1x8_2

	addic.		L,	L,	-1
	bgt		ZTRMM_L1x8_LOOP

ZTRMM_L1x8_LOOP_END:

	dcbt		AO,	PRE
	KERNEL1x8_1
	dcbt		AO,	PRE
	KERNEL1x8_2
	dcbt		AO,	PRE
	KERNEL1x8_1
	dcbt		AO,	PRE
	KERNEL1x8_2

	dcbt		AO,	PRE
	KERNEL1x8_1
	dcbt		AO,	PRE
	KERNEL1x8_2
	dcbt		AO,	PRE
	KERNEL1x8_1
	KERNEL1x8_E2

	b		ZTRMM_L1x8_SUB1

ZTRMM_L1x8_SUB4:

	dcbt		AO,	PRE
	KERNEL1x8_SUBI1
	dcbt		AO,	PRE
	KERNEL1x8_SUB1
	dcbt		AO,	PRE
	KERNEL1x8_SUB1
	dcbt		AO,	PRE
	KERNEL1x8_SUB1

	KERNEL1x8_SUB1
	KERNEL1x8_SUB1
	KERNEL1x8_SUB1
	KERNEL1x8_SUB1

	b		ZTRMM_L1x8_SUB1

ZTRMM_L1x8_SUB0:

	andi.		L,	K1,	7						// K1 & 7 -> L

	KERNEL1x8_SUBI1

	addic.		L,	L,	-1
	ble		ZTRMM_L1x8_SAVE
	b		ZTRMM_L1x8_SUB2

ZTRMM_L1x8_SUB1:

	andi.		L,	K1,	7						// K1 & 7 -> L
	ble		ZTRMM_L1x8_SAVE

ZTRMM_L1x8_SUB2:

	KERNEL1x8_SUB1

	addic.		L,	L,	-1
	bgt		ZTRMM_L1x8_SUB2

ZTRMM_L1x8_SAVE:

	SAVE1x8

#if  (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub		T1,	K,	KKK					// K - KKK -> TEMP1
	slwi		T2,	T1,	4			// TEMP1 * Number of values in B shifted -> TEMP2
	slwi		T1,	T1,	7			// TEMP1 * Number of values in A shifted -> TEMP1
	add		BO,	BO,	T2					// BO += TEMP2 * number of values in B shifted
	add		AO,	AO,	T1					// AO += TEMP1 * number of values in A shifted
#endif

#if defined(LEFT)
	addi		KK,	KK,	8				// KK += Number of values in A
#endif


	addic.		I,	I,	-1
	bgt		ZTRMM_L1x8_BEGIN

ZTRMM_L1x8_END:

ZTRMM_L1x4_BEGIN:
	andi.		T2,	M,	7
	ble		ZTRMM_L1x1_END

	andi.		T1,	M,	4
	ble		ZTRMM_L1x4_END

#if  (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mr		BO,	B					// B -> BO
#else
	mr		BO,	B					// B -> BO
	slwi		T1,	KK,	4				// Number of values in B shifted
	slwi		T2,	KK,	6				// Number of values in A shifted
	add		BO,	BO,	T1				// Add values to BO
	add		AO,	AO,	T2				// Add values to AO
#endif

#if (defined(LEFT) &&  !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub		T1,	K,	KK				// K - KK -> TEMP1
#else
	mr		T1,	KK					// KK -> KTEMP
#ifdef LEFT
	addi		T1,	T1,	4				// KTEMP + Number of values in A -> KTEMP
#else
	addi		T1,	T1,	1				// KTEMP + Number of values in B -> KTEMP
#endif
#endif

	mr		KKK,	T1
	mr		K1,	T1
	srawi.		L,	K1,	3				// KTEMP / 8 -> L
	ble		ZTRMM_L1x4_SUB0
	cmpwi		cr0,	L,	1
	ble		ZTRMM_L1x4_SUB4

ZTRMM_L1x4_LOOP_START:

	LOAD1x4_1
	KERNEL1x4_I1
	KERNEL1x4_2
	KERNEL1x4_1
	KERNEL1x4_2

	KERNEL1x4_1
	KERNEL1x4_2
	KERNEL1x4_1
	KERNEL1x4_2

	addic.		L,	L,	-2
	ble		ZTRMM_L1x4_LOOP_END

	.align 5

ZTRMM_L1x4_LOOP:

	KERNEL1x4_1
	KERNEL1x4_2
	KERNEL1x4_1
	KERNEL1x4_2

	KERNEL1x4_1
	KERNEL1x4_2
	KERNEL1x4_1
	KERNEL1x4_2

	addic.		L,	L,	-1
	bgt		ZTRMM_L1x4_LOOP

ZTRMM_L1x4_LOOP_END:

	KERNEL1x4_1
	KERNEL1x4_2
	KERNEL1x4_1
	KERNEL1x4_2

	KERNEL1x4_1
	KERNEL1x4_2
	KERNEL1x4_1
	KERNEL1x4_E2

	b		ZTRMM_L1x4_SUB1

ZTRMM_L1x4_SUB4:

	KERNEL1x4_SUBI1
	KERNEL1x4_SUB1
	KERNEL1x4_SUB1
	KERNEL1x4_SUB1

	KERNEL1x4_SUB1
	KERNEL1x4_SUB1
	KERNEL1x4_SUB1
	KERNEL1x4_SUB1

	b		ZTRMM_L1x4_SUB1

ZTRMM_L1x4_SUB0:

	andi.		L,	K1,	7						// K1 & 7 -> L

	KERNEL1x4_SUBI1

	addic.		L,	L,	-1
	ble		ZTRMM_L1x4_SAVE
	b		ZTRMM_L1x4_SUB2

ZTRMM_L1x4_SUB1:

	andi.		L,	K1,	7						// K1 & 7 -> L
	ble		ZTRMM_L1x4_SAVE

ZTRMM_L1x4_SUB2:

	KERNEL1x4_SUB1

	addic.		L,	L,	-1
	bgt		ZTRMM_L1x4_SUB2

ZTRMM_L1x4_SAVE:

	SAVE1x4

#if  (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub		T1,	K,	KKK					// K - KKK -> TEMP1
	slwi		T2,	T1,	4			// TEMP1 * Number of values in B shifted -> TEMP2
	slwi		T1,	T1,	6			// TEMP1 * Number of values in A shifted -> TEMP1
	add		BO,	BO,	T2					// BO += TEMP2 * number of values in B shifted
	add		AO,	AO,	T1					// AO += TEMP1 * number of values in A shifted
#endif

#if defined(LEFT)
	addi		KK,	KK,	4				// KK += Number of values in A
#endif


ZTRMM_L1x4_END:

ZTRMM_L1x2_BEGIN:

	andi.		T1,	M,	2
	ble		ZTRMM_L1x2_END

#if  (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mr		BO,	B					// B -> BO
#else
	mr		BO,	B					// B -> BO
	slwi		T1,	KK,	4				// Number of values in B shifted
	slwi		T2,	KK,	5				// Number of values in A shifted
	add		BO,	BO,	T1				// Add values to BO
	add		AO,	AO,	T2				// Add values to AO
#endif

#if (defined(LEFT) &&  !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub		T1,	K,	KK				// K - KK -> TEMP1
#else
	mr		T1,	KK					// KK -> KTEMP
#ifdef LEFT
	addi		T1,	T1,	2				// KTEMP + Number of values in A -> KTEMP
#else
	addi		T1,	T1,	1				// KTEMP + Number of values in B -> KTEMP
#endif
#endif

	mr		KKK,	T1
	mr		K1,	T1
	srawi.		L,	K1,	3				// KTEMP / 8 -> L
	ble		ZTRMM_L1x2_SUB0
	cmpwi		cr0,	L,	1
	ble		ZTRMM_L1x2_SUB4

ZTRMM_L1x2_LOOP_START:

	LOAD1x2_1
	KERNEL1x2_I1
	KERNEL1x2_2
	KERNEL1x2_1
	KERNEL1x2_2

	KERNEL1x2_1
	KERNEL1x2_2
	KERNEL1x2_1
	KERNEL1x2_2

	addic.		L,	L,	-2
	ble		ZTRMM_L1x2_LOOP_END

	.align 5

ZTRMM_L1x2_LOOP:

	KERNEL1x2_1
	KERNEL1x2_2
	KERNEL1x2_1
	KERNEL1x2_2

	KERNEL1x2_1
	KERNEL1x2_2
	KERNEL1x2_1
	KERNEL1x2_2

	addic.		L,	L,	-1
	bgt		ZTRMM_L1x2_LOOP

ZTRMM_L1x2_LOOP_END:

	KERNEL1x2_1
	KERNEL1x2_2
	KERNEL1x2_1
	KERNEL1x2_2

	KERNEL1x2_1
	KERNEL1x2_2
	KERNEL1x2_1
	KERNEL1x2_E2

	b		ZTRMM_L1x2_SUB1

ZTRMM_L1x2_SUB4:

	KERNEL1x2_SUBI1
	KERNEL1x2_SUB1
	KERNEL1x2_SUB1
	KERNEL1x2_SUB1

	KERNEL1x2_SUB1
	KERNEL1x2_SUB1
	KERNEL1x2_SUB1
	KERNEL1x2_SUB1

	b		ZTRMM_L1x2_SUB1

ZTRMM_L1x2_SUB0:

	andi.		L,	K1,	7						// K1 & 7 -> L

	KERNEL1x2_SUBI1

	addic.		L,	L,	-1
	ble		ZTRMM_L1x2_SAVE
	b		ZTRMM_L1x2_SUB2

ZTRMM_L1x2_SUB1:

	andi.		L,	K1,	7						// K1 & 7 -> L
	ble		ZTRMM_L1x2_SAVE

ZTRMM_L1x2_SUB2:

	KERNEL1x2_SUB1

	addic.		L,	L,	-1
	bgt		ZTRMM_L1x2_SUB2

ZTRMM_L1x2_SAVE:

	SAVE1x2

#if  (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub		T1,	K,	KKK					// K - KKK -> TEMP1
	slwi		T2,	T1,	4			// TEMP1 * Number of values in B shifted -> TEMP2
	slwi		T1,	T1,	5			// TEMP1 * Number of values in A shifted -> TEMP1
	add		BO,	BO,	T2					// BO += TEMP2 * number of values in B shifted
	add		AO,	AO,	T1					// AO += TEMP1 * number of values in A shifted
#endif

#if defined(LEFT)
	addi		KK,	KK,	2				// KK += Number of values in A
#endif


ZTRMM_L1x2_END:

ZTRMM_L1x1_BEGIN:

	andi.		T1,	M,	1
	ble		ZTRMM_L1x1_END

#if  (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	mr		BO,	B					// B -> BO
#else
	mr		BO,	B					// B -> BO
	slwi		T1,	KK,	4				// Number of values in B shifted
	slwi		T2,	KK,	4				// Number of values in A shifted
	add		BO,	BO,	T1				// Add values to BO
	add		AO,	AO,	T2				// Add values to AO
#endif

#if (defined(LEFT) &&  !defined(TRANSA)) || (!defined(LEFT) && defined(TRANSA))
	sub		T1,	K,	KK				// K - KK -> TEMP1
#else
	mr		T1,	KK					// KK -> KTEMP
#ifdef LEFT
	addi		T1,	T1,	1				// KTEMP + Number of values in A -> KTEMP
#else
	addi		T1,	T1,	1				// KTEMP + Number of values in B -> KTEMP
#endif
#endif

	mr		KKK,	T1
	mr		K1,	T1
	srawi.		L,	K1,	3				// KTEMP / 8 -> L
	ble		ZTRMM_L1x1_SUB0
	cmpwi		cr0,	L,	1
	ble		ZTRMM_L1x1_SUB4

ZTRMM_L1x1_LOOP_START:

	LOAD1x1_1
	KERNEL1x1_I1
	KERNEL1x1_2
	KERNEL1x1_1
	KERNEL1x1_2

	KERNEL1x1_1
	KERNEL1x1_2
	KERNEL1x1_1
	KERNEL1x1_2

	addic.		L,	L,	-2
	ble		ZTRMM_L1x1_LOOP_END

	.align 5

ZTRMM_L1x1_LOOP:

	KERNEL1x1_1
	KERNEL1x1_2
	KERNEL1x1_1
	KERNEL1x1_2

	KERNEL1x1_1
	KERNEL1x1_2
	KERNEL1x1_1
	KERNEL1x1_2

	addic.		L,	L,	-1
	bgt		ZTRMM_L1x1_LOOP

ZTRMM_L1x1_LOOP_END:

	KERNEL1x1_1
	KERNEL1x1_2
	KERNEL1x1_1
	KERNEL1x1_2

	KERNEL1x1_1
	KERNEL1x1_2
	KERNEL1x1_1
	KERNEL1x1_E2

	b		ZTRMM_L1x1_SUB1

ZTRMM_L1x1_SUB4:

	KERNEL1x1_SUBI1
	KERNEL1x1_SUB1
	KERNEL1x1_SUB1
	KERNEL1x1_SUB1

	KERNEL1x1_SUB1
	KERNEL1x1_SUB1
	KERNEL1x1_SUB1
	KERNEL1x1_SUB1

	b		ZTRMM_L1x1_SUB1

ZTRMM_L1x1_SUB0:

	andi.		L,	K1,	7						// K1 & 7 -> L

	KERNEL1x1_SUBI1

	addic.		L,	L,	-1
	ble		ZTRMM_L1x1_SAVE
	b		ZTRMM_L1x1_SUB2

ZTRMM_L1x1_SUB1:

	andi.		L,	K1,	7						// K1 & 7 -> L
	ble		ZTRMM_L1x1_SAVE

ZTRMM_L1x1_SUB2:

	KERNEL1x1_SUB1

	addic.		L,	L,	-1
	bgt		ZTRMM_L1x1_SUB2

ZTRMM_L1x1_SAVE:

	SAVE1x1

#if  (defined(LEFT) &&  defined(TRANSA)) || (!defined(LEFT) && !defined(TRANSA))
	sub		T1,	K,	KKK					// K - KKK -> TEMP1
	slwi		T2,	T1,	4			// TEMP1 * Number of values in B shifted -> TEMP2
	slwi		T1,	T1,	4			// TEMP1 * Number of values in A shifted -> TEMP1
	add		BO,	BO,	T2					// BO += TEMP2 * number of values in B shifted
	add		AO,	AO,	T1					// AO += TEMP1 * number of values in A shifted
#endif

#if defined(LEFT)
	addi		KK,	KK,	1				// KK += Number of values in A
#endif


ZTRMM_L1x1_END:

#if !defined(LEFT)
	addi		KK,	KK,	1					// KK += Number of values in B
#endif


ZTRMM_L1_END:
