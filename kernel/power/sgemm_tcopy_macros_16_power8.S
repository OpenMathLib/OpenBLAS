
/**********************************************************************************************
* Macros for N=4 and M=16
**********************************************************************************************/

.macro COPY_4x16

	lxvw4x		vs32,	o0,	A0
	lxvw4x		vs33,	o16,	A0
	lxvw4x		vs34,	o32,	A0
	lxvw4x		vs35,	o48,	A0

	lxvw4x		vs36,	o0,	A1
	lxvw4x		vs37,	o16,	A1
	lxvw4x		vs38,	o32,	A1
	lxvw4x		vs39,	o48,	A1

	lxvw4x		vs40,	o0,	A2
	lxvw4x		vs41,	o16,	A2
	lxvw4x		vs42,	o32,	A2
	lxvw4x		vs43,	o48,	A2

	lxvw4x		vs44,	o0,	A3
	lxvw4x		vs45,	o16,	A3
	lxvw4x		vs46,	o32,	A3
	lxvw4x		vs47,	o48,	A3

	mr		T1,	BO

	stxvw4x		vs32,	o0,	T1
	stxvw4x		vs33,	o16,	T1
	stxvw4x		vs34,	o32,	T1
	stxvw4x		vs35,	o48,	T1

	addi		T1,	T1,	64

	stxvw4x		vs36,	o0,	T1
	stxvw4x		vs37,	o16,	T1
	stxvw4x		vs38,	o32,	T1
	stxvw4x		vs39,	o48,	T1

	addi		T1,	T1,	64

	stxvw4x		vs40,	o0,	T1
	stxvw4x		vs41,	o16,	T1
	stxvw4x		vs42,	o32,	T1
	stxvw4x		vs43,	o48,	T1

	addi		T1,	T1,	64

	stxvw4x		vs44,	o0,	T1
	stxvw4x		vs45,	o16,	T1
	stxvw4x		vs46,	o32,	T1
	stxvw4x		vs47,	o48,	T1

.endm

/**********************************************************************************************
* Macros for N=4 and M=8
**********************************************************************************************/

.macro COPY_4x8

	lxvw4x		vs32,	o0,	A0
	lxvw4x		vs33,	o16,	A0

	lxvw4x		vs34,	o0,	A1
	lxvw4x		vs35,	o16,	A1

	lxvw4x		vs36,	o0,	A2
	lxvw4x		vs37,	o16,	A2

	lxvw4x		vs38,	o0,	A3
	lxvw4x		vs39,	o16,	A3

	mr		T1,	BO

	stxvw4x		vs32,	o0,	T1
	stxvw4x		vs33,	o16,	T1

	stxvw4x		vs34,	o32,	T1
	stxvw4x		vs35,	o48,	T1

	addi		T1,	T1,	64

	stxvw4x		vs36,	o0,	T1
	stxvw4x		vs37,	o16,	T1

	stxvw4x		vs38,	o32,	T1
	stxvw4x		vs39,	o48,	T1

.endm

/**********************************************************************************************
* Macros for N=4 and M=4
**********************************************************************************************/

.macro COPY_4x4

	lxvw4x		vs32,	o0,	A0

	lxvw4x		vs33,	o0,	A1

	lxvw4x		vs34,	o0,	A2

	lxvw4x		vs35,	o0,	A3

	mr		T1,	BO

	stxvw4x		vs32,	o0,	T1

	stxvw4x		vs33,	o16,	T1

	stxvw4x		vs34,	o32,	T1

	stxvw4x		vs35,	o48,	T1

.endm

/**********************************************************************************************
* Macros for N=4 and M=2
**********************************************************************************************/

.macro COPY_4x2

	lxsspx		vs32,	o0,	A0
	lxsspx		vs33,	o4,	A0

	lxsspx		vs34,	o0,	A1
	lxsspx		vs35,	o4,	A1

	lxsspx		vs36,	o0,	A2
	lxsspx		vs37,	o4,	A2

	lxsspx		vs38,	o0,	A3
	lxsspx		vs39,	o4,	A3

	mr		T1,	BO

	stxsspx		vs32,	o0,	T1
	stxsspx		vs33,	o4,	T1

	addi		T1,	T1,	8

	stxsspx		vs34,	o0,	T1
	stxsspx		vs35,	o4,	T1

	addi		T1,	T1,	8

	stxsspx		vs36,	o0,	T1
	stxsspx		vs37,	o4,	T1

	addi		T1,	T1,	8

	stxsspx		vs38,	o0,	T1
	stxsspx		vs39,	o4,	T1

.endm

/**********************************************************************************************
* Macros for N=4 and M=1
**********************************************************************************************/

.macro COPY_4x1

	lxsspx		vs32,	o0,	A0

	lxsspx		vs33,	o0,	A1

	lxsspx		vs34,	o0,	A2

	lxsspx		vs35,	o0,	A3

	mr		T1,	BO

	stxsspx		vs32,	o0,	T1

	stxsspx		vs33,	o4,	T1

	addi		T1,	T1,	8

	stxsspx		vs34,	o0,	T1

	stxsspx		vs35,	o4,	T1

.endm

/**********************************************************************************************
* Macros for N=2 and M=16
**********************************************************************************************/

.macro COPY_2x16

	lxvw4x		vs32,	o0,	A0
	lxvw4x		vs33,	o16,	A0
	lxvw4x		vs34,	o32,	A0
	lxvw4x		vs35,	o48,	A0

	lxvw4x		vs36,	o0,	A1
	lxvw4x		vs37,	o16,	A1
	lxvw4x		vs38,	o32,	A1
	lxvw4x		vs39,	o48,	A1

	mr		T1,	BO

	stxvw4x		vs32,	o0,	T1
	stxvw4x		vs33,	o16,	T1
	stxvw4x		vs34,	o32,	T1
	stxvw4x		vs35,	o48,	T1

	addi		T1,	T1,	64

	stxvw4x		vs36,	o0,	T1
	stxvw4x		vs37,	o16,	T1
	stxvw4x		vs38,	o32,	T1
	stxvw4x		vs39,	o48,	T1

.endm

/**********************************************************************************************
* Macros for N=2 and M=8
**********************************************************************************************/

.macro COPY_2x8

	lxvw4x		vs32,	o0,	A0
	lxvw4x		vs33,	o16,	A0

	lxvw4x		vs34,	o0,	A1
	lxvw4x		vs35,	o16,	A1

	mr		T1,	BO

	stxvw4x		vs32,	o0,	T1
	stxvw4x		vs33,	o16,	T1

	stxvw4x		vs34,	o32,	T1
	stxvw4x		vs35,	o48,	T1

.endm

/**********************************************************************************************
* Macros for N=2 and M=4
**********************************************************************************************/

.macro COPY_2x4

	lxvw4x		vs32,	o0,	A0

	lxvw4x		vs33,	o0,	A1

	mr		T1,	BO

	stxvw4x		vs32,	o0,	T1

	stxvw4x		vs33,	o16,	T1

.endm

/**********************************************************************************************
* Macros for N=2 and M=2
**********************************************************************************************/

.macro COPY_2x2

	lxsspx		vs32,	o0,	A0
	lxsspx		vs33,	o4,	A0

	lxsspx		vs34,	o0,	A1
	lxsspx		vs35,	o4,	A1

	mr		T1,	BO

	stxsspx		vs32,	o0,	T1
	stxsspx		vs33,	o4,	T1

	addi		T1,	T1,	8

	stxsspx		vs34,	o0,	T1
	stxsspx		vs35,	o4,	T1

.endm

/**********************************************************************************************
* Macros for N=2 and M=1
**********************************************************************************************/

.macro COPY_2x1

	lxsspx		vs32,	o0,	A0

	lxsspx		vs33,	o0,	A1

	mr		T1,	BO

	stxsspx		vs32,	o0,	T1

	stxsspx		vs33,	o4,	T1

.endm

/**********************************************************************************************
* Macros for N=1 and M=16
**********************************************************************************************/

.macro COPY_1x16

	lxvw4x		vs32,	o0,	A0
	lxvw4x		vs33,	o16,	A0
	lxvw4x		vs34,	o32,	A0
	lxvw4x		vs35,	o48,	A0

	mr		T1,	BO

	stxvw4x		vs32,	o0,	T1
	stxvw4x		vs33,	o16,	T1
	stxvw4x		vs34,	o32,	T1
	stxvw4x		vs35,	o48,	T1

.endm

/**********************************************************************************************
* Macros for N=1 and M=8
**********************************************************************************************/

.macro COPY_1x8

	lxvw4x		vs32,	o0,	A0
	lxvw4x		vs33,	o16,	A0

	mr		T1,	BO

	stxvw4x		vs32,	o0,	T1
	stxvw4x		vs33,	o16,	T1

.endm

/**********************************************************************************************
* Macros for N=1 and M=4
**********************************************************************************************/

.macro COPY_1x4

	lxvw4x		vs32,	o0,	A0

	mr		T1,	BO

	stxvw4x		vs32,	o0,	T1

.endm

/**********************************************************************************************
* Macros for N=1 and M=2
**********************************************************************************************/

.macro COPY_1x2

	lxsspx		vs32,	o0,	A0
	lxsspx		vs33,	o4,	A0

	mr		T1,	BO

	stxsspx		vs32,	o0,	T1
	stxsspx		vs33,	o4,	T1

.endm

/**********************************************************************************************
* Macros for N=1 and M=1
**********************************************************************************************/

.macro COPY_1x1

	lxsspx		vs32,	o0,	A0

	mr		T1,	BO

	stxsspx		vs32,	o0,	T1

.endm

